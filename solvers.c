#include "variables.h"
#include "petsctime.h"

extern PetscInt   block_number, ti, tiout, NumberOfBodies, nbody;
extern PetscInt   movefsi, rotatefsi, immersed, STRONG_COUPLING;
extern PetscInt   implicit, TwoD, fish_c, sediment, wing,rheology,duplicate,turbine;
extern PetscInt   cop, regime, fish,pizza, fishcyl, MHV, BHV, beam, inv_flg, LV, moveframe,rotateframe;
extern PetscReal  max_angle, Flux_in, St_exp, FluxInSum;
extern PetscInt   les, dynamic_freq, tistart, averaging, poisson;

Cmpnts Cross(Cmpnts v1,Cmpnts v2);
Cmpnts MINUS(Cmpnts v1,Cmpnts v2);
Cmpnts rotate_y(Cmpnts v0, Cmpnts v_c, Cmpnts rot);
Cmpnts rotateT_y(Cmpnts v0, Cmpnts v_c, Cmpnts rot);
Cmpnts rotate_xyz(Cmpnts v,Cmpnts v_c, Cmpnts rot);
Cmpnts rotateT_xyz(Cmpnts v, Cmpnts rot);


PetscErrorCode Struc_Solver(UserMG *usermg,IBMNodes *ibm, 
			    FSInfo *fsi, Cstart *cstart,
			    PetscInt itr_sc,
			    PetscInt tistart, 
			    PetscBool *DoSCLoop, FE *fem) {

  PetscReal     dS_sc, dS_MIN=1e-5, dSmax,dS_rel=0.0;
  UserCtx	*user;
  PetscInt	i,j,bi,ibi, level, Add_dUndt=1,MHV_stuck=0,P_No ;
  Cmpnts	u_c,omega_c, a_c, rr;
  PetscReal     rx,ry,rz;

  double PI = 3.141592653589793;
  level = usermg->mglevels-1;
  user = usermg->mgctx[level].user;
  
  /* ==================================================================================             */
  /*     Store old values to determine SC convergence */
  /* ==================================================================================             */
  
  if (movefsi || rotatefsi || MHV || fish || cop || fish_c || rheology || pizza) {
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      for (i=0;i<6;i++){
	fsi[ibi].S_old[i] = fsi[ibi].S_new[i];
	fsi[ibi].S_ang_o[i]=fsi[ibi].S_ang_n[i];      
	if (itr_sc==1) {
	  fsi[ibi].dS[i]=0.;
	  fsi[ibi].atk=0.3;
	}
	fsi[ibi].dS_o[i]=fsi[ibi].dS[i];
	fsi[ibi].atk_o=fsi[ibi].atk;
      }
      if (itr_sc==2){
	fsi[ibi].atk_o=0.298;
      }
      
      fsi[ibi].F_x_old=fsi[ibi].F_x;
      fsi[ibi].F_y_old=fsi[ibi].F_y;
      fsi[ibi].F_z_old=fsi[ibi].F_z;
      
      fsi[ibi].M_x_old=fsi[ibi].M_x;
      fsi[ibi].M_y_old=fsi[ibi].M_y;
      fsi[ibi].M_z_old=fsi[ibi].M_z;
      
      fsi[ibi].L_o[0]=fsi[ibi].L_n[0];
      fsi[ibi].L_o[1]=fsi[ibi].L_n[1];
      fsi[ibi].L_o[2]=fsi[ibi].L_n[2];
      
    }
  }
  
  /* ==================================================================================             */
  /*     Calculating Forces! */
  /* ==================================================================================             */
  if (MHV) Add_dUndt=0;
  
  if (immersed && (!BHV) && (!beam) && (!inv_flg)) {   
    for (bi=0; bi<block_number; bi++) {      
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	
	Calc_forces_SI(&fsi[ibi],&(user[bi]),&ibm[ibi], ti, ibi, bi);
	
	if (rheology) {
	  
	  if (fsi[ibi].clone > 1.){
	    P_No=(int)((fsi[ibi].clone-(int)fsi[ibi].clone)*1001)-1;
	    fsi[P_No].M_x +=fsi[ibi].M_x;
	    fsi[P_No].M_y +=fsi[ibi].M_y;
	    fsi[P_No].M_z +=fsi[ibi].M_z;
	    fsi[P_No].F_x +=fsi[ibi].F_x;
	    fsi[P_No].F_y +=fsi[ibi].F_y;
	    fsi[P_No].F_z +=fsi[ibi].F_z;
	    fsi[ibi].M_x = fsi[P_No].M_x;
	    fsi[ibi].M_y = fsi[P_No].M_y;
	    fsi[ibi].M_z = fsi[P_No].M_z;
	    fsi[ibi].F_x = fsi[P_No].F_x;
	    fsi[ibi].F_y = fsi[P_No].F_y;
	    fsi[ibi].F_z = fsi[P_No].F_z;
	  }
	  
	  PetscPrintf(PETSC_COMM_WORLD, "Momentum on surface body %d : M_x %le M_y %le M_z %le F_x %le F_y %le F_z %le \n",ibi,fsi[ibi].M_x,fsi[ibi].M_y,fsi[ibi].M_z,fsi[ibi].F_x,fsi[ibi].F_y,fsi[ibi].F_z);
	  
	}
      }//ibi
      
      /* ==================================================================================             */
      /*       Ucat is copied here before it is changed by the flow solver */
      /*       it is needed for calculating the forces by CV method */
      /*       Ucat_o shouldn't be equal to Ucat  */
      /* ==================================================================================             */
      
      // Copy Ucat_o here!
      if (itr_sc==1) 
	VecCopy(user[bi].Ucat, user[bi].Ucat_o);
      
      /* Corrector step! start from the same solution  */
      
      if ((MHV || movefsi || rotatefsi || cop || fish || pizza || rheology) && itr_sc>1) {
	PetscPrintf(PETSC_COMM_WORLD, "Corrector Step itr # %d\n", itr_sc);
	
	VecCopy(user[bi].Ucont_o, user[bi].Ucont);
	VecCopy(user[bi].P_o, user[bi].P);
	
	DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
	
	DMGlobalToLocalBegin(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);
	DMGlobalToLocalEnd(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);
	
	Contra2Cart(&(user[bi]));
      }
      PetscBarrier(PETSC_NULL);	  
      
    } // bi
  } //immersed

  /* ==================================================================================             */
  /*     Find The new Position & Move the BODY */
  /* ==================================================================================             */
  if (rheology){
    bi=0;
    
    level = usermg->mglevels-1;
    user = usermg->mgctx[level].user;
    
    for (ibi=0; ibi< NumberOfBodies; ibi++){
      
      if (fsi[ibi].clone < 1.0){
	calc_quarternion(&fsi[ibi],user[bi].dt,ibi);
	Calc_FSI_vel_particle(&fsi[ibi],&ibm[ibi],user[bi].dt);
      }

      else{
	
	P_No=(int)((fsi[ibi].clone-(int)fsi[ibi].clone)*1001)-1;
	
	for (i=0;i<3;i++){
	  
	  fsi[ibi].S_ang_n[2*i+1]= fsi[P_No].S_ang_n[2*i+1];
	  fsi[ibi].alpha[i]=fsi[P_No].alpha[i];
	  fsi[ibi].acc[i]=fsi[P_No].acc[i];
	  fsi[ibi].S_new[2*i+1]=fsi[P_No].S_new[2*i+1];
	  fsi[ibi].L_n[i]=fsi[P_No].L_n[i];
	  
	}
	
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    fsi[ibi].R[i][j]=fsi[P_No].R[i][j];
	  }
	}
      }
    }   
    
    if (NumberOfBodies>1)   CollisionDetectionOfParticles(&(user[bi]),ibm,fsi);
    
    for (ibi=0; ibi< NumberOfBodies; ibi++){
      
      rotate_quarternion(&ibm[ibi],&fsi[ibi],user[bi].dt);
      
      calc_ibm_normal(&ibm[ibi]);
      
      omega_c.x=fsi[ibi].S_ang_n[1];
      omega_c.y=fsi[ibi].S_ang_n[3];
      omega_c.z=fsi[ibi].S_ang_n[5];
      
      // axis of frame rotation
      a_c.x=fsi[ibi].x_c;
      a_c.y=fsi[ibi].y_c;
      a_c.z=fsi[ibi].z_c;
      
      for (i=0; i<ibm[ibi].n_v; i++) {
	
	rx = ibm[ibi].x_bp[i]-a_c.x;
	ry = ibm[ibi].y_bp[i]-a_c.y;
	rz = ibm[ibi].z_bp[i]-a_c.z;
	ibm[ibi].u[i].x =   (rz*omega_c.y-omega_c.z*ry)+fsi[ibi].S_new[1]  ;
	ibm[ibi].u[i].y =   (rx*omega_c.z-omega_c.x*rz)+fsi[ibi].S_new[3]  ;
	ibm[ibi].u[i].z =   (ry*omega_c.x-omega_c.y*rx)+fsi[ibi].S_new[5]  ;
	
      }
     
      ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
      
      if (ti == (ti/tiout)*tiout)  ibm_surface_VTKOut(&ibm[ibi],ibi,ti);
    }
    
  }

  if (movefsi && moveframe==0){
    for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--) {
      user = usermg->mgctx[level].user;
      for (bi=0; bi<block_number; bi++) {
	if (immersed) {	 
	  
	  for (ibi=0;ibi<NumberOfBodies;ibi++) {
	   
	    Calc_FSI_pos_intg(&fsi[ibi], &ibm[ibi], user[bi].dt) ;
	   	  
	  }	
	  
	  CollisionDetectionOfCylinders(fsi,ibm,NumberOfBodies);
	  
	  for (ibi=0;ibi<NumberOfBodies;ibi++) {
	    Elmt_Move_FSI_TRANS(&fsi[ibi], &ibm[ibi],ibi);
	  }
	  PetscBarrier(PETSC_NULL);
	  
	  VecSet(user[bi].Nvert,0.);
	  VecSet(user[bi].lNvert,0.);
	  for (ibi=0;ibi<NumberOfBodies;ibi++) {
	    PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
	    ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
	    PetscBarrier(PETSC_NULL);	  
	  }
	}
      }
    }  
  }
  else if (movefsi && moveframe==1){
    for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--) {
      user = usermg->mgctx[level].user;
      for (bi=0; bi<block_number; bi++) {
       
	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	 
	  Forced_Motion(&fsi[ibi], user->dt);	  
	  
	  // u_c is frame translational speed
	  u_c.x=fsi[ibi].S_new[1];
	  u_c.y=fsi[ibi].S_new[3];
	  u_c.z=fsi[ibi].S_new[5];
	  
	  // omega_c is frame angular speed
	  omega_c.x=fsi[ibi].S_ang_n[1];
	  omega_c.y=fsi[ibi].S_ang_n[3];
	  omega_c.z=fsi[ibi].S_ang_n[5];
	  
	  // axis of frame rotation
	  a_c.x=fsi[ibi].x_c;
	  a_c.y=fsi[ibi].y_c;
	  a_c.z=fsi[ibi].z_c;
	  
	  CalcFrameVelocity(&user[bi],u_c,omega_c,a_c);
	  
	  if (immersed) {
	    
	    for (i=0; i<ibm[ibi].n_v; i++) {
	      ibm[ibi].u[i].x = fsi[ibi].S_new[1];
	      ibm[ibi].u[i].y = fsi[ibi].S_new[3];
	      ibm[ibi].u[i].z = fsi[ibi].S_new[5];
	    }
	    
	  }
	}
      } //ibi
    } //level
    PetscPrintf(PETSC_COMM_WORLD, "MV_FRAME!!!! %le %le %le\n",fsi[0].S_new[3],fsi[0].S_new[5],St_exp);
  }
  
  if (rotatefsi && rotateframe==0) {
    for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--) {
      user = usermg->mgctx[level].user;
      for (bi=0; bi<block_number; bi++) {
	if (immersed) {
	  for (ibi=0;ibi<NumberOfBodies;ibi++) {
	    
	    if (ibi==0) {
	      fsi[ibi].S_ang_n[5]=St_exp;//this is for steady rotation
	      fsi[ibi].S_ang_n[4]=fsi[ibi].S_ang_n[5]*user[bi].dt*ti;
	      
	      rotate_ibm(&ibm[ibi],&fsi[ibi]);
	      calc_ibm_normal(&ibm[ibi]);
	      
	      // omega_c is frame angular speed
	      omega_c.x=fsi[ibi].S_ang_n[1];
	      omega_c.y=fsi[ibi].S_ang_n[3];
	      omega_c.z=fsi[ibi].S_ang_n[5];
	      
	      // axis of frame rotation
	      a_c.x=0.;
	      a_c.y=0.;
	      a_c.z=0.;
	      
	      for (i=0; i<ibm[ibi].n_v; i++) {
		ibm[ibi].u[i].x = fsi[ibi].S_new[1];
		ibm[ibi].u[i].y = fsi[ibi].S_new[3];
		ibm[ibi].u[i].z = fsi[ibi].S_new[5];
		
		rx = ibm[ibi].x_bp[i]-a_c.x;
		ry = ibm[ibi].y_bp[i]-a_c.y;
		rz = ibm[ibi].z_bp[i]-a_c.z;      
		ibm[ibi].u[i].x -=   ry*omega_c.z-omega_c.y*rz  ;
		ibm[ibi].u[i].y += ( rx*omega_c.z-omega_c.x*rz );
		ibm[ibi].u[i].z -=   rx*omega_c.y-omega_c.x*ry  ;      	    
	      }
	      if (ti == (ti/tiout)*tiout) 
		ibm_surface_out(&ibm[ibi],ti,ibi);
	    }
	    
	    PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
	    ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);    
	  }
	}
      }
    }      
  } else if (rotatefsi && rotateframe==1){
    for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--) {
      user = usermg->mgctx[level].user;
      for (bi=0; bi<block_number; bi++) {

	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	  	 
	  if (ti==tistart) {
	    fsi[ibi].S_ang_n[5]=St_exp;//this is for steady rotation
	    // u_c is frame translational speed
	    
	    u_c.x=fsi[ibi].S_new[1];
	    u_c.y=fsi[ibi].S_new[3];
	    u_c.z=fsi[ibi].S_new[5];
	    
	    // omega_c is frame angular speed
	    omega_c.x=fsi[ibi].S_ang_n[1];
	    omega_c.y=fsi[ibi].S_ang_n[3];
	    omega_c.z=fsi[ibi].S_ang_n[5];
	    
	    // axis of frame rotation
	    a_c.x=0.;
	    a_c.y=0.;
	    a_c.z=0.;
	    
	    if (ibi==0) {	    
	      CalcFrameVelocity(&user[bi],u_c,omega_c,a_c);
	      GridDivergence(&user[bi]);
	    }
	    
	    if (immersed) {
	      
	      for (i=0; i<ibm[ibi].n_v; i++) {
		ibm[ibi].u[i].x = fsi[ibi].S_new[1];
		ibm[ibi].u[i].y = fsi[ibi].S_new[3];
		ibm[ibi].u[i].z = fsi[ibi].S_new[5];
		
		rx = ibm[ibi].x_bp[i]-a_c.x;
		ry = ibm[ibi].y_bp[i]-a_c.y;
		rz = ibm[ibi].z_bp[i]-a_c.z;      
		ibm[ibi].u[i].x -=   ry*omega_c.z-omega_c.y*rz  ;
		ibm[ibi].u[i].y += ( rx*omega_c.z-omega_c.x*rz );
		ibm[ibi].u[i].z -=   rx*omega_c.y-omega_c.x*ry  ;      	    
	      }
	      ibm_surface_out(&ibm[ibi],0,ibi);
	      
	    }
	  } // if ti==tistart
	}
      } //ibi
    } //level
    PetscPrintf(PETSC_COMM_WORLD, "MV_FRAME ROTATE!!!! %le %le %le %le %le\n",fsi[0].S_new[3],fsi[0].S_new[5],fsi[0].S_ang_n[0],fsi[0].S_ang_n[5],St_exp);
  }

  /* ==================================================================================             */
  /*   For MHV & LV */
  /* ==================================================================================             */
  if (MHV && ti>-1) {
    for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--) {
      user = usermg->mgctx[level].user;
      PetscReal dir=1.;
      PetscInt  ibi_st=1;
      PetscReal dw=0;
      for (bi=0; bi<block_number; bi++) {
  	if (immersed) {
	  
	  for (ibi=ibi_st;ibi<NumberOfBodies;ibi++) {
	    
	    Calc_FSI_Ang_intg(&fsi[ibi], &ibm[ibi], user[bi].dt,itr_sc,ibi,&user[bi],0.0) ;
	    //
	    if (fabs(fsi[ibi].S_ang_n[0])== fabs(max_angle)) {
	      fsi[ibi].S_ang_n[0]= max_angle;
	      fsi[ibi].S_ang_n[1]= 0.;
	      MHV_stuck=1;
	    }
	    
	    if (fsi[ibi].S_ang_n[0]== 0.0) {
	      fsi[ibi].S_ang_n[0]= 0.0;
	      fsi[ibi].S_ang_n[1]= 0.0;
	    
	    }
	  }
	  //
	  for (ibi=ibi_st;ibi<NumberOfBodies;ibi++) {
	    Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi], user[bi].dt, ibi);
	    if (ti == (ti/tiout)*tiout)
	      ibm_surface_VTKOut(&ibm[ibi],ibi,ti);
	  }
	  //
	  PetscBarrier(PETSC_NULL);
	  for (ibi=0;ibi<NumberOfBodies;ibi++) {
	    PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
	    ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
	    PetscBarrier(PETSC_NULL);
	  }
  	}
      }
    }
  }
  
  /* ==================================================================================             */
  /*   For Copepod  */
  /* ==================================================================================             */
  
  if (cop) {
    level = usermg->mglevels-1;
    user = usermg->mgctx[level].user;
    for (bi=0; bi<block_number; bi++) {
      PetscPrintf(PETSC_COMM_WORLD, "cop swim!!\n");
      cop_swim(ibm, ti*user[bi].dt, user[bi].dt);
      PetscPrintf(PETSC_COMM_WORLD, "cop swim!!\n");
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
	ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
	PetscBarrier(PETSC_NULL);
	
      }//ibi
      
    }
  } 
  
  /* ==================================================================================             */
  /*   For fish  */
  /* ==================================================================================             */
  if (fish_c){
    for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--) {
      user = usermg->mgctx[level].user;
      
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	
	//this is for steady rotation
	form_cstart(cstart, &ibm[ibi], &fsi[ibi], user[0].dt, ibi,rotateframe);		  

	if (rotateframe) {
	  // u_c is frame translational speed
	  u_c.x=fsi[ibi].S_new[1];
	  u_c.y=fsi[ibi].S_new[3];
	  u_c.z=fsi[ibi].S_new[5];
	  
	  // omega_c is frame angular speed
	  omega_c.x=-fsi[ibi].S_ang_n[1];
	  omega_c.y=-fsi[ibi].S_ang_n[3];
	  omega_c.z=-fsi[ibi].S_ang_n[5];
	  
	  // axis of frame rotation
	  a_c.x=0.2;
	  a_c.y=0.;
	  a_c.z=-0.6;
	  
	  if (ibi==0)
	    for (bi=0; bi<block_number; bi++) 
	      CalcFrameVelocity(&user[bi],u_c,omega_c,a_c);
	} // if rotateframe
	
	if (immersed) {
	  if (rotateframe) {
	    
	    for (i=0; i<ibm[ibi].n_v; i++) {
	      ibm[ibi].u[i].x += fsi[ibi].S_new[1];
	      ibm[ibi].u[i].y += fsi[ibi].S_new[3];
	      ibm[ibi].u[i].z += fsi[ibi].S_new[5];
	      
	      rx = ibm[ibi].x_bp[i]-a_c.x;
	      ry = ibm[ibi].y_bp[i]-a_c.y;
	      rz = ibm[ibi].z_bp[i]-a_c.z;
	      ibm[ibi].u[i].x -=   ry*omega_c.z-omega_c.y*rz  ;
	      ibm[ibi].u[i].y += ( rx*omega_c.z-omega_c.x*rz );
	      ibm[ibi].u[i].z -=   rx*omega_c.y-omega_c.x*ry  ;
	    }
	  }//if rotate frame	  	  
	  
	  ibm_surface_VTKOut(&ibm[ibi],ibi,ti);
	  
	  for (bi=block_number-1; bi<block_number; bi++) {
	   
	    PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA %d %d\n", ibi, bi);
	    PetscBarrier(PETSC_NULL);
	    ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
	  }//bi
	}// immersed
	
      } //ibi
    } //level
    
    if (rotateframe) 
      PetscPrintf(PETSC_COMM_WORLD, "MV_FRAME ROTATE!!!! %le %le %le %le %le %le\n",u_c.x, u_c.y,u_c.z,omega_c.x, omega_c.y, omega_c.z);
  }
  
  if ((fish==1)  && (moveframe==0) && (St_exp>0.01)) {
    level = usermg->mglevels-1;
    user = usermg->mgctx[level].user;
    bi=0;
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
    
      if (!fishcyl || ibi<1){
	PetscPrintf(PETSC_COMM_WORLD, "fish swim! \n");
	fish_swim(&ibm[ibi],(ti)*user[bi].dt, user[bi].dt);
	if (ti == (ti/tiout)*tiout)
	  //amir 
	  ibm_surface_VTKOut(&ibm[ibi],ibi,ti);
        
      }
   
    }
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      for (bi=0; bi<block_number; bi++) { 
	PetscBarrier(PETSC_NULL);
	if (bi==ibi+1 && NumberOfBodies>1) {
	  VecSet(user[bi].Nvert,0.);
	  VecSet(user[bi].lNvert,0.);
	  PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA bi %d ibi %d \n", bi, ibi);
	  ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
	} else if (NumberOfBodies==1) {
	  PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
	  ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
	}
	
      }
    }
  } else if ( (fish || pizza) && moveframe ) {
    level = usermg->mglevels-1;
    user = usermg->mgctx[level].user;
    
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      if (block_number>1) bi=1;
      else bi=0;
      
      Calc_FSI_pos_intg(&fsi[ibi], &ibm[ibi], user[bi].dt) ;
      
      if (ibi==0) {
	PetscPrintf(PETSC_COMM_WORLD, "fish swim! \n");
	if (pizza)
	  pizza_box_swim(&ibm[ibi],(ti)*user[bi].dt, user[bi].dt);
	else
	  fish_swim(&ibm[ibi],(ti)*user[bi].dt, user[bi].dt);
      }
      
      for (i=0; i<ibm[ibi].n_v; i++) {
	ibm[ibi].u[i].x += fsi[ibi].S_new[1];
	ibm[ibi].u[i].y += fsi[ibi].S_new[3];
	ibm[ibi].u[i].z += fsi[ibi].S_new[5];
      }
      if (ti == (ti/tiout)*tiout) 
	ibm_surface_VTKOut(&ibm[ibi],ibi,ti);
      
      PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
      ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
      
      for (bi=0; bi<block_number; bi++) {
	u_c.x=fsi[ibi].S_new[1];
	u_c.y=fsi[ibi].S_new[3];
	u_c.z=fsi[ibi].S_new[5];
	
	// omega_c is frame angular speed
	omega_c.x=fsi[ibi].S_ang_n[1];
	omega_c.y=fsi[ibi].S_ang_n[3];
	omega_c.z=fsi[ibi].S_ang_n[5];
	
	// axis of frame rotation
	a_c.x=0.;
	a_c.y=0.;
	a_c.z=0.;
	
	if (ibi==0) {	    
	  CalcFrameVelocity(&user[bi],u_c,omega_c,a_c);
	}
      }
      
    }
  }    
  
  /* ==================================================================================             */
  /*  Sediment  */
  /* ==================================================================================             */
  if (sediment) {
    level = usermg->mglevels-1;
    user = usermg->mgctx[level].user;
    for (bi=0; bi<block_number; bi++) {
      if (immersed) {
	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	  Calc_FSI_pos_intg(&fsi[ibi], &ibm[ibi], user[bi].dt) ;
	  for (i=0; i<ibm[ibi].n_v; i++) {
	    ibm[ibi].u[i].x = fsi[ibi].S_new[1];
	    ibm[ibi].u[i].y = fsi[ibi].S_new[3];
	    ibm[ibi].u[i].z = fsi[ibi].S_new[5];
	  }
	  if (ti == (ti/tiout)*tiout) ibm_surface_VTKOut(&ibm[ibi],ibi,ti);
	  PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
	  ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi); 
	  
	}
      }
    }// bi
  }
  
  /* ==================================================================================             */
  /*  FEM */
  /* ==================================================================================             */
  PetscInt  rank;
  
  if (BHV || inv_flg) {
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    level = usermg->mglevels - 1;
    user = usermg->mgctx[level].user;
    
    bi = 0; 
    for (i=0; i<nbody; i++) {
      // Find Fluid forces for FEM
      PetscPrintf(PETSC_COMM_WORLD, "Begin Force FEM time %d ibi %d!\n", ti, i);
      Projecting(&user[bi], &ibm[i], i);
      PetscBarrier(PETSC_NULL);
      Finalizing_Projecting(&user[bi], &ibm[i]);
      if (ti==(ti/tiout)*tiout) Projecting_VTKOut(&ibm[i], i, ti);
      PetscPrintf(PETSC_COMM_WORLD, "End Force FEM time %d ibi %d!\n", ti, i);
      fem[i].dt = user[bi].dt;
    }
    PetscBarrier(PETSC_NULL);
    
    if(!rank){
      fem_solve(ibm, fem, itr_sc);
    }
    
    PetscBarrier(PETSC_NULL);
    
    for (i=0; i<nbody; i++) { 
      // scatter new ibm to all processes
      BroadcastIBM(&ibm[i]);
      PetscBarrier(PETSC_NULL);
      calc_ibm_normal(&ibm[i]); //clac normals and areas
      calc_ibm_velocity(&ibm[i], user[bi].dt);
    }
    
    PetscBarrier(PETSC_NULL);
    
    for (bi=0; bi<block_number; bi++) {
      VecSet(user[bi].Nvert, 0.);
      VecSet(user[bi].lNvert, 0.);
      
      for (ibi=0; ibi<NumberOfBodies; ibi++) {
    	if (ibi==6) {
    	  ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
    	} else if(ibi<6) {
    	  ibm_search_thin(&(user[bi]), &ibm[ibi], ibi); // search
    	}
      }
    } //bi
  } // if fem
  
  if (beam) {
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    level = usermg->mglevels-1;
    user = usermg->mgctx[level].user;
    
    //for (bi=0; bi<block_number; bi++) {
    bi = 0; //leaflets are in block number one
    for (i=0; i<nbody; i++) {
      // Find Fluid forces for FEM
      PetscPrintf(PETSC_COMM_WORLD, "Begin Force FEM time %d ibi %d!\n", ti, i);
      Projecting(&user[bi], &ibm[i], i);
      PetscBarrier(PETSC_NULL);
      Finalizing_Projecting(&user[bi], &ibm[i]);
      Projecting_VTKOut(&ibm[i], i, ti);
      PetscPrintf(PETSC_COMM_WORLD, "End Force FEM time %d ibi %d!\n", ti, i);
      fem[i].dt = user[bi].dt;
    }
    PetscBarrier(PETSC_NULL);
    
    if (!rank) {fem_solve(ibm, fem, itr_sc);}
    //  if(!rank){rigid_move(ibm,fem,itr_sc);}
    
    PetscBarrier(PETSC_NULL);
    
    for (i=0; i<nbody; i++) {
      // scatter new ibm to all processes
      BroadcastIBM(&ibm[i]);
      PetscBarrier(PETSC_NULL);
      calc_ibm_normal(&ibm[i]); //clac normals and areas
      calc_ibm_velocity(&ibm[i], user[bi].dt);
    }
    
    PetscBarrier(PETSC_NULL);
    
    for (bi=0; bi<block_number; bi++) {
      VecSet(user[bi].Nvert, 0.);
      VecSet(user[bi].lNvert, 0.);
      
      for (ibi=0; ibi<NumberOfBodies; ibi++) {
  	ibm_search_thin(&(user[bi]), &ibm[ibi], ibi); // search
      }
    }//bi
    
    if (!rank){
      if (ti==(ti/tiout)*tiout){
  	ibm_surface_VTKOut(&ibm[0],0,ti);
      }
    }
    
  } // if beam
  
  if (MHV){
    if (!rank){
      if (ti==(ti/tiout)*tiout){
  	for (ibi=0;ibi<NumberOfBodies;ibi++) {
  	  ibm_surface_VTKOut(&ibm[ibi],ibi,ti);
  	}
      }
    }
  }
  /* ==================================================================================             */
  /*  Wing */
  /* ==================================================================================             */
  if (wing) {
    bi=0;
    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      wing_motion(&ibm[ibi], &fsi[ibi], ti, &user[bi].dt);
      if (block_number>1)
	for (bi=1; bi<block_number; bi++)  user[bi].dt=user[0].dt;
      
      if (rotateframe) {
	// u_c is frame translational speed
	u_c.x=fsi[ibi].S_new[1];
	u_c.y=fsi[ibi].S_new[3];
	u_c.z=fsi[ibi].S_new[5];
	
	// omega_c is frame angular speed
	omega_c.x=fsi[ibi].S_ang_n[1];
	omega_c.y=fsi[ibi].S_ang_n[3];
	omega_c.z=fsi[ibi].S_ang_n[5];
	
	// axis of frame rotation
	a_c.x=fsi[ibi].x_c;
	a_c.y=fsi[ibi].y_c;
	a_c.z=fsi[ibi].z_c;	    
	
	// angle of rotation
	Cmpnts alfa_c;
	alfa_c.x=fsi[ibi].S_ang_n[0];
	alfa_c.y=fsi[ibi].S_ang_n[2];
	alfa_c.z=fsi[ibi].S_ang_n[4];
	
	omega_c = rotateT_y(omega_c,a_c, alfa_c);
	PetscPrintf(PETSC_COMM_WORLD, "frame omega %le %le %le wing %le %le %le\n", omega_c.x, omega_c.y, omega_c.z,fsi[ibi].S_ang_n[1],fsi[ibi].S_ang_n[3],fsi[ibi].S_ang_n[5] );
	
	if (ibi==0)
	  for (bi=0; bi<block_number; bi++) 
	    CalcFrameVelocity(&user[bi],u_c,omega_c,a_c);
	
	// wing velocity
	for (i=0; i<ibm[ibi].n_v; i++) {	  
	  rr.x = ibm[ibi].x_bp0[i];
	  rr.y = ibm[ibi].y_bp0[i];
	  rr.z = ibm[ibi].z_bp0[i];
	  rr   = MINUS(rr, a_c);
	  
	  ibm[ibi].u[i] = Cross(omega_c, rr);
	  
	}
	
	// write the surface file to check vel
	ibm_surface_out(&ibm[ibi],ti,ibi);
	
	// rotate the forces from non-inertial to inertial
	for (bi=0; bi<block_number; bi++) 
	  rotate_force(&fsi[ibi], ibi, bi);
	
      } else {
	rotate_ibm(&ibm[ibi],&fsi[ibi]);
	calc_ibm_normal(&ibm[ibi]);	
	if (ti>0) calc_ibm_velocity(&ibm[ibi], user[bi].dt);
	ibm_surface_out(&ibm[ibi],ti,ibi);
	
	PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
	ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
      }
    }// for ibi
   
  } // if wing
  
  /* ==================================================================================             */
  /*   Convergence of the SC Loop */
  /* ==================================================================================             */
  
  *DoSCLoop = PETSC_FALSE;
  dSmax=1e-10;
  for (ibi=0;ibi<NumberOfBodies;ibi++) {
    
    if ((movefsi || fish || pizza || cop) && STRONG_COUPLING) {
      for (i=0;i<6;i++) {
	dS_sc = fabs(fsi[ibi].S_new[i]-fsi[ibi].S_old[i]);
	if (dS_sc > dS_MIN) *DoSCLoop = PETSC_TRUE;
	if (dS_sc > dSmax) dSmax=dS_sc;
      }      
    }
    
    if ((rotatefsi || MHV) && STRONG_COUPLING) {
      dS_sc = fabs(fsi[ibi].S_ang_n[0]-fsi[ibi].S_ang_o[0]);
      if (dS_sc > dS_MIN) *DoSCLoop = PETSC_TRUE;
      dS_sc = fabs(fsi[ibi].S_ang_n[1]-fsi[ibi].S_ang_o[1]);
      if(fabs(fsi[ibi].S_ang_n[1]+fsi[ibi].S_ang_o[1])>2.) 
	dS_sc /= 0.5*fabs(fsi[ibi].S_ang_n[1]+fsi[ibi].S_ang_o[1]);
      
      if (dS_sc > dS_MIN) *DoSCLoop = PETSC_TRUE;
      if (dS_sc > dSmax) dSmax=dS_sc;
      
    }
    if (rheology && STRONG_COUPLING) {
      dS_MIN=0.005;
      dS_sc = fabs(fsi[ibi].S_ang_n[1]-fsi[ibi].S_ang_o[1]);
      dS_rel= fabs(fsi[ibi].S_ang_n[1]/fsi[ibi].S_ang_o[1]);
      
      if (dS_rel<0.05)  goto nextp;  
      
      if (dS_sc > dS_MIN) *DoSCLoop = PETSC_TRUE;   
    nextp: PetscPrintf(PETSC_COMM_WORLD, "S-C Convergence dS_rel %le dS_sc %le \n",dS_rel,dS_sc); 
    }
    if (STRONG_COUPLING && (BHV || inv_flg)) {
      
      if (ibi<3) {
	if (!rank) {
	  VecNorm(fem[ibi].dxsc, NORM_INFINITY, &dS_sc);
	  MPI_Bcast(&dS_sc, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
	} else {
	  MPI_Bcast(&dS_sc, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
	}
	PetscBarrier(PETSC_NULL);
	if (dS_sc>dS_MIN) *DoSCLoop = PETSC_TRUE;
	if (dS_sc>dSmax) dSmax = dS_sc;
	PetscPrintf(PETSC_COMM_WORLD, "S-C Convergence %d %le %le\n", itr_sc, dSmax, dS_MIN);
      }
    }
    if ((movefsi || rotatefsi || MHV || BHV || inv_flg || cop || fish || pizza || rheology) && STRONG_COUPLING && itr_sc<2) 
      *DoSCLoop = PETSC_TRUE;  
    
  } // ibi  
  
  if (itr_sc>30) *DoSCLoop = PETSC_FALSE;
  if (MHV_stuck && itr_sc==1) *DoSCLoop = PETSC_FALSE;
  
  if (cop)
    PetscPrintf(PETSC_COMM_WORLD, "S-C Convergence %d %le %le %le\n", itr_sc, dSmax,fsi[0].S_new[5],fsi[0].S_old[5]);
  else
    PetscPrintf(PETSC_COMM_WORLD, "S-C Convergence %d %le %le %le\n", itr_sc, dSmax,fsi[0].S_ang_n[0],fsi[0].S_ang_o[0]);
  
  for (ibi=0;ibi<NumberOfBodies;ibi++) {
    
    if ((movefsi || rotatefsi || regime || MHV) && !(*DoSCLoop)) {
      FSI_DATA_Output(&fsi[ibi], ibi);
      if (MHV==5 && ibi==0) MomentumJet(&user[0], 189);
    }
    
  }
  /* ==================================================================================             */
  
  return(0);
}
/* ==================================================================================             */

PetscErrorCode Struc_predictor(UserMG *usermg,IBMNodes *ibm, 
			       FSInfo *fsi, PetscInt itr_sc,
			       PetscInt tistart, 
			       PetscBool *DoSCLoop)
{
  UserCtx	*user;
  PetscInt	bi,ibi, level,MHV_stuck=0 ;

  level = usermg->mglevels-1;
  user = usermg->mgctx[level].user;

  if (MHV && ti>-1) {
    for (level = usermg->mglevels-1; level>=usermg->mglevels-1; level--) {
      user = usermg->mgctx[level].user;
      PetscReal dir=1.;
    

      for (bi=0; bi<block_number; bi++) {
	if (immersed) {
	  // for leaflets ibi = 1 & 2
	for (ibi=1;ibi<NumberOfBodies;ibi++) {
	  dir = -1*dir;
	  fsi[ibi].S_ang_n[1]=2.*fsi[ibi].S_ang_r[1]-fsi[ibi].S_ang_rm1[1];
	  fsi[ibi].S_ang_n[0]=fsi[ibi].S_ang_r[0]+0.5*(fsi[ibi].S_ang_n[1]+fsi[ibi].S_ang_r[1])*user->dt;

	  if ((dir*fsi[ibi].S_ang_n[0])> -max_angle) {
	    fsi[ibi].S_ang_n[0]= -dir*max_angle;
	    fsi[ibi].S_ang_n[1]= 0.;
	    MHV_stuck=1;
	  }
	  if ((dir*fsi[ibi].S_ang_n[0])< 0.0) {
	    fsi[ibi].S_ang_n[0]= dir*0.0;
	    fsi[ibi].S_ang_n[1]= 0.;
	    MHV_stuck=1;
	  }


	}
	//PetscBarrier(PETSC_NULL);
	for (ibi=1;ibi<NumberOfBodies;ibi++) {
	  Elmt_Move_FSI_ROT(&fsi[ibi], &ibm[ibi], user[bi].dt, ibi);
	}
	PetscBarrier(PETSC_NULL);
	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	  PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA\n");
	  ibm_search_advanced(&(user[bi]), &ibm[ibi], ibi);
	  PetscBarrier(PETSC_NULL);	  
	}
	}
      }
    }
  }    
  return(0);
}
/* ==================================================================================             */


/* ==================================================================================             */
/*     Flow Solver! */
/* ==================================================================================             */
PetscErrorCode Flow_Solver(UserMG *usermg,IBMNodes *ibm, 
			   FSInfo *fsi)
{

  UserCtx	*user;
  PetscInt	bi, ibi, level;
  PetscReal     tm_s, tm_e, tp_e, tpr_e;

  level = usermg->mglevels-1;
  user = usermg->mgctx[level].user;


/* ==================================================================================             */
/*   Turbulence Models! */
/* ==================================================================================             */

  if(rans){ 	/* rans == 1: Wilcox Low Re; 2: Wilcox High Re, 3: SST Menter */
     for (bi=0; bi<block_number; bi++) { 
       K_Omega_Set_Constant(&user[bi]);
	
       if (ti==tistart)
	 if (immersed) 
	   for (ibi=0;ibi<NumberOfBodies;ibi++) {
	     ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
	   }

       if(ti==0) {
	 PetscPrintf(PETSC_COMM_WORLD, "\nInitializing K-omega ... \n\n");
	 K_Omega_IC(&user[bi]);
       } 
       else if(ti>0 && block_number<2)  {
	 //Solve_K_Omega_DMMG(user);
	 //Solve_K_Omega(&user[bi]);
       }       
     } // for bi

     /* Set from options for DMMG solver:
     IMPORTANT: Always use this option -dmmg_jacobian_mf_fd or it will be
     very slow!!!
     -snes_mf
     -snes_type tr
     -snes_max_it 10
     -snes_monitor
     -snes_ksp_ew
     -snes_ksp_ew_version 3
     */
    /*  if(ti>0 && block_number>1) { */
/* 	 Solve_K_Omega_DMMG(user); */
/*      } */

     for (bi=0; bi<block_number; bi++) {// set constants and calulate nu_t
       K_Omega_Set_Constant(&user[bi]);
       PetscReal kmax,kmin,emax,emin;
       VecStrideMax(user[bi].K_Omega, 0, PETSC_NULL, &kmax);
       VecStrideMin(user[bi].K_Omega, 0, PETSC_NULL, &kmin);
       VecStrideMax(user[bi].K_Omega, 1, PETSC_NULL, &emax);
       VecStrideMin(user[bi].K_Omega, 1, PETSC_NULL, &emin);
       PetscPrintf(PETSC_COMM_WORLD, "Block %d K max/min %le %le and Omega max/min %le %le\n",bi,kmax,kmin, emax,emin);
     }
  } // if rans


  if(les){
    for (bi=0; bi<block_number; bi++) {
      PetscPrintf(PETSC_COMM_WORLD, "LES\n");    
      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
      Contra2Cart(&user[bi]);
      if(ti%dynamic_freq==0 || ti==tistart) 
	Compute_Smagorinsky_Constant_1(&user[bi], user[bi].lUcont, user[bi].lUcat);
      Compute_eddy_viscosity_LES(&user[bi]);
    }
  }

/* ==================================================================================             */
/*   Momentum Solver! */
/* ==================================================================================             */

  PetscTime(&tm_s);
  if (implicit==1) {
    ImplicitMomentumSolver(user, ibm, fsi);
  }else if (implicit==2) {
    ImplicitMomentumSolver1(user, ibm, fsi);
  }else if (implicit==3) {
    ImpRK(user, ibm, fsi);
  }else if (implicit==4) {
    extern PetscErrorCode Implicit_MatrixFree(UserCtx *user, IBMNodes *ibm,
					      FSInfo *fsi);
    Implicit_MatrixFree(user, ibm, fsi);
  }else if (implicit==5) {
     Implicit_SNES(user, ibm, fsi);
  }else if (implicit==7) {
    ImpRK(user, ibm, fsi);
    // Implicit_DMMG(usermg);
  }else {
    RungeKutta(user, ibm, fsi);
  }
  PetscTime(&tm_e);

/* ================================================================================== */
/*    Poisson Solver! */
/* ================================================================================== */


  PetscBarrier(PETSC_NULL);
  PetscInt iPress;
  for (iPress=0; iPress<1; iPress++) {
   
    if (poisson==0) PoissonSolver_MG(usermg);
    else if(poisson==1){ // if need symmetric hypre run with -hypre_symmetric
      /* if (ti==tistart || ti==0) Initialize_Hypre_Var(block_number); */
      /* for (bi=0; bi<block_number; bi++) */
      /* 	PoissonSolver_Hypre(&user[bi], ibm, user[bi].ibm_intp); */
    }
  }
  PetscTime(&tp_e);
/* ================================================================================== */
/*    Velocity Correction! */
/* ================================================================================== */

  for (bi=0; bi<block_number; bi++) {
    //      VecAXPY(user[bi].P, 1., user[bi].Phi);
    UpdatePressure(&user[bi]);
    //Projection(&(user[bi]), user[bi].Ucont);
    Projection(&(user[bi]));
    //UpdateBC(&user[bi]);
    DMGlobalToLocalBegin(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);
    DMGlobalToLocalEnd(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP);
    
  }
  PetscTime(&tpr_e);
  
  
  
  for (bi=0; bi<block_number; bi++) {
    
/* ==================================================================================             *
   Checking Convergence!
   /* ==================================================================================             */
    
    Divergence(&(user[bi]));
    
    //PetscPrintf(PETSC_COMM_WORLD,"divergence done");
    
    /* ==================================================================================             *
       BC!
/* ==================================================================================             */
    /* Update Ucont at the IB nodes, which will be stored as Ucont_o!!
       This is important if an IB node becomes a fluid node at the
       next time step */
    if (immersed) {
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	PetscPrintf(PETSC_COMM_WORLD,"INTP Start");
	ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi,1);
      	PetscPrintf(PETSC_COMM_WORLD,"INTP End");
	}
    }

    /* ==================================================================================             */
    /*     OUTPUT Values! */
/* ==================================================================================             */
    
    if (ti == (ti/tiout) * tiout)// || ti<10)
      Ucont_P_Binary_Output(&(user[bi]));
    
    if (les) 
      KE_Output(&user[bi]);
    
    if(averaging) {	// seokkoo
      extern PetscErrorCode Do_averaging(UserCtx *user);
      Do_averaging(&user[bi]);
      
      /*      
	      VecAXPY(user[0].P_sum, 1.0, user[0].P);
	      PetscScalar *p2sum, *p;
	      VecGetLocalSize(user[0].P_square_sum, &N);
	      VecGetArray(user[0].P_square_sum, &p2sum);
	      VecGetArray(user[0].P, &p);
	      for(v=0; v<N; v++) p2sum[v] += p[v] * p[v];
	      VecRestoreArray(user[0].P_square_sum, &p2sum);
	      VecRestoreArray(user[0].P, &p);
      */
    }
    
  } // for bi
  
  PetscInt profiling=0, rank, size;
  PetscOptionsGetInt(PETSC_NULL, "-profiling", &profiling, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (profiling && !rank){
    FILE *f;
    char filen[80];
    //sprintf(filen, "Converge_dU%1.1d",dir);
    sprintf(filen, "CPU_TIME%4.4d.dat", size);
    f = fopen(filen, "a");
    PetscFPrintf(PETSC_COMM_WORLD, f, "t_mom t_p t_pr t_tot %d %le %le %le %le\n",ti, tm_e-tm_s, tp_e-tm_e, tpr_e-tp_e, tpr_e-tm_s);
    fclose(f);    
  }
/* ==================================================================================             *
    BC!
/* ==================================================================================             */ 
/*   if (block_number>1) { */
/*     Block_Interface_U(user); */
/*   } */

/* ==================================================================================             */
/*     End of Flow Solver! */
/* ==================================================================================             */

  return(0);

  }
