#include "variables.h"
#include "petscvec.h"

PetscErrorCode  FormFunctionFEM(SNES snes, Vec x, Vec R, void *ctx);

extern PetscInt   tistart, ti, tiout, period;
extern PetscReal  E, mu, rho, h0, dampfactor, Flux_in;
extern PetscInt   dof, twod, damping, membrane, bending, outghost, ConstitutiveLawNonLinear;
extern PetscInt   timeinteg, nbody, contact, explicit, STRONG_COUPLING;
PetscReal  teta, tetan, tetanm1, Fnormal;
PetscInt   curvature=6;


PetscErrorCode  fem_Initial(IBMNodes *ibm, FE *fem, PetscReal dt) {

  PetscInt  ibi;
  
  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);

  PetscOptionsGetReal(PETSC_NULL, "-E", &E, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-mu", &mu, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-rho", &rho, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-h0", &h0, PETSC_NULL);
 
  PetscOptionsGetInt(PETSC_NULL, "-timeinteg", &timeinteg, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-explicit", &explicit, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-damping", &damping, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-dampfactor", &dampfactor, PETSC_NULL);
 
  PetscOptionsGetInt(PETSC_NULL, "-twod", &twod, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-membrane", &membrane, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-bending", &bending, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-contact", &contact, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-ConstitutiveLawNonLinear", &ConstitutiveLawNonLinear, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-outghost", &outghost, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-curvature", &curvature, PETSC_NULL);

  for (ibi=0; ibi<nbody; ibi++) {
   
    //IO
    Dimension(&ibm[ibi], &fem[ibi], ibi);
    Create(&fem[ibi], ibi);
    Input(&ibm[ibi], &fem[ibi], ibi);
    ContactZ(&fem[ibi]);
    fem_InitGuess(&fem[ibi], ibi);

    if (explicit) {MassDamp(&fem[ibi]);}

    if (tistart) {
      PetscReal  *xx, *xxd, *xxn, *xxnm1;
      PetscInt   nv;

      LocationIn(&fem[ibi], tistart, ibi);

      VecGetArray(fem[ibi].x, &xx);
      VecGetArray(fem[ibi].xd, &xxd);
      VecGetArray(fem[ibi].xn, &xxn);
      VecGetArray(fem[ibi].xnm1, &xxnm1);

      for (nv=0; nv<fem[ibi].n_v; nv++) {
	ibm[ibi].x_bp[nv] = xx[nv*dof  ];
	ibm[ibi].y_bp[nv] = xx[nv*dof+1];
	ibm[ibi].z_bp[nv] = xx[nv*dof+2];
	
	ibm[ibi].x_bp_o[nv] = xxn[nv*dof  ];
	ibm[ibi].y_bp_o[nv] = xxn[nv*dof+1];
	ibm[ibi].z_bp_o[nv] = xxn[nv*dof+2];
	
	ibm[ibi].u[nv].x = xxd[nv*dof  ];
	ibm[ibi].u[nv].y = xxd[nv*dof+1];
	ibm[ibi].u[nv].z = xxd[nv*dof+2];
	
	ibm[ibi].urm1[nv].x = (xxn[nv*dof  ]- xxnm1[nv*dof  ])/dt;
	ibm[ibi].urm1[nv].y = (xxn[nv*dof+1]- xxnm1[nv*dof+1])/dt;
	ibm[ibi].urm1[nv].z = (xxn[nv*dof+2]- xxnm1[nv*dof+2])/dt;

      }//nv
      VecRestoreArray(fem[ibi].x, &xx);
      VecRestoreArray(fem[ibi].xd, &xxd);
      VecRestoreArray(fem[ibi].xn, &xxn);
      VecRestoreArray(fem[ibi].xnm1, &xxnm1);
    } //tistart
    
  }//ibi

  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode  fem_solve (IBMNodes *ibm, FE *fem, PetscInt itr_sc) {

  PetscInt   ibi, k, ec, nv;
  PetscReal  alpha[4];  alpha[0] = 0.25;  alpha[1] = 1./3.;  alpha[2] = 0.5;  alpha[3] = 1.0; 
  SNES       fem_snes;  
  Mat        J;
  //  
  //PetscInt tt; fem->dt=fem->dt/50;
  //
  PetscReal  time=0.0, dt=fem->dt;
  PetscReal  AvePf=0.0, AveTau=0.0;
  //
  //for (tt=0; tt<50; tt++) {
  //
  for (ibi=0; ibi<nbody; ibi++) {
    time = ti*dt;
    PetscPrintf(PETSC_COMM_WORLD, "body:%d Time(%d) =%le\n", ibi, ti, time);

    //if (itr_sc==1) xAccVel(&fem[ibi], itr_sc);

    //-----Copy Pressures from ibm to fem
    AvePf = 0.0;  AveTau = 0.0;  Fnormal = 0.0;
    for (ec=0; ec<fem[ibi].n_elmt; ec++){
      
      fem[ibi].Pf[ec] = -1*(ibm[ibi].pres_p[ec] - ibm[ibi].pres_n[ec]);// + (ibm[ibi].tauN_p[ec] - ibm[ibi].tauN_n[ec]);
      Fnormal += fem[ibi].dA[ec]*fem[ibi].Pf[ec];
      
      AvePf += fem[ibi].Pf[ec];
      AveTau += -(ibm[ibi].tauN_p[ec] - ibm[ibi].tauN_n[ec]);
    }
    AvePf = AvePf/fem[ibi].n_elmt;
    AveTau = AveTau/fem[ibi].n_elmt;
    
    PetscPrintf(PETSC_COMM_WORLD, "Force Normal (%d)=%le\n", ibi, Fnormal);   
      
    //------------------Explicit RK Solver
    if (explicit) {
      for (k=0; k<4; k++) {
    	FormRK(fem[ibi].Res, fem[ibi].x, fem[ibi].xn, fem[ibi].y, fem[ibi].yn, alpha[k], &fem[ibi]);
      }
      
    }else{
      Vec  U;
      VecDuplicate(fem[ibi].x, &U);
      VecCopy(fem[ibi].x, U);
      //SNES
      SNESCreate(PETSC_COMM_SELF, &fem_snes);
      SNESSetFunction(fem_snes, fem[ibi].Res, FormFunctionFEM, (void *)&fem[ibi]);
      SNESAppendOptionsPrefix(fem_snes, "fem_");
      SNESSetFromOptions(fem_snes);
      MatCreateSNESMF(fem_snes, &J);  //MatrixFree
      SNESSetJacobian(fem_snes, J, J, MatMFFDComputeJacobian, (void *)&fem[ibi]);
      SNESSolve(fem_snes, PETSC_NULL, U); //Cannot pass fem[ibi].x
      
      VecCopy(U, fem[ibi].x);
      SNESDestroy(&fem_snes); VecDestroy(&U); MatDestroy(&J);      
    }
        
    if (STRONG_COUPLING) {FormUnderRelaxedSolution(&fem[ibi], itr_sc);}
    
  } //ibi
  //
  //}
  //fem->dt=fem->dt*50.;
  //

  //---------Update the location in Fluid at each SC iterations
  PetscReal  *xx, *xxd;
  for (ibi=0; ibi<nbody; ibi++) {   
    xAccVel(&fem[ibi], itr_sc);
    VecGetArray(fem[ibi].x, &xx);
    VecGetArray(fem[ibi].xd, &xxd);

    for (nv=0; nv<fem[ibi].n_v + fem[ibi].n_ghosts; nv++) {
      fem->x_bp[nv] = xx[nv*dof  ];
      fem->y_bp[nv] = xx[nv*dof+1];
      fem->z_bp[nv] = xx[nv*dof+2];
    }

    for (nv=0; nv<fem[ibi].n_v; nv++) {
      ibm[ibi].x_bp[nv] = xx[nv*dof  ];
      ibm[ibi].y_bp[nv] = xx[nv*dof+1];
      ibm[ibi].z_bp[nv] = xx[nv*dof+2];

      ibm[ibi].u[nv].x = xxd[nv*dof  ];
      ibm[ibi].u[nv].y = xxd[nv*dof+1];
      ibm[ibi].u[nv].z = xxd[nv*dof+2];
    }

    VecRestoreArray(fem[ibi].x, &xx);
    VecRestoreArray(fem[ibi].xd, &xxd);
  }  

  PetscInt     tcyc=period, t;
  t = ti - ((PetscInt)(ti/tcyc))*tcyc;
  if (contact && t>750) {Fcontact(fem);} //dynamic_bhv
  
  // Printout the results
  PetscInt rank=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    if (ti==(ti/tiout)*tiout){
      for (ibi=0; ibi<nbody; ibi++) {
	Output(&fem[ibi], ti, ibi);  
	if (outghost) {OutputGhost(&fem[ibi], ti, ibi);}
	LocationOut(&fem[ibi], ti, ibi); 
      }
    }
  }

  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode  fem_Finalize(FE *fem) {

  PetscInt  ibi;

  for (ibi=0; ibi<nbody; ibi++) { 
    Free(&fem[ibi]);
  }

  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode FormRK(Vec R,Vec x,Vec xn,Vec y,Vec yn,PetscReal alpha,FE *fem) {

  //y=dx/dt
  PetscReal dt=fem->dt;
  PetscReal *xx,*RR, *RRes,*FF;
  PetscInt nv;
  Vec w1,w2;
  VecDuplicate(R,&w1);
  VecDuplicate(R,&w2);
  
  //---------Update the location
  
  VecGetArray(x, &xx);
  for (nv=0;nv<fem->n_v;nv++) {
    fem->x_bp[nv] = xx[nv*dof  ];
    fem->y_bp[nv] = xx[nv*dof+1];
    fem->z_bp[nv] = xx[nv*dof+2];
    if(twod==3){fem->z_bp[nv] = fem->z_bp0[nv];}  //2d case
    if(twod==2){fem->y_bp[nv] = fem->y_bp0[nv];}  //2d case
  }
  VecRestoreArray(x, &xx);
  AreaNormal(fem);
  if (bending){
    PatchLoc(fem); 
    GhostLoc(fem);
  }
  
  //---------Compute Forces then Residual
  VecSet(fem->Fext,0.0);  VecSet(fem->Fint,0.0); 
  VecSet(R,0.0);  VecSet(fem->FJ,0.0); 
  
  FInternal(fem);
  FExternal(fem);
  
  VecWAXPY(R,-1.,fem->Fint,fem->Fext);
  if (contact){VecAXPY(R,-1., fem->Fcnt);}
  
  VecPointwiseMult(w1, fem->Dissip,y);
  VecAXPY(R,-1.,w1);
  VecPointwiseDivide(w2,R,fem->Mass);
  VecCopy(w2,R);
  VecWAXPY(x,alpha*dt,y,xn);
  VecWAXPY(y,alpha*dt,R,yn);
  
  //---------2d case
  if(twod==3){
    VecGetArray(x, &xx);
    VecGetArray(R, &RR); 
    for (nv=0;nv<fem->n_v;nv++) {
      xx[nv*dof+2]=fem->z_bp0[nv];  
      RR[nv*dof+2]=0.0;
    }
    VecRestoreArray(R, &RR);
    VecRestoreArray(x, &xx);
  }
  if(twod==2){
    VecGetArray(x, &xx);
    VecGetArray(R, &RR); 
    for (nv=0;nv<fem->n_v;nv++) {
      xx[nv*dof+1]=fem->y_bp0[nv];  
      RR[nv*dof+1]=0.0;
    }
    VecRestoreArray(R, &RR);
    VecRestoreArray(x, &xx);
  }
  
  VecDestroy(&w1);VecDestroy(&w2);
  
  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode  FormFunctionFEM(SNES snes, Vec x, Vec R, void *ctx) {

  FE         *fem=(FE *)ctx;
  PetscReal  *xx,*RR, *RRes,*FF;
  PetscInt   nv;

  VecCopy(x, fem->x);
  //---------Update the location
  VecGetArray(fem->x, &xx);
  for (nv=0; nv<fem->n_v + fem->n_ghosts; nv++) {
    fem->x_bp[nv] = xx[nv*dof  ];
    fem->y_bp[nv] = xx[nv*dof+1];
    fem->z_bp[nv] = xx[nv*dof+2];
    if (twod) {fem->z_bp[nv] = fem->z_bp0[nv];}  //2d case
  }
  VecRestoreArray(fem->x, &xx);

  //for displacement BCs
  //MoveBoundary(3, fem);
  //

  AreaNormal(fem);
  if (bending){
    if (curvature==1) {
      PatchLoc(fem);
      GhostLoc(fem);
    } else if (curvature==6) {
      //GlobalGhost(fem);
    }
  }

  //---------Compute Forces then Residual
  VecSet(fem->Fext, 0.0);  VecSet(fem->Fint, 0.0);  VecSet(fem->Fdyn, 0.0);
  VecSet(R, 0.0);  VecSet(fem->FJ, 0.0);

  FInternal(fem);
  FDynamic(fem);
  FExternal(fem);

  VecWAXPY(R,-1., fem->Fext, fem->Fint);
  VecAXPY(R,1., fem->Fdyn);
  if (contact) VecAXPY(R,1., fem->Fcnt);

  //---------2d case
  if(twod){
    VecGetArray(fem->x, &xx);
    VecGetArray(R, &RR);
    for (nv=0; nv<fem->n_v; nv++) {
      xx[nv*dof+2]=fem->z_bp0[nv];
      RR[nv*dof+2]=0.0;
    }
    VecRestoreArray(R, &RR);
    VecRestoreArray(fem->x, &xx);
  }

  //for displacement BCs
  //MoveBoundary(3, fem);
  //
    
  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode  FInternal(FE *fem) {
 
  PetscInt       i;
  PetscReal      Fm[9], Fb[42], Fint[9];
  Cmpnts         x1, x2, x3, X1, X2, X3;
  PetscInt       ec, n1e, n2e, n3e, n4e, n5e, n6e, j;

  for (i=0; i<9; i++) {Fm[i] = 0.0;  Fint[i] = 0.0;}
  for (i=0; i<42; i++) {Fb[i] = 0.0;}

  PetscReal  *FF,*FFJ;
  VecGetArray(fem->Fint, &FF);
  VecGetArray(fem->FJ, &FFJ);

  for (ec=0; ec<fem->n_elmt; ec++) {
    n1e=fem->nv1[ec];n2e=fem->nv2[ec];n3e=fem->nv3[ec];
    n4e=fem->nv4[ec];n5e=fem->nv5[ec];n6e=fem->nv6[ec];
    
    //current location
    x1.x=fem->x_bp[n1e]; x1.y=fem->y_bp[n1e]; x1.z=fem->z_bp[n1e];
    x2.x=fem->x_bp[n2e]; x2.y=fem->y_bp[n2e]; x2.z=fem->z_bp[n2e];
    x3.x=fem->x_bp[n3e]; x3.y=fem->y_bp[n3e]; x3.z=fem->z_bp[n3e];
    //initial location
    X1.x=fem->x_bp0[n1e]; X1.y=fem->y_bp0[n1e]; X1.z=fem->z_bp0[n1e];
    X2.x=fem->x_bp0[n2e]; X2.y=fem->y_bp0[n2e]; X2.z=fem->z_bp0[n2e];
    X3.x=fem->x_bp0[n3e]; X3.y=fem->y_bp0[n3e]; X3.z=fem->z_bp0[n3e];

    if (curvature==1) {
      if (membrane) {
	Fmembrane(ec, X1, X2, X3, x1, x2, x3, Fm, fem);
      }
    }

    if (bending) {
      Fbending(ec, X1, X2, X3, x1, x2, x3, Fb, fem);
    }

    if (curvature==1) { //Only bending     
      for (i=0; i<9; i++) {
	Fint[i]=Fm[i]+Fb[i];
      }  //end loop over nodes of each element
      
      FF[n1e*dof] +=Fint[0];
      FF[n1e*dof+1] +=Fint[1];
      FF[n1e*dof+2] +=Fint[2];
      
      FF[n2e*dof] +=Fint[3];
      FF[n2e*dof+1] +=Fint[4];
      FF[n2e*dof+2] +=Fint[5];
      
      FF[n3e*dof] +=Fint[6];
      FF[n3e*dof+1] +=Fint[7];
      FF[n3e*dof+2] +=Fint[8];     
  
      if(n4e !=1000000){ // Front node is inside domain
	FF[n4e*dof] +=Fb[9];
	FF[n4e*dof+1] +=Fb[10];
	FF[n4e*dof+2] +=Fb[11];
      }else{ // Front node is ghost
	FFJ[n1e*dof] +=Fb[9];
	FFJ[n1e*dof+1] +=Fb[10];
	FFJ[n1e*dof+2] +=Fb[11];
      }
      
      if(n5e !=1000000){// Front node is inside domain
	FF[n5e*dof] +=Fb[12];
	FF[n5e*dof+1] +=Fb[13];
	FF[n5e*dof+2] +=Fb[14];
      }else{ // Front node is ghost
	FFJ[n2e*dof] +=Fb[12];
	FFJ[n2e*dof+1] +=Fb[13];
	FFJ[n2e*dof+2] +=Fb[14];
      }
      
      if(n6e !=1000000){// Front node is inside domain
	FF[n6e*dof] +=Fb[15];
	FF[n6e*dof+1] +=Fb[16];
	FF[n6e*dof+2] +=Fb[17];
      }else{ // Front node is ghost
	FFJ[n3e*dof] +=Fb[15];
	FFJ[n3e*dof+1] +=Fb[16];
	FFJ[n3e*dof+2] +=Fb[17];
      }

    } else if (curvature==6) {
      
      PetscInt  node, v = fem->val[ec];

      for (i=0; i<(v+6); i++) {
      	if (fem->patch[16*ec+i]!=1000000) {
      	  node = fem->patch[16*ec+i];

      	  FF[dof*node] += Fb[dof*i];
      	  FF[dof*node+1] += Fb[dof*i+1];
      	  FF[dof*node+2] += Fb[dof*i+2];
      	 
      	}
      }      
    }
  }//end loop over elements
  
  VecRestoreArray(fem->FJ, &FFJ); 
  VecRestoreArray(fem->Fint, &FF); 
  
  if (bending && curvature==1) {ModifyFbending(fem);}
  
  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode  FDynamic(FE *fem) {

  PetscInt   ec, i;
  PetscReal  M[9], C[9], Fd[9], x[9], xn[9], xnm1[9], xd[9], xdd[9];
  PetscInt   n1e, n2e, n3e;
  PetscReal  Gama=0.5, Beta=0.25, dt=fem->dt;
  PetscReal  M1, C1, M2, C2, M3, C3;

  for (i=0; i<9; i++) {M[i] = 0.0;  C[i] = 0.0;  Fd[i] = 0.0;}
  M1=1./(Beta*pow(dt,2));  M2=1./(Beta*dt);  M3=(1./(2*Beta))-1.;
  C1=Gama/(Beta*dt);  C2=(Gama/Beta)-1;  C3=dt*(Gama/(2.*Beta)-1);

  PetscReal  *xx,*xxn,*xxnm1,*xxd,*xxdd,*FF,facc,fvel;
  VecGetArray(fem->xd, &xxd);
  VecGetArray(fem->xdd, &xxdd);
  VecGetArray(fem->xn, &xxn);
  VecGetArray(fem->xnm1, &xxnm1);
  VecGetArray(fem->x, &xx);
  VecGetArray(fem->Fdyn, &FF);

  for (ec=0; ec<fem->n_elmt+2*fem->n_ghosts; ec++) {
    n1e=fem->nv1[ec];n2e=fem->nv2[ec];n3e=fem->nv3[ec];
  
    x[0]=fem->x_bp[n1e];  x[1]=fem->y_bp[n1e]; x[2]=fem->z_bp[n1e];
    x[3]=fem->x_bp[n2e];  x[4]=fem->y_bp[n2e]; x[5]=fem->z_bp[n2e];
    x[6]=fem->x_bp[n3e];  x[7]=fem->y_bp[n3e]; x[8]=fem->z_bp[n3e];
           
    xn[0]=xxn[n1e*dof];  xn[1]=xxn[n1e*dof+1];  xn[2]=xxn[n1e*dof+2];
    xn[3]=xxn[n2e*dof];  xn[4]=xxn[n2e*dof+1];  xn[5]=xxn[n2e*dof+2];
    xn[6]=xxn[n3e*dof];  xn[7]=xxn[n3e*dof+1];  xn[8]=xxn[n3e*dof+2];

    xnm1[0]=xxnm1[n1e*dof];  xnm1[1]=xxnm1[n1e*dof+1];  xnm1[2]=xxnm1[n1e*dof+2];
    xnm1[3]=xxnm1[n2e*dof];  xnm1[4]=xxnm1[n2e*dof+1];  xnm1[5]=xxnm1[n2e*dof+2];
    xnm1[6]=xxnm1[n3e*dof];  xnm1[7]=xxnm1[n3e*dof+1];  xnm1[8]=xxnm1[n3e*dof+2];

    xd[0]=xxd[n1e*dof];  xd[1]=xxd[n1e*dof+1];  xd[2]=xxd[n1e*dof+2];
    xd[3]=xxd[n2e*dof];  xd[4]=xxd[n2e*dof+1];  xd[5]=xxd[n2e*dof+2];
    xd[6]=xxd[n3e*dof];  xd[7]=xxd[n3e*dof+1];  xd[8]=xxd[n3e*dof+2];

    xdd[0]=xxdd[n1e*dof];  xdd[1]=xxdd[n1e*dof+1];  xdd[2]=xxdd[n1e*dof+2];
    xdd[3]=xxdd[n2e*dof];  xdd[4]=xxdd[n2e*dof+1];  xdd[5]=xxdd[n2e*dof+2];
    xdd[6]=xxdd[n3e*dof];  xdd[7]=xxdd[n3e*dof+1];  xdd[8]=xxdd[n3e*dof+2];
     
    Mass(fem, ec, M);
    if(damping) {Damp(M,C);}

    if(timeinteg==0){ //Newmark constant average acceleration
      for (i=0; i<9; i++) {
	facc = M[i]*(M1*xn[i] + M2*xd[i] + M3*xdd[i]);
	  fvel = C[i]*(C1*xn[i] + C2*xd[i] + C3*xdd[i]);
	  Fd[i] = M[i]*M1*x[i] + C[i]*C1*x[i] - facc - fvel;
      }

    }else if (timeinteg==1){ //central
      for (i=0; i<9; i++) {
	Fd[i] = M[i]*(x[i] - 2*xn[i] + xnm1[i])/(dt*dt) + C[i]*(x[i] - xn[i])/dt;
      }
    }
    //AddToVector
    
    FF[n1e*dof] +=Fd[0];
    FF[n1e*dof+1] +=Fd[1];
    FF[n1e*dof+2] +=Fd[2];
    
    FF[n2e*dof] +=Fd[3];
    FF[n2e*dof+1] +=Fd[4];
    FF[n2e*dof+2] +=Fd[5];

    FF[n3e*dof] +=Fd[6];
    FF[n3e*dof+1] +=Fd[7];
    FF[n3e*dof+2] +=Fd[8];
   
  }//end loop over elements

  VecRestoreArray(fem->x, &xx);
  VecRestoreArray(fem->xd, &xxd);
  VecRestoreArray(fem->xdd, &xxdd);
  VecRestoreArray(fem->xnm1, &xxnm1);
  VecRestoreArray(fem->xn, &xxn);
  VecRestoreArray(fem->Fdyn, &FF);
  
  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode  xAccVel(FE *fem, PetscInt itr_sc) {

  PetscReal      M1, C1, M2, C2, M3, C3;
  PetscReal      *xx, *xxn, *xxd, *xxdsc, *xxdd, *xxddsc, *fdyn, *fcnt;
  PetscReal      Gama=0.5, Beta=0.25, dt=fem->dt;
  PetscInt       nv;
  Cmpnts         u, a;
  
  M1 = 1./(Beta*pow(dt,2));  M2 = 1./(Beta*dt);  M3 = (1./(2*Beta)) - 1.;
  C1 = Gama/(Beta*dt);  C2 = (Gama/Beta) - 1;  C3 = dt*(Gama/(2.*Beta) - 1);
  
  VecGetArray(fem->xd, &xxd);
  VecGetArray(fem->xdsc, &xxdsc);
  VecGetArray(fem->xdd, &xxdd);
  VecGetArray(fem->xddsc, &xxddsc);
  VecGetArray(fem->xn, &xxn);
  VecGetArray(fem->x, &xx);
  VecGetArray(fem->Fdyn, &fdyn);
  VecGetArray(fem->Fcnt, &fcnt);

  for (nv=0; nv<fem->n_v+fem->n_ghosts; nv++) {

    a.x = M1*(xx[nv*dof] - xxn[nv*dof]) - M2*xxd[nv*dof] - M3*xxdd[nv*dof];
    a.y = M1*(xx[nv*dof+1] - xxn[nv*dof+1]) - M2*xxd[nv*dof+1] - M3*xxdd[nv*dof+1];
    a.z = M1*(xx[nv*dof+2] - xxn[nv*dof+2]) - M2*xxd[nv*dof+2] - M3*xxdd[nv*dof+2];

    u.x = C1*(xx[nv*dof]-xxn[nv*dof]) - C2*xxd[nv*dof] - C3*xxdd[nv*dof];
    u.y = C1*(xx[nv*dof+1]-xxn[nv*dof+1]) - C2*xxd[nv*dof+1] - C3*xxdd[nv*dof+1];
    u.z = C1*(xx[nv*dof+2]-xxn[nv*dof+2]) - C2*xxd[nv*dof+2] - C3*xxdd[nv*dof+2];   
      
    if (fem->contact[nv]==1) {
      xxdd[dof*nv] = 0.;  xxdd[dof*nv+1] = 0.;  xxdd[dof*nv+2] = 0.;
      xxd[dof*nv] = 0.;   xxd[dof*nv+1] = 0.;   xxd[dof*nv+2] = 0.;
      fcnt[dof*nv] += fdyn[dof*nv];  fcnt[dof*nv+1] += fdyn[dof*nv+1];  fcnt[dof*nv+2] += fdyn[dof*nv+2];

    } else {
      xxdd[dof*nv] = a.x;  xxdd[dof*nv+1] = a.y;  xxdd[dof*nv+2] = a.z;
      xxd[dof*nv] = u.x;   xxd[dof*nv+1] = u.y;   xxd[dof*nv+2] = u.z;
      fcnt[dof*nv] = 0.;  fcnt[dof*nv+1] = 0.;  fcnt[dof*nv+2] = 0.;
    }    
  }

  //if (itr_sc==1) {
    VecCopy(fem->xn, fem->xnm1);
    VecCopy(fem->x, fem->xn);
    if (explicit) {VecCopy(fem->y, fem->yn);}
    //}

  VecRestoreArray(fem->x, &xx);
  VecRestoreArray(fem->xd, &xxd);
  VecRestoreArray(fem->xdsc, &xxdsc);
  VecRestoreArray(fem->xdd, &xxdd);
  VecRestoreArray(fem->xddsc, &xxddsc);
  VecRestoreArray(fem->xn, &xxn);
  VecRestoreArray(fem->Fdyn, &fdyn);
  VecRestoreArray(fem->Fcnt, &fcnt);

  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode  fem_InitGuess(FE *fem, PetscInt ibi) {
 
  PetscReal  *xx;
  PetscInt   nv, ec;

  if (curvature==6) {GlobalGhostInit(fem);}

  AreaNormal(fem);
  for (ec=0; ec<fem->n_elmt + 2*fem->n_ghosts; ec++) {
    fem->dA0[ec] = fem->dA[ec]; 
    fem->Nf_x[ec] = fem->nf_x[ec];  fem->Nf_y[ec] = fem->nf_y[ec];  fem->Nf_z[ec] = fem->nf_z[ec]; 
  }
  
  if (bending){InitGhost(fem);}

  VecGetArray(fem->x, &xx);
  for (nv=0; nv<fem->n_v+fem->n_ghosts; nv++) {
    xx[nv*dof] = fem->x_bp0[nv];
    xx[nv*dof+1] = fem->y_bp0[nv];
    xx[nv*dof+2] = fem->z_bp0[nv];
  }

  VecRestoreArray(fem->x, &xx);
  VecCopy(fem->x, fem->xn);
  VecCopy(fem->x, fem->xnm1);
  
  if (ConstitutiveLawNonLinear) {
    InitMaterial(fem);
  }
  
  if (!tistart)  Output(fem, 0, ibi);
  if (!tistart && outghost)  OutputGhost(fem, 0, ibi);
  
  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode  FormUnderRelaxedSolution (FE *fem, PetscInt itr_sc) {

  PetscReal  wAit, wAitnsc=fem->wAit, wAitCalc, norm2=0, norm=0;
  Vec        Ddx;
  
  //x~:after scaling
  VecCopy(fem->dx, fem->dxn); //dx(l)=dx(l+1) 
  VecWAXPY(fem->dx,-1., fem->xnsc, fem->x); //dx(l+1)=x(l+1)-x~(l)
    
  if (itr_sc==1) {
    wAit = 0.7;  //1.0 
    VecCopy(fem->xn, fem->xnsc);
  }else{  
    VecDuplicate(fem->x, &Ddx);
    VecWAXPY(Ddx, -1., fem->dxn, fem->dx); //Ddx=dx(l+1)-dx(l)
    
    VecNorm(fem->dx, NORM_1, &norm);
    VecNorm(Ddx, NORM_1, &norm2);
    
    if (norm2>1.e-8) {
      wAit = wAitnsc*(1. - (norm/norm2));
    }else{
      wAit = wAitnsc;
    }
    VecDestroy(&Ddx);
  }
  
  wAitCalc = wAit;
  
  if (wAit>1.0)  wAit = 1.0;
  if (wAit<0.2)  wAit = 0.2;  // 0.75 0.2
  
  fem->wAit = wAit;
  PetscPrintf(PETSC_COMM_WORLD, "it-SC= %d Aitken = %le old-SC = %le Calculated=%le--- norm = %le norm2 = %le\n", itr_sc, wAit, wAitnsc, wAitCalc, norm, norm2);
  
  //Scale x~(l+1)=w x(l+1)+(1-w) x~(l)
  VecScale(fem->x, wAit);
  VecAXPY(fem->x, (1.-wAit), fem->xnsc);
  
  VecWAXPY(fem->dxsc, -1., fem->xnsc, fem->x); //dx~(l+1)=x~(l+1)-x~(l)
  VecCopy(fem->x, fem->xnsc); //x~(l)=x~(l+1)
  
  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode Free(FE *fem){

//Vec
 VecDestroy(&(fem->Res)); VecDestroy(&(fem->x));  VecDestroy(&(fem->xn));  VecDestroy(&(fem->xnm1)); VecDestroy(&(fem->V));
 VecDestroy(&(fem->xd));  VecDestroy(&(fem->xdd));  VecDestroy(&(fem->y));  VecDestroy(&(fem->yn));
 VecDestroy(&(fem->Fext)); VecDestroy(&(fem->Fint)); VecDestroy(&(fem->Fdyn));
 VecDestroy(&(fem->disp)); VecDestroy(&(fem->FJ)); VecDestroy(&(fem->Fcnt));
 VecDestroy(&(fem->Mass)); VecDestroy(&(fem->Dissip));  VecDestroy(&(fem->dx)); VecDestroy(&(fem->dxn));
 VecDestroy(&(fem->dxsc)); VecDestroy(&(fem->xnsc));

 //Malloc
 PetscFree(fem->x_bp); PetscFree(fem->y_bp); PetscFree(fem->z_bp);
 PetscFree(fem->x_bp0); PetscFree(fem->y_bp0); PetscFree(fem->z_bp0);
 PetscFree(fem->nv1); PetscFree(fem->nv2); PetscFree(fem->nv3);
 PetscFree(fem->nv4); PetscFree(fem->nv5); PetscFree(fem->nv6);
 PetscFree(fem->n_bnodes); PetscFree(fem->bnodes); PetscFree(fem->kve0); 
 PetscFree(fem->kve); PetscFree(fem->StressM); PetscFree(fem->StrainM);
PetscFree(fem->StressB); PetscFree(fem->StrainB);PetscFree(fem->Pf);PetscFree(fem->Pfn);


 if(bending){
 PetscFree(fem->belmtsedge);  PetscFree(fem->belmts);  PetscFree(fem->edgefrontnodes);
 PetscFree(fem->edgefrontnodesI); 
 PetscFree(fem->p4x);  PetscFree(fem->p4y);  PetscFree(fem->p4z);
 PetscFree(fem->p5x);  PetscFree(fem->p5y);  PetscFree(fem->p5z);
 PetscFree(fem->p6x);  PetscFree(fem->p6y);  PetscFree(fem->p6z);

 PetscFree(fem->p4x0);  PetscFree(fem->p4y0);  PetscFree(fem->p4z0);
 PetscFree(fem->p5x0);  PetscFree(fem->p5y0);  PetscFree(fem->p5z0);
 PetscFree(fem->p6x0);  PetscFree(fem->p6y0);  PetscFree(fem->p6z0);
 }

 return(0);
}

//-------------------------------------------------------------------------------------
PetscErrorCode ContactZ(FE *fem) {

  PetscInt  nv;

  for (nv=0; nv<fem->n_v; nv++) {
    fem->contact[nv] = 0;
  }

  return(0);
}
