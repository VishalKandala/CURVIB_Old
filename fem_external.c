#include "variables.h"

extern const PetscInt dof,nbody;
extern const PetscReal h0;
extern PetscInt initrest, BHV, inv_flg;
extern Cmpnts PLUS(Cmpnts v1, Cmpnts v2);
extern Cmpnts MINUS(Cmpnts v1, Cmpnts v2);
extern Cmpnts CROSS(Cmpnts v1, Cmpnts v2);
extern PetscReal DOT(Cmpnts v1, Cmpnts v2);
extern Cmpnts UNIT(Cmpnts v1);
extern PetscReal SIZE(Cmpnts v1);
extern Cmpnts AMULT(PetscReal alpha, Cmpnts v1);
extern PetscErrorCode INV(PetscReal T[3][3],PetscReal _Tinv[3][3]);
extern PetscErrorCode MATMULT(PetscReal A[][2],PetscReal B[][2], PetscReal C[][2]);
extern PetscErrorCode TRANS(PetscReal A[][2], PetscReal AT[][2]);

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode NodeForce(PetscInt nv,PetscReal F,PetscInt dir,FE *fem){

  PetscReal *FF;
  
  VecGetArray(fem->Fext, &FF);
  FF[nv*dof+dir]=F;
  VecRestoreArray(fem->Fext, &FF);
 return(0); 
}

//------------------------------------------------------------------------------------------------------------
PetscErrorCode EdgeClamp(PetscInt edge_n, FE *fem) {

  PetscReal *FFint, *FFext, *FFdyn;
  PetscInt start=0, end=0, edge, nbc, nb, ec;

  /* for (edge=0; edge<edge_n+1; edge++) { */
  /*   end += ibm->n_bnodes[edge]; */
  /* } */
  /* start = end - ibm->n_bnodes[edge_n]; */

  VecGetArray(fem->Fint, &FFint);
  VecGetArray(fem->Fext, &FFext);
  VecGetArray(fem->Fdyn, &FFdyn);

  /* for (nbc=start; nbc<end; nbc++) { //fix boundary nodes */
  /*   nb=ibm->bnodes[nbc]; */

  /*   FFint[nb*dof] =0.0; */
  /*   FFint[nb*dof+1] =0.0; */
  /*   FFint[nb*dof+2] =0.0; */

  /*   FFext[nb*dof] =0.0; */
  /*   FFext[nb*dof+1] =0.0; */
  /*   FFext[nb*dof+2] =0.0; */

  /*   FFdyn[nb*dof] =0.0; */
  /*   FFdyn[nb*dof+1] =0.0; */
  /*   FFdyn[nb*dof+2] =0.0; */

  /* } */

  /* for (nb=ibm->n_v; nb<ibm->n_v+ibm->n_ghosts; nb++) { //fix all ghost nodes */
  /*   FFint[nb*dof] =0.0; */
  /*   FFint[nb*dof+1] =0.0; */
  /*   FFint[nb*dof+2] =0.0; */
/*   FFext[nb*dof] =0.0; */
  /*   FFext[nb*dof+1] =0.0; */
  /*   FFext[nb*dof+2] =0.0; */

  /*   FFdyn[nb*dof] =0.0; */
  /*   FFdyn[nb*dof+1] =0.0; */
  /*   FFdyn[nb*dof+2] =0.0; */

  /* } */

  /* for (nb=ibm->n_v+ibm->n_bnodes[edge_n]-1; nb<ibm->n_v+ibm->n_ghosts; nb++) { //fix ghost nodes of BHV */
  /*   FFint[nb*dof] =0.0; */
  /*   FFint[nb*dof+1] =0.0; */
  /*   FFint[nb*dof+2] =0.0; */

  /*   FFext[nb*dof] =0.0; */
  /*   FFext[nb*dof+1] =0.0; */
  /*   FFext[nb*dof+2] =0.0; */

  /*   FFdyn[nb*dof] =0.0; */
  /*   FFdyn[nb*dof+1] =0.0; */
  /*   FFdyn[nb*dof+2] =0.0; */

  /* } */

  /* if (curvature ==6) { //clamped edge */
  /*   for (ec=0; ec<ibm->n_ghosts; ec++) { */
  /*     if (ibm->belmtsedge[ec]==edge_n) { */
  /*    nb = ibm->edgefrontnodes[ec];       */

  /*    FFint[nb*dof] = 0.0; */
  /*    FFint[nb*dof+1] = 0.0; */
  /*    FFint[nb*dof+2] = 0.0; */

  /*    FFext[nb*dof] = 0.0; */
  /*    FFext[nb*dof+1] = 0.0; */
  /*    FFext[nb*dof+2] = 0.0; */

  /*    FFdyn[nb*dof] = 0.0; */
 /*    FFdyn[nb*dof+1] = 0.0; */
  /*    FFdyn[nb*dof+2] = 0.0; */

  /*     } */
  /*   } */
  /* } */

  /* if (curvature ==6) { //general clamped edge for curvature 6 */
  /*   for (nbc=0; nbc<ibm->n_ghosts; nbc++) { */
  /*     nb = ibm->edgefrontnodes[nbc]; */

  /*     ibm->x_bp[nb] = ibm->x_bp0[nb]; */
  /*     ibm->y_bp[nb] = ibm->y_bp0[nb]; */
  /*     ibm->z_bp[nb] = ibm->z_bp0[nb]; */

  /*     FFint[nb*dof] =0.0; */
  /*     FFint[nb*dof+1] =0.0; */
  /*     FFint[nb*dof+2] =0.0; */

  /*     FFext[nb*dof] =0.0; */
  /*     FFext[nb*dof+1] =0.0; */
  /*     FFext[nb*dof+2] =0.0; */

  /*     FFdyn[nb*dof] =0.0; */
  /*     FFdyn[nb*dof+1] =0.0; */
  /*     FFdyn[nb*dof+2] =0.0; */

  /*   } */
  /* } */

    for (nb=0; nb<fem->n_v+fem->n_ghosts; nb++) {
      if (fem->z_bp0[nb]>2.955) {//for canti
        //if (fem->x_bp0[nb]<0.042) {//for plate
        fem->x_bp[nb] = fem->x_bp0[nb];
        fem->y_bp[nb] = fem->y_bp0[nb];
        fem->z_bp[nb] = fem->z_bp0[nb];

        FFint[nb*dof] =0.0;
        FFint[nb*dof+1] =0.0;
        FFint[nb*dof+2] =0.0;

        FFext[nb*dof] =0.0;
        FFext[nb*dof+1] =0.0;
        FFext[nb*dof+2] =0.0;

        FFdyn[nb*dof] =0.0;
        FFdyn[nb*dof+1] =0.0;
        FFdyn[nb*dof+2] =0.0;

      }
    }

  VecRestoreArray(fem->Fdyn, &FFdyn);
  VecRestoreArray(fem->Fint, &FFint);
  VecRestoreArray(fem->Fext, &FFext);

  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode BuoyantForce(PetscReal F,PetscInt dir,FE *fem){
 PetscInt ec,n1e,n2e,n3e;
  PetscReal *FF;

  VecGetArray(fem->Fext, &FF);
  for (ec=0; ec<fem->n_elmt; ec++) {
    n1e=fem->nv1[ec];n2e=fem->nv2[ec];n3e=fem->nv3[ec];

    FF[n1e*dof+dir]+=F/3.;
    FF[n2e*dof+dir]+=F/3.;
    FF[n3e*dof+dir]+=F/3.;
      
  }//end loop over elements
  VecRestoreArray(fem->Fext, &FF);

  return(0);

}
  
//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode  EdgeFix(PetscInt edge_n, FE *fem){
 
  PetscReal  *FFint, *FFext, *FFdyn, *xx, *xdd;
  PetscInt   start=0, end=0, edge, nbc, nb;

  for (edge=0; edge<edge_n+1; edge++) {
    end += fem->n_bnodes[edge];
  }
  start = end - fem->n_bnodes[edge_n];
  
  VecGetArray(fem->Fint, &FFint);
  VecGetArray(fem->Fext, &FFext);
  VecGetArray(fem->Fdyn, &FFdyn);
  VecGetArray(fem->x, &xx);
  VecGetArray(fem->xd, &xdd);

  for (nbc=start; nbc<end; nbc++) { //fix boundary nodes
    nb=fem->bnodes[nbc];

    fem->x_bp[nb] = fem->x_bp0[nb]; //for kinematic contact
    fem->y_bp[nb] = fem->y_bp0[nb];
    fem->z_bp[nb] = fem->z_bp0[nb];
    xx[nb*dof] = fem->x_bp0[nb];
    xx[nb*dof+1] = fem->y_bp0[nb];
    xx[nb*dof+2] = fem->z_bp0[nb];
    fem->contact[nb] = 0;
    xdd[nb*dof] = 0.;
    xdd[nb*dof+1] = 0.;
    xdd[nb*dof+2] = 0.;

    FFint[nb*dof] =0.0;
    FFint[nb*dof+1] =0.0;
    FFint[nb*dof+2] =0.0;
    
    FFext[nb*dof] =0.0;
    FFext[nb*dof+1] =0.0;
    FFext[nb*dof+2] =0.0;
    
    FFdyn[nb*dof] =0.0;
    FFdyn[nb*dof+1] =0.0;
    FFdyn[nb*dof+2] =0.0;

  }

  /* for (nb=fem->n_v; nb<fem->n_v+fem->n_ghosts; nb++) { //fix all ghost nodes */
  /*   FFint[nb*dof] =0.0; */
  /*   FFint[nb*dof+1] =0.0; */
  /*   FFint[nb*dof+2] =0.0; */
    
  /*   FFext[nb*dof] =0.0; */
  /*   FFext[nb*dof+1] =0.0; */
  /*   FFext[nb*dof+2] =0.0; */
    
  /*   FFdyn[nb*dof] =0.0; */
  /*   FFdyn[nb*dof+1] =0.0; */
  /*   FFdyn[nb*dof+2] =0.0; */
    
  /* } */
  if (BHV) {
    for (nb=fem->n_v+fem->n_bnodes[edge_n-1]-1; nb<fem->n_v+fem->n_ghosts; nb++) { //fix ghost nodes of BHV Note:number of ghost nodes are one shorter than boundary nodes on edges 
      FFint[nb*dof] = 0.0;
      FFint[nb*dof+1] = 0.0;
      FFint[nb*dof+2] = 0.0;
      
      FFext[nb*dof] = 0.0;
      FFext[nb*dof+1] = 0.0;
      FFext[nb*dof+2] = 0.0;
      
      FFdyn[nb*dof] = 0.0;
      FFdyn[nb*dof+1] = 0.0;
      FFdyn[nb*dof+2] = 0.0;
      
    }
  } else if (inv_flg) {
    for (nb=0; nb<fem->n_v+fem->n_ghosts; nb++) {
      if (fem->z_bp0[nb]>2.9575) { //2.955 2.0383
	
	FFint[nb*dof] = 0.0;
	FFint[nb*dof+1] = 0.0;
	FFint[nb*dof+2] = 0.0;
	
	FFext[nb*dof] = 0.0;
	FFext[nb*dof+1] = 0.0;
	FFext[nb*dof+2] = 0.0;
	
	FFdyn[nb*dof] = 0.0;
	FFdyn[nb*dof+1] = 0.0;
	FFdyn[nb*dof+2] = 0.0;
	
      }
    }
  }

  VecRestoreArray(fem->Fdyn, &FFdyn);
  VecRestoreArray(fem->Fint, &FFint);
  VecRestoreArray(fem->Fext, &FFext);
  VecRestoreArray(fem->x, &xx);
  VecRestoreArray(fem->xd, &xdd);
 
  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode  SurfaceNormalPressure(FE *fem) {

  PetscInt  i, ec, n1e, n2e, n3e;
  PetscReal A0, *FF, P;
  Cmpnts    x1, x2, x3, dx21, dx31, X1, X2, X3, dX21, dX31, n;
  PetscInt  start, end, edge_n, edge, bn;

  VecGetArray(fem->Fext, &FF);

  for (ec=0; ec<fem->n_elmt; ec++) {
    n1e = fem->nv1[ec];  n2e = fem->nv2[ec];  n3e = fem->nv3[ec];
       
    A0 = fem->dA[ec];
    n.x = fem->nf_x[ec];  n.y = fem->nf_y[ec];  n.z = fem->nf_z[ec];
    P = fem->Pf[ec];

    //P = 4.7;
    
    FF[n1e*dof] += n.x*P*A0/3.;
    FF[n1e*dof+1] += n.y*P*A0/3.;
    FF[n1e*dof+2] += n.z*P*A0/3.;
    
    FF[n2e*dof] += n.x*P*A0/3.;
    FF[n2e*dof+1] += n.y*P*A0/3.;
    FF[n2e*dof+2] += n.z*P*A0/3.;
    
    FF[n3e*dof] += n.x*P*A0/3.;
    FF[n3e*dof+1] += n.y*P*A0/3.;
    FF[n3e*dof+2] += n.z*P*A0/3.;
    
  }//end loop over elements

  VecRestoreArray(fem->Fext, &FF);
 
  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode  SurfaceNormalPressure2(FE *fem) {

  PetscInt  i, ec, n1e, n2e, n3e, nc;
  PetscReal A, *FF, P, *dr_v;
  Cmpnts    n, cent, x1, x2, x3;
  PetscInt  start, end, edge_n, edge, bn;

  PetscMalloc(fem->n_v*sizeof(PetscReal), &dr_v);
  for (nc=0; nc<fem->n_v; nc++)  dr_v[nc] = 0.;

  VecGetArray(fem->Fext, &FF);

  for (ec=0; ec<fem->n_elmt; ec++) {
    n1e = fem->nv1[ec];  n2e = fem->nv2[ec];  n3e = fem->nv3[ec];
       
    A = fem->dA[ec];
    n.x = fem->nf_x[ec];  n.y = fem->nf_y[ec];  n.z = fem->nf_z[ec];
    //P = fem->Pf[ec];
    P = 0.2;
    cent.x = (fem->x_bp[n1e] + fem->x_bp[n2e] + fem->x_bp[n3e])/3.; 
    cent.y = (fem->y_bp[n1e] + fem->y_bp[n2e] + fem->y_bp[n3e])/3.;
    cent.z = (fem->z_bp[n1e] + fem->z_bp[n2e] + fem->z_bp[n3e])/3.;
 
    x1.x = fem->x_bp[n1e];  x1.y = fem->y_bp[n1e];  x1.z = fem->z_bp[n1e];
    x2.x = fem->x_bp[n2e];  x2.y = fem->y_bp[n2e];  x2.z = fem->z_bp[n2e];
    x3.x = fem->x_bp[n3e];  x3.y = fem->y_bp[n3e];  x3.z = fem->z_bp[n3e];

    dr_v[n1e] += 1./SIZE(MINUS(x1, cent));
    dr_v[n2e] += 1./SIZE(MINUS(x2, cent));
    dr_v[n3e] += 1./SIZE(MINUS(x3, cent));

    FF[n1e*dof] += n.x*P*A/SIZE(MINUS(x1, cent));
    FF[n1e*dof+1] += n.y*P*A/SIZE(MINUS(x1, cent));
    FF[n1e*dof+2] += n.z*P*A/SIZE(MINUS(x1, cent));
    
    FF[n2e*dof] += n.x*P*A/SIZE(MINUS(x2, cent));
    FF[n2e*dof+1] += n.y*P*A/SIZE(MINUS(x2, cent));
    FF[n2e*dof+2] += n.z*P*A/SIZE(MINUS(x2, cent));
    
    FF[n3e*dof] += n.x*P*A/SIZE(MINUS(x3, cent));
    FF[n3e*dof+1] += n.y*P*A/SIZE(MINUS(x3, cent));
    FF[n3e*dof+2] += n.z*P*A/SIZE(MINUS(x3, cent));
    
  }

  for (nc=0; nc<fem->n_v; nc++) {
    FF[nc*dof] /= dr_v[nc];
    FF[nc*dof+1] /= dr_v[nc];
    FF[nc*dof+2] /= dr_v[nc];
  }


  PetscFree(dr_v);
  VecRestoreArray(fem->Fext, &FF); 

  return(0);
}

//------------------------------------------------------------------------------------------------------------
PetscErrorCode GhostFix(PetscInt edge_n,FE *fem){

  
  PetscInt n1e,n2e,n3e,ec,be,i;

  Cmpnts N,X1,X2,X3,x1,x2,x3,x4,x5,x6,a1,a2,a3,A1,A2,A3,e1,e2,e3,n1,n2,n3,h1,h2,h3;
  for (ec=0; ec<fem->n_ghosts; ec++) { //Go through ghost nodes
    be=fem->belmts[ec];
    n1e=fem->nv1[be];n2e=fem->nv2[be];n3e=fem->nv3[be];

    //currentlocation
    x1.x=fem->x_bp[n1e]; x1.y=fem->y_bp[n1e]; x1.z=fem->z_bp[n1e];
    x2.x=fem->x_bp[n2e]; x2.y=fem->y_bp[n2e]; x2.z=fem->z_bp[n2e];
    x3.x=fem->x_bp[n3e]; x3.y=fem->y_bp[n3e]; x3.z=fem->z_bp[n3e];
    //initial location
    //  if(!initrest){
     X1.x=fem->x_bp0[n1e]; X1.y=fem->y_bp0[n1e]; X1.z=fem->z_bp0[n1e];
     X2.x=fem->x_bp0[n2e]; X2.y=fem->y_bp0[n2e]; X2.z=fem->z_bp0[n2e];
     X3.x=fem->x_bp0[n3e]; X3.y=fem->y_bp0[n3e]; X3.z=fem->z_bp0[n3e];
  /*   }else{ */
/*       X1.x=fem->x_bpr[n1e]; X1.y=fem->y_bpr[n1e]; X1.z=fem->z_bpr[n1e]; */
/*       X2.x=fem->x_bpr[n2e]; X2.y=fem->y_bpr[n2e]; X2.z=fem->z_bpr[n2e]; */
/*       X3.x=fem->x_bpr[n3e]; X3.y=fem->y_bpr[n3e]; X3.z=fem->z_bpr[n3e]; */
/*     } */


    A1=MINUS(X3,X2) , A2=MINUS(X1,X3), A3=MINUS(X2,X1);
    a1=MINUS(x3,x2) , a2=MINUS(x1,x3), a3=MINUS(x2,x1);
    //  N=UNIT(CROSS(A3,AMULT(-1.,A2)));
     N.x=fem->Nf_x[be]; N.y=fem->Nf_y[be]; N.z=fem->Nf_z[be];

    if(fem->belmtsedge[ec]==edge_n){ //check if it is on right edge
      if(fem->edgefrontnodesI[ec]==1){
	//I=1, J=4 , P=3
	e1=UNIT(A1);
	n1=CROSS(e1,N);
	h1=AMULT(DOT(n1,a3),n1);
	x4=PLUS(x1,AMULT(2.,h1));
	fem->p4x[be]=x4.x; fem->p4y[be]=x4.y; fem->p4z[be]=x4.z; 
      }else if(fem->edgefrontnodesI[ec]==2){
	//I=2, J=5 , P=1
  	e2=UNIT(A2);
	n2=CROSS(e2,N);
	h2=AMULT(DOT(n2,a1),n2);
	x5=PLUS(x2,AMULT(2.,h2));
	fem->p5x[be]=x5.x; fem->p5y[be]=x5.y; fem->p5z[be]=x5.z; 
      }else if(fem->edgefrontnodesI[ec]==3){
	//I=3, J=6 , P=2
	e3=UNIT(A3);
	n3=CROSS(e3,N);
	h3=AMULT(DOT(n3,a2),n3);
	x6=PLUS(x3,AMULT(2.,h3));
	fem->p6x[be]=x6.x; fem->p6y[be]=x6.y; fem->p6z[be]=x6.z; 
      }
    }
  }
 return(0);
}


//------------------------------------------------------------------------------------------------------------
PetscErrorCode ModifyGhostFix(PetscInt edge_n,FE *fem){
  
  
  PetscInt n1e,n2e,n3e,ec,be,i,p;

  Cmpnts N,X1,X2,X3,x1,x2,x3,x4,x5,x6,a1,a2,a3,A1,A2,A3,e1,e2,e3,n1,n2,n3;
  PetscReal E[3][3],Fadd[3],FJ[3],sum;

  PetscReal *FF,*FFJ;
  VecGetArray(fem->Fint, &FF);
  VecGetArray(fem->FJ, &FFJ);

  for (ec=0; ec<fem->n_ghosts; ec++) { //Go through ghost nodes
    be=fem->belmts[ec];
    n1e=fem->nv1[be];n2e=fem->nv2[be];n3e=fem->nv3[be];

    //currentlocation
    x1.x=fem->x_bp[n1e]; x1.y=fem->y_bp[n1e]; x1.z=fem->z_bp[n1e];
    x2.x=fem->x_bp[n2e]; x2.y=fem->y_bp[n2e]; x2.z=fem->z_bp[n2e];
    x3.x=fem->x_bp[n3e]; x3.y=fem->y_bp[n3e]; x3.z=fem->z_bp[n3e];
     //initial location
    //   if(!initrest){
       X1.x=fem->x_bp0[n1e]; X1.y=fem->y_bp0[n1e]; X1.z=fem->z_bp0[n1e];
       X2.x=fem->x_bp0[n2e]; X2.y=fem->y_bp0[n2e]; X2.z=fem->z_bp0[n2e];
       X3.x=fem->x_bp0[n3e]; X3.y=fem->y_bp0[n3e]; X3.z=fem->z_bp0[n3e];
  /*   }else{ */
/*        X1.x=fem->x_bpr[n1e]; X1.y=fem->y_bpr[n1e]; X1.z=fem->z_bpr[n1e]; */
/*        X2.x=fem->x_bpr[n2e]; X2.y=fem->y_bpr[n2e]; X2.z=fem->z_bpr[n2e]; */
/*        X3.x=fem->x_bpr[n3e]; X3.y=fem->y_bpr[n3e]; X3.z=fem->z_bpr[n3e]; */
/*      } */



    A1=MINUS(X3,X2) , A2=MINUS(X1,X3), A3=MINUS(X2,X1);
    a1=MINUS(x3,x2) , a2=MINUS(x1,x3), a3=MINUS(x2,x1);
    // N=UNIT(CROSS(A3,AMULT(-1.,A2)));
    N.x=fem->Nf_x[be]; N.y=fem->Nf_y[be]; N.z=fem->Nf_z[be];

    if(fem->belmtsedge[ec]==edge_n){ //check if it is on right edge
      if(fem->edgefrontnodesI[ec]==1){
	//I=1, J=4 , P=3
	e1=UNIT(A1);
	n1=CROSS(e1,N);

	E[0][0]=1-2*n1.x*n1.x; E[0][1]=-2*n1.x*n1.y;  E[0][2]=-2*n1.x*n1.z;
	E[1][0]=-2*n1.y*n1.x;  E[1][1]=1-2*n1.y*n1.y; E[1][2]=-2*n1.y*n1.z;
	E[2][0]=-2*n1.z*n1.x;  E[2][1]=-2*n1.z*n1.y;  E[2][2]=1-2*n1.z*n1.z;

	FJ[0]=FFJ[n1e*dof]; FJ[1]=FFJ[n1e*dof+1]; FJ[2]=FFJ[n1e*dof+2];

	for (i=0;i<3;i++){
	  sum=0.0;
	  for (p=0;p<3;p++){
	    sum+=E[p][i]*FJ[p];
	  }
	  Fadd[i]=sum;
	}

	FF[n1e*dof] +=Fadd[0]; 	FF[n1e*dof+1] +=Fadd[1];   FF[n1e*dof+2] +=Fadd[2];
	
      }else if(fem->edgefrontnodesI[ec]==2){
	//I=2, J=5 , P=1
  	e2=UNIT(A2);
	n2=CROSS(e2,N);
	
	E[0][0]=1-2*n2.x*n2.x; E[0][1]=-2*n2.x*n2.y;  E[0][2]=-2*n2.x*n2.z;
	E[1][0]=-2*n2.y*n2.x;  E[1][1]=1-2*n2.y*n2.y; E[1][2]=-2*n2.y*n2.z;
	E[2][0]=-2*n2.z*n2.x;  E[2][1]=-2*n2.z*n2.y;  E[2][2]=1-2*n2.z*n2.z;

	FJ[0]=FFJ[n2e*dof]; FJ[1]=FFJ[n2e*dof+1]; FJ[2]=FFJ[n2e*dof+2];

	for (i=0;i<3;i++){
	  sum=0.0;
	  for (p=0;p<3;p++){
	    sum+=E[p][i]*FJ[p];
	  }
	  Fadd[i]=sum;
	}

	FF[n2e*dof] +=Fadd[0]; 	FF[n2e*dof+1] +=Fadd[1];   FF[n2e*dof+2] +=Fadd[2];


      }else if(fem->edgefrontnodesI[ec]==3){
	//I=3, J=6 , P=2
	e3=UNIT(A3);
	n3=CROSS(e3,N);

	E[0][0]=1-2*n3.x*n3.x; E[0][1]=-2*n3.x*n3.y;  E[0][2]=-2*n3.x*n3.z;
	E[1][0]=-2*n3.y*n3.x;  E[1][1]=1-2*n3.y*n3.y; E[1][2]=-2*n3.y*n3.z;
	E[2][0]=-2*n3.z*n3.x;  E[2][1]=-2*n3.z*n3.y;  E[2][2]=1-2*n3.z*n3.z;

	FJ[0]=FFJ[n3e*dof]; FJ[1]=FFJ[n3e*dof+1]; FJ[2]=FFJ[n3e*dof+2];

	for (i=0;i<3;i++){
	  sum=0.0;
	  for (p=0;p<3;p++){
	    sum+=E[p][i]*FJ[p];
	  }
	  Fadd[i]=sum;
	}

	FF[n3e*dof] +=Fadd[0]; 	FF[n3e*dof+1] +=Fadd[1]; 	FF[n3e*dof+2] +=Fadd[2];
      }
    }
  }
  

  VecRestoreArray(fem->FJ, &FFJ); 
  VecRestoreArray(fem->Fint, &FF); 
 return(0);
}

//------------------------------------------------------------------------------------------------------------
PetscErrorCode  GhostFree(PetscInt edge_n, FE *fem) {
  
  PetscInt  n1e, n2e, n3e, ec, be, i;  
  Cmpnts    x1, x2, x3, x4, x5, x6, a1, a2, a3, e1, e2, e3, n1, n2, n3, h1, h2, h3, n;
  for (ec=0; ec<fem->n_ghosts; ec++) { //Go through ghost nodes
    be = fem->belmts[ec];
    n1e = fem->nv1[be];  n2e = fem->nv2[be];  n3e = fem->nv3[be];
    
    //currentlocation
    x1.x = fem->x_bp[n1e];  x1.y = fem->y_bp[n1e];  x1.z = fem->z_bp[n1e];
    x2.x = fem->x_bp[n2e];  x2.y = fem->y_bp[n2e];  x2.z = fem->z_bp[n2e];
    x3.x = fem->x_bp[n3e];  x3.y = fem->y_bp[n3e];  x3.z = fem->z_bp[n3e];
    
    a1 = MINUS(x3,x2);  a2 = MINUS(x1,x3);  a3 = MINUS(x2,x1);
    
    Cmpnts  nt;
    //  n=UNIT(CROSS(a3,AMULT(-1.,a2)));
    n.x = fem->nf_x[be];  n.y = fem->nf_y[be];  n.z = fem->nf_z[be];
    //  PetscPrintf(PETSC_COMM_SELF, "ec:%d n=(%le,%le,%le) nt=(%le,%le,%le)\n",be,n.x,n.y,n.z,nt.x,nt.y,nt.z);
    
    if(fem->belmtsedge[ec]==edge_n){ //check if it is on right edge
      if(fem->edgefrontnodesI[ec]==1){
	//I=1, J=4 , P=3
	e1 = UNIT(a1);
	n1 = CROSS(e1, n);
	h1 = AMULT(DOT(n1, a3), n1);
	x4 = PLUS(x1, AMULT(2., h1));
	fem->p4x[be] = x4.x;  fem->p4y[be] = x4.y;  fem->p4z[be] = x4.z; 
	
      }else if(fem->edgefrontnodesI[ec]==2){
	//I=2, J=5 , P=1
  	e2 = UNIT(a2);
	n2 = CROSS(e2, n);
	h2 = AMULT(DOT(n2, a1), n2);
	x5 = PLUS(x2, AMULT(2., h2));
	fem->p5x[be] = x5.x;  fem->p5y[be] = x5.y;  fem->p5z[be] = x5.z; 
      }else if(fem->edgefrontnodesI[ec]==3){
	//I=3, J=6 , P=2
	e3 = UNIT(a3);
	n3 = CROSS(e3, n);
	h3 = AMULT(DOT(n3, a2), n3);
	x6 = PLUS(x3,AMULT(2., h3));
	fem->p6x[be] = x6.x;  fem->p6y[be] = x6.y;  fem->p6z[be] = x6.z; 
	
      }
    }
    
  }

  return(0);
}

//------------------------------------------------------------------------------------------------------------
PetscErrorCode ModifyGhostFree(PetscInt edge_n,FE *fem) {
  
  PetscInt n1e,n2e,n3e,ec,be,i,j,m,p;
  
  Cmpnts x1,x2,x3,a1,a2,a3,e1,e2,e3,n,n1,n2,n3;
  PetscReal A,c,nIDaP,sum,sum1,sum2;
  PetscReal EI[3][3],EP[3][3],EF[3][3],Er1[3][3],Er2[3][3],Er3[3][3],FJ[3],Fadd[3];
  
  PetscReal NM1[3][3],NM2[3][3],nITPaP[3][3],eINM1[3][3],eINM2[3][3],CI[3][3],nCI[3][3],nITPnI[3][3];
  PetscReal IMnTPn[3][3],M1[3][3],M2[3][3];
  Cmpnts dx21,dx31;
  
  PetscReal *FF,*FFJ;
  VecGetArray(fem->Fint, &FF);
  VecGetArray(fem->FJ, &FFJ);
  for (ec=0; ec<fem->n_ghosts; ec++) { //Go through ghost nodes
    be=fem->belmts[ec];
    n1e=fem->nv1[be];n2e=fem->nv2[be];n3e=fem->nv3[be];
    
    // currentlocation
    x1.x=fem->x_bp[n1e]; x1.y=fem->y_bp[n1e]; x1.z=fem->z_bp[n1e];
    x2.x=fem->x_bp[n2e]; x2.y=fem->y_bp[n2e]; x2.z=fem->z_bp[n2e];
    x3.x=fem->x_bp[n3e]; x3.y=fem->y_bp[n3e]; x3.z=fem->z_bp[n3e];

    dx21=MINUS(x2,x1); dx31=MINUS(x3,x1); //dx21:g1 , dx31:g2
    a1=MINUS(x3,x2) , a2=MINUS(x1,x3), a3=MINUS(x2,x1);
    
    A=fem->dA[be]; 
    A=0.5*SIZE(CROSS(dx21,dx31)); 
    n.x=fem->nf_x[be]; n.y=fem->nf_y[be]; n.z=fem->nf_z[be];
    //  n=UNIT(CROSS(dx21,dx31));
    
    IMnTPn[0][0]=1.-n.x*n.x;  IMnTPn[0][1]=-n.x*n.y;      IMnTPn[0][2]=-n.x*n.z;
    IMnTPn[1][0]=-n.y*n.x;    IMnTPn[1][1]=1.-n.y*n.y;    IMnTPn[1][2]=-n.y*n.z;
    IMnTPn[2][0]=-n.z*n.x;    IMnTPn[2][1]=-n.z*n.y;      IMnTPn[2][2]=1.-n.z*n.z;  
    
    M1[0][0]=0.;        M1[0][1]=-dx21.z;    M1[0][2]=dx21.y;
    M1[1][0]=dx21.z;    M1[1][1]=0.;         M1[1][2]=-dx21.x;
    M1[2][0]=-dx21.y;   M1[2][1]=dx21.x;     M1[2][2]=0.;
    
    M2[0][0]=0.;        M2[0][1]=-dx31.z;  M2[0][2]=dx31.y;
    M2[1][0]=dx31.z;    M2[1][1]=0.;       M2[1][2]=-dx31.x;
    M2[2][0]=-dx31.y;   M2[2][1]=dx31.x;   M2[2][2]=0.;
    
    for (i=0;i<3;i++){
      for (j=0;j<3;j++){
	sum1=0.0; sum2=0.0;
	for (m=0;m<3;m++){
	  sum1+=IMnTPn[i][m]*M1[m][j]/(2.*A);
	  sum2+=IMnTPn[i][m]*M2[m][j]/(2.*A);
	}
	NM1[i][j]=sum1; NM2[i][j]=sum2;
      }
    }   
    
    if(fem->belmtsedge[ec]==edge_n){ //check if it is on right edge
      //--------------------------------------I=1,  P=3 , F=2 , J=4 ----------------------------------------------------------------------
      if(fem->edgefrontnodesI[ec]==1){
	//I=1, P=3 , F=2 , J=4 
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    EI[i][j]=0.0; EP[i][j]=0.0; EF[i][j]=0.0;
	  }
	}
	e1=UNIT(a1);
	n1=CROSS(e1,n);
	
	nIDaP=DOT(n1,a3); //nI.ap
	
	nITPaP[0][0]=n1.x*a3.x; nITPaP[0][1]=n1.x*a3.y; nITPaP[0][2]=n1.x*a3.z; //nI Tensor product ap
	nITPaP[1][0]=n1.y*a3.x; nITPaP[1][1]=n1.y*a3.y; nITPaP[1][2]=n1.y*a3.z;
	nITPaP[2][0]=n1.z*a3.x; nITPaP[2][1]=n1.z*a3.y; nITPaP[2][2]=n1.z*a3.z;
	
	eINM1[0][0]=e1.y*NM1[2][0]-e1.z*NM1[1][0]; eINM1[0][1]=e1.y*NM1[2][1]-e1.z*NM1[1][1];  eINM1[0][2]=e1.y*NM1[2][2]-e1.z*NM1[1][2];
	eINM1[1][0]=e1.z*NM1[0][0]-e1.x*NM1[2][0]; eINM1[1][1]=e1.z*NM1[0][1]-e1.x*NM1[2][1];  eINM1[1][2]=e1.z*NM1[0][2]-e1.x*NM1[2][2];
	eINM1[2][0]=e1.x*NM1[1][0]-e1.y*NM1[0][0]; eINM1[2][1]=e1.x*NM1[1][1]-e1.y*NM1[0][1];  eINM1[2][2]=e1.x*NM1[1][2]-e1.y*NM1[0][2];
	
	eINM2[0][0]=e1.y*NM2[2][0]-e1.z*NM2[1][0]; eINM2[0][1]=e1.y*NM2[2][1]-e1.z*NM2[1][1];  eINM2[0][2]=e1.y*NM2[2][2]-e1.z*NM2[1][2];
	eINM2[1][0]=e1.z*NM2[0][0]-e1.x*NM2[2][0]; eINM2[1][1]=e1.z*NM2[0][1]-e1.x*NM2[2][1];  eINM2[1][2]=e1.z*NM2[0][2]-e1.x*NM2[2][2];
	eINM2[2][0]=e1.x*NM2[1][0]-e1.y*NM2[0][0]; eINM2[2][1]=e1.x*NM2[1][1]-e1.y*NM2[0][1];  eINM2[2][2]=e1.x*NM2[1][2]-e1.y*NM2[0][2];
	
	nITPnI[0][0]=n1.x*n1.x;  nITPnI[0][1]=n1.x*n1.y;  nITPnI[0][2]=n1.x*n1.z;
	nITPnI[1][0]=n1.y*n1.x;  nITPnI[1][1]=n1.y*n1.y;  nITPnI[1][2]=n1.y*n1.z;
	nITPnI[2][0]=n1.z*n1.x;  nITPnI[2][1]=n1.z*n1.y;  nITPnI[2][2]=n1.z*n1.z;
	
	c=1./pow(SIZE(a1),2.); //function(eI,aI)
	CI[0][0]=c*(SIZE(a1)-e1.x*a1.x);   CI[0][1]=c*(-e1.x*a1.y);         CI[0][2]=c*(-e1.x*a1.z);
	CI[1][0]=c*(-e1.y*a1.x);           CI[1][1]=c*(SIZE(a1)-e1.y*a1.y); CI[1][2]=c*(-e1.y*a1.z);
	CI[2][0]=c*(-e1.z*a1.x);           CI[2][1]=c*(-e1.z*a1.y);         CI[2][2]=c*(SIZE(a1)-e1.z*a1.z);
	
	nCI[0][0]=n.z*CI[1][0]-n.y*CI[2][0]; nCI[0][1]=n.z*CI[1][1]-n.y*CI[2][1]; nCI[0][2]=n.z*CI[1][2]-n.y*CI[2][2];
	nCI[1][0]=n.x*CI[2][0]-n.z*CI[0][0]; nCI[1][1]=n.x*CI[2][1]-n.z*CI[0][1]; nCI[1][2]=n.x*CI[2][2]-n.z*CI[0][2];
	nCI[2][0]=n.y*CI[0][0]-n.x*CI[1][0]; nCI[2][1]=n.y*CI[0][1]-n.x*CI[1][1]; nCI[2][2]=n.y*CI[0][2]-n.x*CI[1][2];
	
	//EI
	EI[0][0]=1.-2.*nITPnI[0][0]; EI[0][1]=-2.*nITPnI[0][1];  EI[0][2]=-2.*nITPnI[0][2];
	EI[1][0]=-2.*nITPnI[1][0];  EI[1][1]=1.-2.*nITPnI[1][1]; EI[1][2]=-2.*nITPnI[1][2];
	EI[2][0]=-2.*nITPnI[2][0];  EI[2][1]=-2.*nITPnI[2][1];  EI[2][2]=1.-2.*nITPnI[2][2];
	//EP
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=2*nITPaP[i][m]*nCI[m][j];
	    }
	    EP[i][j]=2*nIDaP*nCI[i][j]+sum;
	  }
	}
	//EF
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=-2*nITPaP[i][m]*nCI[m][j];
	    }
	    EF[i][j]=-2*nIDaP*nCI[i][j]+sum+2*nITPnI[i][j];
	  }
	}
	
	//Er1
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=-2*nITPaP[i][m]*eINM1[m][j]+2*nITPaP[i][m]*eINM2[m][j];
	    }
	    Er1[i][j]=sum-2*nIDaP*eINM1[i][j]+2*nIDaP*eINM2[i][j];
	  }
	}
       	//Er2
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=-2*nITPaP[i][m]*eINM2[m][j];
	    }
	    Er2[i][j]=sum-2*nIDaP*eINM2[i][j];
	  }
	}
	//Er3
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=2*nITPaP[i][m]*eINM1[m][j];
	    }
	    Er3[i][j]=sum+2*nIDaP*eINM1[i][j];
	  }
	}
	//I=1, P=3 , F=2
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    EI[i][j]+=Er1[i][j];
	    EF[i][j]+=Er2[i][j];
	    EP[i][j]+=Er3[i][j];
	  }
	}
	
	FJ[0]=FFJ[n1e*dof]; FJ[1]=FFJ[n1e*dof+1]; FJ[2]=FFJ[n1e*dof+2];
	
	//I=1
	for (i=0;i<3;i++){
	  sum=0.0;
	  for (p=0;p<3;p++){
	    sum+=EI[p][i]*FJ[p];
	  }
	  Fadd[i]=sum;
	}
	
	FF[n1e*dof] +=Fadd[0]; 	FF[n1e*dof+1] +=Fadd[1]; FF[n1e*dof+2] +=Fadd[2];
	
	//F=2
	for (i=0;i<3;i++){
	  sum=0.0;
	  for (p=0;p<3;p++){
	    sum+=EF[p][i]*FJ[p];
	  }
	  Fadd[i]=sum;
	}
       	FF[n2e*dof] +=Fadd[0]; 	FF[n2e*dof+1] +=Fadd[1]; FF[n2e*dof+2] +=Fadd[2];
		
	//P=3
	for (i=0;i<3;i++){
	  sum=0.0;
	  for (p=0;p<3;p++){
	    sum+=EP[p][i]*FJ[p];
	  }
	  Fadd[i]=sum;
	}
       	FF[n3e*dof] +=Fadd[0]; 	FF[n3e*dof+1] +=Fadd[1]; FF[n3e*dof+2] +=Fadd[2];
	
	//--------------------------------------I=2, P=1 , F=3 , J=5 ----------------------------------------------------------------------	
      }else if(fem->edgefrontnodesI[ec]==2){
	//I=2, P=1 , F=3 ,J=5 
  	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    EI[i][j]=0.0; EP[i][j]=0.0; EF[i][j]=0.0;
	  }
	}
	e2=UNIT(a2);
	n2=CROSS(e2,n);
	
	nIDaP=DOT(n2,a1); //nI.ap
	
	nITPaP[0][0]=n2.x*a1.x; nITPaP[0][1]=n2.x*a1.y; nITPaP[0][2]=n2.x*a1.z; //nI Tensor product ap
	nITPaP[1][0]=n2.y*a1.x; nITPaP[1][1]=n2.y*a1.y; nITPaP[1][2]=n2.y*a1.z;
	nITPaP[2][0]=n2.z*a1.x; nITPaP[2][1]=n2.z*a1.y; nITPaP[2][2]=n2.z*a1.z;
	
	eINM1[0][0]=e2.y*NM1[2][0]-e2.z*NM1[1][0]; eINM1[0][1]=e2.y*NM1[2][1]-e2.z*NM1[1][1];  eINM1[0][2]=e2.y*NM1[2][2]-e2.z*NM1[1][2];
	eINM1[1][0]=e2.z*NM1[0][0]-e2.x*NM1[2][0]; eINM1[1][1]=e2.z*NM1[0][1]-e2.x*NM1[2][1];  eINM1[1][2]=e2.z*NM1[0][2]-e2.x*NM1[2][2];
	eINM1[2][0]=e2.x*NM1[1][0]-e2.y*NM1[0][0]; eINM1[2][1]=e2.x*NM1[1][1]-e2.y*NM1[0][1];  eINM1[2][2]=e2.x*NM1[1][2]-e2.y*NM1[0][2];
	
	eINM2[0][0]=e2.y*NM2[2][0]-e2.z*NM2[1][0]; eINM2[0][1]=e2.y*NM2[2][1]-e2.z*NM2[1][1];  eINM2[0][2]=e2.y*NM2[2][2]-e2.z*NM2[1][2];
	eINM2[1][0]=e2.z*NM2[0][0]-e2.x*NM2[2][0]; eINM2[1][1]=e2.z*NM2[0][1]-e2.x*NM2[2][1];  eINM2[1][2]=e2.z*NM2[0][2]-e2.x*NM2[2][2];
	eINM2[2][0]=e2.x*NM2[1][0]-e2.y*NM2[0][0]; eINM2[2][1]=e2.x*NM2[1][1]-e2.y*NM2[0][1];  eINM2[2][2]=e2.x*NM2[1][2]-e2.y*NM2[0][2];
	
	nITPnI[0][0]=n2.x*n2.x;  nITPnI[0][1]=n2.x*n2.y;  nITPnI[0][2]=n2.x*n2.z;
	nITPnI[1][0]=n2.y*n2.x;  nITPnI[1][1]=n2.y*n2.y;  nITPnI[1][2]=n2.y*n2.z;
	nITPnI[2][0]=n2.z*n2.x;  nITPnI[2][1]=n2.z*n2.y;  nITPnI[2][2]=n2.z*n2.z;
	
	c=1./pow(SIZE(a2),2); //function(eI,aI)
	CI[0][0]=c*(SIZE(a2)-e2.x*a2.x);   CI[0][1]=c*(-e2.x*a2.y);         CI[0][2]=c*(-e2.x*a2.z);
	CI[1][0]=c*(-e2.y*a2.x);           CI[1][1]=c*(SIZE(a2)-e2.y*a2.y); CI[1][2]=c*(-e2.y*a2.z);
	CI[2][0]=c*(-e2.z*a2.x);           CI[2][1]=c*(-e2.z*a2.y);         CI[2][2]=c*(SIZE(a2)-e2.z*a2.z);
	
	nCI[0][0]=n.z*CI[1][0]-n.y*CI[2][0]; nCI[0][1]=n.z*CI[1][1]-n.y*CI[2][1]; nCI[0][2]=n.z*CI[1][2]-n.y*CI[2][2];
	nCI[1][0]=n.x*CI[2][0]-n.z*CI[0][0]; nCI[1][1]=n.x*CI[2][1]-n.z*CI[0][1]; nCI[1][2]=n.x*CI[2][2]-n.z*CI[0][2];
	nCI[2][0]=n.y*CI[0][0]-n.x*CI[1][0]; nCI[2][1]=n.y*CI[0][1]-n.x*CI[1][1]; nCI[2][2]=n.y*CI[0][2]-n.x*CI[1][2];
	
	//EI
	EI[0][0]=1-2*nITPnI[0][0]; EI[0][1]=-2*nITPnI[0][1];  EI[0][2]=-2*nITPnI[0][2];
	EI[1][0]=-2*nITPnI[1][0];  EI[1][1]=1-2*nITPnI[1][1]; EI[1][2]=-2*nITPnI[1][2];
	EI[2][0]=-2*nITPnI[2][0];  EI[2][1]=-2*nITPnI[2][1];  EI[2][2]=1-2*nITPnI[2][2];
	//EP
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=2*nITPaP[i][m]*nCI[m][j];
	    }
	    EP[i][j]=2*nIDaP*nCI[i][j]+sum;
	  }
	}
	//EF
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=-2*nITPaP[i][m]*nCI[m][j];
	    }
	    EF[i][j]=-2*nIDaP*nCI[i][j]+sum+2*nITPnI[i][j];
	  }
	}
	
	//Er1
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=-2*nITPaP[i][m]*eINM1[m][j]+2*nITPaP[i][m]*eINM2[m][j];
	    }
	    Er1[i][j]=sum-2*nIDaP*eINM1[i][j]+2*nIDaP*eINM2[i][j];
	  }
	}
       	//Er2
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=-2*nITPaP[i][m]*eINM2[m][j];
	    }
	    Er2[i][j]=sum-2*nIDaP*eINM2[i][j];
	  }
	}
	//Er3
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=2*nITPaP[i][m]*eINM1[m][j];
	    }
	    Er3[i][j]=sum+2*nIDaP*eINM1[i][j];
	  }
	}
	//I=2,  P=1 , F=3
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    EI[i][j]+=Er2[i][j];
	    EF[i][j]+=Er3[i][j];
	    EP[i][j]+=Er1[i][j];
	  }
	}
	
	FJ[0]=FFJ[n2e*dof]; FJ[1]=FFJ[n2e*dof+1]; FJ[2]=FFJ[n2e*dof+2];
	//I=2
	for (i=0;i<3;i++){
	  sum=0.0;
	  for (p=0;p<3;p++){
	    sum+=EI[p][i]*FJ[p];
	  }
	  Fadd[i]=sum;
	}
	FF[n2e*dof] +=Fadd[0]; 	FF[n2e*dof+1] +=Fadd[1]; FF[n2e*dof+2] +=Fadd[2];
	
	//F=3
	for (i=0;i<3;i++){
	  sum=0.0;
	  for (p=0;p<3;p++){
	    sum+=EF[p][i]*FJ[p];
	  }
	  Fadd[i]=sum;
	}
       	FF[n3e*dof] +=Fadd[0]; 	FF[n3e*dof+1] +=Fadd[1]; FF[n3e*dof+2] +=Fadd[2];
	
	//P=1
	for (i=0;i<3;i++){
	  sum=0.0;
	  for (p=0;p<3;p++){
	    sum+=EP[p][i]*FJ[p];
	  }
	  Fadd[i]=sum;
	}
       	FF[n1e*dof] +=Fadd[0]; 	FF[n1e*dof+1] +=Fadd[1]; FF[n1e*dof+2] +=Fadd[2];
	
	//--------------------------------------I=3, P=2 , F=1 , J=6 ----------------------------------------------------------------------	
      }else if(fem->edgefrontnodesI[ec]==3){
	//I=3,  P=2 , F=1 , J=6 
	e3=UNIT(a3);
	n3=CROSS(e3,n);
	
	nIDaP=DOT(n3,a2); //nI.ap
	
	nITPaP[0][0]=n3.x*a2.x; nITPaP[0][1]=n3.x*a2.y; nITPaP[0][2]=n3.x*a2.z; //nI Tensor product ap
	nITPaP[1][0]=n3.y*a2.x; nITPaP[1][1]=n3.y*a2.y; nITPaP[1][2]=n3.y*a2.z;
	nITPaP[2][0]=n3.z*a2.x; nITPaP[2][1]=n3.z*a2.y; nITPaP[2][2]=n3.z*a2.z;
	
	eINM1[0][0]=e3.y*NM1[2][0]-e3.z*NM1[1][0]; eINM1[0][1]=e3.y*NM1[2][1]-e3.z*NM1[1][1];  eINM1[0][2]=e3.y*NM1[2][2]-e3.z*NM1[1][2];
	eINM1[1][0]=e3.z*NM1[0][0]-e3.x*NM1[2][0]; eINM1[1][1]=e3.z*NM1[0][1]-e3.x*NM1[2][1];  eINM1[1][2]=e3.z*NM1[0][2]-e3.x*NM1[2][2];
	eINM1[2][0]=e3.x*NM1[1][0]-e3.y*NM1[0][0]; eINM1[2][1]=e3.x*NM1[1][1]-e3.y*NM1[0][1];  eINM1[2][2]=e3.x*NM1[1][2]-e3.y*NM1[0][2];
	
	eINM2[0][0]=e3.y*NM2[2][0]-e3.z*NM2[1][0]; eINM2[0][1]=e3.y*NM2[2][1]-e3.z*NM2[1][1];  eINM2[0][2]=e3.y*NM2[2][2]-e3.z*NM2[1][2];
	eINM2[1][0]=e3.z*NM2[0][0]-e3.x*NM2[2][0]; eINM2[1][1]=e3.z*NM2[0][1]-e3.x*NM2[2][1];  eINM2[1][2]=e3.z*NM2[0][2]-e3.x*NM2[2][2];
	eINM2[2][0]=e3.x*NM2[1][0]-e3.y*NM2[0][0]; eINM2[2][1]=e3.x*NM2[1][1]-e3.y*NM2[0][1];  eINM2[2][2]=e3.x*NM2[1][2]-e3.y*NM2[0][2];
	
	nITPnI[0][0]=n3.x*n3.x;  nITPnI[0][1]=n3.x*n3.y;  nITPnI[0][2]=n3.x*n3.z;
	nITPnI[1][0]=n3.y*n3.x;  nITPnI[1][1]=n3.y*n3.y;  nITPnI[1][2]=n3.y*n3.z;
	nITPnI[2][0]=n3.z*n3.x;  nITPnI[2][1]=n3.z*n3.y;  nITPnI[2][2]=n3.z*n3.z;
	
	c=1./pow(SIZE(a3),2); //function(eI,aI)
	CI[0][0]=c*(SIZE(a3)-e3.x*a3.x);   CI[0][1]=c*(-e3.x*a3.y);         CI[0][2]=c*(-e3.x*a3.z);
	CI[1][0]=c*(-e3.y*a3.x);           CI[1][1]=c*(SIZE(a3)-e3.y*a3.y); CI[1][2]=c*(-e3.y*a3.z);
	CI[2][0]=c*(-e3.z*a3.x);           CI[2][1]=c*(-e3.z*a3.y);         CI[2][2]=c*(SIZE(a3)-e3.z*a3.z);
	
	nCI[0][0]=n.z*CI[1][0]-n.y*CI[2][0]; nCI[0][1]=n.z*CI[1][1]-n.y*CI[2][1]; nCI[0][2]=n.z*CI[1][2]-n.y*CI[2][2];
	nCI[1][0]=n.x*CI[2][0]-n.z*CI[0][0]; nCI[1][1]=n.x*CI[2][1]-n.z*CI[0][1]; nCI[1][2]=n.x*CI[2][2]-n.z*CI[0][2];
	nCI[2][0]=n.y*CI[0][0]-n.x*CI[1][0]; nCI[2][1]=n.y*CI[0][1]-n.x*CI[1][1]; nCI[2][2]=n.y*CI[0][2]-n.x*CI[1][2];
	
	//EI
	EI[0][0]=1-2*nITPnI[0][0]; EI[0][1]=-2*nITPnI[0][1];  EI[0][2]=-2*nITPnI[0][2];
	EI[1][0]=-2*nITPnI[1][0];  EI[1][1]=1-2*nITPnI[1][1]; EI[1][2]=-2*nITPnI[1][2];
	EI[2][0]=-2*nITPnI[2][0];  EI[2][1]=-2*nITPnI[2][1];  EI[2][2]=1-2*nITPnI[2][2];
	//EP
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=2*nITPaP[i][m]*nCI[m][j];
	    }
	    EP[i][j]=2*nIDaP*nCI[i][j]+sum;
	  }
	}
	//EF
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=-2*nITPaP[i][m]*nCI[m][j];
	    }
	    EF[i][j]=-2*nIDaP*nCI[i][j]+sum+2*nITPnI[i][j];
	  }
	}
	
	//Er1
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=-2*nITPaP[i][m]*eINM1[m][j]+2*nITPaP[i][m]*eINM2[m][j];
	    }
	    Er1[i][j]=sum-2*nIDaP*eINM1[i][j]+2*nIDaP*eINM2[i][j];
	  }
	}
       	//Er2
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=-2*nITPaP[i][m]*eINM2[m][j];
	    }
	    Er2[i][j]=sum-2*nIDaP*eINM2[i][j];
	  }
	}
	//Er3
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    sum=0.0; 
	    for (m=0;m<3;m++){
	      sum+=2*nITPaP[i][m]*eINM1[m][j];
	    }
	    Er3[i][j]=sum+2*nIDaP*eINM1[i][j];
	  }
	}
	//I=3, P=2 , F=1
	for (i=0;i<3;i++){
	  for (j=0;j<3;j++){
	    EI[i][j]+=Er3[i][j];
	    EF[i][j]+=Er1[i][j];
	    EP[i][j]+=Er2[i][j];
	  }
	}
	
	FJ[0]=FFJ[n3e*dof]; FJ[1]=FFJ[n3e*dof+1]; FJ[2]=FFJ[n3e*dof+2];
	//I=3
	for (i=0;i<3;i++){
	  sum=0.0;
	  for (p=0;p<3;p++){
	    sum+=EI[p][i]*FJ[p];
	  }
	  Fadd[i]=sum;
	}
	FF[n3e*dof] +=Fadd[0]; 	FF[n3e*dof+1] +=Fadd[1]; FF[n3e*dof+2] +=Fadd[2];
	
	//F=1
	for (i=0;i<3;i++){
	  sum=0.0;
	  for (p=0;p<3;p++){
	    sum+=EF[p][i]*FJ[p];
	  }
	  Fadd[i]=sum;
	}
       	FF[n1e*dof] +=Fadd[0]; 	FF[n1e*dof+1] +=Fadd[1]; FF[n1e*dof+2] +=Fadd[2];
	
	//P=2
	for (i=0;i<3;i++){
	  sum=0.0;
	  for (p=0;p<3;p++){
	    sum+=EP[p][i]*FJ[p];
	  }
	  Fadd[i]=sum;
	}
       	FF[n2e*dof] +=Fadd[0]; 	FF[n2e*dof+1] +=Fadd[1]; FF[n2e*dof+2] +=Fadd[2];
	
      }//else if

    }
  }
  
  VecRestoreArray(fem->FJ, &FFJ);
  VecRestoreArray(fem->Fint, &FF);
 
  return(0);
}
