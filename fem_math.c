#include "variables.h"


Cmpnts  CROSS(Cmpnts v1, Cmpnts v2) {

  // output = v1 x v2
  Cmpnts  v1cv2;
  v1cv2.x  = v1.y * v2.z - v2.y * v1.z;
  v1cv2.y  =-v1.x * v2.z + v2.x * v1.z;
  v1cv2.z  = v1.x * v2.y - v2.x * v1.y;
  
  return(v1cv2);
}

//--------------------------------------------------------------------------------------- 
PetscReal DOT(Cmpnts v1, Cmpnts v2)
{
  // output = v1.v2
  return(v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

//--------------------------------------------------------------------------------------- 
Cmpnts  UNIT(Cmpnts v1) {

  // output = normal
  Cmpnts  v4;  PetscReal  norm;
  norm = PetscSqrtReal(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z); 
  v4.x = v1.x/norm;  v4.y = v1.y/norm;  v4.z = v1.z/norm; 
  return(v4);
}

//--------------------------------------------------------------------------------------- 
PetscReal SIZE(Cmpnts v1) {

  // output = size
  return(PetscSqrtReal(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z));
}

//---------------------------------------------------------------------------------------
Cmpnts AMULT(PetscReal alpha,Cmpnts v1)
{
  Cmpnts v4;
  v4.x  =alpha* v1.x;
  v4.y  =alpha* v1.y;
  v4.z  =alpha* v1.z;
 
  return(v4);
}
//--------------------------------------------------------------------------------------- 
PetscErrorCode INV(PetscReal T[3][3],PetscReal _Tinv[3][3]){

  PetscReal T11,T12,T13,T21,T22,T23,T31,T32,T33,det,c;
  T11=T[0][0]; T12=T[0][1]; T13=T[0][2];    
  T21=T[1][0]; T22=T[1][1]; T23=T[1][2]; 
  T31=T[2][0]; T32=T[2][1]; T33=T[2][2]; 
       
  det=T11*(T33*T22-T32*T23)-T21*(T33*T12-T32*T13)+T31*(T23*T12-T22*T13);
  c=1/det;
  _Tinv[0][0]=c*(T33*T22-T32*T23);   _Tinv[0][1]=-c*(T33*T12-T32*T13);  _Tinv[0][2]=c*(T23*T12-T22*T13);
  _Tinv[1][0]=-c*(T33*T21-T31*T23);  _Tinv[1][1]=c*(T33*T11-T31*T13);   _Tinv[1][2]=-c*(T23*T11-T21*T13);
  _Tinv[2][0]=c*(T32*T21-T31*T22);   _Tinv[2][1]=-c*(T32*T11-T31*T12);  _Tinv[2][2]=c*(T22*T11-T21*T12);
}
//--------------------------------------------------------------------------------------- 
PetscErrorCode MATMULT(PetscReal A[][2],PetscReal B[][2], PetscReal C[][2])
{
  PetscInt i,j,k;
  for (i=0;i<2;i++){
    for (j=0;j<2;j++){
      C[i][j]=0;
      for (k=0;k<2;k++){
	C[i][j]+=A[i][k]*B[k][j];
      }
    }
  }
  return(0);
}
//--------------------------------------------------------------------------------------- 
PetscErrorCode TRANS(PetscReal A[][2], PetscReal AT[][2])
{
  // no change for diagonal
  AT[0][0]=A[0][0];
  AT[1][1]=A[1][1];
  // change the diagonal
  AT[0][1]=A[1][0];
  AT[1][0]=A[0][1];
  return(0);
}
  
//--------------------------------------------------------------------------------------- 
PetscReal SIGN(PetscReal a) {

  if (a>0) return(1.);
  if (a<0) return(-1.);

  return(0.);
}
