#include "variables.h"

extern PetscReal  E, mu, rho, h0, dampfactor;
extern PetscInt   dof, ConstitutiveLawNonLinear, ti;

extern Cmpnts  PLUS(Cmpnts v1, Cmpnts v2);
extern Cmpnts  MINUS(Cmpnts v1, Cmpnts v2);
extern Cmpnts  CROSS(Cmpnts v1, Cmpnts v2);
extern PetscReal  DOT(Cmpnts v1, Cmpnts v2);
extern Cmpnts  UNIT(Cmpnts v1);
extern PetscReal  SIZE(Cmpnts v1);
extern PetscErrorCode  INV(PetscReal T[3][3], PetscReal _Tinv[3][3]);
extern PetscErrorCode  MATMULT(PetscReal A[][2], PetscReal B[][2], PetscReal C[][2]);
extern Cmpnts  AMULT(PetscReal alpha, Cmpnts v1); 
extern PetscErrorCode TRANS(PetscReal A[3][3], PetscReal _AT[3][3]);
extern PetscReal  SIGN(PetscReal a);

 
PetscErrorCode  InitMaterial(FE *fem) {

  PetscInt  ec, i, j;
  
  // For LV
  Cmpnts  Axis, N, Midd, tcyl;  //Axis:Cylinder axis , Midd: middle of element , tcyl:tangent vector to cylinder
  PetscReal      pi=3.14159 , Teta; // Teta: angle of elemet in respect to the y=0 in xy surface
  PetscReal      R[3][3], nfibR[3], ttcyl[3], angle ;  //angle:angle of fiber direction respect to free edge
  PetscInt       n1e, n2e, n3e;
  PetscReal      e=0.02;
  
  Axis.x = 0.;  Axis.y = 0.;  Axis.z = 1;
  
  for (ec=0; ec<fem->n_elmt; ec++) {
    n1e = fem->nv1[ec];  n2e = fem->nv2[ec];  n3e = fem->nv3[ec];
    N.x = fem->nf_x[ec];  N.y = fem->nf_y[ec];  N.z = fem->nf_z[ec];
    //compute the middle of element in xy plane
    Midd.x = (fem->x_bp0[n1e] + fem->x_bp0[n2e] + fem->x_bp0[n3e])/3.;
    Midd.y = (fem->y_bp0[n1e] + fem->y_bp0[n2e] + fem->y_bp0[n3e])/3.;    
    
    Teta = atan(PetscAbsReal(Midd.y/Midd.x));
    if(Midd.x<0. && Midd.y>0.) {Teta = pi - Teta;}
    if(Midd.x<0. && Midd.y<0.) {Teta = pi + Teta;}
    if(Midd.x>0. && Midd.y<0.) {Teta = 2.*pi - Teta;}
    
    if((Teta>0. && Teta<pi/3.) || (Teta>2.*pi/3. && Teta<pi) || (Teta>4.*pi/3. && Teta<5.*pi/3.) ) {angle = pi - (45./180.)*pi;}
    if((Teta>pi/3. && Teta<2*pi/3.) || (Teta>pi && Teta<4.*pi/3.) || (Teta>5.*pi/3. && Teta<2.*pi) ) {angle = (45./180.)*pi;}
    //if((Teta>-e && Teta<e) || (Teta>2.*pi/3. - e && Teta<2.*pi/3. + e) || (Teta>4.*pi/3. - e && Teta<4.*pi/3. + e) ) {angle = -pi/2.;}

    //angle = pi/4.; //constant fiber direction
    //  FormRotationMatrixAround Normal to element (N) in current mesh it is outward
    R[0][0] = cos(angle) + N.x*N.x*(1 - cos(angle));
    R[0][1] = N.x*N.y*(1 - cos(angle)) - N.z*sin(angle);
    R[0][2] = N.x*N.z*(1 - cos(angle)) + N.y*sin(angle);
    
    R[1][0] = N.x*N.y*(1 - cos(angle)) + N.z*sin(angle);
    R[1][1] = cos(angle) + N.y*N.y*(1 - cos(angle));
    R[1][2] = N.y*N.z*(1 - cos(angle)) - N.x*sin(angle);
    
    R[2][0] = N.z*N.x*(1 - cos(angle)) - N.y*sin(angle);
    R[2][1] = N.z*N.y*(1 - cos(angle)) + N.x*sin(angle);
    R[2][2] = cos(angle) + N.z*N.z*(1 - cos(angle));
    
    tcyl = CROSS(N, Axis);
    ttcyl[0] = tcyl.x;  ttcyl[1] = tcyl.y;  ttcyl[2] = tcyl.z;
   
    //Rotate Vector
    for (i=0; i<3; i++){
      nfibR[i] = 0.;
      for (j=0; j<3; j++){
	nfibR[i] += R[i][j]*ttcyl[j];
      }
    }
    
    fem->n_fib[ec].x = nfibR[0];
    fem->n_fib[ec].y = nfibR[1];
    fem->n_fib[ec].z = nfibR[2]; 
    //for biaxial tests
    /* fem->n_fib[ec].x = 0.707; */
    /* fem->n_fib[ec].y = 0.707; */
    /* fem->n_fib[ec].z = 0.0; */
  }
 
  return(0);
}

//---------------------------------------------------------------------------------------  
PetscErrorCode Mass(FE *fem, PetscInt ec, PetscReal _M[9]) {

  PetscReal A0; PetscInt i;
  A0=fem->dA0[ec];
  
  for (i=0; i<9; i++){
    _M[i]=rho*h0*A0/3;
  }      
  return(0);
}   
  
//---------------------------------------------------------------------------------------  
PetscErrorCode Damp(PetscReal M[9],PetscReal _C[9]) {

  PetscInt i;
  for (i=0; i<9; i++){
    _C[i]=dampfactor; 
  }

  return(0);     
} 

//---------------------------------------------------------------------------------------  
PetscErrorCode MassDamp(FE *fem) {

  PetscReal  A0, M, C;
  PetscReal  *DDissip, *MMass;
  PetscInt   n1e, n2e, n3e, ec;
 
  VecGetArray(fem->Mass, &MMass);
  VecGetArray(fem->Dissip, &DDissip);
  
  for (ec=0; ec<fem->n_elmt + 2*fem->n_ghosts; ec++) {
    A0 = fem->dA0[ec];
    n1e = fem->nv1[ec];  n2e = fem->nv2[ec];  n3e = fem->nv3[ec];

    M = rho*h0*A0/3.;
      
    MMass[n1e*dof] += M;
    MMass[n1e*dof+1] += M;
    MMass[n1e*dof+2] += M;
    
    MMass[n2e*dof] += M;
    MMass[n2e*dof+1] += M;
    MMass[n2e*dof+2] += M;

    MMass[n3e*dof] += M;
    MMass[n3e*dof+1] += M;
    MMass[n3e*dof+2] += M;
 
    //C = dampfactor*M;
    C = dampfactor;

    DDissip[n1e*dof] += C;
    DDissip[n1e*dof+1] += C;
    DDissip[n1e*dof+2] += C;
    
    DDissip[n2e*dof] += C;
    DDissip[n2e*dof+1] += C;
    DDissip[n2e*dof+2] += C;

    DDissip[n3e*dof] += C;
    DDissip[n3e*dof+1] += C;
    DDissip[n3e*dof+2] += C;

  }

  VecRestoreArray(fem->Mass, &MMass);
  VecRestoreArray(fem->Dissip, &DDissip);
  
  return(0);     
}  
       
//--------------------------------------------------------------------------------------- 
PetscErrorCode StressLinear(PetscInt ec, Cmpnts X1, Cmpnts X2, Cmpnts X3, PetscReal Strain[3], PetscReal _S[3], PetscInt method, FE *fem) {
 
  PetscReal  Q[3][3], Dloc[3][3], Strainloc[3], c, sum;
  PetscInt   i, j, m, n;
  PetscReal  q[3][3], xz, xe, yz, ye;
  Cmpnts     dX21, dX31, V1, V2, N, diff1, diff2;

  if (method==0) {
    dX21 = MINUS(X2, X1);  dX31 = MINUS(X3, X1);  //dX21:G1 , dX31:G2
    N.x = fem->Nf_x[ec];  N.y = fem->Nf_y[ec];  N.z = fem->Nf_z[ec];
  } else if (method==1) { 
    dX21.x = fem->G1[ec*dof];  dX21.y = fem->G1[ec*dof+1];  dX21.z = fem->G1[ec*dof+2];
    dX31.x = fem->G2[ec*dof];  dX31.y = fem->G2[ec*dof+1];  dX31.z = fem->G2[ec*dof+2];
    N = UNIT(CROSS(dX21, dX31));
  }  
  V1 = UNIT(dX21);
  V2 = CROSS(N, V1);

  xz = DOT(V1, dX21);  xe = DOT(V1, dX31);
  yz = DOT(V2, dX21);  ye = DOT(V2, dX31);

  q[0][0] = xz*xz;  q[0][1] = xe*xe;  q[0][2] = 2*xz*xe;
  q[1][0] = yz*yz;  q[1][1] = ye*ye;  q[1][2] = 2*yz*ye;
  q[2][0] = xz*yz;  q[2][1] = xe*ye;  q[2][2] = xz*ye + xe*yz;

  INV(q, Q);

  c = E/(1. - pow(mu, 2.));
  Dloc[0][0] = c;     Dloc[0][1] = c*mu;  Dloc[0][2] = 0.;
  Dloc[1][0] = c*mu;  Dloc[1][1] = c;     Dloc[1][2] = 0.;
  Dloc[2][0] = 0.;    Dloc[2][1] = 0.;    Dloc[2][2] = c*(1. - mu)/2.;

  //Form local strain
  for (i=0; i<3; i++){
    sum = 0.0;
    for (n=0; n<3; n++){
      sum += Q[n][i]*Strain[n];
    }
    Strainloc[i] = sum;
  }

  for (i=0; i<3; i++){
    sum = 0.0;
    for (m=0; m<3; m++){
      for (n=0; n<3; n++){
  	sum += Q[i][m]*Dloc[m][n]*Strainloc[n];
      }
    }
    _S[i] = sum;
  }

  return(0);
}

//---------------------------------------------------------------------------------------  
PetscErrorCode CalcCurvStressStrainxyz(PetscInt ec, PetscReal k[3], PetscReal strain[3], PetscReal stress[3], PetscInt method, PetscInt mb, FE *fem) {

  Cmpnts     X1, X2, X3, dX21, dX31;
  PetscInt   i, j, m, n1e, n2e, n3e;
  PetscReal  sum, cart[3];
  PetscReal  q[3][3], Q[3][3], xz, xe, yz, ye;
  
  if(method==0) {
    n1e = fem->nv1[ec];  n2e = fem->nv2[ec];  n3e = fem->nv3[ec];
    X1.x = fem->x_bp0[n1e];  X1.y = fem->y_bp0[n1e];  X1.z = fem->z_bp0[n1e];
    X2.x = fem->x_bp0[n2e];  X2.y = fem->y_bp0[n2e];  X2.z = fem->z_bp0[n2e];
    X3.x = fem->x_bp0[n3e];  X3.y = fem->y_bp0[n3e];  X3.z = fem->z_bp0[n3e];
    dX21 = MINUS(X2, X1);
    dX31 = MINUS(X3, X1);

  } else if(method==1) {
    dX21.x = fem->G1[ec*dof];  dX21.y = fem->G1[ec*dof+1];  dX21.z = fem->G1[ec*dof+2];
    dX31.x = fem->G2[ec*dof];  dX31.y = fem->G2[ec*dof+1];  dX31.z = fem->G2[ec*dof+2];
  }

  xz = dX21.x;  xe = dX31.x;
  yz = dX21.y;  ye = dX31.y;
  
  q[0][0] = xz*xz;  q[0][1] = xe*xe;  q[0][2] = 2*xz*xe; //Transformation based on stress
  q[1][0] = yz*yz;  q[1][1] = ye*ye;  q[1][2] = 2*yz*ye;  
  q[2][0] = xz*yz;  q[2][1] = xe*ye;  q[2][2] = xz*ye + xe*yz; 
  INV(q, Q);

  //curvilinear to global curvature
  for (i=0; i<3; i++) {
    sum = 0.;
    for (j=0; j<3; j++) {
      sum += Q[j][i]*k[j];
    }
    cart[i] = sum; 
  }

  if (mb==1) {
    fem->kve[ec*dof] = cart[0];   fem->kve[ec*dof+1] = cart[1];   fem->kve[ec*dof+2] = cart[2];   
  }	

  //curvilinear to global strain
  for (i=0; i<3; i++) {
    sum = 0.;
    for (j=0; j<3; j++) {
      sum += Q[j][i]*strain[j]; 
    }
    cart[i] = sum;
  }

  if (mb==0) {
    fem->StrainM[ec*(dof+2)] = cart[0];   fem->StrainM[ec*(dof+2)+1] = cart[1];   fem->StrainM[ec*(dof+2)+2] = cart[2];  
    fem->StrainM[ec*(dof+2)+3] = (cart[0] + cart[1])/2. + sqrt(pow((cart[0] - cart[1])/2., 2) + pow(cart[2]/2., 2)); //in-plane max Principal Strain 
    fem->StrainM[ec*(dof+2)+4] = (cart[0] + cart[1])/2. - sqrt(pow((cart[0] - cart[1])/2., 2) + pow(cart[2]/2., 2)); //in-plane min Principal Strain 
  } else if (mb==1) {
    fem->StrainB[ec*dof] = cart[0];   fem->StrainB[ec*dof+1] = cart[1];   fem->StrainB[ec*dof+2] = cart[2];  
  }

  //curvilinear to global stress
  for (i=0; i<3; i++) {
    sum = 0.;
    for (j=0; j<3; j++) {
      sum += q[i][j]*stress[j];
    }
    cart[i] = sum;
  }

  if (mb==0) {
    fem->StressM[ec*dof] = cart[0];   fem->StressM[ec*dof+1] = cart[1];   fem->StressM[ec*dof+2] = cart[2]; 
  } else if (mb==1) {
    fem->StressB[ec*dof] = cart[0];   fem->StressB[ec*dof+1] = cart[1];   fem->StressB[ec*dof+2] = cart[2]; 
  }
  
  return(0);
}

//---------------------------------------------------------------------------------------  
PetscErrorCode MembraneNonLinear(PetscInt ec, Cmpnts X1, Cmpnts X2, Cmpnts X3, PetscReal Strain[3], PetscReal _S[3], PetscInt method, FE *fem) {

  PetscReal  q[3][3], xz, xe, ye, yz, Q[3][3], Strainloc[3], sum;
  PetscInt   i, j, m, n, p; 
  Cmpnts     dX21, dX31, V1, V2, N;

  if (method==0) {
    dX21 = MINUS(X2, X1);  dX31 = MINUS(X3, X1);  //dX21:G1 , dX31:G2
    N.x = fem->Nf_x[ec];  N.y = fem->Nf_y[ec];  N.z = fem->Nf_z[ec];
  } else if (method==1) { 
    dX21.x = fem->G1[ec*dof];  dX21.y = fem->G1[ec*dof+1];  dX21.z = fem->G1[ec*dof+2];
    dX31.x = fem->G2[ec*dof];  dX31.y = fem->G2[ec*dof+1];  dX31.z = fem->G2[ec*dof+2];
    N = UNIT(CROSS(dX21, dX31));
  }  

  V1 = UNIT(dX21);
  V2 = CROSS(N, V1);

  xz = DOT(V1, dX21);  xe = DOT(V1, dX31);
  yz = DOT(V2, dX21);  ye = DOT(V2, dX31);

  q[0][0] = xz*xz;  q[0][1] = xe*xe;  q[0][2] = 2*xz*xe;
  q[1][0] = yz*yz;  q[1][1] = ye*ye;  q[1][2] = 2*yz*ye;
  q[2][0] = xz*yz;  q[2][1] = xe*ye;  q[2][2] = xz*ye + xe*yz;

  INV(q, Q); 
 
  // Find theta: angle between the G1 dir and fiber direction
  PetscReal  cosalpha, costheta, theta;
  Cmpnts     n_fib, Nf;
  n_fib.x = fem->n_fib[ec].x;  n_fib.y = fem->n_fib[ec].y;  n_fib.z = fem->n_fib[ec].z;
  costheta = DOT(dX21, n_fib)/(SIZE(dX21)*SIZE(n_fib));
  if (costheta>1.) costheta = 1.;
  if (costheta<-1.) costheta = -1.;
  theta = acos(costheta);
  Nf = CROSS(dX21, n_fib);
  cosalpha = DOT(Nf, N);
  if (cosalpha<0.) theta =- theta;

  // FormRotationMatrix 
  PetscReal  R[3][3];
  R[0][0] = cos(theta)*cos(theta);  R[0][1] = sin(theta)*sin(theta);   R[0][2] = 2*sin(theta)*cos(theta);
  R[1][0] = sin(theta)*sin(theta);  R[1][1] = cos(theta)*cos(theta);   R[1][2] = -2*sin(theta)*cos(theta);
  R[2][0] = -sin(theta)*cos(theta);  R[2][1] = sin(theta)*cos(theta);   R[2][2] = cos(theta)*cos(theta)-sin(theta)*sin(theta);
 
  //Form local strain
  for (i=0; i<3; i++) {
    sum = 0.0;
    for (n=0; n<3; n++) {
      sum += Q[n][i]*Strain[n];
    }
    Strainloc[i] = sum;
  }
  Strainloc[2] = Strainloc[2]/2.;

  //Compute Strain in fiber direction Ef=R Q E
  PetscReal  Ef[3];
  for (i=0; i<3; i++) {
    sum = 0.0;
    for (m=0; m<3; m++) {
  	sum += R[i][m]*Strainloc[m];
    }
    Ef[i] = sum;
  }

  //  Tissue Constitutive Law according to Fung exponential form
  /* PetscReal c=9.7e3, A1=49.558, A2=5.2871, A3=-3.124, A4=16.031 ,A5=-.004, A6=-0.02; // Values for BHV hammer et al, ABME 2011 */
  PetscReal  RhoU2 = 1200*0.877*0.877;
  PetscReal  c=14.42e3/RhoU2, A1=61.27, A2=70.37, A3=5.11, A4=14.2, A5=3.1, A6=2.01;//Kim et al 2008 annuals of biomedical engineering
  PetscReal  Qf, Sf[3]; //f:fiber direction
  Qf = A1*Ef[0]*Ef[0]+A2*Ef[1]*Ef[1]+2*A3*Ef[0]*Ef[1]+
    A4*Ef[2]*Ef[2]+2*A5*Ef[0]*Ef[2]+2*A6*Ef[1]*Ef[2];
  
  Sf[0] = c*(A1*Ef[0]+A3*Ef[1]+A5*Ef[2])*exp(Qf);
  Sf[1] = c*(A3*Ef[0]+A2*Ef[1]+A6*Ef[2])*exp(Qf);
  Sf[2] = c*(A5*Ef[0]+A6*Ef[1]+A4*Ef[2])*exp(Qf);

  // Convert Sf to Scur Scur=QT Rinv Sf
  PetscReal  Rinv[3][3];
  INV(R, Rinv);
  for (i=0; i<3; i++) {
    sum = 0.0;
    for (m=0; m<3; m++) {
      for (p=0; p<3; p++) {
	sum += Q[i][m]*Rinv[m][p]*Sf[p];
      }
    }
    _S[i] = sum;
  }
  
  return(0);
}

//---------------------------------------------------------------------------------------  
PetscErrorCode BendingNonLinear(PetscInt ec, Cmpnts X1, Cmpnts X2, Cmpnts X3, PetscReal Strain[3], PetscReal _S[3], PetscInt method, FE *fem) {
  
  PetscReal  q[3][3], xz, xe, ye, yz, Q[3][3], Strainloc[3], sum;
  PetscInt   i, j, m, n, p; 
  Cmpnts     dX21, dX31, V1, V2, N;

  if (method==0) {
    dX21 = MINUS(X2, X1);  dX31 = MINUS(X3, X1);  //dX21:G1 , dX31:G2
    N.x = fem->Nf_x[ec];  N.y = fem->Nf_y[ec];  N.z = fem->Nf_z[ec];
  } else if (method==1) { 
    dX21.x = fem->G1[ec*dof];  dX21.y = fem->G1[ec*dof+1];  dX21.z = fem->G1[ec*dof+2];
    dX31.x = fem->G2[ec*dof];  dX31.y = fem->G2[ec*dof+1];  dX31.z = fem->G2[ec*dof+2];
    N = UNIT(CROSS(dX21, dX31));
  }  

  V1 = UNIT(dX21);
  V2 = CROSS(N, V1);

  xz = DOT(V1, dX21);  xe = DOT(V1, dX31);
  yz = DOT(V2, dX21);  ye = DOT(V2, dX31);

  q[0][0] = xz*xz;  q[0][1] = xe*xe;  q[0][2] = 2*xz*xe;
  q[1][0] = yz*yz;  q[1][1] = ye*ye;  q[1][2] = 2*yz*ye;
  q[2][0] = xz*yz;  q[2][1] = xe*ye;  q[2][2] = xz*ye + xe*yz;

  INV(q, Q); 
 
  // Find theta: angle between the G1 dir and fiber direction
  PetscReal  cosalpha, costheta, theta;
  Cmpnts     n_fib, Nf;
  n_fib.x = fem->n_fib[ec].x;  n_fib.y = fem->n_fib[ec].y;  n_fib.z = fem->n_fib[ec].z;
  costheta = DOT(dX21, n_fib)/(SIZE(dX21)*SIZE(n_fib));
  if (costheta>1.) costheta = 1.;
  if (costheta<-1.) costheta = -1.;
  theta = acos(costheta);
  Nf = CROSS(dX21, n_fib);
  cosalpha = DOT(Nf, N);
  if (cosalpha<0.) theta =- theta;

  // FormRotationMatrix 
  PetscReal  R[3][3];
  R[0][0] = cos(theta)*cos(theta);  R[0][1] = sin(theta)*sin(theta);   R[0][2] = 2*sin(theta)*cos(theta);
  R[1][0] = sin(theta)*sin(theta);  R[1][1] = cos(theta)*cos(theta);   R[1][2] =-2*sin(theta)*cos(theta);
  R[2][0] = -sin(theta)*cos(theta);  R[2][1] = sin(theta)*cos(theta);   R[2][2] = cos(theta)*cos(theta)-sin(theta)*sin(theta);
 
  //Form local strain
  for (i=0; i<3; i++) {
    sum = 0.0;
    for (n=0; n<3; n++) {
      sum += Q[n][i]*Strain[n];
    }
    Strainloc[i] = sum;
  }
  Strainloc[2] = Strainloc[2]/2.;

  //Compute Strain in fiber direction Ef=R Q E
  PetscReal  Ef[3];
  for (i=0; i<3; i++) {
    sum = 0.0;
    for (m=0; m<3; m++) {
  	sum += R[i][m]*Strainloc[m];
    }
    Ef[i] = sum;
  }

  //  Tissue Constitutive Law according to Fung exponential form
  PetscReal  RhoU2=1200*0.877*0.877, RhoU2D=1200*0.877*0.877*0.022;
  PetscReal  a1=443000/RhoU2, a2=620000/RhoU2, b1=4302/RhoU2D, b2=6023/RhoU2D; //Kim et al 2008 annuals of biomedical engineering
  PetscReal  Sf[3]; //f:fiber direction
  
  Sf[0] = a1*Ef[0] + b1*SIGN(Ef[0])*Ef[0]*Ef[0] + 0.25*(a1 + a2)*Ef[1];
  Sf[1] = 0.25*(a1 + a2)*Ef[0] + a2*Ef[1] + b2*SIGN(Ef[1])*Ef[1]*Ef[1];
  Sf[2] = 0.25*(a1 + a2)*Ef[2];

  // Convert Sf to Scur Scur=QT Rinv Sf
  PetscReal  Rinv[3][3];
  INV(R, Rinv);
  for (i=0; i<3; i++) {
    sum = 0.0;
    for (m=0; m<3; m++) {
      for (p=0; p<3; p++) {
	sum += Q[i][m]*Rinv[m][p]*Sf[p];
      }
    }
    _S[i] = sum;
  }
  
  return(0);
}
