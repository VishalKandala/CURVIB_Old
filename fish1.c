#include "variables.h"
extern PetscReal CMx_c,CMy_c,CMz_c;
extern PetscInt ti,tistart;
extern PetscInt tiout,eel,pizza;
extern PetscReal St_exp,wavelength;
PetscReal ampl,V_slip,omega, kwave,alpha;
PetscReal T_period_fish;
PetscReal c_0,c_1,c_2;
PetscInt  N_period_fish;
PetscReal Y,Z,yy,zz,wn,z2,theta,Wl;

PetscErrorCode fish_init(PetscReal *delti)
{
  PetscReal  pi = 3.141592653589793, St;
  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  PetscOptionsGetInt(PETSC_NULL, "-N_period_fish", &N_period_fish, PETSC_NULL);
  
  PetscPrintf(PETSC_COMM_WORLD, "N_period_fish is : %d \n",N_period_fish);
  PetscOptionsGetReal(PETSC_NULL, "-theta", &theta, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-Croty", &yy, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-Crotz", &zz, PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-Wavelength", &Wl, PETSC_NULL);
  theta=-pi/36;
  wn=2*pi/Wl;
  
  ampl=.1;
  z2=zz-0.3;
  PetscPrintf(PETSC_COMM_WORLD, " tistart %d omega %le z2 %le  wavenumber %le\n",tistart,omega,z2,wn);
  //*amir

  // V_slip=.7;

  // alpha = 2.18;
  // if (eel) 
  //  ampl=0.1;//0.089;
  // else
  //  ampl=0.1;//0.09 // non-dim ampl a/L
  // c_0=0.02;
  // c_1=-0.08;
  // c_2=0.16;
/*   c_0=0.02; */
/*   c_2=(0.3*(ampl-c_0)+c_0)/0.291;  */
/*   c_1=ampl - c_0 - c_2 ;//  ! Ampl. is defined at the end of tail */

  // St=St_exp/(2.*ampl);  // non-dim  frequency St=fL/U
  // omega=2*pi*St;//2.*pi; //non-dimensionalized w
  omega =2*pi*St_exp;
  //kwave=2.*pi*St*V_slip; //non-dimensionalized k
  //wavelength = 0.95;

  V_slip= kwave*0.5/pi/St;

  T_period_fish=2.*pi/omega;
  *delti = T_period_fish/N_period_fish;
 
  PetscPrintf(PETSC_COMM_WORLD, "fish init: St_exp %le Amplitude %le coeff a0 a1 a2 %le %le %le\n",St_exp,ampl,c_0,c_1,c_2);
  PetscPrintf(PETSC_COMM_WORLD, "fish init: dt %le w %le f %le k %le lambda %le T %le N %d  V_slip %le\n",*delti,omega,St,kwave,wavelength,T_period_fish,N_period_fish,V_slip);
 
  return(0);
}

PetscErrorCode fish_swim(IBMNodes *ibm, PetscReal time
			,PetscReal delti)
{

  PetscReal  pi = 3.141592653589793;

  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  PetscOptionsGetInt(PETSC_NULL, "-pizza", &pizza, PETSC_NULL);
  // PetscOptionsGetInt(PETSC_NULL, "-theta", &theta, PETSC_NULL);
  // PetscPrintf(PETSC_COMM_WORLD, "N_period_fish is : %d theta = %le  \n",N_period_fish,theta);
 
  PetscReal  h0,h1,h2,h3,h4;

  PetscInt   i, n_v=ibm->n_v, n_elmt=ibm->n_elmt;
  PetscReal  x0,x1,y0,y1,z0,z1,beta,r;
  PetscReal  az0,dz,k,theta;

  y0=1.e+5;
  y1=-1.e+5;
  x0=y0;
  x1=y1;
  z0=y0;
  z1=y1;
  for (i=0; i<n_v; i++) {
    x0=PetscMin(x0,(ibm->x_bp0[i]));
    x1=PetscMax(x1,(ibm->x_bp0[i]));
    y0=PetscMin(y0,(ibm->y_bp0[i]));
    y1=PetscMax(y1,(ibm->y_bp0[i]));
    z0=PetscMin(z0,(ibm->z_bp0[i]));
    z1=PetscMax(z1,(ibm->z_bp0[i]));
  }

  ampl=0.10;
  PetscPrintf(PETSC_COMM_WORLD, "MAX fish: %d %le %le %le %le %le %le\n",ti,z1,z0,y1,y0,x1,x0);
  

  for (i=0; i<n_v; i++) {
 
    //wave on y direction
    
    ibm->y_bp[i] = ibm->y_bp0[i];
    ibm->z_bp[i] = ibm->z_bp0[i];
    ibm->x_bp[i] = ibm->x_bp0[i];
    if (ibm->z_bp0[i]>z2){
	dz=ibm->z_bp0[i]-z2;
	//	az0=ampl*exp(alpha*(dz-1.));
	az0=ampl*(dz)/(.8);
	//	ibm->y_bp[i]=ibm->y_bp0[i]+ampl*(1-exp(-20*(time-tistart*delti)))*(1-exp(-10*(ibm->z_bp0[i]-z2)))*sin(omega*time-wn*(ibm->z_bp0[i]-z2));
	ibm->y_bp[i]=ibm->y_bp0[i]+az0*sin(wn*dz-omega*time);
	//ibm->y_bp[i]=ibm->y_bp0[i]+az0*sin(-omega*(time-tistart*delti));
    }
  }

    //rotation
  for (i=0; i<n_v; i++) {
      Y=ibm->y_bp[i];
      Z=ibm->z_bp[i];
      ibm->y_bp[i]=yy+(Z-zz)*sin(theta)+(Y-yy)*cos(theta);
      ibm->z_bp[i]=zz+(Z-zz)*cos(theta)-(Y-yy)*sin(theta);
  }

  
  
  /*    Calculate the new normal & velcity */
  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  
  for (i=0; i<n_elmt; i++) {
    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
    dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e]; 
    dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e]; 
    dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e]; 
    
    dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e]; 
    dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e]; 
    dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e]; 
    
    ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
    ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
    ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;
    
    dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] + 
	      ibm->nf_z[i]*ibm->nf_z[i]);
    
    ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;
  }
  
  PetscReal  v_max=0.;
  PetscInt   i_vmax=0;
  if (ti==0) {
    for (i=0; i<n_v; i++) {
      ibm->u[i].x = 0.0;
      ibm->u[i].y = 0.0;
      ibm->u[i].z = 0.0;
    } 

  }else {
   for (i=0; i<n_v; i++) {
     ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / delti;
     ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / delti;
     ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / delti;
    
     if (v_max<fabs(ibm->u[i].z)) {
       i_vmax=i;
       v_max= fabs(ibm->u[i].z);
     }
     if (v_max<fabs(ibm->u[i].y)) {
       i_vmax=i;
       v_max= fabs(ibm->u[i].y);
     }
     if (v_max<fabs(ibm->u[i].x)) {
       i_vmax=i;
       v_max= fabs(ibm->u[i].x);
     }
   }
  }
 

  PetscPrintf(PETSC_COMM_WORLD, "MAX fish Velocity: %d %le %le %le %le\n",ti, v_max, ibm->x_bp[i_vmax],ibm->y_bp[i_vmax],ibm->z_bp[i_vmax]);
  


  return(0);
}
PetscErrorCode pizza_box_swim(IBMNodes *ibm, PetscReal time
			,PetscReal delti)
{

 
  PetscReal  pi = 3.141592653589793;
  PetscReal  h0=0.,h1=0.,h2=0.,h3=0.,h4=0.;

  PetscInt   i, n_v=ibm->n_v, n_elmt=ibm->n_elmt;
  PetscReal  x0,x1,y0,y1,z0,z1;
  PetscReal  az0,dz,dzz;
  PetscInt EXP=0;

  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  PetscOptionsGetInt(PETSC_NULL, "-EXP", &EXP, PETSC_NULL);
  PetscPrintf(PETSC_COMM_WORLD, "EXP value : %d \n",EXP); 

  y0=1.e+5;
  y1=-1.e+5;
  x0=y0;
  x1=y1;
  z0=y0;
  z1=y1;
  for (i=0; i<n_v; i++) {
    x0=PetscMin(x0,(ibm->x_bp0[i]));
    x1=PetscMax(x1,(ibm->x_bp0[i]));
    y0=PetscMin(y0,(ibm->y_bp0[i]));
    y1=PetscMax(y1,(ibm->y_bp0[i]));
    z0=PetscMin(z0,(ibm->z_bp0[i]));
    z1=PetscMax(z1,(ibm->z_bp0[i]));
  }

  PetscPrintf(PETSC_COMM_WORLD, "MAX fish: %d %le %le %le %le %le %le\n",ti,z1,z0,y1,y0,x1,x0);
  
  PetscReal A11=20.4011,A12=-17.3354,A13=-15.8377;
  PetscReal A21=-11.3845,A22=28.9923,A23=32.7113;
  PetscReal A31=-16.2930,A32=-0.9307,A33=-19.1176;
  PetscReal A41=-3.1176,A42=-17.8372,A43=24.1467;
  
  PetscReal B21=9.5921,B22=-4.3296,B23=-6.0511;
  PetscReal B31=-2.154,B32=4.1442,B33=-5.5296;
  PetscReal B41=7.4193,B42=7.4872,B43=-5.0809;
  
  PetscReal C11=-24.2414,C12=20.5414,C13=16.9785;
  PetscReal C21=-28.2899,C22=14.1204,C23=-0.2142;
  PetscReal C31=3.9762,C32=-17.5688,C33=2.4505;
  PetscReal C41=26.5090,C42=18.9211,C43=-3.5616;
 
  PetscReal D21=-3.8377,D22=6.4793,D23=6.7913;
  PetscReal D31=6.0404,D32=-3.5397,D33=5.5417;
  PetscReal D41=7.9833,D42=0.0597,D43=1.3695;

  PetscReal lambda11=35.2902*0.10096,lambda12=62.3155*0.10096,lambda13=87.7894*0.10096;
  PetscReal lambda21=49.3038*0.10096,lambda22=87.0608*0.10096,lambda23=122.6502*0.10096;
  PetscReal lambda31=49.3038*0.10096,lambda32=87.0608*0.10096,lambda33=122.6502*0.10096;
  PetscReal lambda41=35.2902*0.10096,lambda42=62.3155*0.10096,lambda43=87.7894*0.10096;

  PetscReal alpha11=0.0004e-5/0.10096,alpha12=0.0034e-5/0.10096,alpha13=0.3647e-5/0.10096;
  PetscReal alpha21=0.0004e-5/0.10096,alpha22=0.0048e-5/0.10096,alpha23=0.4978e-5/0.10096;

  PetscReal Volt=500.;
  PetscReal Phi1,Phi2,Phi3;
  PetscReal theta11=0.0033,theta12=-3.1281,theta13=-1.5708,theta21=0.0033,theta22=0.0135,theta23=-1.5708;

  if (EXP){
   
    for (i=0; i<n_v; i++) {
      dz=(ibm->z_bp0[i]-z0) ;
      if(dz <=0.25)  {
	dzz=dz;
	Phi1=A11*(cos(lambda11*dzz)+cosh(lambda11*dzz))+C11*(sin(lambda11*dzz)+sinh(lambda11*dzz));
	Phi2=A12*(cos(lambda12*dzz)+cosh(lambda12*dzz))+C12*(sin(lambda12*dzz)+sinh(lambda12*dzz));
	Phi3=A13*(cos(lambda13*dzz)+cosh(lambda13*dzz))+C13*(sin(lambda13*dzz)+sinh(lambda13*dzz));
	az0=(1.347*(alpha11*Phi1*sin(omega*time+theta11)+alpha12*Phi2*sin(omega*time+theta12)+alpha13*Phi3*sin(omega*time+theta13))
	     +alpha21*Phi1*sin(omega*time-pi+theta21)+alpha22*Phi2*sin(omega*time-pi+theta22)+alpha23*Phi3*sin(omega*time-pi+theta23))*Volt;
      }
      else  if(dz <=0.50) {
	dzz=dz-0.25;
	Phi1=A21*cos(lambda21*dzz)+B21*cosh(lambda21*dzz)+C21*sin(lambda21*dzz)+D21*sinh(lambda21*dzz);
	Phi2=A22*cos(lambda22*dzz)+B22*cosh(lambda22*dzz)+C22*sin(lambda22*dzz)+D22*sinh(lambda22*dzz);
	Phi3=A23*cos(lambda23*dzz)+B23*cosh(lambda23*dzz)+C23*sin(lambda23*dzz)+D23*sinh(lambda23*dzz);
	az0=(1.347*(alpha11*Phi1*sin(omega*time+theta11)+alpha12*Phi2*sin(omega*time+theta12)+alpha13*Phi3*sin(omega*time+theta13))
	     +alpha21*Phi1*sin(omega*time-pi+theta21)+alpha22*Phi2*sin(omega*time-pi+theta22)+alpha23*Phi3*sin(omega*time-pi+theta23))*Volt;
      }
      else  if(dz <=0.75) {
	dzz=dz-0.5;
	Phi1=A31*cos(lambda31*dzz)+B31*cosh(lambda31*dzz)+C31*sin(lambda31*dzz)+D31*sinh(lambda31*dzz);
	Phi2=A32*cos(lambda32*dzz)+B32*cosh(lambda32*dzz)+C32*sin(lambda32*dzz)+D32*sinh(lambda32*dzz);
	Phi3=A33*cos(lambda33*dzz)+B33*cosh(lambda33*dzz)+C33*sin(lambda33*dzz)+D33*sinh(lambda33*dzz);
	az0=(1.347*(alpha11*Phi1*sin(omega*time+theta11)+alpha12*Phi2*sin(omega*time+theta12)+alpha13*Phi3*sin(omega*time+theta13))
	     +alpha21*Phi1*sin(omega*time-pi+theta21)+alpha22*Phi2*sin(omega*time-pi+theta22)+alpha23*Phi3*sin(omega*time-pi+theta23))*Volt;
      } 
      else {
	dzz=dz-0.75;
	Phi1=A41*cos(lambda41*dzz)+B41*cosh(lambda41*dzz)+C41*sin(lambda41*dzz)+D41*sinh(lambda41*dzz);
	Phi2=A42*cos(lambda42*dzz)+B42*cosh(lambda42*dzz)+C42*sin(lambda42*dzz)+D42*sinh(lambda42*dzz);
	Phi3=A43*cos(lambda43*dzz)+B43*cosh(lambda43*dzz)+C43*sin(lambda43*dzz)+D43*sinh(lambda43*dzz);
	az0=(1.347*(alpha11*Phi1*sin(omega*time+theta11)+alpha12*Phi2*sin(omega*time+theta12)+alpha13*Phi3*sin(omega*time+theta13))
	     +alpha21*Phi1*sin(omega*time-pi+theta21)+alpha22*Phi2*sin(omega*time-pi+theta22)+alpha23*Phi3*sin(omega*time-pi+theta23))*Volt;
      } 
      if (dz==0.0)  h0=az0;
      if (dz==0.25) h1=az0;
      if (dz==0.5)  h2=az0;
      if (dz==0.75) h3=az0;
      if (dz==1.0)  h4=az0;   

      ibm->x_bp[i]= ibm->x_bp0[i]+az0;
    }
    PetscPrintf(PETSC_COMM_WORLD, " t %le  h0 %le  h1 %le h2 %le h3 %le h4 %le\n",time,h0,h1,h2,h3,h4);
  }
  else{
    if (eel) {
      h0=ampl*exp(alpha*(-1.))*sin(-omega*time);
      h1=ampl*exp(alpha*(0.25-1.))*sin(kwave*0.25-omega*time);
      h2=ampl*exp(alpha*(0.5-1.))*sin(kwave*0.5-omega*time);
      h3=ampl*exp(alpha*(0.75-1.))*sin(kwave*0.75-omega*time);
      h4=ampl*sin(kwave*1.0-omega*time);
    }else{
      h0=c_0*ampl*sin(-omega*time);
      h1=(c_0+c_1*0.25+c_2*0.25*0.25)*sin(kwave*0.25-omega*time);
      h2=(c_0+c_1*0.5+c_2*0.5*0.5)*sin(kwave*0.5-omega*time);
      h3=(c_0+c_1*0.75+c_2*0.75*0.75)*sin(kwave*0.75-omega*time);
      h4=(c_0+c_1+c_2)*sin(kwave*1.0-omega*time);
    }
    PetscPrintf(PETSC_COMM_WORLD, " t %le  h0 %le  h1 %le h2 %le h3 %le h4 %le\n",time,h0,h1,h2,h3,h4);
    
    /*   oscillation for Mackerel in x-dir */
    for (i=0; i<n_v; i++) {
      dz=(ibm->z_bp0[i]-z0);
      
      if(dz <=0.25)        ibm->x_bp[i]=ibm->x_bp0[i]+h1+(dz-0.25)*(h1-h0)/(0.25);
      else  if(dz <=0.50)  ibm->x_bp[i]=ibm->x_bp0[i]+h2+(dz-0.5)*(h2-h1)/(0.25);
      else  if(dz <=0.75)  ibm->x_bp[i]=ibm->x_bp0[i]+h3+(dz-0.75)*(h3-h2)/(0.25);
      else  if(dz <=1.0)   ibm->x_bp[i]=ibm->x_bp0[i]+h4+(dz-1.0)*(h4-h3)/(0.25);
      else                 ibm->x_bp[i]=ibm->x_bp0[i]+0.0;      
    }
  }
  
/*    Calculate the new normal & velcity */
  PetscInt n1e, n2e, n3e;
  PetscReal dx12, dy12, dz12, dx13, dy13, dz13, dr;
  
  for (i=0; i<n_elmt; i++) {
    n1e = ibm->nv1[i]; n2e =ibm->nv2[i]; n3e =ibm->nv3[i];
    dx12 = ibm->x_bp[n2e] - ibm->x_bp[n1e]; 
    dy12 = ibm->y_bp[n2e] - ibm->y_bp[n1e]; 
    dz12 = ibm->z_bp[n2e] - ibm->z_bp[n1e]; 
    
    dx13 = ibm->x_bp[n3e] - ibm->x_bp[n1e]; 
    dy13 = ibm->y_bp[n3e] - ibm->y_bp[n1e]; 
    dz13 = ibm->z_bp[n3e] - ibm->z_bp[n1e]; 
    
    ibm->nf_x[i] = dy12 * dz13 - dz12 * dy13;
    ibm->nf_y[i] = -dx12 * dz13 + dz12 * dx13;
    ibm->nf_z[i] = dx12 * dy13 - dy12 * dx13;
    
    dr = sqrt(ibm->nf_x[i]*ibm->nf_x[i] + ibm->nf_y[i]*ibm->nf_y[i] + 
	      ibm->nf_z[i]*ibm->nf_z[i]);
    
    ibm->nf_x[i] /=dr; ibm->nf_y[i]/=dr; ibm->nf_z[i]/=dr;
  }
  
  PetscReal  v_max=0.;
  PetscInt   i_vmax=0;
  for (i=0; i<n_v; i++) {
    ibm->u[i].x = (ibm->x_bp[i] - ibm->x_bp_o[i]) / delti;
    ibm->u[i].y = (ibm->y_bp[i] - ibm->y_bp_o[i]) / delti;
    ibm->u[i].z = (ibm->z_bp[i] - ibm->z_bp_o[i]) / delti;
    
    if (v_max<fabs(ibm->u[i].z)) {
      i_vmax=i;
      v_max= fabs(ibm->u[i].z);
    }
    if (v_max<fabs(ibm->u[i].y)) {
      i_vmax=i;
      v_max= fabs(ibm->u[i].y);
    }
    if (v_max<fabs(ibm->u[i].x)) {
      i_vmax=i;
      v_max= fabs(ibm->u[i].x);
    }
  }
  PetscPrintf(PETSC_COMM_WORLD, "MAX fish Velocity: %d %le %le %le %le\n",ti, v_max, ibm->x_bp[i_vmax],ibm->y_bp[i_vmax],ibm->z_bp[i_vmax]);
  
  return(0);
}
