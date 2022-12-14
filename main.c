#include "variables.h"

static char help[] = "Testing programming!";

PetscInt ti,tistart=0;
PetscReal	Flux_in = 4.104388e-04, angle = 0;
PetscInt tiout = 10;
PetscInt block_number;
PetscReal FluxInSum, FluxOutSum;
/* immeresed value>1 determines the type of correction
   1      constant velocity correction
   2      proportional (needs FluxOutSum>0)
   3      proportional to flux
   4      proportional to normal velocity (flux/area)
*/
PetscInt immersed = 0; 
PetscInt invicid = 0;
PetscInt movefsi = 0, rotatefsi=0, sediment=0;
PetscBool rstart_flg;
PetscInt implicit = 0, implicit_type=0;
PetscInt imp_MAX_IT = 50; 
PetscInt radi=10;
PetscInt inletprofile=2, InitialGuessOne=0, period=0;
PetscReal CMx_c=0., CMy_c=0., CMz_c=0.;
PetscInt  mg_MAX_IT=30, mg_idx=1, mg_preItr=1, mg_poItr=1;
PetscReal imp_atol=1e-7, imp_rtol=1e-4, imp_stol=1.e-8;
PetscInt TwoD = 0;
PetscInt STRONG_COUPLING=0;
PetscInt rstart_fsi=0;
PetscInt cop=0, regime=1; // 1 escape regime --- 2 cruise regime
PetscInt fish=0, fish_c=0, eel=0, fishcyl=0, rheology=0,duplicate=0,turbine=0,Pipe=0,pizza=0;
PetscInt wing=0, hydro=0;
PetscReal St_exp=0.5,wavelength=0.95;
PetscInt MHV=0,LV=0, BHV=0, beam=0, inv_flg=0;
PetscReal max_angle = -0.87266;// depends on the initial angle 
PetscInt thin=0;
PetscInt dgf_z=0,dgf_y=1,dgf_x=0;
PetscInt dgf_az=0,dgf_ay=0,dgf_ax=1 ;
PetscInt NumberOfBodies=1;
PetscInt moveframe=0,rotateframe=0, blank=0;
PetscReal L_dim;
PetscInt les=0, rans=0, channelz=0;
PetscInt wallfunction=0;
PetscInt averaging=0;
PetscInt phave=0;
int mixed=0;
int clark=0;
int dynamic_freq=1;
int poisson=0;
int periodic=0;
PetscReal max_cs=0.5;
int grid1d=0;
int i_periodic=0;
int j_periodic=0;
int k_periodic=0;
PetscInt blkpbc=10;
int pseudo_periodic=0;
int testfilter_ik=0;
int testfilter_1d=0;
int i_homo_filter=0;
int j_homo_filter=0;
int k_homo_filter=0;
double poisson_tol=5.e-9;	// relative tolerance

//S-----------------fem variables------------------------------------------------
PetscReal  E=0.0, mu=0.0, rho=0.0, h0=0.0, dampfactor=0.0;
PetscInt   dof=3, twod=0, damping=0, membrane=0, bending=0, outghost=0, ConstitutiveLawNonLinear=0;
PetscInt   nbody=1, timeinteg=0, contact=0, explicit=0, initrest=0, innerloop=1, rstart_fem=0;
//E-----------------fem variables------------------------------------------------

double time_flow[3000],flux_flow[3000];
//IBMNodes	*ibm_ptr;
double PI = 3.141592653589793;

PetscErrorCode Ucont_P_Binary_Input(UserCtx *user)
{

  PetscViewer	viewer;

  char filen[90];

  
  sprintf(filen, "vfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  PetscInt N;

  VecGetSize(user->Ucont, &N);
  PetscPrintf(PETSC_COMM_WORLD, "PPP %d\n", N);
  VecLoad((user->Ucont),viewer);
  
  PetscViewerDestroy(&viewer);

  PetscBarrier(PETSC_NULL);

  PetscOptionsClearValue("-vecload_block_size");

  sprintf(filen, "pfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad((user->P),viewer);
  PetscViewerDestroy(&viewer);

  sprintf(filen, "nvfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
  VecLoad((user->Nvert_o),viewer);
  PetscViewerDestroy(&viewer);

  DMGlobalToLocalBegin(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);
  DMGlobalToLocalEnd(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o);

  if(averaging) {	// Seokkoo Kang
    sprintf(filen, "su0_%06d_%1d.dat",ti, user->_this);
    FILE *fp=fopen(filen, "r");
    
    VecSet(user->Ucat_sum, 0);
    VecSet(user->Ucat_cross_sum, 0);
    VecSet(user->Ucat_square_sum, 0);
    VecSet(user->P_sum, 0);
    VecSet(user->lUstar_sum, 0);

 

    if(fp==NULL) {
      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, setting the statistical quantities to zero and contiues the computation ... ***\n\n", filen);
    }
    else {
      PetscOptionsClearValue("-vecload_block_size");
      fclose(fp);
      PetscBarrier(PETSC_NULL);
      sprintf(filen, "su0_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->Ucat_sum,viewer);
      PetscViewerDestroy(&viewer);
      
      sprintf(filen, "su1_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad(user->Ucat_cross_sum,viewer);
      PetscViewerDestroy(&viewer);
      
      sprintf(filen, "su2_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad((user->Ucat_square_sum),viewer);
      PetscViewerDestroy(&viewer);
      
      sprintf(filen, "sp_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->P_sum,viewer);
      PetscViewerDestroy(&viewer);
 
      sprintf(filen, "sustar_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->lUstar_sum,viewer);
      PetscViewerDestroy(&viewer);

      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Read %s, continuing averaging ... ***\n\n", filen);
    }
  }


    if(averaging && phave) {	// Amir
    sprintf(filen, "suph0_%06d_%1d.dat",ti, user->_this);
    FILE *fp=fopen(filen, "r");

     VecSet(user->Ucat_sum0, 0); 
     VecSet(user->Ucat_sum1, 0);
     VecSet(user->Ucat_sum2, 0); 
     VecSet(user->Ucat_sum3, 0);
     VecSet(user->P_sum0, 0);
     VecSet(user->P_sum1, 0);
     VecSet(user->P_sum2, 0);
     VecSet(user->P_sum3, 0);
 //   } 

     if(fp==NULL) {
      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, setting the statistical quantities to zero and contiues the computation ... ***\n\n", filen);
    }
    else {
      PetscOptionsClearValue("-vecload_block_size");
      fclose(fp);
      PetscBarrier(PETSC_NULL);
      
      sprintf(filen, "suph0_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->Ucat_sum0,viewer);
      PetscViewerDestroy(&viewer);
 
      sprintf(filen, "suph1_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->Ucat_sum1,viewer);
      PetscViewerDestroy(&viewer);
 
      sprintf(filen, "suph2_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->Ucat_sum2,viewer);
      PetscViewerDestroy(&viewer);
 
      sprintf(filen, "suph3_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->Ucat_sum3,viewer);
      PetscViewerDestroy(&viewer);

      sprintf(filen, "spph0_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->P_sum0,viewer);
      PetscViewerDestroy(&viewer);
 
      sprintf(filen, "spph1_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->P_sum1,viewer);
      PetscViewerDestroy(&viewer);
 
      sprintf(filen, "spph2_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->P_sum2,viewer);
      PetscViewerDestroy(&viewer);
 
      sprintf(filen, "spph3_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( user->P_sum3,viewer);
      PetscViewerDestroy(&viewer);
     // }

      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Read %s, continuing averaging ... ***\n\n", filen);
    }
  }
 

  
  if(les) {
    Vec Cs;
    
    VecDuplicate(user->P, &Cs);
    
    sprintf(filen, "cs_%06d_%1d.dat", ti, user->_this);
    FILE *fp=fopen(filen, "r");
    
    if(fp==NULL) {
      PetscPrintf(PETSC_COMM_WORLD,"\n\n*** Cannot open %s, setting Cs to 0 and contiues the computation ... ***\n\n", filen);
      VecSet(Cs, 0);
    }
    else {
      fclose(fp);
      
      PetscBarrier(PETSC_NULL);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad( Cs,viewer);
      PetscViewerDestroy(&viewer);
    }
    
    DMGlobalToLocalBegin(user->da, Cs, INSERT_VALUES, user->lCs);
    DMGlobalToLocalEnd(user->da, Cs, INSERT_VALUES, user->lCs);
    
    VecDestroy(&Cs);
  }

  
   if(rans) {
    // K-Omega
    sprintf(filen, "kfield%06d_%1d.dat", ti, user->_this);
    FILE *fp=fopen(filen, "r");
    
    if(fp!=NULL) {
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
      VecLoad(user->K_Omega,viewer);
      PetscViewerDestroy(&viewer);
    }
    else {
      K_Omega_IC(user);
    }
    
    VecCopy(user->K_Omega, user->K_Omega_o);
    
    DMGlobalToLocalBegin(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega, INSERT_VALUES, user->lK_Omega);
    
    DMGlobalToLocalBegin(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);
    DMGlobalToLocalEnd(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o);    
  }

  return 0;
}

PetscErrorCode Ucont_P_Binary_Output(UserCtx *user)
{
  PetscViewer	viewer;
  char filen[80];

  int rank;
  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
 

  sprintf(filen, "nvfield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  
  VecView(user->Nvert, viewer);
  PetscViewerDestroy(&viewer);
 

   PetscBarrier(PETSC_NULL);

 
  sprintf(filen, "vfield%5.5d_%1.1d.dat", ti, user->_this);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  
  VecView(user->Ucont, viewer);
  
  PetscViewerDestroy(&viewer);
  
  
  sprintf(filen, "ufield%5.5d_%1.1d.dat", ti, user->_this);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  
  VecView(user->Ucat, viewer);
  
  PetscViewerDestroy(&viewer);
  
  sprintf(filen, "pfield%5.5d_%1.1d.dat", ti, user->_this);
  
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
  
  VecView(user->P, viewer);
  PetscViewerDestroy(&viewer);
   
  if(averaging && ti!=0) {	// Seokkoo Kang
    sprintf(filen, "su0_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_sum, viewer);
    PetscViewerDestroy(&viewer);
    sprintf(filen, "su0_%06d_%1d.dat.info", ti, user->_this);	if(!rank) unlink(filen);
    
   PetscBarrier(PETSC_NULL);
    
    sprintf(filen, "su1_%06d_%1d.dat",  ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_cross_sum, viewer);
    PetscViewerDestroy(&viewer);  
    sprintf(filen, "su1_%06d_%1d.dat.info", ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(PETSC_NULL);
    
    sprintf(filen, "su2_%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->Ucat_square_sum, viewer);
    PetscViewerDestroy(&viewer);
    sprintf(filen, "su2_%06d_%1d.dat.info",ti, user->_this);	if(!rank) unlink(filen);
    
    PetscBarrier(PETSC_NULL);
    
    sprintf(filen, "sp_%06d_%1d.dat",ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->P_sum, viewer);
    PetscViewerDestroy(&viewer);
    sprintf(filen, "sp_%06d_%1d.dat.info",ti, user->_this);	if(!rank) unlink(filen);
   
    sprintf(filen, "sustar_%06d_%1d.dat",ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->lUstar_sum, viewer);
    PetscViewerDestroy(&viewer);
    sprintf(filen, "sustar_%06d_%1d.dat.info",ti, user->_this);	if(!rank) unlink(filen);
       
    PetscBarrier(PETSC_NULL);

  if (phave){
      sprintf(filen, "suph0_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
      VecView(user->Ucat_sum0, viewer);
      PetscViewerDestroy(&viewer);
      sprintf(filen, "suph0_%06d_%1d.dat.info", ti, user->_this);	if(!rank) unlink(filen);
   
       PetscBarrier(PETSC_NULL);
 
      sprintf(filen, "suph1_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
      VecView( user->Ucat_sum1,viewer);
      PetscViewerDestroy(&viewer);
      sprintf(filen, "suph1_%06d_%1d.dat.info", ti, user->_this);	if(!rank) unlink(filen);

        PetscBarrier(PETSC_NULL);
 
      sprintf(filen, "suph2_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
      VecView(user->Ucat_sum2, viewer);
      PetscViewerDestroy(&viewer);
      sprintf(filen, "suph2_%06d_%1d.dat.info", ti, user->_this);	if(!rank) unlink(filen);
   
      PetscBarrier(PETSC_NULL);
 
      sprintf(filen, "suph3_%06d_%1d.dat", ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
      VecView(user->Ucat_sum3, viewer);
      PetscViewerDestroy(&viewer);
      sprintf(filen, "suph3_%06d_%1d.dat.info", ti, user->_this);	if(!rank) unlink(filen);
   
      PetscBarrier(PETSC_NULL);
 
      sprintf(filen, "spph0_%06d_%1d.dat",ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
      VecView(user->P_sum0, viewer);
      PetscViewerDestroy(&viewer);
      sprintf(filen, "spph0_%06d_%1d.dat.info",ti, user->_this);	if(!rank) unlink(filen);
 
      PetscBarrier(PETSC_NULL);

      sprintf(filen, "spph1_%06d_%1d.dat",ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
      VecView(user->P_sum1, viewer);
      PetscViewerDestroy(&viewer);
      sprintf(filen, "spph1_%06d_%1d.dat.info",ti, user->_this);	if(!rank) unlink(filen);
 
      PetscBarrier(PETSC_NULL);

      sprintf(filen, "spph2_%06d_%1d.dat",ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
      VecView(user->P_sum2, viewer);
      PetscViewerDestroy(&viewer);
      sprintf(filen, "spph2_%06d_%1d.dat.info",ti, user->_this);	if(!rank) unlink(filen);
 
      PetscBarrier(PETSC_NULL);

      sprintf(filen, "spph3_%06d_%1d.dat",ti, user->_this);
      PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
      VecView(user->P_sum3, viewer);
      PetscViewerDestroy(&viewer);
      sprintf(filen, "spph3_%06d_%1d.dat.info",ti, user->_this);	if(!rank) unlink(filen);
 
      PetscBarrier(PETSC_NULL);

     }




  }
  
  if(les) {
    Vec Cs;
    
    VecDuplicate(user->P, &Cs);
    DMLocalToGlobalBegin(user->da, user->lCs, INSERT_VALUES, Cs);
    DMLocalToGlobalEnd(user->da, user->lCs, INSERT_VALUES, Cs);
    
    sprintf(filen, "cs_%06d_%1d.dat",  ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(Cs, viewer);
    PetscViewerDestroy(&viewer);
    
    PetscBarrier(PETSC_NULL);
    VecDestroy(&Cs);
  }
  
  if(rans) {
    sprintf(filen, "kfield%06d_%1d.dat", ti, user->_this);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_WRITE, &viewer);
    VecView(user->K_Omega, viewer);
    PetscViewerDestroy(&viewer);
    
    PetscBarrier(PETSC_NULL);
  }
  
  return 0;
}

PetscErrorCode Divergence(UserCtx *user )
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  Vec		Div;
  PetscReal	***div, ***aj, ***nvert,***p;
  Cmpnts	***ucont;
  PetscReal	maxdiv;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda,user->lUcont, &ucont);
  DMDAVecGetArray(da, user->lAj, &aj);
  VecDuplicate(user->P, &Div);
  DMDAVecGetArray(da, Div, &div);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	maxdiv = fabs((ucont[k][j][i].x - ucont[k][j][i-1].x +
		       ucont[k][j][i].y - ucont[k][j-1][i].y +
		       ucont[k][j][i].z - ucont[k-1][j][i].z)*aj[k][j][i]);
	if (nvert[k][j][i] + nvert[k+1][j][i] + nvert[k-1][j][i] +
	    nvert[k][j+1][i] + nvert[k][j-1][i] +
	    nvert[k][j][i+1] + nvert[k][j][i-1] > 0.1) maxdiv = 0.;
	div[k][j][i] = maxdiv;

      }
    }
  }

  if (zs==0) {
    k=0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ze == mz) {
    k=mz-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xs==0) {
    i=0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xe==mx) {
    i=mx-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0;
      }
    }
  }

  if (ys==0) {
    j=0;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ye==my) {
    j=my-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }
  DMDAVecRestoreArray(da, Div, &div);
  VecMax(Div, &i, &maxdiv);
  PetscPrintf(PETSC_COMM_WORLD, "Maxdiv %d %d %e\n", ti, i, maxdiv);
  PetscInt mi;
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (mi=xs; mi<xe; mi++) {
	if (Gidx(mi,j,k,user) ==i) {
	  PetscPrintf(PETSC_COMM_SELF, "MMa %d %d %d\n", mi,j, k);
	}
      }
    }
  }
 

 DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lUcont, &ucont);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  VecDestroy(&Div);
  return(0);
}


PetscErrorCode force(UserCtx *user, double psum[97], double plsum[97], double pu17sum[97], double pl17sum[97], double pu24sum[97], double pl24sum[97], double puavesum[97], double plavesum[97])
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  PetscReal	 ***nvert,***p;
 
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  int rank;
  
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
 

        DMDAVecGetArray(da, user->lNvert, &nvert);
        DMDAVecGetArray(da, user->P, &p);
	 double plSum[97]; double pSum[97]; double pu[97]; double pl[97]; double plave[97]; double puave[97]; double plaveSum[97]; double puaveSum[97];
	 double pl17[97]; double pu17[97]; double pu17Sum[97]; double pl17Sum[97];  
	 double pl24[97]; double pu24[97]; double pu24Sum[97]; double pl24Sum[97];  
	 for (i=0;i<97;i++) 
		{
 			       plSum[i]=0.0; pSum[i] = 0.0; pl[i]=0.; pu[i]=0.;puave[i]=0.; plave[i]=0.; puaveSum[i]=0.; plaveSum[i]=0.;
				pu24[i]=0.; pl24[i]=0.; pu24Sum[i]=0.; pl24Sum[i]=0.;
				pu17[i]=0.; pl17[i]=0.; pu17Sum[i]=0.; pl17Sum[i]=0.;
		}

	 
	int ks=71,ke=168,count;
	for (k=zs;k<ze;k++){
	 count=0;
	 for (j=ys;j<ye;j++){
	  for (i=xs;i<xe;i++){
	 
	   if ( k>ks && k<ke ){
		if (i==41 && nvert[k][j][i]==1.0  && nvert[k][j+1][i]<1.0 )
		  pu[k-ks-1]+=p[k][j+1][i];

	
		if (i==41 && nvert[k][j][i]==1. && nvert[k][j-1][i]<1.0 )
		  pl[k-ks-1]+=p[k][j-1][i];
         	
		if (nvert[k][j][i]==1.0 && nvert[k][j-1][i]>2.0 && nvert[k][j+1][i]<2.0 )
	 	   puave[k-ks-1]+=p[k][j][i];

		if ( nvert[k][j][i]==1. && nvert[k][j-1][i]<2.0 && nvert[k][j+1][i]>2.0)
		  { plave[k-ks-1]+=p[k][j][i];  }


		if (i==24 && nvert[k][j][i]==1.0 && nvert[k][j-1][i]>2.0 && nvert[k][j+1][i]<2.0 )
		  pu24[k-ks-1]+=p[k][j][i];

		if (i==24 && nvert[k][j][i]==1. && nvert[k][j-1][i]<2.0 && nvert[k][j+1][i]>2.0)
		  pl24[k-ks-1]+=p[k][j][i];


		if (i==17 && nvert[k][j][i]==1.0 && nvert[k][j-1][i]>2.0 && nvert[k][j+1][i]<2.0 )
		  pu17[k-ks-1]+=p[k][j][i];

		if (i==17 && nvert[k][j][i]==1. && nvert[k][j-1][i]<2.0 && nvert[k][j+1][i]>2.0)
		  pl17[k-ks-1]+=p[k][j][i];

            }
	    
	   }   
	  }
		
	 }
	
	  MPI_Allreduce(&pu, &pSum,97,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
	  MPI_Allreduce(&pl, &plSum,97,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

	  MPI_Allreduce(&pu24, &pu24Sum,97,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
	  MPI_Allreduce(&pl24, &pl24Sum,97,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
	 
	  MPI_Allreduce(&pu17, &pu17Sum,97,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
	  MPI_Allreduce(&pl17, &pl17Sum,97,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
	 
	  MPI_Allreduce(&puave, &puaveSum,97,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
	  MPI_Allreduce(&plave, &plaveSum,97,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
	 


 
	    PetscPrintf(PETSC_COMM_WORLD, " instantanous forcces\n");
	  

	for (i=0;i<97;i++) 
	{
	   psum[i]+=pSum[i];
	   plsum[i]+=plSum[i];
	   pu24sum[i]+=pu24Sum[i];
	   pl24sum[i]+=pl24Sum[i];

	   pu17sum[i]+=pu17Sum[i];
	   pl17sum[i]+=pl17Sum[i];

	   puavesum[i]+=puaveSum[i]/54;
	   plavesum[i]+=plaveSum[i]/54;

	   PetscPrintf(PETSC_COMM_WORLD, " %d  %6f  %6f \n",i,pSum[i],plSum[i]-pSum[i]);
	}

 
  DMDAVecRestoreArray(da, user->P, &p);
  DMDAVecGetArray(da, user->lNvert, &nvert);

}






PetscErrorCode GridDivergence(UserCtx *user)
{
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  Vec		Div;
  PetscReal	***div, ***aj, ***nvert;
  Cmpnts	***ucont;
  PetscReal	maxdiv;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal    norm;
  VecNorm(user->Vcont, NORM_INFINITY, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "Grid Flux norm  %le\n", norm);

  DMDAVecGetArray(fda,user->lVcont, &ucont);
  DMDAVecGetArray(da, user->lAj, &aj);
  VecDuplicate(user->P, &Div);
  DMDAVecGetArray(da, Div, &div);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	maxdiv = fabs((ucont[k][j][i].x - ucont[k][j][i-1].x +
		       ucont[k][j][i].y - ucont[k][j-1][i].y +
		       ucont[k][j][i].z - ucont[k-1][j][i].z)*aj[k][j][i]);
	if (nvert[k][j][i] + nvert[k+1][j][i] + nvert[k-1][j][i] +
	    nvert[k][j+1][i] + nvert[k][j-1][i] +
	    nvert[k][j][i+1] + nvert[k][j][i-1] > 0.1) maxdiv = 0.;
	div[k][j][i] = maxdiv;
      }
    }
  }

  if (zs==0) {
    k=0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ze == mz) {
    k=mz-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xs==0) {
    i=0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (xe==mx) {
    i=mx-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	div[k][j][i] = 0;
      }
    }
  }

  if (ys==0) {
    j=0;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }

  if (ye==my) {
    j=my-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	div[k][j][i] = 0.;
      }
    }
  }
  DMDAVecRestoreArray(da, Div, &div);
  VecMax(Div, &i, &maxdiv);
  PetscPrintf(PETSC_COMM_WORLD, "Maxdiv Grid %d %d %e\n", ti, i, maxdiv);
  PetscInt mi;
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (mi=xs; mi<xe; mi++) {
	if (Gidx(mi,j,k,user) ==i) {
	  PetscPrintf(PETSC_COMM_SELF, "Maxdiv Grid location %d %d %d\n", mi,j, k);
	}
      }
    }
  }
    
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lVcont, &ucont);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  VecDestroy(&Div);
  return(0);
}

//-----------------------------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv) {

  UserCtx    *user;
  PetscInt   i, bi, ibi;  
  IBMNodes   *ibm;
  IBMVNodes  *ibmv;
  FSInfo     *fsi;
  PetscBool  DoSCLoop;
  PetscInt   itr_sc;
  Cstart     cstart;
  PetscInt   level;
  UserMG     usermg;
  PetscBool  flg;
  FE         *fem;

  PetscInitialize(&argc, &argv, (char *)0, help);

  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);
  PetscOptionsGetInt(PETSC_NULL, "-tio", &tiout, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-phave", &phave, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-imm", &immersed, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-inv", &invicid, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rstart", &tistart, &rstart_flg);
  PetscOptionsGetInt(PETSC_NULL, "-rstart_fem", &rstart_fem, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-imp", &implicit, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-imp_type", &implicit_type, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-imp_MAX_IT", &imp_MAX_IT, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fsi", &movefsi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rfsi", &rotatefsi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-radi", &radi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-inlet", &inletprofile, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-str", &STRONG_COUPLING, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rs_fsi", &rstart_fsi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-cop", &cop, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fish", &fish, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-pizza", &pizza, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rheology", &rheology, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-duplicate", &duplicate, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-Pipe", &Pipe, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-turbine", &turbine, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fishcyl", &fishcyl, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-eel", &eel, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-cstart", &fish_c, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-wing", &wing, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-sediment", &sediment, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-mhv", &MHV, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-bhv", &BHV, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-beam", &beam, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-inv_flg", &inv_flg, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-hydro", &hydro, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-lv", &LV, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-reg", &regime, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-twoD", &TwoD, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-thin", &thin, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_z", &dgf_z, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_y", &dgf_y, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_x", &dgf_x, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_az", &dgf_az, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_ay", &dgf_ay, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-dgf_ax", &dgf_ax, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-body", &NumberOfBodies, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-nbody", &nbody, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-mframe", &moveframe, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rframe", &rotateframe, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-blk", &blank, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-init1", &InitialGuessOne, PETSC_NULL);

  PetscOptionsGetReal(PETSC_NULL, "-x_c", &(CMx_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-y_c", &(CMy_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-z_c", &(CMz_c), PETSC_NULL);

  PetscOptionsGetReal(PETSC_NULL, "-imp_atol", &(imp_atol), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-imp_rtol", &(imp_rtol), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-imp_stol", &(imp_stol), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-poisson_tol", &poisson_tol, PETSC_NULL);		// Seokkoo Kang: tolerance of implicit matrix free solver. 1.e-4 is enough for most cases.

  PetscOptionsGetInt(PETSC_NULL, "-les", &les, PETSC_NULL);				// Seokkoo Kang: if 1 Smagorinsky with Cs=0.1, if 2 Dynamic model
  PetscOptionsGetInt(PETSC_NULL, "-rans", &rans, PETSC_NULL);			// Seokkoo Kang
  PetscOptionsGetInt(PETSC_NULL, "-wallfunction", &wallfunction, PETSC_NULL);	// Seokkoo Kang: 1 or 2
  PetscOptionsGetInt(PETSC_NULL, "-mixed", &mixed, PETSC_NULL);			// Seokkoo Kang: mixed model option for LES
  PetscOptionsGetInt(PETSC_NULL, "-clark", &clark, PETSC_NULL);			// Seokkoo Kang: mixed model option for LES
  PetscOptionsGetInt(PETSC_NULL, "-dynamic_freq", &dynamic_freq, PETSC_NULL);		// Seokkoo Kang: LES dynamic compute frequency 
  PetscOptionsGetInt(PETSC_NULL, "-averaging", &averaging, PETSC_NULL);	// Seokkoo Kang: if 1 do averaging; always begin with -rstart 0

  PetscOptionsGetInt(PETSC_NULL, "-grid1d", &grid1d, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-i_periodic", &i_periodic, PETSC_NULL);	
  PetscOptionsGetInt(PETSC_NULL, "-j_periodic", &j_periodic, PETSC_NULL);	
  PetscOptionsGetInt(PETSC_NULL, "-k_periodic", &k_periodic, PETSC_NULL);	
  PetscOptionsGetInt(PETSC_NULL, "-period", &period, PETSC_NULL);

  PetscOptionsGetInt(PETSC_NULL, "-pbc_domain", &blkpbc, PETSC_NULL);

  PetscOptionsGetInt(PETSC_NULL, "-testfilter_ik", &testfilter_ik, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-testfilter_1d", &testfilter_1d, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-poisson", &poisson, PETSC_NULL);
  /*PetscOptionsGetInt(PETSC_NULL, "-i_homo_filter", &i_homo_filter, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-j_homo_filter", &j_homo_filter, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-k_homo_filter", &k_homo_filter, PETSC_NULL);
  */
  PetscOptionsGetReal(PETSC_NULL, "-max_cs", &max_cs, PETSC_NULL);
  
  PetscOptionsGetReal(PETSC_NULL, "-St_exp", &(St_exp), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-wlngth", &(wavelength), PETSC_NULL);
    
  PetscPrintf(PETSC_COMM_WORLD, "tiout %d %le %le thin %d!\n",tiout, imp_atol,imp_rtol,thin);

  if (immersed) {
    ibm=(IBMNodes *)calloc(NumberOfBodies,sizeof(IBMNodes));
    ibmv=(IBMVNodes *)calloc(NumberOfBodies,sizeof(IBMVNodes));
    L_dim = 1; // Channel flow external 
 }

MG_Initial(&usermg, ibm);
PetscPrintf(PETSC_COMM_WORLD,"block number: %d \n",block_number);
  if (immersed) {
    level = usermg.mglevels-1;
    user = usermg.mgctx[level].user;
    for (bi=0; bi<block_number; bi++) {
      user[bi].ibmlist=(IBMList *)malloc(NumberOfBodies*sizeof(IBMList));
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	InitIBMList(&(user[bi].ibmlist[ibi]));
      }
      }
  for(i=0;i<NumberOfBodies;i++){
	PetscPrintf(PETSC_COMM_WORLD,"IBM Read\n");
	ibm_read_Icem(&ibm[i],i);
	ibm_surface_VTKOut(&ibm[i],i,0);
	PetscBarrier(PETSC_NULL);
}
} // Immersed loop ends
  
  if (immersed) {
    ti = 0;
    if (rstart_flg) ti = tistart;
  }
  
  level = usermg.mglevels-1;
  user = usermg.mgctx[level].user;
  if (rstart_flg) {
    ti = tistart; tistart++;
    
    for (bi=0; bi<block_number; bi++) {
      Ucont_P_Binary_Input(&(user[bi]));
      
      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont,INSERT_VALUES, user[bi].lUcont);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont,INSERT_VALUES, user[bi].lUcont);
      
      DMGlobalToLocalBegin(user[bi].da, user[bi].P,INSERT_VALUES, user[bi].lP);
      DMGlobalToLocalEnd(user[bi].da, user[bi].P,INSERT_VALUES, user[bi].lP);
      
      Contra2Cart(&(user[bi]));
  	}
   } //rstart_flg
    
  // do the search once if elmt is not moving!
  if (immersed) {
    for (level = usermg.mglevels-1; level>=usermg.mglevels-1; level--) {
      user = usermg.mgctx[level].user;
      for (bi=0; bi<block_number; bi++) {
// SEARCH
	user[bi].ibm=ibm;
	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	      PetscPrintf(PETSC_COMM_WORLD, "IBM_SERA_REV ibi %d bi %d\n", ibi, bi);
	      ibm_search_advanced_rev(&(user[bi]), &ibm[ibi], ibi);  
	    }//ibi


// INTERPOLATION

	PetscBarrier(PETSC_NULL);
	for (ibi=0;ibi<NumberOfBodies;ibi++) {
	  PetscPrintf(PETSC_COMM_WORLD, "IBM_INTP\n");
	  ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi, 1);
	  PetscPrintf(PETSC_COMM_WORLD, "IBM_INTP End\n");
	  
	} //ibi

      } //bi
    } //mglevels
 } //immersed
  
  ti = 0;
  if (rstart_flg) ti = tistart;
 // PetscPrintf(PETSC_COMM_WORLD," init1: %d\n",InitialGuessOne);
 // PetscPrintf(PETSC_COMM_WORLD,"rstart_flg: %d\n",rstart_flg);
//  PetscPrintf(PETSC_COMM_WORLD," ti: %d\n",ti);
  if (ti==0) {
    PetscPrintf(PETSC_COMM_WORLD," init1: %d\n",InitialGuessOne);
    if (InitialGuessOne) {
      for (bi=0; bi<block_number; bi++) {
	SetInitialGuessToOne(&(user[bi]));
	PetscPrintf(PETSC_COMM_WORLD,"Initial Guess Set to One\n");
	InflowFlux(&(user[bi]));
	OutflowFlux(&(user[bi]));
	FormBCS(&(user[bi]));
      }
      for (bi=0; bi<block_number; bi++) {
	
	PetscReal normZ, normX, normY;
	VecStrideMax(user[bi].Ucat, 0, PETSC_NULL, &normX);
	VecStrideMax(user[bi].Ucat, 1, PETSC_NULL, &normY);
	VecStrideMax(user[bi].Ucat, 2, PETSC_NULL, &normZ);
	PetscPrintf(PETSC_COMM_WORLD, "Initial Eq 11111111111! %le %le %le\n",normX, normY, normZ);
      }
    }
  }
  
  // Copy Ucont to Ucont_o for the finest level
  for (bi=0; bi<block_number; bi++) {
      VecCopy(user[bi].Ucont, user[bi].Ucont_o);
      VecCopy(user[bi].Ucont, user[bi].Ucont_rm1);
      VecCopy(user[bi].Ucat, user[bi].Ucat_o);
      VecCopy(user[bi].P, user[bi].P_o);
    
    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont);
    
    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
    
    DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
    DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
  }
  
  PetscInt tisteps = 2;
  PetscOptionsGetInt(PETSC_NULL, "-totalsteps", &tisteps, &flg);
  
  if (tistart==0) tisteps ++;
  /*  put the time accuracy coefficient 1 for the 1st real-time step */
  /*   COEF_TIME_ACCURACY=1.; */ 
  /* ==================================================================================             */
  /*   physical time Step Loop */
  
  double psum[97]; double plsum[97]; double pu17sum[97]; double pl17sum[97]; double pu24sum[97]; double pl24sum[97]; double puavesum[97]; double plavesum[97];
  for (i=0;i<97;i++) { 
    psum[i]=0.; plsum[i]=0.; pu17sum[i]=0.; pl17sum[i]=0.; pu24sum[i]=0.; pl24sum[i]=0.;puavesum[i]=0.; plavesum[i]=0.;
  }
  for (ti = tistart; ti<tistart + tisteps; ti++) {
    
    PetscPrintf(PETSC_COMM_WORLD, "Time %d\n", ti); 
    /* ==================================================================================             */
    /*     Strong-Coupling (SC) Loop */
    DoSCLoop= PETSC_TRUE ; itr_sc = 0;
    while (DoSCLoop) {
      
      itr_sc++;
      PetscPrintf(PETSC_COMM_WORLD, "SC LOOP itr # %d\n", itr_sc);
      
      //Structral Solver!
      if (immersed){
     //	Struc_Solver(&usermg,ibm,fsi,&cstart, itr_sc,tistart, &DoSCLoop, fem);
    // else
	DoSCLoop = PETSC_FALSE;      
      /*     Flow Solver! */
       Flow_Solver(&usermg, ibm, fsi);
	}
    }// End of while SC loop
    /* ==================================================================================             */
    /*  put the time accuracy coefficient back to 1.5
	after the 1st real-time step */
    /*     COEF_TIME_ACCURACY=1.5; */            
    /* ==================================================================================             */
    /*     Save the old values (at ti) for later */
    PetscPrintf(PETSC_COMM_WORLD," Time loop ends\n");    
    level = usermg.mglevels-1;
    user = usermg.mgctx[level].user;
    for (bi=0; bi<block_number; bi++) {  
     
     if (immersed) {
	VecCopy(user[bi].Nvert, user[bi].Nvert_o);
	DMGlobalToLocalBegin(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
	DMGlobalToLocalEnd(user[bi].da, user[bi].Nvert_o, INSERT_VALUES, user[bi].lNvert_o);
      }
            
      // Copy Ucont to Ucont_o for the finest level
      VecCopy(user[bi].Ucat, user[bi].Ucat_o);
      
      VecCopy(user[bi].Ucont_o, user[bi].Ucont_rm1);
      VecCopy(user[bi].Ucont, user[bi].Ucont_o);
      VecCopy(user[bi].P, user[bi].P_o);
      
      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_o, INSERT_VALUES, user[bi].lUcont_o);
      
      DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
      DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont_rm1, INSERT_VALUES, user[bi].lUcont_rm1);
      

    } // for bi    
    
    if (immersed) {
      
      for (ibi=0;ibi<NumberOfBodies;ibi++) {
	
	for (i=0; i<ibm[ibi].n_v; i++) {
	  ibm[ibi].x_bp_o[i] = ibm[ibi].x_bp[i];
	  ibm[ibi].y_bp_o[i] = ibm[ibi].y_bp[i];
	  ibm[ibi].z_bp_o[i] = ibm[ibi].z_bp[i];
	  
	  ibm[ibi].urm1[i].x = ibm[ibi].uold[i].x;
	  ibm[ibi].urm1[i].y = ibm[ibi].uold[i].y;
	  ibm[ibi].urm1[i].z = ibm[ibi].uold[i].z;
	  
	  ibm[ibi].uold[i].x = ibm[ibi].u[i].x;
	  ibm[ibi].uold[i].y = ibm[ibi].u[i].y;
	  ibm[ibi].uold[i].z = ibm[ibi].u[i].z;
	}
	
	      } //ibi
    } // if immersed
    
    /* ==================================================================================             */
    
  } // ti (physical time) loop
  /* ==================================================================================             */
  //  nextp:   PetscPrintf(PETSC_COMM_WORLD, "Number of iteration is  %d \n",ti);
  
  MG_Finalize(&usermg);
  PetscFinalize();
  
  /* ==================================================================================             */
  return(0); 
}
