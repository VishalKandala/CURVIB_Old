static char help[] = "Testing programming!";

#include "variables.h"

#include "petscda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"
#include "petscksp.h"
#include "petscvec.h"

PetscInt ti,tistart=0;
PetscReal	Flux_in = 4.104388e-04, angle = 0;
PetscInt tiout = 10;
PetscInt block_number;
PetscReal FluxInSum, FluxOutSum;
PetscInt immersed = 0;
PetscInt invicid = 0;
PetscInt movefsi = 0, rotatefsi=0;
PetscInt implicit = 0;
PetscInt imp_MAX_IT = 50; 
PetscInt radi=10;
PetscInt inletprofile=2;
PetscReal CMx_c=0., CMy_c=0., CMz_c=0.;
PetscInt  mg_MAX_IT=30, mg_idx=1, mg_preItr=1, mg_poItr=1;
PetscReal imp_atol=1e-7, imp_rtol=1e-4, imp_stol=1.e-8;
PetscInt TwoD = 1;
PetscInt STRONG_COUPLING=0;
PetscInt rstart_fsi=0;
PetscInt cop=0, regime=1; // 1 escape regime --- 2 cruise regime
PetscInt fish=0, eel=0;
PetscReal St_exp=0.5,wavelength=0.95;
PetscInt MHV=0;
PetscReal max_angle = -54.*3.1415926/180.;// for MHV; min=0
PetscInt thin=0;
PetscInt dgf_z=0,dgf_y=1,dgf_x=0;
PetscInt dgf_az=0,dgf_ay=0,dgf_ax=1 ;
PetscInt NumberOfBodies=1;
PetscInt moveframe=0, wing=0;
PetscReal L_dim=1;

//PetscErrorCode FSI_DATA_Input(FSInfo *FSinf, PetscInt ti)
PetscErrorCode FSI_DATA_Input(FSInfo *FSinf, PetscInt ibi, PetscInt ti)
{
  PetscInt  i;
  PetscReal t;

  FILE *f;
  char filen[80];  
  sprintf(filen, "DATA_FSI%5.5d_%2.2d.dat",ti, ibi);
  f = fopen(filen, "r");
  if (!f) {
    SETERRQ(PETSC_COMM_WORLD,1, "Cannot open FSI DATA file");
    PetscPrintf(PETSC_COMM_WORLD, "FSI_data cannot open file !!!!!!!!!!!!\n");
  }
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input begin %d %s\n",ti,filen);

  fscanf(f, "%le %le %le", &t, &t, &t);	  
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input red vel damp mu %le %le %le \n",FSinf->red_vel,FSinf->damp,FSinf->mu_s);
  fscanf(f, "%le %le %le", &(FSinf->x_c), &(FSinf->y_c), &(FSinf->z_c));	  
  fscanf(f, "%le %le %le \n", &(FSinf->F_x),&(FSinf->F_y), &(FSinf->F_z));	 
  fscanf(f, "%le %le %le \n", &(FSinf->M_x),&(FSinf->M_y), &(FSinf->M_z));	  
 

  for (i=0; i<6; i++) {
    fscanf(f, "%le %le %le %le", &(FSinf->S_new[i]),&(FSinf->S_old[i]), &(FSinf->S_real[i]), &(FSinf->S_realm1[i]));
    fscanf(f, "%le %le %le %le", &(FSinf->S_ang_n[i]),&(FSinf->S_ang_o[i]), &(FSinf->S_ang_r[i]), &(FSinf->S_ang_rm1[i]));
  }
  if (rheology){
    fscanf(f, "%le %le %le ", &FSinf->R[0][0],&FSinf->R[0][1],&FSinf->R[0][2]);
    fscanf(f, "%le %le %le ", &FSinf->R[1][0],&FSinf->R[1][1],&FSinf->R[1][2]);
    fscanf(f, "%le %le %le ", &FSinf->R[2][0],&FSinf->R[2][1],&FSinf->R[2][2]);
    fscanf(f, "%le %le %le ", &FSinf->L_n[0],&FSinf->L_n[1],&FSinf->L_n[2]);
  }
  fclose(f);
  if (rheology)  PetscPrintf(PETSC_COMM_WORLD,"%le %le %le %le %le %le \n",FSinf->R[0][0],FSinf->R[0][1],FSinf->R[0][2],FSinf->R[1][1],FSinf->R[1][2],FSinf->R[2][2]); 
  else{
    PetscPrintf(PETSC_COMM_WORLD, "FSI_data input z, dz/dt  %le %le %le %le\n",FSinf->S_new[4],FSinf->S_new[5],FSinf->red_vel,FSinf->damp);
    PetscPrintf(PETSC_COMM_WORLD, "FSI_data input y, dy/dt  %le %le %le %le\n",FSinf->S_new[2],FSinf->S_new[3],FSinf->red_vel,FSinf->damp);
  }
 
  return(0);
}

//PetscErrorCode FSI_DATA_Input(FSInfo *FSinf, PetscInt ti)
PetscErrorCode FSI_DATA_Input2(FSInfo *FSinf, PetscInt ibi)
{
  PetscInt  i;
  //PetscReal 

  FILE *f;
  char filen[80];  
  //  sprintf(filen, "DATA_FSI%5.5d.dat",ti);
  sprintf(filen, "../TEST_lamprel_R4000St0.6/DATA_FSI%5.5d_%2.2d.dat",ti, ibi);

  f = fopen(filen, "r");
  if (!f) {
    SETERRQ(1, "Cannot open DATA_FSI file");
    PetscPrintf(PETSC_COMM_WORLD, "FSI_data cannot open file !!!!!!!!!!!!\n");
  }
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input begin %d %s\n",ti,filen);
  fscanf(f, "%le %le %le", &(FSinf->red_vel), &(FSinf->damp), &(FSinf->mu_s));	  
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input red vel damp  %le %le\n",FSinf->red_vel,FSinf->damp);
  fscanf(f, "%le %le %le", &(FSinf->x_c), &(FSinf->y_c), &(FSinf->z_c));	  
  fscanf(f, "%le %le %le \n", &(FSinf->F_x),&(FSinf->F_y), &(FSinf->F_z));	 
  fscanf(f, "%le %le %le \n", &(FSinf->M_x),&(FSinf->M_y), &(FSinf->M_z));	  

  for (i=0; i<6; i++) {
    fscanf(f, "%le %le %le %le", &(FSinf->S_new[i]),&(FSinf->S_old[i]), &(FSinf->S_real[i]), &(FSinf->S_realm1[i]));
    fscanf(f, "%le %le %le %le", &(FSinf->S_ang_n[i]),&(FSinf->S_ang_o[i]), &(FSinf->S_ang_r[i]), &(FSinf->S_ang_rm1[i]));
  }
  fclose(f);
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input z, dz/dt  %le %le %le %le\n",FSinf->S_new[4],FSinf->S_new[5],FSinf->red_vel,FSinf->damp);
  PetscPrintf(PETSC_COMM_WORLD, "FSI_data input angle_x, dang_x/dt  %le %le %le %le\n",FSinf->S_ang_n[0],FSinf->S_ang_n[1],FSinf->red_vel,FSinf->damp);
  return(0);
}

PetscErrorCode surface_VTKOut(IBMNodes *ibm, PetscInt ti)
{
    // vtk file name
  PetscInt n_cells=3;
  PetscInt rank,i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "surf_tail%5.5d.vtk",ti);
    f = fopen(filen, "w"); // open file

    PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 2.0\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "Surface Grid\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "ASCII\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET UNSTRUCTURED_GRID\n");
   
    //    PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  5523993 double\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  %d float\n",ibm->n_v);
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->x_bp[i],ibm->y_bp[i],ibm->z_bp[i]);
    }

    PetscFPrintf(PETSC_COMM_WORLD, f, "CELLS %d %d\n",ibm->n_elmt, (n_cells+1)*ibm->n_elmt);
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD,f, "%d  %d %d %d\n",n_cells, ibm->nv1[i], ibm->nv2[i], ibm->nv3[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_TYPES %d\n",ibm->n_elmt);
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD,f, "%d\n",5);
    }
  }
  return(0);
}

int main(int argc, char **argv)
{
  PetscInt      tiend,tistep=10,i;//ti,tistart,
  FSInfo        fsi,fsi2;
  IBMNodes      ibm;
  UserCtx	user;
 
  PetscInitialize(&argc, &argv, (char *)0, help);

  PetscErrorCode flg;

  //  PetscOptionsInsertFile("control.dat");
  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_NULL);

  PetscOptionsGetInt(PETSC_NULL, "-tis", &tistart, &flg);
  PetscOptionsGetInt(PETSC_NULL, "-tie", &tiend, &flg);
  PetscOptionsGetInt(PETSC_NULL, "-ts", &tistep, &flg);
  PetscOptionsGetInt(PETSC_NULL, "-tio", &tiout, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-imm", &immersed, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-inv", &invicid, PETSC_NULL);
  //  PetscOptionsGetInt(PETSC_NULL, "-rstart", &tistart, &rstart_flg);
  PetscOptionsGetInt(PETSC_NULL, "-imp", &implicit, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-imp_MAX_IT", &imp_MAX_IT, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fsi", &movefsi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rfsi", &rotatefsi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-radi", &radi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-inlet", &inletprofile, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-str", &STRONG_COUPLING, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-rs_fsi", &rstart_fsi, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-cop", &cop, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-fish", &fish, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-eel", &eel, PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL, "-mhv", &MHV, PETSC_NULL);
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
  PetscOptionsGetInt(PETSC_NULL, "-mframe", &moveframe, PETSC_NULL);

  PetscOptionsGetReal(PETSC_NULL, "-x_c", &(CMx_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-y_c", &(CMy_c), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-z_c", &(CMz_c), PETSC_NULL);

  PetscOptionsGetReal(PETSC_NULL, "-imp_atol", &(imp_atol), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-imp_rtol", &(imp_rtol), PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL, "-imp_stol", &(imp_stol), PETSC_NULL);

  if (fish) {
    PetscOptionsGetReal(PETSC_NULL, "-St_exp", &(St_exp), PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, "-wlngth", &(wavelength), PETSC_NULL);
  }


  //ibm_read_fish(&ibm,15.);
  ibm_read_Ansys(&ibm,0, 28.8);
  fish_init(&(user.dt));

  for (ti=tistart;ti<tiend;ti=ti+tistep) {
  if (ti == (ti/tiout)*tiout) {

    FSI_DATA_Input(&fsi, 0);
    fish_swim(&ibm,(ti)*user.dt, user.dt);
    //    FSI_DATA_Input2(&fsi2, 0);
    for (i=0; i<ibm.n_v; i++) {
      // ibm.x_bp[i] = ibm.x_bp[i] + 1.;
      ibm.z_bp[i] = ibm.z_bp0[i] + fsi.S_new[4];
    }

    surface_VTKOut(&ibm, ti);
  }
  }
  PetscFinalize();

  return(0);
}
