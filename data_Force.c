static char help[] = "Testing programming!";

#include "petscda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"
#include "petscksp.h"
#include "petscvec.h"

//  fscanf(f, "%le %le %le \n", &(FSinf->M_x),&(FSinf->M_y), &(FSinf->M_z));	  
typedef struct {
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];
  PetscReal    S_ang_n[6],S_ang_o[6],S_ang_r[6],S_ang_rm1[6];
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z, A_tot; //Forces & Area
  PetscReal    M_x,M_y,M_z; // Moments
  PetscReal    x_c,y_c,z_c;
} FSInfo;


PetscErrorCode Force_DATA_Input(FSInfo *fsi, PetscInt ibi, 
				PetscInt tisteps)
{
  PetscInt rank,bi=0,tii;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Force_Coeff_SI%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "r");
    for (tii=0;tii<tisteps;tii++) {
      fscanf(f, "%d %le %le %le %le %le %le %le %le %le\n",ti, F_xSum, F_ySum, F_zSum,Cp_xSum,Cp_ySum,Cp_zSum, A_xSum,A_ySum,A_zSum);     
    }
    fclose(f);

    sprintf(filen, "Momt_Coeff_SI%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "r");
    fscanf(f, "%d %le %le %le %le %le %le\n",ti, M_xSum, M_ySum, M_zSum,FSinfo->M_x_old,FSinfo->M_x_real,FSinfo->Mdpdn_x);
    fclose(f);

    sprintf(filen, "Fish_Energetic_SI%2.2d_%2.2d",ibi,bi);
    f = fopen(filen, "r");
    fscanf(f," %d %le %le %le %le %le %le %le\n",ti, efficiency,Pw_sideSum,Pw_xSum,Pw_ySum,Pw_zSum,F_zSum,A_tSum);
    fclose(f);
  }
}
PetscErrorCode Force_DATA_output(FSInfo *fsi, PetscInt ibi)
{
  PetscInt rank,bi=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "Force_Coeff_SI%2.2d_%2.2d_cleaned.dat",ibi,bi);
    f = fopen(filen, "w");
    
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le\n",ti, F_xSum, F_ySum, F_zSum,Cp_xSum,Cp_ySum,Cp_zSum, A_xSum,A_ySum,A_zSum);
    fclose(f);

    sprintf(filen, "Momt_Coeff_SI%2.2d_%2.2d_cleaned.dat",ibi,bi);
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le\n",ti, M_xSum, M_ySum, M_zSum,FSinfo->M_x_old,FSinfo->M_x_real,FSinfo->Mdpdn_x);
    fclose(f);

    sprintf(filen, "Fish_Energetic_SI%2.2d_%2.2d_cleaned.dat",ibi,bi);
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f," %d %le %le %le %le %le %le %le\n",ti, efficiency,Pw_sideSum,Pw_xSum,Pw_ySum,Pw_zSum,F_zSum,A_tSum);
    fclose(f);
  }
}

int main(int argc, char **argv)
{
  PetscInt      ibi,tisteps=2000;
  FSInfo        *fsi;
  
  PetscInitialize(&argc, &argv, (char *)0, help);

  PetscErrorCode flg;

  PetscOptionsGetInt(PETSC_NULL, "-body", &numberOfbodies, &flg);
  if (!flg)
    numberOfbodies=0;

  PetscOptionsGetInt(PETSC_NULL, "-tm", &tisteps, &flg);

  PetscMalloc(tisteps*sizeof(FSInfo), &fsi);

  //  PetscOptionsGetInt(PETSC_NULL, "-tie", &tiend, &flg);

  for (ibi=0;ibi<=numbreOfbodies;ibi++) {
    Force_DATA_Input(fsi, ibi);
    Force_DATA_output(fsi, ibi);
  }

  PetscFinalize();

  return(0);
}
