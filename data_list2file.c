static char help[] = "Testing programming!";

#include "petscdmda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"
#include "petscksp.h"
#include "petscvec.h"

typedef struct {
  PetscReal    S_new[6],S_old[6],S_real[6],S_realm1[6];
  PetscReal    S_ang_n[6],S_ang_o[6],S_ang_r[6],S_ang_rm1[6];
  PetscReal    red_vel, damp, mu_s; // reduced vel, damping coeff, mass coeff
  PetscReal    F_x,F_y,F_z, A_tot; //Forces & Area
  PetscReal    M_x,M_y,M_z; // Moments
  PetscReal    x_c,y_c,z_c;
} FSInfo;

typedef struct {
  PetscReal    Pw_x,Pw_y,Pw_z,F_z,A_tot;
} ENGinfo;

typedef struct {
  PetscReal    F_x, F_y, F_z,Cp_x,Cp_y,Cp_z, A_x,A_y,A_z;
} FRCinfo;

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

PetscErrorCode FSI_DATA_output(FSInfo *FSinfo, PetscInt ti,
			       PetscInt ibi)
{
  PetscInt rank, i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  PetscBarrier(PETSC_NULL);
  if (!rank) {
    FILE *f;
    char filen[80];
/*     sprintf(filen, "FSI_position%2.2d",ibi); */
/*     f = fopen(filen, "a"); */
/* /\*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, S_new[4],S_new[5],F_z, S_real[4], S_real[5], S_realm1[4], S_realm1[5]); *\/ */
/* /\*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, FSinfo->S_new[2],FSinfo->S_new[3],FSinfo->F_y, FSinfo->S_real[2], FSinfo->S_real[3], FSinfo->S_realm1[2], FSinfo->S_realm1[3]); *\/ */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le\n",ti, FSinfo->S_new[2],FSinfo->S_new[3],FSinfo->F_y, FSinfo->S_new[4],FSinfo->S_new[5],FSinfo->F_z,FSinfo->S_new[0],FSinfo->S_new[1],FSinfo->F_x); */
/*     fclose(f); */

/* /\*     if (dgf_z) { *\/ */
/* /\*     sprintf(filen, "FSI_position_x"); *\/ */
/* /\*     f = fopen(filen, "a"); *\/ */
/* /\*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, FSinfo->S_new[4],FSinfo->S_new[5],FSinfo->F_z, FSinfo->S_real[4], FSinfo->S_real[5], FSinfo->S_realm1[4], FSinfo->S_realm1[5]); *\/ */
/* /\* /\\*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, FSinfo->S_new[2],FSinfo->S_new[3],FSinfo->F_y, FSinfo->S_real[2], FSinfo->S_real[3], FSinfo->S_realm1[2], FSinfo->S_realm1[3]); *\\/ *\/ */
/* /\*     fclose(f); *\/ */
/* /\*     } *\/ */

/*     sprintf(filen, "FSI_Agnle%2.2d",ibi); */
/*     f = fopen(filen, "a"); */
/*     PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le\n",ti, FSinfo->S_ang_n[0],FSinfo->S_ang_n[1],FSinfo->M_x,FSinfo->S_ang_r[0],FSinfo->S_ang_r[1],FSinfo->S_ang_rm1[0],FSinfo->S_ang_rm1[1]); */
/*     fclose(f); */

    sprintf(filen, "DATA_FSI%5.5d_%2.2d.dat",ti, ibi);
    f = fopen(filen, "w");
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->red_vel, FSinfo->damp, FSinfo->mu_s);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->x_c, FSinfo->y_c, FSinfo->z_c);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->F_x, FSinfo->F_y, FSinfo->F_z);	  
    PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le \n", FSinfo->M_x, FSinfo->M_y, FSinfo->M_z);	  
    for (i=0; i<6; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_new[i],FSinfo->S_old[i], FSinfo->S_real[i], FSinfo->S_realm1[i]);
      PetscFPrintf(PETSC_COMM_WORLD, f, "%le %le %le %le \n", FSinfo->S_ang_n[i],FSinfo->S_ang_o[i], FSinfo->S_ang_r[i], FSinfo->S_ang_rm1[i]);
    }
    fclose(f);
  }
  return(0);
}


/* PetscErrorCode FSI_DATA_output(FSInfo *FSinf, PetscInt ti) */
/* { */
/*   FILE *f; */
/*   char filen[80]; */
/*   sprintf(filen, "FSI_data"); */
/*   f = fopen(filen, "a"); */
/*   PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le\n",ti,FSinf->S_new[4],FSinf->S_new[5],FSinf->F_z, FSinf->S_ang_n[0], FSinf->S_ang_n[1],FSinf->M_x); */
/*   fclose(f); */
/*   return(0); */
/* } */

PetscErrorCode ENG_DATA_output(ENGinfo *eng, PetscInt ti)
{
  FILE *f;
  char filen[80];
  //  sprintf(filen, "ENG_data");
  sprintf(filen, "DATA_ENG%5.5d.dat",ti);
  f = fopen(filen, "w");
  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le\n",ti,eng->Pw_x,eng->Pw_y,eng->Pw_z,eng->F_z,eng->A_tot);
  //  PetscFPrintf(PETSC_COMM_WORLD, f, "END of File\n");

  fclose(f);
  return(0);
}

PetscErrorCode FRC_DATA_output(FRCinfo *frc, PetscInt ti)
{
  FILE *f;
  char filen[80];
  //  sprintf(filen, "ENG_data");
  sprintf(filen, "DATA_FRC%5.5d.dat",ti);
  f = fopen(filen, "w");
  PetscFPrintf(PETSC_COMM_WORLD, f, "%d %le %le %le %le %le %le %le %le %le\n",ti,frc->F_x,frc->F_y,frc->F_z,frc->Cp_x,frc->Cp_y,frc->Cp_z,frc->A_x,frc->A_y,frc->A_z);
  //  PetscFPrintf(PETSC_COMM_WORLD, f, "END of File\n");

  fclose(f);
  return(0);
}

int main(int argc, char **argv)
{
  PetscInt      ti,tistart,tiend,tii;
  FSInfo        FSI;
  ENGinfo       eng;
  FRCinfo       frc;
  PetscReal     xxx;
         
  PetscInitialize(&argc, &argv, (char *)0, help);

  PetscBool flg;

  PetscOptionsGetInt(PETSC_NULL, "-tis", &tistart, &flg);
  PetscOptionsGetInt(PETSC_NULL, "-tie", &tiend, &flg);
  PetscBool ENG = PETSC_FALSE;
  PetscOptionsGetBool(PETSC_NULL, "-eng", &ENG, PETSC_NULL);
  PetscBool FRC = PETSC_FALSE;
  PetscOptionsGetBool(PETSC_NULL, "-frc", &FRC, PETSC_NULL);

  PetscInt rank,bi=0,ibi=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];

    if (FRC) {
      sprintf(filen, "Force_Coeff_SI%2.2d_%2.2d",ibi,bi);
    } else if (ENG) {
    
      sprintf(filen, "Power_SI%2.2d_%2.2d",ibi,bi);
    } else {
      sprintf(filen, "FSI_position%2.2d",ibi);
    }
   
    f = fopen(filen, "r");
    for (ti=tistart;ti<=tiend;ti++) {
      //    FSI_DATA_Input(&FSI, ti);
      if (FRC) {
	//	fscanf(f, "%d %le %le %le %le %le %le %le %le %le\n",tii, F_xSum, F_ySum, F_zSum,Cp_xSum,Cp_ySum,Cp_zSum, A_xSum,A_ySum,A_zSum);
	fscanf(f, "%d %le %le %le %le %le %le %le %le %le\n",&tii,&frc.F_x,&frc.F_y,&frc.F_z,&frc.Cp_x,&frc.Cp_y,&frc.Cp_z,&frc.A_x,&frc.A_y,&frc.A_z);
	PetscPrintf(PETSC_COMM_WORLD, "ti, %d %e \n",tii, frc.F_z);
	FRC_DATA_output(&frc, tii);	
      } else if (ENG) {
	fscanf(f," %d %le %le %le %le %le %le %le\n",&tii, &xxx,&xxx,&eng.Pw_x,&eng.Pw_y,&eng.Pw_z,&eng.F_z,&eng.A_tot);
	PetscPrintf(PETSC_COMM_WORLD, "ti, %d %e \n",tii, eng.Pw_x);
	ENG_DATA_output(&eng, tii);
      } else {
	fscanf(f, "%d %le %le %le %le %le %le %le %le %le\n",&tii, &FSI.S_new[2],&FSI.S_new[3],&FSI.F_y, &FSI.S_new[4],&FSI.S_new[5],&FSI.F_z,&FSI.S_new[0],&FSI.S_new[1],&FSI.F_x);
	PetscPrintf(PETSC_COMM_WORLD, "ti, %d %e \n",tii, FSI.S_new[4]);
	FSI_DATA_output(&FSI, tii, ibi);
      }

    }
    fclose(f);
  }
  PetscFinalize();

  return(0);
}
