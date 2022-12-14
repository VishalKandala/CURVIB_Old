static char help[] = "Testing programming!";

#include "variables.h"

#ifdef TECIO
#include "TECIO.h"
#endif

PetscInt NumberOfBodies=1;
PetscReal CMx_c=0., CMy_c=0., CMz_c=0.;
PetscErrorCode ibm_read_tecplot(IBMNodes *ibm, PetscInt nt, PetscInt ti)
{
PetscInt	rank;
PetscInt	n_v , n_elmt ;
PetscInt	i;

PetscPrintf(PETSC_COMM_WORLD, "IBM Reads TECPLOT print out \n ");

  //MPI_Comm_size(PETSC_COMM_WORLD, &size);
MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

// n_v=8; n_elmt=12;

 if(!rank) 
  { // root processor read in the data
    FILE *fd;
    char filen[128],ch[128];

    sprintf(filen, "surface_nf%3.3d_%2.2d.dat",ti,nt);
    fd = fopen(filen, "r");
    // fgets(ch, 128, fd);

    //    fscanf(fd, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z\n");
    fscanf(fd, "%54s\n",ch);//%54s
    //fscanf(fd, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", &n_v, &n_elmt);
    fscanf(fd, "%4s %14s %2s %d %s %2s %d %s %10s %12s %45s\n",ch, ch, ch, &n_v, ch,ch, &n_elmt
	   , ch, ch, ch, ch);
    PetscPrintf(PETSC_COMM_WORLD," nv = %d ne = %d ...\n",n_v, n_elmt);

    if (!fd) 
    {
      PetscPrintf(PETSC_COMM_WORLD, "Cannot open surface file %d\n", nt);
      SETERRQ(PETSC_COMM_WORLD,1, "Cannot open IBM node file surfacexxx_xx.dat!");
     }

    if (fd) 
     {
     
       ibm->n_v=n_v;
       ibm->n_elmt=n_elmt;

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));

      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));
      
      PetscReal cl = 1.,p;
      PetscOptionsGetReal(PETSC_NULL, "-chact_leng_valve", &cl, PETSC_NULL);
/*       cl=1./L_dim; */

      PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->nv3);

      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nf_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nf_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nf_z);

      // Added 4/1/06 iman
      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->dA); //Area

      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nt_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nt_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nt_z);

/*       PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x); */
/*       PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y); */
/*       PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z); */
      
      // Added 6/4/06 iman
      //PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->cent));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));
      // end added


      for (i=0; i<n_v; i++) {
	fscanf(fd, "%le\n", &ibm->x_bp[i]);
      }
      for (i=0; i<n_v; i++) {
	fscanf(fd, "%le\n", &ibm->y_bp[i]);
      }
      for (i=0; i<n_v; i++) {	
	fscanf(fd, "%le\n", &ibm->z_bp[i]);
      }
      for (i=0; i<n_v; i++) {
	fscanf(fd, "%le\n", &ibm->u[i].x);
      }
      for (i=0; i<n_v; i++) {
	fscanf(fd, "%le\n", &ibm->u[i].y);
      }
      for (i=0; i<n_v; i++) {
	fscanf(fd, "%le\n", &ibm->u[i].z);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->nf_x[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->nf_y[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->nf_z[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->nt_x[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->nt_y[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->nt_z[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%d %d %d\n", &ibm->nv1[i], &ibm->nv2[i], &ibm->nv3[i]);
      }
      fclose(fd);
      
      for (i=0; i<n_elmt; i++) {
	ibm->nv1[i]--;
	ibm->nv2[i]--;
	ibm->nv3[i]--;
      }
      
     }
  }


  PetscPrintf(PETSC_COMM_WORLD, "IBM Reads TECPLOT successfully \n ");

  return(0);

}



PetscErrorCode ibm_read_tecplot2(IBMNodes *ibm, PetscInt nt, PetscInt ti)
{
PetscInt	rank;
PetscInt	n_v , n_elmt ;
PetscInt	i;

PetscPrintf(PETSC_COMM_WORLD, "IBM Reads TECPLOT print out \n ");

  //MPI_Comm_size(PETSC_COMM_WORLD, &size);
MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

// n_v=8; n_elmt=12;

 if(!rank) 
  { // root processor read in the data
    FILE *fd;
    char filen[128],ch[128];

    sprintf(filen, "surface%3.3d_%2.2d.dat",ti,nt);
    // sprintf(filen, "surface%3.3d.dat",ti);
    fd = fopen(filen, "r");
    // fgets(ch, 128, fd);

    //    fscanf(fd, "Variables=x,y,z,n_x,n_y,n_z,nt_x,nt_y,nt_z,ns_x,ns_y,ns_z\n");
    fscanf(fd, "%57s\n",ch);//%54s
    //fscanf(fd, "ZONE T='TRIANGLES', N=%d, E=%d, F=FEBLOCK, ET=TRIANGLE, VARLOCATION=([1-3]=NODAL,[4-12]=CELLCENTERED)\n", &n_v, &n_elmt);
    fscanf(fd, "%4s %14s %2s %d %s %2s %d %s %10s %12s %45s\n",ch, ch, ch, &n_v, ch,ch, &n_elmt
	   , ch, ch, ch, ch);
    PetscPrintf(PETSC_COMM_WORLD," nv = %d ne = %d ...\n",n_v, n_elmt);

    if (!fd) 
    {
      PetscPrintf(PETSC_COMM_WORLD, "Cannot open surface file %d\n", nt);
      SETERRQ(PETSC_COMM_WORLD,1, "Cannot open IBM node file surfacexxx_xx.dat!");
     }

    if (fd) 
     {
     
       ibm->n_v=n_v;
       ibm->n_elmt=n_elmt;

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp));

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp_o));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp_o));

      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->x_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->y_bp0));
      PetscMalloc(n_v*sizeof(PetscReal), &(ibm->z_bp0));

      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->u));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->uold));
      PetscMalloc(n_v*sizeof(Cmpnts), &(ibm->urm1));
      
      PetscReal cl = 1.,p;
      PetscOptionsGetReal(PETSC_NULL, "-chact_leng_valve", &cl, PETSC_NULL);
/*       cl=1./L_dim; */

      PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->nv1);
      PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->nv2);
      PetscMalloc(n_elmt*sizeof(PetscInt), &ibm->nv3);

      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nf_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nf_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nf_z);

      // Added 4/1/06 iman
      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->dA); //Area

      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nt_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nt_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->nt_z);

      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->ns_x);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->ns_y);
      PetscMalloc(n_elmt*sizeof(PetscReal), &ibm->ns_z);
/*       PetscMalloc(n_elmt*sizeof(PetscReal), &ns_x); */
/*       PetscMalloc(n_elmt*sizeof(PetscReal), &ns_y); */
/*       PetscMalloc(n_elmt*sizeof(PetscReal), &ns_z); */
      
      // Added 6/4/06 iman
      //PetscMalloc(n_elmt*sizeof(Cmpnts), &(ibm->cent));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_x));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_y));
      PetscMalloc(n_elmt*sizeof(PetscReal), &(ibm->cent_z));
      // end added


      for (i=0; i<n_v; i++) {
	fscanf(fd, "%le\n", &ibm->x_bp[i]);
      }
      for (i=0; i<n_v; i++) {
	fscanf(fd, "%le\n", &ibm->y_bp[i]);
      }
      for (i=0; i<n_v; i++) {	
	fscanf(fd, "%le\n", &ibm->z_bp[i]);
      }
/*       for (i=0; i<n_v; i++) { */
/* 	fscanf(fd, "%le\n", &ibm->u[i].x); */
/*       } */
/*       for (i=0; i<n_v; i++) { */
/* 	fscanf(fd, "%le\n", &ibm->u[i].y); */
/*       } */
/*       for (i=0; i<n_v; i++) {	 */
/* 	fscanf(fd, "%le\n", &ibm->u[i].z); */
/*       } */
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->nf_x[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->nf_y[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->nf_z[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->nt_x[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->nt_y[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->nt_z[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->ns_x[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->ns_y[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%le\n", &ibm->ns_z[i]);
      }
      for (i=0; i<n_elmt; i++) {
	fscanf(fd, "%d %d %d\n", &ibm->nv1[i], &ibm->nv2[i], &ibm->nv3[i]);
      }
      fclose(fd);
      
      for (i=0; i<n_elmt; i++) {
	ibm->nv1[i]--;
	ibm->nv2[i]--;
	ibm->nv3[i]--;
      }
      
     }
  }


  PetscPrintf(PETSC_COMM_WORLD, "IBM Reads TECPLOT successfully \n ");

  return(0);

}

PetscErrorCode PolyVTKOut(IBMNodes *ibm, PetscInt ibi, PetscInt ti)
{
    // vtk file name
  PetscInt n_cells=3;
  PetscInt rank,i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "surface%2.2d_%5.5d.vtk", ibi,ti);
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

    PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n", ibm->n_v);
    PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS u float\n");
    for (i=0; i<ibm->n_v; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->u[i].x,ibm->u[i].y,ibm->u[i].z);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_DATA %d\n", ibm->n_elmt);
    PetscFPrintf(PETSC_COMM_WORLD, f, "NORMALS nf float\n");
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->nf_x[i],ibm->nf_y[i],ibm->nf_z[i]);
    }
    PetscFPrintf(PETSC_COMM_WORLD, f, "NORMALS nt float\n");
    for (i=0; i<ibm->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", ibm->nt_x[i],ibm->nt_y[i],ibm->nt_z[i]);
    }

  }
  return(0);
}

PetscErrorCode PolyVTKOut2(IBMNodes *ibm, PetscInt ibi, PetscInt ti)
{
    // vtk file name
  PetscInt n_cells=3;
  PetscInt rank,i;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (!rank) {
    FILE *f;
    char filen[80];
    sprintf(filen, "surface_nf%2.2d_%5.5d.vtk", ibi,ti);
    f = fopen(filen, "w"); // open file

    PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 2.0\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "Surface Grid\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "ASCII\n");
    PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET UNSTRUCTURED_GRID\n");
   
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

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
  PetscInt      ibi;
  IBMNodes	*ibm;

  PetscInitialize(&argc, &argv, (char *)0, help);
  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);

  PetscBool flag;
  // PetscErrorCode flag;
  PetscInt ti, tis, tie, tsteps=5;
  PetscOptionsGetInt(PETSC_NULL, "-tis", &tis, &flag);
  if (!flag) {
    PetscPrintf(PETSC_COMM_WORLD, "Need the starting number!\n");
  }

  PetscOptionsGetInt(PETSC_NULL, "-tie", &tie,&flag);
  if (!flag) {
    tie = tis;
    //    PetscPrintf(PETSC_COMM_WORLD, "Need the ending number!\n");
  }

  PetscOptionsGetInt(PETSC_NULL, "-ts", &tsteps, &flag);
  if (!flag) {
    tsteps = 10; /* Default increasement is 5 */
    //    PetscPrintf(PETSC_COMM_WORLD, "Need the ending number!\n");
  }

  PetscMalloc(NumberOfBodies*sizeof(IBMNodes), &ibm);

  for (ti=tis; ti<=tie; ti+=tsteps) {

    for (ibi=0;ibi<NumberOfBodies;ibi++) {
      
      ibm_read_tecplot2(&ibm[ibi], ibi, ti);   
      
      PetscPrintf(PETSC_COMM_WORLD, "ibm read tec %d\n", ibi); 
      
      PolyVTKOut2(&ibm[ibi], ibi, ti);   
      // calc_p_force(&ibm[ibi]);
    }
  }

  PetscFinalize();
}
