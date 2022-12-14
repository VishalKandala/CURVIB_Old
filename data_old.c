static char help[] = "Testing programming!";

#include "petscda.h"
#include "petscts.h"
#include "petscpc.h"
#include "petscsnes.h"

PetscInt ti, block_number;

typedef struct {
  PassiveScalar u, v, w, p;
} PassiveField;

typedef struct {
  PetscScalar u, v, w;
} Field;

typedef struct {
  PetscScalar x, y, z;
} Cmpnts;

typedef struct {
  PassiveScalar csi[3], eta[3], zet[3], aj;
} Metrics;

typedef struct {
  Vec	Ubcs; // An ugly hack, waste of memory
} BCS;



typedef struct {
  PetscInt	IM, JM, KM; // dimensions of grid
  DA da;	/* Data structure for scalars (include the grid geometry
		   informaion, to obtain the grid information, 
		   use DAGetCoordinates) */
  DA fda;	// Data Structure for vectors
  DALocalInfo info;

  Vec	Cent;	// Coordinates of cell centers
  Vec 	Csi, Eta, Zet, Aj;
  Vec 	ICsi, IEta, IZet, IAj;
  Vec 	JCsi, JEta, JZet, JAj;
  Vec 	KCsi, KEta, KZet, KAj;
  Vec 	Ucont;	// Contravariant velocity components
  Vec 	Ucat;	// Cartesian velocity components
  Vec	Ucat_o;
  Vec 	P;
  Vec	Phi;
  Vec	GridSpace;
  Vec	Nvert;
  Vec	Nvert_o;
  BCS	Bcs;

  PetscInt	*nvert;//ody property

  PetscReal	ren;	// Reynolds number
  PetscReal	dt; 	// time step
  PetscReal	st;	// Strouhal number

  PetscReal	r[101], tin[101], uinr[101][1001];

  Vec	lUcont, lUcat, lP, lPhi;
  Vec	lCsi, lEta, lZet, lAj;
  Vec	lICsi, lIEta, lIZet, lIAj;
  Vec	lJCsi, lJEta, lJZet, lJAj;
  Vec	lKCsi, lKEta, lKZet, lKAj;
  Vec	lGridSpace;
  Vec	lNvert, lNvert_o;
  Vec	lCent;

  PetscInt this;
} UserCtx;

PetscErrorCode SetCoordinates3d(UserCtx *user)
{
  /* Set 3d grid coordinates. */
  DA	cda = user->fda, da = user->da;
  Vec		Coor;
  Cmpnts	***coor;
  PetscInt	IM, JM, KM;
  PetscInt	rank;
  
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	i, j, k;
  PetscReal	*gc;
  PetscReal	d0 = 1.;

  DASetUniformCoordinates(da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

  DAGetGhostedCoordinates(da, &Coor);
  DAVecGetArray(cda, Coor, &coor);

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if(!rank) {
    FILE *fd;
    fd = fopen("grid.dat", "r");

    fscanf(fd, "%i\n", &block_number);
    MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    fclose(fd);
  }
  else {
    MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
  }
  PetscMalloc(block_number*sizeof(UserCtx), &user);

  ReadCoordinates(user);

/*   DASetUniformCoordinates(da, -0.5, 0.5, -0.5, 0.5, 0.0, 6.0); */

/*   for (k=zs+1; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	coor[k][j][i].x = coor[zs][j][i].x;//1. / (mx-2) * i; */
/* 	coor[k][j][i].y = coor[zs][j][i].y;//1. / (my-2) * j; */
/* 	coor[k][j][i].z = coor[zs][j][i].z - 10. / (mz-2) * k;// 1. / (mz-2) * k; */
/*       } */
/*     } */
/*   } */


  Vec	gCoor;
  DAGetCoordinates(da, &gCoor);
  DALocalToGlobal(cda, Coor, INSERT_VALUES, gCoor);

  DAGlobalToLocalBegin(cda, gCoor, INSERT_VALUES, Coor);
  DAGlobalToLocalEnd(cda, gCoor, INSERT_VALUES, Coor);

  return(0);
}

PetscErrorCode Contra2Cart(UserCtx *user)
{
  DA		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal	mat[3][3], det, det0, det1, det2;

  Vec		Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;
  Vec		Aj = user->lAj;
  Vec		Ucont = user->lUcont, Ucat = user->Ucat;
  Vec		Ubcs = user->Bcs.Ubcs;
  Cmpnts	***csi, ***eta, ***zet;
  PetscReal	***aj;
  Cmpnts	***ucont, ***ucat, ***ubcs;

  PetscReal	***nvert;
  PetscReal	coef = 0.125;

  PetscReal	q[3]; //local working array
  PetscInt	i, j, k;
  PetscReal	c1, c2, c3, c4, g1, g2, g3, g4, ucon;
  Cmpnts	***gs;
  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;

  if (lxs==0) lxs++;
  if (lxe==mx) lxe--;
  if (lys==0) lys++;
  if (lye==my) lye--;
  if (lzs==0) lzs++;
  if (lze==mz) lze--;

  DAVecGetArray(fda, Csi, &csi);
  DAVecGetArray(fda, Eta, &eta);
  DAVecGetArray(fda, Zet, &zet);
  DAVecGetArray(da,  Aj,  &aj);

  DAVecGetArray(fda, Ucont, &ucont);
  DAVecGetArray(fda, Ucat,  &ucat);
  DAVecGetArray(fda, Ubcs,  &ubcs);
/*   DAVecGetArray(fda, user->lGridSpace, &gs); */

  DAVecGetArray(da, user->lNvert, &nvert);

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((int)(nvert[k][j][i]+0.5)==0) {
	  mat[0][0] = 0.5 * (csi[k][j][i-1].x + csi[k][j][i].x);
	  mat[0][1] = 0.5 * (csi[k][j][i-1].y + csi[k][j][i].y);
	  mat[0][2] = 0.5 * (csi[k][j][i-1].z + csi[k][j][i].z);
	         	      
	  mat[1][0] = 0.5 * (eta[k][j-1][i].x + eta[k][j][i].x);
	  mat[1][1] = 0.5 * (eta[k][j-1][i].y + eta[k][j][i].y);
	  mat[1][2] = 0.5 * (eta[k][j-1][i].z + eta[k][j][i].z);
	         	      
	  mat[2][0] = 0.5 * (zet[k-1][j][i].x + zet[k][j][i].x);
	  mat[2][1] = 0.5 * (zet[k-1][j][i].y + zet[k][j][i].y);
	  mat[2][2] = 0.5 * (zet[k-1][j][i].z + zet[k][j][i].z);

/* 	  q[0] = 0.5 * (ucont[k][j][i-1].x + ucont[k][j][i].x); */
/* 	  q[1] = 0.5 * (ucont[k][j-1][i].y + ucont[k][j][i].y); */
/* 	  q[2] = 0.5 * (ucont[k-1][j][i].z + ucont[k][j][i].z); */

	  ucon = 0.5 * (ucont[k][j][i-1].x + ucont[k][j][i].x);
	  if (ucon  > 0) {
	    if (i>1 && (int)(nvert[k][j][i-2]+0.5)==0) {
/* 	      c1 = 0.5 * gs[k][j][i].x; */
/* 	      c2 = 0.5 * gs[k][j][i].x + gs[k][j][i-1].x; */
/* 	      c3 = gs[k][j][i].x; */
/* 	      c4 = gs[k][j][i].x + gs[k][j][i-1].x; */

/* 	      g1 = c1 * c2 / (c3 * c4); */

/* 	      c2 = gs[k][j][i].x * 0.5; */
/* 	      c3 = gs[k][j][i-1].x; */

/* 	      g2 = c1 * c2 / (c3 * c4); */

/* 	      q[0] = (ucont[k][j][i-1].x + */
/* 		      (ucont[k][j][i].x - ucont[k][j][i-1].x) * g1 + */
/* 		      (ucont[k][j][i-1].x - ucont[k][j][i-2].x) * g2); */

	      q[0] = coef * (-    ucont[k][j][i-2].x -
			     2. * ucont[k][j][i-1].x +
			     3. * ucont[k][j][i  ].x) + ucont[k][j][i-1].x;
	    }
	    else {
	      q[0] = 0.5 * (ucont[k][j][i-1].x + ucont[k][j][i].x);
	    }
	  }
	  else {
	    if (i < mx-2 &&(int)(nvert[k][j][i+1]+0.5)==0) {
/* 	      c1 = -0.5 * gs[k][j][i].x; */
/* 	      c2 = -(0.5 * gs[k][j][i].x + gs[k][j][i+1].x); */
/* 	      c3 = -gs[k][j][i+1].x; */
/* 	      c4 = -(gs[k][j][i].x + gs[k][j][i+1].x); */

/* 	      g3 = c1 * c2 / (c3 * c4); */

/* 	      c2 = -(0.5 * gs[k][j][i].x); */
/* 	      c3 = -gs[k][j][i+1].x; */

/* 	      g4 = c1 * c2 / (c3 * c4); */

/* 	      q[0] = (ucont[k][j][i].x + */
/* 		      (ucont[k][j][i-1].x - ucont[k][j][i].x) * g3 + */
/* 		      (ucont[k][j][i].x   - ucont[k][j][i+1].x) * g4); */
	      q[0] = coef * (-    ucont[k][j][i+1].x -
			     2. * ucont[k][j][i  ].x +
			     3. * ucont[k][j][i-1].x) + ucont[k][j][i].x;
	    }
	    else {
	      q[0] = 0.5 * (ucont[k][j][i-1].x + ucont[k][j][i].x);
	    }
	  }

/* 	  q[0] = ucon; */
	  ucon = 0.5 * (ucont[k][j-1][i].y + ucont[k][j][i].y);
	  if (ucon > 0) {
	    if (j>1&&(int)(nvert[k][j-2][i]+0.5)==0) {
/* 	      c1 = 0.5 * gs[k][j][i].y; */
/* 	      c2 = 0.5 * gs[k][j][i].y + gs[k][j-1][i].y; */
/* 	      c3 = gs[k][j][i].y; */
/* 	      c4 = gs[k][j][i].y + gs[k][j-1][i].y; */

/* 	      g1 = c1 * c2 / (c3 * c4); */

/* 	      c2 = gs[k][j][i].y * 0.5; */
/* 	      c3 = gs[k][j-1][i].y; */

/* 	      g2 = c1 * c2 / (c3 * c4); */

/* 	      q[1] = (ucont[k][j-1][i].y + */
/* 		      (ucont[k][j][i].y - ucont[k][j-1][i].y) * g1 + */
/* 		      (ucont[k][j-1][i].y - ucont[k][j-2][i].y) * g2); */
	      q[1] = coef * (-    ucont[k][j-2][i].y -
			     2. * ucont[k][j-1][i].y +
			     3. * ucont[k][j  ][i].y) + ucont[k][j-1][i].y;
	    }
	    else {
	      q[1] = 0.5 * (ucont[k][j-1][i].y + ucont[k][j][i].y);
	    }
	  }
	  else {
	    if (j < my-2 &&(int)(nvert[k][j+1][i]+0.5)==0) {
/* 	      c1 = -0.5 * gs[k][j][i].y; */
/* 	      c2 = -(0.5 * gs[k][j][i].y + gs[k][j+1][i].y); */
/* 	      c3 = -gs[k][j+1][i].y; */
/* 	      c4 = -(gs[k][j][i].y + gs[k][j+1][i].y); */

/* 	      g3 = c1 * c2 / (c3 * c4); */

/* 	      c2 = -(0.5 * gs[k][j][i].y); */
/* 	      c3 = -gs[k][j+1][i].y; */

/* 	      g4 = c1 * c2 / (c3 * c4); */

/* 	      q[1] = (ucont[k][j][i].y + */
/* 		      (ucont[k][j-1][i].y - ucont[k][j][i].y) * g3 + */
/* 		      (ucont[k][j][i].y   - ucont[k][j+1][i].y) * g4); */
	      q[1] = coef * (-    ucont[k][j+1][i].y -
			     2. * ucont[k][j  ][i].y +
			     3. * ucont[k][j-1][i].y) + ucont[k][j][i].y;
	    }
	    else {
	      q[1] = 0.5 * (ucont[k][j-1][i].y + ucont[k][j][i].y);
	    }
	  }
/* 	  q[1]=ucon; */
	  ucon = 0.5 * (ucont[k-1][j][i].z + ucont[k][j][i].z);
	  if (ucon > 0) {
	    if (k>1&&(int)(nvert[k-2][j][i]+0.5)==0) {
/*  	      c1 = 0.5 * gs[k][j][i].z; */
/* 	      c2 = 0.5 * gs[k][j][i].z + gs[k-1][j][i].z; */
/* 	      c3 = gs[k][j][i].z; */
/* 	      c4 = gs[k][j][i].z + gs[k-1][j][i].z; */

/* 	      g1 = c1 * c2 / (c3 * c4); */

/* 	      c2 = gs[k][j][i].z * 0.5; */
/* 	      c3 = gs[k-1][j][i].z; */

/* 	      g2 = c1 * c2 / (c3 * c4); */

/* 	      q[2] = (ucont[k-1][j][i].z + */
/* 		      (ucont[k][j][i].z - ucont[k-1][j][i].z) * g1 + */
/* 		      (ucont[k-1][j][i].z - ucont[k-2][j][i].z) * g2); */
	      q[2] = coef * (-    ucont[k-2][j][i].z -
			     2. * ucont[k-1][j][i].z +
			     3. * ucont[k  ][j][i].z) + ucont[k-1][j][i].z;
	    }
	    else {
	      q[2] = 0.5 * (ucont[k-1][j][i].z + ucont[k][j][i].z);
	    }
	  }
	  else {
	    if (k < mz-2 &&(int)(nvert[k+1][j][i]+0.5)==0) {
/* 	      c1 = -0.5 * gs[k][j][i].z; */
/* 	      c2 = -(0.5 * gs[k][j][i].z + gs[k+1][j][i].z); */
/* 	      c3 = -gs[k+1][j][i].z; */
/* 	      c4 = -(gs[k][j][i].z + gs[k+1][j][i].z); */

/* 	      g3 = c1 * c2 / (c3 * c4); */

/* 	      c2 = -(0.5 * gs[k][j][i].z); */
/* 	      c3 = -gs[k+1][j][i].z; */

/* 	      g4 = c1 * c2 / (c3 * c4); */

/* 	      q[2] = (ucont[k][j][i].z + */
/* 		      (ucont[k-1][j][i].z - ucont[k][j][i].z) * g3 + */
/* 		      (ucont[k][j][i].z   - ucont[k+1][j][i].z) * g4); */
	      q[2] = coef * (-    ucont[k+1][j][i].z -
			     2. * ucont[k  ][j][i].z +
			     3. * ucont[k-1][j][i].z) + ucont[k][j][i].z;
	    }
	    else {
	      q[2] = 0.5 * (ucont[k-1][j][i].z + ucont[k][j][i].z);
	    }
	  }
/* 	  q[2]= ucon; */

	  det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
	    mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
	    mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

	  det0 = q[0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
	    q[1] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) +
	    q[2] * (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]);

	  det1 = -q[0] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
	    q[1] * (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) -
	    q[2] * (mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]);

	  det2 = q[0] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) -
	    q[1] * (mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]) +
	    q[2] * (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);

/* 	  q[2] = q[2] / mat[2][2]; */

/* 	  q[1] = (q[1] - q[2] * mat[1][2]) / mat[1][1]; */

/* 	  q[0] = (q[0] - q[2] * mat[0][2] - q[1] * mat[0][1]) / mat[0][0]; */

/* 	  if (j==6 && (i==1 || i==2)) { */
/* 	    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e %d\n", q[0], q[1], q[2], i); */
/* 	  } */
	
	  /*	  if ((i==21)&& k==1) {
	    PetscPrintf(PETSC_COMM_WORLD, "Ve %e %e %e\n", q[0], q[1], q[2]);
	    }*/
	  ucat[k][j][i].x = det0 / det;
	  ucat[k][j][i].y = det1 / det;
	  ucat[k][j][i].z = det2 / det;

/* 	  if (k==40 && j==17 && (i==39|| i==xe-40)) { */
/* 	    PetscPrintf(PETSC_COMM_SELF, "Ucat %le %le %le\n", ucat[k][j][i].x, ucat[k][j][i].y, ucat[k][j][i].z); */
/* 	    PetscPrintf(PETSC_COMM_SELF, "Ucont %le %le\n", ucont[k][j][i].z, ucont[k-1][j][i].z); */
/* 	  } */

/* 	  if (k==1) PetscPrintf(PETSC_COMM_WORLD, "%le %le\n",ucat[k][j][i].x * mat[2][0] + ucat[k][j][i].y * mat[2][1] + ucat[k][j][i].z* mat[2][2], q[2]); */
/* 	  if (rank) PetscPrintf(PETSC_COMM_SELF, "%d %d %d %le %le %le\n", i, j, k, det0/det, det1/det, det2/det); */
	}
      }
    }
  }

  /*  for (k=xs; k<xe; k++) {
    for (i=xs; i<xe; i++) {
      if (ucont[k][j][i].y>0)
	PetscPrintf(PETSC_COMM_WORLD, "%d %d %d %e vbc\n", i,j,k,ucont[k][j][i].y);
    }
  }
  PetscPrintf(PETSC_COMM_WORLD, "Cart End\n");*/

/*   for (k=zs; k<ze; k++) { */
/*     for (i=xs; i<xe; i++) { */
/*       j=0; */
/*       ubcs[k][j][i].x = ucat[k][ye-2][i].x; */
/*       ubcs[k][j][i].y = ucat[k][ye-2][i].y; */
/*       ubcs[k][j][i].z = ucat[k][ye-2][i].z; */
/*       j=ye-1; */
/*       ubcs[k][j][i].x = ucat[k][1][i].x; */
/*       ubcs[k][j][i].y = ucat[k][1][i].y; */
/*       ubcs[k][j][i].z = ucat[k][1][i].z; */
/*     } */
/*   } */

/*   for (k=zs; k<ze; k++) { */
/*     for (i=ys; j<ye; j++) { */
/*       i=0; */
/*       ubcs[k][j][i].x = ucat[k][j][xe-2].x; */
/*       ubcs[k][j][i].y = ucat[k][j][xe-2].y; */
/*       ubcs[k][j][i].z = ucat[k][j][xe-2].z; */
/*       i=xe-1; */
/*       ubcs[k][j][i].x = ucat[k][j][1].x; */
/*       ubcs[k][j][i].y = ucat[k][j][1].y; */
/*       ubcs[k][j][i].z = ucat[k][j][1].z; */
/*     } */
/*   } */

  DAVecRestoreArray(fda, Csi, &csi);
  DAVecRestoreArray(fda, Eta, &eta);
  DAVecRestoreArray(fda, Zet, &zet);
  DAVecRestoreArray(da,  Aj,  &aj);

  DAVecRestoreArray(fda, Ucont, &ucont);
  DAVecRestoreArray(fda, Ucat,  &ucat);
  DAVecRestoreArray(fda, Ubcs,  &ubcs);

  VecAssemblyBegin(user->Ucat);
  VecAssemblyEnd(user->Ucat);

  DAGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DAGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  DAVecRestoreArray(da, user->lNvert, &nvert);
  //  PetscPrintf(PETSC_COMM_WORLD, "Ve End\n");
  return(0);
}

PetscErrorCode TecOut(UserCtx *user)
{
  DA		da = user->da, fda = user->fda;
  DALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	lxs, lys, lzs, lxe, lye, lze;
  PetscInt	i, j, k;

  PetscReal	***aj;
  Cmpnts	***ucat, ***coor, ***ucat_o;
  PetscReal	***p, ***nvert;
  Vec		Coor, zCoor, nCoor;
  VecScatter	ctx;

  FILE	*f;
  char filen[80];

  PetscInt	tt=0;


  DAGetCoordinates(da, &Coor);

/*   DACreateNaturalVector(fda, &nCoor); */
/*   DAGlobalToNaturalBegin(fda, Coor, INSERT_VALUES, nCoor); */
/*   DAGlobalToNaturalEnd(fda, Coor, INSERT_VALUES, nCoor); */

/*   VecScatterCreateToZero(nCoor, &ctx, &zCoor); */
/*   VecScatterBegin(nCoor, zCoor, INSERT_VALUES, SCATTER_FORWARD, ctx); */
/*   VecScatterEnd(nCoor, zCoor, INSERT_VALUES, SCATTER_FORWARD, ctx); */

/*   VecScatterDestroy(ctx); */
/*   VecDestroy(nCoor); */

/*   PetscPrintf(PETSC_COMM_WORLD, "tecout coor\n "); */

/*   Vec	nUcat, zUcat; */

/*   DACreateNaturalVector(fda, &nUcat); */
/*   DAGlobalToNaturalBegin(fda, user->Ucat, INSERT_VALUES, nUcat); */
/*   DAGlobalToNaturalEnd(fda, user->Ucat, INSERT_VALUES, nUcat); */

/*   VecScatterCreateToZero(nUcat, &ctx, &zUcat); */

/*   VecScatterBegin(nUcat, zUcat, INSERT_VALUES, SCATTER_FORWARD, ctx); */
/*   VecScatterEnd(nUcat, zUcat, INSERT_VALUES, SCATTER_FORWARD, ctx); */

/*   VecScatterDestroy(ctx); */

/*   Vec	zUcat_o; */
/*   VecDuplicate(zUcat, &zUcat_o); */
/*   VecDestroy(nUcat); */


/*   PetscPrintf(PETSC_COMM_WORLD, "tecout zucat\n "); */


/*   Vec	nP, zP; */


/*   DACreateNaturalVector(da, &nP); */
/*   DAGlobalToNaturalBegin(da, user->P, INSERT_VALUES, nP); */
/*   DAGlobalToNaturalEnd(da, user->P, INSERT_VALUES, nP); */

/*   VecScatterCreateToZero(nP, &ctx, &zP); */

/*   VecScatterBegin(nP, zP, INSERT_VALUES, SCATTER_FORWARD, ctx); */
/*   VecScatterEnd(nP, zP, INSERT_VALUES, SCATTER_FORWARD, ctx); */

/*   VecDestroy(nP); */
/*   VecScatterDestroy(ctx); */

/*   PetscPrintf(PETSC_COMM_WORLD, "tecout nP\n "); */

/*   Vec	nNvert, zNvert; */

/*   DACreateNaturalVector(da, &nNvert); */
/*   DAGlobalToNaturalBegin(da, user->Nvert, INSERT_VALUES, nNvert); */
/*   DAGlobalToNaturalEnd(da, user->Nvert, INSERT_VALUES, nNvert); */

/*   VecScatterCreateToZero(nNvert, &ctx, &zNvert); */
  
/*   VecScatterBegin(nNvert, zNvert, INSERT_VALUES, SCATTER_FORWARD, ctx); */
/*   VecScatterEnd(nNvert, zNvert, INSERT_VALUES, SCATTER_FORWARD, ctx); */
/*   VecDestroy(nNvert); */

/*   VecScatterDestroy(ctx); */

/*   PetscPrintf(PETSC_COMM_WORLD, "tecout nvert\n "); */

/* /\*   DAVecGetArray(fda, zCoor, &coor); *\/ */

/*   //  PetscPrintf(PETSC_COMM_WORLD, "TecOutput\n"); */


/* /\*   Cmpnts ***ucont; *\/ */
/* /\*   DAVecGetArray(fda, user->Ucont, &ucont); *\/ */
/* /\*   PetscPrintf(PETSC_COMM_WORLD, "Ucont %e %e %e\n", ucont[100][90][2].x, ucont[100][90][2].y, ucont[100][90][2].z); *\/ */
/* /\*   PetscPrintf(PETSC_COMM_WORLD, "Ucont %e %e %e\n", ucont[100][90][1].x, ucont[100][89][2].y, ucont[99][90][2].z); *\/ */
/* /\*   DAVecRestoreArray(fda, user->Ucont, &ucont); *\/ */
  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) {
/*   VecGetArray3d(zUcat, mz, my, mx*3, 0, 0, 0, &ucat); */
/*   VecGetArray3d(zUcat_o, mz, my, mx*3, 0, 0, 0, &ucat_o); */
  DAVecGetArray(fda, user->Ucat, &ucat);
  DAVecGetArray(fda, user->Ucat_o, &ucat_o);

  // change!!! k=1->k=0
  for (k=0; k<mz-1; k++) {
    for (j=0; j<my-1; j++) {
      for (i=0; i<mx-1; i++) {
/* 	if (i==mx-2) { */
/* 	  ucat_o[k][j][i].x = (ucat_o[k][j][i-1].x * 4 - */
/* 			       ucat_o[k][j][i-2].x) / 3.; */
/* 	  ucat_o[k][j][i].y = (ucat_o[k][j][i-1].y * 4 - */
/* 			       ucat_o[k][j][i-2].y) / 3.; */
/* 	  ucat_o[k][j][i].z = (ucat_o[k][j][i-1].z * 4 - */
/* 			       ucat_o[k][j][i-2].z) / 3.; */
/* 	} */
/* 	else { */
	  ucat_o[k][j][i].x = 0.125 *
	    (ucat[k][j][i].x + ucat[k][j][i+1].x +
	     ucat[k][j+1][i].x + ucat[k][j+1][i+1].x +
	     ucat[k+1][j][i].x + ucat[k+1][j][i+1].x +
	     ucat[k+1][j+1][i].x + ucat[k+1][j+1][i+1].x);
	  ucat_o[k][j][i].y = 0.125 *
	    (ucat[k][j][i].y + ucat[k][j][i+1].y +
	     ucat[k][j+1][i].y + ucat[k][j+1][i+1].y +
	     ucat[k+1][j][i].y + ucat[k+1][j][i+1].y +
	     ucat[k+1][j+1][i].y + ucat[k+1][j+1][i+1].y);
	  ucat_o[k][j][i].z = 0.125 *
	    (ucat[k][j][i].z + ucat[k][j][i+1].z +
	     ucat[k][j+1][i].z + ucat[k][j+1][i+1].z +
	     ucat[k+1][j][i].z + ucat[k+1][j][i+1].z +
	     ucat[k+1][j+1][i].z + ucat[k+1][j+1][i+1].z);
/* 	} */
      }
    }
  }

/*   VecRestoreArray3d(zUcat, mz, my, mx*3, 0, 0, 0, &ucat); */
/*   VecRestoreArray3d(zUcat_o, mz, my, mx*3, 0, 0, 0, &ucat_o); */
  DAVecRestoreArray(fda, user->Ucat, &ucat);
  DAVecRestoreArray(fda, user->Ucat_o, &ucat_o);


  sprintf(filen, "Result%3.3d.dat", ti);
  f = fopen(filen, "a");

/*   VecGetArray3d(zCoor, mz, my, mx*3, 0, 0, 0, &coor); */
  DAVecGetArray(fda, Coor, &coor);

  for (k=1; k<mz; k++) {
    for (j=1; j<my; j++) {
      for (i=1; i<mx; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",
		     coor[k-1][j-1][i-1].x);
      }
    }
  }
  PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
  for (k=1; k<mz; k++) {
    for (j=1; j<my; j++) {
      for (i=1; i<mx; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",
		     coor[k-1][j-1][i-1].y);
      }
    }
  }

  PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
  for (k=1; k<mz; k++) {
    for (j=1; j<my; j++) {
      for (i=1; i<mx; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",
		     coor[k-1][j-1][i-1].z);
      }
    }
  }

/*   VecRestoreArray3d(zCoor, mz, my, mx*3, 0, 0, 0, &coor); */
  DAVecRestoreArray(fda, Coor, &coor);

  PetscFPrintf(PETSC_COMM_WORLD, f, "\n");


/*   VecGetArray3d(zUcat_o, mz, my, mx*3, 0, 0, 0, &ucat_o); */
  DAVecGetArray(fda, user->Ucat_o, &ucat_o);


  for (k=1; k<mz; k++) {
    for (j=1; j<my; j++) {
      for (i=1; i<mx; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",
		     ucat_o[k-1][j-1][i-1].x);
      }
    }
  }
  PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
  for (k=1; k<mz; k++) {
    for (j=1; j<my; j++) {
      for (i=1; i<mx; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",
		     ucat_o[k-1][j-1][i-1].y);
      }
    }
  }
  PetscFPrintf(PETSC_COMM_WORLD, f, "\n"); 
  for (k=1; k<mz; k++) {
    for (j=1; j<my; j++) {
      for (i=1; i<mx; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",
		     ucat_o[k-1][j-1][i-1].z);
      }
    }
  }

  PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
  DAVecRestoreArray(fda, user->Ucat_o, &ucat_o);

/*   VecRestoreArray3d(zUcat_o, mz, my, mx*3, 0, 0, 0, &ucat_o); */


/*   VecGetArray3d(zP, mz, my, mx, 0, 0, 0, &p); */

  DAVecGetArray(da, user->P, &p);

  for (k=1; k<mz-1; k++) {
    for (j=1; j<my-1; j++) {
      for (i=1; i<mx-1; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%f\n",
		     p[k][j][i]);
      }
    }
  }

/*   VecRestoreArray3d(zP, mz, my, mx, 0, 0, 0, &p); */
  DAVecRestoreArray(da, user->P, &p);

  PetscFPrintf(PETSC_COMM_WORLD, f, "\n");
/*   for (k=1; k<mz; k++) { */
/*     for (j=1; j<my; j++) { */
/*       for (i=1; i<mx; i++) { */
/* 	PetscFPrintf(PETSC_COMM_WORLD, f, "%d\n", */
/* 		     (int)(nvert[k-1][j-1][i-1] + 0.5)); */
/*       } */
/*     } */
/*   } */

/*   VecGetArray3d(zNvert, mz, my, mx, 0, 0, 0, &nvert); */
  DAVecGetArray(da, user->Nvert, &nvert);

  for (k=1; k<mz-1; k++) {
    for (j=1; j<my-1; j++) {
      for (i=1; i<mx-1; i++) {
	PetscFPrintf(PETSC_COMM_WORLD, f, "%d\n",
		     (int)(nvert[k][j][i]+0.5));
      }
    }
  }

  fclose(f);
/*   VecRestoreArray3d(zNvert, mz, my, mx, 0, 0, 0, &nvert); */
  DAVecRestoreArray(da, user->Nvert, &nvert);

  }
/*   DAVecRestoreArray(fda, Coor, &coor); */

/*   VecDestroy(zCoor); */
/*   VecDestroy(zUcat_o); */

/*   VecDestroy(zUcat); */
/*   VecDestroy(zP); */
/*   VecDestroy(zNvert); */
  return(0);
}

PetscErrorCode FormMetrics(UserCtx *user)
{
  DA		cda;
  Cmpnts	***csi, ***eta, ***zet;
  PetscScalar	***aj;
  Vec		coords;
  Cmpnts	***coor;

  DA		da = user->da, fda = user->fda;
  Vec		Csi = user->Csi, Eta = user->Eta, Zet = user->Zet;
  Vec		Aj = user->Aj;
  Vec		ICsi = user->ICsi, IEta = user->IEta, IZet = user->IZet;
  Vec		JCsi = user->JCsi, JEta = user->JEta, JZet = user->JZet;
  Vec		KCsi = user->KCsi, KEta = user->KEta, KZet = user->KZet;
  Vec		IAj = user->IAj, JAj = user->JAj, KAj = user->KAj;

  
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;
  Cmpnts	***gs;
  PetscReal	***iaj, ***jaj, ***kaj;

  Vec		Cent = user->Cent; //local working array for storing cell center geometry

  Vec		Centx, Centy, Centz, lCoor;
  Cmpnts	***cent, ***centx, ***centy, ***centz;

  PetscInt	xs, ys, zs, xe, ye, ze;
  DALocalInfo	info;

  PetscInt	mx, my, mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscScalar	dxdc, dydc, dzdc, dxde, dyde, dzde, dxdz, dydz, dzdz;
  PetscInt	i, j, k, ia, ja, ka, ib, jb, kb;
  PetscInt	gxs, gxe, gys, gye, gzs, gze;
  PetscErrorCode	ierr;

  PetscReal	xcp, ycp, zcp, xcm, ycm, zcm;
  DAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  DAGetCoordinateDA(da, &cda);
  DAVecGetArray(cda, Csi, &csi);
  DAVecGetArray(cda, Eta, &eta);
  DAVecGetArray(cda, Zet, &zet);
  ierr = DAVecGetArray(da, Aj,  &aj); CHKERRQ(ierr);

  DAGetGhostedCoordinates(da, &coords);
  DAVecGetArray(fda, coords, &coor);


  //  VecDuplicate(coords, &Cent);
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  /* Calculating transformation metrics in i direction */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=xs; i<lxe; i++){
	/* csi = de X dz */
	dxde = 0.5 * (coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
		      coor[k  ][j-1][i  ].x - coor[k-1][j-1][i  ].x);
	dyde = 0.5 * (coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
		      coor[k  ][j-1][i  ].y - coor[k-1][j-1][i  ].y);
	dzde = 0.5 * (coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
		      coor[k  ][j-1][i  ].z - coor[k-1][j-1][i  ].z);
				       		    	    	 
	dxdz = 0.5 * (coor[k  ][j-1][i  ].x + coor[k  ][j  ][i  ].x -
		      coor[k-1][j-1][i  ].x - coor[k-1][j  ][i  ].x);
	dydz = 0.5 * (coor[k  ][j-1][i  ].y + coor[k  ][j  ][i  ].y -
		      coor[k-1][j-1][i  ].y - coor[k-1][j  ][i  ].y);
	dzdz = 0.5 * (coor[k  ][j-1][i  ].z + coor[k  ][j  ][i  ].z -
		      coor[k-1][j-1][i  ].z - coor[k-1][j  ][i  ].z);
	  
	csi[k][j][i].x = dyde * dzdz - dzde * dydz;
	csi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
	csi[k][j][i].z = dxde * dydz - dyde * dxdz;

	
      }
    }
  }

  // Need more work -- lg65
  /* calculating j direction metrics */
  for (k=lzs; k<lze; k++){
    for (j=ys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

	/* eta = dz X de */
	dxdc = 0.5 * (coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
		      coor[k  ][j  ][i-1].x - coor[k-1][j  ][i-1].x);
	dydc = 0.5 * (coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
		      coor[k  ][j  ][i-1].y - coor[k-1][j  ][i-1].y);
	dzdc = 0.5 * (coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
		      coor[k  ][j  ][i-1].z - coor[k-1][j  ][i-1].z);
			    		         	 		   	 
	dxdz = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
		      coor[k-1][j  ][i  ].x - coor[k-1][j  ][i-1].x);
	dydz = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
		      coor[k-1][j  ][i  ].y - coor[k-1][j  ][i-1].y);
	dzdz = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
		      coor[k-1][j  ][i  ].z - coor[k-1][j  ][i-1].z);
	  
	eta[k][j][i].x = dydz * dzdc - dzdz * dydc;
	eta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
	eta[k][j][i].z = dxdz * dydc - dydz * dxdc;

      }
    }
  }

  /* calculating k direction metrics */
  for (k=zs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){
	dxdc = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x -
		      coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x);
	dydc = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y -
		      coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y);
	dzdc = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z -
		      coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z);
			    		    	     	     	 
	dxde = 0.5 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
		      coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x);
	dyde = 0.5 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
		      coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y);
	dzde = 0.5 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
		      coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z);
	  
	zet[k][j][i].x = dydc * dzde - dzdc * dyde;
	zet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
	zet[k][j][i].z = dxdc * dyde - dydc * dxde;

	/*	if ((i==1 || i==mx-2) && j==1 && (k==1 || k==0)) {
	  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", dxdc * dyde, dydc * dxde, dzdc);
	  PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", dxde, dyde, dzde);
	  PetscPrintf(PETSC_COMM_WORLD, "Met %e %e %e\n", zet[k][j][i].x, zet[k][j][i].y, zet[k][j][i].z);
	  }*/
	
      }
    }
  }

  /* calculating Jacobian of the transformation */
  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs; i<lxe; i++){

	if (i>0 && j>0 && k>0) {
	  dxdc = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
			 coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x -
			 coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x -
			 coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x);
	  dydc = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
			 coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y -
			 coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y -
			 coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y);
	  dzdc = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
			 coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z -
			 coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z -
			 coor[k-1][j  ][i-1].z - coor[k-1][j-1][i-1].z);

	  dxde = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x +
			 coor[k-1][j  ][i  ].x + coor[k-1][j  ][i-1].x - 
			 coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x -
			 coor[k-1][j-1][i  ].x - coor[k-1][j-1][i-1].x);
	  dyde = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y +
			 coor[k-1][j  ][i  ].y + coor[k-1][j  ][i-1].y - 
			 coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y -
			 coor[k-1][j-1][i  ].y - coor[k-1][j-1][i-1].y);
	  dzde = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z +
			 coor[k-1][j  ][i  ].z + coor[k-1][j  ][i-1].z - 
			 coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z -
			 coor[k-1][j-1][i  ].z - coor[k-1][j-1][i-1].z);

	  dxdz = 0.25 * (coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
			 coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x -
			 coor[k-1][j  ][i  ].x - coor[k-1][j-1][i  ].x -
			 coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x);
	  dydz = 0.25 * (coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
			 coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y -
			 coor[k-1][j  ][i  ].y - coor[k-1][j-1][i  ].y -
			 coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y);
	  dzdz = 0.25 * (coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
			 coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z -
			 coor[k-1][j  ][i  ].z - coor[k-1][j-1][i  ].z -
			 coor[k-1][j  ][i-1].z - coor[k-1][j-1][i-1].z);
	  
	  aj[k][j][i] = dxdc * (dyde * dzdz - dzde * dydz) -
	    dydc * (dxde * dzdz - dzde * dxdz) +
	    dzdc * (dxde * dydz - dyde * dxdz);
	  aj[k][j][i] = 1./aj[k][j][i];
	}
      }
    }
  }

  // mirror grid outside the boundary
  if (xs==0) {
    i = xs;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	eta[k][j][i].x = eta[k][j][i+1].x;
	eta[k][j][i].y = eta[k][j][i+1].y;
	eta[k][j][i].z = eta[k][j][i+1].z;
	  
	zet[k][j][i].x = zet[k][j][i+1].x;
	zet[k][j][i].y = zet[k][j][i+1].y;
	zet[k][j][i].z = zet[k][j][i+1].z;

	aj[k][j][i] = aj[k][j][i+1];
      }
    }
  }

  if (xe==mx) {
    i = xe-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	eta[k][j][i].x = eta[k][j][i-1].x;
	eta[k][j][i].y = eta[k][j][i-1].y;
	eta[k][j][i].z = eta[k][j][i-1].z;
	  
	zet[k][j][i].x = zet[k][j][i-1].x;
	zet[k][j][i].y = zet[k][j][i-1].y;
	zet[k][j][i].z = zet[k][j][i-1].z;

	aj[k][j][i] = aj[k][j][i-1];
      }
    }
  }

  if (ys==0) {
    j = ys;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	csi[k][j][i].x = csi[k][j+1][i].x;
	csi[k][j][i].y = csi[k][j+1][i].y;
	csi[k][j][i].z = csi[k][j+1][i].z;
	  
	zet[k][j][i].x = zet[k][j+1][i].x;
	zet[k][j][i].y = zet[k][j+1][i].y;
	zet[k][j][i].z = zet[k][j+1][i].z;

	aj[k][j][i] = aj[k][j+1][i];
      }
    }
  }

  if (ye==my) {
    j = ye-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	csi[k][j][i].x = csi[k][j-1][i].x;
	csi[k][j][i].y = csi[k][j-1][i].y;
	csi[k][j][i].z = csi[k][j-1][i].z;
	  
	zet[k][j][i].x = zet[k][j-1][i].x;
	zet[k][j][i].y = zet[k][j-1][i].y;
	zet[k][j][i].z = zet[k][j-1][i].z;

	aj[k][j][i] = aj[k][j-1][i];
      }
    }
  }

  if (zs==0) {
    k = zs;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	eta[k][j][i].x = eta[k+1][j][i].x;
	eta[k][j][i].y = eta[k+1][j][i].y;
	eta[k][j][i].z = eta[k+1][j][i].z;
	  		      			   
	csi[k][j][i].x = csi[k+1][j][i].x;
	csi[k][j][i].y = csi[k+1][j][i].y;
	csi[k][j][i].z = csi[k+1][j][i].z;

	aj[k][j][i] = aj[k+1][j][i];
      }
    }
  }

  if (ze==mz){
    k = ze-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	eta[k][j][i].x = eta[k-1][j][i].x;
	eta[k][j][i].y = eta[k-1][j][i].y;
	eta[k][j][i].z = eta[k-1][j][i].z;
	  		      			   
	csi[k][j][i].x = csi[k-1][j][i].x;
	csi[k][j][i].y = csi[k-1][j][i].y;
	csi[k][j][i].z = csi[k-1][j][i].z;

	aj[k][j][i] = aj[k-1][j][i];
      }
    }
  }



  //  PetscPrintf(PETSC_COMM_WORLD, "Local info: %d", info.mx);



  DAVecRestoreArray(cda, Csi, &csi);
  DAVecRestoreArray(cda, Eta, &eta);
  DAVecRestoreArray(cda, Zet, &zet);
  DAVecRestoreArray(da, Aj,  &aj);


  DAVecRestoreArray(cda, coords, &coor);


  VecAssemblyBegin(Csi);
  VecAssemblyEnd(Csi);
  VecAssemblyBegin(Eta);
  VecAssemblyEnd(Eta);
  VecAssemblyBegin(Zet);
  VecAssemblyEnd(Zet);
  VecAssemblyBegin(Aj);
  VecAssemblyEnd(Aj);



  DAGlobalToLocalBegin(fda, user->Csi, INSERT_VALUES, user->lCsi);
  DAGlobalToLocalEnd(fda, user->Csi, INSERT_VALUES, user->lCsi);

  DAGlobalToLocalBegin(fda, user->Eta, INSERT_VALUES, user->lEta);
  DAGlobalToLocalEnd(fda, user->Eta, INSERT_VALUES, user->lEta);

  DAGlobalToLocalBegin(fda, user->Zet, INSERT_VALUES, user->lZet);
  DAGlobalToLocalEnd(fda, user->Zet, INSERT_VALUES, user->lZet);


  DAGlobalToLocalBegin(da, user->Aj, INSERT_VALUES, user->lAj);
  DAGlobalToLocalEnd(da, user->Aj, INSERT_VALUES, user->lAj);


  PetscBarrier(PETSC_NULL);
  return 0;
}

PetscErrorCode Ucont_P_Binary_Input(UserCtx *user)
{
/*   PetscReal ***p; */
/*   Cmpnts	***ucont; */
/*   PetscInt	k; */
  PetscViewer	viewer;



  char filen2[90];

  PetscOptionsClearValue("-vecload_block_size");
  sprintf(filen2, "pfield%5.5d_%1.1d.dat", ti, user->this);

  PetscViewer	pviewer;
  Vec temp;
/*   VecDuplicate(user->Ucont, &temp); */
  PetscInt rank;
  DACreateNaturalVector(user->da, &temp);

//#ifdef PETSC_FILE_RDONLY  
//  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2,PETSC_FILE_RDONLY, &pviewer);
//#else
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &pviewer);
//#endif

  VecLoadIntoVector(pviewer, (user->P));
  PetscReal norm;
  VecNorm(user->P, NORM_INFINITY, &norm);
  PetscPrintf(PETSC_COMM_WORLD, "PIn %le\n", norm);
  PetscViewerDestroy(pviewer);

  VecDestroy(temp);

  sprintf(filen2, "nvfield%5.5d_%1.1d.dat", ti, user->this);
//#ifdef PETSC_FILE_RDONLY  
//  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2,  PETSC_FILE_RDONLY, &pviewer);
//#else
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, FILE_MODE_READ, &pviewer);
//#endif
  //  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen2, PETSC_FILE_RDONLY, &pviewer);
  VecLoadIntoVector(pviewer, (user->Nvert));
  PetscViewerDestroy(pviewer);
}

PetscErrorCode Ucont_P_Binary_Input1(UserCtx *user)
{
  PetscViewer viewer;
  char filen[90];

  
  sprintf(filen, "ufield%5.5d_%1.1d.dat", ti, user->this);
//#ifdef PETSC_FILE_RDONLY  
//  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, PETSC_FILE_RDONLY, &viewer);
//#else
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
//#endif

  //  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, PETSC_FILE_RDONLY, &viewer);
  PetscInt N;

  VecGetSize(user->Ucat, &N);
  PetscPrintf(PETSC_COMM_WORLD, "PPP %d\n", N);
  VecLoadIntoVector(viewer, (user->Ucat));
  //VecLoad(viewer, PETSC_NULL, &(user->P));
  PetscViewerDestroy(viewer);

  PetscBarrier(PETSC_NULL);


  sprintf(filen, "vfield%5.5d_%1.1d.dat", ti, user->this);
//#ifdef PETSC_FILE_RDONLY  
//  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, PETSC_FILE_RDONLY, &viewer);
//#else
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);
//#endif

  //  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, PETSC_FILE_RDONLY, &viewer);

  VecGetSize(user->Ucont, &N);
  PetscPrintf(PETSC_COMM_WORLD, "PPP %d\n", N);
  VecLoadIntoVector(viewer, (user->Ucont));
  //VecLoad(viewer, PETSC_NULL, &(user->P));
  PetscViewerDestroy(viewer);

  PetscBarrier(PETSC_NULL);


/*   VecGetSize(temp, &N); */
/*   PetscInt N1; */
/*   VecGetLocalSize(temp, &N1); */
/*   PetscPrintf(PETSC_COMM_WORLD, "PPP %d %d\n", N, N1); */

/*   VecGetSize(user->P, &N); */

/*   VecGetLocalSize(user->P, &N1); */
/*   PetscPrintf(PETSC_COMM_WORLD, "PPP %d %d\n", N, N1); */

/* /\*   VecCopy(temp, user->P); *\/ */

/* /\*   PetscInt ni[3]; *\/ */
/* /\*   PetscReal y[3]; *\/ */
/* /\*   ni[0] = 0; ni[1] = 1; ni[2] = 2; *\/ */
/* /\*   VecGetValues(temp, 3, ni, y); *\/ */
/* /\*   PetscPrintf(PETSC_COMM_SELF, "p %le %le %le\n", y[0], y[1], y[2]); *\/ */
/*   //VecLoadIntoVector(pviewer, (user->P)); */
/*   PetscViewerDestroy(pviewer); */

/*   DAVecGetArray(user->da, user->P, &p); */
/*   DAVecGetArray(user->fda, user->Ucont, &ucont); */
/*   for (k=1; k<70; k++) { */
/*     PetscPrintf(PETSC_COMM_WORLD, "p%le %le %le %le\n", ucont[k][11][10].x, ucont[k][11][10].y, ucont[k][11][10].z, p[k][11][10]); */
/*   } */
/*   DAVecRestoreArray(user->da, user->P, &p); */
/*   DAVecRestoreArray(user->fda, user->Ucont, &ucont); */

}

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **argv)
{
  DA	da, fda;
  Vec	qn, qnm;
  Vec	c;
  UserCtx	*user;

  PetscErrorCode ierr;


  PetscInitialize(&argc, &argv, (char *)0, help);

  PetscInt rank, bi;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if(!rank) {
    FILE *fd;
    fd = fopen("grid.dat", "r");

    fscanf(fd, "%i\n", &block_number);
    MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    fclose(fd);
  }
  else {
    MPI_Bcast(&block_number, 1, MPI_INT, 0, PETSC_COMM_WORLD);
  }
  PetscMalloc(block_number*sizeof(UserCtx), &user);

  ReadCoordinates(user);

  
  for (bi=0; bi<block_number; bi++) {
    ierr = DACreateGlobalVector(user[bi].fda, &user[bi].Csi); CHKERRQ(ierr);
    ierr = VecDuplicate(user[bi].Csi, &user[bi].Eta);
    ierr = VecDuplicate(user[bi].Csi, &user[bi].Zet);

    VecDuplicate(user[bi].Csi, &user[bi].Ucont);
    VecDuplicate(user[bi].Csi, &user[bi].Ucat);
    VecDuplicate(user[bi].Csi, &user[bi].Ucat_o);
    VecDuplicate(user[bi].Csi, &user[bi].Bcs.Ubcs);

    ierr = DACreateGlobalVector(user[bi].da, &user[bi].Aj); CHKERRQ(ierr);
    VecDuplicate(user[bi].Aj, &user[bi].P);

    VecDuplicate(user[bi].Aj, &user[bi].Nvert);

/*     DACreateLocalVector(user[bi].fda, &user[bi].lCsi); */

/*     VecDuplicate(user[bi].lCsi, &user[bi].lEta); */
/*     VecDuplicate(user[bi].lCsi, &user[bi].lZet); */

/*     VecDuplicate(user[bi].lCsi, &user[bi].lUcont); */
/*     VecDuplicate(user[bi].lCsi, &user[bi].lUcat); */

/*     DACreateLocalVector(user[bi].da, &user[bi].lAj); */
/*     VecDuplicate(user[bi].lAj, &user[bi].lP); */

/*     VecDuplicate(user[bi].lAj, &user[bi].lNvert); */

  /*   ierr = FormMetrics(&(user[bi])); */
  }

  PetscInt tis, tie, tistep;
  PetscErrorCode flag;
  PetscOptionsGetInt(PETSC_NULL, "-tis", &tis, &flag);
  if (!flag) {
    PetscPrintf(PETSC_COMM_WORLD, "Need the starting number!\n");    
  }

  PetscOptionsGetInt(PETSC_NULL, "-tie", &tie, PETSC_NULL);
  if (!flag) {
    tie = tis;
    //    PetscPrintf(PETSC_COMM_WORLD, "Need the ending number!\n");
  }

  PetscOptionsGetInt(PETSC_NULL, "-tio", &tistep, &flag);
  if (!flag) {
    tistep=5;
    //PetscPrintf(PETSC_COMM_WORLD, "Need the starting number!\n");    
  }
  
  if (flag) {
    FILE *f;
    char filen[80];
    for (ti=tis; ti<=tie; ti+=tistep) {

      sprintf(filen, "Result%3.3d.dat", ti);
      f = fopen(filen, "w");

  
      PetscFPrintf(PETSC_COMM_WORLD, f, "Variables = X, Y, Z, U, V, W, P, Nv\n");
      fclose(f);
      for (bi=0; bi<block_number; bi++) {
	Ucont_P_Binary_Input(&user[bi]);
	Ucont_P_Binary_Input1(&user[bi]);
      }
      
/*       DAGlobalToLocalBegin(user.fda, user.Ucont, INSERT_VALUES, user.lUcont); */
/*       DAGlobalToLocalEnd(user.fda, user.Ucont, INSERT_VALUES, user.lUcont); */

/*       Contra2Cart(&user); */
      for (bi=0; bi<block_number; bi++) {
	f = fopen(filen, "a");

	PetscFPrintf(PETSC_COMM_WORLD, f, "Zone I=%d, J=%d, K=%d F=BLOCK VARLOCATION=([7-8]=CELLCENTERED)\n", user[bi].info.mx-1, user[bi].info.my-1, user[bi].info.mz-1);
	fclose(f);
	TecOut(&(user[bi]));
      }
    }
  }
  PetscFinalize();
}



PetscErrorCode ReadCoordinates(UserCtx *user)
{
  Cmpnts ***coor;

  Vec Coor;
  PetscInt bi, i, j, k, rank, IM, JM, KM;
  PetscReal *gc;
  FILE *fd;
  PetscReal	d0 = 1.;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

/*   for (bi=0; bi<block_number; bi++) { */
/*     DASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0); */
/*   } */
  if (!rank) {

    fd = fopen("grid.dat", "r");
    fscanf(fd, "%i\n", &i);
  }

/*   for (bi=block_number-1; bi>=0; bi--) { */
  for (bi=0;bi<block_number; bi++) {
    if (!rank) {
      fscanf(fd, "%i %i %i\n", &(user[bi].IM), &(user[bi].JM), &(user[bi].KM));
      IM = user[bi].IM; JM = user[bi].JM; KM = user[bi].KM;

      MPI_Bcast(&(user[bi].IM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].JM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].KM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(&(user[bi].IM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].JM), 1, MPI_INT, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].KM), 1, MPI_INT, 0, PETSC_COMM_WORLD);

      IM = user[bi].IM; JM = user[bi].JM; KM = user[bi].KM;
    }
    DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX,
	       user[bi].IM+1, user[bi].JM+1, user[bi].KM+1, 1,1,
	       PETSC_DECIDE, 1, 2, PETSC_NULL, PETSC_NULL, PETSC_NULL,
	       &(user[bi].da));


    DASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    DAGetCoordinateDA(user[bi].da, &(user[bi].fda));
    DAGetLocalInfo(user[bi].da, &(user[bi].info));

    DALocalInfo	info = user[bi].info;
    PetscInt	xs = info.xs, xe = info.xs + info.xm;
    PetscInt  	ys = info.ys, ye = info.ys + info.ym;
    PetscInt	zs = info.zs, ze = info.zs + info.zm;
    PetscInt	mx = info.mx, my = info.my, mz = info.mz;
    
    PetscMalloc(3*(IM*JM*KM)*sizeof(PetscReal), &gc);
    DAGetGhostedCoordinates(user[bi].da, &Coor);
    DAVecGetArray(user[bi].fda, Coor, &coor);

    if (!rank) {
      for (k=0; k<KM; k++) {
	for (j=0; j<JM; j++) {
	  for (i=0; i<IM; i++) {
	    fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3);
	  }
	}
      }
      
      for (k=0; k<KM; k++) {
	for (j=0; j<JM; j++) {
	  for (i=0; i<IM; i++) {
	    fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 1);
	  }
	}
      }

      for (k=0; k<KM; k++) {
	for (j=0; j<JM; j++) {
	  for (i=0; i<IM; i++) {
	    fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 2);
	  }
	}
      }
    

      MPI_Bcast(gc, 3*(IM*JM*KM), MPIU_REAL, 0, PETSC_COMM_WORLD);
      
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if (k<KM && j<JM && i<IM) {
	      coor[k][j][i].x = *(gc + (k * (IM*JM) + j * IM + i) * 3  );
	      coor[k][j][i].y = *(gc + (k * (IM*JM) + j * IM + i) * 3+1);
	      coor[k][j][i].z = *(gc + (k * (IM*JM) + j * IM + i) * 3+2);
	    }
	  }
	}
      }
      
    }
    else {


      MPI_Bcast(gc, 3*(IM*JM*KM), MPIU_REAL, 0, PETSC_COMM_WORLD);

      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if (k<KM && j<JM && i<IM) {
	      coor[k][j][i].x = *(gc + (k * (IM*JM) + j * IM + i) * 3  );
	      coor[k][j][i].y = *(gc + (k * (IM*JM) + j * IM + i) * 3+1);
	      coor[k][j][i].z = *(gc + (k * (IM*JM) + j * IM + i) * 3+2);
	    }
	  }
	}
      }

    }
    PetscFree(gc);
    DAVecRestoreArray(user[bi].fda, Coor, &coor);

    Vec	gCoor;
    DAGetCoordinates(user[bi].da, &gCoor);
    DALocalToGlobal(user[bi].fda, Coor, INSERT_VALUES, gCoor);

    DAGlobalToLocalBegin(user[bi].fda, gCoor, INSERT_VALUES, Coor);
    DAGlobalToLocalEnd(user[bi].fda, gCoor, INSERT_VALUES, Coor);

  }
  if (!rank) {

    fclose(fd);
  }

  for (bi=0; bi<block_number; bi++) {
    user[bi].this = bi;
  }
  return(0);
}
