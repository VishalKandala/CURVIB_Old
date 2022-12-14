#include "variables.h"
#include <petscvec.h>

extern PetscInt   dof, outghost, ConstitutiveLawNonLinear, initrest, NumberOfBodies;
extern PetscReal  dt;

extern Cmpnts  PLUS(Cmpnts v1,Cmpnts v2);
extern Cmpnts  MINUS(Cmpnts v1,Cmpnts v2);
extern Cmpnts  CROSS(Cmpnts v1, Cmpnts v2);
extern PetscReal  DOT(Cmpnts v1, Cmpnts v2);
extern Cmpnts  UNIT(Cmpnts v1);
extern PetscReal  SIZE(Cmpnts v1);


PetscErrorCode  Dimension(IBMNodes *ibm, FE *fem, PetscInt ibi) {

  char  string[128];
  FILE  *fd;
  char  filen[80];

  fem->n_v = ibm->n_v;
  fem->n_elmt = ibm->n_elmt;
  
  sprintf(filen, "blist%2.2d" , ibi);
  fd = fopen(filen, "r");
  fscanf(fd, "%i", &(fem->n_edge));
  
  fgets(string, 128, fd);// skip line one
  fscanf(fd, "%i", &fem->n_ghosts);

  return(0);
}

//-----------------------------------------------------------------------------------------------------------------------------------
PetscErrorCode  Create(FE *fem, PetscInt ibi) {
 
  fem->ibi = ibi; 

  PetscMalloc((fem->n_v + fem->n_ghosts)*sizeof(PetscReal), &(fem->x_bp));
  PetscMalloc((fem->n_v + fem->n_ghosts)*sizeof(PetscReal), &(fem->y_bp));
  PetscMalloc((fem->n_v + fem->n_ghosts)*sizeof(PetscReal), &(fem->z_bp));
  
  PetscMalloc((fem->n_v + fem->n_ghosts)*sizeof(PetscReal), &(fem->x_bp0));
  PetscMalloc((fem->n_v + fem->n_ghosts)*sizeof(PetscReal), &(fem->y_bp0));
  PetscMalloc((fem->n_v + fem->n_ghosts)*sizeof(PetscReal), &(fem->z_bp0));
 
  PetscMalloc((fem->n_elmt + 2*fem->n_ghosts)*sizeof(PetscInt), &(fem->nv1));
  PetscMalloc((fem->n_elmt + 2*fem->n_ghosts)*sizeof(PetscInt), &(fem->nv2));
  PetscMalloc((fem->n_elmt + 2*fem->n_ghosts)*sizeof(PetscInt), &(fem->nv3));

  PetscMalloc(fem->n_elmt*sizeof(Cmpnts), &(fem->n_fib));

  PetscMalloc(fem->n_edge*sizeof(PetscInt), &(fem->n_bnodes));
  
  PetscMalloc(fem->n_elmt*sizeof(PetscInt), &(fem->nv4));
  PetscMalloc(fem->n_elmt*sizeof(PetscInt), &(fem->nv5));
  PetscMalloc(fem->n_elmt*sizeof(PetscInt), &(fem->nv6));
  
  PetscMalloc(dof*fem->n_elmt*sizeof(PetscReal), &(fem->kve0));
  PetscMalloc(dof*fem->n_elmt*sizeof(PetscReal), &(fem->kve));
 
  PetscMalloc(dof*fem->n_elmt*sizeof(PetscReal), &(fem->StressM)); 
  PetscMalloc((dof + 2)*fem->n_elmt*sizeof(PetscReal), &(fem->StrainM)); 
  PetscMalloc(dof*fem->n_elmt*sizeof(PetscReal), &(fem->StressB)); 
  PetscMalloc(dof*fem->n_elmt*sizeof(PetscReal), &(fem->StrainB));
  
  PetscMalloc(fem->n_elmt*sizeof(PetscReal), &(fem->Pf));
  PetscMalloc(fem->n_elmt*sizeof(PetscReal), &(fem->Pfn));

  PetscMalloc((fem->n_elmt + 2*fem->n_ghosts)*sizeof(PetscReal), &(fem->dA0));
  PetscMalloc((fem->n_elmt + 2*fem->n_ghosts)*sizeof(PetscReal), &(fem->dA)); 

  PetscMalloc((fem->n_elmt + 2*fem->n_ghosts)*sizeof(PetscReal), &(fem->nf_x));
  PetscMalloc((fem->n_elmt + 2*fem->n_ghosts)*sizeof(PetscReal), &(fem->Nf_x)); 
  PetscMalloc((fem->n_elmt + 2*fem->n_ghosts)*sizeof(PetscReal), &(fem->nf_y));
  PetscMalloc((fem->n_elmt + 2*fem->n_ghosts)*sizeof(PetscReal), &(fem->Nf_y)); 
  PetscMalloc((fem->n_elmt + 2*fem->n_ghosts)*sizeof(PetscReal), &(fem->nf_z));
  PetscMalloc((fem->n_elmt + 2*fem->n_ghosts)*sizeof(PetscReal), &(fem->Nf_z));

  PetscMalloc(fem->n_elmt*sizeof(PetscInt), &(fem->ire));
  PetscMalloc(fem->n_elmt*sizeof(PetscInt), &(fem->irv));
  PetscMalloc(fem->n_elmt*sizeof(PetscInt), &(fem->val));

  PetscMalloc(16*fem->n_elmt*sizeof(PetscInt), &(fem->patch));
  PetscMalloc(dof*fem->n_elmt*sizeof(PetscReal), &(fem->G));
  PetscMalloc(dof*fem->n_elmt*sizeof(PetscReal), &(fem->G1));
  PetscMalloc(dof*fem->n_elmt*sizeof(PetscReal), &(fem->G2));
 
  //contact
  PetscMalloc((fem->n_v + fem->n_ghosts)*sizeof(PetscInt), &(fem->contact));
  PetscMalloc(fem->n_elmt*sizeof(Cmpnts), &(fem->qvec));
  PetscMalloc(fem->n_elmt*sizeof(PetscReal), &(fem->radvec));

  VecCreateSeq(PETSC_COMM_SELF, dof*(fem->n_v + fem->n_ghosts), &(fem->Res));
  VecDuplicate(fem->Res, &(fem->x));
  VecDuplicate(fem->Res, &(fem->V));
  VecDuplicate(fem->Res, &(fem->xn));
  VecDuplicate(fem->Res, &(fem->xnm1));
  VecDuplicate(fem->Res, &(fem->xd));
  VecDuplicate(fem->Res, &(fem->xdsc));
  VecDuplicate(fem->Res, &(fem->dx));
  VecDuplicate(fem->Res, &(fem->dxn));
  VecDuplicate(fem->Res, &(fem->xdd));
  VecDuplicate(fem->Res, &(fem->xddsc));
  VecDuplicate(fem->Res, &(fem->y));
  VecDuplicate(fem->Res, &(fem->yn));
  VecDuplicate(fem->Res, &(fem->Fint));
  VecDuplicate(fem->Res, &(fem->Fext));
  VecDuplicate(fem->Res, &(fem->Fdyn));
  VecDuplicate(fem->Res, &(fem->disp));
  VecDuplicate(fem->Res, &(fem->FJ));
  VecDuplicate(fem->Res, &(fem->Fcnt));  //contact force
  VecDuplicate(fem->Res, &(fem->dxsc));
  VecDuplicate(fem->Res, &(fem->xnsc));
  VecDuplicate(fem->Res, &(fem->xninner));
  VecDuplicate(fem->Res, &(fem->yninner));
  VecDuplicate(fem->Res, &(fem->Mass));  
  VecDuplicate(fem->Res, &(fem->Dissip));

  //Initialize
  VecSet(fem->Res, 0.0);  VecSet(fem->x, 0.0);  VecSet(fem->xn, 0.0);  VecSet(fem->xnm1, 0.0);
  VecSet(fem->xd, 0.0);  VecSet(fem->xdd, 0.0);  VecSet(fem->y, 0.0);  VecSet(fem->yn, 0.0);
  VecSet(fem->Fint, 0.0);  VecSet(fem->Fext, 0.0);  VecSet(fem->Fdyn, 0.0);  VecSet(fem->disp, 0.0); 
  VecSet(fem->FJ, 0.0);  VecSet(fem->Fcnt, 0.0);
  VecSet(fem->Mass, 0.0);  VecSet(fem->Dissip, 0.0);

  return (0); 
}

//-----------------------------------------------------------------------------------------------------------------------------------
PetscErrorCode Input(IBMNodes *ibm, FE *fem, PetscInt ibi) {

  PetscInt  i, ii, nc, ec, n_elmt, n_v;
  n_elmt = fem->n_elmt;  n_v = fem->n_v;
  
  char  string[128];
  FILE  *fd;
  char  filen[80];

  //--------------------------------------------Reading nodes list  
  for (nc=0; nc<n_v; nc++) {
    fem->x_bp[nc] = ibm->x_bp[nc];  fem->y_bp[nc] = ibm->y_bp[nc];  fem->z_bp[nc]=ibm->z_bp[nc];
    fem->x_bp0[nc] = ibm->x_bp[nc];  fem->y_bp0[nc] = ibm->y_bp[nc];  fem->z_bp0[nc]=ibm->z_bp[nc];
  }    

  PetscPrintf(PETSC_COMM_WORLD, "number of nodes of list (body:%d) %d \n", ibi, fem->n_v);

  //------------------------------------------Reading elements list 
  PetscInt  *nv1, *nv2, *nv3;   
  PetscMalloc(n_elmt*sizeof(PetscInt), &nv1);
  PetscMalloc(n_elmt*sizeof(PetscInt), &nv2);
  PetscMalloc(n_elmt*sizeof(PetscInt), &nv3);
  
  for (ec=0; ec<n_elmt; ec++) {
    nv1[ec] = ibm->nv1[ec];  nv2[ec] = ibm->nv2[ec];  nv3[ec] = ibm->nv3[ec]; 
    fem->nv1[ec] = nv1[ec];  fem->nv2[ec] = nv2[ec];  fem->nv3[ec] = nv3[ec];
  }

  PetscPrintf(PETSC_COMM_WORLD, "number of element of list(body:%d) %d \n", ibi, fem->n_elmt);

  //--------------------------------------Reading Boundary nodes
  sprintf(filen,"blist%2.2d", ibi);
  fd = fopen(filen, "r");
  PetscInt  n_edge, *bnodes;
  fscanf(fd, "%i", &n_edge);
  
  PetscInt  *n_bnodes, sum_n_bnodes=0;
  
  PetscMalloc(n_edge*sizeof(PetscInt), &n_bnodes);
  i = -1;
  fgets(string, 128, fd);// skip line one
  fscanf(fd, "%i", &(fem->n_ghosts));
  fgets(string, 128, fd);// skip line two
  PetscPrintf(PETSC_COMM_WORLD, "Number of ghost nodes %d\n", fem->n_ghosts);

  while (i+1<n_edge) {
    i++;
    fscanf(fd, "%d \n", &(n_bnodes[i]));
    sum_n_bnodes += n_bnodes[i];
  }

  PetscMalloc(sum_n_bnodes*sizeof(PetscInt), &bnodes);
        
  i = -1;
  while (!feof(fd)) {
    i++;
    fscanf(fd, "%d", &bnodes[i]);
  }
  fclose(fd);
  //Transfer data to ctx
  fem->n_edge = n_edge;
  fem->sum_n_bnodes = sum_n_bnodes;
  for (i=0; i<n_edge; i++) {fem->n_bnodes[i] = n_bnodes[i];}

  PetscMalloc(fem->sum_n_bnodes*sizeof(PetscInt), &(fem->bnodes));

  for (i=0; i<sum_n_bnodes; i++) {fem->bnodes[i] = bnodes[i] - 1;}
    
  //------------------------------Form the patch nodes
  PetscInt  n1e, n2e, n3e, *nv4, *nv5, *nv6;

  PetscMalloc(n_elmt*sizeof(PetscInt), &nv4);
  PetscMalloc(n_elmt*sizeof(PetscInt), &nv5);
  PetscMalloc(n_elmt*sizeof(PetscInt), &nv6);
  for (i=0; i<n_elmt; i++) { // A milion means it does not have patch node (it is on boundary)
    nv4[i] = 1000000;  nv5[i] = 1000000;  nv6[i] = 1000000;
  }
  
  PetscInt   j=0, n1pe, n2pe, n3pe;
  PetscInt   mn, npe; //mn: mutual nodes counter , cn:column number
  PetscReal  cn;

  for (i=0; i<n_elmt; i++) {
    n1e = nv1[i];  n2e = nv2[i];  n3e = nv3[i];
    
    for (j=0; j<n_elmt; j++) {
      n1pe = nv1[j];  n2pe = nv2[j];  n3pe = nv3[j];
      
      mn = 0; cn = 0; npe = 0;
      if(n1e==n1pe || n1e==n2pe || n1e==n3pe){mn = mn+1;  cn = cn + 3.5;}
      if(n2e==n1pe || n2e==n2pe || n2e==n3pe){mn = mn+1;  cn = cn + 2.5;}
      if(n3e==n1pe || n3e==n2pe || n3e==n3pe){mn = mn+1;  cn = cn + 1.5;}
      
      if(mn==2){ //we catch the patch, now find the patch element number
  	if(n1pe!=n1e && n1pe!=n2e && n1pe!=n3e){
  	  npe = n1pe;
  	}else if(n2pe!=n1e && n2pe!=n2e && n2pe!=n3e){
  	  npe = n2pe;
  	}else{
  	  npe = n3pe;
  	}
  	// put it in right location
	if(cn==4.){
  	  nv4[i] = npe;
  	}else if(cn==5.){
  	  nv5[i] = npe;
  	}else{
          nv6[i] = npe;
  	}
      } //end if catch
    }// end neighbor elements check
  }// end patch find
  
  //Transfer data to FE
  for (ec=0; ec<n_elmt; ec++) {
    fem->nv4[ec] = nv4[ec];  fem->nv5[ec] = nv5[ec];  fem->nv6[ec] = nv6[ec];
  }

  PetscFree(nv1);  PetscFree(nv2);  PetscFree(nv3);
  PetscFree(n_bnodes);  PetscFree(bnodes); 
  PetscFree(nv4);  PetscFree(nv5);  PetscFree(nv6);  
  
  return(0);   
}

//-----------------------------------------------------------------------------------------------------------------------------------
PetscErrorCode Output(FE *fem, PetscInt ti, PetscInt ibi){

  PetscInt n_cells=3,i;  
  FILE     *f;
  char     filen[80];

  sprintf(filen, "fem_surface%2.2d_%5.5d.vtk", ibi,ti);
  f = fopen(filen, "w"); // open file
  
  PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 2.0\n");
  PetscFPrintf(PETSC_COMM_WORLD, f, "Surface Grid\n");
  PetscFPrintf(PETSC_COMM_WORLD, f, "ASCII\n");
  PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET UNSTRUCTURED_GRID\n");
  
  PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  %d float\n",(fem->n_v));
  for (i=0; i<fem->n_v; i++) {
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", (fem->x_bp[i]),(fem->y_bp[i]),(fem->z_bp[i]));
  }

  PetscFPrintf(PETSC_COMM_WORLD, f, "CELLS %d %d\n",fem->n_elmt, (n_cells+1)*(fem->n_elmt));
  for (i=0; i<fem->n_elmt; i++) {
    PetscFPrintf(PETSC_COMM_WORLD,f, "%d  %d %d %d\n",n_cells,(fem->nv1[i]),(fem->nv2[i]),(fem->nv3[i]));
  }
  
  PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_TYPES %d\n",fem->n_elmt);
  for (i=0; i<fem->n_elmt; i++) {
    PetscFPrintf(PETSC_COMM_WORLD,f, "%d\n",5);
  }
  
  PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n", fem->n_v);

  PetscFPrintf(PETSC_COMM_WORLD, f,  "SCALARS contact integer\n");
  PetscFPrintf(PETSC_COMM_WORLD, f,  "LOOKUP_TABLE default\n");
  for (i=0; i<fem->n_v; i++){
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d\n", fem->contact[i]);
  }

  //compute displacement
  PetscReal  *dd,*FF; 
  PetscInt   nv;
  VecGetArray(fem->disp, &dd);
  for (nv=0; nv<fem->n_v; nv++) {
    dd[nv*dof] = fem->x_bp[nv]-fem->x_bp0[nv];
    dd[nv*dof+1] = fem->y_bp[nv]-fem->y_bp0[nv];
    dd[nv*dof+2] = fem->z_bp[nv]-fem->z_bp0[nv];
  }
  PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS disp float\n");
  for (i=0; i<fem->n_v; i++){
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", dd[i*dof], dd[i*dof+1], dd[i*dof+2]);
  }
  VecRestoreArray(fem->disp, &dd);
  
  VecGetArray(fem->Fint, &FF);
  PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS Fint float\n");
  for (i=0; i<fem->n_v; i++){
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", FF[i*dof], FF[i*dof+1], FF[i*dof+2]);
  }
  VecRestoreArray(fem->Fint, &FF);
  
  VecGetArray(fem->Fext, &FF);
  PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS Fext float\n");
  for (i=0; i<fem->n_v; i++){
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", FF[i*dof], FF[i*dof+1], FF[i*dof+2]);
  }
  VecRestoreArray(fem->Fext, &FF);

  VecGetArray(fem->Fdyn, &FF);
  PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS Fdyn float\n");
  for (i=0; i<fem->n_v; i++){
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", FF[i*dof], FF[i*dof+1], FF[i*dof+2]);
  }
  VecRestoreArray(fem->Fdyn, &FF);

  VecGetArray(fem->Fcnt, &FF);
  PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS Fcnt float\n");
  for (i=0; i<fem->n_v; i++){
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", FF[i*dof], FF[i*dof+1], FF[i*dof+2]);
  }
  VecRestoreArray(fem->Fcnt, &FF);

  VecGetArray(fem->xd, &FF);
  PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS u float\n");
  for (i=0; i<fem->n_v; i++){
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", FF[i*dof], FF[i*dof+1], FF[i*dof+2]);
  }
  VecRestoreArray(fem->xd, &FF);

  PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_DATA %d\n",fem->n_elmt);
  PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS kve float\n");
  for (i=0; i<fem->n_elmt; i++) {
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n",fem->kve[i*dof], fem->kve[i*dof+1], fem->kve[i*dof+2]);
  }

  PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS nf float\n");
  for (i=0; i<fem->n_elmt; i++) {
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n",fem->nf_x[i], fem->nf_y[i], fem->nf_z[i]);
  }
  
  if(ConstitutiveLawNonLinear){
    PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS nfib float\n");
    for (i=0; i<fem->n_elmt; i++) {
      PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n",fem->n_fib[i].x, fem->n_fib[i].y, fem->n_fib[i].z);
    }
  }

  PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS StrainM float\n");
  for (i=0; i<fem->n_elmt; i++) {
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n",fem->StrainM[i*(dof+2)], fem->StrainM[i*(dof+2)+1], fem->StrainM[i*(dof+2)+2]);
  }

  PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS PStrainM float\n");
  for (i=0; i<fem->n_elmt; i++) {
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n",fem->StrainM[i*(dof+2)+3], fem->StrainM[i*(dof+2)+4], 0.);
  }
  
  PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS StrainB float\n");
  for (i=0; i<fem->n_elmt; i++) {
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n",fem->StrainB[i*dof], fem->StrainB[i*dof+1], fem->StrainB[i*dof+2]);
  }

  PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS StressM float\n");
  for (i=0; i<fem->n_elmt; i++) {
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n",fem->StressM[i*dof], fem->StressM[i*dof+1], fem->StressM[i*dof+2]);
  }

  PetscFPrintf(PETSC_COMM_WORLD, f,  "VECTORS StressB float\n");
  for (i=0; i<fem->n_elmt; i++) {
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n",fem->StressB[i*dof], fem->StressB[i*dof+1], fem->StressB[i*dof+2]);
  }

  PetscFPrintf(PETSC_COMM_WORLD, f,  "SCALARS dA float\n");
  PetscFPrintf(PETSC_COMM_WORLD, f,  "LOOKUP_TABLE default\n");
  for (i=0; i<fem->n_elmt; i++) {
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f \n", fem->dA[i]);
  }
   
  fclose(f);

  /* sprintf(filen, "TipDisp.dat"); */
  /* f = fopen(filen, "a"); */
  /* //PetscFPrintf(PETSC_COMM_WORLD, f, "%le  %le\n", ti*dt, fem->z_bp[56-1]/0.04); */
  /* PetscFPrintf(PETSC_COMM_WORLD, f, "%le  %le\n", ti*dt, fem->z_bp[52-1]); */
  /* fclose(f); */
   
  return(0);
}

//-----------------------------------------------------------------------------------------------------------------------------------
PetscErrorCode OutputGhost(FE *fem, PetscInt ti, PetscInt ibi){
 
  PetscInt   n_cells=3, i;
  PetscReal  x, y, z;
  PetscInt   ec, be, n1e, n2e, n3e, n_ghosts=0;

  n_ghosts = fem->n_ghosts;

  FILE  *f;
  char  filen[80];

  sprintf(filen, "surfaceghost%2.2d_%5.5d.vtk", ibi, ti);
  f = fopen(filen, "w"); // open file

  PetscFPrintf(PETSC_COMM_WORLD, f, "# vtk DataFile Version 2.0\n");
  PetscFPrintf(PETSC_COMM_WORLD, f, "Surface Grid\n");
  PetscFPrintf(PETSC_COMM_WORLD, f, "ASCII\n");
  PetscFPrintf(PETSC_COMM_WORLD, f, "DATASET UNSTRUCTURED_GRID\n");
   
  PetscFPrintf(PETSC_COMM_WORLD, f, "POINTS  %d float\n",(fem->n_v+n_ghosts));
  for (i=0; i<fem->n_v+n_ghosts; i++) {
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", (fem->x_bp[i]), (fem->y_bp[i]), (fem->z_bp[i]));
  }

  //add ghost nodes location
  /* for (i=fem->n_v; i<(fem->n_v+n_ghosts); i++) { */
  /*   ec = i - fem->n_v; */
  /*   be = fem->belmts[ec]; */
  /*   if(fem->edgefrontnodesI[ec]==1){ */
  /*     x = fem->p4x[be];  y = fem->p4y[be];  z = fem->p4z[be]; */
  /*   }else if(fem->edgefrontnodesI[ec]==2){ */
  /*     x = fem->p5x[be];  y = fem->p5y[be];  z = fem->p5z[be]; */
  /*   }else if(fem->edgefrontnodesI[ec]==3){ */
  /*     x = fem->p6x[be];  y = fem->p6y[be];  z = fem->p6z[be]; */
  /*   } */
    
  /*   PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", x, y, z); */
  /* } */

  PetscFPrintf(PETSC_COMM_WORLD, f, "CELLS %d %d\n", (fem->n_elmt+2*n_ghosts), (n_cells+1)*(fem->n_elmt+2*n_ghosts));
  for (i=0; i<fem->n_elmt+2*n_ghosts; i++) {
    PetscFPrintf(PETSC_COMM_WORLD,f, "%d  %d %d %d\n", n_cells, (fem->nv1[i]), (fem->nv2[i]), (fem->nv3[i]));
  }

  /* //add ghost elements  */
  /* for (i=fem->n_elmt; i<(fem->n_elmt+n_ghosts); i++) { */
  /*   ec = i-fem->n_elmt; */
  /*   be = fem->belmts[ec]; */
  /*   if(fem->edgefrontnodesI[ec]==1){ */
  /*     n1e = ec+fem->n_v;  n2e = fem->nv2[be];  n3e = fem->nv3[be]; */
  /*   }else if(fem->edgefrontnodesI[ec]==2){ */
  /*     n1e = fem->nv1[be];  n2e = ec+fem->n_v;  n3e = fem->nv3[be]; */
  /*   }else if(fem->edgefrontnodesI[ec]==3){ */
  /*     n1e = fem->nv1[be];  n2e = fem->nv2[be];  n3e = ec+fem->n_v; */
  /*   } */
  /*   PetscFPrintf(PETSC_COMM_WORLD,f, "%d  %d %d %d\n", n_cells, n1e, n2e, n3e); */
  /* } */

  PetscFPrintf(PETSC_COMM_WORLD, f, "CELL_TYPES %d\n",(fem->n_elmt+2*n_ghosts));
  for (i=0; i<(fem->n_elmt+2*n_ghosts); i++) {
    PetscFPrintf(PETSC_COMM_WORLD, f, "%d\n", 5);
  }

  PetscReal  *dd,*FF;
  PetscInt   nv; 

  VecGetArray(fem->disp, &dd);
  for (nv=0; nv<fem->n_v+n_ghosts; nv++) {
    dd[nv*dof] = fem->x_bp[nv]-fem->x_bp0[nv];
    dd[nv*dof+1] = fem->y_bp[nv]-fem->y_bp0[nv];
    dd[nv*dof+2] = fem->z_bp[nv]-fem->z_bp0[nv];
  }
  
  PetscFPrintf(PETSC_COMM_WORLD, f, "POINT_DATA %d\n", fem->n_v+n_ghosts);
  
  PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS disp float\n");
  for (i=0; i<fem->n_v+n_ghosts; i++){
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", dd[i*dof], dd[i*dof+1], dd[i*dof+2]);
  }
  VecRestoreArray(fem->disp, &dd);
  
  VecGetArray(fem->Fint, &FF);
  PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS Fint float\n");
  for (i=0; i<fem->n_v+n_ghosts; i++){
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", FF[i*dof], FF[i*dof+1], FF[i*dof+2]);
  }
  VecRestoreArray(fem->Fint, &FF);

  VecGetArray(fem->Fext, &FF);
  PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS Fext float\n");
  for (i=0; i<fem->n_v+n_ghosts; i++){
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", FF[i*dof], FF[i*dof+1], FF[i*dof+2]);
  }
  VecRestoreArray(fem->Fext, &FF);

  VecGetArray(fem->Fdyn, &FF);
  PetscFPrintf(PETSC_COMM_WORLD, f, "VECTORS Fdyn float\n");
  for (i=0; i<fem->n_v+n_ghosts; i++){
    PetscFPrintf(PETSC_COMM_WORLD, f, "%f %f %f\n", FF[i*dof], FF[i*dof+1], FF[i*dof+2]);
  }
  VecRestoreArray(fem->Fdyn, &FF);

  fclose(f);

  return(0);
}

//-----------------------------------------------------------------------------------------------------------------------------------
PetscErrorCode  LocationOut(FE *fem, PetscInt ti, PetscInt ibi) {
  
  PetscViewer   viewer;
  char          filen[80];
  PetscInt     fd;

  sprintf(filen, "x%1.1d_%5.5d.dat", ibi, ti);
  PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_WRITE, &viewer);
  VecView(fem->x, viewer);
  PetscViewerDestroy(&viewer);

  sprintf(filen, "xn%1.1d_%5.5d.dat", ibi, ti);
  PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_WRITE, &viewer);
  VecView(fem->xn, viewer);
  PetscViewerDestroy(&viewer);

  sprintf(filen, "xnm1%1.1d_%5.5d.dat", ibi, ti);
  PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_WRITE, &viewer);
  VecView(fem->xnm1, viewer);
  PetscViewerDestroy(&viewer);

  sprintf(filen, "Fcnt1%1.1d_%5.5d.dat", ibi, ti);
  PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_WRITE, &viewer);
  VecView(fem->Fcnt, viewer);
  PetscViewerDestroy(&viewer);

  sprintf(filen, "cnt%1.1d_%5.5d.dat", ibi, ti);
  PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_WRITE, &viewer);
  PetscViewerBinaryGetDescriptor(viewer,&fd);
  PetscBinaryWrite(fd, fem->contact, fem->n_v, PETSC_INT, PETSC_FALSE);
  PetscViewerDestroy(&viewer); 

  sprintf(filen, "xd1%1.1d_%5.5d.dat", ibi, ti);
  PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_WRITE, &viewer);
  VecView(fem->xd, viewer);
  PetscViewerDestroy(&viewer);

  sprintf(filen, "xdd1%1.1d_%5.5d.dat", ibi, ti);
  PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_WRITE, &viewer);
  VecView(fem->xdd, viewer);
  PetscViewerDestroy(&viewer);

  return(0);
}

//-----------------------------------------------------------------------------------------------------------------------------------
PetscErrorCode LocationIn(FE *fem, PetscInt ti, PetscInt ibi) {

  PetscViewer  viewer;
  char         filen[80];
  PetscInt     fd;
  PetscInt rank=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) {
    sprintf(filen, "x%1.1d_%5.5d.dat", ibi, ti);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_READ, &viewer);
    VecLoad(fem->x, viewer);
    
    sprintf(filen, "xn%1.1d_%5.5d.dat", ibi, ti);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_READ, &viewer);
    VecLoad(fem->xn,viewer);
    
    sprintf(filen, "xnm1%1.1d_%5.5d.dat", ibi, ti);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_READ, &viewer);
    VecLoad(fem->xnm1, viewer);
    
    sprintf(filen, "xd1%1.1d_%5.5d.dat", ibi, ti);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_READ, &viewer);
    VecLoad(fem->xd, viewer);
    
    sprintf(filen, "xdd1%1.1d_%5.5d.dat", ibi, ti);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_READ, &viewer);
    VecLoad(fem->xdd, viewer);
    
    sprintf(filen, "Fcnt1%1.1d_%5.5d.dat", ibi, ti);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_READ, &viewer);
    VecLoad(fem->Fcnt, viewer);
    
    sprintf(filen, "cnt%1.1d_%5.5d.dat", ibi, ti);
    PetscViewerBinaryOpen(PETSC_COMM_SELF, filen, FILE_MODE_READ, &viewer);
    PetscViewerBinaryGetDescriptor(viewer, &fd);
    PetscBinaryRead(fd, fem->contact, fem->n_v, PETSC_INT);  
    
    PetscViewerDestroy(&viewer);
  }
  return(0);
}

//-----------------------------------------------------------------------------------------------------------------------------------
PetscErrorCode  AreaNormal(FE *fem) {

  Cmpnts  x1, x2, x3, dx21, dx31, n, cross;
  PetscInt       ec, n1e, n2e, n3e;
  
  for (ec=0; ec<fem->n_elmt + 2*fem->n_ghosts; ec++) {
    n1e = fem->nv1[ec];  n2e = fem->nv2[ec];  n3e = fem->nv3[ec];
    
    //current location
    x1.x = fem->x_bp[n1e];  x1.y = fem->y_bp[n1e];  x1.z = fem->z_bp[n1e];
    x2.x = fem->x_bp[n2e];  x2.y = fem->y_bp[n2e];  x2.z = fem->z_bp[n2e];
    x3.x = fem->x_bp[n3e];  x3.y = fem->y_bp[n3e];  x3.z = fem->z_bp[n3e];
    
    dx21 = MINUS(x2, x1);  dx31 = MINUS(x3, x1);
    cross = CROSS(dx21, dx31);
    
    n = UNIT(cross);
    fem->dA[ec] = 0.5*SIZE(cross); 
    
    fem->nf_x[ec] = n.x;  fem->nf_y[ec] = n.y;  fem->nf_z[ec] = n.z;
  }
 
  return(0);
}

//-----------------------------------------------------------------------------------------------------------------------------------
PetscErrorCode RotateBHVleaflet(IBMNodes *ibm, PetscReal angle) {
  
  PetscInt   nv;
  PetscReal  x,y;
  for (nv=0; nv<ibm->n_v; nv++) {
    ibm->y_bp[nv] = ibm->y_bp_o[nv]*cos(angle) - ibm->x_bp_o[nv]*sin(angle);
    ibm->x_bp[nv] = ibm->y_bp_o[nv]*sin(angle) + ibm->x_bp_o[nv]*cos(angle);
  } 
  
  for (nv=0; nv<ibm->n_v; nv++) {
    ibm->x_bp_o[nv] = ibm->x_bp[nv];
    ibm->y_bp_o[nv] = ibm->y_bp[nv];
  }
  
  return(0);
}

