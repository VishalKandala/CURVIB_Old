#include "variables.h"

extern PetscInt  dof, nbody, ti;

extern Cmpnts  MINUS(Cmpnts v1,Cmpnts v2);
extern Cmpnts  CROSS(Cmpnts v1,Cmpnts v2);
extern PetscReal  SIZE(Cmpnts v1);
extern Cmpnts  UNIT(Cmpnts v1);
extern void initlist(List *ilist);
extern void insertnode(List *ilist, PetscInt Node);


PetscErrorCode  Fcontact(FE *fem) {

  Fcontactij(&fem[0], &fem[1], 0, 1);  Fcontactij(&fem[1], &fem[2], 1, 2);  Fcontactij(&fem[2], &fem[0], 2, 0);
  Fcontactij(&fem[1], &fem[0], 1, 0);  Fcontactij(&fem[2], &fem[1], 2, 1);  Fcontactij(&fem[0], &fem[2], 0, 2);

  return(0);
}

//------------------------------------------------------------------------------------------------------------
PetscErrorCode Fcontactij(FE *fem1, FE *fem2, PetscInt ii, PetscInt jj) {
 
  PetscInt   i, j, k, nv, ncx=15, ncy=15, ncz=15;
  PetscReal  xbp_min, ybp_min, zbp_min, xbp_max, ybp_max, zbp_max, dcx, dcy, dcz;
  PetscReal  *x_bp = fem2->x_bp, *y_bp = fem2->y_bp, *z_bp = fem2->z_bp;
  PetscInt   n_v = fem2->n_v, ln_v;
  Cmpnts     pnv;
  List       *cell_trg;
  PetscReal  xv_min, yv_min, zv_min, xv_max, yv_max, zv_max;
  PetscInt   iv_min, iv_max, jv_min, jv_max, kv_min, kv_max, n1e, n2e, n3e, ic, jc, kc;

  // Bounding Box for fem2
  xbp_min = 1.e23;  xbp_max = -1.e23;
  ybp_min = 1.e23;  ybp_max = -1.e23;
  zbp_min = 1.e23;  zbp_max = -1.e23;
  
  for(i=0; i<n_v; i++) {
    xbp_min = PetscMin(xbp_min, x_bp[i]);
    xbp_max = PetscMax(xbp_max, x_bp[i]);
    
    ybp_min = PetscMin(ybp_min, y_bp[i]);
    ybp_max = PetscMax(ybp_max, y_bp[i]);
    
    zbp_min = PetscMin(zbp_min, z_bp[i]);
    zbp_max = PetscMax(zbp_max, z_bp[i]);
  }
 
  xbp_min -= 0.05;  xbp_max += 0.05;
  ybp_min -= 0.05;  ybp_max += 0.05;
  zbp_min -= 0.05;  zbp_max += 0.05;
 
  //Control cell for fem2
  dcx = (xbp_max - xbp_min)/(ncx - 1.);
  dcy = (ybp_max - ybp_min)/(ncy - 1.);
  dcz = (zbp_max - zbp_min)/(ncz - 1.);

  PetscMalloc(ncz*ncy*ncx*sizeof(List), &cell_trg);
 
  for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {
  	initlist(&cell_trg[k*ncx*ncy + j*ncx + i]);
      }
    }
  }
  
  for (ln_v=0; ln_v<fem2->n_elmt; ln_v++) {

    n1e = fem2->nv1[ln_v];  n2e = fem2->nv2[ln_v];  n3e = fem2->nv3[ln_v];

    xv_min = PetscMin(PetscMin(x_bp[n1e], x_bp[n2e]), x_bp[n3e]);
    xv_max = PetscMax(PetscMax(x_bp[n1e], x_bp[n2e]), x_bp[n3e]);

    yv_min = PetscMin(PetscMin(y_bp[n1e], y_bp[n2e]), y_bp[n3e]);
    yv_max = PetscMax(PetscMax(y_bp[n1e], y_bp[n2e]), y_bp[n3e]);

    zv_min = PetscMin(PetscMin(z_bp[n1e], z_bp[n2e]), z_bp[n3e]);
    zv_max = PetscMax(PetscMax(z_bp[n1e], z_bp[n2e]), z_bp[n3e]);
    
    iv_min = floor((xv_min - xbp_min)/dcx);
    iv_max = floor((xv_max - xbp_min)/dcx) + 1;

    jv_min = floor((yv_min - ybp_min)/dcy);
    jv_max = floor((yv_max - ybp_min)/dcy) + 1;

    kv_min = floor((zv_min - zbp_min)/dcz);
    kv_max = floor((zv_max - zbp_min)/dcz) + 1;

    iv_min = (iv_min<0) ? 0:iv_min;
    iv_max = (iv_max>ncx) ? ncx:iv_max;

    jv_min = (jv_min<0) ? 0:jv_min;
    jv_max = (jv_max>ncy) ? ncy:jv_max;

    kv_min = (kv_min<0) ? 0:kv_min;
    kv_max = (kv_max>ncz) ? ncz:kv_max;
   
    // Insert FEM node information into a list
    for (k=kv_min; k<kv_max; k++) {
      for (j=jv_min; j<jv_max; j++) {
  	for (i=iv_min; i<iv_max; i++) {
  	  insertnode(&(cell_trg[k *ncx*ncy + j*ncx +i]), ln_v);
  	}
      }
    }
  }
  
  // search if fem1 intersects fem2 and find contact forces for fem1
  PetscReal  *xx1, dist, dmin=0.0007/0.022, dir[3], c=0.1;
  PetscInt   elmt2;
  
  VecGetArray(fem1->x, &xx1);  
  
  for (nv=0; nv <fem1->n_v; nv++) {
    pnv.x = fem1->x_bp[nv];
    pnv.y = fem1->y_bp[nv];
    pnv.z = fem1->z_bp[nv];

    if (pnv.x > xbp_min && pnv.x < xbp_max &&
    	pnv.y > ybp_min && pnv.y < ybp_max &&
    	pnv.z > zbp_min && pnv.z < zbp_max) { // if in bounding box
      
      ic = floor((pnv.x - xbp_min)/dcx);
      jc = floor((pnv.y - ybp_min)/dcy);
      kc = floor((pnv.z - zbp_min)/dcz);

      // find the closest triangle in control cells
      fem_nearestcell(pnv, fem2, &dist, &elmt2, dir, ic, jc, kc, ncx, ncy, ncz, cell_trg);
      
      if (elmt2>0) {      
 
	if (PetscAbsReal(dist)<dmin) {
	  
	  fem1->contact[nv] = 1;
	  
	  if (dist>0.) {
	    
	    xx1[dof*nv] = xx1[dof*nv] + c*(dmin - dist)*dir[0];
	    xx1[dof*nv+1] = xx1[dof*nv+1] + c*(dmin - dist)*dir[1];
	    xx1[dof*nv+2] = xx1[dof*nv+2] + c*(dmin - dist)*dir[2];
	    fem1->x_bp[nv] = xx1[dof*nv];
	    fem1->y_bp[nv] = xx1[dof*nv+1];
	    fem1->z_bp[nv] = xx1[dof*nv+2];	   
	    
	  } else {
	    
	    xx1[dof*nv] = xx1[dof*nv] - c*(dmin - dist)*dir[0];
	    xx1[dof*nv+1] = xx1[dof*nv+1] - c*(dmin - dist)*dir[1];
	    xx1[dof*nv+2] = xx1[dof*nv+2] - c*(dmin - dist)*dir[2];
	    fem1->x_bp[nv] = xx1[dof*nv];
	    fem1->y_bp[nv] = xx1[dof*nv+1];
	    fem1->z_bp[nv] = xx1[dof*nv+2];	   
	  }	  
	} // if dist<dmin
      }
    } // if in bounding box  
  } //nv

  VecRestoreArray(fem1->x, &xx1);
  EdgeFix(1, fem1); 
  
  return(0);
}

//------------------------------------------------------------------------------------------------------------
PetscErrorCode fem_nearestcell(Cmpnts p, FE *fem, PetscReal *dmin, PetscInt *cell_min, PetscReal dir[3],
			      PetscInt ic, PetscInt jc, PetscInt kc, PetscInt ncx, PetscInt ncy, PetscInt ncz, List *cell_trg) {

  PetscInt       *nv1=fem->nv1, *nv2=fem->nv2, *nv3=fem->nv3;
  PetscReal      *nf_x=fem->nf_x, *nf_y=fem->nf_y, *nf_z=fem->nf_z;
  PetscReal      *x_bp=fem->x_bp, *y_bp=fem->y_bp, *z_bp=fem->z_bp;
  PetscInt       n_elmt=fem->n_elmt;  
  PetscInt       ln_v;
  Cmpnts  p1, p2, p3, pc;
  PetscReal      tf;
  PetscInt       n1e, n2e, n3e;
  PetscReal      nfx, nfy, nfz;
  Cmpnts  pj; //projection point
  Cmpnts  pmin, po, direct;
  PetscReal      d, d_center;
  node           *current;
  PetscInt       i, j, k, im, jm, km;

  *dmin = 1.e10;
  *cell_min = -100;
  km = ncz;
  jm = ncy;
  im = ncx;
  //additionally search one before and one after the original control cell
  if (kc<ncz-2)  km = kc + 2;
  if (jc<ncy-2)  jm = jc + 2;
  if (ic<ncx-2)  im = ic + 2;

  if (kc>0) kc--;
  if (jc>0) jc--;
  if (ic>0) ic--;

  for (k=kc; k<km; k++) {
    for (j=jc; j<jm; j++) {
      for (i=ic; i<im; i++) {
  	current = cell_trg[k*ncx*ncy+j*ncx+i].head;

  	while (current) {
  	  ln_v = current->Node;
	  fem_BoundingSphere(fem, ln_v);

	  d_center = SIZE(MINUS(p, fem->qvec[ln_v]));
	  
	  if (PetscAbsReal(d_center - fem->radvec[ln_v])<PetscAbsReal(*dmin)) {
	    n1e = nv1[ln_v];  n2e = nv2[ln_v];  n3e = nv3[ln_v];
	    nfx = nf_x[ln_v];  nfy = nf_y[ln_v];  nfz = nf_z[ln_v];
	    
	    p1.x = x_bp[n1e];  p1.y = y_bp[n1e];  p1.z = z_bp[n1e];
	    p2.x = x_bp[n2e];  p2.y = y_bp[n2e];  p2.z = z_bp[n2e];
	    p3.x = x_bp[n3e];  p3.y = y_bp[n3e];  p3.z = z_bp[n3e];
	    
	    tf = ((p.x - x_bp[n1e])*nfx +
		  (p.y - y_bp[n1e])*nfy +
		  (p.z - z_bp[n1e])*nfz);
	    
	    pj.x = p.x - tf*nfx;
	    pj.y = p.y - tf*nfy;
	    pj.z = p.z - tf*nfz;
	    
	    if (fem_ISPointInTriangle(pj, p1, p2, p3, nfx, nfy, nfz) == 1) { /* The projected point is inside the  triangle */
	      
	      if (PetscAbsReal(tf)<PetscAbsReal(*dmin)) {
		*dmin = tf;	 
		pmin.x = pj.x;
		pmin.y = pj.y;
		pmin.z = pj.z;
		*cell_min = ln_v;	  	  
	      }
	      
	    } else {
	      Sign_Dis_P_Line(p, p1, p2, nfx, nfy, nfz, &po, &d);
	      if (PetscAbsReal(d)<PetscAbsReal(*dmin)) {
		*dmin = d;
		pmin.x = po.x;
		pmin.y = po.y;
		pmin.z = po.z;
		*cell_min = ln_v;
	      }
	      Sign_Dis_P_Line(p, p2, p3, nfx, nfy, nfz, &po, &d);
	      if (PetscAbsReal(d)<PetscAbsReal(*dmin)) {
		*dmin = d;
		pmin.x = po.x;
		pmin.y = po.y;
		pmin.z = po.z;
		*cell_min = ln_v;
	      }
	      Sign_Dis_P_Line(p, p3, p1, nfx, nfy, nfz, &po, &d);
	      if (PetscAbsReal(d)<PetscAbsReal(*dmin)) {
		*dmin = d;
		pmin.x = po.x;
		pmin.y = po.y;
		pmin.z = po.z;
		*cell_min = ln_v;
	      }
	    } // is in triangle
	  } // in the bounding sphere
	  current = current->next;
	} // while current
      }
    }
  }

  direct = UNIT(MINUS(p, pmin));  
  dir[0] = direct.x;  dir[1] = direct.y;  dir[2] = direct.z;

  return(0);
}

//------------------------------------------------------------------------------------------------------------
PetscErrorCode fem_BoundingSphere(FE *fem, PetscInt ln_v) {

  PetscInt       *nv1 = fem->nv1, *nv2 = fem->nv2, *nv3 = fem->nv3;
  PetscReal      *x_bp = fem->x_bp, *y_bp = fem->y_bp, *z_bp = fem->z_bp;
  PetscInt       n_elmt = fem->n_elmt;
  //PetscInt       ln_v;
  Cmpnts  p1, p2, p3;
  PetscInt       n1e, n2e, n3e;
  Cmpnts  pa, pb, pc, pu, pv, pf, pd, pt;
  PetscReal      l12, l23, l31;
  PetscReal      gama, lamda;

  n1e = nv1[ln_v];  n2e = nv2[ln_v];  n3e = nv3[ln_v];
  
  p1.x = x_bp[n1e];  p1.y = y_bp[n1e];  p1.z = z_bp[n1e];
  p2.x = x_bp[n2e];  p2.y = y_bp[n2e];  p2.z = z_bp[n2e];
  p3.x = x_bp[n3e];  p3.y = y_bp[n3e];  p3.z = z_bp[n3e];
  
  l12 = SIZE(MINUS(p1, p2));  l23 = SIZE(MINUS(p2, p3));  l31 = SIZE(MINUS(p3, p1));
  
  /* Find the longest edge and assign the corresponding two vertices
     to pa and pb */
  if (l12>l23) {
    if (l12>l31) {
      pa = p1;  pb = p2;  pc = p3;
    }
    else {
      pa = p3;  pb = p1;  pc = p2;
    }
  }
  else {
    if (l31<l23) {
      pa = p2;  pb = p3;  pc = p1;
    }
    else {
      pa = p3;  pb = p1;  pc = p2;
    }
  }
  
  pf.x = 0.5*(pa.x + pb.x);
  pf.y = 0.5*(pa.y + pb.y);
  pf.z = 0.5*(pa.z + pb.z);
  
  // u = a - f; v = c - f;
  pu = MINUS(pa, pf);
  pv = MINUS(pc, pf);
  
  // d = (u X v) X u;
  pt = CROSS(pu, pv);
  pd = CROSS(pt, pu);
  
  // gama = (v^2 - u^2) / (2 d \dot (v - u)); // this is the correct form of point_pair pdf
  gama = -(SIZE(pu)*SIZE(pu) - SIZE(pv)*SIZE(pv));
  
  pt = MINUS(pv, pu);
  lamda = 2*(pd.x*pt.x + pd.y*pt.y + pd.z*pt.z);
  
  gama /= lamda;
  
  if (gama<0) {
    lamda = 0;
  }
  else {
    lamda = gama;
  }
  
  fem->qvec[ln_v].x = pf.x + lamda*pd.x;
  fem->qvec[ln_v].y = pf.y + lamda*pd.y;
  fem->qvec[ln_v].z = pf.z + lamda*pd.z;

  fem->radvec[ln_v] = SIZE(MINUS(fem->qvec[ln_v], pa));   
 
  return 0;
}

//------------------------------------------------------------------------------------------------------------
PetscInt  fem_ISPointInTriangle(Cmpnts p,Cmpnts p1,Cmpnts p2,Cmpnts p3,
			       PetscReal nfx, PetscReal nfy, PetscReal nfz) {
  
  PetscInt  flag;
  Cpt2D	    pj, pj1, pj2, pj3;
  if (fabs(nfz)>=fabs(nfx) && fabs(nfz)>=fabs(nfy)) {
    pj.x =  p.x;  pj.y = p.y;
    pj1.x = p1.x;  pj1.y = p1.y;
    pj2.x = p2.x;  pj2.y = p2.y;
    pj3.x = p3.x;  pj3.y = p3.y;
  }
  else if (fabs(nfx)>=fabs(nfy) && fabs(nfx)>=fabs(nfz)) {
    pj.x = p.z;  pj.y = p.y;
    pj1.x = p1.z;  pj1.y = p1.y;
    pj2.x = p2.z;  pj2.y = p2.y;
    pj3.x = p3.z;  pj3.y = p3.y;
  }
  else {
    pj.x = p.x;  pj.y = p.z;
    pj1.x = p1.x;  pj1.y = p1.z;
    pj2.x = p2.x;  pj2.y = p2.z;
    pj3.x = p3.x;  pj3.y = p3.z;
  }
  flag = fem_ISInsideTriangle2D(pj, pj1, pj2, pj3);
  
  return(flag);
}

//------------------------------------------------------------------------------------------------------------
PetscErrorCode  Sign_Dis_P_Line(Cmpnts p, Cmpnts p1, Cmpnts p2, PetscReal nfx, PetscReal nfy, PetscReal nfz, Cmpnts *po, PetscReal *d) {

  PetscReal  dmin;
  PetscReal  dx21, dy21, dz21, dx31, dy31, dz31, t;

  *d = 1.e6;
  
  dx21 = p2.x - p1.x;  dy21 = p2.y - p1.y;  dz21 = p2.z - p1.z;
  dx31 = p.x  - p1.x;  dy31 = p.y  - p1.y;  dz31 = p.z  - p1.z;
  
  t = (dx31*dx21 + dy31*dy21 + dz31*dz21)/(dx21*dx21 + dy21*dy21 +
						   dz21*dz21);
  if (t<0) { // The closet point is p1
    po->x = p1.x;  po->y = p1.y;  po->z = p1.z;
    *d = sqrt(dx31*dx31 + dy31*dy31 + dz31*dz31);
    
  } else if (t>1) { // The closet point is p2
    po->x = p2.x;  po->y = p2.y;  po->z = p2.z;
    *d = sqrt((p.x - po->x)*(p.x - po->x)+(p.y - po->y)*(p.y - po->y) +
  	      (p.z - po->z)*(p.z - po->z));
  
  } else { // The closet point lies between p1 & p2
    po->x = p1.x + t*dx21;  po->y = p1.y + t*dy21;  po->z = p1.z + t*dz21;
    *d = sqrt((p.x - po->x)*(p.x - po->x) + (p.y - po->y)*(p.y - po->y) +
  	      (p.z - po->z)*(p.z - po->z));
  }

  if (nfx*(p.x-po->x)+nfy*(p.y-po->y)+nfz*(p.z-po->z)<0) {*d = -*d;}

  return(0);
}

//------------------------------------------------------------------------------------------------------------
PetscInt fem_ISInsideTriangle2D(Cpt2D p, Cpt2D pa, Cpt2D pb, Cpt2D pc) {

  // Check if point p and p3 is located on the same side of line p1p2
  PetscInt  ls;
  
  ls = fem_ISSameSide2D(p, pa, pb, pc);
  //  if (flagprint) PetscPrintf(PETSC_COMM_WORLD, "aaa, %d\n", ls);
  if (ls < 0) {
    return (ls);
  }
  ls = fem_ISSameSide2D(p, pb, pc, pa);
  //  if (flagprint) PetscPrintf(PETSC_COMM_WORLD, "bbb, %d\n", ls);
  if (ls < 0) {
    return (ls);
  }
  ls = fem_ISSameSide2D(p, pc, pa, pb);
  //  if (flagprint) PetscPrintf(PETSC_COMM_WORLD, "ccc, %d\n", ls);
  if (ls <0) {
    return(ls);
  }

  return (ls);
}

//------------------------------------------------------------------------------------------------------------
PetscInt  fem_ISSameSide2D(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3) {

 /* Check whether 2D point p is located on the same side of line p1p2
     with point p3. Returns:
     -1	different side
     1	same side (including the case when p is located
     right on the line)
     If p and p3 is located on the same side to line p1p2, then
     the (p-p1) X (p2-p1) and (p3-p1) X (p2-p1) should have the same sign
  */
  PetscReal  t1, t2, t3;
  PetscReal  epsilon=1.e-10;
  PetscReal  A, B, C;

  A = p2.y - p1.y;
  B = -(p2.x - p1.x);
  C = (p2.x - p1.x) * p1.y - (p2.y - p1.y) * p1.x;
  
  t3 = fabs(A * p.x + B * p.y + C) / sqrt(A*A + B*B);
  
  /*   if (t3<1.e-3) return(1); */
  if (t3 < 1.e-3) {
    t1 = A * p.x + B * p.y + C;
    t2 = A * p3.x + B * p3.y + C;
    //    if (flagprint) PetscPrintf(PETSC_COMM_WORLD, "%le %le %le %le %le %le\n", t1, t2, t3, A, B, C);
  }
  else {
    t1 = (p.x - p1.x) * (p2.y - p1.y) - (p.y - p1.y) * (p2.x - p1.x);
    t2 = (p3.x - p1.x) * (p2.y - p1.y) - (p3.y - p1.y) * (p2.x - p1.x);
  }
  
  //!!!!!!!!!!!!1 Change t1, t2 & lt !!!!!!!
  t1 = (p.x - p1.x) * (p2.y - p1.y) - (p.y - p1.y) * (p2.x - p1.x);
  t2 = (p3.x - p1.x) * (p2.y - p1.y) - (p3.y - p1.y) * (p2.x - p1.x);
  PetscReal  lt;
  lt = sqrt((p1.x - p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y));
  //  if(flagprint) PetscPrintf(PETSC_COMM_WORLD, "%le %le %le %le %le %le\n", p1.x, p2.x, p3.x, p1.y, p2.y, p3.y);
  //if (fabs(t1) < epsilon) { // Point is located along the line of p1p2
  if (fabs(t1/lt) < epsilon) { // Point is located along the line of p1p2
    return(1);
  }
  // End of change !!!!!!!!!!!!!1
  
  if (t1 > 0) {
    if (t2 > 0) return (1); // same side
    else return(-1);  // not
  }
  else {
    if (t2 < 0) return(1); // same side
    else return(-1);
  }
}
