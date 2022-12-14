PetscErrorCode ibm_search_advanced_rev(UserCtx *user, IBMNodes *ibm, 
				       PetscInt ibi)

/*      Note : Always go from ibi (immersed body number) 0 -> NumberOfBodies  */
/*             Nvert should be set to zero before any new search  */
/*              */
{
// Read in the Distributed Array; da for variables and fda for co-ordinates.
  DM	da = user->da, fda = user->fda;
// Preprocessing for PETSC.
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscInt	ncx = 40, ncy = 40, ncz = 40; // Setting up the number of control cells.
  List          *cell_trg;
  PetscReal	xbp_min, ybp_min, zbp_min, xbp_max, ybp_max, zbp_max; // These are max and min co-ordinate values of the bounding box in all three directions.
  PetscReal	*x_bp = ibm->x_bp, *y_bp = ibm->y_bp, *z_bp = ibm->z_bp; // Reading in the co-ordinate arrays of the immersed boundary nodes.

  PetscInt 	ln_v, n_v = ibm->n_v; // Reading in the number of nodes in the immersed boundary.

  PetscInt	i, j, k;

  PetscReal	dcx, dcy, dcz;
  PetscInt	n1e, n2e, n3e;
  PetscReal	xv_min, yv_min, zv_min, xv_max, yv_max, zv_max,xg_min,xg_max,yg_min,yg_max,zg_min,zg_max;
  PetscInt	iv_min, iv_max, jv_min, jv_max, kv_min, kv_max;
  PetscReal	***nvert;
  PetscInt	ic, jc, kc;
  Cmpnts 	***coor;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscPrintf(PETSC_COMM_WORLD, "mz,my,mx: %d, %d, %d \n",mz,my,mx);
  xbp_min = 1.e23; xbp_max = -1.e23;
  ybp_min = 1.e23; ybp_max = -1.e23;
  zbp_min = 1.e23; zbp_max = -1.e23;
// Finding the minimum co-ordinate in all three directions.
  PetscPrintf(PETSC_COMM_WORLD, "%i\n", n_v); 
    for(i=0; i<n_v; i++) {
     //PetscPrintf(PETSC_COMM_WORLD, "%i\n", i); 
     //PetscPrintf(PETSC_COMM_WORLD, "%e\n", x_bp[i]); 
    xbp_min = PetscMin(xbp_min, x_bp[i]);
    xbp_max = PetscMax(xbp_max, x_bp[i]);

    ybp_min = PetscMin(ybp_min, y_bp[i]);
    ybp_max = PetscMax(ybp_max, y_bp[i]);

    zbp_min = PetscMin(zbp_min, z_bp[i]);
    zbp_max = PetscMax(zbp_max, z_bp[i]);
  }
//--------------------------------------------------
//  DMDAVecGetArray(fda, user->Cent, &coor); // Reading fluid co-ordinates local vector into an array.
//-----------------------------------------------------------------------------------
  xg_min = 1.e23; xg_max = -1.e23;
  yg_min = 1.e23; yg_max = -1.e23;
  zg_min = 1.e23; zg_max = -1.e23;
/*
for(k=zs;k<ze;k++){
   for(j=ys;j<ye;j++){
      for(i=xs;i<xe;i++){
	                
	   xg_min = PetscMin(xg_min,coor[k][j][i].x);
	   xg_max = PetscMax(xg_max,coor[k][j][i].x);

	   yg_min = PetscMin(yg_min,coor[k][j][i].y);
	   yg_max = PetscMax(yg_max,coor[k][j][i].y);

	   zg_min = PetscMin(zg_min,coor[k][j][i].z);
	   zg_max = PetscMax(zg_max,coor[k][j][i].z);

	 }
      }
 	   if(k == mz-1){
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"k,j,i,z(k),zg_max: %d,%d,%d,%le, %le \n",k,j,i,coor[k][j][i].z,zg_max);
	}
 
 }
PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
*/
//PetscPrintf(PETSC_COMM_WORLD,"zg_max: %le \n",zg_max);
// Adding/subtracting a buffer value to determine the dimensions of the bounding box.
  xbp_min -= 0.01; xbp_max += 0.01;
  ybp_min -= 0.01; ybp_max += 0.01;
  zbp_min -= 0.01; zbp_max += 0.01;
//zbp_min -= 0.09; zbp_max += 0.09;

// Printing the bounding box dimensions.
/*
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "xbp min max  %le %le \n",xbp_min,xbp_max);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "ybp min max  %le %le \n",ybp_min,ybp_max);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "zbp min max  %le %le \n",zbp_min,zbp_max);
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
//
*/
/* 
 PetscSynchronizedPrintf(PETSC_COMM_WORLD, "xg min max  %le %le \n",xg_min,xg_max);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "yg min max  %le %le \n",yg_min,yg_max);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD, "zg min max %le %le \n",zg_min,zg_max);
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

  xbp_min = PetscMax(xg_min,xbp_min);
  xbp_max = PetscMin(xg_max,xbp_max);

  ybp_min = PetscMax(yg_min,ybp_min);
  ybp_max = PetscMin(yg_max,ybp_max);

  zbp_min = PetscMax(zg_min,zbp_min);
  zbp_max = PetscMin(zg_max,zbp_max);
*/
// Determining the dimensions of each control cell.
  dcx = (xbp_max - xbp_min) / (ncx - 1.);
  dcy = (ybp_max - ybp_min) / (ncy - 1.);
  dcz = (zbp_max - zbp_min) / (ncz - 1.);

// Printing the bounding box dimensions
  PetscPrintf(PETSC_COMM_WORLD," After correcting with fluid grid dimensions \n");
  PetscPrintf(PETSC_COMM_WORLD, "xbp min max dcx %le %le %le\n",xbp_min,xbp_max,dcx);
  PetscPrintf(PETSC_COMM_WORLD, "ybp min max dcy %le %le %le\n",ybp_min,ybp_max,dcy);
  PetscPrintf(PETSC_COMM_WORLD, "zbp min max dcz %le %le %le\n",zbp_min,zbp_max,dcz);
// Memory allocation for an array containing one "List" object per control cell.
  PetscMalloc(ncz * ncy * ncx * sizeof(List), &cell_trg);
  //PetscPrintf(PETSC_COMM_SELF, "test00\n");
// Initializing the list for each control cell accessed through a 3D index, although the list is actually 1D. 
 for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {
	initlist(&cell_trg[k*ncx*ncy + j*ncx + i]);
      }
    }
  }
  PetscPrintf(PETSC_COMM_WORLD, "test0\n");
// For each individual element. Finding it's minimum and maximum co-ordinate among all three nodes it is made up of.
    for (ln_v=0; ln_v < ibm->n_elmt; ln_v++) {

    n1e = ibm->nv1[ln_v]; n2e = ibm->nv2[ln_v]; n3e = ibm->nv3[ln_v]; // Reading the three nodes of an element.

    xv_min = PetscMin(PetscMin(x_bp[n1e], x_bp[n2e]), x_bp[n3e]); // Comparision for min along x
    xv_max = PetscMax(PetscMax(x_bp[n1e], x_bp[n2e]), x_bp[n3e]); // "" ""   ""  ""  max ""   ""

    yv_min = PetscMin(PetscMin(y_bp[n1e], y_bp[n2e]), y_bp[n3e]);
    yv_max = PetscMax(PetscMax(y_bp[n1e], y_bp[n2e]), y_bp[n3e]);

    zv_min = PetscMin(PetscMin(z_bp[n1e], z_bp[n2e]), z_bp[n3e]);
    zv_max = PetscMax(PetscMax(z_bp[n1e], z_bp[n2e]), z_bp[n3e]);
    
// Calculating how many control cells are behind the minimum and maximum co-ordinate along each axis.
    iv_min = floor((xv_min - xbp_min) / dcx); //  +1???
    iv_max = floor((xv_max - xbp_min) / dcx) +1;

    jv_min = floor((yv_min - ybp_min) / dcy); //  +1???
    jv_max = floor((yv_max - ybp_min) / dcy) +1;

    kv_min = floor((zv_min - zbp_min) / dcz); //  +1???
    kv_max = floor((zv_max - zbp_min) / dcz) +1;
// Handling exceptions (if the element is the one at the edge end of the immersed object, then ivmin would be 0 but may be calculated as -1. 
    iv_min = (iv_min<0) ? 0:iv_min;
    iv_max = (iv_max>ncx) ? ncx:iv_max;

    jv_min = (jv_min<0) ? 0:jv_min;
    jv_max = (jv_max>ncx) ? ncy:jv_max;

    kv_min = (kv_min<0) ? 0:kv_min;
    kv_max = (kv_max>ncz) ? ncz:kv_max;

//        if (ln_v==25) {
//      PetscPrintf(PETSC_COMM_WORLD, "id25, %d %d %d %d %d %d\n", iv_min, iv_max, jv_min, jv_max, kv_min, kv_max);
//      }
    // Insert IBM node information into a list dedicated for each control cell.
    for (k=kv_min; k<kv_max; k++) {
      for (j=jv_min; j<jv_max; j++) {
	for (i=iv_min; i<iv_max; i++) {
	  insertnode(&(cell_trg[k *ncx*ncy + j*ncx +i]), ln_v);
	}
      }
    }
  } // Element loop closed
 // DMDAVecRestoreArray(fda,user->Cent,&coor);
  PetscPrintf(PETSC_COMM_WORLD, "test001\n");
/*   List test; */
/*   insertnode(&test, 11); */
  PetscInt rank, flg=0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
// Reading the Nvert vector into arrays.
  DMDAVecGetArray(da, user->Nvert, &nvert);
  DMDAVecGetArray(fda,user->lCent,&coor); // Values of cell centers. 0 and IM-1 are ghost nodes.
// For Points INSIDE---------------------- 
  // Initially for this body nvert 4 is inside, 2 is near bndry
  // for previous bodies nvert 3 inside, 1 near bndry

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

//	if (k == 0){
//	PetscPrintf(PETSC_COMM_WORLD,"zbp_min,zbp_max, z(0,j,i),j,i : %le,%le, %le %d,%d \n",zbp_min,zbp_max,coor[k][j][i].z,j,i);
// 	PetscPrintf(PETSC_COMM_WORLD,"ybp_min,ybp_max, y(0,j,i),j,i : %le,%le, %le %d,%d \n",ybp_min,ybp_max,coor[k][j][i].y,j,i);
//	PetscPrintf(PETSC_COMM_WORLD,"xbp_min,xbp_max, x(0,j,i),j,i : %le,%le, %le %d,%d \n",xbp_min,xbp_max,coor[k][j][i].x,j,i);
//		}
//	if(k == mz-1){
// 	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"zbp_max, z(zm-1,j,i) : %le, %le \n",zbp_max,coor[k][j][i].z);
//	} 

	if (ibi==0) nvert[k][j][i] = 0; //reset nvert if new search.

// if the fluid grid point is inside the bounding box (is inside).
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Check whether boundaries are inside the bounding box.
	if (coor[k][j][i].x >= xbp_min && coor[k][j][i].x <= xbp_max &&
	    coor[k][j][i].y >= ybp_min && coor[k][j][i].y <= ybp_max &&
	    coor[k][j][i].z >= zbp_min && coor[k][j][i].z <= zbp_max) {

// ic,jc,kc represent the number of control cells behind the particular grid point in the three directions.
	  ic = floor((coor[k][j][i].x - xbp_min )/ dcx);
	  jc = floor((coor[k][j][i].y - ybp_min )/ dcy);
	  kc = floor((coor[k][j][i].z - zbp_min )/ dcz);

	if(k==0){
	//PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Inside bb: k,j,i,zcor,dcz,kc: %d, %d,%d,%le,%le,%d \n",k,j,i,coor[k][j][i].z,dcz,kc);
}
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD," Inside bounding box: k,j,i = %d, %d, %d \n",k,j,i);
	 
	  /* search only if the node is not inside another body 
	     already! */

	  nvert[k][j][i] =  PetscMax(nvert[k][j][i],
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Check point cell advanced behaviour when ic,jc,kc are zeros.
				     point_cell_advanced2(coor[k][j][i], ic, jc, kc, ibm, ncx, ncy, ncz, dcx, dcy, xbp_min, ybp_min, zbp_max, cell_trg, flg)); // Do the search (Ray Tracing) for these points, the output can be 4 (point is inside the immersed boundary) or 0.
	  nvert[k][j][i] -=4; // Subtract 4 from nvert -- This would make sure interior points have an nvert of 0 and exterior points have an invert <0.
	  if (nvert[k][j][i] < 0) nvert[k][j][i] = 4; // All points outside the immersed boundary but inside the bounding box are set to 4.
	} else 
	  nvert[k][j][i] = 4; // Everywhere outside the bounding box, set the nvert to be four.

      }
    }
  }
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

  //PetscPrintf(PETSC_COMM_SELF, "test01 %d\n",rank);
  DMDAVecRestoreArray(da, user->Nvert, &nvert);

/*   if (user->thislevel < user->mglevels-1) { */
/*     MyNvertRestriction(user->user_f, user); */

  //PetscPrintf(PETSC_COMM_SELF, "test010\n");

  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);

  DMDAVecGetArray(da, user->lNvert, &nvert);

  PetscPrintf(PETSC_COMM_WORLD, "test010\n");
  //PetscPrintf(PETSC_COMM_SELF, "test010\n");
  PetscInt ip, im, jp, jm, kp, km; // These are created to store points neighboring a point. (for i,j,k these are i+1,i-1,j+1,j-1,k+1,k-1).
  PetscInt ii, jj, kk;
//---------------------------------------------------
  for(k=zs;k<ze;k++){
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if(k==1){
//  	PetscSynchronizedPrintf(PETSC_COMM_WORLD, " nvert, k,j,i = %le, %d, %d, %d \n",nvert[k][j][i],k,j,i);
}

	if(k == mz-2){
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "zmax,nvert,i,j = %d,%le, %d, %d \n",k,nvert[k][j][i],i,j);
}
}
}
}
//PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);	

  //For points near boundary-------------------------
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
//	PetscPrintf(PETSC_COMM_WORLD, "test011\n");
	if (nvert[k][j][i] <0) nvert[k][j][i] = 0;
//		if(k == zm -2){
//		PetscSynchronizedPrintf(PETSC_COMM_WORLD, "nvert, i,j,k: %le, %d, %d, %d \n",nvert[k][j][i],i,j,k);
//}
//		nvert[k][j][i] = 0;
//	} // Setting the points inside bounding box but outside immersed boundary to 4.


// ip,jp and kp values are assigned to neughbouring nodes of i,j,k.
	ip = (i<mx-1?(i+1):(i));
	im = (i>0   ?(i-1):(i));

	jp = (j<my-1?(j+1):(j));
	jm = (j>0   ?(j-1):(j));

	kp = (k<mz-1?(k+1):(k));
	km = (k>0   ?(k-1):(k));
 //Check whether the node is inside the bounding box.
//	if (coor[k][j][i].x > xbp_min && coor[k][j][i].x < xbp_max &&
//	    coor[k][j][i].y > ybp_min && coor[k][j][i].y < ybp_max &&
//	    coor[k][j][i].z > zbp_min && coor[k][j][i].z < zbp_max) {

	if ((int)((nvert[k][j][i]) + 0.5) != 4) { // If the node is a fluid node. // 10/11 -> removed + 0.5 from inside the if condition.// 10/25 ->added 0.5 back.
	  for (kk=km; kk<kp+1; kk++) {
	    for (jj=jm; jj<jp+1; jj++) {
	      for (ii=im; ii<ip+1; ii++) {
		if ((int)(nvert[kk][jj][ii]+0.5) == 4) { // If any of the neighbouring points are solid.  // 10/11 -> removed + 0.5 from inside the if condition. // 10/25 ->added 0.5 back. 
	//	PetscPrintf(PETSC_COMM_WORLD,"i,j,k,nvert[k][j][i],nvert[kk][jj][ii]:%d,%d,%d,%f,%f \n",k,j,i,nvert[k][j][i],nvert[kk][jj][ii]);  
		nvert[k][j][i] = PetscMax(2, nvert[k][j][i]); //Set the point to be a boundary point.
		}
	      }
	    }
	  }

// 	  if (nvert[k][j][ip] == 3 || nvert[k][j][im]==3 || 
// 	      nvert[k][jp][i] == 3 || nvert[k][jm][i]==3 || 
// 	      nvert[kp][j][i] == 3 || nvert[km][j][i]==3) { 
// 	    nvert[k][j][i] = 1; }
	}
	}
	}	
}
//	PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
//}
//
//  PetscPrintf(PETSC_COMM_WORLD, "test111: zs,ze,ys,ye,xs,xe,mz-2 = %d,%d,%d,%d,%d,%d,%d\n",zs,ze,ys,ye,xs,xe,mz-2);
  for(k=zs;k<ze;k++){
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
//  	PetscPrintf(PETSC_COMM_WORLD, " nvert, k,j,i = %le, %d, %d, %d \n",nvert[k][j][i],k,j,i);


	if(k == 1){
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "nvert,k,j,i = %le, %d,%d, %d \n",nvert[k][j][i],k,j,i);
}
}
}
}
PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);	


  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if((int)(nvert[k][j][i]+0.5)==2){
//	PetscPrintf(PETSC_COMM_WORLD,"i,j,k,nvert:%d,%d,%d,%le \n",i,j,k,nvert[k][j][i]);
}
}
}
}

PetscPrintf(PETSC_COMM_WORLD,"****************************************\n");

//--------------------------------------------------
//Phase change
//---------------------------------------------------
  PetscBarrier(PETSC_NULL);
//  PetscPrintf(PETSC_COMM_WORLD, "test11\n");


  PetscReal	***nvert_o;
  DMDAVecGetArray(da, user->lNvert_o, &nvert_o);
  if (ibi==NumberOfBodies-1)
  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (nvert_o[k][j][i] >2.5 && nvert[k][j][i] < 0.5) {
	  PetscPrintf(PETSC_COMM_SELF, "Phase Change at %d, %d, %d!\n", i, j, k);
	  nvert[k][j][i]=2;
	}
      }
    }
  }
  DMDAVecRestoreArray(da, user->lNvert_o, &nvert_o);
  //  PetscPrintf(PETSC_COMM_WORLD, "test21\n");

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert);
  DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert);

// WHen multiple bodies are present!
//------------------------------------------------

  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);

  if (user->ibmlist[ibi].head) DestroyIBMList(&(user->ibmlist[ibi])); //List generated to note down immersed boundary mesh elements coinciding with each control cell is deleted.

  InitIBMList(&user->ibmlist[ibi]); // New list initialized.
  PetscInt number;
  number = 0;

  IBMInfo ibm_intp;

  DMDAVecGetArray(da, user->lNvert, &nvert);

  BoundingSphere(ibm); 

   PetscBarrier(PETSC_NULL);
   PetscPrintf(PETSC_COMM_WORLD, "test31\n"); 

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((int)(nvert[k][j][i]+0.5) == 2) {
	  PetscPrintf(PETSC_COMM_WORLD, " nvert: k,j,i, %le, %d, %d ,%d \n",nvert[k][j][i],k,j,i);	

	  number ++;
	  ic = (int)((coor[k][j][i].x - xbp_min) / dcx);
	  jc = (int)((coor[k][j][i].y - ybp_min) / dcy);
	  kc = (int)((coor[k][j][i].z - zbp_min) / dcz);

	  if (ic<0) ic=0;
	  else if (ic>=ncx) ic=ncx-1;

	  if (jc<0) jc=0;
	  else if (jc>=ncy) jc = ncy-1;

	  if (kc<0) kc=0;
	  else if (kc>=ncz) kc = ncz-1;

	  ibm_intp.ni = i;
	  ibm_intp.nj = j;
	  ibm_intp.nk = k;
//----------------------


// 	  if (ibi==0 && i==39 && j==41 && k==63) {
// 	  if (i==1 && j==100 && k==119) { 
// 	    PetscPrintf(PETSC_COMM_SELF, "IJK %d %d %d\n", i, j, k); 
// 	    flagprint =1; 
// 	  } 
	   nearestcell(coor[k][j][i], ibm, &ibm_intp);

//	 nearestcell1(coor[k][j][i], ibm, &ibm_intp);
	  //	  PetscPrintf(PETSC_COMM_WORLD, "nearest cell\n");
	  InterceptionPoint(coor[k][j][i], i, j, k, &ibm_intp, user);
	  //	  PetscPrintf(PETSC_COMM_WORLD, "inteception point\n");

//	  InterceptionPoint2(coor[k][j][i], i, j, k, &ibm_intp, user);

	  if (ibm_intp.imode<0) {
	    PetscInt cell;
	    Cmpnts ptmp;
	    if (i==1 || i==mx-2 ||
		j==1 || j==my-2) {

	      cell = ibm_intp.cell;
	      if (ibm->nf_z[cell] > 0) {
		ptmp = coor[k+1][j][i];
		ibm_intp.d_i = sqrt((coor[k][j][i].x - ptmp.x) *
				    (coor[k][j][i].x - ptmp.x) +
				    (coor[k][j][i].y - ptmp.y) *
				    (coor[k][j][i].y - ptmp.y) +
				    (coor[k][j][i].z - ptmp.z) *
				    (coor[k][j][i].z - ptmp.z));
		ibm_intp.cr1 = 1.;
		ibm_intp.cr2 = 0.;
		ibm_intp.cr3 = 0.;
		ibm_intp.i1 = i;
		ibm_intp.j1 = j;
		ibm_intp.k1 = k+1;

		ibm_intp.i2 = i;
		ibm_intp.j2 = j;
		ibm_intp.k2 = k+1;
		ibm_intp.i3 = i;
		ibm_intp.j3 = j;
		ibm_intp.k3 = k+1;
	      }
	      else {
		ptmp = coor[k-1][j][i];
		ibm_intp.d_i = sqrt((coor[k][j][i].x - ptmp.x) *
				    (coor[k][j][i].x - ptmp.x) +
				    (coor[k][j][i].y - ptmp.y) *
				    (coor[k][j][i].y - ptmp.y) +
				    (coor[k][j][i].z - ptmp.z) *
				    (coor[k][j][i].z - ptmp.z));
		ibm_intp.cr1 = 1.;
		ibm_intp.cr2 = 0.;
		ibm_intp.cr3 = 0.;
		ibm_intp.i1 = i;
		ibm_intp.j1 = j;
		ibm_intp.k1 = k-1;

		ibm_intp.i2 = i;
		ibm_intp.j2 = j;
		ibm_intp.k2 = k-1;
		ibm_intp.i3 = i;
		ibm_intp.j3 = j;
		ibm_intp.k3 = k-1;
	      }
	     }
	    else if (k==1 || k==mz-2) {
	      cell = ibm_intp.cell;
	      ptmp = coor[k][j+1][i];
	      ibm_intp.d_i = sqrt((coor[k][j][i].x - ptmp.x) *
				  (coor[k][j][i].x - ptmp.x) +
				  (coor[k][j][i].y - ptmp.y) *
				  (coor[k][j][i].y - ptmp.y) +
				  (coor[k][j][i].z - ptmp.z) *
				  (coor[k][j][i].z - ptmp.z));
	      ibm_intp.cr1 = 1.;
	      ibm_intp.cr2 = 0.;
	      ibm_intp.cr3 = 0.;
	      ibm_intp.i1 = i;
	      ibm_intp.j1 = j+1;
	      ibm_intp.k1 = k;
	      
	      ibm_intp.i2 = i;
	      ibm_intp.j2 = j+1;
	      ibm_intp.k2 = k;
	      ibm_intp.i3 = i;
	      ibm_intp.j3 = j+1;
	      ibm_intp.k3 = k;
	    }
	    else {
	      PetscPrintf(PETSC_COMM_SELF, "%%%%IBM Searching Fail! %d %d %d\n", i, j, k);
	           }

	  }
	  AddIBMNode(&user->ibmlist[ibi], ibm_intp);
      }
    }
  }
}

   PetscPrintf(PETSC_COMM_WORLD, "test 41 \n");


 //PetscBarrier(PETSC_NULL);

  //  PetscPrintf(PETSC_COMM_WORLD, "test22\n");
  //  PetscBarrier(PETSC_NULL);

  // Back to the old nvert 3 and 1 

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((int)(nvert[k][j][i]+0.5) == 2){ 
		nvert[k][j][i]=1.4;
		  PetscSynchronizedPrintf(PETSC_COMM_WORLD, " nvert: k,j,i, %le, %d, %d ,%d \n",nvert[k][j][i],k,j,i);	
		}
	if ((int)(nvert[k][j][i]+0.5) == 4) nvert[k][j][i]=3.4;
      }
    }
  }
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);	
  
   PetscPrintf(PETSC_COMM_WORLD, "test 51 \n");

  DMDAVecRestoreArray(fda, user->lCent,&coor);
   PetscPrintf(PETSC_COMM_WORLD, "test 61 \n");

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
    PetscPrintf(PETSC_COMM_WORLD, "test 71 \n");
 
  DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert);
    PetscPrintf(PETSC_COMM_WORLD, "test 81 \n");

 DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert);
   PetscPrintf(PETSC_COMM_WORLD, "test 91 \n");

//   PetscPrintf(PETSC_COMM_WORLD, "test23\n"); 

//   PetscBarrier(PETSC_NULL); 

//   PetscPrintf(PETSC_COMM_WORLD, "test24\n"); 


  DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
  DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);

  for (k=0; k<ncz; k++) {
    for (j=0; j<ncy; j++) {
      for (i=0; i<ncx; i++) {
	destroy(&cell_trg[k*ncx*ncy+j*ncx+i]);
      }
    }
  }

  PetscFree(cell_trg) ;
  PetscFree(ibm->qvec);
  PetscFree(ibm->radvec); 

  //PetscPrintf(PETSC_COMM_WORLD, "Interface pts blanked block!!!! %i,\n", block_number);


//   if (block_number>1)
//     Blank_Interface(user);

//   PetscBarrier(PETSC_NULL); 
//   PetscPrintf(PETSC_COMM_WORLD, "test25\n");

  //  PetscBarrier(PETSC_NULL);
 
 return 0;
}
