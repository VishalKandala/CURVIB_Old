#include "variables.h"

extern const PetscInt   dof;
extern const PetscReal  h0;
extern PetscInt         ti, BHV, inv_flg;
extern PetscReal        dt, Fnormal;
extern PetscErrorCode  EdgeFix(PetscInt edge_n,FE *fem);
extern PetscErrorCode  SurfaceNormalPressure(FE *fem);
extern PetscErrorCode  SurfaceNormalPressure2(FE *fem);
extern PetscErrorCode  GhostFix(PetscInt edge_n,FE *fem);
extern PetscErrorCode  GhostFree(PetscInt edge_n,FE *fem);
extern PetscErrorCode  ModifyGhostFix(PetscInt edge_n,FE *fem);
extern PetscErrorCode  ModifyGhostFree(PetscInt edge_n,FE *fem);
extern PetscErrorCode  NodeForce(PetscInt nv,PetscReal F,PetscInt dir,FE *fem);
extern PetscErrorCode  BuoyantForce(PetscReal F,PetscInt dir,FE *fem);


PetscErrorCode  GhostLoc(FE *fem) {

  /* GhostFree(0,fem); */
  /* GhostFree(1,fem); */

  return(0);
}

//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode ModifyFbending(FE *fem) {

  /* ModifyGhostFree(0,fem); */
  /* ModifyGhostFree(1,fem);   */

 return(0);
}
//------------------------------------------------------------------------------------------------------------ 
PetscErrorCode FExternal(FE *fem) {

  //----------------ExternalForce-Body--------------------------------
  SurfaceNormalPressure(fem);  
  
  //----------------ExternalForce-Boundary-----------------------------
  if (BHV) EdgeFix(1,fem);
  if (inv_flg) EdgeFix(1,fem); 
  
  //----------------ExternalForce-Point-----------------------------
  
  return(0);
}

