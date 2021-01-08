/*!
 * \file CFilmVariable.cpp
 * \brief Source of generic thin film variable
 * \author J. Li Volsi
 * \version 7.0.0 "Blackbird"
*/

#include "../include/variables/CFilmVariable.hpp"

/*--- CFILMVARIABLE METHODS ---*/

CFilmVariable::CFilmVariable(su2double density, su2double thickness, const su2double* velocity, su2double val_layer_toll, 
                            unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config, unsigned short val_iLayer):
                             CVariable(npoint,ndim,nvar,config), layer_toll(val_layer_toll) {

/*--- Initializing simple members ---*/
nADim = nDim - 1;
nPrimVar = nADim + 3; // + 4 if we include w
nPrimGrad = nADim; // nDim if we include w
LayerLabel = val_iLayer;
bottom = false;
kth = false;
top = false;
nLayer = config->GetnLayer();

if(LayerLabel == 1)
 bottom = true;
if(LayerLabel == nLayer)
 top = true; 
if(LayerLabel != 1 && LayerLabel != nLayer)
 kth = true;

/*--- Solution initialization (p, h, q_x, q_y) ---*/ 
su2double val_sol[4] = {su2double(0.5),su2double(1.0),velocity[0],velocity[0]};
if(nADim == 2) val_sol[3] = velocity[1];

for(unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
 for (unsigned long iVar = 0; iVar < nVar; iVar++){
  Solution(iPoint,iVar)     = density*val_sol[iVar];
  Solution(iPoint,iVar)     = Solution(iPoint,iVar)*thickness;
 }
 Solution_Old = Solution;

bool classical_rk4 = (config->GetKind_TimeIntScheme_Flow() == CLASSICAL_RK4_EXPLICIT);
if(classical_rk4){
 Solution_New.resize(nPoint,nVar) = su2double(0.0);
 Solution_New = Solution;
}

/*--- Numerics initialization ---*/
Res_TruncError.resize(nPoint,nVar) = su2double(0.0);
Delta_Time.resize(nPoint)          = su2double(0.0);
LocalCFL.resize(nPoint)            = su2double(0.0);
Max_Lambda_Inv.resize(nPoint)      = su2double(0.0);
Max_Lambda_Visc.resize(nPoint)     = su2double(0.0);

/*--- Incompressible flow, primitive variables nADim+3, (p,u,v,rho,delta) ---*/
Primitive.resize(nPoint,nPrimVar) = su2double(0.0);
Viscosity.resize(nPoint)          = su2double(0.0);
for(unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
 Viscosity(iPoint) = config->GetViscosity_FreeStreamND();

BottomTopography.resize(nPoint)   = su2double(0.0);
Bottom_BC.resize(nPoint,nPrimVar) = su2double(0.0);

} // end function 

CFilmVariable::CFilmVariable(su2double density, su2double thickness,const su2double* velocity, unsigned long npoint, 
               unsigned long ndim, unsigned long nvar, CConfig *config, unsigned short val_iLayer): CVariable(npoint,ndim,nvar,config){

/*--- Initializing simple members ---*/
nADim      = nDim - 1;
nPrimVar   = nADim + 3; // + 4 if we include w
nPrimGrad  = nADim; // ndim if we include w
LayerLabel = val_iLayer;
bottom     = false;
kth        = false;
top	   = false;
nLayer     = config->GetnLayer();

if(LayerLabel == 1)
 bottom = true;
if(LayerLabel == nLayer)
 top = true;
if(LayerLabel != 1 && LayerLabel != nLayer)
 kth = true;

/*--- Solution initialization (p, h, q_x, q_y) ---*/ 
su2double val_sol[4] = {su2double(0.5),su2double(1.0),velocity[0],velocity[0]};
if(nADim == 2) val_sol[3] = velocity[1];

for(unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
 for (unsigned long iVar = 0; iVar < nVar; iVar++){
  Solution(iPoint,iVar)     = density*val_sol[iVar];
  Solution(iPoint,iVar)     = Solution(iPoint,iVar)*thickness;
 }
 Solution_Old = Solution;

bool classical_rk4 = (config->GetKind_TimeIntScheme_Flow() == CLASSICAL_RK4_EXPLICIT);
if(classical_rk4){
 Solution_New.resize(nPoint,nVar) = su2double(0.0);
 Solution_New = Solution;
}

/*--- Numerics initialization ---*/
Res_TruncError.resize(nPoint,nVar) = su2double(0.0);
Delta_Time.resize(nPoint)          = su2double(0.0);
LocalCFL.resize(nPoint)            = su2double(0.0);
Max_Lambda_Inv.resize(nPoint);
Max_Lambda_Visc.resize(nPoint);

/*--- Incompressible flow, primitive variables nDim+3, (p,u,v,rho,delta) ---*/
Primitive.resize(nPoint,nPrimVar)  = su2double(0.0);
Viscosity.resize(nPoint)          = su2double(0.0);
for(unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
 Viscosity(iPoint) = config->GetViscosity_FreeStreamND();

BottomTopography.resize(nPoint)   = su2double(0.0);
Bottom_BC.resize(nPoint,nPrimVar) = su2double(0.0);

}


CFilmVariable::~CFilmVariable(){} // end function


su2double* CFilmVariable::GetLayerPrimitive(unsigned long iPoint, unsigned short iLayer){
if(iLayer == (nLayer-1) )
 return Primitive[iPoint];
else
 return NULL;
} // end function

su2double CFilmVariable::GetLayerPrimitive(unsigned long iPoint, unsigned short iVar, unsigned short iLayer){
 if( iLayer == (nLayer-1) )
  return Primitive(iPoint,iVar);
 else
  return 0.0;
} // end function

su2double CFilmVariable::GetLayerSolution_New(unsigned long iPoint, unsigned long iVar, unsigned short iLayer ) const{
 if( iLayer == (nLayer - 1) )
  return Solution_New(iPoint,iVar);
 else
  return 0.0;
} // end function

void CFilmVariable::ReconstructVelocity(unsigned long iPoint){
  for(unsigned short iDim = 0; iDim < nADim; iDim++){
   if(bottom){
    Primitive(iPoint,iDim+1) = GetBottom_BC(iPoint, iDim+1); cout<<"Reconstructed iPoint "<<iPoint<<" for Layer "<<LayerLabel<<endl;
   } else {
    Primitive(iPoint,iDim+1) = GetLayerPrimitive(iPoint, iDim+1, LayerLabel-1); cout<<"Reconstructed iPoint "<<iPoint<<" for Layer "<<LayerLabel<<endl;
   }
  }
} // end function

void CFilmVariable::SetLayerSolution_Old(unsigned long iPoint, su2double* solution, unsigned short iLayer){
 if(iLayer == (nLayer-1)){
  for(unsigned long iVar = 0; iVar < nVar; iVar++)
   this->SetSolution_Old(iPoint,iVar,solution[iVar]);
 }
} // end function

void CFilmVariable::AddSolution_New(unsigned long iPoint, unsigned long iVar, su2double solution){
  Solution_New(iPoint, iVar) += solution;
} // end function

void CFilmVariable::SetLayerSolution(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer){
 if(iLayer == (nLayer-1))
  this->SetSolution(iPoint,iVar,solution);
} // end function

void CFilmVariable::SetLayerSolution_New(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer){
 if(iLayer == (nLayer-1))
  CFilmVariable::SetSolution_New(iPoint, iVar, solution);
} // end function

void CFilmVariable::AddLayerSolution(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer){
 if(iLayer == (nLayer-1))
  AddSolution(iPoint,iVar,solution);
} // end function

void CFilmVariable::AddLayerSolution_New(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer){
 if(iLayer == (nLayer-1))
  CFilmVariable::AddSolution_New(iPoint,iVar,solution);
} // end function

void CFilmVariable::SetLayerPrimitive(unsigned long iPoint, unsigned short iVar, su2double prim_val, unsigned short iLayer){
 if(iLayer == (nLayer-1))
  SetPrimitive(iPoint, iVar, prim_val);
}

void CFilmVariable::ImportBottom_Topography(su2double** topography){
 
 unsigned long iPoint;
  for(iPoint = 0; iPoint < nPoint; iPoint++)
   BottomTopography(iPoint) = topography[iPoint][nADim];

} // end function


void CFilmVariable::InitializeBottom_BC(CConfig *config){

 unsigned long iPoint;
 unsigned short iDim;
 
// su2double press      = config->GetPressure_FreeStreamND();
 su2double density    = config->GetDensity_FreeStreamND();
 su2double bott_press = su2double(0.0);

 for(iPoint = 0; iPoint < nPoint; iPoint++){
 /*--- Initialize bottom pressure ---*/
//  bott_press = press - 0.5*density*BottomTopography(iPoint)/config->GetLength_Ref();
  bott_press = config->GetPressure_FreeStreamND();
  SetBottom_BC(iPoint, 0, bott_press);

 /*--- Initialize bottom velocity ---*/
 /* Hp: no-slip conditions, keep initialization to zero. */
 for(iDim = 0; iDim < nADim; iDim++)
   SetBottom_BC(iPoint, iDim+1, 0.0);

 /*--- Initialize bottom density ---*/
 /* Hp: no normal variation */
  SetBottom_BC(iPoint, nADim+1, density);
 }

} // end function


bool CFilmVariable::SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel){

 unsigned long iVar;
 bool check_thick = false, check_dens = false, physical = true;

 /*--- Setting Primitive Variables ---*/
 SetPressure(iPoint);

 check_dens = SetDensity(iPoint, FluidModel->GetDensity());
 SetLaminarViscosity(iPoint, FluidModel->GetLaminarViscosity());

 su2double thick = Solution(iPoint,1)/GetDensity(iPoint);
 check_thick = SetThickness(iPoint,thick);
 if(!check_thick)
  SetVelocity(iPoint);
 else
  ReconstructVelocity(iPoint);

/*--- Non-physical solution found. Revert to old values. ---*/
 if (check_dens) {
  /*--- Copy the old solution ---*/
  for (iVar = 0; iVar < nVar; iVar++){
//   Solution_New(iPoint, iVar) = Solution(iPoint, iVar);
   Solution(iPoint, iVar)     = Solution_Old(iPoint, iVar);
   }
  /*--- Recompute the primitive variables ---*/
  check_dens = SetDensity(iPoint, FluidModel->GetDensity());
  /*--- Flag this point as non-physical. ---*/ 
  physical = false; 
 }

 return physical;
} // end function


/*-------------- CMULTILAYERFILMVARIABLE METHODS -----------------*/

CMultilayerFilmVariable::CMultilayerFilmVariable(su2double* density, su2double* thickness, const su2double* velocity, 
     su2double val_layer_toll,unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config):
     CFilmVariable( density[config->GetnLayer()-1],thickness[config->GetnLayer()-1],velocity,val_layer_toll,npoint,ndim,nvar,
                    config,config->GetnLayer() ){

KthLayer = NULL;
/*--- Initializing bottom and kth layers ---*/
if((nLayer - 1) > 0 ){
 KthLayer = new CFilmVariable*[nLayer-1];
 for(unsigned short iLayer = 0; iLayer < (nLayer-1); iLayer++){
  KthLayer[iLayer] = new CFilmVariable(density[iLayer],thickness[iLayer],velocity,npoint,ndim,nvar,config,iLayer+1);
  KthLayer[iLayer]->SetLayerToll(val_layer_toll);
 }
}

} // function end


CMultilayerFilmVariable::~CMultilayerFilmVariable(){

 for(unsigned short iLayer = 0; iLayer < (nLayer-1); iLayer++)
  delete KthLayer[iLayer];

 if(KthLayer != NULL)
  delete [] KthLayer;

} // end function 


su2double* CMultilayerFilmVariable::GetLayerPrimitive(unsigned long iPoint, unsigned short iLayer){
 if( iLayer == (nLayer - 1) )
  return this->CFilmVariable::GetPrimitive(iPoint);
 else if( iLayer < (nLayer - 1) )
   return KthLayer[iLayer]->CFilmVariable::GetPrimitive(iPoint);
 else
  return NULL;
} // end function


su2double CMultilayerFilmVariable::GetLayerPrimitive(unsigned long iPoint, unsigned short iVar, unsigned short iLayer){
 if(iLayer == (nLayer - 1) )
  return this->CFilmVariable::GetPrimitive(iPoint, iVar);
 else if( iLayer < (nLayer - 1) )
  return KthLayer[iLayer]->CFilmVariable::GetPrimitive(iPoint, iVar);
 else
  return 0.0;
} // end function


su2double* CMultilayerFilmVariable::GetLayerSolution(unsigned long iPoint, unsigned short iLayer){
 if(iLayer == (nLayer - 1))
  return this->CVariable::GetSolution(iPoint);
 else
  return KthLayer[iLayer]->CVariable::GetSolution(iPoint);
} // end function


su2double CMultilayerFilmVariable::GetLayerSolution_New(unsigned long iPoint, unsigned long iVar, unsigned short iLayer ) const{
 if(iLayer == (nLayer-1) )
  return this->CFilmVariable::GetSolution_New(iPoint,iVar);
 else
  return KthLayer[iLayer]->CFilmVariable::GetSolution_New(iPoint, iVar);
} // end function


bool CMultilayerFilmVariable::SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel){
 
 unsigned long iVar;
 bool check_thick = false, check_dens = false, physical = true;

 /*--- Setting Primitive Variables ---*/
 SetPressure(iPoint);
 check_dens = SetDensity(iPoint, FluidModel->GetDensity());
 SetLaminarViscosity(iPoint, FluidModel->GetLaminarViscosity());

 su2double thick = Solution(iPoint,1)/GetDensity(iPoint);
 check_thick = SetThickness(iPoint,thick);
 if(!check_thick)
  SetVelocity(iPoint);
 else
  ReconstructVelocity(iPoint);

/*--- Non-physical solution found. Revert to old values. ---*/
 if (check_dens) {
  /*--- Copy the old solution ---*/
  for (iVar = 0; iVar < nVar; iVar++)
//   Solution_New(iPoint, iVar) = Solution(iPoint, iVar);
   Solution(iPoint, iVar) = Solution_Old(iPoint, iVar);
  /*--- Recompute the primitive variables ---*/
   check_dens = SetDensity(iPoint, FluidModel->GetDensity());
  /*--- Flag this point as non-physical. ---*/ 
   physical = false; 
 }


 bool chk = true;
 for(unsigned short iLayer = 0; iLayer < (nLayer-1); iLayer++)
   chk = KthLayer[iLayer]->CFilmVariable::SetPrimVar(iPoint,FluidModel);

 return physical;
} // end function


void CMultilayerFilmVariable::SetLayerSolution_Old(unsigned long iPoint, su2double* solution, unsigned short iLayer){
 if(iLayer == (nLayer-1)){
  for(unsigned long iVar = 0; iVar < nVar; iVar++)
   this->CVariable::SetSolution_Old(iPoint,iVar,solution[iVar]);
 } else {
  for(unsigned long iVar = 0; iVar < nVar; iVar++)
   KthLayer[iLayer]->CVariable::SetSolution_Old(iPoint,iVar,solution[iVar]);
 }
} // end function

void CMultilayerFilmVariable::AddLayerSolution(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer){
 if(iLayer == (nLayer -1))
   this->CVariable::AddSolution(iPoint,iVar,solution);
 else
  KthLayer[iLayer]->CVariable::AddSolution(iPoint,iVar,solution);
} // end function

void CMultilayerFilmVariable::AddLayerSolution_New(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer){
 if(iLayer == (nLayer-1) )
  this->CFilmVariable::AddSolution_New(iPoint,iVar,solution);
 else
  KthLayer[iLayer]->CFilmVariable::AddSolution_New(iPoint,iVar,solution);
} // end function

void CMultilayerFilmVariable::SetSolution_New() {
 Solution_New = Solution;
 for(unsigned short iLayer = 0; iLayer < (nLayer-1); iLayer++)
  for(unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
   su2double* sol = KthLayer[iLayer]->GetLayerSolution(iPoint,iLayer);
   for(unsigned short iVar = 0; iVar < nVar; iVar++){
    KthLayer[iLayer]->CFilmVariable::SetSolution_New(iPoint,iVar,sol[iVar]);
   }
  }
} // end function

void CMultilayerFilmVariable::SetLayerPrimitive(unsigned long iPoint, unsigned short iVar, su2double prim_val, unsigned short iLayer) {
 if(iLayer == (nLayer-1) )
  this->CFilmVariable::SetPrimitive(iPoint, iVar, prim_val);
 else
  KthLayer[iLayer]->CFilmVariable::SetPrimitive(iPoint, iVar, prim_val);
} // end function

void CMultilayerFilmVariable::SetLayerSolution(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer){
 if(iLayer == (nLayer-1))
  this->CVariable::SetSolution(iPoint, iVar, solution);
 else
  KthLayer[iLayer]->CVariable::SetSolution(iPoint, iVar, solution);
} // end function

void CMultilayerFilmVariable::SetLayerSolution_New(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer){
 if(iLayer == (nLayer-1))
  this->CFilmVariable::SetSolution_New(iPoint,iVar,solution);
 else
  KthLayer[iLayer]->CFilmVariable::SetSolution_New(iPoint,iVar,solution);
} // end function

CFilmVariable* CMultilayerFilmVariable::GetLayer(unsigned short iLayer){
 if(iLayer == (nLayer-1))
  return this;
 else 
  return KthLayer[iLayer];
} // end function



