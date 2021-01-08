/*!
 * \file CNumerics_Film.cpp
 * \brief Source for the numerical scheme for thin film eq. 
 * \author J. Li Volsi
 * \version 7.0.0 "Blackbird"
*/

#include "../../include/numerics/CFilmNumerics.hpp"

/* ----- Definition of SW_Flux methods ----- */

su2double SW_Flux::operator() (su2double* U, unsigned short iEqs, unsigned short iDim){

 su2double result = 0.0;

 if(iEqs == 0){
  result = U[iDim+1];
 }

 if(iEqs == 1){
  if(iDim == 0){
    if(U[0] != 0.0) result = U[iEqs]*U[iEqs]/U[0] + 0.5*U[0]*U[0]/density;
  }
  if(iDim == 1){
   if(U[0] != 0.0) result = U[iEqs]*U[iEqs+1]/U[0]; 
  }
 }

 if(iEqs == 2){
  if(iDim == 0){ 
   if(U[0] != 0.0) result =  U[iEqs-1]*U[iEqs]/U[0];
  }
  if(iDim == 1){
   if(U[0] != 0.0) result = U[iEqs]*U[iEqs]/U[0] + 0.5*U[0]*U[0]/density;
  }
 }

 return result;

} // end function


/* ----- Definition of CCentUpw_SW methods ----- */

CCentUpw_SW::CCentUpw_SW(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config){

nADim = val_nDim - 1;
nEqs  = nADim + 1; // nADim Momentum eq + Continuty eq
dynamic_grid = config->GetDynamic_Grid();
Coord_i = NULL;
Coord_j = NULL;

Delta           = new su2double [nADim];
Velocity_i      = new su2double [nADim];
Velocity_j      = new su2double [nADim];
Conservatives_i = new su2double [nEqs];
Conservatives_j = new su2double [nEqs];
U_mean          = new su2double [nEqs];
U_i_card        = new su2double*[nEqs];
U_star          = new su2double*[nEqs];
delta_U         = new su2double*[nEqs];
numerical_flux  = new su2double*[nEqs];

for(iEqs = 0; iEqs < nEqs; iEqs++ ){
 U_i_card[iEqs] = new su2double[2*nADim];
 U_star[iEqs]   = new su2double[nADim];
 delta_U[iEqs]  = new su2double[nADim];
 numerical_flux[iEqs] = new su2double[nADim];
}

U_card = new su2double*[2*nADim];
for(iDim = 0; iDim < (2*nADim); iDim++)
 U_card[iDim] = new su2double[nEqs];

adv_coeff = new su2double[2*nADim];

g_const = STANDARD_GRAVITY*config->GetTime_Ref()*config->GetTime_Ref()/config->GetLength_Ref();
gravity_force = config->GetGravityForce();
if(gravity_force){
 body_force_vector = config->GetBody_Force_Vector();
 g_const = abs(body_force_vector[nADim]);
}
//cout<<"Upwind g_const= "<<g_const<<endl;

} // end function


CCentUpw_SW::~CCentUpw_SW(void){
 delete [] Velocity_i;
 delete [] Velocity_j;
 delete [] Conservatives_i;
 delete [] Conservatives_j;
 delete [] Delta;
 delete [] U_mean;
 delete [] adv_coeff;
 for(iEqs = 0; iEqs < nEqs; iEqs++ ){
  if(U_i_card[iEqs] != NULL) delete [] U_i_card[iEqs];
  if(U_star[iEqs]   != NULL) delete [] U_star[iEqs];
  if(delta_U[iEqs]  != NULL) delete [] delta_U[iEqs];
  if(numerical_flux[iEqs] != NULL) delete [] numerical_flux[iEqs];
 }
 delete [] U_i_card;
 delete [] U_star;
 delete [] delta_U;
 delete [] numerical_flux; 

 for(iDim = 0; iDim < 2*nADim; iDim++ )
  if(U_card[iDim] != NULL) delete [] U_card[iDim];
 delete [] U_card;

} // end function


void  CCentUpw_SW::SetCardinals(){

 Delta[0] = (Coord_j[0] - Coord_i[0])/2.0;
 if( nDim == 3 ) Delta[1] = (Coord_j[1] - Coord_i[1])/2.0;

 for(iEqs = 0; iEqs < nEqs; iEqs++){
   // Zero-order Polynomial
   U_mean[iEqs] = Conservatives_i[iEqs]; //(Conservatives_j[iEqs] + Conservatives_i[iEqs])/2.0; 
  }
 
 su2double U_derivative;
 for(iDim = 0; iDim < nADim; iDim++){
  for(iEqs = 0; iEqs < nEqs; iEqs++){
   // Forward derivative
   if( Delta[iDim] > epsy ){ // if( Delta[iDim] != 0.0 ) 
    U_derivative = (Conservatives_j[iEqs] - Conservatives_i[iEqs])/Delta[iDim];
    U_card[iDim][iEqs] = U_mean[iEqs] + Delta[iDim]*U_derivative;
    U_card[iDim+nADim][iEqs] = U_mean[iEqs] - Delta[iDim]*U_derivative;
   } else {
    U_card[iDim][iEqs] = U_mean[iEqs] + Conservatives_j[iEqs] - Conservatives_i[iEqs];
    U_card[iDim+nADim][iEqs] = U_mean[iEqs] - Conservatives_j[iEqs] + Conservatives_i[iEqs];
   }
  }
 }

/*--- Exchange from U_Card to U_i_card ---*/
 for(iEqs = 0; iEqs < nEqs; iEqs++){
  for(iDim = 0; iDim < nADim; iDim++){
   U_i_card[iEqs][iDim] = U_card[iDim][iEqs];
   U_i_card[iEqs][iDim+nADim] = U_card[iDim+nADim][iEqs];
  }
 }

//cout << "U_mean "; for(iEqs=0;iEqs<nEqs;iEqs++)cout<<U_mean[iEqs]<<" ";cout<<endl;
//cout << "U_card "; for(iEqs=0;iEqs<nEqs;iEqs++)cout<<U_card[0][iEqs]<<" ";cout<<endl;

} // end function


void CCentUpw_SW::Compute_adv_coeff(su2double** U){

 for(iDim = 0; iDim < nADim; iDim++){
  adv_coeff[iDim] = max(0.0,max(
   max(U[iDim][iDim+1]/U[iDim][0] + sqrt(g_const*U[iDim][0]/Density_i),U[iDim][iDim+1]/U[iDim][0] - sqrt(g_const*U[iDim][0]/Density_i)),
   max(U[iDim+nADim][iDim+1]/U[iDim+nADim][0] + sqrt(g_const*U[iDim+nADim][0]/Density_i),U[iDim+nADim][iDim+1]/U[iDim+nADim][0] - 
   sqrt(g_const*U[iDim+nADim][0]/Density_i)) ));

  adv_coeff[iDim+nADim] = min(0.0,min(
   min(U[iDim][iDim+1]/U[iDim][0] + sqrt(g_const*U[iDim][0]/Density_i),U[iDim][iDim+1]/U[iDim][0] - sqrt(g_const*U[iDim][0]/Density_i)),
   min(U[iDim+nADim][iDim+1]/U[iDim+nADim][0] + sqrt(g_const*U[iDim+nADim][0]/Density_i),U[iDim+nADim][iDim+1]/U[iDim+nADim][0]-
   sqrt(g_const*U[iDim+nADim][0]/Density_i)) ));
 } 

//cout << "Adv_Coeff ";for(iDim = 0; iDim < 2*nADim; iDim++)cout <<adv_coeff[iDim]<<" ";cout<<endl;
} // end function


void CCentUpw_SW::Compute_U_star(){

/*--- Auxiliary variables --*/
 su2double adv_diff;

 for(iEqs = 0; iEqs < nEqs; iEqs++)
  for(iDim = 0; iDim < nADim; iDim++){
   adv_diff = adv_coeff[iDim] - adv_coeff[iDim + nADim];
   U_star[iEqs][iDim]  = adv_coeff[iDim]*U_i_card[iEqs][iDim+nADim] - adv_coeff[iDim+nADim]*U_i_card[iEqs][iDim];
   U_star[iEqs][iDim] -= (real_flux(U_card[iDim+nADim], iEqs, iDim) - real_flux(U_card[iDim], iEqs, iDim));
   if(adv_diff != 0.0) U_star[iEqs][iDim] /= adv_diff;
  }

//cout << "U_star ";for(iEqs = 0; iEqs < nEqs; iEqs++)cout <<U_star[iEqs][0]<<" ";cout<<endl;
} // end function


void CCentUpw_SW::Compute_delta_U(){

/*--- Auxiliary variables --*/
 su2double aux_U_1, aux_U_2;
 bool positive, same_sign;

 for(iEqs = 0; iEqs < nEqs; iEqs++){
  for(iDim = 0; iDim < nADim; iDim++){
   aux_U_1 = U_i_card[iEqs][iDim+nADim] - U_star[iEqs][iDim];
   aux_U_2 = U_star[iEqs][iDim] - U_i_card[iEqs][iDim];
   positive  =  ( aux_U_1 > 0.0 ) && ( aux_U_2 > 0.0 );
   same_sign = (aux_U_1*aux_U_2 > 0.0);
   if(positive && same_sign) 
    delta_U[iEqs][iDim] = min(aux_U_1,aux_U_2);
   else if(!positive && same_sign)
    delta_U[iEqs][iDim] = max(aux_U_1,aux_U_2);
   else
    delta_U[iEqs][iDim] = 0.0;
  }
 }

//cout << "Delta_U ";for(iEqs = 0; iEqs < nEqs; iEqs++)cout <<delta_U[iEqs][0]<<" ";cout<<endl;
} // end function


void CCentUpw_SW::Compute_numerical_flux(){

/*--- Auxiliary variables --*/
 su2double adv_diff = 0.0;
 su2double adv_prod = 0.0;

/*--- Compute numerical fluxes ---*/
 for(iEqs = 0; iEqs < nEqs; iEqs++){
  for(iDim = 0; iDim < nADim; iDim++){
   adv_diff = adv_coeff[iDim] - adv_coeff[iDim + nADim];
   adv_prod = adv_coeff[iDim]*adv_coeff[iDim + nADim];
   numerical_flux[iEqs][iDim]  = adv_coeff[iDim]*real_flux(U_card[iDim],iEqs,iDim) -
                                 adv_coeff[iDim+nADim]*real_flux(U_card[iDim+nADim],iEqs,iDim);
   numerical_flux[iEqs][iDim] += adv_prod*(U_i_card[iEqs][iDim+nADim] - U_i_card[iEqs][iDim] + delta_U[iEqs][iDim]);
   if(adv_diff != 0.0) numerical_flux[iEqs][iDim] /= adv_diff; else numerical_flux[iEqs][iDim] = 0.0;
 }
}

//cout << "numerical_flux ";for(iEqs = 0; iEqs < nEqs; iEqs++)cout <<numerical_flux[iEqs][0]<<" ";cout<<endl;
} // end function


void CCentUpw_SW::ComputeResidual(su2double *val_residual, CConfig *config){

 AD::StartPreacc(); 
 AD::SetPreaccIn(V_i, nADim+3); AD::SetPreaccIn(V_j, nADim+3); // (p,u,v,rho,delta)
 AD::SetPreaccIn(Normal, nDim);
 if (dynamic_grid) {
  AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
 }

 unsigned short iVar;

 for (iVar = 0; iVar < nVar; iVar++)
   val_residual[iVar] = 0.0;

 /*--- Face area (norm or the normal vector) and unit normal ---*/
 Area = 0.0;
 for (iDim = 0; iDim < nDim; iDim++)
   Area += Normal[iDim]*Normal[iDim];
 Area = sqrt(Area);

 for (iDim = 0; iDim < nDim; iDim++)
   UnitNormal[iDim] = Normal[iDim]/Area;

/*--- Primitive variables at point i ---*/
 Pressure_i  = V_i[0];
 for (iDim = 0; iDim < nADim; iDim++)
   Velocity_i[iDim] = V_i[iDim+1];
 Density_i   = V_i[nADim+1];
 Thickness_i = V_i[nADim+2]; 
 
 real_flux.SetDensity(Density_i);

/*--- Primitive variables at point j ---*/
 Pressure_j  = V_j[0];
 for (iDim = 0; iDim < nADim; iDim++)
  Velocity_j[iDim] = V_j[iDim+1];
 Density_j   = V_j[nADim+1];
 Thickness_j = V_j[nADim+2]; 

/*--- Setting conservative variables ---*/
 Conservatives_i[0] = Density_i*Thickness_i;
 Conservatives_j[0] = Density_j*Thickness_j;
 
 for (iDim = 0; iDim < nADim; iDim++) {
   Conservatives_i[iDim+1] = Density_i*Thickness_i*Velocity_i[iDim];
   Conservatives_j[iDim+1] = Density_j*Thickness_j*Velocity_j[iDim];
 }

//cout << "Conservatives_i "; for(iEqs=0;iEqs<nEqs;iEqs++)cout<<Conservatives_i[iEqs]<<" ";cout<<endl;
//cout << "Conservatives_j "; for(iEqs=0;iEqs<nEqs;iEqs++)cout<<Conservatives_j[iEqs]<<" ";cout<<endl;
/*--- Setting cardinals conservative variables, order important ---*/
 SetCardinals();
 Compute_adv_coeff(U_card);
 Compute_U_star();
 Compute_delta_U();

/*--- Computing numerical flux ---*/
 Compute_numerical_flux();

/*--- Computing pressure under PE hypothesis ---*/
 val_residual[0] = 0.0; 

/*--- Computing residual ---*/
 for(iVar = 1; iVar < nVar; iVar++){   // Oss. nEqs = nVar-1
  for(iDim = 0; iDim < nADim; iDim++){
   val_residual[iVar] = numerical_flux[iVar-1][iDim]*UnitNormal[iDim]*Area;
  }
 }

 AD::SetPreaccOut(val_residual, nVar);
 AD::EndPreacc();

} // end function


/* ----- DEFINITION OF CFILMSOURCE METHODS ----- */

CFilmSource::CFilmSource(unsigned short val_nDim, unsigned short val_nVar, CConfig *config): CNumerics(val_nDim, val_nVar, config){

 unsigned short iOrder, iDim;

 nADim = val_nDim - 1;
 nEqs  = nADim + 1;     // nADim Momentum eq + Continuty eq
 Pev2  = (config->GetFilm_Hp() == PEV2);

 Velocity_i      = new su2double[nADim];
 Conservatives_i = new su2double[nEqs];

 Bottom_Value = new su2double[nADim];
 for(iDim = 0; iDim < nADim; iDim++)
  Bottom_Value[iDim] = config->GetBottom_Value(iDim+1)*PI_NUMBER/180.0;

 dThick_k   = new su2double[nADim];
 dThick_km1 = new su2double[nADim];
 dxk   = new su2double[nADim]; 
 dxkm1 = new su2double[nADim]; 
 
 Coord_i = NULL;
 Coord_j = NULL;

 gravity_force = config->GetGravityForce();
 if(gravity_force){
  body_force_vector = config->GetBody_Force_Vector();
//  for(iDim = 0; iDim < nDim; iDim++)cout<<"Body_Force["<<iDim<<"]= "<<body_force_vector[iDim]<<", ";cout<<endl;
 }

} // end function


CFilmSource::~CFilmSource(){

 delete [] Velocity_i;
 delete [] Conservatives_i;
 delete [] Bottom_Value;
 delete [] dThick_k;
 delete [] dThick_km1;
 delete [] dxk;
 delete [] dxkm1;

} // end function


void CFilmSource::SetdThick(su2double* coord_km1, su2double* coord_k, su2double delta_km1, su2double delta_k){

 unsigned short iDim;
 for(iDim = 0; iDim < nADim; iDim++){
  dThick_k[iDim]   = (coord_k[nADim] - delta_k);
  dThick_km1[iDim] = (coord_km1[nADim] - delta_km1);
  dxk[iDim]   = coord_k[iDim];
  dxkm1[iDim] = coord_km1[iDim];
 }

} // end function


su2double CFilmSource::Hermite(su2double coord, unsigned short order){

 if(order == 0){
  return 1.0;
 } else if(order == 1) {
  return 2.0*coord;
 } else {
  return 2.0*coord*Hermite(coord, order-1) + 2.0*(order-1)*Hermite(coord, order-2);
 }

} // end function


void CFilmSource::SetPressureBC(su2double *coeff, su2double delta_km1, su2double delta_k){

 unsigned short iOrder;

 BC_Pressure_km1 = 0.0;
 BC_Pressure_k   = 0.0;

//
 for(iOrder = 0; iOrder < ProfileOrder; iOrder++){
  BC_Pressure_km1 += coeff[iOrder]*Hermite(delta_km1, iOrder);
  BC_Pressure_k   += coeff[iOrder]*Hermite(delta_k,   iOrder);
 }
//
/*
 BC_Pressure_km1 = Pressure_i*(delta_km1-h_b);
 BC_Pressure_k   = Pressure_i*(delta_k-h_b);
*/
//cout<<"BC_Pressure at d_km1 "<<BC_Pressure_km1<<" and d_k "<<BC_Pressure_k<<endl;
} // end function


void CFilmSource::ComputeResidual(su2double *val_residual, CConfig *config){

 AD::StartPreacc(); 
 AD::SetPreaccIn(V_i, nADim+3); AD::SetPreaccIn(V_j, nADim+3); // (p,u,v,rho,delta)

 Density_i   = V_i[0];
 Thickness_i = V_i[nADim+2]; 

 unsigned short iDim;
 for(iDim = 0; iDim < nADim; iDim++){
  dxk[iDim]   -= Coord_j[iDim];
  dxkm1[iDim] -= Coord_i[iDim];
  if(abs(dxk[iDim])   > epsy) dThick_k[iDim]   /= dxk[iDim]; else dThick_k[iDim]     = 0.0;
  if(abs(dxkm1[iDim]) > epsy) dThick_km1[iDim] /= dxk[iDim]; else dThick_km1[iDim]   = 0.0;
 }

//cout<<"dThick_km1= "<<dThick_km1[0]<<", dThick_k= "<<dThick_k[0]<<endl;
 /*--- Source contribution to the PE hypothesis ---*/
 val_residual[0] = 0.0;

 /*--- Source contribution to the continuity equation ---*/
 val_residual[1] = 0.0;

 /*--- Source contribution to the momentum equations ---*/
 // NB: i segni sono gli opposti al modello perché portato al LHS
 for(iDim = 0; iDim < nADim; iDim++){
  val_residual[iDim+2] = -1.0*Volume*(BC_Pressure_k*dThick_k[iDim] - BC_Pressure_km1*dThick_km1[iDim] );
//cout<<"P contribution=  "<<-1.0*Volume*(BC_Pressure_k*dThick_k[iDim] - BC_Pressure_km1*dThick_km1[iDim] )<<endl;
//cout<<"Mu contribution= "<<-1.0*Volume*(0.5*Laminar_Viscosity_i*BC_Stress_k[iDim] - 0.5*Laminar_Viscosity_i*BC_Stress_km1[iDim])<<endl;
  if(gravity_force){ val_residual[iDim+2] -= Volume*Density_i*body_force_vector[iDim]*Thickness_i; 
//cout<<"g contribution=  "<<-1.0*Volume*body_force_vector[iDim]*Thickness_i<<endl;
  }
 }

 AD::SetPreaccOut(val_residual, nVar);
 AD::EndPreacc();

} // end function


/* ----- DEFINITION OF CFILMVISCOSITY METHODS ----- */

CFilmViscosity::CFilmViscosity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config): CNumerics(val_nDim, val_nVar, config){

 unsigned short iOrder, iDim;

 nADim = val_nDim - 1;
 nEqs  = nADim + 1;     // nADim Momentum eq + Continuty eq
 Pev2  = (config->GetFilm_Hp() == PEV2);

 Velocity_i      = new su2double[nADim];
 Conservatives_i = new su2double[nEqs];

 Bottom_Value = new su2double[nADim];
 for(iDim = 0; iDim < nADim; iDim++)
  Bottom_Value[iDim] = config->GetBottom_Value(iDim+1)*PI_NUMBER/180.0;

 BC_Stress_k   = new su2double[nADim];
 BC_Stress_km1 = new su2double[nADim];
 
 Coord_i = NULL;
 Coord_j = NULL;

} // end function


CFilmViscosity::~CFilmViscosity(){

 delete [] Velocity_i;
 delete [] Conservatives_i;
 delete [] Bottom_Value;
 delete [] BC_Stress_k;
 delete [] BC_Stress_km1;

} // end function


su2double CFilmViscosity::Hermite(su2double coord, unsigned short order){

 if(order == 0){
  return 1.0;
 } else if(order == 1) {
  return 2.0*coord;
 } else {
  return 2.0*coord*Hermite(coord, order-1) + 2.0*(order-1)*Hermite(coord, order-2);
 }

} // end function


void CFilmViscosity::SetStressBC(su2double **coeff, su2double delta_km1, su2double delta_k){

 unsigned short iDim, iOrder;

 for(iDim = 0; iDim < nADim; iDim++){
  BC_Stress_km1[iDim] = 0.0;
  BC_Stress_k[iDim]   = 0.0;
 }

 for(iDim = 0; iDim < nADim; iDim++){
  for(iOrder = 1; iOrder < ProfileOrder; iOrder++){
   //--- Oss on Hermite polynomial: H'_{n}(x) = 2*n*H_{n-1}(x) ---//
   BC_Stress_km1[iDim] += coeff[iDim][iOrder]*2.0*iOrder*Hermite(delta_km1, iOrder - 1);
   BC_Stress_k[iDim]   += coeff[iDim][iOrder]*2.0*iOrder*Hermite(delta_k,   iOrder - 1);
  }
  if(Bottom_Value[iDim] != 0.0){
   BC_Stress_km1[iDim] *= cos(PI_NUMBER - Bottom_Value[iDim]); 
   BC_Stress_k[iDim]   *= cos(PI_NUMBER - Bottom_Value[iDim]);
  }
 }

cout<<"Stress : "<<BC_Stress_k[0]<<" - "<<BC_Stress_km1[0]<<", visc "<<Laminar_Viscosity_i<<endl;
} // end function


void CFilmViscosity::ComputeResidual(su2double *val_residual, CConfig *config){

 AD::StartPreacc(); 
 AD::SetPreaccIn(V_i, nADim+3); AD::SetPreaccIn(V_j, nADim+3); // (p,u,v,rho,delta)

 unsigned short iDim;
 Density_i   = V_i[0];
 Thickness_i = V_i[nADim+2]; 

 /*--- Source contribution to the PE hypothesis ---*/
 val_residual[0] = 0.0;
 if(Pev2) { }
 /*--- Source contribution to the continuity equation ---*/
 val_residual[1] = 0.0;

 /*--- Source contribution to the momentum equations ---*/
 for(iDim = 0; iDim < nADim; iDim++){
  val_residual[iDim+2] = -Volume*(0.5*Laminar_Viscosity_i*BC_Stress_k[iDim] - 0.5*Laminar_Viscosity_i*BC_Stress_km1[iDim]);
//cout<<" = "<<-1.0*Volume*(0.5*Laminar_Viscosity_i*BC_Stress_k[iDim] - 0.5*Laminar_Viscosity_i*BC_Stress_km1[iDim])<<endl;
  if(Pev2) { }
 }

 AD::SetPreaccOut(val_residual, nVar);
 AD::EndPreacc();

} // end function


/* ----- DEFINITION OF CMULTILAYERVISCOSITY METHODS ----- */

CMultilayerViscosity::CMultilayerViscosity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config): CNumerics(val_nDim, val_nVar, config){

 unsigned short iOrder, iDim;

 nADim  = val_nDim - 1;
 nEqs   = nADim + 1;     // nADim Momentum eq + Continuty eq
 nLayer = config->GetnLayer();

 Bottom = false;
 Top    = false;

 Velocity_km1  = new su2double[nADim];
 Velocity_k    = new su2double[nADim];
 Velocity_kp1  = new su2double[nADim];
 
 su2double ModVel_FreeStream, Density_FreeStream, Viscosity_FreeStream, Length_Ref;

 Length_Ref           = config->GetLength_Ref();
 Density_FreeStream   = config->GetDensity_FreeStream();
 Viscosity_FreeStream = config->GetViscosity_FreeStream();
 for(iDim = 0; iDim < nDim; iDim++)
  ModVel_FreeStream += config->GetInc_Velocity_Init()[iDim]*config->GetInc_Velocity_Init()[iDim];
 ModVel_FreeStream   = sqrt(ModVel_FreeStream);                               
 Reynolds_L = Density_FreeStream*ModVel_FreeStream/Viscosity_FreeStream*Length_Ref;
 Friction_Coeff = 0.0592*pow(Reynolds_L, -0.2);
cout<<"1/7 power law friction coefficient = "<<Friction_Coeff<<endl;
 Coord_i = NULL;
 Coord_j = NULL;

} // end function


CMultilayerViscosity::~CMultilayerViscosity(){

 delete [] Velocity_km1;
 delete [] Velocity_k;
 delete [] Velocity_kp1;

} // end function


void CMultilayerViscosity::SetMultilayer_Interaction(CVariable* nodes, unsigned long iPoint, unsigned short iLayer){

 unsigned short iDim;

 if(iLayer == 0)
  Bottom = true;
 else
  Bottom = false;

 if(iLayer == nLayer-1)
  Top = true;
 else
  Top = false;

/*----- Setting Thickness -----*/
  if(!Bottom){
   Thickness_km1 = nodes->GetLayerPrimitive(iPoint, nADim+2, iLayer-1);
  } else {
   Thickness_km1 = 0.0;
  }

  Thickness_k = nodes->GetLayerPrimitive(iPoint, nADim+2, iLayer);
  
  if(!Top) {
   Thickness_kp1 = nodes->GetLayerPrimitive(iPoint, nADim+2, iLayer+1);
  } else {
   Thickness_kp1 = 0.0;
  }

/*----- Setting Velocity -----*/
 for(iDim = 0; iDim < nADim; iDim++){

  if(!Bottom) {
   Velocity_km1[iDim] = nodes->GetLayerPrimitive(iPoint, iDim+1, iLayer-1);
  } else {
   Velocity_km1[iDim] = 0.0;
  }

  Velocity_k[iDim] = nodes->GetLayerPrimitive(iPoint, iDim+1, iLayer);
  
  if(!Top) {
   Velocity_kp1[iDim] = nodes->GetLayerPrimitive(iPoint, iDim+1, iLayer+1);
  } else {
   Velocity_kp1[iDim] = 0.0;
  }

 }

} // end function


void CMultilayerViscosity::ComputeResidual(su2double *val_residual, CConfig *config){

 AD::StartPreacc(); 
 AD::SetPreaccIn(V_i, nADim+3); AD::SetPreaccIn(V_j, nADim+3); // (p,u,v,rho,delta)

 unsigned short iDim;

 val_residual[0] = 0.0;
 /*--- Source contribution to the continuity equation ---*/
 val_residual[1] = 0.0;

 /*--- Source contribution to the momentum equations ---*/
 for(iDim = 0; iDim < nADim; iDim++){

  if(Bottom)
   val_residual[iDim+2]  = Friction_Coeff*Velocity_k[iDim];
  else
   val_residual[iDim+2]  = 2.0*Laminar_Viscosity_i*(Velocity_k[iDim] - Velocity_km1[iDim])/(Thickness_k + Thickness_km1);
  
  if(!Top)
   val_residual[iDim+2] -= 2.0*Laminar_Viscosity_i*(Velocity_kp1[iDim] - Velocity_k[iDim])/(Thickness_kp1 + Thickness_k);

//cout<<"Mu contribution= "<< <<endl;
 }

 AD::SetPreaccOut(val_residual, nVar);
 AD::EndPreacc();

} // end function


