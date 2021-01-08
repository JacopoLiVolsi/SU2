/*!
 * \file CNumerics_Film.hpp
 * \brief Headers for the numerical scheme for thin film eq. 
 *       The subroutines and functions are in the <i>CNumerics_Film.cpp</i>
 * \author J. Li Volsi
 * \version 7.0.0 "Blackbird"
*/

#pragma once

#include "../numerics_structure.hpp"

/*!
*  \struct SW_Flux
*  \brief Functor for 2-D flux of Shallow Water eq 
*/
struct SW_Flux{
su2double density = 0.0;
/*! \brief Constructor of the struct*/
SW_Flux() = default;

inline void SetDensity(su2double _density){ density = _density;}
inline su2double GetDensity(void) const { return density;}

/*! \brief Overloaded access operator*/
su2double operator() (su2double* U, unsigned short iEqs, unsigned short iDim);
}; // end struct SW_Flux


/*!
*  \class CCentUpw_SW
*  \brief Class implementing Kurganov central-upwind method for the convective term of Shallow Water eq.
*  \ingroup Thin_Film_Equations  
*  \author J. Li Volsi
*/
class CCentUpw_SW : public CNumerics{

protected:

 unsigned long  nADim;    //!< \brief Number of Dimension of the Average problem. 
 unsigned short nEqs;     //!< \brief Number of eqs solved with Kurganov central upwind.
 unsigned short iEqs,iDim;
 bool dynamic_grid, gravity_force;
 const su2double epsy = 1.0E-10;
 su2double  g_const;                  // [m/s^2]
 su2double* body_force_vector;        // [m/s^2]
 su2double* Delta;                    // [m]
 su2double  Thickness_i, Thickness_j; // [m]
 su2double* Velocity_i, *Velocity_j;  // [m/s]
 su2double* Conservatives_i, *Conservatives_j;

// U is the conserved vector of mean variables, i.e. U = (rho*delta, rho*u*delta, rho*v*delta)
// The first * is for the component of U, the second for the direction
 su2double*  U_mean;          // [ U1, U2, U3] <-- Hp nEqs == 3
 su2double** U_card;          // [ U1_E, U2_E, U3_E][ U1_N, U2_N, ...][U_W][U_S]
 su2double** U_i_card;        // [iEqs][ U_E, U_N, U_W, U_S ]
 su2double** U_star;          // [iEqs][ U*_(j+1/2,i), U*_(j,i+1/2) ]
 su2double** delta_U;         // [iEqs][ dU_(j+1/2,i), dU_(j,i+1/2) ]
 su2double*  adv_coeff;       // [ a+_(j+1/2,i), b+_(j,i+1/2), a-_(j+1/2,i), b-_(j,i+1/2)]
 su2double** numerical_flux;  // [iEqs][ F_(j+1/2,i), G_(j,i+1/2) ] 

 SW_Flux real_flux;           //!< \brief Functor for real flux


public:

 CCentUpw_SW(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

 ~CCentUpw_SW(void);

 void SetCardinals();

 void Compute_adv_coeff(su2double** U);

 void Compute_U_star();
 
 void Compute_delta_U();

 void Compute_numerical_flux();
 
 void ComputeResidual(su2double *val_residual, CConfig *config) override; 

}; // end class CCentUpw_SW



/*!
*  \class CFilmSource
*  \brief Class implementing multilayer source interaction term.
*  \ingroup Thin_Film_Equations  
*  \author J. Li Volsi
*/
class CFilmSource: public CNumerics{

protected:

 bool Pev2;
 unsigned long  nADim;        //!< \brief Number of Dimension of the Average problem. 
 unsigned short nEqs;         //!< \brief Number of eqs solved with Kurganov central upwind.
 unsigned short ProfileOrder; //!< \brief Nr of point needed for profile reconstruction in source term.

 bool  gravity_force;
 const su2double epsy = 1.0E-6;
 su2double  Thickness_i;
 su2double *Velocity_i;
 su2double *Conservatives_i;
 su2double *dThick_k, *dThick_km1;
 su2double *dxk, *dxkm1;
 su2double *body_force_vector; 
 su2double *Bottom_Value;

 su2double  BC_Pressure_k,  BC_Pressure_km1;

public:

 CFilmSource(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

 ~CFilmSource();

 inline unsigned short GetProfileOrder() {return ProfileOrder;}
 inline void SetProfileOrder(unsigned short ord) override { ProfileOrder = ord;}

 su2double Hermite(su2double coord, unsigned short order);

 void SetdThick(su2double* coord_km1, su2double* coord_k, su2double delta_km1, su2double delta_k) override;
 
 void SetPressureBC(su2double* coeff, su2double delta_km1, su2double delta_k);
 
 void ComputeResidual(su2double *val_residual, CConfig *config) override; 

}; // end class 


/*!
*  \class CFilmViscosity
*  \brief Class implementing PE/PEV2 viscosity contribution.
*  \ingroup Thin_Film_Equations  
*  \author J. Li Volsi
*/
class CFilmViscosity: public CNumerics{

protected:

 bool Pev2;
 unsigned long  nADim;        //!< \brief Number of Dimension of the Average problem. 
 unsigned short nEqs;         //!< \brief Number of eqs solved with Kurganov central upwind.
 unsigned short ProfileOrder; //!< \brief Nr of point needed for profile reconstruction in source term.

 const su2double epsy = 1.0E-6;
 su2double  Thickness_i;
 su2double *Velocity_i;
 su2double *Conservatives_i;
 su2double *Bottom_Value;
 su2double* BC_Stress_k, *BC_Stress_km1;

public:

 CFilmViscosity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

 ~CFilmViscosity();

 inline unsigned short GetProfileOrder() {return ProfileOrder;}
 inline void SetProfileOrder(unsigned short ord) override { ProfileOrder = ord;}

/*!
 * \brief Recursive definition of Hermite polynomial.
*/
 su2double Hermite(su2double coord, unsigned short order);

 void SetStressBC(su2double**  coeff, su2double delta_km1, su2double delta_k);
 
 void ComputeResidual(su2double *val_residual, CConfig *config) override; 

}; // end class 



/*!
*  \class CMultiLayerViscosity
*  \brief Class implementing Audesse PE viscosity contribution.
*  \ingroup Thin_Film_Equations  
*  \author J. Li Volsi
*/
class CMultilayerViscosity: public CNumerics{

protected:

 unsigned long  nADim;        //!< \brief Number of Dimension of the Average problem. 
 unsigned short nLayer;       //!< \brief Number of layers.
 unsigned short nEqs;         //!< \brief Number of eqs solved with Kurganov central upwind.
 unsigned short ProfileOrder; //!< \brief Nr of point needed for profile reconstruction in source term.
 bool Bottom;
 bool Top;

 const su2double epsy = 1.0E-6;
 su2double  Reynolds_L;
 su2double  Friction_Coeff;
 su2double  Thickness_km1, Thickness_k, Thickness_kp1;
 su2double *Velocity_km1, *Velocity_k, *Velocity_kp1;

public:

 CMultilayerViscosity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

 ~CMultilayerViscosity();

 void SetMultilayer_Interaction(CVariable* nodes, unsigned long iPoint, unsigned short iLayer) override;
 
 void ComputeResidual(su2double *val_residual, CConfig *config) override; 

}; // end class 




















