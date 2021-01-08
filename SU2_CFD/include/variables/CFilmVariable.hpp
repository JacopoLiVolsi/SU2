/*!
 * \file CFilmVariable.hpp
 * \brief Headers of generic thin film variable
 * \author J. Li Volsi
 * \version 7.0.0 "Blackbird"
*/

#pragma once

#include "CVariable.hpp"

/*! 
*  \class CFilmVariable 
*  \brief Class for defining variable of thin film problem. 
*  \ingroup Thin_Film_Equations  
*  \author J. Li Volsi
*/
class CFilmVariable: public CVariable{

protected:
 unsigned long nPrimGrad;   //!< \brief Dimension of the gradient of primitive variables. 
 unsigned long nADim;       //!< \brief Number of Dimension of the Average problem. 
 bool bottom, kth, top;     //!< \brief Bool to indicate particular layer.
 unsigned short LayerLabel; //!< \brief Index of the current layer, from 1 to nLayer.
 su2double layer_toll;      //!< \brief Tollerance value for the layer thickness (1% of initial thickness).
 unsigned short nLayer;     //!< \brief Index of max number of layers.
 MatrixType Primitive;      //!< \brief Class containing of primitive variables (p_k,u_k,v_k,rho,delta_(k-1,k)).
 MatrixType Solution_New;   //!< \brief Container for the new solution.
 MatrixType Bottom_BC;      //!< \brief Container of the boundary condition of the bottom topography.
 VectorType BottomTopography; //!< \brief Container of the bottom topography.
 VectorType Viscosity;        //!< \brief Container of fluid viscosity.
 // CVariable::Solution contains (p_k, h_(k-1,k), q_x, q_y)
 
public:


/*! 
 *\brief Constructor of the class 
*/
 CFilmVariable(su2double density, su2double thickness,const su2double* velocity, su2double val_layer_toll, unsigned long npoint, 
               unsigned long ndim, unsigned long nvar, CConfig *config, unsigned short val_iLayer);
/*! 
 *\brief Constructor of the class 
*/
 CFilmVariable(su2double density, su2double thickness,const su2double* velocity, unsigned long npoint, 
               unsigned long ndim, unsigned long nvar, CConfig *config, unsigned short val_iLayer);

/*! 
  *\brief Basic destructor of the class 
 */
 ~CFilmVariable();

 inline const bool isBottom() { return bottom;};
 inline const bool isKth()    { return kth;};
 inline const bool isTop()    { return top;};

/* --- Setter methods ---*/
 virtual bool SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel);
 inline void SetPrimitive(unsigned long iPoint, const su2double *val_prim) override{
      for (unsigned long iVar = 0; iVar < nPrimVar; iVar++)
       Primitive(iPoint,iVar) = val_prim[iVar];}
 inline void SetPrimitive(unsigned long iPoint, unsigned long iVar, su2double val_prim) override{
       Primitive(iPoint,iVar) = val_prim;}
 inline void SetBottom_BC(unsigned long iPoint, const su2double *val_bc){
      for (unsigned long iVar = 0; iVar < nPrimVar; iVar++)
       Bottom_BC(iPoint,iVar) = val_bc[iVar];}
 inline void SetBottom_BC(unsigned long iPoint, unsigned short iVar, const su2double val_bc){
       Bottom_BC(iPoint,iVar) = val_bc;}
 void ImportBottom_Topography(su2double** topography);
 void InitializeBottom_BC(CConfig *config) override;
 inline void SetSolution(unsigned long iPoint, const su2double *val_sol){
      for (unsigned long iVar = 0; iVar < nVar; iVar++)
       Solution(iPoint,iVar) = val_sol[iVar];}
 inline void SetSolution(unsigned long iPoint, unsigned short iVar, su2double val_sol){
       Solution(iPoint,iVar) = val_sol;}
 inline void SetSolution_New(unsigned long iPoint, unsigned short iVar, su2double val_sol){
       Solution_New(iPoint,iVar) = val_sol;}
 inline void SetPressure(unsigned long iPoint) {
     Primitive(iPoint,0) = Solution(iPoint, 0);}
 inline bool SetDensity(unsigned long iPoint,su2double density){
     Primitive(iPoint,nADim+1) = density;
     return density <= 0.0;}
 inline void SetLaminarViscosity(unsigned long iPoint, su2double laminarViscosity) override {
    Viscosity(iPoint) = laminarViscosity;
  }
 inline void SetVelocity(unsigned long iPoint){
    for (unsigned long iDim = 0; iDim < nADim; iDim++)
      Primitive(iPoint,iDim+1) = Solution(iPoint,iDim+2) / Solution(iPoint,1);
  }
 inline bool SetThickness(unsigned long iPoint, su2double thickness){
      if((thickness < layer_toll) || (thickness <= 0.0))
       Primitive(iPoint,nADim+2) = 0.0;
      else
       Primitive(iPoint,nADim+2) = thickness;
      return (thickness < layer_toll) || (thickness <= 0.0);
 }
 inline void SetLayerToll(su2double val_toll){ layer_toll = val_toll;}
 void ReconstructVelocity(unsigned long iPoint);
 void AddSolution_New(unsigned long iPoint, unsigned long iVar, su2double solution) final;
 virtual void SetLayerSolution_Old(unsigned long iPoint, su2double* solution, unsigned short iLayer);
 virtual void SetLayerPrimitive(unsigned long iPoint, unsigned short iVar, su2double prim_val, unsigned short iLayer);
 virtual void AddLayerSolution(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer);
 virtual void AddLayerSolution_New(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer);
 virtual void SetLayerSolution(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer);
 virtual void SetLayerSolution_New(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer);

/* --- Getter methods ---*/
 inline const unsigned short GetnLayer() { return nLayer;};
 inline  su2double* GetPrimitive(unsigned long iPoint) override {return Primitive[iPoint];}
 inline  su2double* GetBottom_BC(unsigned long iPoint) {return Bottom_BC[iPoint];}
 inline  su2double  GetPrimitive(unsigned long iPoint, unsigned long iVar) const override {return Primitive(iPoint,iVar);}
 inline  su2double  GetBottom_BC(unsigned long iPoint, unsigned short iVar) {return Bottom_BC(iPoint,iVar);}
 virtual su2double* GetLayerPrimitive(unsigned long iPoint, unsigned short iLayer);
 virtual su2double  GetLayerPrimitive(unsigned long iPoint, unsigned short iVar, unsigned short iLayer);
 inline  su2double* GetSolution(unsigned long iPoint) {return Solution[iPoint];}
 inline virtual su2double* GetLayerSolution(unsigned long iPoint, unsigned short iLayer) {
                           if(iLayer == (nLayer - 1))
                            return Solution[iPoint];
                           else return NULL;}
 inline su2double  GetSolution_New(unsigned long iPoint, unsigned long iVar) const final { return Solution_New(iPoint,iVar); }
 virtual su2double GetLayerSolution_New(unsigned long iPoint, unsigned long iVar, unsigned short iLayer ) const;
 inline su2double  GetPressure(unsigned long iPoint) const final { return Primitive(iPoint,0); }
 inline su2double  GetDensity(unsigned long iPoint) const final { return Primitive(iPoint,nADim+1);}
 inline su2double  GetLaminarViscosity(unsigned long iPoint) const override { return Viscosity(iPoint);}
 inline su2double* GetVelocity(unsigned long iPoint) const { 
        su2double* vel = new su2double[nADim];
        for(unsigned long iDim = 0; iDim < nADim; iDim++)
         vel[iDim] = Primitive(iPoint,iDim+1);
        return vel;}
 inline su2double GetProjVel(unsigned long iPoint, const su2double *val_vector) const final {
        su2double ProjVel = 0.0;
        for (unsigned long iDim = 0; iDim < nADim; iDim++)
         ProjVel += Primitive(iPoint,iDim+1)*val_vector[iDim];
        return ProjVel;}
 inline virtual CFilmVariable* GetLayer(unsigned short iLayer) { return this; } 

};




/*!
*  \class CMultilayerFilmVariable 
*  \brief Class collecting thin film variable for all layers.
*  The idea is that the inherited members stores the values for the top
*  layer, i.e. the one communicating with the other fluid, while for the other
*  layers the base class  is used.
*  \ingroup Thin_Film_Equations  
*  \author J. Li Volsi
*/
class CMultilayerFilmVariable final: public CFilmVariable{

private:
 CFilmVariable **KthLayer;

public:
/*! 
 *\brief Basic constructor of the class 
*/
 CMultilayerFilmVariable(su2double* density, su2double* thickness,const su2double* velocity, su2double val_layer_toll, 
                         unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);
/*! 
 *\brief Basic destructor of the class 
*/
 ~CMultilayerFilmVariable();

 su2double* GetLayerPrimitive(unsigned long iPoint, unsigned short iLayer) final;
 su2double  GetLayerPrimitive(unsigned long iPoint, unsigned short iVar, unsigned short iLayer) final;
 su2double* GetLayerSolution(unsigned long iPoint, unsigned short iLayer) final;
 su2double  GetLayerSolution_New(unsigned long iPoint, unsigned long iVar, unsigned short iLayer ) const final;

 void SetSolution_New() final;
 bool SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) final;
 void SetLayerSolution_Old(unsigned long iPoint, su2double* solution, unsigned short iLayer) final;
 void SetLayerPrimitive(unsigned long iPoint, unsigned short iVar, su2double prim_val, unsigned short iLayer) final;
 void AddLayerSolution(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer) final;
 void AddLayerSolution_New(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer) final;
 void SetLayerSolution(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer) final;
 void SetLayerSolution_New(unsigned long iPoint, unsigned long iVar, su2double solution, unsigned short iLayer) final;

 CFilmVariable* GetLayer(unsigned short iLayer) final; 

 
};

































