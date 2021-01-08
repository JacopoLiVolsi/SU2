/*!
 * \file CFilmSolver.hpp
 * \brief Headers of the main subroutines for solving thin film equations.
 *        The subroutines and functions are in the <i>CFilmSolver.cpp</i>
 * \author J. Li Volsi
 * \version 7.0.0 "Blackbird"
*/

#pragma once

#include "../solver_structure.hpp"
#include "../variables/CFilmVariable.hpp"
#include "../../../Common/include/geometry_structure.hpp"
#include "../../../Common/include/CMultilayerGeometry.hpp"
#include "../../../Common/include/blas_structure.hpp"

/*! 
*  \class CFilmSolver 
*  \brief Class for defining the solver of thin film problem. 
*  \ingroup Thin_Film_Equations  
*  \author J. Li Volsi
*/
class CFilmSolver : public CSolver {

protected:

 unsigned short nZone;         //!< \brief Number of zones. 
 unsigned short nADim;         //!< \brief Number of Dimension of the Average problem. 
 unsigned long  nEdge;
 unsigned short nPrimVarFlow;
 CFilmVariable* nodes = nullptr;
 unsigned short nLayer;        //!< \brief Total number of layers.
 bool multilayer;              //!< \brief Boolean for multilayer check.
 bool viscous;                 //!< \brief Enable PEV2 hypothesis.
 bool first_bc;
 unsigned long Inlet_Layer_Count  = 0;
 unsigned long Inlet_Vertex_Count = 0;
 su2double g_const = -STANDARD_GRAVITY; // [m/s^2]
 su2double layer_toll;         //!< \brief Tollerance value for the layer thickness (1% of initial thickness).
 su2double *body_force_vector;  
 su2double *delta_layer;       //!< \brief Uniform distance between two layers at each point.
 su2double **total_thickness;  //!< \brief Stores the cumulative thickness for each layer.
 su2double Density_Inf, Pressure_Inf, Pressure_Init;
 su2double Viscosity_Inf, Thickness_Inf;
 su2double *Velocity_Inf, *Velocity_Init;
 CFluidModel  *FluidModel;     //!< \brief Pointer to a modellization of the fluid, useful for density.
 su2double *Primitive;

/*--- Residual variables for each layer ---*/
 su2double **Res_Conv_Layer;         //!< \brief Convective residual for each layer.
 CSysVector<su2double> Layer_SysRes; //!< \brief System of residual for all layers  

/*--- Multilayer geometry variables ---*/
 CMultilayerEdge   ***layer_edge;
 CMultilayerPoint  ***layer_node;
 CMultilayerVertex ****layer_vertex;
 su2double **Bottom_Topography;
 su2double  *Bottom_Value;

/*--- Boundary variables ---*/
 su2double ***CharacPrimVar;
 su2double ****CharacPrimVarLayer;
 su2double ***Inlet_Ptotal;    
 su2double ***Inlet_Ttotal;    
 su2double ****Inlet_FlowDir;

/*--- Sliding meshes variables ---*/
 su2double ****SlidingState;
 int **SlidingStateNodes;

/*--- Profile reconstruction variables ---*/
 int ProfileOrder;
 int RowOrder;

public:

/*! 
 *\brief Void constructor of the class. 
*/
 CFilmSolver();

/*! 
 *\brief Basic constructor of the class. 
*/
 CFilmSolver(CGeometry *geometry, CConfig *config, unsigned short iMGlevel);

/*!
 * \brief Return nodes to allow CSolver::base_nodes to be set.
*/
 inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

 inline CFluidModel* GetFluidModel() override { return FluidModel; }

/*! 
 *\brief Check non-dimensionalization of the problem. 
*/
 void SetNondimensionalization(CConfig *config, unsigned short iMesh);

/*! 
 *\brief Construct the bottom topography. 
*/
 void ReadTopography(su2double **topography, CGeometry* geometry, CConfig* config) override;

 void ComputeProfile(CConfig *config, su2double* Coeff, unsigned long iPoint, unsigned short ivar);

/*!
 * \brief Destructor of the class.
*/
 ~CFilmSolver(void);

/*!
 * \brief Set initial conditions.
*/
 void SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter);

/*!
 * \brief Set the old solution variables to the current solution value for Runge-Kutta iteration.
 * \param[in] geometry - Geometrical definition of the problem.
 */
 void Set_OldSolution(CGeometry *geometry);

/*!
 * \brief Set the new solution variables to the current solution value for classical RK.
 * \param[in] geometry - Geometrical definition of the problem.
 */
 void Set_NewSolution(CGeometry *geometry);

/*!
 * \brief Preprocessing solver container.
*/
 void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);

/*!
 * \brief Update all the geometry structures.
*/
 void Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh);

/*!
 * \brief Set primitive variables.
*/
 unsigned long SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output);

 su2double* GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex);

 su2double* GetCharacPrimVar(unsigned short iLayer, unsigned short val_marker, unsigned long val_vertex);

 int GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex);

 su2double GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, unsigned long donor_index); 

 void SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex);

 su2double GetInterfaceVar(unsigned short val_marker, unsigned long iPoint, unsigned short iVar) override;

 su2double GetTotal_thickness(unsigned short iLayer, unsigned long iPoint) override;

 inline CMultilayerPoint*** GetLayerNode() override { return layer_node; }

 inline CMultilayerPoint* GetLayerNode(unsigned short iLayer, unsigned long iPoint) override{ 
         if(iLayer < (nLayer-1) ) 
          return layer_node[iLayer][iPoint];
         else return NULL;
         }  

/*!
 * \brief Recursive definition of Hermite polynomial.
*/
 su2double Hermite(su2double coord, unsigned short order);

/*!
 * \brief Integrate between a and b Hermite polynomial.
*/
 su2double IntegrateHermite(su2double a, su2double b, unsigned short order);
 
/*!
 * \brief Set time step.
*/
 void SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration);    

/* ----- METHODS FOR SPACE INTERGATION -----*/

 void Centered_Residual (CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh, unsigned short iRKstep);
 void Convective_Residual (CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh, unsigned short iRKstep);
/*!
 * \brief Compute the spatial integration using a upwind scheme.
*/
 void Upwind_Residual (CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) override;

/*!
 * \brief Compute viscous contribution.
*/
 void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
                       CConfig *config, unsigned short iMesh, unsigned short iRKStep) override;

/*!
 * \brief Compute source contribution.
*/
 void Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, 
                      CConfig *config, unsigned short iMesh) override;

 void ComputeResidual_Multizone(CGeometry *geometry, CConfig *config) override;

/* ----- METHODS FOR TIME INTERGATION -----*/

 void ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep); 	
 void ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config); 	
/*! 
 *\brief Explicit Euler time integration.
*/
 void ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config); 
/*! 
 *\brief Classical 4th order Runge-Kutta time integration.
*/	
 void ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep); 	

/* ----- METHODS FOR BOUNDARY CONDITIONS -----*/

 void BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, 
                         CNumerics *visc_numerics, CConfig *config) override;

 void BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
                    CConfig *config, unsigned short val_marker) override;

 void BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                   CConfig *config, unsigned short val_marker) override;
/*!
 * \brief Set uniform inlet boundary conditions.
*/
 void SetUniformInlet(CConfig* config, unsigned short iMarker);

/*!
 * \brief Set a specific profile for inlet boundary conditions.
*/
 void SetInletAtVertex(su2double *val_inlet, unsigned short iMarker, unsigned long iVertex);

 void BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
              CConfig *config, unsigned short val_marker) override;

 void BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                           CConfig *config, unsigned short val_marker) override;

 void BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
                CConfig *config, unsigned short val_marker);

}; // end class CFilmSolver





