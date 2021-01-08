/*!
 * \file CFilmSolver.cpp
 * \brief Source of the main subroutines for solving thin film equations.
 * \author J. Li Volsi
 * \version 7.0.0 "Blackbird"
*/

#include "../../include/solvers/CFilmSolver.hpp"


CFilmSolver::CFilmSolver(): CSolver(){

/*--- Basic initialization to NULL --- */
 Res_Conv_Layer     = NULL;
 FluidModel         = NULL;
 body_force_vector  = NULL;
 delta_layer        = NULL;
 total_thickness    = NULL;
 layer_edge         = NULL;
 layer_node         = NULL;
 Bottom_Topography  = NULL;
 Bottom_Value       = NULL;
 Velocity_Init      = NULL;
 SlidingState       = NULL;
 SlidingStateNodes  = NULL;
 CharacPrimVar      = NULL;
 CharacPrimVarLayer = NULL;
 Inlet_Ptotal       = NULL;
 Inlet_Ttotal       = NULL;
 Inlet_FlowDir      = NULL;

} // end function


CFilmSolver::CFilmSolver(CGeometry *geometry, CConfig *config, unsigned short iMGlevel) : CSolver(){ // oppure -> : CSolver(true)

 unsigned short iVar, iLayer, iDim, iMarker;
 unsigned long  iEdge, iVertex, iPoint;

/*--- Basic initialization to NULL --- */
 Res_Conv_Layer     = NULL;
 FluidModel         = NULL;
 body_force_vector  = NULL;
 delta_layer        = NULL;
 total_thickness    = NULL;
 layer_edge         = NULL;
 layer_node         = NULL;
 Bottom_Topography  = NULL;
 Bottom_Value       = NULL;
 Velocity_Init      = NULL;
 SlidingState       = NULL;
 SlidingStateNodes  = NULL;
 CharacPrimVar      = NULL;
 CharacPrimVarLayer = NULL;
 Inlet_Ptotal       = NULL;
 Inlet_Ttotal       = NULL;
 Inlet_FlowDir      = NULL;
 multilayer         = false; 
 viscous            = false;
 first_bc           = true;

/*--- Define geometry constants in the solver structure ---*/
 nZone = config->GetnZone();   
 nDim  = geometry->GetnDim();
 nADim = nDim  - 1;
 nVar  = nADim + 2; // Unknowns (p, h, q_x, q_y) -- h=rho*delta -- q_x=u*h -- q_y=v*h 

 nPrimVar      = nADim+3; // Primitive variables (p, vx, vy, rho, delta) 
 nPrimVarGrad  = nADim;
 nPrimVarFlow  = nDim+9;
 nSecondaryVar = 0; nSecondaryVarGrad = 0;
 nMarker       = config->GetnMarker_All();
 
 nPoint       = geometry->GetnPoint();
 nPointDomain = geometry->GetnPointDomain();
 nEdge        = geometry->GetnEdge();
 nLayer       = config->GetnLayer();
 if( (config->GetKind_Film_Solver() == MULTI_LAYER_ASYMP) && (nLayer > 1) )
  multilayer = true;
 else if( (config->GetKind_Film_Solver() == MULTI_LAYER_ASYMP) && (nLayer == 1) ){
  multilayer = true;
  nLayer = 3;
 }

 cout << "Solving THIN_FILM problem using ";
 if(multilayer)
  cout << "multilayer asymptotic method with " << nLayer << " layers.";
 else
  cout << "a generic method.";
 cout << endl;
   
/* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
 dynamic_grid = config->GetDynamic_Grid();
 viscous = ( config->GetFilm_Hp() == PEV2 );

 if(config->GetGravityForce()){
  body_force_vector = config->GetBody_Force_Vector();
 }

/*--- Store the number of vertices on each marker for deallocation later ---*/ 
 nVertex = new unsigned long[nMarker];
 for (iMarker = 0; iMarker < nMarker; iMarker++) 
  nVertex[iMarker] = geometry->nVertex[iMarker];

/*--- Adjusting the normal for 2D problems ---*/
 if(nDim == 2){
  su2double *n2d = new su2double[nDim];
  unsigned long jPoint;
  for(iEdge = 0; iEdge < nEdge; iEdge++){
   iPoint = geometry->edge[iEdge]->GetNode(0);
   jPoint = geometry->edge[iEdge]->GetNode(1);
   n2d[0] = geometry->node[iPoint]->GetCoord(0) - geometry->node[jPoint]->GetCoord(0);
   n2d[1] = geometry->node[iPoint]->GetCoord(1) - geometry->node[jPoint]->GetCoord(1);
   geometry->edge[iEdge]->SetNormal(n2d);
  }
  for(iMarker = 0; iMarker < nMarker; iMarker++) {
    for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex++ ) {
     iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
     iEdge  = geometry->node[iPoint]->GetEdge(0);
     jPoint = (geometry->edge[iEdge]->GetNode(1) == iPoint) ? geometry->edge[iEdge]->GetNode(0) : geometry->edge[iEdge]->GetNode(1);
     n2d[0] = geometry->node[iPoint]->GetCoord(0) - geometry->node[jPoint]->GetCoord(0);
     n2d[1] = geometry->node[iPoint]->GetCoord(1) - geometry->node[jPoint]->GetCoord(1);
     geometry->vertex[iMarker][iVertex]->SetNormal(n2d);
    }
  }
 delete [] n2d;
 }

 SetNondimensionalization(config, MESH_0);

/*--- Multilayer variables ---*/
 delta_layer       = new su2double [nPoint];
 total_thickness   = new su2double*[nLayer];
 Bottom_Topography = new su2double*[nPoint];
 Bottom_Value      = new su2double [3];
 Bottom_Value[0]   = config->GetBottom_Value(0);
 Bottom_Value[1]   = config->GetBottom_Value(1)*PI_NUMBER/180.0;
 Bottom_Value[2]   = config->GetBottom_Value(2)*PI_NUMBER/180.0;
 for(iPoint = 0; iPoint < nPoint; iPoint++)
  Bottom_Topography[iPoint] = new su2double[nDim];
 ReadTopography(Bottom_Topography, geometry, config);
 /* Initialize tollerance */  
 su2double Delta_Tot = 0.0;
 for(iDim = 0; iDim < nDim; iDim++)
  Delta_Tot += pow((geometry->node[0]->GetCoord(iDim) - Bottom_Topography[0][iDim]), 2);
 Delta_Tot = sqrt(Delta_Tot);
 layer_toll = 0.01*Delta_Tot/su2double(nLayer); 

 for(iLayer = 0; iLayer < nLayer; iLayer++ ){
  total_thickness[iLayer] = new su2double[nPoint];
  for(iPoint = 0; iPoint < nPoint; iPoint++){
   Delta_Tot = 0.0;
   for(iDim = 0; iDim < nDim; iDim++)
    Delta_Tot += pow((geometry->node[iPoint]->GetCoord(iDim) - Bottom_Topography[iPoint][iDim]), 2);
   Delta_Tot = sqrt(Delta_Tot);
   delta_layer[iPoint] = Delta_Tot/su2double(nLayer);
   layer_toll = min(layer_toll, 0.01*delta_layer[iPoint]);
   total_thickness[iLayer][iPoint] = (iLayer+1)*delta_layer[iPoint];
  } 
 }

 /*--- Check if we are executing a verification case. If so, the
  VerificationSolution object will be instantiated for a particular
  option from the available library of verification solutions. Note
  that this is done after SetNondim(), as problem-specific initial
  parameters are needed by the solution constructors. ---*/
 SetVerificationSolution(nDim, nVar, config);

/*--- Setting layer geometry. ---*/
 if(multilayer) {
  su2double* aux_norm = new su2double [nDim];
  layer_edge   = new CMultilayerEdge**[nLayer-1];
  layer_node   = new CMultilayerPoint**[nLayer-1];
  layer_vertex = new CMultilayerVertex***[nLayer-1];
  
  for(iLayer = 0; iLayer < (nLayer-1); iLayer++){
/*--- Setting multilayer points. ---*/
   layer_node[iLayer] = new CMultilayerPoint*[nPoint];
   for(iPoint = 0; iPoint < nPoint; iPoint++){
    su2double dx, dy, dz;
    dx = (geometry->node[iPoint]->GetCoord(0) - Bottom_Topography[iPoint][0])/su2double(nLayer);
    dy = (geometry->node[iPoint]->GetCoord(1) - Bottom_Topography[iPoint][1])/su2double(nLayer);
    if(nDim == 2){
     layer_node[iLayer][iPoint] = new CMultilayerPoint(geometry->node[iPoint]->GetCoord(0) - (nLayer-1-iLayer)*dx, 
                                                       geometry->node[iPoint]->GetCoord(1) - (nLayer-1-iLayer)*dy,iPoint);
    }
    if(nDim == 3){
     dz = (geometry->node[iPoint]->GetCoord(2) - Bottom_Topography[iPoint][2])/su2double(nLayer);
     layer_node[iLayer][iPoint] = new CMultilayerPoint(geometry->node[iPoint]->GetCoord(0) - (nLayer-1-iLayer)*dx, 
           geometry->node[iPoint]->GetCoord(1) - (nLayer-1-iLayer)*dy, geometry->node[iPoint]->GetCoord(2) - (nLayer-1-iLayer)*dz,iPoint);
    }
    layer_node[iLayer][iPoint]->SetVolume(geometry->node[iPoint]->GetVolume());
   }

/*--- Setting multilayer edges. ---*/
   layer_edge[iLayer] = new CMultilayerEdge*[nEdge];
   for(iEdge = 0; iEdge < nEdge; iEdge++){
    layer_edge[iLayer][iEdge] = new CMultilayerEdge(geometry->edge[iEdge]->GetNode(0),geometry->edge[iEdge]->GetNode(1),nDim);
    geometry->edge[iEdge]->GetNormal(aux_norm);
    layer_edge[iLayer][iEdge]->SetNormal(aux_norm);
   }

/*--- Setting multilayer vertices. ---*/
  layer_vertex[iLayer] = new CMultilayerVertex**[nMarker];
   for(iMarker = 0; iMarker < nMarker; iMarker++){
    layer_vertex[iLayer][iMarker] = NULL;
    if(config->GetMarker_All_KindBC(iMarker) != FLUID_INTERFACE){
     layer_vertex[iLayer][iMarker] = new CMultilayerVertex*[nVertex[iMarker]];
     for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex++ ){
      layer_vertex[iLayer][iMarker][iVertex] = new CMultilayerVertex(geometry->vertex[iMarker][iVertex]->GetNode(), nDim);
      geometry->vertex[iMarker][iVertex]->GetNormal(aux_norm);
      layer_vertex[iLayer][iMarker][iVertex]->SetNormal(aux_norm);
     }
    }
   } // end for iMarker

  } // end for iLayer
  delete [] aux_norm;
 } // end if multilayer

/*--- Res_Conv is for the top layer, Res_Conv_Layer[k] for the kth layer. ---*/
 Residual      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
 Res_Conv      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Conv[iVar]      = 0.0;
 Res_Visc      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Visc[iVar]      = 0.0;
 Res_Sour      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Res_Sour[iVar]      = 0.0; 
 Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
 Residual_BGS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar]  = 0.0;
 Residual_Max_BGS = new su2double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar] = 0.0;
 Point_Max_BGS = new unsigned long[nVar];     for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar] = 0;

 //LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
 unsigned long nlayernpoint    = nLayer*nPoint;
 unsigned long nlayernpointdom = nLayer*nPointDomain;

 Layer_SysRes.Initialize(nlayernpoint, nlayernpointdom, nVar, 0.0);
// LinSysRes.Initialize(nlayernpoint, nlayernpointdom, nVar, 0.0);

 if(multilayer){
 Res_Conv_Layer = new su2double*[nLayer-1];
  for(iLayer = 0; iLayer < (nLayer-1); iLayer++){
   Res_Conv_Layer[iLayer] = new su2double[nVar];
   for(iVar = 0; iVar < nVar; iVar++)
    Res_Conv_Layer[iLayer][iVar] = 0.0;
  }
 }

 // Maybe useless
 Primitive = new su2double[nPrimVar]; 
 for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;

/*--- Read non dimensional initial conditions ---*/
 Velocity_Init = new su2double[nDim];
 Pressure_Init = config->GetPressure_FreeStreamND();
 Velocity_Init = config->GetInc_Velocity_Init();
 for(iDim = 0; iDim < nDim; iDim++){
  Velocity_Init[iDim] /= config->GetVelocity_Ref();
 }
/*--- Read farfield conditions from config ---*/
 Density_Inf     = config->GetDensity_FreeStreamND();
 Pressure_Inf    = config->GetPressure_FreeStreamND();
 Velocity_Inf    = config->GetVelocity_FreeStreamND();
 Viscosity_Inf   = config->GetViscosity_FreeStreamND();
 Thickness_Inf = total_thickness[nLayer-1][0];
 for (iPoint = 1; iPoint < nPoint; iPoint++)
  Thickness_Inf = max(Thickness_Inf, total_thickness[nLayer-1][iPoint]);
 Thickness_Inf  /= config->GetLength_Ref();
/*--- Get physical quantities ---*/
 su2double *Density_Inf_Layer, *Thickness_Inf_Layer;
 Density_Inf_Layer = NULL; Thickness_Inf_Layer = NULL;
 if(multilayer) {
  Density_Inf_Layer   = new su2double[nLayer];
  Thickness_Inf_Layer = new su2double[nLayer];
  for(iLayer = 0; iLayer < nLayer; iLayer++){
   Density_Inf_Layer[iLayer]   = Density_Inf;
   Thickness_Inf_Layer[iLayer] = Thickness_Inf/su2double(nLayer);
  }
  nodes = new CMultilayerFilmVariable(Density_Inf_Layer,Thickness_Inf_Layer,Velocity_Inf,layer_toll,nPoint,nDim,nVar,config);
 }
 else {
  nodes = new CFilmVariable(Density_Inf,Thickness_Inf,Velocity_Inf,layer_toll,nPoint,nDim,nVar,config,1);
 }

 SetBaseClassPointerToNodes();

 nodes->ImportBottom_Topography(Bottom_Topography);
 nodes->InitializeBottom_BC(config);

 if(Density_Inf_Layer   != NULL) delete [] Density_Inf_Layer;
 if(Thickness_Inf_Layer != NULL) delete [] Thickness_Inf_Layer;

 /* Store the initial CFL number for all grid points. */  
 const su2double CFL = config->GetCFL(MGLevel);
 for (iPoint = 0; iPoint < nPoint; iPoint++) {
  nodes->SetLocalCFL(iPoint, CFL);
 }
 Min_CFL_Local = CFL;
 Max_CFL_Local = CFL;
 Avg_CFL_Local = CFL;
  
 /*--- Store the value of the characteristic primitive variables at the boundaries, except for fluid interface ---*/
 CharacPrimVar = new su2double** [nMarker];
 for(iMarker = 0; iMarker < nMarker; iMarker++) {
  CharacPrimVar[iMarker] = new su2double* [geometry->nVertex[iMarker]];
  for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
   CharacPrimVar[iMarker][iVertex] = new su2double [nPrimVar];
   for(iVar = 0; iVar < nPrimVar; iVar++) {
    CharacPrimVar[iMarker][iVertex][iVar] = 0.0;
   }
  }
 }
 
 CharacPrimVarLayer = new su2double*** [nLayer-1];
 for(iLayer = 0; iLayer < nLayer -1; iLayer++){
  CharacPrimVarLayer[iLayer] = new su2double** [nMarker];
  for(iMarker = 0; iMarker < nMarker; iMarker++) {
   CharacPrimVarLayer[iLayer][iMarker] = new su2double* [geometry->nVertex[iMarker]];
   for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
    CharacPrimVarLayer[iLayer][iMarker][iVertex] = new su2double [nPrimVar];
    for(iVar = 0; iVar < nPrimVar; iVar++) {
     CharacPrimVarLayer[iLayer][iMarker][iVertex][iVar] = 0.0;
    }
   }
  }
 }
 /*--- Store variables necessary for inlet profile ---*/
 Inlet_Ptotal  = new su2double**[nLayer];    
 Inlet_Ttotal  = new su2double**[nLayer];    
 Inlet_FlowDir = new su2double***[nLayer];    

 for(iLayer = 0; iLayer < nLayer; iLayer++) {
  Inlet_Ptotal[iLayer]  = new su2double*[nMarker];
  Inlet_Ttotal[iLayer]  = new su2double*[nMarker];
  Inlet_FlowDir[iLayer] = new su2double**[nMarker];
  for(iMarker = 0; iMarker < nMarker; iMarker++) {
   Inlet_Ptotal[iLayer][iMarker]  = new su2double [geometry->nVertex[iMarker]];
   Inlet_Ttotal[iLayer][iMarker]  = new su2double [geometry->nVertex[iMarker]];
   Inlet_FlowDir[iLayer][iMarker] = new su2double* [geometry->nVertex[iMarker]];
   for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
    Inlet_Ptotal[iLayer][iMarker][iVertex]  = 0.0;
    Inlet_Ttotal[iLayer][iMarker][iVertex]  = 0.0;
    Inlet_FlowDir[iLayer][iMarker][iVertex] = new su2double[nDim];
    for(iDim = 0; iDim < nDim; iDim++)
     Inlet_FlowDir[iLayer][iMarker][iVertex][iDim] = 0.0;
   }
  }
 }

 /*--- Initializate quantities for SlidingMesh Interface ---*/
 SlidingState       = new su2double*** [nMarker];
 SlidingStateNodes  = new int*         [nMarker];
   
 for(iMarker = 0; iMarker < nMarker; iMarker++){
  SlidingState[iMarker]      = NULL;
  SlidingStateNodes[iMarker] = NULL;

  if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){
   SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)];
   SlidingStateNodes[iMarker]  = new int [geometry->GetnVertex(iMarker)];
   for(iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++){
    SlidingState[iMarker][iPoint] = new su2double*[nPrimVarFlow+1];
    SlidingStateNodes[iMarker][iPoint] = 0;
    for(iVar = 0; iVar < nPrimVarFlow+1; iVar++)
     SlidingState[iMarker][iPoint][iVar] = NULL;
   }
  }
 }


/*--- Add the solver name (max 8 characters) ---*/
 SolverName = "FILM.SOL";

} // end function


CFilmSolver::~CFilmSolver(void){

 unsigned short iDim, iMarker, iVar, iLayer;
 unsigned long  iPoint, iEdge, iVertex;

 delete [] delta_layer;
 
 if(multilayer){
  for(iLayer = 0; iLayer < (nLayer-1); iLayer++){
   if(layer_node[iLayer] != NULL){
    for(iPoint = 0; iPoint < nPoint; iPoint++)
      if(layer_node[iLayer][iPoint] != NULL) delete layer_node[iLayer][iPoint];
    delete [] layer_node[iLayer];
   }
   if(layer_edge[iLayer] != NULL){
    for(iEdge = 0; iEdge < nEdge; iEdge++)
      if(layer_edge[iLayer][iEdge] != NULL) delete layer_edge[iLayer][iEdge];
    delete [] layer_edge[iLayer];
   }
   if(layer_vertex[iLayer] != NULL){
    for(iMarker = 0; iMarker < nMarker; iMarker++){
     if(layer_vertex[iLayer][iMarker] != NULL){
      for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex++ )
       if(layer_vertex[iLayer][iMarker][iVertex] != NULL) delete layer_vertex[iLayer][iMarker][iVertex];
      delete [] layer_vertex[iLayer][iMarker];
     }
    }
    delete [] layer_vertex[iLayer];
   } 
  } // end for iLayer
  delete [] layer_node;
  delete [] layer_edge;
  delete [] layer_vertex;
 } // end if multilayer

 if(Bottom_Value != NULL) delete [] Bottom_Value;
 if(Bottom_Topography != NULL){ 
  for(iPoint = 0; iPoint < nPoint; iPoint++)
   delete Bottom_Topography[iPoint];
  delete [] Bottom_Topography;
 }

 if(total_thickness != NULL){
  for(iLayer = 0; iLayer < nLayer; iLayer++){
   if(total_thickness[iLayer] != NULL) delete [] total_thickness[iLayer];
  }
  delete [] total_thickness;
 }

 if(multilayer){
  if(Res_Conv_Layer != NULL){
   for(iLayer = 0; iLayer < (nLayer-1); iLayer++){
    if(Res_Conv_Layer[iLayer] != NULL) delete [] Res_Conv_Layer[iLayer];
   }
  delete [] Res_Conv_Layer;
  }
 }
 
 if(FluidModel != NULL) delete FluidModel;
 if (nodes != nullptr)  delete nodes;

 if(SlidingState != NULL) {
  for(iMarker = 0; iMarker < nMarker; iMarker++) {
   if(SlidingState[iMarker] != NULL) {
    for(iVertex = 0; iVertex < nPoint; iVertex++)
     if(SlidingState[iMarker][iVertex] != NULL){
      for (iVar = 0; iVar < nPrimVar+1; iVar++)
       delete [] SlidingState[iMarker][iVertex][iVar];
      delete [] SlidingState[iMarker][iVertex];
     }
    delete [] SlidingState[iMarker];
   }
  }
  delete [] SlidingState;
 }
 if (SlidingStateNodes != NULL){
  for(iMarker = 0; iMarker < nMarker; iMarker++){
   if (SlidingStateNodes[iMarker] != NULL)
    delete [] SlidingStateNodes[iMarker];  
  }
  delete [] SlidingStateNodes;
 }

 if (CharacPrimVar != NULL) {
  for(iMarker = 0; iMarker < nMarker; iMarker++) {
   for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
    delete [] CharacPrimVar[iMarker][iVertex];
   delete [] CharacPrimVar[iMarker];
  }
  delete [] CharacPrimVar;
 }

 if (CharacPrimVarLayer != NULL) {
  for(iLayer = 0; iLayer < (nLayer - 1); iLayer++){
   for(iMarker = 0; iMarker < nMarker; iMarker++) {
    for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
     delete [] CharacPrimVarLayer[iLayer][iMarker][iVertex];
    delete [] CharacPrimVarLayer[iLayer][iMarker];
   }
   delete [] CharacPrimVarLayer[iLayer];
  }
  delete [] CharacPrimVarLayer;
 }

 if (Inlet_Ptotal != NULL) {
  for(iLayer = 0; iLayer < nLayer; iLayer++) {
   if(Inlet_Ptotal[iLayer] != NULL) {
    for(iMarker = 0; iMarker < nMarker; iMarker++){
     if (Inlet_Ptotal[iLayer][iMarker] != NULL)
      delete [] Inlet_Ptotal[iLayer][iMarker];
    }
    delete [] Inlet_Ptotal[iLayer];
   }
  }
  delete [] Inlet_Ptotal;
 }

 if (Inlet_Ttotal != NULL) {
  for(iLayer = 0; iLayer < nLayer; iLayer++) {
   if(Inlet_Ttotal[iLayer] != NULL) {
    for(iMarker = 0; iMarker < nMarker; iMarker++){
     if (Inlet_Ttotal[iLayer][iMarker] != NULL)
      delete [] Inlet_Ttotal[iLayer][iMarker];
    }
    delete [] Inlet_Ttotal[iLayer];
   }
  }
  delete [] Inlet_Ttotal;
 }
   
 if (Inlet_FlowDir != NULL) {
  for(iLayer = 0; iLayer < nLayer; iLayer++) {
   if(Inlet_FlowDir[iLayer] != NULL) {
    for(iMarker = 0; iMarker < nMarker; iMarker++){
     if (Inlet_FlowDir[iLayer][iMarker] != NULL)
      for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex ++){
       if(Inlet_FlowDir[iLayer][iMarker][iVertex] != NULL)
        delete [] Inlet_FlowDir[iLayer][iMarker][iVertex];
      }
      delete [] Inlet_FlowDir[iLayer][iMarker];
    }
    delete [] Inlet_FlowDir[iLayer];
   }  
  }
  delete [] Inlet_FlowDir;
 }

return;

} // end function


void CFilmSolver::SetNondimensionalization(CConfig *config, unsigned short iMesh){

 unsigned short iDim, iVar;

 su2double Density_FreeStream = 0.0, Temperature_FreeStream = 0.0, Pressure_FreeStream;
 su2double Viscosity_FreeStream, ModVel_FreeStream, Thickness_Ref, Force_Ref = 0.0, Viscosity_Ref = 0.0;
 su2double Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Temperature_Ref = 0.0, Velocity_Ref = 0.0;
 su2double Pressure_FreeStreamND, Density_FreeStreamND, Temperature_FreeStreamND;
 su2double *Velocity_FreeStreamND = new su2double[nDim];
 su2double Total_UnstTimeND, Delta_UnstTimeND, Time_Ref, Delta_Ref;
 su2double Mach     = config->GetMach();
 su2double Reynolds_Delta, Reynolds_L;
 
/*--- Compute dimensional free-stream values. ---*/
 Density_FreeStream     = config->GetInc_Density_Init();      config->SetDensity_FreeStream(Density_FreeStream);
 Temperature_FreeStream = config->GetInc_Temperature_Init();  config->SetTemperature_FreeStream(Temperature_FreeStream);
 ModVel_FreeStream      = 0.0;
 if(multilayer)
  Thickness_Ref = Thickness_Inf/su2double(nLayer);
 else
  Thickness_Ref = Thickness_Inf;

 for(iDim = 0; iDim < nDim; iDim++){
  ModVel_FreeStream += config->GetInc_Velocity_Init()[iDim]*config->GetInc_Velocity_Init()[iDim];
  config->SetVelocity_FreeStream(config->GetInc_Velocity_Init()[iDim],iDim);
 }
 ModVel_FreeStream   = sqrt(ModVel_FreeStream);                               
 config->SetModVel_FreeStream(ModVel_FreeStream);
 Pressure_FreeStream = ModVel_FreeStream*ModVel_FreeStream*Density_FreeStream;
 config->SetPressure_FreeStream(Pressure_FreeStream); 

 config->SetMu_RefND(config->GetMu_Ref());
 config->SetMu_ConstantND(config->GetMu_Constant());

 switch(config->GetKind_FluidModel()) { 
 case CONSTANT_DENSITY:
  FluidModel = new CConstantDensity(Density_FreeStream, config->GetSpecific_Heat_Cp());
  FluidModel->SetTDState_T(Temperature_FreeStream);
  FluidModel->SetLaminarViscosityModel(config);
  break;
 }

 Viscosity_FreeStream = FluidModel->GetLaminarViscosity();
 config->SetViscosity_FreeStream(Viscosity_FreeStream);
 
 Reynolds_Delta = config->GetReynolds();
 Reynolds_L = Density_FreeStream*ModVel_FreeStream/Viscosity_FreeStream; 
 Delta_Ref = Reynolds_Delta/Reynolds_L;

/*--- The non-dim. scheme for incompressible flows uses the following ref. values:
      Reference length      = 1 m (fixed by default, grid in meters)
      Reference density     = liquid density or freestream (input)
      Reference velocity    = liquid velocity or freestream (input)
      Reference temperature = liquid temperature or freestream (input)
      Reference pressure    = Reference density * Reference velocity * Reference velocity
      Reference viscosity   = Reference Density * Reference velocity * Reference length
      This is the same non-dim. scheme as in the compressible solver.
      Note that the Re and Re Length are not used as part of initialization. ---*/
 
 if(config->GetRef_Inc_NonDim() == DIMENSIONAL) {
  Density_Ref     = 1.0;
  Length_Ref      = 1.0;
  Velocity_Ref    = 1.0;
  Time_Ref        = 1.0;
  Temperature_Ref = 1.0;
  Pressure_Ref    = 1.0;
 }
 else if(config->GetRef_Inc_NonDim() == INITIAL_VALUES) {
  Length_Ref      = 1.0;
  Density_Ref     = Density_FreeStream;
  Velocity_Ref    = ModVel_FreeStream;
  Time_Ref        = Length_Ref / Velocity_Ref;
  Temperature_Ref = Temperature_FreeStream;
  Pressure_Ref    = -0.5*g_const*Density_Ref*Delta_Ref;
 } 
 else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
  Length_Ref      = config->GetLength_Ref();
  Density_Ref     = config->GetInc_Density_Ref();
  Velocity_Ref    = config->GetInc_Velocity_Ref();
  Time_Ref        = Length_Ref / Velocity_Ref;
  Temperature_Ref = config->GetInc_Temperature_Ref();
  Pressure_Ref    = -0.5*g_const*Density_Ref*Delta_Ref;
 }

 if( (Density_Ref == 0.0) || (Velocity_Ref == 0.0))
  SU2_MPI::Error(string("Invalid reference value for the velocity or the density."), CURRENT_FUNCTION);

 g_const       = g_const*Time_Ref*Time_Ref/Length_Ref;
 Viscosity_Ref = Density_Ref*Delta_Ref*Velocity_Ref; //  Velocity_Ref*Length_Ref;
 Force_Ref     = Density_Ref*Velocity_Ref*Velocity_Ref*Length_Ref*Length_Ref;//*Delta_Ref*Delta_Ref;

 config->SetTime_Ref(Time_Ref); 
 config->SetLength_Ref(Length_Ref); 
 config->SetDensity_Ref(Density_Ref);
 config->SetVelocity_Ref(Velocity_Ref);
 config->SetTemperature_Ref(Temperature_Ref);
 config->SetPressure_Ref(Pressure_Ref);
 config->SetForce_Ref(Force_Ref);
 config->SetViscosity_Ref(Viscosity_Ref);

 Density_FreeStreamND  = Density_FreeStream/config->GetDensity_Ref();   
 Pressure_FreeStreamND = Pressure_FreeStream/config->GetPressure_Ref(); 
 config->SetDensity_FreeStreamND(Density_FreeStreamND);
 config->SetPressure_FreeStreamND(Pressure_FreeStreamND);
 config->SetViscosity_FreeStreamND(Viscosity_FreeStream/Viscosity_Ref);
 Reynolds_L *= Length_Ref;
//cout<<"Viscosity Ref "<<Viscosity_Ref<<" and Viscosity ND "<<Viscosity_FreeStream/Viscosity_Ref<<endl;
 for (iDim = 0; iDim < nDim; iDim++) {
  Velocity_FreeStreamND[iDim] = config->GetVelocity_FreeStream()[iDim]/Velocity_Ref; 
  config->SetVelocity_FreeStreamND(Velocity_FreeStreamND[iDim], iDim);
 }
 
 Temperature_FreeStreamND = Temperature_FreeStream/config->GetTemperature_Ref(); 
 config->SetTemperature_FreeStreamND(Temperature_FreeStreamND);

 Total_UnstTimeND = config->GetTotal_UnstTime() / Time_Ref;    config->SetTotal_UnstTimeND(Total_UnstTimeND);
 Delta_UnstTimeND = config->GetDelta_UnstTime() / Time_Ref;    config->SetDelta_UnstTimeND(Delta_UnstTimeND);

 config->SetMu_ConstantND(config->GetMu_Constant()/Viscosity_Ref);
 config->SetMu_RefND(config->GetMu_Ref()/Viscosity_Ref);
//cout<<"Mu_ConstantND = "<<config->GetMu_Constant()/Viscosity_Ref<<" e SetMu_RefND "<<config->GetMu_Ref()/Viscosity_Ref<<endl; 
 bool gravity_force = config->GetGravityForce();
 if(gravity_force){
  body_force_vector = config->GetBody_Force_Vector();
  for(iDim = 0; iDim < nDim; iDim++) body_force_vector[iDim] = body_force_vector[iDim]*Time_Ref*Time_Ref/Length_Ref;
  g_const = body_force_vector[nADim];
 }

 delete FluidModel;
 
 switch (config->GetKind_FluidModel()) {
 case CONSTANT_DENSITY:
  FluidModel = new CConstantDensity(Density_FreeStreamND, config->GetSpecific_Heat_CpND() );
  FluidModel->SetTDState_T(Temperature_FreeStreamND);
  FluidModel->SetLaminarViscosityModel(config);
  break;
 }

// delete Velocity_FreeStreamND;

} // end function


void CFilmSolver::ReadTopography(su2double **topography, CGeometry* geometry, CConfig* config){
 
 unsigned long iPoint;
 unsigned short iDim;

 if(config->GetBottom_Top_Type() == UNIFORM){
  for(iPoint = 0; iPoint < nPoint; iPoint++){
   for(iDim = 0; iDim < nADim; iDim++)
    topography[iPoint][iDim] = geometry->node[iPoint]->GetCoord(iDim);
   topography[iPoint][nADim] = config->GetUniform_Bottom_Value();
  }
 }

 if(config->GetBottom_Top_Type() == STRAIGHT && nDim == 2){
  su2double m,q,xx;
  q = config->GetBottom_Value(0);
  m = tan(config->GetBottom_Value(1)*PI_NUMBER/180.0);
  for(iPoint = 0; iPoint < nPoint; iPoint++){
   xx = (m*geometry->node[iPoint]->GetCoord(1) - m*q + geometry->node[iPoint]->GetCoord(0))/(1.0+m*m);
   topography[iPoint][0] = xx;   
   topography[iPoint][1] = q + m*xx;
  }
 } else if(config->GetBottom_Top_Type() == STRAIGHT && nDim != 2) {
  SU2_MPI::Error(string("Invalid topography for the current dimension."), CURRENT_FUNCTION);
 }

 if(config->GetBottom_Top_Type() == PLANE && nDim == 3){
  su2double a,b,c,xx,yy;
  a = config->GetBottom_Value(0);
  b = tan(config->GetBottom_Value(1)*PI_NUMBER/180.0);
  c = tan(config->GetBottom_Value(2)*PI_NUMBER/180.0);
  for(iPoint = 0; iPoint < nPoint; iPoint++){
   xx = (a*b + geometry->node[iPoint]->GetCoord(2)*b + geometry->node[iPoint]->GetCoord(0)*b)/(1.0 +b*b); 
   topography[iPoint][0] = xx;
   yy = c*(b*b +1.0)/b/(c*c + 1.0)*(xx - geometry->node[iPoint]->GetCoord(0)) + geometry->node[iPoint]->GetCoord(1);
   topography[iPoint][1] = yy;
   topography[iPoint][nADim] = a + b*xx + c*yy;
  }
 } else if(config->GetBottom_Top_Type() == PLANE && nDim != 3) {
  SU2_MPI::Error(string("Invalid topography for the current dimension."), CURRENT_FUNCTION);
 }

 if(config->GetBottom_Top_Type() == REAL_DATA){

  const string bottom_name = config->GetBottom_FileName();
  ifstream bottom_file(bottom_name, ifstream::in);
  
  if(bottom_file.fail()){
   SU2_MPI::Error(string("Bottom topography file doesn't exist."), CURRENT_FUNCTION);
  } else {
   string text_line;
   string::size_type position;
   while( getline (bottom_file, text_line) ){

    position = text_line.find ("NPOIN=",0);
    if(position != string::npos) {

     while( getline (bottom_file, text_line) ){
      istringstream point_line(text_line);
      passivedouble Coords[3] = {0.0,0.0,0.0};
      if(nDim == 2){
       point_line >> Coords[0];
       point_line >> Coords[1];
      } else {
       point_line >> Coords[0];
       point_line >> Coords[1];
       point_line >> Coords[2];
      }
      point_line >> iPoint;
      for(iDim = 0; iDim < nDim; iDim++)
       topography[iPoint][iDim] = Coords[iDim];
     }
    }
   }
  } // end else

 } // end if REAL_DATA

 if(topography == NULL)
  SU2_MPI::Error(string("Impossibile to intitialize Bottom Topography."), CURRENT_FUNCTION);

} // end function


int CFilmSolver::GetnSlidingStates(unsigned short val_marker, unsigned long val_vertex){return SlidingStateNodes[val_marker][val_vertex];}


su2double CFilmSolver::GetSlidingState(unsigned short val_marker, unsigned long val_vertex, unsigned short val_state, 
                                       unsigned long donor_index) { 
 return SlidingState[val_marker][val_vertex][val_state][donor_index];
} // end function
 

void CFilmSolver::SetSlidingStateStructure(unsigned short val_marker, unsigned long val_vertex){

 int iVar;
 for(iVar = 0; iVar < nPrimVar+1; iVar++){
  if(SlidingState[val_marker][val_vertex][iVar] != NULL )
   delete [] SlidingState[val_marker][val_vertex][iVar];
 }

 for(iVar = 0; iVar < nPrimVar+1; iVar++)
  SlidingState[val_marker][val_vertex][iVar] = new su2double[GetnSlidingStates(val_marker, val_vertex)];

} // end function


su2double CFilmSolver::Hermite(su2double coord, unsigned short order){

 if(order == 0){
  return 1.0;
 } else if(order == 1) {
  return 2.0*coord;
 } else {
  return 2.0*coord*Hermite(coord, order-1) + 2.0*(order-1)*Hermite(coord, order-2);
 }

} // end function


su2double CFilmSolver::IntegrateHermite(su2double a, su2double b, unsigned short order){

 su2double I = 0.0;
 su2double l = (b-a);
 su2double m = (a+b)/2.0;
 
 I = (Hermite(a, order) + 4*Hermite(m, order) + Hermite(b, order))*l/6.0;
 
 return I;

} // end function


void CFilmSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long TimeIter){
 
 unsigned long iPoint;
 unsigned short iMesh, iLayer;
     
 bool restart = (config->GetRestart() || config->GetRestart_Flow());
   
 /*--- Check if a verification solution is to be computed. ---*/
 if( (VerificationSolution) && (TimeIter == 0) && !restart) {  

 /*--- Loop over the multigrid levels. ---*/
  for(iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
       
  /*--- Loop over all grid points. ---*/
   for(iLayer = 0; iLayer < nLayer; iLayer++){
    for(iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
   
    /* Set the pointers to the coordinates and solution of this DOF. */
     const su2double *coor = NULL;
     if(iLayer == (nLayer-1))
      coor = geometry[iMesh]->node[iPoint]->GetCoord();
     else
      coor = layer_node[iLayer][iPoint]->GetCoord();
     
     su2double *solDOF = solver_container[iMesh][FLOW_SOL]->GetNodes()->GetLayerSolution(iPoint, iLayer);
     
     /* Set the solution in this DOF to the initial condition provided by
      the verification solution class. This can be the exact solution,
      but this is not necessary. */
     VerificationSolution->GetInitialCondition(coor, solDOF);
    }
   }
  } // end for nMGLevels
 } // end if !restart

 if(restart)
  SU2_MPI::Error(string("Restart method not yet implemented for THIN_FILM solver."), CURRENT_FUNCTION);
   
} // end function


void CFilmSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output){

 unsigned long iPoint, jPoint, iEdge, iVertex;
 unsigned short iLayer, iVar, iMarker, iDim;
 unsigned long ErrorCounter = 0;
 bool outlet = ((config->GetnMarker_Outlet() != 0)); //  contains also SuperonicOutlet
 bool dynamic_grid = config->GetDynamic_Grid();
 su2double sum_thickness;

/*--- Set the primitive variables ---*/   
 ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);

/*--- Update the beta value based on the maximum velocity / viscosity. ---*/
 SetBeta_Parameter(geometry, solver_container, config, iMesh);

/*--- Compute properties needed for mass flow BCs. ---*/
 if(outlet) GetOutlet_Properties(geometry, config, iMesh, Output);

/*--- Error message ---*/   
 if (config->GetComm_Level() == COMM_FULL) {
  if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
 }


 for(iPoint = 0; iPoint < nPoint; iPoint++){
  sum_thickness = 0.0;
  for(iLayer = 0; iLayer < nLayer; iLayer++){
   sum_thickness += nodes->GetLayerPrimitive(iPoint, nADim + 2, iLayer);
   total_thickness[iLayer][iPoint] = sum_thickness; 
  }
 } // end loop iPoint

 su2double  dx, dy;
 su2double  TimeStep;
 su2double  new_coord, Displacement;
 su2double* aux_normal = new su2double[nDim];
 TimeStep = config->GetDelta_UnstTimeND();

 /*--- Loop updating Point coordinates, grid velocity and volume ---*/
 for(iPoint = 0; iPoint < nPoint; iPoint++){
  for(iLayer = 0; iLayer < (nLayer-1); iLayer++){
   new_coord = pow(total_thickness[iLayer][iPoint], 2);
   dx = layer_node[iLayer][iPoint]->GetCoord(0) - Bottom_Topography[iPoint][0];
   new_coord = sqrt(new_coord - dx*dx);
   new_coord += Bottom_Topography[iPoint][nADim];
   Displacement = new_coord - layer_node[iLayer][iPoint]->GetCoord(nADim);
   layer_node[iLayer][iPoint]->SetCoord(nADim, new_coord);
   if(dynamic_grid)
    layer_node[iLayer][iPoint]->SetGridVel(nADim, Displacement/TimeStep);
  }
  if(nZone == 1){
   new_coord = pow(total_thickness[nLayer-1][iPoint], 2);
   dx = geometry->node[iPoint]->GetCoord(0) - Bottom_Topography[iPoint][0];
   new_coord = sqrt(new_coord - dx*dx);
   new_coord += Bottom_Topography[iPoint][nADim];
   Displacement = new_coord - geometry->node[iPoint]->GetCoord(nADim);
   geometry->node[iPoint]->SetCoord(nADim, new_coord);
   if(dynamic_grid)
    geometry->node[iPoint]->SetGridVel(nADim, Displacement/TimeStep);
  }
 }

 /*--- Loop updating Edge coordinates and normal ---*/
// O mettere if(nDim==2) qui o toglierlo nel constructor e qui
 for(iLayer = 0; iLayer < (nLayer-1); iLayer++) {
  for(iEdge = 0; iEdge < nEdge; iEdge++) {
   iPoint = layer_edge[iLayer][iEdge]->GetNode(0);
   jPoint = layer_edge[iLayer][iEdge]->GetNode(1);
   for(iDim = 0; iDim < nDim; iDim++){
    if(nDim==2)
     aux_normal[iDim] = layer_node[iLayer][iPoint]->GetCoord(iDim) - layer_node[iLayer][jPoint]->GetCoord(iDim);
    else // da verificare!!!
     aux_normal[iDim] = layer_node[iLayer][jPoint]->GetCoord(iDim) - layer_node[iLayer][iPoint]->GetCoord(iDim);
   }
   layer_edge[iLayer][iEdge]->SetNormal(aux_normal);
  }
 }
 if(nZone == 1){
  for(iEdge = 0; iEdge < nEdge; iEdge++) {
   iPoint = geometry->edge[iEdge]->GetNode(0);
   jPoint = geometry->edge[iEdge]->GetNode(1);
   for(iDim = 0; iDim < nDim; iDim++){
    if(nDim==2)
     aux_normal[iDim] = geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim);
    else // da verificare!!!
     aux_normal[iDim] = geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim); 
   }
   geometry->edge[iEdge]->SetNormal(aux_normal);
  }
 }

 /*--- Loop updating Vertex coordinates and volume ---*/
 for(iLayer = 0; iLayer < (nLayer-1); iLayer++) {
  for(iMarker = 0; iMarker < nMarker; iMarker++) {
   if(config->GetMarker_All_KindBC(iMarker) != FLUID_INTERFACE){
    for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
     layer_vertex[iLayer][iMarker][iVertex]->SetNormal(geometry->vertex[iMarker][iVertex]->GetNormal());
    }
   }
  }
 }

 delete [] aux_normal;

 if(Output){
  for(iPoint = 0; iPoint < nPoint; iPoint++){
 /*--- Applying PE hypothesis ---*/
   su2double press_aux;
   for(iLayer = 0; iLayer < nLayer; iLayer++){
    press_aux = -0.5*nodes->GetLayerPrimitive(iPoint, nADim+1, iLayer)*g_const*nodes->GetLayerPrimitive(iPoint, nADim+2, iLayer); 
    nodes->SetLayerSolution(iPoint,  0, press_aux, iLayer);
    nodes->SetLayerPrimitive(iPoint, 0, press_aux, iLayer);
   }
 /*--- Applying PEV2 contribution ---*/
   if(viscous){
   }
  } // end loop iPoint
 } // end if Output

/* --- Resoconto delle variabili fisiche --- //
if(Output){
cout << "----- Physical variables ------" << endl;
cout << "--- Velocity ----"<< endl;
iLayer = nLayer;
while(iLayer != 0){
 cout << '\t' << nodes->GetLayerPrimitive(0,1,iLayer-1)<<" "<< nodes->GetLayerPrimitive(1,1,iLayer-1)<<" "<< nodes->GetLayerPrimitive(nPoint/2,1,iLayer-1)<<" "<< nodes->GetLayerPrimitive(nPoint/2+1,1,iLayer-1)<<" "<< nodes->GetLayerPrimitive(nPoint-2,1,iLayer-1)<<" "<< nodes->GetLayerPrimitive(nPoint-1,1,iLayer-1)  << endl;
iLayer--;
} 
cout << '\t' << "Point 0   Point 1   Point " <<nPoint/2+1<< "   Point " <<nPoint/2+2<< "   Point " <<nPoint-1<<"   Point " <<nPoint<< endl;
cout << endl;
cout << "--- Pressure ----"<< endl;
iLayer = nLayer;
while(iLayer != 0){
 cout << '\t' << nodes->GetLayerPrimitive(0,0,iLayer-1)<<" "<< nodes->GetLayerPrimitive(1,0,iLayer-1)<<" "<< nodes->GetLayerPrimitive(nPoint/2,0,iLayer-1)<<" "<< nodes->GetLayerPrimitive(nPoint/2+1,0,iLayer-1)<<" "<< nodes->GetLayerPrimitive(nPoint-2,0,iLayer-1)<<" "<< nodes->GetLayerPrimitive(nPoint-1,0,iLayer-1)  << endl;
iLayer--;
} 
cout << '\t' << "Point 0   Point 1   Point " <<nPoint/2+1<< "   Point " <<nPoint/2+2<< "   Point " <<nPoint-1<<"   Point " <<nPoint<< endl;
cout << endl;
cout << "--- Delta k-1,k ----"<< endl;
iLayer = nLayer;
while(iLayer != 0){
 cout << '\t' << nodes->GetLayerPrimitive(0,nADim+2,iLayer-1)<<" "<< nodes->GetLayerPrimitive(1,nADim+2,iLayer-1)<<" "<< nodes->GetLayerPrimitive(nPoint/2,nADim+2,iLayer-1)<<" "<< nodes->GetLayerPrimitive(nPoint/2+1,nADim+2,iLayer-1)<<" "<< nodes->GetLayerPrimitive(nPoint-2,nADim+2,iLayer-1)<<" "<< nodes->GetLayerPrimitive(nPoint-1,nADim+2,iLayer-1)  << endl;
iLayer--;
} 
cout << '\t' << "Point 0   Point 1   Point " <<nPoint/2+1<< "   Point " <<nPoint/2+2<< "   Point " <<nPoint-1<<"   Point " <<nPoint<< endl;
cout << endl;
cout << "--- Total Thickness ----"<< endl;
iLayer = nLayer;
while(iLayer != 0){
 cout << '\t' << total_thickness[iLayer-1][0]<<" "<<total_thickness[iLayer-1][1]<<" "<<total_thickness[iLayer-1][nPoint/2]<<" "<<total_thickness[iLayer-1][nPoint/2+1]<<" "<<total_thickness[iLayer-1][nPoint-2]<<" "<<total_thickness[iLayer-1][nPoint-1]<< endl;
iLayer--;
} 
cout << '\t' << "Point 0   Point 1   Point " <<nPoint/2+1<< "   Point " <<nPoint/2+2<< "   Point " <<nPoint-1<<"   Point " <<nPoint<< endl;
cout << endl;
}*/

} // end function


unsigned long CFilmSolver::SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output) {
   
 unsigned long iPoint, nonPhysicalPoints = 0;
 bool physical = true;

 for (iPoint = 0; iPoint < nPoint; iPoint++) {    
 /*--- Incompressible flow, primitive variables nADim+3, (p, vx, vy, rho, delta) ---*/ 
  if(multilayer)
   physical = static_cast<CMultilayerFilmVariable*> (nodes)->SetPrimVar(iPoint, solver_container[FLOW_SOL]->GetFluidModel());
  else
   physical = nodes->SetPrimVar(iPoint, solver_container[FLOW_SOL]->GetFluidModel());
      
 /* Check for non-realizable states for reporting. */   
  if(!physical) nonPhysicalPoints++;     

 /*--- Initialize the convective, source and viscous residual vector ---*/    
  if(!Output){
  for(unsigned short iLayer = 0; iLayer < nLayer; iLayer++)
   Layer_SysRes.SetBlock_Zero(iPoint + iLayer*nPoint);
  } 
 }

 return nonPhysicalPoints;
} // end function


void CFilmSolver::Set_OldSolution(CGeometry *geometry){

 unsigned long  iPoint;
 unsigned short iLayer;
 su2double* solution;
 
 for(iLayer = 0;  iLayer < nLayer; iLayer++){
  for(iPoint = 0; iPoint < nPoint; iPoint++){
   solution = nodes->GetLayerSolution(iPoint, iLayer);
   nodes->SetLayerSolution_Old(iPoint, solution, iLayer);
  } // end for iPoint
 } // end for iLayer

} // end function


void CFilmSolver::Set_NewSolution(CGeometry *geometry){

 unsigned long  iPoint;
 unsigned short iLayer, iVar;
 su2double* solution;
 
 for(iLayer = 0;  iLayer < nLayer; iLayer++){
  for(iPoint = 0; iPoint < nPoint; iPoint++){
   solution = nodes->GetLayerSolution(iPoint, iLayer);
    for(iVar = 0; iVar < nVar; iVar++){
     nodes->SetLayerSolution_New(iPoint, iVar, solution[iVar], iLayer);
    }
  } // end for iPoint
 } // end for iLayer

} // end function


void CFilmSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                 unsigned short iMesh, unsigned long Iteration) {
   
 su2double *Normal, Area, Vol, Mean_ProjVel = 0.0, Mean_Density, Lambda, Local_Delta_Time,
	    Local_Delta_Time_Visc, Global_Delta_Time = 1E6, K_v = 0.25, ProjVel, ProjVel_i, ProjVel_j;
   
 unsigned long iEdge, iVertex, iPoint, jPoint;
 unsigned short iDim, iMarker;

 Min_Delta_Time = 1.E30; Max_Delta_Time = 0.0;

/*--- Set maximum inviscid eigenvalue to zero, and compute sound speed and viscosity ---*/  
  for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
    nodes->SetMax_Lambda_Inv(iPoint,  config->GetVelocity_Ref());
    nodes->SetMax_Lambda_Visc(iPoint, 0.0);
  }

 for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
     
 /*--- Point identification, Normal vector and area ---*/     
  iPoint = geometry->edge[iEdge]->GetNode(0);
  jPoint = geometry->edge[iEdge]->GetNode(1);    
  Normal = geometry->edge[iEdge]->GetNormal();
  Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; 
  Area = sqrt(Area);     
 /*--- Mean Values ---*/     
  Mean_ProjVel    = 0.5 * (nodes->GetProjVel(iPoint, Normal) + nodes->GetProjVel(jPoint, Normal));
  Mean_Density    = 0.5 * (nodes->GetDensity(iPoint) + nodes->GetDensity(jPoint));
 /*--- Adjustment for grid movement ---*/     
  if(dynamic_grid) {
   su2double *GridVel_i = geometry->node[iPoint]->GetGridVel();
   su2double *GridVel_j = geometry->node[jPoint]->GetGridVel();
   ProjVel_i = 0.0; ProjVel_j = 0.0;
   // nDim because we reconstruct also the normal grid velocity component
   for(iDim = 0; iDim < nDim; iDim++) { 
    ProjVel_i += GridVel_i[iDim]*Normal[iDim];
    ProjVel_j += GridVel_j[iDim]*Normal[iDim];
   }
   Mean_ProjVel -= 0.5 * (ProjVel_i + ProjVel_j);
  }
  Lambda = fabs(Mean_ProjVel);
  if (geometry->node[iPoint]->GetDomain()) nodes->AddMax_Lambda_Inv(iPoint,Lambda);
  if (geometry->node[jPoint]->GetDomain()) nodes->AddMax_Lambda_Inv(jPoint,Lambda);
 }

/*--- Loop boundary edges ---*/   
 for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
  if((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
     (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
   for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {       
   /*--- Point identification, Normal vector and area ---*/ 
    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
    Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
    Area = 0.0; for(iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);   
    /*--- Mean Values ---*/   
    Mean_ProjVel = nodes->GetProjVel(iPoint,Normal);
    Mean_Density = nodes->GetDensity(iPoint);
    /*--- Adjustment for grid movement ---*/       
    if (dynamic_grid) {
     su2double *GridVel = geometry->node[iPoint]->GetGridVel();
     ProjVel = 0.0;
     for (iDim = 0; iDim < nDim; iDim++)
      ProjVel += GridVel[iDim]*Normal[iDim];
     Mean_ProjVel -= ProjVel;
    }
   Lambda = fabs(Mean_ProjVel);
   if (geometry->node[iPoint]->GetDomain()) {nodes->AddMax_Lambda_Inv(iPoint,Lambda);}
   }
  }
 }

/*--- Each element uses their own speed, steady state simulation ---*/  
 for (iPoint = 0; iPoint < nPointDomain; iPoint++) {   
  Vol = geometry->node[iPoint]->GetVolume();    
  if (Vol != 0.0) {
   Local_Delta_Time = nodes->GetLocalCFL(iPoint)*Vol / nodes->GetMax_Lambda_Inv(iPoint);
   Local_Delta_Time_Visc = nodes->GetLocalCFL(iPoint)*K_v*Vol*Vol/ nodes->GetMax_Lambda_Visc(iPoint);
   Local_Delta_Time = min(Local_Delta_Time, Local_Delta_Time_Visc);
   Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
   Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
   Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
   if (Local_Delta_Time > config->GetMax_DeltaTime())
    Local_Delta_Time = config->GetMax_DeltaTime();
    nodes->SetDelta_Time(iPoint,Local_Delta_Time);
   }
   else {
    nodes->SetDelta_Time(iPoint,0.0);
   }   
  }

 /*--- For exact time solution use the minimum delta time of the whole mesh ---*/
 if (config->GetTime_Marching() == TIME_STEPPING) {
 /*--- If the unsteady CFL is set to zero, it uses the defined
       unsteady time step, otherwise it computes the time step based
       on the unsteady CFL ---*/    
  if (config->GetUnst_CFL() == 0.0) {
   Global_Delta_Time = config->GetDelta_UnstTime();
  }
  config->SetDelta_UnstTimeND(Global_Delta_Time);

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){     
  /*--- Sets the regular CFL equal to the unsteady CFL ---*/    
   nodes->SetLocalCFL(iPoint, config->GetUnst_CFL());
   nodes->SetDelta_Time(iPoint, Global_Delta_Time);
   Min_Delta_Time = Global_Delta_Time;
   Max_Delta_Time = Global_Delta_Time;
  }
 }


} // end function

void CFilmSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh){

/* bool classical_rk4 = (config->GetKind_TimeIntScheme_Flow() == CLASSICAL_RK4_EXPLICIT);
 unsigned long  iPoint;
 unsigned short iLayer, iVar, iDim;
 su2double* solution;
 
 for(iLayer = 0;  iLayer < nLayer; iLayer++){
  for(iPoint = 0; iPoint < nPoint; iPoint++){
   solution = nodes->GetLayerSolution(iPoint, iLayer);
   nodes->SetLayerSolution_Old(iPoint, solution, iLayer);
   if(classical_rk4){
    for(iVar = 0; iVar < nVar; iVar++){
     nodes->SetLayerSolution_New(iPoint, iVar, solution[iVar], iLayer);
    }
   }
  } // end for iPoint
 } // end for iLayer
*/
/*-- Updating interface velocity for single zone problems --
 if(nZone == 1){
  for(iPoint = 0; iPoint < nPoint; iPoint++){
   for(iDim = 0; iDim <nADim; iDim++)
    Velocity_Init[iDim] = nodes->GetLayerPrimitive(iPoint, iDim+1, nLayer-1);
  }
 }*/

} // end function


void CFilmSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, 
                                   CConfig *config, unsigned short iMesh){

 su2double *Normal = new su2double[nDim];
 su2double *V_i, *V_j;
 su2double **V_i_Layer, **V_j_Layer;
 unsigned long iEdge, iPoint, jPoint;
 unsigned short iDim, iVar, iLayer;
 
 if(multilayer)
 { V_i_Layer = new su2double*[nLayer-1]; V_j_Layer = new su2double*[nLayer-1];}
 
 for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

  if(multilayer){    
   for(iLayer = 0; iLayer < (nLayer-1); iLayer++){

    iPoint = layer_edge[iLayer][iEdge]->GetNode(0); jPoint = layer_edge[iLayer][iEdge]->GetNode(1);
    numerics->SetCoord(layer_node[iLayer][iPoint]->GetCoord(), layer_node[iLayer][jPoint]->GetCoord());
    layer_edge[iLayer][iEdge]->GetNormal(Normal);
    numerics->SetNormal(Normal);
//cout<<"Check for Edge "<<iEdge<<", between points "<<iPoint<<" & "<<jPoint;
//cout<<", with normal ";for(iDim = 0; iDim < nDim; iDim++)cout<<Normal[iDim]<<" ";cout<<endl;

   /*--- Grid movement ---*/     
    if (dynamic_grid)
     numerics->SetGridVel(layer_node[iLayer][iPoint]->GetGridVel(), layer_node[iLayer][jPoint]->GetGridVel());

   /*--- Get primitive variables ---*/     
    V_i_Layer[iLayer] = nodes->GetLayerPrimitive(iPoint, iLayer);
    V_j_Layer[iLayer] = nodes->GetLayerPrimitive(jPoint, iLayer);
    numerics->SetPrimitive(V_i_Layer[iLayer], V_j_Layer[iLayer]);
 
    numerics->ComputeResidual(Res_Conv_Layer[iLayer], config);

    Layer_SysRes.AddBlock(iPoint + iLayer*nPoint, Res_Conv_Layer[iLayer]);
    Layer_SysRes.SubtractBlock(jPoint + iLayer*nPoint, Res_Conv_Layer[iLayer]);

//cout << "Upwind Residual: "; for(iVar =  0; iVar<nVar; iVar++) cout << Res_Conv_Layer[iLayer][iVar] << " "; cout << endl;
   }
  } // end if(multilayer)

  iPoint = geometry->edge[iEdge]->GetNode(0); jPoint = geometry->edge[iEdge]->GetNode(1);
  numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
  geometry->edge[iEdge]->GetNormal(Normal);
  numerics->SetNormal(Normal);
//cout<<"Check for Edge "<<iEdge<<", between points "<<iPoint<<" & "<<jPoint;
//cout<<", with normal ";for(iDim = 0; iDim < nDim; iDim++)cout<<Normal[iDim]<<" ";cout<<endl;
 /*--- Grid movement ---*/     
  if (dynamic_grid)
   numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
 
 /*--- Get primitive variables ---*/     
  V_i = nodes->GetPrimitive(iPoint); V_j = nodes->GetPrimitive(jPoint);

 /*--- Set conservative variables without reconstruction ---*/   
  numerics->SetPrimitive(V_i, V_j);
 /*--- Compute the residual ---*/    
  numerics->ComputeResidual(Res_Conv,config);

  Layer_SysRes.AddBlock(iPoint + (nLayer-1)*nPoint, Res_Conv);
  Layer_SysRes.SubtractBlock(jPoint + (nLayer-1)*nPoint, Res_Conv);

//cout<< "Top Upwind Residual: "; for(iVar = 0; iVar<nVar; iVar++) cout << Res_Conv[iVar] << " "; cout << endl;

 } // end loop iEdge

 delete [] Normal;

} // end function


void CFilmSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, 
                  unsigned short iMesh, unsigned short iRKStep){ return;} // end function

su2double CFilmSolver::GetInterfaceVar(unsigned short val_marker, unsigned long iPoint, unsigned short iVar){

 unsigned short iState;
 su2double interface_val;
 interface_val = 0.0;

 if(nZone > 1){
  for(iState = 0.0; iState < GetnSlidingStates(val_marker, iPoint); iState++)
   interface_val += GetSlidingState(val_marker,iPoint,0,iState);

 } else {

  switch(iVar){
   case 0:
     interface_val = Pressure_Init;
     break;
   case 1:
     interface_val = Velocity_Init[0];
     break;
   case 2:
    if(nADim == 2)
     interface_val = Velocity_Init[1];
    else
     interface_val = Density_Inf;
    break;      
   case 3:
    interface_val = Density_Inf;
    break;      
  }
 }

 return interface_val;

} // end function


void CFilmSolver::ComputeProfile(CConfig *config, su2double* Coeff, unsigned long iPoint, unsigned short ivar){

 unsigned short iOrder, iLayer;
 unsigned short intMarker, iMarker;
 for(iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
  if(config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE)
   intMarker = iMarker;

 CBlasStructure BlasSolver;
 su2double*  Hermite_Matrix = new su2double[RowOrder*ProfileOrder];

 su2double a_max, b_max, ah, bh, bh_ah;

  /*--- Profile reconstruction ---*/ 
  a_max = 0.0;
  b_max = total_thickness[nLayer-1][iPoint];


 /*-- Setting boundary entries of known term ---*/
  Coeff[0]  = nodes->GetBottom_BC(iPoint, ivar); 
  Coeff[RowOrder-1] = GetInterfaceVar(intMarker, iPoint, ivar);
 
 /*-- Setting boundary entries of the matrix ---*/
  for(iOrder = 0; iOrder < ProfileOrder; iOrder++){
   Hermite_Matrix[RowOrder*iOrder] = Hermite(a_max, iOrder);
   Hermite_Matrix[(RowOrder - 1) + RowOrder*iOrder] = Hermite(b_max, iOrder);
  }

   for(iLayer = 0; iLayer < nLayer; iLayer++ ){
 
    if(iLayer == 0)
     ah = 0.0;
    else
     ah = total_thickness[iLayer-1][iPoint];
    bh = total_thickness[iLayer][iPoint];
    bh_ah = bh - ah;
 
  /*-- Setting entries of known term ---*/
    Coeff[iLayer+1] = nodes->GetLayerPrimitive(iPoint, ivar, iLayer)*bh_ah;

  /*-- Setting integral entries of the matrix ---*/
    for(iOrder = 0; iOrder < ProfileOrder; iOrder++){
     Hermite_Matrix[(iLayer + 1) + RowOrder*iOrder] = IntegrateHermite(ah, bh, iOrder);
    }

   } // end for iLayer
   
 /*------ Solving the systems ------*/ 
 BlasSolver.dgetrs(ProfileOrder, RowOrder, Hermite_Matrix, Coeff);
if(ivar==1){
//cout<<"H = "<<endl;for(int ii=0;ii<ProfileOrder;ii++){cout<<"[ ";for(int jj=0;jj<RowOrder;jj++)cout<<Hermite_Matrix[ii+jj*RowOrder]<<" ";cout<<"]"<<endl;}
cout<<"Coeff = ";for(iOrder=0;iOrder<RowOrder;iOrder++)cout<<Coeff[iOrder]<<" ";cout<<endl;
}

 delete [] Hermite_Matrix;

} // end function


void CFilmSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, 
                   CConfig *config, unsigned short iMesh){

 unsigned short iDim, iOrder;
 unsigned short iVar, iLayer, jLayer;
 unsigned long  iPoint, jPoint, iEdge;
 
/*--- Get the order for profile reconstruction ---*/
 RowOrder     = nLayer + 2; // Nr. of rows = nLayer + 2 bc
 ProfileOrder = nLayer + 2; // Nr. of columns = max Hermite degree
 numerics->SetProfileOrder(ProfileOrder);
 second_numerics->SetProfileOrder(ProfileOrder);

 su2double*  P_Coeff = new su2double[ProfileOrder];
 su2double** V_Coeff = new su2double*[nADim];
 for(iDim = 0; iDim < nADim; iDim++) {
  V_Coeff[iDim] = new su2double[ProfileOrder];
 }

/*--- Loop over all points ---*/
 for(iPoint = 0; iPoint < nPoint; iPoint++) {

/*--- Comupet the coefficient of the profile reconstruction ---*/
 ComputeProfile(config, P_Coeff, iPoint, 0);
 for(iDim = 0; iDim < nADim; iDim++)
  ComputeProfile(config, V_Coeff[iDim], iPoint, iDim+1);

 /*--- Loop over all layers ---*/
 for(iLayer = 0; iLayer < nLayer; iLayer++ ){
 
 /*--- Set incompressible density  ---*/
  numerics->SetDensity(nodes->GetDensity(iPoint), nodes->GetDensity(iPoint));
  second_numerics->SetDensity(nodes->GetDensity(iPoint), nodes->GetDensity(iPoint));

 /*--- Set laminar viscosity  ---*/
  second_numerics->SetLaminarViscosity(nodes->GetLaminarViscosity(iPoint), nodes->GetLaminarViscosity(iPoint)); 

 /*--- Load the volume of the dual mesh cell ---*/
  su2double vol;
  if( iLayer != (nLayer-1) )
   vol = layer_node[iLayer][iPoint]->GetVolume(); 
  else
   vol = geometry->node[iPoint]->GetVolume(); 
  numerics->SetVolume(vol); 
  second_numerics->SetVolume(vol); 
//cout<<"Internal CV Volume = "<<vol<<endl;

 /*--- Set coord ---*/
  iEdge  = geometry->node[iPoint]->GetEdge(0);
  jPoint = geometry->edge[iEdge]->GetNode(1);
  jPoint = (jPoint == iPoint) ? geometry->edge[iEdge]->GetNode(0) : jPoint;
  su2double Aux_Coordjm1[nDim], Aux_Coordj[nDim];

 /*--- Setting Thickness derivative ---*/ 
  for(iDim = 0; iDim < nADim; iDim++){
   if(iLayer == 0){
    Aux_Coordjm1[iDim] = Bottom_Topography[jPoint][iDim];
    Aux_Coordj[iDim]   = layer_node[iLayer][jPoint]->GetCoord(iDim);    
   } else if(iLayer != 0 && iLayer < (nLayer-1)){
    Aux_Coordjm1[iDim] = layer_node[iLayer-1][jPoint]->GetCoord(iDim);
    Aux_Coordj[iDim]   = layer_node[iLayer][jPoint]->GetCoord(iDim);
   } else {
    Aux_Coordjm1[iDim] = layer_node[nLayer-2][jPoint]->GetCoord(iDim);
    Aux_Coordj[iDim]   = geometry->node[jPoint]->GetCoord(iDim);
   }
  }

   if(iLayer == 0){
    Aux_Coordjm1[nADim] = 0.0;
    Aux_Coordj[nADim]   = total_thickness[iLayer][jPoint];    
   } else {
    Aux_Coordjm1[nADim] = total_thickness[iLayer-1][jPoint];
    Aux_Coordj[nADim]   = total_thickness[iLayer][jPoint];
   }

  if(iLayer == 0){
   numerics->SetdThick(Aux_Coordjm1, Aux_Coordj, 0.0, total_thickness[iLayer][iPoint]);
  }
  else {
   numerics->SetdThick(Aux_Coordjm1, Aux_Coordj, total_thickness[iLayer-1][iPoint], total_thickness[iLayer][iPoint]);
  }

 /*--- Set normal pressure and velocity profiles ---*/
//cout<<"K-th Layer "<<iLayer+1<<","<<endl;
  if(iLayer == 0){
   numerics->SetPressureBC(P_Coeff, 0.0, total_thickness[0][iPoint]);
   second_numerics->SetStressBC(V_Coeff,   0.0, total_thickness[0][iPoint]);
  }
  else {
   numerics->SetPressureBC(P_Coeff, total_thickness[iLayer-1][iPoint], total_thickness[iLayer][iPoint]);
   second_numerics->SetStressBC(V_Coeff,   total_thickness[iLayer-1][iPoint], total_thickness[iLayer][iPoint]);
  }

  if(iLayer == 0){
   numerics->SetCoord(Bottom_Topography[iPoint], layer_node[iLayer][iPoint]->GetCoord());
  }
  else if(iLayer != 0 && iLayer < (nLayer-1)){
   numerics->SetCoord(layer_node[iLayer-1][iPoint]->GetCoord(),layer_node[iLayer][iPoint]->GetCoord());
  }
  else {
   numerics->SetCoord(layer_node[nLayer-2][iPoint]->GetCoord(),geometry->node[iPoint]->GetCoord());
  }

  numerics->SetPrimitive(nodes->GetLayerPrimitive(iPoint, iLayer), nodes->GetLayerPrimitive(iPoint, iLayer));
  second_numerics->SetPrimitive(nodes->GetLayerPrimitive(iPoint, iLayer), nodes->GetLayerPrimitive(iPoint, iLayer));

 /*--- Compute source residual ---*/
  numerics->ComputeResidual(Res_Sour, config); 
 /*--- Add the source residual to the total ---*/
  Layer_SysRes.AddBlock(iPoint + iLayer*nPoint, Res_Sour);

 /*--- Compute second source residual ---*/
  second_numerics->ComputeResidual(Res_Sour, config); 
  Layer_SysRes.AddBlock(iPoint + iLayer*nPoint, Res_Sour);

//cout << "Source residual: "; for(iVar =  0; iVar<nVar; iVar++) cout << Res_Sour[iVar] << " "; cout << endl;
  } // end for iLayer
//cout<<"--------------------"<<endl;
 } // end for iPoint

 /*--- Free allocated memory --*/
 delete [] P_Coeff;
 for(iDim = 0; iDim < nADim; iDim++){
  delete [] V_Coeff[iDim];
 }
 delete [] V_Coeff;
 
 return;
} // end function


void CFilmSolver::Centered_Residual (CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh, unsigned short iRKstep){
 SU2_MPI::Error(string("Method not yet implemented for THIN_FILM problems."), CURRENT_FUNCTION);
}
void CFilmSolver::Convective_Residual (CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh, unsigned short iRKstep){
 SU2_MPI::Error(string("Method not yet implemented for THIN_FILM problems."), CURRENT_FUNCTION);
}

void CFilmSolver::ComputeResidual_Multizone(CGeometry *geometry, CConfig *config){

  unsigned short iVar;
  unsigned long iPoint;
  su2double residual;

  /*--- Set Residuals to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++){
    SetRes_BGS(iVar,0.0);
    SetRes_Max_BGS(iVar,0.0,0);
  }

  /*--- Set the residuals ---
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      residual = base_nodes->GetSolution(iPoint,iVar) - base_nodes->Get_BGSSolution_k(iPoint,iVar);
      AddRes_BGS(iVar,residual*residual);
      AddRes_Max_BGS(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
    }
  }

  SetResidual_BGS(geometry, config);
   */

} // end function


void CFilmSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config){ 

 su2double *local_Residual, Vol, Delta, Res;
 unsigned short iVar, jVar, iLayer;
 unsigned long iPoint;

 for(iVar = 0; iVar < nVar; iVar++)
  SetRes_RMS(iVar, 0.0);

/*--- Update the solution  for each layer ---*/
for(iLayer = 0; iLayer < (nLayer - 1); iLayer++){
 for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
  Vol = layer_node[iLayer][iPoint]->GetVolume();
  Delta = nodes->GetDelta_Time(iPoint) / Vol;
     
  local_Residual = Layer_SysRes.GetBlock(iPoint + iLayer*nPoint);

 /*--- Setting the solution for the pressure, no residual required ---*/ 
  Res = local_Residual[0];
  nodes->AddLayerSolution(iPoint, 0, Res, iLayer); 
  for(iVar = 1; iVar < nVar; iVar++) {
   Res = 0.0;
   Res = local_Residual[iVar];
   nodes->AddLayerSolution(iPoint,iVar, -Res*Delta, iLayer);
   }
 }
}

 for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
  Vol = (geometry->node[iPoint]->GetVolume() + geometry->node[iPoint]->GetPeriodicVolume());
  Delta = nodes->GetDelta_Time(iPoint) / Vol;
     
  local_Residual = Layer_SysRes.GetBlock(iPoint + (nLayer - 1)*nPoint);

 /*--- Setting the solution for the pressure, no residual required ---*/ 
  Res = local_Residual[0];
  nodes->AddSolution(iPoint, 0, Res);  
  for(iVar = 1; iVar < nVar; iVar ++ ) {
   Res = 0.0;
   Res = local_Residual[iVar];
   nodes->AddSolution(iPoint, iVar, -Res*Delta);
   AddRes_RMS(iVar, Res*Res);
  }
 }

/*--- MPI solution ---*/   
 InitiateComms(geometry, config, SOLUTION);
 CompleteComms(geometry, config, SOLUTION);

 SetResidual_RMS(geometry, config);
   
/*--- For verification cases, compute the global error metrics. ---*/   
 ComputeVerificationError(geometry, config);

 return;
} // end function 


void CFilmSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep){

 su2double *local_Residual, Vol, Delta, Res;
 unsigned short iVar, iLayer;
 unsigned long iPoint;
 su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  
 for(iVar = 0; iVar < nVar; iVar++)
   SetRes_RMS(iVar, 0.0);
  
 /*--- Update the solution for the layers ---*/  
for(iLayer = 0; iLayer < (nLayer -1); iLayer++){
 for(iPoint = 0; iPoint < nPointDomain; iPoint++){
  Vol = layer_node[iLayer][iPoint]->GetVolume();
  Delta = nodes->GetDelta_Time(iPoint)/Vol;
  local_Residual = Layer_SysRes.GetBlock(iPoint + iLayer*nPoint);
  for(iVar = 1; iVar < nVar; iVar++){
   Res = local_Residual[iVar];
   nodes->AddLayerSolution(iPoint,iVar, -Res*Delta*RK_AlphaCoeff, iLayer);
  }  
 }
}
  
 /*--- Update the solution for the top layer ---*/  
 for(iPoint = 0; iPoint < nPointDomain; iPoint++){
  Vol = geometry->node[iPoint]->GetVolume();
  Delta = nodes->GetDelta_Time(iPoint)/Vol;
  local_Residual = Layer_SysRes.GetBlock(iPoint + (nLayer-1)*nPoint);
  for(iVar = 1; iVar < nVar; iVar++){
   Res = local_Residual[iVar];
   nodes->AddSolution(iPoint,iVar, -Res*Delta*RK_AlphaCoeff);
   AddRes_RMS(iVar, Res*Res);
  }
 }

 /*--- MPI solution ---*/ 
 InitiateComms(geometry, config, SOLUTION);
 CompleteComms(geometry, config, SOLUTION);
//cout<<"---------------------end iRK "<<iRKStep<<"-------------------------"<<endl;  
 /*--- Compute the root mean square residual ---*/ 
 SetResidual_RMS(geometry, config);
  
 /*--- For verification cases, compute the global error metrics. ---*/  
 ComputeVerificationError(geometry, config);

} 
 	
void CFilmSolver::ClassicalRK4_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep){

 su2double *local_Residual, Vol, Delta, Res, tmp_time, tmp_func;
 unsigned short iVar, iLayer;
 unsigned long iPoint;
 
 /*--- Classical RK4 coefficients. ---*/
 su2double RK_FuncCoeff[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
 su2double RK_TimeCoeff[4] = {0.5, 0.5, 1.0, 1.0};
 
 bool adjoint = config->GetContinuous_Adjoint();
 for(iVar = 0; iVar < nVar; iVar++)
  SetRes_RMS(iVar, 0.0);

 /*--- Update the solution for each layer---*/
 for(iLayer = 0; iLayer < (nLayer - 1); iLayer++){
  for(iPoint = 0; iPoint < nPointDomain; iPoint++) {

   Vol = layer_node[iLayer][iPoint]->GetVolume();
   Delta = nodes->GetDelta_Time(iPoint) / Vol;
   local_Residual = Layer_SysRes.GetBlock(iPoint + iLayer*nPoint);
//cout<<"Total iLayer "<<iLayer+1<<" residual ";for(iVar = 0; iVar < nVar; iVar++)cout<<local_Residual[iVar]<<" ";cout<<endl;
 
   tmp_time = -1.0*RK_TimeCoeff[iRKStep]*Delta;
   tmp_func = -1.0*RK_FuncCoeff[iRKStep]*Delta;

    for (iVar = 1; iVar < nVar; iVar++) {
     Res = local_Residual[iVar];
     if (iRKStep < 3) {
     /* Base Solution Update */
     nodes->AddLayerSolution(iPoint, iVar, tmp_time*Res, iLayer); 
     /* New Solution Update */
     nodes->AddLayerSolution_New(iPoint, iVar, tmp_func*Res, iLayer);
     } else {
      nodes->SetLayerSolution(iPoint, iVar, nodes->GetLayerSolution_New(iPoint, iVar, iLayer) + tmp_func*Res, iLayer);
     } 
    }

  } // end for iPoint
 } // end for iLayer

 /*--- Update the solution for the free surface ---*/
 for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
 
  Vol = geometry->node[iPoint]->GetVolume();
  Delta = nodes->GetDelta_Time(iPoint)/Vol;
  local_Residual = Layer_SysRes.GetBlock(iPoint + (nLayer-1)*nPoint);
//cout<<"Total top-layer residual ";for(iVar = 0; iVar < nVar; iVar++)cout<<local_Residual[iVar]<<" ";cout<<endl;

  tmp_time = -1.0*RK_TimeCoeff[iRKStep]*Delta;
  tmp_func = -1.0*RK_FuncCoeff[iRKStep]*Delta;
 
   for (iVar = 1; iVar < nVar; iVar++) {
    Res = local_Residual[iVar];
    if (iRKStep < 3) {
    /* Base Solution Update */
    nodes->AddSolution(iPoint, iVar, tmp_time*Res); 
    /* New Solution Update */
    nodes->AddSolution_New(iPoint, iVar, tmp_func*Res);
    } else {
     nodes->SetSolution(iPoint, iVar, nodes->GetSolution_New(iPoint, iVar) + tmp_func*Res);
     AddRes_RMS(iVar, Res*Res);
    }
   }
 } // end if iPoint

 /*--- MPI solution ---*/
 InitiateComms(geometry, config, SOLUTION);
 CompleteComms(geometry, config, SOLUTION);
cout<<"---------------------end iRK "<<iRKStep<<"-------------------------"<<endl;
 SetResidual_RMS(geometry, config);

 /*--- For verification cases, compute the global error metrics. ---*/   
 ComputeVerificationError(geometry, config);

} // end function

	
void CFilmSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config){
 SU2_MPI::Error(string("Method not yet implemented for THIN_FILM problems."), CURRENT_FUNCTION);
} // end function	


su2double* CFilmSolver::GetCharacPrimVar(unsigned short val_marker, unsigned long val_vertex){ 
 return CharacPrimVar[val_marker][val_vertex];
} // end function

su2double* CFilmSolver::GetCharacPrimVar(unsigned short iLayer, unsigned short val_marker, unsigned long val_vertex){
 if(iLayer == (nLayer-1))
  return GetCharacPrimVar(val_marker, val_vertex);
 else if(iLayer < (nLayer-1))
 return CharacPrimVarLayer[iLayer][val_marker][val_vertex];
 else
  return NULL;
}

su2double CFilmSolver::GetTotal_thickness(unsigned short iLayer, unsigned long iPoint){
 return total_thickness[iLayer][iPoint];
}

void CFilmSolver::BC_Fluid_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, 
                                     CNumerics *visc_numerics, CConfig *config){

 unsigned long iVertex, jVertex, iPoint, Point_Normal = 0;
 unsigned short iDim, iVar, iMarker, nDonorVertex;
 su2double weight;
 bool interaction = config->GetEnergy_Equation();
 su2double *Normal    = new su2double[nDim];
 su2double *PrimVar_i = new su2double[nPrimVar];
 su2double *PrimVar_j = new su2double[nPrimVar];

 for (iDim = 0; iDim < nDim; iDim++) 
  Normal[iDim] = 0.0;

 for(iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
  if(config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) { 
   for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

    for(iVar = 0; iVar < nVar; iVar++)
     Res_Conv[iVar] = 0.0;

    if(geometry->node[iPoint]->GetDomain()) { 
     nDonorVertex = GetnSlidingStates(iMarker, iVertex);
    /*--- Loop over the nDonorVertexes and compute the averaged flux ---*/
     for(jVertex = 0; jVertex < nDonorVertex; jVertex++) {

      for(iVar = 0; iVar < nPrimVar; iVar++) {
       PrimVar_i[iVar] = nodes->GetPrimitive(iPoint,iVar);
       PrimVar_j[iVar] = GetSlidingState(iMarker, iVertex, iVar, jVertex);
      }

     /*--- Get the weight computed in the interpolator class for the j-th donor vertex ---*/
      weight = GetSlidingState(iMarker, iVertex, nPrimVar, jVertex);

     /*--- Set the normal vector ---*/
      geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) 
       Normal[iDim] = -Normal[iDim];


      if(interaction){

       for(iDim = 0; iDim < nADim; iDim++)
        Res_Conv[1] = PrimVar_i[nADim+1]*(PrimVar_i[iDim+1] - PrimVar_j[iDim+1])*Normal[iDim];
      
       for(iVar = 2; iVar < nVar; iVar++)
        for(iDim = 0; iDim < nADim; iDim++)
         Res_Conv[iVar] = PrimVar_i[nADim+1]*PrimVar_i[iVar-1]*(PrimVar_i[iDim+1] - PrimVar_j[iDim+1])*Normal[iDim];

       for(iVar = 0; iVar < nVar; iVar++)
        Res_Conv[iVar] *= weight;

       Layer_SysRes.AddBlock(iPoint + (nLayer-1)*nPoint, Res_Conv);
      }

     }
    }
   }
  }
 }

 delete [] Normal;
 delete [] PrimVar_i;
 delete [] PrimVar_j;
} // end function


void CFilmSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
                   CConfig *config, unsigned short val_marker){ 
cout<<"-----------------WALL--------------------"<<endl;
 /*--- Call the equivalent symmetry plane boundary condition. ---*/
 BC_Sym_Plane(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

 return;
} // end function

void CFilmSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                   CConfig *config, unsigned short val_marker){

 unsigned short iDim, iVar, iLayer;
 unsigned long  iPoint, jPoint;
 unsigned long  iVertex, iEdge;
 
 su2double  Area, ProjVelocity_i, *V_reflected, *V_domain;
 su2double* Normal     = new su2double[nDim];
 su2double* UnitNormal = new su2double[nDim];

/*--- Loop over all layers ---*/
for(iLayer = 0; iLayer < nLayer; iLayer++){
 /*--- Loop over all the vertices on this boundary marker. ---*/
 for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

  iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

 /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
  if (geometry->node[iPoint]->GetDomain()) {
   
   /*--- Allocate the reflected state at the symmetry boundary. ---*/
   V_reflected = GetCharacPrimVar(iLayer, val_marker, iVertex);
       
   /*--- Get current solution at this boundary node ---*/
   V_domain = nodes->GetLayerPrimitive(iPoint, iLayer);

   if(iLayer != (nLayer -1))
    layer_vertex[iLayer][val_marker][iVertex]->GetNormal(Normal);
   else 
    geometry->vertex[val_marker][iVertex]->GetNormal(Normal); 

   for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
//cout<<"Check for point "<<iPoint<<", with normal ";for(iDim = 0; iDim < nDim; iDim++)cout<<Normal[iDim]<<" ";cout<<endl;

   /*--- Compute unit normal, to be used for unit tangential, projected velocity and velocity 
         component gradients. ---*/
   Area = 0.0;
   for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
   Area = sqrt (Area);

   for (iDim = 0; iDim < nDim; iDim++)
    if(Area != 0.0)
     UnitNormal[iDim] = Normal[iDim]/Area;
   conv_numerics->SetNormal(Normal);

   iEdge  = geometry->node[iPoint]->GetEdge(0);
   jPoint = geometry->edge[iEdge]->GetNode(1);
   jPoint = (jPoint == iPoint) ? geometry->edge[iEdge]->GetNode(0) : jPoint;

   if(iLayer != (nLayer -1) )
    conv_numerics->SetCoord(layer_node[iLayer][iPoint]->GetCoord(), layer_node[iLayer][jPoint]->GetCoord());
   else
    conv_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
        
   /*--- Set the reflected state based on the boundary node. Scalars are copied and 
         the velocity is mirrored along the symmetry boundary, i.e. the velocity in 
         normal direction is substracted twice. ---*/
   for(iVar = 0; iVar < nPrimVar; iVar++)
     V_reflected[iVar] = nodes->GetLayerPrimitive(iPoint, iVar, iLayer);
   /*--- Compute velocity in normal direction (ProjVelcity_i=(v*n)) und substract twice from
         velocity in normal direction: v_r = v - 2 (v*n)n ---*/
   ProjVelocity_i = 0.0;
   for(iDim = 0; iDim < nADim; iDim++) 
    ProjVelocity_i += V_domain[iDim+1]*UnitNormal[iDim];     
   for(iDim = 0; iDim < nADim; iDim++)
    V_reflected[iDim+1] = V_domain[iDim+1] - 2.0 * ProjVelocity_i*UnitNormal[iDim];

  // Audesse hp: prescribed thickness
   V_reflected[nADim+2] = delta_layer[iPoint];
/*
cout<<"V_domain  variables: "; for(iVar=0;iVar<nPrimVar;iVar++)cout<<V_domain[iVar]    <<" ";cout<<endl;   
cout<<"V_reflect variables: "; for(iVar=0;iVar<nPrimVar;iVar++)cout<<V_reflected[iVar] <<" ";cout<<endl;   
*/

   /*--- Set Primitive and Secondary for numerics class. ---*/
   conv_numerics->SetPrimitive(V_domain, V_reflected);

   /*--- Compute the residual using an upwind scheme. ---*/
   conv_numerics->ComputeResidual(Residual, config);
       
   /*--- Update residual value ---*/
   Layer_SysRes.AddBlock(iPoint + iLayer*nPoint, Residual);
//cout<<"Euler Wall residual: ";for(iVar=0;iVar<nVar;iVar++)cout<<Residual[iVar] <<" ";cout<<endl; 
  } // end if Domain                  
 } // end for iVertex
} // end for iLayer

 /*--- Free locally allocated memory ---*/
 delete [] Normal;
 delete [] UnitNormal;

} // end function


void CFilmSolver::SetUniformInlet(CConfig* config, unsigned short iMarker){

 unsigned short iLayer;

 if(config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) {
     
  string Marker_Tag   = config->GetMarker_All_TagBound(iMarker);
  su2double  p_total  = config->GetInlet_Ptotal(Marker_Tag);
  su2double  t_total  = config->GetInlet_Ttotal(Marker_Tag);
  su2double* flow_dir = config->GetInlet_FlowDir(Marker_Tag);
   
  for(iLayer = 0; iLayer < nLayer; iLayer++)  
   for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
    Inlet_Ttotal[iLayer][iMarker][iVertex] = t_total;
    Inlet_Ptotal[iLayer][iMarker][iVertex] = p_total;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
     Inlet_FlowDir[iLayer][iMarker][iVertex][iDim] = flow_dir[iDim];
    }     
  } else{
     
   /*--- For now, non-inlets just get set to zero. ---*/    
   for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
    Inlet_Ttotal[iLayer][iMarker][iVertex] = 0.0;
    Inlet_Ptotal[iLayer][iMarker][iVertex] = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
     Inlet_FlowDir[iLayer][iMarker][iVertex][iDim] = 0.0;
   }
  }

} // end function

void CFilmSolver::SetInletAtVertex(su2double *val_inlet, unsigned short iMarker, unsigned long iVertex){

 /*--- IMPORTANT ASSUMPTION on val_inlet: inlet input files stores the profile for the
       complete geometry using properly the nr of columnes of the inlet.dat input file, 
       starting from the bottom of the geometry. ---*/

 /*--- Alias positions within inlet file for readability ---*/  
 unsigned short T_position       = nDim;
 unsigned short P_position       = nDim+1;
 unsigned short FlowDir_position = nDim+2;

 for(unsigned short iLayer = 0; iLayer < nLayer; iLayer++){
  Inlet_Ttotal[iLayer][iMarker][iVertex] = val_inlet[T_position + (nADim+2)*iLayer];
  Inlet_Ptotal[iLayer][iMarker][iVertex] = val_inlet[P_position + (nADim+2)*iLayer];
  for (unsigned short iDim = 0; iDim < nADim; iDim++) {
   Inlet_FlowDir[iLayer][iMarker][iVertex][iDim] =  val_inlet[FlowDir_position + iDim + (nADim+2)*iLayer];
  }
 }

} // end function


void CFilmSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
              CConfig *config, unsigned short val_marker){
//cout<<"-----------------INLET--------------------"<<endl;
 unsigned short iDim, iLayer, iVar;
 unsigned long iVertex, iPoint, jPoint, iEdge, Point_Normal;
 unsigned short iOrder, jLayer;
 su2double  Area;
 su2double* V_inlet, *V_domain;
 su2double* Normal      = new su2double[nDim];
 su2double* UnitFlowDir = new su2double[nDim];
 su2double* Flow_Dir    = new su2double[nDim];
 su2double  P_domain, P_inlet, Vel_Mag, Flow_Dir_Mag;

 string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
 unsigned short Kind_Inlet = config->GetKind_Inc_Inlet(Marker_Tag);
 
 /*--- Loop over all layers ---*/   
 for(iLayer = 0; iLayer < nLayer; iLayer++){
 /*--- Loop over all the vertices on this boundary marker ---*/   
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
     
   if(iLayer != (nLayer -1))
    iPoint = layer_vertex[iLayer][val_marker][iVertex]->GetNode();
   else 
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

  /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/     
  if(geometry->node[iPoint]->GetDomain()) {

   /*--- Setting coordinates ---*/ 
   iEdge  = geometry->node[iPoint]->GetEdge(0);
   jPoint = geometry->edge[iEdge]->GetNode(1);
   jPoint = (jPoint == iPoint) ? geometry->edge[iEdge]->GetNode(0) : jPoint;
   if(iLayer != (nLayer -1) )
    conv_numerics->SetCoord(layer_node[iLayer][iPoint]->GetCoord(), layer_node[iLayer][jPoint]->GetCoord());
   else
    conv_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
     
  /*--- Normal vector for this vertex (negate for outward convention) ---*/
   if(iLayer != (nLayer -1))
    layer_vertex[iLayer][val_marker][iVertex]->GetNormal(Normal);
   else 
    geometry->vertex[val_marker][iVertex]->GetNormal(Normal);

   for(iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
   conv_numerics->SetNormal(Normal);
//cout<<"Check for point "<<iPoint<<", with coord ";for(iDim = 0; iDim < nDim; iDim++)cout<<geometry->node[iPoint]->GetCoord()[iDim]<<" ";cout<<endl;

   Area = 0.0;
   for(iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
   Area = sqrt (Area);

   /*--- Allocate the value at the inlet ---*/   
   V_inlet = GetCharacPrimVar(iLayer, val_marker, iVertex);
   /*--- Retrieve solution at this boundary node. ---*/
   V_domain = nodes->GetLayerPrimitive(iPoint, iLayer);     

   Flow_Dir = Inlet_FlowDir[iLayer][val_marker][iVertex];
   Flow_Dir_Mag = 0.0;
   for (iDim = 0; iDim < nDim; iDim++)
    Flow_Dir_Mag += Flow_Dir[iDim]*Flow_Dir[iDim];
   Flow_Dir_Mag = sqrt(Flow_Dir_Mag);

   /*--- Store the unit flow direction vector. ---*/
   for (iDim = 0; iDim < nDim; iDim++)
    UnitFlowDir[iDim] = Flow_Dir[iDim]/Flow_Dir_Mag;

   switch(Kind_Inlet){
   
   case VELOCITY_INLET:
    
    V_inlet[0] = V_domain[0];
  /*--- Store the velocity in the primitive variable vector. ---*/
    Vel_Mag  = Inlet_Ptotal[iLayer][val_marker][iVertex]/config->GetVelocity_Ref();
    for (iDim = 0; iDim < nADim; iDim++)
     V_inlet[iDim+1] = Vel_Mag*UnitFlowDir[iDim];    

    break;

   case PRESSURE_INLET:
     
    P_inlet = Inlet_Ptotal[iLayer][val_marker][iVertex]/config->GetPressure_Ref();        
   /*--- Store the current static pressure for clarity. ---*/          
    P_domain = nodes->GetPressure(iPoint);
    if (P_domain > P_inlet) {
             
    /*--- Back flow: use the prescribed P_inlet as static pressure. ---*/
     V_inlet[0] = P_inlet;
             
    /*--- Neumann condition for velocity. ---*/ 
     for (iDim = 0; iDim < nADim; iDim++)
      V_inlet[iDim+1] = V_domain[iDim+1];

    } else {
  
      V_inlet[0] = P_domain;
     /*--- Update the velocity magnitude using the total pressure. ---*/
      Vel_Mag = sqrt( 2*(P_inlet - P_domain)/(nodes->GetDensity(iPoint)) );

      for(iDim = 0; iDim < nDim; iDim++){
       if(Area != 0.0) UnitFlowDir[iDim] = -Normal[iDim]/Area;
      }

      for(iDim = 0; iDim < nADim; iDim++)
       V_inlet[iDim+1] = Vel_Mag*UnitFlowDir[iDim];    
    }

    break; 
 
   }

    V_inlet[nADim+1] = nodes->GetDensity(iPoint);
    V_inlet[nADim+2] = V_domain[nADim+2];    

   /*--- Set various quantities in the solver class ---*/ 
/*
cout<<"V_domain variables: "; for(iVar=0;iVar<nPrimVar;iVar++)cout<<V_domain[iVar]<<" ";cout<<endl;   
cout<<"V_inlet  variables: "; for(iVar=0;iVar<nPrimVar;iVar++)cout<<V_inlet[iVar] <<" ";cout<<endl;   
*/
   conv_numerics->SetPrimitive(V_domain, V_inlet);

   if (dynamic_grid){
    if(iLayer != (nLayer -1))
     conv_numerics->SetGridVel(layer_node[iLayer][iPoint]->GetGridVel(), layer_node[iLayer][iPoint]->GetGridVel());
    else
     conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
   }

   /*--- Compute the residual using an upwind scheme ---*/
   conv_numerics->ComputeResidual(Residual, config);

//cout<<"Inlet Conv residual: ";for(iVar=0;iVar<nVar;iVar++)cout<<Residual[iVar] <<" ";cout<<endl; 
   /*--- Update residual value ---*/
   Layer_SysRes.AddBlock(iPoint + iLayer*nPoint, Residual);
  } // end if Domain                  
  } // end for iVertex
 } // end for iLayer


////////////Controllo///////////////
/*for(iLayer = 0; iLayer < nLayer; iLayer++ ){
 for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++){
  iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  V_domain = nodes->GetLayerPrimitive(iPoint, iLayer);     
  Flow_Dir = Inlet_FlowDir[iLayer][val_marker][iVertex];
  Flow_Dir_Mag = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
   Flow_Dir_Mag += Flow_Dir[iDim]*Flow_Dir[iDim];
  Flow_Dir_Mag = sqrt(Flow_Dir_Mag);
  for (iDim = 0; iDim < nDim; iDim++)
   UnitFlowDir[iDim] = Flow_Dir[iDim]/Flow_Dir_Mag;
  Vel_Mag  = Inlet_Ptotal[iLayer][val_marker][iVertex]/config->GetVelocity_Ref();
  for (iDim = 0; iDim < nADim; iDim++)
   V_inlet[iDim+1] = Vel_Mag*UnitFlowDir[iDim];    
  V_inlet[nADim+1] = nodes->GetDensity(iPoint);
  V_inlet[nADim+2] = V_domain[nADim+2];    
  Layer_SysRes.SetBlock_Zero(iPoint + iLayer*nPoint);
  for(iVar=0;iVar<nPrimVar;iVar++)
   nodes->SetLayerPrimitive(iPoint,iVar,V_inlet[iVar],iLayer);
  nodes->SetLayerSolution_New(iPoint,0,V_inlet[0],iLayer);
  nodes->SetLayerSolution_New(iPoint,1,V_inlet[nADim+1]*V_inlet[nADim+2],iLayer);
  for(iDim=0;iDim<nADim;iDim++)
   nodes->SetLayerSolution_New(iPoint,iDim+2,V_inlet[iDim+1]*V_inlet[nADim+2],iLayer);
 }
}
*////////////////////////////////

 /*--- Free allocated memory --*/
 delete [] Normal;
 delete [] UnitFlowDir;

} // end function


void CFilmSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
                       CConfig *config, unsigned short val_marker) { 

/*--- Temporary call to supersonic outlet boundary condition. ---*/
 if(first_bc){
  cout << "Warning: using Supersonic_Outlet, simple Outlet not implemented yet." << endl;
  first_bc = false;
 }
 BC_Supersonic_Outlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
} // end function

void CFilmSolver::BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, 
              CConfig *config, unsigned short val_marker){ 

//cout<<"----------------OUTLET--------------------"<<endl;
 unsigned short iDim, iLayer, iVar;
 unsigned long  iPoint, jPoint;
 unsigned long  iVertex, iEdge;
 unsigned short iOrder, jLayer;
 string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
 su2double *V_outlet, *V_domain;
 su2double  P_Outlet;
 su2double *Normal = new su2double[nDim];

 /*--- Supersonic outlet flow: there are no ingoing characteristics,
       so all flow variables can should be interpolated from the domain. ---*/
 
 for(iLayer = 0; iLayer < nLayer; iLayer++ ){  
 /*--- Loop over all the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
   
   iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
     
   /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/     
   if (geometry->node[iPoint]->GetDomain()) {
       
   /*--- Current solution at this boundary node ---*/       
    V_domain = nodes->GetLayerPrimitive(iPoint, iLayer);
       
   /*--- Allocate the value at the outlet ---*/
    V_outlet = GetCharacPrimVar(iLayer, val_marker, iVertex);

   /*--- Primitive variables, using the derived quantities ---*/       
    P_Outlet    = config->GetOutlet_Pressure(Marker_Tag)/config->GetPressure_Ref();       
    V_outlet[0] = P_Outlet/su2double(nLayer);  

    for (iDim = 0; iDim < nADim; iDim++)
     V_outlet[iDim+1] = V_domain[iDim+1];

    V_outlet[nADim+1] = V_domain[nADim+1];
    V_outlet[nADim+2] = V_domain[nADim+2];
 
   /*--- Normal vector for this vertex (negate for outward convention) ---*/ 
    if(iLayer != (nLayer-1))    
     layer_vertex[iLayer][val_marker][iVertex]->GetNormal(Normal);
    else
     geometry->vertex[val_marker][iVertex]->GetNormal(Normal);

    for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
    conv_numerics->SetNormal(Normal);
//cout<<"Check for point "<<iPoint<<", with coord ";for(iDim = 0; iDim < nDim; iDim++)cout<<geometry->node[iPoint]->GetCoord()[iDim]<<" ";cout<<endl;
      
   /*--- Set various quantities in the solver class ---*/      
    iEdge  = geometry->node[iPoint]->GetEdge(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    jPoint = (jPoint == iPoint) ? geometry->edge[iEdge]->GetNode(0) : jPoint;
    if(iLayer != (nLayer-1) )
     conv_numerics->SetCoord(layer_node[iLayer][iPoint]->GetCoord(), layer_node[iLayer][jPoint]->GetCoord());
    else
     conv_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
/*
cout << "V_domain variables: "; for(iVar=0;iVar<nPrimVar;iVar++)cout<<V_domain[iVar]<<" ";cout<<endl;   
cout << "V_outlet variables: "; for(iVar=0;iVar<nPrimVar;iVar++)cout<<V_outlet[iVar]<<" ";cout<<endl;   
*/
     conv_numerics->SetPrimitive(V_domain, V_outlet);
       
    if (dynamic_grid){
     if(iLayer != (nLayer -1))
      conv_numerics->SetGridVel(layer_node[iLayer][iPoint]->GetGridVel(), layer_node[iLayer][iPoint]->GetGridVel());
     else
      conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
    }

   /*--- Compute the residual using an upwind scheme ---*/      
    conv_numerics->ComputeResidual(Residual, config);
    Layer_SysRes.AddBlock(iPoint + iLayer*nPoint, Residual);
//cout<<"Outlet Conv residual: ";for(iVar=0;iVar<nVar;iVar++)cout<<Residual[iVar] <<" ";cout<<endl; 

   } // end if(iPoint)
  } // end for
 } // end for iLayer


////////////Controllo///////////////
/*for(iLayer = 0; iLayer < nLayer; iLayer++ ){
 for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++){
  iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  V_domain = nodes->GetLayerPrimitive(iPoint, iLayer);       
  V_outlet = GetCharacPrimVar(iLayer, val_marker, iVertex);
  P_Outlet    = config->GetOutlet_Pressure(Marker_Tag)/config->GetPressure_Ref();       
  V_outlet[0] = P_Outlet/su2double(nLayer);  
  for (iDim = 0; iDim < nADim; iDim++)
   V_outlet[iDim+1] = V_domain[iDim+1];
  V_outlet[nADim+1] = V_domain[nADim+1];
  V_outlet[nADim+2] = V_domain[nADim+2];
  Layer_SysRes.SetBlock_Zero(iPoint + iLayer*nPoint);
  for(iVar=0;iVar<nPrimVar;iVar++)
   nodes->SetLayerPrimitive(iPoint,iVar,V_outlet[iVar],iLayer);
  nodes->SetLayerSolution_New(iPoint,0,V_outlet[0],iLayer);
  nodes->SetLayerSolution_New(iPoint,1,V_outlet[nADim+1]*V_outlet[nADim+2],iLayer);
  for(iDim=0;iDim<nADim;iDim++)
   nodes->SetLayerSolution_New(iPoint,iDim+2,V_outlet[iDim+1]*V_outlet[nADim+2],iLayer);
 }
}
*////////////////////////////////


 /*--- Free allocated memory --*/
 delete [] Normal;

} // end function





