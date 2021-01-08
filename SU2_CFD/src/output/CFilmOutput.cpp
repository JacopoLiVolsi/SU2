
#include "../../include/output/CFilmOutput.hpp"

CFilmOutput::CFilmOutput(CConfig *config, unsigned short nDim, bool femOutput) : COutput(config,nDim,femOutput){

 nADim = nDim - 1;
 initialize_bottom = true;
 nLayer = config->GetnLayer();

 /*--- Set the volume filename --- */
 volumeFilename = config->GetVolume_FileName();
 LayerVolumeFilename = new string[nLayer-1];

 for(unsigned short iLayer = 1; iLayer < nLayer; iLayer++){
  stringstream ilay;
  ilay << iLayer;
  LayerVolumeFilename[iLayer-1] = volumeFilename + "Layer" + ilay.str();
 }
 /*--- Set the default history fields if nothing is set in the config file ---*/

 if (nRequestedHistoryFields == 0){
  requestedHistoryFields.emplace_back("ITER");
  requestedHistoryFields.emplace_back("RMS_RES");
  nRequestedHistoryFields = requestedHistoryFields.size();
 }

 if (nRequestedScreenFields == 0){
  if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
  requestedScreenFields.emplace_back("TIME_ITER");
  requestedScreenFields.emplace_back("INNER_ITER");
  requestedScreenFields.emplace_back("RMS_PRESSURE");
  requestedScreenFields.emplace_back("RMS_VELOCITY-X");
  if (nADim == 2) requestedScreenFields.emplace_back("RMS_VELOCITY-Y");
  nRequestedScreenFields = requestedScreenFields.size();
 }

 if (nRequestedVolumeFields == 0){
  requestedVolumeFields.emplace_back("COORDINATES");
  requestedVolumeFields.emplace_back("SOLUTION");
  requestedVolumeFields.emplace_back("PRIMITIVE");
  nRequestedVolumeFields = requestedVolumeFields.size();
 }


 if (convFields.empty() ) convFields.emplace_back("RMS_PRESSURE");
 
} // end function

 
CFilmOutput::~CFilmOutput() {
 if(Bottom_Topography != NULL){
  for(unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
   delete [] Bottom_Topography[iPoint];
  delete [] Bottom_Topography;
 }
} // end function


bool CFilmOutput::WriteScreen_Header(CConfig *config){

 unsigned long RestartIter = 0;
 if (config->GetRestart() && config->GetTime_Domain()){
  RestartIter = config->GetRestart_Iter();
 }
 unsigned long ScreenWrt_Freq_Inner = config->GetScreen_Wrt_Freq(2);
 unsigned long ScreenWrt_Freq_Outer = config->GetScreen_Wrt_Freq(1);
 unsigned long ScreenWrt_Freq_Time  = config->GetScreen_Wrt_Freq(0);

/*--- Header is always disabled for multizone problems unless explicitely requested --- */
 if (config->GetMultizone_Problem() && !config->GetWrt_ZoneConv()){
   return false;
 }

 if( (curInnerIter == 0) && (curOuterIter == 0) && (curTimeIter == RestartIter)){
  return true;
 }
 if( !PrintOutput(curTimeIter, ScreenWrt_Freq_Time) && !(curTimeIter == config->GetnTime_Iter() - 1) ){
  return false;
 }

/*--- If there is no inner or outer iteration, don't print header ---*/
 if(ScreenWrt_Freq_Outer == 0 && ScreenWrt_Freq_Inner == 0){
  return false;
 }

/*--- Print header if we are at the first inner iteration ---*/
 if (curInnerIter == 0){
   return true;
 }

 return false;

} // end function


void CFilmOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

 unsigned short iLayer;
 nPoint = geometry->GetnPoint();
 
 if(initialize_bottom)
  Initialize_Bottom(nPoint, solver[FLOW_SOL], geometry, config);

 CVariable* Node_Flow = solver[FLOW_SOL]->GetNodes();

 CPoint* Node_Geo  = geometry->node[iPoint];
 CMultilayerPoint** Layer_Node_Geo = new CMultilayerPoint*[(nLayer-1)];
  
 for(iLayer = 0; iLayer < (nLayer-1); iLayer++)
  Layer_Node_Geo[iLayer] = solver[FLOW_SOL]->GetLayerNode(iLayer, iPoint);
 
 /*--- Adding coordinates to volume output file ---*/

 su2double mean_coord; /* The normal coordinate is averaged because we assume the average unkowns
                        * are evaluated in the middle verical point. */

 mean_coord = (Layer_Node_Geo[0]->GetCoord(nADim) + Bottom_Topography[iPoint][nADim] ) / 2.0;
 /*--- Bottom Layer ---*/
 SetVolumeOutputValue("COORD-X", iPoint,  Layer_Node_Geo[0]->GetCoord(0));
 if(nDim == 2){
    SetVolumeOutputValue("COORD-Y", iPoint, mean_coord);
 } else {
  SetVolumeOutputValue("COORD-Y", iPoint,  Layer_Node_Geo[0]->GetCoord(1));
  SetVolumeOutputValue("COORD-Z", iPoint, mean_coord);
 }

 /*--- Kth Layers ---*/
 for(iLayer = 1; iLayer < (nLayer-1); iLayer++){
  mean_coord = (Layer_Node_Geo[iLayer-1]->GetCoord(nADim) + Layer_Node_Geo[iLayer]->GetCoord(nADim) ) / 2.0;
  SetVolumeOutputValue("COORD-X", iPoint + iLayer*nPoint,  Layer_Node_Geo[iLayer]->GetCoord(0));
  if(nDim == 2){
   SetVolumeOutputValue("COORD-Y", iPoint + iLayer*nPoint, mean_coord);
  } else {
   SetVolumeOutputValue("COORD-Y", iPoint + iLayer*nPoint,  Layer_Node_Geo[iLayer]->GetCoord(1));
   SetVolumeOutputValue("COORD-Z", iPoint + iLayer*nPoint, mean_coord);
  } 
 }

 /*--- Top Layer ---*/
 mean_coord = (Layer_Node_Geo[nLayer-2]->GetCoord(nADim) + Node_Geo->GetCoord(nADim) ) / 2.0;
 SetVolumeOutputValue("COORD-X", iPoint + (nLayer-1)*nPoint,  Node_Geo->GetCoord(0));
 if(nDim == 2){
    SetVolumeOutputValue("COORD-Y", iPoint + (nLayer-1)*nPoint, mean_coord);
 } else {
  SetVolumeOutputValue("COORD-Y", iPoint + (nLayer-1)*nPoint,  Node_Geo->GetCoord(1));
  SetVolumeOutputValue("COORD-Z", iPoint + (nLayer-1)*nPoint, mean_coord);
 }

//cout << "Pressure and velocity for iPoint " << iPoint << " "; 
/*--- Adding physical quantities to volume output file ---*/
 for(iLayer = 0; iLayer < nLayer; iLayer++){ 
  SetVolumeOutputValue("MEAN_PRESSURE",   iPoint + iLayer*nPoint, Node_Flow->GetLayerPrimitive(iPoint, 0, iLayer) );
  SetVolumeOutputValue("MEAN_VELOCITY-X", iPoint + iLayer*nPoint, Node_Flow->GetLayerPrimitive(iPoint, 1, iLayer) );
/*
if(iLayer != (nLayer-1))
cout<<Layer_Node_Geo[iLayer]->GetCoord(0) <<", "<<Layer_Node_Geo[iLayer]->GetCoord(1)<<endl;
cout<< "Primitive for Layer "<< iLayer <<" "<<Node_Flow->GetLayerPrimitive(iPoint, 0, iLayer)<<" "<<Node_Flow->GetLayerPrimitive(iPoint, 1, iLayer)<< endl;
*/
  if(nADim == 1){
    SetVolumeOutputValue("MEAN_VELOCITY-Y", iPoint + iLayer*nPoint, 0.0);
  } else {
    SetVolumeOutputValue("MEAN_VELOCITY-Y", iPoint + iLayer*nPoint, Node_Flow->GetLayerPrimitive(iPoint, 2, iLayer) );
    SetVolumeOutputValue("MEAN_VELOCITY-Z", iPoint + iLayer*nPoint, 0.0 );
  }
 }

/*--- Adding thickness to volume output file ---*/
 for(iLayer = 0; iLayer < (nLayer-1); iLayer++)
  SetVolumeOutputValue("LAYER_HEIGHT", iPoint + iLayer*nPoint,  Layer_Node_Geo[iLayer]->GetCoord(nADim));
 SetVolumeOutputValue("LAYER_HEIGHT", iPoint + (nLayer-1)*nPoint,  Node_Geo->GetCoord(nADim));

 for(iLayer = 0; iLayer < (nLayer-1); iLayer++)
  SetVolumeOutputValue("THICKNESS", iPoint + iLayer*nPoint,  solver[FLOW_SOL]->GetTotal_thickness(iLayer, iPoint));
 SetVolumeOutputValue("THICKNESS", iPoint + (nLayer-1)*nPoint,  solver[FLOW_SOL]->GetTotal_thickness((nLayer-1), iPoint));

} // end function

void CFilmOutput::SetHistoryOutputFields(CConfig *config){ 

 AddHistoryOutput("RMS_PRESSURE",   "rms[P]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the pressure.", HistoryFieldType::RESIDUAL);
 AddHistoryOutput("RMS_VELOCITY-X", "rms[U]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the velocity x-component.", HistoryFieldType::RESIDUAL);

if(nDim == 3)
 AddHistoryOutput("RMS_VELOCITY-Y", "rms[V]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the velocity y-component.", HistoryFieldType::RESIDUAL);

 AddHistoryOutput("CUR_TIME", "Cur_Time", ScreenOutputFormat::SCIENTIFIC, "TIME_DOMAIN", "Current physical time (s)");
 AddHistoryOutput("TIME_STEP", "Time_Step", ScreenOutputFormat::SCIENTIFIC, "TIME_DOMAIN", "Current time step (s)");

} // end function


void CFilmOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {
 
 CSolver* flow_solver = solver[FLOW_SOL];
 
 SetHistoryOutputValue("RMS_PRESSURE", log10(flow_solver->GetRes_RMS(0)));
 SetHistoryOutputValue("RMS_VELOCITY-X", log10(flow_solver->GetRes_RMS(1)));
 if(nDim == 3)  
  SetHistoryOutputValue("RMS_VELOCITY-Y", log10(flow_solver->GetRes_RMS(2)));

 if (config->GetTime_Domain()){
  SetHistoryOutputValue("TIME_STEP", config->GetDelta_UnstTimeND()*config->GetTime_Ref());
  if (curInnerIter == 0){
   SetHistoryOutputValue("CUR_TIME",  GetHistoryFieldValue("CUR_TIME") + GetHistoryFieldValue("TIME_STEP"));
  }
 }

} // end function

void CFilmOutput::Initialize_Bottom(unsigned long nPoint, CSolver *solver, CGeometry *geometry, CConfig *config){

 nPoint = geometry->GetnPoint();
 Bottom_Topography = new su2double*[nPoint];
 for(unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
  Bottom_Topography[iPoint] = new su2double[nDim];

 solver->ReadTopography(Bottom_Topography, geometry, config);
 initialize_bottom = false;

} // end function


void CFilmOutput::SetVolumeOutputFields(CConfig *config){

  // Grid coordinates
  AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
  AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
  if (nDim == 3)
    AddVolumeOutput("COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");

  // SOLUTION variables
  AddVolumeOutput("MEAN_PRESSURE",   "Pressure",   "SOLUTION", "Mean Pressure");
  AddVolumeOutput("MEAN_VELOCITY-X", "Velocity_x", "SOLUTION", "mean x-component of the velocity vector");
  AddVolumeOutput("MEAN_VELOCITY-Y", "Velocity_y", "SOLUTION", "mean y-component of the velocity vector");
  if (nDim == 3)
   AddVolumeOutput("MEAN_VELOCITY-Z", "Velocity_z", "SOLUTION", "mean z-component of the velocity vector");

  AddVolumeOutput("LAYER_HEIGHT", "Layer_height", "SOLUTION", "Layer height");
  AddVolumeOutput("THICKNESS", "Layer_thickness", "SOLUTION", "Layer thickness");

  // Grid velocity
  if (config->GetGrid_Movement()){
   AddVolumeOutput("GRID_VELOCITY-X", "Grid_Velocity_x", "GRID_VELOCITY", "x-component of the grid velocity vector");
   if (nDim == 3 )
    AddVolumeOutput("GRID_VELOCITY-Y", "Grid_Velocity_y", "GRID_VELOCITY", "y-component of the grid velocity vector");
  }

  // Primitive variables
//  AddVolumeOutput("PRESSURE_COEFF", "Pressure_Coefficient", "PRIMITIVE", "Pressure coefficient");
  AddVolumeOutput("DENSITY",        "Density",              "PRIMITIVE", "Density");

  //Residuals
  AddVolumeOutput("RES_PRESSURE",   "Residual_Pressure", "RESIDUAL", "Residual of the pressure");
  AddVolumeOutput("RES_VELOCITY-X", "Residual_Velocity_x", "RESIDUAL", "Residual of the x-velocity component");
  if (nDim == 3)
   AddVolumeOutput("RES_VELOCITY-Y", "Residual_Velocity_y", "RESIDUAL", "Residual of the y-velocity component");
  AddVolumeOutput("RES_TEMPERATURE", "Residual_Temperature", "RESIDUAL", "Residual of the temperature");


} // end function


void CFilmOutput::WriteAdditionalFiles(CConfig *config, CGeometry* geometry, CSolver** solver){

 unsigned short iLayer, iMarker, iVar;
 unsigned long  iPoint;
 su2double mean_coord;
 su2double time = config->GetTimeIter()*config->GetDelta_UnstTime();
 su2double Length_Ref = 0.0, Density_Ref = 0.0, Pressure_Ref = 0.0, Temperature_Ref = 0.0, Velocity_Ref = 0.0;
 Length_Ref   = config->GetLength_Ref(); 
 Density_Ref  = config->GetDensity_Ref();
 Velocity_Ref = config->GetVelocity_Ref();
 Pressure_Ref = config->GetPressure_Ref();

 nPoint = geometry->GetnPoint();
 if(initialize_bottom)
  Initialize_Bottom(nPoint, solver[FLOW_SOL], geometry, config);
 CVariable* Node_Flow = solver[FLOW_SOL]->GetNodes();
 CMultilayerPoint*** layer_node = solver[FLOW_SOL]->GetLayerNode(); 

 ofstream resul_file;
 resul_file.precision(15);
/*--- Retrieving the name of the solution file ---*/
 string config_filename = config->GetSolution_FileName();
 string resul_filename;
 resul_filename.reserve(config_filename.size());
 unsigned short ii = 0;
 char punto = '.';
 while(config_filename[ii] != punto  && ii < config_filename.size()){
  resul_filename.push_back(config_filename[ii]);
  ii++;
 }
 stringstream time_it;
 time_it << config->GetTimeIter();
 resul_filename += time_it.str();
 resul_filename = resul_filename + ".txt"; // + ".csv"

 if(rank == MASTER_NODE) {
  (*fileWritingTable) << "CSV file" << resul_filename;
 }

 resul_file.open(resul_filename,ios_base::out);
 resul_file << "----- THIN FILM RESULT ------" << endl;
 resul_file << "nPoint" << '\t' << nPoint << '\t' << "nDim" << '\t' << nDim << '\t'
            << "nLayer" << '\t' << nLayer << '\t' << "Time" << '\t' << time << endl;
 resul_file << "---- X ----"<< endl;
 for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << scientific << geometry->node[iPoint]->GetCoord(0)<<'\t'; resul_file << endl;
 iLayer = nLayer-1;
 while(iLayer != 0){
  for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << layer_node[iLayer-1][iPoint]->GetCoord(0)<<'\t'; resul_file << endl;
  iLayer--;
 }  
 resul_file << endl;
 resul_file << "-- Bottom Topography X--" << endl;
 for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << Bottom_Topography[iPoint][0]<<'\t'; resul_file << endl;
 resul_file << endl;

 resul_file << "---- Y ----"<< endl;
 for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << geometry->node[iPoint]->GetCoord(1)<<'\t'; resul_file << endl;
 iLayer = nLayer-1;
 while(iLayer != 0){
  for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << layer_node[iLayer-1][iPoint]->GetCoord(1)<<'\t'; resul_file << endl;
  iLayer--;
 } 

 if(nDim == 2){
  resul_file << "---- Y_Mean ----"<< endl;
  for(iPoint = 0; iPoint < nPoint; iPoint++){
  mean_coord = (layer_node[nLayer-2][iPoint]->GetCoord(nADim) + geometry->node[iPoint]->GetCoord(nADim) ) / 2.0;
   resul_file << mean_coord<<'\t';
  }
  resul_file << endl;
  iLayer = nLayer-1;
  while(iLayer != 0){
   for(iPoint = 0; iPoint < nPoint; iPoint++){
    if(iLayer-1 == 0)
     mean_coord = (layer_node[0][iPoint]->GetCoord(nADim) + Bottom_Topography[iPoint][nADim] ) / 2.0;
    else
     mean_coord = (layer_node[iLayer-2][iPoint]->GetCoord(nADim) + layer_node[iLayer-1][iPoint]->GetCoord(nADim) ) / 2.0;
    resul_file << mean_coord <<'\t'; 
   }
   resul_file << endl;
   iLayer--;
  } 
 }

 resul_file << endl;
 resul_file << "-- Bottom Topography Y--" << endl;
 for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << Bottom_Topography[iPoint][1]<<'\t'; resul_file << endl;
 resul_file << endl;


 if(nDim == 3){
  resul_file << "---- Z ----"<< endl;
  for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << geometry->node[iPoint]->GetCoord(2)<<'\t'; resul_file << endl;
  iLayer = nLayer-1;
  while(iLayer != 0){
   for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << layer_node[iLayer-1][iPoint]->GetCoord(2)<<'\t'; resul_file << endl;
   iLayer--;
  } 
  resul_file << endl;
  resul_file << "---- Z_Mean ----"<< endl;
  for(iPoint = 0; iPoint < nPoint; iPoint++){
   mean_coord = (layer_node[nLayer-2][iPoint]->GetCoord(nADim) + geometry->node[iPoint]->GetCoord(nADim) ) / 2.0;
   resul_file << mean_coord<<'\t';
  }
  resul_file << endl;
  iLayer = nLayer-1;
  while(iLayer != 0){
   for(iPoint = 0; iPoint < nPoint; iPoint++){
    if(iLayer-1 == 0)
     mean_coord = (layer_node[0][iPoint]->GetCoord(nADim) + Bottom_Topography[iPoint][nADim] ) / 2.0;
    else
     mean_coord = (layer_node[iLayer-2][iPoint]->GetCoord(nADim) + layer_node[iLayer-1][iPoint]->GetCoord(nADim) ) / 2.0;
    resul_file << mean_coord <<'\t'; 
   }
   resul_file << endl;
   iLayer--;
  } 
  resul_file << endl;
  resul_file << "-- Bottom Topography Z--" << endl;
   for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << Bottom_Topography[iPoint][2]<<'\t'; resul_file << endl;
  resul_file << endl;
 }

 for(iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
  if(config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){
  resul_file << "--- Interface X-Velocity ----" << endl;
  for(iPoint = 0; iPoint < nPoint; iPoint++)
   resul_file << solver[FLOW_SOL]->GetInterfaceVar(iMarker, iPoint, 1)*Velocity_Ref<<'\t';
  resul_file << endl;
 }

 resul_file << "--- Mean X-Velocity ----" << endl;
 iLayer = nLayer;
 while(iLayer != 0){
  for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << Node_Flow->GetLayerPrimitive(iPoint,1,iLayer-1)*Velocity_Ref<<'\t'; resul_file << endl;
  iLayer--;
 } 
 resul_file << endl;

 if(nDim == 3){

 for(iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
  if(config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){
   resul_file << "--- Interface Y-Velocity ----" << endl;
   for(iPoint = 0; iPoint < nPoint; iPoint++)
    resul_file << solver[FLOW_SOL]->GetInterfaceVar(iMarker, iPoint, nADim)*Velocity_Ref<<'\t';
   resul_file << endl;
  }

 resul_file << "--- Mean Y-Velocity ----" << endl;
 iLayer = nLayer;
 while(iLayer != 0){
  for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << Node_Flow->GetLayerPrimitive(iPoint,2,iLayer-1)*Velocity_Ref<<'\t'; resul_file << endl;
  iLayer--;
 }  
 resul_file << endl;
}

 resul_file << "--- Pressure ----"<< endl;
 iLayer = nLayer;
 while(iLayer != 0){
  for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << Node_Flow->GetLayerPrimitive(iPoint,0,iLayer-1)*Pressure_Ref<<'\t'; resul_file << endl;
  iLayer--;
 } 
 resul_file << endl;

 resul_file << "--- Total Thickness ----"<< endl; 
 iLayer = nLayer;
 while(iLayer != 0){
  for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << solver[FLOW_SOL]->GetTotal_thickness(iLayer-1, iPoint)<<'\t'; resul_file << endl;
  iLayer--;
 } 
 resul_file << endl;

 resul_file << "--- Layer Thickness ----"<< endl;
 iLayer = nLayer;
 while(iLayer != 0){
  for(iPoint = 0; iPoint < nPoint; iPoint++) resul_file << Node_Flow->GetLayerPrimitive(iPoint,nADim+2,iLayer-1) <<'\t'; resul_file << endl;
  iLayer--;
 } 
 resul_file << endl;
 resul_file.close();

 fileWritingTable->PrintFooter();

} // end function

