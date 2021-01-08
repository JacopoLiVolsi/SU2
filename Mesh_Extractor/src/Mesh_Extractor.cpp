/*!
 * \file Utilities.hpp
 * \brief Useful free functions
 * \author J. Li Volsi
*/

// __INCLUDE.HPP__
#include "../include/Mesh_Extractor.hpp"

void split_string_by_space(std::string const long_str, std::vector<std::string>& vec_str){
  // vector assumed to be clear
  bool prev_space = true;
  size_t j = 0; // index of vec_str

  for(size_t i = 0; i < long_str.size(); i++ ){
    if(!isspace(long_str[i])){
      vec_str[j] += long_str[i];
      prev_space = false;
    }
    if(isspace(long_str[i]) && !prev_space){
      prev_space = true;
      j++;
    }
  } // end for i

} // end function split_string_by_space


void parse_options(std::string const long_str, std::vector<std::string>& vec_str){
 // vec_str assumed to be clear
 bool wrong_char = true;
 char colon    = ',';
 char Obracket = '(';
 char Cbracket = ')';
 size_t j = 0; // index of vec_str

 for(size_t i = 0; i < long_str.size(); i++ ){
  if(!isspace(long_str[i]) && long_str[i] != colon && long_str[i] != Obracket && long_str[i] != Cbracket){
   vec_str[j] += long_str[i];
   wrong_char = false;
  }
  if(isspace(long_str[i]) && long_str[i] != colon && long_str[i] != Obracket && long_str[i] != Cbracket &&  !wrong_char){
   wrong_char = true;
  j++;
  }
 } // end for i

} // end function


Mesh_Extractor::Mesh_Extractor(std::string configname): config_name(configname){
} // end function


bool Mesh_Extractor::preliminary_read(){

 bool check = false;
 std::ifstream config_file(config_name, std::ifstream::in);
 std::string line;
 std::vector<std::string> line_vec;


 while(getline(config_file, line) && config_file.good()) {
  line_vec.clear();
  line_vec.resize(line.size());
  parse_options(line, line_vec);
  if(line_vec[0].compare("MESH_FILENAME=") == 0)
   mesh_name = line_vec[1];
  if(line_vec[0].compare("MARKER_FLUID_INTERFACE=") == 0){
   interface_tag = line_vec[1];
   new_interface = line_vec[2];
   check = true;
  }
  if(line_vec[0].compare("MESH_OUT_FILENAME=") == 0)
  out_mesh_name = line_vec[1];
 }
 config_file.close();

 return check;
} // end function preliminary_read


bool Mesh_Extractor::get_surf_mesh_su2(){

/*--- Local variables ---*/
int wNdime, Nelem, Npoin;
Point null_point;
char chk_m = 'M';
char chk_n = 'N';
char chk_comment = '%';
int new_idx = 0;
bool marker_found = false;

/*  --------- READING FROM INPUT FILE ---------*/
std::ifstream mesh_file(mesh_name, std::ifstream::in);
std::string line;
/* --------- WRITING IN OUTPUT FILE ---------*/
std::ofstream surface_mesh; // creates the new surface mesh
surface_mesh.open(out_mesh_name,std::ios_base::out | std::ios_base::trunc);

while (getline(mesh_file, line) && mesh_file.good()) {
// Reading problem dimension
  if ( line.compare("NDIME= 3") == 0 ){
    Ndime = 3;
    wNdime = 3;}
  if ( line.compare("NDIME= 2") == 0 ){
    Ndime = 2;
    wNdime = 2;}
}

// Writing problem dimension
surface_mesh << "%" << endl;
surface_mesh << "% Problem dimension" << endl;
surface_mesh << "%" << endl;
surface_mesh << "NDIME= " << wNdime << endl;

// Reading number of element from the element of the free surface boundary
mesh_file.clear();
mesh_file.seekg(0);
std::string tag_searched;
tag_searched = "MARKER_TAG= " + interface_tag;

while (getline(mesh_file, line) && mesh_file.good()) {
if(line.compare(tag_searched) == 0){
    marker_found = true;
    std::string nelem_str;
    getline(mesh_file,nelem_str);/* go to the next Line */
    std::string aux_str ="MARKER_ELEMS= ";
    nelem_str.erase(0,aux_str.size()-1);/* read only the number after MARKER_ELEMS= (Oss. imp space after =) */
    Nelem = stoi(nelem_str);
    surface_mesh << "%" << endl;
    surface_mesh << "% Inner element connectivity" << endl;
    surface_mesh << "%" << endl;
    surface_mesh << "NELEM= " << Nelem << endl;
    /* Reading number of points and points doing so:
        - read all the elements correspond free-surface
        - put them in a map<key=old index,value= Point> so that they are unique
        - the size of the set is NPOIN
    */
  // Auxiliary variables
  int elem_idx = 0;

  unsigned short n_vtk = 1;
  bool vtk_set = false;

  std::string elem_str;
  std::vector<std::string> vaux_elem;// contains all point indices of the chosen MARKER_TAG

  while (getline(mesh_file,elem_str) && mesh_file.good()) {
    if(elem_str[0] == chk_m) {break;} /* if we encounter the line MARKER_TAG= other_tag, it uniquely
                                         means that the mesh for the free surface is complete */
    vaux_elem.clear();
    vaux_elem.resize(elem_str.size());
    split_string_by_space(elem_str, vaux_elem);
    surface_mesh << vaux_elem[0] << " ";
    unsigned long aux2_int;

    if(!vtk_set) {
     std::string test_vtk = vaux_elem[0];
     //-- Line element --//
     if(test_vtk == "3")
      n_vtk = 3;
     //-- Triangle element --//
     if(test_vtk == "5")
      n_vtk = 4;
     //-- Quadrilateral element --//
     if(test_vtk == "9")
      n_vtk = 5;
     if(n_vtk == 1){
      std::cerr << "Incountered error in the vtk format for the elements." << endl;
      std::cerr << "------- STOPPING EXECUTION -------" << endl;
      return false;
     }
     vtk_set = true;    
    } 

    for (size_t kk = 1; kk < n_vtk; kk++) { //vaux_elem.size()
      if(!vaux_elem[kk].empty()){
      aux2_int = stoul(vaux_elem[kk], NULL, 10);
      if(old_to_new_idx.find(aux2_int) == old_to_new_idx.end()){
        old_to_new_idx.insert(std::pair<unsigned long,unsigned long>(aux2_int,new_idx) );
        new_idx++;
      }
      surface_mesh << old_to_new_idx.at(aux2_int) << " ";
      map_idx_point.insert(std::pair<unsigned long,Point>(old_to_new_idx.at(aux2_int),null_point));
      } // if !empty
    }
    surface_mesh << elem_idx << endl;
    elem_idx++;
   } // end while elem_str
 } // end if tag_searched ok  maybe useful a break; here
} // end while line

if(!marker_found){
  std::cerr << "It does not exist a boundary called " << interface_tag << "." << endl;
  std::cerr << "------- STOPPING EXECUTION -------" << endl;
  return false;
}

Npoin = map_idx_point.size();
surface_mesh << "%" << endl;
surface_mesh << "% Node coordinate" << endl;
surface_mesh << "%" << endl;
surface_mesh << "NPOIN= " << Npoin << endl;
// Rewinding mesh_file read
mesh_file.clear();
mesh_file.seekg(0);

/* Reading physical point, but only the ones present in old_to_new_idx */
while (getline(mesh_file, line) && mesh_file.good()) {
std::string chk_npoin;
for(size_t i=0; i < line.size(); i++){ /* loop on each char of the string that extract the nr */
  if(!isspace(line[i]) )
     chk_npoin += line[i];
  else
    break;}
if(chk_npoin.compare("NPOIN=") == 0 ){
std::string poin_str;
std::vector<std::string> vaux_poin;
while (getline(mesh_file,poin_str) && mesh_file.good() ) {
  if(poin_str[0]==chk_n || poin_str[0]==chk_m || poin_str[0] == chk_comment)
    {break;}
  vaux_poin.clear(); // reset for vaux_poin
  vaux_poin.resize(Ndime+1); // contains coordinate and index
  split_string_by_space(poin_str, vaux_poin);
  unsigned long aux3_int, aux4_int;
  aux3_int = stoul(vaux_poin[Ndime]);
  if(old_to_new_idx.find(aux3_int) != old_to_new_idx.end()){
  aux4_int = old_to_new_idx.at(aux3_int);
  if (map_idx_point.find(aux4_int) != map_idx_point.end()) {
   Point aux3_point(Ndime); // IMP!!!
   for (int idx_aux_poin = 0; idx_aux_poin < Ndime; idx_aux_poin++){
   if(!vaux_poin[idx_aux_poin].empty()){
    aux3_point[idx_aux_poin] = stod(vaux_poin[idx_aux_poin]);
//      surface_mesh << aux3_point[idx_aux_poin] << " ";
      }
    } // end for idx_aux_poin
  map_idx_point.at(aux4_int) = aux3_point;
//  surface_mesh << aux4_int << endl;
} //end if find mapidx_point
} // end if find old_to_new_idx
 } // end while poin_str
 } // end if chk_npoin -  maybe a break; here
} // end while

/*-- Loop writing points --*/
for(int kk = 0; kk < Npoin; kk++){

 for(unsigned short idim = 0; idim < Ndime; idim++)
  surface_mesh << map_idx_point[kk][idim] << " ";

 surface_mesh << kk << endl;
}


mesh_file.close();
surface_mesh.close();

return true;
} // end get_surf_mesh_su2


void Mesh_Extractor::set_boundary_coditions_su2(){

 if(Ndime == 2){
   std::ofstream surface_mesh;
   surface_mesh.open(out_mesh_name, std::ios_base::out | std::ios_base::app);
   surface_mesh << "NMARK= 3" << endl;
   surface_mesh.close();
   write_fluid_interface();
   set_bc_1D();
  }
 if(Ndime == 3){
   set_bc_2D();
  }
 return;
} // end set_boundary_coditions_su2


void Mesh_Extractor::write_fluid_interface(){


std::string line;
char chk_m = 'M';
std::string tag_searched;
tag_searched = "MARKER_TAG= " + interface_tag;
int Nelem;

/*  --------- READING FROM INPUT FILE ---------*/
std::ifstream mesh_file(mesh_name, std::ifstream::in);

std::ofstream surface_mesh;
surface_mesh.open(out_mesh_name,std::ios_base::out | std::ios_base::app);

surface_mesh << "MARKER_TAG= " << new_interface << endl;

while (getline(mesh_file, line) && mesh_file.good()) {
  if(line.compare(tag_searched) == 0){
    std::string nelem_str;
    getline(mesh_file,nelem_str);/* go to the next Line */
    std::string aux_str ="MARKER_ELEMS= ";
    nelem_str.erase(0,aux_str.size()-1);/* read only the number after MARKER_ELEMS= (oss imp space after =) */
    Nelem = stoi(nelem_str);
    surface_mesh << "MARKER_ELEM= " << Nelem << endl;
    std::string elem_str;
    std::vector<std::string> vaux_elem;// contains all point indices of the chosen MARKER_TAG
    while (getline(mesh_file,elem_str) && mesh_file.good()) {
     if(elem_str[0] == chk_m) {break;}
     vaux_elem.clear();
     vaux_elem.resize(elem_str.size());
     split_string_by_space(elem_str, vaux_elem);
     surface_mesh << vaux_elem[0] << " ";
     unsigned long aux2_int;
     unsigned long aux3_int;
     for (size_t kk = 1; kk < vaux_elem.size(); kk++) {
      if(!vaux_elem[kk].empty()){
        aux2_int = stoul(vaux_elem[kk], NULL, 10);
        aux3_int = old_to_new_idx.at(aux2_int);
        surface_mesh << aux3_int << " ";
        } // if !empty
      }
      surface_mesh << endl;
     } // end while elem_str
   } // end if tag_searched ok  maybe useful a break; here
  } // end while line
 surface_mesh.close();
} // end function


void Mesh_Extractor::set_bc_1D(){
 int uno = 1;
 double x_max, x_min;
 int idx_max, idx_min;
 idx_max = 0; idx_min = 0;
 x_max = map_idx_point.at(idx_max).get_x();
 x_min = map_idx_point.at(idx_min).get_x();
for(auto it = map_idx_point.cbegin(); it != map_idx_point.cend(); it++){
  if ( it->second.get_x() > x_max ){
    x_max = it->second.get_x();
    idx_max = it->first;
  }
  if ( it->second.get_x() < x_min ){
    x_min = it->second.get_x();
    idx_min = it->first;
  }
}
std::ofstream surface_mesh;
surface_mesh.open(out_mesh_name, std::ios_base::out | std::ios_base::app);

surface_mesh << "MARKER_TAG= left" << endl;
surface_mesh << "MARKER_ELEMS= 1" << endl;
surface_mesh << uno << " " << idx_min << endl;

surface_mesh << "MARKER_TAG= right" << endl;
surface_mesh << "MARKER_ELEMS= 1" << endl;
surface_mesh << uno << " " << idx_max << endl;

std::cout << "-----EXIT SUCCESSFUL-----" << endl;   
return;
} // end set_bc_1D



void Mesh_Extractor::set_bc_2D(){
 
  /*--- CConfig object initialization --*/
  char config_file_name[MAX_STRING_SIZE];
  strcpy(config_file_name, config_name.c_str());
  CConfig *config = new CConfig(config_file_name, SU2_MSH, false);

  /*--- CGeometry object initialization --*/
  std::cout << "Creating an auxiliary geometry object." << endl;
  CGeometry *geometry;
  geometry = NULL;
  /*--- Definition of the geometry class to store the primal grid in the
     partitioning process. ---*/
  CGeometry *geometry_aux = NULL;
  geometry_aux = new CPhysicalGeometry(config, ZONE_0, 1);


  /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/
  geometry_aux->SetColorGrid_Parallel(config);

  geometry = new CPhysicalGeometry(geometry_aux, config);

  delete geometry_aux;

  geometry->SetPoint_Connectivity();
  geometry->SetElement_Connectivity();
  geometry->SetBoundVolume();

  /*--- Create the edge structure ---*/
  geometry->SetEdges();
  geometry->SetVertex(config);


  unsigned long  iPoint, jPoint, aux_Point, iElem;
  size_t *VecSize;
  unsigned short iMarker, nMarker, nCommonMarker, interface_Marker = -1;
  unsigned short iNode, iNode_Conn;
  bool *isCommonMarker;
  std::string vtk_line = "3";
  std::vector<std::string> marker_names;
  std::map<unsigned short,std::vector<unsigned long>> Marker_common_points;
  std::vector<std::pair<unsigned long, unsigned long>> *boundary_lines;

  nMarker = config->GetnMarker_All();
  nCommonMarker = 1; // surely FLUID_INTERFACE is common

  for(iMarker = 0; iMarker < nMarker; iMarker++){
   if(config->GetMarker_All_TagBound(iMarker) == interface_tag){
     interface_Marker = iMarker; 
     std::cout << "Interface marker is tagged as: " << interface_Marker << "." << endl;
   }
  }

  if(interface_Marker == -1){
   std::cout << "Error, no Marker corresponds to FLUID_INTERFACE." << endl;   
   delete config;
   delete geometry;
   std::cout << "-----EXIT UNSUCCESSFUL-----" << endl;   
   return;
  }

  marker_names.resize(nMarker);
  isCommonMarker = new bool[nMarker];
  VecSize = new size_t[nMarker];
  for(iMarker = 0; iMarker < nMarker; iMarker++){
   marker_names[iMarker] = config->GetMarker_All_TagBound(iMarker);
   isCommonMarker[iMarker] = false;
  }

 //--- Loop identifying boundary nodes ---//
  for(iMarker = 0; iMarker < nMarker; iMarker++) {
   if( iMarker != interface_Marker ){
    for(iElem = 0; iElem < geometry->nElem_Bound[iMarker]; iElem++){
     for(iNode = 0; iNode < geometry->bound[iMarker][iElem]->GetnNodes(); iNode++){
      iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
      if(old_to_new_idx.find(iPoint) != old_to_new_idx.end() ){
       if(map_idx_point.find(old_to_new_idx.at(iPoint)) != map_idx_point.end())
        Marker_common_points[iMarker].push_back(iPoint);
      }
     }
    }
    if(!Marker_common_points[iMarker].empty()){
     nCommonMarker++;
     isCommonMarker[iMarker] = true;
    }
   } // end iMarker != interface_Marker
  } // end for iMarker

 // Eliminating duplicates
  for(iMarker = 0; iMarker < nMarker; iMarker++){
   std::sort(Marker_common_points[iMarker].begin(), Marker_common_points[iMarker].end());
   auto it = std::unique(Marker_common_points[iMarker].begin(), Marker_common_points[iMarker].end());
   Marker_common_points[iMarker].erase(it, Marker_common_points[iMarker].end());
  }

  for(iMarker = 0; iMarker < nMarker; iMarker++)
   VecSize[iMarker] = Marker_common_points[iMarker].size();

  boundary_lines = new std::vector<std::pair<unsigned long, unsigned long>>[nMarker];

  for(iMarker = 0; iMarker < nMarker; iMarker++) {
   if(iMarker != interface_Marker){
    for(size_t ii_vec = 0; ii_vec < VecSize[iMarker]; ii_vec++){
     iPoint = Marker_common_points[iMarker][ii_vec];
     for(size_t jj_vec = ii_vec; jj_vec < VecSize[iMarker]; jj_vec++) { // forward check to avoid duplicates couple
      jPoint = Marker_common_points[iMarker][jj_vec];
      if(iPoint != jPoint){
       for(iElem = 0; iElem < geometry->nElem_Bound[iMarker]; iElem++){
         for(iNode = 0; iNode < geometry->bound[iMarker][iElem]->GetnNodes(); iNode++){
          if(iPoint == geometry->bound[iMarker][iElem]->GetNode(iNode) ){
           for(iNode_Conn = 0; iNode_Conn < geometry->bound[iMarker][iElem]->GetnNeighbor_Nodes(iNode); iNode_Conn++){
            aux_Point = geometry->bound[iMarker][iElem]->GetNeighbor_Nodes(iNode,iNode_Conn);
            if(jPoint == geometry->bound[iMarker][iElem]->GetNode(aux_Point) )
             boundary_lines[iMarker].push_back(std::pair<unsigned long,unsigned long>(iPoint, jPoint));
          }
         }
        }
       } // end iElem
      } // end iPoint != jPoint
     } // end for jj_vec
    } // end for ii_vec
   } // end if iMarker != interface_marker
  } // end for iMarker

 // Opening output mesh file
  std::ofstream surface_mesh;
  surface_mesh.open(out_mesh_name, std::ios_base::out | std::ios_base::app);
  surface_mesh << "NMARK= " << nCommonMarker << endl;
  surface_mesh.close();
  write_fluid_interface();
  surface_mesh.open(out_mesh_name, std::ios_base::out | std::ios_base::app);

 // Writing markers
  for(iMarker = 0; iMarker < nMarker; iMarker++) {
   if(isCommonMarker[iMarker]){
    surface_mesh << "MARKER_TAG= " << marker_names[iMarker] << endl;
    surface_mesh << "MARKER_ELEMS= " << boundary_lines[iMarker].size() << endl;
    for(auto it = boundary_lines[iMarker].begin(); it != boundary_lines[iMarker].end(); it++){
     surface_mesh << vtk_line << " " << old_to_new_idx.at(it->first)<< " "<< old_to_new_idx.at(it->second) <<endl;
    }
   } // end if isCommonMarker
  } // end for iMarker

  //--- Free allocated memory ---//
  for(iMarker = 0; iMarker < nMarker; iMarker++)
   boundary_lines[iMarker].clear();
  if( boundary_lines != NULL) delete [] boundary_lines;
  if( isCommonMarker != NULL) delete [] isCommonMarker;
  if( VecSize != NULL) delete [] VecSize;
  std::cout << "Deleted auxiliary geometrical structures." << endl;

  delete config;
  std::cout << "Deleted auxiliary CConfig object." << endl;

  delete geometry;
  std::cout << "Deleted auxiliary CGeometry object." << endl;

  std::cout << "-----EXIT SUCCESSFUL-----" << endl;   

} // end set_bc_2D


void Mesh_Extractor::Excavate(){

 bool interface_exist;
 interface_exist = preliminary_read();
 if(!interface_exist){
  std::cerr << "Interface boundary does not exist, please check "<< mesh_name <<" file." << '\n';
  return;
} else {
  std::cout << "Creating "<< Get_out_mesh_name() << " mesh file, extracting " <<
                Get_interface_tag() << " from " << Get_mesh_name() << " file."<< '\n';
}


 interface_exist = false;

 std::cout << endl;
 std::cout << "------- EXTRACTING MESH FROM " << interface_tag << " BOUNDARY -------" << endl;
 std::cout << endl;
 interface_exist = get_surf_mesh_su2();

if(interface_exist){
 std::cout << "------- SETTING BOUNDARY CONDITIONS -------" << endl;
 std::cout << endl;
 set_boundary_coditions_su2();
}

} // end function
