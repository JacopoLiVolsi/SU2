/*!
 * \file Utilities.hpp
 * \brief Useful free functions
 * \author J. Li Volsi
*/

#pragma once

#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<map>
#include<string>
#include<cstring>
#include<cmath>
#include <cstdlib>
#include <iomanip>

#include "../../Common/include/mpi_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/option_structure.hpp"
#include"Point.hpp"

// __USING__
using namespace Geometry;
using std::endl;
using std::size_t;
using std::stoul;

//! \brief Function to split words/data separated by a separator
void split_string(std::string const long_str, char separator, std::vector<std::string>& vec_str);

//! \brief A function to split words/data separated by "\space"
/*! This function does not consider in any case the spaces
* that can be find in a row and it stores in a vector any character
* that is not a space but in between of at least 2 spaces.
*/
void split_string_by_space(std::string const long_str, std::vector<std::string>& vec_str);


/*! \class Mesh_Extractor
 *  \brief Extrac a mesh from the interface boundary.
 *  \ingroup Thin_Film_Equations 
 *  \author J. Li Volsi
*/
class Mesh_Extractor{

private:
 std::string config_name;
 std::string mesh_name;
 std::string out_mesh_name;
 int Ndime;
 std::string interface_tag;
 std::string new_interface;
 std::map<unsigned long,Point> map_idx_point;
 std::map<unsigned long,unsigned long> old_to_new_idx;

public:

//! \brief Constructor of the class.
 Mesh_Extractor(std::string configname);

//! \brief Preliminary read to the configuration file.
 bool preliminary_read();

inline std::string Get_mesh_name() const { return mesh_name;}
inline std::string Get_out_mesh_name() const { return out_mesh_name;}
inline std::string Get_interface_tag() const { return interface_tag;}

//! Function to extract the mesh on a given surface tag in su2 format
 bool get_surf_mesh_su2();

//! Function to set BC for a thin film mesh in su2 format
 void set_boundary_coditions_su2();

//! Function to write the free surface interface
 void write_fluid_interface();

//! Function to set BC for a 1D thin film problem
 void set_bc_1D();

//! Function to set BC for a 2D thin film problem
 void set_bc_2D();

//! Function to lead the process.
 void Excavate();

};


