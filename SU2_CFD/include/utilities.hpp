/*!
 * \file utilities.hpp
 * \brief Useful free functions
 * \author J. Li Volsi
*/

#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__
// __INCLUDE__
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<cmath>
#include<map>
// __USING__
using namespace std;

//! \brief A function to split words/data separated by "\space"
/*! This function does not consider in any case the spaces
* that can be find in a row and it stores in a vector any character
* that is not a space but in between of at least 2 spaces.
*/
void split_string_by_space(string const long_str, vector<string>& vec_str);

//! \brief 
void execute_thin_film(string const config_filename);

//! \brief Function to extract the mesh on a given surface tag in su2 format
void get_surf_mesh_su2(const string configflow_name,const string interfacetag);

#endif
