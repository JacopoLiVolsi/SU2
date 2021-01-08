/*!
 * \file main.cpp
 * \brief Start computation to extract a mesh for average problems.
 * \author J. Li Volsi
*/

#include<iostream>
#include<cstring>

#include "../include/Point.hpp"
#include "../include/Mesh_Extractor.hpp"

using std::cout;
using std::endl;

int main(int argc, char* argv[]){

 if(argc < 2){
  std::cerr << "No configuration file passed!" << '\n';
  return 0;
 }
 if(argc > 3){
  std::cerr << "Too many input arguments!" << '\n';
  return 0;
 }

 const char *config_name;
 if(argc == 2){
  config_name = argv[1];
  std::cout << "Extracting informations from " << config_name << " file." << '\n';

  Mesh_Extractor Archeologist(config_name);

  Archeologist.Excavate();

 }

 return 0;
}
