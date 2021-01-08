/*!
 * \file CMultilayerGeometry.cpp
 * \brief Source of multilayer geometry.
 * \author J. Li Volsi
 * \version 7.0.0 "Blackbird"
*/

#include "../include/CMultilayerGeometry.hpp"

CMultilayerGeometry::CMultilayerGeometry(){
 Normal = NULL;
}

CMultilayerGeometry::CMultilayerGeometry(unsigned short ndim): nDim(ndim) {

 Normal = new su2double[nDim];
 Volume = 0.0;
}

CMultilayerGeometry::~CMultilayerGeometry(){};

CMultilayerPoint::CMultilayerPoint(su2double coord0, su2double coord1, unsigned long ipoint): CMultilayerGeometry(2), iPoint(ipoint){

 nDim    = 2;
 Coord   = new su2double[nDim];
 GridVel = new su2double[nDim];

 Coord[0] = coord0;
 Coord[1] = coord1;
 
 for(unsigned short iDim = 0; iDim < nDim; iDim++)
  GridVel[iDim] = 0.0;
};


CMultilayerPoint::CMultilayerPoint(su2double coord0, su2double coord1, su2double coord2, unsigned long ipoint): CMultilayerGeometry(3), iPoint(ipoint){

 nDim    = 3;
 Coord   = new su2double[nDim];
 GridVel = new su2double[nDim];

 Coord[0] = coord0;
 Coord[1] = coord1;
 Coord[2] = coord2;
 
 for(unsigned short iDim = 0; iDim < nDim; iDim++)
  GridVel[iDim] = 0.0;
};


CMultilayerPoint::~CMultilayerPoint(){
 
 delete [] Normal;
 delete [] Coord;
 delete [] GridVel; 

}


CMultilayerEdge::CMultilayerEdge(unsigned long node0, unsigned long node1, unsigned short ndim): CMultilayerGeometry(ndim){

 Node[0] = node0;
 Node[1] = node1;

}
 
unsigned long CMultilayerEdge::GetNode(unsigned short idx){
 if(idx == 0 || idx == 1)
  return Node[idx];
 else
  return 0;
}


CMultilayerEdge::~CMultilayerEdge(){ 
 delete [] Normal;
}

CMultilayerVertex::CMultilayerVertex(unsigned long node, unsigned short ndim): CMultilayerGeometry(ndim), Node(node) {};
 
CMultilayerVertex::~CMultilayerVertex(){
 delete [] Normal;
}



