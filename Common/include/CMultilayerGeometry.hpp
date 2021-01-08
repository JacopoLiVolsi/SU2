/*!
 * \file CMultilayerGeometry.hpp
 * \brief Headers of the main subroutines for storing multilayer geometry.
 *        The subroutines and functions are in the <i>CMultilayerGeometry.cpp</i>
 * \author J. Li Volsi
 * \version 7.0.0 "Blackbird"
*/

#pragma once
#include "geometry_structure.hpp"

/*! 
*  \class CMultilayerGeometry 
*  \brief Base class for defining the multilayer geometry. 
*  \ingroup Thin_Film_Equations  
*  \author J. Li Volsi
*/
class CMultilayerGeometry{

protected:
 unsigned short nDim;
 su2double* Normal;
 su2double  Volume;

public:

 CMultilayerGeometry();

 CMultilayerGeometry(unsigned short ndim);
 
 virtual ~CMultilayerGeometry() = 0;


 inline void SetVolume(su2double vol) { Volume = vol;}

 inline su2double GetVolume(void) { return Volume;}

 inline void SetNormal(su2double *norm) { for(unsigned short iDim = 0; iDim < nDim; iDim++) Normal[iDim]=norm[iDim];}
 
 inline su2double* GetNormal(void) { return Normal;}
};


/*! 
*  \class CMultilayerPoint 
*  \ingroup Thin_Film_Equations  
*  \author J. Li Volsi
*/
class CMultilayerPoint: public CMultilayerGeometry{

private:
 unsigned long iPoint; 
 su2double *Coord;
 su2double *GridVel;

public:

 CMultilayerPoint(su2double coord0, su2double coord1,  unsigned long ipoint);

 CMultilayerPoint(su2double coord0, su2double coord1, su2double coord2,  unsigned long ipoint);
 
 ~CMultilayerPoint();

 inline void SetCoord(su2double *coord) { for(unsigned short iDim = 0; iDim < nDim; iDim++) Coord[iDim] = coord[iDim]; }

 inline void SetCoord(unsigned short iDim, su2double val_coord) { Coord[iDim] = val_coord; }

 inline su2double* GetCoord(void) { return Coord;}

 inline su2double GetCoord(unsigned short iDim) { return Coord[iDim];}

 inline void SetGridVel(su2double *gridvel) { for(unsigned short iDim = 0; iDim < nDim; iDim++) GridVel[iDim] = gridvel[iDim]; }

 inline void SetGridVel(unsigned short iDim, su2double gridvel) { GridVel[iDim] = gridvel; }

 inline su2double* GetGridVel(void) { return GridVel;}

};


/*! 
*  \class CMultilayerEdge 
*  \ingroup Thin_Film_Equations  
*  \author J. Li Volsi
*/
class CMultilayerEdge: public CMultilayerGeometry{

private:
 unsigned long Node[2];

public:

 CMultilayerEdge(unsigned long node0, unsigned long node1, unsigned short ndim);
 
 ~CMultilayerEdge();

 unsigned long GetNode(unsigned short idx);

 inline void GetNormal(su2double *norm) { for(unsigned short iDim=0; iDim<nDim; iDim++) norm[iDim] = Normal[iDim];}

};


/*! 
*  \class CMultilayerVertex 
*  \ingroup Thin_Film_Equations  
*  \author J. Li Volsi
*/
class CMultilayerVertex: public CMultilayerGeometry{

private:
 unsigned long Node;

public:

 CMultilayerVertex(unsigned long node, unsigned short ndim);
 
 ~CMultilayerVertex();

 inline unsigned long GetNode(void) {return Node;}

 inline void GetNormal(su2double *norm) { for(unsigned short iDim=0; iDim<nDim; iDim++) norm[iDim] = Normal[iDim];}
};





