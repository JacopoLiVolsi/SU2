/*!
 * \file CLINE1.cpp
 * \brief Definition of the fictitious 2D line element with one Gauss point.
 * \author J. Li Volsi
 * \version 7.0.0 "Blackbird"
*/

#include "../../../include/geometry/elements/CElement.hpp"

CLINE1::CLINE1() : CElementWithKnownSizes<NGAUSS,NNODE,NDIM>() {

  /*--- Gauss coordinates and weights ---*/

  GaussCoord[0][0] = 0.0;  GaussCoord[0][1] = 1.0/2.0;  GaussWeight(0) = 1;

  /*--- Store the values of the shape functions and their derivatives ---*/

  unsigned short iGauss;
  su2double Xi, Eta, val_Ni;

  for (iGauss = 0; iGauss < NGAUSS; iGauss++) {

    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];

    val_Ni = Xi;        GaussPoint[iGauss].SetNi(val_Ni,0);
    val_Ni = Eta;       GaussPoint[iGauss].SetNi(val_Ni,1);
    val_Ni = 1-Xi-Eta;  GaussPoint[iGauss].SetNi(val_Ni,2);

    /*--- dN/d xi, dN/d eta ---*/

    dNiXj[iGauss][0][0] =  1.0;  dNiXj[iGauss][0][1] =  0.0;
    dNiXj[iGauss][1][0] =  0.0;  dNiXj[iGauss][1][1] =  1.0;
    dNiXj[iGauss][2][0] = -1.0;  dNiXj[iGauss][2][1] = -1.0;

  }

  /*--- Shape functions evaluated at the nodes for extrapolation of the stresses at the Gaussian Points ---*/
  /*--- The stress is constant over a TRIA1 element ---*/

  NodalExtrap[0][0] = 1.0;
  NodalExtrap[1][0] = 1.0;
  NodalExtrap[2][0] = 1.0;

}

su2double CLINE1::ComputeArea(const FrameType mode) const {

  unsigned short iDim;
  su2double a[2] = {0.0,0.0};
  su2double Area = 0.0;

  /*--- Select the appropriate source for the nodal coordinates depending on the frame requested
        for the gradient computation, REFERENCE (undeformed) or CURRENT (deformed) ---*/
  const su2activematrix& Coord = (mode==REFERENCE) ? RefCoord : CurrentCoord;

  for (iDim = 0; iDim < NDIM; iDim++) {
    a[iDim] = Coord[0][iDim]-Coord[1][iDim];
  }

  Area = sqrt( fabs(a[0]*a[0]+a[1]*a[1]) );

  return Area;

}

