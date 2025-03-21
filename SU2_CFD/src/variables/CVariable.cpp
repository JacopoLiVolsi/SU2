/*!
 * \file CVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */


#include "../../include/variables/CVariable.hpp"

CVariable::CVariable(unsigned long ndim){
 nDim = ndim;
}

CVariable::CVariable(unsigned long npoint, unsigned long nvar, CConfig *config) {

  /*--- Initialize the number of solution variables. This version
   of the constructor will be used primarily for converting the
   restart files into solution files (SU2_SOL). ---*/
  nPoint = npoint;
  nVar = nvar;

  /*--- Allocate the solution array. ---*/
  Solution.resize(nPoint,nVar) = su2double(0.0);

  if (config->GetMultizone_Problem())
    Solution_BGS_k.resize(nPoint,nVar) = su2double(0.0);

}

CVariable::CVariable(unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config) {

  /*--- Initializate the number of dimension and number of variables ---*/
  nPoint = npoint;
  nDim = ndim;
  nVar = nvar;

  /*--- Allocate fields common to all problems. Do not allocate fields
   that are specific to one solver, i.e. not common, in this class. ---*/
  Solution.resize(nPoint,nVar) = su2double(0.0);

  Solution_Old.resize(nPoint,nVar) = su2double(0.0);

  if (config->GetTime_Marching() != NO) {
    Solution_time_n.resize(nPoint, nVar);
    Solution_time_n1.resize(nPoint,nVar);
  }
  else if (config->GetTime_Domain()) {
    Solution_time_n.resize(nPoint,nVar) = su2double(0.0);
  }

  if (config->GetFSI_Simulation() && config->GetDiscrete_Adjoint()) {
    Solution_Adj_Old.resize(nPoint,nVar);
  }

  if(config->GetMultizone_Problem() && config->GetAD_Mode()) {
    AD_InputIndex.resize(nPoint,nVar)  = -1;
    AD_OutputIndex.resize(nPoint,nVar) = -1;
  }

  if (config->GetMultizone_Problem())
    Solution_BGS_k.resize(nPoint,nVar) = su2double(0.0);
}

void CVariable::Set_OldSolution() { Solution_Old = Solution; }

void CVariable::Set_Solution() { Solution = Solution_Old; }

void CVariable::Set_Solution_time_n() { Solution_time_n = Solution; }

void CVariable::Set_Solution_time_n1() { Solution_time_n1 = Solution_time_n; }

void CVariable::Set_BGSSolution_k() { Solution_BGS_k = Solution; }

void CVariable::Restore_BGSSolution_k() { Solution = Solution_BGS_k; }

void CVariable::SetResidualSumZero() { Residual_Sum.setConstant(0.0); }

void CVariable::SetAuxVarGradientZero() { Grad_AuxVar.setConstant(0.0); }

void CVariable::SetGradientZero() { Gradient.storage.setConstant(0.0); }

void CVariable::SetRmatrixZero() { Rmatrix.storage.setConstant(0.0); }

void CVariable::SetUnd_LaplZero() { Undivided_Laplacian.setConstant(0.0); }

void CVariable::SetExternalZero() { External.setConstant(0.0); }

void CVariable::RegisterSolution(bool input, bool push_index) {
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    for(unsigned long iVar=0; iVar<nVar; ++iVar) {
      if(input) {
        if(push_index) {
          AD::RegisterInput(Solution(iPoint,iVar));
        }
        else {
          AD::RegisterInput(Solution(iPoint,iVar), false);
          AD::SetIndex(AD_InputIndex(iPoint,iVar), Solution(iPoint,iVar));
        }
      }
      else {
        AD::RegisterOutput(Solution(iPoint,iVar));
        if(!push_index)
          AD::SetIndex(AD_OutputIndex(iPoint,iVar), Solution(iPoint,iVar));
      }
    }
  }
}

void CVariable::RegisterSolution_time_n() {
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint)
    for(unsigned long iVar=0; iVar<nVar; ++iVar)
      AD::RegisterInput(Solution_time_n(iPoint,iVar));
}

void CVariable::RegisterSolution_time_n1() {
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint)
    for(unsigned long iVar=0; iVar<nVar; ++iVar)
      AD::RegisterInput(Solution_time_n1(iPoint,iVar));
}
