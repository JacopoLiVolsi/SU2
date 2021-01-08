#pragma once
 
#include "COutput.hpp"
#include "../variables/CVariable.hpp"
#include "../../../Common/include/geometry_structure.hpp"
#include "../solver_structure.hpp"

/*!
 *  \class CFilmOutput
 *  \brief Output class for thin film problem.
 *  \ingroup Thin_Film_Equations
 *  \author J. Li Volsi
*/
class CFilmOutput : public COutput {
protected:
 
 unsigned short nADim, nLayer;
 unsigned long  nPoint;
 bool initialize_bottom;
 su2double** Bottom_Topography;
 string* LayerVolumeFilename;

public:
 
/*!
 *  \brief Constructor of the class.
*/ 
 CFilmOutput(CConfig *config, unsigned short nDim, bool femOutput);

/*!
 *  \brief Destructor of the class.
*/ 
 ~CFilmOutput() override;

/*!
 *  \brief Initialize bottom topography.
*/ 
 void Initialize_Bottom(unsigned long nPoint, CSolver *solver, CGeometry *geometry, CConfig *config);

/*!
  * \brief Determines if the screen header should be written.
  * \param[in] config - Definition of the particular problem.
  */
 bool WriteScreen_Header(CConfig *config) override;

/*!
 *  \brief Sets and load the output fields.
*/  
 void LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint) override;

/*!
 *  \brief Sets the history output fields.
*/  
 void SetHistoryOutputFields(CConfig *config) override;

/*!
 *  \brief Sets the volume output fields.
*/  
 void SetVolumeOutputFields(CConfig *config) override;

/*!
 *  \brief Sets the history output fields.
*/  
 void LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) override;

/*!
 * \brief Write any additional files defined for the current solver.
 * \param[in] config - Definition of the particular problem per zone.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] solver_container - The container holding all solution data.
 */
 void WriteAdditionalFiles(CConfig *config, CGeometry* geometry, CSolver** solver_container) override;


}; // end CFilmOutput class








