//---------------------------------------------------------------------------//
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                         MAJORANA Simulation                               //
//                                                                           //
//      This code implementation is the intellectual property of the         //
//      MAJORANA Collaboration. It is based on Geant4, an intellectual       //
//      property of the RD44 GEANT4 collaboration.                           //
//                                                                           //
//                        *********************                              //
//                                                                           //
//    Neither the authors of this software system, nor their employing       //
//    institutes, nor the agencies providing financial support for this      //
//    work  make  any representation or  warranty, express or implied,       //
//    regarding this software system or assume any liability for its use.    //
//    By copying, distributing or modifying the Program (or any work based   //
//    on on the Program) you indicate your acceptance of this statement,     //
//    and all its terms.                                                     //
//                                                                           //
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//---------------------------------------------------------------------------//
/**                                                            
 *      
 * CLASS DECLARATION:  MGWFCFStartTime.hh
 *
 * DESCRIPTION: 
 *
 * Estimates the start time (t0) of a trapezoidal filter using Constant Fraction method
 *
 * AUTHOR: B. Zhu
 * FIRST SUBMISSION: 
 */

#ifndef _MGWFCFStartTime_HH
#define _MGWFCFStartTime_HH

#include "MGVWaveformParameterCalculator.hh"

class MGWFCFStartTime : public MGVWaveformParameterCalculator
{
  public:
    MGWFCFStartTime(const std::string& aName= "MGWFCFStartTime");
    virtual ~MGWFCFStartTime() {}

    /// Set delay offset
    virtual inline void SetOffset(double offset) {fOffset = offset;}

    /// Set ratio for evaluation
    virtual inline void SetRatio(double ratio) {fRatio = ratio;}

    virtual inline void RestrictToRegion(size_t beg, size_t end) { fRegion.SetBeginning(beg); fRegion.SetEnd(end); }

    virtual void CalculateParameters(const MGWaveform& inputWF);

  protected:
    double fRatio; // Ratio to evaluate fraction
    double fOffset; // Time delay in units of samples
    std::vector<double> fScaledInput; // Internal vector for inversion
    std::vector<double> fSummedVector; // Internal vector for summing

    MGWaveformRegion fRegion;
    MGWaveform fWF;
};

#endif
