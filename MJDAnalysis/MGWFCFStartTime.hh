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
 * Simple start time (t0) estimator using constant fraction (CF) timing method:
 * The premise of CF timing is to always evaluate the t0 at a fixed time
 * after the leading edge of the pulse has reached a constant fraction of the peak pulse amplitude
 *
 * The original signal is first inverted and multiplied by a fraction (representing the evaluation point)
 * The original signal is also delayed (typically by a time larger than the rise-time) and the two signals
 * are summed together; the zero point crossing represents the time at which the pulse reaches the fraction. 
 *
 * As the rise-time and shape of the waveform can vary greatly, trigger walk is reduced by 
 * reducing the delay time (also known as amplitude and rise-time compensated timing)
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
    double fRatio; // Ratio of the delayed waveform
    double fOffset; // Time delay
    std::vector<double> fScaledInput; // Internal vector for inverted + scaled waveform
    std::vector<double> fSummedVector; // Internal vector for summing and evaluating t0

    MGWaveformRegion fRegion;
    MGWaveform fWF;
};

#endif
