#include "MGWFCFStartTime.hh"

using namespace std;

MGWFCFStartTime::MGWFCFStartTime(const string& aName) :
  MGVWaveformParameterCalculator(aName), fRatio(0.5), fOffset(150)
{
  AddParameter("CFt0");
}

void MGWFCFStartTime::CalculateParameters(const MGWaveform& anInput)
{
  size_t offSet = static_cast<size_t>(fOffset*anInput.GetSamplingFrequency());

  // Create the delayed + inverted pulse and add it to the main waveform
  // Scale down by some fraction and invert, start with the fOffset'th sample (shifts scaled pulse backwards in time)
  // Vector length is reduced by fOffset samples
  for(size_t i = offSet; i < anInput.GetLength(); i++) fScaledInput.push_back(-fRatio*anInput[i]);

  // Sum together the inverted and delayed waveform with the original
  for(size_t i = 0; i < anInput.GetLength()-offSet; i++) fSummedVector.push_back( anInput[i] + fScaledInput[i] );
  
  // Start at offset sample and then walk foward until reaching 0 crossing
  size_t iRef = offSet;
  while(anInput[iRef] < 0.) iRef++;
  SetParameterValue(0, anInput.InterpolateForTLocal(0., iRef-1, iRef));

}
