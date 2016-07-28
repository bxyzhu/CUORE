#include "MGWFCFStartTime.hh"

using namespace std;

MGWFCFStartTime::MGWFCFStartTime(const string& aName) :
  MGVWaveformParameterCalculator(aName), fRatio(0.5), fOffset(1.5*us)
{
  AddParameter("CFt0");
}

void MGWFCFStartTime::CalculateParameters(const MGWaveform& anInput)
{
  // Convert delay from units of time to sample
  size_t iOffset = static_cast<size_t>(fOffset*anInput.GetSamplingFrequency());

  // Create the delayed + inverted pulse
  // Scale down by some fraction and invert, start with the iOffset'th sample
  for(size_t i = iOffset; i < anInput.GetLength(); i++) fScaledInput.push_back(-fRatio*anInput[i]);

  // Sum together the inverted and delayed waveform with the original, vector length is reduced by iOffset samples
  for(size_t i = 0; i < anInput.GetLength()-iOffset; i++) fSummedVector.push_back( anInput[i] + fScaledInput[i] );
  
  // Start at offset sample and then walk foward until reaching 0 crossing
  size_t iRef = iOffset;
  while(fSummedVector[iRef] < 0.) iRef++;

  // Interpolate with previous sample for t0
  SetParameterValue(0, anInput.InterpolateForTLocal(0., iRef-1, iRef));
  // SetParameterValue(0, iRef/anInput.GetSamplingFrequency());
}
