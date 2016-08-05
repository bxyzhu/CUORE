#include "MGWFCFStartTime.hh"

using namespace std;

MGWFCFStartTime::MGWFCFStartTime(const string& aName) :
  MGVWaveformParameterCalculator(aName), fThreshold(0)
{
  AddParameter("CFt0");
}

void MGWFCFStartTime::CalculateParameters(const MGWaveform& anInput)
{
  // Create the delayed + inverted pulse
  // Scale down by some fraction and invert
  for(size_t i = 0; i < anInput.GetLength(); i++) fScaledInput.push_back(-fRatio*anInput[i]);

  // Sum together the inverted and delayed waveform with the original, vector length is reduced by iOffset samples
  for(size_t i = fOffset; i < anInput.GetLength(); i++) fSummedVector.push_back( anInput[i - fOffset] + fScaledInput[i] );

  double threshold = fThreshold;
  
  // Start at offset sample and then walk foward until reaching 0 crossing
  // size_t iRef = fOffset;
  // while(fSummedVector[iRef] < 0.) iRef++;

  // Interpolate with previous sample for t0
  // SetParameterValue(0, anInput.InterpolateForTLocal(0., iRef-1, iRef));

  // Start at maximum of trapezoidal filter
  fExtFinder.FindExtremum(anInput);
  size_t iRef = fExtFinder.GetTheExtremumPoint();
  
  // cout << "Extremum point: " << iRef << endl; 
  size_t iStart = fRegion.GetBeginning();
  size_t iStop = fRegion.GetEnd();
  if(iStop == 0 || iStop > anInput.GetLength()) iStop = anInput.GetLength();


  // Start at maximum of trapezoidal filter and walk backwards until crossing 0
  if(iRef >= iStop) iRef = iStop-1;
  while(iRef > iStart && fSummedVector[iRef] > threshold) iRef--;
  
  // Interpolate to find zero crossing (if within good region)
  if(iRef > iStart) SetParameterValue(0, anInput.InterpolateForTLocal(threshold, iRef, iRef+1));
  else SetParameterValue(0, anInput.GetTLocal(iRef));
}
