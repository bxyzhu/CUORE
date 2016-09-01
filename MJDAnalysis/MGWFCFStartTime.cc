#include "MGWFCFStartTime.hh"

using namespace std;

MGWFCFStartTime::MGWFCFStartTime(const string& aName) :
  MGVWaveformParameterCalculator(aName), fThreshold(0)
{
  AddParameter("CFt0");
}

void MGWFCFStartTime::CalculateParameters(const MGWaveform& anInput)
{
  // Create the inverted pulse
  // Scale down by some fraction and invert
  for(size_t i = 0; i < anInput.GetLength(); i++) fScaledInput.push_back(-fRatio*anInput[i]);

  // Sum together the inverted with the original waveform delayed by fOffset, vector length is reduced by fOffset samples
  for(size_t i = fOffset; i < anInput.GetLength(); i++) fSummedVector.push_back( anInput[i - fOffset] + fScaledInput[i] );

  double threshold = fThreshold;

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
