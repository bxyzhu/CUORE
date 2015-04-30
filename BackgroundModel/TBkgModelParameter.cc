#include "TBkgModelParameter.hh"

ClassImp(TBkgModelParameter)

TBkgModelParameter::TBkgModelParameter()
{
}

TBkgModelParameter::TBkgModelParameter(const char *fParName, int dParIndex, double fInitialValue, double fMinLimit, double fMaxLimit, TH1D *fHist)
{
	dParName = fParName;
	dParIndex = fParIndex;
	dInitialValue = fInitialValue;
	dMinLimit = fMinLimit;
	dMaxLimit = fMaxLimit;	
	dHist = fHist;
}

TBkgModelParameter::~TBkgModelParameter()
{}

int TBkgModelParameter::GetParIndex()
{
	return dParIndex;
}

const char *TBkgModelParameter::GetParName()
{
	return dParName;
}

double TBkgModelParameter::GetParInitial()
{
	return dInitialValue;
}

double TBkgModelParameter::GetParMin()
{
	return dParMin;
}

double TBkgModelParameter::GetParMax()
{
	return dParMax;
}

TH1D *TBkgModelParameter::GetHist()
{
	return dHist;
}

void TBkgModelParameter::SetInitValue(double fInitialValue)
{
	dInitialValue = fInitialValue;
}
