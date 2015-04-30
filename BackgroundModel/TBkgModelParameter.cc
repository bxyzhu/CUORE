#include "TBkgModelParameter.hh"

// using namespace std;

ClassImp(TBkgModelParameter)

TBkgModelParameter::TBkgModelParameter()
{
}

TBkgModelParameter::TBkgModelParameter(const char *fParName, int fParIndex, double fInitialValue, double fMinLimit, double fMaxLimit, TH1D *fHist)
{
	// if(fHist == NULL) std::cout << Form(" \"%s\" histogram not created", fParName) << std::endl;
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
	return dMinLimit;
}

double TBkgModelParameter::GetParMax()
{
	return dMaxLimit;
}

TH1D *TBkgModelParameter::GetHist()
{
	return dHist;
}

void TBkgModelParameter::SetInitValue(double fInitialValue)
{
	dInitialValue = fInitialValue;
}
