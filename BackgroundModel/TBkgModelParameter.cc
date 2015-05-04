#include "TBkgModelParameter.hh"

// using namespace std;

ClassImp(TBkgModelParameter)

TBkgModelParameter::TBkgModelParameter()
{
}

TBkgModelParameter::TBkgModelParameter(const char *fParName, int fParIndex, double fInitialValue, double fInitErr, double fMinLimit, double fMaxLimit, TH1D *&fHistM1, TH1D *&fHistM2, TH1D *&fHistM2Sum)
{
	dParName = fParName;
	dParIndex = fParIndex;
	dInitialValue = fInitialValue;
	dInitErr = fInitErr; 
	dMinLimit = fMinLimit;
	dMaxLimit = fMaxLimit;	
	dHistM1 = fHistM1;
	dHistM2 = fHistM2;
	dHistM2Sum = fHistM2Sum;
}

TBkgModelParameter::TBkgModelParameter(const char *fParName, int fParIndex, double fInitialValue, double fInitErr, double fMinLimit, double fMaxLimit, TH1D *&fHist)
{
	dParName = fParName;
	dParIndex = fParIndex;
	dInitialValue = fInitialValue;
	dInitErr = fInitErr; 
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

double TBkgModelParameter::GetParInitErr()
{
	return dInitErr;
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

TH1D *TBkgModelParameter::GetHistM1()
{
	return dHistM1;
}

TH1D *TBkgModelParameter::GetHistM2()
{
	return dHistM2;
}

TH1D *TBkgModelParameter::GetHistM2Sum()
{
	return dHistM2Sum;
}

void TBkgModelParameter::SetInitValue(double fInitialValue)
{
	dInitialValue = fInitialValue;
}
