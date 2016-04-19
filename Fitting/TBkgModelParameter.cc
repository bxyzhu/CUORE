#include "TBkgModelParameter.hh"

// using namespace std;

ClassImp(TBkgModelParameter)

TBkgModelParameter::TBkgModelParameter()
{
}

TBkgModelParameter::TBkgModelParameter(const char *fParName, int fParIndex, double fInitialValue, double fInitErr, double fMinLimit, double fMaxLimit, TH1D *&fHistM1Adap, TH1D *&fHistM2Adap, TH1D *&fHistM1, TH1D *&fHistM2)
{
	dParName = fParName;
	dParIndex = fParIndex;
	dInitialValue = fInitialValue;
	dInitErr = fInitErr; 
	dMinLimit = fMinLimit;
	dMaxLimit = fMaxLimit;	
	dHistM1Adap = fHistM1Adap;
	dHistM2Adap = fHistM2Adap;
	dHistM1 = fHistM1;
	dHistM2 = fHistM2;	
}

TBkgModelParameter::TBkgModelParameter(const char *fParName, int fParIndex, double fInitialValue, double fInitErr, double fMinLimit, double fMaxLimit, TH1D *&fHistM1Adap, TH1D *&fHistM2Adap)
{
	dParName = fParName;
	dParIndex = fParIndex;
	dInitialValue = fInitialValue;
	dInitErr = fInitErr; 
	dMinLimit = fMinLimit;
	dMaxLimit = fMaxLimit;	
	dHistM1Adap = fHistM1Adap;
	dHistM2Adap = fHistM2Adap;
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

TH1D *TBkgModelParameter::GetHistM1(int fIndex)
{
	if(dHistM1Adap != NULL && fIndex == 1) return dHistM1Adap;
	else if(dHistM1 != NULL && fIndex == 2) return dHistM1;
	else return NULL;
}

TH1D *TBkgModelParameter::GetHistM2(int fIndex)
{
	if(dHistM2Adap != NULL && fIndex == 1) return dHistM2Adap;
	else if(dHistM2 != NULL && fIndex == 2) return dHistM2;
	else return NULL;
}


void TBkgModelParameter::SetInitValue(double fInitialValue)
{
	dInitialValue = fInitialValue;
}
