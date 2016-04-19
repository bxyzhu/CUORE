// Class for the parameters of the fit
// Includes name, parameter number, initial value, etc
// Inherits only from TObject, not with other classes (maybe this is stupid) ...

#ifndef __TBkgModelParameter__
#define __TBkgModelParameter__
#include "TObject.h"
#include "TH1D.h"
#include "TStyle.h"

class TBkgModelParameter : public TObject {

public:
	TBkgModelParameter();

	TBkgModelParameter(const char *fParName, int fParIndex, double fInitialValue, double fInitErr, double fMinLimit, double fMaxLimit, TH1D *&fHistM1Adap, TH1D *&fHistM2Adap, TH1D *&fHistM1, TH1D *&fHistM2);

	TBkgModelParameter(const char *fParName, int fParIndex, double fInitialValue, double fInitErr, double fMinLimit, double fMaxLimit, TH1D *&fHistM1, TH1D *&fHistM2);

	TBkgModelParameter(const char *fParName, int fParIndex, double fInitialValue, double fInitErr, double fMinLimit, double fMaxLimit, TH1D *&fHist);

	virtual ~TBkgModelParameter();

	int GetParIndex();

	const char *GetParName();

	double GetParInitial();

	double GetParInitErr();

	double GetParMin();

	double GetParMax();

	TH1D *GetHist();

	// 1 for Adap, 2 for normal
	TH1D *GetHistM1(int fIndex);

	TH1D *GetHistM2(int fIndex);

	void SetInitValue(double fInitialValue);

	bool 			bFixed;
	const char 	 	*dParName;
	int 	 		dParIndex;
	double	 		dInitialValue;
	double			dInitErr;
	double	 		dMinLimit;
	double	 		dMaxLimit;
	TH1D 			*dHist;
	TH1D 			*dHistM1;
	TH1D 			*dHistM2;
	TH1D 			*dHistM1Adap;
	TH1D 			*dHistM2Adap;

	ClassDef(TBkgModelParameter, 1)
};

#endif // __TBkgModelParameter__
