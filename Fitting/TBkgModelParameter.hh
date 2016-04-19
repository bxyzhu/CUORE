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

	// Constructor that includes the normal binned histograms
	TBkgModelParameter(const char *fParName, int fParIndex, double fInitialValue, double fInitErr, double fMinLimit, double fMaxLimit, TH1D *&fHistM1Adap, TH1D *&fHistM2Adap, TH1D *&fHistM1, TH1D *&fHistM2);

	// Constructor with only adaptive binned histogrmas
	TBkgModelParameter(const char *fParName, int fParIndex, double fInitialValue, double fInitErr, double fMinLimit, double fMaxLimit, TH1D *&fHistM1Adap, TH1D *&fHistM2Adap);

	virtual ~TBkgModelParameter();

	// Returns parameter index
	int GetParIndex();

	// Returns parameter name
	const char *GetParName();

	// Returns initial normalization value
	double GetParInitial();

	// Returns initial error on normalization value
	double GetParInitErr();

	// Returns parameter normalization minimum
	double GetParMin();

	// Returns parameter normalization maximum
	double GetParMax();

	// Returns histogram
	// 1 for Adap, 2 for normal
	TH1D *GetHistM1(int fIndex);
	TH1D *GetHistM2(int fIndex);

	// Sets initial normalization value for parameter
	void SetInitValue(double fInitialValue);

private:
	const char 	 	*dParName;
	int 	 		dParIndex;
	double	 		dInitialValue;
	double			dInitErr;
	double	 		dMinLimit;
	double	 		dMaxLimit;
	TH1D 			*dHistM1;
	TH1D 			*dHistM2;
	TH1D 			*dHistM1Adap;
	TH1D 			*dHistM2Adap;

	ClassDef(TBkgModelParameter, 1)
};

#endif // __TBkgModelParameter__
