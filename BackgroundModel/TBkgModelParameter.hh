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

	TBkgModelParameter(const char *fParName, int fParIndex, double fInitialValue, double fMinLimit, double fMaxLimit, TH1D *&fHist);

	virtual ~TBkgModelParameter();

	int GetParIndex();

	const char *GetParName();

	double GetParInitial();

	double GetParMin();

	double GetParMax();

	TH1D *GetHist();

	void SetInitValue(double fInitialValue);

	bool 			bFixed;
	const char 	 	*dParName;
	int 	 		dParIndex;
	double	 		dInitialValue;
	double	 		dMinLimit;
	double	 		dMaxLimit;
	TH1D 			*dHist;

	ClassDef(TBkgModelParameter, 1)
};

#endif // __TBkgModelParameter__
