// Class for the parameters of the fit
// Includes name, parameter number, initial value, etc
// Inherits only from TObject, not with other classes (maybe this is stupid) ...

#ifndef __TBkgModelParameter__
#define __TBkgModelParameter__
#include "TObject.h"
#include "TH1D.h"

class TBkgModelParameter : public TObject {

public:
	TBkgModelParameter();

	virtual ~TBkgModelParameter();

	int GetParIndex();

	const char *GetParName();

	double GetParInitial();

	double GetParMin();

	double GetParMax();

	TH1D *GetHist();

	void SetInitValue();

protected:
	bool bFixed;
	const char 	 	*dParName;
	int 	 		dParIndex;
	double	 		dInitalValue;
	double	 		dMinLimit;
	double	 		dMaxLimit;
	TH1D 			*fHist

	ClassDef(TBkgModelParameter, 1)
};

#endif // __TBkgModelParameter__
