// Class for the parameters of the fit
// Should include changing the initial parameters, etc

#ifndef __TBkgModelParameter__
#define __TBkgModelParameter__
#include "TObject.h"

class TBkgModelParameter : public TObject {

public:
	TBkgModelParameter();

	// TBkgModelParameter(bool bFixParameter);

	virtual ~TBkgModelParameter();

private:

	bool bFixed;

	char 	 dParName;
	double	 dInitalValue;
	double	 dMinLimit;
	double	 dMaxLimit;


	ClassDef(TBkgModelParameter, 1)
};

#endif // __TBkgModelParameter__
