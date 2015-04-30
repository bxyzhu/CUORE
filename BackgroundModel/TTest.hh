#ifndef __TTest__
#define __TTest__
#include "TObject.h"
#include "TBkgModelSource.hh"
#include "TBkgModelParameter.hh"

class TTest : public TBkgModelSource {

public:
	TTest();

	// Call using base class header
	// TTest(double fFitMin, double fFitMax, int fBinBase, int fDataset);

	virtual ~TTest();

	void TestStuff();

	ClassDef(TTest, 1)
};

#endif // __TTest__