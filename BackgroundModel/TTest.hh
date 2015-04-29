#ifndef __TTest__
#define __TTest__
#include "TObject.h"
#include "TBackgroundModel.hh"

class TTest : public TBackgroundModel {

public:
	TTest();

	TTest(double fFitMin, double fFitMax, int fBinBase, int fDataset, bool fSave);

	virtual ~TTest();

	double GetChiSquareM1();

	double GetChiSquareM2();



	ClassDef(TTest, 1)
};

#endif // __TTest__