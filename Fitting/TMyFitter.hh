#ifndef __TMyFitter__
#define __TMyFitter__
#include "TObject.h"
#include "TH1D.h"

class TMyFitter : public TObject {

public:

	TMyFitter();
	virtual ~TMyFitter();
  
	double GetChiSquare();

	void SetParameters(int index, double value);

	bool DoTheFit();

	TH1D*	fDataHisto;
	TH1D* 	fModelHisto;

	void FillData();
void PrintParameters();
	void UpdateModel();

protected:
  
  double fParameters[2]; 


//  ClassDef(TMyFitter,1) // 
    };

#endif // __TMyFitter__
