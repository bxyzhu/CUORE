#ifndef __TBackgroundModel__
#define __TBackgroundModel__
#include "TObject.h"
#include "TH1D.h"
#include "TChain.h"
#include "TCut.h"
#include "TRandom1.h"
class TBackgroundModel : public TObject {

public:

	TBackgroundModel();
	virtual ~TBackgroundModel();
  
	double GetChiSquare();

	void SetParameters(int index, double value);

	bool DoTheFit();

	void LoadData();
	void PrintParameters();
	void UpdateModel();

private:

	// Data
	TChain		*qtree;
	TCut 		base_cut;

	TH1D		*fDataHistoTot;
	TH1D		*fDataHistoM1;
	TH1D		*fDataHistoM2;

	// Model
	TH1D		*fModel50mKTh;
	TH1D		*fModelMixingTh;
	TH1D		*fModel600mKTh;
	TH1D		*fModelIVCTh;
	TH1D		*fModelOVCTh;

	// Smearing
	TRandom3	*fRandomGenerator;	
	
	double	fParameters[2] = {0}; 



//  ClassDef(TMyFitter,1) // 
    };

#endif // __TBackgroundModel__
