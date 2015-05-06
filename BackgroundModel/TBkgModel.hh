// Fitter class of the background model

#ifndef __TBkgModel__
#define __TBkgModel__
#include "TObject.h"
#include "TBkgModelSource.hh"
#include "TBkgModelParameter.hh"
#include "TFile.h"
#include "TH1.h"
#include "TH2C.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TChain.h"
#include "TCut.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "TMatrixT.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include <TMinuitMinimizer.h>
#include "TROOT.h"
#include "TStyle.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

class TBkgModel : public TBkgModelSource {

public:
	TBkgModel();

	TBkgModel(double fFitMin, double fFitMax, int fBinBase, int fDataset, bool fSave):TBkgModelSource(fFitMin, fFitMax, fBinBase, fDataset)
	{
		bSave = fSave;
  		dNParam = 28; // number of fitting parameters
  		dNumCalls = 0;
  		minuit = new TMinuit(dNParam);
  		nLoop = 0;
  		// Need to set save directory here...
		dSaveDir = "/Users/brian/Dropbox/code/BackgroundModel/";


  		// Generates parameter list
  		TBkgModel::GenerateParameters();
	}

	virtual ~TBkgModel();

	void CorrectForEfficiency();

	bool DoTheFit();

	void FixPar(int fParIndex);

	void ReleasePar(int fParIndex);

	double GetChiSquare();

	void GenerateParameters();

	void PrintParameters();

	void PrintParActivity();

	void ResetParameters();

	void SetParameters(int index, double value);
	
	int ShowNParameters();

	void UpdateModel();


  	TBkgModelParameter *BkgParM1[100];
  	TBkgModelParameter *BkgParM2[100];
  	TBkgModelParameter *BkgParM2Sum[100];

  	bool 	bFixedArray[100];

	double 	dChiSquare;

private:

	TMinuit			*minuit;

	TH1D 			*hChiSquaredProgressM1;
	TH1D 			*hChiSquaredProgressM2;
	TH1D 			*hChiSquaredProgressM2Sum;

	TDatime 		*tTime;

	TMatrixT<double> 	*mCorrMatrix;
	TPaveText 			*pPave; 

	// Cut Efficiency correction
	TH1D 			*hEfficiency;

	int 			nLoop;

	TFile *fSaveResult;
	std::string 	dSaveDir;


	bool			bFixedRes;
	bool			bAdaptiveBinning;
	bool 			bSave;

	int 			dNumCalls;
	int 			dMult;

	// Parameters

	int 				dNParam;
	int 				dNumFreeParameters;
	int 				dNDF;

	double				fParameters[139];
	double				fParError[139];
	double 				fParActivity[139]; // Integral of all events
	double 				fParActivityErr[139];
	double 				fParMass[139]; // Mass of all elements
	double 				fParSurfaceArea[139]; // Surface area of all elements
	double				fResolution[52];
	double 				fParEfficiencyM1[139]; // Efficiency of the parameters 
	double				dSecToYears;
	double				fMCEff[62];

 ClassDef(TBkgModel,1) // 
 };

#endif // __TBkgModel__
