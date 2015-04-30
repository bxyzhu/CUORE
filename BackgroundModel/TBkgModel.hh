// Fitter class of the background model

#ifndef __TBkgModel__
#define __TBkgModel__
#include "TObject.h"
#include "TBkgModelSource.hh"
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
#include "TMatrixDEigen.h"
#include "TGaxis.h"
#include <TMinuitMinimizer.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TBkgModelParameter.hh"

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

	TBkgModel(double fFitMin, double fFitMax, int fBinBase, int fDataset, bool fSave);
	virtual ~TBkgModel();

	std::vector<double> AdaptiveBinning(TH1D *h1, int dBinBase);

	TH1D* CalculateResidualsAdaptive(TH1D *h1, TH1D *h2, TH1D *hResid, int binMin, int binMax, int dMult);

	bool DoTheFitAdaptive(double f2nuValue, double fVariableValue);

	void DrawBkg();

	double GetChiSquareAdaptive();

	void GenerateParameters();

	void Initialize();

	TH1D *Kernal(TH1D *hMC, TH1D *hSMC);

	void LatexResultTable(double fValue);

	void PrintParameters();

	void PrintParActivity();

	void ProfileNLL(double fBestFitInit, double fBestFitChiSq);

	void ProfileNLL2D(double fBestFitInit, double fBestFitInit2, double fBestFitChiSq);

	void ResetParameters();

	void SanityCheck();

	void SetParameters(int index, double value);

	void SetParEfficiency();
	
	int ShowNParameters();

	void ToyFit(int fNumFits);

	void UpdateModelAdaptive();

	int 	dNParam;
	int 	dBinSize;
	int 	dBaseBinSize;
	int 	dNBins;
	int 	dBinBase;

	double	dMinEnergy;
	double 	dMaxEnergy;
	double	dFitMin;
	double	dFitMax;
	double 	dNorm;
	int 	dFitMinBinM1;
	int 	dFitMaxBinM1;
	int 	dFitMinBinM2;
	int 	dFitMaxBinM2;
	int 	dFitMinBinM2Sum;
	int 	dFitMaxBinM2Sum;	
	int 	dNumFreeParameters;
	int 	dNDF;

	double  	dDataIntegralTot;
	double 	dDataIntegralM1;
	double 	dDataIntegralM2;
	double 	dDataIntegralM2Sum;

	int 	dDataSet;
	double 	dLivetime;
	double 	dLivetimeYr;

	double 	dChiSquare;
	double 	dResidualRMSM1;
	double 	dResidualRMSM2;
	double 	dResidualRMSM2Sum;
	double 	dResidualRMSTot;


	int 	dAdaptiveBinsM1;
	int 	dAdaptiveBinsM2;
	int 	dAdaptiveBinsM2Sum;
	std::vector<double> dAdaptiveVectorM1;
	std::vector<double> dAdaptiveVectorM2;
	std::vector<double> dAdaptiveVectorM2Sum;
	double 	*dAdaptiveArrayM1;
	double 	*dAdaptiveArrayM2;
	double 	*dAdaptiveArrayM2Sum;

	std::map<std::string, int> dParMap;

protected:
	TMinuit			*minuit;

	TH1D 			*hChiSquaredProgressM1;
	TH1D 			*hChiSquaredProgressM2;
	TH1D 			*hChiSquaredProgressM2Sum;

	TDatime 		*tTime;

	// Cut Efficiency
	// TF1 			*fEfficiency;
	TH1D 			*hEfficiency;
	TH1D 			*hEfficiencyM2;
	// TH1 			*hEfficiencyM1;

	ofstream 		OutFile;
	ofstream 	 	OutPNLL;
	ofstream 		OutToy;

	int 			nLoop;
	std::vector<double> 	fInitValues;
	std::vector<double> 	fInitValues2;

	TFile *fBulkInner;
	TFile *fBulkInnerOld;
	TFile *fBulkInnerM2Sum;
	TFile *fBulkOuter;
	TFile *fBulkOuterOld;
	TFile *fBulkOuterM2Sum;
	TFile *fSurfaceCrystal;
	TFile *fSurfaceCrystalOld;
	TFile *fSurfaceOther;
	TFile *fSurfaceOtherOld;

	TFile *fFudge;

	TFile *fBulkSmeared;
	TFile *fSurfaceSmeared;

	// For accidental coincidence test
	TFile *fFileCoin;
	TFile *fFileCorrection;
	TH1D *fCorrectionM2; // Correction spectra for M2 (for accidental coincidences)
	TH1D *fCorrectionM2Tot;
	TH1D *fTotCorrection;

	TFile *fSaveResult;
	TFile *fToyData;

	std::string		dDataDir;
	std::string 	dMCDir;
	std::string 	dSaveDir;


	// Error Matrix
	// TMatrixT<double> 	*mCorrMatrix;

	bool			bFixedRes;
	bool			bAdaptiveBinning;
	bool 			bSave;
	bool 			bToyData;

	int 			dNumCalls;
	int 			dMult;
	double			dBestChiSq;

	// Parameters
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
