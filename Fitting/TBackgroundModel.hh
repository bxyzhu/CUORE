/* 
The fitter uses two classes:
TBkgModelParameter -- container for background model parameters. 
	Primarily used to keep track of the histograms of each background model parameter.

TBackgroundModel -- Class with the actual fitter, fit is performed using TMinuit. 
	If I were to re-write the code, I would probably try use TMinuit2 or RooFit.


I usually compile the two classes and then load the objects in a ROOT macro and run whatever functions I need

eg:
	.L TBkgModelParameter.cc+
	.L TBackgroundModel.cc+

Then:
	gSystem->Load("TBkgModelParameter_cc.so");
	gSystem->Load("TBackgroundModel_cc.so");
	TBackgroundModel *f1 = new TBackgroundModel(500, 7000, 50, 0, false);
	f1->DoTheFit();
	f1->PrintParameters();
*/

#ifndef __TBackgroundModel__
#define __TBackgroundModel__
#include "TObject.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2C.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TChain.h"
#include "TCut.h"
#include "TTree.h"
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
#include "TArrayD.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

class TBackgroundModel : public TObject {

public:
	TBackgroundModel();

	// Takes in arguments: Minimum energy for fit (fFitMin), maximum energy for fit (fFitMax)
	// base number of counts in each bin (fBinBase), Dataset (fDataset) -- only use 0, Save plots (fSave)
	TBackgroundModel(double fFitMin, double fFitMax, int fBinBase, int fDataset, bool fSave);
	virtual ~TBackgroundModel();

	// Binning for M1 and M2 spectra are separated for tests
	// In reality they can probably be binned the same
	std::vector<double> AdaptiveBinningM1(TH1D *h1, int dBinBase);
	std::vector<double> AdaptiveBinningM2(TH1D *h1, int dBinBase);

	// Adds Poisson constraint on ChiSquared, outputs additional ChiSq to be added to ChiSq
	// *** Not working properly right now ***
	double AddConstraint(int fParInput);

	// Calculates and outputs the normalized residual spectrum, also fills the residual distribution into fResidual
	TH1D* CalculateResiduals(TH1D *fData, TH1D *fModel, TH1D *fResidual, int binMin, int binMax, int dMult);

	// Creates the histograms for all of the MC
	// Must be used AFTER adaptive binning as the adaptive binn
	void CreateModelHistograms();

	// Main function that calls Minuit
	bool DoTheFit();

	// Calculates and returns the ChiSq
	// You can modify how the ChiSq is calculated as well as add/remove constraints
	double GetChiSquare();

	// Initializes an array of TBkgModelParameters
	void GenerateParameters();

	// Loads all of the MC histograms from file as well as rebins them
	void Initialize();

	// Loads CUORE-0 background data and applies cuts
	void LoadData();

	// Prints the parameter name, normalization, and err on normalization
	void PrintParameters();

	// Prints Parameter activity in units of Bq/kg
	void PrintParActivity();

	// ProfileNLL of parameter based off of its best fit value
	// Varies a fixed parameter around best fit value and writes fits to a ROOT file
	void ProfileNLL(int fParFixed);

	// Creates contour plot of 2nbb and fixed parameter value 
	void ProfileNLL2D(int fParFixed);

	// Resets all parameters to 0
	// Resets model histograms (only total model histograms)
	void ResetParameters();

	// Adds additional percentage of 2nbb events
	// Was used to check if fitter can evaluate additional 2nbb events without bias
	void SanityCheck();

	// Set parameter normalization
	// Changed from SetParameters as it overlaps with TF1 to avoid confusion
	void SetParValue(int index, double value);
	
	// Sets MC efficiency of each element, mass of each element, and prior (limit)
	void SetParEfficiency();

	// Used to calculate 90% C.L limit for a parameter
	// Does the fit (gets best fit ChiSq), calculates change in ChiSq in limited parameter
	//    by varying the parameter normalization (amount of varying depends on parameter and error)
	// Outputs into a text file
	void SetLimit(int fParFixed);

	// Performs Toy fits and writes to a ROOT file
	// Calculates change in ChiSq, Pull, etc
	void ToyFit(int fStart, int fStop);

	// Creates and updates the model histogram
	// Total Model PDF is first reset and then components are added to the total model as a linear combination
	void UpdateModel();

	int 	dNParam; // Number of fitting parameters

private:
	int 	dBinSize; // Base Bin size
	int 	dBaseBinSize; // Rebinning the base bin size (if I want to modify the base bin size for testing)
	int 	dNBins; // 
	int 	dBinBase;

	double 	dMass;	
	double	dMinEnergy;
	double 	dMaxEnergy;
	double	dFitMin;
	double	dFitMax;
	double 	dNorm;

	int 	dFitMinBinM1;
	int 	dFitMaxBinM1;
	int 	dFitMinBinM2;
	int 	dFitMaxBinM2;
	int 	dNumFreeParameters;
	int 	dNDF;

	double  dDataIntegralTot;
	double 	dDataIntegralM1;
	double 	dDataIntegralM2;

	int 	dDataSet;
	double 	dLivetime;
	double 	dLivetimeYr;
	double 	dExposure;

	double 	dChiSquare;
	double 	dResidualRMSM1;
	double 	dResidualRMSM2;
	double 	dResidualRMSTot;

	int 	dAdaptiveBinsM1;
	int 	dAdaptiveBinsM2;
	std::vector<double> dAdaptiveVectorM1;
	std::vector<double> dAdaptiveVectorM2;
	double 	*dAdaptiveArrayM1;
	double 	*dAdaptiveArrayM2;

  	TBkgModelParameter *BkgPar[100];
  	bool 	bFixedArray[100];


	// Data
	TChain			*qtree;
	TCut 			base_cut;

	TMinuit			*minuit;

	TH1D			*fDataHistoTot;
	TH1D			*fDataHistoM1;
	TH1D			*fDataHistoM2;

	TH1D			*fAdapDataHistoM1;
	TH1D			*fAdapDataHistoM2;


	// Updated 01-20-2015
	// Total PDFs M1
	TH1D			*fModelTotM1;
	TH1D			*fModelTotAdapM1;
	
	TH1D			*fModelTotM2;
	TH1D			*fModelTotAdapM2;
	

//////////// Residual distributions
	TGraph			*gResidualM1;
	TGraph 			*gResidualM2;

	TH1D			*hResidualDistM1;
	TH1D			*hResidualDistM2;

	TH1D			*hResidualGausM1;
	TH1D			*hResidualGausM2;

//////////////////// MC Histograms
////////// Crystal M1 and M2
	TH1D			*hTeO20nuM1;
	TH1D			*hTeO22nuM1;
	TH1D			*hTeO22nuIKM1;
	TH1D			*hTeO2co60M1;
	TH1D			*hTeO2k40M1;
	TH1D			*hTeO2pb210M1;
	TH1D			*hTeO2po210M1;
	TH1D			*hTeO2te125M1;
	TH1D			*hTeO2th232M1;
	TH1D			*hTeO2th228M1;
	TH1D			*hTeO2ra226M1;
	TH1D			*hTeO2rn222M1;
	TH1D			*hTeO2u238M1;
	TH1D			*hTeO2th230M1;
	TH1D			*hTeO2u234M1;
	TH1D 			*hTeO2sb125M1;

	TH1D			*hTeO2th232onlyM1;
	TH1D			*hTeO2ra228pb208M1;
	TH1D			*hTeO2th230onlyM1;
	TH1D 			*hTeO2u238th230M1;
	TH1D 			*hTeO2ra226pb210M1;

	TH1D			*hTeO2Sxth232onlyM1_001;
	TH1D			*hTeO2Sxra228pb208M1_001;
	TH1D			*hTeO2Sxu238th230M1_001;
	TH1D			*hTeO2Sxth230onlyM1_001;
	TH1D			*hTeO2Sxra226pb210M1_001;
	TH1D			*hTeO2Sxpb210M1_0001;

	TH1D			*hTeO2Sxth232onlyM1_01;
	TH1D			*hTeO2Sxra228pb208M1_01;
	TH1D			*hTeO2Sxu238th230M1_01;
	TH1D			*hTeO2Sxth230onlyM1_01;
	TH1D			*hTeO2Sxra226pb210M1_01;

	TH1D			*hTeO2Sxth232onlyM1_0001;
	TH1D			*hTeO2Sxra228pb208M1_0001;
	TH1D			*hTeO2Sxu238th230M1_0001;
	TH1D			*hTeO2Sxth230onlyM1_0001;
	TH1D			*hTeO2Sxra226pb210M1_0001;

	TH1D			*hTeO2Spb210M1_01;
	TH1D			*hTeO2Spo210M1_001;
	TH1D			*hTeO2Spo210M1_01;
	TH1D			*hTeO2Sth232M1_01;
	TH1D			*hTeO2Su238M1_01;
	TH1D			*hTeO2Sxpb210M1_001;
	TH1D			*hTeO2Sxpb210M1_01;
	TH1D			*hTeO2Sxpb210M1_1;
	TH1D			*hTeO2Sxpb210M1_10;
	TH1D			*hTeO2Sxpo210M1_001;
	TH1D			*hTeO2Sxpo210M1_01;
	TH1D			*hTeO2Sxpo210M1_1;
	TH1D			*hTeO2Sxth232M1_001;
	TH1D			*hTeO2Sxth232M1_01;
	TH1D			*hTeO2Sxth232M1_1;
	TH1D			*hTeO2Sxth232M1_10;
	TH1D			*hTeO2Sxu238M1_001;
	TH1D			*hTeO2Sxu238M1_01;
	TH1D			*hTeO2Sxu238M1_1;
	TH1D			*hTeO2Sxu238M1_10;

	TH1D			*hTeO2Sxu238M1_100;
	TH1D			*hTeO2Sxth232M1_100;
	TH1D			*hTeO2Sxpb210M1_100;

	TH1D			*hTeO20nuM2;
	TH1D			*hTeO22nuM2;
	TH1D			*hTeO22nuIKM2;	
	TH1D			*hTeO2co60M2;
	TH1D			*hTeO2k40M2;
	TH1D			*hTeO2pb210M2;
	TH1D			*hTeO2po210M2;
	TH1D			*hTeO2te125M2;
	TH1D			*hTeO2th232M2;
	TH1D			*hTeO2th228M2;
	TH1D			*hTeO2ra226M2;
	TH1D			*hTeO2rn222M2;
	TH1D			*hTeO2u238M2;
	TH1D			*hTeO2th230M2;
	TH1D			*hTeO2u234M2;
	TH1D 			*hTeO2sb125M2;

	TH1D			*hTeO2th232onlyM2;
	TH1D			*hTeO2ra228pb208M2;
	TH1D			*hTeO2th230onlyM2;
	TH1D 			*hTeO2u238th230M2;
	TH1D 			*hTeO2ra226pb210M2;

	TH1D			*hTeO2Sxth232onlyM2_001;
	TH1D			*hTeO2Sxra228pb208M2_001;
	TH1D			*hTeO2Sxu238th230M2_001;
	TH1D			*hTeO2Sxth230onlyM2_001;
	TH1D			*hTeO2Sxra226pb210M2_001;
	TH1D			*hTeO2Sxpb210M2_0001;

	TH1D			*hTeO2Sxth232onlyM2_01;
	TH1D			*hTeO2Sxra228pb208M2_01;
	TH1D			*hTeO2Sxu238th230M2_01;
	TH1D			*hTeO2Sxth230onlyM2_01;
	TH1D			*hTeO2Sxra226pb210M2_01;

	TH1D			*hTeO2Sxth232onlyM2_0001;
	TH1D			*hTeO2Sxra228pb208M2_0001;
	TH1D			*hTeO2Sxu238th230M2_0001;
	TH1D			*hTeO2Sxth230onlyM2_0001;
	TH1D			*hTeO2Sxra226pb210M2_0001;

	TH1D			*hTeO2Spb210M2_01;
	TH1D			*hTeO2Spo210M2_001;
	TH1D			*hTeO2Spo210M2_01;
	TH1D			*hTeO2Sth232M2_01;
	TH1D			*hTeO2Su238M2_01;
	TH1D			*hTeO2Sxpb210M2_001;
	TH1D			*hTeO2Sxpb210M2_01;
	TH1D			*hTeO2Sxpb210M2_1;
	TH1D			*hTeO2Sxpb210M2_10;
	TH1D			*hTeO2Sxpo210M2_001;
	TH1D			*hTeO2Sxpo210M2_01;
	TH1D			*hTeO2Sxpo210M2_1;
	TH1D			*hTeO2Sxth232M2_001;
	TH1D			*hTeO2Sxth232M2_01;
	TH1D			*hTeO2Sxth232M2_1;
	TH1D			*hTeO2Sxth232M2_10;
	TH1D			*hTeO2Sxu238M2_001;
	TH1D			*hTeO2Sxu238M2_01;
	TH1D			*hTeO2Sxu238M2_1;
	TH1D			*hTeO2Sxu238M2_10;

	TH1D			*hTeO2Sxu238M2_100;
	TH1D			*hTeO2Sxth232M2_100;
	TH1D			*hTeO2Sxpb210M2_100;

/////////// Frame + Box
	TH1D			*hCuBox_CuFrameco60M1;
	TH1D			*hCuBox_CuFramek40M1;
	TH1D			*hCuBox_CuFrameth232M1;
	TH1D			*hCuBox_CuFrameu238M1;

	TH1D			*hCuBox_CuFramemn54M1;
	TH1D			*hCuBox_CuFramebi207M1;

	TH1D			*hCuBox_CuFrameth232M1_10;
	TH1D			*hCuBox_CuFrameu238M1_10;
	TH1D			*hCuBox_CuFramepb210M1_10;
	TH1D			*hCuBox_CuFramepb210M1_1;
	TH1D			*hCuBox_CuFramepb210M1_01;
	TH1D			*hCuBox_CuFramepb210M1_001;

	TH1D			*hCuBox_CuFrameth232M1_1;
	TH1D			*hCuBox_CuFrameu238M1_1;
	TH1D			*hCuBox_CuFrameth232M1_01;
	TH1D			*hCuBox_CuFrameu238M1_01;
	TH1D			*hCuBox_CuFrameth232M1_001;
	TH1D			*hCuBox_CuFrameu238M1_001;

	TH1D			*hCuBox_CuFramepb210M1_100;
	TH1D			*hCuBox_CuFrameth232M1_100;
	TH1D			*hCuBox_CuFrameu238M1_100;
	TH1D			*hCuBox_CuFramepb210M1_50;
	TH1D			*hCuBox_CuFrameth232M1_50;
	TH1D			*hCuBox_CuFrameu238M1_50;
	TH1D			*hCuBox_CuFramepb210M1_5;
	TH1D			*hCuBox_CuFrameth232M1_5;
	TH1D			*hCuBox_CuFrameu238M1_5;

	TH1D			*hCuBox_CuFrameco60M2;
	TH1D			*hCuBox_CuFramek40M2;
	TH1D			*hCuBox_CuFrameth232M2;
	TH1D			*hCuBox_CuFrameu238M2;

	TH1D			*hCuBox_CuFramemn54M2;
	TH1D			*hCuBox_CuFramebi207M2;

	TH1D			*hCuBox_CuFrameth232M2_10;
	TH1D			*hCuBox_CuFrameu238M2_10;
	TH1D			*hCuBox_CuFramepb210M2_10;
	TH1D			*hCuBox_CuFramepb210M2_1;
	TH1D			*hCuBox_CuFramepb210M2_01;
	TH1D			*hCuBox_CuFramepb210M2_001;

	TH1D			*hCuBox_CuFrameth232M2_1;
	TH1D			*hCuBox_CuFrameu238M2_1;
	TH1D			*hCuBox_CuFrameth232M2_01;
	TH1D			*hCuBox_CuFrameu238M2_01;
	TH1D			*hCuBox_CuFrameth232M2_001;
	TH1D			*hCuBox_CuFrameu238M2_001;

	TH1D			*hCuBox_CuFramepb210M2_100;
	TH1D			*hCuBox_CuFrameth232M2_100;
	TH1D			*hCuBox_CuFrameu238M2_100;
	TH1D			*hCuBox_CuFramepb210M2_50;
	TH1D			*hCuBox_CuFrameth232M2_50;
	TH1D			*hCuBox_CuFrameu238M2_50;
	TH1D			*hCuBox_CuFramepb210M2_5;
	TH1D			*hCuBox_CuFrameth232M2_5;
	TH1D			*hCuBox_CuFrameu238M2_5;

///////////// 50mK M1 and M2
	TH1D			*h50mKcs137M1;
	TH1D			*h50mKcs137M2;

//////////// Internal Shields M1 and M2
	TH1D			*hInternalco60M1;
	TH1D			*hInternalk40M1;
	TH1D			*hInternalth232M1;
	TH1D			*hInternalu238M1;

	TH1D			*hInternalco60M2;
	TH1D			*hInternalk40M2;
	TH1D			*hInternalth232M2;
	TH1D			*hInternalu238M2;



//////////// (PbRom) Roman Lead M1 and M2
	TH1D			*hPbRombi207M1;
	TH1D			*hPbRomco60M1;
	TH1D			*hPbRomcs137M1;
	TH1D			*hPbRomk40M1;
	TH1D			*hPbRompb210M1;
	TH1D			*hPbRomth232M1;
	TH1D			*hPbRomu238M1;		

	TH1D			*hPbRombi207M2;
	TH1D			*hPbRomco60M2;
	TH1D			*hPbRomcs137M2;
	TH1D			*hPbRomk40M2;
	TH1D			*hPbRompb210M2;
	TH1D			*hPbRomth232M2;
	TH1D			*hPbRomu238M2;		

//////////// External Shield M1 and M2
	TH1D			*hExtPbbi210M1;
	TH1D			*hExtPbk40M1;
	TH1D			*hExtPbth232M1;
	TH1D			*hExtPbu238M1;
	TH1D			*hExtPbpb210M1;

	TH1D			*hExtPbbi210M2;
	TH1D			*hExtPbk40M2;
	TH1D			*hExtPbth232M2;
	TH1D			*hExtPbu238M2;
	TH1D			*hExtPbpb210M2;

	TH1D 			*hExtMuonM1;
	TH1D			*hCuBox_th232spotM1;
	TH1D			*hCuBox_k40spotM1;
	TH1D 			*hBotExtPb_k40spotM1;

	TH1D 			*hExtMuonM2;
	TH1D			*hCuBox_th232spotM2;
	TH1D			*hCuBox_k40spotM2;
	TH1D 			*hBotExtPb_k40spotM2;

/////////// OVC M1 and M2
	TH1D			*hOVCco60M1;
	TH1D			*hOVCk40M1;
	TH1D			*hOVCth232M1;
	TH1D			*hOVCu238M1;		
	TH1D			*hOVCbi207M1;

	TH1D			*hOVCco60M2;
	TH1D			*hOVCk40M2;
	TH1D			*hOVCth232M2;
	TH1D			*hOVCu238M2;	
	TH1D			*hOVCbi207M2;


///////////////// Adaptive binned histograms
/////////// Crystal M1 and M2
	TH1D			*hAdapTeO20nuM1;
	TH1D			*hAdapTeO22nuM1;
	TH1D			*hAdapTeO22nuIKM1;
	TH1D			*hAdapTeO2co60M1;
	TH1D			*hAdapTeO2k40M1;
	TH1D			*hAdapTeO2pb210M1;
	TH1D			*hAdapTeO2po210M1;
	TH1D			*hAdapTeO2te125M1;
	TH1D			*hAdapTeO2th232M1;
	TH1D			*hAdapTeO2th228M1;
	TH1D			*hAdapTeO2ra226M1;
	TH1D			*hAdapTeO2rn222M1;
	TH1D			*hAdapTeO2u238M1;
	TH1D			*hAdapTeO2th230M1;
	TH1D			*hAdapTeO2u234M1;
	TH1D 			*hAdapTeO2sb125M1;

	TH1D			*hAdapTeO2th232onlyM1;
	TH1D			*hAdapTeO2ra228pb208M1;
	TH1D			*hAdapTeO2th230onlyM1;
	TH1D 			*hAdapTeO2u238th230M1;
	TH1D 			*hAdapTeO2ra226pb210M1;

	TH1D			*hAdapTeO2Sxth232onlyM1_001;
	TH1D			*hAdapTeO2Sxra228pb208M1_001;
	TH1D			*hAdapTeO2Sxu238th230M1_001;
	TH1D			*hAdapTeO2Sxth230onlyM1_001;
	TH1D			*hAdapTeO2Sxra226pb210M1_001;
	TH1D			*hAdapTeO2Sxpb210M1_0001;

	TH1D			*hAdapTeO2Sxth232onlyM1_01;
	TH1D			*hAdapTeO2Sxra228pb208M1_01;
	TH1D			*hAdapTeO2Sxu238th230M1_01;
	TH1D			*hAdapTeO2Sxth230onlyM1_01;
	TH1D			*hAdapTeO2Sxra226pb210M1_01;

	TH1D			*hAdapTeO2Sxth232onlyM1_0001;
	TH1D			*hAdapTeO2Sxra228pb208M1_0001;
	TH1D			*hAdapTeO2Sxu238th230M1_0001;
	TH1D			*hAdapTeO2Sxth230onlyM1_0001;
	TH1D			*hAdapTeO2Sxra226pb210M1_0001;

	TH1D			*hAdapTeO2Spb210M1_01;
	TH1D			*hAdapTeO2Spo210M1_001;
	TH1D			*hAdapTeO2Spo210M1_01;
	TH1D			*hAdapTeO2Sth232M1_01;
	TH1D			*hAdapTeO2Su238M1_01;
	TH1D			*hAdapTeO2Sxpb210M1_001;
	TH1D			*hAdapTeO2Sxpb210M1_01;
	TH1D			*hAdapTeO2Sxpb210M1_1;
	TH1D			*hAdapTeO2Sxpb210M1_10;
	TH1D			*hAdapTeO2Sxpo210M1_001;
	TH1D			*hAdapTeO2Sxpo210M1_01;
	TH1D			*hAdapTeO2Sxpo210M1_1;
	TH1D			*hAdapTeO2Sxth232M1_001;
	TH1D			*hAdapTeO2Sxth232M1_01;
	TH1D			*hAdapTeO2Sxth232M1_1;
	TH1D			*hAdapTeO2Sxth232M1_10;
	TH1D			*hAdapTeO2Sxu238M1_001;
	TH1D			*hAdapTeO2Sxu238M1_01;
	TH1D			*hAdapTeO2Sxu238M1_1;
	TH1D			*hAdapTeO2Sxu238M1_10;

	TH1D			*hAdapTeO2Sxu238M1_100;
	TH1D			*hAdapTeO2Sxth232M1_100;
	TH1D			*hAdapTeO2Sxpb210M1_100;

	TH1D			*hAdapTeO20nuM2;
	TH1D			*hAdapTeO22nuM2;
	TH1D			*hAdapTeO22nuIKM2;	
	TH1D			*hAdapTeO2co60M2;
	TH1D			*hAdapTeO2k40M2;
	TH1D			*hAdapTeO2pb210M2;
	TH1D			*hAdapTeO2po210M2;
	TH1D			*hAdapTeO2te125M2;
	TH1D			*hAdapTeO2th232M2;
	TH1D			*hAdapTeO2th228M2;
	TH1D			*hAdapTeO2ra226M2;
	TH1D			*hAdapTeO2rn222M2;
	TH1D			*hAdapTeO2u238M2;
	TH1D			*hAdapTeO2th230M2;
	TH1D			*hAdapTeO2u234M2;
	TH1D 			*hAdapTeO2sb125M2;

	TH1D			*hAdapTeO2th232onlyM2;
	TH1D			*hAdapTeO2ra228pb208M2;
	TH1D			*hAdapTeO2th230onlyM2;
	TH1D 			*hAdapTeO2u238th230M2;
	TH1D 			*hAdapTeO2ra226pb210M2;

	TH1D			*hAdapTeO2Sxth232onlyM2_001;
	TH1D			*hAdapTeO2Sxra228pb208M2_001;
	TH1D			*hAdapTeO2Sxu238th230M2_001;
	TH1D			*hAdapTeO2Sxth230onlyM2_001;
	TH1D			*hAdapTeO2Sxra226pb210M2_001;
	TH1D			*hAdapTeO2Sxpb210M2_0001;

	TH1D			*hAdapTeO2Sxth232onlyM2_01;
	TH1D			*hAdapTeO2Sxra228pb208M2_01;
	TH1D			*hAdapTeO2Sxu238th230M2_01;
	TH1D			*hAdapTeO2Sxth230onlyM2_01;
	TH1D			*hAdapTeO2Sxra226pb210M2_01;

	TH1D			*hAdapTeO2Sxth232onlyM2_0001;
	TH1D			*hAdapTeO2Sxra228pb208M2_0001;
	TH1D			*hAdapTeO2Sxu238th230M2_0001;
	TH1D			*hAdapTeO2Sxth230onlyM2_0001;
	TH1D			*hAdapTeO2Sxra226pb210M2_0001;

	TH1D			*hAdapTeO2Spb210M2_01;
	TH1D			*hAdapTeO2Spo210M2_001;
	TH1D			*hAdapTeO2Spo210M2_01;
	TH1D			*hAdapTeO2Sth232M2_01;
	TH1D			*hAdapTeO2Su238M2_01;
	TH1D			*hAdapTeO2Sxpb210M2_001;
	TH1D			*hAdapTeO2Sxpb210M2_01;
	TH1D			*hAdapTeO2Sxpb210M2_1;
	TH1D			*hAdapTeO2Sxpb210M2_10;
	TH1D			*hAdapTeO2Sxpo210M2_001;
	TH1D			*hAdapTeO2Sxpo210M2_01;
	TH1D			*hAdapTeO2Sxpo210M2_1;
	TH1D			*hAdapTeO2Sxth232M2_001;
	TH1D			*hAdapTeO2Sxth232M2_01;
	TH1D			*hAdapTeO2Sxth232M2_1;
	TH1D			*hAdapTeO2Sxth232M2_10;
	TH1D			*hAdapTeO2Sxu238M2_001;
	TH1D			*hAdapTeO2Sxu238M2_01;
	TH1D			*hAdapTeO2Sxu238M2_1;
	TH1D			*hAdapTeO2Sxu238M2_10;

	TH1D			*hAdapTeO2Sxu238M2_100;
	TH1D			*hAdapTeO2Sxth232M2_100;
	TH1D			*hAdapTeO2Sxpb210M2_100;

////////////// Frame + Box
	TH1D			*hAdapCuBox_CuFrameco60M1;
	TH1D			*hAdapCuBox_CuFramek40M1;
	TH1D			*hAdapCuBox_CuFrameth232M1;
	TH1D			*hAdapCuBox_CuFrameu238M1;

	TH1D			*hAdapCuBox_CuFramemn54M1;
	TH1D			*hAdapCuBox_CuFramebi207M1;

	TH1D			*hAdapCuBox_CuFrameth232M1_10;
	TH1D			*hAdapCuBox_CuFrameu238M1_10;

	TH1D			*hAdapCuBox_CuFramepb210M1_10;
	TH1D			*hAdapCuBox_CuFramepb210M1_1;
	TH1D			*hAdapCuBox_CuFramepb210M1_01;
	TH1D			*hAdapCuBox_CuFramepb210M1_001;

	TH1D			*hAdapCuBox_CuFrameth232M1_1;
	TH1D			*hAdapCuBox_CuFrameu238M1_1;
	TH1D			*hAdapCuBox_CuFrameth232M1_01;
	TH1D			*hAdapCuBox_CuFrameu238M1_01;
	TH1D			*hAdapCuBox_CuFrameth232M1_001;
	TH1D			*hAdapCuBox_CuFrameu238M1_001;

	TH1D			*hAdapCuBox_CuFramepb210M1_100;
	TH1D			*hAdapCuBox_CuFrameth232M1_100;
	TH1D			*hAdapCuBox_CuFrameu238M1_100;
	TH1D			*hAdapCuBox_CuFramepb210M1_50;
	TH1D			*hAdapCuBox_CuFrameth232M1_50;
	TH1D			*hAdapCuBox_CuFrameu238M1_50;
	TH1D			*hAdapCuBox_CuFramepb210M1_5;
	TH1D			*hAdapCuBox_CuFrameth232M1_5;
	TH1D			*hAdapCuBox_CuFrameu238M1_5;

	TH1D			*hAdapCuBox_CuFrameco60M2;
	TH1D			*hAdapCuBox_CuFramek40M2;
	TH1D			*hAdapCuBox_CuFrameth232M2;
	TH1D			*hAdapCuBox_CuFrameu238M2;

	TH1D			*hAdapCuBox_CuFramemn54M2;
	TH1D			*hAdapCuBox_CuFramebi207M2;

	TH1D			*hAdapCuBox_CuFrameth232M2_10;
	TH1D			*hAdapCuBox_CuFrameu238M2_10;
	TH1D			*hAdapCuBox_CuFramepb210M2_10;
	TH1D			*hAdapCuBox_CuFramepb210M2_1;
	TH1D			*hAdapCuBox_CuFramepb210M2_01;
	TH1D			*hAdapCuBox_CuFramepb210M2_001;

	TH1D			*hAdapCuBox_CuFrameth232M2_1;
	TH1D			*hAdapCuBox_CuFrameu238M2_1;
	TH1D			*hAdapCuBox_CuFrameth232M2_01;
	TH1D			*hAdapCuBox_CuFrameu238M2_01;
	TH1D			*hAdapCuBox_CuFrameth232M2_001;
	TH1D			*hAdapCuBox_CuFrameu238M2_001;

	TH1D			*hAdapCuBox_CuFramepb210M2_100;
	TH1D			*hAdapCuBox_CuFrameth232M2_100;
	TH1D			*hAdapCuBox_CuFrameu238M2_100;
	TH1D			*hAdapCuBox_CuFramepb210M2_50;
	TH1D			*hAdapCuBox_CuFrameth232M2_50;
	TH1D			*hAdapCuBox_CuFrameu238M2_50;
	TH1D			*hAdapCuBox_CuFramepb210M2_5;
	TH1D			*hAdapCuBox_CuFrameth232M2_5;
	TH1D			*hAdapCuBox_CuFrameu238M2_5;

/////////// 50mK M1 and M2
	TH1D			*hAdap50mKcs137M1;
	TH1D			*hAdap50mKcs137M2;
	
///////// Roman Lead M1 and M2
	TH1D			*hAdapPbRombi207M1;
	TH1D			*hAdapPbRomco60M1;
	TH1D			*hAdapPbRomcs137M1;
	TH1D			*hAdapPbRomk40M1;
	TH1D			*hAdapPbRompb210M1;
	TH1D			*hAdapPbRomth232M1;
	TH1D			*hAdapPbRomu238M1;		

	TH1D			*hAdapPbRombi207M2;
	TH1D			*hAdapPbRomco60M2;
	TH1D			*hAdapPbRomcs137M2;
	TH1D			*hAdapPbRomk40M2;
	TH1D			*hAdapPbRompb210M2;
	TH1D			*hAdapPbRomth232M2;
	TH1D			*hAdapPbRomu238M2;		

////////// Internal Shields M1 and M2
	TH1D			*hAdapInternalco60M1;
	TH1D			*hAdapInternalk40M1;
	TH1D			*hAdapInternalth232M1;
	TH1D			*hAdapInternalu238M1;

	TH1D			*hAdapInternalco60M2;
	TH1D			*hAdapInternalk40M2;
	TH1D			*hAdapInternalth232M2;
	TH1D			*hAdapInternalu238M2;


///////// OVC M1 and M2
	TH1D			*hAdapOVCco60M1;
	TH1D			*hAdapOVCk40M1;
	TH1D			*hAdapOVCth232M1;
	TH1D			*hAdapOVCu238M1;		
	TH1D			*hAdapOVCbi207M1;

	TH1D			*hAdapOVCco60M2;
	TH1D			*hAdapOVCk40M2;
	TH1D			*hAdapOVCth232M2;
	TH1D			*hAdapOVCu238M2;	
	TH1D			*hAdapOVCbi207M2;

///////// External Shield M1 and M2
	TH1D			*hAdapExtPbbi210M1;
	TH1D			*hAdapExtPbk40M1;
	TH1D			*hAdapExtPbth232M1;
	TH1D			*hAdapExtPbu238M1;
	TH1D			*hAdapExtPbpb210M1;

	TH1D			*hAdapExtPbbi210M2;
	TH1D			*hAdapExtPbk40M2;
	TH1D			*hAdapExtPbth232M2;
	TH1D			*hAdapExtPbu238M2;
	TH1D			*hAdapExtPbpb210M2;

	TH1D 			*hAdapExtMuonM1;
	TH1D			*hAdapCuBox_th232spotM1;
	TH1D			*hAdapCuBox_k40spotM1;
	TH1D 			*hAdapBotExtPb_k40spotM1;

	TH1D 			*hAdapExtMuonM2;
	TH1D			*hAdapCuBox_th232spotM2;
	TH1D			*hAdapCuBox_k40spotM2;
	TH1D 			*hAdapBotExtPb_k40spotM2;


////////////// Dummy histograms
	TH1 			*hnewM1;
	TH1				*hnewM2;

////////// Crystal M1 and M2
	TH1			*hnewTeO20nuM1;
	TH1			*hnewTeO22nuM1;
	TH1			*hnewTeO22nuIKM1;
	TH1			*hnewTeO2co60M1;
	TH1			*hnewTeO2k40M1;
	TH1			*hnewTeO2pb210M1;
	TH1			*hnewTeO2po210M1;
	TH1			*hnewTeO2te125M1;
	TH1			*hnewTeO2th232M1;
	TH1			*hnewTeO2th228M1;
	TH1			*hnewTeO2ra226M1;
	TH1			*hnewTeO2rn222M1;
	TH1			*hnewTeO2u238M1;
	TH1			*hnewTeO2th230M1;
	TH1			*hnewTeO2u234M1;
	TH1 		*hnewTeO2sb125M1;

	TH1			*hnewTeO2th232onlyM1;
	TH1			*hnewTeO2ra228pb208M1;
	TH1			*hnewTeO2th230onlyM1;
	TH1 		*hnewTeO2u238th230M1;
	TH1 		*hnewTeO2ra226pb210M1;

	TH1			*hnewTeO2Sxth232onlyM1_001;
	TH1			*hnewTeO2Sxra228pb208M1_001;
	TH1			*hnewTeO2Sxu238th230M1_001;
	TH1			*hnewTeO2Sxth230onlyM1_001;
	TH1			*hnewTeO2Sxra226pb210M1_001;
	TH1			*hnewTeO2Sxpb210M1_0001;

	TH1			*hnewTeO2Sxth232onlyM1_01;
	TH1			*hnewTeO2Sxra228pb208M1_01;
	TH1			*hnewTeO2Sxu238th230M1_01;
	TH1			*hnewTeO2Sxth230onlyM1_01;
	TH1			*hnewTeO2Sxra226pb210M1_01;

	TH1			*hnewTeO2Sxth232onlyM1_0001;
	TH1			*hnewTeO2Sxra228pb208M1_0001;
	TH1			*hnewTeO2Sxu238th230M1_0001;
	TH1			*hnewTeO2Sxth230onlyM1_0001;
	TH1			*hnewTeO2Sxra226pb210M1_0001;	

	TH1			*hnewTeO2Spb210M1_01;
	TH1			*hnewTeO2Spo210M1_001;
	TH1			*hnewTeO2Spo210M1_01;
	TH1			*hnewTeO2Sth232M1_01;
	TH1			*hnewTeO2Su238M1_01;
	TH1			*hnewTeO2Sxpb210M1_001;
	TH1			*hnewTeO2Sxpb210M1_01;
	TH1			*hnewTeO2Sxpb210M1_1;
	TH1			*hnewTeO2Sxpb210M1_10;
	TH1			*hnewTeO2Sxpo210M1_001;
	TH1			*hnewTeO2Sxpo210M1_01;
	TH1			*hnewTeO2Sxpo210M1_1;
	TH1			*hnewTeO2Sxth232M1_001;
	TH1			*hnewTeO2Sxth232M1_01;
	TH1			*hnewTeO2Sxth232M1_1;
	TH1			*hnewTeO2Sxth232M1_10;
	TH1			*hnewTeO2Sxu238M1_001;
	TH1			*hnewTeO2Sxu238M1_01;
	TH1			*hnewTeO2Sxu238M1_1;
	TH1			*hnewTeO2Sxu238M1_10;

	TH1			*hnewTeO2Sxu238M1_100;
	TH1			*hnewTeO2Sxth232M1_100;
	TH1			*hnewTeO2Sxpb210M1_100;


	TH1			*hnewTeO20nuM2;
	TH1			*hnewTeO22nuM2;
	TH1			*hnewTeO22nuIKM2;	
	TH1			*hnewTeO2co60M2;
	TH1			*hnewTeO2k40M2;
	TH1			*hnewTeO2pb210M2;
	TH1			*hnewTeO2po210M2;
	TH1			*hnewTeO2te125M2;
	TH1			*hnewTeO2th232M2;
	TH1			*hnewTeO2th228M2;
	TH1			*hnewTeO2ra226M2;
	TH1			*hnewTeO2rn222M2;
	TH1			*hnewTeO2u238M2;
	TH1			*hnewTeO2th230M2;
	TH1			*hnewTeO2u234M2;
	TH1 		*hnewTeO2sb125M2;

	TH1			*hnewTeO2th232onlyM2;
	TH1			*hnewTeO2ra228pb208M2;
	TH1			*hnewTeO2th230onlyM2;
	TH1 		*hnewTeO2u238th230M2;
	TH1 		*hnewTeO2ra226pb210M2;

	TH1			*hnewTeO2Sxth232onlyM2_001;
	TH1			*hnewTeO2Sxra228pb208M2_001;
	TH1			*hnewTeO2Sxu238th230M2_001;
	TH1			*hnewTeO2Sxth230onlyM2_001;
	TH1			*hnewTeO2Sxra226pb210M2_001;
	TH1			*hnewTeO2Sxpb210M2_0001;

	TH1			*hnewTeO2Sxth232onlyM2_01;
	TH1			*hnewTeO2Sxra228pb208M2_01;
	TH1			*hnewTeO2Sxu238th230M2_01;
	TH1			*hnewTeO2Sxth230onlyM2_01;
	TH1			*hnewTeO2Sxra226pb210M2_01;

	TH1			*hnewTeO2Sxth232onlyM2_0001;
	TH1			*hnewTeO2Sxra228pb208M2_0001;
	TH1			*hnewTeO2Sxu238th230M2_0001;
	TH1			*hnewTeO2Sxth230onlyM2_0001;
	TH1			*hnewTeO2Sxra226pb210M2_0001;	

	TH1			*hnewTeO2Spb210M2_01;
	TH1			*hnewTeO2Spo210M2_001;
	TH1			*hnewTeO2Spo210M2_01;
	TH1			*hnewTeO2Sth232M2_01;
	TH1			*hnewTeO2Su238M2_01;
	TH1			*hnewTeO2Sxpb210M2_001;
	TH1			*hnewTeO2Sxpb210M2_01;
	TH1			*hnewTeO2Sxpb210M2_1;
	TH1			*hnewTeO2Sxpb210M2_10;
	TH1			*hnewTeO2Sxpo210M2_001;
	TH1			*hnewTeO2Sxpo210M2_01;
	TH1			*hnewTeO2Sxpo210M2_1;
	TH1			*hnewTeO2Sxth232M2_001;
	TH1			*hnewTeO2Sxth232M2_01;
	TH1			*hnewTeO2Sxth232M2_1;
	TH1			*hnewTeO2Sxth232M2_10;
	TH1			*hnewTeO2Sxu238M2_001;
	TH1			*hnewTeO2Sxu238M2_01;
	TH1			*hnewTeO2Sxu238M2_1;
	TH1			*hnewTeO2Sxu238M2_10;

	TH1			*hnewTeO2Sxu238M2_100;
	TH1			*hnewTeO2Sxth232M2_100;
	TH1			*hnewTeO2Sxpb210M2_100;

/////////// Frame + Box
	TH1			*hnewCuBox_CuFrameco60M1;
	TH1			*hnewCuBox_CuFramek40M1;
	TH1			*hnewCuBox_CuFrameth232M1;
	TH1			*hnewCuBox_CuFrameu238M1;

	TH1 		*hnewCuBox_CuFramemn54M1;
	TH1 		*hnewCuBox_CuFramebi207M1;

	TH1			*hnewCuBox_CuFrameth232M1_10;
	TH1			*hnewCuBox_CuFrameu238M1_10;
	TH1			*hnewCuBox_CuFramepb210M1_10;
	TH1			*hnewCuBox_CuFramepb210M1_1;
	TH1			*hnewCuBox_CuFramepb210M1_01;
	TH1			*hnewCuBox_CuFramepb210M1_001;

	TH1			*hnewCuBox_CuFrameth232M1_1;
	TH1			*hnewCuBox_CuFrameu238M1_1;
	TH1			*hnewCuBox_CuFrameth232M1_01;
	TH1			*hnewCuBox_CuFrameu238M1_01;
	TH1			*hnewCuBox_CuFrameth232M1_001;
	TH1			*hnewCuBox_CuFrameu238M1_001;

	TH1			*hnewCuBox_CuFramepb210M1_100;
	TH1			*hnewCuBox_CuFrameth232M1_100;
	TH1			*hnewCuBox_CuFrameu238M1_100;
	TH1			*hnewCuBox_CuFramepb210M1_50;
	TH1			*hnewCuBox_CuFrameth232M1_50;
	TH1			*hnewCuBox_CuFrameu238M1_50;
	TH1			*hnewCuBox_CuFramepb210M1_5;
	TH1			*hnewCuBox_CuFrameth232M1_5;
	TH1			*hnewCuBox_CuFrameu238M1_5;

	TH1			*hnewCuBox_CuFrameco60M2;
	TH1			*hnewCuBox_CuFramek40M2;
	TH1			*hnewCuBox_CuFrameth232M2;
	TH1			*hnewCuBox_CuFrameu238M2;

	TH1 		*hnewCuBox_CuFramemn54M2;
	TH1 		*hnewCuBox_CuFramebi207M2;

	TH1			*hnewCuBox_CuFrameth232M2_10;
	TH1			*hnewCuBox_CuFrameu238M2_10;
	TH1			*hnewCuBox_CuFramepb210M2_10;
	TH1			*hnewCuBox_CuFramepb210M2_1;
	TH1			*hnewCuBox_CuFramepb210M2_01;
	TH1			*hnewCuBox_CuFramepb210M2_001;

	TH1			*hnewCuBox_CuFrameth232M2_1;
	TH1			*hnewCuBox_CuFrameu238M2_1;
	TH1			*hnewCuBox_CuFrameth232M2_01;
	TH1			*hnewCuBox_CuFrameu238M2_01;
	TH1			*hnewCuBox_CuFrameth232M2_001;
	TH1			*hnewCuBox_CuFrameu238M2_001;

	TH1			*hnewCuBox_CuFramepb210M2_100;
	TH1			*hnewCuBox_CuFrameth232M2_100;
	TH1			*hnewCuBox_CuFrameu238M2_100;
	TH1			*hnewCuBox_CuFramepb210M2_50;
	TH1			*hnewCuBox_CuFrameth232M2_50;
	TH1			*hnewCuBox_CuFrameu238M2_50;
	TH1			*hnewCuBox_CuFramepb210M2_5;
	TH1			*hnewCuBox_CuFrameth232M2_5;
	TH1			*hnewCuBox_CuFrameu238M2_5;


///////////// 50mK M1 and M2
	TH1			*hnew50mKcs137M1;
	TH1			*hnew50mKcs137M2;

//////////// Internal Shields M1 and M2
	TH1			*hnewInternalco60M1;
	TH1			*hnewInternalk40M1;
	TH1			*hnewInternalth232M1;
	TH1			*hnewInternalu238M1;

	TH1			*hnewInternalco60M2;
	TH1			*hnewInternalk40M2;
	TH1			*hnewInternalth232M2;
	TH1			*hnewInternalu238M2;


//////////// (PbRom) Roman Lead M1 and M2
	TH1			*hnewPbRombi207M1;
	TH1			*hnewPbRomco60M1;
	TH1			*hnewPbRomcs137M1;
	TH1			*hnewPbRomk40M1;
	TH1			*hnewPbRompb210M1;
	TH1			*hnewPbRomth232M1;
	TH1			*hnewPbRomu238M1;		

	TH1			*hnewPbRombi207M2;
	TH1			*hnewPbRomco60M2;
	TH1			*hnewPbRomcs137M2;
	TH1			*hnewPbRomk40M2;
	TH1			*hnewPbRompb210M2;
	TH1			*hnewPbRomth232M2;
	TH1			*hnewPbRomu238M2;		

//////////// External Shield M1 and M2
	TH1			*hnewExtPbbi210M1;
	TH1			*hnewExtPbk40M1;
	TH1			*hnewExtPbth232M1;
	TH1			*hnewExtPbu238M1;
	TH1			*hnewExtPbpb210M1;

	TH1			*hnewExtPbbi210M2;
	TH1			*hnewExtPbk40M2;
	TH1			*hnewExtPbth232M2;
	TH1			*hnewExtPbu238M2;
	TH1			*hnewExtPbpb210M2;
	
	TH1 		*hnewExtMuonM1;
	TH1			*hnewCuBox_th232spotM1;
	TH1			*hnewCuBox_k40spotM1;
	TH1 		*hnewBotExtPb_k40spotM1;

	TH1 		*hnewExtMuonM2;
	TH1			*hnewCuBox_th232spotM2;
	TH1			*hnewCuBox_k40spotM2;
	TH1 		*hnewBotExtPb_k40spotM2;


/////////// OVC M1 and M2
	TH1			*hnewOVCco60M1;
	TH1			*hnewOVCk40M1;
	TH1			*hnewOVCth232M1;
	TH1			*hnewOVCu238M1;		
	TH1 		*hnewOVCbi207M1;

	TH1			*hnewOVCco60M2;
	TH1			*hnewOVCk40M2;
	TH1			*hnewOVCth232M2;
	TH1			*hnewOVCu238M2;	
	TH1 		*hnewOVCbi207M2;

	TH1D 			*hOut;

	TH1D 			*hChiSquaredProgressM1;
	TH1D 			*hChiSquaredProgressM2;

	TDatime 		*tTime;

	// Cut Efficiency
	TH1D 			*hEfficiency;
	TH1D 			*hEfficiencyM2;

	std::ofstream 		OutFile;

	TTree 		*ProfileTree;
	TTree 		*ToyTree;
	TFile 		*ProfileFile;
	TFile 		*ToyTreeFile;

	std::vector<double> 	fInitValues;
	std::vector<double> 	fInitValues2;	

	TFile *fBulk;
	TFile *fSurface;
	TFile *fBulk_CDR;
	TFile *fBulk_CDRInternal;

	// TFile *fSaveResult;
	TFile *fToyData;

	std::string		dDataDir;
	std::string 	dMCDir;
	std::string 	dSaveDir;

	TVectorD *fParArray;
	TVectorD *fParArrayErr;

	bool 			bSave;

	int 			dNumCalls;
	int 			dMult;
	double			dBestChiSq;

	int 			dStartAlpha;
	int 			dEndAlpha;

	// Parameters
	double				fParameters[50];
	double				fParError[50];
	double				fParCountsM1[50]; // Integral of events in M1 spectrum
	double 				fParCountsM2[50];
	double 				fParActivityM1[50]; // Integral of events in M1 spectrum
	double 				fParActivityM2[50];
	double 				fParActivityErr[50];
	double 				fParMass[50]; // Mass of all elements
	double 				fParSurfaceArea[50]; // Surface area of all elements
	double				fResolution[52];
	double 				fParEfficiencyM1[50]; // Efficiency of the parameters 
	double 				fParPrior[50];
	double				dSecToYears;
	double				fMCEff[62];

 ClassDef(TBackgroundModel,1) // 
    };

#endif // __TBackgroundModel__
