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

class TBackgroundModel : public TObject {

public:
	TBackgroundModel();

	TBackgroundModel(double fFitMin, double fFitMax, int fBinBase, int fDataset, bool fSave);
	virtual ~TBackgroundModel();

	std::vector<double> AdaptiveBinning(TH1D *h1, int dBinBase);

	TH1D* CalculateResidualsAdaptive(TH1D *h1, TH1D *h2, TH1D *hResid, int binMin, int binMax, int dMult);

	bool DoTheFitAdaptive();

	void DrawBkg();

	double GetChiSquareAdaptive();

	void Initialize();

	TH1D *Kernal(TH1D *hMC, TH1D *hSMC);

	void LatexResultTable(double fValue);

	// Dumb to have all of these but w/e
	void LoadData();

	void PrintParameters();

	void PrintParActivity();

	void ProfileNLL();

	void ProfileNLL2D();

	void ResetParameters();

	void SanityCheck();

	void SetParameters(int index, double value);

	void SetParEfficiency();
	
	void ToyFit(int fNumFits);

	void UpdateModelAdaptive();

	int 	dNParam;
	int 	dBinSize;
	int 	dBaseBinSize;
	int 	dNBins;
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

  	TBkgModelParameter *BkgPar[100];

  	bool 	bFixedArray[100];

private:

	// Data
	TChain			*qtree;
	TCut 			base_cut;
	TCut			ener_cut;

	double			dDataIntegral;

	TMinuit			*minuit;

	TH1D			*fDataHistoTot;
	TH1D			*fDataHistoM1;
	TH1D			*fDataHistoM2;
	TH1D			*fDataHistoM2Sum;

	TH1D			*fAdapDataHistoM1;
	TH1D			*fAdapDataHistoM2;
	TH1D			*fAdapDataHistoM2Sum;


	// Updated 01-20-2015
	// Total PDFs M1
	TH1D			*fModelTotM1;
	TH1D			*fModelTotAdapM1;
	
	TH1D			*fModelTotM2;
	TH1D			*fModelTotAdapM2;
	
	TH1D			*fModelTotM2Sum;
	TH1D			*fModelTotAdapM2Sum;
	
//////////// Residual distributions
	TGraph			*gResidualM1;
	TGraph 			*gResidualM2;
	TGraph 			*gResidualM2Sum;

	TH1D			*hResidualDistM1;
	TH1D			*hResidualDistM2;
	TH1D			*hResidualDistM2Sum;

	TH1D			*hResidualGausM1;
	TH1D			*hResidualGausM2;
	TH1D			*hResidualGausM2Sum;


//////////////////// MC Histograms
////////// Crystal M1 and M2
	TH1D			*hTeO20nuM1;
	TH1D			*hTeO22nuM1;
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
	TH1D 			*hTeO2u238th230M2;

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

	TH1D			*hTeO20nuM2Sum;
	TH1D			*hTeO22nuM2Sum;
	TH1D			*hTeO2co60M2Sum;
	TH1D			*hTeO2k40M2Sum;
	TH1D			*hTeO2pb210M2Sum;
	TH1D			*hTeO2po210M2Sum;
	TH1D			*hTeO2te125M2Sum;
	TH1D			*hTeO2th232M2Sum;
	TH1D			*hTeO2th228M2Sum;
	TH1D			*hTeO2ra226M2Sum;
	TH1D			*hTeO2rn222M2Sum;
	TH1D			*hTeO2u238M2Sum;
	TH1D			*hTeO2th230M2Sum;
	TH1D			*hTeO2u234M2Sum;
	TH1D 			*hTeO2sb125M2Sum;

	TH1D			*hTeO2th232onlyM2Sum;
	TH1D			*hTeO2ra228pb208M2Sum;
	TH1D			*hTeO2th230onlyM2Sum;
	TH1D 			*hTeO2u238th230M2Sum;
	TH1D 			*hTeO2ra226pb210M2Sum;

	TH1D			*hTeO2Sxth232onlyM2Sum_001;
	TH1D			*hTeO2Sxra228pb208M2Sum_001;
	TH1D			*hTeO2Sxu238th230M2Sum_001;
	TH1D			*hTeO2Sxth230onlyM2Sum_001;
	TH1D			*hTeO2Sxra226pb210M2Sum_001;
	TH1D			*hTeO2Sxpb210M2Sum_0001;

	TH1D			*hTeO2Sxth232onlyM2Sum_01;
	TH1D			*hTeO2Sxra228pb208M2Sum_01;
	TH1D			*hTeO2Sxu238th230M2Sum_01;
	TH1D			*hTeO2Sxth230onlyM2Sum_01;
	TH1D			*hTeO2Sxra226pb210M2Sum_01;

	TH1D			*hTeO2Sxth232onlyM2Sum_0001;
	TH1D			*hTeO2Sxra228pb208M2Sum_0001;
	TH1D			*hTeO2Sxu238th230M2Sum_0001;
	TH1D			*hTeO2Sxth230onlyM2Sum_0001;
	TH1D			*hTeO2Sxra226pb210M2Sum_0001;

	TH1D			*hTeO2Spb210M2Sum_01;
	TH1D			*hTeO2Spo210M2Sum_001;
	TH1D			*hTeO2Spo210M2Sum_01;
	TH1D			*hTeO2Sth232M2Sum_01;
	TH1D			*hTeO2Su238M2Sum_01;
	TH1D			*hTeO2Sxpb210M2Sum_001;
	TH1D			*hTeO2Sxpb210M2Sum_01;
	TH1D			*hTeO2Sxpb210M2Sum_1;
	TH1D			*hTeO2Sxpb210M2Sum_10;
	TH1D			*hTeO2Sxpo210M2Sum_001;
	TH1D			*hTeO2Sxpo210M2Sum_01;
	TH1D			*hTeO2Sxpo210M2Sum_1;
	TH1D			*hTeO2Sxth232M2Sum_001;
	TH1D			*hTeO2Sxth232M2Sum_01;
	TH1D			*hTeO2Sxth232M2Sum_1;
	TH1D			*hTeO2Sxth232M2Sum_10;
	TH1D			*hTeO2Sxu238M2Sum_001;
	TH1D			*hTeO2Sxu238M2Sum_01;
	TH1D			*hTeO2Sxu238M2Sum_1;
	TH1D			*hTeO2Sxu238M2Sum_10;

	TH1D			*hTeO2Sxu238M2Sum_100;
	TH1D			*hTeO2Sxth232M2Sum_100;
	TH1D			*hTeO2Sxpb210M2Sum_100;


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

	TH1D			*hCuBox_CuFrameco60M2Sum;
	TH1D			*hCuBox_CuFramek40M2Sum;
	TH1D			*hCuBox_CuFrameth232M2Sum;
	TH1D			*hCuBox_CuFrameu238M2Sum;

	TH1D			*hCuBox_CuFrameth232M2Sum_10;
	TH1D			*hCuBox_CuFrameu238M2Sum_10;
	TH1D			*hCuBox_CuFramepb210M2Sum_10;
	TH1D			*hCuBox_CuFramepb210M2Sum_1;
	TH1D			*hCuBox_CuFramepb210M2Sum_01;
	TH1D			*hCuBox_CuFramepb210M2Sum_001;

	TH1D			*hCuBox_CuFrameth232M2Sum_1;
	TH1D			*hCuBox_CuFrameu238M2Sum_1;
	TH1D			*hCuBox_CuFrameth232M2Sum_01;
	TH1D			*hCuBox_CuFrameu238M2Sum_01;
	TH1D			*hCuBox_CuFrameth232M2Sum_001;
	TH1D			*hCuBox_CuFrameu238M2Sum_001;

	TH1D			*hCuBox_CuFramepb210M2Sum_100;
	TH1D			*hCuBox_CuFrameth232M2Sum_100;
	TH1D			*hCuBox_CuFrameu238M2Sum_100;
	TH1D			*hCuBox_CuFramepb210M2Sum_50;
	TH1D			*hCuBox_CuFrameth232M2Sum_50;
	TH1D			*hCuBox_CuFrameu238M2Sum_50;
	TH1D			*hCuBox_CuFramepb210M2Sum_5;
	TH1D			*hCuBox_CuFrameth232M2Sum_5;
	TH1D			*hCuBox_CuFrameu238M2Sum_5;

///////////// 50mK M1 and M2
	TH1D			*h50mKcs137M1;
	TH1D			*h50mKcs137M2;
	TH1D			*h50mKcs137M2Sum;

//////////// Internal Shields M1 and M2
	TH1D			*hInternalco60M1;
	TH1D			*hInternalk40M1;
	TH1D			*hInternalth232M1;
	TH1D			*hInternalu238M1;

	TH1D			*hInternalco60M2;
	TH1D			*hInternalk40M2;
	TH1D			*hInternalth232M2;
	TH1D			*hInternalu238M2;

	TH1D			*hInternalco60M2Sum;
	TH1D			*hInternalk40M2Sum;
	TH1D			*hInternalth232M2Sum;
	TH1D			*hInternalu238M2Sum;	


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

	TH1D			*hPbRombi207M2Sum;
	TH1D			*hPbRomco60M2Sum;
	TH1D			*hPbRomcs137M2Sum;
	TH1D			*hPbRomk40M2Sum;
	TH1D			*hPbRompb210M2Sum;
	TH1D			*hPbRomth232M2Sum;
	TH1D			*hPbRomu238M2Sum;	

//////////// External Shield M1 and M2
	TH1D			*hExtPbbi210M1;

	TH1D			*hExtPbbi210M2;
	
	TH1D			*hExtPbbi210M2Sum;

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

	TH1D			*hOVCco60M2Sum;
	TH1D			*hOVCk40M2Sum;
	TH1D			*hOVCth232M2Sum;
	TH1D			*hOVCu238M2Sum;	

////////// Fudge Factor Models
	TH1D 		*hOVC804M1;
	TH1D 		*hOVC835M1;
	TH1D 		*hOVC1063M1;

	TH1D 		*hOVC804M2;
	TH1D 		*hOVC835M2;
	TH1D 		*hOVC1063M2;

	TH1D 		*hOVC804M2Sum;
	TH1D 		*hOVC835M2Sum;
	TH1D 		*hOVC1063M2Sum;

///////////////// Adaptive binned histograms
/////////// Crystal M1 and M2
	TH1D			*hAdapTeO20nuM1;
	TH1D			*hAdapTeO22nuM1;
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
	TH1D 			*hAdapTeO2u238th230M1;

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
	TH1D 			*hAdapTeO2u238th230M2;

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


	TH1D			*hAdapTeO20nuM2Sum;
	TH1D			*hAdapTeO22nuM2Sum;
	TH1D			*hAdapTeO2co60M2Sum;
	TH1D			*hAdapTeO2k40M2Sum;
	TH1D			*hAdapTeO2pb210M2Sum;
	TH1D			*hAdapTeO2po210M2Sum;
	TH1D			*hAdapTeO2te125M2Sum;
	TH1D			*hAdapTeO2th232M2Sum;
	TH1D			*hAdapTeO2th228M2Sum;
	TH1D			*hAdapTeO2ra226M2Sum;
	TH1D			*hAdapTeO2rn222M2Sum;
	TH1D			*hAdapTeO2u238M2Sum;
	TH1D			*hAdapTeO2th230M2Sum;
	TH1D			*hAdapTeO2u234M2Sum;
	TH1D 			*hAdapTeO2sb125M2Sum;


	TH1D			*hAdapTeO2th232onlyM2Sum;
	TH1D			*hAdapTeO2ra228pb208M2Sum;
	TH1D			*hAdapTeO2th230onlyM2Sum;
	TH1D 			*hAdapTeO2u238th230M2Sum;
	TH1D 			*hAdapTeO2ra226pb210M2Sum;

	TH1D			*hAdapTeO2Sxth232onlyM2Sum_001;
	TH1D			*hAdapTeO2Sxra228pb208M2Sum_001;
	TH1D			*hAdapTeO2Sxu238th230M2Sum_001;
	TH1D			*hAdapTeO2Sxth230onlyM2Sum_001;
	TH1D			*hAdapTeO2Sxra226pb210M2Sum_001;
	TH1D			*hAdapTeO2Sxpb210M2Sum_0001;

	TH1D			*hAdapTeO2Sxth232onlyM2Sum_01;
	TH1D			*hAdapTeO2Sxra228pb208M2Sum_01;
	TH1D			*hAdapTeO2Sxu238th230M2Sum_01;
	TH1D			*hAdapTeO2Sxth230onlyM2Sum_01;
	TH1D			*hAdapTeO2Sxra226pb210M2Sum_01;

	TH1D			*hAdapTeO2Sxth232onlyM2Sum_0001;
	TH1D			*hAdapTeO2Sxra228pb208M2Sum_0001;
	TH1D			*hAdapTeO2Sxu238th230M2Sum_0001;
	TH1D			*hAdapTeO2Sxth230onlyM2Sum_0001;
	TH1D			*hAdapTeO2Sxra226pb210M2Sum_0001;

	TH1D			*hAdapTeO2Spb210M2Sum_01;
	TH1D			*hAdapTeO2Spo210M2Sum_001;
	TH1D			*hAdapTeO2Spo210M2Sum_01;
	TH1D			*hAdapTeO2Sth232M2Sum_01;
	TH1D			*hAdapTeO2Su238M2Sum_01;
	TH1D			*hAdapTeO2Sxpb210M2Sum_001;
	TH1D			*hAdapTeO2Sxpb210M2Sum_01;
	TH1D			*hAdapTeO2Sxpb210M2Sum_1;
	TH1D			*hAdapTeO2Sxpb210M2Sum_10;
	TH1D			*hAdapTeO2Sxpo210M2Sum_001;
	TH1D			*hAdapTeO2Sxpo210M2Sum_01;
	TH1D			*hAdapTeO2Sxpo210M2Sum_1;
	TH1D			*hAdapTeO2Sxth232M2Sum_001;
	TH1D			*hAdapTeO2Sxth232M2Sum_01;
	TH1D			*hAdapTeO2Sxth232M2Sum_1;
	TH1D			*hAdapTeO2Sxth232M2Sum_10;
	TH1D			*hAdapTeO2Sxu238M2Sum_001;
	TH1D			*hAdapTeO2Sxu238M2Sum_01;
	TH1D			*hAdapTeO2Sxu238M2Sum_1;
	TH1D			*hAdapTeO2Sxu238M2Sum_10;

	TH1D			*hAdapTeO2Sxu238M2Sum_100;
	TH1D			*hAdapTeO2Sxth232M2Sum_100;
	TH1D			*hAdapTeO2Sxpb210M2Sum_100;

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

	TH1D			*hAdapCuBox_CuFrameco60M2Sum;
	TH1D			*hAdapCuBox_CuFramek40M2Sum;
	TH1D			*hAdapCuBox_CuFrameth232M2Sum;
	TH1D			*hAdapCuBox_CuFrameu238M2Sum;

	TH1D			*hAdapCuBox_CuFrameth232M2Sum_10;
	TH1D			*hAdapCuBox_CuFrameu238M2Sum_10;
	TH1D			*hAdapCuBox_CuFramepb210M2Sum_10;
	TH1D			*hAdapCuBox_CuFramepb210M2Sum_1;
	TH1D			*hAdapCuBox_CuFramepb210M2Sum_01;
	TH1D			*hAdapCuBox_CuFramepb210M2Sum_001;

	TH1D			*hAdapCuBox_CuFrameth232M2Sum_1;
	TH1D			*hAdapCuBox_CuFrameu238M2Sum_1;
	TH1D			*hAdapCuBox_CuFrameth232M2Sum_01;
	TH1D			*hAdapCuBox_CuFrameu238M2Sum_01;
	TH1D			*hAdapCuBox_CuFrameth232M2Sum_001;
	TH1D			*hAdapCuBox_CuFrameu238M2Sum_001;

	TH1D			*hAdapCuBox_CuFramepb210M2Sum_100;
	TH1D			*hAdapCuBox_CuFrameth232M2Sum_100;
	TH1D			*hAdapCuBox_CuFrameu238M2Sum_100;
	TH1D			*hAdapCuBox_CuFramepb210M2Sum_50;
	TH1D			*hAdapCuBox_CuFrameth232M2Sum_50;
	TH1D			*hAdapCuBox_CuFrameu238M2Sum_50;
	TH1D			*hAdapCuBox_CuFramepb210M2Sum_5;
	TH1D			*hAdapCuBox_CuFrameth232M2Sum_5;
	TH1D			*hAdapCuBox_CuFrameu238M2Sum_5;

/////////// 50mK M1 and M2
	TH1D			*hAdap50mKcs137M1;
	TH1D			*hAdap50mKcs137M2;
	TH1D			*hAdap50mKcs137M2Sum;
	
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

	TH1D			*hAdapPbRombi207M2Sum;
	TH1D			*hAdapPbRomco60M2Sum;
	TH1D			*hAdapPbRomcs137M2Sum;
	TH1D			*hAdapPbRomk40M2Sum;
	TH1D			*hAdapPbRompb210M2Sum;
	TH1D			*hAdapPbRomth232M2Sum;
	TH1D			*hAdapPbRomu238M2Sum;	


////////// Internal Shields M1 and M2
	TH1D			*hAdapInternalco60M1;
	TH1D			*hAdapInternalk40M1;
	TH1D			*hAdapInternalth232M1;
	TH1D			*hAdapInternalu238M1;

	TH1D			*hAdapInternalco60M2;
	TH1D			*hAdapInternalk40M2;
	TH1D			*hAdapInternalth232M2;
	TH1D			*hAdapInternalu238M2;

	TH1D			*hAdapInternalco60M2Sum;
	TH1D			*hAdapInternalk40M2Sum;
	TH1D			*hAdapInternalth232M2Sum;
	TH1D			*hAdapInternalu238M2Sum;	

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

	TH1D			*hAdapOVCco60M2Sum;
	TH1D			*hAdapOVCk40M2Sum;
	TH1D			*hAdapOVCth232M2Sum;
	TH1D			*hAdapOVCu238M2Sum;

///////// External Shield M1 and M2
	TH1D			*hAdapExtPbbi210M1;

	TH1D			*hAdapExtPbbi210M2;
	
	TH1D			*hAdapExtPbbi210M2Sum;

////////// Fudge Factor Models
	TH1D 		*hAdapOVC804M1;
	TH1D 		*hAdapOVC835M1;
	TH1D 		*hAdapOVC1063M1;

	TH1D 		*hAdapOVC804M2;
	TH1D 		*hAdapOVC835M2;
	TH1D 		*hAdapOVC1063M2;

	TH1D 		*hAdapOVC804M2Sum;
	TH1D 		*hAdapOVC835M2Sum;
	TH1D 		*hAdapOVC1063M2Sum;

////////////// Dummy histograms
	TH1 			*hnewM1;
	TH1				*hnewM2;
	TH1				*hnewM2Sum;

////////// Crystal M1 and M2
	TH1			*hnewTeO20nuM1;
	TH1			*hnewTeO22nuM1;
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
	TH1 		*hnewTeO2u238th230M1;

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
	TH1 		*hnewTeO2u238th230M2;

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


	TH1			*hnewTeO20nuM2Sum;
	TH1			*hnewTeO22nuM2Sum;
	TH1			*hnewTeO2co60M2Sum;
	TH1			*hnewTeO2k40M2Sum;
	TH1			*hnewTeO2pb210M2Sum;
	TH1			*hnewTeO2po210M2Sum;
	TH1			*hnewTeO2te125M2Sum;
	TH1			*hnewTeO2th232M2Sum;
	TH1			*hnewTeO2th228M2Sum;
	TH1			*hnewTeO2ra226M2Sum;
	TH1			*hnewTeO2rn222M2Sum;
	TH1			*hnewTeO2u238M2Sum;
	TH1			*hnewTeO2th230M2Sum;
	TH1			*hnewTeO2u234M2Sum;
	TH1 		*hnewTeO2sb125M2Sum;

	TH1			*hnewTeO2th232onlyM2Sum;
	TH1			*hnewTeO2ra228pb208M2Sum;
	TH1			*hnewTeO2th230onlyM2Sum;
	TH1 		*hnewTeO2u238th230M2Sum;
	TH1 		*hnewTeO2ra226pb210M2Sum;

	TH1			*hnewTeO2Sxth232onlyM2Sum_001;
	TH1			*hnewTeO2Sxra228pb208M2Sum_001;
	TH1			*hnewTeO2Sxu238th230M2Sum_001;
	TH1			*hnewTeO2Sxth230onlyM2Sum_001;
	TH1			*hnewTeO2Sxra226pb210M2Sum_001;
	TH1			*hnewTeO2Sxpb210M2Sum_0001;

	TH1			*hnewTeO2Sxth232onlyM2Sum_01;
	TH1			*hnewTeO2Sxra228pb208M2Sum_01;
	TH1			*hnewTeO2Sxu238th230M2Sum_01;
	TH1			*hnewTeO2Sxth230onlyM2Sum_01;
	TH1			*hnewTeO2Sxra226pb210M2Sum_01;

	TH1			*hnewTeO2Sxth232onlyM2Sum_0001;
	TH1			*hnewTeO2Sxra228pb208M2Sum_0001;
	TH1			*hnewTeO2Sxu238th230M2Sum_0001;
	TH1			*hnewTeO2Sxth230onlyM2Sum_0001;
	TH1			*hnewTeO2Sxra226pb210M2Sum_0001;	

	TH1			*hnewTeO2Spb210M2Sum_01;
	TH1			*hnewTeO2Spo210M2Sum_001;
	TH1			*hnewTeO2Spo210M2Sum_01;
	TH1			*hnewTeO2Sth232M2Sum_01;
	TH1			*hnewTeO2Su238M2Sum_01;
	TH1			*hnewTeO2Sxpb210M2Sum_001;
	TH1			*hnewTeO2Sxpb210M2Sum_01;
	TH1			*hnewTeO2Sxpb210M2Sum_1;
	TH1			*hnewTeO2Sxpb210M2Sum_10;
	TH1			*hnewTeO2Sxpo210M2Sum_001;
	TH1			*hnewTeO2Sxpo210M2Sum_01;
	TH1			*hnewTeO2Sxpo210M2Sum_1;
	TH1			*hnewTeO2Sxth232M2Sum_001;
	TH1			*hnewTeO2Sxth232M2Sum_01;
	TH1			*hnewTeO2Sxth232M2Sum_1;
	TH1			*hnewTeO2Sxth232M2Sum_10;
	TH1			*hnewTeO2Sxu238M2Sum_001;
	TH1			*hnewTeO2Sxu238M2Sum_01;
	TH1			*hnewTeO2Sxu238M2Sum_1;
	TH1			*hnewTeO2Sxu238M2Sum_10;

	TH1			*hnewTeO2Sxu238M2Sum_100;
	TH1			*hnewTeO2Sxth232M2Sum_100;
	TH1			*hnewTeO2Sxpb210M2Sum_100;

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

	TH1			*hnewCuBox_CuFrameco60M2Sum;
	TH1			*hnewCuBox_CuFramek40M2Sum;
	TH1			*hnewCuBox_CuFrameth232M2Sum;
	TH1			*hnewCuBox_CuFrameu238M2Sum;

	TH1			*hnewCuBox_CuFrameth232M2Sum_10;
	TH1			*hnewCuBox_CuFrameu238M2Sum_10;
	TH1			*hnewCuBox_CuFramepb210M2Sum_10;
	TH1			*hnewCuBox_CuFramepb210M2Sum_1;
	TH1			*hnewCuBox_CuFramepb210M2Sum_01;
	TH1			*hnewCuBox_CuFramepb210M2Sum_001;

	TH1			*hnewCuBox_CuFrameth232M2Sum_1;
	TH1			*hnewCuBox_CuFrameu238M2Sum_1;
	TH1			*hnewCuBox_CuFrameth232M2Sum_01;
	TH1			*hnewCuBox_CuFrameu238M2Sum_01;
	TH1			*hnewCuBox_CuFrameth232M2Sum_001;
	TH1			*hnewCuBox_CuFrameu238M2Sum_001;

	TH1			*hnewCuBox_CuFramepb210M2Sum_100;
	TH1			*hnewCuBox_CuFrameth232M2Sum_100;
	TH1			*hnewCuBox_CuFrameu238M2Sum_100;
	TH1			*hnewCuBox_CuFramepb210M2Sum_50;
	TH1			*hnewCuBox_CuFrameth232M2Sum_50;
	TH1			*hnewCuBox_CuFrameu238M2Sum_50;
	TH1			*hnewCuBox_CuFramepb210M2Sum_5;
	TH1			*hnewCuBox_CuFrameth232M2Sum_5;
	TH1			*hnewCuBox_CuFrameu238M2Sum_5;

///////////// 50mK M1 and M2
	TH1			*hnew50mKcs137M1;
	TH1			*hnew50mKcs137M2;
	TH1			*hnew50mKcs137M2Sum;

//////////// Internal Shields M1 and M2
	TH1			*hnewInternalco60M1;
	TH1			*hnewInternalk40M1;
	TH1			*hnewInternalth232M1;
	TH1			*hnewInternalu238M1;

	TH1			*hnewInternalco60M2;
	TH1			*hnewInternalk40M2;
	TH1			*hnewInternalth232M2;
	TH1			*hnewInternalu238M2;

	TH1			*hnewInternalco60M2Sum;
	TH1			*hnewInternalk40M2Sum;
	TH1			*hnewInternalth232M2Sum;
	TH1			*hnewInternalu238M2Sum;	


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

	TH1			*hnewPbRombi207M2Sum;
	TH1			*hnewPbRomco60M2Sum;
	TH1			*hnewPbRomcs137M2Sum;
	TH1			*hnewPbRomk40M2Sum;
	TH1			*hnewPbRompb210M2Sum;
	TH1			*hnewPbRomth232M2Sum;
	TH1			*hnewPbRomu238M2Sum;	

//////////// External Shield M1 and M2
	TH1			*hnewExtPbbi210M1;

	TH1			*hnewExtPbbi210M2;
	
	TH1			*hnewExtPbbi210M2Sum;

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

	TH1			*hnewOVCco60M2Sum;
	TH1			*hnewOVCk40M2Sum;
	TH1			*hnewOVCth232M2Sum;
	TH1			*hnewOVCu238M2Sum;	


////////// Fudge Factor Models
	TH1 		*hnewOVC804M1;
	TH1 		*hnewOVC835M1;
	TH1 		*hnewOVC1063M1;

	TH1 		*hnewOVC804M2;
	TH1 		*hnewOVC835M2;
	TH1 		*hnewOVC1063M2;

	TH1 		*hnewOVC804M2Sum;
	TH1 		*hnewOVC835M2Sum;
	TH1 		*hnewOVC1063M2Sum;

	TH1D 			*hOut;

	TH1D 			*hChiSquaredProgressM1;
	TH1D 			*hChiSquaredProgressM2;
	TH1D 			*hChiSquaredProgressM2Sum;

	TDatime 		*tTime;

	// Smearing
	TF1				*gaus;

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

	TFile *fBulk;
	TFile *fSurface;
	TFile *fBulk_CDR;

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

	// TFile *fSaveResult;
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
	double 				fParActivityM1[139]; // Integral of events in M1 spectrum
	double 				fParActivityM2[139];
	double 				fParActivityErr[139];
	double 				fParMass[139]; // Mass of all elements
	double 				fParSurfaceArea[139]; // Surface area of all elements
	double				fResolution[52];
	double 				fParEfficiencyM1[139]; // Efficiency of the parameters 
	double				dSecToYears;
	double				fMCEff[62];

 ClassDef(TBackgroundModel,1) // 
    };

#endif // __TBackgroundModel__
