#ifndef __TBackgroundModel__
#define __TBackgroundModel__
#include "TObject.h"
#include "TFile.h"
#include "TH1.h"
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

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <vector>
#include <map>
// #include <iostream>
#include <fstream>

class TBackgroundModel : public TObject {

public:
	TBackgroundModel();

	TBackgroundModel(double fFitMin, double fFitMax, int fBinBase, int fDataset);
	virtual ~TBackgroundModel();

	std::vector<double> AdaptiveBinning(TH1D *h1, int dBinBase);

	TH1D* CalculateResidualsAdaptive(TH1D *h1, TH1D *h2, TH1D *hResid, int binMin, int binMax, int dMult);

	bool DoTheFitAdaptive(double f2nuValue, double fk40Value);

	void DrawBkg();

	void DrawMC();

	TH1D *EnergyScale(TH1D *hIn, TH1D *hDummy, double dConst, double dSlope);

	double GetChiSquareAdaptive();

	void Initialize();

	void Initialize2();

	TH1D *Kernal(TH1D *hMC, TH1D *hSMC);

	void LatexResultTable(double fValue);

	// Dumb to have all of these but w/e
	void LoadData();

	void LoadPDFs();

	void PrintParameters();

	void PrintParActivity();

	void ProfileNLL(double fBestFitInit, double fBestFitChiSq);

	void ProfileNLL2D(double fBestFitInit, double fBestFitInit2, double fBestFitChiSq);

	void ReadMC();

	void ResetParameters();

	void SanityCheck();

	void SetParameters(int index, double value);

	void SetParEfficiency();
	
	void Test();

	void UpdateModelAdaptive();

	void WriteParameters();	


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

	int  	dDataIntegralTot;
	int 	dDataIntegralM1;
	int 	dDataIntegralM2;
	int 	dDataIntegralM2Sum;

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
	TH1D			*fModelTotthM1;
	TH1D			*fModelTotuM1;
	TH1D			*fModelTotkM1;
	TH1D			*fModelTotcoM1;
	TH1D			*fModelTotmnM1;
	TH1D			*fModelTotNDBDM1;
	TH1D			*fModelTot2NDBDM1;
	TH1D			*fModelTotbiM1;
	TH1D			*fModelTotbi2M1;	
	TH1D			*fModelTotptM1;
	TH1D			*fModelTotpbM1;
	TH1D			*fModelTotcsM1;
	TH1D			*fModelTotco2M1;
	TH1D			*fModelTotteo2M1;

	TH1D			*fModelTotSthM1;
	TH1D			*fModelTotSuM1;
	TH1D 			*fModelTotSpoM1;
	TH1D			*fModelTotSpbM1;
	TH1D			*fModelTotExtM1;

	TH1D			*fModelTotAdapM1;
	TH1D			*fModelTotAdapthM1;
	TH1D			*fModelTotAdapuM1;
	TH1D			*fModelTotAdapkM1;
	TH1D			*fModelTotAdapcoM1;
	TH1D			*fModelTotAdapmnM1;
	TH1D			*fModelTotAdapNDBDM1;
	TH1D			*fModelTotAdap2NDBDM1;
	TH1D			*fModelTotAdapbiM1;
	TH1D			*fModelTotAdapbi2M1;
	TH1D			*fModelTotAdapptM1;
	TH1D			*fModelTotAdappbM1;
	TH1D			*fModelTotAdapcsM1;
	TH1D			*fModelTotAdapco2M1;
	TH1D			*fModelTotAdapteo2M1;

	TH1D			*fModelTotAdapSthM1;
	TH1D			*fModelTotAdapSuM1;
	TH1D 			*fModelTotAdapSpoM1;
	TH1D			*fModelTotAdapSpbM1;
	TH1D			*fModelTotAdapExtM1;


	TH1D			*fModelTotAdapAlphaM1;
	TH1D			*fModelTotAdapAlphaHighM1;
	TH1D			*fModelTotAdapAlphaLowM1;

	TH1D 			*fModelTotAdapFudgeM1;

	// Total PDFs M2
	TH1D			*fModelTotM2;
	TH1D			*fModelTotthM2;
	TH1D			*fModelTotuM2;
	TH1D			*fModelTotkM2;
	TH1D			*fModelTotcoM2;
	TH1D			*fModelTotmnM2;
	TH1D			*fModelTotNDBDM2;
	TH1D			*fModelTot2NDBDM2;
	TH1D			*fModelTotbiM2;
	TH1D			*fModelTotbi2M2;
	TH1D			*fModelTotptM2;
	TH1D			*fModelTotpbM2;
	TH1D			*fModelTotcsM2;
	TH1D			*fModelTotco2M2;
	TH1D			*fModelTotteo2M2;

	TH1D			*fModelTotSthM2;
	TH1D			*fModelTotSuM2;
	TH1D 			*fModelTotSpoM2;
	TH1D			*fModelTotSpbM2;
	TH1D			*fModelTotExtM2;

	TH1D			*fModelTotAdapM2;
	TH1D			*fModelTotAdapthM2;
	TH1D			*fModelTotAdapuM2;
	TH1D			*fModelTotAdapkM2;
	TH1D			*fModelTotAdapcoM2;
	TH1D			*fModelTotAdapmnM2;
	TH1D			*fModelTotAdapNDBDM2;
	TH1D			*fModelTotAdap2NDBDM2;
	TH1D			*fModelTotAdapbiM2;
	TH1D			*fModelTotAdapbi2M2;
	TH1D			*fModelTotAdapptM2;
	TH1D			*fModelTotAdappbM2;
	TH1D			*fModelTotAdapcsM2;
	TH1D			*fModelTotAdapco2M2;
	TH1D			*fModelTotAdapteo2M2;

	TH1D			*fModelTotAdapSthM2;
	TH1D			*fModelTotAdapSuM2;
	TH1D 			*fModelTotAdapSpoM2;
	TH1D			*fModelTotAdapSpbM2;
	TH1D			*fModelTotAdapExtM2;


	// TH1D			*fModelTotAdapAlphaM2;
	// TH1D			*fModelTotAdapAlphaHighM2;
	// TH1D			*fModelTotAdapAlphaLowM2;

	// Total PDFs M2Sum
	TH1D			*fModelTotM2Sum;
	TH1D			*fModelTotthM2Sum;
	TH1D			*fModelTotuM2Sum;
	TH1D			*fModelTotkM2Sum;
	TH1D			*fModelTotcoM2Sum;
	TH1D			*fModelTotmnM2Sum;
	TH1D			*fModelTotNDBDM2Sum;
	TH1D			*fModelTot2NDBDM2Sum;
	TH1D			*fModelTotbiM2Sum;
	TH1D			*fModelTotbi2M2Sum;
	TH1D			*fModelTotptM2Sum;
	TH1D			*fModelTotpbM2Sum;
	TH1D			*fModelTotcsM2Sum;
	TH1D			*fModelTotco2M2Sum;
	TH1D			*fModelTotteo2M2Sum;

	TH1D			*fModelTotSthM2Sum;
	TH1D			*fModelTotSuM2Sum;
	TH1D 			*fModelTotSpoM2Sum;
	TH1D			*fModelTotSpbM2Sum;
	TH1D			*fModelTotExtM2Sum;

	TH1D			*fModelTotAdapM2Sum;
	TH1D			*fModelTotAdapthM2Sum;
	TH1D			*fModelTotAdapuM2Sum;
	TH1D			*fModelTotAdapkM2Sum;
	TH1D			*fModelTotAdapcoM2Sum;
	TH1D			*fModelTotAdapmnM2Sum;
	TH1D			*fModelTotAdapNDBDM2Sum;
	TH1D			*fModelTotAdap2NDBDM2Sum;
	TH1D			*fModelTotAdapbiM2Sum;
	TH1D			*fModelTotAdapbi2M2Sum;
	TH1D			*fModelTotAdapptM2Sum;
	TH1D			*fModelTotAdappbM2Sum;
	TH1D			*fModelTotAdapcsM2Sum;
	TH1D			*fModelTotAdapco2M2Sum;
	TH1D			*fModelTotAdapteo2M2Sum;

	TH1D			*fModelTotAdapSthM2Sum;
	TH1D			*fModelTotAdapSuM2Sum;
	TH1D 			*fModelTotAdapSpoM2Sum;
	TH1D			*fModelTotAdapSpbM2Sum;
	TH1D			*fModelTotAdapExtM2Sum;

	TH1D			*fModelTotAdapAlphaM2Sum;
	TH1D			*fModelTotAdapAlphaHighM2Sum;
	TH1D			*fModelTotAdapAlphaLowM2Sum;

	TH1D 			*fModelTotAdapFudgeM2;


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

	TH1D			*hTeO2th232onlyM1;
	TH1D			*hTeO2ra228pb208M1;
	TH1D			*hTeO2th230onlyM1;

	TH1D			*hTeO2Sxth232onlyM1_001;
	TH1D			*hTeO2Sxra228pb208M1_001;
	TH1D			*hTeO2Sxu238th230M1_001;
	TH1D			*hTeO2Sxth230onlyM1_001;
	TH1D			*hTeO2Sxra226pb210M1_001;
	TH1D			*hTeO2Sxpb210M1_0001;

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

	TH1D			*hTeO2th232onlyM2;
	TH1D			*hTeO2ra228pb208M2;
	TH1D			*hTeO2th230onlyM2;

	TH1D			*hTeO2Sxth232onlyM2_001;
	TH1D			*hTeO2Sxra228pb208M2_001;
	TH1D			*hTeO2Sxu238th230M2_001;
	TH1D			*hTeO2Sxth230onlyM2_001;
	TH1D			*hTeO2Sxra226pb210M2_001;
	TH1D			*hTeO2Sxpb210M2_0001;


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

	TH1D			*hTeO2th232onlyM2Sum;
	TH1D			*hTeO2ra228pb208M2Sum;
	TH1D			*hTeO2th230onlyM2Sum;

	TH1D			*hTeO2Sxth232onlyM2Sum_001;
	TH1D			*hTeO2Sxra228pb208M2Sum_001;
	TH1D			*hTeO2Sxu238th230M2Sum_001;
	TH1D			*hTeO2Sxth230onlyM2Sum_001;
	TH1D			*hTeO2Sxra226pb210M2Sum_001;
	TH1D			*hTeO2Sxpb210M2Sum_0001;

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

///////// Frame M1 and M2
	TH1D			*hCuFrameco58M1;
	TH1D			*hCuFrameco60M1;
	TH1D			*hCuFramecs137M1;
	TH1D			*hCuFramek40M1;
	TH1D			*hCuFramemn54M1;
	TH1D			*hCuFramepb210M1;
	TH1D			*hCuFrameth232M1;
	TH1D			*hCuFrameu238M1;

	TH1D			*hCuFrameSth232M1_1;
	TH1D			*hCuFrameSu238M1_1;
	TH1D			*hCuFrameSxpb210M1_001;
	TH1D			*hCuFrameSxpb210M1_01;
	TH1D			*hCuFrameSxpb210M1_1;
	TH1D			*hCuFrameSxpb210M1_10;
	TH1D			*hCuFrameSxth232M1_001;
	TH1D			*hCuFrameSxth232M1_01;
	TH1D			*hCuFrameSxth232M1_1;
	TH1D			*hCuFrameSxth232M1_10;
	TH1D			*hCuFrameSxu238M1_001;
	TH1D			*hCuFrameSxu238M1_01;
	TH1D			*hCuFrameSxu238M1_1;
	TH1D			*hCuFrameSxu238M1_10;

	TH1D			*hCuFrameco58M2;
	TH1D			*hCuFrameco60M2;
	TH1D			*hCuFramecs137M2;
	TH1D			*hCuFramek40M2;
	TH1D			*hCuFramemn54M2;
	TH1D			*hCuFramepb210M2;
	TH1D			*hCuFrameth232M2;
	TH1D			*hCuFrameu238M2;

	TH1D			*hCuFrameSth232M2_1;
	TH1D			*hCuFrameSu238M2_1;
	TH1D			*hCuFrameSxpb210M2_001;
	TH1D			*hCuFrameSxpb210M2_01;
	TH1D			*hCuFrameSxpb210M2_1;
	TH1D			*hCuFrameSxpb210M2_10;
	TH1D			*hCuFrameSxth232M2_001;
	TH1D			*hCuFrameSxth232M2_01;
	TH1D			*hCuFrameSxth232M2_1;
	TH1D			*hCuFrameSxth232M2_10;
	TH1D			*hCuFrameSxu238M2_001;
	TH1D			*hCuFrameSxu238M2_01;
	TH1D			*hCuFrameSxu238M2_1;
	TH1D			*hCuFrameSxu238M2_10;

	TH1D			*hCuFrameco58M2Sum;
	TH1D			*hCuFrameco60M2Sum;
	TH1D			*hCuFramecs137M2Sum;
	TH1D			*hCuFramek40M2Sum;
	TH1D			*hCuFramemn54M2Sum;
	TH1D			*hCuFramepb210M2Sum;
	TH1D			*hCuFrameth232M2Sum;
	TH1D			*hCuFrameu238M2Sum;

	TH1D			*hCuFrameSth232M2Sum_1;
	TH1D			*hCuFrameSu238M2Sum_1;
	TH1D			*hCuFrameSxpb210M2Sum_001;
	TH1D			*hCuFrameSxpb210M2Sum_01;
	TH1D			*hCuFrameSxpb210M2Sum_1;
	TH1D			*hCuFrameSxpb210M2Sum_10;
	TH1D			*hCuFrameSxth232M2Sum_001;
	TH1D			*hCuFrameSxth232M2Sum_01;
	TH1D			*hCuFrameSxth232M2Sum_1;
	TH1D			*hCuFrameSxth232M2Sum_10;
	TH1D			*hCuFrameSxu238M2Sum_001;
	TH1D			*hCuFrameSxu238M2Sum_01;
	TH1D			*hCuFrameSxu238M2Sum_1;
	TH1D			*hCuFrameSxu238M2Sum_10;

/////////// CuBox (TShield) M1 and M2
	TH1D			*hCuBoxco58M1;
	TH1D			*hCuBoxco60M1;
	TH1D			*hCuBoxcs137M1;
	TH1D			*hCuBoxk40M1;
	TH1D			*hCuBoxmn54M1;
	TH1D			*hCuBoxpb210M1;
	TH1D			*hCuBoxth232M1;
	TH1D			*hCuBoxu238M1;	

	TH1D			*hCuBoxSth232M1_1;
	TH1D			*hCuBoxSu238M1_1;
	TH1D			*hCuBoxSxpb210M1_001;
	TH1D			*hCuBoxSxpb210M1_01;
	TH1D			*hCuBoxSxpb210M1_1;
	TH1D			*hCuBoxSxpb210M1_10;
	TH1D			*hCuBoxSxth232M1_001;
	TH1D			*hCuBoxSxth232M1_01;
	TH1D			*hCuBoxSxth232M1_1;
	TH1D			*hCuBoxSxth232M1_10;
	TH1D			*hCuBoxSxu238M1_001;
	TH1D			*hCuBoxSxu238M1_01;
	TH1D			*hCuBoxSxu238M1_1;
	TH1D			*hCuBoxSxu238M1_10;

	TH1D			*hCuBoxco58M2;
	TH1D			*hCuBoxco60M2;
	TH1D			*hCuBoxcs137M2;
	TH1D			*hCuBoxk40M2;
	TH1D			*hCuBoxmn54M2;
	TH1D			*hCuBoxpb210M2;
	TH1D			*hCuBoxth232M2;
	TH1D			*hCuBoxu238M2;	

	TH1D			*hCuBoxSth232M2_1;
	TH1D			*hCuBoxSu238M2_1;
	TH1D			*hCuBoxSxpb210M2_001;
	TH1D			*hCuBoxSxpb210M2_01;
	TH1D			*hCuBoxSxpb210M2_1;
	TH1D			*hCuBoxSxpb210M2_10;
	TH1D			*hCuBoxSxth232M2_001;
	TH1D			*hCuBoxSxth232M2_01;
	TH1D			*hCuBoxSxth232M2_1;
	TH1D			*hCuBoxSxth232M2_10;
	TH1D			*hCuBoxSxu238M2_001;
	TH1D			*hCuBoxSxu238M2_01;
	TH1D			*hCuBoxSxu238M2_1;
	TH1D			*hCuBoxSxu238M2_10;

	TH1D			*hCuBoxco58M2Sum;
	TH1D			*hCuBoxco60M2Sum;
	TH1D			*hCuBoxcs137M2Sum;
	TH1D			*hCuBoxk40M2Sum;
	TH1D			*hCuBoxmn54M2Sum;
	TH1D			*hCuBoxpb210M2Sum;
	TH1D			*hCuBoxth232M2Sum;
	TH1D			*hCuBoxu238M2Sum;	

	TH1D			*hCuBoxSth232M2Sum_1;
	TH1D			*hCuBoxSu238M2Sum_1;
	TH1D			*hCuBoxSxpb210M2Sum_001;
	TH1D			*hCuBoxSxpb210M2Sum_01;
	TH1D			*hCuBoxSxpb210M2Sum_1;
	TH1D			*hCuBoxSxpb210M2Sum_10;
	TH1D			*hCuBoxSxth232M2Sum_001;
	TH1D			*hCuBoxSxth232M2Sum_01;
	TH1D			*hCuBoxSxth232M2Sum_1;
	TH1D			*hCuBoxSxth232M2Sum_10;
	TH1D			*hCuBoxSxu238M2Sum_001;
	TH1D			*hCuBoxSxu238M2Sum_01;
	TH1D			*hCuBoxSxu238M2Sum_1;
	TH1D			*hCuBoxSxu238M2Sum_10;

/////////// Frame + Box
	TH1D			*hCuBox_CuFrameco60M1;
	TH1D			*hCuBox_CuFramek40M1;
	TH1D			*hCuBox_CuFrameth232M1;
	TH1D			*hCuBox_CuFrameu238M1;

	TH1D			*hCuBox_CuFrameth232M1_10;
	TH1D			*hCuBox_CuFrameu238M1_10;
	TH1D			*hCuBox_CuFramepb210M1_10;
	TH1D			*hCuBox_CuFramepb210M1_1;
	TH1D			*hCuBox_CuFramepb210M1_01;
	TH1D			*hCuBox_CuFramepb210M1_001;

	TH1D			*hCuBox_CuFrameco60M2;
	TH1D			*hCuBox_CuFramek40M2;
	TH1D			*hCuBox_CuFrameth232M2;
	TH1D			*hCuBox_CuFrameu238M2;

	TH1D			*hCuBox_CuFrameth232M2_10;
	TH1D			*hCuBox_CuFrameu238M2_10;
	TH1D			*hCuBox_CuFramepb210M2_10;
	TH1D			*hCuBox_CuFramepb210M2_1;
	TH1D			*hCuBox_CuFramepb210M2_01;
	TH1D			*hCuBox_CuFramepb210M2_001;

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

///////////// 50mK M1 and M2
	TH1D			*h50mKco58M1;
	TH1D			*h50mKco60M1;
	TH1D			*h50mKcs137M1;
	TH1D			*h50mKk40M1;
	TH1D			*h50mKmn54M1;
	TH1D			*h50mKpb210M1;
	TH1D			*h50mKth232M1;
	TH1D			*h50mKu238M1;		

	TH1D			*h50mKco58M2;
	TH1D			*h50mKco60M2;
	TH1D			*h50mKcs137M2;
	TH1D			*h50mKk40M2;
	TH1D			*h50mKmn54M2;
	TH1D			*h50mKpb210M2;
	TH1D			*h50mKth232M2;
	TH1D			*h50mKu238M2;	

	TH1D			*h50mKco58M2Sum;
	TH1D			*h50mKco60M2Sum;
	TH1D			*h50mKcs137M2Sum;
	TH1D			*h50mKk40M2Sum;
	TH1D			*h50mKmn54M2Sum;
	TH1D			*h50mKpb210M2Sum;
	TH1D			*h50mKth232M2Sum;
	TH1D			*h50mKu238M2Sum;

//////////// 600mK M1 and M2
	TH1D			*h600mKco60M1;
	TH1D			*h600mKk40M1;
	TH1D			*h600mKth232M1;
	TH1D			*h600mKu238M1;		

	TH1D			*h600mKco60M2;
	TH1D			*h600mKk40M2;
	TH1D			*h600mKth232M2;
	TH1D			*h600mKu238M2;	

	TH1D			*h600mKco60M2Sum;
	TH1D			*h600mKk40M2Sum;
	TH1D			*h600mKth232M2Sum;
	TH1D			*h600mKu238M2Sum;	

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

/////////// Main Bath M1 and M2
	TH1D			*hMBco60M1;
	TH1D			*hMBk40M1;
	TH1D			*hMBth232M1;
	TH1D			*hMBu238M1;		

	TH1D			*hMBco60M2;
	TH1D			*hMBk40M2;
	TH1D			*hMBth232M2;
	TH1D			*hMBu238M2;	

	TH1D			*hMBco60M2Sum;
	TH1D			*hMBk40M2Sum;
	TH1D			*hMBth232M2Sum;
	TH1D			*hMBu238M2Sum;	

///////// Super Insulation M1 and M2
	TH1D			*hSIk40M1;
	TH1D			*hSIth232M1;
	TH1D			*hSIu238M1;

	TH1D			*hSIk40M2;
	TH1D			*hSIth232M2;
	TH1D			*hSIu238M2;

	TH1D			*hSIk40M2Sum;
	TH1D			*hSIth232M2Sum;
	TH1D			*hSIu238M2Sum;

//////////// External Shield M1 and M2
	TH1D			*hExtPbbi210M1;

	TH1D			*hExtPbbi210M2;
	
	TH1D			*hExtPbbi210M2Sum;

/////////// IVC M1 and M2
	TH1D			*hIVCco60M1;
	TH1D			*hIVCk40M1;
	TH1D			*hIVCth232M1;
	TH1D			*hIVCu238M1;		

	TH1D			*hIVCco60M2;
	TH1D			*hIVCk40M2;
	TH1D			*hIVCth232M2;
	TH1D			*hIVCu238M2;	

	TH1D			*hIVCco60M2Sum;
	TH1D			*hIVCk40M2Sum;
	TH1D			*hIVCth232M2Sum;
	TH1D			*hIVCu238M2Sum;	

/////////// OVC M1 and M2
	TH1D			*hOVCco60M1;
	TH1D			*hOVCk40M1;
	TH1D			*hOVCth232M1;
	TH1D			*hOVCu238M1;		

	TH1D			*hOVCco60M2;
	TH1D			*hOVCk40M2;
	TH1D			*hOVCth232M2;
	TH1D			*hOVCu238M2;	

	TH1D			*hOVCco60M2Sum;
	TH1D			*hOVCk40M2Sum;
	TH1D			*hOVCth232M2Sum;
	TH1D			*hOVCu238M2Sum;	

////////// Fudge Factors
	TH1D 		*hFudge661M1;
	TH1D 		*hFudge803M1;
	TH1D 		*hFudge1063M1;

	TH1D 		*hFudge661M1;
	TH1D 		*hFudge803M1;
	TH1D 		*hFudge1063M1;

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

	TH1D			*hAdapTeO2th232onlyM1;
	TH1D			*hAdapTeO2ra228pb208M1;
	TH1D			*hAdapTeO2th230onlyM1;

	TH1D			*hAdapTeO2Sxth232onlyM1_001;
	TH1D			*hAdapTeO2Sxra228pb208M1_001;
	TH1D			*hAdapTeO2Sxu238th230M1_001;
	TH1D			*hAdapTeO2Sxth230onlyM1_001;
	TH1D			*hAdapTeO2Sxra226pb210M1_001;
	TH1D			*hAdapTeO2Sxpb210M1_0001;

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

	TH1D			*hAdapTeO2th232onlyM2;
	TH1D			*hAdapTeO2ra228pb208M2;
	TH1D			*hAdapTeO2th230onlyM2;

	TH1D			*hAdapTeO2Sxth232onlyM2_001;
	TH1D			*hAdapTeO2Sxra228pb208M2_001;
	TH1D			*hAdapTeO2Sxu238th230M2_001;
	TH1D			*hAdapTeO2Sxth230onlyM2_001;
	TH1D			*hAdapTeO2Sxra226pb210M2_001;
	TH1D			*hAdapTeO2Sxpb210M2_0001;

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

	TH1D			*hAdapTeO2th232onlyM2Sum;
	TH1D			*hAdapTeO2ra228pb208M2Sum;
	TH1D			*hAdapTeO2th230onlyM2Sum;

	TH1D			*hAdapTeO2Sxth232onlyM2Sum_001;
	TH1D			*hAdapTeO2Sxra228pb208M2Sum_001;
	TH1D			*hAdapTeO2Sxu238th230M2Sum_001;
	TH1D			*hAdapTeO2Sxth230onlyM2Sum_001;
	TH1D			*hAdapTeO2Sxra226pb210M2Sum_001;
	TH1D			*hAdapTeO2Sxpb210M2Sum_0001;

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

//////////// Frame M1 and M2
	TH1D			*hAdapCuFrameco58M1;
	TH1D			*hAdapCuFrameco60M1;
	TH1D			*hAdapCuFramecs137M1;
	TH1D			*hAdapCuFramek40M1;
	TH1D			*hAdapCuFramemn54M1;
	TH1D			*hAdapCuFramepb210M1;
	TH1D			*hAdapCuFrameth232M1;
	TH1D			*hAdapCuFrameu238M1;

	TH1D			*hAdapCuFrameSth232M1_1;
	TH1D			*hAdapCuFrameSu238M1_1;
	TH1D			*hAdapCuFrameSxpb210M1_001;
	TH1D			*hAdapCuFrameSxpb210M1_01;
	TH1D			*hAdapCuFrameSxpb210M1_1;
	TH1D			*hAdapCuFrameSxpb210M1_10;
	TH1D			*hAdapCuFrameSxth232M1_001;
	TH1D			*hAdapCuFrameSxth232M1_01;
	TH1D			*hAdapCuFrameSxth232M1_1;
	TH1D			*hAdapCuFrameSxth232M1_10;
	TH1D			*hAdapCuFrameSxu238M1_001;
	TH1D			*hAdapCuFrameSxu238M1_01;
	TH1D			*hAdapCuFrameSxu238M1_1;
	TH1D			*hAdapCuFrameSxu238M1_10;

	TH1D			*hAdapCuFrameco58M2;
	TH1D			*hAdapCuFrameco60M2;
	TH1D			*hAdapCuFramecs137M2;
	TH1D			*hAdapCuFramek40M2;
	TH1D			*hAdapCuFramemn54M2;
	TH1D			*hAdapCuFramepb210M2;
	TH1D			*hAdapCuFrameth232M2;
	TH1D			*hAdapCuFrameu238M2;

	TH1D			*hAdapCuFrameSth232M2_1;
	TH1D			*hAdapCuFrameSu238M2_1;
	TH1D			*hAdapCuFrameSxpb210M2_001;
	TH1D			*hAdapCuFrameSxpb210M2_01;
	TH1D			*hAdapCuFrameSxpb210M2_1;
	TH1D			*hAdapCuFrameSxpb210M2_10;
	TH1D			*hAdapCuFrameSxth232M2_001;
	TH1D			*hAdapCuFrameSxth232M2_01;
	TH1D			*hAdapCuFrameSxth232M2_1;
	TH1D			*hAdapCuFrameSxth232M2_10;
	TH1D			*hAdapCuFrameSxu238M2_001;
	TH1D			*hAdapCuFrameSxu238M2_01;
	TH1D			*hAdapCuFrameSxu238M2_1;
	TH1D			*hAdapCuFrameSxu238M2_10;

	TH1D			*hAdapCuFrameco58M2Sum;
	TH1D			*hAdapCuFrameco60M2Sum;
	TH1D			*hAdapCuFramecs137M2Sum;
	TH1D			*hAdapCuFramek40M2Sum;
	TH1D			*hAdapCuFramemn54M2Sum;
	TH1D			*hAdapCuFramepb210M2Sum;
	TH1D			*hAdapCuFrameth232M2Sum;
	TH1D			*hAdapCuFrameu238M2Sum;

	TH1D			*hAdapCuFrameSth232M2Sum_1;
	TH1D			*hAdapCuFrameSu238M2Sum_1;
	TH1D			*hAdapCuFrameSxpb210M2Sum_001;
	TH1D			*hAdapCuFrameSxpb210M2Sum_01;
	TH1D			*hAdapCuFrameSxpb210M2Sum_1;
	TH1D			*hAdapCuFrameSxpb210M2Sum_10;
	TH1D			*hAdapCuFrameSxth232M2Sum_001;
	TH1D			*hAdapCuFrameSxth232M2Sum_01;
	TH1D			*hAdapCuFrameSxth232M2Sum_1;
	TH1D			*hAdapCuFrameSxth232M2Sum_10;
	TH1D			*hAdapCuFrameSxu238M2Sum_001;
	TH1D			*hAdapCuFrameSxu238M2Sum_01;
	TH1D			*hAdapCuFrameSxu238M2Sum_1;
	TH1D			*hAdapCuFrameSxu238M2Sum_10;

///////////// CuBox (TShield) M1 and M2
	TH1D			*hAdapCuBoxco58M1;
	TH1D			*hAdapCuBoxco60M1;
	TH1D			*hAdapCuBoxcs137M1;
	TH1D			*hAdapCuBoxk40M1;
	TH1D			*hAdapCuBoxmn54M1;
	TH1D			*hAdapCuBoxpb210M1;
	TH1D			*hAdapCuBoxth232M1;
	TH1D			*hAdapCuBoxu238M1;	

	TH1D			*hAdapCuBoxSth232M1_1;
	TH1D			*hAdapCuBoxSu238M1_1;
	TH1D			*hAdapCuBoxSxpb210M1_001;
	TH1D			*hAdapCuBoxSxpb210M1_01;
	TH1D			*hAdapCuBoxSxpb210M1_1;
	TH1D			*hAdapCuBoxSxpb210M1_10;
	TH1D			*hAdapCuBoxSxth232M1_001;
	TH1D			*hAdapCuBoxSxth232M1_01;
	TH1D			*hAdapCuBoxSxth232M1_1;
	TH1D			*hAdapCuBoxSxth232M1_10;
	TH1D			*hAdapCuBoxSxu238M1_001;
	TH1D			*hAdapCuBoxSxu238M1_01;
	TH1D			*hAdapCuBoxSxu238M1_1;
	TH1D			*hAdapCuBoxSxu238M1_10;

	TH1D			*hAdapCuBoxco58M2;
	TH1D			*hAdapCuBoxco60M2;
	TH1D			*hAdapCuBoxcs137M2;
	TH1D			*hAdapCuBoxk40M2;
	TH1D			*hAdapCuBoxmn54M2;
	TH1D			*hAdapCuBoxpb210M2;
	TH1D			*hAdapCuBoxth232M2;
	TH1D			*hAdapCuBoxu238M2;	

	TH1D			*hAdapCuBoxSth232M2_1;
	TH1D			*hAdapCuBoxSu238M2_1;
	TH1D			*hAdapCuBoxSxpb210M2_001;
	TH1D			*hAdapCuBoxSxpb210M2_01;
	TH1D			*hAdapCuBoxSxpb210M2_1;
	TH1D			*hAdapCuBoxSxpb210M2_10;
	TH1D			*hAdapCuBoxSxth232M2_001;
	TH1D			*hAdapCuBoxSxth232M2_01;
	TH1D			*hAdapCuBoxSxth232M2_1;
	TH1D			*hAdapCuBoxSxth232M2_10;
	TH1D			*hAdapCuBoxSxu238M2_001;
	TH1D			*hAdapCuBoxSxu238M2_01;
	TH1D			*hAdapCuBoxSxu238M2_1;
	TH1D			*hAdapCuBoxSxu238M2_10;

	TH1D			*hAdapCuBoxco58M2Sum;
	TH1D			*hAdapCuBoxco60M2Sum;
	TH1D			*hAdapCuBoxcs137M2Sum;
	TH1D			*hAdapCuBoxk40M2Sum;
	TH1D			*hAdapCuBoxmn54M2Sum;
	TH1D			*hAdapCuBoxpb210M2Sum;
	TH1D			*hAdapCuBoxth232M2Sum;
	TH1D			*hAdapCuBoxu238M2Sum;	

	TH1D			*hAdapCuBoxSth232M2Sum_1;
	TH1D			*hAdapCuBoxSu238M2Sum_1;
	TH1D			*hAdapCuBoxSxpb210M2Sum_001;
	TH1D			*hAdapCuBoxSxpb210M2Sum_01;
	TH1D			*hAdapCuBoxSxpb210M2Sum_1;
	TH1D			*hAdapCuBoxSxpb210M2Sum_10;
	TH1D			*hAdapCuBoxSxth232M2Sum_001;
	TH1D			*hAdapCuBoxSxth232M2Sum_01;
	TH1D			*hAdapCuBoxSxth232M2Sum_1;
	TH1D			*hAdapCuBoxSxth232M2Sum_10;
	TH1D			*hAdapCuBoxSxu238M2Sum_001;
	TH1D			*hAdapCuBoxSxu238M2Sum_01;
	TH1D			*hAdapCuBoxSxu238M2Sum_1;
	TH1D			*hAdapCuBoxSxu238M2Sum_10;

////////////// Frame + Box
	TH1D			*hAdapCuBox_CuFrameco60M1;
	TH1D			*hAdapCuBox_CuFramek40M1;
	TH1D			*hAdapCuBox_CuFrameth232M1;
	TH1D			*hAdapCuBox_CuFrameu238M1;

	TH1D			*hAdapCuBox_CuFrameth232M1_10;
	TH1D			*hAdapCuBox_CuFrameu238M1_10;

	TH1D			*hAdapCuBox_CuFramepb210M1_10;
	TH1D			*hAdapCuBox_CuFramepb210M1_1;
	TH1D			*hAdapCuBox_CuFramepb210M1_01;
	TH1D			*hAdapCuBox_CuFramepb210M1_001;

	TH1D			*hAdapCuBox_CuFrameco60M2;
	TH1D			*hAdapCuBox_CuFramek40M2;
	TH1D			*hAdapCuBox_CuFrameth232M2;
	TH1D			*hAdapCuBox_CuFrameu238M2;

	TH1D			*hAdapCuBox_CuFrameth232M2_10;
	TH1D			*hAdapCuBox_CuFrameu238M2_10;
	TH1D			*hAdapCuBox_CuFramepb210M2_10;
	TH1D			*hAdapCuBox_CuFramepb210M2_1;
	TH1D			*hAdapCuBox_CuFramepb210M2_01;
	TH1D			*hAdapCuBox_CuFramepb210M2_001;

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


/////////// 50mK M1 and M2
	TH1D			*hAdap50mKco58M1;
	TH1D			*hAdap50mKco60M1;
	TH1D			*hAdap50mKcs137M1;
	TH1D			*hAdap50mKk40M1;
	TH1D			*hAdap50mKmn54M1;
	TH1D			*hAdap50mKpb210M1;
	TH1D			*hAdap50mKth232M1;
	TH1D			*hAdap50mKu238M1;		

	TH1D			*hAdap50mKco58M2;
	TH1D			*hAdap50mKco60M2;
	TH1D			*hAdap50mKcs137M2;
	TH1D			*hAdap50mKk40M2;
	TH1D			*hAdap50mKmn54M2;
	TH1D			*hAdap50mKpb210M2;
	TH1D			*hAdap50mKth232M2;
	TH1D			*hAdap50mKu238M2;	

	TH1D			*hAdap50mKco58M2Sum;
	TH1D			*hAdap50mKco60M2Sum;
	TH1D			*hAdap50mKcs137M2Sum;
	TH1D			*hAdap50mKk40M2Sum;
	TH1D			*hAdap50mKmn54M2Sum;
	TH1D			*hAdap50mKpb210M2Sum;
	TH1D			*hAdap50mKth232M2Sum;
	TH1D			*hAdap50mKu238M2Sum;	

///////// 600mK M1 and M2
	TH1D			*hAdap600mKco60M1;
	TH1D			*hAdap600mKk40M1;
	TH1D			*hAdap600mKth232M1;
	TH1D			*hAdap600mKu238M1;		

	TH1D			*hAdap600mKco60M2;
	TH1D			*hAdap600mKk40M2;
	TH1D			*hAdap600mKth232M2;
	TH1D			*hAdap600mKu238M2;	

	TH1D			*hAdap600mKco60M2Sum;
	TH1D			*hAdap600mKk40M2Sum;
	TH1D			*hAdap600mKth232M2Sum;
	TH1D			*hAdap600mKu238M2Sum;	

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

///////// Main Bath M1 and M2
	TH1D			*hAdapMBco60M1;
	TH1D			*hAdapMBk40M1;
	TH1D			*hAdapMBth232M1;
	TH1D			*hAdapMBu238M1;		

	TH1D			*hAdapMBco60M2;
	TH1D			*hAdapMBk40M2;
	TH1D			*hAdapMBth232M2;
	TH1D			*hAdapMBu238M2;	

	TH1D			*hAdapMBco60M2Sum;
	TH1D			*hAdapMBk40M2Sum;
	TH1D			*hAdapMBth232M2Sum;
	TH1D			*hAdapMBu238M2Sum;	

///////// Super Insulation M1 and M2
	TH1D			*hAdapSIk40M1;
	TH1D			*hAdapSIth232M1;
	TH1D			*hAdapSIu238M1;

	TH1D			*hAdapSIk40M2;
	TH1D			*hAdapSIth232M2;
	TH1D			*hAdapSIu238M2;

	TH1D			*hAdapSIk40M2Sum;
	TH1D			*hAdapSIth232M2Sum;
	TH1D			*hAdapSIu238M2Sum;

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

/////////// IVC M1 and M2
	TH1D			*hAdapIVCco60M1;
	TH1D			*hAdapIVCk40M1;
	TH1D			*hAdapIVCth232M1;
	TH1D			*hAdapIVCu238M1;		

	TH1D			*hAdapIVCco60M2;
	TH1D			*hAdapIVCk40M2;
	TH1D			*hAdapIVCth232M2;
	TH1D			*hAdapIVCu238M2;	

	TH1D			*hAdapIVCco60M2Sum;
	TH1D			*hAdapIVCk40M2Sum;
	TH1D			*hAdapIVCth232M2Sum;
	TH1D			*hAdapIVCu238M2Sum;	

///////// OVC M1 and M2
	TH1D			*hAdapOVCco60M1;
	TH1D			*hAdapOVCk40M1;
	TH1D			*hAdapOVCth232M1;
	TH1D			*hAdapOVCu238M1;		

	TH1D			*hAdapOVCco60M2;
	TH1D			*hAdapOVCk40M2;
	TH1D			*hAdapOVCth232M2;
	TH1D			*hAdapOVCu238M2;	

	TH1D			*hAdapOVCco60M2Sum;
	TH1D			*hAdapOVCk40M2Sum;
	TH1D			*hAdapOVCth232M2Sum;
	TH1D			*hAdapOVCu238M2Sum;

///////// External Shield M1 and M2
	TH1D			*hAdapExtPbbi210M1;

	TH1D			*hAdapExtPbbi210M2;
	
	TH1D			*hAdapExtPbbi210M2Sum;

////////// Fudge Factors
	TH1D 		*hAdapFudge661M1;
	TH1D 		*hAdapFudge803M1;
	TH1D 		*hAdapFudge1063M1;

	TH1 		*hAdapFudge661M2;
	TH1 		*hAdapFudge803M2;
	TH1 		*hAdapFudge1063M2;


	TH1D			*hEnergyScaleDummyM1;
	TH1D			*hEnergyScaleDummyM2;
	TH1D			*hEnergyScaleDummyM2Sum;


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

	TH1			*hnewTeO2th232onlyM1;
	TH1			*hnewTeO2ra228pb208M1;
	TH1			*hnewTeO2th230onlyM1;

	TH1			*hnewTeO2Sxth232onlyM1_001;
	TH1			*hnewTeO2Sxra228pb208M1_001;
	TH1			*hnewTeO2Sxu238th230M1_001;
	TH1			*hnewTeO2Sxth230onlyM1_001;
	TH1			*hnewTeO2Sxra226pb210M1_001;
	TH1			*hnewTeO2Sxpb210M1_0001;

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

	TH1			*hnewTeO2th232onlyM2;
	TH1			*hnewTeO2ra228pb208M2;
	TH1			*hnewTeO2th230onlyM2;

	TH1			*hnewTeO2Sxth232onlyM2_001;
	TH1			*hnewTeO2Sxra228pb208M2_001;
	TH1			*hnewTeO2Sxu238th230M2_001;
	TH1			*hnewTeO2Sxth230onlyM2_001;
	TH1			*hnewTeO2Sxra226pb210M2_001;
	TH1			*hnewTeO2Sxpb210M2_0001;


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

	TH1			*hnewTeO2th232onlyM2Sum;
	TH1			*hnewTeO2ra228pb208M2Sum;
	TH1			*hnewTeO2th230onlyM2Sum;

	TH1			*hnewTeO2Sxth232onlyM2Sum_001;
	TH1			*hnewTeO2Sxra228pb208M2Sum_001;
	TH1			*hnewTeO2Sxu238th230M2Sum_001;
	TH1			*hnewTeO2Sxth230onlyM2Sum_001;
	TH1			*hnewTeO2Sxra226pb210M2Sum_001;
	TH1			*hnewTeO2Sxpb210M2Sum_0001;

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

///////// Frame M1 and M2
	TH1			*hnewCuFrameco58M1;
	TH1			*hnewCuFrameco60M1;
	TH1			*hnewCuFramecs137M1;
	TH1			*hnewCuFramek40M1;
	TH1			*hnewCuFramemn54M1;
	TH1			*hnewCuFramepb210M1;
	TH1			*hnewCuFrameth232M1;
	TH1			*hnewCuFrameu238M1;

	TH1			*hnewCuFrameSth232M1_1;
	TH1			*hnewCuFrameSu238M1_1;
	TH1			*hnewCuFrameSxpb210M1_001;
	TH1			*hnewCuFrameSxpb210M1_01;
	TH1			*hnewCuFrameSxpb210M1_1;
	TH1			*hnewCuFrameSxpb210M1_10;
	TH1			*hnewCuFrameSxth232M1_001;
	TH1			*hnewCuFrameSxth232M1_01;
	TH1			*hnewCuFrameSxth232M1_1;
	TH1			*hnewCuFrameSxth232M1_10;
	TH1			*hnewCuFrameSxu238M1_001;
	TH1			*hnewCuFrameSxu238M1_01;
	TH1			*hnewCuFrameSxu238M1_1;
	TH1			*hnewCuFrameSxu238M1_10;

	TH1			*hnewCuFrameco58M2;
	TH1			*hnewCuFrameco60M2;
	TH1			*hnewCuFramecs137M2;
	TH1			*hnewCuFramek40M2;
	TH1			*hnewCuFramemn54M2;
	TH1			*hnewCuFramepb210M2;
	TH1			*hnewCuFrameth232M2;
	TH1			*hnewCuFrameu238M2;

	TH1			*hnewCuFrameSth232M2_1;
	TH1			*hnewCuFrameSu238M2_1;
	TH1			*hnewCuFrameSxpb210M2_001;
	TH1			*hnewCuFrameSxpb210M2_01;
	TH1			*hnewCuFrameSxpb210M2_1;
	TH1			*hnewCuFrameSxpb210M2_10;
	TH1			*hnewCuFrameSxth232M2_001;
	TH1			*hnewCuFrameSxth232M2_01;
	TH1			*hnewCuFrameSxth232M2_1;
	TH1			*hnewCuFrameSxth232M2_10;
	TH1			*hnewCuFrameSxu238M2_001;
	TH1			*hnewCuFrameSxu238M2_01;
	TH1			*hnewCuFrameSxu238M2_1;
	TH1			*hnewCuFrameSxu238M2_10;

	TH1			*hnewCuFrameco58M2Sum;
	TH1			*hnewCuFrameco60M2Sum;
	TH1			*hnewCuFramecs137M2Sum;
	TH1			*hnewCuFramek40M2Sum;
	TH1			*hnewCuFramemn54M2Sum;
	TH1			*hnewCuFramepb210M2Sum;
	TH1			*hnewCuFrameth232M2Sum;
	TH1			*hnewCuFrameu238M2Sum;

	TH1			*hnewCuFrameSth232M2Sum_1;
	TH1			*hnewCuFrameSu238M2Sum_1;
	TH1			*hnewCuFrameSxpb210M2Sum_001;
	TH1			*hnewCuFrameSxpb210M2Sum_01;
	TH1			*hnewCuFrameSxpb210M2Sum_1;
	TH1			*hnewCuFrameSxpb210M2Sum_10;
	TH1			*hnewCuFrameSxth232M2Sum_001;
	TH1			*hnewCuFrameSxth232M2Sum_01;
	TH1			*hnewCuFrameSxth232M2Sum_1;
	TH1			*hnewCuFrameSxth232M2Sum_10;
	TH1			*hnewCuFrameSxu238M2Sum_001;
	TH1			*hnewCuFrameSxu238M2Sum_01;
	TH1			*hnewCuFrameSxu238M2Sum_1;
	TH1			*hnewCuFrameSxu238M2Sum_10;

/////////// CuBox (TShield) M1 and M2
	TH1			*hnewCuBoxco58M1;
	TH1			*hnewCuBoxco60M1;
	TH1			*hnewCuBoxcs137M1;
	TH1			*hnewCuBoxk40M1;
	TH1			*hnewCuBoxmn54M1;
	TH1			*hnewCuBoxpb210M1;
	TH1			*hnewCuBoxth232M1;
	TH1			*hnewCuBoxu238M1;	

	TH1			*hnewCuBoxSth232M1_1;
	TH1			*hnewCuBoxSu238M1_1;
	TH1			*hnewCuBoxSxpb210M1_001;
	TH1			*hnewCuBoxSxpb210M1_01;
	TH1			*hnewCuBoxSxpb210M1_1;
	TH1			*hnewCuBoxSxpb210M1_10;
	TH1			*hnewCuBoxSxth232M1_001;
	TH1			*hnewCuBoxSxth232M1_01;
	TH1			*hnewCuBoxSxth232M1_1;
	TH1			*hnewCuBoxSxth232M1_10;
	TH1			*hnewCuBoxSxu238M1_001;
	TH1			*hnewCuBoxSxu238M1_01;
	TH1			*hnewCuBoxSxu238M1_1;
	TH1			*hnewCuBoxSxu238M1_10;

	TH1			*hnewCuBoxco58M2;
	TH1			*hnewCuBoxco60M2;
	TH1			*hnewCuBoxcs137M2;
	TH1			*hnewCuBoxk40M2;
	TH1			*hnewCuBoxmn54M2;
	TH1			*hnewCuBoxpb210M2;
	TH1			*hnewCuBoxth232M2;
	TH1			*hnewCuBoxu238M2;	

	TH1			*hnewCuBoxSth232M2_1;
	TH1			*hnewCuBoxSu238M2_1;
	TH1			*hnewCuBoxSxpb210M2_001;
	TH1			*hnewCuBoxSxpb210M2_01;
	TH1			*hnewCuBoxSxpb210M2_1;
	TH1			*hnewCuBoxSxpb210M2_10;
	TH1			*hnewCuBoxSxth232M2_001;
	TH1			*hnewCuBoxSxth232M2_01;
	TH1			*hnewCuBoxSxth232M2_1;
	TH1			*hnewCuBoxSxth232M2_10;
	TH1			*hnewCuBoxSxu238M2_001;
	TH1			*hnewCuBoxSxu238M2_01;
	TH1			*hnewCuBoxSxu238M2_1;
	TH1			*hnewCuBoxSxu238M2_10;

	TH1			*hnewCuBoxco58M2Sum;
	TH1			*hnewCuBoxco60M2Sum;
	TH1			*hnewCuBoxcs137M2Sum;
	TH1			*hnewCuBoxk40M2Sum;
	TH1			*hnewCuBoxmn54M2Sum;
	TH1			*hnewCuBoxpb210M2Sum;
	TH1			*hnewCuBoxth232M2Sum;
	TH1			*hnewCuBoxu238M2Sum;	

	TH1			*hnewCuBoxSth232M2Sum_1;
	TH1			*hnewCuBoxSu238M2Sum_1;
	TH1			*hnewCuBoxSxpb210M2Sum_001;
	TH1			*hnewCuBoxSxpb210M2Sum_01;
	TH1			*hnewCuBoxSxpb210M2Sum_1;
	TH1			*hnewCuBoxSxpb210M2Sum_10;
	TH1			*hnewCuBoxSxth232M2Sum_001;
	TH1			*hnewCuBoxSxth232M2Sum_01;
	TH1			*hnewCuBoxSxth232M2Sum_1;
	TH1			*hnewCuBoxSxth232M2Sum_10;
	TH1			*hnewCuBoxSxu238M2Sum_001;
	TH1			*hnewCuBoxSxu238M2Sum_01;
	TH1			*hnewCuBoxSxu238M2Sum_1;
	TH1			*hnewCuBoxSxu238M2Sum_10;

/////////// Frame + Box
	TH1			*hnewCuBox_CuFrameco60M1;
	TH1			*hnewCuBox_CuFramek40M1;
	TH1			*hnewCuBox_CuFrameth232M1;
	TH1			*hnewCuBox_CuFrameu238M1;

	TH1			*hnewCuBox_CuFrameth232M1_10;
	TH1			*hnewCuBox_CuFrameu238M1_10;
	TH1			*hnewCuBox_CuFramepb210M1_10;
	TH1			*hnewCuBox_CuFramepb210M1_1;
	TH1			*hnewCuBox_CuFramepb210M1_01;
	TH1			*hnewCuBox_CuFramepb210M1_001;

	TH1			*hnewCuBox_CuFrameco60M2;
	TH1			*hnewCuBox_CuFramek40M2;
	TH1			*hnewCuBox_CuFrameth232M2;
	TH1			*hnewCuBox_CuFrameu238M2;

	TH1			*hnewCuBox_CuFrameth232M2_10;
	TH1			*hnewCuBox_CuFrameu238M2_10;
	TH1			*hnewCuBox_CuFramepb210M2_10;
	TH1			*hnewCuBox_CuFramepb210M2_1;
	TH1			*hnewCuBox_CuFramepb210M2_01;
	TH1			*hnewCuBox_CuFramepb210M2_001;

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

///////////// 50mK M1 and M2
	TH1			*hnew50mKco58M1;
	TH1			*hnew50mKco60M1;
	TH1			*hnew50mKcs137M1;
	TH1			*hnew50mKk40M1;
	TH1			*hnew50mKmn54M1;
	TH1			*hnew50mKpb210M1;
	TH1			*hnew50mKth232M1;
	TH1			*hnew50mKu238M1;		

	TH1			*hnew50mKco58M2;
	TH1			*hnew50mKco60M2;
	TH1			*hnew50mKcs137M2;
	TH1			*hnew50mKk40M2;
	TH1			*hnew50mKmn54M2;
	TH1			*hnew50mKpb210M2;
	TH1			*hnew50mKth232M2;
	TH1			*hnew50mKu238M2;	

	TH1			*hnew50mKco58M2Sum;
	TH1			*hnew50mKco60M2Sum;
	TH1			*hnew50mKcs137M2Sum;
	TH1			*hnew50mKk40M2Sum;
	TH1			*hnew50mKmn54M2Sum;
	TH1			*hnew50mKpb210M2Sum;
	TH1			*hnew50mKth232M2Sum;
	TH1			*hnew50mKu238M2Sum;

//////////// 600mK M1 and M2
	TH1			*hnew600mKco60M1;
	TH1			*hnew600mKk40M1;
	TH1			*hnew600mKth232M1;
	TH1			*hnew600mKu238M1;		

	TH1			*hnew600mKco60M2;
	TH1			*hnew600mKk40M2;
	TH1			*hnew600mKth232M2;
	TH1			*hnew600mKu238M2;	

	TH1			*hnew600mKco60M2Sum;
	TH1			*hnew600mKk40M2Sum;
	TH1			*hnew600mKth232M2Sum;
	TH1			*hnew600mKu238M2Sum;	

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

/////////// Main Bath M1 and M2
	TH1			*hnewMBco60M1;
	TH1			*hnewMBk40M1;
	TH1			*hnewMBth232M1;
	TH1			*hnewMBu238M1;		

	TH1			*hnewMBco60M2;
	TH1			*hnewMBk40M2;
	TH1			*hnewMBth232M2;
	TH1			*hnewMBu238M2;	

	TH1			*hnewMBco60M2Sum;
	TH1			*hnewMBk40M2Sum;
	TH1			*hnewMBth232M2Sum;
	TH1			*hnewMBu238M2Sum;	

///////// Super Insulation M1 and M2
	TH1			*hnewSIk40M1;
	TH1			*hnewSIth232M1;
	TH1			*hnewSIu238M1;

	TH1			*hnewSIk40M2;
	TH1			*hnewSIth232M2;
	TH1			*hnewSIu238M2;

	TH1			*hnewSIk40M2Sum;
	TH1			*hnewSIth232M2Sum;
	TH1			*hnewSIu238M2Sum;

//////////// External Shield M1 and M2
	TH1			*hnewExtPbbi210M1;

	TH1			*hnewExtPbbi210M2;
	
	TH1			*hnewExtPbbi210M2Sum;

/////////// IVC M1 and M2
	TH1			*hnewIVCco60M1;
	TH1			*hnewIVCk40M1;
	TH1			*hnewIVCth232M1;
	TH1			*hnewIVCu238M1;		

	TH1			*hnewIVCco60M2;
	TH1			*hnewIVCk40M2;
	TH1			*hnewIVCth232M2;
	TH1			*hnewIVCu238M2;	

	TH1			*hnewIVCco60M2Sum;
	TH1			*hnewIVCk40M2Sum;
	TH1			*hnewIVCth232M2Sum;
	TH1			*hnewIVCu238M2Sum;	

/////////// OVC M1 and M2
	TH1			*hnewOVCco60M1;
	TH1			*hnewOVCk40M1;
	TH1			*hnewOVCth232M1;
	TH1			*hnewOVCu238M1;		

	TH1			*hnewOVCco60M2;
	TH1			*hnewOVCk40M2;
	TH1			*hnewOVCth232M2;
	TH1			*hnewOVCu238M2;	

	TH1			*hnewOVCco60M2Sum;
	TH1			*hnewOVCk40M2Sum;
	TH1			*hnewOVCth232M2Sum;
	TH1			*hnewOVCu238M2Sum;	


////////// Fudge Factors
	TH1 		*hnewFudge661M1;
	TH1 		*hnewFudge803M1;
	TH1 		*hnewFudge1063M1;

	TH1 		*hnewFudge661M2;
	TH1 		*hnewFudge803M2;
	TH1 		*hnewFudge1063M2;


	TH1D 			*hOut;


	TDatime 		*tTime;

	// Smearing
	TF1				*gaus;

	// Cut Efficiency
	TF1 			*fEfficiency;
	TH1 			*hEfficiencyM1;
	TH1 			*hEfficiencyM2;
	// TH1 			*hEfficiencyM1;

	ofstream 		OutFile;
	ofstream 	 	OutPNLL;

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

	std::string		dDataDir;
	std::string 	dMCDir;
	std::string 	dSaveDir;


	// Error Matrix
	// TMatrixT<double> 	*mCorrMatrix;

	bool			bFixedRes;
	bool			bAdaptiveBinning;
	bool 			bSave;

	int 			dNumCalls;
	int 			dMult;
	double			dBestChiSq;

	// Parameters
	double				fParameters[139];
	double				fParError[139];
	double 				fParActivity[139]; // Integral of all events
	double				fResolution[52];
	double 				fParEfficiencyM1[139]; // Efficiency of the parameters 
	double				dSecToYears;
	double				fMCEff[62];

 ClassDef(TBackgroundModel,1) // 
    };

#endif // __TBackgroundModel__
