#ifndef __TBackgroundModel__
#define __TBackgroundModel__
#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TChain.h"
#include "TCut.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include <vector>

class TBackgroundModel : public TObject {

public:

	TBackgroundModel(double fFitMin, double fFitMax);
	virtual ~TBackgroundModel();

	vector<double> AdaptiveBinning(TH1D *h1);

	TGraphErrors *CalculateResiduals(TH1D *h1, TH1D *h2, TH1D *hResid);

	TGraphErrors *CalculateResidualsAdaptive(TH1D *h1, TH1D *h2, TH1D *hResid, int binMin, int binMax);
  
	double GetChiSquare();

	double GetChiSquareAdaptive();

	bool DoTheFit();

	bool DoTheFitAdaptive();

	void DrawBkg();

	void DrawMC();

	void Initialize();

	// Dumb to have all of these but w/e
	void LoadData();

	void LoadPDFs();

	TChain *LoadMC(std::string dDir, std::string dLocation, std::string dSource, std::string dSType, int dMult);

	void NormalizePDF(TH1D *h1, int minE, int maxE);

	void NormalizePDFPair(TH1D *h1, TH1D *h2, int minE, int maxE);

	void PrintParameters();

	void ReadMC();

	void SetParameters(int index, double value);

	TH1D *SmearMC(TH1D *hMC, TH1D *hSMC, double resolution1, double resolution2);

	TH1D *SmearMCOld(TH1D *hMC, TH1D *hSMC, double resolution1);
	
	void Test();

	void UpdateModel();

	void UpdateModelAdaptive();


	int 	dBinSize;
	int 	dNBins;
	double	dMinEnergy;
	double 	dMaxEnergy;
	double	dFitMin;
	double	dFitMax;
	int 	dFitMinBinM1;
	int 	dFitMaxBinM1;
	int 	dFitMinBinM2;
	int 	dFitMaxBinM2;

	int 	dAdaptiveBinsM1;
	int 	dAdaptiveBinsM2;
	vector<double> dAdaptiveVectorM1;
	vector<double> dAdaptiveVectorM2;
	double 	*dAdaptiveArrayM1;
	double 	*dAdaptiveArrayM2;

private:

	// Data
	TChain			*qtree;
	TCut 			base_cut;
	TCut			ener_cut;

	double			dDataIntegral;

	TH1D			*fDataHistoTot;
	TH1D			*fDataHistoM1;
	TH1D			*fDataHistoM2;

	TH1D			*fAdapDataHistoM1;
	TH1D			*fAdapDataHistoM2;


	// Total PDFs M1
	TH1D			*fModelTotM1;
	TH1D			*fModelTotThM1;
	TH1D			*fModelTotRaM1;
	TH1D			*fModelTotKM1;
	TH1D			*fModelTotCoM1;
	TH1D			*fModelTotMnM1;
	TH1D			*fModelTotNDBDM1;
	TH1D			*fModelTot2NDBDM1;
	TH1D			*fModelTotBiM1;
	TH1D			*fModelTotPtM1;
	TH1D			*fModelTotPbM1;

	TH1D			*fModelTotAdapM1;
	TH1D			*fModelTotAdapThM1;
	TH1D			*fModelTotAdapRaM1;
	TH1D			*fModelTotAdapKM1;
	TH1D			*fModelTotAdapCoM1;
	TH1D			*fModelTotAdapMnM1;
	TH1D			*fModelTotAdapNDBDM1;
	TH1D			*fModelTotAdap2NDBDM1;
	TH1D			*fModelTotAdapBiM1;
	TH1D			*fModelTotAdapPtM1;
	TH1D			*fModelTotAdapPbM1;

	// Total PDFs M2
	TH1D			*fModelTotM2;
	TH1D			*fModelTotThM2;
	TH1D			*fModelTotRaM2;
	TH1D			*fModelTotKM2;
	TH1D			*fModelTotCoM2;
	TH1D			*fModelTotMnM2;
	TH1D			*fModelTotNDBDM2;
	TH1D			*fModelTot2NDBDM2;
	TH1D			*fModelTotBiM2;
	TH1D			*fModelTotPtM2;
	TH1D			*fModelTotPbM2;

	TH1D			*fModelTotAdapM2;
	TH1D			*fModelTotAdapThM2;
	TH1D			*fModelTotAdapRaM2;
	TH1D			*fModelTotAdapKM2;
	TH1D			*fModelTotAdapCoM2;
	TH1D			*fModelTotAdapMnM2;
	TH1D			*fModelTotAdapNDBDM2;
	TH1D			*fModelTotAdap2NDBDM2;
	TH1D			*fModelTotAdapBiM2;
	TH1D			*fModelTotAdapPtM2;
	TH1D			*fModelTotAdapPbM2;


	// Residuals
	TGraph			*gResidualM1;
	TGraph 			*gResidualM2;

	TH1D			*hResidualDistM1;
	TH1D			*hResidualDistM2;

	TH1D			*hResidualGausM1;
	TH1D			*hResidualGausM2;

	// Smeared PDFs M1
	TH1D			*fModelDummyM1;

	TH1D			*fModelNDBDM1;
	TH1D			*fModel2NDBDM1;
	TH1D			*fModelBiM1;

	TH1D			*fModelFrameThS01M1;
	TH1D			*fModelFrameThS1M1;
	TH1D			*fModelFrameThS10M1;
	TH1D			*fModelFrameThS100M1;

	TH1D			*fModelFrameRaS01M1;
	TH1D			*fModelFrameRaS1M1;
	TH1D			*fModelFrameRaS10M1;
	TH1D			*fModelFrameRaS100M1;

	TH1D			*fModelTShieldThS01M1;
	TH1D			*fModelTShieldThS1M1;
	TH1D			*fModelTShieldThS10M1;
	TH1D			*fModelTShieldThS100M1;


	TH1D			*fModelFrameThM1;
	TH1D			*fModelTShieldThM1;
	TH1D			*fModel50mKThM1;
	TH1D			*fModel600mKThM1;	
	TH1D			*fModelIVCThM1;
	TH1D			*fModelOVCThM1;

	TH1D			*fModelFrameRaM1;
	TH1D			*fModelTShieldRaM1;
	TH1D			*fModel50mKRaM1;
	TH1D			*fModel600mKRaM1;
	TH1D			*fModelIVCRaM1;
	TH1D			*fModelOVCRaM1;

	TH1D			*fModelFrameKM1;
	TH1D			*fModelTShieldKM1;
	TH1D			*fModel50mKKM1;
	TH1D			*fModel600mKKM1;
	TH1D			*fModelIVCKM1;
	TH1D			*fModelOVCKM1;


	TH1D			*fModelFrameCoM1;
	TH1D			*fModelTShieldCoM1;
	TH1D			*fModel50mKCoM1;
	TH1D			*fModel600mKCoM1;
	TH1D			*fModelIVCCoM1;
	TH1D			*fModelOVCCoM1;

	TH1D			*fModelTShieldMnM1;
	TH1D			*fModelIVCMnM1;

	TH1D			*fModelCrystalBi2M1;
	TH1D			*fModelFrameBi2M1;

	TH1D			*fModelCrystalPtM1;
	TH1D			*fModelCrystalPbBM1;
	TH1D			*fModelCrystalPbS01M1;
	TH1D			*fModelCrystalPbS1M1;
	TH1D			*fModelCrystalPbS10M1;
	TH1D			*fModelCrystalPbS100M1;
	TH1D			*fModelFramePbBM1;
	TH1D			*fModelFramePbS01M1;
	TH1D			*fModelFramePbS1M1;
	TH1D			*fModelFramePbS10M1;
	TH1D			*fModelFramePbS100M1;


	// Modeled PDFs M2
	TH1D			*fModelDummyM2;

	TH1D			*fModelNDBDM2;
	TH1D			*fModel2NDBDM2;
	TH1D			*fModelBiM2;

	TH1D			*fModelFrameThS01M2;
	TH1D			*fModelFrameThS1M2;
	TH1D			*fModelFrameThS10M2;
	TH1D			*fModelFrameThS100M2;

	TH1D			*fModelFrameRaS01M2;
	TH1D			*fModelFrameRaS1M2;
	TH1D			*fModelFrameRaS10M2;
	TH1D			*fModelFrameRaS100M2;

	TH1D			*fModelTShieldThS01M2;
	TH1D			*fModelTShieldThS1M2;
	TH1D			*fModelTShieldThS10M2;
	TH1D			*fModelTShieldThS100M2;


	TH1D			*fModelFrameThM2;
	TH1D			*fModelTShieldThM2;
	TH1D			*fModel50mKThM2;
	TH1D			*fModel600mKThM2;	
	TH1D			*fModelIVCThM2;
	TH1D			*fModelOVCThM2;

	TH1D			*fModelFrameRaM2;
	TH1D			*fModelTShieldRaM2;
	TH1D			*fModel50mKRaM2;
	TH1D			*fModel600mKRaM2;
	TH1D			*fModelIVCRaM2;
	TH1D			*fModelOVCRaM2;

	TH1D			*fModelFrameKM2;
	TH1D			*fModelTShieldKM2;
	TH1D			*fModel50mKKM2;
	TH1D			*fModel600mKKM2;
	TH1D			*fModelIVCKM2;
	TH1D			*fModelOVCKM2;


	TH1D			*fModelFrameCoM2;
	TH1D			*fModelTShieldCoM2;
	TH1D			*fModel50mKCoM2;
	TH1D			*fModel600mKCoM2;
	TH1D			*fModelIVCCoM2;
	TH1D			*fModelOVCCoM2;

	TH1D			*fModelTShieldMnM2;
	TH1D			*fModelIVCMnM2;

	TH1D			*fModelCrystalBi2M2;
	TH1D			*fModelFrameBi2M2;

	TH1D			*fModelCrystalPtM2;
	TH1D			*fModelCrystalPbBM2;
	TH1D			*fModelCrystalPbS01M2;
	TH1D			*fModelCrystalPbS1M2;
	TH1D			*fModelCrystalPbS10M2;
	TH1D			*fModelCrystalPbS100M2;
	TH1D			*fModelFramePbBM2;
	TH1D			*fModelFramePbS01M2;
	TH1D			*fModelFramePbS1M2;
	TH1D			*fModelFramePbS10M2;
	TH1D			*fModelFramePbS100M2;



	// M1 Alphas
	TH1D			*fModelCrystalPt190BM1;
	TH1D			*fModelCrystalPt190S01M1;
	TH1D			*fModelCrystalPt190S1M1;
	TH1D			*fModelCrystalPt190S10M1;
	TH1D			*fModelCrystalPt190S100M1;

	TH1D			*fModelCrystalPb210BM1;
	TH1D			*fModelCrystalPb210S01M1;
	TH1D			*fModelCrystalPb210S1M1;
	TH1D			*fModelCrystalPb210S10M1;
	TH1D			*fModelCrystalPb210S100M1;

	TH1D			*fModelCrystalRa226BM1;
	TH1D			*fModelCrystalRa226S01M1;
	TH1D			*fModelCrystalRa226S1M1;
	TH1D			*fModelCrystalRa226S10M1;
	TH1D			*fModelCrystalRa226S100M1;

	TH1D			*fModelFramePb210BM1;
	TH1D			*fModelFramePb210S01M1;
	TH1D			*fModelFramePb210S1M1;
	TH1D			*fModelFramePb210S10M1;
	TH1D			*fModelFramePb210S100M1;

	// TH1D			*fModelFrameRa226BM1;
	TH1D			*fModelFrameRa226S01M1;
	TH1D			*fModelFrameRa226S1M1;
	TH1D			*fModelFrameRa226S10M1;
	TH1D			*fModelFrameRa226S100M1;

	// M2 Alphas
	TH1D			*fModelCrystalPt190BM2;
	TH1D			*fModelCrystalPt190S01M2;
	TH1D			*fModelCrystalPt190S1M2;
	TH1D			*fModelCrystalPt190S10M2;
	TH1D			*fModelCrystalPt190S100M2;

	TH1D			*fModelCrystalPb210BM2;
	TH1D			*fModelCrystalPb210S01M2;
	TH1D			*fModelCrystalPb210S1M2;
	TH1D			*fModelCrystalPb210S10M2;
	TH1D			*fModelCrystalPb210S100M2;

	TH1D			*fModelCrystalRa226BM2;
	TH1D			*fModelCrystalRa226S01M2;
	TH1D			*fModelCrystalRa226S1M2;
	TH1D			*fModelCrystalRa226S10M2;
	TH1D			*fModelCrystalRa226S100M2;

	TH1D			*fModelFramePb210BM2;
	TH1D			*fModelFramePb210S01M2;
	TH1D			*fModelFramePb210S1M2;
	TH1D			*fModelFramePb210S10M2;
	TH1D			*fModelFramePb210S100M2;

	// TH1D			*fModelFrameRa226BM2;
	TH1D			*fModelFrameRa226S01M2;
	TH1D			*fModelFrameRa226S1M2;
	TH1D			*fModelFrameRa226S10M2;
	TH1D			*fModelFrameRa226S100M2;


	// Adaptive binned histograms

	TH1D *fAdapFrameThM1;
	TH1D *fAdapTShieldThM1;
	TH1D *fAdap50mKThM1;
	TH1D *fAdap600mKThM1;
	TH1D *fAdapIVCThM1;
	TH1D *fAdapOVCThM1;

	TH1D *fAdapFrameRaM1;
	TH1D *fAdapTShieldRaM1;
	TH1D *fAdap50mKRaM1;
	TH1D *fAdap600mKRaM1;
	TH1D *fAdapIVCRaM1;
	TH1D *fAdapOVCRaM1;
	
	TH1D *fAdapFrameKM1;
	TH1D *fAdapTShieldKM1;
	TH1D *fAdap50mKKM1;
	TH1D *fAdap600mKKM1;
	TH1D *fAdapIVCKM1;
	TH1D *fAdapOVCKM1;
	
	TH1D *fAdapFrameCoM1;
	TH1D *fAdapTShieldCoM1;
	TH1D *fAdap50mKCoM1;
	TH1D *fAdap600mKCoM1;
	TH1D *fAdapIVCCoM1;
	TH1D *fAdapOVCCoM1;
	
	TH1D *fAdapTShieldMnM1;
	TH1D *fAdapIVCMnM1;
	
	TH1D *fAdapCrystalBi2M1;
	TH1D *fAdapFrameBi2M1;

	TH1D *fAdapNDBDM1;
	TH1D *fAdap2NDBDM1;
	TH1D *fAdapBiM1;



	TH1D *fAdapFrameThM2;
	TH1D *fAdapTShieldThM2;
	TH1D *fAdap50mKThM2;
	TH1D *fAdap600mKThM2;
	TH1D *fAdapIVCThM2;
	TH1D *fAdapOVCThM2;

	TH1D *fAdapFrameRaM2;
	TH1D *fAdapTShieldRaM2;
	TH1D *fAdap50mKRaM2;
	TH1D *fAdap600mKRaM2;
	TH1D *fAdapIVCRaM2;
	TH1D *fAdapOVCRaM2;
	
	TH1D *fAdapFrameKM2;
	TH1D *fAdapTShieldKM2;
	TH1D *fAdap50mKKM2;
	TH1D *fAdap600mKKM2;
	TH1D *fAdapIVCKM2;
	TH1D *fAdapOVCKM2;
	
	TH1D *fAdapFrameCoM2;
	TH1D *fAdapTShieldCoM2;
	TH1D *fAdap50mKCoM2;
	TH1D *fAdap600mKCoM2;
	TH1D *fAdapIVCCoM2;
	TH1D *fAdapOVCCoM2;
	
	TH1D *fAdapTShieldMnM2;
	TH1D *fAdapIVCMnM2;
	
	TH1D *fAdapCrystalBi2M2;
	TH1D *fAdapFrameBi2M2;


	// M1 Alphas
	TH1D			*fAdapCrystalPt190BM1;
	TH1D			*fAdapCrystalPt190S01M1;
	TH1D			*fAdapCrystalPt190S1M1;
	TH1D			*fAdapCrystalPt190S10M1;
	TH1D			*fAdapCrystalPt190S100M1;

	TH1D			*fAdapCrystalPb210BM1;
	TH1D			*fAdapCrystalPb210S01M1;
	TH1D			*fAdapCrystalPb210S1M1;
	TH1D			*fAdapCrystalPb210S10M1;
	TH1D			*fAdapCrystalPb210S100M1;

	TH1D			*fAdapCrystalRa226BM1;
	TH1D			*fAdapCrystalRa226S01M1;
	TH1D			*fAdapCrystalRa226S1M1;
	TH1D			*fAdapCrystalRa226S10M1;
	TH1D			*fAdapCrystalRa226S100M1;

	TH1D			*fAdapFramePb210BM1;
	TH1D			*fAdapFramePb210S01M1;
	TH1D			*fAdapFramePb210S1M1;
	TH1D			*fAdapFramePb210S10M1;
	TH1D			*fAdapFramePb210S100M1;

	// TH1D			*fAdapFrameRa226BM1;
	TH1D			*fAdapFrameRa226S01M1;
	TH1D			*fAdapFrameRa226S1M1;
	TH1D			*fAdapFrameRa226S10M1;
	TH1D			*fAdapFrameRa226S100M1;

	// M2 Alphas
	TH1D			*fAdapCrystalPt190BM2;
	TH1D			*fAdapCrystalPt190S01M2;
	TH1D			*fAdapCrystalPt190S1M2;
	TH1D			*fAdapCrystalPt190S10M2;
	TH1D			*fAdapCrystalPt190S100M2;

	TH1D			*fAdapCrystalPb210BM2;
	TH1D			*fAdapCrystalPb210S01M2;
	TH1D			*fAdapCrystalPb210S1M2;
	TH1D			*fAdapCrystalPb210S10M2;
	TH1D			*fAdapCrystalPb210S100M2;

	TH1D			*fAdapCrystalRa226BM2;
	TH1D			*fAdapCrystalRa226S01M2;
	TH1D			*fAdapCrystalRa226S1M2;
	TH1D			*fAdapCrystalRa226S10M2;
	TH1D			*fAdapCrystalRa226S100M2;

	TH1D			*fAdapFramePb210BM2;
	TH1D			*fAdapFramePb210S01M2;
	TH1D			*fAdapFramePb210S1M2;
	TH1D			*fAdapFramePb210S10M2;
	TH1D			*fAdapFramePb210S100M2;

	// TH1D			*fAdapFrameRa226BM2;
	TH1D			*fAdapFrameRa226S01M2;
	TH1D			*fAdapFrameRa226S1M2;
	TH1D			*fAdapFrameRa226S10M2;
	TH1D			*fAdapFrameRa226S100M2;



	// For accidental coincidence test
	TFile *fFileCoin;
	TFile *fFileCorrection;
	TH1D *fCorrectionM2; // Correction spectra for M2 (for accidental coincidences)
	TH1D *fCorrectionM2Tot;
	TH1D *fTotCorrection;


	std::string		dDataDir;
	bool			bFixedRes;
	bool			bAdaptiveBinning;
	
	int 			dNumCalls;
	int 			dMult;

	// Parameters
	double				fParameters[26];
	double				fParError[26];
	double				fResolution[52];
	double				dSecToYears;
	double				fMCEff[26];

//  ClassDef(TMyFitter,1) // 
    };

#endif // __TBackgroundModel__
