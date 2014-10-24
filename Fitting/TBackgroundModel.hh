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
  
	double GetChiSquare();

	double GetChiSquareAdaptive();

	double GetMCEff(TH1D *h1);

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

	void SaveAlpha();

	void SaveGamma();

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


	// Toy data
	TH1D			*fToyData;
	TH1D			*fToyDataThTot;
	TH1D			*fToyDataRaTot;
	TH1D			*fToyDataCoTot;
	TH1D			*fToyDataKTot;

	TH1D			*fToyDataTh;
	TH1D			*fToyDataRa;
	TH1D			*fToyDataCo;
	TH1D			*fToyDataK;
	TH1D			*fToyDataNDBD;
	TH1D			*fToyDataBi;

	TChain			*outTreeToyTh;
	TChain			*outTreeToyRa;
	TChain			*outTreeToyCo;
	TChain			*outTreeToyK;
	TChain			*outTreeToyNDBD;
	TChain			*outTreeToyBi;

	// Model M1
	TChain			*outTreeFrameThS01M1;
	TChain			*outTreeFrameThS1M1;
	TChain			*outTreeFrameThS10M1;
	TChain			*outTreeFrameThS100M1;

	TChain			*outTreeFrameRaS01M1;
	TChain			*outTreeFrameRaS1M1;
	TChain			*outTreeFrameRaS10M1;
	TChain			*outTreeFrameRaS100M1;

	TChain			*outTreeTShieldThS01M1;
	TChain			*outTreeTShieldThS1M1;
	TChain			*outTreeTShieldThS10M1;
	TChain			*outTreeTShieldThS100M1;

	TChain			*outTreeFrameThM1;
	TChain			*outTreeTShieldThM1;	
	TChain			*outTree50mKThM1;
	TChain			*outTree600mKThM1;
	TChain			*outTreeIVCThM1;
	TChain			*outTreeOVCThM1;

	TChain			*outTreeFrameRaM1;
	TChain			*outTreeTShieldRaM1;
	TChain			*outTree50mKRaM1;
	TChain			*outTree600mKRaM1;
	TChain			*outTreeIVCRaM1;
	TChain			*outTreeOVCRaM1;

	TChain			*outTreeFrameKM1;
	TChain			*outTreeTShieldKM1;
	TChain			*outTree50mKKM1;
	TChain			*outTree600mKKM1;
	TChain			*outTreeIVCKM1;
	TChain			*outTreeOVCKM1;

	TChain			*outTreeFrameCoM1;
	TChain			*outTreeTShieldCoM1;
	TChain			*outTree50mKCoM1;
	TChain			*outTree600mKCoM1;
	TChain			*outTreeIVCCoM1;
	TChain			*outTreeOVCCoM1;

	TChain			*outTreeTShieldMnM1;
	TChain			*outTreeIVCMnM1;

	TChain			*outTreeNDBDM1;
	TChain			*outTree2NDBDM1;
	TChain			*outTreeBiM1;


	// Model M2
	TChain			*outTreeFrameThS01M2;
	TChain			*outTreeFrameThS1M2;
	TChain			*outTreeFrameThS10M2;
	TChain			*outTreeFrameThS100M2;

	TChain			*outTreeFrameRaS01M2;
	TChain			*outTreeFrameRaS1M2;
	TChain			*outTreeFrameRaS10M2;
	TChain			*outTreeFrameRaS100M2;

	TChain			*outTreeTShieldThS01M2;
	TChain			*outTreeTShieldThS1M2;
	TChain			*outTreeTShieldThS10M2;
	TChain			*outTreeTShieldThS100M2;

	TChain			*outTreeFrameThM2;
	TChain			*outTreeTShieldThM2;	
	TChain			*outTree50mKThM2;
	TChain			*outTree600mKThM2;
	TChain			*outTreeIVCThM2;
	TChain			*outTreeOVCThM2;

	TChain			*outTreeFrameRaM2;
	TChain			*outTreeTShieldRaM2;
	TChain			*outTree50mKRaM2;
	TChain			*outTree600mKRaM2;
	TChain			*outTreeIVCRaM2;
	TChain			*outTreeOVCRaM2;

	TChain			*outTreeFrameKM2;
	TChain			*outTreeTShieldKM2;
	TChain			*outTree50mKKM2;
	TChain			*outTree600mKKM2;
	TChain			*outTreeIVCKM2;
	TChain			*outTreeOVCKM2;

	TChain			*outTreeFrameCoM2;
	TChain			*outTreeTShieldCoM2;
	TChain			*outTree50mKCoM2;
	TChain			*outTree600mKCoM2;
	TChain			*outTreeIVCCoM2;
	TChain			*outTreeOVCCoM2;

	TChain			*outTreeTShieldMnM2;
	TChain			*outTreeIVCMnM2;

	TChain			*outTreeNDBDM2;
	TChain			*outTree2NDBDM2;	
	TChain			*outTreeBiM2;


	// Alphas
	TChain			*outTreeCrystalPt190BM1;
	TChain			*outTreeCrystalPt190S01M1;
	TChain			*outTreeCrystalPt190S1M1;
	TChain			*outTreeCrystalPt190S10M1;
	TChain			*outTreeCrystalPt190S100M1;

	TChain			*outTreeCrystalPb210BM1;
	TChain			*outTreeCrystalPb210S01M1;
	TChain			*outTreeCrystalPb210S1M1;
	TChain			*outTreeCrystalPb210S10M1;
	TChain			*outTreeCrystalPb210S100M1;

	TChain			*outTreeCrystalRa226BM1;
	TChain			*outTreeCrystalRa226S01M1;
	TChain			*outTreeCrystalRa226S1M1;
	TChain			*outTreeCrystalRa226S10M1;
	TChain			*outTreeCrystalRa226S100M1;

	TChain			*outTreeFramePb210BM1;
	TChain			*outTreeFramePb210S01M1;
	TChain			*outTreeFramePb210S1M1;
	TChain			*outTreeFramePb210S10M1;
	TChain			*outTreeFramePb210S100M1;

	TChain			*outTreeFrameRa226BM1;
	TChain			*outTreeFrameRa226S01M1;
	TChain			*outTreeFrameRa226S1M1;
	TChain			*outTreeFrameRa226S10M1;
	TChain			*outTreeFrameRa226S100M1;


	TChain			*outTreeCrystalPt190BM2;
	TChain			*outTreeCrystalPt190S01M2;
	TChain			*outTreeCrystalPt190S1M2;
	TChain			*outTreeCrystalPt190S10M2;
	TChain			*outTreeCrystalPt190S100M2;

	TChain			*outTreeCrystalPb210BM2;
	TChain			*outTreeCrystalPb210S01M2;
	TChain			*outTreeCrystalPb210S1M2;
	TChain			*outTreeCrystalPb210S10M2;
	TChain			*outTreeCrystalPb210S100M2;

	TChain			*outTreeCrystalRa226BM2;
	TChain			*outTreeCrystalRa226S01M2;
	TChain			*outTreeCrystalRa226S1M2;
	TChain			*outTreeCrystalRa226S10M2;
	TChain			*outTreeCrystalRa226S100M2;

	TChain			*outTreeFramePb210BM2;
	TChain			*outTreeFramePb210S01M2;
	TChain			*outTreeFramePb210S1M2;
	TChain			*outTreeFramePb210S10M2;
	TChain			*outTreeFramePb210S100M2;

	TChain			*outTreeFrameRa226BM2;
	TChain			*outTreeFrameRa226S01M2;
	TChain			*outTreeFrameRa226S1M2;
	TChain			*outTreeFrameRa226S10M2;
	TChain			*outTreeFrameRa226S100M2;


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

	TH1D			*fModelNDBDM1;
	TH1D			*fModel2NDBDM1;
	TH1D			*fModelBiM1;
	TH1D			*fModelBi2M1;

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

	TH1D			*fModelNDBDM2;
	TH1D			*fModel2NDBDM2;
	TH1D			*fModelBiM2;
	TH1D			*fModelBi2M2;

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



	// Residuals
	TGraph			*gResidualM1;
	TGraph 			*gResidualM2;

	TH1D			*hResidualDistM1;
	TH1D			*hResidualDistM2;

	TH1D			*hResidualGausM1;
	TH1D			*hResidualGausM2;

	// Smearing
	TRandom3		*fRandomGenerator;	
	TF1				*gaus;
	TF1 			*gaus2;



	// Smeared PDFs M1
	TH1D			*fSmearDummyM1;

	TH1D			*fSmearNDBDM1;
	TH1D			*fSmear2NDBDM1;
	TH1D			*fSmearBiM1;

	TH1D			*fSmearFrameThS01M1;
	TH1D			*fSmearFrameThS1M1;
	TH1D			*fSmearFrameThS10M1;
	TH1D			*fSmearFrameThS100M1;

	TH1D			*fSmearFrameRaS01M1;
	TH1D			*fSmearFrameRaS1M1;
	TH1D			*fSmearFrameRaS10M1;
	TH1D			*fSmearFrameRaS100M1;

	TH1D			*fSmearTShieldThS01M1;
	TH1D			*fSmearTShieldThS1M1;
	TH1D			*fSmearTShieldThS10M1;
	TH1D			*fSmearTShieldThS100M1;


	TH1D			*fSmearFrameThM1;
	TH1D			*fSmearTShieldThM1;
	TH1D			*fSmear50mKThM1;
	TH1D			*fSmear600mKThM1;	
	TH1D			*fSmearIVCThM1;
	TH1D			*fSmearOVCThM1;

	TH1D			*fSmearFrameRaM1;
	TH1D			*fSmearTShieldRaM1;
	TH1D			*fSmear50mKRaM1;
	TH1D			*fSmear600mKRaM1;
	TH1D			*fSmearIVCRaM1;
	TH1D			*fSmearOVCRaM1;

	TH1D			*fSmearFrameKM1;
	TH1D			*fSmearTShieldKM1;
	TH1D			*fSmear50mKKM1;
	TH1D			*fSmear600mKKM1;
	TH1D			*fSmearIVCKM1;
	TH1D			*fSmearOVCKM1;


	TH1D			*fSmearFrameCoM1;
	TH1D			*fSmearTShieldCoM1;
	TH1D			*fSmear50mKCoM1;
	TH1D			*fSmear600mKCoM1;
	TH1D			*fSmearIVCCoM1;
	TH1D			*fSmearOVCCoM1;

	TH1D			*fSmearTShieldMnM1;
	TH1D			*fSmearIVCMnM1;

	TH1D			*fSmearCrystalBi2M1;
	TH1D			*fSmearFrameBi2M1;

	TH1D			*fSmearCrystalPtM1;
	TH1D			*fSmearCrystalPbBM1;
	TH1D			*fSmearCrystalPbS01M1;
	TH1D			*fSmearCrystalPbS1M1;
	TH1D			*fSmearCrystalPbS10M1;
	TH1D			*fSmearCrystalPbS100M1;
	TH1D			*fSmearFramePbBM1;
	TH1D			*fSmearFramePbS01M1;
	TH1D			*fSmearFramePbS1M1;
	TH1D			*fSmearFramePbS10M1;
	TH1D			*fSmearFramePbS100M1;


	// Smeared PDFs M2
	TH1D			*fSmearDummyM2;

	TH1D			*fSmearNDBDM2;
	TH1D			*fSmear2NDBDM2;
	TH1D			*fSmearBiM2;

	TH1D			*fSmearFrameThS01M2;
	TH1D			*fSmearFrameThS1M2;
	TH1D			*fSmearFrameThS10M2;
	TH1D			*fSmearFrameThS100M2;

	TH1D			*fSmearFrameRaS01M2;
	TH1D			*fSmearFrameRaS1M2;
	TH1D			*fSmearFrameRaS10M2;
	TH1D			*fSmearFrameRaS100M2;

	TH1D			*fSmearTShieldThS01M2;
	TH1D			*fSmearTShieldThS1M2;
	TH1D			*fSmearTShieldThS10M2;
	TH1D			*fSmearTShieldThS100M2;


	TH1D			*fSmearFrameThM2;
	TH1D			*fSmearTShieldThM2;
	TH1D			*fSmear50mKThM2;
	TH1D			*fSmear600mKThM2;	
	TH1D			*fSmearIVCThM2;
	TH1D			*fSmearOVCThM2;

	TH1D			*fSmearFrameRaM2;
	TH1D			*fSmearTShieldRaM2;
	TH1D			*fSmear50mKRaM2;
	TH1D			*fSmear600mKRaM2;
	TH1D			*fSmearIVCRaM2;
	TH1D			*fSmearOVCRaM2;

	TH1D			*fSmearFrameKM2;
	TH1D			*fSmearTShieldKM2;
	TH1D			*fSmear50mKKM2;
	TH1D			*fSmear600mKKM2;
	TH1D			*fSmearIVCKM2;
	TH1D			*fSmearOVCKM2;


	TH1D			*fSmearFrameCoM2;
	TH1D			*fSmearTShieldCoM2;
	TH1D			*fSmear50mKCoM2;
	TH1D			*fSmear600mKCoM2;
	TH1D			*fSmearIVCCoM2;
	TH1D			*fSmearOVCCoM2;

	TH1D			*fSmearTShieldMnM2;
	TH1D			*fSmearIVCMnM2;

	TH1D			*fSmearCrystalBi2M2;
	TH1D			*fSmearFrameBi2M2;

	TH1D			*fSmearCrystalPtM2;
	TH1D			*fSmearCrystalPbBM2;
	TH1D			*fSmearCrystalPbS01M2;
	TH1D			*fSmearCrystalPbS1M2;
	TH1D			*fSmearCrystalPbS10M2;
	TH1D			*fSmearCrystalPbS100M2;
	TH1D			*fSmearFramePbBM2;
	TH1D			*fSmearFramePbS01M2;
	TH1D			*fSmearFramePbS1M2;
	TH1D			*fSmearFramePbS10M2;
	TH1D			*fSmearFramePbS100M2;



	// M1 Alphas
	TH1D			*fSmearCrystalPt190BM1;
	TH1D			*fSmearCrystalPt190S01M1;
	TH1D			*fSmearCrystalPt190S1M1;
	TH1D			*fSmearCrystalPt190S10M1;
	TH1D			*fSmearCrystalPt190S100M1;

	TH1D			*fSmearCrystalPb210BM1;
	TH1D			*fSmearCrystalPb210S01M1;
	TH1D			*fSmearCrystalPb210S1M1;
	TH1D			*fSmearCrystalPb210S10M1;
	TH1D			*fSmearCrystalPb210S100M1;

	TH1D			*fSmearCrystalRa226BM1;
	TH1D			*fSmearCrystalRa226S01M1;
	TH1D			*fSmearCrystalRa226S1M1;
	TH1D			*fSmearCrystalRa226S10M1;
	TH1D			*fSmearCrystalRa226S100M1;

	TH1D			*fSmearFramePb210BM1;
	TH1D			*fSmearFramePb210S01M1;
	TH1D			*fSmearFramePb210S1M1;
	TH1D			*fSmearFramePb210S10M1;
	TH1D			*fSmearFramePb210S100M1;

	// TH1D			*fSmearFrameRa226BM1;
	TH1D			*fSmearFrameRa226S01M1;
	TH1D			*fSmearFrameRa226S1M1;
	TH1D			*fSmearFrameRa226S10M1;
	TH1D			*fSmearFrameRa226S100M1;

	// M2 Alphas
	TH1D			*fSmearCrystalPt190BM2;
	TH1D			*fSmearCrystalPt190S01M2;
	TH1D			*fSmearCrystalPt190S1M2;
	TH1D			*fSmearCrystalPt190S10M2;
	TH1D			*fSmearCrystalPt190S100M2;

	TH1D			*fSmearCrystalPb210BM2;
	TH1D			*fSmearCrystalPb210S01M2;
	TH1D			*fSmearCrystalPb210S1M2;
	TH1D			*fSmearCrystalPb210S10M2;
	TH1D			*fSmearCrystalPb210S100M2;

	TH1D			*fSmearCrystalRa226BM2;
	TH1D			*fSmearCrystalRa226S01M2;
	TH1D			*fSmearCrystalRa226S1M2;
	TH1D			*fSmearCrystalRa226S10M2;
	TH1D			*fSmearCrystalRa226S100M2;

	TH1D			*fSmearFramePb210BM2;
	TH1D			*fSmearFramePb210S01M2;
	TH1D			*fSmearFramePb210S1M2;
	TH1D			*fSmearFramePb210S10M2;
	TH1D			*fSmearFramePb210S100M2;

	// TH1D			*fSmearFrameRa226BM2;
	TH1D			*fSmearFrameRa226S01M2;
	TH1D			*fSmearFrameRa226S1M2;
	TH1D			*fSmearFrameRa226S10M2;
	TH1D			*fSmearFrameRa226S100M2;


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

	TH1D *fModelTestM1;
	TH1D *fModelTestM2;
	TH1D *fModelTest1;	
	TH1D *fModelTest2;

	TH1D *fFrameThM1; 
	TH1D *fTShieldThM1;
	TH1D *f50mKThM1;
	TH1D *f600mKThM1;
	TH1D *fIVCThM1;
	TH1D *fOVCThM1;

	TH1D *fFrameRaM1; 
	TH1D *fTShieldRaM1;
	TH1D *f50mKRaM1;
	TH1D *f600mKRaM1;
	TH1D *fIVCRaM1;
	TH1D *fOVCRaM1;

	TH1D *fFrameKM1; 
	TH1D *fTShieldKM1;
	TH1D *f50mKKM1;
	TH1D *f600mKKM1;
	TH1D *fIVCKM1;
	TH1D *fOVCKM1;

	TH1D *fFrameCoM1; 
	TH1D *fTShieldCoM1;
	TH1D *f50mKCoM1;
	TH1D *f600mKCoM1;
	TH1D *fIVCCoM1;
	TH1D *fOVCCoM1;

	TH1D *fTShieldMnM1;
	TH1D *fIVCMnM1;

	TH1D *fCrystalBi2M1;
	TH1D *fFrameBi2M1;

	TH1D *fNDBDM1;
	TH1D *f2NDBDM1;
	TH1D *fBiM1;


	TH1D *fFrameThM1; 
	TH1D *fTShieldThM1;
	TH1D *f50mKThM1;
	TH1D *f600mKThM1;
	TH1D *fIVCThM1;
	TH1D *fOVCThM1;

	TH1D *fFrameRaM1; 
	TH1D *fTShieldRaM1;
	TH1D *f50mKRaM1;
	TH1D *f600mKRaM1;
	TH1D *fIVCRaM1;
	TH1D *fOVCRaM1;

	TH1D *fFrameKM1; 
	TH1D *fTShieldKM1;
	TH1D *f50mKKM1;
	TH1D *f600mKKM1;
	TH1D *fIVCKM1;
	TH1D *fOVCKM1;

	TH1D *fFrameCoM1; 
	TH1D *fTShieldCoM1;
	TH1D *f50mKCoM1;
	TH1D *f600mKCoM1;
	TH1D *fIVCCoM1;
	TH1D *fOVCCoM1;

	TH1D *fTShieldMnM1;
	TH1D *fIVCMnM1;

	TH1D *fCrystalBi2M1;
	TH1D *fFrameBi2M1;

	TH1D *fNDBDM1;
	TH1D *f2NDBDM1;
	TH1D *fBiM1;

	TFile *fFileCorrection;
	TH1D *fCorrectionM2; // Correction spectra for M2 (for accidental coincidences)
	TH1D *fCorrectionM2Tot;
	TH1D *fTotCorrection;


	std::string		dDataDir;
	bool			bToyFit;
	bool			bFixedRes;
	bool			bUnSmeared;
	
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
