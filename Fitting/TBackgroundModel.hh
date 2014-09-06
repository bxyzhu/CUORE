#ifndef __TBackgroundModel__
#define __TBackgroundModel__
#include "TObject.h"
#include "TH1D.h"
#include "TF1.h"
#include "TChain.h"
#include "TCut.h"
#include "TGraphErrors.h"
#include "TRandom3.h"

class TBackgroundModel : public TObject {

public:

	TBackgroundModel(double fFitMin, double fFitMax, int fMult, bool fFixedRes);
	virtual ~TBackgroundModel();

	TGraphErrors *CalculateResiduals(TH1D *h1, TH1D *h2);
  
	double GetChiSquare();

	double GetMCEff(TH1D *h1);

	bool DoTheFit();

	void DrawBkg();

	void DrawMC();

	void Initialize();

	void LoadData();

	TChain *LoadMC(std::string dDir, std::string dLocation, std::string dSource, std::string dSType, int dMult);

	void NormalizePDF(TH1D *h1, TChain *hChain, int minE, int maxE);

	void PrintParameters();

	void ReadMC();

	void SetParameters(int index, double value);

	TH1D *SmearMC(TH1D *hMC, TH1D *hSMC, double resolution);

	TH1D *SmearChannel(TChain *fChain, TH1D *fMC, TH1D *fSMC, int channel, double resolution);

	void TestSmear();

	void UpdateModel();


	int 	dBinSize;
	int 	dNBins;
	double	dMinEnergy;
	double 	dMaxEnergy;
	double	dFitMin;
	double	dFitMax;

private:

	// Data
	TChain			*qtree;
	TCut 			base_cut;
	TCut			ener_cut;

	double			dDataIntegral;

	TH1D			*fDataHistoTot;
	TH1D			*fDataHistoM1;
	TH1D			*fDataHistoM2;

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

	TChain			*outTreeNDBDM1;
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

	TChain			*outTreeNDBDM2;
	TChain			*outTreeBiM2;


	// Total PDFs M1
	TH1D			*fModelTotM1;
	TH1D			*fModelTotThM1;
	TH1D			*fModelTotRaM1;
	TH1D			*fModelTotKM1;
	TH1D			*fModelTotCoM1;

	TH1D			*fModelTotNDBDM1;
	TH1D			*fModelTotBiM1;
	TH1D			*fModelNDBDM1;
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

	// Total PDFs M2
	TH1D			*fModelTotM2;
	TH1D			*fModelTotThM2;
	TH1D			*fModelTotRaM2;
	TH1D			*fModelTotKM2;
	TH1D			*fModelTotCoM2;

	TH1D			*fModelTotNDBDM2;
	TH1D			*fModelTotBiM2;
	TH1D			*fModelNDBDM2;
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


	// Smeared PDFs M1
	TH1D			*fSmearDummyM1;

	TH1D			*fSmearNDBDM1;
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

	// Smeared PDFs M2
	TH1D			*fSmearDummyM2;

	TH1D			*fSmearNDBDM2;
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



	TGraph			*gResidualM1;
	TGraph 			*gResidualM2;

	TH1D			*hResidualDistM1;
	TH1D			*hResidualDistM2;

	// Smearing
	TRandom3		*fRandomGenerator;	
	TF1				*gaus;

	std::string		dDataDir;
	bool			bToyFit;
	bool			bFixedRes;
	
	int 			dNumCalls;
	int 			dMult;

	// Parameters
	double				fParameters[19];
	double				fParError[19];
	double				fResolution[52];
	double				dSecToYears;
	double				fMCEff[26];

//  ClassDef(TMyFitter,1) // 
    };

#endif // __TBackgroundModel__
