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

	// Model
	TChain			*outTreeFrameThS01;
	TChain			*outTreeFrameThS1;
	TChain			*outTreeFrameThS10;
	TChain			*outTreeFrameThS100;

	TChain			*outTreeFrameRaS01;
	TChain			*outTreeFrameRaS1;
	TChain			*outTreeFrameRaS10;
	TChain			*outTreeFrameRaS100;

	TChain			*outTreeTShieldThS01;
	TChain			*outTreeTShieldThS1;
	TChain			*outTreeTShieldThS10;
	TChain			*outTreeTShieldThS100;

	TChain			*outTreeFrameTh;
	TChain			*outTreeTShieldTh;	
	TChain			*outTree50mKTh;
	TChain			*outTree600mKTh;
	TChain			*outTreeIVCTh;
	TChain			*outTreeOVCTh;

	TChain			*outTreeFrameRa;
	TChain			*outTreeTShieldRa;
	TChain			*outTree50mKRa;
	TChain			*outTree600mKRa;
	TChain			*outTreeIVCRa;
	TChain			*outTreeOVCRa;

	TChain			*outTreeFrameK;
	TChain			*outTreeTShieldK;
	TChain			*outTree50mKK;
	TChain			*outTree600mKK;
	TChain			*outTreeIVCK;
	TChain			*outTreeOVCK;

	TChain			*outTreeFrameCo;
	TChain			*outTreeTShieldCo;
	TChain			*outTree50mKCo;
	TChain			*outTree600mKCo;
	TChain			*outTreeIVCCo;
	TChain			*outTreeOVCCo;

	TChain			*outTreeNDBD;
	TChain			*outTreeBi;

	// Total PDFs
	TH1D			*fModelTot;
	TH1D			*fModelTotTh;
	TH1D			*fModelTotRa;
	TH1D			*fModelTotK;
	TH1D			*fModelTotCo;

	TH1D			*fModelTotNDBD;
	TH1D			*fModelTotBi;
	TH1D			*fModelNDBD;
	TH1D			*fModelBi;

	TH1D			*fModelFrameThS01;
	TH1D			*fModelFrameThS1;
	TH1D			*fModelFrameThS10;
	TH1D			*fModelFrameThS100;

	TH1D			*fModelFrameRaS01;
	TH1D			*fModelFrameRaS1;
	TH1D			*fModelFrameRaS10;
	TH1D			*fModelFrameRaS100;

	TH1D			*fModelTShieldThS01;
	TH1D			*fModelTShieldThS1;
	TH1D			*fModelTShieldThS10;
	TH1D			*fModelTShieldThS100;


	TH1D			*fModelFrameTh;
	TH1D			*fModelTShieldTh;
	TH1D			*fModel50mKTh;
	TH1D			*fModel600mKTh;	
	TH1D			*fModelIVCTh;
	TH1D			*fModelOVCTh;

	TH1D			*fModelFrameRa;
	TH1D			*fModelTShieldRa;
	TH1D			*fModel50mKRa;
	TH1D			*fModel600mKRa;
	TH1D			*fModelIVCRa;
	TH1D			*fModelOVCRa;

	TH1D			*fModelFrameK;
	TH1D			*fModelTShieldK;
	TH1D			*fModel50mKK;
	TH1D			*fModel600mKK;
	TH1D			*fModelIVCK;
	TH1D			*fModelOVCK;


	TH1D			*fModelFrameCo;
	TH1D			*fModelTShieldCo;
	TH1D			*fModel50mKCo;
	TH1D			*fModel600mKCo;
	TH1D			*fModelIVCCo;
	TH1D			*fModelOVCCo;


	// Smeared PDFs
	TH1D			*fSmearDummy;

	TH1D			*fSmearNDBD;
	TH1D			*fSmearBi;

	TH1D			*fSmearFrameThS01;
	TH1D			*fSmearFrameThS1;
	TH1D			*fSmearFrameThS10;
	TH1D			*fSmearFrameThS100;

	TH1D			*fSmearFrameRaS01;
	TH1D			*fSmearFrameRaS1;
	TH1D			*fSmearFrameRaS10;
	TH1D			*fSmearFrameRaS100;

	TH1D			*fSmearTShieldThS01;
	TH1D			*fSmearTShieldThS1;
	TH1D			*fSmearTShieldThS10;
	TH1D			*fSmearTShieldThS100;


	TH1D			*fSmearFrameTh;
	TH1D			*fSmearTShieldTh;
	TH1D			*fSmear50mKTh;
	TH1D			*fSmear600mKTh;	
	TH1D			*fSmearIVCTh;
	TH1D			*fSmearOVCTh;

	TH1D			*fSmearFrameRa;
	TH1D			*fSmearTShieldRa;
	TH1D			*fSmear50mKRa;
	TH1D			*fSmear600mKRa;
	TH1D			*fSmearIVCRa;
	TH1D			*fSmearOVCRa;

	TH1D			*fSmearFrameK;
	TH1D			*fSmearTShieldK;
	TH1D			*fSmear50mKK;
	TH1D			*fSmear600mKK;
	TH1D			*fSmearIVCK;
	TH1D			*fSmearOVCK;


	TH1D			*fSmearFrameCo;
	TH1D			*fSmearTShieldCo;
	TH1D			*fSmear50mKCo;
	TH1D			*fSmear600mKCo;
	TH1D			*fSmearIVCCo;
	TH1D			*fSmearOVCCo;

	TGraph			*gResidual;
	TH1D			*hResidualDist;

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
