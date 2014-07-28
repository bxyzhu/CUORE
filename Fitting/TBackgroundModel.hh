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

	TBackgroundModel(double fFitMin, double fFitMax, int fMult);
	virtual ~TBackgroundModel();

	TGraphErrors *CalculateResiduals(TH1D *h1, TH1D *h2);
  
	double GetChiSquare();

	bool DoTheFit();

	void DrawBkg();
	void DrawMC();

	void Initialize();

	void LoadData();

	TChain *LoadMC(std::string dDir, std::string dLocation, std::string dSource, int dMult);

	void NormalizePDF(TH1D *h1, TChain *hChain, int minE, int maxE);

	void PrintParameters();

	void ReadMC();

	void SetParameters(int index, double value);

	TH1D *SmearMC(TH1D *hMC, TH1D *hSMC, double resolution);

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

	TChain			*outTreeToyTh;
	TChain			*outTreeToyRa;
	TChain			*outTreeToyCo;
	TChain			*outTreeToyK;
	TChain			*outTreeToyNDBD;


	// Model
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

	// Total PDFs
	TH1D			*fModelTot;
	TH1D			*fModelTotTh;
	TH1D			*fModelTotRa;
	TH1D			*fModelTotK;
	TH1D			*fModelTotCo;

	TH1D			*fModelTotNDBD;
	TH1D			*fModelNDBD;

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
	
	int 			dNumCalls;
	int 			dMult;

	// Parameters
	double				fParameters[10];
	double				dSecToYears;

//  ClassDef(TMyFitter,1) // 
    };

#endif // __TBackgroundModel__
