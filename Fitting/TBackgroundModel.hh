#ifndef __TBackgroundModel__
#define __TBackgroundModel__
#include "TObject.h"
#include "TH1D.h"
#include "TChain.h"
#include "TCut.h"
#include "TGraphErrors.h"
#include "TRandom3.h"

class TBackgroundModel : public TObject {

public:

	TBackgroundModel();
	virtual ~TBackgroundModel();

	TGraphErrors *CalculateResiduals(TH1D *h1, TH1D *h2);
  
	double GetChiSquare();

	bool DoTheFit();

	void DrawMC();

	void Initialize();

	void LoadData();

	TChain *LoadMC(std::string dLocation, std::string dSource, int dMult);

	void NormalizePDF(TH1D *h1);

	void PrintParameters();

	void ReadMC();

	void SetParameters(int index, double value);

	void UpdateModel();


	int 	dBinSize;
	int 	dNBins;
	double	dMinEnergy;
	double 	dMaxEnergy;

private:

	// Data
	TChain			*qtree;
	TCut 			base_cut;
	TCut			ener_cut;

	double			dDataIntegral;

	TH1D			*fDataHistoTot;
	TH1D			*fDataHistoM1;
	TH1D			*fDataHistoM2;

	// Model
	TChain			*outTree50mKTh;
	TChain			*outTree600mKTh;
	TChain			*outTreeIVCTh;
	TChain			*outTreeOVCTh;
	TChain			*outTreeFrameTh;

	TChain			*outTree50mKRa;
	TChain			*outTree600mKRa;
	TChain			*outTreeIVCRa;
	TChain			*outTreeOVCRa;
	TChain			*outTreeFrameRa;

	TChain			*outTreeFrameK;
	TChain			*outTreeFrameCo;

	// Total PDFs
	TH1D			*fModelTot;
	TH1D			*fModelTotTh;
	TH1D			*fModelTotRa;
	TH1D			*fModelTotK;
	TH1D			*fModelTotCo;

	TH1D			*fModel50mKTh;
	TH1D			*fModel600mKTh;
	TH1D			*fModelIVCTh;
	TH1D			*fModelOVCTh;
	TH1D			*fModelFrameTh;

	TH1D			*fModel50mKRa;
	TH1D			*fModel600mKRa;
	TH1D			*fModelIVCRa;
	TH1D			*fModelOVCRa;
	TH1D			*fModelFrameRa;

	TH1D			*fModelFrameK;
	TH1D			*fModelFrameCo;

	TGraph			*gResidual;
	// TH1D			*hResidualDist;

	// Smearing
	TRandom3		*fRandomGenerator;	
	
	// Parameters
	double			fParameters[12]; 



//  ClassDef(TMyFitter,1) // 
    };

#endif // __TBackgroundModel__
