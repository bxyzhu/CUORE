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

	void Initialize();

	// Dumb to have all of these but w/e
	void LoadData();

	void LoadPDFs();

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


	// To be updated 11-05-2014
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
	TH1D			*fModelTotptM1;
	TH1D			*fModelTotpbM1;
	TH1D			*fModelTotcsM1;
	TH1D			*fModelTotco2M1;
	TH1D			*fModelTotteo2M1;


	TH1D			*fModelTotAdapM1;
	TH1D			*fModelTotAdapthM1;
	TH1D			*fModelTotAdapuM1;
	TH1D			*fModelTotAdapkM1;
	TH1D			*fModelTotAdapcoM1;
	TH1D			*fModelTotAdapmnM1;
	TH1D			*fModelTotAdapNDBDM1;
	TH1D			*fModelTotAdap2NDBDM1;
	TH1D			*fModelTotAdapbiM1;
	TH1D			*fModelTotAdapptM1;
	TH1D			*fModelTotAdappbM1;
	TH1D			*fModelTotAdapcsM1;
	TH1D			*fModelTotAdapco2M1;
	TH1D			*fModelTotAdapteo2M1;

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
	TH1D			*fModelTotptM2;
	TH1D			*fModelTotpbM2;
	TH1D			*fModelTotcsM2;
	TH1D			*fModelTotco2M2;
	TH1D			*fModelTotteo2M2;


	TH1D			*fModelTotAdapM2;
	TH1D			*fModelTotAdapthM2;
	TH1D			*fModelTotAdapuM2;
	TH1D			*fModelTotAdapkM2;
	TH1D			*fModelTotAdapcoM2;
	TH1D			*fModelTotAdapmnM2;
	TH1D			*fModelTotAdapNDBDM2;
	TH1D			*fModelTotAdap2NDBDM2;
	TH1D			*fModelTotAdapbiM2;
	TH1D			*fModelTotAdapptM2;
	TH1D			*fModelTotAdappbM2;
	TH1D			*fModelTotAdapcsM2;
	TH1D			*fModelTotAdapco2M2;
	TH1D			*fModelTotAdapteo2M2;


	// Residuals
	TGraph			*gResidualM1;
	TGraph 			*gResidualM2;

	TH1D			*hResidualDistM1;
	TH1D			*hResidualDistM2;

	TH1D			*hResidualGausM1;
	TH1D			*hResidualGausM2;



	TH1D			*fModelDummyM1;


//////////// Bulk Histograms
	// Crystal M1 and M2
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

	// Frame M1 and M2
	TH1D			*hCuFrameco58M1;
	TH1D			*hCuFrameco60M1;
	TH1D			*hCuFramecs137M1;
	TH1D			*hCuFramek40M1;
	TH1D			*hCuFramemn54M1;
	TH1D			*hCuFramepb210M1;
	TH1D			*hCuFrameth232M1;
	TH1D			*hCuFrameu238M1;

	TH1D			*hCuFrameco58M2;
	TH1D			*hCuFrameco60M2;
	TH1D			*hCuFramecs137M2;
	TH1D			*hCuFramek40M2;
	TH1D			*hCuFramemn54M2;
	TH1D			*hCuFramepb210M2;
	TH1D			*hCuFrameth232M2;
	TH1D			*hCuFrameu238M2;

	// CuBox (TShield) M1 and M2
	TH1D			*hCuBoxco58M1;
	TH1D			*hCuBoxco60M1;
	TH1D			*hCuBoxcs137M1;
	TH1D			*hCuBoxk40M1;
	TH1D			*hCuBoxmn54M1;
	TH1D			*hCuBoxpb210M1;
	TH1D			*hCuBoxth232M1;
	TH1D			*hCuBoxu238M1;	

	TH1D			*hCuBoxco58M2;
	TH1D			*hCuBoxco60M2;
	TH1D			*hCuBoxcs137M2;
	TH1D			*hCuBoxk40M2;
	TH1D			*hCuBoxmn54M2;
	TH1D			*hCuBoxpb210M2;
	TH1D			*hCuBoxth232M2;
	TH1D			*hCuBoxu238M2;	

	// 50mK M1 and M2
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

	// 600mK M1 and M2
	TH1D			*h600mKco60M1;
	TH1D			*h600mKk40M1;
	TH1D			*h600mKth232M1;
	TH1D			*h600mKu238M1;		

	TH1D			*h600mKco60M2;
	TH1D			*h600mKk40M2;
	TH1D			*h600mKth232M2;
	TH1D			*h600mKu238M2;	

	// Roman Lead M1 and M2
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


	// Main Bath M1 and M2
	TH1D			*hMBco60M1;
	TH1D			*hMBk40M1;
	TH1D			*hMBth232M1;
	TH1D			*hMBu238M1;		

	TH1D			*hMBco60M2;
	TH1D			*hMBk40M2;
	TH1D			*hMBth232M2;
	TH1D			*hMBu238M2;	

	// IVC M1 and M2
	TH1D			*hIVCco60M1;
	TH1D			*hIVCk40M1;
	TH1D			*hIVCth232M1;
	TH1D			*hIVCu238M1;		

	TH1D			*hIVCco60M2;
	TH1D			*hIVCk40M2;
	TH1D			*hIVCth232M2;
	TH1D			*hIVCu238M2;	

	// OVC M1 and M2
	TH1D			*hOVCco60M1;
	TH1D			*hOVCk40M1;
	TH1D			*hOVCth232M1;
	TH1D			*hOVCu238M1;		

	TH1D			*hOVCco60M2;
	TH1D			*hOVCk40M2;
	TH1D			*hOVCth232M2;
	TH1D			*hOVCu238M2;	


//////// Adaptive binned histograms
	// Crystal M1 and M2
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

	// Frame M1 and M2
	TH1D			*hAdapCuFrameco58M1;
	TH1D			*hAdapCuFrameco60M1;
	TH1D			*hAdapCuFramecs137M1;
	TH1D			*hAdapCuFramek40M1;
	TH1D			*hAdapCuFramemn54M1;
	TH1D			*hAdapCuFramepb210M1;
	TH1D			*hAdapCuFrameth232M1;
	TH1D			*hAdapCuFrameu238M1;

	TH1D			*hAdapCuFrameco58M2;
	TH1D			*hAdapCuFrameco60M2;
	TH1D			*hAdapCuFramecs137M2;
	TH1D			*hAdapCuFramek40M2;
	TH1D			*hAdapCuFramemn54M2;
	TH1D			*hAdapCuFramepb210M2;
	TH1D			*hAdapCuFrameth232M2;
	TH1D			*hAdapCuFrameu238M2;

	// CuBox (TShield) M1 and M2
	TH1D			*hAdapCuBoxco58M1;
	TH1D			*hAdapCuBoxco60M1;
	TH1D			*hAdapCuBoxcs137M1;
	TH1D			*hAdapCuBoxk40M1;
	TH1D			*hAdapCuBoxmn54M1;
	TH1D			*hAdapCuBoxpb210M1;
	TH1D			*hAdapCuBoxth232M1;
	TH1D			*hAdapCuBoxu238M1;	

	TH1D			*hAdapCuBoxco58M2;
	TH1D			*hAdapCuBoxco60M2;
	TH1D			*hAdapCuBoxcs137M2;
	TH1D			*hAdapCuBoxk40M2;
	TH1D			*hAdapCuBoxmn54M2;
	TH1D			*hAdapCuBoxpb210M2;
	TH1D			*hAdapCuBoxth232M2;
	TH1D			*hAdapCuBoxu238M2;	

	// 50mK M1 and M2
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

	// 600mK M1 and M2
	TH1D			*hAdap600mKco60M1;
	TH1D			*hAdap600mKk40M1;
	TH1D			*hAdap600mKth232M1;
	TH1D			*hAdap600mKu238M1;		

	TH1D			*hAdap600mKco60M2;
	TH1D			*hAdap600mKk40M2;
	TH1D			*hAdap600mKth232M2;
	TH1D			*hAdap600mKu238M2;	

	// Roman Lead M1 and M2
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


	// Main Bath M1 and M2
	TH1D			*hAdapMBco60M1;
	TH1D			*hAdapMBk40M1;
	TH1D			*hAdapMBth232M1;
	TH1D			*hAdapMBu238M1;		

	TH1D			*hAdapMBco60M2;
	TH1D			*hAdapMBk40M2;
	TH1D			*hAdapMBth232M2;
	TH1D			*hAdapMBu238M2;	

	// IVC M1 and M2
	TH1D			*hAdapIVCco60M1;
	TH1D			*hAdapIVCk40M1;
	TH1D			*hAdapIVCth232M1;
	TH1D			*hAdapIVCu238M1;		

	TH1D			*hAdapIVCco60M2;
	TH1D			*hAdapIVCk40M2;
	TH1D			*hAdapIVCth232M2;
	TH1D			*hAdapIVCu238M2;	

	// OVC M1 and M2
	TH1D			*hAdapOVCco60M1;
	TH1D			*hAdapOVCk40M1;
	TH1D			*hAdapOVCth232M1;
	TH1D			*hAdapOVCu238M1;		

	TH1D			*hAdapOVCco60M2;
	TH1D			*hAdapOVCk40M2;
	TH1D			*hAdapOVCth232M2;
	TH1D			*hAdapOVCu238M2;	



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
	double				fParameters[62];
	double				fParError[62];
	double				fResolution[52];
	double				dSecToYears;
	double				fMCEff[62];

//  ClassDef(TMyFitter,1) // 
    };

#endif // __TBackgroundModel__
