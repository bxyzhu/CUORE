#include "TMinuit.h"
#include "TBackgroundModel.hh"
#include "TRandom3.h"
#include "TPaveText.h"
#include "TAxis.h"
#include <cmath>
#include <iostream>
#include <string>


using namespace std;

//ClassImp(TBackgroundModel)
  
//first set up a global function that calls your classes method that calculates the quantity to minimise
void myExternal_FCN(int &n, double *grad, double &fval, double x[], int code)
{
	// Required External Wrapper for function to be minimized by Minuit 
 
 	// This gets called for each value of the parameters minuit tries
	// here the x array contains the parameters you are trying to fit
  
	// here myClass should inherit from TObject
	TBackgroundModel* Obj = (TBackgroundModel*)gMinuit->GetObjectFit(); 

	// implement a method in your class for setting the parameters and thus update the parameters of your fitter class 
	Obj->SetParameters(0,	x[0]);   
	Obj->SetParameters(1,	x[1]);  
	Obj->SetParameters(2,	x[2]);
	Obj->SetParameters(3,	x[3]);
	Obj->SetParameters(4,	x[4]);
	Obj->SetParameters(5,	x[5]);   
	Obj->SetParameters(6,	x[6]);  
	Obj->SetParameters(7,	x[7]);
	Obj->SetParameters(8,	x[8]);
	Obj->SetParameters(9,	x[9]);
  Obj->SetParameters(10, x[10]);
  // Obj->SetParameters(11, x[11]);
  // Obj->SetParameters(12, x[12]);
  // Obj->SetParameters(13, x[13]);
  // Obj->SetParameters(14, x[14]);
  // Obj->SetParameters(15, x[15]);
  // Obj->SetParameters(16, x[16]);
  // Obj->SetParameters(17, x[17]);
  // Obj->SetParameters(18, x[18]);


	Obj->UpdateModel();

	//implement a method in your class that calculates the quantity you want to minimise, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimise fval
	fval = Obj->GetChiSquare();
}


TBackgroundModel::TBackgroundModel(double fFitMin, double fFitMax, int fMult, bool fFixedRes)
{

  dNumCalls = 0;
  dSecToYears = 1./(60*60*24*365);

  dDataDir =  "/Users/brian/macros/Simulations/Bkg/Unsmeared/";
  dDataIntegral = 0;
  bToyFit = false;
  bFixedRes = fFixedRes;

  // Bin size (keV)
  dBinSize = 5;
  // Histogram range
  dMinEnergy = 0.;
  dMaxEnergy = 3500.;

  if(fFitMin >= fFitMax)
  {
    cout << "Trouble's a brewing: Fit Min >= Fit Max!" << endl;
  }

  // Fitting range
  dFitMin = fFitMin;
  dFitMax = fFitMax;

  dMult = fMult;


  dNBins = (dMaxEnergy - dMinEnergy)/ dBinSize;
  // Data
  qtree = new TChain("qtree");

 base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
 base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
 base_cut = base_cut && "abs(BaselineSlope)<0.1";
 base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";
//  base_cut = base_cut && Form("Energy > %f && Energy < %f", dMinEnergy, dMaxEnergy);

//  ener_cut = ener_cut && Form("Ener1 > %f && Ener1 < %f", dMinEnergy, dMaxEnergy);

  fDataHistoTot  = new TH1D("fDataHistoTot",  "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM1   = new TH1D("fDataHistoM1",   "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM2   = new TH1D("fDataHistoM2",   "", dNBins, dMinEnergy, dMaxEnergy);

  // Toy Data
  fToyData       = new TH1D("fToyData",        "", dNBins, dMinEnergy, dMaxEnergy);
  fToyDataThTot  = new TH1D("fToyDataThTot",   "", dNBins, dMinEnergy, dMaxEnergy);
  fToyDataRaTot  = new TH1D("fToyDataRaTot",   "", dNBins, dMinEnergy, dMaxEnergy);
  fToyDataCoTot  = new TH1D("fToyDataCoTot",   "", dNBins, dMinEnergy, dMaxEnergy);
  fToyDataKTot   = new TH1D("fToyDataKTot",    "", dNBins, dMinEnergy, dMaxEnergy);

  fToyDataTh     = new TH1D("fToyDataTh",   "", dNBins, dMinEnergy, dMaxEnergy);
  fToyDataRa     = new TH1D("fToyDataRa",   "", dNBins, dMinEnergy, dMaxEnergy);
  fToyDataCo     = new TH1D("fToyDataCo",   "", dNBins, dMinEnergy, dMaxEnergy);
  fToyDataK      = new TH1D("fToyDataK",    "", dNBins, dMinEnergy, dMaxEnergy);
  fToyDataNDBD   = new TH1D("fToyDataNDBD", "", dNBins, dMinEnergy, dMaxEnergy);
  fToyDataBi     = new TH1D("fToyDataBi",   "", dNBins, dMinEnergy, dMaxEnergy);


  // Model histograms
  fModelFrameThS01   = new TH1D("fModelFrameThS01",  "Frame Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameThS1    = new TH1D("fModelFrameThS1",  "Frame Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameThS10   = new TH1D("fModelFrameThS10",  "Frame Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameThS100  = new TH1D("fModelFrameThS100",  "Frame Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameRaS01   = new TH1D("fModelFrameRaS01",  "Frame Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRaS1    = new TH1D("fModelFrameRaS1",  "Frame Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRaS10   = new TH1D("fModelFrameRaS10",  "Frame Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRaS100  = new TH1D("fModelFrameRaS100",  "Frame Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameTh    = new TH1D("fModelFrameTh",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldTh  = new TH1D("fModelTShieldTh","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fModel50mKTh     = new TH1D("fModel50mKTh",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKTh    = new TH1D("fModel600mKTh",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCTh      = new TH1D("fModelIVCTh",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCTh      = new TH1D("fModelOVCTh",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameRa    = new TH1D("fModelFrameRa",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldRa  = new TH1D("fModelTShieldRa","TShield",  dNBins, dMinEnergy, dMaxEnergy);  
  fModel50mKRa     = new TH1D("fModel50mKRa",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKRa    = new TH1D("fModel600mKRa",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCRa      = new TH1D("fModelIVCRa",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCRa      = new TH1D("fModelOVCRa",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameK     = new TH1D("fModelFrameK",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldK   = new TH1D("fModelTShieldK", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fModel50mKK      = new TH1D("fModel50mKK",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKK     = new TH1D("fModel600mKK",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCK       = new TH1D("fModelIVCK",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCK       = new TH1D("fModelOVCK",   "OVC",        dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameCo    = new TH1D("fModelFrameCo",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldCo  = new TH1D("fModelTShieldCo","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fModel50mKCo     = new TH1D("fModel50mKCo",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKCo    = new TH1D("fModel600mKCo",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCCo      = new TH1D("fModelIVCCo",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCCo      = new TH1D("fModelOVCCo",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelNDBD       = new TH1D("fModelNDBD",     "NDBD",     dNBins, dMinEnergy, dMaxEnergy);
  fModelBi         = new TH1D("fModelBi",       "Bi207",    dNBins, dMinEnergy, dMaxEnergy);


  // Total model histograms
  fModelTot      = new TH1D("fModelTot",      "Frame",        dNBins, dMinEnergy, dMaxEnergy);  
  fModelTotTh    = new TH1D("fModelTotTh",    "Total Th232",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotRa    = new TH1D("fModelTotRa",    "Total Ra226",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotK     = new TH1D("fModelTotK",     "Total K40",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTotCo    = new TH1D("fModelTotCo",    "Total Co60",   dNBins, dMinEnergy, dMaxEnergy);

  fModelTotNDBD  = new TH1D("fModelTotNDBD",  "Total NDBD",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotBi    = new TH1D("fModelTotBi",   "Total Bi207",   dNBins, dMinEnergy, dMaxEnergy);


  // Smearing
  gaus = new TF1("gaus","gaus(0)", dMinEnergy, dMaxEnergy);

  // Smeared Histograms
  fSmearDummy      = new TH1D("fSmearDummy",  "Dummy smeared",  dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameThS01   = new TH1D("fSmearFrameThS01",  "Frame Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameThS1    = new TH1D("fSmearFrameThS1",  "Frame Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameThS10   = new TH1D("fSmearFrameThS10",  "Frame Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameThS100  = new TH1D("fSmearFrameThS100",  "Frame Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameRaS01   = new TH1D("fSmearFrameRaS01",  "Frame Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRaS1    = new TH1D("fSmearFrameRaS1",  "Frame Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRaS10   = new TH1D("fSmearFrameRaS10",  "Frame Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRaS100  = new TH1D("fSmearFrameRaS100",  "Frame Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameTh    = new TH1D("fSmearFrameTh",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldTh  = new TH1D("fSmearTShieldTh","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fSmear50mKTh     = new TH1D("fSmear50mKTh",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fSmear600mKTh    = new TH1D("fSmear600mKTh",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCTh      = new TH1D("fSmearIVCTh",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearOVCTh      = new TH1D("fSmearOVCTh",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameRa    = new TH1D("fSmearFrameRa",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldRa  = new TH1D("fSmearTShieldRa","TShield",  dNBins, dMinEnergy, dMaxEnergy);  
  fSmear50mKRa     = new TH1D("fSmear50mKRa",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fSmear600mKRa    = new TH1D("fSmear600mKRa",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCRa      = new TH1D("fSmearIVCRa",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearOVCRa      = new TH1D("fSmearOVCRa",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameK     = new TH1D("fSmearFrameK",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldK   = new TH1D("fSmearTShieldK", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fSmear50mKK      = new TH1D("fSmear50mKK",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fSmear600mKK     = new TH1D("fSmear600mKK",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCK       = new TH1D("fSmearIVCK",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearOVCK       = new TH1D("fSmearOVCK",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameCo    = new TH1D("fSmearFrameCo",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldCo  = new TH1D("fSmearTShieldCo","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fSmear50mKCo     = new TH1D("fSmear50mKCo",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fSmear600mKCo    = new TH1D("fSmear600mKCo",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCCo      = new TH1D("fSmearIVCCo",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearOVCCo      = new TH1D("fSmearOVCCo",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearNDBD       = new TH1D("fSmearNDBD",   "NDBD",       dNBins, dMinEnergy, dMaxEnergy);
  fSmearBi         = new TH1D("fSmearBi",     "Bi",         dNBins, dMinEnergy, dMaxEnergy);


  // Set Initial Parameters/Errors to 0
  fParameters[0]  = 0.;
  fParameters[1]  = 0.;
  fParameters[2]  = 0.;
  fParameters[3]  = 0.;
  fParameters[4]  = 0.;
  fParameters[5]  = 0.;
  fParameters[6]  = 0.;
  fParameters[7]  = 0.;
  fParameters[8]  = 0.;
  fParameters[9]  = 0.;
  fParameters[10] = 0.;
  fParameters[11] = 0.;
  fParameters[12] = 0.;
  fParameters[13] = 0.;
  fParameters[14] = 0.;
  fParameters[15] = 0.;
  fParameters[16] = 0.;
  fParameters[17] = 0.;
  fParameters[18] = 0.;


  fParError[0]  = 0.;
  fParError[1]  = 0.;
  fParError[2]  = 0.;
  fParError[3]  = 0.;
  fParError[4]  = 0.;
  fParError[5]  = 0.;
  fParError[6]  = 0.;
  fParError[7]  = 0.;
  fParError[8]  = 0.;
  fParError[9]  = 0.;
  fParError[10] = 0.;
  fParError[11] = 0.;
  fParError[12] = 0.;
  fParError[13] = 0.;
  fParError[14] = 0.;
  fParError[15] = 0.;
  fParError[16] = 0.;
  fParError[17] = 0.;
  fParError[18] = 0.;



  // Resolutions of individual channels
  fResolution[0] = 5.0/2.355;
  fResolution[1] = 6.0/2.355;
  fResolution[2] = 9.8/2.355;
  fResolution[3] = 5.3/2.355;

  fResolution[4] = 5.6/2.355;
  fResolution[5] = 10.0/2.355;
  fResolution[6] = 4.8/2.355;
  fResolution[7] = 4.5/2.355;

  fResolution[8] = 4.3/2.355;
  fResolution[9] = 7.5/2.355;
  fResolution[10] = 6.8/2.355;
  fResolution[11] = 5.2/2.355;

  fResolution[12] = 5.4/2.355;
  fResolution[13] = 4.3/2.355;
  fResolution[14] = 6.3/2.355;
  fResolution[15] = 5.9/2.355;

  fResolution[16] = 10.3/2.355;
  fResolution[17] = 5.0/2.355;
  fResolution[18] = 7.4/2.355;
  fResolution[19] = 5.2/2.355;

  fResolution[20] = 4.5/2.355;
  fResolution[21] = 4.5/2.355;
  fResolution[22] = 5.5/2.355;
  fResolution[23] = 4.7/2.355;

  fResolution[24] = 5.9/2.355;
  fResolution[25] = 8.0/2.355;
  fResolution[26] = 8.4/2.355;
  fResolution[27] = 6.4/2.355;

  fResolution[28] = 6.4/2.355;
  fResolution[29] = 4.7/2.355;
  fResolution[30] = 7.2/2.355;
  fResolution[31] = 6.1/2.355;

  fResolution[32] = 4.5/2.355;
  fResolution[33] = 11.9/2.355;
  fResolution[34] = 13.0/2.355;
  fResolution[35] = 6.7/2.355;

  fResolution[36] = 5.0/2.355;
  fResolution[37] = 4.9/2.355;
  fResolution[38] = 5.3/2.355;
  fResolution[39] = 5.0/2.355;

  fResolution[40] = 5.1/2.355;
  fResolution[41] = 6.5/2.355;
  fResolution[42] = 9.2/2.355;
  fResolution[43] = 12.0/2.355;

  fResolution[44] = 5.4/2.355;
  fResolution[45] = 5.5/2.355;
  fResolution[46] = 6.3/2.355;
  fResolution[47] = 6.6/2.355;

  fResolution[48] = 5.2/2.355;
  fResolution[49] = 6.4/2.355;
  fResolution[50] = 6.8/2.355;
  fResolution[51] = 7.8/2.355;

}
  

TBackgroundModel::~TBackgroundModel()
{
	delete	fDataHistoTot;
	delete	fDataHistoM1;
	delete	fDataHistoM2;
	delete	fToyData;

	delete 	fModelFrameTh;
	delete	fModelTShieldTh;
	delete	fModel50mKTh;
	delete	fModel600mKTh;
	delete 	fModelIVCTh;
	delete	fModelOVCTh;

	delete	fModelFrameRa
	delete 	fModelTShieldRa;
	delete	fModel50mKRa;
	delete	fModel600mKRa;
	delete 	fModelIVCRa;
	delete	fModelOVCRa;

	delete	fModelFrameK;
	delete 	fModelTShieldK;
	delete	fModel50mKK;
	delete	fModel600mKK;
	delete 	fModelIVCK;
	delete	fModelOVCK;

	delete	fModelFrameCo;
	delete 	fModelTShieldCo;
	delete	fModel50mKCo;
	delete	fModel600mKCo;
	delete 	fModelIVCCo;
	delete	fModelOVCCo;

	delete 	fModelTotTh;
	delete	fModelTotRa;
	delete	fModelTotK;
	delete	fModelTotCo;
  delete  fModelTotNDBD;
}



TH1D *TBackgroundModel::CalculateResiduals(TH1D *h1, TH1D *h2)
{

	// Clone histograms for rebinning
	TH1D 	*hCloneBkg 		= (TH1D*)h1->Clone("hCloneBkg");
	TH1D 	*hCloneMC		= (TH1D*)h2->Clone("hCloneMC");
	TH1D	*h1				= new TH1D("h1", "Fit Residuals", dNBins, dMinEnergy, dMaxEnergy);



	// Variables used in Residual calculations
	double dResidualX, dResidualY, dResidualXErr = 0, dResidualYErr;

	// Residual plot and distribution
	for (int j = dFitMin/dBinSize+1; j <= dFitMax/dBinSize; j++)
	{
		dResidualX 		= hCloneBkg->GetBinCenter(j);
		dResidualY 		= (hCloneBkg->GetBinContent(j) - hCloneMC->GetBinContent(j)) /
							TMath::Sqrt(hCloneBkg->GetBinContent(j)); // Sqrt(MC + data) = sigma for poisson distribution

		// g1->SetPoint(j, dResidualX, dResidualY);
		h1->SetBinContent(j, dResidualY);
		h1->SetBinError(j, 1);
	}

	return h1;
}



bool TBackgroundModel::DoTheFit()
{
  if(!bToyFit)
  {
    Initialize();
  }
	gStyle->SetOptStat(0);
   // This method actually sets up minuit and does the fit

 
   TMinuit minuit(11); //initialize minuit, n is the number of parameters

   // Reduce Minuit Output
   minuit.SetPrintLevel(1);
//   minuit.Command("SET MINImize 1000 0.001");
   minuit.Command("SET STRategy 2");
  //minuit.Command("SET IMProve 1000 ");
   minuit.SetMaxIterations(1000);
   minuit.SetObjectFit(this); //see the external FCN  above
   
   //define the parameters and set the ranges and initial guesses see ROOTs TMinuit documentation
   // Range is from 0 to integral of the data
   // Around 60000 events in background spectrum



   ////////////////////////////////////////////////
   // Using less parameters
   ////////////////////////////////////////////////
   ///// Maximum values are typically 1.5x to 2x measured rate of background just to be sure
   minuit.DefineParameter(0, "Close Th",  3000., 100.0, 0., 60000);
   // minuit.DefineParameter(0, "Close Th",  20000, 100.0, 0., 60000);   
   // minuit.DefineParameter(1, "Far Th",	 	2830., 50.0, 0., 4000);
   minuit.DefineParameter(1, "Far Th",   35000., 100.0, 0., 100000);
   minuit.DefineParameter(2, "Close Ra",  100., 100.0, 0., 50000);   
   minuit.DefineParameter(3, "Far Ra",    50000., 100.0, 0., 80000);
   // minuit.DefineParameter(4, "Close K", 	0., 100.0, 0., 500000);
   minuit.DefineParameter(4, "Close K",   100., 100.0, 0., 500000);
   // minuit.DefineParameter(5, "Far K",     0., 100.0, 0., 500000);
   minuit.DefineParameter(5, "Far K", 		25000., 100.0, 0., 500000);
   minuit.DefineParameter(6, "Close Co", 	3000., 100.0, 0., 50000); 
   minuit.DefineParameter(7, "Far Co",	 	100., 100.0, 0., 50000);  
   minuit.DefineParameter(8, "Resolution",	6., 1, 3, 10);  
   minuit.DefineParameter(9, "NDBD",      91.7., 100.0, 0., 1000);     
   minuit.DefineParameter(10, "Lead Bi",	 	5000., 100.0, 0., 200000);  
   // minuit.DefineParameter(10, "Lead Bi",    0., 100.0, 0., 200000);  

   // minuit.DefineParameter(11, "Surface Th",  100, 100.0, 0., 60000);   




/*
   // Close Th and Close Ra now split into its various sections, far Th and Ra still the same
   ////////////////////////////////////////////////
   // Using more parameters
   ////////////////////////////////////////////////
   ///// Maximum values are typically 1.5x to 2x measured rate of background just to be sure
   // minuit.DefineParameter(0, "Close Th",  0., 100.0, 0., 60000);
   minuit.DefineParameter(0, "Frame Th",  3380, 100.0, 0., 100000);   
   // minuit.DefineParameter(1, "Far Th",   2830., 50.0, 0., 4000);
   minuit.DefineParameter(1, "IVC Th",   40000., 100.0, 0., 100000);
   minuit.DefineParameter(2, "Frame Ra",  100., 100.0, 0., 50000);   
   minuit.DefineParameter(3, "IVC Ra",    55000., 100.0, 0., 80000);
   minuit.DefineParameter(4, "Close K",   100., 100.0, 0., 500000);
   minuit.DefineParameter(5, "Far K",     30000., 100.0, 0., 500000);
   minuit.DefineParameter(6, "Close Co",  3000., 100.0, 0., 50000); 
   minuit.DefineParameter(7, "Far Co",    100., 100.0, 0., 50000);  
   minuit.DefineParameter(8, "Resolution",  5., 1, 3, 10);  
   minuit.DefineParameter(9, "NDBD",      91.7., 100.0, 0., 1000);     
   minuit.DefineParameter(10, "Lead Bi",    5000., 100.0, 0., 200000);  
   minuit.DefineParameter(11, "TShield Th",  3380, 100.0, 0., 100000);   
   minuit.DefineParameter(12, "50mK Th",  3380, 100.0, 0., 100000);   
   minuit.DefineParameter(13, "600mK Th",  3380, 100.0, 0., 100000);   
   minuit.DefineParameter(14, "TShield Ra",  100, 100.0, 0., 60000);   
   minuit.DefineParameter(15, "50mK Ra",  100, 100.0, 0., 60000);   
   minuit.DefineParameter(16, "600mK Ra",  100, 100.0, 0., 60000);   
   minuit.DefineParameter(17, "OVC Th",  40000, 100.0, 0., 150000);   
   minuit.DefineParameter(18, "OVC Ra",  55000, 100.0, 0., 150000);   
*/




   // Fix parameters for testing
   // minuit.FixParameter(0); // Close Th
   // minuit.FixParameter(1); // Far Th
   // minuit.FixParameter(2); // Close Ra
   // minuit.FixParameter(3); // Far Ra
   // minuit.FixParameter(4); // Close K
   // minuit.FixParameter(5); // Far K
   // minuit.FixParameter(6); // Close Co
   // minuit.FixParameter(7); // Far Co
   if(bFixedRes)
   {
    minuit.FixParameter(8); // Resolution
   }
   // minuit.FixParameter(9); // NDBD
   // minuit.FixParameter(10); // Bi207

  // Number of Parameters! (for Chi-squared/NDF calculation)
  int dNumParameters = 10;




   //Tell minuit what external function to use 
   minuit.SetFCN(myExternal_FCN);
   
   int status = minuit.Migrad(); // this actually does the minimisation


    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->SetLogy();


    if(bToyFit)
    {
		////// Draw Toy Data
   		fToyData->SetLineColor(1);
   		fToyData->SetLineWidth(2);
   		fToyData->GetXaxis()->SetTitle("Energy (keV)");
   		fToyData->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));
		  fToyData->Draw();
	 }
	 else
	 {

    if(dMult == 1)
    {
   		///// Draw Data
  	 	fDataHistoM1->SetLineColor(1);
  	 	fDataHistoM1->SetLineWidth(2);
  	 	fDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
   		fDataHistoM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)/yr", dBinSize));
      fDataHistoM1->SetMaximum(90000);
      fDataHistoM1->GetXaxis()->SetRange(1, 2650/dBinSize+5);
		  fDataHistoM1->Draw();
    }
    else if(dMult == 2)
    {
      ///// Draw Data
      fDataHistoM2->SetLineColor(1);
      fDataHistoM2->SetLineWidth(2);
      fDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
      fDataHistoM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)/yr", dBinSize));
      fDataHistoM2->GetYaxis()->SetRange(0, 50000);
      fDataHistoM2->GetXaxis()->SetRange(2000/dBinSize-5, 2650/dBinSize+5);
      fDataHistoM2->Draw();

    }
	 }


	minuit.GetParameter(0,	fParameters[0],		fParError[0]);
	minuit.GetParameter(1,	fParameters[1],		fParError[1]);
	minuit.GetParameter(2,	fParameters[2],		fParError[2]);
	minuit.GetParameter(3,	fParameters[3],		fParError[3]);
	minuit.GetParameter(4,	fParameters[4],		fParError[4]);
	minuit.GetParameter(5,	fParameters[5],		fParError[5]);
	minuit.GetParameter(6,	fParameters[6],		fParError[6]);
	minuit.GetParameter(7,	fParameters[7],		fParError[7]);
	minuit.GetParameter(8,	fParameters[8],		fParError[8]);
	minuit.GetParameter(9,	fParameters[9],		fParError[9]);
  minuit.GetParameter(10,  fParameters[10],   fParError[10]);
  // minuit.GetParameter(11,  fParameters[11],   fParError[11]);
  // minuit.GetParameter(12,  fParameters[12],   fParError[12]);
  // minuit.GetParameter(13,  fParameters[13],   fParError[13]);
  // minuit.GetParameter(14,  fParameters[14],   fParError[14]);
  // minuit.GetParameter(15,  fParameters[15],   fParError[15]);
  // minuit.GetParameter(16,  fParameters[16],   fParError[16]);
  // minuit.GetParameter(17,  fParameters[17],   fParError[17]);
  // minuit.GetParameter(18,  fParameters[18],   fParError[18]);


	UpdateModel();
	
	cout << "At the end; ChiSq/NDF = " << GetChiSquare()/((dFitMax-dFitMin)/dBinSize - dNumParameters) << endl;
  cout << "Total number of calls = " << dNumCalls << endl;

  ///////////////////////////////////////////
  //// Few Parameters
  ///////////////////////////////////////////
  /// Add Histograms after chi-squared minimization

  // Surface....
  // fModelTotTh->Add(fSmearFrameThS1,   fParameters[11]);

  fModelTotTh->Add(fSmearFrameTh,   fParameters[0]);
  fModelTotTh->Add(fSmearTShieldTh, fParameters[0]);
  fModelTotTh->Add(fSmear50mKTh,    fParameters[0]);
  fModelTotTh->Add(fSmear600mKTh,   fParameters[0]);
  fModelTotTh->Add(fSmearIVCTh,     fParameters[1]);
  fModelTotTh->Add(fSmearOVCTh,     fParameters[1]);

  fModelTotRa->Add(fSmearFrameRa,   fParameters[2]);
  fModelTotRa->Add(fSmearTShieldRa, fParameters[2]);
  fModelTotRa->Add(fSmear50mKRa,    fParameters[2]);
  fModelTotRa->Add(fSmear600mKRa,   fParameters[2]);
  fModelTotRa->Add(fSmearIVCRa,     fParameters[3]);
  fModelTotRa->Add(fSmearOVCRa,     fParameters[3]);

  fModelTotK->Add(fSmearFrameK,     fParameters[4]);
  fModelTotK->Add(fSmearTShieldK,   fParameters[4]);
  fModelTotK->Add(fSmear50mKK,      fParameters[4]);
  fModelTotK->Add(fSmear600mKK,     fParameters[4]);
  fModelTotK->Add(fSmearIVCK,       fParameters[5]);
  fModelTotK->Add(fSmearOVCK,       fParameters[5]);

  fModelTotCo->Add(fSmearFrameCo,   fParameters[6]);
  fModelTotCo->Add(fSmearTShieldCo, fParameters[6]);
  fModelTotCo->Add(fSmear50mKCo,    fParameters[6]);
  fModelTotCo->Add(fSmear600mKCo,   fParameters[6]);
  fModelTotCo->Add(fSmearIVCCo,     fParameters[7]);
  fModelTotCo->Add(fSmearOVCCo,     fParameters[7]);

  fModelTotNDBD->Add(fSmearNDBD,    fParameters[9]);
  fModelTotBi->Add(fSmearBi,      fParameters[10]);



/*
  ///////////////////////////////////////////
  //// Many Parameters
  ///////////////////////////////////////////
  /// Add Histograms after chi-squared minimization

  // Surface....
  // fModelTotTh->Add(fSmearFrameThS1,   fParameters[11]);

  fModelTotTh->Add(fSmearFrameTh,   fParameters[0]);
  fModelTotTh->Add(fSmearTShieldTh, fParameters[11]);
  fModelTotTh->Add(fSmear50mKTh,    fParameters[12]);
  fModelTotTh->Add(fSmear600mKTh,   fParameters[13]);
  fModelTotTh->Add(fSmearIVCTh,     fParameters[1]);
  fModelTotTh->Add(fSmearOVCTh,     fParameters[17]);

  fModelTotRa->Add(fSmearFrameRa,   fParameters[2]);
  fModelTotRa->Add(fSmearTShieldRa, fParameters[14]);
  fModelTotRa->Add(fSmear50mKRa,    fParameters[15]);
  fModelTotRa->Add(fSmear600mKRa,   fParameters[16]);
  fModelTotRa->Add(fSmearIVCRa,     fParameters[3]);
  fModelTotRa->Add(fSmearOVCRa,     fParameters[18]);

  fModelTotK->Add(fSmearFrameK,     fParameters[4]);
  fModelTotK->Add(fSmearTShieldK,   fParameters[4]);
  fModelTotK->Add(fSmear50mKK,      fParameters[4]);
  fModelTotK->Add(fSmear600mKK,     fParameters[4]);
  fModelTotK->Add(fSmearIVCK,       fParameters[5]);
  fModelTotK->Add(fSmearOVCK,       fParameters[5]);

  fModelTotCo->Add(fSmearFrameCo,   fParameters[6]);
  fModelTotCo->Add(fSmearTShieldCo, fParameters[6]);
  fModelTotCo->Add(fSmear50mKCo,    fParameters[6]);
  fModelTotCo->Add(fSmear600mKCo,   fParameters[6]);
  fModelTotCo->Add(fSmearIVCCo,     fParameters[7]);
  fModelTotCo->Add(fSmearOVCCo,     fParameters[7]);

  fModelTotNDBD->Add(fSmearNDBD,    fParameters[9]);
  fModelTotBi->Add(fSmearBi,      fParameters[10]);

  ////////////////////////////////////////////////////////
*/

	
	fModelTot->SetLineColor(2);
	fModelTot->Draw("SAME");

	fModelTotTh->SetLineColor(3);
	fModelTotTh->SetLineStyle(2);
	fModelTotRa->SetLineColor(4);
	fModelTotRa->SetLineStyle(2);
	fModelTotK->SetLineColor(6);
	fModelTotK->SetLineStyle(2);
	fModelTotCo->SetLineColor(7);
	fModelTotCo->SetLineStyle(2);
	fModelTotNDBD->SetLineColor(42);
	fModelTotNDBD->SetLineStyle(2);
  fModelTotBi->SetLineColor(5);
  fModelTotBi->SetLineStyle(2);

	fModelTotTh->Draw("SAME");
	fModelTotRa->Draw("SAME");
	fModelTotK->Draw("SAME");
	fModelTotCo->Draw("SAME");
	fModelTotNDBD->Draw("SAME");
  fModelTotBi->Draw("SAME");



  // Few Parameters
  TPaveText *pt = new TPaveText(0.35,0.78,0.70,0.98,"NB NDC");
  pt->AddText(Form("Fit Range: %.0f to %.0f keV -- #chi^{2}/NDF: %0.3f", dFitMin, dFitMax, (GetChiSquare()/((dFitMax-dFitMin)/dBinSize - dNumParameters)) ));
  pt->AddText(Form("Close Th: %0.2E#pm%0.2E --- Far Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
  pt->AddText(Form("Close Ra: %0.2E#pm%0.2E --- Far Ra: %0.2E#pm%0.2E", fParameters[2], fParError[2], fParameters[3], fParError[3] ));
  pt->AddText(Form("Close K: %0.2E#pm%0.2E --- Far K: %0.2E#pm%0.2E", fParameters[4], fParError[4], fParameters[5], fParError[5] ));
  pt->AddText(Form("Close Co: %0.2E#pm%0.2E --- Far Co: %0.2E#pm%0.2E", fParameters[6], fParError[6], fParameters[7], fParError[7] ));
  pt->AddText(Form("Bi-207: %0.2E#pm%0.2E --- NDBD: %0.2E#pm%0.2E", fParameters[10], fParError[10], fParameters[9], fParError[9] ));

/*
  // Many Parameters
  TPaveText *pt = new TPaveText(0.35,0.77,0.70,0.99,"NB NDC");
  pt->AddText(Form("Fit Range: %.0f to %.0f keV -- #chi^{2}/NDF: %0.3f", dFitMin, dFitMax, (GetChiSquare()/((dFitMax-dFitMin)/dBinSize - dNumParameters)) ));
  // pt->AddText(Form("#chi^{2}/NDF: %0.3f --  Resolution %0.4f", (GetChiSquare()/((dFitMax-dFitMin)/dBinSize - dNumParameters)) ,fParameters[8]));
  pt->AddText(Form("Frame Th: %0.2E#pm%0.2E --- TShield Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[11], fParError[11] ));
  pt->AddText(Form("50mK Th: %0.2E#pm%0.2E --- 600mK Th: %0.2E#pm%0.2E", fParameters[12], fParError[12], fParameters[13], fParError[13] ));
  pt->AddText(Form("IVC Th: %0.2E#pm%0.2E --- OVC Th: %0.2E#pm%0.2E", fParameters[1], fParError[1], fParameters[17], fParError[17] ));
  pt->AddText(Form("Frame Ra: %0.2E#pm%0.2E --- TShield Ra: %0.2E#pm%0.2E", fParameters[2], fParError[2], fParameters[14], fParError[14] ));
  pt->AddText(Form("50mK Ra: %0.2E#pm%0.2E --- 600mK Ra: %0.2E#pm%0.2E", fParameters[15], fParError[15], fParameters[16], fParError[16] ));
  pt->AddText(Form("IVC Ra: %0.2E#pm%0.2E --- OVC Ra: %0.2E#pm%0.2E", fParameters[3], fParError[3], fParameters[18], fParError[18] ));
  pt->AddText(Form("Close K: %0.2E#pm%0.2E --- Far K: %0.2E#pm%0.2E", fParameters[4], fParError[4], fParameters[5], fParError[5] ));
  pt->AddText(Form("Close Co: %0.2E#pm%0.2E --- Far Co: %0.2E#pm%0.2E", fParameters[6], fParError[6], fParameters[7], fParError[7] ));
  pt->AddText(Form("Bi-207: %0.2E#pm%0.2E --- NDBD: %0.2E#pm%0.2E", fParameters[10], fParError[10], fParameters[9], fParError[9] ));
*/




  // pt->AddText(Form(,  ));
  pt->Draw();

 	TLegend *legfit = new TLegend(0.82,0.82,0.95,0.95);
 	legfit->AddEntry(fModelTot, "Total PDF", "l");
 	legfit->AddEntry(fModelTotTh, "Total Th-232", "l");
  legfit->AddEntry(fModelTotRa, "Total Ra-226", "l");
 	legfit->AddEntry(fModelTotK, "Total K-40", "l");
 	legfit->AddEntry(fModelTotCo, "Total Co-60", "l");
 	legfit->AddEntry(fModelTotNDBD, "NDBD", "l");
  legfit->AddEntry(fModelTotBi, "Bi-207", "l");

 	legfit->Draw();


	// Residuals
	TCanvas *cResidual = new TCanvas("cResidual", "cResidual", 1200, 800);
	hResidualDist = CalculateResiduals(fModelTot, fDataHistoM1);

	hResidualDist->SetName("Residuals");
	hResidualDist->GetXaxis()->SetTitle("Energy (keV)");
	// hResidualDist->GetXaxis()->SetTitleSize(0.04);
	// hResidualDist->GetXaxis()->SetLabelSize(0.05);
	// hResidualDist->GetYaxis()->SetLabelSize(0.05);
	// hResidualDist->GetYaxis()->SetTitleSize(0.04);	
	hResidualDist->GetYaxis()->SetTitle("Residuals (#sigma)");

	hResidualDist->GetXaxis()->SetRange(dFitMin/dBinSize-5, dFitMax/dBinSize+5);
	hResidualDist->Draw("E");

	return true;
   
 }

// Draws background data, must Initialize first!
void TBackgroundModel::DrawBkg()
{

 	gStyle->SetOptStat(0);
  gStyle->SetOptFit();
 	// gStyle->SetOptTitle(0);	
  TCanvas *cBkg = new TCanvas("cBkg", "cBkg", 1200, 800);
  cBkg->SetLogy();

  if(dMult == 1)
  {
    fDataHistoM1->SetLineColor(1);
    fDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
    fDataHistoM1->GetYaxis()->SetTitle("Counts/(10 keV)/yr");
    fDataHistoM1->Draw();
  }
  else if(dMult == 2)
  {
    fDataHistoM2->SetLineColor(1);
    fDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
    fDataHistoM2->GetYaxis()->SetTitle("Counts/(10 keV)/yr");
    fDataHistoM2->Draw();
  }

}

void TBackgroundModel::DrawMC()
{
// Draws all MC spectra, must Initialize first!

 	gStyle->SetOptStat(0);
 	gStyle->SetOptTitle(0);

 	TLegend *legth = new TLegend(0.75,0.80,0.925,0.925);
 	TLegend *legra = new TLegend(0.75,0.80,0.925,0.925);
 	TLegend *legk = new TLegend(0.75,0.80,0.925,0.925);
 	TLegend *legco = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legndbd = new TLegend(0.8,0.890,0.925,0.925);
  TLegend *legbi = new TLegend(0.8,0.890,0.925,0.925);
  TLegend *legths = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legras = new TLegend(0.75,0.80,0.925,0.925);


    TCanvas *cThSurf = new TCanvas("cThSurf", "cThSurf", 1200, 800);
    cThSurf->SetLogy();

    fModelFrameTh->SetLineColor(1);
    fModelFrameThS01->SetLineColor(2);
    fModelFrameThS1->SetLineColor(3);
    fModelFrameThS10->SetLineColor(4);
    fModelFrameThS100->SetLineColor(6);

    fModelFrameTh->DrawNormalized();
    fModelFrameThS01->DrawNormalized("SAME");
    fModelFrameThS1->DrawNormalized("SAME");
    fModelFrameThS10->DrawNormalized("SAME");
    fModelFrameThS100->DrawNormalized("SAME");
    fModelFrameTh->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);


    legths->AddEntry(fModelFrameTh, "Frame Bulk" ,"l");
    legths->AddEntry(fModelFrameThS01, "Frame Surface 0.1 #mum", "l");
    legths->AddEntry(fModelFrameThS1, "Frame Surface 1 #mum" ,"l");
    legths->AddEntry(fModelFrameThS10, "Frame Surface 10 #mum" ,"l");
    legths->AddEntry(fModelFrameThS100, "Frame Surface 100 #mum" ,"l");
    legths->Draw();



    TCanvas *cRaSurf = new TCanvas("cRaSurf", "cRaSurf", 1200, 800);
    cRaSurf->SetLogy();

    fModelFrameRa->SetLineColor(1);
    fModelFrameRaS01->SetLineColor(2);
    fModelFrameRaS1->SetLineColor(3);
    fModelFrameRaS10->SetLineColor(4);
    fModelFrameRaS100->SetLineColor(6);

    fModelFrameRa->DrawNormalized();
    fModelFrameRaS01->DrawNormalized("SAME");
    fModelFrameRaS1->DrawNormalized("SAME");
    fModelFrameRaS10->DrawNormalized("SAME");
    fModelFrameRaS100->DrawNormalized("SAME");
    fModelFrameRa->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);


    legras->AddEntry(fModelFrameRa, "Frame Bulk" ,"l");
    legras->AddEntry(fModelFrameRaS01, "Frame Surface 0.1 #mum", "l");
    legras->AddEntry(fModelFrameRaS1, "Frame Surface 1 #mum" ,"l");
    legras->AddEntry(fModelFrameRaS10, "Frame Surface 10 #mum" ,"l");
    legras->AddEntry(fModelFrameRaS100, "Frame Surface 100 #mum" ,"l");
    legras->Draw();



    TCanvas *cTh232 = new TCanvas("cTh232", "cTh232", 1200, 800);
    cTh232->SetLogy();

    fModelFrameTh->SetLineColor(1);
    fModelTShieldTh->SetLineColor(2);
    fModel50mKTh->SetLineColor(3);
    fModel600mKTh->SetLineColor(4);
    fModelIVCTh->SetLineColor(6);
    fModelOVCTh->SetLineColor(7);


    fModelFrameTh->DrawNormalized();
    fModelTShieldTh->DrawNormalized("SAME");
    fModel50mKTh->DrawNormalized("SAME");
    fModel600mKTh->DrawNormalized("SAME");
    fModelIVCTh->DrawNormalized("SAME");
    fModelOVCTh->DrawNormalized("SAME");
    fModelFrameTh->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);


    legth->AddEntry(fModelFrameTh, "Frame" ,"l");
    legth->AddEntry(fModelTShieldTh, "TShield", "l");
    legth->AddEntry(fModel50mKTh, "50mK" ,"l");
    legth->AddEntry(fModel600mKTh, "600mK" ,"l");
    legth->AddEntry(fModelIVCTh, "IVC" ,"l");
    legth->AddEntry(fModelOVCTh, "OVC" ,"l");
    legth->Draw();

    TCanvas *cRa226 = new TCanvas("cRa226", "cRa226", 1200, 800);
    cRa226->SetLogy();

    fModelFrameRa->SetLineColor(1);
    fModelTShieldRa->SetLineColor(2);
    fModel50mKRa->SetLineColor(3);
    fModel600mKRa->SetLineColor(4);
    fModelIVCRa->SetLineColor(6);
    fModelOVCRa->SetLineColor(7);

    fModelFrameRa->DrawNormalized();
    fModelTShieldRa->DrawNormalized("SAME");
    fModel50mKRa->DrawNormalized("SAME");
    fModel600mKRa->DrawNormalized("SAME");
    fModelIVCRa->DrawNormalized("SAME");
    fModelOVCRa->DrawNormalized("SAME");
    fModelFrameRa->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);


    legra->AddEntry(fModelFrameRa, "Frame" ,"l");
    legra->AddEntry(fModelTShieldRa, "TShield", "l");
    legra->AddEntry(fModel50mKRa, "50mK" ,"l");
    legra->AddEntry(fModel600mKRa, "600mK" ,"l");
    legra->AddEntry(fModelIVCRa, "IVC" ,"l");
    legra->AddEntry(fModelOVCRa, "OVC" ,"l");
    legra->Draw();

    TCanvas *cK40 = new TCanvas("cK40", "cK40", 1200, 800);
    cK40->SetLogy();

    fModelFrameK->SetLineColor(1);
    fModelTShieldK->SetLineColor(2);
    fModel50mKK->SetLineColor(3);
    fModel600mKK->SetLineColor(4);
    fModelIVCK->SetLineColor(6);
    fModelOVCK->SetLineColor(7);

    fModelFrameK->DrawNormalized();
    fModelTShieldK->DrawNormalized("SAME");
    fModel50mKK->DrawNormalized("SAME");
    fModel600mKK->DrawNormalized("SAME");
    fModelIVCK->DrawNormalized("SAME");
    fModelOVCK->DrawNormalized("SAME");
    fModelFrameK->GetXaxis()->SetRange(0/dBinSize, 1600/dBinSize);


    legk->AddEntry(fModelFrameK, "Frame" ,"l");
    legk->AddEntry(fModelTShieldK, "TShield", "l");
    legk->AddEntry(fModel50mKK, "50mK" ,"l");
    legk->AddEntry(fModel600mKK, "600mK" ,"l");
    legk->AddEntry(fModelIVCK, "IVC" ,"l");
    legk->AddEntry(fModelOVCK, "OVC" ,"l");    
    legk->Draw();

    TCanvas *cCo60 = new TCanvas("cCo60", "cCo60", 1200, 800);
    cCo60->SetLogy();

    fModelFrameCo->SetLineColor(1);
    fModelTShieldCo->SetLineColor(2);
    fModel50mKCo->SetLineColor(3);
    fModel600mKCo->SetLineColor(4);
    fModelIVCCo->SetLineColor(6);
    fModelOVCCo->SetLineColor(7);

    fModelFrameCo->Draw();
    fModelTShieldCo->Draw("SAME");
    fModel50mKCo->Draw("SAME");
    fModel600mKCo->Draw("SAME");
    fModelIVCCo->Draw("SAME");
    fModelOVCCo->Draw("SAME");
    fModelFrameCo->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);


    legco->AddEntry(fModelFrameCo, "Frame" ,"l");
    legco->AddEntry(fModelTShieldCo, "TShield", "l");
    legco->AddEntry(fModel50mKCo, "50mK" ,"l");
    legco->AddEntry(fModel600mKCo, "600mK" ,"l");
    legco->AddEntry(fModelIVCCo, "IVC" ,"l");
    legco->AddEntry(fModelOVCCo, "OVC" ,"l");
    legco->Draw();



    TCanvas *cNDBD = new TCanvas("cNDBD", "cNDBD", 1200, 800);
    cNDBD->SetLogy();
    fModelNDBD->SetLineColor(1);
    fModelNDBD->Draw();
    fModelNDBD->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);
    legndbd->AddEntry(fModelNDBD, "0#nu#beta#beta" ,"l");
    legndbd->Draw();

    TCanvas *cBi = new TCanvas("cBi", "cBi", 1200, 800);
    cBi->SetLogy();
    fModelBi->SetLineColor(1);
    fModelBi->Draw();
    fModelBi->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);
    legbi->AddEntry(fModelBi, "Bi-207" ,"l");
    legbi->Draw();

}


// Generates Toy Data using MC histograms
void TBackgroundModel::GenerateToyData()
{
	bToyFit = true;
  Initialize();
	// Create some RNG for weights?
	// Currently just put in by hand

  fToyDataThTot->Add( SmearMC(fModelFrameTh, fToyDataTh, 5),   3000);
  fToyDataThTot->Add( SmearMC(fModelTShieldTh, fToyDataTh, 5), 3000);  
  fToyDataThTot->Add( SmearMC(fModel50mKTh, fToyDataTh, 5),    3000);
  fToyDataThTot->Add( SmearMC(fModel600mKTh, fToyDataTh, 5),   3000);
  fToyDataThTot->Add( SmearMC(fModelIVCTh, fToyDataTh, 5),     1500);
  fToyDataThTot->Add( SmearMC(fModelOVCTh, fToyDataTh, 5),     1500);

  fToyDataRaTot->Add( SmearMC(fModelFrameRa, fToyDataRa, 5),    500);
  fToyDataRaTot->Add( SmearMC(fModelTShieldRa, fToyDataRa, 5),  500);  
  fToyDataRaTot->Add( SmearMC(fModel50mKRa, fToyDataRa, 5),     500);
  fToyDataRaTot->Add( SmearMC(fModel600mKRa, fToyDataRa, 5),    500);
  fToyDataRaTot->Add( SmearMC(fModelIVCRa, fToyDataRa, 5),      100);
  fToyDataRaTot->Add( SmearMC(fModelOVCRa, fToyDataRa, 5),      100);

  fToyDataCoTot->Add( SmearMC(fModelFrameK, fToyDataCo, 5),     900);
  fToyDataCoTot->Add( SmearMC(fModelTShieldK, fToyDataCo, 5),   900);
  fToyDataCoTot->Add( SmearMC(fModel50mKK, fToyDataCo, 5),      900);
  fToyDataCoTot->Add( SmearMC(fModel600mKK, fToyDataCo, 5),     900);
  fToyDataCoTot->Add( SmearMC(fModelIVCK, fToyDataCo, 5),       50);
  fToyDataCoTot->Add( SmearMC(fModelOVCK, fToyDataCo, 5),       50); 

  fToyDataKTot->Add( SmearMC(fModelFrameCo, fToyDataK, 5),      400);
  fToyDataKTot->Add( SmearMC(fModelTShieldCo, fToyDataK, 5),    400);
  fToyDataKTot->Add( SmearMC(fModel50mKCo, fToyDataK, 5),       400);
  fToyDataKTot->Add( SmearMC(fModel600mKCo, fToyDataK, 5),      400);
  fToyDataKTot->Add( SmearMC(fModelIVCCo, fToyDataK, 5),        800);
  fToyDataKTot->Add( SmearMC(fModelOVCCo, fToyDataK, 5),        800);  

	cout << "Loaded Toy Histograms" << endl;


	fToyData->Add(fToyDataThTot);
	fToyData->Add(fToyDataRaTot);
	// fToyData->Add(fToyDataCo,	1);
	// fToyData->Add(fToyDataK,	1);

}


// Calculates ChiSquare... model parameters not set here!
double TBackgroundModel::GetChiSquare()
{
	//cout<<"Calling GetChiSquare()"<<endl;
	double chiSquare = 0.;
	double datam1_i, errm1_i;
  double datam2_i, errm2_i;
  double modelm1_i, modelm2_i;	

	for(int i = dFitMin/dBinSize+1; i <= dFitMax/dBinSize; i++)
	{
		if(bToyFit)
		{
			datam1_i = fToyData->GetBinContent(i); // For Toy data
		}
		else 
		{
      if(dMult == 1)
      {
			 datam1_i = fDataHistoM1->GetBinContent(i); // For real data
      }
      else if(dMult == 2)
      {
       datam2_i = fDataHistoM2->GetBinContent(i); // For real data
      }      
		}

		modelm1_i = fModelTot->GetBinContent(i);


		// Log-likelihood Chi-Squared
    // Avoiding 0's... correct or no?
		if(modelm1_i != 0 && datam1_i != 0)
		{
			chiSquare += 2 * (modelm1_i - datam1_i + datam1_i * TMath::Log(datam1_i/modelm1_i));
      // Adding on M2 portion

		}


		// Neyman chi-squared
/*		
		err_i = sqrt(data_i);	// Assuming no data bins are 0!
		if(err_i>0)
		{
			chiSquare += pow(data_i - model_i,2)/pow(err_i,2);
		}
		else
		{

			chiSquare+=pow(data_i - model_i,2);	
		}
*/


	}

	return chiSquare;
}

// Gets efficiency of MC as fraction of normalization
double TBackgroundModel::GetMCEff(TH1D *h1)
{
  double fIntTotal;
  double fIntRange;
  double dEff;

  fIntTotal = h1->GetIntegral(500/dBinSize, 2650/dBinSize);
  fIntRange = h1->GetIntegral(dFitMin/dBinSize, dFitMax/dBinSize);

  dEff = fIntRange/fIntTotal;

  return dEff;
}


void TBackgroundModel::Initialize()
{	
	// Loading background data
	LoadData();	


  // Fills and Loads MC data
  outTreeFrameTh 		= LoadMC(dDataDir.c_str(),	"Frame", 	"Th232", "B", dMult);
  outTreeTShieldTh 	= LoadMC(dDataDir.c_str(),	"TShield","Th232", "B", dMult);
  outTree50mKTh 		= LoadMC(dDataDir.c_str(),	"50mK",		"Th232", "B", dMult);
  outTree600mKTh 		= LoadMC(dDataDir.c_str(),	"600mK", 	"Th232", "B", dMult);
  outTreeIVCTh 	  	= LoadMC(dDataDir.c_str(),	"IVC", 		"Th232", "B", dMult);
  outTreeOVCTh 	  	= LoadMC(dDataDir.c_str(),	"OVC", 		"Th232", "B", dMult);

  outTreeFrameRa	 	= LoadMC(dDataDir.c_str(),	"Frame", 	"Ra226", "B", dMult);
  outTreeTShieldRa 	= LoadMC(dDataDir.c_str(),	"TShield","Ra226", "B", dMult);    
  outTree50mKRa	  	= LoadMC(dDataDir.c_str(),	"50mK", 	"Ra226", "B", dMult);
  outTree600mKRa		= LoadMC(dDataDir.c_str(),	"600mK", 	"Ra226", "B", dMult);
  outTreeIVCRa	   	= LoadMC(dDataDir.c_str(),	"IVC", 		"Ra226", "B", dMult);
  outTreeOVCRa 	  	= LoadMC(dDataDir.c_str(),	"OVC", 		"Ra226", "B", dMult);

  outTreeFrameK 		= LoadMC(dDataDir.c_str(),	"Frame", 	"K40", "B", dMult);
  outTreeTShieldK 	= LoadMC(dDataDir.c_str(),	"TShield","K40", "B", dMult);    
  outTree50mKK	   	= LoadMC(dDataDir.c_str(),	"50mK", 	"K40", "B", dMult);
  outTree600mKK	  	= LoadMC(dDataDir.c_str(),	"600mK", 	"K40", "B", dMult);
  outTreeIVCK		   	= LoadMC(dDataDir.c_str(),	"IVC", 		"K40", "B", dMult);
  outTreeOVCK 	   	= LoadMC(dDataDir.c_str(),	"OVC", 		"K40", "B", dMult);


  outTreeFrameCo 		= LoadMC(dDataDir.c_str(),	"Frame", 	"Co60",	"B", dMult);
  outTreeTShieldCo 	= LoadMC(dDataDir.c_str(),	"TShield","Co60", "B", dMult);    
  outTree50mKCo	  	= LoadMC(dDataDir.c_str(),	"50mK", 	"Co60", "B", dMult);
  outTree600mKCo		= LoadMC(dDataDir.c_str(),	"600mK", 	"Co60", "B", dMult);
  outTreeIVCCo	   	= LoadMC(dDataDir.c_str(),	"IVC", 		"Co60", "B", dMult);
  outTreeOVCCo 	  	= LoadMC(dDataDir.c_str(),	"OVC", 		"Co60", "B", dMult);

  outTreeNDBD 	   	= LoadMC(dDataDir.c_str(),	"Crystal", "0NDBD", "B", dMult);
  outTreeBi         = LoadMC(dDataDir.c_str(),  "RLead",   "Bi207", "B", dMult);

  outTreeFrameThS01   = LoadMC(dDataDir.c_str(),  "Frame",  "Th232", "S01", dMult);
  outTreeFrameThS1    = LoadMC(dDataDir.c_str(),  "Frame",  "Th232", "S1", dMult);
  outTreeFrameThS10   = LoadMC(dDataDir.c_str(),  "Frame",  "Th232", "S10", dMult);
  outTreeFrameThS100  = LoadMC(dDataDir.c_str(),  "Frame",  "Th232", "S100", dMult);

  outTreeFrameRaS01   = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S01", dMult);
  outTreeFrameRaS1    = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S1", dMult);
  outTreeFrameRaS10   = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S10", dMult);
  outTreeFrameRaS100  = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S100", dMult);



  // Surface
  outTreeFrameThS01->Project("fModelFrameThS01", "Ener1", ener_cut);
  outTreeFrameThS1->Project("fModelFrameThS1", "Ener1", ener_cut);
  outTreeFrameThS10->Project("fModelFrameThS10", "Ener1", ener_cut);
  outTreeFrameThS100->Project("fModelFrameThS100", "Ener1", ener_cut);

  outTreeFrameRaS01->Project("fModelFrameRaS01", "Ener1", ener_cut);
  outTreeFrameRaS1->Project("fModelFrameRaS1", "Ener1", ener_cut);
  outTreeFrameRaS10->Project("fModelFrameRaS10", "Ener1", ener_cut);
  outTreeFrameRaS100->Project("fModelFrameRaS100", "Ener1", ener_cut);

  outTreeNDBD->Project("fModelNDBD",				"Ener1", ener_cut);
  outTreeBi->Project("fModelBi",            "Ener1", ener_cut);  

	outTreeFrameTh->Project("fModelFrameTh", 		"Ener1", ener_cut);
	outTreeTShieldTh->Project("fModelTShieldTh",	"Ener1", ener_cut);
  outTree50mKTh->Project("fModel50mKTh", 			"Ener1", ener_cut);
  outTree600mKTh->Project("fModel600mKTh", 		"Ener1", ener_cut);
  outTreeIVCTh->Project("fModelIVCTh", 			"Ener1", ener_cut);
  outTreeOVCTh->Project("fModelOVCTh", 			"Ener1", ener_cut);

	outTreeFrameRa->Project("fModelFrameRa", 		"Ener1", ener_cut);
	outTreeTShieldRa->Project("fModelTShieldRa",	"Ener1", ener_cut);	
  outTree50mKRa->Project("fModel50mKRa", 			"Ener1", ener_cut);
  outTree600mKRa->Project("fModel600mKRa", 		"Ener1", ener_cut);
  outTreeIVCRa->Project("fModelIVCRa", 			"Ener1", ener_cut);
  outTreeOVCRa->Project("fModelOVCRa", 			"Ener1", ener_cut);

	outTreeFrameK->Project("fModelFrameK", 			"Ener1", ener_cut);
	outTreeTShieldK->Project("fModelTShieldK",		"Ener1", ener_cut);	
  outTree50mKK->Project("fModel50mKK", 			"Ener1", ener_cut);
  outTree600mKK->Project("fModel600mKK", 			"Ener1", ener_cut);
  outTreeIVCK->Project("fModelIVCK", 				"Ener1", ener_cut);
  outTreeOVCK->Project("fModelOVCK", 				"Ener1", ener_cut);	

	outTreeFrameCo->Project("fModelFrameCo", 		"Ener1", ener_cut);
	outTreeTShieldCo->Project("fModelTShieldCo",	"Ener1", ener_cut);	
  outTree50mKCo->Project("fModel50mKCo", 			"Ener1", ener_cut);
  outTree600mKCo->Project("fModel600mKCo", 		"Ener1", ener_cut);
  outTreeIVCCo->Project("fModelIVCCo", 			"Ener1", ener_cut);
  outTreeOVCCo->Project("fModelOVCCo", 			"Ener1", ener_cut);




	cout << "Loaded MC" << endl;


  // Will be for detector efficiencies... not added yet
/*  
  fMCEff[0] = GetMCEff(fModelFrameTh);
  fMCEff[1] = GetMCEff(fModelTShieldTh);
  fMCEff[2] = GetMCEff(fModel50mKTh);
  fMCEff[3] = GetMCEff(fModel600mKTh);
  fMCEff[4] = GetMCEff(fModelIVCTh);
  fMCEff[5] = GetMCEff(fModelOVCTh);

  fMCEff[6] = GetMCEff(fModelFrameRa);
  fMCEff[7] = GetMCEff(fModelTShieldRa);
  fMCEff[8] = GetMCEff(fModel50mKRa);
  fMCEff[9] = GetMCEff(fModel600mKRa);
  fMCEff[10] = GetMCEff(fModelIVCRa);
  fMCEff[11] = GetMCEff(fModelOVCRa);

  fMCEff[12] = GetMCEff(fModelFrameK);
  fMCEff[13] = GetMCEff(fModelTShieldK);
  fMCEff[14] = GetMCEff(fModel50mKK);
  fMCEff[15] = GetMCEff(fModel600mKK);
  fMCEff[16] = GetMCEff(fModelIVCK);
  fMCEff[17] = GetMCEff(fModelOVCK);

  fMCEff[18] = GetMCEff(fModelFrameCo);
  fMCEff[19] = GetMCEff(fModelTShieldCo);
  fMCEff[20] = GetMCEff(fModel50mKCo);
  fMCEff[21] = GetMCEff(fModel600mKCo);
  fMCEff[22] = GetMCEff(fModelIVCCo);
  fMCEff[23] = GetMCEff(fModelOVCCo);

  fMCEff[24] = GetMCEff(fModelNDBD);
  fMCEff[25] = GetMCEff(fModelBi);
*/

	// Normalize all MC histograms
  // Fixing normalization of NDBD from 2000 to 2650
	NormalizePDF(fModelFrameTh, outTreeFrameTh, 	50, 2700);
	NormalizePDF(fModelTShieldTh, outTreeTShieldTh,	50, 2700);
	NormalizePDF(fModel50mKTh, outTree50mKTh,		50, 2700);
	NormalizePDF(fModel600mKTh, outTree600mKTh,		50, 2700);
	NormalizePDF(fModelIVCTh, outTreeIVCTh,			50, 2700);
	NormalizePDF(fModelOVCTh, outTreeOVCTh,			50, 2700);

	NormalizePDF(fModelFrameRa,  outTreeFrameRa,	50, 2700);
	NormalizePDF(fModelTShieldRa, outTreeTShieldRa,	50, 2700);	
	NormalizePDF(fModel50mKRa, outTree50mKRa,		50, 2700);
	NormalizePDF(fModel600mKRa, outTree600mKRa,		50, 2700);
	NormalizePDF(fModelIVCRa, outTreeIVCRa,			50, 2700);
	NormalizePDF(fModelOVCRa, outTreeOVCRa,			50, 2700);

	NormalizePDF(fModelFrameK, 	outTreeFrameK,		50, 2700);
	NormalizePDF(fModelTShieldK, outTreeTShieldK,	50, 2700);	
	NormalizePDF(fModel50mKK, outTree50mKK,			50, 2700);
	NormalizePDF(fModel600mKK, outTree600mKK,		50, 2700);
	NormalizePDF(fModelIVCK, outTreeIVCK,			50, 2700);
	NormalizePDF(fModelOVCK, outTreeOVCK,			50, 2700);

	NormalizePDF(fModelFrameCo, outTreeFrameCo,		50, 2700);
	NormalizePDF(fModelTShieldCo, outTreeTShieldCo,	50, 2700);	
	NormalizePDF(fModel50mKCo, outTree50mKCo,		50, 2700);
	NormalizePDF(fModel600mKCo, outTree600mKCo,		50, 2700);
	NormalizePDF(fModelIVCCo, outTreeIVCCo,			50, 2700);
	NormalizePDF(fModelOVCCo, outTreeOVCCo,			50, 2700);

  NormalizePDF(fModelNDBD, outTreeNDBD,     50, 2700);

  NormalizePDF(fModelBi, outTreeBi,         50, 2700);


  NormalizePDF(fModelFrameThS01,   outTreeFrameThS01, 50, 2700);
  NormalizePDF(fModelFrameThS1,    outTreeFrameThS1, 50, 2700);
  NormalizePDF(fModelFrameThS10,   outTreeFrameThS10, 50, 2700);
  NormalizePDF(fModelFrameThS100,  outTreeFrameThS100, 50, 2700);

  NormalizePDF(fModelFrameRaS01,   outTreeFrameRaS01, 50, 2700);
  NormalizePDF(fModelFrameRaS1,    outTreeFrameRaS1, 50, 2700);
  NormalizePDF(fModelFrameRaS10,   outTreeFrameRaS10, 50, 2700);
  NormalizePDF(fModelFrameRaS100,  outTreeFrameRaS100, 50, 2700);



	cout << "Normalized MC PDFs" << endl;

  if(bFixedRes)
  {
    // cout << "Fixed resolution: " << endl;
    cout << "Smearing histograms with constant 6 keV resolution" << endl;

    // double dRes = 4;
    double dRes = 6.0/2.355;

    // Adding the 10 micron distribution for now...
    SmearMC(fModelFrameThS10, fSmearFrameThS10, dRes);
    // SmearMC();

    SmearMC(fModelFrameTh, fSmearFrameTh, dRes);
    SmearMC(fModelTShieldTh, fSmearTShieldTh, dRes);  
    SmearMC(fModel50mKTh, fSmear50mKTh, dRes);
    SmearMC(fModel600mKTh, fSmear600mKTh, dRes);
    SmearMC(fModelIVCTh, fSmearIVCTh, dRes);
    SmearMC(fModelOVCTh, fSmearOVCTh, dRes);

    SmearMC(fModelFrameRa, fSmearFrameRa, dRes);
    SmearMC(fModelTShieldRa, fSmearTShieldRa, dRes);  
    SmearMC(fModel50mKRa, fSmear50mKRa, dRes);
    SmearMC(fModel600mKRa, fSmear600mKRa, dRes);
    SmearMC(fModelIVCRa, fSmearIVCRa, dRes);
    SmearMC(fModelOVCRa, fSmearOVCRa, dRes);

    SmearMC(fModelFrameK, fSmearFrameK, dRes);
    SmearMC(fModelTShieldK, fSmearTShieldK, dRes);
    SmearMC(fModel50mKK, fSmear50mKK, dRes);
    SmearMC(fModel600mKK, fSmear600mKK, dRes);
    SmearMC(fModelIVCK, fSmearIVCK, dRes);
    SmearMC(fModelOVCK, fSmearOVCK, dRes); 

    SmearMC(fModelFrameCo, fSmearFrameCo, dRes);
    SmearMC(fModelTShieldCo, fSmearTShieldCo, dRes);
    SmearMC(fModel50mKCo, fSmear50mKCo, dRes);
    SmearMC(fModel600mKCo, fSmear600mKCo, dRes);
    SmearMC(fModelIVCCo, fSmearIVCCo, dRes);
    SmearMC(fModelOVCCo, fSmearOVCCo, dRes);  

    SmearMC(fModelNDBD, fSmearNDBD, dRes);  

    SmearMC(fModelBi, fSmearBi, dRes);  

    cout << "Finished smearing MC histograms" << endl;

  }

}


// Loads the data
void TBackgroundModel::LoadData()
{
	if(fDataHistoTot == NULL) 
	{
		cout << "Data Histograms Not Created" << endl;
		return;
	}
	else
	{
		cout << "Data Histograms Created" << endl;
	}

  qtree->Add("/Users/brian/macros/CUOREZ/Bkg/ReducedBkg-ds*.root");	
  qtree->Project("fDataHistoTot", "Energy", base_cut);
  qtree->Project("fDataHistoM1", 	"Energy", base_cut && "Multiplicity_OFTime==1");
  qtree->Project("fDataHistoM2", 	"Energy", base_cut && "Multiplicity_OFTime==2");

	cout << "Loaded Data" << endl;

	// Normalizing data (don't!)
	// bin 0 = underflow, bin dNBins = last bin with upper-edge xup Excluded
  if(dMult == 1)
  {
	  dDataIntegral = fDataHistoM1->Integral(1, dNBins);
  }
  else if (dMult == 2)
  {
    dDataIntegral = fDataHistoM2->Integral(1, dNBins);
  }  

  int dDataIntegralTot = qtree->GetEntries();

  cout << "Total Events in background spectrum: " << dDataIntegralTot << endl; 
	cout << "Events in background spectrum: " << dDataIntegral << " Multiplicity: " << dMult << endl;

  // Scale by Live-time (ds 2061 - 2100) 14647393.0 seconds
  fDataHistoM1->Scale(1/(14647393.0 * dSecToYears));
  fDataHistoM2->Scale(1/(14647393.0 * dSecToYears));  

  cout << "Normalized Data" << endl;

}


// Loads MC files into Trees
TChain *TBackgroundModel::LoadMC(std::string dDir, std::string dLocation, std::string dSource, std::string dSType, int dMult)
{
    TChain *outTree = new TChain("outTree");
    outTree->Add(Form("%s%s-%s-%s-M%d-T50-r0.0425.root", dDir.c_str(), dLocation.c_str(), dSource.c_str(), dSType.c_str(), dMult));

    return outTree;
}


// Normalize histogram
void TBackgroundModel::NormalizePDF(TH1D *h1, TChain *hChain, int minE, int maxE)
{
	double dIntegral = 0;
	double Time = 0;

	// hChain->SetBranchAddress("Time", &Time);

	// int dEvents = hChain->GetEntries();
	// hChain->GetEntry(dEvents - 1);

	// bin 0 = underflow, bin dNBins = last bin with upper-edge xup Excluded
	dIntegral = h1->Integral(minE/dBinSize, maxE/dBinSize);
	// cout << "Integral for " << h1->GetTitle() << " :" << dIntegral << endl;

  // How to apply as efficiency...
  // Integrate full range... 

	// Make sure integral isn't 0
  // If it is 0, clear model... 
  if(dIntegral == 0)
  {
    cout << Form("Integral of %s is 0, resetting histogram", h1->GetName()) << endl;
    h1->Reset();
  }

	if(dIntegral != 0)
	{
		h1->Scale(1/dIntegral);
	}

}

// Prints parameters, make sure to update
void TBackgroundModel::PrintParameters()
{	
	cout<< "Par0 = "	<< fParameters[0]	<< " +/- " << fParError[0] << endl;
	cout<< "Par1 = "	<< fParameters[1]	<< " +/- " << fParError[1] << endl;
	cout<< "Par2 = "	<< fParameters[2]	<< " +/- " << fParError[2] << endl;
	cout<< "Par3 = "	<< fParameters[3]	<< " +/- " << fParError[3] << endl;
	cout<< "Par4 = "	<< fParameters[4]	<< " +/- " << fParError[4] << endl;
	cout<< "Par5 = "	<< fParameters[5]	<< " +/- " << fParError[5] << endl;
	cout<< "Par6 = "	<< fParameters[6]	<< " +/- " << fParError[6] << endl;
	cout<< "Par7 = "	<< fParameters[7]	<< " +/- " << fParError[7] << endl;
	cout<< "Par8 = "	<< fParameters[8]	<< " +/- " << fParError[8] << endl;
	cout<< "Par9 = "	<< fParameters[9]	<< " +/- " << fParError[9] << endl;
  cout<< "Par10 = "  << fParameters[10] << " +/- " << fParError[10] << endl;
  cout<< "Par11 = "  << fParameters[11] << " +/- " << fParError[11] << endl;
  cout<< "Par11 = "  << fParameters[12] << " +/- " << fParError[12] << endl;
  cout<< "Par11 = "  << fParameters[13] << " +/- " << fParError[13] << endl;
  cout<< "Par11 = "  << fParameters[14] << " +/- " << fParError[14] << endl;
  cout<< "Par11 = "  << fParameters[15] << " +/- " << fParError[15] << endl;
  cout<< "Par11 = "  << fParameters[16] << " +/- " << fParError[16] << endl;
  cout<< "Par11 = "  << fParameters[17] << " +/- " << fParError[17] << endl;
  cout<< "Par11 = "  << fParameters[18] << " +/- " << fParError[18] << endl;


//	double dSum = fParameters[0] + fParameters[1] + fParameters[2] + fParameters[3]
//					+ fParameters[4] + fParameters[5] + fParameters[6] + fParameters[7] + fParameters[9];

	// cout << "Par11 (1 - Sum) = " << 1 - dSum << endl;
//	cout << "Sum = " << dSum << endl; 


}


// Set Parameters in Model
void TBackgroundModel::SetParameters(int index, double value)
{
	// Change the index max depending on model
	if(index > 24) cout << "Index too large" << endl;
	else fParameters[index] = value;

}


// For custom smearing with resolution, currently constant resolution 
TH1D *TBackgroundModel::SmearMC(TH1D *hMC, TH1D *hSMC, double resolution)
{
	// Reset previously smeared histogram
	hSMC->Reset();

	double dArea;
	double dSmearedValue;

	for(int i = 0; i<dNBins; i++)
	{
		for(int j = 0; j<dNBins; j++)
		{
			// Normalization of gaussian = (bsin size * Area of bin j in MC) / Sigma of bin j (fit function evaluated at bin center)
			dArea = dBinSize*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*resolution);

			// Set parameters of gaussian ... resolution floating in fit
			gaus->SetParameters(dArea, hMC->GetBinCenter(j), resolution);

			// Smeared contribution from gaussian centered at bin j for bin i 
			dSmearedValue = gaus->Eval(hSMC->GetBinCenter(i));

			// Fill bin i with contribution from gaussian centered at bin j
			hSMC->Fill(hSMC->GetBinCenter(i), dSmearedValue);
		}
	}
	
	return hSMC;
}



// Creates/updates the background model
void TBackgroundModel::UpdateModel()
{
	if(fModelTot == NULL) 
	{
		cout << "Model Histogram Not Created" << endl;
		return;
	}

	// Reset all bins in model histogram(s)
	fModelTot->Reset();

	// Create model


	////////////////////////////////////////
	// Few parameters
	////////////////////////////////////////
  dNumCalls++;
  if(dNumCalls%50==0)
  {
    cout << "Call #: "<< dNumCalls << endl;
  }


  // Efficiency = Integral over fit range/integral over entire range -> Normalize 
  if(bFixedRes)
  {



  /////////////////////////////////////
  //// Few Parameters ////////////////
  ////////////////////////////////////
  // Testing 10 micron Frame for now
  fModelTot->Add( fSmearFrameThS1,    fParameters[11]);

  fModelTot->Add( fSmearFrameTh,    fParameters[0]);
  fModelTot->Add( fSmearTShieldTh,  fParameters[0]);  
  fModelTot->Add( fSmear50mKTh,     fParameters[0]);
  fModelTot->Add( fSmear600mKTh,    fParameters[0]);
  fModelTot->Add( fSmearIVCTh,      fParameters[1]);
  fModelTot->Add( fSmearOVCTh,      fParameters[1]);

  fModelTot->Add( fSmearFrameRa,    fParameters[2]);
  fModelTot->Add( fSmearTShieldRa,  fParameters[2]);  
  fModelTot->Add( fSmear50mKRa,     fParameters[2]);
  fModelTot->Add( fSmear600mKRa,    fParameters[2]);
  fModelTot->Add( fSmearIVCRa,      fParameters[3]);
  fModelTot->Add( fSmearOVCRa,      fParameters[3]);

  fModelTot->Add( fSmearFrameK,    fParameters[4]);
  fModelTot->Add( fSmearTShieldK,  fParameters[4]);
  fModelTot->Add( fSmear50mKK,     fParameters[4]);
  fModelTot->Add( fSmear600mKK,    fParameters[4]);
  fModelTot->Add( fSmearIVCK,      fParameters[5]);
  fModelTot->Add( fSmearOVCK,      fParameters[5]); 

  fModelTot->Add( fSmearFrameCo,    fParameters[6]);
  fModelTot->Add( fSmearTShieldCo,  fParameters[6]);
  fModelTot->Add( fSmear50mKCo,     fParameters[6]);
  fModelTot->Add( fSmear600mKCo,    fParameters[6]);
  fModelTot->Add( fSmearIVCCo,      fParameters[7]);
  fModelTot->Add( fSmearOVCCo,      fParameters[7]);  

  fModelTot->Add( fSmearNDBD,      fParameters[9]);  
  fModelTot->Add( fSmearBi,      fParameters[10]);  


/*
  /////////////////////////////////////
  //// Many parameters
  ////////////////////////////////////

  fModelTot->Add( fSmearFrameTh,    fParameters[0]);
  fModelTot->Add( fSmearTShieldTh,  fParameters[11]);  
  fModelTot->Add( fSmear50mKTh,     fParameters[12]);
  fModelTot->Add( fSmear600mKTh,    fParameters[13]);
  fModelTot->Add( fSmearIVCTh,      fParameters[1]);
  fModelTot->Add( fSmearOVCTh,      fParameters[17]);

  fModelTot->Add( fSmearFrameRa,    fParameters[2]);
  fModelTot->Add( fSmearTShieldRa,  fParameters[14]);  
  fModelTot->Add( fSmear50mKRa,     fParameters[15]);
  fModelTot->Add( fSmear600mKRa,    fParameters[16]);
  fModelTot->Add( fSmearIVCRa,      fParameters[3]);
  fModelTot->Add( fSmearOVCRa,      fParameters[18]);

  fModelTot->Add( fSmearFrameK,    fParameters[4]);
  fModelTot->Add( fSmearTShieldK,  fParameters[4]);
  fModelTot->Add( fSmear50mKK,     fParameters[4]);
  fModelTot->Add( fSmear600mKK,    fParameters[4]);
  fModelTot->Add( fSmearIVCK,      fParameters[5]);
  fModelTot->Add( fSmearOVCK,      fParameters[5]); 

  fModelTot->Add( fSmearFrameCo,    fParameters[6]);
  fModelTot->Add( fSmearTShieldCo,  fParameters[6]);
  fModelTot->Add( fSmear50mKCo,     fParameters[6]);
  fModelTot->Add( fSmear600mKCo,    fParameters[6]);
  fModelTot->Add( fSmearIVCCo,      fParameters[7]);
  fModelTot->Add( fSmearOVCCo,      fParameters[7]);  

  fModelTot->Add( fSmearNDBD,      fParameters[9]);  

  fModelTot->Add( fSmearBi,      fParameters[10]);  
*/


  }






  else
  {

  // Not updated! (as of 8/27/2014)

	fModelTot->Add( SmearMC(fModelFrameTh, fSmearFrameTh, fParameters[8]), 		fParameters[0]);
	fModelTot->Add( SmearMC(fModelTShieldTh, fSmearTShieldTh, fParameters[8]), 	fParameters[0]);	
	fModelTot->Add( SmearMC(fModel50mKTh, fSmear50mKTh, fParameters[8]), 		fParameters[0]);
	fModelTot->Add( SmearMC(fModel600mKTh, fSmear600mKTh,fParameters[8]), 		fParameters[0]);
	fModelTot->Add( SmearMC(fModelIVCTh, fSmearIVCTh, fParameters[8]), 			fParameters[1]);
	fModelTot->Add( SmearMC(fModelOVCTh, fSmearOVCTh, fParameters[8]), 			fParameters[1]);

	fModelTot->Add( SmearMC(fModelFrameRa, fSmearFrameRa, fParameters[8]), 		fParameters[2]);
	fModelTot->Add( SmearMC(fModelTShieldRa, fSmearTShieldRa, fParameters[8]), 	fParameters[2]);	
	fModelTot->Add( SmearMC(fModel50mKRa, fSmear50mKRa, fParameters[8]), 		fParameters[2]);
	fModelTot->Add( SmearMC(fModel600mKRa, fSmear600mKRa, fParameters[8]), 		fParameters[2]);
	fModelTot->Add( SmearMC(fModelIVCRa, fSmearIVCRa, fParameters[8]), 			fParameters[3]);
	fModelTot->Add( SmearMC(fModelOVCRa, fSmearOVCRa, fParameters[8]), 			fParameters[3]);

	fModelTot->Add( SmearMC(fModelFrameK, fSmearFrameK, fParameters[8]), 		fParameters[4]);
	fModelTot->Add( SmearMC(fModelTShieldK, fSmearTShieldK, fParameters[8]), 	fParameters[4]);
	fModelTot->Add( SmearMC(fModel50mKK, fSmear50mKK, fParameters[8]), 			fParameters[4]);
	fModelTot->Add( SmearMC(fModel600mKK, fSmear600mKK, fParameters[8]), 		fParameters[4]);
	fModelTot->Add( SmearMC(fModelIVCK, fSmearIVCK, fParameters[8]), 			fParameters[5]);
	fModelTot->Add( SmearMC(fModelOVCK, fSmearOVCK, fParameters[8]), 			fParameters[5]); 

	fModelTot->Add( SmearMC(fModelFrameCo, fSmearFrameCo, fParameters[8]), 		fParameters[6]);
	fModelTot->Add( SmearMC(fModelTShieldCo, fSmearTShieldCo, fParameters[8]), 	fParameters[6]);
	fModelTot->Add( SmearMC(fModel50mKCo, fSmear50mKCo, fParameters[8]), 		fParameters[6]);
	fModelTot->Add( SmearMC(fModel600mKCo, fSmear600mKCo, fParameters[8]), 		fParameters[6]);
	fModelTot->Add( SmearMC(fModelIVCCo, fSmearIVCCo, fParameters[8]), 			fParameters[7]);
	fModelTot->Add( SmearMC(fModelOVCCo, fSmearOVCCo, fParameters[8]), 			fParameters[7]);	

	fModelTot->Add( SmearMC(fModelNDBD, fSmearNDBD, fParameters[8]), 			fParameters[9]);	

  fModelTot->Add( SmearMC(fModelBi, fSmearBi, fParameters[8]),      fParameters[10]);  
  }



}















/*


////////////// Just testing
TH1D *TBackgroundModel::SmearChannel(TChain *fChain, TH1D *fMC, TH1D *fSMC, int channel, double resolution)
{
  // Reset previously smeared histogram
  // fMC->Reset(); 
  fSMC->Reset();

  fChain->Project(fMC->GetName() , "Ener1", Form("Channel==%d", channel));

  double dArea = 0;
  double dSmearedValue = 0;

  for(int i = 0; i<dBin; i++)
  {
    for(int j = 0; j<dBin; j++)
    {
      // Normalization of gaussian = (bsin size * Area of bin j in MC) / Sigma of bin j (fit function evaluated at bin center)
      dArea = dBinSize*fMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*resolution);

      // Set parameters of gaussian ... resolution floating in fit
      gaus->SetParameters(dArea, fMC->GetBinCenter(j), resolution);

      // Smeared contribution from gaussian centered at bin j for bin i 
      dSmearedValue = gaus->Eval(fSMC->GetBinCenter(i));

      // Fill bin i with contribution from gaussian centered at bin j
      fSMC->Fill(fSMC->GetBinCenter(i), dSmearedValue);
    }
  }

  return fSMC;
}

// Test to see how long the initial smearing takes
void TBackgroundModel::TestSmear()
{
  // int i = 1;
  for (int i = 0; i<52; i++)
  {
    cout << "Channel: " << i << endl;
    hDummy = new TH1D(Form("hC%d", i), "", dBin, dEMin, dEMax);

    fSmearFrameTh->Add( SmearChannel(outTreeFrameTh, hDummy, hSmearDummy, i, fResolution[i]), 1 );
    fSmearFrameRa->Add( SmearChannel(outTreeFrameRa, hDummy, hSmearDummy, i, fResolution[i]), 1 );
    fSmearFrameCo->Add( SmearChannel(outTreeFrameCo, hDummy, hSmearDummy, i, fResolution[i]), 1 );
    fSmearFrameK->Add( SmearChannel(outTreeFrameK, hDummy, hSmearDummy, i, fResolution[i]), 1 );




  } 

  // TCanvas *ctest2 = new TCanvas ("ctest2");
  hMC->Draw();

}
*/





