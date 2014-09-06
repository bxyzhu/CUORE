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
  dMaxEnergy = 8000.;

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


  // Model histograms M1
  fModelFrameThS01M1   = new TH1D("fModelFrameThS01M1",  "Frame Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameThS1M1    = new TH1D("fModelFrameThS1M1",  "Frame Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameThS10M1   = new TH1D("fModelFrameThS10M1",  "Frame Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameThS100M1  = new TH1D("fModelFrameThS100M1",  "Frame Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameRaS01M1   = new TH1D("fModelFrameRaS01M1",  "Frame Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRaS1M1    = new TH1D("fModelFrameRaS1M1",  "Frame Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRaS10M1   = new TH1D("fModelFrameRaS10M1",  "Frame Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRaS100M1  = new TH1D("fModelFrameRaS100M1",  "Frame Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fModelTShieldThS01M1   = new TH1D("fModelTShieldThS01M1",  "TShield Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldThS1M1    = new TH1D("fModelTShieldThS1M1",  "TShield Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldThS10M1   = new TH1D("fModelTShieldThS10M1",  "TShield Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldThS100M1  = new TH1D("fModelTShieldThS100M1",  "TShield Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameThM1    = new TH1D("fModelFrameThM1",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldThM1  = new TH1D("fModelTShieldThM1","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fModel50mKThM1     = new TH1D("fModel50mKThM1",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKThM1    = new TH1D("fModel600mKThM1",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCThM1      = new TH1D("fModelIVCThM1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCThM1      = new TH1D("fModelOVCThM1",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameRaM1    = new TH1D("fModelFrameRaM1",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldRaM1  = new TH1D("fModelTShieldRaM1","TShield",  dNBins, dMinEnergy, dMaxEnergy);  
  fModel50mKRaM1     = new TH1D("fModel50mKRaM1",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKRaM1    = new TH1D("fModel600mKRaM1",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCRaM1      = new TH1D("fModelIVCRaM1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCRaM1      = new TH1D("fModelOVCRaM1",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameKM1     = new TH1D("fModelFrameKM1",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldKM1   = new TH1D("fModelTShieldKM1", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fModel50mKKM1      = new TH1D("fModel50mKKM1",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKKM1     = new TH1D("fModel600mKKM1",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCKM1       = new TH1D("fModelIVCKM1",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCKM1       = new TH1D("fModelOVCKM1",   "OVC",        dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameCoM1    = new TH1D("fModelFrameCoM1",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldCoM1  = new TH1D("fModelTShieldCoM1","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fModel50mKCoM1     = new TH1D("fModel50mKCoM1",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKCoM1    = new TH1D("fModel600mKCoM1",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCCoM1      = new TH1D("fModelIVCCoM1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCCoM1      = new TH1D("fModelOVCCoM1",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelNDBDM1       = new TH1D("fModelNDBDM1",     "NDBD",     dNBins, dMinEnergy, dMaxEnergy);
  fModelBiM1         = new TH1D("fModelBiM1",       "Bi207",    dNBins, dMinEnergy, dMaxEnergy);

  // Model histograms M2
  fModelFrameThS01M2   = new TH1D("fModelFrameThS01M2",  "Frame Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameThS1M2    = new TH1D("fModelFrameThS1M2",  "Frame Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameThS10M2   = new TH1D("fModelFrameThS10M2",  "Frame Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameThS100M2  = new TH1D("fModelFrameThS100M2",  "Frame Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameRaS01M2   = new TH1D("fModelFrameRaS01M2",  "Frame Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRaS1M2    = new TH1D("fModelFrameRaS1M2",  "Frame Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRaS10M2   = new TH1D("fModelFrameRaS10M2",  "Frame Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRaS100M2  = new TH1D("fModelFrameRaS100M2",  "Frame Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fModelTShieldThS01M2   = new TH1D("fModelTShieldThS01M2",  "TShield Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldThS1M2    = new TH1D("fModelTShieldThS1M2",  "TShield Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldThS10M2   = new TH1D("fModelTShieldThS10M2",  "TShield Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldThS100M2  = new TH1D("fModelTShieldThS100M2",  "TShield Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameThM2    = new TH1D("fModelFrameThM2",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldThM2  = new TH1D("fModelTShieldThM2","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fModel50mKThM2     = new TH1D("fModel50mKThM2",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKThM2    = new TH1D("fModel600mKThM2",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCThM2      = new TH1D("fModelIVCThM2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCThM2      = new TH1D("fModelOVCThM2",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameRaM2    = new TH1D("fModelFrameRaM2",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldRaM2  = new TH1D("fModelTShieldRaM2","TShield",  dNBins, dMinEnergy, dMaxEnergy);  
  fModel50mKRaM2     = new TH1D("fModel50mKRaM2",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKRaM2    = new TH1D("fModel600mKRaM2",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCRaM2      = new TH1D("fModelIVCRaM2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCRaM2      = new TH1D("fModelOVCRaM2",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameKM2     = new TH1D("fModelFrameKM2",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldKM2   = new TH1D("fModelTShieldKM2", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fModel50mKKM2      = new TH1D("fModel50mKKM2",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKKM2     = new TH1D("fModel600mKKM2",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCKM2       = new TH1D("fModelIVCKM2",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCKM2       = new TH1D("fModelOVCKM2",   "OVC",        dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameCoM2    = new TH1D("fModelFrameCoM2",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldCoM2  = new TH1D("fModelTShieldCoM2","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fModel50mKCoM2     = new TH1D("fModel50mKCoM2",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKCoM2    = new TH1D("fModel600mKCoM2",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCCoM2      = new TH1D("fModelIVCCoM2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCCoM2      = new TH1D("fModelOVCCoM2",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelNDBDM2       = new TH1D("fModelNDBDM2",     "NDBD",     dNBins, dMinEnergy, dMaxEnergy);
  fModelBiM2         = new TH1D("fModelBiM2",       "Bi207",    dNBins, dMinEnergy, dMaxEnergy);



  // Total model histograms M1
  fModelTotM1      = new TH1D("fModelTotM1",      "Frame",        dNBins, dMinEnergy, dMaxEnergy);  
  fModelTotThM1    = new TH1D("fModelTotThM1",    "Total Th232",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotRaM1    = new TH1D("fModelTotRaM1",    "Total Ra226",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotKM1     = new TH1D("fModelTotKM1",     "Total K40",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTotCoM1    = new TH1D("fModelTotCoM1",    "Total Co60",   dNBins, dMinEnergy, dMaxEnergy);

  fModelTotNDBDM1  = new TH1D("fModelTotNDBDM1",  "Total NDBD",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotBiM1    = new TH1D("fModelTotBiM1",   "Total Bi207",   dNBins, dMinEnergy, dMaxEnergy);

  // Total model histograms M2
  fModelTotM2      = new TH1D("fModelTotM2",      "Frame",        dNBins, dMinEnergy, dMaxEnergy);  
  fModelTotThM2    = new TH1D("fModelTotThM2",    "Total Th232",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotRaM2    = new TH1D("fModelTotRaM2",    "Total Ra226",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotKM2     = new TH1D("fModelTotKM2",     "Total K40",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTotCoM2    = new TH1D("fModelTotCoM2",    "Total Co60",   dNBins, dMinEnergy, dMaxEnergy);

  fModelTotNDBDM2  = new TH1D("fModelTotNDBDM2",  "Total NDBD",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotBiM2    = new TH1D("fModelTotBiM2",   "Total Bi207",   dNBins, dMinEnergy, dMaxEnergy);



  // Smearing
  gaus = new TF1("gaus","gaus(0)", dMinEnergy, dMaxEnergy);



  // Smeared Histograms M1
  fSmearDummyM1      = new TH1D("fSmearDummyM1",  "Dummy smeared",  dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameThS01M1   = new TH1D("fSmearFrameThS01M1",  "Frame Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameThS1M1    = new TH1D("fSmearFrameThS1M1",  "Frame Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameThS10M1   = new TH1D("fSmearFrameThS10M1",  "Frame Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameThS100M1  = new TH1D("fSmearFrameThS100M1",  "Frame Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameRaS01M1   = new TH1D("fSmearFrameRaS01M1",  "Frame Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRaS1M1    = new TH1D("fSmearFrameRaS1M1",  "Frame Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRaS10M1   = new TH1D("fSmearFrameRaS10M1",  "Frame Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRaS100M1  = new TH1D("fSmearFrameRaS100M1",  "Frame Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearTShieldThS01M1   = new TH1D("fSmearTShieldThS01M1",  "TShield Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldThS1M1    = new TH1D("fSmearTShieldThS1M1",  "TShield Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldThS10M1   = new TH1D("fSmearTShieldThS10M1",  "TShield Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldThS100M1  = new TH1D("fSmearTShieldThS100M1",  "TShield Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameThM1    = new TH1D("fSmearFrameThM1",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldThM1  = new TH1D("fSmearTShieldThM1","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fSmear50mKThM1     = new TH1D("fSmear50mKThM1",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fSmear600mKThM1    = new TH1D("fSmear600mKThM1",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCThM1      = new TH1D("fSmearIVCThM1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearOVCThM1      = new TH1D("fSmearOVCThM1",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameRaM1    = new TH1D("fSmearFrameRaM1",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldRaM1  = new TH1D("fSmearTShieldRaM1","TShield",  dNBins, dMinEnergy, dMaxEnergy);  
  fSmear50mKRaM1     = new TH1D("fSmear50mKRaM1",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fSmear600mKRaM1    = new TH1D("fSmear600mKRaM1",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCRaM1      = new TH1D("fSmearIVCRaM1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearOVCRaM1      = new TH1D("fSmearOVCRaM1",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameKM1     = new TH1D("fSmearFrameKM1",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldKM1   = new TH1D("fSmearTShieldKM1", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fSmear50mKKM1      = new TH1D("fSmear50mKKM1",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fSmear600mKKM1     = new TH1D("fSmear600mKKM1",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCKM1       = new TH1D("fSmearIVCKM1",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearOVCKM1       = new TH1D("fSmearOVCKM1",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameCoM1    = new TH1D("fSmearFrameCoM1",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldCoM1  = new TH1D("fSmearTShieldCoM1","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fSmear50mKCoM1     = new TH1D("fSmear50mKCoM1",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fSmear600mKCoM1    = new TH1D("fSmear600mKCoM1",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCCoM1      = new TH1D("fSmearIVCCoM1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearOVCCoM1      = new TH1D("fSmearOVCCoM1",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearNDBDM1       = new TH1D("fSmearNDBDM1",   "NDBD",       dNBins, dMinEnergy, dMaxEnergy);
  fSmearBiM1         = new TH1D("fSmearBiM1",     "Bi",         dNBins, dMinEnergy, dMaxEnergy);


  // Smeared Histograms M2
  fSmearDummyM2      = new TH1D("fSmearDummyM2",  "Dummy smeared",  dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameThS01M2   = new TH1D("fSmearFrameThS01M2",  "Frame Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameThS1M2    = new TH1D("fSmearFrameThS1M2",  "Frame Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameThS10M2   = new TH1D("fSmearFrameThS10M2",  "Frame Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameThS100M2  = new TH1D("fSmearFrameThS100M2",  "Frame Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameRaS01M2   = new TH1D("fSmearFrameRaS01M2",  "Frame Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRaS1M2    = new TH1D("fSmearFrameRaS1M2",  "Frame Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRaS10M2   = new TH1D("fSmearFrameRaS10M2",  "Frame Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRaS100M2  = new TH1D("fSmearFrameRaS100M2",  "Frame Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearTShieldThS01M2   = new TH1D("fSmearTShieldThS01M2",  "TShield Surface 0.1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldThS1M2    = new TH1D("fSmearTShieldThS1M2",  "TShield Surface 1 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldThS10M2   = new TH1D("fSmearTShieldThS10M2",  "TShield Surface 10 #mum",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldThS100M2  = new TH1D("fSmearTShieldThS100M2",  "TShield Surface 100 #mum",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameThM2    = new TH1D("fSmearFrameThM2",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldThM2  = new TH1D("fSmearTShieldThM2","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fSmear50mKThM2     = new TH1D("fSmear50mKThM2",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fSmear600mKThM2    = new TH1D("fSmear600mKThM2",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCThM2      = new TH1D("fSmearIVCThM2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearOVCThM2      = new TH1D("fSmearOVCThM2",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameRaM2    = new TH1D("fSmearFrameRaM2",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldRaM2  = new TH1D("fSmearTShieldRaM2","TShield",  dNBins, dMinEnergy, dMaxEnergy);  
  fSmear50mKRaM2     = new TH1D("fSmear50mKRaM2",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fSmear600mKRaM2    = new TH1D("fSmear600mKRaM2",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCRaM2      = new TH1D("fSmearIVCRaM2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearOVCRaM2      = new TH1D("fSmearOVCRaM2",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameKM2     = new TH1D("fSmearFrameKM2",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldKM2   = new TH1D("fSmearTShieldKM2", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fSmear50mKKM2      = new TH1D("fSmear50mKKM2",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fSmear600mKKM2     = new TH1D("fSmear600mKKM2",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCKM2       = new TH1D("fSmearIVCKM2",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearOVCKM2       = new TH1D("fSmearOVCKM2",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameCoM2    = new TH1D("fSmearFrameCoM2",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearTShieldCoM2  = new TH1D("fSmearTShieldCoM2","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fSmear50mKCoM2     = new TH1D("fSmear50mKCoM2",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fSmear600mKCoM2    = new TH1D("fSmear600mKCoM2",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCCoM2      = new TH1D("fSmearIVCCoM2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearOVCCoM2      = new TH1D("fSmearOVCCoM2",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearNDBDM2       = new TH1D("fSmearNDBDM2",   "NDBD",       dNBins, dMinEnergy, dMaxEnergy);
  fSmearBiM2         = new TH1D("fSmearBiM2",     "Bi",         dNBins, dMinEnergy, dMaxEnergy);




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
   // minuit.DefineParameter(0, "Close Th",  100., 10.0, 0., 60000);
   minuit.DefineParameter(0, "Surface Th",  3000, 100.0, 0., 500000);   
   // minuit.DefineParameter(1, "Far Th",	 	100., 50.0, 0., 100000);
   minuit.DefineParameter(1, "Far Th",    45000., 10.0, 0., 100000);
   minuit.DefineParameter(2, "Close Ra",  600., 10.0, 0., 80000);   
   // minuit.DefineParameter(2, "Close Ra",  30000., 100.0, 0., 80000);   
   minuit.DefineParameter(3, "Far Ra",    55000., 10.0, 0., 80000);
   // minuit.DefineParameter(3, "Far Ra",    100., 100.0, 0., 80000);
   minuit.DefineParameter(4, "Close K", 	0., 100.0, 0., 500000);
   // minuit.DefineParameter(4, "Close K",   100., 10.0, 0., 500000);
   minuit.DefineParameter(5, "Far K",     0., 100.0, 0., 500000);
   // minuit.DefineParameter(5, "Far K", 		38000., 100.0, 0., 500000);
   minuit.DefineParameter(6, "Close Co", 	100., 100.0, 0., 80000); 
   minuit.DefineParameter(7, "Far Co",    11000, 100.0, 0., 80000);  
   // minuit.DefineParameter(7, "Far Co",	 	100., 100.0, 0., 50000);  
   minuit.DefineParameter(8, "Resolution",	6., 1, 3, 10);  
   minuit.DefineParameter(9, "NDBD",       65., 10.0, 0., 1000);     
   // minuit.DefineParameter(10, "Lead Bi",	 	6900., 100.0, 0., 100000);  
   minuit.DefineParameter(10, "Lead Bi",    0., 100.0, 0., 100000);  


/*
   ////////////////////////////////////////////////
   // Using less parameters WITH Surface 
   ////////////////////////////////////////////////
   // minuit.DefineParameter(0, "Close Th",  100., 10.0, 0., 60000);
   minuit.DefineParameter(0, "Surface Th",  220000, 100.0, 0., 500000);   
   // minuit.DefineParameter(1, "Far Th",   100., 50.0, 0., 100000);
   minuit.DefineParameter(1, "Far Th",    1000., 10.0, 0., 100000);
   minuit.DefineParameter(2, "Close Ra",  600., 10.0, 0., 80000);   
   // minuit.DefineParameter(2, "Close Ra",  30000., 100.0, 0., 80000);   
   minuit.DefineParameter(3, "Far Ra",    55000., 10.0, 0., 80000);
   // minuit.DefineParameter(3, "Far Ra",    100., 100.0, 0., 80000);
   // minuit.DefineParameter(4, "Close K",   0., 100.0, 0., 500000);
   minuit.DefineParameter(4, "Close K",   100., 10.0, 0., 500000);
   // minuit.DefineParameter(5, "Far K",     0., 100.0, 0., 500000);
   minuit.DefineParameter(5, "Far K",    38000., 100.0, 0., 500000);
   minuit.DefineParameter(6, "Close Co",  100., 100.0, 0., 80000); 
   minuit.DefineParameter(7, "Far Co",    11000, 100.0, 0., 80000);  
   // minuit.DefineParameter(7, "Far Co",   100., 100.0, 0., 50000);  
   minuit.DefineParameter(8, "Resolution",  6., 1, 3, 10);  
   minuit.DefineParameter(9, "NDBD",       54., 10.0, 0., 1000);     
   minuit.DefineParameter(10, "Lead Bi",   6900., 100.0, 0., 100000);  
   // minuit.DefineParameter(10, "Lead Bi",    0., 100.0, 0., 100000);  
*/



   // Close Th and Close Ra now split into its various sections, far Th and Ra still the same
   // This step after previous step full fit converges, just to see if any differences show up
   ////////////////////////////////////////////////
   // Using more parameters
   ////////////////////////////////////////////////
/*
   minuit.DefineParameter(0, "Frame Th",  100, 10.0, 0., 100000);   
   minuit.DefineParameter(1, "IVC Th",    100., 10.0, 0., 100000);
   minuit.DefineParameter(2, "Frame Ra",  100., 10.0, 0., 50000);   
   minuit.DefineParameter(3, "IVC Ra",    100., 10.0, 0., 80000);
   minuit.DefineParameter(4, "Close K",   100., 10.0, 0., 500000);
   minuit.DefineParameter(5, "Far K",     38000., 100.0, 0., 500000);
   minuit.DefineParameter(6, "Close Co",  100., 10.0, 0., 50000); 
   minuit.DefineParameter(7, "Far Co",    11000., 100.0, 0., 50000);  
   minuit.DefineParameter(8, "Resolution",  6., 1, 3, 10);  
   minuit.DefineParameter(9, "NDBD",      65., 100.0, 0., 1000);     
   minuit.DefineParameter(10, "Lead Bi",    6700., 100.0, 0., 200000);  
   minuit.DefineParameter(11, "TShield Th",  3200, 100.0, 0., 100000);   
   minuit.DefineParameter(12, "50mK Th",  100, 10.0, 0., 100000);   
   minuit.DefineParameter(13, "600mK Th",  20000, 100.0, 0., 100000);   
   minuit.DefineParameter(14, "TShield Ra",  3600, 100.0, 0., 60000);   
   minuit.DefineParameter(15, "50mK Ra",  100, 10.0, 0., 60000);   
   minuit.DefineParameter(16, "600mK Ra",  100, 10.0, 0., 60000);   
   minuit.DefineParameter(17, "OVC Th",  83000, 100.0, 0., 150000);   
   minuit.DefineParameter(18, "OVC Ra",  110000, 100.0, 0., 150000);   
*/




   // Fix parameters for testing
   // minuit.FixParameter(0); // Close Th (or Surface)
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
		  fDataHistoM1->Draw("E");
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
  fModelTotTh->Add(fSmearTShieldThS10,   fParameters[0]);

  // fModelTotTh->Add(fSmearFrameTh,   fParameters[0]);
  // fModelTotTh->Add(fSmearTShieldTh, fParameters[0]);
  // fModelTotTh->Add(fSmear50mKTh,    fParameters[0]);
  // fModelTotTh->Add(fSmear600mKTh,   fParameters[0]);
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




  ///////////////////////////////////////////
  //// Many Parameters
  ///////////////////////////////////////////
  /// Use only after previous step converges!

  // Surface....
  // fModelTotTh->Add(fSmearFrameThS1,   fParameters[11]);
/*
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
*/
  ////////////////////////////////////////////////////////


	
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
  // pt->AddText(Form("Close Th: %0.2E#pm%0.2E --- Far Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
  pt->AddText(Form("Surface Th: %0.2E#pm%0.2E --- Far Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
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

 	TLegend *legfit = new TLegend(0.8,0.8,0.97,0.97);
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
    fDataHistoM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)/yr", dBinSize));
    fDataHistoM1->Draw();
  }
  else if(dMult == 2)
  {
    fDataHistoM2->SetLineColor(1);
    fDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
    fDataHistoM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)/yr", dBinSize));
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
    fModelFrameTh->GetXaxis()->SetTitle("Energy (keV)");
    fModelFrameTh->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));    
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
    fModelFrameRa->GetXaxis()->SetTitle("Energy (keV)");
    fModelFrameRa->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));      
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
    fModelFrameTh->GetXaxis()->SetTitle("Energy (keV)");
    fModelFrameTh->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));      
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
    fModelFrameRa->GetXaxis()->SetTitle("Energy (keV)");
    fModelFrameRa->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));      
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
    fModelFrameK->GetXaxis()->SetTitle("Energy (keV)");
    fModelFrameK->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));  
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
    fModelFrameCo->GetXaxis()->SetTitle("Energy (keV)");
    fModelFrameCo->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));  
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
    fModelNDBD->GetXaxis()->SetTitle("Energy (keV)");
    fModelNDBD->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));  
    fModelNDBD->Draw();
    fModelNDBD->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);
    legndbd->AddEntry(fModelNDBD, "0#nu#beta#beta" ,"l");
    legndbd->Draw();

    TCanvas *cBi = new TCanvas("cBi", "cBi", 1200, 800);
    cBi->SetLogy();
    fModelBi->SetLineColor(1);
    fModelBi->GetXaxis()->SetTitle("Energy (keV)");
    fModelBi->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));      
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
  // Bulk M1
  outTreeFrameThM1 		= LoadMC(dDataDir.c_str(),	"Frame", 	"Th232", "B", 1);
  outTreeTShieldThM1 	= LoadMC(dDataDir.c_str(),	"TShield","Th232", "B", 1);
  outTree50mKThM1 		= LoadMC(dDataDir.c_str(),	"50mK",		"Th232", "B", 1);
  outTree600mKThM1 		= LoadMC(dDataDir.c_str(),	"600mK", 	"Th232", "B", 1);
  outTreeIVCThM1 	  	= LoadMC(dDataDir.c_str(),	"IVC", 		"Th232", "B", 1);
  outTreeOVCThM1 	  	= LoadMC(dDataDir.c_str(),	"OVC", 		"Th232", "B", 1);

  outTreeFrameRaM1	 	= LoadMC(dDataDir.c_str(),	"Frame", 	"Ra226", "B", 1);
  outTreeTShieldRaM1 	= LoadMC(dDataDir.c_str(),	"TShield","Ra226", "B", 1);    
  outTree50mKRaM1	  	= LoadMC(dDataDir.c_str(),	"50mK", 	"Ra226", "B", 1);
  outTree600mKRaM1		= LoadMC(dDataDir.c_str(),	"600mK", 	"Ra226", "B", 1);
  outTreeIVCRaM1	   	= LoadMC(dDataDir.c_str(),	"IVC", 		"Ra226", "B", 1);
  outTreeOVCRaM1 	  	= LoadMC(dDataDir.c_str(),	"OVC", 		"Ra226", "B", 1);

  outTreeFrameKM1 		= LoadMC(dDataDir.c_str(),	"Frame", 	"K40", "B", 1);
  outTreeTShieldKM1 	= LoadMC(dDataDir.c_str(),	"TShield","K40", "B", 1);    
  outTree50mKKM1	   	= LoadMC(dDataDir.c_str(),	"50mK", 	"K40", "B", 1);
  outTree600mKKM1	  	= LoadMC(dDataDir.c_str(),	"600mK", 	"K40", "B", 1);
  outTreeIVCKM1		   	= LoadMC(dDataDir.c_str(),	"IVC", 		"K40", "B", 1);
  outTreeOVCKM1 	   	= LoadMC(dDataDir.c_str(),	"OVC", 		"K40", "B", 1);


  outTreeFrameCoM1 		= LoadMC(dDataDir.c_str(),	"Frame", 	"Co60",	"B", 1);
  outTreeTShieldCoM1 	= LoadMC(dDataDir.c_str(),	"TShield","Co60", "B", 1);    
  outTree50mKCoM1	  	= LoadMC(dDataDir.c_str(),	"50mK", 	"Co60", "B", 1);
  outTree600mKCoM1		= LoadMC(dDataDir.c_str(),	"600mK", 	"Co60", "B", 1);
  outTreeIVCCoM1	   	= LoadMC(dDataDir.c_str(),	"IVC", 		"Co60", "B", 1);
  outTreeOVCCoM1 	  	= LoadMC(dDataDir.c_str(),	"OVC", 		"Co60", "B", 1);

  outTreeNDBDM1 	   	= LoadMC(dDataDir.c_str(),	"Crystal", "0NDBD", "B", 1);
  outTreeBiM1         = LoadMC(dDataDir.c_str(),  "RLead",   "Bi207", "B", 1);

  outTreeFrameThS01M1   = LoadMC(dDataDir.c_str(),  "Frame",  "Th232", "S01", 1);
  outTreeFrameThS1M1    = LoadMC(dDataDir.c_str(),  "Frame",  "Th232", "S1", 1);
  outTreeFrameThS10M1   = LoadMC(dDataDir.c_str(),  "Frame",  "Th232", "S10", 1);
  outTreeFrameThS100M1  = LoadMC(dDataDir.c_str(),  "Frame",  "Th232", "S100", 1);

  outTreeFrameRaS01M1   = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S01", 1);
  outTreeFrameRaS1M1    = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S1", 1);
  outTreeFrameRaS10M1   = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S10", 1);
  outTreeFrameRaS100M1  = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S100", 1);

  outTreeTShieldThS01M1   = LoadMC(dDataDir.c_str(),  "TShield",  "Th232", "S01", 1);
  outTreeTShieldThS1M1    = LoadMC(dDataDir.c_str(),  "TShield",  "Th232", "S1", 1);
  outTreeTShieldThS10M1   = LoadMC(dDataDir.c_str(),  "TShield",  "Th232", "S10", 1);
  outTreeTShieldThS100M1  = LoadMC(dDataDir.c_str(),  "TShield",  "Th232", "S100", 1);


  // Bulk M2
  outTreeFrameThM2    = LoadMC(dDataDir.c_str(),  "Frame",  "Th232", "B", 2);
  outTreeTShieldThM2  = LoadMC(dDataDir.c_str(),  "TShield","Th232", "B", 2);
  outTree50mKThM2     = LoadMC(dDataDir.c_str(),  "50mK",   "Th232", "B", 2);
  outTree600mKThM2    = LoadMC(dDataDir.c_str(),  "600mK",  "Th232", "B", 2);
  outTreeIVCThM2      = LoadMC(dDataDir.c_str(),  "IVC",    "Th232", "B", 2);
  outTreeOVCThM2      = LoadMC(dDataDir.c_str(),  "OVC",    "Th232", "B", 2);

  outTreeFrameRaM2    = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "B", 2);
  outTreeTShieldRaM2  = LoadMC(dDataDir.c_str(),  "TShield","Ra226", "B", 2);    
  outTree50mKRaM2     = LoadMC(dDataDir.c_str(),  "50mK",   "Ra226", "B", 2);
  outTree600mKRaM2    = LoadMC(dDataDir.c_str(),  "600mK",  "Ra226", "B", 2);
  outTreeIVCRaM2      = LoadMC(dDataDir.c_str(),  "IVC",    "Ra226", "B", 2);
  outTreeOVCRaM2      = LoadMC(dDataDir.c_str(),  "OVC",    "Ra226", "B", 2);

  outTreeFrameKM2     = LoadMC(dDataDir.c_str(),  "Frame",  "K40", "B", 2);
  outTreeTShieldKM2   = LoadMC(dDataDir.c_str(),  "TShield","K40", "B", 2);    
  outTree50mKKM2      = LoadMC(dDataDir.c_str(),  "50mK",   "K40", "B", 2);
  outTree600mKKM2     = LoadMC(dDataDir.c_str(),  "600mK",  "K40", "B", 2);
  outTreeIVCKM2       = LoadMC(dDataDir.c_str(),  "IVC",    "K40", "B", 2);
  outTreeOVCKM2       = LoadMC(dDataDir.c_str(),  "OVC",    "K40", "B", 2);

  outTreeFrameCoM2    = LoadMC(dDataDir.c_str(),  "Frame",  "Co60", "B", 2);
  outTreeTShieldCoM2  = LoadMC(dDataDir.c_str(),  "TShield","Co60", "B", 2);    
  outTree50mKCoM2     = LoadMC(dDataDir.c_str(),  "50mK",   "Co60", "B", 2);
  outTree600mKCoM2    = LoadMC(dDataDir.c_str(),  "600mK",  "Co60", "B", 2);
  outTreeIVCCoM2      = LoadMC(dDataDir.c_str(),  "IVC",    "Co60", "B", 2);
  outTreeOVCCoM2      = LoadMC(dDataDir.c_str(),  "OVC",    "Co60", "B", 2);

  outTreeNDBDM2       = LoadMC(dDataDir.c_str(),  "Crystal", "0NDBD", "B", 2);
  outTreeBiM2         = LoadMC(dDataDir.c_str(),  "RLead",   "Bi207", "B", 2);

  outTreeFrameThS01M2   = LoadMC(dDataDir.c_str(),  "Frame",  "Th232", "S01", 2);
  outTreeFrameThS1M2    = LoadMC(dDataDir.c_str(),  "Frame",  "Th232", "S1", 2);
  outTreeFrameThS10M2   = LoadMC(dDataDir.c_str(),  "Frame",  "Th232", "S10", 2);
  outTreeFrameThS100M2  = LoadMC(dDataDir.c_str(),  "Frame",  "Th232", "S100", 2);

  outTreeFrameRaS01M2   = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S01", 2);
  outTreeFrameRaS1M2    = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S1", 2);
  outTreeFrameRaS10M2   = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S10", 2);
  outTreeFrameRaS100M2  = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S100", 2);

  outTreeTShieldThS01M2   = LoadMC(dDataDir.c_str(),  "TShield",  "Th232", "S01", 2);
  outTreeTShieldThS1M2    = LoadMC(dDataDir.c_str(),  "TShield",  "Th232", "S1", 2);
  outTreeTShieldThS10M2   = LoadMC(dDataDir.c_str(),  "TShield",  "Th232", "S10", 2);
  outTreeTShieldThS100M2  = LoadMC(dDataDir.c_str(),  "TShield",  "Th232", "S100", 2);


  // Projecting to histograms
  // M1
  outTreeNDBDM1->Project("fModelNDBDM1",				"Ener1", ener_cut);
  outTreeBiM1->Project("fModelBiM1",            "Ener1", ener_cut);  

	outTreeFrameThM1->Project("fModelFrameThM1", 	   	"Ener1", ener_cut);
	outTreeTShieldThM1->Project("fModelTShieldThM1",	"Ener1", ener_cut);
  outTree50mKThM1->Project("fModel50mKThM1", 		   	"Ener1", ener_cut);
  outTree600mKThM1->Project("fModel600mKThM1", 	  	"Ener1", ener_cut);
  outTreeIVCThM1->Project("fModelIVCThM1", 		    	"Ener1", ener_cut);
  outTreeOVCThM1->Project("fModelOVCThM1", 		     	"Ener1", ener_cut);

	outTreeFrameRaM1->Project("fModelFrameRaM1", 	   	"Ener1", ener_cut);
	outTreeTShieldRaM1->Project("fModelTShieldRaM1",	"Ener1", ener_cut);	
  outTree50mKRaM1->Project("fModel50mKRaM1", 			  "Ener1", ener_cut);
  outTree600mKRaM1->Project("fModel600mKRaM1", 		  "Ener1", ener_cut);
  outTreeIVCRaM1->Project("fModelIVCRaM1", 		     	"Ener1", ener_cut);
  outTreeOVCRaM1->Project("fModelOVCRaM1", 		     	"Ener1", ener_cut);

	outTreeFrameKM1->Project("fModelFrameKM1", 		  	"Ener1", ener_cut);
	outTreeTShieldKM1->Project("fModelTShieldKM1",		"Ener1", ener_cut);	
  outTree50mKKM1->Project("fModel50mKKM1", 		     	"Ener1", ener_cut);
  outTree600mKKM1->Project("fModel600mKKM1", 		   	"Ener1", ener_cut);
  outTreeIVCKM1->Project("fModelIVCKM1", 			    	"Ener1", ener_cut);
  outTreeOVCKM1->Project("fModelOVCKM1", 			    	"Ener1", ener_cut);	

	outTreeFrameCoM1->Project("fModelFrameCoM1", 	   	"Ener1", ener_cut);
	outTreeTShieldCoM1->Project("fModelTShieldCoM1",	"Ener1", ener_cut);	
  outTree50mKCoM1->Project("fModel50mKCoM1", 		   	"Ener1", ener_cut);
  outTree600mKCoM1->Project("fModel600mKCoM1", 	  	"Ener1", ener_cut);
  outTreeIVCCoM1->Project("fModelIVCCoM1", 		     	"Ener1", ener_cut);
  outTreeOVCCoM1->Project("fModelOVCCoM1", 		     	"Ener1", ener_cut);

  outTreeFrameThS01M1->Project("fModelFrameThS01M1",   "Ener1", ener_cut);
  outTreeFrameThS1M1->Project("fModelFrameThS1M1",     "Ener1", ener_cut);
  outTreeFrameThS10M1->Project("fModelFrameThS10M1",   "Ener1", ener_cut);
  outTreeFrameThS100M1->Project("fModelFrameThS100M1", "Ener1", ener_cut);

  outTreeFrameRaS01M1->Project("fModelFrameRaS01M1",   "Ener1", ener_cut);
  outTreeFrameRaS1M1->Project("fModelFrameRaS1M1",     "Ener1", ener_cut);
  outTreeFrameRaS10M1->Project("fModelFrameRaS10M1",   "Ener1", ener_cut);
  outTreeFrameRaS100M1->Project("fModelFrameRaS100M1", "Ener1", ener_cut);

  outTreeTShieldThS01M1->Project("fModelTShieldThS01M1",   "Ener1", ener_cut);
  outTreeTShieldThS1M1->Project("fModelTShieldThS1M1",     "Ener1", ener_cut);
  outTreeTShieldThS10M1->Project("fModelTShieldThS10M1",   "Ener1", ener_cut);
  outTreeTShieldThS100M1->Project("fModelTShieldThS100M1", "Ener1", ener_cut);
  
  // M1
  outTreeNDBDM2->Project("fModelNDBDM2",        "Ener1", ener_cut);
  outTreeBiM2->Project("fModelBiM2",            "Ener1", ener_cut);  

  outTreeFrameThM2->Project("fModelFrameThM2",      "Ener1", ener_cut);
  outTreeTShieldThM2->Project("fModelTShieldThM2",  "Ener1", ener_cut);
  outTree50mKThM2->Project("fModel50mKThM2",        "Ener1", ener_cut);
  outTree600mKThM2->Project("fModel600mKThM2",      "Ener1", ener_cut);
  outTreeIVCThM2->Project("fModelIVCThM2",          "Ener1", ener_cut);
  outTreeOVCThM2->Project("fModelOVCThM2",          "Ener1", ener_cut);

  outTreeFrameRaM2->Project("fModelFrameRaM2",      "Ener1", ener_cut);
  outTreeTShieldRaM2->Project("fModelTShieldRaM2",  "Ener1", ener_cut); 
  outTree50mKRaM2->Project("fModel50mKRaM2",        "Ener1", ener_cut);
  outTree600mKRaM2->Project("fModel600mKRaM2",      "Ener1", ener_cut);
  outTreeIVCRaM2->Project("fModelIVCRaM2",          "Ener1", ener_cut);
  outTreeOVCRaM2->Project("fModelOVCRaM2",          "Ener1", ener_cut);

  outTreeFrameKM2->Project("fModelFrameKM2",        "Ener1", ener_cut);
  outTreeTShieldKM2->Project("fModelTShieldKM2",    "Ener1", ener_cut); 
  outTree50mKKM2->Project("fModel50mKKM2",          "Ener1", ener_cut);
  outTree600mKKM2->Project("fModel600mKKM2",        "Ener1", ener_cut);
  outTreeIVCKM2->Project("fModelIVCKM2",            "Ener1", ener_cut);
  outTreeOVCKM2->Project("fModelOVCKM2",            "Ener1", ener_cut); 

  outTreeFrameCoM2->Project("fModelFrameCoM2",      "Ener1", ener_cut);
  outTreeTShieldCoM2->Project("fModelTShieldCoM2",  "Ener1", ener_cut); 
  outTree50mKCoM2->Project("fModel50mKCoM2",        "Ener1", ener_cut);
  outTree600mKCoM2->Project("fModel600mKCoM2",      "Ener1", ener_cut);
  outTreeIVCCoM2->Project("fModelIVCCoM2",          "Ener1", ener_cut);
  outTreeOVCCoM2->Project("fModelOVCCoM2",          "Ener1", ener_cut);

  outTreeFrameThS01M2->Project("fModelFrameThS01M2",   "Ener1", ener_cut);
  outTreeFrameThS1M2->Project("fModelFrameThS1M2",     "Ener1", ener_cut);
  outTreeFrameThS10M2->Project("fModelFrameThS10M2",   "Ener1", ener_cut);
  outTreeFrameThS100M2->Project("fModelFrameThS100M2", "Ener1", ener_cut);

  outTreeFrameRaS01M2->Project("fModelFrameRaS01M2",   "Ener1", ener_cut);
  outTreeFrameRaS1M2->Project("fModelFrameRaS1M2",     "Ener1", ener_cut);
  outTreeFrameRaS10M2->Project("fModelFrameRaS10M2",   "Ener1", ener_cut);
  outTreeFrameRaS100M2->Project("fModelFrameRaS100M2", "Ener1", ener_cut);

  outTreeTShieldThS01M2->Project("fModelTShieldThS01M2",   "Ener1", ener_cut);
  outTreeTShieldThS1M2->Project("fModelTShieldThS1M2",     "Ener1", ener_cut);
  outTreeTShieldThS10M2->Project("fModelTShieldThS10M2",   "Ener1", ener_cut);
  outTreeTShieldThS100M2->Project("fModelTShieldThS100M2", "Ener1", ener_cut);


	cout << "Loaded MC" << endl;


	// Normalize all MC histograms
  // Fixing normalization of NDBD from 2000 to 2650
  // M1
	NormalizePDF(fModelFrameThM1, outTreeFrameThM1, 	50, 2700);
	NormalizePDF(fModelTShieldThM1, outTreeTShieldThM1,	50, 2700);
	NormalizePDF(fModel50mKThM1, outTree50mKThM1,		50, 2700);
	NormalizePDF(fModel600mKThM1, outTree600mKThM1,		50, 2700);
	NormalizePDF(fModelIVCThM1, outTreeIVCThM1,			50, 2700);
	NormalizePDF(fModelOVCThM1, outTreeOVCThM1,			50, 2700);

	NormalizePDF(fModelFrameRaM1,  outTreeFrameRaM1,	50, 2700);
	NormalizePDF(fModelTShieldRaM1, outTreeTShieldRaM1,	50, 2700);	
	NormalizePDF(fModel50mKRaM1, outTree50mKRaM1,		50, 2700);
	NormalizePDF(fModel600mKRaM1, outTree600mKRaM1,		50, 2700);
	NormalizePDF(fModelIVCRaM1, outTreeIVCRaM1,			50, 2700);
	NormalizePDF(fModelOVCRaM1, outTreeOVCRaM1,			50, 2700);

	NormalizePDF(fModelFrameKM1, 	outTreeFrameKM1,		50, 2700);
	NormalizePDF(fModelTShieldKM1, outTreeTShieldKM1,	50, 2700);	
	NormalizePDF(fModel50mKKM1, outTree50mKKM1,			50, 2700);
	NormalizePDF(fModel600mKKM1, outTree600mKKM1,		50, 2700);
	NormalizePDF(fModelIVCKM1, outTreeIVCKM1,			50, 2700);
	NormalizePDF(fModelOVCKM1, outTreeOVCKM1,			50, 2700);

	NormalizePDF(fModelFrameCoM1, outTreeFrameCoM1,		50, 2700);
	NormalizePDF(fModelTShieldCoM1, outTreeTShieldCoM1,	50, 2700);	
	NormalizePDF(fModel50mKCoM1, outTree50mKCoM1,		50, 2700);
	NormalizePDF(fModel600mKCoM1, outTree600mKCoM1,		50, 2700);
	NormalizePDF(fModelIVCCoM1, outTreeIVCCoM1,			50, 2700);
	NormalizePDF(fModelOVCCoM1, outTreeOVCCoM1,			50, 2700);

  NormalizePDF(fModelNDBDM1, outTreeNDBDM1,     50, 2700);
  NormalizePDF(fModelBiM1, outTreeBiM1,         50, 2700);


  NormalizePDF(fModelFrameThS01M1,   outTreeFrameThS01M1, 50, 2700);
  NormalizePDF(fModelFrameThS1M1,    outTreeFrameThS1M1, 50, 2700);
  NormalizePDF(fModelFrameThS10M1,   outTreeFrameThS10M1, 50, 2700);
  NormalizePDF(fModelFrameThS100M1,  outTreeFrameThS100M1, 50, 2700);

  NormalizePDF(fModelFrameRaS01M1,   outTreeFrameRaS01M1, 50, 2700);
  NormalizePDF(fModelFrameRaS1M1,    outTreeFrameRaS1M1, 50, 2700);
  NormalizePDF(fModelFrameRaS10M1,   outTreeFrameRaS10M1, 50, 2700);
  NormalizePDF(fModelFrameRaS100M1,  outTreeFrameRaS100M1, 50, 2700);

  NormalizePDF(fModelTShieldThS01M1,   outTreeTShieldThS01M1, 50, 2700);
  NormalizePDF(fModelTShieldThS1M1,    outTreeTShieldThS1M1, 50, 2700);
  NormalizePDF(fModelTShieldThS10M1,   outTreeTShieldThS10M1, 50, 2700);
  NormalizePDF(fModelTShieldThS100M1,  outTreeTShieldThS100M1, 50, 2700);

  // M2
  NormalizePDF(fModelFrameThM2, outTreeFrameThM2,   50, 2700);
  NormalizePDF(fModelTShieldThM2, outTreeTShieldThM2, 50, 2700);
  NormalizePDF(fModel50mKThM2, outTree50mKThM2,   50, 2700);
  NormalizePDF(fModel600mKThM2, outTree600mKThM2,   50, 2700);
  NormalizePDF(fModelIVCThM2, outTreeIVCThM2,     50, 2700);
  NormalizePDF(fModelOVCThM2, outTreeOVCThM2,     50, 2700);

  NormalizePDF(fModelFrameRaM2,  outTreeFrameRaM2,  50, 2700);
  NormalizePDF(fModelTShieldRaM2, outTreeTShieldRaM2, 50, 2700);  
  NormalizePDF(fModel50mKRaM2, outTree50mKRaM2,   50, 2700);
  NormalizePDF(fModel600mKRaM2, outTree600mKRaM2,   50, 2700);
  NormalizePDF(fModelIVCRaM2, outTreeIVCRaM2,     50, 2700);
  NormalizePDF(fModelOVCRaM2, outTreeOVCRaM2,     50, 2700);

  NormalizePDF(fModelFrameKM2,  outTreeFrameKM2,    50, 2700);
  NormalizePDF(fModelTShieldKM2, outTreeTShieldKM2, 50, 2700);  
  NormalizePDF(fModel50mKKM2, outTree50mKKM2,     50, 2700);
  NormalizePDF(fModel600mKKM2, outTree600mKKM2,   50, 2700);
  NormalizePDF(fModelIVCKM2, outTreeIVCKM2,     50, 2700);
  NormalizePDF(fModelOVCKM2, outTreeOVCKM2,     50, 2700);

  NormalizePDF(fModelFrameCoM2, outTreeFrameCoM2,   50, 2700);
  NormalizePDF(fModelTShieldCoM2, outTreeTShieldCoM2, 50, 2700);  
  NormalizePDF(fModel50mKCoM2, outTree50mKCoM2,   50, 2700);
  NormalizePDF(fModel600mKCoM2, outTree600mKCoM2,   50, 2700);
  NormalizePDF(fModelIVCCoM2, outTreeIVCCoM2,     50, 2700);
  NormalizePDF(fModelOVCCoM2, outTreeOVCCoM2,     50, 2700);

  NormalizePDF(fModelNDBDM2, outTreeNDBDM2,     50, 2700);
  NormalizePDF(fModelBiM2, outTreeBiM2,         50, 2700);


  NormalizePDF(fModelFrameThS01M2,   outTreeFrameThS01M1, 50, 2700);
  NormalizePDF(fModelFrameThS1M2,    outTreeFrameThS1M1, 50, 2700);
  NormalizePDF(fModelFrameThS10M2,   outTreeFrameThS10M1, 50, 2700);
  NormalizePDF(fModelFrameThS100M2,  outTreeFrameThS100M1, 50, 2700);

  NormalizePDF(fModelFrameRaS01M2,   outTreeFrameRaS01M1, 50, 2700);
  NormalizePDF(fModelFrameRaS1M2,    outTreeFrameRaS1M1, 50, 2700);
  NormalizePDF(fModelFrameRaS10M2,   outTreeFrameRaS10M1, 50, 2700);
  NormalizePDF(fModelFrameRaS100M2,  outTreeFrameRaS100M1, 50, 2700);

  NormalizePDF(fModelTShieldThS01M2,   outTreeTShieldThS01M2, 50, 2700);
  NormalizePDF(fModelTShieldThS1M2,    outTreeTShieldThS1M2, 50, 2700);
  NormalizePDF(fModelTShieldThS10M2,   outTreeTShieldThS10M2, 50, 2700);
  NormalizePDF(fModelTShieldThS100M2,  outTreeTShieldThS100M2, 50, 2700);




	cout << "Normalized MC PDFs" << endl;

  if(bFixedRes)
  {
    // cout << "Fixed resolution: " << endl;
    cout << "Smearing histograms with constant 6 keV resolution" << endl;

    // double dRes = 4;
    double dRes = 6.0/2.355;

    // Adding the 10 micron distribution for now...
    SmearMC(fModelFrameThS10M1, fSmearFrameThS10M1, dRes);
    SmearMC(fModelTShieldThS10M1, fSmearTShieldThS10M1, dRes);

    // M1
    SmearMC(fModelFrameThM1, fSmearFrameThM1, dRes);
    SmearMC(fModelTShieldThM1, fSmearTShieldThM1, dRes);  
    SmearMC(fModel50mKThM1, fSmear50mKThM1, dRes);
    SmearMC(fModel600mKThM1, fSmear600mKThM1, dRes);
    SmearMC(fModelIVCThM1, fSmearIVCThM1, dRes);
    SmearMC(fModelOVCThM1, fSmearOVCThM1, dRes);

    SmearMC(fModelFrameRaM1, fSmearFrameRaM1, dRes);
    SmearMC(fModelTShieldRaM1, fSmearTShieldRaM1, dRes);  
    SmearMC(fModel50mKRaM1, fSmear50mKRaM1, dRes);
    SmearMC(fModel600mKRaM1, fSmear600mKRaM1, dRes);
    SmearMC(fModelIVCRaM1, fSmearIVCRaM1, dRes);
    SmearMC(fModelOVCRaM1, fSmearOVCRaM1, dRes);

    SmearMC(fModelFrameKM1, fSmearFrameKM1, dRes);
    SmearMC(fModelTShieldKM1, fSmearTShieldKM1, dRes);
    SmearMC(fModel50mKKM1, fSmear50mKKM1, dRes);
    SmearMC(fModel600mKKM1, fSmear600mKKM1, dRes);
    SmearMC(fModelIVCKM1, fSmearIVCKM1, dRes);
    SmearMC(fModelOVCKM1, fSmearOVCKM1, dRes); 

    SmearMC(fModelFrameCoM1, fSmearFrameCoM1, dRes);
    SmearMC(fModelTShieldCoM1, fSmearTShieldCoM1, dRes);
    SmearMC(fModel50mKCoM1, fSmear50mKCoM1, dRes);
    SmearMC(fModel600mKCoM1, fSmear600mKCoM1, dRes);
    SmearMC(fModelIVCCoM1, fSmearIVCCoM1, dRes);
    SmearMC(fModelOVCCoM1, fSmearOVCCoM1, dRes);  

    SmearMC(fModelNDBDM1, fSmearNDBDM1, dRes);  
    SmearMC(fModelBiM1, fSmearBiM1, dRes);  

    // M2
    SmearMC(fModelFrameThM2, fSmearFrameThM2, dRes);
    SmearMC(fModelTShieldThM2, fSmearTShieldThM2, dRes);  
    SmearMC(fModel50mKThM2, fSmear50mKThM2, dRes);
    SmearMC(fModel600mKThM2, fSmear600mKThM2, dRes);
    SmearMC(fModelIVCThM2, fSmearIVCThM2, dRes);
    SmearMC(fModelOVCThM2, fSmearOVCThM2, dRes);

    SmearMC(fModelFrameRaM2, fSmearFrameRaM2, dRes);
    SmearMC(fModelTShieldRaM2, fSmearTShieldRaM2, dRes);  
    SmearMC(fModel50mKRaM2, fSmear50mKRaM2, dRes);
    SmearMC(fModel600mKRaM2, fSmear600mKRaM2, dRes);
    SmearMC(fModelIVCRaM2, fSmearIVCRaM2, dRes);
    SmearMC(fModelOVCRaM2, fSmearOVCRaM2, dRes);

    SmearMC(fModelFrameKM2, fSmearFrameKM2, dRes);
    SmearMC(fModelTShieldKM2, fSmearTShieldKM2, dRes);
    SmearMC(fModel50mKKM2, fSmear50mKKM2, dRes);
    SmearMC(fModel600mKKM2, fSmear600mKKM2, dRes);
    SmearMC(fModelIVCKM2, fSmearIVCKM2, dRes);
    SmearMC(fModelOVCKM2, fSmearOVCKM2, dRes); 

    SmearMC(fModelFrameCoM2, fSmearFrameCoM2, dRes);
    SmearMC(fModelTShieldCoM2, fSmearTShieldCoM2, dRes);
    SmearMC(fModel50mKCoM2, fSmear50mKCoM2, dRes);
    SmearMC(fModel600mKCoM2, fSmear600mKCoM2, dRes);
    SmearMC(fModelIVCCoM2, fSmearIVCCoM2, dRes);
    SmearMC(fModelOVCCoM2, fSmearOVCCoM2, dRes);  

    SmearMC(fModelNDBDM2, fSmearNDBDM2, dRes);  
    SmearMC(fModelBiM2, fSmearBiM2, dRes);  

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
  // cout<< "Par12 = "  << fParameters[12] << " +/- " << fParError[12] << endl;
  // cout<< "Par13 = "  << fParameters[13] << " +/- " << fParError[13] << endl;
  // cout<< "Par14 = "  << fParameters[14] << " +/- " << fParError[14] << endl;
  // cout<< "Par15 = "  << fParameters[15] << " +/- " << fParError[15] << endl;
  // cout<< "Par16 = "  << fParameters[16] << " +/- " << fParError[16] << endl;
  // cout<< "Par17 = "  << fParameters[17] << " +/- " << fParError[17] << endl;
  // cout<< "Par18 = "  << fParameters[18] << " +/- " << fParError[18] << endl;


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
	fModelTotM1->Reset();
  fModelTotM2->Reset();

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
  // fModelTot->Add( fSmearTShieldThS10,    fParameters[0]);
  // M1
  fModelTotM1->Add( fSmearFrameThM1,    fParameters[0]);
  fModelTotM1->Add( fSmearTShieldThM1,  fParameters[0]);  
  fModelTotM1->Add( fSmear50mKThM1,     fParameters[0]);
  fModelTotM1->Add( fSmear600mKThM1,    fParameters[0]);
  fModelTotM1->Add( fSmearIVCThM1,      fParameters[1]);
  fModelTotM1->Add( fSmearOVCThM1,      fParameters[1]);

  fModelTotM1->Add( fSmearFrameRaM1,    fParameters[2]);
  fModelTotM1->Add( fSmearTShieldRaM1,  fParameters[2]);  
  fModelTotM1->Add( fSmear50mKRaM1,     fParameters[2]);
  fModelTotM1->Add( fSmear600mKRaM1,    fParameters[2]);
  fModelTotM1->Add( fSmearIVCRaM1,      fParameters[3]);
  fModelTotM1->Add( fSmearOVCRaM1,      fParameters[3]);

  fModelTotM1->Add( fSmearFrameKM1,    fParameters[4]);
  fModelTotM1->Add( fSmearTShieldKM1,  fParameters[4]);
  fModelTotM1->Add( fSmear50mKKM1,     fParameters[4]);
  fModelTotM1->Add( fSmear600mKKM1,    fParameters[4]);
  fModelTotM1->Add( fSmearIVCKM1,      fParameters[5]);
  fModelTotM1->Add( fSmearOVCKM1,      fParameters[5]); 

  fModelTotM1->Add( fSmearFrameCoM1,    fParameters[6]);
  fModelTotM1->Add( fSmearTShieldCoM1,  fParameters[6]);
  fModelTotM1->Add( fSmear50mKCoM1,     fParameters[6]);
  fModelTotM1->Add( fSmear600mKCoM1,    fParameters[6]);
  fModelTotM1->Add( fSmearIVCCoM1,      fParameters[7]);
  fModelTotM1->Add( fSmearOVCCoM1,      fParameters[7]);  

  fModelTotM1->Add( fSmearNDBDM1,      fParameters[9]);  
  fModelTotM1->Add( fSmearBiM1,        fParameters[10]);  



  // M2
  fModelTotM2->Add( fSmearFrameThM2,    fParameters[0]);
  fModelTotM2->Add( fSmearTShieldThM2,  fParameters[0]);  
  fModelTotM2->Add( fSmear50mKThM2,     fParameters[0]);
  fModelTotM2->Add( fSmear600mKThM2,    fParameters[0]);
  fModelTotM2->Add( fSmearIVCThM2,      fParameters[1]);
  fModelTotM2->Add( fSmearOVCThM2,      fParameters[1]);

  fModelTotM2->Add( fSmearFrameRaM2,    fParameters[2]);
  fModelTotM2->Add( fSmearTShieldRaM2,  fParameters[2]);  
  fModelTotM2->Add( fSmear50mKRaM2,     fParameters[2]);
  fModelTotM2->Add( fSmear600mKRaM2,    fParameters[2]);
  fModelTotM2->Add( fSmearIVCRaM2,      fParameters[3]);
  fModelTotM2->Add( fSmearOVCRaM2,      fParameters[3]);

  fModelTotM2->Add( fSmearFrameKM2,    fParameters[4]);
  fModelTotM2->Add( fSmearTShieldKM2,  fParameters[4]);
  fModelTotM2->Add( fSmear50mKKM2,     fParameters[4]);
  fModelTotM2->Add( fSmear600mKKM2,    fParameters[4]);
  fModelTotM2->Add( fSmearIVCKM2,      fParameters[5]);
  fModelTotM2->Add( fSmearOVCKM2,      fParameters[5]); 

  fModelTotM2->Add( fSmearFrameCoM2,    fParameters[6]);
  fModelTotM2->Add( fSmearTShieldCoM2,  fParameters[6]);
  fModelTotM2->Add( fSmear50mKCoM2,     fParameters[6]);
  fModelTotM2->Add( fSmear600mKCoM2,    fParameters[6]);
  fModelTotM2->Add( fSmearIVCCoM2,      fParameters[7]);
  fModelTotM2->Add( fSmearOVCCoM2,      fParameters[7]);  

  fModelTotM2->Add( fSmearNDBDM2,      fParameters[9]);  
  fModelTotM2->Add( fSmearBiM2,        fParameters[10]);  


  /////////////////////////////////////
  //// Many parameters
  ////////////////////////////////////
/*
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

  fModelTot->Add( fSmearFrameK,     fParameters[4]);
  fModelTot->Add( fSmearTShieldK,   fParameters[4]);
  fModelTot->Add( fSmear50mKK,      fParameters[4]);
  fModelTot->Add( fSmear600mKK,     fParameters[4]);
  fModelTot->Add( fSmearIVCK,       fParameters[5]);
  fModelTot->Add( fSmearOVCK,       fParameters[5]); 

  fModelTot->Add( fSmearFrameCo,    fParameters[6]);
  fModelTot->Add( fSmearTShieldCo,  fParameters[6]);
  fModelTot->Add( fSmear50mKCo,     fParameters[6]);
  fModelTot->Add( fSmear600mKCo,    fParameters[6]);
  fModelTot->Add( fSmearIVCCo,      fParameters[7]);
  fModelTot->Add( fSmearOVCCo,      fParameters[7]);  

  fModelTot->Add( fSmearNDBD,       fParameters[9]);  

  fModelTot->Add( fSmearBi,         fParameters[10]);  
*/

  }





/*
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
*/


}




