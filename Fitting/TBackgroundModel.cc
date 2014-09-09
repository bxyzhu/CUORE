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


TBackgroundModel::TBackgroundModel(double fFitMin, double fFitMax)
{

  dNumCalls = 0;
  dSecToYears = 1./(60*60*24*365);

  dDataDir =  "/Users/brian/macros/Simulations/Bkg/";
  dDataIntegral = 0;
  bToyFit = false;

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
  fModelTotBiM2    = new TH1D("fModelTotBiM2",    "Total Bi207",   dNBins, dMinEnergy, dMaxEnergy);



  // Modeling
  gaus = new TF1("gaus","gaus(0)", dMinEnergy, dMaxEnergy);


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

}
  
// Needs to be updated
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


   ////////////////////////////////////////////////
   // Using less parameters
   ////////////////////////////////////////////////
   minuit.DefineParameter(0, "Close Th",  3000., 10.0, 0., 100000);
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


   // Fix parameters here
   // minuit.FixParameter(0); // Close Th (or Surface)
   // minuit.FixParameter(1); // Far Th
   // minuit.FixParameter(2); // Close Ra
   // minuit.FixParameter(3); // Far Ra
   minuit.FixParameter(4); // Close K
   minuit.FixParameter(5); // Far K
   // minuit.FixParameter(6); // Close Co
   // minuit.FixParameter(7); // Far Co
    minuit.FixParameter(8); // Resolution
   // minuit.FixParameter(9); // NDBD
   minuit.FixParameter(10); // Bi207

  // Number of Parameters! (for Chi-squared/NDF calculation)
  int dNumParameters = 7;




   //Tell minuit what external function to use 
   minuit.SetFCN(myExternal_FCN);
   
   int status = minuit.Migrad(); // this actually does the minimisation


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
	
	cout << "At the end; ChiSq/NDF = " << GetChiSquare()/(2*(dFitMax-dFitMin)/dBinSize - dNumParameters) << endl;
  cout << "Total number of calls = " << dNumCalls << endl;


  ///////////////////////////////////////////
  //// Few Parameters
  ///////////////////////////////////////////
  /// Add Modeled Histograms after chi-squared minimization

  // Surface....
  // fModelTotThM1->Add(fModelTShieldThS10,   fParameters[0]);

  // M1 Parameters
  fModelTotThM1->Add(fModelFrameThM1,   fParameters[0]);
  fModelTotThM1->Add(fModelTShieldThM1, fParameters[0]);
  fModelTotThM1->Add(fModel50mKThM1,    fParameters[0]);
  fModelTotThM1->Add(fModel600mKThM1,   fParameters[0]);
  fModelTotThM1->Add(fModelIVCThM1,     fParameters[1]);
  fModelTotThM1->Add(fModelOVCThM1,     fParameters[1]);

  fModelTotRaM1->Add(fModelFrameRaM1,   fParameters[2]);
  fModelTotRaM1->Add(fModelTShieldRaM1, fParameters[2]);
  fModelTotRaM1->Add(fModel50mKRaM1,    fParameters[2]);
  fModelTotRaM1->Add(fModel600mKRaM1,   fParameters[2]);
  fModelTotRaM1->Add(fModelIVCRaM1,     fParameters[3]);
  fModelTotRaM1->Add(fModelOVCRaM1,     fParameters[3]);

  fModelTotKM1->Add(fModelFrameKM1,     fParameters[4]);
  fModelTotKM1->Add(fModelTShieldKM1,   fParameters[4]);
  fModelTotKM1->Add(fModel50mKKM1,      fParameters[4]);
  fModelTotKM1->Add(fModel600mKKM1,     fParameters[4]);
  fModelTotKM1->Add(fModelIVCKM1,       fParameters[5]);
  fModelTotKM1->Add(fModelOVCKM1,       fParameters[5]);

  fModelTotCoM1->Add(fModelFrameCoM1,   fParameters[6]);
  fModelTotCoM1->Add(fModelTShieldCoM1, fParameters[6]);
  fModelTotCoM1->Add(fModel50mKCoM1,    fParameters[6]);
  fModelTotCoM1->Add(fModel600mKCoM1,   fParameters[6]);
  fModelTotCoM1->Add(fModelIVCCoM1,     fParameters[7]);
  fModelTotCoM1->Add(fModelOVCCoM1,     fParameters[7]);

  fModelTotNDBDM1->Add(fModelNDBDM1,    fParameters[9]);
  fModelTotBiM1->Add(fModelBiM1,        fParameters[10]);

  // M2 Parameters
  fModelTotThM2->Add(fModelFrameThM2,   fParameters[0]);
  fModelTotThM2->Add(fModelTShieldThM2, fParameters[0]);
  fModelTotThM2->Add(fModel50mKThM2,    fParameters[0]);
  fModelTotThM2->Add(fModel600mKThM2,   fParameters[0]);
  fModelTotThM2->Add(fModelIVCThM2,     fParameters[1]);
  fModelTotThM2->Add(fModelOVCThM2,     fParameters[1]);

  fModelTotRaM2->Add(fModelFrameRaM2,   fParameters[2]);
  fModelTotRaM2->Add(fModelTShieldRaM2, fParameters[2]);
  fModelTotRaM2->Add(fModel50mKRaM2,    fParameters[2]);
  fModelTotRaM2->Add(fModel600mKRaM2,   fParameters[2]);
  fModelTotRaM2->Add(fModelIVCRaM2,     fParameters[3]);
  fModelTotRaM2->Add(fModelOVCRaM2,     fParameters[3]);

  fModelTotKM2->Add(fModelFrameKM2,     fParameters[4]);
  fModelTotKM2->Add(fModelTShieldKM2,   fParameters[4]);
  fModelTotKM2->Add(fModel50mKKM2,      fParameters[4]);
  fModelTotKM2->Add(fModel600mKKM2,     fParameters[4]);
  fModelTotKM2->Add(fModelIVCKM2,       fParameters[5]);
  fModelTotKM2->Add(fModelOVCKM2,       fParameters[5]);

  fModelTotCoM2->Add(fModelFrameCoM2,   fParameters[6]);
  fModelTotCoM2->Add(fModelTShieldCoM2, fParameters[6]);
  fModelTotCoM2->Add(fModel50mKCoM2,    fParameters[6]);
  fModelTotCoM2->Add(fModel600mKCoM2,   fParameters[6]);
  fModelTotCoM2->Add(fModelIVCCoM2,     fParameters[7]);
  fModelTotCoM2->Add(fModelOVCCoM2,     fParameters[7]);

  fModelTotNDBDM2->Add(fModelNDBDM2,    fParameters[9]);
  fModelTotBiM2->Add(fModelBiM2,      fParameters[10]);



  ///////////////////////////////////////////
  //// Many Parameters
  ///////////////////////////////////////////
  /// Use only after previous step converges!

  // Surface....
  // fModelTotTh->Add(fModelFrameThS1,   fParameters[11]);
/*
  // M1
  fModelTotThM1->Add(fModelFrameThM1,   fParameters[0]);
  fModelTotThM1->Add(fModelTShieldThM1, fParameters[11]);
  fModelTotThM1->Add(fModel50mKThM1,    fParameters[12]);
  fModelTotThM1->Add(fModel600mKThM1,   fParameters[13]);
  fModelTotThM1->Add(fModelIVCThM1,     fParameters[1]);
  fModelTotThM1->Add(fModelOVCThM1,     fParameters[17]);

  fModelTotRaM1->Add(fModelFrameRaM1,   fParameters[2]);
  fModelTotRaM1->Add(fModelTShieldRaM1, fParameters[14]);
  fModelTotRaM1->Add(fModel50mKRaM1,    fParameters[15]);
  fModelTotRaM1->Add(fModel600mKRaM1,   fParameters[16]);
  fModelTotRaM1->Add(fModelIVCRaM1,     fParameters[3]);
  fModelTotRaM1->Add(fModelOVCRaM1,     fParameters[18]);

  fModelTotKM1->Add(fModelFrameKM1,     fParameters[4]);
  fModelTotKM1->Add(fModelTShieldKM1,   fParameters[4]);
  fModelTotKM1->Add(fModel50mKKM1,      fParameters[4]);
  fModelTotKM1->Add(fModel600mKKM1,     fParameters[4]);
  fModelTotKM1->Add(fModelIVCKM1,       fParameters[5]);
  fModelTotKM1->Add(fModelOVCKM1,       fParameters[5]);

  fModelTotCoM1->Add(fModelFrameCoM1,   fParameters[6]);
  fModelTotCoM1->Add(fModelTShieldCoM1, fParameters[6]);
  fModelTotCoM1->Add(fModel50mKCoM1,    fParameters[6]);
  fModelTotCoM1->Add(fModel600mKCoM1,   fParameters[6]);
  fModelTotCoM1->Add(fModelIVCCoM1,     fParameters[7]);
  fModelTotCoM1->Add(fModelOVCCoM1,     fParameters[7]);

  fModelTotNDBDM1->Add(fModelNDBDM1,    fParameters[9]);
  fModelTotBiM1->Add(fModelBiM1,      fParameters[10]);

  // M2
  fModelTotThM2->Add(fModelFrameThM2,   fParameters[0]);
  fModelTotThM2->Add(fModelTShieldThM2, fParameters[11]);
  fModelTotThM2->Add(fModel50mKThM2,    fParameters[12]);
  fModelTotThM2->Add(fModel600mKThM2,   fParameters[13]);
  fModelTotThM2->Add(fModelIVCThM2,     fParameters[1]);
  fModelTotThM2->Add(fModelOVCThM2,     fParameters[17]);

  fModelTotRaM2->Add(fModelFrameRaM2,   fParameters[2]);
  fModelTotRaM2->Add(fModelTShieldRaM2, fParameters[14]);
  fModelTotRaM2->Add(fModel50mKRaM2,    fParameters[15]);
  fModelTotRaM2->Add(fModel600mKRaM2,   fParameters[16]);
  fModelTotRaM2->Add(fModelIVCRaM2,     fParameters[3]);
  fModelTotRaM2->Add(fModelOVCRaM2,     fParameters[18]);

  fModelTotKM2->Add(fModelFrameKM2,     fParameters[4]);
  fModelTotKM2->Add(fModelTShieldKM2,   fParameters[4]);
  fModelTotKM2->Add(fModel50mKKM2,      fParameters[4]);
  fModelTotKM2->Add(fModel600mKKM2,     fParameters[4]);
  fModelTotKM2->Add(fModelIVCKM2,       fParameters[5]);
  fModelTotKM2->Add(fModelOVCKM2,       fParameters[5]);

  fModelTotCoM2->Add(fModelFrameCoM2,   fParameters[6]);
  fModelTotCoM2->Add(fModelTShieldCoM2, fParameters[6]);
  fModelTotCoM2->Add(fModel50mKCoM2,    fParameters[6]);
  fModelTotCoM2->Add(fModel600mKCoM2,   fParameters[6]);
  fModelTotCoM2->Add(fModelIVCCoM2,     fParameters[7]);
  fModelTotCoM2->Add(fModelOVCCoM2,     fParameters[7]);

  fModelTotNDBDM2->Add(fModelNDBDM2,    fParameters[9]);
  fModelTotBiM2->Add(fModelBiM2,      fParameters[10]);

*/
  ////////////////////////////////////////////////////////



  // Toy fit plots currently aren't complete...
  if(bToyFit)
  {

    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->SetLogy();    
  ////// Draw Toy Data
    fToyData->SetLineColor(1);
    fToyData->SetLineWidth(2);
    fToyData->GetXaxis()->SetTitle("Energy (keV)");
    fToyData->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));
    fToyData->Draw();

  }
  else
  {

    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->SetLogy();

    ///// Draw Data M1
    fDataHistoM1->SetLineColor(1);
    fDataHistoM1->SetLineWidth(2);
    fDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
    fDataHistoM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)/yr", dBinSize));
    fDataHistoM1->SetMaximum(90000);
    fDataHistoM1->GetXaxis()->SetRange(1, 2650/dBinSize+5);
    fDataHistoM1->Draw("E");

  
    fModelTotM1->SetLineColor(2);
    fModelTotM1->Draw("SAME");

    fModelTotThM1->SetLineColor(3);
    fModelTotThM1->SetLineStyle(2);
    fModelTotRaM1->SetLineColor(4);
    fModelTotRaM1->SetLineStyle(2);
    fModelTotKM1->SetLineColor(6);
    fModelTotKM1->SetLineStyle(2);
    fModelTotCoM1->SetLineColor(7);
    fModelTotCoM1->SetLineStyle(2);
    fModelTotNDBDM1->SetLineColor(42);
    fModelTotNDBDM1->SetLineStyle(2);
    fModelTotBiM1->SetLineColor(5);
    fModelTotBiM1->SetLineStyle(2);

    fModelTotThM1->Draw("SAME");
    fModelTotRaM1->Draw("SAME");
    fModelTotKM1->Draw("SAME");
    fModelTotCoM1->Draw("SAME");
    fModelTotNDBDM1->Draw("SAME");
    fModelTotBiM1->Draw("SAME");

    // Few Parameters
    TPaveText *pt1 = new TPaveText(0.35,0.78,0.70,0.98,"NB NDC");
    pt1->AddText(Form("Fit Range (M1): %.0f to %.0f keV -- #chi^{2}/NDF: %0.3f", dFitMin, dFitMax, (GetChiSquare()/((dFitMax-dFitMin)/dBinSize - dNumParameters)) ));
    pt1->AddText(Form("Close Th: %0.2E#pm%0.2E --- Far Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
    // pt1->AddText(Form("Surface Th: %0.2E#pm%0.2E --- Far Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
    pt1->AddText(Form("Close Ra: %0.2E#pm%0.2E --- Far Ra: %0.2E#pm%0.2E", fParameters[2], fParError[2], fParameters[3], fParError[3] ));
    pt1->AddText(Form("Close K: %0.2E#pm%0.2E --- Far K: %0.2E#pm%0.2E", fParameters[4], fParError[4], fParameters[5], fParError[5] ));
    pt1->AddText(Form("Close Co: %0.2E#pm%0.2E --- Far Co: %0.2E#pm%0.2E", fParameters[6], fParError[6], fParameters[7], fParError[7] ));
    pt1->AddText(Form("Bi-207: %0.2E#pm%0.2E --- NDBD: %0.2E#pm%0.2E", fParameters[10], fParError[10], fParameters[9], fParError[9] ));


/*
    // Many Parameters
    TPaveText *pt1 = new TPaveText(0.35,0.77,0.70,0.99,"NB NDC");
    pt1->AddText(Form("Fit Range: %.0f to %.0f keV -- #chi^{2}/NDF: %0.3f", dFitMin, dFitMax, (GetChiSquare()/((dFitMax-dFitMin)/dBinSize - dNumParameters)) ));
    // pt1->AddText(Form("#chi^{2}/NDF: %0.3f --  Resolution %0.4f", (GetChiSquare()/((dFitMax-dFitMin)/dBinSize - dNumParameters)) ,fParameters[8]));
    pt1->AddText(Form("Frame Th: %0.2E#pm%0.2E --- TShield Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[11], fParError[11] ));
    pt1->AddText(Form("50mK Th: %0.2E#pm%0.2E --- 600mK Th: %0.2E#pm%0.2E", fParameters[12], fParError[12], fParameters[13], fParError[13] ));
    pt1->AddText(Form("IVC Th: %0.2E#pm%0.2E --- OVC Th: %0.2E#pm%0.2E", fParameters[1], fParError[1], fParameters[17], fParError[17] ));
    pt1->AddText(Form("Frame Ra: %0.2E#pm%0.2E --- TShield Ra: %0.2E#pm%0.2E", fParameters[2], fParError[2], fParameters[14], fParError[14] ));
    pt1->AddText(Form("50mK Ra: %0.2E#pm%0.2E --- 600mK Ra: %0.2E#pm%0.2E", fParameters[15], fParError[15], fParameters[16], fParError[16] ));
    pt1->AddText(Form("IVC Ra: %0.2E#pm%0.2E --- OVC Ra: %0.2E#pm%0.2E", fParameters[3], fParError[3], fParameters[18], fParError[18] ));
    pt1->AddText(Form("Close K: %0.2E#pm%0.2E --- Far K: %0.2E#pm%0.2E", fParameters[4], fParError[4], fParameters[5], fParError[5] ));
    pt1->AddText(Form("Close Co: %0.2E#pm%0.2E --- Far Co: %0.2E#pm%0.2E", fParameters[6], fParError[6], fParameters[7], fParError[7] ));
    pt1->AddText(Form("Bi-207: %0.2E#pm%0.2E --- NDBD: %0.2E#pm%0.2E", fParameters[10], fParError[10], fParameters[9], fParError[9] ));
*/
    pt1->Draw();

    TLegend *legfit1 = new TLegend(0.8,0.8,0.97,0.97);
    legfit1->AddEntry(fModelTotM1, "Total PDF", "l");
    legfit1->AddEntry(fModelTotThM1, "Total Th-232", "l");
    legfit1->AddEntry(fModelTotRaM1, "Total Ra-226", "l");
    legfit1->AddEntry(fModelTotKM1, "Total K-40", "l");
    legfit1->AddEntry(fModelTotCoM1, "Total Co-60", "l");
    legfit1->AddEntry(fModelTotNDBDM1, "NDBD", "l");
    legfit1->AddEntry(fModelTotBiM1, "Bi-207", "l");
    legfit1->Draw();





    TCanvas *c2 = new TCanvas("c2", "c2", 1200, 800);
    c2->SetLogy();

    ///// Draw Data M2
    fDataHistoM2->SetLineColor(1);
    fDataHistoM2->SetLineWidth(2);
    fDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
    fDataHistoM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)/yr", dBinSize));
    fDataHistoM2->SetMaximum(9000);
    fDataHistoM2->GetXaxis()->SetRange(1/dBinSize-5, 2650/dBinSize+5);
    fDataHistoM2->Draw("E");

  
    fModelTotM2->SetLineColor(2);
    fModelTotM2->Draw("SAME");

    fModelTotThM2->SetLineColor(3);
    fModelTotThM2->SetLineStyle(2);
    fModelTotRaM2->SetLineColor(4);
    fModelTotRaM2->SetLineStyle(2);
    fModelTotKM2->SetLineColor(6);
    fModelTotKM2->SetLineStyle(2);
    fModelTotCoM2->SetLineColor(7);
    fModelTotCoM2->SetLineStyle(2);
    fModelTotNDBDM2->SetLineColor(42);
    fModelTotNDBDM2->SetLineStyle(2);
    fModelTotBiM2->SetLineColor(5);
    fModelTotBiM2->SetLineStyle(2);

    fModelTotThM2->Draw("SAME");
    fModelTotRaM2->Draw("SAME");
    fModelTotKM2->Draw("SAME");
    fModelTotCoM2->Draw("SAME");
    fModelTotNDBDM2->Draw("SAME");
    fModelTotBiM2->Draw("SAME");    

    // Few Parameters
    TPaveText *pt2 = new TPaveText(0.35,0.78,0.70,0.98,"NB NDC");
    pt2->AddText(Form("Fit Range (M2): %.0f to %.0f keV -- #chi^{2}/NDF: %0.3f", dFitMin, dFitMax, (GetChiSquare()/((dFitMax-dFitMin)/dBinSize - dNumParameters)) ));
    pt2->AddText(Form("Close Th: %0.2E#pm%0.2E --- Far Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
    // pt2->AddText(Form("Surface Th: %0.2E#pm%0.2E --- Far Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
    pt2->AddText(Form("Close Ra: %0.2E#pm%0.2E --- Far Ra: %0.2E#pm%0.2E", fParameters[2], fParError[2], fParameters[3], fParError[3] ));
    pt2->AddText(Form("Close K: %0.2E#pm%0.2E --- Far K: %0.2E#pm%0.2E", fParameters[4], fParError[4], fParameters[5], fParError[5] ));
    pt2->AddText(Form("Close Co: %0.2E#pm%0.2E --- Far Co: %0.2E#pm%0.2E", fParameters[6], fParError[6], fParameters[7], fParError[7] ));
    pt2->AddText(Form("Bi-207: %0.2E#pm%0.2E --- NDBD: %0.2E#pm%0.2E", fParameters[10], fParError[10], fParameters[9], fParError[9] ));


/*
    // Many Parameters
    TPaveText *pt2 = new TPaveText(0.35,0.77,0.70,0.99,"NB NDC");
    pt2->AddText(Form("Fit Range: %.0f to %.0f keV -- #chi^{2}/NDF: %0.3f", dFitMin, dFitMax, (GetChiSquare()/((dFitMax-dFitMin)/dBinSize - dNumParameters)) ));
    // pt2->AddText(Form("#chi^{2}/NDF: %0.3f --  Resolution %0.4f", (GetChiSquare()/((dFitMax-dFitMin)/dBinSize - dNumParameters)) ,fParameters[8]));
    pt2->AddText(Form("Frame Th: %0.2E#pm%0.2E --- TShield Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[11], fParError[11] ));
    pt2->AddText(Form("50mK Th: %0.2E#pm%0.2E --- 600mK Th: %0.2E#pm%0.2E", fParameters[12], fParError[12], fParameters[13], fParError[13] ));
    pt2->AddText(Form("IVC Th: %0.2E#pm%0.2E --- OVC Th: %0.2E#pm%0.2E", fParameters[1], fParError[1], fParameters[17], fParError[17] ));
    pt2->AddText(Form("Frame Ra: %0.2E#pm%0.2E --- TShield Ra: %0.2E#pm%0.2E", fParameters[2], fParError[2], fParameters[14], fParError[14] ));
    pt2->AddText(Form("50mK Ra: %0.2E#pm%0.2E --- 600mK Ra: %0.2E#pm%0.2E", fParameters[15], fParError[15], fParameters[16], fParError[16] ));
    pt2->AddText(Form("IVC Ra: %0.2E#pm%0.2E --- OVC Ra: %0.2E#pm%0.2E", fParameters[3], fParError[3], fParameters[18], fParError[18] ));
    pt2->AddText(Form("Close K: %0.2E#pm%0.2E --- Far K: %0.2E#pm%0.2E", fParameters[4], fParError[4], fParameters[5], fParError[5] ));
    pt2->AddText(Form("Close Co: %0.2E#pm%0.2E --- Far Co: %0.2E#pm%0.2E", fParameters[6], fParError[6], fParameters[7], fParError[7] ));
    pt2->AddText(Form("Bi-207: %0.2E#pm%0.2E --- NDBD: %0.2E#pm%0.2E", fParameters[10], fParError[10], fParameters[9], fParError[9] ));
*/
    pt2->Draw();


    TLegend *legfit2 = new TLegend(0.8,0.8,0.97,0.97);
    legfit2->AddEntry(fModelTotM2, "Total PDF", "l");
    legfit2->AddEntry(fModelTotThM2, "Total Th-232", "l");
    legfit2->AddEntry(fModelTotRaM2, "Total Ra-226", "l");
    legfit2->AddEntry(fModelTotKM2, "Total K-40", "l");
    legfit2->AddEntry(fModelTotCoM2, "Total Co-60", "l");
    legfit2->AddEntry(fModelTotNDBDM2, "NDBD", "l");
    legfit2->AddEntry(fModelTotBiM2, "Bi-207", "l");

    legfit2->Draw();


  }


	// Residuals
	TCanvas *cResidual1 = new TCanvas("cResidual1", "cResidual1", 1200, 800);
	hResidualDistM1 = CalculateResiduals(fModelTotM1, fDataHistoM1);
	hResidualDistM1->SetName("Residuals");
	hResidualDistM1->GetXaxis()->SetTitle("Energy (keV)");
	// hResidualDistM1->GetXaxis()->SetTitleSize(0.04);
	// hResidualDistM1->GetXaxis()->SetLabelSize(0.05);
	// hResidualDistM1->GetYaxis()->SetLabelSize(0.05);
	// hResidualDistM1->GetYaxis()->SetTitleSize(0.04);	
	hResidualDistM1->GetYaxis()->SetTitle("Residuals (#sigma)");

	hResidualDistM1->GetXaxis()->SetRange(dFitMin/dBinSize-5, dFitMax/dBinSize+5);
	hResidualDistM1->Draw("E");

  TCanvas *cResidual2 = new TCanvas("cResidual2", "cResidual2", 1200, 800);
  hResidualDistM2 = CalculateResiduals(fModelTotM2, fDataHistoM2);
  hResidualDistM2->SetName("Residuals");
  hResidualDistM2->GetXaxis()->SetTitle("Energy (keV)");
  // hResidualDistM2->GetXaxis()->SetTitleSize(0.04);
  // hResidualDistM2->GetXaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetTitleSize(0.04); 
  hResidualDistM2->GetYaxis()->SetTitle("Residuals (#sigma)");

  hResidualDistM2->GetXaxis()->SetRange(dFitMin/dBinSize-5, dFitMax/dBinSize+5);
  hResidualDistM2->Draw("E");


  // Prints out covariance and correlation matrices
  // minuit.mnmatu(1);

	return true;
   
 }

// Draws background data, must Initialize first!
void TBackgroundModel::DrawBkg()
{

 	gStyle->SetOptStat(0);
  gStyle->SetOptFit();
 	// gStyle->SetOptTitle(0);	
  TCanvas *cBkg1 = new TCanvas("cBkg1", "cBkg1", 1200, 800);
  cBkg1->SetLogy();
  fDataHistoM1->SetLineColor(1);
  fDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
  fDataHistoM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)/yr", dBinSize));
  fDataHistoM1->Draw();

  TCanvas *cBkg2 = new TCanvas("cBkg2", "cBkg2", 1200, 800);
  cBkg2->SetLogy();
  fDataHistoM2->SetLineColor(1);
  fDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
  fDataHistoM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)/yr", dBinSize));
  fDataHistoM2->Draw();

}

void TBackgroundModel::DrawMC()
{
// Draws all MC spectra, must Initialize first!

 	gStyle->SetOptStat(0);
 	gStyle->SetOptTitle(0);

 	TLegend *legth1 = new TLegend(0.75,0.80,0.925,0.925);
 	TLegend *legra1 = new TLegend(0.75,0.80,0.925,0.925);
 	TLegend *legk1 = new TLegend(0.75,0.80,0.925,0.925);
 	TLegend *legco1 = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legndbd1 = new TLegend(0.8,0.890,0.925,0.925);
  TLegend *legbi1 = new TLegend(0.8,0.890,0.925,0.925);
  TLegend *legths1 = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legras1 = new TLegend(0.75,0.80,0.925,0.925);

  TLegend *legth1 = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legra1 = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legk1 = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legco1 = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legndbd1 = new TLegend(0.8,0.890,0.925,0.925);
  TLegend *legbi1 = new TLegend(0.8,0.890,0.925,0.925);
  TLegend *legths1 = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legras1 = new TLegend(0.75,0.80,0.925,0.925);


  TCanvas *cThSurf1 = new TCanvas("cThSurf1", "cThSurf1", 1200, 800);
  cThSurf1->SetLogy();

  fModelFrameThM1->SetLineColor(1);
  fModelFrameThM1->GetXaxis()->SetTitle("Energy (keV)");
  fModelFrameThM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));    
  fModelFrameThS01M1->SetLineColor(2);
  fModelFrameThS1M1->SetLineColor(3);
  fModelFrameThS10M1->SetLineColor(4);
  fModelFrameThS100M1->SetLineColor(6);

  fModelFrameThM1->DrawNormalized();
  fModelFrameThS01M1->DrawNormalized("SAME");
  fModelFrameThS1M1->DrawNormalized("SAME");
  fModelFrameThS10M1->DrawNormalized("SAME");
  fModelFrameThS100M1->DrawNormalized("SAME");
  fModelFrameThM1->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);


  legths1->AddEntry(fModelFrameThM1, "Frame Bulk (M2)" ,"l");
  legths1->AddEntry(fModelFrameThS01M1, "Frame Surface 0.1 #mum (M2)", "l");
  legths1->AddEntry(fModelFrameThS1M1, "Frame Surface 1 #mum (M2)" ,"l");
  legths1->AddEntry(fModelFrameThS10M1, "Frame Surface 10 #mum (M2)" ,"l");
  legths1->AddEntry(fModelFrameThS100M1, "Frame Surface 100 #mum (M2)" ,"l");
  legths1->Draw();

  TCanvas *cThSurf2 = new TCanvas("cThSurf2", "cThSurf2", 1200, 800);
  cThSurf2->SetLogy();

  fModelFrameThM2->SetLineColor(1);
  fModelFrameThM2->GetXaxis()->SetTitle("Energy (keV)");
  fModelFrameThM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));    
  fModelFrameThS01M2->SetLineColor(2);
  fModelFrameThS1M2->SetLineColor(3);
  fModelFrameThS10M2->SetLineColor(4);
  fModelFrameThS100M2->SetLineColor(6);

  fModelFrameThM2->DrawNormalized();
  fModelFrameThS01M2->DrawNormalized("SAME");
  fModelFrameThS1M2->DrawNormalized("SAME");
  fModelFrameThS10M2->DrawNormalized("SAME");
  fModelFrameThS100M2->DrawNormalized("SAME");
  fModelFrameThM2->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);


  legths2->AddEntry(fModelFrameThM2, "Frame Bulk (M2)" ,"l");
  legths2->AddEntry(fModelFrameThS01M2, "Frame Surface 0.1 #mum (M2)", "l");
  legths2->AddEntry(fModelFrameThS1M2, "Frame Surface 1 #mum (M2)" ,"l");
  legths2->AddEntry(fModelFrameThS10M2, "Frame Surface 10 #mum (M2)" ,"l");
  legths2->AddEntry(fModelFrameThS100M2, "Frame Surface 100 #mum (M2)" ,"l");
  legths2->Draw();



  TCanvas *cRaSurf1 = new TCanvas("cRaSurf1", "cRaSurf1", 1200, 800);
  cRaSurf1->SetLogy();

  fModelFrameRaM1->SetLineColor(1);
  fModelFrameRaM1->GetXaxis()->SetTitle("Energy (keV)");
  fModelFrameRaM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));      
  fModelFrameRaS01M1->SetLineColor(2);
  fModelFrameRaS1M1->SetLineColor(3);
  fModelFrameRaS10M1->SetLineColor(4);
  fModelFrameRaS100M1->SetLineColor(6);

  fModelFrameRaM1->DrawNormalized();
  fModelFrameRaS01M1->DrawNormalized("SAME");
  fModelFrameRaS1M1->DrawNormalized("SAME");
  fModelFrameRaS10M1->DrawNormalized("SAME");
  fModelFrameRaS100M1->DrawNormalized("SAME");
  fModelFrameRaM1->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);


  legras1->AddEntry(fModelFrameRaM1, "Frame Bulk (M1)" ,"l");
  legras1->AddEntry(fModelFrameRaS01M1, "Frame Surface 0.1 #mum (M1)", "l");
  legras1->AddEntry(fModelFrameRaS1M1, "Frame Surface 1 #mum (M1)" ,"l");
  legras1->AddEntry(fModelFrameRaS10M1, "Frame Surface 10 #mum (M1)" ,"l");
  legras1->AddEntry(fModelFrameRaS100M1, "Frame Surface 100 #mum (M1)" ,"l");
  legras1->Draw();

  TCanvas *cRaSurf2 = new TCanvas("cRaSurf2", "cRaSurf2", 1200, 800);
  cRaSurf2->SetLogy();

  fModelFrameRaM2->SetLineColor(1);
  fModelFrameRaM2->GetXaxis()->SetTitle("Energy (keV)");
  fModelFrameRaM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));      
  fModelFrameRaS01M2->SetLineColor(2);
  fModelFrameRaS1M2->SetLineColor(3);
  fModelFrameRaS10M2->SetLineColor(4);
  fModelFrameRaS100M2->SetLineColor(6);

  fModelFrameRaM2->DrawNormalized();
  fModelFrameRaS01M2->DrawNormalized("SAME");
  fModelFrameRaS1M2->DrawNormalized("SAME");
  fModelFrameRaS10M2->DrawNormalized("SAME");
  fModelFrameRaS100M2->DrawNormalized("SAME");
  fModelFrameRaM2->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);


  legras2->AddEntry(fModelFrameRaM2, "Frame Bulk (M2)" ,"l");
  legras2->AddEntry(fModelFrameRaS01M2, "Frame Surface 0.1 #mum (M2)", "l");
  legras2->AddEntry(fModelFrameRaS1M2, "Frame Surface 1 #mum (M2)" ,"l");
  legras2->AddEntry(fModelFrameRaS10M2, "Frame Surface 10 #mum (M2)" ,"l");
  legras2->AddEntry(fModelFrameRaS100M2, "Frame Surface 100 #mum (M2)" ,"l");
  legras2->Draw();



  TCanvas *cTh2321 = new TCanvas("cTh2321", "cTh2321", 1200, 800);
  cTh2321->SetLogy();

  fModelFrameThM1->SetLineColor(1);
  fModelFrameThM1->GetXaxis()->SetTitle("Energy (keV)");
  fModelFrameThM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));      
  fModelTShieldThM1->SetLineColor(2);
  fModel50mKThM1->SetLineColor(3);
  fModel600mKThM1->SetLineColor(4);
  fModelIVCThM1->SetLineColor(6);
  fModelOVCThM1->SetLineColor(7);

  fModelFrameThM1->DrawNormalized();
  fModelTShieldThM1->DrawNormalized("SAME");
  fModel50mKThM1->DrawNormalized("SAME");
  fModel600mKThM1->DrawNormalized("SAME");
  fModelIVCThM1->DrawNormalized("SAME");
  fModelOVCThM1->DrawNormalized("SAME");
  fModelFrameThM1->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);

  legth1->AddEntry(fModelFrameThM1, "Frame (M1)" ,"l");
  legth1->AddEntry(fModelTShieldThM1, "TShield (M1)", "l");
  legth1->AddEntry(fModel50mKThM1, "50mK (M1)" ,"l");
  legth1->AddEntry(fModel600mKThM1, "600mK (M1)" ,"l");
  legth1->AddEntry(fModelIVCThM1, "IVC (M1)" ,"l");
  legth1->AddEntry(fModelOVCThM1, "OVC (M1)" ,"l");
  legth1->Draw();

  TCanvas *cTh2322 = new TCanvas("cTh2322", "cTh2322", 1200, 800);
  cTh2322->SetLogy();

  fModelFrameThM2->SetLineColor(1);
  fModelFrameThM2->GetXaxis()->SetTitle("Energy (keV)");
  fModelFrameThM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));      
  fModelTShieldThM2->SetLineColor(2);
  fModel50mKThM2->SetLineColor(3);
  fModel600mKThM2->SetLineColor(4);
  fModelIVCThM2->SetLineColor(6);
  fModelOVCThM2->SetLineColor(7);

  fModelFrameThM2->DrawNormalized();
  fModelTShieldThM2->DrawNormalized("SAME");
  fModel50mKThM2->DrawNormalized("SAME");
  fModel600mKThM2->DrawNormalized("SAME");
  fModelIVCThM2->DrawNormalized("SAME");
  fModelOVCThM2->DrawNormalized("SAME");
  fModelFrameThM2->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);

  legth2->AddEntry(fModelFrameThM2, "Frame (M2)" ,"l");
  legth2->AddEntry(fModelTShieldThM2, "TShield (M2)", "l");
  legth2->AddEntry(fModel50mKThM2, "50mK (M2)" ,"l");
  legth2->AddEntry(fModel600mKThM2, "600mK (M2)" ,"l");
  legth2->AddEntry(fModelIVCThM2, "IVC (M2)" ,"l");
  legth2->AddEntry(fModelOVCThM2, "OVC (M2)" ,"l");
  legth2->Draw();



  TCanvas *cRa2261 = new TCanvas("cRa2261", "cRa2261", 1200, 800);
  cRa2261->SetLogy();

  fModelFrameRaM1->SetLineColor(1);
  fModelFrameRaM1->GetXaxis()->SetTitle("Energy (keV)");
  fModelFrameRaM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));      
  fModelTShieldRaM1->SetLineColor(2);
  fModel50mKRaM1->SetLineColor(3);
  fModel600mKRaM1->SetLineColor(4);
  fModelIVCRaM1->SetLineColor(6);
  fModelOVCRaM1->SetLineColor(7);

  fModelFrameRaM1->DrawNormalized();
  fModelTShieldRaM1->DrawNormalized("SAME");
  fModel50mKRaM1->DrawNormalized("SAME");
  fModel600mKRaM1->DrawNormalized("SAME");
  fModelIVCRaM1->DrawNormalized("SAME");
  fModelOVCRaM1->DrawNormalized("SAME");
  fModelFrameRaM1->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);

  legra1->AddEntry(fModelFrameRaM1, "Frame (M1)" ,"l");
  legra1->AddEntry(fModelTShieldRaM1, "TShield (M1)", "l");
  legra1->AddEntry(fModel50mKRaM1, "50mK (M1)" ,"l");
  legra1->AddEntry(fModel600mKRaM1, "600mK (M1)" ,"l");
  legra1->AddEntry(fModelIVCRaM1, "IVC (M1)" ,"l");
  legra1->AddEntry(fModelOVCRaM1, "OVC (M1)" ,"l");
  legra1->Draw();

  TCanvas *cRa2262 = new TCanvas("cRa2262", "cRa2262", 1200, 800);
  cRa2262->SetLogy();

  fModelFrameRaM2->SetLineColor(1);
  fModelFrameRaM2->GetXaxis()->SetTitle("Energy (keV)");
  fModelFrameRaM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));      
  fModelTShieldRaM2->SetLineColor(2);
  fModel50mKRaM2->SetLineColor(3);
  fModel600mKRaM2->SetLineColor(4);
  fModelIVCRaM2->SetLineColor(6);
  fModelOVCRaM2->SetLineColor(7);

  fModelFrameRaM2->DrawNormalized();
  fModelTShieldRaM2->DrawNormalized("SAME");
  fModel50mKRaM2->DrawNormalized("SAME");
  fModel600mKRaM2->DrawNormalized("SAME");
  fModelIVCRaM2->DrawNormalized("SAME");
  fModelOVCRaM2->DrawNormalized("SAME");
  fModelFrameRaM2->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);

  legra2->AddEntry(fModelFrameRaM2, "Frame (M2)" ,"l");
  legra2->AddEntry(fModelTShieldRaM2, "TShield (M2)", "l");
  legra2->AddEntry(fModel50mKRaM2, "50mK (M2)" ,"l");
  legra2->AddEntry(fModel600mKRaM2, "600mK (M2)" ,"l");
  legra2->AddEntry(fModelIVCRaM2, "IVC (M2)" ,"l");
  legra2->AddEntry(fModelOVCRaM2, "OVC (M2)" ,"l");
  legra2->Draw();



  TCanvas *cK401 = new TCanvas("cK401", "cK401", 1200, 800);
  cK401->SetLogy();

  fModelFrameKM1->SetLineColor(1);
  fModelFrameKM1->GetXaxis()->SetTitle("Energy (keV)");
  fModelFrameKM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));  
  fModelTShieldKM1->SetLineColor(2);
  fModel50mKKM1->SetLineColor(3);
  fModel600mKKM1->SetLineColor(4);
  fModelIVCKM1->SetLineColor(6);
  fModelOVCKM1->SetLineColor(7);

  fModelFrameKM1->DrawNormalized();
  fModelTShieldKM1->DrawNormalized("SAME");
  fModel50mKKM1->DrawNormalized("SAME");
  fModel600mKKM1->DrawNormalized("SAME");
  fModelIVCKM1->DrawNormalized("SAME");
  fModelOVCKM1->DrawNormalized("SAME");
  fModelFrameKM1->GetXaxis()->SetRange(0/dBinSize, 1600/dBinSize);

  legk1->AddEntry(fModelFrameKM1, "Frame (M1)" ,"l");
  legk1->AddEntry(fModelTShieldKM1, "TShield (M1)", "l");
  legk1->AddEntry(fModel50mKKM1, "50mK (M1)" ,"l");
  legk1->AddEntry(fModel600mKKM1, "600mK (M1)" ,"l");
  legk1->AddEntry(fModelIVCKM1, "IVC (M1)" ,"l");
  legk1->AddEntry(fModelOVCKM1, "OVC (M1)" ,"l");    
  legk1->Draw();

  TCanvas *cK402 = new TCanvas("cK402", "cK402", 1200, 800);
  cK402->SetLogy();

  fModelFrameKM2->SetLineColor(1);
  fModelFrameKM2->GetXaxis()->SetTitle("Energy (keV)");
  fModelFrameKM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));  
  fModelTShieldKM2->SetLineColor(2);
  fModel50mKKM2->SetLineColor(3);
  fModel600mKKM2->SetLineColor(4);
  fModelIVCKM2->SetLineColor(6);
  fModelOVCKM2->SetLineColor(7);

  fModelFrameKM2->DrawNormalized();
  fModelTShieldKM2->DrawNormalized("SAME");
  fModel50mKKM2->DrawNormalized("SAME");
  fModel600mKKM2->DrawNormalized("SAME");
  fModelIVCKM2->DrawNormalized("SAME");
  fModelOVCKM2->DrawNormalized("SAME");
  fModelFrameKM2->GetXaxis()->SetRange(0/dBinSize, 1600/dBinSize);

  legk2->AddEntry(fModelFrameKM2, "Frame (M2)" ,"l");
  legk2->AddEntry(fModelTShieldKM2, "TShield (M2)", "l");
  legk2->AddEntry(fModel50mKKM2, "50mK (M2)" ,"l");
  legk2->AddEntry(fModel600mKKM2, "600mK (M2)" ,"l");
  legk2->AddEntry(fModelIVCKM2, "IVC (M2)" ,"l");
  legk2->AddEntry(fModelOVCKM2, "OVC (M2)" ,"l");    
  legk2->Draw();



  TCanvas *cCo601 = new TCanvas("cCo601", "cCo601", 1200, 800);
  cCo601->SetLogy();

  fModelFrameCoM1->SetLineColor(1);
  fModelFrameCoM1->GetXaxis()->SetTitle("Energy (keV)");
  fModelFrameCoM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));  
  fModelTShieldCoM1->SetLineColor(2);
  fModel50mKCoM1->SetLineColor(3);
  fModel600mKCoM1->SetLineColor(4);
  fModelIVCCoM1->SetLineColor(6);
  fModelOVCCoM1->SetLineColor(7);

  fModelFrameCoM1->Draw();
  fModelTShieldCoM1->Draw("SAME");
  fModel50mKCoM1->Draw("SAME");
  fModel600mKCoM1->Draw("SAME");
  fModelIVCCoM1->Draw("SAME");
  fModelOVCCoM1->Draw("SAME");
  fModelFrameCoM1->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);

  legco1->AddEntry(fModelFrameCoM1, "Frame (M1)" ,"l");
  legco1->AddEntry(fModelTShieldCoM1, "TShield (M1)", "l");
  legco1->AddEntry(fModel50mKCoM1, "50mK (M1)" ,"l");
  legco1->AddEntry(fModel600mKCoM1, "600mK (M1)" ,"l");
  legco1->AddEntry(fModelIVCCoM1, "IVC (M1)" ,"l");
  legco1->AddEntry(fModelOVCCoM1, "OVC (M1)" ,"l");
  legco1->Draw();

  TCanvas *cCo602 = new TCanvas("cCo602", "cCo602", 1200, 800);
  cCo602->SetLogy();

  fModelFrameCoM2->SetLineColor(1);
  fModelFrameCoM2->GetXaxis()->SetTitle("Energy (keV)");
  fModelFrameCoM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));  
  fModelTShieldCoM2->SetLineColor(2);
  fModel50mKCoM2->SetLineColor(3);
  fModel600mKCoM2->SetLineColor(4);
  fModelIVCCoM2->SetLineColor(6);
  fModelOVCCoM2->SetLineColor(7);

  fModelFrameCoM2->Draw();
  fModelTShieldCoM2->Draw("SAME");
  fModel50mKCoM2->Draw("SAME");
  fModel600mKCoM2->Draw("SAME");
  fModelIVCCoM2->Draw("SAME");
  fModelOVCCoM2->Draw("SAME");
  fModelFrameCoM2->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);

  legco2->AddEntry(fModelFrameCoM2, "Frame (M2)" ,"l");
  legco2->AddEntry(fModelTShieldCoM2, "TShield (M2)", "l");
  legco2->AddEntry(fModel50mKCoM2, "50mK (M2)" ,"l");
  legco2->AddEntry(fModel600mKCoM2, "600mK (M2)" ,"l");
  legco2->AddEntry(fModelIVCCoM2, "IVC (M2)" ,"l");
  legco2->AddEntry(fModelOVCCoM2, "OVC (M2)" ,"l");
  legco2->Draw();



  TCanvas *cNDBD1 = new TCanvas("cNDBD1", "cNDBD1", 1200, 800);
  cNDBD1->SetLogy();
  fModelNDBDM1->SetLineColor(1);
  fModelNDBDM1->GetXaxis()->SetTitle("Energy (keV)");
  fModelNDBDM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));  
  fModelNDBDM1->Draw();
  fModelNDBDM1->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);
  legndbd1->AddEntry(fModelNDBDM1, "0#nu#beta#beta (M1)" ,"l");
  legndbd1->Draw();

  TCanvas *cNDBD2 = new TCanvas("cNDBD2", "cNDBD2", 1200, 800);
  cNDBD2->SetLogy();
  fModelNDBDM2->SetLineColor(1);
  fModelNDBDM2->GetXaxis()->SetTitle("Energy (keV)");
  fModelNDBDM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));  
  fModelNDBDM2->Draw();
  fModelNDBDM2->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);
  legndbd2->AddEntry(fModelNDBDM2, "0#nu#beta#beta (M2)" ,"l");
  legndbd2->Draw();


  TCanvas *cBi1 = new TCanvas("cBi1", "cBi1", 1200, 800);
  cBi->SetLogy();
  fModelBiM1->SetLineColor(1);
  fModelBiM1->GetXaxis()->SetTitle("Energy (keV)");
  fModelBiM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));      
  fModelBiM1->Draw();
  fModelBiM1->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);
  legbi1->AddEntry(fModelBiM1, "Bi-207 (M1)" ,"l");
  legbi1->Draw();

  TCanvas *cBi2 = new TCanvas("cBi2", "cBi2", 1200, 800);
  cBi2->SetLogy();
  fModelBiM2->SetLineColor(1);
  fModelBiM2->GetXaxis()->SetTitle("Energy (keV)");
  fModelBiM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));      
  fModelBiM2->Draw();
  fModelBiM2->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);
  legbi2->AddEntry(fModelBiM2, "Bi-207 (M2)" ,"l");
  legbi2->Draw();

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
			datam1_i = fDataHistoM1->GetBinContent(i); // For real data
      datam2_i = fDataHistoM2->GetBinContent(i); // For real data
		}

    // From MC
		modelm1_i = fModelTotM1->GetBinContent(i);
    modelm2_i = fModelTotM2->GetBinContent(i);


		// Log-likelihood Chi-Squared
    // Avoiding 0's... correct or no?
		if(modelm1_i != 0 && datam1_i != 0)
		{
      // M1 portion
			chiSquare += 2 * (modelm1_i - datam1_i + datam1_i * TMath::Log(datam1_i/modelm1_i));
		}

    if(modelm2_i != 0 && datam2_i != 0)
    {
      // M2 portion
      chiSquare += 2 * (modelm2_i - datam2_i + datam2_i * TMath::Log(datam2_i/modelm2_i));
    }
	}

	return chiSquare;
}

// Gets efficiency of MC as fraction of normalization -- not used at all right now
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
  // Load M1
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



  // Load M2
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

	dDataIntegral = fDataHistoM1->Integral(1, dNBins);
  int dDataIntegralTot = qtree->GetEntries();

  cout << "Total Events in background spectrum: " << dDataIntegralTot << endl; 
	cout << "Events in background spectrum (M1): " << dDataIntegral << endl;
  cout << "Events in background spectrum (M2): " << fDataHistoM2->Integral(1, dNBins) << endl;

  // Scale by Live-time (ds 2061 - 2100) 14647393.0 seconds
  fDataHistoM1->Scale(1/(14647393.0 * dSecToYears));
  fDataHistoM2->Scale(1/(14647393.0 * dSecToYears));  

  cout << "Normalized Data using Livetime of: " << 14647393.0 * dSecToYears << " years" <<endl;

}


// Loads MC files into Trees
TChain *TBackgroundModel::LoadMC(std::string dDir, std::string dLocation, std::string dSource, std::string dSType, int dMult)
{
    TChain *outTree = new TChain("outTree");
    outTree->Add(Form("%s%s-%s-%s-M%d-r0.0500.root", dDir.c_str(), dLocation.c_str(), dSource.c_str(), dSType.c_str(), dMult));

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
  // cout<< "Par11 = "  << fParameters[11] << " +/- " << fParError[11] << endl;
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


// For custom Smearing with resolution, currently constant resolution 
TH1D *TBackgroundModel::SmearMC(TH1D *hMC, TH1D *hSMC, double resolution)
{
	// Reset previously Modeled histogram
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
	if(fModelTotM1 == NULL) 
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

  /////////////////////////////////////
  //// Few Parameters ////////////////
  ////////////////////////////////////
  // fModelTot->Add( fModelTShieldThS10,    fParameters[0]);
  // M1
  fModelTotM1->Add( fModelFrameThM1,    fParameters[0]);
  fModelTotM1->Add( fModelTShieldThM1,  fParameters[0]);  
  fModelTotM1->Add( fModel50mKThM1,     fParameters[0]);
  fModelTotM1->Add( fModel600mKThM1,    fParameters[0]);
  fModelTotM1->Add( fModelIVCThM1,      fParameters[1]);
  fModelTotM1->Add( fModelOVCThM1,      fParameters[1]);

  fModelTotM1->Add( fModelFrameRaM1,    fParameters[2]);
  fModelTotM1->Add( fModelTShieldRaM1,  fParameters[2]);  
  fModelTotM1->Add( fModel50mKRaM1,     fParameters[2]);
  fModelTotM1->Add( fModel600mKRaM1,    fParameters[2]);
  fModelTotM1->Add( fModelIVCRaM1,      fParameters[3]);
  fModelTotM1->Add( fModelOVCRaM1,      fParameters[3]);

  fModelTotM1->Add( fModelFrameKM1,    fParameters[4]);
  fModelTotM1->Add( fModelTShieldKM1,  fParameters[4]);
  fModelTotM1->Add( fModel50mKKM1,     fParameters[4]);
  fModelTotM1->Add( fModel600mKKM1,    fParameters[4]);
  fModelTotM1->Add( fModelIVCKM1,      fParameters[5]);
  fModelTotM1->Add( fModelOVCKM1,      fParameters[5]); 

  fModelTotM1->Add( fModelFrameCoM1,    fParameters[6]);
  fModelTotM1->Add( fModelTShieldCoM1,  fParameters[6]);
  fModelTotM1->Add( fModel50mKCoM1,     fParameters[6]);
  fModelTotM1->Add( fModel600mKCoM1,    fParameters[6]);
  fModelTotM1->Add( fModelIVCCoM1,      fParameters[7]);
  fModelTotM1->Add( fModelOVCCoM1,      fParameters[7]);  

  fModelTotM1->Add( fModelNDBDM1,      fParameters[9]);  
  fModelTotM1->Add( fModelBiM1,        fParameters[10]);  



  // M2
  fModelTotM2->Add( fModelFrameThM2,    fParameters[0]);
  fModelTotM2->Add( fModelTShieldThM2,  fParameters[0]);  
  fModelTotM2->Add( fModel50mKThM2,     fParameters[0]);
  fModelTotM2->Add( fModel600mKThM2,    fParameters[0]);
  fModelTotM2->Add( fModelIVCThM2,      fParameters[1]);
  fModelTotM2->Add( fModelOVCThM2,      fParameters[1]);

  fModelTotM2->Add( fModelFrameRaM2,    fParameters[2]);
  fModelTotM2->Add( fModelTShieldRaM2,  fParameters[2]);  
  fModelTotM2->Add( fModel50mKRaM2,     fParameters[2]);
  fModelTotM2->Add( fModel600mKRaM2,    fParameters[2]);
  fModelTotM2->Add( fModelIVCRaM2,      fParameters[3]);
  fModelTotM2->Add( fModelOVCRaM2,      fParameters[3]);

  fModelTotM2->Add( fModelFrameKM2,    fParameters[4]);
  fModelTotM2->Add( fModelTShieldKM2,  fParameters[4]);
  fModelTotM2->Add( fModel50mKKM2,     fParameters[4]);
  fModelTotM2->Add( fModel600mKKM2,    fParameters[4]);
  fModelTotM2->Add( fModelIVCKM2,      fParameters[5]);
  fModelTotM2->Add( fModelOVCKM2,      fParameters[5]); 

  fModelTotM2->Add( fModelFrameCoM2,    fParameters[6]);
  fModelTotM2->Add( fModelTShieldCoM2,  fParameters[6]);
  fModelTotM2->Add( fModel50mKCoM2,     fParameters[6]);
  fModelTotM2->Add( fModel600mKCoM2,    fParameters[6]);
  fModelTotM2->Add( fModelIVCCoM2,      fParameters[7]);
  fModelTotM2->Add( fModelOVCCoM2,      fParameters[7]);  

  fModelTotM2->Add( fModelNDBDM2,      fParameters[9]);  
  fModelTotM2->Add( fModelBiM2,        fParameters[10]);  


  /////////////////////////////////////
  //// Many parameters
  ////////////////////////////////////
/*
  fModelTot->Add( fModelFrameTh,    fParameters[0]);
  fModelTot->Add( fModelTShieldTh,  fParameters[11]);  
  fModelTot->Add( fModel50mKTh,     fParameters[12]);
  fModelTot->Add( fModel600mKTh,    fParameters[13]);
  fModelTot->Add( fModelIVCTh,      fParameters[1]);
  fModelTot->Add( fModelOVCTh,      fParameters[17]);

  fModelTot->Add( fModelFrameRa,    fParameters[2]);
  fModelTot->Add( fModelTShieldRa,  fParameters[14]);  
  fModelTot->Add( fModel50mKRa,     fParameters[15]);
  fModelTot->Add( fModel600mKRa,    fParameters[16]);
  fModelTot->Add( fModelIVCRa,      fParameters[3]);
  fModelTot->Add( fModelOVCRa,      fParameters[18]);

  fModelTot->Add( fModelFrameK,     fParameters[4]);
  fModelTot->Add( fModelTShieldK,   fParameters[4]);
  fModelTot->Add( fModel50mKK,      fParameters[4]);
  fModelTot->Add( fModel600mKK,     fParameters[4]);
  fModelTot->Add( fModelIVCK,       fParameters[5]);
  fModelTot->Add( fModelOVCK,       fParameters[5]); 

  fModelTot->Add( fModelFrameCo,    fParameters[6]);
  fModelTot->Add( fModelTShieldCo,  fParameters[6]);
  fModelTot->Add( fModel50mKCo,     fParameters[6]);
  fModelTot->Add( fModel600mKCo,    fParameters[6]);
  fModelTot->Add( fModelIVCCo,      fParameters[7]);
  fModelTot->Add( fModelOVCCo,      fParameters[7]);  

  fModelTot->Add( fModelNDBD,       fParameters[9]);  

  fModelTot->Add( fModelBi,         fParameters[10]);  
*/
}




