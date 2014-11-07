#include "TMinuit.h"
#include "TBackgroundModel.hh"
#include "TRandom3.h"
#include "TPaveText.h"
#include "TAxis.h"
#include "TLine.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>


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
  Obj->SetParameters(11, x[11]);
  Obj->SetParameters(12, x[12]);
  Obj->SetParameters(13, x[13]);
  Obj->SetParameters(14, x[14]);
  Obj->SetParameters(15, x[15]);
  Obj->SetParameters(16, x[16]);
  Obj->SetParameters(17, x[17]);
  Obj->SetParameters(18, x[18]);
  Obj->SetParameters(19, x[19]);
  Obj->SetParameters(20, x[20]);
  Obj->SetParameters(21, x[21]);
  Obj->SetParameters(22, x[22]);
  Obj->SetParameters(23, x[23]);
  Obj->SetParameters(24, x[24]);
  Obj->SetParameters(25, x[25]);


	// Obj->UpdateModel();
  Obj->UpdateModelAdaptive();

	//implement a method in your class that calculates the quantity you want to minimise, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimise fval
	// fval = Obj->GetChiSquare();
  fval = Obj->GetChiSquareAdaptive();

}


TBackgroundModel::TBackgroundModel(double fFitMin, double fFitMax)
{

  dNumCalls = 0;
  dSecToYears = 1./(60*60*24*365);

  dDataDir =  "/Users/brian/macros/Simulations/Bkg/Unsmeared/";
  dDataIntegral = 0;
  bToyFit = false;
  bUnSmeared = false;

  // Bin size (keV)
  dBinSize = 2;
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
  fDataHistoTot  = new TH1D("fDataHistoTot",  "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM1   = new TH1D("fDataHistoM1",   "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM2   = new TH1D("fDataHistoM2",   "", dNBins, dMinEnergy, dMaxEnergy);

  qtree = new TChain("qredtree");
  // qtree = new TChain("qtree");

  // base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
  // base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
  // base_cut = base_cut && "abs(BaselineSlope)<0.1";
  // base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";

  // Load data here
  qtree->Add("/Users/brian/macros/CUOREZ/Bkg/Q0_DR2_BackgroundSignalData.root"); 
  // qtree->Add("/Users/brian/macros/CUOREZ/Bkg/ReducedBkg-ds*.root");
  qtree->Project("fDataHistoTot", "Energy", base_cut);
  qtree->Project("fDataHistoM1",  "Energy", base_cut && "Multiplicity == 1");
  qtree->Project("fDataHistoM2",  "Energy", base_cut && "Multiplicity == 2");

  cout << "Loaded Data" << endl;

  // Normalizing data (don't!)
  // bin 0 = underflow, bin dNBins = last bin with upper-edge xup Excluded

  dDataIntegral = fDataHistoM1->Integral(1, dNBins);
  int dDataIntegralTot = qtree->GetEntries();

  cout << "Total Events in background spectrum: " << dDataIntegralTot << endl; 
  cout << "Events in background spectrum (M1): " << fDataHistoM1->Integral(1, 3000/dBinSize) << endl;
  cout << "Events in background spectrum (M2): " << fDataHistoM2->Integral(1, 3000/dBinSize) << endl;

  // Scale by Live-time (ds 2061 - 2100) 14647393.0 seconds
  fDataHistoM1->Scale(1/((936398+14647393.0) * dSecToYears));
  fDataHistoM2->Scale(1/((936398+14647393.0) * dSecToYears));  

  cout << "Normalized Data using Livetime of: " << (936398+14647393.0) * dSecToYears << " years" <<endl;

  dAdaptiveVectorM1 = AdaptiveBinning(fDataHistoM1);
  dAdaptiveBinsM1 = dAdaptiveVectorM1.size() - 1;
  dAdaptiveArrayM1 = &dAdaptiveVectorM1[0];
  dAdaptiveVectorM2 = AdaptiveBinning(fDataHistoM2);
  dAdaptiveBinsM2 = dAdaptiveVectorM2.size() - 1;
  dAdaptiveArrayM2 = &dAdaptiveVectorM2[0];

  // Adaptive binning data
  fAdapDataHistoM1   = new TH1D("fAdapDataHistoM1",   "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapDataHistoM2   = new TH1D("fAdapDataHistoM2",   "", dAdaptiveBinsM2, dAdaptiveArrayM2);

  fDataHistoM1->Rebin(dAdaptiveBinsM1, "hnewM1", dAdaptiveArrayM1);
  fDataHistoM2->Rebin(dAdaptiveBinsM2, "hnewM2", dAdaptiveArrayM2);

  for(int i = 1; i <= dAdaptiveBinsM1; i++)
  {
    fAdapDataHistoM1->SetBinContent(i, 2*hnewM1->GetBinContent(i)/hnewM1->GetBinWidth(i));
  }

  for(int i = 1; i <= dAdaptiveBinsM2; i++)
  {
    fAdapDataHistoM2->SetBinContent(i, 2*hnewM2->GetBinContent(i)/hnewM2->GetBinWidth(i));
  }


  dFitMinBinM1 = fAdapDataHistoM1->FindBin(dFitMin);
  dFitMinBinM2 = fAdapDataHistoM2->FindBin(dFitMin);
  dFitMaxBinM1 = fAdapDataHistoM1->FindBin(dFitMax);
  dFitMaxBinM2 = fAdapDataHistoM2->FindBin(dFitMax);

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
  fModelOVCKM1       = new TH1D("fModelOVCKM1",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameCoM1    = new TH1D("fModelFrameCoM1",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldCoM1  = new TH1D("fModelTShieldCoM1","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fModel50mKCoM1     = new TH1D("fModel50mKCoM1",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKCoM1    = new TH1D("fModel600mKCoM1",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCCoM1      = new TH1D("fModelIVCCoM1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCCoM1      = new TH1D("fModelOVCCoM1",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelNDBDM1       = new TH1D("fModelNDBDM1",     "NDBD",     dNBins, dMinEnergy, dMaxEnergy);
  fModel2NDBDM1      = new TH1D("fModel2NDBDM1",    "2NDBD",    dNBins, dMinEnergy, dMaxEnergy);
  fModelBiM1         = new TH1D("fModelBiM1",       "Bi207",    dNBins, dMinEnergy, dMaxEnergy);

  fModelTShieldMnM1  = new TH1D("fModelTShieldMnM1","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCMnM1      = new TH1D("fModelIVCMnM1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelCrystalBi2M1  = new TH1D("fModelCrystalBi2M1","Crystal",  dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameBi2M1    = new TH1D("fModelFrameBi2M1","Frame",  dNBins, dMinEnergy, dMaxEnergy);


  fModelCrystalPt190BM1        = new TH1D("fModelCrystalPt190BM1",    "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPt190S01M1      = new TH1D("fModelCrystalPt190S01M1",  "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPt190S1M1       = new TH1D("fModelCrystalPt190S1M1",   "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPt190S10M1      = new TH1D("fModelCrystalPt190S10M1",  "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPt190S100M1     = new TH1D("fModelCrystalPt190S100M1", "Pt190",    dNBins, dMinEnergy, dMaxEnergy);

  fModelCrystalPb210BM1        = new TH1D("fModelCrystalPb210BM1",    "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPb210S01M1      = new TH1D("fModelCrystalPb210S01M1",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPb210S1M1       = new TH1D("fModelCrystalPb210S1M1",   "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPb210S10M1      = new TH1D("fModelCrystalPb210S10M1",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPb210S100M1     = new TH1D("fModelCrystalPb210S100M1", "Pb210",    dNBins, dMinEnergy, dMaxEnergy);

  fModelCrystalRa226BM1        = new TH1D("fModelCrystalRa226BM1",    "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalRa226S01M1      = new TH1D("fModelCrystalRa226S01M1",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalRa226S1M1       = new TH1D("fModelCrystalRa226S1M1",   "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalRa226S10M1      = new TH1D("fModelCrystalRa226S10M1",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalRa226S100M1     = new TH1D("fModelCrystalRa226S100M1", "Ra226",    dNBins, dMinEnergy, dMaxEnergy);

  fModelFramePb210BM1        = new TH1D("fModelFramePb210BM1",    "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFramePb210S01M1      = new TH1D("fModelFramePb210S01M1",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFramePb210S1M1       = new TH1D("fModelFramePb210S1M1",   "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFramePb210S10M1      = new TH1D("fModelFramePb210S10M1",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFramePb210S100M1     = new TH1D("fModelFramePb210S100M1", "Pb210",    dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameRa226BM1        = new TH1D("fModelFrameRa226BM1",    "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRa226S01M1      = new TH1D("fModelFrameRa226S01M1",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRa226S1M1       = new TH1D("fModelFrameRa226S1M1",   "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRa226S10M1      = new TH1D("fModelFrameRa226S10M1",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRa226S100M1     = new TH1D("fModelFrameRa226S100M1", "Ra226",    dNBins, dMinEnergy, dMaxEnergy);

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
  fModelOVCKM2       = new TH1D("fModelOVCKM2",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameCoM2    = new TH1D("fModelFrameCoM2",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTShieldCoM2  = new TH1D("fModelTShieldCoM2","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fModel50mKCoM2     = new TH1D("fModel50mKCoM2",   "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  fModel600mKCoM2    = new TH1D("fModel600mKCoM2",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCCoM2      = new TH1D("fModelIVCCoM2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fModelOVCCoM2      = new TH1D("fModelOVCCoM2",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelNDBDM2       = new TH1D("fModelNDBDM2",     "NDBD",     dNBins, dMinEnergy, dMaxEnergy);
  fModel2NDBDM2      = new TH1D("fModel2NDBDM2",    "2NDBD",    dNBins, dMinEnergy, dMaxEnergy);
  fModelBiM2         = new TH1D("fModelBiM2",       "Bi207",    dNBins, dMinEnergy, dMaxEnergy);

  fModelTShieldMnM2  = new TH1D("fModelTShieldMnM2","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fModelIVCMnM2      = new TH1D("fModelIVCMnM2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);

  fModelCrystalBi2M2  = new TH1D("fModelCrystalBi2M2","Crystal",  dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameBi2M2    = new TH1D("fModelFrameBi2M2","Frame",  dNBins, dMinEnergy, dMaxEnergy);

  fModelCrystalPt190BM2        = new TH1D("fModelCrystalPtBM2",       "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPt190S01M2      = new TH1D("fModelCrystalPt190S01M2",  "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPt190S1M2       = new TH1D("fModelCrystalPt190S1M2",   "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPt190S10M2      = new TH1D("fModelCrystalPt190S10M2",  "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPt190S100M2     = new TH1D("fModelCrystalPt190S100M2", "Pt190",    dNBins, dMinEnergy, dMaxEnergy);

  fModelCrystalPb210BM2        = new TH1D("fModelCrystalPb210BM2",    "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPb210S01M2      = new TH1D("fModelCrystalPb210S01M2",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPb210S1M2       = new TH1D("fModelCrystalPb210S1M2",   "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPb210S10M2      = new TH1D("fModelCrystalPb210S10M2",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalPb210S100M2     = new TH1D("fModelCrystalPb210S100M2", "Pb210",    dNBins, dMinEnergy, dMaxEnergy);

  fModelCrystalRa226BM2        = new TH1D("fModelCrystalRa226BM2",    "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalRa226S01M2      = new TH1D("fModelCrystalRa226S01M2",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalRa226S1M2       = new TH1D("fModelCrystalRa226S1M2",   "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalRa226S10M2      = new TH1D("fModelCrystalRa226S10M2",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelCrystalRa226S100M2     = new TH1D("fModelCrystalRa226S100M2", "Ra226",    dNBins, dMinEnergy, dMaxEnergy);

  fModelFramePb210BM2        = new TH1D("fModelFramePb210BM2",    "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFramePb210S01M2      = new TH1D("fModelFramePb210S01M2",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFramePb210S1M2       = new TH1D("fModelFramePb210S1M2",   "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFramePb210S10M2      = new TH1D("fModelFramePb210S10M2",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFramePb210S100M2     = new TH1D("fModelFramePb210S100M2", "Pb210",    dNBins, dMinEnergy, dMaxEnergy);

  fModelFrameRa226BM2        = new TH1D("fModelFrameRa226BM2",    "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRa226S01M2      = new TH1D("fModelFrameRa226S01M2",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRa226S1M2       = new TH1D("fModelFrameRa226S1M2",   "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRa226S10M2      = new TH1D("fModelFrameRa226S10M2",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fModelFrameRa226S100M2     = new TH1D("fModelFrameRa226S100M2", "Ra226",    dNBins, dMinEnergy, dMaxEnergy);


  // Total model histograms M1
  fModelTotM1      = new TH1D("fModelTotM1",      "Frame",        dNBins, dMinEnergy, dMaxEnergy);  
  fModelTotThM1    = new TH1D("fModelTotThM1",    "Total Th232",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotRaM1    = new TH1D("fModelTotRaM1",    "Total Ra226",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotKM1     = new TH1D("fModelTotKM1",     "Total K40",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTotCoM1    = new TH1D("fModelTotCoM1",    "Total Co60",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotMnM1    = new TH1D("fModelTotMnM1",    "Total Mn54",   dNBins, dMinEnergy, dMaxEnergy);

  fModelTotNDBDM1  = new TH1D("fModelTotNDBDM1",  "Total NDBD",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTot2NDBDM1 = new TH1D("fModelTot2NDBDM1", "Total 2NDBD",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotBiM1    = new TH1D("fModelTotBiM1",   "Total Bi207",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotBi2M1   = new TH1D("fModelTotBi2M1",   "Total Bi210",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotPtM1    = new TH1D("fModelTotPtM1",   "Total Pt190",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotPbM1    = new TH1D("fModelTotPbM1",   "Total Pb210",   dNBins, dMinEnergy, dMaxEnergy);


  // Total model histograms M2
  fModelTotM2      = new TH1D("fModelTotM2",      "Frame",        dNBins, dMinEnergy, dMaxEnergy);  
  fModelTotThM2    = new TH1D("fModelTotThM2",    "Total Th232",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotRaM2    = new TH1D("fModelTotRaM2",    "Total Ra226",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotKM2     = new TH1D("fModelTotKM2",     "Total K40",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTotCoM2    = new TH1D("fModelTotCoM2",    "Total Co60",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotMnM2    = new TH1D("fModelTotMnM2",    "Total Mn54",   dNBins, dMinEnergy, dMaxEnergy);

  fModelTotNDBDM2  = new TH1D("fModelTotNDBDM2",  "Total NDBD",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTot2NDBDM2 = new TH1D("fModelTot2NDBDM2", "Total 2NDBD",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotBiM2    = new TH1D("fModelTotBiM2",    "Total Bi207",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotBi2M2   = new TH1D("fModelTotBi2M2",   "Total Bi210",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotPtM2    = new TH1D("fModelTotPtM2",   "Total Pt190",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotPbM2    = new TH1D("fModelTotPbM2",   "Total Pb210",   dNBins, dMinEnergy, dMaxEnergy);


  // Total adaptive binning created later

  // Modeling
  gaus = new TF1("gaus","gaus(0)", dMinEnergy, dMaxEnergy);
  gaus2 = new TF1("gaus2","gaus(0)", dMinEnergy, dMaxEnergy);



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

  fSmearTShieldMnM1  = new TH1D("fSmearTShieldMnM1","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCMnM1      = new TH1D("fSmearIVCMnM1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearCrystalBi2M1  = new TH1D("fSmearCrystalBi2M1","Crystal",  dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameBi2M1    = new TH1D("fSmearFrameBi2M1","Frame",  dNBins, dMinEnergy, dMaxEnergy);


  fSmearNDBDM1       = new TH1D("fSmearNDBDM1",   "NDBD",       dNBins, dMinEnergy, dMaxEnergy);
  fSmear2NDBDM1      = new TH1D("fSmear2NDBDM1",  "2NDBD",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearBiM1         = new TH1D("fSmearBiM1",     "Bi",         dNBins, dMinEnergy, dMaxEnergy);

  fSmearCrystalPt190BM1        = new TH1D("fSmearCrystalPtBM1",       "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPt190S01M1      = new TH1D("fSmearCrystalPt190S01M1",  "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPt190S1M1       = new TH1D("fSmearCrystalPt190S1M1",   "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPt190S10M1      = new TH1D("fSmearCrystalPt190S10M1",  "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPt190S100M1     = new TH1D("fSmearCrystalPt190S100M1", "Pt190",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearCrystalPb210BM1        = new TH1D("fSmearCrystalPb210BM1",    "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPb210S01M1      = new TH1D("fSmearCrystalPb210S01M1",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPb210S1M1       = new TH1D("fSmearCrystalPb210S1M1",   "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPb210S10M1      = new TH1D("fSmearCrystalPb210S10M1",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPb210S100M1     = new TH1D("fSmearCrystalPb210S100M1", "Pb210",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearCrystalRa226BM1        = new TH1D("fSmearCrystalRa226BM1",    "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalRa226S01M1      = new TH1D("fSmearCrystalRa226S01M1",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalRa226S1M1       = new TH1D("fSmearCrystalRa226S1M1",   "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalRa226S10M1      = new TH1D("fSmearCrystalRa226S10M1",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalRa226S100M1     = new TH1D("fSmearCrystalRa226S100M1", "Ra226",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearFramePb210BM1        = new TH1D("fSmearFramePb210BM1",    "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFramePb210S01M1      = new TH1D("fSmearFramePb210S01M1",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFramePb210S1M1       = new TH1D("fSmearFramePb210S1M1",   "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFramePb210S10M1      = new TH1D("fSmearFramePb210S10M1",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFramePb210S100M1     = new TH1D("fSmearFramePb210S100M1", "Pb210",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameRa226BM1        = new TH1D("fSmearFrameRa226BM1",    "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRa226S01M1      = new TH1D("fSmearFrameRa226S01M1",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRa226S1M1       = new TH1D("fSmearFrameRa226S1M1",   "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRa226S10M1      = new TH1D("fSmearFrameRa226S10M1",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRa226S100M1     = new TH1D("fSmearFrameRa226S100M1", "Ra226",    dNBins, dMinEnergy, dMaxEnergy);


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

  fSmearTShieldMnM2  = new TH1D("fSmearTShieldMnM2","TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fSmearIVCMnM2      = new TH1D("fSmearIVCMnM2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);

  fSmearCrystalBi2M2  = new TH1D("fSmearCrystalBi2M2","Crystal",  dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameBi2M2    = new TH1D("fSmearFrameBi2M2","Frame",  dNBins, dMinEnergy, dMaxEnergy);


  fSmearNDBDM2       = new TH1D("fSmearNDBDM2",   "NDBD",       dNBins, dMinEnergy, dMaxEnergy);
  fSmear2NDBDM2      = new TH1D("fSmear2NDBDM2",  "2NDBD",      dNBins, dMinEnergy, dMaxEnergy);
  fSmearBiM2         = new TH1D("fSmearBiM2",     "Bi",         dNBins, dMinEnergy, dMaxEnergy);

  fSmearCrystalPt190BM2        = new TH1D("fSmearCrystalPtBM2",       "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPt190S01M2      = new TH1D("fSmearCrystalPt190S01M2",  "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPt190S1M2       = new TH1D("fSmearCrystalPt190S1M2",   "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPt190S10M2      = new TH1D("fSmearCrystalPt190S10M2",  "Pt190",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPt190S100M2     = new TH1D("fSmearCrystalPt190S100M2", "Pt190",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearCrystalPb210BM2        = new TH1D("fSmearCrystalPb210BM2",    "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPb210S01M2      = new TH1D("fSmearCrystalPb210S01M2",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPb210S1M2       = new TH1D("fSmearCrystalPb210S1M2",   "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPb210S10M2      = new TH1D("fSmearCrystalPb210S10M2",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalPb210S100M2     = new TH1D("fSmearCrystalPb210S100M2", "Pb210",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearCrystalRa226BM2        = new TH1D("fSmearCrystalRa226BM2",    "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalRa226S01M2      = new TH1D("fSmearCrystalRa226S01M2",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalRa226S1M2       = new TH1D("fSmearCrystalRa226S1M2",   "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalRa226S10M2      = new TH1D("fSmearCrystalRa226S10M2",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearCrystalRa226S100M2     = new TH1D("fSmearCrystalRa226S100M2", "Ra226",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearFramePb210BM2        = new TH1D("fSmearFramePb210BM2",    "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFramePb210S01M2      = new TH1D("fSmearFramePb210S01M2",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFramePb210S1M2       = new TH1D("fSmearFramePb210S1M2",   "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFramePb210S10M2      = new TH1D("fSmearFramePb210S10M2",  "Pb210",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFramePb210S100M2     = new TH1D("fSmearFramePb210S100M2", "Pb210",    dNBins, dMinEnergy, dMaxEnergy);

  fSmearFrameRa226BM2        = new TH1D("fSmearFrameRa226BM2",    "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRa226S01M2      = new TH1D("fSmearFrameRa226S01M2",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRa226S1M2       = new TH1D("fSmearFrameRa226S1M2",   "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRa226S10M2      = new TH1D("fSmearFrameRa226S10M2",  "Ra226",    dNBins, dMinEnergy, dMaxEnergy);
  fSmearFrameRa226S100M2     = new TH1D("fSmearFrameRa226S100M2", "Ra226",    dNBins, dMinEnergy, dMaxEnergy);



//////////////// Adaptive binned histograms
  fAdapFrameThM1    = new TH1D("fAdapFrameThM1",  "Frame",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapTShieldThM1  = new TH1D("fAdapTShieldThM1",  "TShield",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdap50mKThM1     = new TH1D("fAdap50mKThM1",  "50mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdap600mKThM1    = new TH1D("fAdap600mKThM1",  "600mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapIVCThM1      = new TH1D("fAdapIVCThM1",    "IVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapOVCThM1      = new TH1D("fAdapOVCThM1",    "OVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);

  fAdapFrameRaM1    = new TH1D("fAdapFrameRaM1",  "Frame",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapTShieldRaM1  = new TH1D("fAdapTShieldRaM1",  "TShield",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdap50mKRaM1     = new TH1D("fAdap50mKRaM1",  "50mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdap600mKRaM1    = new TH1D("fAdap600mKRaM1",  "600mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapIVCRaM1      = new TH1D("fAdapIVCRaM1",    "IVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapOVCRaM1      = new TH1D("fAdapOVCRaM1",    "OVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);

  fAdapFrameKM1     = new TH1D("fAdapFrameKM1",   "Frame",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapTShieldKM1   = new TH1D("fAdapTShieldKM1", "TShield",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdap50mKKM1      = new TH1D("fAdap50mKKM1",    "50mK",     dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdap600mKKM1     = new TH1D("fAdap600mKKM1",   "600mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapIVCKM1       = new TH1D("fAdapIVCKM1",     "IVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapOVCKM1       = new TH1D("fAdapOVCKM1",     "OVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);

  fAdapFrameCoM1    = new TH1D("fAdapFrameCoM1",  "Frame",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapTShieldCoM1  = new TH1D("fAdapTShieldCoM1", "TShield",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdap50mKCoM1     = new TH1D("fAdap50mKCoM1",  "50mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdap600mKCoM1    = new TH1D("fAdap600mKCoM1",  "600mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapIVCCoM1      = new TH1D("fAdapIVCCoM1",    "IVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapOVCCoM1      = new TH1D("fAdapOVCCoM1",    "OVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);

  fAdapTShieldMnM1  = new TH1D("fAdapTShieldMnM1",    "TShield",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapIVCMnM1      = new TH1D("fAdapIVCMnM1",    "IVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);

  fAdapCrystalBi2M1 = new TH1D("fAdapCrystalBi2M1",    "Crystal",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapFrameBi2M1   = new TH1D("fAdapFrameBi2M1",    "Frame",      dAdaptiveBinsM1, dAdaptiveArrayM1);


  fAdapNDBDM1       = new TH1D("fAdapNDBDM1",   "NDBD",       dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdap2NDBDM1      = new TH1D("fAdap2NDBDM1",  "2NDBD",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapBiM1         = new TH1D("fAdapBiM1",     "Bi",         dAdaptiveBinsM1, dAdaptiveArrayM1);

  fAdapCrystalPt190BM1        = new TH1D("fAdapCrystalPtBM1",       "Pt190",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapCrystalPt190S01M1      = new TH1D("fAdapCrystalPt190S01M1",  "Pt190",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapCrystalPt190S1M1       = new TH1D("fAdapCrystalPt190S1M1",   "Pt190",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapCrystalPt190S10M1      = new TH1D("fAdapCrystalPt190S10M1",  "Pt190",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapCrystalPt190S100M1     = new TH1D("fAdapCrystalPt190S100M1", "Pt190",    dAdaptiveBinsM1, dAdaptiveArrayM1);

  fAdapCrystalPb210BM1        = new TH1D("fAdapCrystalPb210BM1",    "Pb210",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapCrystalPb210S01M1      = new TH1D("fAdapCrystalPb210S01M1",  "Pb210",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapCrystalPb210S1M1       = new TH1D("fAdapCrystalPb210S1M1",   "Pb210",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapCrystalPb210S10M1      = new TH1D("fAdapCrystalPb210S10M1",  "Pb210",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapCrystalPb210S100M1     = new TH1D("fAdapCrystalPb210S100M1", "Pb210",    dAdaptiveBinsM1, dAdaptiveArrayM1);

  fAdapCrystalRa226BM1        = new TH1D("fAdapCrystalRa226BM1",    "Ra226",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapCrystalRa226S01M1      = new TH1D("fAdapCrystalRa226S01M1",  "Ra226",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapCrystalRa226S1M1       = new TH1D("fAdapCrystalRa226S1M1",   "Ra226",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapCrystalRa226S10M1      = new TH1D("fAdapCrystalRa226S10M1",  "Ra226",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapCrystalRa226S100M1     = new TH1D("fAdapCrystalRa226S100M1", "Ra226",    dAdaptiveBinsM1, dAdaptiveArrayM1);

  fAdapFramePb210BM1        = new TH1D("fAdapFramePb210BM1",    "Pb210",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapFramePb210S01M1      = new TH1D("fAdapFramePb210S01M1",  "Pb210",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapFramePb210S1M1       = new TH1D("fAdapFramePb210S1M1",   "Pb210",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapFramePb210S10M1      = new TH1D("fAdapFramePb210S10M1",  "Pb210",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapFramePb210S100M1     = new TH1D("fAdapFramePb210S100M1", "Pb210",    dAdaptiveBinsM1, dAdaptiveArrayM1);

  fAdapFrameRa226BM1        = new TH1D("fAdapFrameRa226BM1",    "Ra226",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapFrameRa226S01M1      = new TH1D("fAdapFrameRa226S01M1",  "Ra226",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapFrameRa226S1M1       = new TH1D("fAdapFrameRa226S1M1",   "Ra226",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapFrameRa226S10M1      = new TH1D("fAdapFrameRa226S10M1",  "Ra226",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapFrameRa226S100M1     = new TH1D("fAdapFrameRa226S100M1", "Ra226",    dAdaptiveBinsM1, dAdaptiveArrayM1);



  fAdapFrameThM2    = new TH1D("fAdapFrameThM2",  "Frame",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapTShieldThM2  = new TH1D("fAdapTShieldThM2",  "TShield",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdap50mKThM2     = new TH1D("fAdap50mKThM2",  "50mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdap600mKThM2    = new TH1D("fAdap600mKThM2",  "600mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapIVCThM2      = new TH1D("fAdapIVCThM2",    "IVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapOVCThM2      = new TH1D("fAdapOVCThM2",    "OVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);

  fAdapFrameRaM2    = new TH1D("fAdapFrameRaM2",  "Frame",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapTShieldRaM2  = new TH1D("fAdapTShieldRaM2",  "TShield",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdap50mKRaM2     = new TH1D("fAdap50mKRaM2",  "50mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdap600mKRaM2    = new TH1D("fAdap600mKRaM2",  "600mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapIVCRaM2      = new TH1D("fAdapIVCRaM2",    "IVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapOVCRaM2      = new TH1D("fAdapOVCRaM2",    "OVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);

  fAdapFrameKM2     = new TH1D("fAdapFrameKM2",   "Frame",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapTShieldKM2   = new TH1D("fAdapTShieldKM2", "TShield",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdap50mKKM2      = new TH1D("fAdap50mKKM2",    "50mK",     dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdap600mKKM2     = new TH1D("fAdap600mKKM2",   "600mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapIVCKM2       = new TH1D("fAdapIVCKM2",     "IVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapOVCKM2       = new TH1D("fAdapOVCKM2",     "OVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);

  fAdapFrameCoM2    = new TH1D("fAdapFrameCoM2",  "Frame",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapTShieldCoM2  = new TH1D("fAdapTShieldCoM2",  "TShield",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdap50mKCoM2     = new TH1D("fAdap50mKCoM2",  "50mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdap600mKCoM2    = new TH1D("fAdap600mKCoM2",  "600mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapIVCCoM2      = new TH1D("fAdapIVCCoM2",    "IVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapOVCCoM2      = new TH1D("fAdapOVCCoM2",    "OVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);

  fAdapTShieldMnM2  = new TH1D("fAdapTShieldMnM2",    "TShield",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapIVCMnM2      = new TH1D("fAdapIVCMnM2",    "IVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);

  fAdapCrystalBi2M2 = new TH1D("fAdapCrystalBi2M2",    "Crystal",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapFrameBi2M2   = new TH1D("fAdapFrameBi2M2",    "Frame",      dAdaptiveBinsM2, dAdaptiveArrayM2);

  fAdapNDBDM2       = new TH1D("fAdapNDBDM2",   "NDBD",       dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdap2NDBDM2      = new TH1D("fAdap2NDBDM2",  "2NDBD",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapBiM2         = new TH1D("fAdapBiM2",     "Bi",         dAdaptiveBinsM2, dAdaptiveArrayM2);


  fAdapCrystalPt190BM2        = new TH1D("fAdapCrystalPtBM2",       "Pt190",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapCrystalPt190S01M2      = new TH1D("fAdapCrystalPt190S01M2",  "Pt190",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapCrystalPt190S1M2       = new TH1D("fAdapCrystalPt190S1M2",   "Pt190",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapCrystalPt190S10M2      = new TH1D("fAdapCrystalPt190S10M2",  "Pt190",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapCrystalPt190S100M2     = new TH1D("fAdapCrystalPt190S100M2", "Pt190",    dAdaptiveBinsM2, dAdaptiveArrayM2);

  fAdapCrystalPb210BM2        = new TH1D("fAdapCrystalPb210BM2",    "Pb210",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapCrystalPb210S01M2      = new TH1D("fAdapCrystalPb210S01M2",  "Pb210",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapCrystalPb210S1M2       = new TH1D("fAdapCrystalPb210S1M2",   "Pb210",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapCrystalPb210S10M2      = new TH1D("fAdapCrystalPb210S10M2",  "Pb210",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapCrystalPb210S100M2     = new TH1D("fAdapCrystalPb210S100M2", "Pb210",    dAdaptiveBinsM2, dAdaptiveArrayM2);

  fAdapCrystalRa226BM2        = new TH1D("fAdapCrystalRa226BM2",    "Ra226",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapCrystalRa226S01M2      = new TH1D("fAdapCrystalRa226S01M2",  "Ra226",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapCrystalRa226S1M2       = new TH1D("fAdapCrystalRa226S1M2",   "Ra226",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapCrystalRa226S10M2      = new TH1D("fAdapCrystalRa226S10M2",  "Ra226",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapCrystalRa226S100M2     = new TH1D("fAdapCrystalRa226S100M2", "Ra226",    dAdaptiveBinsM2, dAdaptiveArrayM2);

  fAdapFramePb210BM2        = new TH1D("fAdapFramePb210BM2",    "Pb210",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapFramePb210S01M2      = new TH1D("fAdapFramePb210S01M2",  "Pb210",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapFramePb210S1M2       = new TH1D("fAdapFramePb210S1M2",   "Pb210",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapFramePb210S10M2      = new TH1D("fAdapFramePb210S10M2",  "Pb210",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapFramePb210S100M2     = new TH1D("fAdapFramePb210S100M2", "Pb210",    dAdaptiveBinsM2, dAdaptiveArrayM2);

  fAdapFrameRa226BM2        = new TH1D("fAdapFrameRa226BM2",    "Ra226",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapFrameRa226S01M2      = new TH1D("fAdapFrameRa226S01M2",  "Ra226",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapFrameRa226S1M2       = new TH1D("fAdapFrameRa226S1M2",   "Ra226",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapFrameRa226S10M2      = new TH1D("fAdapFrameRa226S10M2",  "Ra226",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapFrameRa226S100M2     = new TH1D("fAdapFrameRa226S100M2", "Ra226",    dAdaptiveBinsM2, dAdaptiveArrayM2);


  // Total Adaptive binning histograms M1
  fModelTotAdapM1      = new TH1D("fModelTotAdapM1",      "Total PDF M1", dAdaptiveBinsM1, dAdaptiveArrayM1);  
  fModelTotAdapThM1    = new TH1D("fModelTotAdapThM1",    "Total Th232",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapRaM1    = new TH1D("fModelTotAdapRaM1",    "Total Ra226",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapKM1     = new TH1D("fModelTotAdapKM1",     "Total K40",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapCoM1    = new TH1D("fModelTotAdapCoM1",    "Total Co60",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapMnM1    = new TH1D("fModelTotAdapMnM1",    "Total Mn54",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  fModelTotAdapNDBDM1  = new TH1D("fModelTotAdapNDBDM1",  "Total NDBD",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdap2NDBDM1 = new TH1D("fModelTotAdap2NDBDM1", "Total 2NDBD",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapBiM1    = new TH1D("fModelTotAdapBiM1",   "Total Bi207",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapBi2M1   = new TH1D("fModelTotAdapBi2M1",   "Total Bi210",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapPtM1    = new TH1D("fModelTotAdapPtM1",   "Total Pt190",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapPbM1    = new TH1D("fModelTotAdapPbM1",   "Total Pb210",   dAdaptiveBinsM1, dAdaptiveArrayM1);


  // Total Adaptive binning histograms M2
  fModelTotAdapM2      = new TH1D("fModelTotAdapM2",      "Total PDF M2", dAdaptiveBinsM2, dAdaptiveArrayM2);  
  fModelTotAdapThM2    = new TH1D("fModelTotAdapThM2",    "Total Th232",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapRaM2    = new TH1D("fModelTotAdapRaM2",    "Total Ra226",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapKM2     = new TH1D("fModelTotAdapKM2",     "Total K40",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapCoM2    = new TH1D("fModelTotAdapCoM2",    "Total Co60",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapMnM2    = new TH1D("fModelTotAdapMnM2",    "Total Mn54",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  fModelTotAdapNDBDM2  = new TH1D("fModelTotAdapNDBDM2",  "Total NDBD",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdap2NDBDM2 = new TH1D("fModelTotAdap2NDBDM2", "Total 2NDBD",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapBiM2    = new TH1D("fModelTotAdapBiM2",    "Total Bi207",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapBi2M2   = new TH1D("fModelTotAdapBi2M2",   "Total Bi210",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapPtM2    = new TH1D("fModelTotAdapPtM2",   "Total Pt190",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapPbM2    = new TH1D("fModelTotAdapPbM2",   "Total Pb210",   dAdaptiveBinsM2, dAdaptiveArrayM2);


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
  fParameters[19] = 0.;
  fParameters[20] = 0.;
  fParameters[21] = 0.;
  fParameters[22] = 0.;
  fParameters[23] = 0.;
  fParameters[24] = 0.;
  fParameters[25] = 0.;


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
  fParError[19] = 0.;
  fParError[20] = 0.;
  fParError[21] = 0.;
  fParError[22] = 0.;
  fParError[23] = 0.;
  fParError[24] = 0.;
  fParError[25] = 0.;


  ////////////////////////////// Histograms for accidental coincidence test
  fCorrectionM2     = new TH1D("fCorrectionM2",      "Correction Spectra",        dNBins, dMinEnergy, dMaxEnergy);  
  fCorrectionM2Tot     = new TH1D("fCorrectionM2Tot",      "Total Correction Spectra",        dNBins, dMinEnergy, dMaxEnergy);  
  fTotCorrection     = new TH1D("fTotCorrection",      "Total Correction Spectra",        dNBins, dMinEnergy, dMaxEnergy);  


  fModelTestM1      = new TH1D("fModelTestM1",    "Test",        dNBins, dMinEnergy, dMaxEnergy);  
  fModelTest2      = new TH1D("fModelTest2",      "Test",        dNBins, dMinEnergy, dMaxEnergy);  

  fFrameThM1     = new TH1D("fFrameThM1",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fTShieldThM1   = new TH1D("fTShieldThM1", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  f50mKThM1      = new TH1D("f50mKThM1",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  f600mKThM1     = new TH1D("f600mKThM1",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fIVCThM1       = new TH1D("fIVCThM1",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fOVCThM1       = new TH1D("fOVCThM1",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fFrameRaM1     = new TH1D("fFrameRaM1",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fTShieldRaM1   = new TH1D("fTShieldRaM1", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  f50mKRaM1      = new TH1D("f50mKRaM1",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  f600mKRaM1     = new TH1D("f600mKRaM1",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fIVCRaM1       = new TH1D("fIVCRaM1",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fOVCRaM1       = new TH1D("fOVCRaM1",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fFrameKM1     = new TH1D("fFrameKM1",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fTShieldKM1   = new TH1D("fTShieldKM1", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  f50mKKM1      = new TH1D("f50mKKM1",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  f600mKKM1     = new TH1D("f600mKKM1",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fIVCKM1       = new TH1D("fIVCKM1",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fOVCKM1       = new TH1D("fOVCKM1",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fFrameCoM1     = new TH1D("fFrameCoM1",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fTShieldCoM1   = new TH1D("fTShieldCoM1", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  f50mKCoM1      = new TH1D("f50mKCoM1",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  f600mKCoM1     = new TH1D("f600mKCoM1",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fIVCCoM1       = new TH1D("fIVCCoM1",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fOVCCoM1       = new TH1D("fOVCCoM1",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);


  fTShieldMnM1     = new TH1D("fTShieldMnM1", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  fIVCMnM1         = new TH1D("fIVCMnM1",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);

  fCrystalBi2M1    = new TH1D("fCrystalBi2M1",  "Crystal",    dNBins, dMinEnergy, dMaxEnergy);
  fFrameBi2M1      = new TH1D("fFrameBi2M1",    "Frame",      dNBins, dMinEnergy, dMaxEnergy);

  fNDBDM1       = new TH1D("fNDBDM1",   "NDBD",       dNBins, dMinEnergy, dMaxEnergy);
  f2NDBDM1      = new TH1D("f2NDBDM1",  "2NDBD",      dNBins, dMinEnergy, dMaxEnergy);
  fBiM1         = new TH1D("fBiM1",     "Bi",         dNBins, dMinEnergy, dMaxEnergy);


  fFrameThM2     = new TH1D("fFrameThM2",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fTShieldThM2   = new TH1D("fTShieldThM2", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  f50mKThM2      = new TH1D("f50mKThM2",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  f600mKThM2     = new TH1D("f600mKThM2",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fIVCThM2       = new TH1D("fIVCThM2",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fOVCThM2       = new TH1D("fOVCThM2",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fFrameRaM2     = new TH1D("fFrameRaM2",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fTShieldRaM2   = new TH1D("fTShieldRaM2", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  f50mKRaM2      = new TH1D("f50mKRaM2",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  f600mKRaM2     = new TH1D("f600mKRaM2",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fIVCRaM2       = new TH1D("fIVCRaM2",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fOVCRaM2       = new TH1D("fOVCRaM2",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fFrameKM2     = new TH1D("fFrameKM2",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fTShieldKM2   = new TH1D("fTShieldKM2", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  f50mKKM2      = new TH1D("f50mKKM2",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  f600mKKM2     = new TH1D("f600mKKM2",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fIVCKM2       = new TH1D("fIVCKM2",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fOVCKM2       = new TH1D("fOVCKM2",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  fFrameCoM2     = new TH1D("fFrameCoM2",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  fTShieldCoM2   = new TH1D("fTShieldCoM2", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  f50mKCoM2      = new TH1D("f50mKCoM2",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  f600mKCoM2     = new TH1D("f600mKCoM2",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  fIVCCoM2       = new TH1D("fIVCCoM2",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  fOVCCoM2       = new TH1D("fOVCCoM2",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);


  fTShieldMnM2   = new TH1D("fTShieldMnM2", "TShield",      dNBins, dMinEnergy, dMaxEnergy);
  fIVCMnM2       = new TH1D("fIVCMnM2",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);

  fCrystalBi2M2    = new TH1D("fCrystalBi2M2",  "Crystal",      dNBins, dMinEnergy, dMaxEnergy);
  fFrameBi2M2      = new TH1D("fFrameBi2M2",    "Frame",      dNBins, dMinEnergy, dMaxEnergy);

  fNDBDM2       = new TH1D("fNDBDM2",   "NDBD",       dNBins, dMinEnergy, dMaxEnergy);
  f2NDBDM2      = new TH1D("f2NDBDM2",  "2NDBD",      dNBins, dMinEnergy, dMaxEnergy);
  fBiM2         = new TH1D("fBiM2",     "Bi",         dNBins, dMinEnergy, dMaxEnergy);


  fFileCorrection = new TFile(Form("MCCorrection-%dkeV.root", dBinSize));
  fCorrectionM2 = (TH1D*)fFileCorrection->Get("fModelTotM1");
  // fCorrectionM2Tot = (TH1D*)fFileCorrection->Get("fCorrectionM2Tot");

}
  
// Needs to be updated
TBackgroundModel::~TBackgroundModel()
{
	delete	fDataHistoTot;
	delete	fDataHistoM1;
	delete	fDataHistoM2;
	delete	fToyData;

	delete 	fModelFrameThM1;
	delete	fModelTShieldThM1;
	delete	fModel50mKThM1;
	delete	fModel600mKThM1;
	delete 	fModelIVCThM1;
	delete	fModelOVCThM1;

  delete  fModelFrameThM2;
  delete  fModelTShieldThM2;
  delete  fModel50mKThM2;
  delete  fModel600mKThM2;
  delete  fModelIVCThM2;
  delete  fModelOVCThM2;


	delete	fModelFrameRaM1;
	delete 	fModelTShieldRaM1;
	delete	fModel50mKRaM1;
	delete	fModel600mKRaM1;
	delete 	fModelIVCRaM1;
	delete	fModelOVCRaM1;

  delete  fModelFrameRaM2;
  delete  fModelTShieldRaM2;
  delete  fModel50mKRaM2;
  delete  fModel600mKRaM2;
  delete  fModelIVCRaM2;
  delete  fModelOVCRaM2;

	delete	fModelFrameKM1;
	delete 	fModelTShieldKM1;
	delete	fModel50mKKM1;
	delete	fModel600mKKM1;
	delete 	fModelIVCKM1;
	delete	fModelOVCKM1;

  delete  fModelFrameKM2;
  delete  fModelTShieldKM2;
  delete  fModel50mKKM2;
  delete  fModel600mKKM2;
  delete  fModelIVCKM2;
  delete  fModelOVCKM2;

	delete	fModelFrameCoM1;
	delete 	fModelTShieldCoM1;
	delete	fModel50mKCoM1;
	delete	fModel600mKCoM1;
	delete 	fModelIVCCoM1;
	delete	fModelOVCCoM1;

  delete  fModelFrameCoM2;
  delete  fModelTShieldCoM2;
  delete  fModel50mKCoM2;
  delete  fModel600mKCoM2;
  delete  fModelIVCCoM2;
  delete  fModelOVCCoM2;

	delete 	fModelTotThM1;
	delete	fModelTotRaM1;
	delete	fModelTotKM1;
	delete	fModelTotCoM1;
  delete  fModelTotBiM1;
  delete  fModelTotNDBDM1;
  delete  fModelTot2NDBDM1;

  delete  fModelTotThM2;
  delete  fModelTotRaM2;
  delete  fModelTotKM2;
  delete  fModelTotCoM2;
  delete  fModelTotBiM2;
  delete  fModelTotNDBDM2;
  delete  fModelTot2NDBDM2;

}

vector<double> TBackgroundModel::AdaptiveBinning(TH1D *h1)
{

  vector<double> dBinArrayThing;

  double dDummy = 0;
  double dDummyFill = 0;
  int j = 0;

  for(int i = 1; i < 25; i++)
  {
    dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(i));
  }

  for(int i = 25; i < dNBins; i++)
  {
    dDummy = h1->GetBinContent(i);
    dDummyFill += dDummy;

    if(dDummyFill >= 50)
    {
      dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(i-j));
      dDummyFill = 0;
      j = 0;
    }
    else if(i == dNBins-1) // for the very end if it doesn't reach 50 events (which it won't)
    {
      dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(i-j));
    }
    else 
    {
      j++;
    }
  }

return dBinArrayThing;
}


TH1D *TBackgroundModel::CalculateResiduals(TH1D *h1, TH1D *h2, TH1D *hResid)
{

	// Clone histograms for rebinning
	TH1D 	*hCloneBkg 		= (TH1D*)h1->Clone("hCloneBkg");
	TH1D 	*hCloneMC		= (TH1D*)h2->Clone("hCloneMC");

  TH1D  *hOut       = new TH1D("hOut", "Fit Residuals", dNBins, dMinEnergy, dMaxEnergy);


	// Variables used in Residual calculations
	double dResidualX, dResidualY, dResidualXErr = 0, dResidualYErr;

	// Residual plot and distribution
	for (int j = dFitMin/dBinSize+1; j <= dFitMax/dBinSize; j++)
	{
		dResidualX 		= hCloneBkg->GetBinCenter(j);
		dResidualY 		= (hCloneBkg->GetBinContent(j) - hCloneMC->GetBinContent(j)) /
							TMath::Sqrt(hCloneBkg->GetBinContent(j)); // Sqrt(MC + data) = sigma for poisson distribution

		// g1->SetPoint(j, dResidualX, dResidualY);
		hOut->SetBinContent(j, dResidualY);
		hOut->SetBinError(j, 0.1);
    hResid->Fill(dResidualY);
	}


	return hOut;
}

TH1D *TBackgroundModel::CalculateResidualsAdaptive(TH1D *h1, TH1D *h2, TH1D *hResid, int binMin, int binMax)
{

  // Clone histograms for rebinning
  TH1D  *hCloneBkg    = (TH1D*)h1->Clone("hCloneBkg");
  TH1D  *hCloneMC   = (TH1D*)h2->Clone("hCloneMC");

  TH1D  *hOut       = new TH1D("hOut", "Fit Residuals", );


  // Variables used in Residual calculations
  double dResidualX, dResidualY, dResidualXErr = 0, dResidualYErr;

  // Residual plot and distribution
  for (int j = binMin ; j <= binMax ; j++)
  {
    dResidualX    = hCloneBkg->GetBinCenter(j);
    dResidualY    = (hCloneBkg->GetBinContent(j) - hCloneMC->GetBinContent(j)) /
              TMath::Sqrt(hCloneBkg->GetBinContent(j)); // Sqrt(MC + data) = sigma for poisson distribution

    // g1->SetPoint(j, dResidualX, dResidualY);
    hOut->SetBinContent(j, dResidualY);
    hOut->SetBinError(j, 0.1);
    hResid->Fill(dResidualY);
  }


  return hOut;
}

bool TBackgroundModel::DoTheFit()
{
  if(!bToyFit)
  {
    Initialize();
  }
	gStyle->SetOptStat(0);
   // This method actually sets up minuit and does the fit


   // TMinuit minuit(14); //initialize minuit, n is the number of parameters
   TMinuit minuit(26); // for more parameters

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
/*
   minuit.DefineParameter(0, "Close Th",  20., 10.0, 0., 100000);
   // minuit.DefineParameter(0, "Close Th",  100., 1.0, 0., 100000);
   // minuit.DefineParameter(1, "Far Th",	 	100., 50.0, 0., 100000);
   minuit.DefineParameter(1, "Far Th",    130000., 10.0, 0., 500000);
   // minuit.DefineParameter(2, "Close Ra",  37800., 10.0, 0., 80000);   
   minuit.DefineParameter(2, "Close Ra",  22000., 10.0, 0., 80000);   
   // minuit.DefineParameter(3, "Far Ra",    55000., 10.0, 0., 80000);
   minuit.DefineParameter(3, "Far Ra",    100., 10.0, 0., 80000);
   // minuit.DefineParameter(4, "Close K", 	0., 100.0, 0., 500000);
   minuit.DefineParameter(4, "Close K",   100., 1.0, 0., 500000);
   // minuit.DefineParameter(5, "Far K",     0., 100.0, 0., 500000);
   minuit.DefineParameter(5, "Far K", 		90000., 10.0, 0., 500000);
   // minuit.DefineParameter(6, "Close Co", 	0., 10.0, 0., 80000); 
   minuit.DefineParameter(6, "Close Co",  30., 1.0, 0., 80000);    
   minuit.DefineParameter(7, "Far Co",    30000, 10.0, 0., 500000);  
   // minuit.DefineParameter(7, "Far Co",	 	0., 100.0, 0., 50000);  
   // minuit.DefineParameter(8, "Resolution",	6., 1, 3, 10);  
   minuit.DefineParameter(8, "2NDBD",    33394., 10.0, 0., 100000);     
   // minuit.DefineParameter(8, "2NDBD",    0., 10.0, 0., 100000);  
   minuit.DefineParameter(9, "NDBD",       87.4, 1.0, 0., 500);     
   minuit.DefineParameter(10, "Lead Bi",	 	10000., 10.0, 0., 100000);  
   // minuit.DefineParameter(10, "Lead Bi",    0., 100.0, 0., 100000);  
   minuit.DefineParameter(11, "TShield Mn",    100., 10.0, 0., 100000);  
   minuit.DefineParameter(12, "IVC Mn",    100., 10.0, 0., 100000);  
   // minuit.DefineParameter(13, "IVC Mn",    0., 10.0, 0., 100000);  
*/

   // Close Th and Close Ra now split into its various sections, far Th and Ra still the same
   // This step after previous step full fit converges, just to see if any differences show up
   ////////////////////////////////////////////////
   // Using more parameters
   ////////////////////////////////////////////////

   minuit.DefineParameter(0, "Frame Th",  0., 10.0, 0., 100000);
   minuit.DefineParameter(1, "TShield Th",    0., 10.0, 0., 500000);
   minuit.DefineParameter(2, "Frame Ra",  0., 10.0, 0., 80000);   
   minuit.DefineParameter(3, "TShield Ra", 0., 10.0, 0., 80000);
   minuit.DefineParameter(4, "Close K",   21732.5, 1.0, 0., 500000);
   minuit.DefineParameter(5, "Far K",     3463.03, 10.0, 0., 500000);
   minuit.DefineParameter(6, "Frame Co",  1828.31, 1.0, 0., 80000);    
   minuit.DefineParameter(7, "TShield Co",    0, 10.0, 0., 80000);  
   minuit.DefineParameter(8, "2NDBD",    33394., 10.0, 0., 100000);        
   // minuit.DefineParameter(8, "2NDBD",    53000.6., 10.0, 0., 100000);   
   minuit.DefineParameter(9, "NDBD",       92.45, 1.0, 0., 500);     
   minuit.DefineParameter(10, "Lead Bi",    7723.64, 10.0, 0., 100000);  
   minuit.DefineParameter(11, "TShield Mn",    0., 10.0, 0., 100000);  
   minuit.DefineParameter(12, "IVC Mn",    5435.33, 10.0, 0., 100000);  
   minuit.DefineParameter(13, "50mK Th",    0., 10.0, 0., 500000);
   minuit.DefineParameter(14, "600mK Th",    37952.8, 10.0, 0., 500000);
   minuit.DefineParameter(15, "IVC Th",    42461.1, 10.0, 0., 500000);
   minuit.DefineParameter(16, "OVC Th",    16575.1, 10.0, 0., 500000);
   minuit.DefineParameter(17, "50mK Ra", 0., 10.0, 0., 80000);
   minuit.DefineParameter(18, "600mK Ra", 9666.94, 10.0, 0., 80000);
   minuit.DefineParameter(19, "IVC Ra",    0., 10.0, 0., 500000);
   minuit.DefineParameter(20, "OVC Ra",    106770., 10.0, 0., 500000);
   minuit.DefineParameter(21, "50mK Co",  0., 1.0, 0., 80000);    
   minuit.DefineParameter(22, "600mK Co",    0., 10.0, 0., 80000);
   minuit.DefineParameter(23, "IVC Co",    0, 10.0, 0., 500000);  
   minuit.DefineParameter(24, "OVC Co",    20815.6, 10.0, 0., 500000);  
   // minuit.DefineParameter(25, "Constant",    500000, 1, 0., 1000000);  
   minuit.DefineParameter(25, "Bi-210",    110000, 1, 0., 1000000);  

   



   // Fix parameters here
   // minuit.FixParameter(0); // Close Th
   // minuit.FixParameter(1); // Far Th
   // minuit.FixParameter(2); // Close Ra
   // minuit.FixParameter(3); // Far Ra
   // minuit.FixParameter(4); // Close K
   // minuit.FixParameter(5); // Far K
   // minuit.FixParameter(6); // Close Co
   // minuit.FixParameter(7); // Far Co
   minuit.FixParameter(8); // 2NDBD
   // minuit.FixParameter(9); // NDBD
   // minuit.FixParameter(10); // Bi207
   // minuit.FixParameter(11); // Close Mn
   // minuit.FixParameter(12); // Far Mn
   // minuit.FixParameter(13); // 
   // minuit.FixParameter(14); // 
   // minuit.FixParameter(15); // 
   // minuit.FixParameter(16); // 
   // minuit.FixParameter(17); // 
   // minuit.FixParameter(18); // 
   // minuit.FixParameter(19); // 
   // minuit.FixParameter(20); // 
   // minuit.FixParameter(21); // 
   // minuit.FixParameter(22); // 
   // minuit.FixParameter(23); // 
   // minuit.FixParameter(24); // 
   // minuit.FixParameter(25); // 

  // Number of Parameters! (for Chi-squared/NDF calculation)
  int dNumParameters = 26;




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
  minuit.GetParameter(11,  fParameters[11],   fParError[11]);
  minuit.GetParameter(12,  fParameters[12],   fParError[12]);
  minuit.GetParameter(13,  fParameters[13],   fParError[13]);
  minuit.GetParameter(14,  fParameters[14],   fParError[14]);
  minuit.GetParameter(15,  fParameters[15],   fParError[15]);
  minuit.GetParameter(16,  fParameters[16],   fParError[16]);
  minuit.GetParameter(17,  fParameters[17],   fParError[17]);
  minuit.GetParameter(18,  fParameters[18],   fParError[18]);
  minuit.GetParameter(19,  fParameters[19],   fParError[19]);
  minuit.GetParameter(20,  fParameters[20],   fParError[20]);
  minuit.GetParameter(21,  fParameters[21],   fParError[21]);
  minuit.GetParameter(22,  fParameters[22],   fParError[22]);
  minuit.GetParameter(23,  fParameters[23],   fParError[23]);
  minuit.GetParameter(24,  fParameters[24],   fParError[24]);
  minuit.GetParameter(25,  fParameters[25],   fParError[25]);
	UpdateModel();
	
	cout << "At the end; ChiSq/NDF = " << GetChiSquare()/(2*(dFitMax-dFitMin)/dBinSize - dNumParameters) << endl; // for M1 and M2
  // cout << "At the end; ChiSq/NDF = " << GetChiSquare()/((dFitMax-dFitMin)/dBinSize - dNumParameters) << endl;  // for M1
  cout << "Total number of calls = " << dNumCalls << endl;


  ///////////////////////////////////////////
  //// Few Parameters - Pre-smeared
  ///////////////////////////////////////////
  /// Add Modeled Histograms after chi-squared minimization
/*
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
  fModelTot2NDBDM1->Add(fModel2NDBDM1,  fParameters[11]);
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
  fModelTot2NDBDM2->Add(fModel2NDBDM2,  fParameters[11]);
  fModelTotBiM2->Add(fModelBiM2,        fParameters[10]);
*/


  ///////////////////////////////////////////
  //// Few Parameters - no smear
  ///////////////////////////////////////////
  /// Add Smeared Histograms after chi-squared minimization
/*
  // M1 Parameters
  fModelTotThM1->Add(fSmearFrameThM1,   fParameters[0]);
  fModelTotThM1->Add(fSmearTShieldThM1, fParameters[0]);
  fModelTotThM1->Add(fSmear50mKThM1,    fParameters[0]);
  fModelTotThM1->Add(fSmear600mKThM1,   fParameters[0]);
  fModelTotThM1->Add(fSmearIVCThM1,     fParameters[1]);
  fModelTotThM1->Add(fSmearOVCThM1,     fParameters[1]);

  fModelTotRaM1->Add(fSmearFrameRaM1,   fParameters[2]);
  fModelTotRaM1->Add(fSmearTShieldRaM1, fParameters[2]);
  fModelTotRaM1->Add(fSmear50mKRaM1,    fParameters[2]);
  fModelTotRaM1->Add(fSmear600mKRaM1,   fParameters[2]);
  fModelTotRaM1->Add(fSmearIVCRaM1,     fParameters[3]);
  // fModelTotRaM1->Add(fSmearOVCRaM1,     fParameters[3]);

  // fModelTotKM1->Add(fSmearFrameKM1,     fParameters[4]);
  fModelTotKM1->Add(fSmearTShieldKM1,   fParameters[4]);
  fModelTotKM1->Add(fSmear50mKKM1,      fParameters[4]);
  fModelTotKM1->Add(fSmear600mKKM1,     fParameters[4]);
  fModelTotKM1->Add(fSmearIVCKM1,       fParameters[5]);
  // fModelTotKM1->Add(fSmearOVCKM1,       fParameters[5]);

  // fModelTotCoM1->Add(fSmearFrameCoM1,   fParameters[6]);
  fModelTotCoM1->Add(fSmearTShieldCoM1, fParameters[6]);
  fModelTotCoM1->Add(fSmear50mKCoM1,    fParameters[6]);
  fModelTotCoM1->Add(fSmear600mKCoM1,   fParameters[6]);
  fModelTotCoM1->Add(fSmearIVCCoM1,     fParameters[7]);
  // fModelTotCoM1->Add(fSmearOVCCoM1,     fParameters[7]);

  fModelTotMnM1->Add(fSmearTShieldMnM1, fParameters[11]);
  fModelTotMnM1->Add(fSmearIVCMnM1,     fParameters[12]);

  fModelTotNDBDM1->Add(fSmearNDBDM1,    fParameters[9]);
  fModelTot2NDBDM1->Add(fSmear2NDBDM1,  fParameters[8]);
  fModelTotBiM1->Add(fSmearBiM1,        fParameters[10]);


  // M2 Parameters
  fModelTotThM2->Add(fSmearFrameThM2,   fParameters[0]);
  fModelTotThM2->Add(fSmearTShieldThM2, fParameters[0]);
  fModelTotThM2->Add(fSmear50mKThM2,    fParameters[0]);
  fModelTotThM2->Add(fSmear600mKThM2,   fParameters[0]);
  fModelTotThM2->Add(fSmearIVCThM2,     fParameters[1]);
  fModelTotThM2->Add(fSmearOVCThM2,     fParameters[1]);

  fModelTotRaM2->Add(fSmearFrameRaM2,   fParameters[2]);
  fModelTotRaM2->Add(fSmearTShieldRaM2, fParameters[2]);
  fModelTotRaM2->Add(fSmear50mKRaM2,    fParameters[2]);
  fModelTotRaM2->Add(fSmear600mKRaM2,   fParameters[2]);
  fModelTotRaM2->Add(fSmearIVCRaM2,     fParameters[3]);
  // fModelTotRaM2->Add(fSmearOVCRaM2,     fParameters[3]);

  // fModelTotKM2->Add(fSmearFrameKM2,     fParameters[4]);
  fModelTotKM2->Add(fSmearTShieldKM2,   fParameters[4]);
  fModelTotKM2->Add(fSmear50mKKM2,      fParameters[4]);
  fModelTotKM2->Add(fSmear600mKKM2,     fParameters[4]);
  fModelTotKM2->Add(fSmearIVCKM2,       fParameters[5]);
  // fModelTotKM2->Add(fSmearOVCKM2,       fParameters[5]);

  // fModelTotCoM2->Add(fSmearFrameCoM2,   fParameters[6]);
  fModelTotCoM2->Add(fSmearTShieldCoM2, fParameters[6]);
  fModelTotCoM2->Add(fSmear50mKCoM2,    fParameters[6]);
  fModelTotCoM2->Add(fSmear600mKCoM2,   fParameters[6]);
  fModelTotCoM2->Add(fSmearIVCCoM2,     fParameters[7]);
  // fModelTotCoM2->Add(fSmearOVCCoM2,     fParameters[7]);

  fModelTotMnM2->Add(fSmearTShieldMnM2, fParameters[11]);
  fModelTotMnM2->Add(fSmearIVCMnM2,     fParameters[12]);

  fModelTotNDBDM2->Add(fSmearNDBDM2,    fParameters[9]);
  fModelTot2NDBDM2->Add(fSmear2NDBDM2,  fParameters[8]);
  fModelTotBiM2->Add(fSmearBiM2,      fParameters[10]);
*/

  ///////////////////////////////////////////
  //// Many Parameters
  ///////////////////////////////////////////
  /// Use only after previous step converges!
  // 

  // M1 Parameters
  fModelTotThM1->Add(fSmearFrameThM1,   fParameters[0]);
  fModelTotThM1->Add(fSmearTShieldThM1, fParameters[1]);
  fModelTotThM1->Add(fSmear50mKThM1,    fParameters[13]);
  fModelTotThM1->Add(fSmear600mKThM1,   fParameters[14]);
  fModelTotThM1->Add(fSmearIVCThM1,     fParameters[15]);
  fModelTotThM1->Add(fSmearOVCThM1,     fParameters[16]);

  fModelTotRaM1->Add(fSmearFrameRaM1,   fParameters[2]);
  fModelTotRaM1->Add(fSmearTShieldRaM1, fParameters[3]);
  fModelTotRaM1->Add(fSmear50mKRaM1,    fParameters[17]);
  fModelTotRaM1->Add(fSmear600mKRaM1,   fParameters[18]);
  fModelTotRaM1->Add(fSmearIVCRaM1,     fParameters[19]);
  fModelTotRaM1->Add(fSmearOVCRaM1,     fParameters[20]);

  fModelTotKM1->Add(fSmearFrameKM1,     fParameters[4]);
  fModelTotKM1->Add(fSmearTShieldKM1,   fParameters[4]);
  fModelTotKM1->Add(fSmear50mKKM1,      fParameters[4]);
  fModelTotKM1->Add(fSmear600mKKM1,     fParameters[4]);
  fModelTotKM1->Add(fSmearIVCKM1,       fParameters[5]);
  fModelTotKM1->Add(fSmearOVCKM1,       fParameters[5]);

  fModelTotCoM1->Add(fSmearFrameCoM1,   fParameters[6]);
  fModelTotCoM1->Add(fSmearTShieldCoM1, fParameters[7]);
  fModelTotCoM1->Add(fSmear50mKCoM1,    fParameters[21]);
  fModelTotCoM1->Add(fSmear600mKCoM1,   fParameters[22]);
  fModelTotCoM1->Add(fSmearIVCCoM1,     fParameters[23]);
  fModelTotCoM1->Add(fSmearOVCCoM1,     fParameters[24]);

  fModelTotMnM1->Add(fSmearTShieldMnM1, fParameters[11]);
  fModelTotMnM1->Add(fSmearIVCMnM1,     fParameters[12]);


  fModelTotNDBDM1->Add(fSmearNDBDM1,    fParameters[9]);
  fModelTot2NDBDM1->Add(fSmear2NDBDM1,  fParameters[8]);
  fModelTotBiM1->Add(fSmearBiM1,        fParameters[10]);

  fModelTotPbM1->Add(fSmearCrystalBi2M1, fParameters[25]);


  // M2 Parameters
  fModelTotThM2->Add(fSmearFrameThM2,   fParameters[0]);
  fModelTotThM2->Add(fSmearTShieldThM2, fParameters[1]);
  fModelTotThM2->Add(fSmear50mKThM2,    fParameters[13]);
  fModelTotThM2->Add(fSmear600mKThM2,   fParameters[14]);
  fModelTotThM2->Add(fSmearIVCThM2,     fParameters[15]);
  fModelTotThM2->Add(fSmearOVCThM2,     fParameters[16]);

  fModelTotRaM2->Add(fSmearFrameRaM2,   fParameters[2]);
  fModelTotRaM2->Add(fSmearTShieldRaM2, fParameters[3]);
  fModelTotRaM2->Add(fSmear50mKRaM2,    fParameters[17]);
  fModelTotRaM2->Add(fSmear600mKRaM2,   fParameters[18]);
  fModelTotRaM2->Add(fSmearIVCRaM2,     fParameters[19]);
  fModelTotRaM2->Add(fSmearOVCRaM2,     fParameters[20]);

  fModelTotKM2->Add(fSmearFrameKM2,     fParameters[4]);
  fModelTotKM2->Add(fSmearTShieldKM2,   fParameters[4]);
  fModelTotKM2->Add(fSmear50mKKM2,      fParameters[4]);
  fModelTotKM2->Add(fSmear600mKKM2,     fParameters[4]);
  fModelTotKM2->Add(fSmearIVCKM2,       fParameters[5]);
  fModelTotKM2->Add(fSmearOVCKM2,       fParameters[5]);

  fModelTotCoM2->Add(fSmearFrameCoM2,   fParameters[6]);
  fModelTotCoM2->Add(fSmearTShieldCoM2, fParameters[7]);
  fModelTotCoM2->Add(fSmear50mKCoM2,    fParameters[21]);
  fModelTotCoM2->Add(fSmear600mKCoM2,   fParameters[22]);
  fModelTotCoM2->Add(fSmearIVCCoM2,     fParameters[23]);
  fModelTotCoM2->Add(fSmearOVCCoM2,     fParameters[24]);

  fModelTotMnM2->Add(fSmearTShieldMnM2, fParameters[11]);
  fModelTotMnM2->Add(fSmearIVCMnM2,     fParameters[12]);

  fModelTotNDBDM2->Add(fSmearNDBDM2,    fParameters[9]);
  fModelTot2NDBDM2->Add(fSmear2NDBDM2,  fParameters[8]);
  fModelTotBiM2->Add(fSmearBiM2,      fParameters[10]);

  fModelTotPbM2->Add( fSmearCrystalBi2M2, fParameters[25]);

  ////////// Only for testing
  // Correction for M2 spectra, it's the M1 spectra but scaled down by N_M1*1-Exp(2*R*T)
  // 2 because of double counting in M2 spectrum...
  fTotCorrection->Add(fCorrectionM2, 180197*(1-TMath::Exp(-2*0.05*0.1)));



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
    fModelTotM1->SetLineWidth(1);
    fModelTotM1->Draw("SAME");

  
//////////////////////////
    /*
    fModelTest2->SetLineColor(4);
    fModelTest2->SetLineWidth(1);
    fModelTest2->Draw("SAME");

    TLegend *legfit1 = new TLegend(0.67,0.87,0.97,0.97);
    legfit1->AddEntry(fModelTotM1, "Standard PDF", "l");
    legfit1->AddEntry(fModelTest2, "Accidental Coincidence corrected PDF", "l");
    legfit1->Draw();
    */
//////////////////////////////////////////

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
    fModelTot2NDBDM1->SetLineColor(46);
    fModelTot2NDBDM1->SetLineStyle(2);
    fModelTotBiM1->SetLineColor(5);
    fModelTotBiM1->SetLineStyle(2);
    fModelTotMnM1->SetLineColor(40);
    fModelTotMnM1->SetLineStyle(2);


    fModelTotPbM1->SetLineStyle(2);
    fModelTotPbM1->SetLineColor(38);

    fModelTotThM1->Draw("SAME");
    fModelTotRaM1->Draw("SAME");
    fModelTotKM1->Draw("SAME");
    fModelTotCoM1->Draw("SAME");
    fModelTotNDBDM1->Draw("SAME");
    fModelTot2NDBDM1->Draw("SAME");
    fModelTotBiM1->Draw("SAME");
    fModelTotMnM1->Draw("SAME");

    fModelTotPbM1->Draw("SAME");

/*
    // Few Parameters
    TPaveText *pt1 = new TPaveText(0.35,0.78,0.70,0.98,"NB NDC");
    pt1->AddText(Form("Fit Range (M1): %.0f to %.0f keV -- #chi^{2}/NDF: %0.3f", dFitMin, dFitMax, (GetChiSquare()/(2*(dFitMax-dFitMin)/dBinSize - dNumParameters)) ));
    pt1->AddText(Form("Close Th: %0.2E#pm%0.2E --- Far Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
    // pt1->AddText(Form("Surface Th: %0.2E#pm%0.2E --- Far Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
    pt1->AddText(Form("Close Ra: %0.2E#pm%0.2E --- Far Ra: %0.2E#pm%0.2E", fParameters[2], fParError[2], fParameters[3], fParError[3] ));
    pt1->AddText(Form("Close K: %0.2E#pm%0.2E --- Far K: %0.2E#pm%0.2E", fParameters[4], fParError[4], fParameters[5], fParError[5] ));
    pt1->AddText(Form("Close Co: %0.2E#pm%0.2E --- Far Co: %0.2E#pm%0.2E", fParameters[6], fParError[6], fParameters[7], fParError[7] ));
    pt1->AddText(Form("Bi-207: %0.2E#pm%0.2E --- NDBD: %0.2E#pm%0.2E", fParameters[10], fParError[10], fParameters[9], fParError[9] ));
    pt1->AddText(Form("Close Mn-54: %0.2E#pm%0.2E --- Far Mn-54: %0.2E#pm%0.2E", fParameters[11], fParError[11], fParameters[12], fParError[12] ));
    pt1->AddText(Form("2NDBD: %0.2E#pm%0.2E", fParameters[8], fParError[8] ));
*/


    // Many Parameters
    TPaveText *pt1 = new TPaveText(0.35,0.77,0.70,0.99,"NB NDC");
    pt1->AddText(Form("Fit Range (M1): %.0f to %.0f keV -- #chi^{2}/NDF: %0.3f", dFitMin, dFitMax, (GetChiSquare()/(2*(dFitMax-dFitMin)/dBinSize - dNumParameters)) ));
    pt1->AddText(Form("Frame Th: %0.2E#pm%0.2E --- TShield Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
    pt1->AddText(Form("50mK Th: %0.2E#pm%0.2E --- 600mK Th: %0.2E#pm%0.2E", fParameters[13], fParError[13], fParameters[14], fParError[14] ));
    pt1->AddText(Form("IVC Th: %0.2E#pm%0.2E --- OVC Th: %0.2E#pm%0.2E", fParameters[15], fParError[15], fParameters[16], fParError[16] ));
    pt1->AddText(Form("Frame Ra: %0.2E#pm%0.2E --- TShield Ra: %0.2E#pm%0.2E", fParameters[2], fParError[2], fParameters[3], fParError[3] ));
    pt1->AddText(Form("50mK Ra: %0.2E#pm%0.2E --- 600mK Ra: %0.2E#pm%0.2E", fParameters[17], fParError[17], fParameters[18], fParError[18] ));
    pt1->AddText(Form("IVC Ra: %0.2E#pm%0.2E --- OVC Ra: %0.2E#pm%0.2E", fParameters[19], fParError[19], fParameters[20], fParError[20] ));
    pt1->AddText(Form("Close K: %0.2E#pm%0.2E --- Far K: %0.2E#pm%0.2E", fParameters[4], fParError[4], fParameters[5], fParError[5] ));
    pt1->AddText(Form("Close Co: %0.2E#pm%0.2E --- Far Co: %0.2E#pm%0.2E", fParameters[6], fParError[6], fParameters[7], fParError[7] ));
    pt1->AddText(Form("Bi-207: %0.2E#pm%0.2E --- NDBD: %0.2E#pm%0.2E", fParameters[10], fParError[10], fParameters[9], fParError[9] ));
    pt1->AddText(Form("Close Mn-54: %0.2E#pm%0.2E --- Far Mn-54: %0.2E#pm%0.2E", fParameters[11], fParError[11], fParameters[12], fParError[12] ));
    pt1->AddText(Form("2NDBD: %0.2E#pm%0.2E -- Bi-210: %0.2E#pm%0.2E" , fParameters[8], fParError[8], fParameters[25], fParError[25] ));
    pt1->Draw();

    TLegend *legfit1 = new TLegend(0.8,0.8,0.97,0.97);
    legfit1->AddEntry(fModelTotM1, "Total PDF", "l");
    legfit1->AddEntry(fModelTotThM1, "Total Th-232", "l");
    legfit1->AddEntry(fModelTotRaM1, "Total Ra-226", "l");
    legfit1->AddEntry(fModelTotKM1, "Total K-40", "l");
    legfit1->AddEntry(fModelTotCoM1, "Total Co-60", "l");
    legfit1->AddEntry(fModelTotNDBDM1, "NDBD", "l");
    legfit1->AddEntry(fModelTot2NDBDM1, "2NDBD", "l");
    legfit1->AddEntry(fModelTotBiM1, "Bi-207", "l");
    legfit1->AddEntry(fModelTotMnM1, "Mn-54", "l");
    legfit1->AddEntry(fModelTotPbM1 , "Bi-210", "l");    
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
    fModelTotM2->SetLineWidth(1);
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
    fModelTot2NDBDM2->SetLineColor(46);
    fModelTot2NDBDM2->SetLineStyle(2);
    fModelTotBiM2->SetLineColor(5);
    fModelTotBiM2->SetLineStyle(2);
    fModelTotMnM2->SetLineColor(40);
    fModelTotMnM2->SetLineStyle(2);

    fModelTotPbM2->SetLineStyle(2);
    fModelTotPbM2->SetLineColor(38);
    fTotCorrection->SetLineStyle(2);
    fTotCorrection->SetLineColor(kBlue+2);

    fModelTotThM2->Draw("SAME");
    fModelTotRaM2->Draw("SAME");
    fModelTotKM2->Draw("SAME");
    fModelTotCoM2->Draw("SAME");
    fModelTotNDBDM2->Draw("SAME");
    fModelTot2NDBDM2->Draw("SAME");
    fModelTotBiM2->Draw("SAME");    
    fModelTotMnM2->Draw("SAME");

    fModelTotPbM2->Draw("SAME");
    fTotCorrection->Draw("SAME");    
/*
    // Few Parameters
    TPaveText *pt2 = new TPaveText(0.35,0.78,0.70,0.98,"NB NDC");
    pt2->AddText(Form("Fit Range (M2): %.0f to %.0f keV -- #chi^{2}/NDF: %0.3f", dFitMin, dFitMax, (GetChiSquare()/(2*(dFitMax-dFitMin)/dBinSize - dNumParameters)) ));
    pt2->AddText(Form("Close Th: %0.2E#pm%0.2E --- Far Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
    // pt2->AddText(Form("Surface Th: %0.2E#pm%0.2E --- Far Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
    pt2->AddText(Form("Close Ra: %0.2E#pm%0.2E --- Far Ra: %0.2E#pm%0.2E", fParameters[2], fParError[2], fParameters[3], fParError[3] ));
    pt2->AddText(Form("Close K: %0.2E#pm%0.2E --- Far K: %0.2E#pm%0.2E", fParameters[4], fParError[4], fParameters[5], fParError[5] ));
    pt2->AddText(Form("Close Co: %0.2E#pm%0.2E --- Far Co: %0.2E#pm%0.2E", fParameters[6], fParError[6], fParameters[7], fParError[7] ));
    pt2->AddText(Form("Bi-207: %0.2E#pm%0.2E --- NDBD: %0.2E#pm%0.2E", fParameters[10], fParError[10], fParameters[9], fParError[9] ));
    pt2->AddText(Form("Close Mn-54: %0.2E#pm%0.2E --- Far Mn-54: %0.2E#pm%0.2E", fParameters[11], fParError[11], fParameters[12], fParError[12] ));
    pt2->AddText(Form("2NDBD: %0.2E#pm%0.2E", fParameters[8], fParError[8] ));
*/

    // Many Parameters
    TPaveText *pt2 = new TPaveText(0.35,0.77,0.70,0.99,"NB NDC");
    pt2->AddText(Form("Fit Range (M2): %.0f to %.0f keV -- #chi^{2}/NDF: %0.3f", dFitMin, dFitMax, (GetChiSquare()/(2*(dFitMax-dFitMin)/dBinSize - dNumParameters)) ));
    pt2->AddText(Form("Frame Th: %0.2E#pm%0.2E --- TShield Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
    pt2->AddText(Form("50mK Th: %0.2E#pm%0.2E --- 600mK Th: %0.2E#pm%0.2E", fParameters[13], fParError[13], fParameters[14], fParError[14] ));
    pt2->AddText(Form("IVC Th: %0.2E#pm%0.2E --- OVC Th: %0.2E#pm%0.2E", fParameters[15], fParError[15], fParameters[16], fParError[16] ));
    pt2->AddText(Form("Frame Ra: %0.2E#pm%0.2E --- TShield Ra: %0.2E#pm%0.2E", fParameters[2], fParError[2], fParameters[3], fParError[3] ));
    pt2->AddText(Form("50mK Ra: %0.2E#pm%0.2E --- 600mK Ra: %0.2E#pm%0.2E", fParameters[17], fParError[17], fParameters[18], fParError[18] ));
    pt2->AddText(Form("IVC Ra: %0.2E#pm%0.2E --- OVC Ra: %0.2E#pm%0.2E", fParameters[19], fParError[19], fParameters[20], fParError[20] ));
    pt2->AddText(Form("Close K: %0.2E#pm%0.2E --- Far K: %0.2E#pm%0.2E", fParameters[4], fParError[4], fParameters[5], fParError[5] ));
    pt2->AddText(Form("Close Co: %0.2E#pm%0.2E --- Far Co: %0.2E#pm%0.2E", fParameters[6], fParError[6], fParameters[7], fParError[7] ));
    pt2->AddText(Form("Bi-207: %0.2E#pm%0.2E --- NDBD: %0.2E#pm%0.2E", fParameters[10], fParError[10], fParameters[9], fParError[9] ));
    pt2->AddText(Form("Close Mn-54: %0.2E#pm%0.2E --- Far Mn-54: %0.2E#pm%0.2E", fParameters[11], fParError[11], fParameters[12], fParError[12] ));
    pt2->AddText(Form("2NDBD: %0.2E#pm%0.2E -- Bi-210: %0.2E#pm%0.2E" , fParameters[8], fParError[8], fParameters[25], fParError[25] ));
    pt2->Draw();


    TLegend *legfit2 = new TLegend(0.8,0.8,0.97,0.97);
    legfit2->AddEntry(fModelTotM2, "Total PDF", "l");
    legfit2->AddEntry(fModelTotThM2, "Total Th-232", "l");
    legfit2->AddEntry(fModelTotRaM2, "Total Ra-226", "l");
    legfit2->AddEntry(fModelTotKM2, "Total K-40", "l");
    legfit2->AddEntry(fModelTotCoM2, "Total Co-60", "l");
    legfit2->AddEntry(fModelTotNDBDM2, "NDBD", "l");
    legfit2->AddEntry(fModelTot2NDBDM2, "2NDBD", "l");
    legfit2->AddEntry(fModelTotBiM2, "Bi-207", "l");
    legfit2->AddEntry(fModelTotMnM2, "Mn-54", "l");
    legfit2->AddEntry(fModelTotPbM2 , "Bi-210", "l");
    legfit2->AddEntry(fTotCorrection, "Accidental coincidence correction", "l");
    legfit2->Draw();


  }


	// Residuals
	TCanvas *cResidual1 = new TCanvas("cResidual1", "cResidual1", 1200, 800);
  hResidualGausM1 = new TH1D("hResidualGausM1", "Residual Distribution (M1)", 100, -50, 50);
	hResidualDistM1 = CalculateResiduals(fModelTotM1, fDataHistoM1, hResidualGausM1);
  hResidualDistM1->SetLineColor(kBlack);
	hResidualDistM1->SetName("Residuals");
  hResidualDistM1->SetTitle("Fit Residuals (M1)");
  hResidualDistM1->SetMarkerStyle(25);
	hResidualDistM1->GetXaxis()->SetTitle("Energy (keV)");
	// hResidualDistM1->GetXaxis()->SetTitleSize(0.04);
	// hResidualDistM1->GetXaxis()->SetLabelSize(0.05);
	// hResidualDistM1->GetYaxis()->SetLabelSize(0.05);
	// hResidualDistM1->GetYaxis()->SetTitleSize(0.04);	
	hResidualDistM1->GetYaxis()->SetTitle("Residuals (#sigma)");

	hResidualDistM1->GetXaxis()->SetRange(dFitMin/dBinSize-5, dFitMax/dBinSize+5);
	hResidualDistM1->Draw("E");

/*
  TLine *lineth = new TLine();
  lineth->SetLineStyle(9);
  lineth->SetLineWidth(1);
  lineth->SetLineColor(3);
  lineth->DrawLine(2615, -15, 2615, 15);
  lineth->DrawLine(2104, -15, 2104, 15);
  lineth->DrawLine(1593, -15, 1593, 15);
  lineth->DrawLine(968, -15, 968, 15);
  lineth->DrawLine(911, -15, 911, 15);
  lineth->DrawLine(583, -15, 583, 15);


  TLine *linera = new TLine();
  linera->SetLineStyle(9);
  linera->SetLineWidth(1);
  linera->SetLineColor(4);
  linera->DrawLine(609, -13, 609, 15);
  linera->DrawLine(665, -13, 665, 15);
  linera->DrawLine(768, -13, 768, 15);
  linera->DrawLine(806, -13, 806, 15);
  linera->DrawLine(934, -13, 934, 15);
  linera->DrawLine(1120, -13, 1120, 15);
  linera->DrawLine(1155, -13, 1155, 15);
  linera->DrawLine(1238, -13, 1238, 15);
  linera->DrawLine(1377, -13, 1377, 15);
  linera->DrawLine(1408, -13, 1408, 15);
  linera->DrawLine(1729, -13, 1729, 15);
  linera->DrawLine(1764, -13, 1764, 15);
  linera->DrawLine(1847, -13, 1847, 15);
  linera->DrawLine(2204, -13, 2204, 15);
  linera->DrawLine(2447, -13, 2447, 15);

  TLine *linek = new TLine();
  linek->SetLineStyle(9);
  linek->SetLineWidth(1);
  linek->SetLineColor(6);
  linek->DrawLine(1461, -13, 1461, 15);

  TLine *lineco = new TLine();
  lineco->SetLineStyle(9);
  lineco->SetLineWidth(1);
  lineco->SetLineColor(7);
  lineco->DrawLine(1173, -13, 1173, 15);
  lineco->DrawLine(1332, -13, 1332, 15);
*/



  TCanvas *cres1 = new TCanvas();
  hResidualGausM1->Draw();

  TCanvas *cResidual2 = new TCanvas("cResidual2", "cResidual2", 1200, 800);
  hResidualGausM2 = new TH1D("hResidualGausM2", "Residual Distribution (M2)", 100, -50, 50);  
  hResidualDistM2 = CalculateResiduals(fModelTotM2, fDataHistoM2, hResidualGausM2);
  hResidualDistM2->SetLineColor(kBlack);
  hResidualDistM2->SetMarkerStyle(25);
  hResidualDistM2->SetName("Residuals");
  hResidualDistM2->SetTitle("Fit Residuals (M2)");
  hResidualDistM2->GetXaxis()->SetTitle("Energy (keV)");
  // hResidualDistM2->GetXaxis()->SetTitleSize(0.04);
  // hResidualDistM2->GetXaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetTitleSize(0.04); 
  hResidualDistM2->GetYaxis()->SetTitle("Residuals (#sigma)");

  hResidualDistM2->GetXaxis()->SetRange(dFitMin/dBinSize-5, dFitMax/dBinSize+5);
  hResidualDistM2->Draw("E");

  TCanvas *cres2 = new TCanvas();
  hResidualGausM2->Draw();

  // Output integrals of stuff for limits
  cout << "Integral Data in ROI: " << fDataHistoM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fDataHistoM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;
  cout << "Integral Total PDF in ROI: " << fModelTotM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTotM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;
  cout << "Integral Total Th PDF in ROI: " << fModelTotThM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTotThM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;
  cout << "Integral Total Ra PDF in ROI: " << fModelTotRaM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTotRaM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;
  cout << "Integral Total Co PDF in ROI: " << fModelTotCoM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTotCoM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;
  cout << "Integral Total K PDF in ROI: " << fModelTotKM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTotKM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;
  cout << "Integral Total Bi-207 PDF in ROI: " << fModelTotBiM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTotBiM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;  
  cout << "Integral Total 2NDBD PDF in ROI: " << fModelTot2NDBDM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTot2NDBDM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;
  cout << "Integral Total 0NDBD PDF in ROI: " << fModelTotNDBDM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTotNDBDM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;
  cout << "Integral Frame Th PDF in ROI: " << fParameters[0]*fSmearFrameThM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << endl;
  cout << "Integral TShield Th PDF in ROI: " << fParameters[1]*fSmearTShieldThM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << endl;
  cout << "Integral 50mK Th PDF in ROI: " << fParameters[13]*fSmear50mKThM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << endl;
  cout << "Integral 600mK Th PDF in ROI: " << fParameters[14]*fSmear600mKThM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << endl;
  cout << "Integral IVC Th PDF in ROI: " << fParameters[15]*fSmearIVCThM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << endl;
  cout << "Integral OVC Th PDF in ROI: " << fParameters[16]*fSmearOVCThM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << endl;
  cout << "Integral Bi-210 (300 keV to 1000 keV): " << fParameters[25]*fSmearCrystalBi2M1->Integral(300/dBinSize, 1000/dBinSize) << " +/- " << endl;

  cout << "M2/(M1+M2) = " << (double)fModelTotM2->Integral(300/dBinSize, 3000/dBinSize)/(fModelTotM1->Integral(300/dBinSize, 3000/dBinSize)+fModelTotM2->Integral(300/dBinSize, 3000/dBinSize)) << endl;


  // cout << fTotCorrection->Integral(1300/dBinSize, 1400/dBinSize) << endl;

  // Write
/*
  TH1D  *hCloneResultM1    = (TH1D*)fModelTotM1->Clone("fModelTotM1");
  // TH1D  *hCloneResultM2    = (TH1D*)fModelTotM2->Clone("fModelTotM2");
  NormalizePDF(hCloneResultM1, 0, 2700);
  TFile *fFileResult = new TFile("Result-2keV.root", "RECREATE");
  hCloneResultM1->Write();
  // hCloneResultM2->Write();
  fFileResult->Write();
*/
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

  TLegend *legth2 = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legra2 = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legk2 = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legco2 = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legndbd2 = new TLegend(0.8,0.890,0.925,0.925);
  TLegend *legbi2 = new TLegend(0.8,0.890,0.925,0.925);
  TLegend *legths2 = new TLegend(0.75,0.80,0.925,0.925);
  TLegend *legras2 = new TLegend(0.75,0.80,0.925,0.925);


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

  fModelFrameKM1->GetXaxis()->SetRange(0/dBinSize, 1600/dBinSize);
  fModelFrameKM1->DrawNormalized();
  fModelTShieldKM1->DrawNormalized("SAME");
  fModel50mKKM1->DrawNormalized("SAME");
  fModel600mKKM1->DrawNormalized("SAME");
  fModelIVCKM1->DrawNormalized("SAME");
  fModelOVCKM1->DrawNormalized("SAME");

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

  fModelFrameKM2->GetXaxis()->SetRange(0/dBinSize, 1600/dBinSize);
  fModelFrameKM2->DrawNormalized();
  fModelTShieldKM2->DrawNormalized("SAME");
  fModel50mKKM2->DrawNormalized("SAME");
  fModel600mKKM2->DrawNormalized("SAME");
  fModelIVCKM2->DrawNormalized("SAME");
  fModelOVCKM2->DrawNormalized("SAME");

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
  fModel2NDBDM1->SetLineColor(2);  
  fModelNDBDM1->GetXaxis()->SetTitle("Energy (keV)");
  fModelNDBDM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));  
  fModelNDBDM1->Draw();
  fModelNDBDM1->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);
  fModel2NDBDM1->Draw("SAME");  
  legndbd1->AddEntry(fModelNDBDM1, "0#nu#beta#beta (M1)" ,"l");
  legndbd1->AddEntry(fModel2NDBDM1, "2#nu#beta#beta (M1)" ,"l");
  legndbd1->Draw();

  TCanvas *cNDBD2 = new TCanvas("cNDBD2", "cNDBD2", 1200, 800);
  cNDBD2->SetLogy();
  fModelNDBDM2->SetLineColor(1);
  fModel2NDBDM2->SetLineColor(2);
  fModelNDBDM2->GetXaxis()->SetTitle("Energy (keV)");
  fModelNDBDM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));  
  fModelNDBDM2->Draw();
  fModel2NDBDM2->Draw("SAME");
  fModelNDBDM2->GetXaxis()->SetRange(0/dBinSize, 2700/dBinSize);
  legndbd2->AddEntry(fModelNDBDM2, "0#nu#beta#beta (M2)" ,"l");
  legndbd2->AddEntry(fModel2NDBDM2, "2#nu#beta#beta (M2)" ,"l");
  legndbd2->Draw();


  TCanvas *cBi1 = new TCanvas("cBi1", "cBi1", 1200, 800);
  cBi1->SetLogy();
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


// Generates Toy Data using MC histograms // Needs to be updated
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
    // modelm1_i = fModelTest2->GetBinContent(i); // For testing

		// Log-likelihood Chi-Squared
    // Avoiding 0's... correct or no?
    // This only doens't work if model_i = 0, need to make sure statistics are high enough in every bin for model
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

double TBackgroundModel::GetChiSquareAdaptive()
{
  double chiSquare = 0.;
  double datam1_i, errm1_i;
  double datam2_i, errm2_i;
  double modelm1_i, modelm2_i;  

  for(int i = dFitMinBinM1 ; i <= dFitMaxBinM1; i++)
  {

    datam1_i = fAdapDataHistoM1->GetBinContent(i); // For real data

    modelm1_i = fModelTotAdapM1->GetBinContent(i);

    if(modelm1_i != 0 && datam1_i != 0)
    {
      chiSquare += 2 * (modelm1_i - datam1_i + datam1_i * TMath::Log(datam1_i/modelm1_i));
    }
  }

  for(int i = dFitMinBinM2; i <= dFitMaxBinM2; i++)
  {
    datam2_i = fAdapDataHistoM2->GetBinContent(i); // For real data

    modelm2_i = fModelTotAdapM2->GetBinContent(i);

    if(modelm2_i != 0 && datam2_i != 0)
    {
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
	// LoadData();	

  // If using unsmeared
  if(bUnSmeared)
  {
    LoadPDFs();
  }
  // If using smeared custom
  else
  {
    cout << "Loading Smeared Histograms from file" << endl;
    // fFile = new TFile(Form("MCHist-%dkeV.root", dBinSize));
    fFile = new TFile("MCAdaptiveGamma.root");    
    // fFile = new TFile(Form("MCData-%dkeV.root", dBinSize));

    fSmearFrameThM1   = (TH1D*)fFile->Get("fSmearFrameThM1");
    fSmearTShieldThM1 = (TH1D*)fFile->Get("fSmearTShieldThM1");
    fSmear50mKThM1    = (TH1D*)fFile->Get("fSmear50mKThM1");
    fSmear600mKThM1   = (TH1D*)fFile->Get("fSmear600mKThM1");
    fSmearIVCThM1     = (TH1D*)fFile->Get("fSmearIVCThM1");
    fSmearOVCThM1     = (TH1D*)fFile->Get("fSmearOVCThM1");

    fSmearFrameRaM1   = (TH1D*)fFile->Get("fSmearFrameRaM1");
    fSmearTShieldRaM1 = (TH1D*)fFile->Get("fSmearTShieldRaM1");
    fSmear50mKRaM1    = (TH1D*)fFile->Get("fSmear50mKRaM1");
    fSmear600mKRaM1   = (TH1D*)fFile->Get("fSmear600mKRaM1");
    fSmearIVCRaM1     = (TH1D*)fFile->Get("fSmearIVCRaM1");
    fSmearOVCRaM1     = (TH1D*)fFile->Get("fSmearOVCRaM1");

    fSmearFrameCoM1   = (TH1D*)fFile->Get("fSmearFrameCoM1");
    fSmearTShieldCoM1 = (TH1D*)fFile->Get("fSmearTShieldCoM1");
    fSmear50mKCoM1    = (TH1D*)fFile->Get("fSmear50mKCoM1");
    fSmear600mKCoM1   = (TH1D*)fFile->Get("fSmear600mKCoM1");
    fSmearIVCCoM1     = (TH1D*)fFile->Get("fSmearIVCCoM1");
    fSmearOVCCoM1     = (TH1D*)fFile->Get("fSmearOVCCoM1");    

    fSmearFrameKM1    = (TH1D*)fFile->Get("fSmearFrameKM1");
    fSmearTShieldKM1  = (TH1D*)fFile->Get("fSmearTShieldKM1");
    fSmear50mKKM1     = (TH1D*)fFile->Get("fSmear50mKKM1");
    fSmear600mKKM1    = (TH1D*)fFile->Get("fSmear600mKKM1");
    fSmearIVCKM1      = (TH1D*)fFile->Get("fSmearIVCKM1");
    fSmearOVCKM1      = (TH1D*)fFile->Get("fSmearOVCKM1");    

    fSmearTShieldMnM1 = (TH1D*)fFile->Get("fSmearTShieldMnM1");
    fSmearIVCMnM1     = (TH1D*)fFile->Get("fSmearIVCMnM1");

    fSmearCrystalBi2M1 = (TH1D*)fFile->Get("fSmearCrystalBi2M1");
    fSmearFrameBi2M1 = (TH1D*)fFile->Get("fSmearFrameBi2M1");


    fSmear2NDBDM1     = (TH1D*)fFile->Get("fSmear2NDBDM1");
    fSmearNDBDM1      = (TH1D*)fFile->Get("fSmearNDBDM1");
    fSmearBiM1        = (TH1D*)fFile->Get("fSmearBiM1");

    // fSmearCrystalPtM1 = (TH1D*)fFile->Get("fSmearCrystalPtM1");
    // fSmearCrystalPbBM1 = (TH1D*)fFile->Get("fSmearCrystalPbBM1");
    // fSmearCrystalPbS01M1 = (TH1D*)fFile->Get("fSmearCrystalPbS01M1");
    // fSmearCrystalPbS1M1 = (TH1D*)fFile->Get("fSmearCrystalPbS1M1");
    // fSmearCrystalPbS10M1 = (TH1D*)fFile->Get("fSmearCrystalPbS10M1");
    // fSmearCrystalPbS100M1 = (TH1D*)fFile->Get("fSmearCrystalPbS100M1");
    // fSmearFramePbBM1 = (TH1D*)fFile->Get("fSmearFramePbBM1");
    // fSmearFramePbS01M1 = (TH1D*)fFile->Get("fSmearFramePbS01M1");
    // fSmearFramePbS1M1 = (TH1D*)fFile->Get("fSmearFramePbS1M1");
    // fSmearFramePbS10M1 = (TH1D*)fFile->Get("fSmearFramePbS10M1");
    // fSmearFramePbS100M1 = (TH1D*)fFile->Get("fSmearFramePbS100M1");

    fSmearFrameThM2   = (TH1D*)fFile->Get("fSmearFrameThM2");
    fSmearTShieldThM2 = (TH1D*)fFile->Get("fSmearTShieldThM2");
    fSmear50mKThM2    = (TH1D*)fFile->Get("fSmear50mKThM2");
    fSmear600mKThM2   = (TH1D*)fFile->Get("fSmear600mKThM2");
    fSmearIVCThM2     = (TH1D*)fFile->Get("fSmearIVCThM2");
    fSmearOVCThM2     = (TH1D*)fFile->Get("fSmearOVCThM2");

    fSmearFrameRaM2   = (TH1D*)fFile->Get("fSmearFrameRaM2");
    fSmearTShieldRaM2 = (TH1D*)fFile->Get("fSmearTShieldRaM2");
    fSmear50mKRaM2    = (TH1D*)fFile->Get("fSmear50mKRaM2");
    fSmear600mKRaM2   = (TH1D*)fFile->Get("fSmear600mKRaM2");
    fSmearIVCRaM2     = (TH1D*)fFile->Get("fSmearIVCRaM2");
    fSmearOVCRaM2     = (TH1D*)fFile->Get("fSmearOVCRaM2");

    fSmearFrameCoM2   = (TH1D*)fFile->Get("fSmearFrameCoM2");
    fSmearTShieldCoM2 = (TH1D*)fFile->Get("fSmearTShieldCoM2");
    fSmear50mKCoM2    = (TH1D*)fFile->Get("fSmear50mKCoM2");
    fSmear600mKCoM2   = (TH1D*)fFile->Get("fSmear600mKCoM2");
    fSmearIVCCoM2     = (TH1D*)fFile->Get("fSmearIVCCoM2");
    fSmearOVCCoM2     = (TH1D*)fFile->Get("fSmearOVCCoM2");    

    fSmearFrameKM2    = (TH1D*)fFile->Get("fSmearFrameKM2");
    fSmearTShieldKM2  = (TH1D*)fFile->Get("fSmearTShieldKM2");
    fSmear50mKKM2     = (TH1D*)fFile->Get("fSmear50mKKM2");
    fSmear600mKKM2    = (TH1D*)fFile->Get("fSmear600mKKM2");
    fSmearIVCKM2      = (TH1D*)fFile->Get("fSmearIVCKM2");
    fSmearOVCKM2      = (TH1D*)fFile->Get("fSmearOVCKM2");    

    fSmearTShieldMnM2 = (TH1D*)fFile->Get("fSmearTShieldMnM2");
    fSmearIVCMnM2     = (TH1D*)fFile->Get("fSmearIVCMnM2");

    fSmearCrystalBi2M2 = (TH1D*)fFile->Get("fSmearCrystalBi2M2");
    fSmearFrameBi2M2 = (TH1D*)fFile->Get("fSmearFrameBi2M2");

    fSmear2NDBDM2     = (TH1D*)fFile->Get("fSmear2NDBDM2");
    fSmearNDBDM2      = (TH1D*)fFile->Get("fSmearNDBDM2");
    fSmearBiM2        = (TH1D*)fFile->Get("fSmearBiM2");    
 
    // fSmearCrystalPtM2 = (TH1D*)fFile->Get("fSmearCrystalPtM2");
    // fSmearCrystalPbBM2 = (TH1D*)fFile->Get("fSmearCrystalPbBM2");
    // fSmearCrystalPbS01M2 = (TH1D*)fFile->Get("fSmearCrystalPbS01M2");
    // fSmearCrystalPbS1M2 = (TH1D*)fFile->Get("fSmearCrystalPbS1M2");
    // fSmearCrystalPbS10M2 = (TH1D*)fFile->Get("fSmearCrystalPbS10M2");
    // fSmearCrystalPbS100M2 = (TH1D*)fFile->Get("fSmearCrystalPbS100M2");
    // fSmearFramePbBM2 = (TH1D*)fFile->Get("fSmearFramePbBM2");
    // fSmearFramePbS01M2 = (TH1D*)fFile->Get("fSmearFramePbS01M2"); 
    // fSmearFramePbS1M2 = (TH1D*)fFile->Get("fSmearFramePbS1M2"); 
    // fSmearFramePbS10M2 = (TH1D*)fFile->Get("fSmearFramePbS10M2"); 
    // fSmearFramePbS100M2 = (TH1D*)fFile->Get("fSmearFramePbS100M2"); 

    fAdap600mKThM1   = (TH1D*)fFile->Get("fAdap600mKThM1");
    fAdapIVCThM1     = (TH1D*)fFile->Get("fAdapIVCThM1");
    fAdapOVCThM1     = (TH1D*)fFile->Get("fAdapOVCThM1");

    fAdap600mKRaM1   = (TH1D*)fFile->Get("fAdap600mKRaM1");
    fAdapOVCRaM1     = (TH1D*)fFile->Get("fAdapOVCRaM1");

    fAdapFrameCoM1   = (TH1D*)fFile->Get("fAdapFrameCoM1");
    fAdapOVCCoM1     = (TH1D*)fFile->Get("fAdapOVCCoM1");    

    fAdapFrameKM1    = (TH1D*)fFile->Get("fAdapFrameKM1");
    fAdapTShieldKM1  = (TH1D*)fFile->Get("fAdapTShieldKM1");
    fAdap50mKKM1     = (TH1D*)fFile->Get("fAdap50mKKM1");
    fAdap600mKKM1    = (TH1D*)fFile->Get("fAdap600mKKM1");
    fAdapIVCKM1      = (TH1D*)fFile->Get("fAdapIVCKM1");
    fAdapOVCKM1      = (TH1D*)fFile->Get("fAdapOVCKM1");    

    fAdapIVCMnM1     = (TH1D*)fFile->Get("fAdapIVCMnM1");

    fAdap2NDBDM1     = (TH1D*)fFile->Get("fAdap2NDBDM1");
    fAdapNDBDM1      = (TH1D*)fFile->Get("fAdapNDBDM1");
    fAdapBiM1        = (TH1D*)fFile->Get("fAdapBiM1");  


    fAdap600mKThM2   = (TH1D*)fFile->Get("fAdap600mKThM2");
    fAdapIVCThM2     = (TH1D*)fFile->Get("fAdapIVCThM2");
    fAdapOVCThM2     = (TH1D*)fFile->Get("fAdapOVCThM2");

    fAdap600mKRaM2   = (TH1D*)fFile->Get("fAdap600mKRaM2");
    fAdapOVCRaM2     = (TH1D*)fFile->Get("fAdapOVCRaM2");

    fAdapFrameCoM2   = (TH1D*)fFile->Get("fAdapFrameCoM2");
    fAdapOVCCoM2     = (TH1D*)fFile->Get("fAdapOVCCoM2");    

    fAdapFrameKM2    = (TH1D*)fFile->Get("fAdapFrameKM2");
    fAdapTShieldKM2  = (TH1D*)fFile->Get("fAdapTShieldKM2");
    fAdap50mKKM2     = (TH1D*)fFile->Get("fAdap50mKKM2");
    fAdap600mKKM2    = (TH1D*)fFile->Get("fAdap600mKKM2");
    fAdapIVCKM2      = (TH1D*)fFile->Get("fAdapIVCKM2");
    fAdapOVCKM2      = (TH1D*)fFile->Get("fAdapOVCKM2");    

    fAdapIVCMnM2     = (TH1D*)fFile->Get("fAdapIVCMnM2");

    fAdap2NDBDM2     = (TH1D*)fFile->Get("fAdap2NDBDM2");
    fAdapNDBDM2      = (TH1D*)fFile->Get("fAdapNDBDM2");
    fAdapBiM2        = (TH1D*)fFile->Get("fAdapBiM2");   

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

  // qtree->Add("/Users/brian/macros/CUOREZ/Bkg/Q0_DR2_BackgroundSignalData.root");	
  qtree->Add("/Users/brian/macros/CUOREZ/Bkg/ReducedBkg-ds*.root");
  qtree->Project("fDataHistoTot", "Energy", base_cut);
  qtree->Project("fDataHistoM1", 	"Energy", base_cut);
  qtree->Project("fDataHistoM2", 	"Energy");

	cout << "Loaded Data" << endl;

	// Normalizing data (don't!)
	// bin 0 = underflow, bin dNBins = last bin with upper-edge xup Excluded

	dDataIntegral = fDataHistoTot->Integral(1, dNBins);
  int dDataIntegralTot = qtree->GetEntries();

  cout << "Total Events in background spectrum: " << dDataIntegralTot << endl; 
	cout << "Events in background spectrum (M1): " << dDataIntegral << endl;
  cout << "Events in background spectrum (M2): " << fDataHistoM2->Integral(1, dNBins) << endl;

  // Scale by Live-time (ds 2061 - 2100) 14647393.0 seconds
  fDataHistoM1->Scale(1/((936398+14647393.0) * dSecToYears));
  fDataHistoM2->Scale(1/((936398+14647393.0) * dSecToYears));  

  cout << "Normalized Data using Livetime of: " << (936398+14647393.0) * dSecToYears << " years" <<endl;

}

void TBackgroundModel::LoadPDFs()
{
  // Fills and Loads MC data
  // Load M1
  outTreeCrystalPt190BM1    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pt190", "B", 1);
  // outTreeCrystalPt190S01M1    = LoadMC(dDataDir.c_str(),  "Crystal", "Pt190", "S01", 1);
  outTreeCrystalPt190S1M1    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pt190", "S1", 1);
  outTreeCrystalPt190S10M1    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pt190", "S10", 1);
  outTreeCrystalPt190S100M1    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pt190", "S100", 1);

  outTreeCrystalPb210BM1    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pb210", "B", 1);
  outTreeCrystalPb210S01M1    = LoadMC(dDataDir.c_str(),  "Crystal", "Pb210", "S01", 1);
  outTreeCrystalPb210S1M1    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pb210", "S1", 1);
  outTreeCrystalPb210S10M1    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pb210", "S10", 1);
  outTreeCrystalPb210S100M1    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pb210", "S100", 1);

  outTreeCrystalRa226BM1    = LoadMC(dDataDir.c_str(),  "Crystal",  "Ra226", "B", 1);
  outTreeCrystalRa226S01M1    = LoadMC(dDataDir.c_str(),  "Crystal", "Ra226", "S01", 1);
  outTreeCrystalRa226S1M1    = LoadMC(dDataDir.c_str(),  "Crystal",  "Ra226", "S1", 1);
  outTreeCrystalRa226S10M1    = LoadMC(dDataDir.c_str(),  "Crystal",  "Ra226", "S10", 1);
  outTreeCrystalRa226S100M1    = LoadMC(dDataDir.c_str(),  "Crystal",  "Ra226", "S100", 1);

  outTreeFramePb210BM1    = LoadMC(dDataDir.c_str(),  "Frame",  "Pb210", "B", 1);
  outTreeFramePb210S01M1    = LoadMC(dDataDir.c_str(),  "Frame", "Pb210", "S01", 1);
  outTreeFramePb210S1M1    = LoadMC(dDataDir.c_str(),  "Frame",  "Pb210", "S1", 1);
  outTreeFramePb210S10M1    = LoadMC(dDataDir.c_str(),  "Frame",  "Pb210", "S10", 1);
  outTreeFramePb210S100M1    = LoadMC(dDataDir.c_str(),  "Frame",  "Pb210", "S100", 1);

  // outTreeFrameRa226BM1    = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "B", 1);
  outTreeFrameRa226S01M1    = LoadMC(dDataDir.c_str(),  "Frame", "Ra226", "S01", 1);
  outTreeFrameRa226S1M1    = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S1", 1);
  outTreeFrameRa226S10M1    = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S10", 1);
  outTreeFrameRa226S100M1    = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S100", 1);


  outTreeCrystalPt190BM2    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pt190", "B", 2);
  // outTreeCrystalPt190S01M2    = LoadMC(dDataDir.c_str(),  "Crystal", "Pt190", "S01", 2);
  outTreeCrystalPt190S1M2    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pt190", "S1", 2);
  outTreeCrystalPt190S10M2    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pt190", "S10", 2);
  outTreeCrystalPt190S100M2    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pt190", "S100", 2);

  outTreeCrystalPb210BM2    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pb210", "B", 2);
  outTreeCrystalPb210S01M2    = LoadMC(dDataDir.c_str(),  "Crystal", "Pb210", "S01", 2);
  outTreeCrystalPb210S1M2    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pb210", "S1", 2);
  outTreeCrystalPb210S10M2    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pb210", "S10", 2);
  outTreeCrystalPb210S100M2    = LoadMC(dDataDir.c_str(),  "Crystal",  "Pb210", "S100", 2);

  outTreeCrystalRa226BM2    = LoadMC(dDataDir.c_str(),  "Crystal",  "Ra226", "B", 2);
  outTreeCrystalRa226S01M2    = LoadMC(dDataDir.c_str(),  "Crystal", "Ra226", "S01", 2);
  outTreeCrystalRa226S1M2    = LoadMC(dDataDir.c_str(),  "Crystal",  "Ra226", "S1", 2);
  outTreeCrystalRa226S10M2    = LoadMC(dDataDir.c_str(),  "Crystal",  "Ra226", "S10", 2);
  outTreeCrystalRa226S100M2    = LoadMC(dDataDir.c_str(),  "Crystal",  "Ra226", "S100", 2);

  outTreeFramePb210BM2    = LoadMC(dDataDir.c_str(),  "Frame",  "Pb210", "B", 2);
  outTreeFramePb210S01M2    = LoadMC(dDataDir.c_str(),  "Frame", "Pb210", "S01", 2);
  outTreeFramePb210S1M2    = LoadMC(dDataDir.c_str(),  "Frame",  "Pb210", "S1", 2);
  outTreeFramePb210S10M2    = LoadMC(dDataDir.c_str(),  "Frame",  "Pb210", "S10", 2);
  outTreeFramePb210S100M2    = LoadMC(dDataDir.c_str(),  "Frame",  "Pb210", "S100", 2);

  // outTreeFrameRa226BM2    = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "B", 2);
  outTreeFrameRa226S01M2    = LoadMC(dDataDir.c_str(),  "Frame", "Ra226", "S01", 2);
  outTreeFrameRa226S1M2    = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S1", 2);
  outTreeFrameRa226S10M2    = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S10", 2);
  outTreeFrameRa226S100M2    = LoadMC(dDataDir.c_str(),  "Frame",  "Ra226", "S100", 2);


  outTreeCrystalPt190BM1->Project("fModelCrystalPt190BM1",        "Ener2");
  // outTreeCrystalPt190S01M1->Project("fModelCrystalPt190S01M1",     "Ener2");
  outTreeCrystalPt190S1M1->Project("fModelCrystalPt190S1M1",        "Ener2");
  outTreeCrystalPt190S10M1->Project("fModelCrystalPt190S10M1",        "Ener2");
  outTreeCrystalPt190S100M1->Project("fModelCrystalPt190S100M1",        "Ener2");

  outTreeCrystalPb210BM1->Project("fModelCrystalPb210BM1",        "Ener2");
  outTreeCrystalPb210S01M1->Project("fModelCrystalPb210S01M1",     "Ener2");
  outTreeCrystalPb210S1M1->Project("fModelCrystalPb210S1M1",        "Ener2");
  outTreeCrystalPb210S10M1->Project("fModelCrystalPb210S10M1",        "Ener2");
  outTreeCrystalPb210S100M1->Project("fModelCrystalPb210S100M1",        "Ener2");

  outTreeCrystalRa226BM1->Project("fModelCrystalRa226BM1",        "Ener2");
  outTreeCrystalRa226S01M1->Project("fModelCrystalRa226S01M1",     "Ener2");
  outTreeCrystalRa226S1M1->Project("fModelCrystalRa226S1M1",        "Ener2");
  outTreeCrystalRa226S10M1->Project("fModelCrystalRa226S10M1",        "Ener2");
  outTreeCrystalRa226S100M1->Project("fModelCrystalRa226S100M1",        "Ener2");


  outTreeFramePb210BM1->Project("fModelFramePb210BM1",        "Ener2");
  outTreeFramePb210S01M1->Project("fModelFramePb210S01M1",     "Ener2");
  outTreeFramePb210S1M1->Project("fModelFramePb210S1M1",        "Ener2");
  outTreeFramePb210S10M1->Project("fModelFramePb210S10M1",        "Ener2");
  outTreeFramePb210S100M1->Project("fModelFramePb210S100M1",        "Ener2");

  // outTreeFrameRa226BM1->Project("fModelFrameRa226BM1",        "Ener2");
  outTreeFrameRa226S01M1->Project("fModelFrameRa226S01M1",     "Ener2");
  outTreeFrameRa226S1M1->Project("fModelFrameRa226S1M1",        "Ener2");
  outTreeFrameRa226S10M1->Project("fModelFrameRa226S10M1",        "Ener2");
  outTreeFrameRa226S100M1->Project("fModelFrameRa226S100M1",        "Ener2");



  outTreeCrystalPt190BM2->Project("fModelCrystalPt190BM2",        "Ener2");
  // outTreeCrystalPt190S01M2->Project("fModelCrystalPt190S01M2",     "Ener2");
  outTreeCrystalPt190S1M2->Project("fModelCrystalPt190S1M2",        "Ener2");
  outTreeCrystalPt190S10M2->Project("fModelCrystalPt190S10M2",        "Ener2");
  outTreeCrystalPt190S100M2->Project("fModelCrystalPt190S100M2",        "Ener2");

  outTreeCrystalPb210BM2->Project("fModelCrystalPb210BM2",        "Ener2");
  outTreeCrystalPb210S01M2->Project("fModelCrystalPb210S01M2",     "Ener2");
  outTreeCrystalPb210S1M2->Project("fModelCrystalPb210S1M2",        "Ener2");
  outTreeCrystalPb210S10M2->Project("fModelCrystalPb210S10M2",        "Ener2");
  outTreeCrystalPb210S100M2->Project("fModelCrystalPb210S100M2",        "Ener2");

  outTreeCrystalRa226BM2->Project("fModelCrystalRa226BM2",        "Ener2");
  outTreeCrystalRa226S01M2->Project("fModelCrystalRa226S01M2",     "Ener2");
  outTreeCrystalRa226S1M2->Project("fModelCrystalRa226S1M2",        "Ener2");
  outTreeCrystalRa226S10M2->Project("fModelCrystalRa226S10M2",        "Ener2");
  outTreeCrystalRa226S100M2->Project("fModelCrystalRa226S100M2",        "Ener2");


  outTreeFramePb210BM2->Project("fModelFramePb210BM2",        "Ener2");
  outTreeFramePb210S01M2->Project("fModelFramePb210S01M2",     "Ener2");
  outTreeFramePb210S1M2->Project("fModelFramePb210S1M2",        "Ener2");
  outTreeFramePb210S10M2->Project("fModelFramePb210S10M2",        "Ener2");
  outTreeFramePb210S100M2->Project("fModelFramePb210S100M2",        "Ener2");

  // outTreeFrameRa226BM2->Project("fModelFrameRa226BM2",        "Ener2");
  outTreeFrameRa226S01M2->Project("fModelFrameRa226S01M2",     "Ener2");
  outTreeFrameRa226S1M2->Project("fModelFrameRa226S1M2",        "Ener2");
  outTreeFrameRa226S10M2->Project("fModelFrameRa226S10M2",        "Ener2");
  outTreeFrameRa226S100M2->Project("fModelFrameRa226S100M2",        "Ener2");

  NormalizePDFPair(fModelCrystalPt190BM1, fModelCrystalPt190BM2,      50, 10000);
  // NormalizePDFPair(fModelCrystalPt190S01M1, fModelCrystalPt190S01M2,      50, 10000);
  NormalizePDFPair(fModelCrystalPt190S1M1, fModelCrystalPt190S1M2,      50, 10000);
  NormalizePDFPair(fModelCrystalPt190S10M1, fModelCrystalPt190S10M2,      50, 10000);
  NormalizePDFPair(fModelCrystalPt190S100M1, fModelCrystalPt190S100M2,      50, 10000);

  NormalizePDFPair(fModelCrystalPb210BM1, fModelCrystalPb210BM2,      50, 10000);
  NormalizePDFPair(fModelCrystalPb210S01M1, fModelCrystalPb210S01M2,      50, 10000);
  NormalizePDFPair(fModelCrystalPb210S1M1, fModelCrystalPb210S1M2,      50, 10000);
  NormalizePDFPair(fModelCrystalPb210S10M1, fModelCrystalPb210S10M2,      50, 10000);
  NormalizePDFPair(fModelCrystalPb210S100M1, fModelCrystalPb210S100M2,      50, 10000);

  NormalizePDFPair(fModelCrystalRa226BM1, fModelCrystalRa226BM2,      50, 10000);
  NormalizePDFPair(fModelCrystalRa226S01M1, fModelCrystalRa226S01M2,      50, 10000);
  NormalizePDFPair(fModelCrystalRa226S1M1, fModelCrystalRa226S1M2,      50, 10000);
  NormalizePDFPair(fModelCrystalRa226S10M1, fModelCrystalRa226S10M2,      50, 10000);
  NormalizePDFPair(fModelCrystalRa226S100M1, fModelCrystalRa226S100M2,      50, 10000);


  NormalizePDFPair(fModelFramePb210BM1, fModelFramePb210BM2,      50, 10000);
  NormalizePDFPair(fModelFramePb210S01M1, fModelFramePb210S01M2,      50, 10000);
  NormalizePDFPair(fModelFramePb210S1M1, fModelFramePb210S1M2,      50, 10000);
  NormalizePDFPair(fModelFramePb210S10M1, fModelFramePb210S10M2,      50, 10000);
  NormalizePDFPair(fModelFramePb210S100M1, fModelFramePb210S100M2,      50, 10000);

  // NormalizePDFPair(fModelFrameRa226BM1, fModelFrameRa226BM2,      50, 10000);
  NormalizePDFPair(fModelFrameRa226S01M1, fModelFrameRa226S01M2,      50, 10000);
  NormalizePDFPair(fModelFrameRa226S1M1, fModelFrameRa226S1M2,      50, 10000);
  NormalizePDFPair(fModelFrameRa226S10M1, fModelFrameRa226S10M2,      50, 10000);
  NormalizePDFPair(fModelFrameRa226S100M1, fModelFrameRa226S100M2,      50, 10000);

  cout << "Normalized Histograms, beginning smearing" << endl;

  // Smearing, needs to be updated
  double dRes1 = 1.986;
  double dRes2 = 5.332;

  // M1
  SmearMC(fModelCrystalPt190BM1, fSmearCrystalPt190BM1, dRes1 , dRes2);
  // SmearMC(fModelCrystalPt190S01M1, fSmearCrystalPt190S01M1, dRes1 , dRes2);
  SmearMC(fModelCrystalPt190S1M1, fSmearCrystalPt190S1M1, dRes1 , dRes2);
  SmearMC(fModelCrystalPt190S10M1, fSmearCrystalPt190S10M1, dRes1 , dRes2);
  SmearMC(fModelCrystalPt190S100M1, fSmearCrystalPt190S100M1, dRes1 , dRes2);

  SmearMC(fModelCrystalPb210BM1, fSmearCrystalPb210BM1, dRes1 , dRes2);
  SmearMC(fModelCrystalPb210S01M1, fSmearCrystalPb210S01M1, dRes1 , dRes2);
  SmearMC(fModelCrystalPb210S1M1, fSmearCrystalPb210S1M1, dRes1 , dRes2);
  SmearMC(fModelCrystalPb210S10M1, fSmearCrystalPb210S10M1, dRes1 , dRes2);
  SmearMC(fModelCrystalPb210S100M1, fSmearCrystalPb210S100M1, dRes1 , dRes2);

  SmearMC(fModelCrystalRa226BM1, fSmearCrystalRa226BM1, dRes1 , dRes2);
  SmearMC(fModelCrystalRa226S01M1, fSmearCrystalRa226S01M1, dRes1 , dRes2);
  SmearMC(fModelCrystalRa226S1M1, fSmearCrystalRa226S1M1, dRes1 , dRes2);
  SmearMC(fModelCrystalRa226S10M1, fSmearCrystalRa226S10M1, dRes1 , dRes2);
  SmearMC(fModelCrystalRa226S100M1, fSmearCrystalRa226S100M1, dRes1 , dRes2);

  SmearMC(fModelFramePb210BM1, fSmearFramePb210BM1, dRes1 , dRes2);
  SmearMC(fModelFramePb210S01M1, fSmearFramePb210S01M1, dRes1 , dRes2);
  SmearMC(fModelFramePb210S1M1, fSmearFramePb210S1M1, dRes1 , dRes2);
  SmearMC(fModelFramePb210S10M1, fSmearFramePb210S10M1, dRes1 , dRes2);
  SmearMC(fModelFramePb210S100M1, fSmearFramePb210S100M1, dRes1 , dRes2);

  // SmearMC(fModelFrameRa226BM1, fSmearFrameRa226BM1, dRes1 , dRes2);
  SmearMC(fModelFrameRa226S01M1, fSmearFrameRa226S01M1, dRes1 , dRes2);
  SmearMC(fModelFrameRa226S1M1, fSmearFrameRa226S1M1, dRes1 , dRes2);
  SmearMC(fModelFrameRa226S10M1, fSmearFrameRa226S10M1, dRes1 , dRes2);
  SmearMC(fModelFrameRa226S100M1, fSmearFrameRa226S100M1, dRes1 , dRes2);

  //M2
  SmearMC(fModelCrystalPt190BM2, fSmearCrystalPt190BM2, dRes1 , dRes2);
  // SmearMC(fModelCrystalPt190S01M2, fSmearCrystalPt190S01M2, dRes1 , dRes2);
  SmearMC(fModelCrystalPt190S1M2, fSmearCrystalPt190S1M2, dRes1 , dRes2);
  SmearMC(fModelCrystalPt190S10M2, fSmearCrystalPt190S10M2, dRes1 , dRes2);
  SmearMC(fModelCrystalPt190S100M2, fSmearCrystalPt190S100M2, dRes1 , dRes2);

  SmearMC(fModelCrystalPb210BM2, fSmearCrystalPb210BM2, dRes1 , dRes2);
  SmearMC(fModelCrystalPb210S01M2, fSmearCrystalPb210S01M2, dRes1 , dRes2);
  SmearMC(fModelCrystalPb210S1M2, fSmearCrystalPb210S1M2, dRes1 , dRes2);
  SmearMC(fModelCrystalPb210S10M2, fSmearCrystalPb210S10M2, dRes1 , dRes2);
  SmearMC(fModelCrystalPb210S100M2, fSmearCrystalPb210S100M2, dRes1 , dRes2);

  SmearMC(fModelCrystalRa226BM2, fSmearCrystalRa226BM2, dRes1 , dRes2);
  SmearMC(fModelCrystalRa226S01M2, fSmearCrystalRa226S01M2, dRes1 , dRes2);
  SmearMC(fModelCrystalRa226S1M2, fSmearCrystalRa226S1M2, dRes1 , dRes2);
  SmearMC(fModelCrystalRa226S10M2, fSmearCrystalRa226S10M2, dRes1 , dRes2);
  SmearMC(fModelCrystalRa226S100M2, fSmearCrystalRa226S100M2, dRes1 , dRes2);

  SmearMC(fModelFramePb210BM2, fSmearFramePb210BM2, dRes1 , dRes2);
  SmearMC(fModelFramePb210S01M2, fSmearFramePb210S01M2, dRes1 , dRes2);
  SmearMC(fModelFramePb210S1M2, fSmearFramePb210S1M2, dRes1 , dRes2);
  SmearMC(fModelFramePb210S10M2, fSmearFramePb210S10M2, dRes1 , dRes2);
  SmearMC(fModelFramePb210S100M2, fSmearFramePb210S100M2, dRes1 , dRes2);

  // SmearMC(fModelFrameRa226BM2, fSmearFrameRa226BM2, dRes1 , dRes2);
  SmearMC(fModelFrameRa226S01M2, fSmearFrameRa226S01M2, dRes1 , dRes2);
  SmearMC(fModelFrameRa226S1M2, fSmearFrameRa226S1M2, dRes1 , dRes2);
  SmearMC(fModelFrameRa226S10M2, fSmearFrameRa226S10M2, dRes1 , dRes2);
  SmearMC(fModelFrameRa226S100M2, fSmearFrameRa226S100M2, dRes1 , dRes2);

    cout << "Finished smearing MC histograms" << endl;

  fSmearCrystalPt190BM1->Rebin(dAdaptiveBinsM1, "fnewCrystalPt190BM1", dAdaptiveArrayM1);
  // fSmearCrystalPt190S01M1->Rebin(dAdaptiveBinsM1, "fnewCrystalPt190S01M1", dAdaptiveArrayM1);
  fSmearCrystalPt190S1M1->Rebin(dAdaptiveBinsM1, "fnewCrystalPt190S1M1", dAdaptiveArrayM1);
  fSmearCrystalPt190S10M1->Rebin(dAdaptiveBinsM1, "fnewCrystalPt190S10M1", dAdaptiveArrayM1);
  fSmearCrystalPt190S100M1->Rebin(dAdaptiveBinsM1, "fnewCrystalPt190S100M1", dAdaptiveArrayM1);

  fSmearCrystalPb210BM1->Rebin(dAdaptiveBinsM1, "fnewCrystalPb210BM1", dAdaptiveArrayM1);
  fSmearCrystalPb210S01M1->Rebin(dAdaptiveBinsM1, "fnewCrystalPb210S01M1", dAdaptiveArrayM1);
  fSmearCrystalPb210S1M1->Rebin(dAdaptiveBinsM1, "fnewCrystalPb210S1M1", dAdaptiveArrayM1);
  fSmearCrystalPb210S10M1->Rebin(dAdaptiveBinsM1, "fnewCrystalPb210S10M1", dAdaptiveArrayM1);
  fSmearCrystalPb210S100M1->Rebin(dAdaptiveBinsM1, "fnewCrystalPb210S100M1", dAdaptiveArrayM1);

  fSmearCrystalRa226BM1->Rebin(dAdaptiveBinsM1, "fnewCrystalRa226BM1", dAdaptiveArrayM1);
  fSmearCrystalRa226S01M1->Rebin(dAdaptiveBinsM1, "fnewCrystalRa226S01M1", dAdaptiveArrayM1);
  fSmearCrystalRa226S1M1->Rebin(dAdaptiveBinsM1, "fnewCrystalRa226S1M1", dAdaptiveArrayM1);
  fSmearCrystalRa226S10M1->Rebin(dAdaptiveBinsM1, "fnewCrystalRa226S10M1", dAdaptiveArrayM1);
  fSmearCrystalRa226S100M1->Rebin(dAdaptiveBinsM1, "fnewCrystalRa226S100M1", dAdaptiveArrayM1);

  fSmearFramePb210BM1->Rebin(dAdaptiveBinsM1, "fnewFramePb210BM1", dAdaptiveArrayM1);
  fSmearFramePb210S01M1->Rebin(dAdaptiveBinsM1, "fnewFramePb210S01M1", dAdaptiveArrayM1);
  fSmearFramePb210S1M1->Rebin(dAdaptiveBinsM1, "fnewFramePb210S1M1", dAdaptiveArrayM1);
  fSmearFramePb210S10M1->Rebin(dAdaptiveBinsM1, "fnewFramePb210S10M1", dAdaptiveArrayM1);
  fSmearFramePb210S100M1->Rebin(dAdaptiveBinsM1, "fnewFramePb210S100M1", dAdaptiveArrayM1);

  // fSmearFrameRa226BM1->Rebin(dAdaptiveBinsM1, "fnewFrameRa226BM1", dAdaptiveArrayM1);
  fSmearFrameRa226S01M1->Rebin(dAdaptiveBinsM1, "fnewFrameRa226S01M1", dAdaptiveArrayM1);
  fSmearFrameRa226S1M1->Rebin(dAdaptiveBinsM1, "fnewFrameRa226S1M1", dAdaptiveArrayM1);
  fSmearFrameRa226S10M1->Rebin(dAdaptiveBinsM1, "fnewFrameRa226S10M1", dAdaptiveArrayM1);
  fSmearFrameRa226S100M1->Rebin(dAdaptiveBinsM1, "fnewFrameRa226S100M1", dAdaptiveArrayM1);  



  fSmearCrystalPt190BM2->Rebin(dAdaptiveBinsM2, "fnewCrystalPt190BM2", dAdaptiveArrayM2);
  // fSmearCrystalPt190S01M2->Rebin(dAdaptiveBinsM2, "fnewCrystalPt190S01M2", dAdaptiveArrayM2);
  fSmearCrystalPt190S1M2->Rebin(dAdaptiveBinsM2, "fnewCrystalPt190S1M2", dAdaptiveArrayM2);
  fSmearCrystalPt190S10M2->Rebin(dAdaptiveBinsM2, "fnewCrystalPt190S10M2", dAdaptiveArrayM2);
  fSmearCrystalPt190S100M2->Rebin(dAdaptiveBinsM2, "fnewCrystalPt190S100M2", dAdaptiveArrayM2);

  fSmearCrystalPb210BM2->Rebin(dAdaptiveBinsM2, "fnewCrystalPb210BM2", dAdaptiveArrayM2);
  fSmearCrystalPb210S01M2->Rebin(dAdaptiveBinsM2, "fnewCrystalPb210S01M2", dAdaptiveArrayM2);
  fSmearCrystalPb210S1M2->Rebin(dAdaptiveBinsM2, "fnewCrystalPb210S1M2", dAdaptiveArrayM2);
  fSmearCrystalPb210S10M2->Rebin(dAdaptiveBinsM2, "fnewCrystalPb210S10M2", dAdaptiveArrayM2);
  fSmearCrystalPb210S100M2->Rebin(dAdaptiveBinsM2, "fnewCrystalPb210S100M2", dAdaptiveArrayM2);

  fSmearCrystalRa226BM2->Rebin(dAdaptiveBinsM2, "fnewCrystalRa226BM2", dAdaptiveArrayM2);
  fSmearCrystalRa226S01M2->Rebin(dAdaptiveBinsM2, "fnewCrystalRa226S01M2", dAdaptiveArrayM2);
  fSmearCrystalRa226S1M2->Rebin(dAdaptiveBinsM2, "fnewCrystalRa226S1M2", dAdaptiveArrayM2);
  fSmearCrystalRa226S10M2->Rebin(dAdaptiveBinsM2, "fnewCrystalRa226S10M2", dAdaptiveArrayM2);
  fSmearCrystalRa226S100M2->Rebin(dAdaptiveBinsM2, "fnewCrystalRa226S100M2", dAdaptiveArrayM2);

  fSmearFramePb210BM2->Rebin(dAdaptiveBinsM2, "fnewFramePb210BM2", dAdaptiveArrayM2);
  fSmearFramePb210S01M2->Rebin(dAdaptiveBinsM2, "fnewFramePb210S01M2", dAdaptiveArrayM2);
  fSmearFramePb210S1M2->Rebin(dAdaptiveBinsM2, "fnewFramePb210S1M2", dAdaptiveArrayM2);
  fSmearFramePb210S10M2->Rebin(dAdaptiveBinsM2, "fnewFramePb210S10M2", dAdaptiveArrayM2);
  fSmearFramePb210S100M2->Rebin(dAdaptiveBinsM2, "fnewFramePb210S100M2", dAdaptiveArrayM2);

  // fSmearFrameRa226BM2->Rebin(dAdaptiveBinsM2, "fnewFrameRa226BM2", dAdaptiveArrayM2);
  fSmearFrameRa226S01M2->Rebin(dAdaptiveBinsM2, "fnewFrameRa226S01M2", dAdaptiveArrayM2);
  fSmearFrameRa226S1M2->Rebin(dAdaptiveBinsM2, "fnewFrameRa226S1M2", dAdaptiveArrayM2);
  fSmearFrameRa226S10M2->Rebin(dAdaptiveBinsM2, "fnewFrameRa226S10M2", dAdaptiveArrayM2);
  fSmearFrameRa226S100M2->Rebin(dAdaptiveBinsM2, "fnewFrameRa226S100M2", dAdaptiveArrayM2); 

  for(int i = 1; i <= dAdaptiveBinsM1; i++)
  {
    fAdapCrystalPt190BM1->SetBinContent(i, 2*fnewCrystalPt190BM1->GetBinContent(i)/fnewCrystalPt190BM1->GetBinWidth(i));
    // fAdapCrystalPt190S01M1->SetBinContent(i, 2*fnewCrystalPt190S01M1->GetBinContent(i)/fnewCrystalPt190S01M1->GetBinWidth(i));
    fAdapCrystalPt190S1M1->SetBinContent(i, 2*fnewCrystalPt190S1M1->GetBinContent(i)/fnewCrystalPt190S1M1->GetBinWidth(i));
    fAdapCrystalPt190S10M1->SetBinContent(i, 2*fnewCrystalPt190S10M1->GetBinContent(i)/fnewCrystalPt190S10M1->GetBinWidth(i));
    fAdapCrystalPt190S100M1->SetBinContent(i, 2*fnewCrystalPt190S100M1->GetBinContent(i)/fnewCrystalPt190S100M1->GetBinWidth(i));

    fAdapCrystalPb210BM1->SetBinContent(i, 2*fnewCrystalPb210BM1->GetBinContent(i)/fnewCrystalPb210BM1->GetBinWidth(i));
    fAdapCrystalPb210S01M1->SetBinContent(i, 2*fnewCrystalPb210S01M1->GetBinContent(i)/fnewCrystalPb210S01M1->GetBinWidth(i));
    fAdapCrystalPb210S1M1->SetBinContent(i, 2*fnewCrystalPb210S1M1->GetBinContent(i)/fnewCrystalPb210S1M1->GetBinWidth(i));
    fAdapCrystalPb210S10M1->SetBinContent(i, 2*fnewCrystalPb210S10M1->GetBinContent(i)/fnewCrystalPb210S10M1->GetBinWidth(i));
    fAdapCrystalPb210S100M1->SetBinContent(i, 2*fnewCrystalPb210S100M1->GetBinContent(i)/fnewCrystalPb210S100M1->GetBinWidth(i));    

    fAdapCrystalRa226BM1->SetBinContent(i, 2*fnewCrystalRa226BM1->GetBinContent(i)/fnewCrystalRa226BM1->GetBinWidth(i));
    fAdapCrystalRa226S01M1->SetBinContent(i, 2*fnewCrystalRa226S01M1->GetBinContent(i)/fnewCrystalRa226S01M1->GetBinWidth(i));
    fAdapCrystalRa226S1M1->SetBinContent(i, 2*fnewCrystalRa226S1M1->GetBinContent(i)/fnewCrystalRa226S1M1->GetBinWidth(i));
    fAdapCrystalRa226S10M1->SetBinContent(i, 2*fnewCrystalRa226S10M1->GetBinContent(i)/fnewCrystalRa226S10M1->GetBinWidth(i));
    fAdapCrystalRa226S100M1->SetBinContent(i, 2*fnewCrystalRa226S100M1->GetBinContent(i)/fnewCrystalRa226S100M1->GetBinWidth(i));   

    fAdapFramePb210BM1->SetBinContent(i, 2*fnewFramePb210BM1->GetBinContent(i)/fnewFramePb210BM1->GetBinWidth(i));
    fAdapFramePb210S01M1->SetBinContent(i, 2*fnewFramePb210S01M1->GetBinContent(i)/fnewFramePb210S01M1->GetBinWidth(i));
    fAdapFramePb210S1M1->SetBinContent(i, 2*fnewFramePb210S1M1->GetBinContent(i)/fnewFramePb210S1M1->GetBinWidth(i));
    fAdapFramePb210S10M1->SetBinContent(i, 2*fnewFramePb210S10M1->GetBinContent(i)/fnewFramePb210S10M1->GetBinWidth(i));
    fAdapFramePb210S100M1->SetBinContent(i, 2*fnewFramePb210S100M1->GetBinContent(i)/fnewFramePb210S100M1->GetBinWidth(i));    

    // fAdapFrameRa226BM1->SetBinContent(i, 2*fnewFrameRa226BM1->GetBinContent(i)/fnewFrameRa226BM1->GetBinWidth(i));
    fAdapFrameRa226S01M1->SetBinContent(i, 2*fnewFrameRa226S01M1->GetBinContent(i)/fnewFrameRa226S01M1->GetBinWidth(i));
    fAdapFrameRa226S1M1->SetBinContent(i, 2*fnewFrameRa226S1M1->GetBinContent(i)/fnewFrameRa226S1M1->GetBinWidth(i));
    fAdapFrameRa226S10M1->SetBinContent(i, 2*fnewFrameRa226S10M1->GetBinContent(i)/fnewFrameRa226S10M1->GetBinWidth(i));
    fAdapFrameRa226S100M1->SetBinContent(i, 2*fnewFrameRa226S100M1->GetBinContent(i)/fnewFrameRa226S100M1->GetBinWidth(i)); 
  }

  for(int i = 1; i <= dAdaptiveBinsM2; i++)
  {
    fAdapCrystalPt190BM2->SetBinContent(i, 2*fnewCrystalPt190BM2->GetBinContent(i)/fnewCrystalPt190BM2->GetBinWidth(i));
    // fAdapCrystalPt190S01M2->SetBinContent(i, 2*fnewCrystalPt190S01M2->GetBinContent(i)/fnewCrystalPt190S01M2->GetBinWidth(i));
    fAdapCrystalPt190S1M2->SetBinContent(i, 2*fnewCrystalPt190S1M2->GetBinContent(i)/fnewCrystalPt190S1M2->GetBinWidth(i));
    fAdapCrystalPt190S10M2->SetBinContent(i, 2*fnewCrystalPt190S10M2->GetBinContent(i)/fnewCrystalPt190S10M2->GetBinWidth(i));
    fAdapCrystalPt190S100M2->SetBinContent(i, 2*fnewCrystalPt190S100M2->GetBinContent(i)/fnewCrystalPt190S100M2->GetBinWidth(i));

    fAdapCrystalPb210BM2->SetBinContent(i, 2*fnewCrystalPb210BM2->GetBinContent(i)/fnewCrystalPb210BM2->GetBinWidth(i));
    fAdapCrystalPb210S01M2->SetBinContent(i, 2*fnewCrystalPb210S01M2->GetBinContent(i)/fnewCrystalPb210S01M2->GetBinWidth(i));
    fAdapCrystalPb210S1M2->SetBinContent(i, 2*fnewCrystalPb210S1M2->GetBinContent(i)/fnewCrystalPb210S1M2->GetBinWidth(i));
    fAdapCrystalPb210S10M2->SetBinContent(i, 2*fnewCrystalPb210S10M2->GetBinContent(i)/fnewCrystalPb210S10M2->GetBinWidth(i));
    fAdapCrystalPb210S100M2->SetBinContent(i, 2*fnewCrystalPb210S100M2->GetBinContent(i)/fnewCrystalPb210S100M2->GetBinWidth(i));    

    fAdapCrystalRa226BM2->SetBinContent(i, 2*fnewCrystalRa226BM2->GetBinContent(i)/fnewCrystalRa226BM2->GetBinWidth(i));
    fAdapCrystalRa226S01M2->SetBinContent(i, 2*fnewCrystalRa226S01M2->GetBinContent(i)/fnewCrystalRa226S01M2->GetBinWidth(i));
    fAdapCrystalRa226S1M2->SetBinContent(i, 2*fnewCrystalRa226S1M2->GetBinContent(i)/fnewCrystalRa226S1M2->GetBinWidth(i));
    fAdapCrystalRa226S10M2->SetBinContent(i, 2*fnewCrystalRa226S10M2->GetBinContent(i)/fnewCrystalRa226S10M2->GetBinWidth(i));
    fAdapCrystalRa226S100M2->SetBinContent(i, 2*fnewCrystalRa226S100M2->GetBinContent(i)/fnewCrystalRa226S100M2->GetBinWidth(i));   

    fAdapFramePb210BM2->SetBinContent(i, 2*fnewFramePb210BM2->GetBinContent(i)/fnewFramePb210BM2->GetBinWidth(i));
    fAdapFramePb210S01M2->SetBinContent(i, 2*fnewFramePb210S01M2->GetBinContent(i)/fnewFramePb210S01M2->GetBinWidth(i));
    fAdapFramePb210S1M2->SetBinContent(i, 2*fnewFramePb210S1M2->GetBinContent(i)/fnewFramePb210S1M2->GetBinWidth(i));
    fAdapFramePb210S10M2->SetBinContent(i, 2*fnewFramePb210S10M2->GetBinContent(i)/fnewFramePb210S10M2->GetBinWidth(i));
    fAdapFramePb210S100M2->SetBinContent(i, 2*fnewFramePb210S100M2->GetBinContent(i)/fnewFramePb210S100M2->GetBinWidth(i));    

    // fAdapFrameRa226BM2->SetBinContent(i, 2*fnewFrameRa226BM2->GetBinContent(i)/fnewFrameRa226BM2->GetBinWidth(i));
    fAdapFrameRa226S01M2->SetBinContent(i, 2*fnewFrameRa226S01M2->GetBinContent(i)/fnewFrameRa226S01M2->GetBinWidth(i));
    fAdapFrameRa226S1M2->SetBinContent(i, 2*fnewFrameRa226S1M2->GetBinContent(i)/fnewFrameRa226S1M2->GetBinWidth(i));
    fAdapFrameRa226S10M2->SetBinContent(i, 2*fnewFrameRa226S10M2->GetBinContent(i)/fnewFrameRa226S10M2->GetBinWidth(i));
    fAdapFrameRa226S100M2->SetBinContent(i, 2*fnewFrameRa226S100M2->GetBinContent(i)/fnewFrameRa226S100M2->GetBinWidth(i)); 
  }

}


// Loads MC files into Trees... Careful about coincidence timing!
TChain *TBackgroundModel::LoadMC(std::string dDir, std::string dLocation, std::string dSource, std::string dSType, int dMult)
{
    TChain *outTree = new TChain("outTree");
    outTree->Add(Form("%s%s-%s-%s-M%d-T50-r0.0425.root", dDir.c_str(), dLocation.c_str(), dSource.c_str(), dSType.c_str(), dMult));
    // outTree->Add(Form("%s%s-%s-%s-M%d-r0.0100.root", dDir.c_str(), dLocation.c_str(), dSource.c_str(), dSType.c_str(), dMult));

    return outTree;
}


// Normalize histogram
void TBackgroundModel::NormalizePDF(TH1D *h1, int minE, int maxE)
{
	double dIntegral = 0;
	double Time = 0;

	// hChain->SetBranchAddress("Time", &Time);

	// int dEvents = hChain->GetEntries();
	// hChain->GetEntry(dEvents - 1);

	// bin 0 = underflow, bin dNBins = last bin with upper-edge xup Excluded
	dIntegral = h1->Integral(minE/dBinSize, maxE/dBinSize);
	// cout << "Integral for " << h1->GetTitle() << " :" << dIntegral << endl;


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

void TBackgroundModel::NormalizePDFPair(TH1D *h1, TH1D *h2, int minE, int maxE)
{
  double dIntegral = 0;

  // bin 0 = underflow, bin dNBins = last bin with upper-edge xup Excluded
  dIntegral = h1->Integral(minE/dBinSize, maxE/dBinSize);
  // cout << "Integral for " << h1->GetTitle() << " :" << dIntegral << endl;


  // Make sure integral isn't 0
  // If it is 0, clear model... 
  if(dIntegral == 0)
  {
    cout << Form("Integral of %s is 0, resetting histogram", h1->GetName()) << endl;
    h1->Reset();
    h2->Reset();
  }

  if(dIntegral != 0)
  {
    h1->Scale(1.0/dIntegral);
    h2->Scale(1.0/dIntegral);
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
  cout<< "Par12 = "  << fParameters[12] << " +/- " << fParError[12] << endl;
  cout<< "Par13 = "  << fParameters[13] << " +/- " << fParError[13] << endl;
  cout<< "Par14 = "  << fParameters[14] << " +/- " << fParError[14] << endl;
  cout<< "Par15 = "  << fParameters[15] << " +/- " << fParError[15] << endl;
  cout<< "Par16 = "  << fParameters[16] << " +/- " << fParError[16] << endl;
  cout<< "Par17 = "  << fParameters[17] << " +/- " << fParError[17] << endl;
  cout<< "Par18 = "  << fParameters[18] << " +/- " << fParError[18] << endl;
  cout<< "Par19 = "  << fParameters[19] << " +/- " << fParError[19] << endl;
  cout<< "Par20 = "  << fParameters[20] << " +/- " << fParError[20] << endl;
  cout<< "Par21 = "  << fParameters[21] << " +/- " << fParError[21] << endl;
  cout<< "Par22 = "  << fParameters[22] << " +/- " << fParError[22] << endl;
  cout<< "Par23 = "  << fParameters[23] << " +/- " << fParError[23] << endl;
  cout<< "Par24 = "  << fParameters[24] << " +/- " << fParError[24] << endl;
  cout<< "Par25 = "  << fParameters[25] << " +/- " << fParError[25] << endl;
}


// Set Parameters in Model
void TBackgroundModel::SetParameters(int index, double value)
{
	// Change the index max depending on model
	if(index > 26) cout << "Index too large" << endl;
	else fParameters[index] = value;
}


// For custom Smearing with resolution, currently constant resolution 
TH1D *TBackgroundModel::SmearMC(TH1D *hMC, TH1D *hSMC, double resolution1, double resolution2)
{
	// Reset previously Modeled histogram
	hSMC->Reset();

	double dNorm;
  double dNorm2;
	double dSmearedValue;

	for(int i = 0; i<dNBins; i++)
	{
		for(int j = 0; j<dNBins; j++)
		{
			// Normalization of gaussian = (bsin size * Area of bin j in MC) / Sigma of bin j (fit function evaluated at bin center)
			dNorm = (853.3/1215.8)*dBinSize*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*resolution1);
      dNorm2 = (362.5/1215.8)*dBinSize*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*resolution2);

			// Set parameters of gaussian ... 2nd gaussian *slightly* shifted... not sure if this works
			gaus->SetParameters(dNorm, hMC->GetBinCenter(j), resolution1);
      gaus2->SetParameters(dNorm2, hMC->GetBinCenter(j)-1, resolution2);

			// Smeared contribution from gaussian centered at bin j for bin i 
			dSmearedValue = gaus->Eval(hSMC->GetBinCenter(i)) + gaus2->Eval(hSMC->GetBinCenter(i));

			// Fill bin i with contribution from gaussian centered at bin j
			hSMC->Fill(hSMC->GetBinCenter(i), dSmearedValue);
		}
	}
	
	return hSMC;
}

// Smearing with single gaussian
TH1D *TBackgroundModel::SmearMCOld(TH1D *hMC, TH1D *hSMC, double resolution1)
{
  // Reset previously Modeled histogram
  hSMC->Reset();

  double dNorm;
  double dSmearedValue;

  for(int i = 0; i<dNBins; i++)
  {
    for(int j = 0; j<dNBins; j++)
    {
      // Normalization of gaussian = (bsin size * Area of bin j in MC) / Sigma of bin j (fit function evaluated at bin center)
      dNorm = dBinSize*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*resolution1);

      // Set parameters of gaussian ... 2nd gaussian *slightly* shifted... not sure if this works
      gaus->SetParameters(dNorm, hMC->GetBinCenter(j), resolution1);

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

  /////////////////////////////////////////////////
  //// Few Parameters //////////////// Smeared 
  ////////////////////////////////////////////////
  /*
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
  fModelTotM1->Add( fModel2NDBDM1,     fParameters[11]);  
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
  fModelTotM2->Add( fModel2NDBDM2,     fParameters[11]);  
  fModelTotM2->Add( fModelBiM2,        fParameters[10]);  
*/




  ///////////////////////////////////////////////////
  //// Few Parameters //////////////// Unsmeared
  //////////////////////////////////////////////////
  // M1
/*
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
  // fModelTotM1->Add( fSmearOVCRaM1,      fParameters[3]);

  // fModelTotM1->Add( fSmearFrameKM1,    fParameters[4]);
  fModelTotM1->Add( fSmearTShieldKM1,  fParameters[4]);
  fModelTotM1->Add( fSmear50mKKM1,     fParameters[4]);
  fModelTotM1->Add( fSmear600mKKM1,    fParameters[4]);
  fModelTotM1->Add( fSmearIVCKM1,      fParameters[5]);
  // fModelTotM1->Add( fSmearOVCKM1,      fParameters[5]); 

  // fModelTotM1->Add( fSmearFrameCoM1,    fParameters[6]);
  fModelTotM1->Add( fSmearTShieldCoM1,  fParameters[6]);
  fModelTotM1->Add( fSmear50mKCoM1,     fParameters[6]);
  fModelTotM1->Add( fSmear600mKCoM1,    fParameters[6]);
  fModelTotM1->Add( fSmearIVCCoM1,      fParameters[7]);
  // fModelTotM1->Add( fSmearOVCCoM1,      fParameters[7]);  

  fModelTotM1->Add( fSmearTShieldMnM1,  fParameters[11]);
  fModelTotM1->Add( fSmearIVCMnM1,      fParameters[12]);

  fModelTotM1->Add( fSmearNDBDM1,      fParameters[9]);  
  fModelTotM1->Add( fSmear2NDBDM1,     fParameters[8]);  
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
  // fModelTotM2->Add( fSmearOVCRaM2,      fParameters[3]);

  // fModelTotM2->Add( fSmearFrameKM2,    fParameters[4]);
  fModelTotM2->Add( fSmearTShieldKM2,  fParameters[4]);
  fModelTotM2->Add( fSmear50mKKM2,     fParameters[4]);
  fModelTotM2->Add( fSmear600mKKM2,    fParameters[4]);
  fModelTotM2->Add( fSmearIVCKM2,      fParameters[5]);
  // fModelTotM2->Add( fSmearOVCKM2,      fParameters[5]); 

  // fModelTotM2->Add( fSmearFrameCoM2,    fParameters[6]);
  fModelTotM2->Add( fSmearTShieldCoM2,  fParameters[6]);
  fModelTotM2->Add( fSmear50mKCoM2,     fParameters[6]);
  fModelTotM2->Add( fSmear600mKCoM2,    fParameters[6]);
  fModelTotM2->Add( fSmearIVCCoM2,      fParameters[7]);
  // fModelTotM2->Add( fSmearOVCCoM2,      fParameters[7]);  

  fModelTotM2->Add( fSmearTShieldMnM2,  fParameters[11]);
  fModelTotM2->Add( fSmearIVCMnM2,      fParameters[12]);

  fModelTotM2->Add( fSmearNDBDM2,      fParameters[9]);  
  fModelTotM2->Add( fSmear2NDBDM2,     fParameters[8]);  
  fModelTotM2->Add( fSmearBiM2,        fParameters[10]);  
*/


  /////////////////////////////////////
  //// Many parameters
  ////////////////////////////////////

  // M1
  fModelTotM1->Add( fSmearFrameThM1,    fParameters[0]);
  fModelTotM1->Add( fSmearTShieldThM1,  fParameters[1]);  
  fModelTotM1->Add( fSmear50mKThM1,     fParameters[13]);
  fModelTotM1->Add( fSmear600mKThM1,    fParameters[14]);
  fModelTotM1->Add( fSmearIVCThM1,      fParameters[15]);
  fModelTotM1->Add( fSmearOVCThM1,      fParameters[16]);

  fModelTotM1->Add( fSmearFrameRaM1,    fParameters[2]);
  fModelTotM1->Add( fSmearTShieldRaM1,  fParameters[3]);  
  fModelTotM1->Add( fSmear50mKRaM1,     fParameters[17]);
  fModelTotM1->Add( fSmear600mKRaM1,    fParameters[18]);
  fModelTotM1->Add( fSmearIVCRaM1,      fParameters[19]);
  fModelTotM1->Add( fSmearOVCRaM1,      fParameters[20]);

  fModelTotM1->Add( fSmearFrameKM1,     fParameters[4]);
  fModelTotM1->Add( fSmearTShieldKM1,   fParameters[4]);
  fModelTotM1->Add( fSmear50mKKM1,      fParameters[4]);
  fModelTotM1->Add( fSmear600mKKM1,     fParameters[4]);
  fModelTotM1->Add( fSmearIVCKM1,       fParameters[5]);
  fModelTotM1->Add( fSmearOVCKM1,       fParameters[5]); 

  fModelTotM1->Add( fSmearFrameCoM1,    fParameters[6]);
  fModelTotM1->Add( fSmearTShieldCoM1,  fParameters[7]);
  fModelTotM1->Add( fSmear50mKCoM1,     fParameters[21]);
  fModelTotM1->Add( fSmear600mKCoM1,    fParameters[22]);
  fModelTotM1->Add( fSmearIVCCoM1,      fParameters[23]);
  fModelTotM1->Add( fSmearOVCCoM1,      fParameters[24]);  

  fModelTotM1->Add( fSmearTShieldMnM1,  fParameters[11]);
  fModelTotM1->Add( fSmearIVCMnM1,      fParameters[12]);

  fModelTotM1->Add( fSmearNDBDM1,      fParameters[9]);  
  fModelTotM1->Add( fSmear2NDBDM1,     fParameters[8]);  
  fModelTotM1->Add( fSmearBiM1,        fParameters[10]);  

  fModelTotM1->Add( fSmearCrystalBi2M1, fParameters[25]);  

  // M2
  fModelTotM2->Add( fSmearFrameThM2,    fParameters[0]);
  fModelTotM2->Add( fSmearTShieldThM2,  fParameters[1]);  
  fModelTotM2->Add( fSmear50mKThM2,     fParameters[13]);
  fModelTotM2->Add( fSmear600mKThM2,    fParameters[14]);
  fModelTotM2->Add( fSmearIVCThM2,      fParameters[15]);
  fModelTotM2->Add( fSmearOVCThM2,      fParameters[16]);

  fModelTotM2->Add( fSmearFrameRaM2,    fParameters[2]);
  fModelTotM2->Add( fSmearTShieldRaM2,  fParameters[3]);  
  fModelTotM2->Add( fSmear50mKRaM2,     fParameters[17]);
  fModelTotM2->Add( fSmear600mKRaM2,    fParameters[18]);
  fModelTotM2->Add( fSmearIVCRaM2,      fParameters[19]);
  fModelTotM2->Add( fSmearOVCRaM2,      fParameters[20]);

  fModelTotM2->Add( fSmearFrameKM2,     fParameters[4]);
  fModelTotM2->Add( fSmearTShieldKM2,   fParameters[4]);
  fModelTotM2->Add( fSmear50mKKM2,      fParameters[4]);
  fModelTotM2->Add( fSmear600mKKM2,     fParameters[4]);
  fModelTotM2->Add( fSmearIVCKM2,       fParameters[5]);
  fModelTotM2->Add( fSmearOVCKM2,       fParameters[5]); 

  fModelTotM2->Add( fSmearFrameCoM2,    fParameters[6]);
  fModelTotM2->Add( fSmearTShieldCoM2,  fParameters[7]);
  fModelTotM2->Add( fSmear50mKCoM2,     fParameters[21]);
  fModelTotM2->Add( fSmear600mKCoM2,    fParameters[22]);
  fModelTotM2->Add( fSmearIVCCoM2,      fParameters[23]);
  fModelTotM2->Add( fSmearOVCCoM2,      fParameters[24]);  

  fModelTotM2->Add( fSmearTShieldMnM2,  fParameters[11]);
  fModelTotM2->Add( fSmearIVCMnM2,      fParameters[12]);

  fModelTotM2->Add( fSmearNDBDM2,      fParameters[9]);  
  fModelTotM2->Add( fSmear2NDBDM2,     fParameters[8]);  
  fModelTotM2->Add( fSmearBiM2,        fParameters[10]);  

  fModelTotM1->Add( fSmearCrystalBi2M2 , fParameters[25]);  


  // Adding on correction for M2.. (just the M1 spectrum)
  fModelTotM2->Add( fCorrectionM2, 180197*(1-TMath::Exp(-2*0.05*0.1)) );



}

void TBackgroundModel::UpdateModelAdaptive()
{
  if(fModelTotM1 == NULL) 
  {
    cout << "Model Histogram Not Created" << endl;
    return;
  }

  // Reset all bins in model histogram(s)
  fModelTotAdapM1->Reset();
  fModelTotAdapM2->Reset();

  dNumCalls++;

  // Create model
  if(dNumCalls%50==0)
  {
    cout << "Call #: "<< dNumCalls << endl;
  }

  /////////////////////////////////////
  //// Adaptive Binning parameters
  ////////////////////////////////////

  // M1
  fModelTotAdapM1->Add( fAdapFrameThM1,    fParameters[0]);
  fModelTotAdapM1->Add( fAdapTShieldThM1,  fParameters[1]);  
  fModelTotAdapM1->Add( fAdap50mKThM1,     fParameters[13]);
  fModelTotAdapM1->Add( fAdap600mKThM1,    fParameters[14]);
  fModelTotAdapM1->Add( fAdapIVCThM1,      fParameters[15]);
  fModelTotAdapM1->Add( fAdapOVCThM1,      fParameters[16]);

  fModelTotAdapM1->Add( fAdapFrameRaM1,    fParameters[2]);
  fModelTotAdapM1->Add( fAdapTShieldRaM1,  fParameters[3]);  
  fModelTotAdapM1->Add( fAdap50mKRaM1,     fParameters[17]);
  fModelTotAdapM1->Add( fAdap600mKRaM1,    fParameters[18]);
  fModelTotAdapM1->Add( fAdapIVCRaM1,      fParameters[19]);
  fModelTotAdapM1->Add( fAdapOVCRaM1,      fParameters[20]);

  fModelTotAdapM1->Add( fAdapFrameKM1,     fParameters[4]);
  fModelTotAdapM1->Add( fAdapTShieldKM1,   fParameters[4]);
  fModelTotAdapM1->Add( fAdap50mKKM1,      fParameters[4]);
  fModelTotAdapM1->Add( fAdap600mKKM1,     fParameters[4]);
  fModelTotAdapM1->Add( fAdapIVCKM1,       fParameters[5]);
  fModelTotAdapM1->Add( fAdapOVCKM1,       fParameters[5]); 

  fModelTotAdapM1->Add( fAdapFrameCoM1,    fParameters[6]);
  fModelTotAdapM1->Add( fAdapTShieldCoM1,  fParameters[7]);
  fModelTotAdapM1->Add( fAdap50mKCoM1,     fParameters[21]);
  fModelTotAdapM1->Add( fAdap600mKCoM1,    fParameters[22]);
  fModelTotAdapM1->Add( fAdapIVCCoM1,      fParameters[23]);
  fModelTotAdapM1->Add( fAdapOVCCoM1,      fParameters[24]);  

  fModelTotAdapM1->Add( fAdapTShieldMnM1,  fParameters[11]);
  fModelTotAdapM1->Add( fAdapIVCMnM1,      fParameters[12]);

  fModelTotAdapM1->Add( fAdapNDBDM1,      fParameters[9]);  
  fModelTotAdapM1->Add( fAdap2NDBDM1,     fParameters[8]);  
  fModelTotAdapM1->Add( fAdapBiM1,        fParameters[10]);  

  // M2
  fModelTotAdapM2->Add( fAdapFrameThM2,    fParameters[0]);
  fModelTotAdapM2->Add( fAdapTShieldThM2,  fParameters[1]);  
  fModelTotAdapM2->Add( fAdap50mKThM2,     fParameters[13]);
  fModelTotAdapM2->Add( fAdap600mKThM2,    fParameters[14]);
  fModelTotAdapM2->Add( fAdapIVCThM2,      fParameters[15]);
  fModelTotAdapM2->Add( fAdapOVCThM2,      fParameters[16]);

  fModelTotAdapM2->Add( fAdapFrameRaM2,    fParameters[2]);
  fModelTotAdapM2->Add( fAdapTShieldRaM2,  fParameters[3]);  
  fModelTotAdapM2->Add( fAdap50mKRaM2,     fParameters[17]);
  fModelTotAdapM2->Add( fAdap600mKRaM2,    fParameters[18]);
  fModelTotAdapM2->Add( fAdapIVCRaM2,      fParameters[19]);
  fModelTotAdapM2->Add( fAdapOVCRaM2,      fParameters[20]);

  fModelTotAdapM2->Add( fAdapFrameKM2,     fParameters[4]);
  fModelTotAdapM2->Add( fAdapTShieldKM2,   fParameters[4]);
  fModelTotAdapM2->Add( fAdap50mKKM2,      fParameters[4]);
  fModelTotAdapM2->Add( fAdap600mKKM2,     fParameters[4]);
  fModelTotAdapM2->Add( fAdapIVCKM2,       fParameters[5]);
  fModelTotAdapM2->Add( fAdapOVCKM2,       fParameters[5]); 

  fModelTotAdapM2->Add( fAdapFrameCoM2,    fParameters[6]);
  fModelTotAdapM2->Add( fAdapTShieldCoM2,  fParameters[7]);
  fModelTotAdapM2->Add( fAdap50mKCoM2,     fParameters[21]);
  fModelTotAdapM2->Add( fAdap600mKCoM2,    fParameters[22]);
  fModelTotAdapM2->Add( fAdapIVCCoM2,      fParameters[23]);
  fModelTotAdapM2->Add( fAdapOVCCoM2,      fParameters[24]);  

  fModelTotAdapM2->Add( fAdapTShieldMnM2,  fParameters[11]);
  fModelTotAdapM2->Add( fAdapIVCMnM2,      fParameters[12]);

  fModelTotAdapM2->Add( fAdapNDBDM2,      fParameters[9]);  
  fModelTotAdapM2->Add( fAdap2NDBDM2,     fParameters[8]);  
  fModelTotAdapM2->Add( fAdapBiM2,        fParameters[10]);  

}


// Saves gamma histograms
void TBackgroundModel::SaveGamma()
{
  double dRes1 = 1.986;
  double dRes2 = 5.332;

  // Bin size (keV)
  int dBinSize = 2;
  // Histogram range
  int dMinEnergy = 0.;
  int dMaxEnergy = 8000.;
  int dNBins = (dMaxEnergy - dMinEnergy)/ dBinSize;


  TH1D *fFrameThM1    = new TH1D("fFrameThM1",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fTShieldThM1  = new TH1D("fTShieldThM1",  "TShield",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f50mKThM1     = new TH1D("f50mKThM1",   "50mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f600mKThM1    = new TH1D("f600mKThM1",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fIVCThM1      = new TH1D("fIVCThM1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fOVCThM1      = new TH1D("fOVCThM1",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fFrameRaM1    = new TH1D("fFrameRaM1",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fTShieldRaM1  = new TH1D("fTShieldRaM1","TShield",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f50mKRaM1     = new TH1D("f50mKRaM1",   "50mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f600mKRaM1    = new TH1D("f600mKRaM1",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fIVCRaM1      = new TH1D("fIVCRaM1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fOVCRaM1      = new TH1D("fOVCRaM1",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fFrameKM1     = new TH1D("fFrameKM1",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fTShieldKM1   = new TH1D("fTShieldKM1", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f50mKKM1      = new TH1D("f50mKKM1",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f600mKKM1     = new TH1D("f600mKKM1",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fIVCKM1       = new TH1D("fIVCKM1",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fOVCKM1       = new TH1D("fOVCKM1",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fFrameCoM1    = new TH1D("fFrameCoM1",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fTShieldCoM1  = new TH1D("fTShieldCoM1","TShield",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f50mKCoM1     = new TH1D("f50mKCoM1",   "50mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f600mKCoM1    = new TH1D("f600mKCoM1",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fIVCCoM1      = new TH1D("fIVCCoM1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fOVCCoM1      = new TH1D("fOVCCoM1",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fTShieldMnM1  = new TH1D("fTShieldMnM1",    "TShield",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fIVCMnM1      = new TH1D("fIVCMnM1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fCrystalBi2M1    = new TH1D("fCrystalBi2M1",    "Crystal",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fFrameBi2M1      = new TH1D("fFrameBi2M1",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fNDBDM1       = new TH1D("fNDBDM1",   "NDBD",       dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f2NDBDM1      = new TH1D("f2NDBDM1",  "2NDBD",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fBiM1         = new TH1D("fBiM1",     "Bi",         dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fCrystalPtM1  = new TH1D("fCrystalPtM1", "Pt",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fCrystalPbBM1 = new TH1D("fCrystalPbBM1", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fCrystalPbS01M1 = new TH1D("fCrystalPbS01M1", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fCrystalPbS1M1 = new TH1D("fCrystalPbS1M1", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fCrystalPbS10M1 = new TH1D("fCrystalPbS10M1", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fCrystalPbS100M1 = new TH1D("fCrystalPbS100M1", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fFramePbBM1 = new TH1D("fFramePbBM1", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fFramePbS01M1 = new TH1D("fFramePbS01M1", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fFramePbS1M1 = new TH1D("fFramePbS1M1", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fFramePbS10M1 = new TH1D("fFramePbS10M1", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fFramePbS100M1 = new TH1D("fFramePbS100M1", "Pb",     dNBins, dMinEnergy, dMaxEnergy);



  TH1D *fFrameThM2    = new TH1D("fFrameThM2",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fTShieldThM2  = new TH1D("fTShieldThM2",  "TShield",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f50mKThM2     = new TH1D("f50mKThM2",   "50mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f600mKThM2    = new TH1D("f600mKThM2",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fIVCThM2      = new TH1D("fIVCThM2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fOVCThM2      = new TH1D("fOVCThM2",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fFrameRaM2    = new TH1D("fFrameRaM2",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fTShieldRaM2  = new TH1D("fTShieldRaM2","TShield",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f50mKRaM2     = new TH1D("f50mKRaM2",   "50mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f600mKRaM2    = new TH1D("f600mKRaM2",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fIVCRaM2      = new TH1D("fIVCRaM2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fOVCRaM2      = new TH1D("fOVCRaM2",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fFrameKM2     = new TH1D("fFrameKM2",   "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fTShieldKM2   = new TH1D("fTShieldKM2", "TShield",  dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f50mKKM2      = new TH1D("f50mKKM2",    "50mK",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f600mKKM2     = new TH1D("f600mKKM2",   "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fIVCKM2       = new TH1D("fIVCKM2",     "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fOVCKM2       = new TH1D("fOVCKM2",     "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fFrameCoM2    = new TH1D("fFrameCoM2",  "Frame",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fTShieldCoM2  = new TH1D("fTShieldCoM2","TShield",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f50mKCoM2     = new TH1D("f50mKCoM2",   "50mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f600mKCoM2    = new TH1D("f600mKCoM2",  "600mK",    dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fIVCCoM2      = new TH1D("fIVCCoM2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fOVCCoM2      = new TH1D("fOVCCoM2",    "OVC",      dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fTShieldMnM2  = new TH1D("fTShieldMnM2",    "TShield",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fIVCMnM2      = new TH1D("fIVCMnM2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fCrystalBi2M2    = new TH1D("fCrystalBi2M2",    "Crystal",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fFrameBi2M2      = new TH1D("fFrameBi2M2",    "IVC",      dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fNDBDM2       = new TH1D("fNDBDM2",   "NDBD",       dNBins, dMinEnergy, dMaxEnergy);
  TH1D *f2NDBDM2      = new TH1D("f2NDBDM2",  "2NDBD",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fBiM2         = new TH1D("fBiM2",     "Bi",         dNBins, dMinEnergy, dMaxEnergy);

  TH1D *fCrystalPtM2  = new TH1D("fCrystalPtM2", "Pt",      dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fCrystalPbBM2 = new TH1D("fCrystalPbBM2", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fCrystalPbS01M2 = new TH1D("fCrystalPbS01M2", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fCrystalPbS1M2 = new TH1D("fCrystalPbS1M2", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fCrystalPbS10M2 = new TH1D("fCrystalPbS10M2", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fCrystalPbS100M2 = new TH1D("fCrystalPbS100M2", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fFramePbBM2 = new TH1D("fFramePbBM2", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fFramePbS01M2 = new TH1D("fFramePbS01M2", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fFramePbS1M2 = new TH1D("fFramePbS1M2", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fFramePbS10M2 = new TH1D("fFramePbS10M2", "Pb",     dNBins, dMinEnergy, dMaxEnergy);
  TH1D *fFramePbS100M2 = new TH1D("fFramePbS100M2", "Pb",     dNBins, dMinEnergy, dMaxEnergy);



  TH1D *fAdapFrameThM1    = new TH1D("fAdapFrameThM1",  "Frame",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapTShieldThM1  = new TH1D("fAdapTShieldThM1", "TShield",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdap50mKThM1     = new TH1D("fAdap50mKThM1",  "50mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdap600mKThM1    = new TH1D("fAdap600mKThM1",  "600mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapIVCThM1      = new TH1D("fAdapIVCThM1",    "IVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapOVCThM1      = new TH1D("fAdapOVCThM1",    "OVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);

  TH1D *fAdapFrameRaM1    = new TH1D("fAdapFrameRaM1",  "Frame",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapTShieldRaM1  = new TH1D("fAdapTShieldRaM1", "TShield",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdap50mKRaM1     = new TH1D("fAdap50mKRaM1",  "50mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdap600mKRaM1    = new TH1D("fAdap600mKRaM1",  "600mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapIVCRaM1      = new TH1D("fAdapIVCRaM1",    "IVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapOVCRaM1      = new TH1D("fAdapOVCRaM1",    "OVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);

  TH1D *fAdapFrameKM1     = new TH1D("fAdapFrameKM1",   "Frame",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapTShieldKM1   = new TH1D("fAdapTShieldKM1", "TShield",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdap50mKKM1      = new TH1D("fAdap50mKKM1",    "50mK",     dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdap600mKKM1     = new TH1D("fAdap600mKKM1",   "600mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapIVCKM1       = new TH1D("fAdapIVCKM1",     "IVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapOVCKM1       = new TH1D("fAdapOVCKM1",     "OVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);

  TH1D *fAdapFrameCoM1    = new TH1D("fAdapFrameCoM1",  "Frame",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapTShieldCoM1  = new TH1D("fAdapTShieldCoM1", "TShield",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdap50mKCoM1     = new TH1D("fAdap50mKCoM1",  "50mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdap600mKCoM1    = new TH1D("fAdap600mKCoM1",  "600mK",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapIVCCoM1      = new TH1D("fAdapIVCCoM1",    "IVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapOVCCoM1      = new TH1D("fAdapOVCCoM1",    "OVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);

  TH1D *fAdapTShieldMnM1  = new TH1D("fAdapTShieldMnM1", "TShield",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapIVCMnM1      = new TH1D("fAdapIVCMnM1",    "IVC",      dAdaptiveBinsM1, dAdaptiveArrayM1);

  TH1D *fAdapCrystalBi2M1  = new TH1D("fAdapCrystalBi2M1", "Crystal",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapFrameBi2M1  = new TH1D("fAdapFrameBi2M1", "Frame",      dAdaptiveBinsM1, dAdaptiveArrayM1);

  TH1D *fAdapNDBDM1       = new TH1D("fAdapNDBDM1",   "NDBD",       dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdap2NDBDM1      = new TH1D("fAdap2NDBDM1",  "2NDBD",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapBiM1         = new TH1D("fAdapBiM1",     "Bi",         dAdaptiveBinsM1, dAdaptiveArrayM1);

  TH1D *fAdapCrystalPtM1     = new TH1D("fAdapCrystalPtM1", "Pt",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapCrystalPbBM1    = new TH1D("fAdapCrystalPbBM1", "Pb",     dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapCrystalPbS01M1  = new TH1D("fAdapCrystalPbS01M1", "Pb",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapCrystalPbS1M1   = new TH1D("fAdapCrystalPbS1M1", "Pb",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapCrystalPbS10M1  = new TH1D("fAdapCrystalPbS10M1", "Pb",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapCrystalPbS100M1 = new TH1D("fAdapCrystalPbS100M1", "Pb",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapFramePbBM1      = new TH1D("fAdapFramePbBM1", "Pb",       dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapFramePbS01M1    = new TH1D("fAdapFramePbS01M1", "Pb",     dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapFramePbS1M1     = new TH1D("fAdapFramePbS1M1", "Pb",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapFramePbS10M1    = new TH1D("fAdapFramePbS10M1", "Pb",     dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *fAdapFramePbS100M1   = new TH1D("fAdapFramePbS100M1", "Pb",    dAdaptiveBinsM1, dAdaptiveArrayM1);



  TH1D *fAdapFrameThM2    = new TH1D("fAdapFrameThM2",  "Frame",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapTShieldThM2  = new TH1D("fAdapTShieldThM2", "TShield",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdap50mKThM2     = new TH1D("fAdap50mKThM2",  "50mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdap600mKThM2    = new TH1D("fAdap600mKThM2",  "600mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapIVCThM2      = new TH1D("fAdapIVCThM2",    "IVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapOVCThM2      = new TH1D("fAdapOVCThM2",    "OVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);

  TH1D *fAdapFrameRaM2    = new TH1D("fAdapFrameRaM2",  "Frame",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapTShieldRaM2  = new TH1D("fAdapTShieldRaM2", "TShield",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdap50mKRaM2     = new TH1D("fAdap50mKRaM2",  "50mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdap600mKRaM2    = new TH1D("fAdap600mKRaM2",  "600mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapIVCRaM2      = new TH1D("fAdapIVCRaM2",    "IVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapOVCRaM2      = new TH1D("fAdapOVCRaM2",    "OVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);

  TH1D *fAdapFrameKM2     = new TH1D("fAdapFrameKM2",   "Frame",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapTShieldKM2   = new TH1D("fAdapTShieldKM2", "TShield",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdap50mKKM2      = new TH1D("fAdap50mKKM2",    "50mK",     dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdap600mKKM2     = new TH1D("fAdap600mKKM2",   "600mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapIVCKM2       = new TH1D("fAdapIVCKM2",     "IVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapOVCKM2       = new TH1D("fAdapOVCKM2",     "OVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);

  TH1D *fAdapFrameCoM2    = new TH1D("fAdapFrameCoM2",  "Frame",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapTShieldCoM2  = new TH1D("fAdapTShieldCoM2", "TShield",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdap50mKCoM2     = new TH1D("fAdap50mKCoM2",  "50mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdap600mKCoM2    = new TH1D("fAdap600mKCoM2",  "600mK",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapIVCCoM2      = new TH1D("fAdapIVCCoM2",    "IVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapOVCCoM2      = new TH1D("fAdapOVCCoM2",    "OVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);

  TH1D *fAdapTShieldMnM2  = new TH1D("fAdapTShieldMnM2", "TShield",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapIVCMnM2      = new TH1D("fAdapIVCMnM2",    "IVC",      dAdaptiveBinsM2, dAdaptiveArrayM2);

  TH1D *fAdapCrystalBi2M2  = new TH1D("fAdapCrystalBi2M2", "Crystal",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapFrameBi2M2  = new TH1D("fAdapFrameBi2M2", "Frame",      dAdaptiveBinsM2, dAdaptiveArrayM2);


  TH1D *fAdapNDBDM2       = new TH1D("fAdapNDBDM2",   "NDBD",       dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdap2NDBDM2      = new TH1D("fAdap2NDBDM2",  "2NDBD",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapBiM2         = new TH1D("fAdapBiM2",     "Bi",         dAdaptiveBinsM2, dAdaptiveArrayM2);

  TH1D *fAdapCrystalPtM2     = new TH1D("fAdapCrystalPtM2", "Pt",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapCrystalPbBM2    = new TH1D("fAdapCrystalPbBM2", "Pb",     dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapCrystalPbS01M2  = new TH1D("fAdapCrystalPbS01M2", "Pb",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapCrystalPbS1M2   = new TH1D("fAdapCrystalPbS1M2", "Pb",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapCrystalPbS10M2  = new TH1D("fAdapCrystalPbS10M2", "Pb",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapCrystalPbS100M2 = new TH1D("fAdapCrystalPbS100M2", "Pb",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapFramePbBM2      = new TH1D("fAdapFramePbBM2", "Pb",       dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapFramePbS01M2    = new TH1D("fAdapFramePbS01M2", "Pb",     dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapFramePbS1M2     = new TH1D("fAdapFramePbS1M2", "Pb",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapFramePbS10M2    = new TH1D("fAdapFramePbS10M2", "Pb",     dAdaptiveBinsM2, dAdaptiveArrayM2);
  TH1D *fAdapFramePbS100M2   = new TH1D("fAdapFramePbS100M2", "Pb",    dAdaptiveBinsM2, dAdaptiveArrayM2);



  TFile *fFrameTh1 = new TFile(Form("%s/Frame-Th232-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *fTShieldTh1 = new TFile(Form("%s/TShield-Th232-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *f50mKTh1 = new TFile(Form("%s/50mK-Th232-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *f600mKTh1 = new TFile(Form("%s/600mK-Th232-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *fIVCTh1 = new TFile(Form("%s/IVC-Th232-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *fOVCTh1 = new TFile(Form("%s/OVC-Th232-B-M1-T50.root", dDataDir.c_str()) );  

  TFile *fFrameRa1 = new TFile(Form("%s/Frame-Ra226-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *fTShieldRa1 = new TFile(Form("%s/TShield-Ra226-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *f50mKRa1 = new TFile(Form("%s/50mK-Ra226-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *f600mKRa1 = new TFile(Form("%s/600mK-Ra226-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *fIVCRa1 = new TFile(Form("%s/IVC-Ra226-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *fOVCRa1 = new TFile(Form("%s/OVC-Ra226-B-M1-T50.root", dDataDir.c_str()) );  

  TFile *fCrystalNDBD1 = new TFile(Form("%s/Crystal-0NDBD-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *fCrystal2NDBD1 = new TFile(Form("%s/Crystal-2NDBD-B-M1-T50.root", dDataDir.c_str()) ); 

  TFile *fFrameCo1 = new TFile(Form("%s/Frame-Co60-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *fTShieldCo1 = new TFile(Form("%s/TShield-Co60-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *f50mKCo1 = new TFile(Form("%s/50mK-Co60-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *f600mKCo1 = new TFile(Form("%s/600mK-Co60-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *fIVCCo1 = new TFile(Form("%s/IVC-Co60-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *fOVCCo1 = new TFile(Form("%s/OVC-Co60-B-M1-T50.root", dDataDir.c_str()) );  

  TFile *fRLeadBi1 = new TFile(Form("%s/RLead-Bi207-B-M1-T50.root", dDataDir.c_str()) );

  TFile *fTShieldMn1 = new TFile(Form("%s/TShield-Mn54-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *fIVCMn1 = new TFile(Form("%s/IVC-Mn54-B-M1-T50.root", dDataDir.c_str()) );

  TFile *fCrystalBi21 = new TFile(Form("%s/Crystal-Bi210-B-M1-T50.root", dDataDir.c_str()) );  
  TFile *fFrameBi21 = new TFile(Form("%s/Frame-Bi210-B-M1-T50.root", dDataDir.c_str()) );  


  TFile *fFrameK1 = new TFile(Form("%s/Frame-K40-B-M1-T50.root", dDataDir.c_str()) );
  TFile *fTShieldK1 = new TFile(Form("%s/TShield-K40-B-M1-T50.root", dDataDir.c_str()) );
  TFile *f50mKK1 = new TFile(Form("%s/50mK-K40-B-M1-T50.root", dDataDir.c_str()) );
  TFile *f600mKK1 = new TFile(Form("%s/600mK-K40-B-M1-T50.root", dDataDir.c_str()) );
  TFile *fIVCK1 = new TFile(Form("%s/IVC-K40-B-M1-T50.root", dDataDir.c_str()) );
  TFile *fOVCK1 = new TFile(Form("%s/OVC-K40-B-M1-T50.root", dDataDir.c_str()) );



  TFile *fFrameTh2 = new TFile(Form("%s/Frame-Th232-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *fTShieldTh2 = new TFile(Form("%s/TShield-Th232-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *f50mKTh2 = new TFile(Form("%s/50mK-Th232-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *f600mKTh2 = new TFile(Form("%s/600mK-Th232-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *fIVCTh2 = new TFile(Form("%s/IVC-Th232-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *fOVCTh2 = new TFile(Form("%s/OVC-Th232-B-M2-T50.root", dDataDir.c_str()) );  

  TFile *fFrameRa2 = new TFile(Form("%s/Frame-Ra226-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *fTShieldRa2 = new TFile(Form("%s/TShield-Ra226-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *f50mKRa2 = new TFile(Form("%s/50mK-Ra226-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *f600mKRa2 = new TFile(Form("%s/600mK-Ra226-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *fIVCRa2 = new TFile(Form("%s/IVC-Ra226-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *fOVCRa2 = new TFile(Form("%s/OVC-Ra226-B-M2-T50.root", dDataDir.c_str()) );  

  TFile *fCrystalNDBD2 = new TFile(Form("%s/Crystal-0NDBD-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *fCrystal2NDBD2 = new TFile(Form("%s/Crystal-2NDBD-B-M2-T50.root", dDataDir.c_str()) ); 

  TFile *fFrameCo2 = new TFile(Form("%s/Frame-Co60-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *fTShieldCo2 = new TFile(Form("%s/TShield-Co60-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *f50mKCo2 = new TFile(Form("%s/50mK-Co60-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *f600mKCo2 = new TFile(Form("%s/600mK-Co60-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *fIVCCo2 = new TFile(Form("%s/IVC-Co60-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *fOVCCo2 = new TFile(Form("%s/OVC-Co60-B-M2-T50.root", dDataDir.c_str()) );  

  TFile *fRLeadBi2 = new TFile(Form("%s/RLead-Bi207-B-M2-T50.root", dDataDir.c_str()) );
  TFile *fTShieldMn2 = new TFile(Form("%s/TShield-Mn54-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *fIVCMn2 = new TFile(Form("%s/IVC-Mn54-B-M2-T50.root", dDataDir.c_str()) );
  
  TFile *fCrystalBi22 = new TFile(Form("%s/Crystal-Bi210-B-M2-T50.root", dDataDir.c_str()) );  
  TFile *fFrameBi22 = new TFile(Form("%s/Frame-Bi210-B-M2-T50.root", dDataDir.c_str()) );  

  TFile *fFrameK2 = new TFile(Form("%s/Frame-K40-B-M2-T50.root", dDataDir.c_str()) );
  TFile *fTShieldK2 = new TFile(Form("%s/TShield-K40-B-M2-T50.root", dDataDir.c_str()) );
  TFile *f50mKK2 = new TFile(Form("%s/50mK-K40-B-M2-T50.root", dDataDir.c_str()) );
  TFile *f600mKK2 = new TFile(Form("%s/600mK-K40-B-M2-T50.root", dDataDir.c_str()) );
  TFile *fIVCK2 = new TFile(Form("%s/IVC-K40-B-M2-T50.root", dDataDir.c_str()) );
  TFile *fOVCK2 = new TFile(Form("%s/OVC-K40-B-M2-T50.root", dDataDir.c_str()) );


  cout << "Loaded Files" << endl;


  // Have to rebin because bin size is 2 but default is 1 keV bins in the output...
  fFrameThM1 = (TH1D*)fFrameTh1->Get("h0"); fFrameThM1->Rebin();
  fTShieldThM1 = (TH1D*)fTShieldTh1->Get("h0"); fTShieldThM1->Rebin();
  f50mKThM1 = (TH1D*)f50mKTh1->Get("h0"); f50mKThM1->Rebin();
  f600mKThM1 = (TH1D*)f600mKTh1->Get("h0"); f600mKThM1->Rebin();
  fIVCThM1 = (TH1D*)fIVCTh1->Get("h0");  fIVCThM1->Rebin();
  fOVCThM1 = (TH1D*)fOVCTh1->Get("h0");  fOVCThM1->Rebin();

  fFrameRaM1 = (TH1D*)fFrameRa1->Get("h0");  fFrameRaM1->Rebin();
  fTShieldRaM1 = (TH1D*)fTShieldRa1->Get("h0");  fTShieldRaM1->Rebin();
  f50mKRaM1 = (TH1D*)f50mKRa1->Get("h0");  f50mKRaM1->Rebin();
  f600mKRaM1 = (TH1D*)f600mKRa1->Get("h0");  f600mKRaM1->Rebin();
  fIVCRaM1 = (TH1D*)fIVCRa1->Get("h0");  fIVCRaM1->Rebin();
  fOVCRaM1 = (TH1D*)fOVCRa1->Get("h0");  fOVCRaM1->Rebin();

  fFrameKM1 = (TH1D*)fFrameK1->Get("h0"); fFrameKM1->Rebin();
  fTShieldKM1 = (TH1D*)fTShieldK1->Get("h0"); fTShieldKM1->Rebin();
  f50mKKM1 = (TH1D*)f50mKK1->Get("h0"); f50mKKM1->Rebin();
  f600mKKM1 = (TH1D*)f600mKK1->Get("h0"); f600mKKM1->Rebin();
  fIVCKM1 = (TH1D*)fIVCK1->Get("h0");  fIVCKM1->Rebin();
  fOVCKM1 = (TH1D*)fOVCK1->Get("h0");  fOVCKM1->Rebin();

  fFrameCoM1 = (TH1D*)fFrameCo1->Get("h0"); fFrameCoM1->Rebin();
  fTShieldCoM1 = (TH1D*)fTShieldCo1->Get("h0"); fTShieldCoM1->Rebin();
  f50mKCoM1 = (TH1D*)f50mKCo1->Get("h0"); f50mKCoM1->Rebin();
  f600mKCoM1 = (TH1D*)f600mKCo1->Get("h0"); f600mKCoM1->Rebin();
  fIVCCoM1 = (TH1D*)fIVCCo1->Get("h0"); fIVCCoM1->Rebin();
  fOVCCoM1 = (TH1D*)fOVCCo1->Get("h0"); fOVCCoM1->Rebin();

  fTShieldMnM1 = (TH1D*)fTShieldMn1->Get("h0"); fTShieldMnM1->Rebin();
  fIVCMnM1 = (TH1D*)fIVCMn1->Get("h0"); fIVCMnM1->Rebin();

  fCrystalBi2M1 = (TH1D*)fCrystalBi21->Get("h0"); fCrystalBi2M1->Rebin();
  fFrameBi2M1 = (TH1D*)fFrameBi21->Get("h0"); fFrameBi2M1->Rebin();

  fNDBDM1 = (TH1D*)fCrystalNDBD1->Get("h0"); fNDBDM1->Rebin();
  f2NDBDM1 = (TH1D*)fCrystal2NDBD1->Get("h0"); f2NDBDM1->Rebin();
  fBiM1 = (TH1D*)fRLeadBi1->Get("h0"); fBiM1->Rebin();



  fFrameThM2 = (TH1D*)fFrameTh2->Get("h0"); fFrameThM2->Rebin();
  fTShieldThM2 = (TH1D*)fTShieldTh2->Get("h0"); fTShieldThM2->Rebin();
  f50mKThM2 = (TH1D*)f50mKTh2->Get("h0"); f50mKThM2->Rebin();
  f600mKThM2 = (TH1D*)f600mKTh2->Get("h0"); f600mKThM2->Rebin();
  fIVCThM2 = (TH1D*)fIVCTh2->Get("h0");  fIVCThM2->Rebin();
  fOVCThM2 = (TH1D*)fOVCTh2->Get("h0");  fOVCThM2->Rebin();

  fFrameRaM2 = (TH1D*)fFrameRa2->Get("h0");  fFrameRaM2->Rebin();
  fTShieldRaM2 = (TH1D*)fTShieldRa2->Get("h0");  fTShieldRaM2->Rebin();
  f50mKRaM2 = (TH1D*)f50mKRa2->Get("h0");  f50mKRaM2->Rebin();
  f600mKRaM2 = (TH1D*)f600mKRa2->Get("h0");  f600mKRaM2->Rebin();
  fIVCRaM2 = (TH1D*)fIVCRa2->Get("h0");  fIVCRaM2->Rebin();
  fOVCRaM2 = (TH1D*)fOVCRa2->Get("h0");  fOVCRaM2->Rebin();

  fFrameKM2 = (TH1D*)fFrameK2->Get("h0"); fFrameKM2->Rebin();
  fTShieldKM2 = (TH1D*)fTShieldK2->Get("h0"); fTShieldKM2->Rebin();
  f50mKKM2 = (TH1D*)f50mKK2->Get("h0"); f50mKKM2->Rebin();
  f600mKKM2 = (TH1D*)f600mKK2->Get("h0"); f600mKKM2->Rebin();
  fIVCKM2 = (TH1D*)fIVCK2->Get("h0");  fIVCKM2->Rebin();
  fOVCKM2 = (TH1D*)fOVCK2->Get("h0");  fOVCKM2->Rebin();

  fFrameCoM2 = (TH1D*)fFrameCo2->Get("h0"); fFrameCoM2->Rebin();
  fTShieldCoM2 = (TH1D*)fTShieldCo2->Get("h0"); fTShieldCoM2->Rebin();
  f50mKCoM2 = (TH1D*)f50mKCo2->Get("h0"); f50mKCoM2->Rebin();
  f600mKCoM2 = (TH1D*)f600mKCo2->Get("h0"); f600mKCoM2->Rebin();
  fIVCCoM2 = (TH1D*)fIVCCo2->Get("h0"); fIVCCoM2->Rebin();
  fOVCCoM2 = (TH1D*)fOVCCo2->Get("h0"); fOVCCoM2->Rebin();

  fTShieldMnM2 = (TH1D*)fTShieldMn2->Get("h0"); fTShieldMnM2->Rebin();
  fIVCMnM2 = (TH1D*)fIVCMn2->Get("h0"); fIVCMnM2->Rebin();

  fCrystalBi2M2 = (TH1D*)fCrystalBi22->Get("h0"); fCrystalBi2M2->Rebin();
  fFrameBi2M2 = (TH1D*)fFrameBi22->Get("h0"); fFrameBi2M2->Rebin();

  fNDBDM2 = (TH1D*)fCrystalNDBD2->Get("h0"); fNDBDM2->Rebin();
  f2NDBDM2 = (TH1D*)fCrystal2NDBD2->Get("h0"); f2NDBDM2->Rebin();
  fBiM2 = (TH1D*)fRLeadBi2->Get("h0"); fBiM2->Rebin();



  NormalizePDFPair(fFrameThM1, fFrameThM2, 50, 2700);
  NormalizePDFPair(fTShieldThM1, fTShieldThM2, 50, 2700);
  NormalizePDFPair(f50mKThM1, f50mKThM2, 50, 2700);
  NormalizePDFPair(f600mKThM1, f600mKThM2, 50, 2700);
  NormalizePDFPair(fIVCThM1, fIVCThM2, 50, 2700);
  NormalizePDFPair(fOVCThM1, fOVCThM2, 50, 2700);

  NormalizePDFPair(fFrameRaM1, fFrameRaM2, 50, 2700);
  NormalizePDFPair(fTShieldRaM1, fTShieldRaM2, 50, 2700);
  NormalizePDFPair(f50mKRaM1, f50mKRaM2, 50, 2700);
  NormalizePDFPair(f600mKRaM1, f600mKRaM2, 50, 2700);
  NormalizePDFPair(fIVCRaM1, fIVCRaM2, 50, 2700);
  NormalizePDFPair(fOVCRaM1, fOVCRaM2, 50, 2700);

  NormalizePDFPair(fFrameKM1, fFrameKM2, 50, 2700);
  NormalizePDFPair(fTShieldKM1, fTShieldKM2, 50, 2700);
  NormalizePDFPair(f50mKKM1, f50mKKM2, 50, 2700);
  NormalizePDFPair(f600mKKM1, f600mKKM2, 50, 2700);
  NormalizePDFPair(fIVCKM1, fIVCKM2, 50, 2700);
  NormalizePDFPair(fOVCKM1, fOVCKM2, 50, 2700);

  NormalizePDFPair(fFrameCoM1, fFrameCoM2, 50, 2700);
  NormalizePDFPair(fTShieldCoM1, fTShieldCoM2, 50, 2700);
  NormalizePDFPair(f50mKCoM1, f50mKCoM2, 50, 2700);
  NormalizePDFPair(f600mKCoM1, f600mKCoM2, 50, 2700);
  NormalizePDFPair(fIVCCoM1, fIVCCoM2, 50, 2700);
  NormalizePDFPair(fOVCCoM1, fOVCCoM2, 50, 2700);

  NormalizePDFPair(fTShieldMnM1, fTShieldMnM2, 50, 2700);
  NormalizePDFPair(fIVCMnM1, fIVCMnM2, 50, 2700);

  NormalizePDFPair(fCrystalBi2M1, fCrystalBi2M2, 50, 2700);
  NormalizePDFPair(fFrameBi2M1, fFrameBi2M2, 50, 2700);

  NormalizePDFPair(fNDBDM1, fNDBDM2, 50, 2700);
  NormalizePDFPair(f2NDBDM1, f2NDBDM2, 50, 2700);
  NormalizePDFPair(fBiM1, fBiM2, 50, 2700);

  cout << "Normalized Histograms" << endl;


  SmearMC(fFrameThM1, fSmearFrameThM1, dRes1, dRes2);
  SmearMC(fTShieldThM1, fSmearTShieldThM1, dRes1, dRes2);
  SmearMC(f50mKThM1, fSmear50mKThM1, dRes1, dRes2);
  SmearMC(f600mKThM1, fSmear600mKThM1, dRes1, dRes2);
  SmearMC(fIVCThM1, fSmearIVCThM1, dRes1, dRes2);
  SmearMC(fOVCThM1, fSmearOVCThM1, dRes1, dRes2);

  SmearMC(fFrameRaM1, fSmearFrameRaM1, dRes1, dRes2);
  SmearMC(fTShieldRaM1, fSmearTShieldRaM1, dRes1, dRes2);
  SmearMC(f50mKRaM1, fSmear50mKRaM1, dRes1, dRes2);
  SmearMC(f600mKRaM1, fSmear600mKRaM1, dRes1, dRes2);
  SmearMC(fIVCRaM1, fSmearIVCRaM1, dRes1, dRes2);
  SmearMC(fOVCRaM1, fSmearOVCRaM1, dRes1, dRes2);

  SmearMC(fFrameKM1, fSmearFrameKM1, dRes1, dRes2);
  SmearMC(fTShieldKM1, fSmearTShieldKM1, dRes1, dRes2);
  SmearMC(f50mKKM1, fSmear50mKKM1, dRes1, dRes2);
  SmearMC(f600mKKM1, fSmear600mKKM1, dRes1, dRes2);
  SmearMC(fIVCKM1, fSmearIVCKM1, dRes1, dRes2);
  SmearMC(fOVCKM1, fSmearOVCKM1, dRes1, dRes2);

  SmearMC(fFrameCoM1, fSmearFrameCoM1, dRes1, dRes2);
  SmearMC(fTShieldCoM1, fSmearTShieldCoM1, dRes1, dRes2);
  SmearMC(f50mKCoM1, fSmear50mKCoM1, dRes1, dRes2);
  SmearMC(f600mKCoM1, fSmear600mKCoM1, dRes1, dRes2);
  SmearMC(fIVCCoM1, fSmearIVCCoM1, dRes1, dRes2);
  SmearMC(fOVCCoM1, fSmearOVCCoM1, dRes1, dRes2);

  SmearMC(fTShieldMnM1, fSmearTShieldMnM1, dRes1, dRes2);
  SmearMC(fIVCMnM1, fSmearIVCMnM1, dRes1, dRes2);

  SmearMC(fCrystalBi2M1, fSmearCrystalBi2M1, dRes1, dRes2);
  SmearMC(fFrameBi2M1, fSmearFrameBi2M1, dRes1, dRes2);

  SmearMC(fNDBDM1, fSmearNDBDM1, dRes1, dRes2);
  SmearMC(f2NDBDM1, fSmear2NDBDM1, dRes1, dRes2);
  SmearMC(fBiM1, fSmearBiM1, dRes1, dRes2);


  SmearMC(fFrameThM2, fSmearFrameThM2, dRes1, dRes2);
  SmearMC(fTShieldThM2, fSmearTShieldThM2, dRes1, dRes2);
  SmearMC(f50mKThM2, fSmear50mKThM2, dRes1, dRes2);
  SmearMC(f600mKThM2, fSmear600mKThM2, dRes1, dRes2);
  SmearMC(fIVCThM2, fSmearIVCThM2, dRes1, dRes2);
  SmearMC(fOVCThM2, fSmearOVCThM2, dRes1, dRes2);

  SmearMC(fFrameRaM2, fSmearFrameRaM2, dRes1, dRes2);
  SmearMC(fTShieldRaM2, fSmearTShieldRaM2, dRes1, dRes2);
  SmearMC(f50mKRaM2, fSmear50mKRaM2, dRes1, dRes2);
  SmearMC(f600mKRaM2, fSmear600mKRaM2, dRes1, dRes2);
  SmearMC(fIVCRaM2, fSmearIVCRaM2, dRes1, dRes2);
  SmearMC(fOVCRaM2, fSmearOVCRaM2, dRes1, dRes2);

  SmearMC(fFrameKM2, fSmearFrameKM2, dRes1, dRes2);
  SmearMC(fTShieldKM2, fSmearTShieldKM2, dRes1, dRes2);
  SmearMC(f50mKKM2, fSmear50mKKM2, dRes1, dRes2);
  SmearMC(f600mKKM2, fSmear600mKKM2, dRes1, dRes2);
  SmearMC(fIVCKM2, fSmearIVCKM2, dRes1, dRes2);
  SmearMC(fOVCKM2, fSmearOVCKM2, dRes1, dRes2);

  SmearMC(fFrameCoM2, fSmearFrameCoM2, dRes1, dRes2);
  SmearMC(fTShieldCoM2, fSmearTShieldCoM2, dRes1, dRes2);
  SmearMC(f50mKCoM2, fSmear50mKCoM2, dRes1, dRes2);
  SmearMC(f600mKCoM2, fSmear600mKCoM2, dRes1, dRes2);
  SmearMC(fIVCCoM2, fSmearIVCCoM2, dRes1, dRes2);
  SmearMC(fOVCCoM2, fSmearOVCCoM2, dRes1, dRes2);

  SmearMC(fTShieldMnM2, fSmearTShieldMnM2, dRes1, dRes2);
  SmearMC(fIVCMnM2, fSmearIVCMnM2, dRes1, dRes2);

  SmearMC(fCrystalBi2M2, fSmearCrystalBi2M2, dRes1, dRes2);
  SmearMC(fFrameBi2M2, fSmearFrameBi2M2, dRes1, dRes2);

  SmearMC(fNDBDM2, fSmearNDBDM2, dRes1, dRes2);
  SmearMC(f2NDBDM2, fSmear2NDBDM2, dRes1, dRes2);
  SmearMC(fBiM2, fSmearBiM2, dRes1, dRes2);



  /////////////// Rebinning for Adaptive


  fSmearFrameThM1->Rebin(dAdaptiveBinsM1, "fnewFrameThM1", dAdaptiveArrayM1);
  fSmearTShieldThM1->Rebin(dAdaptiveBinsM1, "fnewTShieldThM1", dAdaptiveArrayM1);
  fSmear50mKThM1->Rebin(dAdaptiveBinsM1, "fnew50mKThM1", dAdaptiveArrayM1);
  fSmear600mKThM1->Rebin(dAdaptiveBinsM1, "fnew600mKThM1", dAdaptiveArrayM1);
  fSmearIVCThM1->Rebin(dAdaptiveBinsM1, "fnewIVCThM1", dAdaptiveArrayM1);
  fSmearOVCThM1->Rebin(dAdaptiveBinsM1, "fnewOVCThM1", dAdaptiveArrayM1);

  fSmearFrameRaM1->Rebin(dAdaptiveBinsM1, "fnewFrameRaM1", dAdaptiveArrayM1);
  fSmearTShieldRaM1->Rebin(dAdaptiveBinsM1, "fnewTShieldRaM1", dAdaptiveArrayM1);
  fSmear50mKRaM1->Rebin(dAdaptiveBinsM1, "fnew50mKRaM1", dAdaptiveArrayM1);
  fSmear600mKRaM1->Rebin(dAdaptiveBinsM1, "fnew600mKRaM1", dAdaptiveArrayM1);
  fSmearIVCRaM1->Rebin(dAdaptiveBinsM1, "fnewIVCRaM1", dAdaptiveArrayM1);
  fSmearOVCRaM1->Rebin(dAdaptiveBinsM1, "fnewOVCRaM1", dAdaptiveArrayM1);

  fSmearFrameKM1->Rebin(dAdaptiveBinsM1, "fnewFrameKM1", dAdaptiveArrayM1);
  fSmearTShieldKM1->Rebin(dAdaptiveBinsM1, "fnewTShieldKM1", dAdaptiveArrayM1);
  fSmear50mKKM1->Rebin(dAdaptiveBinsM1, "fnew50mKKM1", dAdaptiveArrayM1);
  fSmear600mKKM1->Rebin(dAdaptiveBinsM1, "fnew600mKKM1", dAdaptiveArrayM1);
  fSmearIVCKM1->Rebin(dAdaptiveBinsM1, "fnewIVCKM1", dAdaptiveArrayM1);
  fSmearOVCKM1->Rebin(dAdaptiveBinsM1, "fnewOVCKM1", dAdaptiveArrayM1);

  fSmearFrameCoM1->Rebin(dAdaptiveBinsM1, "fnewFrameCoM1", dAdaptiveArrayM1);
  fSmearTShieldCoM1->Rebin(dAdaptiveBinsM1, "fnewTShieldCoM1", dAdaptiveArrayM1);
  fSmear50mKCoM1->Rebin(dAdaptiveBinsM1, "fnew50mKCoM1", dAdaptiveArrayM1);
  fSmear600mKCoM1->Rebin(dAdaptiveBinsM1, "fnew600mKCoM1", dAdaptiveArrayM1);
  fSmearIVCCoM1->Rebin(dAdaptiveBinsM1, "fnewIVCCoM1", dAdaptiveArrayM1);
  fSmearOVCCoM1->Rebin(dAdaptiveBinsM1, "fnewOVCCoM1", dAdaptiveArrayM1);

  fSmearTShieldMnM1->Rebin(dAdaptiveBinsM1, "fnewTShieldMnM1", dAdaptiveArrayM1);
  fSmearIVCMnM1->Rebin(dAdaptiveBinsM1, "fnewIVCMnM1", dAdaptiveArrayM1);

  fSmearCrystalBi2M1->Rebin(dAdaptiveBinsM1, "fnewCrystalBi2M1", dAdaptiveArrayM1);
  fSmearFrameBi2M1->Rebin(dAdaptiveBinsM1, "fnewFrameBi2M1", dAdaptiveArrayM1);

  fSmearNDBDM1->Rebin(dAdaptiveBinsM1, "fnewNDBDM1", dAdaptiveArrayM1);
  fSmear2NDBDM1->Rebin(dAdaptiveBinsM1, "fnew2NDBDM1", dAdaptiveArrayM1);
  fSmearBiM1->Rebin(dAdaptiveBinsM1, "fnewBiM1", dAdaptiveArrayM1);



  fSmearFrameThM2->Rebin(dAdaptiveBinsM2, "fnewFrameThM2", dAdaptiveArrayM2);
  fSmearTShieldThM2->Rebin(dAdaptiveBinsM2, "fnewTShieldThM2", dAdaptiveArrayM2);
  fSmear50mKThM2->Rebin(dAdaptiveBinsM2, "fnew50mKThM2", dAdaptiveArrayM2);
  fSmear600mKThM2->Rebin(dAdaptiveBinsM2, "fnew600mKThM2", dAdaptiveArrayM2);
  fSmearIVCThM2->Rebin(dAdaptiveBinsM2, "fnewIVCThM2", dAdaptiveArrayM2);
  fSmearOVCThM2->Rebin(dAdaptiveBinsM2, "fnewOVCThM2", dAdaptiveArrayM2);

  fSmearFrameRaM2->Rebin(dAdaptiveBinsM2, "fnewFrameRaM2", dAdaptiveArrayM2);
  fSmearTShieldRaM2->Rebin(dAdaptiveBinsM2, "fnewTShieldRaM2", dAdaptiveArrayM2);
  fSmear50mKRaM2->Rebin(dAdaptiveBinsM2, "fnew50mKRaM2", dAdaptiveArrayM2);
  fSmear600mKRaM2->Rebin(dAdaptiveBinsM2, "fnew600mKRaM2", dAdaptiveArrayM2);
  fSmearIVCRaM2->Rebin(dAdaptiveBinsM2, "fnewIVCRaM2", dAdaptiveArrayM2);
  fSmearOVCRaM2->Rebin(dAdaptiveBinsM2, "fnewOVCRaM2", dAdaptiveArrayM2);

  fSmearFrameKM2->Rebin(dAdaptiveBinsM2, "fnewFrameKM2", dAdaptiveArrayM2);
  fSmearTShieldKM2->Rebin(dAdaptiveBinsM2, "fnewTShieldKM2", dAdaptiveArrayM2);
  fSmear50mKKM2->Rebin(dAdaptiveBinsM2, "fnew50mKKM2", dAdaptiveArrayM2);
  fSmear600mKKM2->Rebin(dAdaptiveBinsM2, "fnew600mKKM2", dAdaptiveArrayM2);
  fSmearIVCKM2->Rebin(dAdaptiveBinsM2, "fnewIVCKM2", dAdaptiveArrayM2);
  fSmearOVCKM2->Rebin(dAdaptiveBinsM2, "fnewOVCKM2", dAdaptiveArrayM2);

  fSmearFrameCoM2->Rebin(dAdaptiveBinsM2, "fnewFrameCoM2", dAdaptiveArrayM2);
  fSmearTShieldCoM2->Rebin(dAdaptiveBinsM2, "fnewTShieldCoM2", dAdaptiveArrayM2);
  fSmear50mKCoM2->Rebin(dAdaptiveBinsM2, "fnew50mKCoM2", dAdaptiveArrayM2);
  fSmear600mKCoM2->Rebin(dAdaptiveBinsM2, "fnew600mKCoM2", dAdaptiveArrayM2);
  fSmearIVCCoM2->Rebin(dAdaptiveBinsM2, "fnewIVCCoM2", dAdaptiveArrayM2);
  fSmearOVCCoM2->Rebin(dAdaptiveBinsM2, "fnewOVCCoM2", dAdaptiveArrayM2);

  fSmearTShieldMnM2->Rebin(dAdaptiveBinsM2, "fnewTShieldMnM2", dAdaptiveArrayM2);
  fSmearIVCMnM2->Rebin(dAdaptiveBinsM2, "fnewIVCMnM2", dAdaptiveArrayM2);

  fSmearCrystalBi2M2->Rebin(dAdaptiveBinsM2, "fnewCrystalBi2M2", dAdaptiveArrayM2);
  fSmearFrameBi2M2->Rebin(dAdaptiveBinsM2, "fnewFrameBi2M2", dAdaptiveArrayM2);

  fSmearNDBDM2->Rebin(dAdaptiveBinsM2, "fnewNDBDM2", dAdaptiveArrayM2);
  fSmear2NDBDM2->Rebin(dAdaptiveBinsM2, "fnew2NDBDM2", dAdaptiveArrayM2);
  fSmearBiM2->Rebin(dAdaptiveBinsM2, "fnewBiM2", dAdaptiveArrayM2);



  for(int i = 1; i <= dAdaptiveBinsM1; i++)
  {
    fAdapFrameThM1->SetBinContent(i, 2*fnewFrameThM1->GetBinContent(i)/fnewFrameThM1->GetBinWidth(i));
    fAdapTShieldThM1->SetBinContent(i, 2*fnewTShieldThM1->GetBinContent(i)/fnewTShieldThM1->GetBinWidth(i));
    fAdap50mKThM1->SetBinContent(i, 2*fnew50mKThM1->GetBinContent(i)/fnew50mKThM1->GetBinWidth(i));
    fAdap600mKThM1->SetBinContent(i, 2*fnew600mKThM1->GetBinContent(i)/fnew600mKThM1->GetBinWidth(i));
    fAdapIVCThM1->SetBinContent(i, 2*fnewIVCThM1->GetBinContent(i)/fnewIVCThM1->GetBinWidth(i));
    fAdapOVCThM1->SetBinContent(i, 2*fnewOVCThM1->GetBinContent(i)/fnewOVCThM1->GetBinWidth(i));

    fAdapFrameRaM1->SetBinContent(i, 2*fnewFrameRaM1->GetBinContent(i)/fnewFrameRaM1->GetBinWidth(i));
    fAdapTShieldRaM1->SetBinContent(i, 2*fnewTShieldRaM1->GetBinContent(i)/fnewTShieldRaM1->GetBinWidth(i));
    fAdap50mKRaM1->SetBinContent(i, 2*fnew50mKRaM1->GetBinContent(i)/fnew50mKRaM1->GetBinWidth(i));
    fAdap600mKRaM1->SetBinContent(i, 2*fnew600mKRaM1->GetBinContent(i)/fnew600mKRaM1->GetBinWidth(i));
    fAdapIVCRaM1->SetBinContent(i, 2*fnewIVCRaM1->GetBinContent(i)/fnewIVCRaM1->GetBinWidth(i));
    fAdapOVCRaM1->SetBinContent(i, 2*fnewOVCRaM1->GetBinContent(i)/fnewOVCRaM1->GetBinWidth(i));

    fAdapFrameKM1->SetBinContent(i, 2*fnewFrameKM1->GetBinContent(i)/fnewFrameKM1->GetBinWidth(i));
    fAdapTShieldKM1->SetBinContent(i, 2*fnewTShieldKM1->GetBinContent(i)/fnewTShieldKM1->GetBinWidth(i));
    fAdap50mKKM1->SetBinContent(i, 2*fnew50mKKM1->GetBinContent(i)/fnew50mKKM1->GetBinWidth(i));
    fAdap600mKKM1->SetBinContent(i, 2*fnew600mKKM1->GetBinContent(i)/fnew600mKKM1->GetBinWidth(i));
    fAdapIVCKM1->SetBinContent(i, 2*fnewIVCKM1->GetBinContent(i)/fnewIVCKM1->GetBinWidth(i));
    fAdapOVCKM1->SetBinContent(i, 2*fnewOVCKM1->GetBinContent(i)/fnewOVCKM1->GetBinWidth(i));

    fAdapFrameCoM1->SetBinContent(i, 2*fnewFrameCoM1->GetBinContent(i)/fnewFrameCoM1->GetBinWidth(i));
    fAdapTShieldCoM1->SetBinContent(i, 2*fnewTShieldCoM1->GetBinContent(i)/fnewTShieldCoM1->GetBinWidth(i));
    fAdap50mKCoM1->SetBinContent(i, 2*fnew50mKCoM1->GetBinContent(i)/fnew50mKCoM1->GetBinWidth(i));
    fAdap600mKCoM1->SetBinContent(i, 2*fnew600mKCoM1->GetBinContent(i)/fnew600mKCoM1->GetBinWidth(i));
    fAdapIVCCoM1->SetBinContent(i, 2*fnewIVCCoM1->GetBinContent(i)/fnewIVCCoM1->GetBinWidth(i));
    fAdapOVCCoM1->SetBinContent(i, 2*fnewOVCCoM1->GetBinContent(i)/fnewOVCCoM1->GetBinWidth(i));

    fAdapTShieldMnM1->SetBinContent(i, 2*fnewTShieldMnM1->GetBinContent(i)/fnewTShieldMnM1->GetBinWidth(i));
    fAdapIVCMnM1->SetBinContent(i, 2*fnewIVCMnM1->GetBinContent(i)/fnewIVCMnM1->GetBinWidth(i));

    fAdapCrystalBi2M1->SetBinContent(i, 2*fnewCrystalBi2M1->GetBinContent(i)/fnewCrystalBi2M1->GetBinWidth(i));
    fAdapFrameBi2M1->SetBinContent(i, 2*fnewFrameBi2M1->GetBinContent(i)/fnewFrameBi2M1->GetBinWidth(i));

    fAdapNDBDM1->SetBinContent(i, 2*fnewNDBDM1->GetBinContent(i)/fnewNDBDM1->GetBinWidth(i));
    fAdap2NDBDM1->SetBinContent(i, 2*fnew2NDBDM1->GetBinContent(i)/fnew600mKThM1->GetBinWidth(i));
    fAdapBiM1->SetBinContent(i, 2*fnewBiM1->GetBinContent(i)/fnewBiM1->GetBinWidth(i));


  }

  for(int i = 1; i <= dAdaptiveBinsM2; i++)
  {
    fAdapFrameThM2->SetBinContent(i, 2*fnewFrameThM2->GetBinContent(i)/fnewFrameThM2->GetBinWidth(i));
    fAdapTShieldThM2->SetBinContent(i, 2*fnewTShieldThM2->GetBinContent(i)/fnewTShieldThM2->GetBinWidth(i));
    fAdap50mKThM2->SetBinContent(i, 2*fnew50mKThM2->GetBinContent(i)/fnew50mKThM2->GetBinWidth(i));
    fAdap600mKThM2->SetBinContent(i, 2*fnew600mKThM2->GetBinContent(i)/fnew600mKThM2->GetBinWidth(i));
    fAdapIVCThM2->SetBinContent(i, 2*fnewIVCThM2->GetBinContent(i)/fnewIVCThM2->GetBinWidth(i));
    fAdapOVCThM2->SetBinContent(i, 2*fnewOVCThM2->GetBinContent(i)/fnewOVCThM2->GetBinWidth(i));

    fAdapFrameRaM2->SetBinContent(i, 2*fnewFrameRaM2->GetBinContent(i)/fnewFrameRaM2->GetBinWidth(i));
    fAdapTShieldRaM2->SetBinContent(i, 2*fnewTShieldRaM2->GetBinContent(i)/fnewTShieldRaM2->GetBinWidth(i));
    fAdap50mKRaM2->SetBinContent(i, 2*fnew50mKRaM2->GetBinContent(i)/fnew50mKRaM2->GetBinWidth(i));
    fAdap600mKRaM2->SetBinContent(i, 2*fnew600mKRaM2->GetBinContent(i)/fnew600mKRaM2->GetBinWidth(i));
    fAdapIVCRaM2->SetBinContent(i, 2*fnewIVCRaM2->GetBinContent(i)/fnewIVCRaM2->GetBinWidth(i));
    fAdapOVCRaM2->SetBinContent(i, 2*fnewOVCRaM2->GetBinContent(i)/fnewOVCRaM2->GetBinWidth(i));

    fAdapFrameKM2->SetBinContent(i, 2*fnewFrameKM2->GetBinContent(i)/fnewFrameKM2->GetBinWidth(i));
    fAdapTShieldKM2->SetBinContent(i, 2*fnewTShieldKM2->GetBinContent(i)/fnewTShieldKM2->GetBinWidth(i));
    fAdap50mKKM2->SetBinContent(i, 2*fnew50mKKM2->GetBinContent(i)/fnew50mKKM2->GetBinWidth(i));
    fAdap600mKKM2->SetBinContent(i, 2*fnew600mKKM2->GetBinContent(i)/fnew600mKKM2->GetBinWidth(i));
    fAdapIVCKM2->SetBinContent(i, 2*fnewIVCKM2->GetBinContent(i)/fnewIVCKM2->GetBinWidth(i));
    fAdapOVCKM2->SetBinContent(i, 2*fnewOVCKM2->GetBinContent(i)/fnewOVCKM2->GetBinWidth(i));

    fAdapFrameCoM2->SetBinContent(i, 2*fnewFrameCoM2->GetBinContent(i)/fnewFrameCoM2->GetBinWidth(i));
    fAdapTShieldCoM2->SetBinContent(i, 2*fnewTShieldCoM2->GetBinContent(i)/fnewTShieldCoM2->GetBinWidth(i));
    fAdap50mKCoM2->SetBinContent(i, 2*fnew50mKCoM2->GetBinContent(i)/fnew50mKCoM2->GetBinWidth(i));
    fAdap600mKCoM2->SetBinContent(i, 2*fnew600mKCoM2->GetBinContent(i)/fnew600mKCoM2->GetBinWidth(i));
    fAdapIVCCoM2->SetBinContent(i, 2*fnewIVCCoM2->GetBinContent(i)/fnewIVCCoM2->GetBinWidth(i));
    fAdapOVCCoM2->SetBinContent(i, 2*fnewOVCCoM2->GetBinContent(i)/fnewOVCCoM2->GetBinWidth(i));

    fAdapTShieldMnM2->SetBinContent(i, 2*fnewTShieldMnM2->GetBinContent(i)/fnewTShieldMnM2->GetBinWidth(i));
    fAdapIVCMnM2->SetBinContent(i, 2*fnewIVCMnM2->GetBinContent(i)/fnewIVCMnM2->GetBinWidth(i));

    fAdapCrystalBi2M2->SetBinContent(i, 2*fnewCrystalBi2M2->GetBinContent(i)/fnewCrystalBi2M2->GetBinWidth(i));
    fAdapFrameBi2M2->SetBinContent(i, 2*fnewFrameBi2M2->GetBinContent(i)/fnewFrameBi2M2->GetBinWidth(i));

    fAdapNDBDM2->SetBinContent(i, 2*fnewNDBDM2->GetBinContent(i)/fnewNDBDM2->GetBinWidth(i));
    fAdap2NDBDM2->SetBinContent(i, 2*fnew2NDBDM2->GetBinContent(i)/fnew600mKThM2->GetBinWidth(i));
    fAdapBiM2->SetBinContent(i, 2*fnewBiM2->GetBinContent(i)/fnewBiM2->GetBinWidth(i));

  }



  TFile *file1 = new TFile("MCAdaptiveGamma-new.root", "RECREATE");


    fSmearFrameThM1->Write();
    fSmearTShieldThM1->Write();  
    fSmear50mKThM1->Write();
    fSmear600mKThM1->Write();
    fSmearIVCThM1->Write();
    fSmearOVCThM1->Write();

    fSmearFrameRaM1->Write();
    fSmearTShieldRaM1->Write();  
    fSmear50mKRaM1->Write();
    fSmear600mKRaM1->Write();
    fSmearIVCRaM1->Write();
    fSmearOVCRaM1->Write();

    fSmearFrameKM1->Write();
    fSmearTShieldKM1->Write();
    fSmear50mKKM1->Write();
    fSmear600mKKM1->Write();
    fSmearIVCKM1->Write();
    fSmearOVCKM1->Write(); 

    fSmearFrameCoM1->Write();
    fSmearTShieldCoM1->Write();
    fSmear50mKCoM1->Write();
    fSmear600mKCoM1->Write();
    fSmearIVCCoM1->Write();
    fSmearOVCCoM1->Write();  

    fSmearTShieldMnM1->Write();
    fSmearIVCMnM1->Write();

    fSmearCrystalBi2M1->Write();
    fSmearFrameBi2M1->Write();

    fSmearNDBDM1->Write();  
    fSmear2NDBDM1->Write();  
    fSmearBiM1->Write();  



    fSmearFrameThM2->Write();
    fSmearTShieldThM2->Write();  
    fSmear50mKThM2->Write();
    fSmear600mKThM2->Write();
    fSmearIVCThM2->Write();
    fSmearOVCThM2->Write();

    fSmearFrameRaM2->Write();
    fSmearTShieldRaM2->Write();  
    fSmear50mKRaM2->Write();
    fSmear600mKRaM2->Write();
    fSmearIVCRaM2->Write();
    fSmearOVCRaM2->Write();

    fSmearFrameKM2->Write();
    fSmearTShieldKM2->Write();
    fSmear50mKKM2->Write();
    fSmear600mKKM2->Write();
    fSmearIVCKM2->Write();
    fSmearOVCKM2->Write(); 

    fSmearFrameCoM2->Write();
    fSmearTShieldCoM2->Write();
    fSmear50mKCoM2->Write();
    fSmear600mKCoM2->Write();
    fSmearIVCCoM2->Write();
    fSmearOVCCoM2->Write();  

    fSmearTShieldMnM2->Write();
    fSmearIVCMnM2->Write();

    fSmearCrystalBi2M2->Write();
    fSmearFrameBi2M2->Write();

    fSmearNDBDM2->Write();  
    fSmear2NDBDM2->Write();  
    fSmearBiM2->Write();  


    fAdapFrameThM1->Write();
    fAdapTShieldThM1->Write();
    fAdap50mKThM1->Write();
    fAdap600mKThM1->Write();
    fAdapIVCThM1->Write();
    fAdapOVCThM1->Write();

    fAdapFrameRaM1->Write();
    fAdapTShieldRaM1->Write();
    fAdap50mKRaM1->Write();
    fAdap600mKRaM1->Write();
    fAdapIVCRaM1->Write();
    fAdapOVCRaM1->Write();

    fAdapFrameKM1->Write();
    fAdapTShieldKM1->Write();
    fAdap50mKKM1->Write();
    fAdap600mKKM1->Write();
    fAdapIVCKM1->Write();
    fAdapOVCKM1->Write();

    fAdapFrameCoM1->Write();
    fAdapTShieldCoM1->Write();
    fAdap50mKCoM1->Write();
    fAdap600mKCoM1->Write();
    fAdapIVCCoM1->Write();
    fAdapOVCCoM1->Write();

    fAdapTShieldMnM1->Write();
    fAdapIVCMnM1->Write();

    fAdapCrystalBi2M1->Write();
    fAdapFrameBi2M1->Write();

    fAdapNDBDM1->Write();
    fAdap2NDBDM1->Write();
    fAdapBiM1->Write();


    fAdapFrameThM2->Write();
    fAdapTShieldThM2->Write();
    fAdap50mKThM2->Write();
    fAdap600mKThM2->Write();
    fAdapIVCThM2->Write();
    fAdapOVCThM2->Write();

    fAdapFrameRaM2->Write();
    fAdapTShieldRaM2->Write();
    fAdap50mKRaM2->Write();
    fAdap600mKRaM2->Write();
    fAdapIVCRaM2->Write();
    fAdapOVCRaM2->Write();

    fAdapFrameKM2->Write();
    fAdapTShieldKM2->Write();
    fAdap50mKKM2->Write();
    fAdap600mKKM2->Write();
    fAdapIVCKM2->Write();
    fAdapOVCKM2->Write();

    fAdapFrameCoM2->Write();
    fAdapTShieldCoM2->Write();
    fAdap50mKCoM2->Write();
    fAdap600mKCoM2->Write();
    fAdapIVCCoM2->Write();
    fAdapOVCCoM2->Write();

    fAdapTShieldMnM2->Write();
    fAdapIVCMnM2->Write();

    fAdapCrystalBi2M2->Write();
    fAdapFrameBi2M2->Write();

    fAdapNDBDM2->Write();
    fAdap2NDBDM2->Write();
    fAdapBiM2->Write();

    file1->Write();

}

// Saves alpha histograms
void TBackgroundModel::SaveAlpha()
{
  bUnSmeared = true;
  // First initialize
  Initialize();

  // Now store data
  TFile *file1 = new TFile("MCAdaptiveAlpha-new.root", "RECREATE");

  // M1
  fSmearCrystalPt190BM1->Write();
  fSmearCrystalPt190S01M1->Write();
  fSmearCrystalPt190S1M1->Write();
  fSmearCrystalPt190S10M1->Write();
  fSmearCrystalPt190S100M1->Write();

  fSmearCrystalPb210BM1->Write();
  fSmearCrystalPb210S01M1->Write();
  fSmearCrystalPb210S1M1->Write();
  fSmearCrystalPb210S10M1->Write();
  fSmearCrystalPb210S100M1->Write();

  fSmearCrystalRa226BM1->Write();
  fSmearCrystalRa226S01M1->Write();
  fSmearCrystalRa226S1M1->Write();
  fSmearCrystalRa226S10M1->Write();
  fSmearCrystalRa226S100M1->Write();

  fSmearFramePb210BM1->Write();
  fSmearFramePb210S01M1->Write();
  fSmearFramePb210S1M1->Write();
  fSmearFramePb210S10M1->Write();
  fSmearFramePb210S100M1->Write();

  fSmearFrameRa226BM1->Write();
  fSmearFrameRa226S01M1->Write();
  fSmearFrameRa226S1M1->Write();
  fSmearFrameRa226S10M1->Write();
  fSmearFrameRa226S100M1->Write();


  // M2
  fSmearCrystalPt190BM2->Write();
  fSmearCrystalPt190S01M2->Write();
  fSmearCrystalPt190S1M2->Write();
  fSmearCrystalPt190S10M2->Write();
  fSmearCrystalPt190S100M2->Write();

  fSmearCrystalPb210BM2->Write();
  fSmearCrystalPb210S01M2->Write();
  fSmearCrystalPb210S1M2->Write();
  fSmearCrystalPb210S10M2->Write();
  fSmearCrystalPb210S100M2->Write();

  fSmearCrystalRa226BM2->Write();
  fSmearCrystalRa226S01M2->Write();
  fSmearCrystalRa226S1M2->Write();
  fSmearCrystalRa226S10M2->Write();
  fSmearCrystalRa226S100M2->Write();

  fSmearFramePb210BM2->Write();
  fSmearFramePb210S01M2->Write();
  fSmearFramePb210S1M2->Write();
  fSmearFramePb210S10M2->Write();
  fSmearFramePb210S100M2->Write();

  fSmearFrameRa226BM2->Write();
  fSmearFrameRa226S01M2->Write();
  fSmearFrameRa226S1M2->Write();
  fSmearFrameRa226S10M2->Write();
  fSmearFrameRa226S100M2->Write();

  // Adaptive
  fAdapCrystalPt190BM1->Write();
  fAdapCrystalPt190S01M1->Write();
  fAdapCrystalPt190S1M1->Write();
  fAdapCrystalPt190S10M1->Write();
  fAdapCrystalPt190S100M1->Write();

  fAdapCrystalPb210BM1->Write();
  fAdapCrystalPb210S01M1->Write();
  fAdapCrystalPb210S1M1->Write();
  fAdapCrystalPb210S10M1->Write();
  fAdapCrystalPb210S100M1->Write();

  fAdapCrystalRa226BM1->Write();
  fAdapCrystalRa226S01M1->Write();
  fAdapCrystalRa226S1M1->Write();
  fAdapCrystalRa226S10M1->Write();
  fAdapCrystalRa226S100M1->Write();

  fAdapFramePb210BM1->Write();
  fAdapFramePb210S01M1->Write();
  fAdapFramePb210S1M1->Write();
  fAdapFramePb210S10M1->Write();
  fAdapFramePb210S100M1->Write();

  fAdapFrameRa226BM1->Write();
  fAdapFrameRa226S01M1->Write();
  fAdapFrameRa226S1M1->Write();
  fAdapFrameRa226S10M1->Write();
  fAdapFrameRa226S100M1->Write();



  fAdapCrystalPt190BM2->Write();
  fAdapCrystalPt190S01M2->Write();
  fAdapCrystalPt190S1M2->Write();
  fAdapCrystalPt190S10M2->Write();
  fAdapCrystalPt190S100M2->Write();

  fAdapCrystalPb210BM2->Write();
  fAdapCrystalPb210S01M2->Write();
  fAdapCrystalPb210S1M2->Write();
  fAdapCrystalPb210S10M2->Write();
  fAdapCrystalPb210S100M2->Write();

  fAdapCrystalRa226BM2->Write();
  fAdapCrystalRa226S01M2->Write();
  fAdapCrystalRa226S1M2->Write();
  fAdapCrystalRa226S10M2->Write();
  fAdapCrystalRa226S100M2->Write();

  fAdapFramePb210BM2->Write();
  fAdapFramePb210S01M2->Write();
  fAdapFramePb210S1M2->Write();
  fAdapFramePb210S10M2->Write();
  fAdapFramePb210S100M2->Write();

  fAdapFrameRa226BM2->Write();
  fAdapFrameRa226S01M2->Write();
  fAdapFrameRa226S1M2->Write();
  fAdapFrameRa226S10M2->Write();
  fAdapFrameRa226S100M2->Write();
    file1->Write();

}

// For whatever tests...
void TBackgroundModel::Test()
{ 
  // vector<double> Test;

  // Test = AdaptiveBinning(fDataHistoM1);

  // int bins = Test.size();
  gStyle->SetOptStat(0);

  cout << "Old Number of bins: " << dNBins << " New number of bins (M1): " << dAdaptiveBinsM1 << " New number of bins (M2): " << dAdaptiveBinsM2 << endl;

  double dFill = 0;

  fDataHistoM1->Rebin(dAdaptiveBinsM1, "hnewM1", dAdaptiveArrayM1);
  fDataHistoM2->Rebin(dAdaptiveBinsM2, "hnewM2", dAdaptiveArrayM2);



  TH1D *hAdjustedM1 = new TH1D("hAdjustedM1", "M1 spectrum with adaptive binning", dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *hAdjustedM2 = new TH1D("hAdjustedM2", "M2 spectrum with adaptive binning", dAdaptiveBinsM2, dAdaptiveArrayM2);



  for(int i = 1; i <= dAdaptiveBinsM1; i++)
  {
    dFill = 2*hnewM1->GetBinContent(i)/hnewM1->GetBinWidth(i);
    hAdjustedM1->SetBinContent(i, dFill);
    hAdjustedM1->SetBinError(i, TMath::Sqrt(2*hnewM1->GetBinContent(i))/hnewM1->GetBinWidth(i));
  }

  for(int i = 1; i <= dAdaptiveBinsM2; i++)
  {
    dFill = 2*hnewM2->GetBinContent(i)/hnewM2->GetBinWidth(i);
    hAdjustedM2->SetBinContent(i, dFill);    
    hAdjustedM2->SetBinError(i, TMath::Sqrt(2*hnewM2->GetBinContent(i))/hnewM2->GetBinWidth(i));
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1200);
  c1->Divide(1, 2);
  c1->cd(1);
  c1->SetLogy();
  // fDataHistoM1->Draw("E");
  // hAdjustedM1->SetLineColor(kRed);
  // hAdjustedM1->SetLineStyle(2);
  hAdjustedM1->Draw("E");
  // hnew->SetLineColor(kRed);
  // hnew->SetLineStyle(2);
  // hnew->Draw("ESAME");

  // TCanvas *c2 = new TCanvas("c2");
  c1->cd(2); 
  c1->SetLogy();
  // fDataHistoM2->Draw("E");
  // hAdjustedM2->SetLineColor(kRed);
  // hAdjustedM2->SetLineStyle(2);
  hAdjustedM2->Draw("E");
}




bool TBackgroundModel::DoTheFitAdaptive()
{
  if(!bToyFit)
  {
    Initialize();
  }
  gStyle->SetOptStat(0);
   // This method actually sets up minuit and does the fit


   // TMinuit minuit(14); //initialize minuit, n is the number of parameters
   TMinuit minuit(26); // for more parameters

   // Reduce Minuit Output
   minuit.SetPrintLevel(1);
//   minuit.Command("SET MINImize 1000 0.001");
   minuit.Command("SET STRategy 2");
  //minuit.Command("SET IMProve 1000 ");
   minuit.SetMaxIterations(1000);
   minuit.SetObjectFit(this); //see the external FCN  above
   
   //define the parameters and set the ranges and initial guesses see ROOTs TMinuit documentation

   // Close Th and Close Ra now split into its various sections, far Th and Ra still the same
   // This step after previous step full fit converges, just to see if any differences show up
   ////////////////////////////////////////////////
   // Using more parameters
   ////////////////////////////////////////////////

   minuit.DefineParameter(0, "Frame Th",  20000., 10.0, 0., 100000);
   minuit.DefineParameter(1, "TShield Th",    380000., 10.0, 0., 500000);
   minuit.DefineParameter(2, "Frame Ra",  1100., 10.0, 0., 80000);   
   minuit.DefineParameter(3, "TShield Ra", 1100., 10.0, 0., 80000);
   minuit.DefineParameter(4, "Close K",   21732.5, 1.0, 0., 500000);
   minuit.DefineParameter(5, "Far K",     3463.03, 10.0, 0., 500000);
   minuit.DefineParameter(6, "Frame Co",  1828.31, 1.0, 0., 80000);    
   minuit.DefineParameter(7, "TShield Co",    1100, 10.0, 0., 80000);  
   // minuit.DefineParameter(8, "2NDBD",    33394., 10.0, 0., 100000);        
   minuit.DefineParameter(8, "2NDBD",    53000.6., 10.0, 0., 100000);   
   minuit.DefineParameter(9, "NDBD",       92.45, 1.0, 0., 500);     
   minuit.DefineParameter(10, "Lead Bi",    7723.64, 10.0, 0., 100000);  
   minuit.DefineParameter(11, "TShield Mn",    20000., 10.0, 0., 100000);  
   minuit.DefineParameter(12, "IVC Mn",    5435.33, 10.0, 0., 100000);  
   minuit.DefineParameter(13, "50mK Th",    380000., 10.0, 0., 500000);
   minuit.DefineParameter(14, "600mK Th",    37952.8, 10.0, 0., 500000);
   minuit.DefineParameter(15, "IVC Th",    42461.1, 10.0, 0., 500000);
   minuit.DefineParameter(16, "OVC Th",    16575.1, 10.0, 0., 500000);
   minuit.DefineParameter(17, "50mK Ra", 1100., 10.0, 0., 80000);
   minuit.DefineParameter(18, "600mK Ra", 9666.94, 10.0, 0., 80000);
   minuit.DefineParameter(19, "IVC Ra",    380000., 10.0, 0., 500000);
   minuit.DefineParameter(20, "OVC Ra",    106770., 10.0, 0., 500000);
   minuit.DefineParameter(21, "50mK Co",  40000., 1.0, 0., 80000);    
   minuit.DefineParameter(22, "600mK Co",    1100., 10.0, 0., 80000);
   minuit.DefineParameter(23, "IVC Co",    380000, 10.0, 0., 500000);  
   minuit.DefineParameter(24, "OVC Co",    20815.6, 10.0, 0., 500000);  
   // minuit.DefineParameter(25, "Constant",    500000, 1, 0., 1000000);  
   minuit.DefineParameter(25, "Pb-210 chain",    0, 1, 0., 1000000);  

   



   // Fix parameters here
   // minuit.FixParameter(0); // Frame Th
   // minuit.FixParameter(1); // TShield Th
   // minuit.FixParameter(2); // Frame Ra
   // minuit.FixParameter(3); // TShield Ra
   // minuit.FixParameter(4); // Close K
   // minuit.FixParameter(5); // Far K
   // minuit.FixParameter(6); // Frame Co
   // minuit.FixParameter(7); // TShield Co
   // minuit.FixParameter(8); // 2NDBD
   // minuit.FixParameter(9); // NDBD
   // minuit.FixParameter(10); // Bi207
   // minuit.FixParameter(11); // TShield Mn
   // minuit.FixParameter(12); // IVC Mn
   // minuit.FixParameter(13); // 50mK Th
   // minuit.FixParameter(14); // 600mK Th
   // minuit.FixParameter(15); // IVC Th
   // minuit.FixParameter(16); // OVC Th
   // minuit.FixParameter(17); // 50mK Ra
   // minuit.FixParameter(18); // 600mK Ra
   // minuit.FixParameter(19); // IVC Ra
   // minuit.FixParameter(20); // OVC Ra
   // minuit.FixParameter(21); // 50mK Co
   // minuit.FixParameter(22); // 600 mK Co
   // minuit.FixParameter(23); // IVC Co
   // minuit.FixParameter(24); // OVC Co
   minuit.FixParameter(25); // 

  // Number of Parameters! (for Chi-squared/NDF calculation)
  int dNumParameters = 25;




   //Tell minuit what external function to use 
   minuit.SetFCN(myExternal_FCN);
   
   int status = minuit.Migrad(); // this actually does the minimisation


  minuit.GetParameter(0,  fParameters[0],   fParError[0]);
  minuit.GetParameter(1,  fParameters[1],   fParError[1]);
  minuit.GetParameter(2,  fParameters[2],   fParError[2]);
  minuit.GetParameter(3,  fParameters[3],   fParError[3]);
  minuit.GetParameter(4,  fParameters[4],   fParError[4]);
  minuit.GetParameter(5,  fParameters[5],   fParError[5]);
  minuit.GetParameter(6,  fParameters[6],   fParError[6]);
  minuit.GetParameter(7,  fParameters[7],   fParError[7]);
  minuit.GetParameter(8,  fParameters[8],   fParError[8]);
  minuit.GetParameter(9,  fParameters[9],   fParError[9]);
  minuit.GetParameter(10,  fParameters[10],   fParError[10]);
  minuit.GetParameter(11,  fParameters[11],   fParError[11]);
  minuit.GetParameter(12,  fParameters[12],   fParError[12]);
  minuit.GetParameter(13,  fParameters[13],   fParError[13]);
  minuit.GetParameter(14,  fParameters[14],   fParError[14]);
  minuit.GetParameter(15,  fParameters[15],   fParError[15]);
  minuit.GetParameter(16,  fParameters[16],   fParError[16]);
  minuit.GetParameter(17,  fParameters[17],   fParError[17]);
  minuit.GetParameter(18,  fParameters[18],   fParError[18]);
  minuit.GetParameter(19,  fParameters[19],   fParError[19]);
  minuit.GetParameter(20,  fParameters[20],   fParError[20]);
  minuit.GetParameter(21,  fParameters[21],   fParError[21]);
  minuit.GetParameter(22,  fParameters[22],   fParError[22]);
  minuit.GetParameter(23,  fParameters[23],   fParError[23]);
  minuit.GetParameter(24,  fParameters[24],   fParError[24]);
  minuit.GetParameter(25,  fParameters[25],   fParError[25]);
  UpdateModelAdaptive();
  
  cout << "At the end; ChiSq/NDF = " << GetChiSquareAdaptive()/(dFitMaxBinM1+dFitMaxBinM2-dFitMinBinM1-dFitMinBinM2-dNumParameters) << endl; // for M1 and M2
  // cout << "At the end; ChiSq/NDF = " << GetChiSquare()/((dFitMax-dFitMin)/dBinSize - dNumParameters) << endl;  // for M1
  cout << "Total number of calls = " << dNumCalls << endl;

  ///////////////////////////////////////////
  //// Many Parameters
  ///////////////////////////////////////////
  /// Use only after previous step converges!
  // 

  // M1 Parameters
  fModelTotAdapThM1->Add(fAdapFrameThM1,   fParameters[0]);
  fModelTotAdapThM1->Add(fAdapTShieldThM1, fParameters[1]);
  fModelTotAdapThM1->Add(fAdap50mKThM1,    fParameters[13]);
  fModelTotAdapThM1->Add(fAdap600mKThM1,   fParameters[14]);
  fModelTotAdapThM1->Add(fAdapIVCThM1,     fParameters[15]);
  fModelTotAdapThM1->Add(fAdapOVCThM1,     fParameters[16]);

  fModelTotAdapRaM1->Add(fAdapFrameRaM1,   fParameters[2]);
  fModelTotAdapRaM1->Add(fAdapTShieldRaM1, fParameters[3]);
  fModelTotAdapRaM1->Add(fAdap50mKRaM1,    fParameters[17]);
  fModelTotAdapRaM1->Add(fAdap600mKRaM1,   fParameters[18]);
  fModelTotAdapRaM1->Add(fAdapIVCRaM1,     fParameters[19]);
  fModelTotAdapRaM1->Add(fAdapOVCRaM1,     fParameters[20]);

  fModelTotAdapKM1->Add(fAdapFrameKM1,     fParameters[4]);
  fModelTotAdapKM1->Add(fAdapTShieldKM1,   fParameters[4]);
  fModelTotAdapKM1->Add(fAdap50mKKM1,      fParameters[4]);
  fModelTotAdapKM1->Add(fAdap600mKKM1,     fParameters[4]);
  fModelTotAdapKM1->Add(fAdapIVCKM1,       fParameters[5]);
  fModelTotAdapKM1->Add(fAdapOVCKM1,       fParameters[5]);

  fModelTotAdapCoM1->Add(fAdapFrameCoM1,   fParameters[6]);
  fModelTotAdapCoM1->Add(fAdapTShieldCoM1, fParameters[7]);
  fModelTotAdapCoM1->Add(fAdap50mKCoM1,    fParameters[21]);
  fModelTotAdapCoM1->Add(fAdap600mKCoM1,   fParameters[22]);
  fModelTotAdapCoM1->Add(fAdapIVCCoM1,     fParameters[23]);
  fModelTotAdapCoM1->Add(fAdapOVCCoM1,     fParameters[24]);

  fModelTotAdapMnM1->Add(fAdapTShieldMnM1, fParameters[11]);
  fModelTotAdapMnM1->Add(fAdapIVCMnM1,     fParameters[12]);


  fModelTotAdapNDBDM1->Add(fAdapNDBDM1,    fParameters[9]);
  fModelTotAdap2NDBDM1->Add(fAdap2NDBDM1,  fParameters[8]);
  fModelTotAdapBiM1->Add(fAdapBiM1,        fParameters[10]);

  // fModelTotAdapPbM1->Add(fAdapCrystalPbBM1, fParameters[25]);


  // M2 Parameters
  fModelTotAdapThM2->Add(fAdapFrameThM2,   fParameters[0]);
  fModelTotAdapThM2->Add(fAdapTShieldThM2, fParameters[1]);
  fModelTotAdapThM2->Add(fAdap50mKThM2,    fParameters[13]);
  fModelTotAdapThM2->Add(fAdap600mKThM2,   fParameters[14]);
  fModelTotAdapThM2->Add(fAdapIVCThM2,     fParameters[15]);
  fModelTotAdapThM2->Add(fAdapOVCThM2,     fParameters[16]);

  fModelTotAdapRaM2->Add(fAdapFrameRaM2,   fParameters[2]);
  fModelTotAdapRaM2->Add(fAdapTShieldRaM2, fParameters[3]);
  fModelTotAdapRaM2->Add(fAdap50mKRaM2,    fParameters[17]);
  fModelTotAdapRaM2->Add(fAdap600mKRaM2,   fParameters[18]);
  fModelTotAdapRaM2->Add(fAdapIVCRaM2,     fParameters[19]);
  fModelTotAdapRaM2->Add(fAdapOVCRaM2,     fParameters[20]);

  fModelTotAdapKM2->Add(fAdapFrameKM2,     fParameters[4]);
  fModelTotAdapKM2->Add(fAdapTShieldKM2,   fParameters[4]);
  fModelTotAdapKM2->Add(fAdap50mKKM2,      fParameters[4]);
  fModelTotAdapKM2->Add(fAdap600mKKM2,     fParameters[4]);
  fModelTotAdapKM2->Add(fAdapIVCKM2,       fParameters[5]);
  fModelTotAdapKM2->Add(fAdapOVCKM2,       fParameters[5]);

  fModelTotAdapCoM2->Add(fAdapFrameCoM2,   fParameters[6]);
  fModelTotAdapCoM2->Add(fAdapTShieldCoM2, fParameters[7]);
  fModelTotAdapCoM2->Add(fAdap50mKCoM2,    fParameters[21]);
  fModelTotAdapCoM2->Add(fAdap600mKCoM2,   fParameters[22]);
  fModelTotAdapCoM2->Add(fAdapIVCCoM2,     fParameters[23]);
  fModelTotAdapCoM2->Add(fAdapOVCCoM2,     fParameters[24]);

  fModelTotAdapMnM2->Add(fAdapTShieldMnM2, fParameters[11]);
  fModelTotAdapMnM2->Add(fAdapIVCMnM2,     fParameters[12]);

  fModelTotAdapNDBDM2->Add(fAdapNDBDM2,    fParameters[9]);
  fModelTotAdap2NDBDM2->Add(fAdap2NDBDM2,  fParameters[8]);
  fModelTotAdapBiM2->Add(fAdapBiM2,      fParameters[10]);

  // fModelTotPbM2->Add(fAdapCrystalPbBM2, fParameters[25]);

  ////////// Only for testing
  // Correction for M2 spectra, it's the M1 spectra but scaled down by N_M1*1-Exp(R*T)
  // fTotCorrection->Add(fCorrectionM2, 180197*(1-TMath::Exp(-0.05*0.1)));


    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->SetLogy();

    ///// Draw Data M1
    fAdapDataHistoM1->SetLineColor(1);
    fAdapDataHistoM1->SetLineWidth(2);
    fAdapDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
    fAdapDataHistoM1->GetYaxis()->SetTitle("Counts/yr");
    fAdapDataHistoM1->SetMaximum(90000);
    fAdapDataHistoM1->GetXaxis()->SetRange(1, 2650/dBinSize+5);
    fAdapDataHistoM1->Draw();


    fModelTotAdapM1->SetLineColor(2);
    fModelTotAdapM1->SetLineWidth(1);
    fModelTotAdapM1->Draw("SAME");

    fModelTotAdapThM1->SetLineColor(3);
    fModelTotAdapThM1->SetLineStyle(2);
    fModelTotAdapRaM1->SetLineColor(4);
    fModelTotAdapRaM1->SetLineStyle(2);
    fModelTotAdapKM1->SetLineColor(6);
    fModelTotAdapKM1->SetLineStyle(2);
    fModelTotAdapCoM1->SetLineColor(7);
    fModelTotAdapCoM1->SetLineStyle(2);
    fModelTotAdapNDBDM1->SetLineColor(42);
    fModelTotAdapNDBDM1->SetLineStyle(2);
    fModelTotAdap2NDBDM1->SetLineColor(46);
    fModelTotAdap2NDBDM1->SetLineStyle(2);
    fModelTotAdapBiM1->SetLineColor(5);
    fModelTotAdapBiM1->SetLineStyle(2);
    fModelTotAdapMnM1->SetLineColor(40);
    fModelTotAdapMnM1->SetLineStyle(2);


    fModelTotAdapPbM1->SetLineStyle(2);
    fModelTotAdapPbM1->SetLineColor(38);

    fModelTotAdapThM1->Draw("SAME");
    fModelTotAdapRaM1->Draw("SAME");
    fModelTotAdapKM1->Draw("SAME");
    fModelTotAdapCoM1->Draw("SAME");
    fModelTotAdapNDBDM1->Draw("SAME");
    fModelTotAdap2NDBDM1->Draw("SAME");
    fModelTotAdapBiM1->Draw("SAME");
    fModelTotAdapMnM1->Draw("SAME");

    // fModelTotAdapPbM1->Draw("SAME");

    // Many Parameters
    TPaveText *pt1 = new TPaveText(0.35,0.77,0.70,0.99,"NB NDC");
    pt1->AddText(Form("Fit Range (M1): %.0f to %.0f keV -- #chi^{2}/NDF: %0.3f", dFitMin, dFitMax, (GetChiSquareAdaptive()/(dFitMaxBinM1+dFitMaxBinM2-dFitMinBinM1-dFitMinBinM2-dNumParameters) )));
    pt1->AddText(Form("Frame Th: %0.2E#pm%0.2E --- TShield Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
    pt1->AddText(Form("50mK Th: %0.2E#pm%0.2E --- 600mK Th: %0.2E#pm%0.2E", fParameters[13], fParError[13], fParameters[14], fParError[14] ));
    pt1->AddText(Form("IVC Th: %0.2E#pm%0.2E --- OVC Th: %0.2E#pm%0.2E", fParameters[15], fParError[15], fParameters[16], fParError[16] ));
    pt1->AddText(Form("Frame Ra: %0.2E#pm%0.2E --- TShield Ra: %0.2E#pm%0.2E", fParameters[2], fParError[2], fParameters[3], fParError[3] ));
    pt1->AddText(Form("50mK Ra: %0.2E#pm%0.2E --- 600mK Ra: %0.2E#pm%0.2E", fParameters[17], fParError[17], fParameters[18], fParError[18] ));
    pt1->AddText(Form("IVC Ra: %0.2E#pm%0.2E --- OVC Ra: %0.2E#pm%0.2E", fParameters[19], fParError[19], fParameters[20], fParError[20] ));
    pt1->AddText(Form("Close K: %0.2E#pm%0.2E --- Far K: %0.2E#pm%0.2E", fParameters[4], fParError[4], fParameters[5], fParError[5] ));
    pt1->AddText(Form("Close Co: %0.2E#pm%0.2E --- Far Co: %0.2E#pm%0.2E", fParameters[6], fParError[6], fParameters[7], fParError[7] ));
    pt1->AddText(Form("Bi-207: %0.2E#pm%0.2E --- NDBD: %0.2E#pm%0.2E", fParameters[10], fParError[10], fParameters[9], fParError[9] ));
    pt1->AddText(Form("Close Mn-54: %0.2E#pm%0.2E --- Far Mn-54: %0.2E#pm%0.2E", fParameters[11], fParError[11], fParameters[12], fParError[12] ));
    pt1->AddText(Form("2NDBD: %0.2E#pm%0.2E" , fParameters[8], fParError[8] ));
    pt1->Draw();

    TLegend *legfit1 = new TLegend(0.8,0.8,0.97,0.97);
    legfit1->AddEntry(fModelTotAdapM1, "Total PDF", "l");
    legfit1->AddEntry(fModelTotAdapThM1, "Total Th-232", "l");
    legfit1->AddEntry(fModelTotAdapRaM1, "Total Ra-226", "l");
    legfit1->AddEntry(fModelTotAdapKM1, "Total K-40", "l");
    legfit1->AddEntry(fModelTotAdapCoM1, "Total Co-60", "l");
    legfit1->AddEntry(fModelTotAdapNDBDM1, "NDBD", "l");
    legfit1->AddEntry(fModelTotAdap2NDBDM1, "2NDBD", "l");
    legfit1->AddEntry(fModelTotAdapBiM1, "Bi-207", "l");
    legfit1->AddEntry(fModelTotAdapMnM1, "Mn-54", "l");
    legfit1->Draw();





    TCanvas *c2 = new TCanvas("c2", "c2", 1200, 800);
    c2->SetLogy();

    ///// Draw Data M2
    fAdapDataHistoM2->SetLineColor(1);
    fAdapDataHistoM2->SetLineWidth(2);
    fAdapDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
    fAdapDataHistoM2->GetYaxis()->SetTitle("Counts/yr");
    fAdapDataHistoM2->SetMaximum(9000);
    fAdapDataHistoM2->GetXaxis()->SetRange(1/dBinSize-5, 2650/dBinSize+5);
    fAdapDataHistoM2->Draw();

  
    fModelTotAdapM2->SetLineColor(2);
    fModelTotAdapM2->SetLineWidth(1);
    fModelTotAdapM2->Draw("SAME");

    fModelTotAdapThM2->SetLineColor(3);
    fModelTotAdapThM2->SetLineStyle(2);
    fModelTotAdapRaM2->SetLineColor(4);
    fModelTotAdapRaM2->SetLineStyle(2);
    fModelTotAdapKM2->SetLineColor(6);
    fModelTotAdapKM2->SetLineStyle(2);
    fModelTotAdapCoM2->SetLineColor(7);
    fModelTotAdapCoM2->SetLineStyle(2);
    fModelTotAdapNDBDM2->SetLineColor(42);
    fModelTotAdapNDBDM2->SetLineStyle(2);
    fModelTotAdap2NDBDM2->SetLineColor(46);
    fModelTotAdap2NDBDM2->SetLineStyle(2);
    fModelTotAdapBiM2->SetLineColor(5);
    fModelTotAdapBiM2->SetLineStyle(2);
    fModelTotAdapMnM2->SetLineColor(40);
    fModelTotAdapMnM2->SetLineStyle(2);

    // fTotCorrection->SetLineStyle(2);
    // fTotCorrection->SetLineColor(38);

    fModelTotAdapThM2->Draw("SAME");
    fModelTotAdapRaM2->Draw("SAME");
    fModelTotAdapKM2->Draw("SAME");
    fModelTotAdapCoM2->Draw("SAME");
    fModelTotAdapNDBDM2->Draw("SAME");
    fModelTotAdap2NDBDM2->Draw("SAME");
    fModelTotAdapBiM2->Draw("SAME");    
    fModelTotAdapMnM2->Draw("SAME"); 

    // Many Parameters
    TPaveText *pt2 = new TPaveText(0.35,0.77,0.70,0.99,"NB NDC");
    pt2->AddText(Form("Fit Range (M2): %.0f to %.0f keV -- #chi^{2}/NDF: %0.3f", dFitMin, dFitMax, (GetChiSquareAdaptive()/(dFitMaxBinM1+dFitMaxBinM2-dFitMinBinM1-dFitMinBinM2-dNumParameters) )));
    pt2->AddText(Form("Frame Th: %0.2E#pm%0.2E --- TShield Th: %0.2E#pm%0.2E", fParameters[0], fParError[0], fParameters[1], fParError[1] ));
    pt2->AddText(Form("50mK Th: %0.2E#pm%0.2E --- 600mK Th: %0.2E#pm%0.2E", fParameters[13], fParError[13], fParameters[14], fParError[14] ));
    pt2->AddText(Form("IVC Th: %0.2E#pm%0.2E --- OVC Th: %0.2E#pm%0.2E", fParameters[15], fParError[15], fParameters[16], fParError[16] ));
    pt2->AddText(Form("Frame Ra: %0.2E#pm%0.2E --- TShield Ra: %0.2E#pm%0.2E", fParameters[2], fParError[2], fParameters[3], fParError[3] ));
    pt2->AddText(Form("50mK Ra: %0.2E#pm%0.2E --- 600mK Ra: %0.2E#pm%0.2E", fParameters[17], fParError[17], fParameters[18], fParError[18] ));
    pt2->AddText(Form("IVC Ra: %0.2E#pm%0.2E --- OVC Ra: %0.2E#pm%0.2E", fParameters[19], fParError[19], fParameters[20], fParError[20] ));
    pt2->AddText(Form("Close K: %0.2E#pm%0.2E --- Far K: %0.2E#pm%0.2E", fParameters[4], fParError[4], fParameters[5], fParError[5] ));
    pt2->AddText(Form("Close Co: %0.2E#pm%0.2E --- Far Co: %0.2E#pm%0.2E", fParameters[6], fParError[6], fParameters[7], fParError[7] ));
    pt2->AddText(Form("Bi-207: %0.2E#pm%0.2E --- NDBD: %0.2E#pm%0.2E", fParameters[10], fParError[10], fParameters[9], fParError[9] ));
    pt2->AddText(Form("Close Mn-54: %0.2E#pm%0.2E --- Far Mn-54: %0.2E#pm%0.2E", fParameters[11], fParError[11], fParameters[12], fParError[12] ));
    pt2->AddText(Form("2NDBD: %0.2E#pm%0.2E" , fParameters[8], fParError[8] ));
    pt2->Draw();


    TLegend *legfit2 = new TLegend(0.8,0.8,0.97,0.97);
    legfit2->AddEntry(fModelTotAdapM2, "Total PDF", "l");
    legfit2->AddEntry(fModelTotAdapThM2, "Total Th-232", "l");
    legfit2->AddEntry(fModelTotAdapRaM2, "Total Ra-226", "l");
    legfit2->AddEntry(fModelTotAdapKM2, "Total K-40", "l");
    legfit2->AddEntry(fModelTotAdapCoM2, "Total Co-60", "l");
    legfit2->AddEntry(fModelTotAdapNDBDM2, "NDBD", "l");
    legfit2->AddEntry(fModelTotAdap2NDBDM2, "2NDBD", "l");
    legfit2->AddEntry(fModelTotAdapBiM2, "Bi-207", "l");
    legfit2->AddEntry(fModelTotAdapMnM2, "Mn-54", "l");

    legfit2->Draw();


  // Residuals
  TCanvas *cResidual1 = new TCanvas("cResidual1", "cResidual1", 1200, 800);
  hResidualGausM1 = new TH1D("hResidualGausM1", "Residual Distribution (M1)", 100, -50, 50);
  hResidualDistM1 = CalculateResidualsAdaptive(fModelTotAdapM1, fDataHistoM1, hResidualGausM1, dFitMinBinM1, dFitMinBinM2);
  hResidualDistM1->SetLineColor(kBlack);
  hResidualDistM1->SetName("Residuals");
  hResidualDistM1->SetTitle("Fit Residuals (M1)");
  hResidualDistM1->SetMarkerStyle(25);
  hResidualDistM1->GetXaxis()->SetTitle("Energy (keV)");
  // hResidualDistM1->GetXaxis()->SetTitleSize(0.04);
  // hResidualDistM1->GetXaxis()->SetLabelSize(0.05);
  // hResidualDistM1->GetYaxis()->SetLabelSize(0.05);
  // hResidualDistM1->GetYaxis()->SetTitleSize(0.04); 
  hResidualDistM1->GetYaxis()->SetTitle("Residuals (#sigma)");

  hResidualDistM1->GetXaxis()->SetRange(dFitMin/dBinSize-5, dFitMax/dBinSize+5);
  hResidualDistM1->Draw("E");


  TCanvas *cres1 = new TCanvas();
  hResidualGausM1->Draw();

  TCanvas *cResidual2 = new TCanvas("cResidual2", "cResidual2", 1200, 800);
  hResidualGausM2 = new TH1D("hResidualGausM2", "Residual Distribution (M2)", 100, -50, 50);  
  hResidualDistM2 = CalculateResidualsAdaptive(fModelTotAdapM2, fDataHistoM2, hResidualGausM2, dFitMinBinM2, dFitMaxBinM2);
  hResidualDistM2->SetLineColor(kBlack);
  hResidualDistM2->SetMarkerStyle(25);
  hResidualDistM2->SetName("Residuals");
  hResidualDistM2->SetTitle("Fit Residuals (M2)");
  hResidualDistM2->GetXaxis()->SetTitle("Energy (keV)");
  // hResidualDistM2->GetXaxis()->SetTitleSize(0.04);
  // hResidualDistM2->GetXaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetTitleSize(0.04); 
  hResidualDistM2->GetYaxis()->SetTitle("Residuals (#sigma)");

  hResidualDistM2->GetXaxis()->SetRange(dFitMin/dBinSize-5, dFitMax/dBinSize+5);
  hResidualDistM2->Draw("E");

  TCanvas *cres2 = new TCanvas();
  hResidualGausM2->Draw();


  // Output integrals of stuff for limits
  cout << "ROI bin: " << fAdapDataHistoM1->FindBin(2470) << " " << fAdapDataHistoM1->FindBin(2570) << endl;
  cout << "Integral Data in ROI: " << fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt( fAdapDataHistoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total PDF in ROI: " << fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470)) << " +/- " << sqrt( fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total Th PDF in ROI: " << fModelTotAdapThM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt( fModelTotAdapThM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total Ra PDF in ROI: " << fModelTotAdapRaM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt(fModelTotAdapRaM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total Co PDF in ROI: " << fModelTotAdapCoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt(fModelTotAdapCoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total K PDF in ROI: " << fModelTotAdapKM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt(fModelTotAdapKM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total Bi PDF in ROI: " << fModelTotAdapBiM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt(fModelTotAdapBiM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;  
  cout << "Integral Total 2NDBD PDF in ROI: " << fModelTotAdap2NDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt(fModelTotAdap2NDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total 0NDBD PDF in ROI: " << fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt(fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Frame Th PDF in ROI: " << fParameters[0]*fSmearFrameThM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << endl;
  cout << "Integral TShield Th PDF in ROI: " << fParameters[1]*fSmearTShieldThM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << endl;
  cout << "Integral 50mK Th PDF in ROI: " << fParameters[13]*fSmear50mKThM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << endl;
  cout << "Integral 600mK Th PDF in ROI: " << fParameters[14]*fSmear600mKThM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << endl;
  cout << "Integral IVC Th PDF in ROI: " << fParameters[15]*fSmearIVCThM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << endl;
  cout << "Integral OVC Th PDF in ROI: " << fParameters[16]*fSmearOVCThM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << endl;


  return true;

}

