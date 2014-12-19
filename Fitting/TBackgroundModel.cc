#include "TMinuit.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TBackgroundModel.hh"
#include "TRandom3.h"
#include "TPaveText.h"
#include "TAxis.h"
#include "TLine.h"
#include "TMath.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>


using namespace std;

// In an effort to make things more legible I'm trying to clean up the code...

//ClassImp(TBackgroundModel)
  
//first set up a global function that calls your classes method that calculates the quantity to minimise
void myExternal_FCNAdap(int &n, double *grad, double &fval, double x[], int code)
{
  // Required External Wrapper for function to be minimized by Minuit 
 
  // This gets called for each value of the parameters minuit tries
  // here the x array contains the parameters you are trying to fit
  
  // here myClass should inherit from TObject
  TBackgroundModel* Obj = (TBackgroundModel*)gMinuit->GetObjectFit(); 

  // implement a method in your class for setting the parameters and thus update the parameters of your fitter class 
  for(int i = 0; i < 116; i++ )
  {
    Obj->SetParameters(i, x[i]);
  }
  // Implement a method in your class that calculates the quantity you want to minimize, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimise fval
    Obj->UpdateModelAdaptive();
    fval = Obj->GetChiSquareAdaptive();
}


TBackgroundModel::TBackgroundModel(double fFitMin, double fFitMax, int dBinBase)
{
  dNParam = 116; // number of fitting parameters
  dNumCalls = 0;
  dSecToYears = 1./(60*60*24*365);

  // dDataDir =  "/Users/brian/macros/Simulations/Production/";
  dDataDir =  "/Users/brian/macros/CUOREZ/Bkg/";
  dDataIntegral = 0;

  // Bin size (keV) -- base binning is 2 keV
  dBinSize = 2; 
  // Histogram range - from 0 to 10 MeV
  dMinEnergy = 0.;
  dMaxEnergy = 10000.;

  if(fFitMin >= fFitMax)
  {
    cout << "Fit Min >= Fit Max!" << endl;
  }

  minuit = new TMinuit(dNParam);

  // Fitting range
  dFitMin = fFitMin;
  dFitMax = fFitMax;

  dNBins = (dMaxEnergy - dMinEnergy)/ dBinSize;

  // Data Histograms
  fDataHistoTot  = new TH1D("fDataHistoTot",  "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM1   = new TH1D("fDataHistoM1",   "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM2   = new TH1D("fDataHistoM2",   "", dNBins, dMinEnergy, dMaxEnergy);

  // Data variables change when switching between Toy/Real data
  // qtree = new TChain("qredtree");
  // qtree = new TChain("CombiTree");
  // Data cuts 
  qtree = new TChain("qtree");
  base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
  base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
  base_cut = base_cut && "abs(BaselineSlope)<0.1";
  base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";

  // Load data here
  LoadData();


  // Scaling by livetime, don't use with Toy data
  dDataIntegral = fDataHistoM1->Integral(1, dNBins);
  int dDataIntegralTot = qtree->GetEntries();

  cout << "Total Events in background spectrum: " << dDataIntegralTot << endl; 
  cout << "Events in background spectrum (M1): " << fDataHistoM1->Integral(1, 10000/dBinSize) << endl;
  cout << "Events in background spectrum (M2): " << fDataHistoM2->Integral(1, 10000/dBinSize) << endl;

  // Scale by Live-time (ds 2061 - 2100) 14647393.0 seconds
  // fDataHistoM1->Scale(1/((936398+14647393.0) * dSecToYears));
  // fDataHistoM2->Scale(1/((936398+14647393.0) * dSecToYears));  

  // cout << "Normalized Data using Livetime of: " << (936398+14647393.0) * dSecToYears << " years" << endl;


  // Modeling
  gaus = new TF1("gaus","gaus(0)", dMinEnergy, dMaxEnergy);


  // Total model histograms M1
  fModelTotM1      = new TH1D("fModelTotM1",      "Frame",        dNBins, dMinEnergy, dMaxEnergy);  
  fModelTotthM1    = new TH1D("fModelTotthM1",    "Total th232",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotuM1     = new TH1D("fModelTotuM1",     "Total u238",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotkM1     = new TH1D("fModelTotkM1",     "Total k40",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTotcoM1    = new TH1D("fModelTotvoM1",    "Total co60",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotmnM1    = new TH1D("fModelTotmnM1",    "Total mn54",   dNBins, dMinEnergy, dMaxEnergy);

  fModelTotNDBDM1  = new TH1D("fModelTotNDBDM1",  "Total NDBD",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTot2NDBDM1 = new TH1D("fModelTot2NDBDM1", "Total 2NDBD",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotbiM1    = new TH1D("fModelTotbiM1",    "Total bi207",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotbi2M1   = new TH1D("fModelTotbi2M1",   "Total bi210",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotptM1    = new TH1D("fModelTotptM1",    "Total pt190",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotpbM1    = new TH1D("fModelTotpbM1",    "Total pb210",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotcsM1    = new TH1D("fModelTotcsM1",    "Total cs137",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotco2M1   = new TH1D("fModelTotco2M1",   "Total co58",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotteo2M1  = new TH1D("fModelTotteo2M1",  "Total TeO2",   dNBins, dMinEnergy, dMaxEnergy);

  fModelTotSthM1   = new TH1D("fModelTotSthM1",   "Total S th232",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotSuM1    = new TH1D("fModelTotSuM1",    "Total S u238",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotSpbM1   = new TH1D("fModelTotSpbM1",   "Total S pb210",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotSpoM1   = new TH1D("fModelTotSpoM1",   "Total S po210",  dNBins, dMinEnergy, dMaxEnergy);



  // Total model histograms M2
  fModelTotM2      = new TH1D("fModelTotM2",      "Frame",        dNBins, dMinEnergy, dMaxEnergy);  
  fModelTotthM2    = new TH1D("fModelTotthM2",    "Total th232",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotuM2     = new TH1D("fModelTotuM2",     "Total u238",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotkM2     = new TH1D("fModelTotkM2",     "Total k40",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTotcoM2    = new TH1D("fModelTotcoM2",    "Total co60",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotmnM2    = new TH1D("fModelTotmnM2",    "Total mn54",   dNBins, dMinEnergy, dMaxEnergy);

  fModelTotNDBDM2  = new TH1D("fModelTotNDBDM2",  "Total NDBD",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTot2NDBDM2 = new TH1D("fModelTot2NDBDM2", "Total 2NDBD",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotbiM2    = new TH1D("fModelTotbiM2",    "Total bi207",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotbi2M2   = new TH1D("fModelTotbi2M2",   "Total bi210",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotptM2    = new TH1D("fModelTotptM2",    "Total pt190",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotpbM2    = new TH1D("fModelTotpbM2",    "Total pb210",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotcsM2    = new TH1D("fModelTotcsM2",    "Total cs137",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotco2M2   = new TH1D("fModelTotco2M2",   "Total co58",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotteo2M2  = new TH1D("fModelTotteo2M2",  "Total TeO2",   dNBins, dMinEnergy, dMaxEnergy);

  fModelTotSthM2   = new TH1D("fModelTotSthM2",   "Total S th232",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotSuM2    = new TH1D("fModelTotSuM2",    "Total S u238",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotSpbM2   = new TH1D("fModelTotSpbM2",   "Total S pb210",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotSpoM2   = new TH1D("fModelTotSpoM2",   "Total S po210",  dNBins, dMinEnergy, dMaxEnergy);

//////////// Bulk model histograms
  // Crystal M1 and M2
  hTeO20nuM1       = new TH1D("hTeO20nuM1",    "hTeO20nuM1",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO22nuM1       = new TH1D("hTeO22nuM1",    "hTeO22nuM1",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO2co60M1      = new TH1D("hTeO2co60M1",   "hTeO2co60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2k40M1       = new TH1D("hTeO2k40M1",    "hTeO2k40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO2pb210M1     = new TH1D("hTeO2pb210M1",  "hTeO2pb210M1",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2po210M1     = new TH1D("hTeO2po210M1",  "hTeO2po210M1",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2te125M1     = new TH1D("hTeO2te125M1",  "hTeO2te125M1",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2th232M1     = new TH1D("hTeO2th232M1",  "hTeO2th232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hTeO2th228M1     = new TH1D("hTeO2th228M1",  "hTeO2th228M1",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2ra226M1     = new TH1D("hTeO2ra226M1",  "hTeO2ra226M1",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2rn222M1     = new TH1D("hTeO2rn222M1",  "hTeO2rn222M1",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2u238M1      = new TH1D("hTeO2u238M1",   "hTeO2u238M1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2th230M1     = new TH1D("hTeO2th230M1",  "hTeO2th230M1",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2u234M1      = new TH1D("hTeO2u234M1",   "hTeO2u234M1",   dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Spb210M1_01      = new TH1D("hTeO2Spb210M1_01",   "hTeO2Spb210M1_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Spo210M1_001     = new TH1D("hTeO2Spo210M1_001",  "hTeO2Spo210M1_001",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Spo210M1_01      = new TH1D("hTeO2Spo210M1_01",   "hTeO2Spo210M1_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sth232M1_01      = new TH1D("hTeO2Sth232M1_01",   "hTeO2Sth232M1_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Su238M1_01       = new TH1D("hTeO2Su238M1_01",    "hTeO2Su238M1_01",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M1_001    = new TH1D("hTeO2Sxpb210M1_001", "hTeO2Sxpb210M1_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M1_01     = new TH1D("hTeO2Sxpb210M1_01",  "hTeO2Sxpb210M1_01",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M1_1      = new TH1D("hTeO2Sxpb210M1_1",   "hTeO2Sxpb210M1_1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M1_10     = new TH1D("hTeO2Sxpb210M1_10",  "hTeO2Sxpb210M1_10",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpo210M1_001    = new TH1D("hTeO2Sxpo210M1_001", "hTeO2Sxpo210M1_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpo210M1_01     = new TH1D("hTeO2Sxpo210M1_01",  "hTeO2Sxpo210M1_01",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpo210M1_1      = new TH1D("hTeO2Sxpo210M1_1",   "hTeO2Sxpo210M1_1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M1_001    = new TH1D("hTeO2Sxth232M1_001", "hTeO2Sxth232M1_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M1_01     = new TH1D("hTeO2Sxth232M1_01",  "hTeO2Sxth232M1_01",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M1_1      = new TH1D("hTeO2Sxth232M1_1",   "hTeO2Sxth232M1_1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M1_10     = new TH1D("hTeO2Sxth232M1_10",  "hTeO2Sxth232M1_10",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238M1_001     = new TH1D("hTeO2Sxu238M1_001",  "hTeO2Sxu238M1_001",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238M1_01      = new TH1D("hTeO2Sxu238M1_01",   "hTeO2Sxu238M1_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238M1_1       = new TH1D("hTeO2Sxu238M1_1",    "hTeO2Sxu238M1_1",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238M1_10      = new TH1D("hTeO2Sxu238M1_10",   "hTeO2Sxu238M1_10",   dNBins, dMinEnergy, dMaxEnergy);


  hTeO20nuM2       = new TH1D("hTeO20nuM2",    "hTeO20nuM2",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO22nuM2       = new TH1D("hTeO22nuM2",    "hTeO22nuM2",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO2co60M2      = new TH1D("hTeO2co60M2",   "hTeO2co60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2k40M2       = new TH1D("hTeO2k40M2",    "hTeO2k40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO2pb210M2     = new TH1D("hTeO2pb210M2",  "hTeO2pb210M2",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2po210M2     = new TH1D("hTeO2po210M2",  "hTeO2po210M2",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2te125M2     = new TH1D("hTeO2te125M2",  "hTeO2te125M2",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2th232M2     = new TH1D("hTeO2th232M2",  "hTeO2th232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hTeO2th228M2     = new TH1D("hTeO2th228M2",  "hTeO2th228M2",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2ra226M2     = new TH1D("hTeO2ra226M2",  "hTeO2ra226M2",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2rn222M2     = new TH1D("hTeO2rn222M2",  "hTeO2rn222M2",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2u238M2      = new TH1D("hTeO2u238M2",   "hTeO2u238M2",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2th230M2     = new TH1D("hTeO2th230M2",  "hTeO2th230M2",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2u234M2      = new TH1D("hTeO2u234M2",   "hTeO2u234M2",   dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Spb210M2_01      = new TH1D("hTeO2Spb210M2_01",   "hTeO2Spb210M2_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Spo210M2_001     = new TH1D("hTeO2Spo210M2_001",  "hTeO2Spo210M2_001",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Spo210M2_01      = new TH1D("hTeO2Spo210M2_01",   "hTeO2Spo210M2_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sth232M2_01      = new TH1D("hTeO2Sth232M2_01",   "hTeO2Sth232M2_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Su238M2_01       = new TH1D("hTeO2Su238M2_01",    "hTeO2Su238M2_01",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2_001    = new TH1D("hTeO2Sxpb210M2_001", "hTeO2Sxpb210M2_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2_01     = new TH1D("hTeO2Sxpb210M2_01",  "hTeO2Sxpb210M2_01",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2_1      = new TH1D("hTeO2Sxpb210M2_1",   "hTeO2Sxpb210M2_1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2_10     = new TH1D("hTeO2Sxpb210M2_10",  "hTeO2Sxpb210M2_10",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpo210M2_001    = new TH1D("hTeO2Sxpo210M2_001", "hTeO2Sxpo210M2_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpo210M2_01     = new TH1D("hTeO2Sxpo210M2_01",  "hTeO2Sxpo210M2_01",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpo210M2_1      = new TH1D("hTeO2Sxpo210M2_1",   "hTeO2Sxpo210M2_1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M2_001    = new TH1D("hTeO2Sxth232M2_001", "hTeO2Sxth232M2_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M2_01     = new TH1D("hTeO2Sxth232M2_01",  "hTeO2Sxth232M2_01",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M2_1      = new TH1D("hTeO2Sxth232M2_1",   "hTeO2Sxth232M2_1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M2_10     = new TH1D("hTeO2Sxth232M2_10",  "hTeO2Sxth232M2_10",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238M2_001     = new TH1D("hTeO2Sxu238M2_001",  "hTeO2Sxu238M2_001",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238M2_01      = new TH1D("hTeO2Sxu238M2_01",   "hTeO2Sxu238M2_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238M2_1       = new TH1D("hTeO2Sxu238M2_1",    "hTeO2Sxu238M2_1",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238M2_10      = new TH1D("hTeO2Sxu238M2_10",   "hTeO2Sxu238M2_10",   dNBins, dMinEnergy, dMaxEnergy);


  // Frame M1 and M2
  hCuFrameco58M1      = new TH1D("hCuFrameco58M1",   "hCuFrameco58M1",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameco60M1      = new TH1D("hCuFrameco60M1",   "hCuFrameco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFramecs137M1     = new TH1D("hCuFramecs137M1",  "hCuFramecs137M1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFramek40M1       = new TH1D("hCuFramek40M1",    "hCuFramek40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFramemn54M1      = new TH1D("hCuFramemn54M1",   "hCuFramemn54M1",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFramepb210M1     = new TH1D("hCuFramepb210M1",  "hCuFramepb210M1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameth232M1     = new TH1D("hCuFrameth232M1",  "hCuFrameth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hCuFrameu238M1      = new TH1D("hCuFrameu238M1",   "hCuFrameu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hCuFrameSth232M1_1     = new TH1D("hCuFrameSth232M1_1",    "hCuFrameSth232M1_1",     dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSu238M1_1      = new TH1D("hCuFrameSu238M1_1",     "hCuFrameSu238M1_1",      dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxpb210M1_001  = new TH1D("hCuFrameSxpb210M1_001", "hCuFrameSxpb210M1_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxpb210M1_01   = new TH1D("hCuFrameSxpb210M1_01",  "hCuFrameSxpb210M1_01",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxpb210M1_1    = new TH1D("hCuFrameSxpb210M1_1",   "hCuFrameSxpb210M1_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxpb210M1_10   = new TH1D("hCuFrameSxpb210M1_10",  "hCuFrameSxpb210M1_10",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxth232M1_001  = new TH1D("hCuFrameSxth232M1_001", "hCuFrameSxth232M1_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxth232M1_01   = new TH1D("hCuFrameSxth232M1_01",  "hCuFrameSxth232M1_01",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxth232M1_1    = new TH1D("hCuFrameSxth232M1_1",   "hCuFrameSxth232M1_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxth232M1_10   = new TH1D("hCuFrameSxth232M1_10",  "hCuFrameSxth232M1_10",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxu238M1_001   = new TH1D("hCuFrameSxu238M1_001",  "hCuFrameSxu238M1_001",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxu238M1_01    = new TH1D("hCuFrameSxu238M1_01",   "hCuFrameSxu238M1_01",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxu238M1_1     = new TH1D("hCuFrameSxu238M1_1",    "hCuFrameSxu238M1_1",     dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxu238M1_10    = new TH1D("hCuFrameSxu238M1_10",   "hCuFrameSxu238M1_10",    dNBins, dMinEnergy, dMaxEnergy);


  hCuFrameco58M2      = new TH1D("hCuFrameco58M2",   "hCuFrameco58M2",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameco60M2      = new TH1D("hCuFrameco60M2",   "hCuFrameco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFramecs137M2     = new TH1D("hCuFramecs137M2",  "hCuFramecs137M2",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFramek40M2       = new TH1D("hCuFramek40M2",    "hCuFramek40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFramemn54M2      = new TH1D("hCuFramemn54M2",   "hCuFramemn54M2",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFramepb210M2     = new TH1D("hCuFramepb210M2",  "hCuFramepb210M2",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameth232M2     = new TH1D("hCuFrameth232M2",  "hCuFrameth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hCuFrameu238M2      = new TH1D("hCuFrameu238M2",   "hCuFrameu238M2",   dNBins, dMinEnergy, dMaxEnergy);

  hCuFrameSth232M2_1     = new TH1D("hCuFrameSth232M2_1",    "hCuFrameSth232M2_1",     dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSu238M2_1      = new TH1D("hCuFrameSu238M2_1",     "hCuFrameSu238M2_1",      dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxpb210M2_001  = new TH1D("hCuFrameSxpb210M2_001", "hCuFrameSxpb210M2_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxpb210M2_01   = new TH1D("hCuFrameSxpb210M2_01",  "hCuFrameSxpb210M2_01",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxpb210M2_1    = new TH1D("hCuFrameSxpb210M2_1",   "hCuFrameSxpb210M2_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxpb210M2_10   = new TH1D("hCuFrameSxpb210M2_10",  "hCuFrameSxpb210M2_10",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxth232M2_001  = new TH1D("hCuFrameSxth232M2_001", "hCuFrameSxth232M2_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxth232M2_01   = new TH1D("hCuFrameSxth232M2_01",  "hCuFrameSxth232M2_01",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxth232M2_1    = new TH1D("hCuFrameSxth232M2_1",   "hCuFrameSxth232M2_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxth232M2_10   = new TH1D("hCuFrameSxth232M2_10",  "hCuFrameSxth232M2_10",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxu238M2_001   = new TH1D("hCuFrameSxu238M2_001",  "hCuFrameSxu238M2_001",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxu238M2_01    = new TH1D("hCuFrameSxu238M2_01",   "hCuFrameSxu238M2_01",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxu238M2_1     = new TH1D("hCuFrameSxu238M2_1",    "hCuFrameSxu238M2_1",     dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxu238M2_10    = new TH1D("hCuFrameSxu238M2_10",   "hCuFrameSxu238M2_10",    dNBins, dMinEnergy, dMaxEnergy);

  // CuBox (TShield) M1 and M2
  hCuBoxco58M1      = new TH1D("hCuBoxco58M1",   "hCuBoxco58M1",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxco60M1      = new TH1D("hCuBoxco60M1",   "hCuBoxco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxcs137M1     = new TH1D("hCuBoxcs137M1",  "hCuBoxcs137M1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxk40M1       = new TH1D("hCuBoxk40M1",    "hCuBoxk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxmn54M1      = new TH1D("hCuBoxmn54M1",   "hCuBoxmn54M1",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxpb210M1     = new TH1D("hCuBoxpb210M1",  "hCuBoxpb210M1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxth232M1     = new TH1D("hCuBoxth232M1",  "hCuBoxth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hCuBoxu238M1      = new TH1D("hCuBoxu238M1",   "hCuBoxu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hCuBoxSth232M1_1     = new TH1D("hCuBoxSth232M1_1",    "hCuBoxSth232M1_1",     dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSu238M1_1      = new TH1D("hCuBoxSu238M1_1",     "hCuBoxSu238M1_1",      dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxpb210M1_001  = new TH1D("hCuBoxSxpb210M1_001", "hCuBoxSxpb210M1_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxpb210M1_01   = new TH1D("hCuBoxSxpb210M1_01",  "hCuBoxSxpb210M1_01",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxpb210M1_1    = new TH1D("hCuBoxSxpb210M1_1",   "hCuBoxSxpb210M1_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxpb210M1_10   = new TH1D("hCuBoxSxpb210M1_10",  "hCuBoxSxpb210M1_10",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxth232M1_001  = new TH1D("hCuBoxSxth232M1_001", "hCuBoxSxth232M1_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxth232M1_01   = new TH1D("hCuBoxSxth232M1_01",  "hCuBoxSxth232M1_01",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxth232M1_1    = new TH1D("hCuBoxSxth232M1_1",   "hCuBoxSxth232M1_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxth232M1_10   = new TH1D("hCuBoxSxth232M1_10",  "hCuBoxSxth232M1_10",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxu238M1_001   = new TH1D("hCuBoxSxu238M1_001",  "hCuBoxSxu238M1_001",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxu238M1_01    = new TH1D("hCuBoxSxu238M1_01",   "hCuBoxSxu238M1_01",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxu238M1_1     = new TH1D("hCuBoxSxu238M1_1",    "hCuBoxSxu238M1_1",     dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxu238M1_10    = new TH1D("hCuBoxSxu238M1_10",   "hCuBoxSxu238M1_10",    dNBins, dMinEnergy, dMaxEnergy);

  hCuBoxco58M2      = new TH1D("hCuBoxco58M2",   "hCuBoxco58M2",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxco60M2      = new TH1D("hCuBoxco60M2",   "hCuBoxco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxcs137M2     = new TH1D("hCuBoxcs137M2",  "hCuBoxcs137M2",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxk40M2       = new TH1D("hCuBoxk40M2",    "hCuBoxk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxmn54M2      = new TH1D("hCuBoxmn54M2",   "hCuBoxmn54M2",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxpb210M2     = new TH1D("hCuBoxpb210M2",  "hCuBoxpb210M2",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxth232M2     = new TH1D("hCuBoxth232M2",  "hCuBoxth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hCuBoxu238M2      = new TH1D("hCuBoxu238M2",   "hCuBoxu238M2",   dNBins, dMinEnergy, dMaxEnergy);

  hCuBoxSth232M2_1     = new TH1D("hCuBoxSth232M2_1",    "hCuBoxSth232M2_1",     dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSu238M2_1      = new TH1D("hCuBoxSu238M2_1",     "hCuBoxSu238M2_1",      dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxpb210M2_001  = new TH1D("hCuBoxSxpb210M2_001", "hCuBoxSxpb210M2_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxpb210M2_01   = new TH1D("hCuBoxSxpb210M2_01",  "hCuBoxSxpb210M2_01",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxpb210M2_1    = new TH1D("hCuBoxSxpb210M2_1",   "hCuBoxSxpb210M2_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxpb210M2_10   = new TH1D("hCuBoxSxpb210M2_10",  "hCuBoxSxpb210M2_10",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxth232M2_001  = new TH1D("hCuBoxSxth232M2_001", "hCuBoxSxth232M2_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxth232M2_01   = new TH1D("hCuBoxSxth232M2_01",  "hCuBoxSxth232M2_01",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxth232M2_1    = new TH1D("hCuBoxSxth232M2_1",   "hCuBoxSxth232M2_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxth232M2_10   = new TH1D("hCuBoxSxth232M2_10",  "hCuBoxSxth232M2_10",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxu238M2_001   = new TH1D("hCuBoxSxu238M2_001",  "hCuBoxSxu238M2_001",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxu238M2_01    = new TH1D("hCuBoxSxu238M2_01",   "hCuBoxSxu238M2_01",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxu238M2_1     = new TH1D("hCuBoxSxu238M2_1",    "hCuBoxSxu238M2_1",     dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxu238M2_10    = new TH1D("hCuBoxSxu238M2_10",   "hCuBoxSxu238M2_10",    dNBins, dMinEnergy, dMaxEnergy);

  // 50mK M1 and M2
  h50mKco58M1      = new TH1D("h50mKco58M1",   "h50mKco58M1",   dNBins, dMinEnergy, dMaxEnergy);
  h50mKco60M1      = new TH1D("h50mKco60M1",   "h50mKco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  h50mKcs137M1     = new TH1D("h50mKcs137M1",  "h50mKcs137M1",  dNBins, dMinEnergy, dMaxEnergy);
  h50mKk40M1       = new TH1D("h50mKk40M1",    "h50mKk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  h50mKmn54M1      = new TH1D("h50mKmn54M1",   "h50mKmn54M1",   dNBins, dMinEnergy, dMaxEnergy);
  h50mKpb210M1     = new TH1D("h50mKpb210M1",  "h50mKpb210M1",  dNBins, dMinEnergy, dMaxEnergy);
  h50mKth232M1     = new TH1D("h50mKth232M1",  "h50mKth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  h50mKu238M1      = new TH1D("h50mKu238M1",   "h50mKu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  h50mKco58M2      = new TH1D("h50mKco58M2",   "h50mKco58M2",   dNBins, dMinEnergy, dMaxEnergy);
  h50mKco60M2      = new TH1D("h50mKco60M2",   "h50mKco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  h50mKcs137M2     = new TH1D("h50mKcs137M2",  "h50mKcs137M2",  dNBins, dMinEnergy, dMaxEnergy);
  h50mKk40M2       = new TH1D("h50mKk40M2",    "h50mKk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  h50mKmn54M2      = new TH1D("h50mKmn54M2",   "h50mKmn54M2",   dNBins, dMinEnergy, dMaxEnergy);
  h50mKpb210M2     = new TH1D("h50mKpb210M2",  "h50mKpb210M2",  dNBins, dMinEnergy, dMaxEnergy);
  h50mKth232M2     = new TH1D("h50mKth232M2",  "h50mKth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  h50mKu238M2      = new TH1D("h50mKu238M2",   "h50mKu238M2",   dNBins, dMinEnergy, dMaxEnergy);

  // 600mK M1 and M2
  h600mKco60M1      = new TH1D("h600mKco60M1",   "h600mKco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  h600mKk40M1       = new TH1D("h600mKk40M1",    "h600mKk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  h600mKth232M1     = new TH1D("h600mKth232M1",  "h600mKth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  h600mKu238M1      = new TH1D("h600mKu238M1",   "h600mKu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  h600mKco60M2      = new TH1D("h600mKco60M2",   "h600mKco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  h600mKk40M2       = new TH1D("h600mKk40M2",    "h600mKk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  h600mKth232M2     = new TH1D("h600mKth232M2",  "h600mKth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  h600mKu238M2      = new TH1D("h600mKu238M2",   "h600mKu238M2",   dNBins, dMinEnergy, dMaxEnergy);  

  // Roman Lead M1 and M2
  hPbRombi207M1     = new TH1D("hPbRombi207M1",  "hPbRombi207M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomco60M1      = new TH1D("hPbRomco60M1",   "hPbRomco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hPbRomcs137M1     = new TH1D("hPbRomcs137M1",  "hPbRomcs137M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomk40M1       = new TH1D("hPbRomk40M1",    "hPbRomk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hPbRompb210M1     = new TH1D("hPbRompb210M1",  "hPbRompb210M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomth232M1     = new TH1D("hPbRomth232M1",  "hPbRomth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomu238M1      = new TH1D("hPbRomu238M1",   "hPbRomu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hPbRombi207M2     = new TH1D("hPbRombi207M2",  "hPbRombi207M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomco60M2      = new TH1D("hPbRomco60M2",   "hPbRomco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hPbRomcs137M2     = new TH1D("hPbRomcs137M2",  "hPbRomcs137M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomk40M2       = new TH1D("hPbRomk40M2",    "hPbRomk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hPbRompb210M2     = new TH1D("hPbRompb210M2",  "hPbRompb210M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomth232M2     = new TH1D("hPbRomth232M2",  "hPbRomth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomu238M2      = new TH1D("hPbRomu238M2",   "hPbRomu238M2",   dNBins, dMinEnergy, dMaxEnergy);

  // Main bath M1 and M2
  hMBco60M1      = new TH1D("hMBco60M1",   "hMBco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hMBk40M1       = new TH1D("hMBk40M1",    "hMBk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hMBth232M1     = new TH1D("hMBth232M1",  "hMBth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hMBu238M1      = new TH1D("hMBu238M1",   "hMBu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hMBco60M2      = new TH1D("hMBco60M2",   "hMBco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hMBk40M2       = new TH1D("hMBk40M2",    "hMBk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hMBth232M2     = new TH1D("hMBth232M2",  "hMBth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hMBu238M2      = new TH1D("hMBu238M2",   "hMBu238M2",   dNBins, dMinEnergy, dMaxEnergy);  

  // Super Insulation M1 and M2
  hSIk40M1       = new TH1D("hSIk40M1",    "hSIk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hSIth232M1     = new TH1D("hSIth232M1",  "hSIth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hSIu238M1      = new TH1D("hSIu238M1",   "hSIu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hSIk40M2       = new TH1D("hSIk40M2",    "hSIk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hSIth232M2     = new TH1D("hSIth232M2",  "hSIth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hSIu238M2      = new TH1D("hSIu238M2",   "hSIu238M2",   dNBins, dMinEnergy, dMaxEnergy);


  // IVC M1 and M2
  hIVCco60M1      = new TH1D("hIVCco60M1",   "hIVCco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hIVCk40M1       = new TH1D("hIVCk40M1",    "hIVCk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hIVCth232M1     = new TH1D("hIVCth232M1",  "hIVCth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hIVCu238M1      = new TH1D("hIVCu238M1",   "hIVCu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hIVCco60M2      = new TH1D("hIVCco60M2",   "hIVCco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hIVCk40M2       = new TH1D("hIVCk40M2",    "hIVCk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hIVCth232M2     = new TH1D("hIVCth232M2",  "hIVCth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hIVCu238M2      = new TH1D("hIVCu238M2",   "hIVCu238M2",   dNBins, dMinEnergy, dMaxEnergy);  

  // OVC M1 and M2
  hOVCco60M1      = new TH1D("hOVCco60M1",   "hOVCco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hOVCk40M1       = new TH1D("hOVCk40M1",    "hOVCk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hOVCth232M1     = new TH1D("hOVCth232M1",  "hOVCth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hOVCu238M1      = new TH1D("hOVCu238M1",   "hOVCu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hOVCco60M2      = new TH1D("hOVCco60M2",   "hOVCco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hOVCk40M2       = new TH1D("hOVCk40M2",    "hOVCk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hOVCth232M2     = new TH1D("hOVCth232M2",  "hOVCth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hOVCu238M2      = new TH1D("hOVCu238M2",   "hOVCu238M2",   dNBins, dMinEnergy, dMaxEnergy);  

/////// Adaptive binning
 // Calculates adaptive binning vectors
  dAdaptiveVectorM1 = AdaptiveBinning(fDataHistoM1, dBinBase);
  dAdaptiveBinsM1 = dAdaptiveVectorM1.size() - 1;
  dAdaptiveArrayM1 = &dAdaptiveVectorM1[0];
  dAdaptiveVectorM2 = AdaptiveBinning(fDataHistoM2, dBinBase);
  dAdaptiveBinsM2 = dAdaptiveVectorM2.size() - 1;
  dAdaptiveArrayM2 = &dAdaptiveVectorM2[0];

  // Adaptive binning data
  fAdapDataHistoM1   = new TH1D("fAdapDataHistoM1",   "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapDataHistoM2   = new TH1D("fAdapDataHistoM2",   "", dAdaptiveBinsM2, dAdaptiveArrayM2);
  
  fDataHistoM1->Rebin(dAdaptiveBinsM1, "hnewM1", dAdaptiveArrayM1);
  fDataHistoM2->Rebin(dAdaptiveBinsM2, "hnewM2", dAdaptiveArrayM2);

  for(int i = 1; i <= dAdaptiveBinsM1; i++)
  {
    fAdapDataHistoM1->SetBinContent(i, dBinSize * hnewM1->GetBinContent(i)/hnewM1->GetBinWidth(i));
  }

  for(int i = 1; i <= dAdaptiveBinsM2; i++)
  {
    fAdapDataHistoM2->SetBinContent(i, dBinSize * hnewM2->GetBinContent(i)/hnewM2->GetBinWidth(i));
  }

  dFitMinBinM1 = fAdapDataHistoM1->FindBin(dFitMin);
  dFitMinBinM2 = fAdapDataHistoM2->FindBin(dFitMin);
  dFitMaxBinM1 = fAdapDataHistoM1->FindBin(dFitMax);
  dFitMaxBinM2 = fAdapDataHistoM2->FindBin(dFitMax);

//////////////// Adaptive binned histograms
  // Total Adaptive binning histograms M1
  fModelTotAdapM1      = new TH1D("fModelTotAdapM1",      "Total PDF M1", dAdaptiveBinsM1, dAdaptiveArrayM1);  
  fModelTotAdapthM1    = new TH1D("fModelTotAdapthM1",    "Total th232",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapuM1     = new TH1D("fModelTotAdapuM1",     "Total u238",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapkM1     = new TH1D("fModelTotAdapkM1",     "Total k40",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapcoM1    = new TH1D("fModelTotAdapcoM1",    "Total co60",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapmnM1    = new TH1D("fModelTotAdapmnM1",    "Total mn54",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  fModelTotAdapNDBDM1  = new TH1D("fModelTotAdapNDBDM1",  "Total NDBD",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdap2NDBDM1 = new TH1D("fModelTotAdap2NDBDM1", "Total 2NDBD",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapbiM1    = new TH1D("fModelTotAdapbiM1",    "Total bi207",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapbi2M1   = new TH1D("fModelTotAdapbi2M1",   "Total bi210",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapptM1    = new TH1D("fModelTotAdapptM1",    "Total pt190",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdappbM1    = new TH1D("fModelTotAdappbM1",    "Total pb210",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapcsM1    = new TH1D("fModelTotAdapcsM1",    "Total cs137",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapco2M1   = new TH1D("fModelTotAdapco2M1",   "Total co58",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapteo2M1  = new TH1D("fModelTotAdapteo2M1",  "Total TeO2",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  fModelTotAdapSthM1   = new TH1D("fModelTotAdapSthM1",   "Total S th232",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapSuM1    = new TH1D("fModelTotAdapSuM1",    "Total S u238",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapSpbM1   = new TH1D("fModelTotAdapSpbM1",   "Total S pb210",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapSpoM1   = new TH1D("fModelTotAdapSpoM1",   "Total S po210",  dAdaptiveBinsM1, dAdaptiveArrayM1);

  fModelTotAdapAlphaM1      = new TH1D("fModelTotAdapAlphaM1",   "Total Alphas",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapAlphaHighM1  = new TH1D("fModelTotAdapAlphaHighM1",   "Total Alphas High",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapAlphaLowM1  = new TH1D("fModelTotAdapAlphaLowM1",   "Total Alphas Low",  dAdaptiveBinsM1, dAdaptiveArrayM1);


  // Total Adaptive binning histograms M2
  fModelTotAdapM2      = new TH1D("fModelTotAdapM2",      "Total PDF M2", dAdaptiveBinsM2, dAdaptiveArrayM2);  
  fModelTotAdapthM2    = new TH1D("fModelTotAdapthM2",    "Total th232",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapuM2     = new TH1D("fModelTotAdapuM2",     "Total u238",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapkM2     = new TH1D("fModelTotAdapkM2",     "Total k40",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapcoM2    = new TH1D("fModelTotAdapcoM2",    "Total co60",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapmnM2    = new TH1D("fModelTotAdapmnM2",    "Total mn54",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  fModelTotAdapNDBDM2  = new TH1D("fModelTotAdapNDBDM2",  "Total NDBD",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdap2NDBDM2 = new TH1D("fModelTotAdap2NDBDM2", "Total 2NDBD",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapbiM2    = new TH1D("fModelTotAdapbiM2",    "Total bi207",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapbi2M2   = new TH1D("fModelTotAdapbi2M2",   "Total bi210",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapptM2    = new TH1D("fModelTotAdapotM2",    "Total pt190",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdappbM2    = new TH1D("fModelTotAdappbM2",    "Total pb210",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapcsM2    = new TH1D("fModelTotAdapcsM2",    "Total cs137",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapco2M2   = new TH1D("fModelTotAdapco2M2",   "Total co58",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapteo2M2  = new TH1D("fModelTotAdapteo2M2",  "Total TeO2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  fModelTotAdapSthM2   = new TH1D("fModelTotAdapSthM2",   "Total S th232",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapSuM2    = new TH1D("fModelTotAdapSuM2",    "Total S u238",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapSpbM2   = new TH1D("fModelTotAdapSpbM2",   "Total S pb210",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapSpoM2   = new TH1D("fModelTotAdapSpoM2",   "Total S po210",  dAdaptiveBinsM2, dAdaptiveArrayM2);

  fModelTotAdapAlphaM2      = new TH1D("fModelTotAdapAlphaM2",   "Total Alphas",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapAlphaHighM2  = new TH1D("fModelTotAdapAlphaHighM2",   "Total Alphas High",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapAlphaLowM2  = new TH1D("fModelTotAdapAlphaLowM2",   "Total Alphas Low",  dAdaptiveBinsM2, dAdaptiveArrayM2);

  // Crystal M1 and M2
  hAdapTeO20nuM1       = new TH1D("hAdapTeO20nuM1",    "hAdapTeO20nuM1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO22nuM1       = new TH1D("hAdapTeO22nuM1",    "hAdapTeO22nuM1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2co60M1      = new TH1D("hAdapTeO2co60M1",   "hAdapTeO2co60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2k40M1       = new TH1D("hAdapTeO2k40M1",    "hAdapTeO2k40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2pb210M1     = new TH1D("hAdapTeO2pb210M1",  "hAdapTeO2pb210M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2po210M1     = new TH1D("hAdapTeO2po210M1",  "hAdapTeO2po210M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2te125M1     = new TH1D("hAdapTeO2te125M1",  "hAdapTeO2te125M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2th232M1     = new TH1D("hAdapTeO2th232M1",  "hAdapTeO2th232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapTeO2th228M1     = new TH1D("hAdapTeO2th228M1",  "hAdapTeO2th228M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2ra226M1     = new TH1D("hAdapTeO2ra226M1",  "hAdapTeO2ra226M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2rn222M1     = new TH1D("hAdapTeO2rn222M1",  "hAdapTeO2rn222M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2u238M1      = new TH1D("hAdapTeO2u238M1",   "hAdapTeO2u238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2th230M1     = new TH1D("hAdapTeO2th230M1",  "hAdapTeO2th230M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2u234M1      = new TH1D("hAdapTeO2u234M1",   "hAdapTeO2u234M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapTeO2Spb210M1_01      = new TH1D("hAdapTeO2Spb210M1_01",   "hAdapTeO2Spb210M1_01",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Spo210M1_001     = new TH1D("hAdapTeO2Spo210M1_001",  "hAdapTeO2Spo210M1_001",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Spo210M1_01      = new TH1D("hAdapTeO2Spo210M1_01",   "hAdapTeO2Spo210M1_01",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sth232M1_01      = new TH1D("hAdapTeO2Sth232M1_01",   "hAdapTeO2Sth232M1_01",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Su238M1_01       = new TH1D("hAdapTeO2Su238M1_01",    "hAdapTeO2Su238M1_01",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpb210M1_001    = new TH1D("hAdapTeO2Sxpb210M1_001", "hAdapTeO2Sxpb210M1_001", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpb210M1_01     = new TH1D("hAdapTeO2Sxpb210M1_01",  "hAdapTeO2Sxpb210M1_01",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpb210M1_1      = new TH1D("hAdapTeO2Sxpb210M1_1",   "hAdapTeO2Sxpb210M1_1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpb210M1_10     = new TH1D("hAdapTeO2Sxpb210M1_10",  "hAdapTeO2Sxpb210M1_10",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpo210M1_001    = new TH1D("hAdapTeO2Sxpo210M1_001", "hAdapTeO2Sxpo210M1_001", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpo210M1_01     = new TH1D("hAdapTeO2Sxpo210M1_01",  "hAdapTeO2Sxpo210M1_01",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpo210M1_1      = new TH1D("hAdapTeO2Sxpo210M1_1",   "hAdapTeO2Sxpo210M1_1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxth232M1_001    = new TH1D("hAdapTeO2Sxth232M1_001", "hAdapTeO2Sxth232M1_001", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxth232M1_01     = new TH1D("hAdapTeO2Sxth232M1_01",  "hAdapTeO2Sxth232M1_01",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxth232M1_1      = new TH1D("hAdapTeO2Sxth232M1_1",   "hAdapTeO2Sxth232M1_1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxth232M1_10     = new TH1D("hAdapTeO2Sxth232M1_10",  "hAdapTeO2Sxth232M1_10",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxu238M1_001     = new TH1D("hAdapTeO2Sxu238M1_001",  "hAdapTeO2Sxu238M1_001",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxu238M1_01      = new TH1D("hAdapTeO2Sxu238M1_01",   "hAdapTeO2Sxu238M1_01",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxu238M1_1       = new TH1D("hAdapTeO2Sxu238M1_1",    "hAdapTeO2Sxu238M1_1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxu238M1_10      = new TH1D("hAdapTeO2Sxu238M1_10",   "hAdapTeO2Sxu238M1_10",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapTeO20nuM2       = new TH1D("hAdapTeO20nuM2",    "hAdapTeO20nuM2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO22nuM2       = new TH1D("hAdapTeO22nuM2",    "hAdapTeO22nuM2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2co60M2      = new TH1D("hAdapTeO2co60M2",   "hAdapTeO2co60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2k40M2       = new TH1D("hAdapTeO2k40M2",    "hAdapTeO2k40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2pb210M2     = new TH1D("hAdapTeO2pb210M2",  "hAdapTeO2pb210M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2po210M2     = new TH1D("hAdapTeO2po210M2",  "hAdapTeO2po210M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2te125M2     = new TH1D("hAdapTeO2te125M2",  "hAdapTeO2te125M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2th232M2     = new TH1D("hAdapTeO2th232M2",  "hAdapTeO2th232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapTeO2th228M2     = new TH1D("hAdapTeO2th228M2",  "hAdapTeO2th228M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2ra226M2     = new TH1D("hAdapTeO2ra226M2",  "hAdapTeO2ra226M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2rn222M2     = new TH1D("hAdapTeO2rn222M2",  "hAdapTeO2rn222M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2u238M2      = new TH1D("hAdapTeO2u238M2",   "hAdapTeO2u238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2th230M2     = new TH1D("hAdapTeO2th230M2",  "hAdapTeO2th230M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2u234M2      = new TH1D("hAdapTeO2u234M2",   "hAdapTeO2u234M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapTeO2Spb210M2_01      = new TH1D("hAdapTeO2Spb210M2_01",   "hAdapTeO2Spb210M2_01",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Spo210M2_001     = new TH1D("hAdapTeO2Spo210M2_001",  "hAdapTeO2Spo210M2_001",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Spo210M2_01      = new TH1D("hAdapTeO2Spo210M2_01",   "hAdapTeO2Spo210M2_01",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sth232M2_01      = new TH1D("hAdapTeO2Sth232M2_01",   "hAdapTeO2Sth232M2_01",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Su238M2_01       = new TH1D("hAdapTeO2Su238M2_01",    "hAdapTeO2Su238M2_01",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpb210M2_001    = new TH1D("hAdapTeO2Sxpb210M2_001", "hAdapTeO2Sxpb210M2_001", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpb210M2_01     = new TH1D("hAdapTeO2Sxpb210M2_01",  "hAdapTeO2Sxpb210M2_01",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpb210M2_1      = new TH1D("hAdapTeO2Sxpb210M2_1",   "hAdapTeO2Sxpb210M2_1",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpb210M2_10     = new TH1D("hAdapTeO2Sxpb210M2_10",  "hAdapTeO2Sxpb210M2_10",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpo210M2_001    = new TH1D("hAdapTeO2Sxpo210M2_001", "hAdapTeO2Sxpo210M2_001", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpo210M2_01     = new TH1D("hAdapTeO2Sxpo210M2_01",  "hAdapTeO2Sxpo210M2_01",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpo210M2_1      = new TH1D("hAdapTeO2Sxpo210M2_1",   "hAdapTeO2Sxpo210M2_1",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxth232M2_001    = new TH1D("hAdapTeO2Sxth232M2_001", "hAdapTeO2Sxth232M2_001", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxth232M2_01     = new TH1D("hAdapTeO2Sxth232M2_01",  "hAdapTeO2Sxth232M2_01",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxth232M2_1      = new TH1D("hAdapTeO2Sxth232M2_1",   "hAdapTeO2Sxth232M2_1",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxth232M2_10     = new TH1D("hAdapTeO2Sxth232M2_10",  "hAdapTeO2Sxth232M2_10",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxu238M2_001     = new TH1D("hAdapTeO2Sxu238M2_001",  "hAdapTeO2Sxu238M2_001",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxu238M2_01      = new TH1D("hAdapTeO2Sxu238M2_01",   "hAdapTeO2Sxu238M2_01",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxu238M2_1       = new TH1D("hAdapTeO2Sxu238M2_1",    "hAdapTeO2Sxu238M2_1",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxu238M2_10      = new TH1D("hAdapTeO2Sxu238M2_10",   "hAdapTeO2Sxu238M2_10",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  // Frame M1 and M2
  hAdapCuFrameco58M1      = new TH1D("hAdapCuFrameco58M1",   "hAdapCuFrameco58M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameco60M1      = new TH1D("hAdapCuFrameco60M1",   "hAdapCuFrameco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFramecs137M1     = new TH1D("hAdapCuFramecs137M1",  "hAdapCuFramecs137M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFramek40M1       = new TH1D("hAdapCuFramek40M1",    "hAdapCuFramek40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFramemn54M1      = new TH1D("hAdapCuFramemn54M1",   "hAdapCuFramemn54M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFramepb210M1     = new TH1D("hAdapCuFramepb210M1",  "hAdapCuFramepb210M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameth232M1     = new TH1D("hAdapCuFrameth232M1",  "hAdapCuFrameth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapCuFrameu238M1      = new TH1D("hAdapCuFrameu238M1",   "hAdapCuFrameu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuFrameSth232M1_1     = new TH1D("hAdapCuFrameSth232M1_1",    "hAdapCuFrameSth232M1_1",     dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSu238M1_1      = new TH1D("hAdapCuFrameSu238M1_1",     "hAdapCuFrameSu238M1_1",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxpb210M1_001  = new TH1D("hAdapCuFrameSxpb210M1_001", "hAdapCuFrameSxpb210M1_001",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxpb210M1_01   = new TH1D("hAdapCuFrameSxpb210M1_01",  "hAdapCuFrameSxpb210M1_01",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxpb210M1_1    = new TH1D("hAdapCuFrameSxpb210M1_1",   "hAdapCuFrameSxpb210M1_1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxpb210M1_10   = new TH1D("hAdapCuFrameSxpb210M1_10",  "hAdapCuFrameSxpb210M1_10",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxth232M1_001  = new TH1D("hAdapCuFrameSxth232M1_001", "hAdapCuFrameSxth232M1_001",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxth232M1_01   = new TH1D("hAdapCuFrameSxth232M1_01",  "hAdapCuFrameSxth232M1_01",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxth232M1_1    = new TH1D("hAdapCuFrameSxth232M1_1",   "hAdapCuFrameSxth232M1_1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxth232M1_10   = new TH1D("hAdapCuFrameSxth232M1_10",  "hAdapCuFrameSxth232M1_10",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxu238M1_001   = new TH1D("hAdapCuFrameSxu238M1_001",  "hAdapCuFrameSxu238M1_001",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxu238M1_01    = new TH1D("hAdapCuFrameSxu238M1_01",   "hAdapCuFrameSxu238M1_01",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxu238M1_1     = new TH1D("hAdapCuFrameSxu238M1_1",    "hAdapCuFrameSxu238M1_1",     dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxu238M1_10    = new TH1D("hAdapCuFrameSxu238M1_10",   "hAdapCuFrameSxu238M1_10",    dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuFrameco58M2      = new TH1D("hAdapCuFrameco58M2",   "hAdapCuFrameco58M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameco60M2      = new TH1D("hAdapCuFrameco60M2",   "hAdapCuFrameco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFramecs137M2     = new TH1D("hAdapCuFramecs137M2",  "hAdapCuFramecs137M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFramek40M2       = new TH1D("hAdapCuFramek40M2",    "hAdapCuFramek40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFramemn54M2      = new TH1D("hAdapCuFramemn54M2",   "hAdapCuFramemn54M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFramepb210M2     = new TH1D("hAdapCuFramepb210M2",  "hAdapCuFramepb210M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameth232M2     = new TH1D("hAdapCuFrameth232M2",  "hAdapCuFrameth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapCuFrameu238M2      = new TH1D("hAdapCuFrameu238M2",   "hAdapCuFrameu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapCuFrameSth232M2_1     = new TH1D("hAdapCuFrameSth232M2_1",    "hAdapCuFrameSth232M2_1",     dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSu238M2_1      = new TH1D("hAdapCuFrameSu238M2_1",     "hAdapCuFrameSu238M2_1",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxpb210M2_001  = new TH1D("hAdapCuFrameSxpb210M2_001", "hAdapCuFrameSxpb210M2_001",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxpb210M2_01   = new TH1D("hAdapCuFrameSxpb210M2_01",  "hAdapCuFrameSxpb210M2_01",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxpb210M2_1    = new TH1D("hAdapCuFrameSxpb210M2_1",   "hAdapCuFrameSxpb210M2_1",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxpb210M2_10   = new TH1D("hAdapCuFrameSxpb210M2_10",  "hAdapCuFrameSxpb210M2_10",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxth232M2_001  = new TH1D("hAdapCuFrameSxth232M2_001", "hAdapCuFrameSxth232M2_001",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxth232M2_01   = new TH1D("hAdapCuFrameSxth232M2_01",  "hAdapCuFrameSxth232M2_01",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxth232M2_1    = new TH1D("hAdapCuFrameSxth232M2_1",   "hAdapCuFrameSxth232M2_1",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxth232M2_10   = new TH1D("hAdapCuFrameSxth232M2_10",  "hAdapCuFrameSxth232M2_10",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxu238M2_001   = new TH1D("hAdapCuFrameSxu238M2_001",  "hAdapCuFrameSxu238M2_001",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxu238M2_01    = new TH1D("hAdapCuFrameSxu238M2_01",   "hAdapCuFrameSxu238M2_01",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxu238M2_1     = new TH1D("hAdapCuFrameSxu238M2_1",    "hAdapCuFrameSxu238M2_1",     dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxu238M2_10    = new TH1D("hAdapCuFrameSxu238M2_10",   "hAdapCuFrameSxu238M2_10",    dAdaptiveBinsM2, dAdaptiveArrayM2);

  // CuBox (TShield) M1 and M2
  hAdapCuBoxco58M1      = new TH1D("hAdapCuBoxco58M1",   "hAdapCuBoxco58M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxco60M1      = new TH1D("hAdapCuBoxco60M1",   "hAdapCuBoxco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxcs137M1     = new TH1D("hAdapCuBoxcs137M1",  "hAdapCuBoxcs137M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxk40M1       = new TH1D("hAdapCuBoxk40M1",    "hAdapCuBoxk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxmn54M1      = new TH1D("hAdapCuBoxmn54M1",   "hAdapCuBoxmn54M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxpb210M1     = new TH1D("hAdapCuBoxpb210M1",  "hAdapCuBoxpb210M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxth232M1     = new TH1D("hAdapCuBoxth232M1",  "hAdapCuBoxth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapCuBoxu238M1      = new TH1D("hAdapCuBoxu238M1",   "hAdapCuBoxu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuBoxSth232M1_1     = new TH1D("hAdapCuBoxSth232M1_1",    "hAdapCuBoxSth232M1_1",     dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSu238M1_1      = new TH1D("hAdapCuBoxSu238M1_1",     "hAdapCuBoxSu238M1_1",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxpb210M1_001  = new TH1D("hAdapCuBoxSxpb210M1_001", "hAdapCuBoxSxpb210M1_001",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxpb210M1_01   = new TH1D("hAdapCuBoxSxpb210M1_01",  "hAdapCuBoxSxpb210M1_01",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxpb210M1_1    = new TH1D("hAdapCuBoxSxpb210M1_1",   "hAdapCuBoxSxpb210M1_1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxpb210M1_10   = new TH1D("hAdapCuBoxSxpb210M1_10",  "hAdapCuBoxSxpb210M1_10",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxth232M1_001  = new TH1D("hAdapCuBoxSxth232M1_001", "hAdapCuBoxSxth232M1_001",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxth232M1_01   = new TH1D("hAdapCuBoxSxth232M1_01",  "hAdapCuBoxSxth232M1_01",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxth232M1_1    = new TH1D("hAdapCuBoxSxth232M1_1",   "hAdapCuBoxSxth232M1_1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxth232M1_10   = new TH1D("hAdapCuBoxSxth232M1_10",  "hAdapCuBoxSxth232M1_10",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxu238M1_001   = new TH1D("hAdapCuBoxSxu238M1_001",  "hAdapCuBoxSxu238M1_001",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxu238M1_01    = new TH1D("hAdapCuBoxSxu238M1_01",   "hAdapCuBoxSxu238M1_01",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxu238M1_1     = new TH1D("hAdapCuBoxSxu238M1_1",    "hAdapCuBoxSxu238M1_1",     dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxu238M1_10    = new TH1D("hAdapCuBoxSxu238M1_10",   "hAdapCuBoxSxu238M1_10",    dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuBoxco58M2      = new TH1D("hAdapCuBoxco58M2",   "hAdapCuBoxco58M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxco60M2      = new TH1D("hAdapCuBoxco60M2",   "hAdapCuBoxco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxcs137M2     = new TH1D("hAdapCuBoxcs137M2",  "hAdapCuBoxcs137M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxk40M2       = new TH1D("hAdapCuBoxk40M2",    "hAdapCuBoxk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxmn54M2      = new TH1D("hAdapCuBoxmn54M2",   "hAdapCuBoxmn54M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxpb210M2     = new TH1D("hAdapCuBoxpb210M2",  "hAdapCuBoxpb210M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxth232M2     = new TH1D("hAdapCuBoxth232M2",  "hAdapCuBoxth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapCuBoxu238M2      = new TH1D("hAdapCuBoxu238M2",   "hAdapCuBoxu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapCuBoxSth232M2_1     = new TH1D("hAdapCuBoxSth232M2_1",    "hAdapCuBoxSth232M2_1",     dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSu238M2_1      = new TH1D("hAdapCuBoxSu238M2_1",     "hAdapCuBoxSu238M2_1",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxpb210M2_001  = new TH1D("hAdapCuBoxSxpb210M2_001", "hAdapCuBoxSxpb210M2_001",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxpb210M2_01   = new TH1D("hAdapCuBoxSxpb210M2_01",  "hAdapCuBoxSxpb210M2_01",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxpb210M2_1    = new TH1D("hAdapCuBoxSxpb210M2_1",   "hAdapCuBoxSxpb210M2_1",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxpb210M2_10   = new TH1D("hAdapCuBoxSxpb210M2_10",  "hAdapCuBoxSxpb210M2_10",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxth232M2_001  = new TH1D("hAdapCuBoxSxth232M2_001", "hAdapCuBoxSxth232M2_001",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxth232M2_01   = new TH1D("hAdapCuBoxSxth232M2_01",  "hAdapCuBoxSxth232M2_01",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxth232M2_1    = new TH1D("hAdapCuBoxSxth232M2_1",   "hAdapCuBoxSxth232M2_1",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxth232M2_10   = new TH1D("hAdapCuBoxSxth232M2_10",  "hAdapCuBoxSxth232M2_10",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxu238M2_001   = new TH1D("hAdapCuBoxSxu238M2_001",  "hAdapCuBoxSxu238M2_001",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxu238M2_01    = new TH1D("hAdapCuBoxSxu238M2_01",   "hAdapCuBoxSxu238M2_01",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxu238M2_1     = new TH1D("hAdapCuBoxSxu238M2_1",    "hAdapCuBoxSxu238M2_1",     dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxu238M2_10    = new TH1D("hAdapCuBoxSxu238M2_10",   "hAdapCuBoxSxu238M2_10",    dAdaptiveBinsM2, dAdaptiveArrayM2);

  // 50mK M1 and M2
  hAdap50mKco58M1      = new TH1D("hAdap50mKco58M1",   "hAdap50mKco58M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKco60M1      = new TH1D("hAdap50mKco60M1",   "hAdap50mKco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKcs137M1     = new TH1D("hAdap50mKcs137M1",  "hAdap50mKcs137M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKk40M1       = new TH1D("hAdap50mKk40M1",    "hAdap50mKk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKmn54M1      = new TH1D("hAdap50mKmn54M1",   "hAdap50mKmn54M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKpb210M1     = new TH1D("hAdap50mKpb210M1",  "hAdap50mKpb210M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKth232M1     = new TH1D("hAdap50mKth232M1",  "hAdap50mKth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdap50mKu238M1      = new TH1D("hAdap50mKu238M1",   "hAdap50mKu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdap50mKco58M2      = new TH1D("hAdap50mKco58M2",   "hAdap50mKco58M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap50mKco60M2      = new TH1D("hAdap50mKco60M2",   "hAdap50mKco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap50mKcs137M2     = new TH1D("hAdap50mKcs137M2",  "hAdap50mKcs137M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap50mKk40M2       = new TH1D("hAdap50mKk40M2",    "hAdap50mKk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap50mKmn54M2      = new TH1D("hAdap50mKmn54M2",   "hAdap50mKmn54M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap50mKpb210M2     = new TH1D("hAdap50mKpb210M2",  "hAdap50mKpb210M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap50mKth232M2     = new TH1D("hAdap50mKth232M2",  "hAdap50mKth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdap50mKu238M2      = new TH1D("hAdap50mKu238M2",   "hAdap50mKu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  // 600mK M1 and M2
  hAdap600mKco60M1      = new TH1D("hAdap600mKco60M1",   "hAdap600mKco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap600mKk40M1       = new TH1D("hAdap600mKk40M1",    "hAdap600mKk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap600mKth232M1     = new TH1D("hAdap600mKth232M1",  "hAdap600mKth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdap600mKu238M1      = new TH1D("hAdap600mKu238M1",   "hAdap600mKu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdap600mKco60M2      = new TH1D("hAdap600mKco60M2",   "hAdap600mKco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap600mKk40M2       = new TH1D("hAdap600mKk40M2",    "hAdap600mKk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap600mKth232M2     = new TH1D("hAdap600mKth232M2",  "hAdap600mKth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdap600mKu238M2      = new TH1D("hAdap600mKu238M2",   "hAdap600mKu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  

  // Roman Lead M1 and M2
  hAdapPbRombi207M1     = new TH1D("hAdapPbRombi207M1",  "hAdapPbRombi207M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapPbRomco60M1      = new TH1D("hAdapPbRomco60M1",   "hAdapPbRomco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapPbRomcs137M1     = new TH1D("hAdapPbRomcs137M1",  "hAdapPbRomcs137M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapPbRomk40M1       = new TH1D("hAdapPbRomk40M1",    "hAdapPbRomk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapPbRompb210M1     = new TH1D("hAdapPbRompb210M1",  "hAdapPbRompb210M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapPbRomth232M1     = new TH1D("hAdapPbRomth232M1",  "hAdapPbRomth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapPbRomu238M1      = new TH1D("hAdapPbRomu238M1",   "hAdapPbRomu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapPbRombi207M2     = new TH1D("hAdapPbRombi207M2",  "hAdapPbRombi207M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapPbRomco60M2      = new TH1D("hAdapPbRomco60M2",   "hAdapPbRomco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapPbRomcs137M2     = new TH1D("hAdapPbRomcs137M2",  "hAdapPbRomcs137M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapPbRomk40M2       = new TH1D("hAdapPbRomk40M2",    "hAdapPbRomk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapPbRompb210M2     = new TH1D("hAdapPbRompb210M2",  "hAdapPbRompb210M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapPbRomth232M2     = new TH1D("hAdapPbRomth232M2",  "hAdapPbRomth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapPbRomu238M2      = new TH1D("hAdapPbRomu238M2",   "hAdapPbRomu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  // Main bath M1 and M2
  hAdapMBco60M1      = new TH1D("hAdapMBco60M1",   "hAdapMBco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapMBk40M1       = new TH1D("hAdapMBk40M1",    "hAdapMBk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapMBth232M1     = new TH1D("hAdapMBth232M1",  "hAdapMBth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapMBu238M1      = new TH1D("hAdapMBu238M1",   "hAdapMBu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapMBco60M2      = new TH1D("hAdapMBco60M2",   "hAdapMBco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapMBk40M2       = new TH1D("hAdapMBk40M2",    "hAdapMBk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapMBth232M2     = new TH1D("hAdapMBth232M2",  "hAdapMBth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapMBu238M2      = new TH1D("hAdapMBu238M2",   "hAdapMBu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  

  // Super Insulation M1 and M2
  hAdapSIk40M1       = new TH1D("hAdapSIk40M1",    "hAdapSIk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1); 
  hAdapSIth232M1     = new TH1D("hAdapSIth232M1",  "hAdapSIth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapSIu238M1      = new TH1D("hAdapSIu238M1",   "hAdapSIu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapSIk40M2       = new TH1D("hAdapSIk40M2",    "hAdapSIk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapSIth232M2     = new TH1D("hAdapSIth232M2",  "hAdapSIth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapSIu238M2      = new TH1D("hAdapSIu238M2",   "hAdapSIu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);


  // IVC M1 and M2
  hAdapIVCco60M1      = new TH1D("hAdapIVCco60M1",   "hAdapIVCco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapIVCk40M1       = new TH1D("hAdapIVCk40M1",    "hAdapIVCk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapIVCth232M1     = new TH1D("hAdapIVCth232M1",  "hAdapIVCth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapIVCu238M1      = new TH1D("hAdapIVCu238M1",   "hAdapIVCu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapIVCco60M2      = new TH1D("hAdapIVCco60M2",   "hAdapIVCco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapIVCk40M2       = new TH1D("hAdapIVCk40M2",    "hAdapIVCk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapIVCth232M2     = new TH1D("hAdapIVCth232M2",  "hAdapIVCth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapIVCu238M2      = new TH1D("hAdapIVCu238M2",   "hAdapIVCu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  

  // OVC M1 and M2
  hAdapOVCco60M1      = new TH1D("hAdapOVCco60M1",   "hAdapOVCco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapOVCk40M1       = new TH1D("hAdapOVCk40M1",    "hAdapOVCk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapOVCth232M1     = new TH1D("hAdapOVCth232M1",  "hAdapOVCth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapOVCu238M1      = new TH1D("hAdapOVCu238M1",   "hAdapOVCu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapOVCco60M2      = new TH1D("hAdapOVCco60M2",   "hAdapOVCco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapOVCk40M2       = new TH1D("hAdapOVCk40M2",    "hAdapOVCk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapOVCth232M2     = new TH1D("hAdapOVCth232M2",  "hAdapOVCth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapOVCu238M2      = new TH1D("hAdapOVCu238M2",   "hAdapOVCu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  

  hEnergyScaleDummyM1 = new TH1D("hEnergyScaleDummyM1",   "Energy Scale M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hEnergyScaleDummyM2 = new TH1D("hEnergyScaleDummyM2",   "Energy Scale M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);



  // Loads all of the PDFs from file
  Initialize();

  // Do the fit now 
  // DoTheFitAdaptive();
}
  
// Needs updating  
TBackgroundModel::~TBackgroundModel()
{
 
  delete fDataHistoTot;
  delete fDataHistoM1;
  delete fDataHistoM2;

  delete fAdapDataHistoM1;
  delete fAdapDataHistoM2;


  // To be updated 11-05-2014
  // Total PDFs M1
  delete fModelTotM1;
  delete fModelTotthM1;
  delete fModelTotraM1;
  delete fModelTotkM1;
  delete fModelTotcoM1;
  delete fModelTotmnM1;
  delete fModelTotNDBDM1;
  delete fModelTot2NDBDM1;
  delete fModelTotbiM1;
  delete fModelTotptM1;
  delete fModelTotpbM1;

  delete fModelTotAdapM1;
  delete fModelTotAdapThM1;
  delete fModelTotAdapRaM1;
  delete fModelTotAdapKM1;
  delete fModelTotAdapCoM1;
  delete fModelTotAdapMnM1;
  delete fModelTotAdapNDBDM1;
  delete fModelTotAdap2NDBDM1;
  delete fModelTotAdapBiM1;
  delete fModelTotAdapPtM1;
  delete fModelTotAdapPbM1;

  // Total PDFs M2
  delete fModelTotM2;
  delete fModelTotThM2;
  delete fModelTotRaM2;
  delete fModelTotKM2;
  delete fModelTotCoM2;
  delete fModelTotMnM2;
  delete fModelTotNDBDM2;
  delete fModelTot2NDBDM2;
  delete fModelTotBiM2;
  delete fModelTotPtM2;
  delete fModelTotPbM2;

  delete fModelTotAdapM2;
  delete fModelTotAdapThM2;
  delete fModelTotAdapRaM2;
  delete fModelTotAdapKM2;
  delete fModelTotAdapCoM2;
  delete fModelTotAdapMnM2;
  delete fModelTotAdapNDBDM2;
  delete fModelTotAdap2NDBDM2;
  delete fModelTotAdapBiM2;
  delete fModelTotAdapPtM2;
  delete fModelTotAdapPbM2;

  delete hResidualDistM1;
  delete hResidualDistM2;

  delete hResidualGausM1;
  delete hResidualGausM2;



  delete fModelDummyM1;

  // Crystal M1 and M2
  delete hTeO20nuM1;
  delete hTeO22nuM1;
  delete hTeO2co60M1;
  delete hTeO2k40M1;
  delete hTeO2pb210M1;
  delete hTeO2po210M1;
  delete hTeO2te125M1;
  delete hTeO2th232M1;
  delete hTeO2th228M1;
  delete hTeO2ra226M1;
  delete hTeO2rn222M1;
  delete hTeO2u238M1;
  delete hTeO2th230M1;
  delete hTeO2u234M1;

  delete hTeO20nuM2;
  delete hTeO22nuM2;
  delete hTeO2co60M2;
  delete hTeO2k40M2;
  delete hTeO2pb210M2;
  delete hTeO2po210M2;
  delete hTeO2te125M2;
  delete hTeO2th232M2;
  delete hTeO2th228M2;
  delete hTeO2ra226M2;
  delete hTeO2rn222M2;
  delete hTeO2u238M2;
  delete hTeO2th230M2;
  delete hTeO2u234M2;

  // Frame M1 and M2
  delete hCuFrameco58M1;
  delete hCuFrameco60M1;
  delete hCuFramecs137M1;
  delete hCuFramek40M1;
  delete hCuFramemn54M1;
  delete hCuFramepb210M1;
  delete hCuFrameth232M1;
  delete hCuFrameu238M1;

  delete hCuFrameco58M2;
  delete hCuFrameco60M2;
  delete hCuFramecs137M2;
  delete hCuFramek40M2;
  delete hCuFramemn54M2;
  delete hCuFramepb210M2;
  delete hCuFrameth232M2;
  delete hCuFrameu238M2;

  // CuBox (TShield) M1 and M2
  delete hCuBoxco58M1;
  delete hCuBoxco60M1;
  delete hCuBoxcs137M1;
  delete hCuBoxk40M1;
  delete hCuBoxmn54M1;
  delete hCuBoxpb210M1;
  delete hCuBoxth232M1;
  delete hCuBoxu238M1;  

  delete hCuBoxco58M2;
  delete hCuBoxco60M2;
  delete hCuBoxcs137M2;
  delete hCuBoxk40M2;
  delete hCuBoxmn54M2;
  delete hCuBoxpb210M2;
  delete hCuBoxth232M2;
  delete hCuBoxu238M2;  

  // 50mK M1 and M2
  delete h50mKco58M1;
  delete h50mKco60M1;
  delete h50mKcs137M1;
  delete h50mKk40M1;
  delete h50mKmn54M1;
  delete h50mKpb210M1;
  delete h50mKth232M1;
  delete h50mKu238M1;   

  delete h50mKco58M2;
  delete h50mKco60M2;
  delete h50mKcs137M2;
  delete h50mKk40M2;
  delete h50mKmn54M2;
  delete h50mKpb210M2;
  delete h50mKth232M2;
  delete h50mKu238M2; 

  // 600mK M1 and M2
  delete h600mKco60M1;
  delete h600mKk40M1;
  delete h600mKth232M1;
  delete h600mKu238M1;    

  delete h600mKco60M2;
  delete h600mKk40M2;
  delete h600mKth232M2;
  delete h600mKu238M2;  

  // Roman Lead M1 and M2
  delete hPbRombi207M1;
  delete hPbRomco60M1;
  delete hPbRomcs137M1;
  delete hPbRomk40M1;
  delete hPbRompb210M1;
  delete hPbRomth232M1;
  delete hPbRomu238M1;    

  delete hPbRombi207M2;
  delete hPbRomco60M2;
  delete hPbRomcs137M2;
  delete hPbRomk40M2;
  delete hPbRompb210M2;
  delete hPbRomth232M2;
  delete hPbRomu238M2;    


  // Main Bath M1 and M2
  delete hMBco60M1;
  delete hMBk40M1;
  delete hMBth232M1;
  delete hMBu238M1;   

  delete hMBco60M2;
  delete hMBk40M2;
  delete hMBth232M2;
  delete hMBu238M2; 

  // IVC M1 and M2
  delete hIVCco60M1;
  delete hIVCk40M1;
  delete hIVCth232M1;
  delete hIVCu238M1;    

  delete hIVCco60M2;
  delete hIVCk40M2;
  delete hIVCth232M2;
  delete hIVCu238M2;  

  // OVC M1 and M2
  delete hOVCco60M1;
  delete hOVCk40M1;
  delete hOVCth232M1;
  delete hOVCu238M1;    

  delete hOVCco60M2;
  delete hOVCk40M2;
  delete hOVCth232M2;
  delete hOVCu238M2;  
  // Crystal M1 and M2
  delete hAdapTeO20nuM1;
  delete hAdapTeO22nuM1;
  delete hAdapTeO2co60M1;
  delete hAdapTeO2k40M1;
  delete hAdapTeO2pb210M1;
  delete hAdapTeO2po210M1;
  delete hAdapTeO2te125M1;
  delete hAdapTeO2th232M1;
  delete hAdapTeO2th228M1;
  delete hAdapTeO2ra226M1;
  delete hAdapTeO2rn222M1;
  delete hAdapTeO2u238M1;
  delete hAdapTeO2th230M1;
  delete hAdapTeO2u234M1;

  delete hAdapTeO20nuM2;
  delete hAdapTeO22nuM2;
  delete hAdapTeO2co60M2;
  delete hAdapTeO2k40M2;
  delete hAdapTeO2pb210M2;
  delete hAdapTeO2po210M2;
  delete hAdapTeO2te125M2;
  delete hAdapTeO2th232M2;
  delete hAdapTeO2th228M2;
  delete hAdapTeO2ra226M2;
  delete hAdapTeO2rn222M2;
  delete hAdapTeO2u238M2;
  delete hAdapTeO2th230M2;
  delete hAdapTeO2u234M2;

  // Frame M1 and M2
  delete hAdapCuFrameco58M1;
  delete hAdapCuFrameco60M1;
  delete hAdapCuFramecs137M1;
  delete hAdapCuFramek40M1;
  delete hAdapCuFramemn54M1;
  delete hAdapCuFramepb210M1;
  delete hAdapCuFrameth232M1;
  delete hAdapCuFrameu238M1;

  delete hAdapCuFrameco58M2;
  delete hAdapCuFrameco60M2;
  delete hAdapCuFramecs137M2;
  delete hAdapCuFramek40M2;
  delete hAdapCuFramemn54M2;
  delete hAdapCuFramepb210M2;
  delete hAdapCuFrameth232M2;
  delete hAdapCuFrameu238M2;

  // CuBox (TShield) M1 and M2
  delete hAdapCuBoxco58M1;
  delete hAdapCuBoxco60M1;
  delete hAdapCuBoxcs137M1;
  delete hAdapCuBoxk40M1;
  delete hAdapCuBoxmn54M1;
  delete hAdapCuBoxpb210M1;
  delete hAdapCuBoxth232M1;
  delete hAdapCuBoxu238M1;  

  delete hAdapCuBoxco58M2;
  delete hAdapCuBoxco60M2;
  delete hAdapCuBoxcs137M2;
  delete hAdapCuBoxk40M2;
  delete hAdapCuBoxmn54M2;
  delete hAdapCuBoxpb210M2;
  delete hAdapCuBoxth232M2;
  delete hAdapCuBoxu238M2;  

  // 50mK M1 and M2
  delete hAdap50mKco58M1;
  delete hAdap50mKco60M1;
  delete hAdap50mKcs137M1;
  delete hAdap50mKk40M1;
  delete hAdap50mKmn54M1;
  delete hAdap50mKpb210M1;
  delete hAdap50mKth232M1;
  delete hAdap50mKu238M1;   

  delete hAdap50mKco58M2;
  delete hAdap50mKco60M2;
  delete hAdap50mKcs137M2;
  delete hAdap50mKk40M2;
  delete hAdap50mKmn54M2;
  delete hAdap50mKpb210M2;
  delete hAdap50mKth232M2;
  delete hAdap50mKu238M2; 

  // 600mK M1 and M2
  delete hAdap600mKco60M1;
  delete hAdap600mKk40M1;
  delete hAdap600mKth232M1;
  delete hAdap600mKu238M1;    

  delete hAdap600mKco60M2;
  delete hAdap600mKk40M2;
  delete hAdap600mKth232M2;
  delete hAdap600mKu238M2;  

  // Roman Lead M1 and M2
  delete hAdapPbRombi207M1;
  delete hAdapPbRomco60M1;
  delete hAdapPbRomcs137M1;
  delete hAdapPbRomk40M1;
  delete hAdapPbRompb210M1;
  delete hAdapPbRomth232M1;
  delete hAdapPbRomu238M1;    

  delete hAdapPbRombi207M2;
  delete hAdapPbRomco60M2;
  delete hAdapPbRomcs137M2;
  delete hAdapPbRomk40M2;
  delete hAdapPbRompb210M2;
  delete hAdapPbRomth232M2;
  delete hAdapPbRomu238M2;    


  // Main Bath M1 and M2
  delete hAdapMBco60M1;
  delete hAdapMBk40M1;
  delete hAdapMBth232M1;
  delete hAdapMBu238M1;   

  delete hAdapMBco60M2;
  delete hAdapMBk40M2;
  delete hAdapMBth232M2;
  delete hAdapMBu238M2; 

  // IVC M1 and M2
  delete hAdapIVCco60M1;
  delete hAdapIVCk40M1;
  delete hAdapIVCth232M1;
  delete hAdapIVCu238M1;    

  delete hAdapIVCco60M2;
  delete hAdapIVCk40M2;
  delete hAdapIVCth232M2;
  delete hAdapIVCu238M2;  

  // OVC M1 and M2
  delete hAdapOVCco60M1;
  delete hAdapOVCk40M1;
  delete hAdapOVCth232M1;
  delete hAdapOVCu238M1;    

  delete hAdapOVCco60M2;
  delete hAdapOVCk40M2;
  delete hAdapOVCth232M2;
  delete hAdapOVCu238M2;  
}

// Returns vector of bin low-edge for adaptive binning
vector<double> TBackgroundModel::AdaptiveBinning(TH1D *h1, int dBinBase)
{
  vector<double> dBinArrayThing;

  double dDummy = 0;
  double dDummyFill = 0;
  int j = 0;
  // 25 since start at 50 keV with 2 keV bin width
  for(int i = 1; i < 25; i++)
  {
    dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(i));
  }

  for(int i = 25; i < dNBins; i++)
  {
    dDummy = h1->GetBinContent(i);
    dDummyFill += dDummy;

    if(dDummyFill >= dBinBase)
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

TH1D *TBackgroundModel::CalculateResidualsAdaptive(TH1D *h1, TH1D *h2, TH1D *hResid, int binMin, int binMax, int dMult)
{

  if(binMin >= binMax)
  {
    cout << " min bin >= max bin" << endl;
    break;
  }

  if(dMult == 1)
  {
  TH1D  *hOut       = new TH1D("hOutResidualM1", "Fit Residuals M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  }

  if(dMult == 2)
  {
  TH1D  *hOut       = new TH1D("hOutResidualM2", "Fit Residuals M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  }

  // Clone histograms for rebinning
  TH1D  *hCloneBkg    = (TH1D*)h1->Clone("hCloneBkg");
  TH1D  *hCloneMC   = (TH1D*)h2->Clone("hCloneMC");

  // Variables used in Residual calculations
  double dResidualX = 0, dResidualY = 0, dResidualXErr = 0, dResidualYErr = 0;

  // binMin = 50;

  // Residual plot and distribution
  for(int j = binMin ; j < binMax ; j++)
  {


    dResidualX    = hCloneBkg->GetBinCenter(j);
    dResidualY    = (hCloneBkg->GetBinContent(j) - hCloneMC->GetBinContent(j))*hCloneBkg->GetBinWidth(j)/dBinSize /
              (TMath::Sqrt((hCloneBkg->GetBinContent(j)+hCloneMC->GetBinContent(j))*hCloneBkg->GetBinWidth(j)/dBinSize) ); // Sqrt(MC + data) = sigma for poisson distribution
    // dResidualY  = (hCloneBkg->GetBinContent(j) - hCloneMC->GetBinContent(j))*hCloneBkg->GetBinWidth(j)/dBinSize; // not normalized
    // g1->SetPoint(j, dResidualX, dResidualY);
    hOut->SetBinContent(j, dResidualY);
    hOut->SetBinError(j, 0.01);
    // cout << "bin: " << j << " Energy: " << hCloneBkg->GetBinCenter(j) << " Bkg: " <<  hCloneBkg->GetBinContent(j) << " MC: " << hCloneMC->GetBinContent(j) << endl;
    // cout << "   Residual: " << dResidualY << endl;
    hResid->Fill(dResidualY);
  }
  return hOut;
}

// Draws background spectrum
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

  // TCanvas *cBkgAdap1 = new TCanvas("cBkgAdap1", "cBkgAdap1", 1200, 800);
  // cBkgAdap1->SetLogy();
  // fAdapDataHistoM1->SetLineColor(1);
  // fAdapDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
  // fAdapDataHistoM1->GetYaxis()->SetTitle("Counts/Bin");
  // fAdapDataHistoM1->Draw();

  // TCanvas *cBkgAdap2 = new TCanvas("cBkgAdap2", "cBkgAdap2", 1200, 800);
  // cBkgAdap2->SetLogy();
  // fAdapDataHistoM2->SetLineColor(1);
  // fAdapDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
  // fAdapDataHistoM2->GetYaxis()->SetTitle("Counts/Bin");
  // fAdapDataHistoM2->Draw();
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

// This function applies an energy scaling of E' = dSlope*E + dConst
TH1D *TBackgroundModel::EnergyScale(TH1D *hIn, TH1D *hDummy, double dConst, double dSlope)
{
  hDummy->Reset();
  int fBins = hDummy->GetNbinsX();

    for(int i = 1; i <= fBins; i++)
    {
      hDummy->Fill(dConst + dSlope*hIn->GetBinCenter(i), hIn->GetBinContent(i));
    }

  return hDummy;
}

double TBackgroundModel::GetChiSquareAdaptive()
{
  double chiSquare = 0.;
  double datam1_i, errm1_i;
  double datam2_i, errm2_i;
  double modelm1_i, modelm2_i;  

  for(int i = dFitMinBinM1 ; i <= dFitMaxBinM1; i++)
  {

    if( fAdapDataHistoM1->GetBinCenter(i) >= 3200 && fAdapDataHistoM1->GetBinCenter(i) <= 3400)continue;
    datam1_i = fAdapDataHistoM1->GetBinContent(i)*fAdapDataHistoM1->GetBinWidth(i)/dBinSize; // For real data

    modelm1_i = fModelTotAdapM1->GetBinContent(i)*fAdapDataHistoM1->GetBinWidth(i)/dBinSize;

    if(modelm1_i != 0 && datam1_i != 0)
    {
      chiSquare += 2 * (modelm1_i - datam1_i + datam1_i * TMath::Log(datam1_i/modelm1_i));
    }
  }

  for(int i = dFitMinBinM2; i <= dFitMaxBinM2; i++)
  {
    // Pt peak can be ignored in M2 since doesn't exist?
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 3200 && fAdapDataHistoM2->GetBinCenter(i) <= 3500)continue;
    datam2_i = fAdapDataHistoM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i)/dBinSize; // For real data

    modelm2_i = fModelTotAdapM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i)/dBinSize;

    if(modelm2_i != 0 && datam2_i != 0)
    {
      chiSquare += 2 * (modelm2_i - datam2_i + datam2_i * TMath::Log(datam2_i/modelm2_i));
    }

  }


  return chiSquare;
}

void TBackgroundModel::Initialize()
{	

  // Loads PDFs from file
  cout << "Loading PDF Histograms from file" << endl;
  // fFile = new TFile("/Users/brian/macros/Simulations/Production/MCProduction_Bulk_nonnormalized.root"); 
  // fFile2 = new TFile("/Users/brian/macros/Simulations/Production/MCProduction_Surface_nonnormalized.root");

  fFile = new TFile("/Users/brian/macros/Simulations/Production/MCProduction_Bulk.root"); // These are normalized PDFs
  fFile2 = new TFile("/Users/brian/macros/Simulations/Production/MCProduction_Surface.root");

///////////// Bulk Histograms
  // Crystal M1 and M2
  hTeO20nuM1     = (TH1D*)fFile->Get("hTeO20nuM1");
  hTeO22nuM1     = (TH1D*)fFile->Get("hTeO22nuM1");
  hTeO2co60M1    = (TH1D*)fFile->Get("hTeO2co60M1");
  hTeO2k40M1     = (TH1D*)fFile->Get("hTeO2k40M1");
  hTeO2pb210M1   = (TH1D*)fFile->Get("hTeO2pb210M1");
  hTeO2po210M1   = (TH1D*)fFile->Get("hTeO2po210M1");
  hTeO2te125M1   = (TH1D*)fFile->Get("hTeO2te125M1");
  hTeO2th232M1   = (TH1D*)fFile->Get("hTeO2th232M1");
  hTeO2th228M1   = (TH1D*)fFile->Get("hTeO2th228M1");
  hTeO2ra226M1   = (TH1D*)fFile->Get("hTeO2ra226M1");
  hTeO2rn222M1   = (TH1D*)fFile->Get("hTeO2rn222M1");
  hTeO2u238M1    = (TH1D*)fFile->Get("hTeO2u238M1");
  hTeO2th230M1   = (TH1D*)fFile->Get("hTeO2th230M1");
  hTeO2u234M1    = (TH1D*)fFile->Get("hTeO2u234M1");

  hTeO20nuM2     = (TH1D*)fFile->Get("hTeO20nuM2");
  hTeO22nuM2     = (TH1D*)fFile->Get("hTeO22nuM2");
  hTeO2co60M2    = (TH1D*)fFile->Get("hTeO2co60M2");
  hTeO2k40M2     = (TH1D*)fFile->Get("hTeO2k40M2");
  hTeO2pb210M2   = (TH1D*)fFile->Get("hTeO2pb210M2");
  hTeO2po210M2   = (TH1D*)fFile->Get("hTeO2po210M2");
  hTeO2te125M2   = (TH1D*)fFile->Get("hTeO2te125M2");
  hTeO2th232M2   = (TH1D*)fFile->Get("hTeO2th232M2");
  hTeO2th228M2   = (TH1D*)fFile->Get("hTeO2th228M2");
  hTeO2ra226M2   = (TH1D*)fFile->Get("hTeO2ra226M2");
  hTeO2rn222M2   = (TH1D*)fFile->Get("hTeO2rn222M2");
  hTeO2u238M2    = (TH1D*)fFile->Get("hTeO2u238M2");
  hTeO2th230M2   = (TH1D*)fFile->Get("hTeO2th230M2");
  hTeO2u234M2    = (TH1D*)fFile->Get("hTeO2u234M2");

  // Frame M1 and M2
  hCuFrameco58M1    = (TH1D*)fFile->Get("hCuFrameco58M1");
  hCuFrameco60M1    = (TH1D*)fFile->Get("hCuFrameco60M1");
  hCuFramecs137M1   = (TH1D*)fFile->Get("hCuFramecs137M1");
  hCuFramek40M1     = (TH1D*)fFile->Get("hCuFramek40M1");
  hCuFramemn54M1    = (TH1D*)fFile->Get("hCuFramemn54M1");
  hCuFramepb210M1   = (TH1D*)fFile->Get("hCuFramepb210M1");
  hCuFrameth232M1   = (TH1D*)fFile->Get("hCuFrameth232M1");
  hCuFrameu238M1    = (TH1D*)fFile->Get("hCuFrameu238M1");
   
  hCuFrameco58M2    = (TH1D*)fFile->Get("hCuFrameco58M2");
  hCuFrameco60M2    = (TH1D*)fFile->Get("hCuFrameco60M2");
  hCuFramecs137M2   = (TH1D*)fFile->Get("hCuFramecs137M2");
  hCuFramek40M2     = (TH1D*)fFile->Get("hCuFramek40M2");
  hCuFramemn54M2    = (TH1D*)fFile->Get("hCuFramemn54M2");
  hCuFramepb210M2   = (TH1D*)fFile->Get("hCuFramepb210M2");
  hCuFrameth232M2   = (TH1D*)fFile->Get("hCuFrameth232M2");
  hCuFrameu238M2    = (TH1D*)fFile->Get("hCuFrameu238M2");

  // CuBox (TShield) M1 and M2
  hCuBoxco58M1    = (TH1D*)fFile->Get("hCuBoxco58M1");
  hCuBoxco60M1    = (TH1D*)fFile->Get("hCuBoxco60M1");
  hCuBoxcs137M1   = (TH1D*)fFile->Get("hCuBoxcs137M1");
  hCuBoxk40M1     = (TH1D*)fFile->Get("hCuBoxk40M1");
  hCuBoxmn54M1    = (TH1D*)fFile->Get("hCuBoxmn54M1");
  hCuBoxpb210M1   = (TH1D*)fFile->Get("hCuBoxpb210M1");
  hCuBoxth232M1   = (TH1D*)fFile->Get("hCuBoxth232M1");
  hCuBoxu238M1    = (TH1D*)fFile->Get("hCuBoxu238M1");
   
  hCuBoxco58M2    = (TH1D*)fFile->Get("hCuBoxco58M2");
  hCuBoxco60M2    = (TH1D*)fFile->Get("hCuBoxco60M2");
  hCuBoxcs137M2   = (TH1D*)fFile->Get("hCuBoxcs137M2");
  hCuBoxk40M2     = (TH1D*)fFile->Get("hCuBoxk40M2");
  hCuBoxmn54M2    = (TH1D*)fFile->Get("hCuBoxmn54M2");
  hCuBoxpb210M2   = (TH1D*)fFile->Get("hCuBoxpb210M2");
  hCuBoxth232M2   = (TH1D*)fFile->Get("hCuBoxth232M2");
  hCuBoxu238M2    = (TH1D*)fFile->Get("hCuBoxu238M2");

  // 50mK M1 and M2
  h50mKco58M1    = (TH1D*)fFile->Get("h50mKco58M1");
  h50mKco60M1    = (TH1D*)fFile->Get("h50mKco60M1");
  h50mKcs137M1   = (TH1D*)fFile->Get("h50mKcs137M1");
  h50mKk40M1     = (TH1D*)fFile->Get("h50mKk40M1");
  h50mKmn54M1    = (TH1D*)fFile->Get("h50mKmn54M1");
  h50mKpb210M1   = (TH1D*)fFile->Get("h50mKpb210M1");
  h50mKth232M1   = (TH1D*)fFile->Get("h50mKth232M1");
  h50mKu238M1    = (TH1D*)fFile->Get("h50mKu238M1");
   
  h50mKco58M2    = (TH1D*)fFile->Get("h50mKco58M2");
  h50mKco60M2    = (TH1D*)fFile->Get("h50mKco60M2");
  h50mKcs137M2   = (TH1D*)fFile->Get("h50mKcs137M2");
  h50mKk40M2     = (TH1D*)fFile->Get("h50mKk40M2");
  h50mKmn54M2    = (TH1D*)fFile->Get("h50mKmn54M2");
  h50mKpb210M2   = (TH1D*)fFile->Get("h50mKpb210M2");
  h50mKth232M2   = (TH1D*)fFile->Get("h50mKth232M2");
  h50mKu238M2    = (TH1D*)fFile->Get("h50mKu238M2");


  // 600mK M1 and M2
  h600mKco60M1    = (TH1D*)fFile->Get("h600mKco60M1");
  h600mKk40M1     = (TH1D*)fFile->Get("h600mKk40M1");
  h600mKth232M1   = (TH1D*)fFile->Get("h600mKth232M1");
  h600mKu238M1    = (TH1D*)fFile->Get("h600mKu238M1");
   
  h600mKco60M2    = (TH1D*)fFile->Get("h600mKco60M2");
  h600mKk40M2     = (TH1D*)fFile->Get("h600mKk40M2");
  h600mKth232M2   = (TH1D*)fFile->Get("h600mKth232M2");
  h600mKu238M2    = (TH1D*)fFile->Get("h600mKu238M2");


  // Roman Lead M1 and M2
  hPbRombi207M1   = (TH1D*)fFile->Get("hPbRombi207M1");
  hPbRomco60M1    = (TH1D*)fFile->Get("hPbRomco60M1");
  hPbRomcs137M1   = (TH1D*)fFile->Get("hPbRomcs137M1");
  hPbRomk40M1     = (TH1D*)fFile->Get("hPbRomk40M1");
  hPbRompb210M1   = (TH1D*)fFile->Get("hPbRompb210M1");
  hPbRomth232M1   = (TH1D*)fFile->Get("hPbRomth232M1");
  hPbRomu238M1    = (TH1D*)fFile->Get("hPbRomu238M1");
   
  hPbRombi207M2   = (TH1D*)fFile->Get("hPbRombi207M2");
  hPbRomco60M2    = (TH1D*)fFile->Get("hPbRomco60M2");
  hPbRomcs137M2   = (TH1D*)fFile->Get("hPbRomcs137M2");
  hPbRomk40M2     = (TH1D*)fFile->Get("hPbRomk40M2");
  hPbRompb210M2   = (TH1D*)fFile->Get("hPbRompb210M2");
  hPbRomth232M2   = (TH1D*)fFile->Get("hPbRomth232M2");
  hPbRomu238M2    = (TH1D*)fFile->Get("hPbRomu238M2");

  // Main Bath M1 and M2
  hMBco60M1    = (TH1D*)fFile->Get("hMBco60M1");
  hMBk40M1     = (TH1D*)fFile->Get("hMBk40M1");
  hMBth232M1   = (TH1D*)fFile->Get("hMBth232M1");
  hMBu238M1    = (TH1D*)fFile->Get("hMBu238M1");
   
  hMBco60M2    = (TH1D*)fFile->Get("hMBco60M2");
  hMBk40M2     = (TH1D*)fFile->Get("hMBk40M2");
  hMBth232M2   = (TH1D*)fFile->Get("hMBth232M2");
  hMBu238M2    = (TH1D*)fFile->Get("hMBu238M2");

  // Super Insulation M1 and M2
  hSIk40M1     = (TH1D*)fFile->Get("hSIk40M1");
  hSIth232M1   = (TH1D*)fFile->Get("hSIth232M1");
  hSIu238M1    = (TH1D*)fFile->Get("hSIu238M1");  

  hSIk40M2     = (TH1D*)fFile->Get("hSIk40M2");
  hSIth232M2   = (TH1D*)fFile->Get("hSIth232M2");
  hSIu238M2    = (TH1D*)fFile->Get("hSIu238M2"); 

  // IVC M1 and M2
  hIVCco60M1    = (TH1D*)fFile->Get("hIVCco60M1");
  hIVCk40M1     = (TH1D*)fFile->Get("hIVCk40M1");
  hIVCth232M1   = (TH1D*)fFile->Get("hIVCth232M1");
  hIVCu238M1    = (TH1D*)fFile->Get("hIVCu238M1");
   
  hIVCco60M2    = (TH1D*)fFile->Get("hIVCco60M2");
  hIVCk40M2     = (TH1D*)fFile->Get("hIVCk40M2");
  hIVCth232M2   = (TH1D*)fFile->Get("hIVCth232M2");
  hIVCu238M2    = (TH1D*)fFile->Get("hIVCu238M2");

  // OVC M1 and M2
  hOVCco60M1    = (TH1D*)fFile->Get("hOVCco60M1");
  hOVCk40M1     = (TH1D*)fFile->Get("hOVCk40M1");
  hOVCth232M1   = (TH1D*)fFile->Get("hOVCth232M1");
  hOVCu238M1    = (TH1D*)fFile->Get("hOVCu238M1");
   
  hOVCco60M2    = (TH1D*)fFile->Get("hOVCco60M2");
  hOVCk40M2     = (TH1D*)fFile->Get("hOVCk40M2");
  hOVCth232M2   = (TH1D*)fFile->Get("hOVCth232M2");
  hOVCu238M2    = (TH1D*)fFile->Get("hOVCu238M2");

//////// Surface PDFs
  // Crystal M1 and M2
  hTeO2Spb210M1_01    = (TH1D*)fFile2->Get("hTeO2Spb210M1_01");
  hTeO2Spo210M1_001   = (TH1D*)fFile2->Get("hTeO2Spo210M1_001");
  hTeO2Spo210M1_01    = (TH1D*)fFile2->Get("hTeO2Spo210M1_01");
  hTeO2Sth232M1_01    = (TH1D*)fFile2->Get("hTeO2Sth232M1_01");
  hTeO2Su238M1_01     = (TH1D*)fFile2->Get("hTeO2Su238M1_01");
  hTeO2Sxpb210M1_001  = (TH1D*)fFile2->Get("hTeO2Sxpb210M1_001");
  hTeO2Sxpb210M1_01   = (TH1D*)fFile2->Get("hTeO2Sxpb210M1_01");
  hTeO2Sxpb210M1_1    = (TH1D*)fFile2->Get("hTeO2Sxpb210M1_1");
  hTeO2Sxpb210M1_10   = (TH1D*)fFile2->Get("hTeO2Sxpb210M1_10");
  hTeO2Sxpo210M1_001  = (TH1D*)fFile2->Get("hTeO2Sxpo210M1_001");
  hTeO2Sxpo210M1_01   = (TH1D*)fFile2->Get("hTeO2Sxpo210M1_01");
  hTeO2Sxpo210M1_1    = (TH1D*)fFile2->Get("hTeO2Sxpo210M1_1");
  hTeO2Sxth232M1_001  = (TH1D*)fFile2->Get("hTeO2Sxth232M1_001");
  hTeO2Sxth232M1_01   = (TH1D*)fFile2->Get("hTeO2Sxth232M1_01");
  hTeO2Sxth232M1_1    = (TH1D*)fFile2->Get("hTeO2Sxth232M1_1");
  hTeO2Sxth232M1_10   = (TH1D*)fFile2->Get("hTeO2Sxth232M1_10");
  hTeO2Sxu238M1_001   = (TH1D*)fFile2->Get("hTeO2Sxu238M1_001");
  hTeO2Sxu238M1_01    = (TH1D*)fFile2->Get("hTeO2Sxu238M1_01");
  hTeO2Sxu238M1_1     = (TH1D*)fFile2->Get("hTeO2Sxu238M1_1");
  hTeO2Sxu238M1_10    = (TH1D*)fFile2->Get("hTeO2Sxu238M1_10");

  hTeO2Spb210M2_01    = (TH1D*)fFile2->Get("hTeO2Spb210M2_01");
  hTeO2Spo210M2_001   = (TH1D*)fFile2->Get("hTeO2Spo210M2_001");
  hTeO2Spo210M2_01    = (TH1D*)fFile2->Get("hTeO2Spo210M2_01");
  hTeO2Sth232M2_01    = (TH1D*)fFile2->Get("hTeO2Sth232M2_01");
  hTeO2Su238M2_01     = (TH1D*)fFile2->Get("hTeO2Su238M2_01");
  hTeO2Sxpb210M2_001  = (TH1D*)fFile2->Get("hTeO2Sxpb210M2_001");
  hTeO2Sxpb210M2_01   = (TH1D*)fFile2->Get("hTeO2Sxpb210M2_01");
  hTeO2Sxpb210M2_1    = (TH1D*)fFile2->Get("hTeO2Sxpb210M2_1");
  hTeO2Sxpb210M2_10   = (TH1D*)fFile2->Get("hTeO2Sxpb210M2_10");
  hTeO2Sxpo210M2_001  = (TH1D*)fFile2->Get("hTeO2Sxpo210M2_001");
  hTeO2Sxpo210M2_01   = (TH1D*)fFile2->Get("hTeO2Sxpo210M2_01");
  hTeO2Sxpo210M2_1    = (TH1D*)fFile2->Get("hTeO2Sxpo210M2_1");
  hTeO2Sxth232M2_001  = (TH1D*)fFile2->Get("hTeO2Sxth232M2_001");
  hTeO2Sxth232M2_01   = (TH1D*)fFile2->Get("hTeO2Sxth232M2_01");
  hTeO2Sxth232M2_1    = (TH1D*)fFile2->Get("hTeO2Sxth232M2_1");
  hTeO2Sxth232M2_10   = (TH1D*)fFile2->Get("hTeO2Sxth232M2_10");
  hTeO2Sxu238M2_001   = (TH1D*)fFile2->Get("hTeO2Sxu238M2_001");
  hTeO2Sxu238M2_01    = (TH1D*)fFile2->Get("hTeO2Sxu238M2_01");
  hTeO2Sxu238M2_1     = (TH1D*)fFile2->Get("hTeO2Sxu238M2_1");
  hTeO2Sxu238M2_10    = (TH1D*)fFile2->Get("hTeO2Sxu238M2_10");

  // Frame M1 and M2
  hCuFrameSth232M1_1    = (TH1D*)fFile2->Get("hCuFrameSth232M1_1");
  hCuFrameSu238M1_1     = (TH1D*)fFile2->Get("hCuFrameSu238M1_1");
  hCuFrameSxpb210M1_001 = (TH1D*)fFile2->Get("hCuFrameSxpb210M1_001");
  hCuFrameSxpb210M1_01  = (TH1D*)fFile2->Get("hCuFrameSxpb210M1_01");
  hCuFrameSxpb210M1_1   = (TH1D*)fFile2->Get("hCuFrameSxpb210M1_1");
  hCuFrameSxpb210M1_10  = (TH1D*)fFile2->Get("hCuFrameSxpb210M1_10");
  hCuFrameSxth232M1_001 = (TH1D*)fFile2->Get("hCuFrameSxth232M1_001");
  hCuFrameSxth232M1_01  = (TH1D*)fFile2->Get("hCuFrameSxth232M1_01");
  hCuFrameSxth232M1_1   = (TH1D*)fFile2->Get("hCuFrameSxth232M1_1");
  hCuFrameSxth232M1_10  = (TH1D*)fFile2->Get("hCuFrameSxth232M1_10");
  hCuFrameSxu238M1_001  = (TH1D*)fFile2->Get("hCuFrameSxu238M1_001");
  hCuFrameSxu238M1_01   = (TH1D*)fFile2->Get("hCuFrameSxu238M1_01");
  hCuFrameSxu238M1_1    = (TH1D*)fFile2->Get("hCuFrameSxu238M1_1");
  hCuFrameSxu238M1_10   = (TH1D*)fFile2->Get("hCuFrameSxu238M1_10");

  hCuFrameSth232M2_1    = (TH1D*)fFile2->Get("hCuFrameSth232M2_1");
  hCuFrameSu238M2_1     = (TH1D*)fFile2->Get("hCuFrameSu238M2_1");
  hCuFrameSxpb210M2_001 = (TH1D*)fFile2->Get("hCuFrameSxpb210M2_001");
  hCuFrameSxpb210M2_01  = (TH1D*)fFile2->Get("hCuFrameSxpb210M2_01");
  hCuFrameSxpb210M2_1   = (TH1D*)fFile2->Get("hCuFrameSxpb210M2_1");
  hCuFrameSxpb210M2_10  = (TH1D*)fFile2->Get("hCuFrameSxpb210M2_10");
  hCuFrameSxth232M2_001 = (TH1D*)fFile2->Get("hCuFrameSxth232M2_001");
  hCuFrameSxth232M2_01  = (TH1D*)fFile2->Get("hCuFrameSxth232M2_01");
  hCuFrameSxth232M2_1   = (TH1D*)fFile2->Get("hCuFrameSxth232M2_1");
  hCuFrameSxth232M2_10  = (TH1D*)fFile2->Get("hCuFrameSxth232M2_10");
  hCuFrameSxu238M2_001  = (TH1D*)fFile2->Get("hCuFrameSxu238M2_001");
  hCuFrameSxu238M2_01   = (TH1D*)fFile2->Get("hCuFrameSxu238M2_01");
  hCuFrameSxu238M2_1    = (TH1D*)fFile2->Get("hCuFrameSxu238M2_1");
  hCuFrameSxu238M2_10   = (TH1D*)fFile2->Get("hCuFrameSxu238M2_10");

  // CuBox M1 and M2
  hCuBoxSth232M1_1    = (TH1D*)fFile2->Get("hCuBoxSth232M1_1");
  hCuBoxSu238M1_1     = (TH1D*)fFile2->Get("hCuBoxSu238M1_1");
  hCuBoxSxpb210M1_001 = (TH1D*)fFile2->Get("hCuBoxSxpb210M1_001");
  hCuBoxSxpb210M1_01  = (TH1D*)fFile2->Get("hCuBoxSxpb210M1_01");
  hCuBoxSxpb210M1_1   = (TH1D*)fFile2->Get("hCuBoxSxpb210M1_1");
  hCuBoxSxpb210M1_10  = (TH1D*)fFile2->Get("hCuBoxSxpb210M1_10");
  hCuBoxSxth232M1_001 = (TH1D*)fFile2->Get("hCuBoxSxth232M1_001");
  hCuBoxSxth232M1_01  = (TH1D*)fFile2->Get("hCuBoxSxth232M1_01");
  hCuBoxSxth232M1_1   = (TH1D*)fFile2->Get("hCuBoxSxth232M1_1");
  hCuBoxSxth232M1_10  = (TH1D*)fFile2->Get("hCuBoxSxth232M1_10");
  hCuBoxSxu238M1_001  = (TH1D*)fFile2->Get("hCuBoxSxu238M1_001");
  hCuBoxSxu238M1_01   = (TH1D*)fFile2->Get("hCuBoxSxu238M1_01");
  hCuBoxSxu238M1_1    = (TH1D*)fFile2->Get("hCuBoxSxu238M1_1");
  hCuBoxSxu238M1_10   = (TH1D*)fFile2->Get("hCuBoxSxu238M1_10");

  hCuBoxSth232M2_1    = (TH1D*)fFile2->Get("hCuBoxSth232M2_1");
  hCuBoxSu238M2_1     = (TH1D*)fFile2->Get("hCuBoxSu238M2_1");
  hCuBoxSxpb210M2_001 = (TH1D*)fFile2->Get("hCuBoxSxpb210M2_001");
  hCuBoxSxpb210M2_01  = (TH1D*)fFile2->Get("hCuBoxSxpb210M2_01");
  hCuBoxSxpb210M2_1   = (TH1D*)fFile2->Get("hCuBoxSxpb210M2_1");
  hCuBoxSxpb210M2_10  = (TH1D*)fFile2->Get("hCuBoxSxpb210M2_10");
  hCuBoxSxth232M2_001 = (TH1D*)fFile2->Get("hCuBoxSxth232M2_001");
  hCuBoxSxth232M2_01  = (TH1D*)fFile2->Get("hCuBoxSxth232M2_01");
  hCuBoxSxth232M2_1   = (TH1D*)fFile2->Get("hCuBoxSxth232M2_1");
  hCuBoxSxth232M2_10  = (TH1D*)fFile2->Get("hCuBoxSxth232M2_10");
  hCuBoxSxu238M2_001  = (TH1D*)fFile2->Get("hCuBoxSxu238M2_001");
  hCuBoxSxu238M2_01   = (TH1D*)fFile2->Get("hCuBoxSxu238M2_01");
  hCuBoxSxu238M2_1    = (TH1D*)fFile2->Get("hCuBoxSxu238M2_1");
  hCuBoxSxu238M2_10   = (TH1D*)fFile2->Get("hCuBoxSxu238M2_10");

//////// Get adaptive binned histograms
  // Crystal M1 and M2
  hTeO20nuM1->Rebin(dAdaptiveBinsM1, "hnewTeO20nuM1", dAdaptiveArrayM1);
  hTeO22nuM1->Rebin(dAdaptiveBinsM1, "hnewTeO22nuM1", dAdaptiveArrayM1);
  hTeO2co60M1->Rebin(dAdaptiveBinsM1, "hnewTeO2co60M1", dAdaptiveArrayM1);
  hTeO2k40M1->Rebin(dAdaptiveBinsM1, "hnewTeO2k40M1", dAdaptiveArrayM1);
  hTeO2pb210M1->Rebin(dAdaptiveBinsM1, "hnewTeO2pb210M1", dAdaptiveArrayM1);
  hTeO2po210M1->Rebin(dAdaptiveBinsM1, "hnewTeO2po210M1", dAdaptiveArrayM1);
  hTeO2te125M1->Rebin(dAdaptiveBinsM1, "hnewTeO2te125M1", dAdaptiveArrayM1);
  hTeO2th232M1->Rebin(dAdaptiveBinsM1, "hnewTeO2th232M1", dAdaptiveArrayM1);
  hTeO2th228M1->Rebin(dAdaptiveBinsM1, "hnewTeO2th228M1", dAdaptiveArrayM1);
  hTeO2ra226M1->Rebin(dAdaptiveBinsM1, "hnewTeO2ra226M1", dAdaptiveArrayM1);
  hTeO2rn222M1->Rebin(dAdaptiveBinsM1, "hnewTeO2rn222M1", dAdaptiveArrayM1);
  hTeO2u238M1->Rebin(dAdaptiveBinsM1, "hnewTeO2u238M1", dAdaptiveArrayM1);
  hTeO2th230M1->Rebin(dAdaptiveBinsM1, "hnewTeO2th230M1", dAdaptiveArrayM1);
  hTeO2u234M1->Rebin(dAdaptiveBinsM1, "hnewTeO2u234M1", dAdaptiveArrayM1);

  hTeO2Spb210M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Spb210M1_01", dAdaptiveArrayM1);
  hTeO2Spo210M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Spo210M1_001", dAdaptiveArrayM1);
  hTeO2Spo210M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Spo210M1_01", dAdaptiveArrayM1);
  hTeO2Sth232M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sth232M1_01", dAdaptiveArrayM1);
  hTeO2Su238M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Su238M1_01", dAdaptiveArrayM1);
  hTeO2Sxpb210M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpb210M1_001", dAdaptiveArrayM1);
  hTeO2Sxpb210M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpb210M1_01", dAdaptiveArrayM1);
  hTeO2Sxpb210M1_1->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpb210M1_1", dAdaptiveArrayM1);
  hTeO2Sxpb210M1_10->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpb210M1_10", dAdaptiveArrayM1);
  hTeO2Sxpo210M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpo210M1_001", dAdaptiveArrayM1);
  hTeO2Sxpo210M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpo210M1_01", dAdaptiveArrayM1);
  hTeO2Sxpo210M1_1->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpo210M1_1", dAdaptiveArrayM1);
  hTeO2Sxth232M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232M1_001", dAdaptiveArrayM1);
  hTeO2Sxth232M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232M1_01", dAdaptiveArrayM1);
  hTeO2Sxth232M1_1->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232M1_1", dAdaptiveArrayM1);
  hTeO2Sxth232M1_10->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232M1_10", dAdaptiveArrayM1);
  hTeO2Sxu238M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238M1_001", dAdaptiveArrayM1);
  hTeO2Sxu238M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238M1_01", dAdaptiveArrayM1);
  hTeO2Sxu238M1_1->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238M1_1", dAdaptiveArrayM1);
  hTeO2Sxu238M1_10->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238M1_10", dAdaptiveArrayM1);

  hTeO20nuM2->Rebin(dAdaptiveBinsM2, "hnewTeO20nuM2", dAdaptiveArrayM2);
  hTeO22nuM2->Rebin(dAdaptiveBinsM2, "hnewTeO22nuM2", dAdaptiveArrayM2);
  hTeO2co60M2->Rebin(dAdaptiveBinsM2, "hnewTeO2co60M2", dAdaptiveArrayM2);
  hTeO2k40M2->Rebin(dAdaptiveBinsM2, "hnewTeO2k40M2", dAdaptiveArrayM2);
  hTeO2pb210M2->Rebin(dAdaptiveBinsM2, "hnewTeO2pb210M2", dAdaptiveArrayM2);
  hTeO2po210M2->Rebin(dAdaptiveBinsM2, "hnewTeO2po210M2", dAdaptiveArrayM2);
  hTeO2te125M2->Rebin(dAdaptiveBinsM2, "hnewTeO2te125M2", dAdaptiveArrayM2);
  hTeO2th232M2->Rebin(dAdaptiveBinsM2, "hnewTeO2th232M2", dAdaptiveArrayM2);
  hTeO2th228M2->Rebin(dAdaptiveBinsM2, "hnewTeO2th228M2", dAdaptiveArrayM2);
  hTeO2ra226M2->Rebin(dAdaptiveBinsM2, "hnewTeO2ra226M2", dAdaptiveArrayM2);
  hTeO2rn222M2->Rebin(dAdaptiveBinsM2, "hnewTeO2rn222M2", dAdaptiveArrayM2);
  hTeO2u238M2->Rebin(dAdaptiveBinsM2, "hnewTeO2u238M2", dAdaptiveArrayM2);
  hTeO2th230M2->Rebin(dAdaptiveBinsM2, "hnewTeO2th230M2", dAdaptiveArrayM2);
  hTeO2u234M2->Rebin(dAdaptiveBinsM2, "hnewTeO2u234M2", dAdaptiveArrayM2);

  hTeO2Spb210M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Spb210M2_01", dAdaptiveArrayM2);
  hTeO2Spo210M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Spo210M2_001", dAdaptiveArrayM2);
  hTeO2Spo210M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Spo210M2_01", dAdaptiveArrayM2);
  hTeO2Sth232M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sth232M2_01", dAdaptiveArrayM2);
  hTeO2Su238M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Su238M2_01", dAdaptiveArrayM2);
  hTeO2Sxpb210M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpb210M2_001", dAdaptiveArrayM2);
  hTeO2Sxpb210M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpb210M2_01", dAdaptiveArrayM2);
  hTeO2Sxpb210M2_1->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpb210M2_1", dAdaptiveArrayM2);
  hTeO2Sxpb210M2_10->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpb210M2_10", dAdaptiveArrayM2);
  hTeO2Sxpo210M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpo210M2_001", dAdaptiveArrayM2);
  hTeO2Sxpo210M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpo210M2_01", dAdaptiveArrayM2);
  hTeO2Sxpo210M2_1->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpo210M2_1", dAdaptiveArrayM2);
  hTeO2Sxth232M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232M2_001", dAdaptiveArrayM2);
  hTeO2Sxth232M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232M2_01", dAdaptiveArrayM2);
  hTeO2Sxth232M2_1->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232M2_1", dAdaptiveArrayM2);
  hTeO2Sxth232M2_10->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232M2_10", dAdaptiveArrayM2);
  hTeO2Sxu238M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238M2_001", dAdaptiveArrayM2);
  hTeO2Sxu238M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238M2_01", dAdaptiveArrayM2);
  hTeO2Sxu238M2_1->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238M2_1", dAdaptiveArrayM2);
  hTeO2Sxu238M2_10->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238M2_10", dAdaptiveArrayM2);

  // Frame M1 and M2
  hCuFrameco58M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameco58M1", dAdaptiveArrayM1);
  hCuFrameco60M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameco60M1", dAdaptiveArrayM1);
  hCuFramecs137M1->Rebin(dAdaptiveBinsM1, "hnewCuFramecs137M1", dAdaptiveArrayM1);
  hCuFramek40M1->Rebin(dAdaptiveBinsM1, "hnewCuFramek40M1", dAdaptiveArrayM1);
  hCuFramemn54M1->Rebin(dAdaptiveBinsM1, "hnewCuFramemn54M1", dAdaptiveArrayM1);
  hCuFramepb210M1->Rebin(dAdaptiveBinsM1, "hnewCuFramepb210M1", dAdaptiveArrayM1);
  hCuFrameth232M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameth232M1", dAdaptiveArrayM1);
  hCuFrameu238M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameu238M1", dAdaptiveArrayM1);

  hCuFrameSth232M1_1->Rebin(dAdaptiveBinsM1, "hnewCuFrameSth232M1_1", dAdaptiveArrayM1);
  hCuFrameSu238M1_1->Rebin(dAdaptiveBinsM1, "hnewCuFrameSu238M1_1", dAdaptiveArrayM1);
  hCuFrameSxpb210M1_001->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxpb210M1_001", dAdaptiveArrayM1);
  hCuFrameSxpb210M1_01->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxpb210M1_01", dAdaptiveArrayM1);
  hCuFrameSxpb210M1_1->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxpb210M1_1", dAdaptiveArrayM1);
  hCuFrameSxpb210M1_10->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxpb210M1_10", dAdaptiveArrayM1);
  hCuFrameSxth232M1_001->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxth232M1_001", dAdaptiveArrayM1);
  hCuFrameSxth232M1_01->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxth232M1_01", dAdaptiveArrayM1);
  hCuFrameSxth232M1_1->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxth232M1_1", dAdaptiveArrayM1);
  hCuFrameSxth232M1_10->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxth232M1_10", dAdaptiveArrayM1);
  hCuFrameSxu238M1_001->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxu238M1_001", dAdaptiveArrayM1);
  hCuFrameSxu238M1_01->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxu238M1_01", dAdaptiveArrayM1);
  hCuFrameSxu238M1_1->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxu238M1_1", dAdaptiveArrayM1);
  hCuFrameSxu238M1_10->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxu238M1_10", dAdaptiveArrayM1);

  hCuFrameco58M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameco58M2", dAdaptiveArrayM2);
  hCuFrameco60M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameco60M2", dAdaptiveArrayM2);
  hCuFramecs137M2->Rebin(dAdaptiveBinsM2, "hnewCuFramecs137M2", dAdaptiveArrayM2);
  hCuFramek40M2->Rebin(dAdaptiveBinsM2, "hnewCuFramek40M2", dAdaptiveArrayM2);
  hCuFramemn54M2->Rebin(dAdaptiveBinsM2, "hnewCuFramemn54M2", dAdaptiveArrayM2);
  hCuFramepb210M2->Rebin(dAdaptiveBinsM2, "hnewCuFramepb210M2", dAdaptiveArrayM2);
  hCuFrameth232M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameth232M2", dAdaptiveArrayM2);
  hCuFrameu238M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameu238M2", dAdaptiveArrayM2);

  hCuFrameSth232M2_1->Rebin(dAdaptiveBinsM2, "hnewCuFrameSth232M2_1", dAdaptiveArrayM2);
  hCuFrameSu238M2_1->Rebin(dAdaptiveBinsM2, "hnewCuFrameSu238M2_1", dAdaptiveArrayM2);
  hCuFrameSxpb210M2_001->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxpb210M2_001", dAdaptiveArrayM2);
  hCuFrameSxpb210M2_01->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxpb210M2_01", dAdaptiveArrayM2);
  hCuFrameSxpb210M2_1->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxpb210M2_1", dAdaptiveArrayM2);
  hCuFrameSxpb210M2_10->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxpb210M2_10", dAdaptiveArrayM2);
  hCuFrameSxth232M2_001->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxth232M2_001", dAdaptiveArrayM2);
  hCuFrameSxth232M2_01->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxth232M2_01", dAdaptiveArrayM2);
  hCuFrameSxth232M2_1->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxth232M2_1", dAdaptiveArrayM2);
  hCuFrameSxth232M2_10->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxth232M2_10", dAdaptiveArrayM2);
  hCuFrameSxu238M2_001->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxu238M2_001", dAdaptiveArrayM2);
  hCuFrameSxu238M2_01->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxu238M2_01", dAdaptiveArrayM2);
  hCuFrameSxu238M2_1->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxu238M2_1", dAdaptiveArrayM2);
  hCuFrameSxu238M2_10->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxu238M2_10", dAdaptiveArrayM2);


  // CuBox (TShield) M1 and M2
  hCuBoxco58M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxco58M1", dAdaptiveArrayM1);
  hCuBoxco60M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxco60M1", dAdaptiveArrayM1);
  hCuBoxcs137M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxcs137M1", dAdaptiveArrayM1);
  hCuBoxk40M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxk40M1", dAdaptiveArrayM1);
  hCuBoxmn54M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxmn54M1", dAdaptiveArrayM1);
  hCuBoxpb210M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxpb210M1", dAdaptiveArrayM1);
  hCuBoxth232M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxth232M1", dAdaptiveArrayM1);
  hCuBoxu238M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxu238M1", dAdaptiveArrayM1);

  hCuBoxSth232M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBoxSth232M1_1", dAdaptiveArrayM1);
  hCuBoxSu238M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBoxSu238M1_1", dAdaptiveArrayM1);
  hCuBoxSxpb210M1_001->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxpb210M1_001", dAdaptiveArrayM1);
  hCuBoxSxpb210M1_01->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxpb210M1_01", dAdaptiveArrayM1);
  hCuBoxSxpb210M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxpb210M1_1", dAdaptiveArrayM1);
  hCuBoxSxpb210M1_10->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxpb210M1_10", dAdaptiveArrayM1);
  hCuBoxSxth232M1_001->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxth232M1_001", dAdaptiveArrayM1);
  hCuBoxSxth232M1_01->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxth232M1_01", dAdaptiveArrayM1);
  hCuBoxSxth232M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxth232M1_1", dAdaptiveArrayM1);
  hCuBoxSxth232M1_10->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxth232M1_10", dAdaptiveArrayM1);
  hCuBoxSxu238M1_001->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxu238M1_001", dAdaptiveArrayM1);
  hCuBoxSxu238M1_01->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxu238M1_01", dAdaptiveArrayM1);
  hCuBoxSxu238M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxu238M1_1", dAdaptiveArrayM1);
  hCuBoxSxu238M1_10->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxu238M1_10", dAdaptiveArrayM1);

  hCuBoxco58M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxco58M2", dAdaptiveArrayM2);
  hCuBoxco60M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxco60M2", dAdaptiveArrayM2);
  hCuBoxcs137M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxcs137M2", dAdaptiveArrayM2);
  hCuBoxk40M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxk40M2", dAdaptiveArrayM2);
  hCuBoxmn54M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxmn54M2", dAdaptiveArrayM2);
  hCuBoxpb210M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxpb210M2", dAdaptiveArrayM2);
  hCuBoxth232M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxth232M2", dAdaptiveArrayM2);
  hCuBoxu238M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxu238M2", dAdaptiveArrayM2);

  hCuBoxSth232M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBoxSth232M2_1", dAdaptiveArrayM2);
  hCuBoxSu238M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBoxSu238M2_1", dAdaptiveArrayM2);
  hCuBoxSxpb210M2_001->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxpb210M2_001", dAdaptiveArrayM2);
  hCuBoxSxpb210M2_01->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxpb210M2_01", dAdaptiveArrayM2);
  hCuBoxSxpb210M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxpb210M2_1", dAdaptiveArrayM2);
  hCuBoxSxpb210M2_10->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxpb210M2_10", dAdaptiveArrayM2);
  hCuBoxSxth232M2_001->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxth232M2_001", dAdaptiveArrayM2);
  hCuBoxSxth232M2_01->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxth232M2_01", dAdaptiveArrayM2);
  hCuBoxSxth232M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxth232M2_1", dAdaptiveArrayM2);
  hCuBoxSxth232M2_10->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxth232M2_10", dAdaptiveArrayM2);
  hCuBoxSxu238M2_001->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxu238M2_001", dAdaptiveArrayM2);
  hCuBoxSxu238M2_01->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxu238M2_01", dAdaptiveArrayM2);
  hCuBoxSxu238M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxu238M2_1", dAdaptiveArrayM2);
  hCuBoxSxu238M2_10->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxu238M2_10", dAdaptiveArrayM2);


  // 50mK M1 and M2
  h50mKco58M1->Rebin(dAdaptiveBinsM1, "hnew50mKco58M1", dAdaptiveArrayM1);
  h50mKco60M1->Rebin(dAdaptiveBinsM1, "hnew50mKco60M1", dAdaptiveArrayM1);
  h50mKcs137M1->Rebin(dAdaptiveBinsM1, "hnew50mKcs137M1", dAdaptiveArrayM1);
  h50mKk40M1->Rebin(dAdaptiveBinsM1, "hnew50mKk40M1", dAdaptiveArrayM1);
  h50mKmn54M1->Rebin(dAdaptiveBinsM1, "hnew50mKmn54M1", dAdaptiveArrayM1);
  h50mKpb210M1->Rebin(dAdaptiveBinsM1, "hnew50mKpb210M1", dAdaptiveArrayM1);
  h50mKth232M1->Rebin(dAdaptiveBinsM1, "hnew50mKth232M1", dAdaptiveArrayM1);
  h50mKu238M1->Rebin(dAdaptiveBinsM1, "hnew50mKu238M1", dAdaptiveArrayM1);

  h50mKco58M2->Rebin(dAdaptiveBinsM2, "hnew50mKco58M2", dAdaptiveArrayM2);
  h50mKco60M2->Rebin(dAdaptiveBinsM2, "hnew50mKco60M2", dAdaptiveArrayM2);
  h50mKcs137M2->Rebin(dAdaptiveBinsM2, "hnew50mKcs137M2", dAdaptiveArrayM2);
  h50mKk40M2->Rebin(dAdaptiveBinsM2, "hnew50mKk40M2", dAdaptiveArrayM2);
  h50mKmn54M2->Rebin(dAdaptiveBinsM2, "hnew50mKmn54M2", dAdaptiveArrayM2);
  h50mKpb210M2->Rebin(dAdaptiveBinsM2, "hnew50mKpb210M2", dAdaptiveArrayM2);
  h50mKth232M2->Rebin(dAdaptiveBinsM2, "hnew50mKth232M2", dAdaptiveArrayM2);
  h50mKu238M2->Rebin(dAdaptiveBinsM2, "hnew50mKu238M2", dAdaptiveArrayM2);

  // 600mK
  h600mKco60M1->Rebin(dAdaptiveBinsM1, "hnew600mKco60M1", dAdaptiveArrayM1);
  h600mKk40M1->Rebin(dAdaptiveBinsM1, "hnew600mKk40M1", dAdaptiveArrayM1);
  h600mKth232M1->Rebin(dAdaptiveBinsM1, "hnew600mKth232M1", dAdaptiveArrayM1);
  h600mKu238M1->Rebin(dAdaptiveBinsM1, "hnew600mKu238M1", dAdaptiveArrayM1);

  h600mKco60M2->Rebin(dAdaptiveBinsM2, "hnew600mKco60M2", dAdaptiveArrayM2);
  h600mKk40M2->Rebin(dAdaptiveBinsM2, "hnew600mKk40M2", dAdaptiveArrayM2);
  h600mKth232M2->Rebin(dAdaptiveBinsM2, "hnew600mKth232M2", dAdaptiveArrayM2);
  h600mKu238M2->Rebin(dAdaptiveBinsM2, "hnew600mKu238M2", dAdaptiveArrayM2);

  // Roman Lead M1 and M2
  hPbRombi207M1->Rebin(dAdaptiveBinsM1, "hnewPbRombi207M1", dAdaptiveArrayM1);
  hPbRomco60M1->Rebin(dAdaptiveBinsM1, "hnewPbRomco60M1", dAdaptiveArrayM1);
  hPbRomcs137M1->Rebin(dAdaptiveBinsM1, "hnewPbRomcs137M1", dAdaptiveArrayM1);
  hPbRomk40M1->Rebin(dAdaptiveBinsM1, "hnewPbRomk40M1", dAdaptiveArrayM1);
  hPbRompb210M1->Rebin(dAdaptiveBinsM1, "hnewPbRompb210M1", dAdaptiveArrayM1);
  hPbRomth232M1->Rebin(dAdaptiveBinsM1, "hnewPbRomth232M1", dAdaptiveArrayM1);
  hPbRomu238M1->Rebin(dAdaptiveBinsM1, "hnewPbRomu238M1", dAdaptiveArrayM1);

  hPbRombi207M2->Rebin(dAdaptiveBinsM2, "hnewPbRombi207M2", dAdaptiveArrayM2);
  hPbRomco60M2->Rebin(dAdaptiveBinsM2, "hnewPbRomco60M2", dAdaptiveArrayM2);
  hPbRomcs137M2->Rebin(dAdaptiveBinsM2, "hnewPbRomcs137M2", dAdaptiveArrayM2);
  hPbRomk40M2->Rebin(dAdaptiveBinsM2, "hnewPbRomk40M2", dAdaptiveArrayM2);
  hPbRompb210M2->Rebin(dAdaptiveBinsM2, "hnewPbRompb210M2", dAdaptiveArrayM2);
  hPbRomth232M2->Rebin(dAdaptiveBinsM2, "hnewPbRomth232M2", dAdaptiveArrayM2);
  hPbRomu238M2->Rebin(dAdaptiveBinsM2, "hnewPbRomu238M2", dAdaptiveArrayM2);

  // Main Bath M1 and M2
  hMBco60M1->Rebin(dAdaptiveBinsM1, "hnewMBco60M1", dAdaptiveArrayM1);
  hMBk40M1->Rebin(dAdaptiveBinsM1, "hnewMBk40M1", dAdaptiveArrayM1);
  hMBth232M1->Rebin(dAdaptiveBinsM1, "hnewMBth232M1", dAdaptiveArrayM1);
  hMBu238M1->Rebin(dAdaptiveBinsM1, "hnewMBu238M1", dAdaptiveArrayM1);

  hMBco60M2->Rebin(dAdaptiveBinsM2, "hnewMBco60M2", dAdaptiveArrayM2);
  hMBk40M2->Rebin(dAdaptiveBinsM2, "hnewMBk40M2", dAdaptiveArrayM2);
  hMBth232M2->Rebin(dAdaptiveBinsM2, "hnewMBth232M2", dAdaptiveArrayM2);
  hMBu238M2->Rebin(dAdaptiveBinsM2, "hnewMBu238M2", dAdaptiveArrayM2);

  // Super Insulation M1 and M2
  hSIk40M1->Rebin(dAdaptiveBinsM1, "hnewSIk40M1", dAdaptiveArrayM1);
  hSIth232M1->Rebin(dAdaptiveBinsM1, "hnewSIth232M1", dAdaptiveArrayM1);
  hSIu238M1->Rebin(dAdaptiveBinsM1, "hnewSIu238M1", dAdaptiveArrayM1);

  hSIk40M2->Rebin(dAdaptiveBinsM2, "hnewSIk40M2", dAdaptiveArrayM2);
  hSIth232M2->Rebin(dAdaptiveBinsM2, "hnewSIth232M2", dAdaptiveArrayM2);
  hSIu238M2->Rebin(dAdaptiveBinsM2, "hnewSIu238M2", dAdaptiveArrayM2);

  // IVC M1 and M2
  hIVCco60M1->Rebin(dAdaptiveBinsM1, "hnewIVCco60M1", dAdaptiveArrayM1);
  hIVCk40M1->Rebin(dAdaptiveBinsM1, "hnewIVCk40M1", dAdaptiveArrayM1);
  hIVCth232M1->Rebin(dAdaptiveBinsM1, "hnewIVCth232M1", dAdaptiveArrayM1);
  hIVCu238M1->Rebin(dAdaptiveBinsM1, "hnewIVCu238M1", dAdaptiveArrayM1);

  hIVCco60M2->Rebin(dAdaptiveBinsM2, "hnewIVCco60M2", dAdaptiveArrayM2);
  hIVCk40M2->Rebin(dAdaptiveBinsM2, "hnewIVCk40M2", dAdaptiveArrayM2);
  hIVCth232M2->Rebin(dAdaptiveBinsM2, "hnewIVCth232M2", dAdaptiveArrayM2);
  hIVCu238M2->Rebin(dAdaptiveBinsM2, "hnewIVCu238M2", dAdaptiveArrayM2);

  // OVC M1 and M2
  hOVCco60M1->Rebin(dAdaptiveBinsM1, "hnewOVCco60M1", dAdaptiveArrayM1);
  hOVCk40M1->Rebin(dAdaptiveBinsM1, "hnewOVCk40M1", dAdaptiveArrayM1);
  hOVCth232M1->Rebin(dAdaptiveBinsM1, "hnewOVCth232M1", dAdaptiveArrayM1);
  hOVCu238M1->Rebin(dAdaptiveBinsM1, "hnewOVCu238M1", dAdaptiveArrayM1);

  hOVCco60M2->Rebin(dAdaptiveBinsM2, "hnewOVCco60M2", dAdaptiveArrayM2);
  hOVCk40M2->Rebin(dAdaptiveBinsM2, "hnewOVCk40M2", dAdaptiveArrayM2);
  hOVCth232M2->Rebin(dAdaptiveBinsM2, "hnewOVCth232M2", dAdaptiveArrayM2);
  hOVCu238M2->Rebin(dAdaptiveBinsM2, "hnewOVCu238M2", dAdaptiveArrayM2);


  // Fill adaptive binning histograms
  for(int i = 1; i <= dAdaptiveBinsM1; i++)
  {
    hAdapTeO20nuM1->SetBinContent(i, dBinSize * hnewTeO20nuM1->GetBinContent(i)/hnewTeO20nuM1->GetBinWidth(i));
    hAdapTeO22nuM1->SetBinContent(i, dBinSize * hnewTeO22nuM1->GetBinContent(i)/hnewTeO22nuM1->GetBinWidth(i));
    hAdapTeO2co60M1->SetBinContent(i, dBinSize * hnewTeO2co60M1->GetBinContent(i)/hnewTeO2co60M1->GetBinWidth(i));
    hAdapTeO2k40M1->SetBinContent(i, dBinSize * hnewTeO2k40M1->GetBinContent(i)/hnewTeO2k40M1->GetBinWidth(i));
    hAdapTeO2pb210M1->SetBinContent(i, dBinSize * hnewTeO2pb210M1->GetBinContent(i)/hnewTeO2pb210M1->GetBinWidth(i));
    hAdapTeO2po210M1->SetBinContent(i, dBinSize * hnewTeO2po210M1->GetBinContent(i)/hnewTeO2po210M1->GetBinWidth(i));
    hAdapTeO2te125M1->SetBinContent(i, dBinSize * hnewTeO2te125M1->GetBinContent(i)/hnewTeO2te125M1->GetBinWidth(i));
    hAdapTeO2th232M1->SetBinContent(i, dBinSize * hnewTeO2th232M1->GetBinContent(i)/hnewTeO2th232M1->GetBinWidth(i));
    hAdapTeO2th228M1->SetBinContent(i, dBinSize * hnewTeO2th228M1->GetBinContent(i)/hnewTeO2th228M1->GetBinWidth(i));
    hAdapTeO2ra226M1->SetBinContent(i, dBinSize * hnewTeO2ra226M1->GetBinContent(i)/hnewTeO2ra226M1->GetBinWidth(i));
    hAdapTeO2rn222M1->SetBinContent(i, dBinSize * hnewTeO2rn222M1->GetBinContent(i)/hnewTeO2rn222M1->GetBinWidth(i));
    hAdapTeO2u238M1->SetBinContent(i, dBinSize * hnewTeO2u238M1->GetBinContent(i)/hnewTeO2u238M1->GetBinWidth(i));
    hAdapTeO2th230M1->SetBinContent(i, dBinSize * hnewTeO2th230M1->GetBinContent(i)/hnewTeO2th230M1->GetBinWidth(i));
    hAdapTeO2u234M1->SetBinContent(i, dBinSize * hnewTeO2u234M1->GetBinContent(i)/hnewTeO2u234M1->GetBinWidth(i));

    hAdapTeO2Spb210M1_01->SetBinContent(i, dBinSize * hnewTeO2Spb210M1_01->GetBinContent(i)/hnewTeO2Spb210M1_01->GetBinWidth(i));
    hAdapTeO2Spo210M1_001->SetBinContent(i, dBinSize * hnewTeO2Spo210M1_001->GetBinContent(i)/hnewTeO2Spo210M1_001->GetBinWidth(i));
    hAdapTeO2Spo210M1_01->SetBinContent(i, dBinSize * hnewTeO2Spo210M1_01->GetBinContent(i)/hnewTeO2Spo210M1_01->GetBinWidth(i));
    hAdapTeO2Sth232M1_01->SetBinContent(i, dBinSize * hnewTeO2Sth232M1_01->GetBinContent(i)/hnewTeO2Sth232M1_01->GetBinWidth(i));
    hAdapTeO2Su238M1_01->SetBinContent(i, dBinSize * hnewTeO2Su238M1_01->GetBinContent(i)/hnewTeO2Su238M1_01->GetBinWidth(i));
    hAdapTeO2Sxpb210M1_001->SetBinContent(i, dBinSize * hnewTeO2Sxpb210M1_001->GetBinContent(i)/hnewTeO2Sxpb210M1_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M1_01->SetBinContent(i, dBinSize * hnewTeO2Sxpb210M1_01->GetBinContent(i)/hnewTeO2Sxpb210M1_01->GetBinWidth(i));
    hAdapTeO2Sxpb210M1_1->SetBinContent(i, dBinSize * hnewTeO2Sxpb210M1_1->GetBinContent(i)/hnewTeO2Sxpb210M1_1->GetBinWidth(i));
    hAdapTeO2Sxpb210M1_10->SetBinContent(i, dBinSize * hnewTeO2Sxpb210M1_10->GetBinContent(i)/hnewTeO2Sxpb210M1_10->GetBinWidth(i));
    hAdapTeO2Sxpo210M1_001->SetBinContent(i, dBinSize * hnewTeO2Sxpo210M1_001->GetBinContent(i)/hnewTeO2Sxpo210M1_001->GetBinWidth(i));
    hAdapTeO2Sxpo210M1_01->SetBinContent(i, dBinSize * hnewTeO2Sxpo210M1_01->GetBinContent(i)/hnewTeO2Sxpo210M1_01->GetBinWidth(i));
    hAdapTeO2Sxpo210M1_1->SetBinContent(i, dBinSize * hnewTeO2Sxpo210M1_1->GetBinContent(i)/hnewTeO2Sxpo210M1_1->GetBinWidth(i));
    hAdapTeO2Sxth232M1_001->SetBinContent(i, dBinSize * hnewTeO2Sxth232M1_001->GetBinContent(i)/hnewTeO2Sxth232M1_001->GetBinWidth(i));
    hAdapTeO2Sxth232M1_01->SetBinContent(i, dBinSize * hnewTeO2Sxth232M1_01->GetBinContent(i)/hnewTeO2Sxth232M1_01->GetBinWidth(i));
    hAdapTeO2Sxth232M1_1->SetBinContent(i, dBinSize * hnewTeO2Sxth232M1_1->GetBinContent(i)/hnewTeO2Sxth232M1_1->GetBinWidth(i));
    hAdapTeO2Sxth232M1_10->SetBinContent(i, dBinSize * hnewTeO2Sxth232M1_10->GetBinContent(i)/hnewTeO2Sxth232M1_10->GetBinWidth(i));
    hAdapTeO2Sxu238M1_001->SetBinContent(i, dBinSize * hnewTeO2Sxu238M1_001->GetBinContent(i)/hnewTeO2Sxu238M1_001->GetBinWidth(i));
    hAdapTeO2Sxu238M1_01->SetBinContent(i, dBinSize * hnewTeO2Sxu238M1_01->GetBinContent(i)/hnewTeO2Sxu238M1_01->GetBinWidth(i));
    hAdapTeO2Sxu238M1_1->SetBinContent(i, dBinSize * hnewTeO2Sxu238M1_1->GetBinContent(i)/hnewTeO2Sxu238M1_1->GetBinWidth(i));
    hAdapTeO2Sxu238M1_10->SetBinContent(i, dBinSize * hnewTeO2Sxu238M1_10->GetBinContent(i)/hnewTeO2Sxu238M1_10->GetBinWidth(i));

    hAdapCuFrameco58M1->SetBinContent(i, dBinSize * hnewCuFrameco58M1->GetBinContent(i)/hnewCuFrameco58M1->GetBinWidth(i));
    hAdapCuFrameco60M1->SetBinContent(i, dBinSize * hnewCuFrameco60M1->GetBinContent(i)/hnewCuFrameco60M1->GetBinWidth(i));
    hAdapCuFramecs137M1->SetBinContent(i, dBinSize * hnewCuFramecs137M1->GetBinContent(i)/hnewCuFramecs137M1->GetBinWidth(i));
    hAdapCuFramek40M1->SetBinContent(i, dBinSize * hnewCuFramek40M1->GetBinContent(i)/hnewCuFramek40M1->GetBinWidth(i));
    hAdapCuFramemn54M1->SetBinContent(i, dBinSize * hnewCuFramemn54M1->GetBinContent(i)/hnewCuFramemn54M1->GetBinWidth(i));
    hAdapCuFramepb210M1->SetBinContent(i, dBinSize * hnewCuFramepb210M1->GetBinContent(i)/hnewCuFramepb210M1->GetBinWidth(i));
    hAdapCuFrameth232M1->SetBinContent(i, dBinSize * hnewCuFrameth232M1->GetBinContent(i)/hnewCuFrameth232M1->GetBinWidth(i));
    hAdapCuFrameu238M1->SetBinContent(i, dBinSize * hnewCuFrameu238M1->GetBinContent(i)/hnewCuFrameu238M1->GetBinWidth(i));

    hAdapCuFrameSth232M1_1->SetBinContent(i, dBinSize * hnewCuFrameSth232M1_1->GetBinContent(i)/hnewCuFrameSth232M1_1->GetBinWidth(i));
    hAdapCuFrameSu238M1_1->SetBinContent(i, dBinSize * hnewCuFrameSu238M1_1->GetBinContent(i)/hnewCuFrameSu238M1_1->GetBinWidth(i));
    hAdapCuFrameSxpb210M1_001->SetBinContent(i, dBinSize * hnewCuFrameSxpb210M1_001->GetBinContent(i)/hnewCuFrameSxpb210M1_001->GetBinWidth(i));
    hAdapCuFrameSxpb210M1_01->SetBinContent(i, dBinSize * hnewCuFrameSxpb210M1_01->GetBinContent(i)/hnewCuFrameSxpb210M1_01->GetBinWidth(i));
    hAdapCuFrameSxpb210M1_1->SetBinContent(i, dBinSize * hnewCuFrameSxpb210M1_1->GetBinContent(i)/hnewCuFrameSxpb210M1_1->GetBinWidth(i));
    hAdapCuFrameSxpb210M1_10->SetBinContent(i, dBinSize * hnewCuFrameSxpb210M1_10->GetBinContent(i)/hnewCuFrameSxpb210M1_10->GetBinWidth(i));
    hAdapCuFrameSxth232M1_001->SetBinContent(i, dBinSize * hnewCuFrameSxth232M1_001->GetBinContent(i)/hnewCuFrameSxth232M1_001->GetBinWidth(i));
    hAdapCuFrameSxth232M1_01->SetBinContent(i, dBinSize * hnewCuFrameSxth232M1_01->GetBinContent(i)/hnewCuFrameSxth232M1_01->GetBinWidth(i));
    hAdapCuFrameSxth232M1_1->SetBinContent(i, dBinSize * hnewCuFrameSxth232M1_1->GetBinContent(i)/hnewCuFrameSxth232M1_1->GetBinWidth(i));
    hAdapCuFrameSxth232M1_10->SetBinContent(i, dBinSize * hnewCuFrameSxth232M1_10->GetBinContent(i)/hnewCuFrameSxth232M1_10->GetBinWidth(i));
    hAdapCuFrameSxu238M1_001->SetBinContent(i, dBinSize * hnewCuFrameSxu238M1_001->GetBinContent(i)/hnewCuFrameSxu238M1_001->GetBinWidth(i));
    hAdapCuFrameSxu238M1_01->SetBinContent(i, dBinSize * hnewCuFrameSxu238M1_01->GetBinContent(i)/hnewCuFrameSxu238M1_01->GetBinWidth(i));
    hAdapCuFrameSxu238M1_1->SetBinContent(i, dBinSize * hnewCuFrameSxu238M1_1->GetBinContent(i)/hnewCuFrameSxu238M1_1->GetBinWidth(i));
    hAdapCuFrameSxu238M1_10->SetBinContent(i, dBinSize * hnewCuFrameSxu238M1_10->GetBinContent(i)/hnewCuFrameSxu238M1_10->GetBinWidth(i));

    hAdapCuBoxco58M1->SetBinContent(i, dBinSize * hnewCuBoxco58M1->GetBinContent(i)/hnewCuBoxco58M1->GetBinWidth(i));
    hAdapCuBoxco60M1->SetBinContent(i, dBinSize * hnewCuBoxco60M1->GetBinContent(i)/hnewCuBoxco60M1->GetBinWidth(i));
    hAdapCuBoxcs137M1->SetBinContent(i, dBinSize * hnewCuBoxcs137M1->GetBinContent(i)/hnewCuBoxcs137M1->GetBinWidth(i));
    hAdapCuBoxk40M1->SetBinContent(i, dBinSize * hnewCuBoxk40M1->GetBinContent(i)/hnewCuBoxk40M1->GetBinWidth(i));
    hAdapCuBoxmn54M1->SetBinContent(i, dBinSize * hnewCuBoxmn54M1->GetBinContent(i)/hnewCuBoxmn54M1->GetBinWidth(i));
    hAdapCuBoxpb210M1->SetBinContent(i, dBinSize * hnewCuBoxpb210M1->GetBinContent(i)/hnewCuBoxpb210M1->GetBinWidth(i));
    hAdapCuBoxth232M1->SetBinContent(i, dBinSize * hnewCuBoxth232M1->GetBinContent(i)/hnewCuBoxth232M1->GetBinWidth(i));
    hAdapCuBoxu238M1->SetBinContent(i, dBinSize * hnewCuBoxu238M1->GetBinContent(i)/hnewCuBoxu238M1->GetBinWidth(i));

    hAdapCuBoxSth232M1_1->SetBinContent(i, dBinSize * hnewCuBoxSth232M1_1->GetBinContent(i)/hnewCuBoxSth232M1_1->GetBinWidth(i));
    hAdapCuBoxSu238M1_1->SetBinContent(i, dBinSize * hnewCuBoxSu238M1_1->GetBinContent(i)/hnewCuBoxSu238M1_1->GetBinWidth(i));
    hAdapCuBoxSxpb210M1_001->SetBinContent(i, dBinSize * hnewCuBoxSxpb210M1_001->GetBinContent(i)/hnewCuBoxSxpb210M1_001->GetBinWidth(i));
    hAdapCuBoxSxpb210M1_01->SetBinContent(i, dBinSize * hnewCuBoxSxpb210M1_01->GetBinContent(i)/hnewCuBoxSxpb210M1_01->GetBinWidth(i));
    hAdapCuBoxSxpb210M1_1->SetBinContent(i, dBinSize * hnewCuBoxSxpb210M1_1->GetBinContent(i)/hnewCuBoxSxpb210M1_1->GetBinWidth(i));
    hAdapCuBoxSxpb210M1_10->SetBinContent(i, dBinSize * hnewCuBoxSxpb210M1_10->GetBinContent(i)/hnewCuBoxSxpb210M1_10->GetBinWidth(i));
    hAdapCuBoxSxth232M1_001->SetBinContent(i, dBinSize * hnewCuBoxSxth232M1_001->GetBinContent(i)/hnewCuBoxSxth232M1_001->GetBinWidth(i));
    hAdapCuBoxSxth232M1_01->SetBinContent(i, dBinSize * hnewCuBoxSxth232M1_01->GetBinContent(i)/hnewCuBoxSxth232M1_01->GetBinWidth(i));
    hAdapCuBoxSxth232M1_1->SetBinContent(i, dBinSize * hnewCuBoxSxth232M1_1->GetBinContent(i)/hnewCuBoxSxth232M1_1->GetBinWidth(i));
    hAdapCuBoxSxth232M1_10->SetBinContent(i, dBinSize * hnewCuBoxSxth232M1_10->GetBinContent(i)/hnewCuBoxSxth232M1_10->GetBinWidth(i));
    hAdapCuBoxSxu238M1_001->SetBinContent(i, dBinSize * hnewCuBoxSxu238M1_001->GetBinContent(i)/hnewCuBoxSxu238M1_001->GetBinWidth(i));
    hAdapCuBoxSxu238M1_01->SetBinContent(i, dBinSize * hnewCuBoxSxu238M1_01->GetBinContent(i)/hnewCuBoxSxu238M1_01->GetBinWidth(i));
    hAdapCuBoxSxu238M1_1->SetBinContent(i, dBinSize * hnewCuBoxSxu238M1_1->GetBinContent(i)/hnewCuBoxSxu238M1_1->GetBinWidth(i));
    hAdapCuBoxSxu238M1_10->SetBinContent(i, dBinSize * hnewCuBoxSxu238M1_10->GetBinContent(i)/hnewCuBoxSxu238M1_10->GetBinWidth(i));

    hAdap50mKco58M1->SetBinContent(i, dBinSize * hnew50mKco58M1->GetBinContent(i)/hnew50mKco58M1->GetBinWidth(i));
    hAdap50mKco60M1->SetBinContent(i, dBinSize * hnew50mKco60M1->GetBinContent(i)/hnew50mKco60M1->GetBinWidth(i));
    hAdap50mKcs137M1->SetBinContent(i, dBinSize * hnew50mKcs137M1->GetBinContent(i)/hnew50mKcs137M1->GetBinWidth(i));
    hAdap50mKk40M1->SetBinContent(i, dBinSize * hnew50mKk40M1->GetBinContent(i)/hnew50mKk40M1->GetBinWidth(i));
    hAdap50mKmn54M1->SetBinContent(i, dBinSize * hnew50mKmn54M1->GetBinContent(i)/hnew50mKmn54M1->GetBinWidth(i));
    hAdap50mKpb210M1->SetBinContent(i, dBinSize * hnew50mKpb210M1->GetBinContent(i)/hnew50mKpb210M1->GetBinWidth(i));
    hAdap50mKth232M1->SetBinContent(i, dBinSize * hnew50mKth232M1->GetBinContent(i)/hnew50mKth232M1->GetBinWidth(i));
    hAdap50mKu238M1->SetBinContent(i, dBinSize * hnew50mKu238M1->GetBinContent(i)/hnew50mKu238M1->GetBinWidth(i));

    hAdap600mKco60M1->SetBinContent(i, dBinSize * hnew600mKco60M1->GetBinContent(i)/hnew600mKco60M1->GetBinWidth(i));
    hAdap600mKk40M1->SetBinContent(i, dBinSize * hnew600mKk40M1->GetBinContent(i)/hnew600mKk40M1->GetBinWidth(i));
    hAdap600mKth232M1->SetBinContent(i, dBinSize * hnew600mKth232M1->GetBinContent(i)/hnew600mKth232M1->GetBinWidth(i));
    hAdap600mKu238M1->SetBinContent(i, dBinSize * hnew600mKu238M1->GetBinContent(i)/hnew600mKu238M1->GetBinWidth(i));

    hAdapPbRombi207M1->SetBinContent(i, dBinSize * hnewPbRombi207M1->GetBinContent(i)/hnewPbRombi207M1->GetBinWidth(i));
    hAdapPbRomco60M1->SetBinContent(i, dBinSize * hnewPbRomco60M1->GetBinContent(i)/hnewPbRomco60M1->GetBinWidth(i));
    hAdapPbRomcs137M1->SetBinContent(i, dBinSize * hnewPbRomcs137M1->GetBinContent(i)/hnewPbRomcs137M1->GetBinWidth(i));
    hAdapPbRomk40M1->SetBinContent(i, dBinSize * hnewPbRomk40M1->GetBinContent(i)/hnewPbRomk40M1->GetBinWidth(i));
    hAdapPbRompb210M1->SetBinContent(i, dBinSize * hnewPbRompb210M1->GetBinContent(i)/hnewPbRompb210M1->GetBinWidth(i));
    hAdapPbRomth232M1->SetBinContent(i, dBinSize * hnewPbRomth232M1->GetBinContent(i)/hnewPbRomth232M1->GetBinWidth(i));
    hAdapPbRomu238M1->SetBinContent(i, dBinSize * hnewPbRomu238M1->GetBinContent(i)/hnewPbRomu238M1->GetBinWidth(i));

    hAdapMBco60M1->SetBinContent(i, dBinSize * hnewMBco60M1->GetBinContent(i)/hnewMBco60M1->GetBinWidth(i));
    hAdapMBk40M1->SetBinContent(i, dBinSize * hnewMBk40M1->GetBinContent(i)/hnewMBk40M1->GetBinWidth(i));
    hAdapMBth232M1->SetBinContent(i, dBinSize * hnewMBth232M1->GetBinContent(i)/hnewMBth232M1->GetBinWidth(i));
    hAdapMBu238M1->SetBinContent(i, dBinSize * hnewMBu238M1->GetBinContent(i)/hnewMBu238M1->GetBinWidth(i));

    hAdapSIk40M1->SetBinContent(i, dBinSize * hnewSIk40M1->GetBinContent(i)/hnewSIk40M1->GetBinWidth(i));
    hAdapSIth232M1->SetBinContent(i, dBinSize * hnewSIth232M1->GetBinContent(i)/hnewSIth232M1->GetBinWidth(i));
    hAdapSIu238M1->SetBinContent(i, dBinSize * hnewSIu238M1->GetBinContent(i)/hnewSIu238M1->GetBinWidth(i));

    hAdapIVCco60M1->SetBinContent(i, dBinSize * hnewIVCco60M1->GetBinContent(i)/hnewIVCco60M1->GetBinWidth(i));
    hAdapIVCk40M1->SetBinContent(i, dBinSize * hnewIVCk40M1->GetBinContent(i)/hnewIVCk40M1->GetBinWidth(i));
    hAdapIVCth232M1->SetBinContent(i, dBinSize * hnewIVCth232M1->GetBinContent(i)/hnewIVCth232M1->GetBinWidth(i));
    hAdapIVCu238M1->SetBinContent(i, dBinSize * hnewIVCu238M1->GetBinContent(i)/hnewIVCu238M1->GetBinWidth(i));

    hAdapOVCco60M1->SetBinContent(i, dBinSize * hnewOVCco60M1->GetBinContent(i)/hnewOVCco60M1->GetBinWidth(i));
    hAdapOVCk40M1->SetBinContent(i, dBinSize * hnewOVCk40M1->GetBinContent(i)/hnewOVCk40M1->GetBinWidth(i));
    hAdapOVCth232M1->SetBinContent(i, dBinSize * hnewOVCth232M1->GetBinContent(i)/hnewOVCth232M1->GetBinWidth(i));
    hAdapOVCu238M1->SetBinContent(i, dBinSize * hnewOVCu238M1->GetBinContent(i)/hnewOVCu238M1->GetBinWidth(i));

  }

  for(int i = 1; i <= dAdaptiveBinsM2; i++)
  {
    hAdapTeO20nuM2->SetBinContent(i, dBinSize * hnewTeO20nuM2->GetBinContent(i)/hnewTeO20nuM2->GetBinWidth(i));
    hAdapTeO22nuM2->SetBinContent(i, dBinSize * hnewTeO22nuM2->GetBinContent(i)/hnewTeO22nuM2->GetBinWidth(i));
    hAdapTeO2co60M2->SetBinContent(i, dBinSize * hnewTeO2co60M2->GetBinContent(i)/hnewTeO2co60M2->GetBinWidth(i));
    hAdapTeO2k40M2->SetBinContent(i, dBinSize * hnewTeO2k40M2->GetBinContent(i)/hnewTeO2k40M2->GetBinWidth(i));
    hAdapTeO2pb210M2->SetBinContent(i, dBinSize * hnewTeO2pb210M2->GetBinContent(i)/hnewTeO2pb210M2->GetBinWidth(i));
    hAdapTeO2po210M2->SetBinContent(i, dBinSize * hnewTeO2po210M2->GetBinContent(i)/hnewTeO2po210M2->GetBinWidth(i));
    hAdapTeO2te125M2->SetBinContent(i, dBinSize * hnewTeO2te125M2->GetBinContent(i)/hnewTeO2te125M2->GetBinWidth(i));
    hAdapTeO2th232M2->SetBinContent(i, dBinSize * hnewTeO2th232M2->GetBinContent(i)/hnewTeO2th232M2->GetBinWidth(i));
    hAdapTeO2th228M2->SetBinContent(i, dBinSize * hnewTeO2th228M2->GetBinContent(i)/hnewTeO2th228M2->GetBinWidth(i));
    hAdapTeO2ra226M2->SetBinContent(i, dBinSize * hnewTeO2ra226M2->GetBinContent(i)/hnewTeO2ra226M2->GetBinWidth(i));
    hAdapTeO2rn222M2->SetBinContent(i, dBinSize * hnewTeO2rn222M2->GetBinContent(i)/hnewTeO2rn222M2->GetBinWidth(i));
    hAdapTeO2u238M2->SetBinContent(i, dBinSize * hnewTeO2u238M2->GetBinContent(i)/hnewTeO2u238M2->GetBinWidth(i));
    hAdapTeO2th230M2->SetBinContent(i, dBinSize * hnewTeO2th230M2->GetBinContent(i)/hnewTeO2th230M2->GetBinWidth(i));
    hAdapTeO2u234M2->SetBinContent(i, dBinSize * hnewTeO2u234M2->GetBinContent(i)/hnewTeO2u234M2->GetBinWidth(i));

    hAdapTeO2Spb210M2_01->SetBinContent(i, dBinSize * hnewTeO2Spb210M2_01->GetBinContent(i)/hnewTeO2Spb210M2_01->GetBinWidth(i));
    hAdapTeO2Spo210M2_001->SetBinContent(i, dBinSize * hnewTeO2Spo210M2_001->GetBinContent(i)/hnewTeO2Spo210M2_001->GetBinWidth(i));
    hAdapTeO2Spo210M2_01->SetBinContent(i, dBinSize * hnewTeO2Spo210M2_01->GetBinContent(i)/hnewTeO2Spo210M2_01->GetBinWidth(i));
    hAdapTeO2Sth232M2_01->SetBinContent(i, dBinSize * hnewTeO2Sth232M2_01->GetBinContent(i)/hnewTeO2Sth232M2_01->GetBinWidth(i));
    hAdapTeO2Su238M2_01->SetBinContent(i, dBinSize * hnewTeO2Su238M2_01->GetBinContent(i)/hnewTeO2Su238M2_01->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_001->SetBinContent(i, dBinSize * hnewTeO2Sxpb210M2_001->GetBinContent(i)/hnewTeO2Sxpb210M2_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_01->SetBinContent(i, dBinSize * hnewTeO2Sxpb210M2_01->GetBinContent(i)/hnewTeO2Sxpb210M2_01->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_1->SetBinContent(i, dBinSize * hnewTeO2Sxpb210M2_1->GetBinContent(i)/hnewTeO2Sxpb210M2_1->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_10->SetBinContent(i, dBinSize * hnewTeO2Sxpb210M2_10->GetBinContent(i)/hnewTeO2Sxpb210M2_10->GetBinWidth(i));
    hAdapTeO2Sxpo210M2_001->SetBinContent(i, dBinSize * hnewTeO2Sxpo210M2_001->GetBinContent(i)/hnewTeO2Sxpo210M2_001->GetBinWidth(i));
    hAdapTeO2Sxpo210M2_01->SetBinContent(i, dBinSize * hnewTeO2Sxpo210M2_01->GetBinContent(i)/hnewTeO2Sxpo210M2_01->GetBinWidth(i));
    hAdapTeO2Sxpo210M2_1->SetBinContent(i, dBinSize * hnewTeO2Sxpo210M2_1->GetBinContent(i)/hnewTeO2Sxpo210M2_1->GetBinWidth(i));
    hAdapTeO2Sxth232M2_001->SetBinContent(i, dBinSize * hnewTeO2Sxth232M2_001->GetBinContent(i)/hnewTeO2Sxth232M2_001->GetBinWidth(i));
    hAdapTeO2Sxth232M2_01->SetBinContent(i, dBinSize * hnewTeO2Sxth232M2_01->GetBinContent(i)/hnewTeO2Sxth232M2_01->GetBinWidth(i));
    hAdapTeO2Sxth232M2_1->SetBinContent(i, dBinSize * hnewTeO2Sxth232M2_1->GetBinContent(i)/hnewTeO2Sxth232M2_1->GetBinWidth(i));
    hAdapTeO2Sxth232M2_10->SetBinContent(i, dBinSize * hnewTeO2Sxth232M2_10->GetBinContent(i)/hnewTeO2Sxth232M2_10->GetBinWidth(i));
    hAdapTeO2Sxu238M2_001->SetBinContent(i, dBinSize * hnewTeO2Sxu238M2_001->GetBinContent(i)/hnewTeO2Sxu238M2_001->GetBinWidth(i));
    hAdapTeO2Sxu238M2_01->SetBinContent(i, dBinSize * hnewTeO2Sxu238M2_01->GetBinContent(i)/hnewTeO2Sxu238M2_01->GetBinWidth(i));
    hAdapTeO2Sxu238M2_1->SetBinContent(i, dBinSize * hnewTeO2Sxu238M2_1->GetBinContent(i)/hnewTeO2Sxu238M2_1->GetBinWidth(i));
    hAdapTeO2Sxu238M2_10->SetBinContent(i, dBinSize * hnewTeO2Sxu238M2_10->GetBinContent(i)/hnewTeO2Sxu238M2_10->GetBinWidth(i));

    hAdapCuFrameco58M2->SetBinContent(i, dBinSize * hnewCuFrameco58M2->GetBinContent(i)/hnewCuFrameco58M2->GetBinWidth(i));
    hAdapCuFrameco60M2->SetBinContent(i, dBinSize * hnewCuFrameco60M2->GetBinContent(i)/hnewCuFrameco60M2->GetBinWidth(i));
    hAdapCuFramecs137M2->SetBinContent(i, dBinSize * hnewCuFramecs137M2->GetBinContent(i)/hnewCuFramecs137M2->GetBinWidth(i));
    hAdapCuFramek40M2->SetBinContent(i, dBinSize * hnewCuFramek40M2->GetBinContent(i)/hnewCuFramek40M2->GetBinWidth(i));
    hAdapCuFramemn54M2->SetBinContent(i, dBinSize * hnewCuFramemn54M2->GetBinContent(i)/hnewCuFramemn54M2->GetBinWidth(i));
    hAdapCuFramepb210M2->SetBinContent(i, dBinSize * hnewCuFramepb210M2->GetBinContent(i)/hnewCuFramepb210M2->GetBinWidth(i));
    hAdapCuFrameth232M2->SetBinContent(i, dBinSize * hnewCuFrameth232M2->GetBinContent(i)/hnewCuFrameth232M2->GetBinWidth(i));
    hAdapCuFrameu238M2->SetBinContent(i, dBinSize * hnewCuFrameu238M2->GetBinContent(i)/hnewCuFrameu238M2->GetBinWidth(i));

    hAdapCuFrameSth232M2_1->SetBinContent(i, dBinSize * hnewCuFrameSth232M2_1->GetBinContent(i)/hnewCuFrameSth232M2_1->GetBinWidth(i));
    hAdapCuFrameSu238M2_1->SetBinContent(i, dBinSize * hnewCuFrameSu238M2_1->GetBinContent(i)/hnewCuFrameSu238M2_1->GetBinWidth(i));
    hAdapCuFrameSxpb210M2_001->SetBinContent(i, dBinSize * hnewCuFrameSxpb210M2_001->GetBinContent(i)/hnewCuFrameSxpb210M2_001->GetBinWidth(i));
    hAdapCuFrameSxpb210M2_01->SetBinContent(i, dBinSize * hnewCuFrameSxpb210M2_01->GetBinContent(i)/hnewCuFrameSxpb210M2_01->GetBinWidth(i));
    hAdapCuFrameSxpb210M2_1->SetBinContent(i, dBinSize * hnewCuFrameSxpb210M2_1->GetBinContent(i)/hnewCuFrameSxpb210M2_1->GetBinWidth(i));
    hAdapCuFrameSxpb210M2_10->SetBinContent(i, dBinSize * hnewCuFrameSxpb210M2_10->GetBinContent(i)/hnewCuFrameSxpb210M2_10->GetBinWidth(i));
    hAdapCuFrameSxth232M2_001->SetBinContent(i, dBinSize * hnewCuFrameSxth232M2_001->GetBinContent(i)/hnewCuFrameSxth232M2_001->GetBinWidth(i));
    hAdapCuFrameSxth232M2_01->SetBinContent(i, dBinSize * hnewCuFrameSxth232M2_01->GetBinContent(i)/hnewCuFrameSxth232M2_01->GetBinWidth(i));
    hAdapCuFrameSxth232M2_1->SetBinContent(i, dBinSize * hnewCuFrameSxth232M2_1->GetBinContent(i)/hnewCuFrameSxth232M2_1->GetBinWidth(i));
    hAdapCuFrameSxth232M2_10->SetBinContent(i, dBinSize * hnewCuFrameSxth232M2_10->GetBinContent(i)/hnewCuFrameSxth232M2_10->GetBinWidth(i));
    hAdapCuFrameSxu238M2_001->SetBinContent(i, dBinSize * hnewCuFrameSxu238M2_001->GetBinContent(i)/hnewCuFrameSxu238M2_001->GetBinWidth(i));
    hAdapCuFrameSxu238M2_01->SetBinContent(i, dBinSize * hnewCuFrameSxu238M2_01->GetBinContent(i)/hnewCuFrameSxu238M2_01->GetBinWidth(i));
    hAdapCuFrameSxu238M2_1->SetBinContent(i, dBinSize * hnewCuFrameSxu238M2_1->GetBinContent(i)/hnewCuFrameSxu238M2_1->GetBinWidth(i));
    hAdapCuFrameSxu238M2_10->SetBinContent(i, dBinSize * hnewCuFrameSxu238M2_10->GetBinContent(i)/hnewCuFrameSxu238M2_10->GetBinWidth(i));

    hAdapCuBoxco58M2->SetBinContent(i, dBinSize * hnewCuBoxco58M2->GetBinContent(i)/hnewCuBoxco58M2->GetBinWidth(i));
    hAdapCuBoxco60M2->SetBinContent(i, dBinSize * hnewCuBoxco60M2->GetBinContent(i)/hnewCuBoxco60M2->GetBinWidth(i));
    hAdapCuBoxcs137M2->SetBinContent(i, dBinSize * hnewCuBoxcs137M2->GetBinContent(i)/hnewCuBoxcs137M2->GetBinWidth(i));
    hAdapCuBoxk40M2->SetBinContent(i, dBinSize * hnewCuBoxk40M2->GetBinContent(i)/hnewCuBoxk40M2->GetBinWidth(i));
    hAdapCuBoxmn54M2->SetBinContent(i, dBinSize * hnewCuBoxmn54M2->GetBinContent(i)/hnewCuBoxmn54M2->GetBinWidth(i));
    hAdapCuBoxpb210M2->SetBinContent(i, dBinSize * hnewCuBoxpb210M2->GetBinContent(i)/hnewCuBoxpb210M2->GetBinWidth(i));
    hAdapCuBoxth232M2->SetBinContent(i, dBinSize * hnewCuBoxth232M2->GetBinContent(i)/hnewCuBoxth232M2->GetBinWidth(i));
    hAdapCuBoxu238M2->SetBinContent(i, dBinSize * hnewCuBoxu238M2->GetBinContent(i)/hnewCuBoxu238M2->GetBinWidth(i));

    hAdapCuBoxSth232M2_1->SetBinContent(i, dBinSize * hnewCuBoxSth232M2_1->GetBinContent(i)/hnewCuBoxSth232M2_1->GetBinWidth(i));
    hAdapCuBoxSu238M2_1->SetBinContent(i, dBinSize * hnewCuBoxSu238M2_1->GetBinContent(i)/hnewCuBoxSu238M2_1->GetBinWidth(i));
    hAdapCuBoxSxpb210M2_001->SetBinContent(i, dBinSize * hnewCuBoxSxpb210M2_001->GetBinContent(i)/hnewCuBoxSxpb210M2_001->GetBinWidth(i));
    hAdapCuBoxSxpb210M2_01->SetBinContent(i, dBinSize * hnewCuBoxSxpb210M2_01->GetBinContent(i)/hnewCuBoxSxpb210M2_01->GetBinWidth(i));
    hAdapCuBoxSxpb210M2_1->SetBinContent(i, dBinSize * hnewCuBoxSxpb210M2_1->GetBinContent(i)/hnewCuBoxSxpb210M2_1->GetBinWidth(i));
    hAdapCuBoxSxpb210M2_10->SetBinContent(i, dBinSize * hnewCuBoxSxpb210M2_10->GetBinContent(i)/hnewCuBoxSxpb210M2_10->GetBinWidth(i));
    hAdapCuBoxSxth232M2_001->SetBinContent(i, dBinSize * hnewCuBoxSxth232M2_001->GetBinContent(i)/hnewCuBoxSxth232M2_001->GetBinWidth(i));
    hAdapCuBoxSxth232M2_01->SetBinContent(i, dBinSize * hnewCuBoxSxth232M2_01->GetBinContent(i)/hnewCuBoxSxth232M2_01->GetBinWidth(i));
    hAdapCuBoxSxth232M2_1->SetBinContent(i, dBinSize * hnewCuBoxSxth232M2_1->GetBinContent(i)/hnewCuBoxSxth232M2_1->GetBinWidth(i));
    hAdapCuBoxSxth232M2_10->SetBinContent(i, dBinSize * hnewCuBoxSxth232M2_10->GetBinContent(i)/hnewCuBoxSxth232M2_10->GetBinWidth(i));
    hAdapCuBoxSxu238M2_001->SetBinContent(i, dBinSize * hnewCuBoxSxu238M2_001->GetBinContent(i)/hnewCuBoxSxu238M2_001->GetBinWidth(i));
    hAdapCuBoxSxu238M2_01->SetBinContent(i, dBinSize * hnewCuBoxSxu238M2_01->GetBinContent(i)/hnewCuBoxSxu238M2_01->GetBinWidth(i));
    hAdapCuBoxSxu238M2_1->SetBinContent(i, dBinSize * hnewCuBoxSxu238M2_1->GetBinContent(i)/hnewCuBoxSxu238M2_1->GetBinWidth(i));
    hAdapCuBoxSxu238M2_10->SetBinContent(i, dBinSize * hnewCuBoxSxu238M2_10->GetBinContent(i)/hnewCuBoxSxu238M2_10->GetBinWidth(i));

    hAdap50mKco58M2->SetBinContent(i, dBinSize * hnew50mKco58M2->GetBinContent(i)/hnew50mKco58M2->GetBinWidth(i));
    hAdap50mKco60M2->SetBinContent(i, dBinSize * hnew50mKco60M2->GetBinContent(i)/hnew50mKco60M2->GetBinWidth(i));
    hAdap50mKcs137M2->SetBinContent(i, dBinSize * hnew50mKcs137M2->GetBinContent(i)/hnew50mKcs137M2->GetBinWidth(i));
    hAdap50mKk40M2->SetBinContent(i, dBinSize * hnew50mKk40M2->GetBinContent(i)/hnew50mKk40M2->GetBinWidth(i));
    hAdap50mKmn54M2->SetBinContent(i, dBinSize * hnew50mKmn54M2->GetBinContent(i)/hnew50mKmn54M2->GetBinWidth(i));
    hAdap50mKpb210M2->SetBinContent(i, dBinSize * hnew50mKpb210M2->GetBinContent(i)/hnew50mKpb210M2->GetBinWidth(i));
    hAdap50mKth232M2->SetBinContent(i, dBinSize * hnew50mKth232M2->GetBinContent(i)/hnew50mKth232M2->GetBinWidth(i));
    hAdap50mKu238M2->SetBinContent(i, dBinSize * hnew50mKu238M2->GetBinContent(i)/hnew50mKu238M2->GetBinWidth(i));

    hAdap600mKco60M2->SetBinContent(i, dBinSize * hnew600mKco60M2->GetBinContent(i)/hnew600mKco60M2->GetBinWidth(i));
    hAdap600mKk40M2->SetBinContent(i, dBinSize * hnew600mKk40M2->GetBinContent(i)/hnew600mKk40M2->GetBinWidth(i));
    hAdap600mKth232M2->SetBinContent(i, dBinSize * hnew600mKth232M2->GetBinContent(i)/hnew600mKth232M2->GetBinWidth(i));
    hAdap600mKu238M2->SetBinContent(i, dBinSize * hnew600mKu238M2->GetBinContent(i)/hnew600mKu238M2->GetBinWidth(i));

    hAdapPbRombi207M2->SetBinContent(i, dBinSize * hnewPbRombi207M2->GetBinContent(i)/hnewPbRombi207M2->GetBinWidth(i));
    hAdapPbRomco60M2->SetBinContent(i, dBinSize * hnewPbRomco60M2->GetBinContent(i)/hnewPbRomco60M2->GetBinWidth(i));
    hAdapPbRomcs137M2->SetBinContent(i, dBinSize * hnewPbRomcs137M2->GetBinContent(i)/hnewPbRomcs137M2->GetBinWidth(i));
    hAdapPbRomk40M2->SetBinContent(i, dBinSize * hnewPbRomk40M2->GetBinContent(i)/hnewPbRomk40M2->GetBinWidth(i));
    hAdapPbRompb210M2->SetBinContent(i, dBinSize * hnewPbRompb210M2->GetBinContent(i)/hnewPbRompb210M2->GetBinWidth(i));
    hAdapPbRomth232M2->SetBinContent(i, dBinSize * hnewPbRomth232M2->GetBinContent(i)/hnewPbRomth232M2->GetBinWidth(i));
    hAdapPbRomu238M2->SetBinContent(i, dBinSize * hnewPbRomu238M2->GetBinContent(i)/hnewPbRomu238M2->GetBinWidth(i));

    hAdapMBco60M2->SetBinContent(i, dBinSize * hnewMBco60M2->GetBinContent(i)/hnewMBco60M2->GetBinWidth(i));
    hAdapMBk40M2->SetBinContent(i, dBinSize * hnewMBk40M2->GetBinContent(i)/hnewMBk40M2->GetBinWidth(i));
    hAdapMBth232M2->SetBinContent(i, dBinSize * hnewMBth232M2->GetBinContent(i)/hnewMBth232M2->GetBinWidth(i));
    hAdapMBu238M2->SetBinContent(i, dBinSize * hnewMBu238M2->GetBinContent(i)/hnewMBu238M2->GetBinWidth(i));

    hAdapSIk40M2->SetBinContent(i, dBinSize * hnewSIk40M2->GetBinContent(i)/hnewSIk40M2->GetBinWidth(i));
    hAdapSIth232M2->SetBinContent(i, dBinSize * hnewSIth232M2->GetBinContent(i)/hnewSIth232M2->GetBinWidth(i));
    hAdapSIu238M2->SetBinContent(i, dBinSize * hnewSIu238M2->GetBinContent(i)/hnewSIu238M2->GetBinWidth(i));


    hAdapIVCco60M2->SetBinContent(i, dBinSize * hnewIVCco60M2->GetBinContent(i)/hnewIVCco60M2->GetBinWidth(i));
    hAdapIVCk40M2->SetBinContent(i, dBinSize * hnewIVCk40M2->GetBinContent(i)/hnewIVCk40M2->GetBinWidth(i));
    hAdapIVCth232M2->SetBinContent(i, dBinSize * hnewIVCth232M2->GetBinContent(i)/hnewIVCth232M2->GetBinWidth(i));
    hAdapIVCu238M2->SetBinContent(i, dBinSize * hnewIVCu238M2->GetBinContent(i)/hnewIVCu238M2->GetBinWidth(i));

    hAdapOVCco60M2->SetBinContent(i, dBinSize * hnewOVCco60M2->GetBinContent(i)/hnewOVCco60M2->GetBinWidth(i));
    hAdapOVCk40M2->SetBinContent(i, dBinSize * hnewOVCk40M2->GetBinContent(i)/hnewOVCk40M2->GetBinWidth(i));
    hAdapOVCth232M2->SetBinContent(i, dBinSize * hnewOVCth232M2->GetBinContent(i)/hnewOVCth232M2->GetBinWidth(i));
    hAdapOVCu238M2->SetBinContent(i, dBinSize * hnewOVCu238M2->GetBinContent(i)/hnewOVCu238M2->GetBinWidth(i));

  }
}

// Loads the background data
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

  // Currently using Jon's reduced file -- change for other input files
/*
  qtree->Add("/Users/brian/macros/Simulations/Toy/combi2/combi2.root"); 
  qtree->Project("fDataHistoTot", "Ener2");
  qtree->Project("fDataHistoM1",  "Ener2", "Multiplicity == 1");
  qtree->Project("fDataHistoM2",  "Ener2", "Multiplicity == 2");
*/
  
  // qtree->Add("/Users/brian/macros/CUOREZ/Bkg/Q0_DR2_BackgroundSignalData.root"); 
  qtree->Add("/Users/brian/macros/CUOREZ/Bkg/ReducedBkgSync-ds*.root");   
  qtree->Project("fDataHistoTot", "Energy", base_cut);
  qtree->Project("fDataHistoM1",  "Energy", base_cut && "Multiplicity_Sync == 1");
  qtree->Project("fDataHistoM2",  "Energy", base_cut && "Multiplicity_Sync == 2");

  cout << "Loaded Data" << endl;
}

// Prints parameters, needs update 11-06-2014
void TBackgroundModel::PrintParameters()
{
  for(int i = 0; i < dNParam; i++)
  {
    cout << Form("Par%d = ", i) << fParameters[i] << " +/- " << fParError[i] << endl;
  }
}


// Set Parameters in Model
void TBackgroundModel::SetParameters(int index, double value)
{
	// Change the index max depending on model
	if(index > dNParam) cout << "Index too large" << endl;
	else fParameters[index] = value;
}


// Creates/updates the background model
void TBackgroundModel::UpdateModelAdaptive()
{
  if(fModelTotAdapM1 == NULL) 
  {
    cout << "Model Histogram Not Created" << endl;
    return;
  }

  // Reset all bins in model histogram(s)
  fModelTotAdapM1->Reset();
  fModelTotAdapM2->Reset();
  fModelTotAdapAlphaM1->Reset();
  fModleTotAdapAlphaM2->Reset();

  dNumCalls++;
  if(dNumCalls%500==0)
  {
    cout << "Call #: "<< dNumCalls << endl;
  }

  // Create model
  /////////////////////////////////////
  //// Adaptive Binning parameters
  ////////////////////////////////////
  // M1
  fModelTotAdapM1->Add( hAdapTeO20nuM1,      fParameters[0]);
  fModelTotAdapM1->Add( hAdapTeO22nuM1,      fParameters[1]);
  fModelTotAdapM1->Add( hAdapTeO2co60M1,     fParameters[2]);
  fModelTotAdapM1->Add( hAdapTeO2k40M1,      fParameters[3]);
  fModelTotAdapM1->Add( hAdapTeO2te125M1,    fParameters[6]);
  // fModelTotAdapM1->Add( hAdapTeO2pb210M1,    fParameters[4]);
  // fModelTotAdapM1->Add( hAdapTeO2po210M1,    fParameters[5]);
  // fModelTotAdapM1->Add( hAdapTeO2th232M1,    fParameters[7]);
  // fModelTotAdapM1->Add( hAdapTeO2th228M1,    fParameters[8]);
  // fModelTotAdapM1->Add( hAdapTeO2ra226M1,    fParameters[9]);
  // fModelTotAdapM1->Add( hAdapTeO2rn222M1,    fParameters[10]);
  // fModelTotAdapM1->Add( hAdapTeO2u238M1,     fParameters[11]);
  // fModelTotAdapM1->Add( hAdapTeO2th230M1,    fParameters[12]);
  // fModelTotAdapM1->Add( hAdapTeO2u234M1,     fParameters[13]);

  fModelTotAdapM1->Add( hAdapCuFrameco58M1,  fParameters[14]);
  fModelTotAdapM1->Add( hAdapCuFrameco60M1,  fParameters[15]);
  fModelTotAdapM1->Add( hAdapCuFramecs137M1, fParameters[16]);
  fModelTotAdapM1->Add( hAdapCuFramek40M1,   fParameters[17]);
  fModelTotAdapM1->Add( hAdapCuFramemn54M1,  fParameters[18]);
  fModelTotAdapM1->Add( hAdapCuFramepb210M1, fParameters[19]);
  fModelTotAdapM1->Add( hAdapCuFrameth232M1, fParameters[20]);
  fModelTotAdapM1->Add( hAdapCuFrameu238M1,  fParameters[21]);

  fModelTotAdapM1->Add( hAdapCuBoxco58M1,    fParameters[22]);
  fModelTotAdapM1->Add( hAdapCuBoxco60M1,    fParameters[23]);
  fModelTotAdapM1->Add( hAdapCuBoxcs137M1,   fParameters[24]);
  fModelTotAdapM1->Add( hAdapCuBoxk40M1,     fParameters[25]);
  fModelTotAdapM1->Add( hAdapCuBoxmn54M1,    fParameters[26]);
  fModelTotAdapM1->Add( hAdapCuBoxpb210M1,   fParameters[27]);
  fModelTotAdapM1->Add( hAdapCuBoxth232M1,   fParameters[28]);
  fModelTotAdapM1->Add( hAdapCuBoxu238M1,    fParameters[29]);

  fModelTotAdapM1->Add( hAdap50mKco58M1,     fParameters[30]);
  fModelTotAdapM1->Add( hAdap50mKco60M1,     fParameters[31]);
  fModelTotAdapM1->Add( hAdap50mKcs137M1,    fParameters[32]);
  fModelTotAdapM1->Add( hAdap50mKk40M1,      fParameters[33]);
  fModelTotAdapM1->Add( hAdap50mKmn54M1,     fParameters[34]);
  fModelTotAdapM1->Add( hAdap50mKpb210M1,    fParameters[35]);
  fModelTotAdapM1->Add( hAdap50mKth232M1,    fParameters[36]);
  fModelTotAdapM1->Add( hAdap50mKu238M1,     fParameters[37]);

  fModelTotAdapM1->Add( hAdap600mKco60M1,    fParameters[38]);
  fModelTotAdapM1->Add( hAdap600mKk40M1,     fParameters[39]);
  fModelTotAdapM1->Add( hAdap600mKth232M1,   fParameters[40]);
  fModelTotAdapM1->Add( hAdap600mKu238M1,    fParameters[41]);

  fModelTotAdapM1->Add( hAdapPbRombi207M1,   fParameters[42]);
  fModelTotAdapM1->Add( hAdapPbRomco60M1,    fParameters[43]);
  fModelTotAdapM1->Add( hAdapPbRomcs137M1,   fParameters[44]);
  fModelTotAdapM1->Add( hAdapPbRomk40M1,     fParameters[45]);
  fModelTotAdapM1->Add( hAdapPbRompb210M1,   fParameters[46]);
  fModelTotAdapM1->Add( hAdapPbRomth232M1,   fParameters[47]);
  fModelTotAdapM1->Add( hAdapPbRomu238M1,    fParameters[48]);

  fModelTotAdapM1->Add( hAdapMBco60M1,       fParameters[49]);
  fModelTotAdapM1->Add( hAdapMBk40M1,        fParameters[50]);
  fModelTotAdapM1->Add( hAdapMBth232M1,      fParameters[51]);
  fModelTotAdapM1->Add( hAdapMBu238M1,       fParameters[52]);

  fModelTotAdapM1->Add( hAdapIVCco60M1,      fParameters[53]);
  fModelTotAdapM1->Add( hAdapIVCk40M1,       fParameters[54]);
  fModelTotAdapM1->Add( hAdapIVCth232M1,     fParameters[55]);
  fModelTotAdapM1->Add( hAdapIVCu238M1,      fParameters[56]);

  fModelTotAdapM1->Add( hAdapOVCco60M1,      fParameters[57]);
  fModelTotAdapM1->Add( hAdapOVCk40M1,       fParameters[58]);
  fModelTotAdapM1->Add( hAdapOVCth232M1,     fParameters[59]);
  fModelTotAdapM1->Add( hAdapOVCu238M1,      fParameters[60]);

  fModelTotAdapM1->Add( hAdapSIk40M1,        fParameters[61]);
  fModelTotAdapM1->Add( hAdapSIth232M1,      fParameters[62]);
  fModelTotAdapM1->Add( hAdapSIu238M1,       fParameters[63]);

/*
  fModelTotAdapM1->Add( hAdapTeO2Spb210M1_01,      fParameters[64]);
  fModelTotAdapM1->Add( hAdapTeO2Spo210M1_001,     fParameters[65]);
  fModelTotAdapM1->Add( hAdapTeO2Spo210M1_01,      fParameters[66]);
  fModelTotAdapM1->Add( hAdapTeO2Sth232M1_01,      fParameters[67]);
  fModelTotAdapM1->Add( hAdapTeO2Su238M1_01,       fParameters[68]);
  fModelTotAdapM1->Add( hAdapTeO2Sxpb210M1_001,    fParameters[69]);
  fModelTotAdapM1->Add( hAdapTeO2Sxpb210M1_01,     fParameters[70]);
  fModelTotAdapM1->Add( hAdapTeO2Sxpb210M1_1,      fParameters[71]);
  fModelTotAdapM1->Add( hAdapTeO2Sxpb210M1_10,     fParameters[72]);
  fModelTotAdapM1->Add( hAdapTeO2Sxpo210M1_001,    fParameters[73]);
  fModelTotAdapM1->Add( hAdapTeO2Sxpo210M1_01,     fParameters[74]);
  fModelTotAdapM1->Add( hAdapTeO2Sxpo210M1_1,      fParameters[75]);
  fModelTotAdapM1->Add( hAdapTeO2Sxth232M1_001,    fParameters[76]);
  fModelTotAdapM1->Add( hAdapTeO2Sxth232M1_01,     fParameters[77]);
  fModelTotAdapM1->Add( hAdapTeO2Sxth232M1_1,      fParameters[78]);
  fModelTotAdapM1->Add( hAdapTeO2Sxth232M1_10,     fParameters[79]);
  fModelTotAdapM1->Add( hAdapTeO2Sxu238M1_001,     fParameters[80]);
  fModelTotAdapM1->Add( hAdapTeO2Sxu238M1_01,      fParameters[81]);
  fModelTotAdapM1->Add( hAdapTeO2Sxu238M1_1,       fParameters[82]);
  fModelTotAdapM1->Add( hAdapTeO2Sxu238M1_10,      fParameters[83]);

  fModelTotAdapM1->Add( hAdapCuFrameSth232M1_1,    fParameters[84]);
  fModelTotAdapM1->Add( hAdapCuFrameSu238M1_1,     fParameters[85]);
  fModelTotAdapM1->Add( hAdapCuFrameSxpb210M1_001, fParameters[86]);
  fModelTotAdapM1->Add( hAdapCuFrameSxpb210M1_01,  fParameters[87]);
  fModelTotAdapM1->Add( hAdapCuFrameSxpb210M1_1,   fParameters[88]);
  fModelTotAdapM1->Add( hAdapCuFrameSxpb210M1_10,  fParameters[89]);
  fModelTotAdapM1->Add( hAdapCuFrameSxth232M1_001, fParameters[90]);
  fModelTotAdapM1->Add( hAdapCuFrameSxth232M1_01,  fParameters[91]);
  fModelTotAdapM1->Add( hAdapCuFrameSxth232M1_1,   fParameters[92]);
  fModelTotAdapM1->Add( hAdapCuFrameSxth232M1_10,  fParameters[93]);
  fModelTotAdapM1->Add( hAdapCuFrameSxu238M1_001,  fParameters[94]);
  fModelTotAdapM1->Add( hAdapCuFrameSxu238M1_01,   fParameters[95]);
  fModelTotAdapM1->Add( hAdapCuFrameSxu238M1_1,    fParameters[96]);
  fModelTotAdapM1->Add( hAdapCuFrameSxu238M1_10,   fParameters[97]);

  fModelTotAdapM1->Add( hAdapCuBoxSth232M1_1,      fParameters[98]);
  fModelTotAdapM1->Add( hAdapCuBoxSu238M1_1,       fParameters[99]);
  fModelTotAdapM1->Add( hAdapCuBoxSxpb210M1_001,   fParameters[100]);
  fModelTotAdapM1->Add( hAdapCuBoxSxpb210M1_01,    fParameters[101]);
  fModelTotAdapM1->Add( hAdapCuBoxSxpb210M1_1,     fParameters[102]);
  fModelTotAdapM1->Add( hAdapCuBoxSxpb210M1_10,    fParameters[103]);
  fModelTotAdapM1->Add( hAdapCuBoxSxth232M1_001,   fParameters[104]);
  fModelTotAdapM1->Add( hAdapCuBoxSxth232M1_01,    fParameters[105]);
  fModelTotAdapM1->Add( hAdapCuBoxSxth232M1_1,     fParameters[106]);
  fModelTotAdapM1->Add( hAdapCuBoxSxth232M1_10,    fParameters[107]);
  fModelTotAdapM1->Add( hAdapCuBoxSxu238M1_001,    fParameters[108]);
  fModelTotAdapM1->Add( hAdapCuBoxSxu238M1_01,     fParameters[109]);
  fModelTotAdapM1->Add( hAdapCuBoxSxu238M1_1,      fParameters[110]);
  fModelTotAdapM1->Add( hAdapCuBoxSxu238M1_10,     fParameters[111]);   
*/

// M2
  fModelTotAdapM2->Add( hAdapTeO20nuM2,      fParameters[0]);
  fModelTotAdapM2->Add( hAdapTeO22nuM2,      fParameters[1]);
  fModelTotAdapM2->Add( hAdapTeO2co60M2,     fParameters[2]);
  fModelTotAdapM2->Add( hAdapTeO2k40M2,      fParameters[3]);
  fModelTotAdapM2->Add( hAdapTeO2te125M2,    fParameters[6]);
  // fModelTotAdapM2->Add( hAdapTeO2pb210M2,    fParameters[4]);
  // fModelTotAdapM2->Add( hAdapTeO2po210M2,    fParameters[5]);
  // fModelTotAdapM2->Add( hAdapTeO2th232M2,    fParameters[7]);
  // fModelTotAdapM2->Add( hAdapTeO2th228M2,    fParameters[8]);
  // fModelTotAdapM2->Add( hAdapTeO2ra226M2,    fParameters[9]);
  // fModelTotAdapM2->Add( hAdapTeO2rn222M2,    fParameters[10]);
  // fModelTotAdapM2->Add( hAdapTeO2u238M2,     fParameters[11]);
  // fModelTotAdapM2->Add( hAdapTeO2th230M2,    fParameters[12]);
  // fModelTotAdapM2->Add( hAdapTeO2u234M2,     fParameters[13]);

  fModelTotAdapM2->Add( hAdapCuFrameco58M2,  fParameters[14]);
  fModelTotAdapM2->Add( hAdapCuFrameco60M2,  fParameters[15]);
  fModelTotAdapM2->Add( hAdapCuFramecs137M2, fParameters[16]);
  fModelTotAdapM2->Add( hAdapCuFramek40M2,   fParameters[17]);
  fModelTotAdapM2->Add( hAdapCuFramemn54M2,  fParameters[18]);
  fModelTotAdapM2->Add( hAdapCuFramepb210M2, fParameters[19]);
  fModelTotAdapM2->Add( hAdapCuFrameth232M2, fParameters[20]);
  fModelTotAdapM2->Add( hAdapCuFrameu238M2,  fParameters[21]);

  fModelTotAdapM2->Add( hAdapCuBoxco58M2,    fParameters[22]);
  fModelTotAdapM2->Add( hAdapCuBoxco60M2,    fParameters[23]);
  fModelTotAdapM2->Add( hAdapCuBoxcs137M2,   fParameters[24]);
  fModelTotAdapM2->Add( hAdapCuBoxk40M2,     fParameters[25]);
  fModelTotAdapM2->Add( hAdapCuBoxmn54M2,    fParameters[26]);
  fModelTotAdapM2->Add( hAdapCuBoxpb210M2,   fParameters[27]);
  fModelTotAdapM2->Add( hAdapCuBoxth232M2,   fParameters[28]);
  fModelTotAdapM2->Add( hAdapCuBoxu238M2,    fParameters[29]);

  fModelTotAdapM2->Add( hAdap50mKco58M2,     fParameters[30]);
  fModelTotAdapM2->Add( hAdap50mKco60M2,     fParameters[31]);
  fModelTotAdapM2->Add( hAdap50mKcs137M2,    fParameters[32]);
  fModelTotAdapM2->Add( hAdap50mKk40M2,      fParameters[33]);
  fModelTotAdapM2->Add( hAdap50mKmn54M2,     fParameters[34]);
  fModelTotAdapM2->Add( hAdap50mKpb210M2,    fParameters[35]);
  fModelTotAdapM2->Add( hAdap50mKth232M2,    fParameters[36]);
  fModelTotAdapM2->Add( hAdap50mKu238M2,     fParameters[37]);

  fModelTotAdapM2->Add( hAdap600mKco60M2,    fParameters[38]);
  fModelTotAdapM2->Add( hAdap600mKk40M2,     fParameters[39]);
  fModelTotAdapM2->Add( hAdap600mKth232M2,   fParameters[40]);
  fModelTotAdapM2->Add( hAdap600mKu238M2,    fParameters[41]);

  fModelTotAdapM2->Add( hAdapPbRombi207M2,   fParameters[42]);
  fModelTotAdapM2->Add( hAdapPbRomco60M2,    fParameters[43]);
  fModelTotAdapM2->Add( hAdapPbRomcs137M2,   fParameters[44]);
  fModelTotAdapM2->Add( hAdapPbRomk40M2,     fParameters[45]);
  fModelTotAdapM2->Add( hAdapPbRompb210M2,   fParameters[46]);
  fModelTotAdapM2->Add( hAdapPbRomth232M2,   fParameters[47]);
  fModelTotAdapM2->Add( hAdapPbRomu238M2,    fParameters[48]);

  fModelTotAdapM2->Add( hAdapMBco60M2,       fParameters[49]);
  fModelTotAdapM2->Add( hAdapMBk40M2,        fParameters[50]);
  fModelTotAdapM2->Add( hAdapMBth232M2,      fParameters[51]);
  fModelTotAdapM2->Add( hAdapMBu238M2,       fParameters[52]);

  fModelTotAdapM2->Add( hAdapIVCco60M2,      fParameters[53]);
  fModelTotAdapM2->Add( hAdapIVCk40M2,       fParameters[54]);
  fModelTotAdapM2->Add( hAdapIVCth232M2,     fParameters[55]);
  fModelTotAdapM2->Add( hAdapIVCu238M2,      fParameters[56]);

  fModelTotAdapM2->Add( hAdapOVCco60M2,      fParameters[57]);
  fModelTotAdapM2->Add( hAdapOVCk40M2,       fParameters[58]);
  fModelTotAdapM2->Add( hAdapOVCth232M2,     fParameters[59]);
  fModelTotAdapM2->Add( hAdapOVCu238M2,      fParameters[60]);  

  fModelTotAdapM2->Add( hAdapSIk40M2,        fParameters[61]);
  fModelTotAdapM2->Add( hAdapSIth232M2,      fParameters[62]);
  fModelTotAdapM2->Add( hAdapSIu238M2,       fParameters[63]);  

/*
  fModelTotAdapM2->Add( hAdapTeO2Spb210M2_01,      fParameters[64]);
  fModelTotAdapM2->Add( hAdapTeO2Spo210M2_001,     fParameters[65]);
  fModelTotAdapM2->Add( hAdapTeO2Spo210M2_01,      fParameters[66]);
  fModelTotAdapM2->Add( hAdapTeO2Sth232M2_01,      fParameters[67]);
  fModelTotAdapM2->Add( hAdapTeO2Su238M2_01,       fParameters[68]);
  fModelTotAdapM2->Add( hAdapTeO2Sxpb210M2_001,    fParameters[69]);
  fModelTotAdapM2->Add( hAdapTeO2Sxpb210M2_01,     fParameters[70]);
  fModelTotAdapM2->Add( hAdapTeO2Sxpb210M2_1,      fParameters[71]);
  fModelTotAdapM2->Add( hAdapTeO2Sxpb210M2_10,     fParameters[72]);
  fModelTotAdapM2->Add( hAdapTeO2Sxpo210M2_001,    fParameters[73]);
  fModelTotAdapM2->Add( hAdapTeO2Sxpo210M2_01,     fParameters[74]);
  fModelTotAdapM2->Add( hAdapTeO2Sxpo210M2_1,      fParameters[75]);
  fModelTotAdapM2->Add( hAdapTeO2Sxth232M2_001,    fParameters[76]);
  fModelTotAdapM2->Add( hAdapTeO2Sxth232M2_01,     fParameters[77]);
  fModelTotAdapM2->Add( hAdapTeO2Sxth232M2_1,      fParameters[78]);
  fModelTotAdapM2->Add( hAdapTeO2Sxth232M2_10,     fParameters[79]);
  fModelTotAdapM2->Add( hAdapTeO2Sxu238M2_001,     fParameters[80]);
  fModelTotAdapM2->Add( hAdapTeO2Sxu238M2_01,      fParameters[81]);
  fModelTotAdapM2->Add( hAdapTeO2Sxu238M2_1,       fParameters[82]);
  fModelTotAdapM2->Add( hAdapTeO2Sxu238M2_10,      fParameters[83]);

  fModelTotAdapM2->Add( hAdapCuFrameSth232M2_1,    fParameters[84]);
  fModelTotAdapM2->Add( hAdapCuFrameSu238M2_1,     fParameters[85]);
  fModelTotAdapM2->Add( hAdapCuFrameSxpb210M2_001, fParameters[86]);
  fModelTotAdapM2->Add( hAdapCuFrameSxpb210M2_01,  fParameters[87]);
  fModelTotAdapM2->Add( hAdapCuFrameSxpb210M2_1,   fParameters[88]);
  fModelTotAdapM2->Add( hAdapCuFrameSxpb210M2_10,  fParameters[89]);
  fModelTotAdapM2->Add( hAdapCuFrameSxth232M2_001, fParameters[90]);
  fModelTotAdapM2->Add( hAdapCuFrameSxth232M2_01,  fParameters[91]);
  fModelTotAdapM2->Add( hAdapCuFrameSxth232M2_1,   fParameters[92]);
  fModelTotAdapM2->Add( hAdapCuFrameSxth232M2_10,  fParameters[93]);
  fModelTotAdapM2->Add( hAdapCuFrameSxu238M2_001,  fParameters[94]);
  fModelTotAdapM2->Add( hAdapCuFrameSxu238M2_01,   fParameters[95]);
  fModelTotAdapM2->Add( hAdapCuFrameSxu238M2_1,    fParameters[96]);
  fModelTotAdapM2->Add( hAdapCuFrameSxu238M2_10,   fParameters[97]);

  fModelTotAdapM2->Add( hAdapCuBoxSth232M2_1,      fParameters[98]);
  fModelTotAdapM2->Add( hAdapCuBoxSu238M2_1,       fParameters[99]);
  fModelTotAdapM2->Add( hAdapCuBoxSxpb210M2_001,   fParameters[100]);
  fModelTotAdapM2->Add( hAdapCuBoxSxpb210M2_01,    fParameters[101]);
  fModelTotAdapM2->Add( hAdapCuBoxSxpb210M2_1,     fParameters[102]);
  fModelTotAdapM2->Add( hAdapCuBoxSxpb210M2_10,    fParameters[103]);
  fModelTotAdapM2->Add( hAdapCuBoxSxth232M2_001,   fParameters[104]);
  fModelTotAdapM2->Add( hAdapCuBoxSxth232M2_01,    fParameters[105]);
  fModelTotAdapM2->Add( hAdapCuBoxSxth232M2_1,     fParameters[106]);
  fModelTotAdapM2->Add( hAdapCuBoxSxth232M2_10,    fParameters[107]);
  fModelTotAdapM2->Add( hAdapCuBoxSxu238M2_001,    fParameters[108]);
  fModelTotAdapM2->Add( hAdapCuBoxSxu238M2_01,     fParameters[109]);
  fModelTotAdapM2->Add( hAdapCuBoxSxu238M2_1,      fParameters[110]);
  fModelTotAdapM2->Add( hAdapCuBoxSxu238M2_10,     fParameters[111]);  
*/

  //// Energy scale changes for alphas
  // Need 1 histogram to add all the alphas before energy scale, 2 more histograms for the energy scale

  // Add together all the alphas into 1 histogram
  fModelTotAdapAlphaM1->Add( hAdapTeO2pb210M1,    fParameters[4]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2po210M1,    fParameters[5]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2th232M1,    fParameters[7]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2th228M1,    fParameters[8]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2ra226M1,    fParameters[9]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2rn222M1,    fParameters[10]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2u238M1,     fParameters[11]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2th230M1,    fParameters[12]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2u234M1,     fParameters[13]);

  fModelTotAdapAlphaM1->Add( hAdapTeO2Spb210M1_01,      fParameters[64]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Spo210M1_001,     fParameters[65]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Spo210M1_01,      fParameters[66]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sth232M1_01,      fParameters[67]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Su238M1_01,       fParameters[68]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxpb210M1_001,    fParameters[69]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxpb210M1_01,     fParameters[70]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxpb210M1_1,      fParameters[71]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxpb210M1_10,     fParameters[72]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxpo210M1_001,    fParameters[73]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxpo210M1_01,     fParameters[74]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxpo210M1_1,      fParameters[75]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxth232M1_001,    fParameters[76]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxth232M1_01,     fParameters[77]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxth232M1_1,      fParameters[78]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxth232M1_10,     fParameters[79]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxu238M1_001,     fParameters[80]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxu238M1_01,      fParameters[81]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxu238M1_1,       fParameters[82]);
  fModelTotAdapAlphaM1->Add( hAdapTeO2Sxu238M1_10,      fParameters[83]);

  fModelTotAdapAlphaM1->Add( hAdapCuFrameSth232M1_1,    fParameters[84]);
  fModelTotAdapAlphaM1->Add( hAdapCuFrameSu238M1_1,     fParameters[85]);
  fModelTotAdapAlphaM1->Add( hAdapCuFrameSxpb210M1_001, fParameters[86]);
  fModelTotAdapAlphaM1->Add( hAdapCuFrameSxpb210M1_01,  fParameters[87]);
  fModelTotAdapAlphaM1->Add( hAdapCuFrameSxpb210M1_1,   fParameters[88]);
  fModelTotAdapAlphaM1->Add( hAdapCuFrameSxpb210M1_10,  fParameters[89]);
  fModelTotAdapAlphaM1->Add( hAdapCuFrameSxth232M1_001, fParameters[90]);
  fModelTotAdapAlphaM1->Add( hAdapCuFrameSxth232M1_01,  fParameters[91]);
  fModelTotAdapAlphaM1->Add( hAdapCuFrameSxth232M1_1,   fParameters[92]);
  fModelTotAdapAlphaM1->Add( hAdapCuFrameSxth232M1_10,  fParameters[93]);
  fModelTotAdapAlphaM1->Add( hAdapCuFrameSxu238M1_001,  fParameters[94]);
  fModelTotAdapAlphaM1->Add( hAdapCuFrameSxu238M1_01,   fParameters[95]);
  fModelTotAdapAlphaM1->Add( hAdapCuFrameSxu238M1_1,    fParameters[96]);
  fModelTotAdapAlphaM1->Add( hAdapCuFrameSxu238M1_10,   fParameters[97]);

  fModelTotAdapAlphaM1->Add( hAdapCuBoxSth232M1_1,      fParameters[98]);
  fModelTotAdapAlphaM1->Add( hAdapCuBoxSu238M1_1,       fParameters[99]);
  fModelTotAdapAlphaM1->Add( hAdapCuBoxSxpb210M1_001,   fParameters[100]);
  fModelTotAdapAlphaM1->Add( hAdapCuBoxSxpb210M1_01,    fParameters[101]);
  fModelTotAdapAlphaM1->Add( hAdapCuBoxSxpb210M1_1,     fParameters[102]);
  fModelTotAdapAlphaM1->Add( hAdapCuBoxSxpb210M1_10,    fParameters[103]);
  fModelTotAdapAlphaM1->Add( hAdapCuBoxSxth232M1_001,   fParameters[104]);
  fModelTotAdapAlphaM1->Add( hAdapCuBoxSxth232M1_01,    fParameters[105]);
  fModelTotAdapAlphaM1->Add( hAdapCuBoxSxth232M1_1,     fParameters[106]);
  fModelTotAdapAlphaM1->Add( hAdapCuBoxSxth232M1_10,    fParameters[107]);
  fModelTotAdapAlphaM1->Add( hAdapCuBoxSxu238M1_001,    fParameters[108]);
  fModelTotAdapAlphaM1->Add( hAdapCuBoxSxu238M1_01,     fParameters[109]);
  fModelTotAdapAlphaM1->Add( hAdapCuBoxSxu238M1_1,      fParameters[110]);
  fModelTotAdapAlphaM1->Add( hAdapCuBoxSxu238M1_10,     fParameters[111]);   


  fModelTotAdapAlphaM2->Add( hAdapTeO2pb210M2,    fParameters[4]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2po210M2,    fParameters[5]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2th232M2,    fParameters[7]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2th228M2,    fParameters[8]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2ra226M2,    fParameters[9]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2rn222M2,    fParameters[10]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2u238M2,     fParameters[11]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2th230M2,    fParameters[12]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2u234M2,     fParameters[13]);

  fModelTotAdapAlphaM2->Add( hAdapTeO2Spb210M2_01,      fParameters[64]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Spo210M2_001,     fParameters[65]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Spo210M2_01,      fParameters[66]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sth232M2_01,      fParameters[67]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Su238M2_01,       fParameters[68]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxpb210M2_001,    fParameters[69]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxpb210M2_01,     fParameters[70]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxpb210M2_1,      fParameters[71]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxpb210M2_10,     fParameters[72]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxpo210M2_001,    fParameters[73]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxpo210M2_01,     fParameters[74]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxpo210M2_1,      fParameters[75]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxth232M2_001,    fParameters[76]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxth232M2_01,     fParameters[77]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxth232M2_1,      fParameters[78]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxth232M2_10,     fParameters[79]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxu238M2_001,     fParameters[80]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxu238M2_01,      fParameters[81]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxu238M2_1,       fParameters[82]);
  fModelTotAdapAlphaM2->Add( hAdapTeO2Sxu238M2_10,      fParameters[83]);

  fModelTotAdapAlphaM2->Add( hAdapCuFrameSth232M2_1,    fParameters[84]);
  fModelTotAdapAlphaM2->Add( hAdapCuFrameSu238M2_1,     fParameters[85]);
  fModelTotAdapAlphaM2->Add( hAdapCuFrameSxpb210M2_001, fParameters[86]);
  fModelTotAdapAlphaM2->Add( hAdapCuFrameSxpb210M2_01,  fParameters[87]);
  fModelTotAdapAlphaM2->Add( hAdapCuFrameSxpb210M2_1,   fParameters[88]);
  fModelTotAdapAlphaM2->Add( hAdapCuFrameSxpb210M2_10,  fParameters[89]);
  fModelTotAdapAlphaM2->Add( hAdapCuFrameSxth232M2_001, fParameters[90]);
  fModelTotAdapAlphaM2->Add( hAdapCuFrameSxth232M2_01,  fParameters[91]);
  fModelTotAdapAlphaM2->Add( hAdapCuFrameSxth232M2_1,   fParameters[92]);
  fModelTotAdapAlphaM2->Add( hAdapCuFrameSxth232M2_10,  fParameters[93]);
  fModelTotAdapAlphaM2->Add( hAdapCuFrameSxu238M2_001,  fParameters[94]);
  fModelTotAdapAlphaM2->Add( hAdapCuFrameSxu238M2_01,   fParameters[95]);
  fModelTotAdapAlphaM2->Add( hAdapCuFrameSxu238M2_1,    fParameters[96]);
  fModelTotAdapAlphaM2->Add( hAdapCuFrameSxu238M2_10,   fParameters[97]);

  fModelTotAdapAlphaM2->Add( hAdapCuBoxSth232M2_1,      fParameters[98]);
  fModelTotAdapAlphaM2->Add( hAdapCuBoxSu238M2_1,       fParameters[99]);
  fModelTotAdapAlphaM2->Add( hAdapCuBoxSxpb210M2_001,   fParameters[100]);
  fModelTotAdapAlphaM2->Add( hAdapCuBoxSxpb210M2_01,    fParameters[101]);
  fModelTotAdapAlphaM2->Add( hAdapCuBoxSxpb210M2_1,     fParameters[102]);
  fModelTotAdapAlphaM2->Add( hAdapCuBoxSxpb210M2_10,    fParameters[103]);
  fModelTotAdapAlphaM2->Add( hAdapCuBoxSxth232M2_001,   fParameters[104]);
  fModelTotAdapAlphaM2->Add( hAdapCuBoxSxth232M2_01,    fParameters[105]);
  fModelTotAdapAlphaM2->Add( hAdapCuBoxSxth232M2_1,     fParameters[106]);
  fModelTotAdapAlphaM2->Add( hAdapCuBoxSxth232M2_10,    fParameters[107]);
  fModelTotAdapAlphaM2->Add( hAdapCuBoxSxu238M2_001,    fParameters[108]);
  fModelTotAdapAlphaM2->Add( hAdapCuBoxSxu238M2_01,     fParameters[109]);
  fModelTotAdapAlphaM2->Add( hAdapCuBoxSxu238M2_1,      fParameters[110]);
  fModelTotAdapAlphaM2->Add( hAdapCuBoxSxu238M2_10,     fParameters[111]);   

  // Create the 2 energyscale alphas...
  fModelTotAdapAlphaHighM1->Add( EnergyScale(fModelTotAdapAlphaM1, hEnergyScaleDummyM1, 0, fParameters[112]) );
  fModelTotAdapAlphaLowM1->Add( EnergyScale(fModelTotAdapAlphaM1, hEnergyScaleDummyM1, 0, fParameters[113]) );

  fModelTotAdapAlphaHighM2->Add( EnergyScale(fModelTotAdapAlphaM2, hEnergyScaleDummyM2, 0, fParameters[112]) );
  fModelTotAdapAlphaLowM2->Add( EnergyScale(fModelTotAdapAlphaM2, hEnergyScaleDummyM2, 0, fParameters[113]) );

  // Add together the 3 energy scale histograms into total histogram.. parameter can become floating in the future..
  fModelTotAdapM1->Add( fModelTotAdapAlphaM1, 0.8 );
  fModelTotAdapM1->Add( fModelTotAdapAlphaHighM1, 0.1 );
  fModelTotAdapM1->Add( fModelTotAdapAlphaLowM1, 0.1 );

  fModelTotAdapM2->Add( fModelTotAdapAlphaM2, 0.8 );
  fModelTotAdapM2->Add( fModelTotAdapAlphaHighM2, 0.1 );
  fModelTotAdapM2->Add( fModelTotAdapAlphaLowM2, 0.1 );


}

// For whatever tests...
void TBackgroundModel::Test()
{ 
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->SetLogy();
  TH1D *hTest1 = new TH1D("hTest1", "hTest1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  TH1D *hTest2 = new TH1D("hTest2", "hTest2", dAdaptiveBinsM1, dAdaptiveArrayM1);


  hTest1->Add(EnergyScale(hAdapTeO2po210M1, hEnergyScaleDummyM1, 0, 0.999),1);
  hTest2->Add(EnergyScale(hAdapTeO2po210M1, hEnergyScaleDummyM1, 0, 1.001),1);
  hAdapTeO2po210M1->SetLineColor(1);
  hTest1->SetLineColor(2);
  // hTest1->SetLineStyle(2);
  hTest2->SetLineColor(3);
  // hTest2->SetLineStyle(2);


  TH1D *hTest3 = new TH1D("hTest3", "hTest3", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hTest3->Add(hAdapTeO2po210M1, 0.8);
  hTest3->Add(hTest1, 0.1);
  hTest3->Add(hTest2, 0.1);

  hTest3->SetLineColor(4);
  // hTest3->SetLineStyle(4);

  hAdapTeO2po210M1->Draw();
  hTest1->Draw("SAME");
  hTest2->Draw("SAME");
  hTest3->Draw("SAME");

  TLegend *leg = new TLegend(0.8,0.8,0.97,0.97);
  leg->AddEntry(hAdapTeO2po210M1, "Original", "l");
  leg->AddEntry(hTest1, "Scaled Low", "l");
  leg->AddEntry(hTest2, "Scaled High", "l");
  leg->AddEntry(hTest3, "Summed", "l");
  leg->Draw();

}
  
bool TBackgroundModel::DoTheFitAdaptive()
{ 
  gStyle->SetOptStat(0);
   // This method actually sets up minuit and does the fit

  // Reduce Minuit Output
  minuit->SetPrintLevel(1);
  minuit->Command("SET STRategy 2");
  minuit->SetMaxIterations(10000);
  minuit->SetObjectFit(this); //see the external FCN  above

      //define the parameters and set the ranges and initial guesses see ROOTs TMinuit documentation

   // Close Th and Close Ra now split into its various sections, far Th and Ra still the same
   // This step after previous step full fit converges, just to see if any differences show up
   ////////////////////////////////////////////////
   // Using more parameters
   ////////////////////////////////////////////////
   minuit->DefineParameter(0, "TeO2 0nu",  80, 0.1, 0., 1000000);
   minuit->DefineParameter(1, "TeO2 2nu",  52000., 0.1, 0., 1000000);
   minuit->DefineParameter(2, "TeO2 co60",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(3, "TeO2 k40",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(4, "TeO2 pb210",  0, 0.1, 0., 1000000);
   minuit->DefineParameter(5, "TeO2 po210",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(6, "TeO2 te125",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(7, "TeO2 th232",  0, 0.1, 0., 1000000);
   minuit->DefineParameter(8, "TeO2 th228",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(9, "TeO2 ra226",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(10, "TeO2 rn222",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(11, "TeO2 u238",  0, 0.1, 0., 1000000);
   minuit->DefineParameter(12, "TeO2 th230",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(13, "TeO2 u234",  0., 0.1, 0., 1000000);

   minuit->DefineParameter(14, "CuFrame co58",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(15, "CuFrame co60",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(16, "CuFrame cs137",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(17, "CuFrame k40",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(18, "CuFrame mn54",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(19, "CuFrame pb210",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(20, "CuFrame th232",  0, 0.1, 0., 1000000);
   minuit->DefineParameter(21, "CuFrame u238",  0, 0.1, 0., 1000000);

   minuit->DefineParameter(22, "CuBoxco58",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(23, "CuBoxco60",  0, 0.1, 0., 1000000);
   minuit->DefineParameter(24, "CuBoxcs137",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(25, "CuBoxk40",  0, 0.1, 0., 1000000);
   minuit->DefineParameter(26, "CuBoxmn54",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(27, "CuBoxpb210",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(28, "CuBoxth232",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(29, "CuBoxu238",  0., 0.1, 0., 1000000);

   minuit->DefineParameter(30, "50mKco58",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(31, "50mKco60",  100, 0.1, 0., 1000000);
   minuit->DefineParameter(32, "50mKcs137",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(33, "50mKk40",  100, 0.1, 0., 1000000);
   minuit->DefineParameter(34, "50mKmn54",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(35, "50mKpb210",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(36, "50mKth232",  100, 0.1, 0., 1000000);
   minuit->DefineParameter(37, "50mKu238",  100, 0.1, 0., 1000000);

   minuit->DefineParameter(38, "600mKco60",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(39, "600mKk40",  0., 0.1, 0., 1000000); 
   minuit->DefineParameter(40, "600mKth232",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(41, "600mKu238",  0., 0.1, 0., 1000000);

   minuit->DefineParameter(42, "PbRombi207",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(43, "PbRomco60",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(44, "PbRomcs137",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(45, "PbRomk40",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(46, "PbRompb210",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(47, "PbRomth232",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(48, "PbRomu238",  0., 0.1, 0., 1000000);

   minuit->DefineParameter(49, "MBco60",  0., 0.1, 0., 1000000); 
   minuit->DefineParameter(50, "MBk40",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(51, "MBth232",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(52, "MBu238",  0., 0.1, 0., 1000000);

   minuit->DefineParameter(53, "IVCco60",  100, 0.1, 0., 1000000);
   minuit->DefineParameter(54, "IVCk40",  100, 0.1, 0., 1000000);
   minuit->DefineParameter(55, "IVCth232",  100, 0.1, 0., 1000000);
   minuit->DefineParameter(56, "IVCu238",  100, 0.1, 0., 1000000);

   minuit->DefineParameter(57, "OVCco60",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(58, "OVCk40",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(59, "OVCth232",  0., 0.1, 0., 1000000);    
   minuit->DefineParameter(60, "OVCu238",  0., 0.1, 0., 1000000);

   minuit->DefineParameter(61, "SIk40",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(62, "SIth232",  0., 0.1, 0., 1000000);    
   minuit->DefineParameter(63, "SIu238",  0., 0.1, 0., 1000000);

//////////////////////////////////////
   minuit->DefineParameter(64, "TeO2Spb210_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(65, "TeO2Spo210_001",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(66, "TeO2Spo210_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(67, "TeO2Sth232_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(68, "TeO2Su238_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(69, "TeO2Sxpb210_001",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(70, "TeO2Sxpb210_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(71, "TeO2Sxpb210_1",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(72, "TeO2Sxpb210_10",  0., 0.1, 0., 1000000);    
   minuit->DefineParameter(73, "TeO2Sxpo210_001",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(74, "TeO2Sxpo210_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(75, "TeO2Sxpo210_1",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(76, "TeO2Sxth232_001",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(77, "TeO2Sxth232_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(78, "TeO2Sxth232_1",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(79, "TeO2Sxth232_10",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(80, "TeO2Sxu238_001",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(81, "TeO2Sxu238_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(82, "TeO2Sxu238_1",  0., 0.1, 0., 1000000);   
   minuit->DefineParameter(83, "TeO2Sxu238_10",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(84, "CuFrameSth232_1",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(85, "CuFrameSu238_1",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(86, "CuFrameSxpb210_001",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(87, "CuFrameSxpb210_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(88, "CuFrameSxpb210_1",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(89, "CuFrameSxpb210_10",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(90, "CuFrameSxth232_001",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(91, "CuFrameSxth232_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(92, "CuFrameSxth232_1",  0., 0.1, 0., 1000000);   
   minuit->DefineParameter(93, "CuFrameSxth232_10",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(94, "CuFrameSxu238_001",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(95, "CuFrameSxu238_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(96, "CuFrameSxu238_1",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(97, "CuFrameSxu238_10",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(98, "CuBoxSth232_1",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(99, "CuBoxSu238_1",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(100, "CuBoxSxpb210_001",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(101, "CuBoxSxpb210_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(102, "CuBoxSxpb210_1",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(103, "CuBoxSxpb210_10",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(104, "CuBoxSxth232_001",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(105, "CuBoxSxth232_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(106, "CuBoxSxth232_1",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(107, "CuBoxSxth232_10",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(108, "CuBoxSxu238_001",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(109, "CuBoxSxu238_01",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(110, "CuBoxSxu238_1",  0., 0.1, 0., 1000000);
   minuit->DefineParameter(111, "CuBoxSxu238_10",  0., 0.1, 0., 1000000);

   minuit->DefineParameter(112, "Energy Scale High",  1.001., 0.00001, 1., 1.1);
   minuit->DefineParameter(113, "Energy Scale Low",  0.999., 0.00001, 0.9, 1);
   // minuit->DefineParameter(114, "Extra Smearing",         0., 0.001, 0., 5); // For alpha contamination...


//////////////////////////////////////

   // Fix parameters here
   // minuit->FixParameter(0); // TeO2 0nu
   // minuit->FixParameter(1); // TeO2 2nu
   // minuit->FixParameter(2); // TeO2 co60
   minuit->FixParameter(3); // TeO2 k40
   minuit->FixParameter(4); // TeO2 pb210
   // minuit->FixParameter(5); // TeO2 po210
   // minuit->FixParameter(6); // TeO2 te125m
   // minuit->FixParameter(7); // TeO2 th232
   minuit->FixParameter(8); // TeO2 th232-th228
   minuit->FixParameter(9); // TeO2 u238-ra226
   minuit->FixParameter(10); // TeO2 u238-rn222
   // minuit->FixParameter(11); // TeO2 u238
   minuit->FixParameter(12); // TeO2 u238-th230
   minuit->FixParameter(13); // TeO2 u238-u234
   minuit->FixParameter(14); // Frame co58
   // minuit->FixParameter(15); // Frame co60
   minuit->FixParameter(16); // Frame cs137
   minuit->FixParameter(17); // Frame k40
   minuit->FixParameter(18); // Frame mn54
   minuit->FixParameter(19); // Frame pb210
   // minuit->FixParameter(20); // Frame th232
   // minuit->FixParameter(21); // Frame u238
   minuit->FixParameter(22); // CuBox co58
   // minuit->FixParameter(23); // CuBox co60
   minuit->FixParameter(24); // CuBox cs137
   minuit->FixParameter(25); // CuBox k40
   minuit->FixParameter(26); // CuBox mn54
   minuit->FixParameter(27); // CuBox pb210
   // minuit->FixParameter(28); // CuBox th232
   // minuit->FixParameter(29); // CuBox u238
   minuit->FixParameter(30); // 50mK co58
   // minuit->FixParameter(31); // 50mK co60
   minuit->FixParameter(32); // 50mK cs137
   // minuit->FixParameter(33); // 50mK k40
   minuit->FixParameter(34); // 50mK mn54
   minuit->FixParameter(35); // 50mK pb210
   // minuit->FixParameter(36); // 50mK th232
   // minuit->FixParameter(37); // 50mK u238
   minuit->FixParameter(38); // 600mK co60
   minuit->FixParameter(39); // 600mK k40
   minuit->FixParameter(40); // 600mK th232   
   minuit->FixParameter(41); // 600mK u238
   minuit->FixParameter(42); // RLead bi207
   minuit->FixParameter(43); // RLead co60
   minuit->FixParameter(44); // RLead cs137
   minuit->FixParameter(45); // RLead k40
   minuit->FixParameter(46); // RLead pb210
   minuit->FixParameter(47); // RLead th232
   minuit->FixParameter(48); // RLead u238
   minuit->FixParameter(49); // MB co60
   minuit->FixParameter(50); // MB k40
   minuit->FixParameter(51); // MB th232
   minuit->FixParameter(52); // MB u238
   // minuit->FixParameter(53); // IVC co60
   // minuit->FixParameter(54); // IVC k40
   // minuit->FixParameter(55); // IVC th232
   // minuit->FixParameter(56); // IVC u238
   minuit->FixParameter(57); // OVC co60
   minuit->FixParameter(58); // OVC k40
   minuit->FixParameter(59); // OVC th232
   minuit->FixParameter(60); // OVC u238
   minuit->FixParameter(61); // SI k40
   minuit->FixParameter(62); // SI th232
   minuit->FixParameter(63); // SI u238   
   
   minuit->FixParameter(64); // TeO2 S pb210 01
   minuit->FixParameter(65); // TeO2 S po210 001
   minuit->FixParameter(66); // TeO2 S po210 01
   minuit->FixParameter(67); // TeO2 S th232 01
   minuit->FixParameter(68); // TeO2 S u238 01
   minuit->FixParameter(69); // TeO2 Sx pb210 001
   minuit->FixParameter(70); // TeO2 Sx pb210 01
   minuit->FixParameter(71); // TeO2 Sx pb210 1
   minuit->FixParameter(72); // TeO2 Sx pb210 10
   // minuit->FixParameter(73); // TeO2 Sx po210 001
   // minuit->FixParameter(74); // TeO2 Sx po210 01
   // minuit->FixParameter(75); // TeO2 Sx po210 1
   minuit->FixParameter(76); // TeO2 Sx th232 001
   minuit->FixParameter(77); // TeO2 Sx th232 01
   minuit->FixParameter(78); // TeO2 Sx th232 1
   minuit->FixParameter(79); // TeO2 Sx th232 10
   minuit->FixParameter(80); // TeO2 Sx u238 001
   minuit->FixParameter(81); // TeO2 Sx u238 01
   minuit->FixParameter(82); // TeO2 Sx u238 1
   minuit->FixParameter(83); // TeO2 Sx u238 10
   minuit->FixParameter(84); // Frame S th232 1
   minuit->FixParameter(85); // Frame S u238 1
   // minuit->FixParameter(86); // Frame Sx pb210 001
   // minuit->FixParameter(87); // Frame Sx pb210 01
   minuit->FixParameter(88); // Frame Sx pb210 1
   minuit->FixParameter(89); // Frame Sx pb210 10
   minuit->FixParameter(90); // Frame Sx th232 001
   minuit->FixParameter(91); // Frame Sx th232 01
   minuit->FixParameter(92); // Frame Sx th232 1
   minuit->FixParameter(93); // Frame Sx th232 10
   minuit->FixParameter(94); // Frame Sx u238 001
   minuit->FixParameter(95); // Frame Sx u238 01
   minuit->FixParameter(96); // Frame Sx u238 1
   minuit->FixParameter(97); // Frame Sx u238 10
   minuit->FixParameter(98); // CuBox S th232 1
   minuit->FixParameter(99); // CuBox S u238 1
   // minuit->FixParameter(100); // CuBox Sx pb210 001
   // minuit->FixParameter(101); // CuBox Sx pb210 01
   minuit->FixParameter(102); // CuBox Sx pb210 1
   minuit->FixParameter(103); // CuBox Sx pb210 10 
   minuit->FixParameter(104); // CuBox Sx th232 001
   minuit->FixParameter(105); // CuBox Sx th232 01
   minuit->FixParameter(106); // CuBox Sx th232 1
   minuit->FixParameter(107); // CuBox Sx th232 10
   minuit->FixParameter(108); // CuBox Sx u238 001
   minuit->FixParameter(109); // CuBox Sx u238 01
   minuit->FixParameter(110); // CuBox Sx u238 1
   minuit->FixParameter(111); // CuBox Sx u238 10
   minuit->FixParameter(112); // Energy scale factor High
   minuit->FixParameter(113); // Energy scale factor Low

   // Number of Parameters (for Chi-squared/NDF calculation)
   int dNumParameters = 28;
   //Tell minuit what external function to use 
   minuit->SetFCN(myExternal_FCNAdap);
   // int status = minuit->Migrad(); // this actually does the minimisation
   int status = minuit->Command("MINImize 100000 0.1"); // better method
   
  // Get final parameters from fit
  for(int i = 0; i < dNParam; i++)
  {
    minuit->GetParameter(i, fParameters[i], fParError[i]);
  }
  // Update model with final parameters
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
  fModelTotAdapthM1->Add( hAdapTeO2th232M1,     fParameters[7]  );
  fModelTotAdapthM1->Add( hAdapCuFrameth232M1,  fParameters[20] );
  fModelTotAdapthM1->Add( hAdapCuBoxth232M1,    fParameters[28] );
  fModelTotAdapthM1->Add( hAdap50mKth232M1,     fParameters[36] );
  fModelTotAdapthM1->Add( hAdap600mKth232M1,    fParameters[40] );
  fModelTotAdapthM1->Add( hAdapPbRomth232M1,    fParameters[47] );
  fModelTotAdapthM1->Add( hAdapMBth232M1,       fParameters[51] );
  fModelTotAdapthM1->Add( hAdapIVCth232M1,      fParameters[55] );
  fModelTotAdapthM1->Add( hAdapOVCth232M1,      fParameters[59] );
  fModelTotAdapthM1->Add( hAdapSIth232M1,       fParameters[62] );

  fModelTotAdapthM1->Add( hAdapTeO2Sxth232M1_001,      fParameters[76] );
  fModelTotAdapthM1->Add( hAdapTeO2Sxth232M1_01,       fParameters[77] );
  fModelTotAdapthM1->Add( hAdapTeO2Sxth232M1_1,        fParameters[78] );
  fModelTotAdapthM1->Add( hAdapTeO2Sxth232M1_10,       fParameters[79] );
  fModelTotAdapthM1->Add( hAdapCuFrameSxth232M1_001,   fParameters[90] );
  fModelTotAdapthM1->Add( hAdapCuFrameSxth232M1_01,    fParameters[91] );
  fModelTotAdapthM1->Add( hAdapCuFrameSxth232M1_1,     fParameters[92] );
  fModelTotAdapthM1->Add( hAdapCuFrameSxth232M1_10,    fParameters[93] );
  fModelTotAdapthM1->Add( hAdapCuBoxSxth232M1_001,     fParameters[104] );
  fModelTotAdapthM1->Add( hAdapCuBoxSxth232M1_01,      fParameters[105] );
  fModelTotAdapthM1->Add( hAdapCuBoxSxth232M1_1,       fParameters[106] );
  fModelTotAdapthM1->Add( hAdapCuBoxSxth232M1_10,      fParameters[107] );

  fModelTotAdapuM1->Add( hAdapTeO2u238M1,       fParameters[11] );
  fModelTotAdapuM1->Add( hAdapCuFrameu238M1,    fParameters[21] );
  fModelTotAdapuM1->Add( hAdapCuBoxu238M1,      fParameters[29] );
  fModelTotAdapuM1->Add( hAdap50mKu238M1,       fParameters[37] );
  fModelTotAdapuM1->Add( hAdap600mKu238M1,      fParameters[41] );
  fModelTotAdapuM1->Add( hAdapPbRomu238M1,      fParameters[48] );
  fModelTotAdapuM1->Add( hAdapMBu238M1,         fParameters[52] );
  fModelTotAdapuM1->Add( hAdapIVCu238M1,        fParameters[56] );
  fModelTotAdapuM1->Add( hAdapOVCu238M1,        fParameters[60] );
  fModelTotAdapuM1->Add( hAdapSIu238M1,         fParameters[63] );

  fModelTotAdapuM1->Add( hAdapTeO2Sxu238M1_001,        fParameters[80] );
  fModelTotAdapuM1->Add( hAdapTeO2Sxu238M1_01,         fParameters[81] );
  fModelTotAdapuM1->Add( hAdapTeO2Sxu238M1_1,          fParameters[82] );
  fModelTotAdapuM1->Add( hAdapTeO2Sxu238M1_10,         fParameters[83] );
  fModelTotAdapuM1->Add( hAdapCuFrameSxu238M1_001,     fParameters[94] );
  fModelTotAdapuM1->Add( hAdapCuFrameSxu238M1_01,      fParameters[95] );
  fModelTotAdapuM1->Add( hAdapCuFrameSxu238M1_1,       fParameters[96] );
  fModelTotAdapuM1->Add( hAdapCuFrameSxu238M1_10,      fParameters[97] );
  fModelTotAdapuM1->Add( hAdapCuBoxSxu238M1_001,       fParameters[108] );
  fModelTotAdapuM1->Add( hAdapCuBoxSxu238M1_01,        fParameters[109] );
  fModelTotAdapuM1->Add( hAdapCuBoxSxu238M1_1,         fParameters[110] );
  fModelTotAdapuM1->Add( hAdapCuBoxSxu238M1_10,        fParameters[111] );


  fModelTotAdapkM1->Add( hAdapTeO2k40M1,        fParameters[3]  );
  fModelTotAdapkM1->Add( hAdapCuFramek40M1,     fParameters[17] );
  fModelTotAdapkM1->Add( hAdapCuBoxk40M1,       fParameters[25] );
  fModelTotAdapkM1->Add( hAdap50mKk40M1,        fParameters[33] );
  fModelTotAdapkM1->Add( hAdap600mKk40M1,       fParameters[39] );
  fModelTotAdapkM1->Add( hAdapPbRomk40M1,       fParameters[45] );
  fModelTotAdapkM1->Add( hAdapMBk40M1,          fParameters[50] );
  fModelTotAdapkM1->Add( hAdapIVCk40M1,         fParameters[54] );
  fModelTotAdapkM1->Add( hAdapOVCk40M1,         fParameters[58] ); 
  fModelTotAdapkM1->Add( hAdapSIk40M1,          fParameters[61] );


  fModelTotAdapcoM1->Add( hAdapTeO2co60M1,      fParameters[2]  );
  fModelTotAdapcoM1->Add( hAdapCuFrameco60M1,   fParameters[15] );
  fModelTotAdapcoM1->Add( hAdapCuBoxco60M1,     fParameters[23] );
  fModelTotAdapcoM1->Add( hAdap50mKco60M1,      fParameters[31] );
  fModelTotAdapcoM1->Add( hAdap600mKco60M1,     fParameters[38] );
  fModelTotAdapcoM1->Add( hAdapPbRomco60M1,     fParameters[43] );
  fModelTotAdapcoM1->Add( hAdapMBco60M1,        fParameters[49] );
  fModelTotAdapcoM1->Add( hAdapIVCco60M1,       fParameters[53] );
  fModelTotAdapcoM1->Add( hAdapOVCco60M1,       fParameters[57] );

  fModelTotAdappbM1->Add( hAdapTeO2pb210M1,     fParameters[4]  );
  fModelTotAdappbM1->Add( hAdapCuFramepb210M1,  fParameters[19] );
  fModelTotAdappbM1->Add( hAdapCuBoxpb210M1,    fParameters[27] );
  fModelTotAdappbM1->Add( hAdap50mKpb210M1,     fParameters[35] );
  fModelTotAdappbM1->Add( hAdapPbRompb210M1,    fParameters[46] );

  fModelTotAdapcsM1->Add( hAdapCuFramecs137M1,  fParameters[16] );
  fModelTotAdapcsM1->Add( hAdapCuBoxcs137M1,    fParameters[24] );
  fModelTotAdapcsM1->Add( hAdap50mKcs137M1,     fParameters[32] );
  fModelTotAdapcsM1->Add( hAdapPbRomcs137M1,    fParameters[44] );

  fModelTotAdapcoM1->Add( hAdapCuFrameco58M1,   fParameters[14] );
  fModelTotAdapcoM1->Add( hAdapCuBoxco58M1,     fParameters[22] );
  fModelTotAdapcoM1->Add( hAdap50mKco58M1,      fParameters[30] );

  // fModelTotAdapteo2M1->Add( hAdapTeO2po210M1,   fParameters[5] );
  fModelTotAdapteo2M1->Add( hAdapTeO2te125M1,   fParameters[6] );
  fModelTotAdapteo2M1->Add( hAdapTeO2th228M1,   fParameters[8] );
  fModelTotAdapteo2M1->Add( hAdapTeO2ra226M1,   fParameters[9] );
  fModelTotAdapteo2M1->Add( hAdapTeO2rn222M1,   fParameters[10] );
  fModelTotAdapteo2M1->Add( hAdapTeO2th230M1,   fParameters[12] );
  fModelTotAdapteo2M1->Add( hAdapTeO2u234M1,    fParameters[13] );

  fModelTotAdapmnM1->Add( hAdapCuFramemn54M1,   fParameters[18] );
  fModelTotAdapmnM1->Add( hAdapCuBoxmn54M1,     fParameters[26] );
  fModelTotAdapmnM1->Add( hAdap50mKmn54M1,      fParameters[34] );

  fModelTotAdapbiM1->Add( hAdapPbRombi207M1,    fParameters[42] );
  fModelTotAdapNDBDM1->Add( hAdapTeO20nuM1,     fParameters[0] );
  fModelTotAdap2NDBDM1->Add( hAdapTeO22nuM1,    fParameters[1] );

  fModelTotAdapSpbM1->Add( hAdapTeO2Spb210M1_01,      fParameters[64] );
  fModelTotAdapSpbM1->Add( hAdapTeO2Sxpb210M1_001,    fParameters[69] );
  fModelTotAdapSpbM1->Add( hAdapTeO2Sxpb210M1_01,     fParameters[70] );
  fModelTotAdapSpbM1->Add( hAdapTeO2Sxpb210M1_1,      fParameters[71] );
  fModelTotAdapSpbM1->Add( hAdapTeO2Sxpb210M1_10,     fParameters[72] );
  fModelTotAdapSpbM1->Add( hAdapCuFrameSxpb210M1_001, fParameters[86] );
  fModelTotAdapSpbM1->Add( hAdapCuFrameSxpb210M1_01,  fParameters[87] );
  fModelTotAdapSpbM1->Add( hAdapCuFrameSxpb210M1_1,   fParameters[88] );
  fModelTotAdapSpbM1->Add( hAdapCuFrameSxpb210M1_10,  fParameters[89] );
  fModelTotAdapSpbM1->Add( hAdapCuBoxSxpb210M1_001,   fParameters[100] );
  fModelTotAdapSpbM1->Add( hAdapCuBoxSxpb210M1_01,    fParameters[101] );
  fModelTotAdapSpbM1->Add( hAdapCuBoxSxpb210M1_1,     fParameters[102] );
  fModelTotAdapSpbM1->Add( hAdapCuBoxSxpb210M1_10,    fParameters[103] );

  // Bulk Po added to here, because it's the same
  fModelTotAdapSpoM1->Add( hAdapTeO2po210M1,          fParameters[5] );
  fModelTotAdapSpoM1->Add( hAdapTeO2Spo210M1_001,     fParameters[65] );
  fModelTotAdapSpoM1->Add( hAdapTeO2Spo210M1_01,      fParameters[66] );
  fModelTotAdapSpoM1->Add( hAdapTeO2Sxpo210M1_001,    fParameters[73] );
  fModelTotAdapSpoM1->Add( hAdapTeO2Sxpo210M1_01,     fParameters[74] );
  fModelTotAdapSpoM1->Add( hAdapTeO2Sxpo210M1_1,      fParameters[75] );

  fModelTotAdapSthM1->Add( hAdapTeO2Sth232M1_01,      fParameters[67] );
  fModelTotAdapSthM1->Add( hAdapTeO2Sxth232M1_001,    fParameters[76] );
  fModelTotAdapSthM1->Add( hAdapTeO2Sxth232M1_01,     fParameters[77] );
  fModelTotAdapSthM1->Add( hAdapTeO2Sxth232M1_1,      fParameters[78] );
  fModelTotAdapSthM1->Add( hAdapTeO2Sxth232M1_10,     fParameters[79] );
  fModelTotAdapSthM1->Add( hAdapCuFrameSth232M1_1,    fParameters[84] );
  fModelTotAdapSthM1->Add( hAdapCuFrameSxth232M1_001, fParameters[90] );
  fModelTotAdapSthM1->Add( hAdapCuFrameSxth232M1_01,  fParameters[91] );
  fModelTotAdapSthM1->Add( hAdapCuFrameSxth232M1_1,   fParameters[92] );
  fModelTotAdapSthM1->Add( hAdapCuFrameSxth232M1_10,  fParameters[93] );
  fModelTotAdapSthM1->Add( hAdapCuBoxSth232M1_1,      fParameters[98] );
  fModelTotAdapSthM1->Add( hAdapCuBoxSxth232M1_001,   fParameters[104] );
  fModelTotAdapSthM1->Add( hAdapCuBoxSxth232M1_01,    fParameters[105] );
  fModelTotAdapSthM1->Add( hAdapCuBoxSxth232M1_1,     fParameters[106] );
  fModelTotAdapSthM1->Add( hAdapCuBoxSxth232M1_10,    fParameters[107] );

  fModelTotAdapSuM1->Add( hAdapTeO2Su238M1_01,        fParameters[68] );
  fModelTotAdapSuM1->Add( hAdapTeO2Sxu238M1_001,      fParameters[80] );
  fModelTotAdapSuM1->Add( hAdapTeO2Sxu238M1_01,       fParameters[81] );
  fModelTotAdapSuM1->Add( hAdapTeO2Sxu238M1_1,        fParameters[82] );
  fModelTotAdapSuM1->Add( hAdapTeO2Sxu238M1_10,       fParameters[83] );
  fModelTotAdapSuM1->Add( hAdapCuFrameSu238M1_1,      fParameters[85] );
  fModelTotAdapSuM1->Add( hAdapCuFrameSxu238M1_001,   fParameters[94] );
  fModelTotAdapSuM1->Add( hAdapCuFrameSxu238M1_01,    fParameters[95] );
  fModelTotAdapSuM1->Add( hAdapCuFrameSxu238M1_1,     fParameters[96] );
  fModelTotAdapSuM1->Add( hAdapCuFrameSxu238M1_10,    fParameters[97] );
  fModelTotAdapSuM1->Add( hAdapCuBoxSu238M1_1,        fParameters[99] );
  fModelTotAdapSuM1->Add( hAdapCuBoxSxu238M1_001,     fParameters[108] );
  fModelTotAdapSuM1->Add( hAdapCuBoxSxu238M1_01,      fParameters[109] );
  fModelTotAdapSuM1->Add( hAdapCuBoxSxu238M1_1,       fParameters[110] );
  fModelTotAdapSuM1->Add( hAdapCuBoxSxu238M1_10,      fParameters[111] );

// M2
  fModelTotAdapthM2->Add( hAdapTeO2th232M2,     fParameters[7]  );
  fModelTotAdapthM2->Add( hAdapCuFrameth232M2,  fParameters[20] );
  fModelTotAdapthM2->Add( hAdapCuBoxth232M2,    fParameters[28] );
  fModelTotAdapthM2->Add( hAdap50mKth232M2,     fParameters[36] );
  fModelTotAdapthM2->Add( hAdap600mKth232M2,    fParameters[40] );
  fModelTotAdapthM2->Add( hAdapPbRomth232M2,    fParameters[47] );
  fModelTotAdapthM2->Add( hAdapMBth232M2,       fParameters[51] );
  fModelTotAdapthM2->Add( hAdapIVCth232M2,      fParameters[55] );
  fModelTotAdapthM2->Add( hAdapOVCth232M2,      fParameters[59] );
  fModelTotAdapthM2->Add( hAdapSIth232M2,       fParameters[62] );

  fModelTotAdapthM2->Add( hAdapTeO2Sxth232M2_001,      fParameters[76] );
  fModelTotAdapthM2->Add( hAdapTeO2Sxth232M2_01,       fParameters[77] );
  fModelTotAdapthM2->Add( hAdapTeO2Sxth232M2_1,        fParameters[78] );
  fModelTotAdapthM2->Add( hAdapTeO2Sxth232M2_10,       fParameters[79] );
  fModelTotAdapthM2->Add( hAdapCuFrameSxth232M2_001,   fParameters[90] );
  fModelTotAdapthM2->Add( hAdapCuFrameSxth232M2_01,    fParameters[91] );
  fModelTotAdapthM2->Add( hAdapCuFrameSxth232M2_1,     fParameters[92] );
  fModelTotAdapthM2->Add( hAdapCuFrameSxth232M2_10,    fParameters[93] );
  fModelTotAdapthM2->Add( hAdapCuBoxSxth232M2_001,     fParameters[104] );
  fModelTotAdapthM2->Add( hAdapCuBoxSxth232M2_01,      fParameters[105] );
  fModelTotAdapthM2->Add( hAdapCuBoxSxth232M2_1,       fParameters[106] );
  fModelTotAdapthM2->Add( hAdapCuBoxSxth232M2_10,      fParameters[107] );


  fModelTotAdapuM2->Add( hAdapTeO2u238M2,       fParameters[11] );
  fModelTotAdapuM2->Add( hAdapCuFrameu238M2,    fParameters[21] );
  fModelTotAdapuM2->Add( hAdapCuBoxu238M2,      fParameters[29] );
  fModelTotAdapuM2->Add( hAdap50mKu238M2,       fParameters[37] );
  fModelTotAdapuM2->Add( hAdap600mKu238M2,      fParameters[41] );
  fModelTotAdapuM2->Add( hAdapPbRomu238M2,      fParameters[48] );
  fModelTotAdapuM2->Add( hAdapMBu238M2,         fParameters[52] );
  fModelTotAdapuM2->Add( hAdapIVCu238M2,        fParameters[56] );
  fModelTotAdapuM2->Add( hAdapOVCu238M2,        fParameters[60] );
  fModelTotAdapuM2->Add( hAdapSIu238M2,         fParameters[63] );

  fModelTotAdapuM2->Add( hAdapTeO2Sxu238M2_001,        fParameters[80] );
  fModelTotAdapuM2->Add( hAdapTeO2Sxu238M2_01,         fParameters[81] );
  fModelTotAdapuM2->Add( hAdapTeO2Sxu238M2_1,          fParameters[82] );
  fModelTotAdapuM2->Add( hAdapTeO2Sxu238M2_10,         fParameters[83] );
  fModelTotAdapuM2->Add( hAdapCuFrameSxu238M2_001,     fParameters[94] );
  fModelTotAdapuM2->Add( hAdapCuFrameSxu238M2_01,      fParameters[95] );
  fModelTotAdapuM2->Add( hAdapCuFrameSxu238M2_1,       fParameters[96] );
  fModelTotAdapuM2->Add( hAdapCuFrameSxu238M2_10,      fParameters[97] );
  fModelTotAdapuM2->Add( hAdapCuBoxSxu238M2_001,       fParameters[108] );
  fModelTotAdapuM2->Add( hAdapCuBoxSxu238M2_01,        fParameters[109] );
  fModelTotAdapuM2->Add( hAdapCuBoxSxu238M2_1,         fParameters[110] );
  fModelTotAdapuM2->Add( hAdapCuBoxSxu238M2_10,        fParameters[111] );


  fModelTotAdapkM2->Add( hAdapTeO2k40M2,        fParameters[3]  );
  fModelTotAdapkM2->Add( hAdapCuFramek40M2,     fParameters[17] );
  fModelTotAdapkM2->Add( hAdapCuBoxk40M2,       fParameters[25] );
  fModelTotAdapkM2->Add( hAdap50mKk40M2,        fParameters[33] );
  fModelTotAdapkM2->Add( hAdap600mKk40M2,       fParameters[39] );
  fModelTotAdapkM2->Add( hAdapPbRomk40M2,       fParameters[45] );
  fModelTotAdapkM2->Add( hAdapMBk40M2,          fParameters[50] );
  fModelTotAdapkM2->Add( hAdapIVCk40M2,         fParameters[54] );
  fModelTotAdapkM2->Add( hAdapOVCk40M2,         fParameters[58] ); 
  fModelTotAdapkM2->Add( hAdapSIk40M2,          fParameters[61] );


  fModelTotAdapcoM2->Add( hAdapTeO2co60M2,      fParameters[2]  );
  fModelTotAdapcoM2->Add( hAdapCuFrameco60M2,   fParameters[15] );
  fModelTotAdapcoM2->Add( hAdapCuBoxco60M2,     fParameters[23] );
  fModelTotAdapcoM2->Add( hAdap50mKco60M2,      fParameters[31] );
  fModelTotAdapcoM2->Add( hAdap600mKco60M2,     fParameters[38] );
  fModelTotAdapcoM2->Add( hAdapPbRomco60M2,     fParameters[43] );
  fModelTotAdapcoM2->Add( hAdapMBco60M2,        fParameters[49] );
  fModelTotAdapcoM2->Add( hAdapIVCco60M2,       fParameters[53] );
  fModelTotAdapcoM2->Add( hAdapOVCco60M2,       fParameters[57] );

  fModelTotAdappbM2->Add( hAdapTeO2pb210M2,     fParameters[4]  );
  fModelTotAdappbM2->Add( hAdapCuFramepb210M2,  fParameters[19] );
  fModelTotAdappbM2->Add( hAdapCuBoxpb210M2,    fParameters[27] );
  fModelTotAdappbM2->Add( hAdap50mKpb210M2,     fParameters[35] );
  fModelTotAdappbM2->Add( hAdapPbRompb210M2,    fParameters[46] );

  fModelTotAdapcsM2->Add( hAdapCuFramecs137M2,  fParameters[16] );
  fModelTotAdapcsM2->Add( hAdapCuBoxcs137M2,    fParameters[24] );
  fModelTotAdapcsM2->Add( hAdap50mKcs137M2,     fParameters[32] );
  fModelTotAdapcsM2->Add( hAdapPbRomcs137M2,    fParameters[44] );

  fModelTotAdapcoM2->Add( hAdapCuFrameco58M2,   fParameters[14] );
  fModelTotAdapcoM2->Add( hAdapCuBoxco58M2,     fParameters[22] );
  fModelTotAdapcoM2->Add( hAdap50mKco58M2,      fParameters[30] );

  // fModelTotAdapteo2M2->Add( hAdapTeO2po210M2,   fParameters[5] );
  fModelTotAdapteo2M2->Add( hAdapTeO2te125M2,   fParameters[6] );
  fModelTotAdapteo2M2->Add( hAdapTeO2th228M2,   fParameters[8] );
  fModelTotAdapteo2M2->Add( hAdapTeO2ra226M2,   fParameters[9] );
  fModelTotAdapteo2M2->Add( hAdapTeO2rn222M2,   fParameters[10] );
  fModelTotAdapteo2M2->Add( hAdapTeO2th230M2,   fParameters[12] );
  fModelTotAdapteo2M2->Add( hAdapTeO2u234M2,    fParameters[13] );

  fModelTotAdapmnM2->Add( hAdapCuFramemn54M2,   fParameters[18] );
  fModelTotAdapmnM2->Add( hAdapCuBoxmn54M2,     fParameters[26] );
  fModelTotAdapmnM2->Add( hAdap50mKmn54M2,      fParameters[34] );

  fModelTotAdapbiM2->Add( hAdapPbRombi207M2,    fParameters[42] );
  fModelTotAdapNDBDM2->Add( hAdapTeO20nuM2,     fParameters[0] );
  fModelTotAdap2NDBDM2->Add( hAdapTeO22nuM2,    fParameters[1] );

  fModelTotAdapSpbM2->Add( hAdapTeO2Spb210M2_01,      fParameters[64] );
  fModelTotAdapSpbM2->Add( hAdapTeO2Sxpb210M2_001,    fParameters[69] );
  fModelTotAdapSpbM2->Add( hAdapTeO2Sxpb210M2_01,     fParameters[70] );
  fModelTotAdapSpbM2->Add( hAdapTeO2Sxpb210M2_1,      fParameters[71] );
  fModelTotAdapSpbM2->Add( hAdapTeO2Sxpb210M2_10,     fParameters[72] );
  fModelTotAdapSpbM2->Add( hAdapCuFrameSxpb210M2_001, fParameters[86] );
  fModelTotAdapSpbM2->Add( hAdapCuFrameSxpb210M2_01,  fParameters[87] );
  fModelTotAdapSpbM2->Add( hAdapCuFrameSxpb210M2_1,   fParameters[88] );
  fModelTotAdapSpbM2->Add( hAdapCuFrameSxpb210M2_10,  fParameters[89] );
  fModelTotAdapSpbM2->Add( hAdapCuBoxSxpb210M2_001,   fParameters[100] );
  fModelTotAdapSpbM2->Add( hAdapCuBoxSxpb210M2_01,    fParameters[101] );
  fModelTotAdapSpbM2->Add( hAdapCuBoxSxpb210M2_1,     fParameters[102] );
  fModelTotAdapSpbM2->Add( hAdapCuBoxSxpb210M2_10,    fParameters[103] );

  fModelTotAdapSpoM2->Add( hAdapTeO2po210M2,          fParameters[5] );
  fModelTotAdapSpoM2->Add( hAdapTeO2Spo210M2_001,     fParameters[65] );
  fModelTotAdapSpoM2->Add( hAdapTeO2Spo210M2_01,      fParameters[66] );
  fModelTotAdapSpoM2->Add( hAdapTeO2Sxpo210M2_001,    fParameters[73] );
  fModelTotAdapSpoM2->Add( hAdapTeO2Sxpo210M2_01,     fParameters[74] );
  fModelTotAdapSpoM2->Add( hAdapTeO2Sxpo210M2_1,      fParameters[75] );

  fModelTotAdapSthM2->Add( hAdapTeO2Sth232M2_01,      fParameters[67] );
  fModelTotAdapSthM2->Add( hAdapTeO2Sxth232M2_001,    fParameters[76] );
  fModelTotAdapSthM2->Add( hAdapTeO2Sxth232M2_01,     fParameters[77] );
  fModelTotAdapSthM2->Add( hAdapTeO2Sxth232M2_1,      fParameters[78] );
  fModelTotAdapSthM2->Add( hAdapTeO2Sxth232M2_10,     fParameters[79] );
  fModelTotAdapSthM2->Add( hAdapCuFrameSth232M2_1,    fParameters[84] );
  fModelTotAdapSthM2->Add( hAdapCuFrameSxth232M2_001, fParameters[90] );
  fModelTotAdapSthM2->Add( hAdapCuFrameSxth232M2_01,  fParameters[91] );
  fModelTotAdapSthM2->Add( hAdapCuFrameSxth232M2_1,   fParameters[92] );
  fModelTotAdapSthM2->Add( hAdapCuFrameSxth232M2_10,  fParameters[93] );
  fModelTotAdapSthM2->Add( hAdapCuBoxSth232M2_1,      fParameters[98] );
  fModelTotAdapSthM2->Add( hAdapCuBoxSxth232M2_001,   fParameters[104] );
  fModelTotAdapSthM2->Add( hAdapCuBoxSxth232M2_01,    fParameters[105] );
  fModelTotAdapSthM2->Add( hAdapCuBoxSxth232M2_1,     fParameters[106] );
  fModelTotAdapSthM2->Add( hAdapCuBoxSxth232M2_10,    fParameters[107] );

  fModelTotAdapSuM2->Add( hAdapTeO2Su238M2_01,        fParameters[68] );
  fModelTotAdapSuM2->Add( hAdapTeO2Sxu238M2_001,      fParameters[80] );
  fModelTotAdapSuM2->Add( hAdapTeO2Sxu238M2_01,       fParameters[81] );
  fModelTotAdapSuM2->Add( hAdapTeO2Sxu238M2_1,        fParameters[82] );
  fModelTotAdapSuM2->Add( hAdapTeO2Sxu238M2_10,       fParameters[83] );
  fModelTotAdapSuM2->Add( hAdapCuFrameSu238M2_1,      fParameters[85] );
  fModelTotAdapSuM2->Add( hAdapCuFrameSxu238M2_001,   fParameters[94] );
  fModelTotAdapSuM2->Add( hAdapCuFrameSxu238M2_01,    fParameters[95] );
  fModelTotAdapSuM2->Add( hAdapCuFrameSxu238M2_1,     fParameters[96] );
  fModelTotAdapSuM2->Add( hAdapCuFrameSxu238M2_10,    fParameters[97] );
  fModelTotAdapSuM2->Add( hAdapCuBoxSu238M2_1,        fParameters[99] );
  fModelTotAdapSuM2->Add( hAdapCuBoxSxu238M2_001,     fParameters[108] );
  fModelTotAdapSuM2->Add( hAdapCuBoxSxu238M2_01,      fParameters[109] );
  fModelTotAdapSuM2->Add( hAdapCuBoxSxu238M2_1,       fParameters[110] );
  fModelTotAdapSuM2->Add( hAdapCuBoxSxu238M2_10,      fParameters[111] );



  TCanvas *cadap1 = new TCanvas("cadap1", "cadap1", 1200, 800);
  cadap1->SetLogy();

  ///// Draw Data M1
  fAdapDataHistoM1->SetLineColor(1);
  fAdapDataHistoM1->SetLineWidth(2);
  fAdapDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM1->GetYaxis()->SetTitle("Counts/Bin");
  fAdapDataHistoM1->SetMaximum(90000);
  // fAdapDataHistoM1->GetXaxis()->SetRange(1, fAdapDataHistoM1->FindBin(3000));
  fAdapDataHistoM1->Draw();


  fModelTotAdapM1->SetLineColor(2);
  fModelTotAdapM1->SetLineWidth(1);
  fModelTotAdapM1->Draw("SAME");

  fModelTotAdapthM1->SetLineColor(3);
  fModelTotAdapthM1->SetLineStyle(2);
  fModelTotAdapuM1->SetLineColor(4);
  fModelTotAdapuM1->SetLineStyle(2);
  fModelTotAdapkM1->SetLineColor(6);
  fModelTotAdapkM1->SetLineStyle(2);
  fModelTotAdapcoM1->SetLineColor(7);
  fModelTotAdapcoM1->SetLineStyle(2);
  fModelTotAdapNDBDM1->SetLineColor(42);
  fModelTotAdapNDBDM1->SetLineStyle(2);
  fModelTotAdap2NDBDM1->SetLineColor(46);
  fModelTotAdap2NDBDM1->SetLineStyle(2);
  // fModelTotAdapbiM1->SetLineColor(5);
  // fModelTotAdapbiM1->SetLineStyle(2);
  // fModelTotAdapmnM1->SetLineColor(40);
  // fModelTotAdapmnM1->SetLineStyle(2);

  fModelTotAdapSpoM1->SetLineStyle(2);
  fModelTotAdapSpoM1->SetLineColor(8);

  fModelTotAdapSpbM1->SetLineStyle(2);
  fModelTotAdapSpbM1->SetLineColor(40);

  fModelTotAdapteo2M1->SetLineStyle(2);
  fModelTotAdapteo2M1->SetLineColor(41);


  fModelTotAdapthM1->Draw("SAME");
  fModelTotAdapuM1->Draw("SAME");
  fModelTotAdapkM1->Draw("SAME");
  fModelTotAdapcoM1->Draw("SAME");
  fModelTotAdapNDBDM1->Draw("SAME");
  fModelTotAdap2NDBDM1->Draw("SAME");
  fModelTotAdapbiM1->Draw("SAME");
  fModelTotAdapmnM1->Draw("SAME");
  fModelTotAdapSpbM1->Draw("SAME");
  fModelTotAdapSpoM1->Draw("SAME");
  fModelTotAdappbM1->Draw("SAME");
  fModelTotAdapteo2M1->Draw("SAME");
  // fModelTotAdapPbM1->Draw("SAME");
/*
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
*/
  TLegend *legfit1 = new TLegend(0.8,0.8,0.97,0.97);
  legfit1->AddEntry(fModelTotAdapM1, "Total PDF", "l");
  legfit1->AddEntry(fModelTotAdapthM1, "Total th-232", "l");
  legfit1->AddEntry(fModelTotAdapuM1, "Total u-238", "l");
  legfit1->AddEntry(fModelTotAdapkM1, "Total k-40", "l");
  legfit1->AddEntry(fModelTotAdapcoM1, "Total co-60", "l");
  legfit1->AddEntry(fModelTotAdapNDBDM1, "NDBD", "l");
  legfit1->AddEntry(fModelTotAdap2NDBDM1, "2NDBD", "l");
  legfit1->AddEntry(fModelTotAdapSpoM1, "Total po-210", "l");
  legfit1->AddEntry(fModelTotAdapSpbM1, "Surface pb-210", "l");
  legfit1->AddEntry(fModelTotAdapteo2M1, "Other Crystal Bulk", "l");
  legfit1->Draw();





  TCanvas *cadap2 = new TCanvas("cadap2", "cadap2", 1200, 800);
  cadap2->SetLogy();

  ///// Draw Data M2
  fAdapDataHistoM2->SetLineColor(1);
  fAdapDataHistoM2->SetLineWidth(2);
  fAdapDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM2->GetYaxis()->SetTitle("Counts/Bin");
  fAdapDataHistoM2->SetMaximum(9000);
  // fAdapDataHistoM2->GetXaxis()->SetRange(1, fAdapDataHistoM2->FindBin(3000));
  fAdapDataHistoM2->Draw();

  
  fModelTotAdapM2->SetLineColor(2);
  fModelTotAdapM2->SetLineWidth(1);
  fModelTotAdapM2->Draw("SAME");

  fModelTotAdapthM2->SetLineColor(3);
  fModelTotAdapthM2->SetLineStyle(2);
  fModelTotAdapuM2->SetLineColor(4);
  fModelTotAdapuM2->SetLineStyle(2);
  fModelTotAdapkM2->SetLineColor(6);
  fModelTotAdapkM2->SetLineStyle(2);
  fModelTotAdapcoM2->SetLineColor(7);
  fModelTotAdapcoM2->SetLineStyle(2);
  fModelTotAdapNDBDM2->SetLineColor(42);
  fModelTotAdapNDBDM2->SetLineStyle(2);
  fModelTotAdap2NDBDM2->SetLineColor(46);
  fModelTotAdap2NDBDM2->SetLineStyle(2);
  fModelTotAdapbiM2->SetLineColor(5);
  fModelTotAdapbiM2->SetLineStyle(2);
  fModelTotAdapmnM2->SetLineColor(40);
  fModelTotAdapmnM2->SetLineStyle(2);

  // fTotCorrection->SetLineStyle(2);
  // fTotCorrection->SetLineColor(38);


  fModelTotAdapSpoM2->SetLineStyle(2);
  fModelTotAdapSpoM2->SetLineColor(8);

  fModelTotAdapSpbM2->SetLineStyle(2);
  fModelTotAdapSpbM2->SetLineColor(40);

  fModelTotAdapteo2M2->SetLineStyle(2);
  fModelTotAdapteo2M2->SetLineColor(41);

  fModelTotAdapthM2->Draw("SAME");
  fModelTotAdapuM2->Draw("SAME");
  fModelTotAdapkM2->Draw("SAME");
  fModelTotAdapcoM2->Draw("SAME");
  fModelTotAdapNDBDM2->Draw("SAME");
  fModelTotAdap2NDBDM2->Draw("SAME");
  // fModelTotAdapbiM2->Draw("SAME");    
  // fModelTotAdapmnM2->Draw("SAME"); 
  fModelTotAdapSpbM2->Draw("SAME");
  fModelTotAdapSpoM2->Draw("SAME");
  fModelTotAdappbM2->Draw("SAME");
  fModelTotAdapteo2M2->Draw("SAME");

/* 
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
*/

  TLegend *legfit2 = new TLegend(0.8,0.8,0.97,0.97);
  legfit2->AddEntry(fModelTotAdapM2, "Total PDF", "l");
  legfit2->AddEntry(fModelTotAdapthM2, "Total th-232", "l");
  legfit2->AddEntry(fModelTotAdapuM2, "Total u-238", "l");
  legfit2->AddEntry(fModelTotAdapkM2, "Total k-40", "l");
  legfit2->AddEntry(fModelTotAdapcoM2, "Total co-60", "l");
  legfit2->AddEntry(fModelTotAdapNDBDM2, "NDBD", "l");
  legfit2->AddEntry(fModelTotAdap2NDBDM2, "2NDBD", "l");
  legfit2->AddEntry(fModelTotAdapSpoM2, "Total po-210", "l");
  legfit2->AddEntry(fModelTotAdapSpbM2, "Surface pb-210", "l");
  legfit2->AddEntry(fModelTotAdapteo2M2, "Other Crystal Bulk", "l");

  legfit2->Draw();


  // Residuals
  TCanvas *cResidual1 = new TCanvas("cResidual1", "cResidual1", 1200, 800);
  hResidualGausM1 = new TH1D("hResidualGausM1", "Residual Distribution (M1)", 100, -50, 50);
  hResidualDistM1 = CalculateResidualsAdaptive(fModelTotAdapM1, fAdapDataHistoM1, hResidualGausM1, dFitMinBinM1, dFitMaxBinM1, 1);
  hResidualDistM1->SetLineColor(kBlack);
  hResidualDistM1->SetName("Residuals");
  hResidualDistM1->SetTitle("Fit Residuals (M1)");
  hResidualDistM1->SetMarkerStyle(25);
  hResidualDistM1->GetXaxis()->SetTitle("Energy (keV)");
  // hResidualDistM1->GetXaxis()->SetTitleSize(0.04);
  // hResidualDistM1->GetXaxis()->SetLabelSize(0.05);
  // hResidualDistM1->GetYaxis()->SetLabelSize(0.05);
  // hResidualDistM1->GetYaxis()->SetTitleSize(0.04);
  // hResidualDistM1->GetXaxis()->SetRange(1, fAdapDataHistoM2->FindBin(3000));
  hResidualDistM1->GetYaxis()->SetTitle("Residuals (#sigma)");
  hResidualDistM1->Draw("E");


  TCanvas *cres1 = new TCanvas("cres1");
  hResidualGausM1->Draw();



  TCanvas *cResidual2 = new TCanvas("cResidual2", "cResidual2", 1200, 800);
  hResidualGausM2 = new TH1D("hResidualGausM2", "Residual Distribution (M2)", 100, -50, 50);  
  hResidualDistM2 = CalculateResidualsAdaptive(fModelTotAdapM2, fAdapDataHistoM2, hResidualGausM2, dFitMinBinM2, dFitMaxBinM2, 2);
  hResidualDistM2->SetLineColor(kBlack);
  hResidualDistM2->SetName("Residuals");
  hResidualDistM2->SetTitle("Fit Residuals (M2)");
  hResidualDistM2->SetMarkerStyle(25);
  hResidualDistM2->GetXaxis()->SetTitle("Energy (keV)");
  // hResidualDistM2->GetXaxis()->SetTitleSize(0.04);
  // hResidualDistM2->GetXaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetTitleSize(0.04); 
  // hResidualDistM2->GetXaxis()->SetRange(1, fAdapDataHistoM2->FindBin(3000));
  hResidualDistM2->GetYaxis()->SetTitle("Residuals (#sigma)");

  hResidualDistM2->Draw("E");

  TCanvas *cres2 = new TCanvas("cres2");
  hResidualGausM2->Draw();



/*
  // Output integrals of stuff for limits
  cout << "ROI bin: " << fAdapDataHistoM1->FindBin(2470) << " " << fAdapDataHistoM1->FindBin(2570) << endl;
  cout << "Integral Data in ROI: " << fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt( fAdapDataHistoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total PDF in ROI: " << fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470)) << " +/- " << sqrt( fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total Th PDF in ROI: " << fModelTotAdapthM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt( fModelTotAdapthM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total Ra PDF in ROI: " << fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt(fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total Co PDF in ROI: " << fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt(fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total K PDF in ROI: " << fModelTotAdapkM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt(fModelTotAdapkM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total Bi PDF in ROI: " << fModelTotAdapbiM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt(fModelTotAdapbiM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;  
  cout << "Integral Total 2NDBD PDF in ROI: " << fModelTotAdap2NDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt(fModelTotAdap2NDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
  cout << "Integral Total 0NDBD PDF in ROI: " << fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) ) << " +/- " << sqrt(fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2470) )) << endl;
*/

  return true;

}


void TBackgroundModel::DrawMC()
{
// Draws all MC spectra, must Initialize first!

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend *legpb1 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legpb2 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legpo1 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legpo2 = new TLegend(0.65,0.7,0.95,0.95);

  TCanvas *cSPb2101 = new TCanvas("cSPb2101", "cSPb2101", 1200, 800);
  cSPb2101->SetLogy();

  hCuFrameSxpb210M1_001->SetLineColor(1);
  hCuFrameSxpb210M1_01->SetLineColor(2);
  hCuFrameSxpb210M1_1->SetLineColor(3);
  hCuFrameSxpb210M1_10->SetLineColor(4);
  hCuBoxSxpb210M1_001->SetLineColor(6);
  hCuBoxSxpb210M1_01->SetLineColor(7);
  hCuBoxSxpb210M1_1->SetLineColor(8);
  hCuBoxSxpb210M1_10->SetLineColor(9);

  hCuFrameSxpb210M1_001->GetXaxis()->SetTitle("Energy (keV)");
  hCuFrameSxpb210M1_001->GetYaxis()->SetTitle("Probability");  
  hCuFrameSxpb210M1_001->DrawNormalized();
  hCuFrameSxpb210M1_01->DrawNormalized("SAME");
  hCuFrameSxpb210M1_1->DrawNormalized("SAME");
  hCuFrameSxpb210M1_10->DrawNormalized("SAME");
  hCuBoxSxpb210M1_001->DrawNormalized("SAME");
  hCuBoxSxpb210M1_01->DrawNormalized("SAME");
  hCuBoxSxpb210M1_1->DrawNormalized("SAME");
  hCuBoxSxpb210M1_10->DrawNormalized("SAME");


  legpb1->AddEntry(hCuFrameSxpb210M1_001, "CuFrame Sx Pb210 0.01", "l");  
  legpb1->AddEntry(hCuFrameSxpb210M1_01, "CuFrame Sx Pb210 0.1", "l");
  legpb1->AddEntry(hCuFrameSxpb210M1_1, "CuFrame Sx Pb210 1", "l");
  legpb1->AddEntry(hCuFrameSxpb210M1_10, "CuFrame Sx Pb210 10", "l");
  legpb1->AddEntry(hCuBoxSxpb210M1_001, "CuBox S Pb210 0.01", "l");  
  legpb1->AddEntry(hCuBoxSxpb210M1_01, "CuBox Sx Pb210 0.1", "l");
  legpb1->AddEntry(hCuBoxSxpb210M1_1, "CuBox Sx Pb210 1", "l");
  legpb1->AddEntry(hCuBoxSxpb210M1_10, "CuBox Sx Pb210 10", "l");
  legpb1->Draw();


  TCanvas *cSPb2102 = new TCanvas("cSPb2102", "cSPb2102", 1200, 800);
  cSPb2102->SetLogy();

  hCuFrameSxpb210M2_001->SetLineColor(1);
  hCuFrameSxpb210M2_01->SetLineColor(2);
  hCuFrameSxpb210M2_1->SetLineColor(3);
  hCuFrameSxpb210M2_10->SetLineColor(4);
  hCuBoxSxpb210M2_001->SetLineColor(6);
  hCuBoxSxpb210M2_01->SetLineColor(7);
  hCuBoxSxpb210M2_1->SetLineColor(8);
  hCuBoxSxpb210M2_10->SetLineColor(9);

  hCuFrameSxpb210M2_001->GetXaxis()->SetTitle("Energy (keV)");
  hCuFrameSxpb210M2_001->GetYaxis()->SetTitle("Probability");  
  hCuFrameSxpb210M2_001->DrawNormalized();
  hCuFrameSxpb210M2_01->DrawNormalized("SAME");
  hCuFrameSxpb210M2_1->DrawNormalized("SAME");
  hCuFrameSxpb210M2_10->DrawNormalized("SAME");
  hCuBoxSxpb210M2_001->DrawNormalized("SAME");
  hCuBoxSxpb210M2_01->DrawNormalized("SAME");
  hCuBoxSxpb210M2_1->DrawNormalized("SAME");
  hCuBoxSxpb210M2_10->DrawNormalized("SAME");


  legpb2->AddEntry(hCuFrameSxpb210M2_001, "CuFrame Sx Pb210 0.01", "l");  
  legpb2->AddEntry(hCuFrameSxpb210M2_01, "CuFrame Sx Pb210 0.1", "l");
  legpb2->AddEntry(hCuFrameSxpb210M2_1, "CuFrame Sx Pb210 1", "l");
  legpb2->AddEntry(hCuFrameSxpb210M2_10, "CuFrame Sx Pb210 10", "l");
  legpb2->AddEntry(hCuBoxSxpb210M2_001, "CuBox S Pb210 0.01", "l");  
  legpb2->AddEntry(hCuBoxSxpb210M2_01, "CuBox Sx Pb210 0.1", "l");
  legpb2->AddEntry(hCuBoxSxpb210M2_1, "CuBox Sx Pb210 1", "l");
  legpb2->AddEntry(hCuBoxSxpb210M2_10, "CuBox Sx Pb210 10", "l");
  legpb2->Draw();

/*
  TCanvas *cPb2101 = new TCanvas("cPb2101", "cPb2101", 1200, 800);
  cPb2101->SetLogy();

  hTeO2Spb210M1_01->SetLineColor(1);
  hTeO2Sxpb210M1_10->SetLineColor(2);
  hTeO2Sxpb210M1_1->SetLineColor(3);
  hTeO2Sxpb210M1_01->SetLineColor(4);
  hTeO2Sxpb210M1_001->SetLineColor(6);

  hTeO2Spb210M1_01->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2Spb210M1_01->GetYaxis()->SetTitle("Probability");  
  hTeO2Spb210M1_01->DrawNormalized();
  hTeO2Sxpb210M1_10->DrawNormalized("SAME");
  hTeO2Sxpb210M1_1->DrawNormalized("SAME");
  hTeO2Sxpb210M1_01->DrawNormalized("SAME");
  hTeO2Sxpb210M1_001->DrawNormalized("SAME");

  legpb1->AddEntry(hTeO2Spb210M1_01, "TeO2 S Pb210 0.1", "l");  
  legpb1->AddEntry(hTeO2Sxpb210M1_10, "TeO2 Sx Pb210 10", "l");
  legpb1->AddEntry(hTeO2Sxpb210M1_1, "TeO2 Sx Pb210 1", "l");
  legpb1->AddEntry(hTeO2Sxpb210M1_01, "TeO2 Sx Pb210 0.1", "l");
  legpb1->AddEntry(hTeO2Sxpb210M1_001, "TeO2 Sx Pb210 0.01", "l");
  legpb1->Draw();

  TCanvas *cPb2102 = new TCanvas("cPb2102", "cPb2102", 1200, 800);
  cPb2102->SetLogy();

  hTeO2Spb210M2_01->SetLineColor(1);
  hTeO2Sxpb210M2_10->SetLineColor(2);
  hTeO2Sxpb210M2_1->SetLineColor(3);
  hTeO2Sxpb210M2_01->SetLineColor(4);
  hTeO2Sxpb210M2_001->SetLineColor(6);

  hTeO2Spb210M2_01->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2Spb210M2_01->GetYaxis()->SetTitle("Probability");  
  hTeO2Spb210M2_01->DrawNormalized();
  hTeO2Sxpb210M2_10->DrawNormalized("SAME");
  hTeO2Sxpb210M2_1->DrawNormalized("SAME");
  hTeO2Sxpb210M2_01->DrawNormalized("SAME");
  hTeO2Sxpb210M2_001->DrawNormalized("SAME");

  legpb2->AddEntry(hTeO2Spb210M2_01, "TeO2 S Pb210 0.1", "l");  
  legpb2->AddEntry(hTeO2Sxpb210M2_10, "TeO2 Sx Pb210 10", "l");
  legpb2->AddEntry(hTeO2Sxpb210M2_1, "TeO2 Sx Pb210 1", "l");
  legpb2->AddEntry(hTeO2Sxpb210M2_01, "TeO2 Sx Pb210 0.1", "l");
  legpb2->AddEntry(hTeO2Sxpb210M2_001, "TeO2 Sx Pb210 0.01", "l");
  legpb2->Draw();



  TCanvas *cPo2101 = new TCanvas("cPo2101", "cPo2101", 1200, 800);
  cPo2101->SetLogy();

  hTeO2Spo210M1_01->SetLineColor(1);
  hTeO2Spo210M1_001->SetLineColor(2);
  hTeO2Sxpo210M1_1->SetLineColor(3);
  hTeO2Sxpo210M1_01->SetLineColor(4);
  hTeO2Sxpo210M1_001->SetLineColor(6);

  hTeO2Spo210M1_01->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2Spo210M1_01->GetYaxis()->SetTitle("Probability");  
  hTeO2Spo210M1_01->DrawNormalized();
  hTeO2Spo210M1_001->DrawNormalized("SAME");
  hTeO2Sxpo210M1_1->DrawNormalized("SAME");
  hTeO2Sxpo210M1_01->DrawNormalized("SAME");
  hTeO2Sxpo210M1_001->DrawNormalized("SAME");

  legpo1->AddEntry(hTeO2Spo210M1_01, "TeO2 S po210 0.1", "l");  
  legpo1->AddEntry(hTeO2Spo210M1_001, "TeO2 S po210 0.01", "l");  
  legpo1->AddEntry(hTeO2Sxpo210M1_1, "TeO2 Sx po210 1", "l");
  legpo1->AddEntry(hTeO2Sxpo210M1_01, "TeO2 Sx po210 0.1", "l");
  legpo1->AddEntry(hTeO2Sxpo210M1_001, "TeO2 Sx po210 0.01", "l");
  legpo1->Draw();


  TCanvas *cPo2102 = new TCanvas("cPo2102", "cPo2102", 1200, 800);
  cPo2102->SetLogy();

  hTeO2Spo210M2_01->SetLineColor(1);
  hTeO2Spo210M2_001->SetLineColor(2);
  hTeO2Sxpo210M2_1->SetLineColor(3);
  hTeO2Sxpo210M2_01->SetLineColor(4);
  hTeO2Sxpo210M2_001->SetLineColor(6);

  hTeO2Spo210M2_01->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2Spo210M2_01->GetYaxis()->SetTitle("Probability");  
  hTeO2Spo210M2_01->DrawNormalized();
  hTeO2Spo210M2_001->DrawNormalized("SAME");
  hTeO2Sxpo210M2_1->DrawNormalized("SAME");
  hTeO2Sxpo210M2_01->DrawNormalized("SAME");
  hTeO2Sxpo210M2_001->DrawNormalized("SAME");

  legpo2->AddEntry(hTeO2Spo210M2_01, "TeO2 S po210 0.1", "l");  
  legpo2->AddEntry(hTeO2Spo210M2_001, "TeO2 S po210 0.01", "l");  
  legpo2->AddEntry(hTeO2Sxpo210M2_1, "TeO2 Sx po210 1", "l");
  legpo2->AddEntry(hTeO2Sxpo210M2_01, "TeO2 Sx po210 0.1", "l");
  legpo2->AddEntry(hTeO2Sxpo210M2_001, "TeO2 Sx po210 0.01", "l");
  legpo2->Draw();
*/

/*
  TLegend *legth1 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legu1 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legk1 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legco1 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legndbd1 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legbi1 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legths1 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legus1 = new TLegend(0.65,0.7,0.95,0.95);

  TLegend *legth2 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legu2 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legk2 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legco2 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legndbd2 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legbi2 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legths2 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legus2 = new TLegend(0.65,0.7,0.95,0.95);


  TCanvas *cTh2321 = new TCanvas("cTh2321", "cTh2321", 1200, 800);
  cTh2321->SetLogy();

  hTeO2th232M1->SetLineColor(1);
  hCuFrameth232M1->SetLineColor(2);
  hCuBoxth232M1->SetLineColor(3);
  h50mKth232M1->SetLineColor(4);
  h600mKth232M1->SetLineColor(5);
  hIVCth232M1->SetLineColor(6);
  hOVCth232M1->SetLineColor(7);
  hPbRomth232M1->SetLineColor(8);
  hMBth232M1->SetLineColor(9);
  hSIth232M1->SetLineColor(11);

  hTeO2th232M1->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2th232M1->GetYaxis()->SetTitle("Probability");  
  hTeO2th232M1->DrawNormalized();
  hCuFrameth232M1->DrawNormalized("SAME");
  hCuBoxth232M1->DrawNormalized("SAME");
  h50mKth232M1->DrawNormalized("SAME");
  h600mKth232M1->DrawNormalized("SAME");
  hIVCth232M1->DrawNormalized("SAME");
  hOVCth232M1->DrawNormalized("SAME");
  hPbRomth232M1->DrawNormalized("SAME");
  hMBth232M1->DrawNormalized("SAME");
  hSIth232M1->DrawNormalized("SAME");

  legth1->AddEntry(hTeO2th232M1, "TeO2", "l");  
  legth1->AddEntry(hCuFrameth232M1, "CuFrame", "l");
  legth1->AddEntry(hCuBoxth232M1, "CuBox", "l");
  legth1->AddEntry(h50mKth232M1, "50mK", "l");
  legth1->AddEntry(h600mKth232M1, "600mK", "l");
  legth1->AddEntry(hIVCth232M1, "IVC", "l");
  legth1->AddEntry(hOVCth232M1, "OVC", "l");
  legth1->AddEntry(hPbRomth232M1, "PbRom", "l");
  legth1->AddEntry(hMBth232M1, "MB", "l");
  legth1->AddEntry(hSIth232M1, "SI", "l");
  legth1->Draw();


  TCanvas *cTh2322 = new TCanvas("cTh2322", "cTh2322", 1200, 800);
  cTh2322->SetLogy();

  hTeO2th232M2->SetLineColor(1);
  hCuFrameth232M2->SetLineColor(2);
  hCuBoxth232M2->SetLineColor(3);
  h50mKth232M2->SetLineColor(4);
  h600mKth232M2->SetLineColor(5);
  hIVCth232M2->SetLineColor(6);
  hOVCth232M2->SetLineColor(7);
  hPbRomth232M2->SetLineColor(8);
  hMBth232M2->SetLineColor(9);
  hSIth232M2->SetLineColor(11);

  hTeO2th232M2->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2th232M2->GetYaxis()->SetTitle("Probability"); 
  hTeO2th232M2->DrawNormalized();   
  hCuFrameth232M2->DrawNormalized("SAME");
  hCuBoxth232M2->DrawNormalized("SAME");
  h50mKth232M2->DrawNormalized("SAME");
  h600mKth232M2->DrawNormalized("SAME");
  hIVCth232M2->DrawNormalized("SAME");
  hOVCth232M2->DrawNormalized("SAME");
  hPbRomth232M2->DrawNormalized("SAME");
  hMBth232M2->DrawNormalized("SAME");
  hSIth232M2->DrawNormalized("SAME");

  legth2->AddEntry(hTeO2th232M2, "TeO2", "l");  
  legth2->AddEntry(hCuFrameth232M2, "CuFrame", "l");
  legth2->AddEntry(hCuBoxth232M2, "CuBox", "l");
  legth2->AddEntry(h50mKth232M2, "50mK", "l");
  legth2->AddEntry(h600mKth232M2, "600mK", "l");
  legth2->AddEntry(hIVCth232M2, "IVC", "l");
  legth2->AddEntry(hOVCth232M2, "OVC", "l");
  legth2->AddEntry(hPbRomth232M2, "PbRom", "l");
  legth2->AddEntry(hMBth232M2, "MB", "l");
  legth2->AddEntry(hSIth232M2, "SI", "l");
  legth2->Draw();

  TCanvas *cU2381 = new TCanvas("cU2381", "cU2381", 1200, 800);
  cU2381->SetLogy();

  hTeO2u238M1->SetLineColor(1);
  hCuFrameu238M1->SetLineColor(2);
  hCuBoxu238M1->SetLineColor(3);
  h50mKu238M1->SetLineColor(4);
  h600mKu238M1->SetLineColor(5);
  hIVCu238M1->SetLineColor(6);
  hOVCu238M1->SetLineColor(7);
  hPbRomu238M1->SetLineColor(8);
  hMBu238M1->SetLineColor(9);
  hSIu238M1->SetLineColor(11);

  hTeO2u238M1->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2u238M1->GetYaxis()->SetTitle("Probability");  
  hTeO2u238M1->DrawNormalized();
  hCuFrameu238M1->DrawNormalized("SAME");
  hCuBoxu238M1->DrawNormalized("SAME");
  h50mKu238M1->DrawNormalized("SAME");
  h600mKu238M1->DrawNormalized("SAME");
  hIVCu238M1->DrawNormalized("SAME");
  hOVCu238M1->DrawNormalized("SAME");
  hPbRomu238M1->DrawNormalized("SAME");
  hMBu238M1->DrawNormalized("SAME");
  hSIu238M1->DrawNormalized("SAME");

  legu1->AddEntry(hTeO2u238M2, "TeO2", "l");  
  legu1->AddEntry(hCuFrameu238M2, "CuFrame", "l");
  legu1->AddEntry(hCuBoxu238M2, "CuBox", "l");
  legu1->AddEntry(h50mKu238M2, "50mK", "l");
  legu1->AddEntry(h600mKu238M2, "600mK", "l");
  legu1->AddEntry(hIVCu238M2, "IVC", "l");
  legu1->AddEntry(hOVCu238M2, "OVC", "l");
  legu1->AddEntry(hPbRomu238M2, "PbRom", "l");
  legu1->AddEntry(hMBu238M2, "MB", "l");
  legu1->AddEntry(hSIu238M2, "SI", "l");
  legu1->Draw();


  TCanvas *cU2382 = new TCanvas("cU2382", "cU2382", 1200, 800);
  cU2382->SetLogy();

  hTeO2u238M2->SetLineColor(1);
  hCuFrameu238M2->SetLineColor(2);
  hCuBoxu238M2->SetLineColor(3);
  h50mKu238M2->SetLineColor(4);
  h600mKu238M2->SetLineColor(5);
  hIVCu238M2->SetLineColor(6);
  hOVCu238M2->SetLineColor(7);
  hPbRomu238M2->SetLineColor(8);
  hMBu238M2->SetLineColor(9);
  hSIu238M2->SetLineColor(11);

  hTeO2u238M2->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2u238M2->GetYaxis()->SetTitle("Probability"); 
  hTeO2u238M2->DrawNormalized();   
  hCuFrameu238M2->DrawNormalized("SAME");
  hCuBoxu238M2->DrawNormalized("SAME");
  h50mKu238M2->DrawNormalized("SAME");
  h600mKu238M2->DrawNormalized("SAME");
  hIVCu238M2->DrawNormalized("SAME");
  hOVCu238M2->DrawNormalized("SAME");
  hPbRomu238M2->DrawNormalized("SAME");
  hMBu238M2->DrawNormalized("SAME");
  hSIu238M2->DrawNormalized("SAME");

  legu2->AddEntry(hTeO2u238M2, "TeO2", "l");  
  legu2->AddEntry(hCuFrameu238M2, "CuFrame", "l");
  legu2->AddEntry(hCuBoxu238M2, "CuBox", "l");
  legu2->AddEntry(h50mKu238M2, "50mK", "l");
  legu2->AddEntry(h600mKu238M2, "600mK", "l");
  legu2->AddEntry(hIVCu238M2, "IVC", "l");
  legu2->AddEntry(hOVCu238M2, "OVC", "l");
  legu2->AddEntry(hPbRomu238M2, "PbRom", "l");
  legu2->AddEntry(hMBu238M2, "MB", "l");
  legu2->AddEntry(hSIu238M2, "SI", "l");
  legu2->Draw();




  TCanvas *cK401 = new TCanvas("cK401", "cK401", 1200, 800);
  cK401->SetLogy();

  hTeO2k40M1->SetLineColor(1);
  hCuFramek40M1->SetLineColor(2);
  hCuBoxk40M1->SetLineColor(3);
  h50mKk40M1->SetLineColor(4);
  h600mKk40M1->SetLineColor(5);
  hIVCk40M1->SetLineColor(6);
  hOVCk40M1->SetLineColor(7);
  hPbRomk40M1->SetLineColor(8);
  hMBk40M1->SetLineColor(9);
  hSIk40M1->SetLineColor(11);

  hTeO2k40M1->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2k40M1->GetYaxis()->SetTitle("Probability");  
  hTeO2k40M1->DrawNormalized();
  hCuFramek40M1->DrawNormalized("SAME");
  hCuBoxk40M1->DrawNormalized("SAME");
  h50mKk40M1->DrawNormalized("SAME");
  h600mKk40M1->DrawNormalized("SAME");
  hIVCk40M1->DrawNormalized("SAME");
  hOVCk40M1->DrawNormalized("SAME");
  hPbRomk40M1->DrawNormalized("SAME");
  hMBk40M1->DrawNormalized("SAME");
  hSIk40M1->DrawNormalized("SAME");

  legk1->AddEntry(hTeO2k40M1, "TeO2", "l");  
  legk1->AddEntry(hCuFramek40M1, "CuFrame", "l");
  legk1->AddEntry(hCuBoxk40M1, "CuBox", "l");
  legk1->AddEntry(h50mKk40M1, "50mK", "l");
  legk1->AddEntry(h600mKk40M1, "600mK", "l");
  legk1->AddEntry(hIVCk40M1, "IVC", "l");
  legk1->AddEntry(hOVCk40M1, "OVC", "l");
  legk1->AddEntry(hPbRomk40M1, "PbRom", "l");
  legk1->AddEntry(hMBk40M1, "MB", "l");
  legk1->AddEntry(hSIk40M1, "SI", "l");
  legk1->Draw();



  TCanvas *cK402 = new TCanvas("cK402", "cK402", 1200, 800);
  cK402->SetLogy();

  hTeO2k40M2->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2k40M2->GetYaxis()->SetTitle("Probability");  
  hTeO2k40M2->SetLineColor(1);
  hCuFramek40M2->SetLineColor(2);
  hCuBoxk40M2->SetLineColor(3);
  h50mKk40M2->SetLineColor(4);
  h600mKk40M2->SetLineColor(5);
  hIVCk40M2->SetLineColor(6);
  hOVCk40M2->SetLineColor(7);
  hPbRomk40M2->SetLineColor(8);
  hMBk40M2->SetLineColor(9);
  hSIk40M2->SetLineColor(11);


  hTeO2k40M2->DrawNormalized();
  hCuFramek40M2->DrawNormalized("SAME");
  hCuBoxk40M2->DrawNormalized("SAME");
  h50mKk40M2->DrawNormalized("SAME");
  h600mKk40M2->DrawNormalized("SAME");
  hIVCk40M2->DrawNormalized("SAME");
  hOVCk40M2->DrawNormalized("SAME");
  hPbRomk40M2->DrawNormalized("SAME");
  hMBk40M2->DrawNormalized("SAME");
  hSIk40M2->DrawNormalized("SAME");

  legk2->AddEntry(hTeO2k40M2, "TeO2", "l");  
  legk2->AddEntry(hCuFramek40M2, "CuFrame", "l");
  legk2->AddEntry(hCuBoxk40M2, "CuBox", "l");
  legk2->AddEntry(h50mKk40M2, "50mK", "l");
  legk2->AddEntry(h600mKk40M2, "600mK", "l");
  legk2->AddEntry(hIVCk40M2, "IVC", "l");
  legk2->AddEntry(hOVCk40M2, "OVC", "l");
  legk2->AddEntry(hPbRomk40M2, "PbRom", "l");
  legk2->AddEntry(hMBk40M2, "MB", "l");
  legk2->AddEntry(hSIk40M2, "SI", "l");
  legk2->Draw();


  TCanvas *cCo601 = new TCanvas("cCo601", "cCo601", 1200, 800);
  cCo601->SetLogy();

  hTeO2co60M1->SetLineColor(1);
  hCuFrameco60M1->SetLineColor(2);
  hCuBoxco60M1->SetLineColor(3);
  h50mKco60M1->SetLineColor(4);
  h600mKco60M1->SetLineColor(5);
  hIVCco60M1->SetLineColor(6);
  hOVCco60M1->SetLineColor(7);
  hPbRomco60M1->SetLineColor(8);
  hMBco60M1->SetLineColor(9);

  hTeO2co60M1->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2co60M1->GetYaxis()->SetTitle("Probability");  
  hTeO2co60M1->DrawNormalized();
  hCuFrameco60M1->DrawNormalized("SAME");
  hCuBoxco60M1->DrawNormalized("SAME");
  h50mKco60M1->DrawNormalized("SAME");
  h600mKco60M1->DrawNormalized("SAME");
  hIVCco60M1->DrawNormalized("SAME");
  hOVCco60M1->DrawNormalized("SAME");
  hPbRomco60M1->DrawNormalized("SAME");
  hMBco60M1->DrawNormalized("SAME");

  legco1->AddEntry(hTeO2co60M1, "TeO2", "l");  
  legco1->AddEntry(hCuFrameco60M1, "CuFrame", "l");
  legco1->AddEntry(hCuBoxco60M1, "CuBox", "l");
  legco1->AddEntry(h50mKco60M1, "50mK", "l");
  legco1->AddEntry(h600mKco60M1, "600mK", "l");
  legco1->AddEntry(hIVCco60M1, "IVC", "l");
  legco1->AddEntry(hOVCco60M1, "OVC", "l");
  legco1->AddEntry(hPbRomco60M1, "PbRom", "l");
  legco1->AddEntry(hMBco60M1, "MB", "l");
  legco1->Draw();



  TCanvas *cCo602 = new TCanvas("cCo602", "cCo602", 1200, 800);
  cCo602->SetLogy();

  hTeO2co60M2->SetLineColor(1);
  hCuFrameco60M2->SetLineColor(2);
  hCuBoxco60M2->SetLineColor(3);
  h50mKco60M2->SetLineColor(4);
  h600mKco60M2->SetLineColor(5);
  hIVCco60M2->SetLineColor(6);
  hOVCco60M2->SetLineColor(7);
  hPbRomco60M2->SetLineColor(8);
  hMBco60M2->SetLineColor(9);

  hTeO2co60M2->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2co60M2->GetYaxis()->SetTitle("Probability");  
  hTeO2co60M2->DrawNormalized();
  hCuFrameco60M2->DrawNormalized("SAME");
  hCuBoxco60M2->DrawNormalized("SAME");
  h50mKco60M2->DrawNormalized("SAME");
  h600mKco60M2->DrawNormalized("SAME");
  hIVCco60M2->DrawNormalized("SAME");
  hOVCco60M2->DrawNormalized("SAME");
  hPbRomco60M2->DrawNormalized("SAME");
  hMBco60M2->DrawNormalized("SAME");

  legco2->AddEntry(hTeO2co60M2, "TeO2", "l");  
  legco2->AddEntry(hCuFrameco60M2, "CuFrame", "l");
  legco2->AddEntry(hCuBoxco60M2, "CuBox", "l");
  legco2->AddEntry(h50mKco60M2, "50mK", "l");
  legco2->AddEntry(h600mKco60M2, "600mK", "l");
  legco2->AddEntry(hIVCco60M2, "IVC", "l");
  legco2->AddEntry(hOVCco60M2, "OVC", "l");
  legco2->AddEntry(hPbRomco60M2, "PbRom", "l");
  legco2->AddEntry(hMBco60M2, "MB", "l");
  legco2->Draw();
*/
}

// For custom Smearing with resolution, currently constant resolution 
TH1D *TBackgroundModel::SmearMC(TH1D *hMC, TH1D *hSMC, double resolution1, int dMult)
{
  // Reset previously Modeled histogram
  hSMC->Reset();

  double dNorm;
  double dSmearedValue;

  for(dMult == 1)
  {
  for(int i = 0; i < dAdaptiveBinsM1 ; i++)
  {
    for(int j = 0; j < dAdaptiveBinsM1; j++)
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
  }

  for(dMult == 2)
  {
  for(int i = 0; i < dAdaptiveBinsM2 ; i++)
  {
    for(int j = 0; j < dAdaptiveBinsM2; j++)
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
  }

  return hSMC;
}

// Only use after fit
// Need maybe an array with all of the names of the variables
void TBackgroundModel::LatexResultTable()
{
  ofstream OutFile;
  OutFile.open("OutputTable.txt");

  for(int i = 0; i < dNParam; i++)
  {
    OutFile << "Name & " << fParameters[i] <<  " $\pm$" << fParError[i] << "\\" << endl;
  }
  OutFile.wcite();

}


