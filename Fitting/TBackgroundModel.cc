// The save option crashes when exiting for some reason...
// #include "TMinuit.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TBackgroundModel.hh"
#include "TApplication.h"
#include "TRandom3.h"
#include "TAxis.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TPave.h"
#include "TMath.h"
#include <cmath>
#include <string>
#include <vector>


using namespace std;
using std::cout;
using std::endl;
using std::map;
using std::vector;

ClassImp(TBackgroundModel)

//first set up a global function that calls your classes method that calculates the quantity to minimize
void myExternal_FCNAdap(int &n, double *grad, double &fval, double x[], int code)
{
  // Required External Wrapper for function to be minimized by Minuit 
 
  // This gets called for each value of the parameters minuit tries
  // here the x array contains the parameters you are trying to fit
  
  // here myClass should inherit from TObject
  TBackgroundModel* Obj = (TBackgroundModel*)gMinuit->GetObjectFit(); 

  // implement a method in your class for setting the parameters and thus update the parameters of your fitter class 
  for(int i = 0; i < Obj->dNParam; i++ )
  {
    Obj->SetParameters(i, x[i]);
  }
  // Implement a method in your class that calculates the quantity you want to minimize, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimize fval
    Obj->UpdateModelAdaptive();
    fval = Obj->GetChiSquareAdaptive();
}

TBackgroundModel::TBackgroundModel()
{}

TBackgroundModel::TBackgroundModel(double fFitMin, double fFitMax, int fBinBase, int fDataSet, bool fSave)
{
  // gSystem->Load("TBkgModelParameter_cc.so");
  bSave = fSave;
  bToyData = false;

  tTime = new TDatime();
  // dNParam = 59; // number of fitting parameters
  dNParam = 28; // number of fitting parameters (reduced)
  dNumCalls = 0;
  dSecToYears = 1./(60*60*24*365);
  dMass = 36.75;

  // Data directories depending on QCC/local
  dDataDir =  "/Users/brian/macros/CUOREZ/Bkg";
   // dDataDir = "/cuore/user/zhubrian/CUORE0/scratch/";
  dMCDir = "/Users/brian/macros/Simulations/Production";
   // dMCDir = "/cuore/user/zhubrian/MC/Bkg";
  dSaveDir = "/Users/brian/Dropbox/code/Scratch";
   // dSaveDir = "/cuore/user/zhubrian/";
  dDataIntegral = 0;

  // Bin size (keV) -- base binning is 1 keV
  dBinSize = 1; 
  // Histogram range - from 0 to 10 MeV
  dMinEnergy = 0.;
  dMaxEnergy = 10000.;
  dBinBase = fBinBase;
  dDataSet = fDataSet;

  dChiSquare = 0;

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
  fDataHistoTot     = new TH1D("fDataHistoTot",  "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM1      = new TH1D("fDataHistoM1",   "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM2      = new TH1D("fDataHistoM2",   "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM2Sum   = new TH1D("fDataHistoM2Sum",   "", dNBins, dMinEnergy, dMaxEnergy);

  // Data cuts 
  qtree = new TChain("qtree");
  base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
  base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
  base_cut = base_cut && "Channel != 1 && Channel != 10"; 

  // Old PSA cuts
  // base_cut = base_cut && "abs(BaselineSlope)<0.1";
  // base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";

  // New PSA cuts
  base_cut = base_cut && "NormBaselineSlope < 4.8 && NormBaselineSlope > -4";
  base_cut = base_cut && "NormRiseTime < 4.8 && NormRiseTime > -4";
  base_cut = base_cut && "NormDecayTime < 4.8 && NormDecayTime > -4";
  base_cut = base_cut && "NormDelay < 4.8 && NormDelay > -4";
  base_cut = base_cut && "NormTVL < 5.3 && NormTVL > -6";
  base_cut = base_cut && "NormTVR < 5.3 && NormTVR > -6";

  LoadData();


  // TF1 *fEff = new TF1("fEff", "[0]+[1]/(1+[2]*exp(-[3]*x)) + [4]/(1+[5]*exp(-[6]*x))", dMinEnergy, dMaxEnergy);
  // fEff->SetParameters(-4.71e-2, 1.12e-1, 2.29, -8.81e-5, 9.68e-1, 2.09, 1.58e-2);

  TF1 *fEff = new TF1("fEff", "[0]+[1]*TMath::Exp([2]*x)", dMinEnergy, dMaxEnergy);
  fEff->SetParameters(0.934, -4.982, -0.083);

  hEfficiency = new TH1D("hEfficiency", "", dNBins, dMinEnergy, dMaxEnergy);

  for(int i = 1; i <= hEfficiency->GetNbinsX(); i++)
  {
    hEfficiency->SetBinContent(i, fEff->Eval(hEfficiency->GetBinCenter(i)));
  }

  fDataHistoM1->Divide( hEfficiency );
  fDataHistoM2->Divide( hEfficiency );
  fDataHistoM2Sum->Divide( hEfficiency );


  // Total model histograms
  fModelTotM1      = new TH1D("fModelTotM1",      "Total Model",        dNBins, dMinEnergy, dMaxEnergy);  
  fModelTotM2      = new TH1D("fModelTotM2",      "Total Model",        dNBins, dMinEnergy, dMaxEnergy);  
  fModelTotM2Sum      = new TH1D("fModelTotM2Sum",      "Total Model",        dNBins, dMinEnergy, dMaxEnergy);  


/////////////////// Model histograms
//////////// Crystal M1 and M2
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
  hTeO2sb125M1     = new TH1D("hTeO2sb125M1",  "hTeO2sb125M1",  dNBins, dMinEnergy, dMaxEnergy);  

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
  hTeO2Sxu238M1_100      = new TH1D("hTeO2Sxu238M1_100",   "hTeO2Sxu238M1_100",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M1_100     = new TH1D("hTeO2Sxth232M1_100",  "hTeO2Sxth232M1_100",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M1_100     = new TH1D("hTeO2Sxpb210M1_100",  "hTeO2Sxpb210M1_100",  dNBins, dMinEnergy, dMaxEnergy);


  hTeO2th232onlyM1      = new TH1D("hTeO2th232onlyM1", "hTeO2th232onlyM1",     dNBins, dMinEnergy, dMaxEnergy);
  hTeO2ra228pb208M1     = new TH1D("hTeO2ra228pb208M1", "hTeO2ra228pb208M1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2th230onlyM1      = new TH1D("hTeO2th230onlyM1", "hTeO2th230onlyM1",     dNBins, dMinEnergy, dMaxEnergy);
  hTeO2u238th230M1      = new TH1D("hTeO2u238th230M1", "hTeO2u238th230M1",     dNBins, dMinEnergy, dMaxEnergy);
  hTeO2ra226pb210M1     = new TH1D("hTeO2ra226pb210M1", "hTeO2ra226pb210M1",   dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Sxth232onlyM1_001  = new TH1D("hTeO2Sxth232onlyM1_001", "hTeO2Sxth232onlyM1_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra228pb208M1_001 = new TH1D("hTeO2Sxra228pb208M1_001", "hTeO2Sxra228pb208M1_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238th230M1_001  = new TH1D("hTeO2Sxu238th230M1_001", "hTeO2Sxu238th230M1_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth230onlyM1_001  = new TH1D("hTeO2Sxth230onlyM1_001", "hTeO2Sxth230onlyM1_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra226pb210M1_001 = new TH1D("hTeO2Sxra226pb210M1_001", "hTeO2Sxra226pb210M1_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M1_0001     = new TH1D("hTeO2Sxpb210M1_0001", "hTeO2Sxpb210M1_0001",         dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Sxth232onlyM1_01  = new TH1D("hTeO2Sxth232onlyM1_01", "hTeO2Sxth232onlyM1_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra228pb208M1_01 = new TH1D("hTeO2Sxra228pb208M1_01", "hTeO2Sxra228pb208M1_01", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238th230M1_01  = new TH1D("hTeO2Sxu238th230M1_01", "hTeO2Sxu238th230M1_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth230onlyM1_01  = new TH1D("hTeO2Sxth230onlyM1_01", "hTeO2Sxth230onlyM1_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra226pb210M1_01 = new TH1D("hTeO2Sxra226pb210M1_01", "hTeO2Sxra226pb210M1_01", dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Sxth232onlyM1_0001  = new TH1D("hTeO2Sxth232onlyM1_0001", "hTeO2Sxth232onlyM1_0001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra228pb208M1_0001 = new TH1D("hTeO2Sxra228pb208M1_0001", "hTeO2Sxra228pb208M1_0001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238th230M1_0001  = new TH1D("hTeO2Sxu238th230M1_0001", "hTeO2Sxu238th230M1_0001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth230onlyM1_0001  = new TH1D("hTeO2Sxth230onlyM1_0001", "hTeO2Sxth230onlyM1_0001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra226pb210M1_0001 = new TH1D("hTeO2Sxra226pb210M1_0001", "hTeO2Sxra226pb210M1_0001", dNBins, dMinEnergy, dMaxEnergy);

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
  hTeO2sb125M2     = new TH1D("hTeO2sb125M2",  "hTeO2sb125M2",  dNBins, dMinEnergy, dMaxEnergy);

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
  hTeO2Sxu238M2_100      = new TH1D("hTeO2Sxu238M2_100",   "hTeO2Sxu238M2_100",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M2_100     = new TH1D("hTeO2Sxth232M2_100",  "hTeO2Sxth232M2_100",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2_100     = new TH1D("hTeO2Sxpb210M2_100",  "hTeO2Sxpb210M2_100",  dNBins, dMinEnergy, dMaxEnergy);

  hTeO2th232onlyM2      = new TH1D("hTeO2th232onlyM2", "hTeO2th232onlyM2",     dNBins, dMinEnergy, dMaxEnergy);
  hTeO2ra228pb208M2     = new TH1D("hTeO2ra228pb208M2", "hTeO2ra228pb208M2",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2u238th230M2      = new TH1D("hTeO2u238th230M2", "hTeO2u238th230M2",     dNBins, dMinEnergy, dMaxEnergy);
  hTeO2th230onlyM2      = new TH1D("hTeO2th230onlyM2", "hTeO2th230onlyM2",     dNBins, dMinEnergy, dMaxEnergy);
  hTeO2ra226pb210M2     = new TH1D("hTeO2ra226pb210M2", "hTeO2ra226pb210M2",   dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Sxth232onlyM2_001  = new TH1D("hTeO2Sxth232onlyM2_001", "hTeO2Sxth232onlyM2_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra228pb208M2_001 = new TH1D("hTeO2Sxra228pb208M2_001", "hTeO2Sxra228pb208M2_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238th230M2_001  = new TH1D("hTeO2Sxu238th230M2_001", "hTeO2Sxu238th230M2_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth230onlyM2_001  = new TH1D("hTeO2Sxth230onlyM2_001", "hTeO2Sxth230onlyM2_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra226pb210M2_001 = new TH1D("hTeO2Sxra226pb210M2_001", "hTeO2Sxra226pb210M2_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2_0001     = new TH1D("hTeO2Sxpb210M2_0001", "hTeO2Sxpb210M2_0001",         dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Sxth232onlyM2_01  = new TH1D("hTeO2Sxth232onlyM2_01", "hTeO2Sxth232onlyM2_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra228pb208M2_01 = new TH1D("hTeO2Sxra228pb208M2_01", "hTeO2Sxra228pb208M2_01", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238th230M2_01  = new TH1D("hTeO2Sxu238th230M2_01", "hTeO2Sxu238th230M2_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth230onlyM2_01  = new TH1D("hTeO2Sxth230onlyM2_01", "hTeO2Sxth230onlyM2_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra226pb210M2_01 = new TH1D("hTeO2Sxra226pb210M2_01", "hTeO2Sxra226pb210M2_01", dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Sxth232onlyM2_0001  = new TH1D("hTeO2Sxth232onlyM2_0001", "hTeO2Sxth232onlyM2_0001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra228pb208M2_0001 = new TH1D("hTeO2Sxra228pb208M2_0001", "hTeO2Sxra228pb208M2_0001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238th230M2_0001  = new TH1D("hTeO2Sxu238th230M2_0001", "hTeO2Sxu238th230M2_0001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth230onlyM2_0001  = new TH1D("hTeO2Sxth230onlyM2_0001", "hTeO2Sxth230onlyM2_0001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra226pb210M2_0001 = new TH1D("hTeO2Sxra226pb210M2_0001", "hTeO2Sxra226pb210M2_0001", dNBins, dMinEnergy, dMaxEnergy);


  hTeO20nuM2Sum       = new TH1D("hTeO20nuM2Sum",    "hTeO20nuM2Sum",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO22nuM2Sum       = new TH1D("hTeO22nuM2Sum",    "hTeO22nuM2Sum",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO2co60M2Sum      = new TH1D("hTeO2co60M2Sum",   "hTeO2co60M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2k40M2Sum       = new TH1D("hTeO2k40M2Sum",    "hTeO2k40M2Sum",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO2pb210M2Sum     = new TH1D("hTeO2pb210M2Sum",  "hTeO2pb210M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2po210M2Sum     = new TH1D("hTeO2po210M2Sum",  "hTeO2po210M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2te125M2Sum     = new TH1D("hTeO2te125M2Sum",  "hTeO2te125M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2th232M2Sum     = new TH1D("hTeO2th232M2Sum",  "hTeO2th232M2Sum",  dNBins, dMinEnergy, dMaxEnergy);  
  hTeO2th228M2Sum     = new TH1D("hTeO2th228M2Sum",  "hTeO2th228M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2ra226M2Sum     = new TH1D("hTeO2ra226M2Sum",  "hTeO2ra226M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2rn222M2Sum     = new TH1D("hTeO2rn222M2Sum",  "hTeO2rn222M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2u238M2Sum      = new TH1D("hTeO2u238M2Sum",   "hTeO2u238M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2th230M2Sum     = new TH1D("hTeO2th230M2Sum",  "hTeO2th230M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2u234M2Sum      = new TH1D("hTeO2u234M2Sum",   "hTeO2u234M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2sb125M2Sum     = new TH1D("hTeO2sb125M2Sum",  "hTeO2sb125M2Sum",  dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Spb210M2Sum_01      = new TH1D("hTeO2Spb210M2Sum_01",   "hTeO2Spb210M2Sum_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Spo210M2Sum_001     = new TH1D("hTeO2Spo210M2Sum_001",  "hTeO2Spo210M2Sum_001",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Spo210M2Sum_01      = new TH1D("hTeO2Spo210M2Sum_01",   "hTeO2Spo210M2Sum_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sth232M2Sum_01      = new TH1D("hTeO2Sth232M2Sum_01",   "hTeO2Sth232M2Sum_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Su238M2Sum_01       = new TH1D("hTeO2Su238M2Sum_01",    "hTeO2Su238M2Sum_01",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2Sum_001    = new TH1D("hTeO2Sxpb210M2Sum_001", "hTeO2Sxpb210M2Sum_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2Sum_01     = new TH1D("hTeO2Sxpb210M2Sum_01",  "hTeO2Sxpb210M2Sum_01",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2Sum_1      = new TH1D("hTeO2Sxpb210M2Sum_1",   "hTeO2Sxpb210M2Sum_1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2Sum_10     = new TH1D("hTeO2Sxpb210M2Sum_10",  "hTeO2Sxpb210M2Sum_10",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpo210M2Sum_001    = new TH1D("hTeO2Sxpo210M2Sum_001", "hTeO2Sxpo210M2Sum_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpo210M2Sum_01     = new TH1D("hTeO2Sxpo210M2Sum_01",  "hTeO2Sxpo210M2Sum_01",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpo210M2Sum_1      = new TH1D("hTeO2Sxpo210M2Sum_1",   "hTeO2Sxpo210M2Sum_1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M2Sum_001    = new TH1D("hTeO2Sxth232M2Sum_001", "hTeO2Sxth232M2Sum_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M2Sum_01     = new TH1D("hTeO2Sxth232M2Sum_01",  "hTeO2Sxth232M2Sum_01",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M2Sum_1      = new TH1D("hTeO2Sxth232M2Sum_1",   "hTeO2Sxth232M2Sum_1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M2Sum_10     = new TH1D("hTeO2Sxth232M2Sum_10",  "hTeO2Sxth232M2Sum_10",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238M2Sum_001     = new TH1D("hTeO2Sxu238M2Sum_001",  "hTeO2Sxu238M2Sum_001",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238M2Sum_01      = new TH1D("hTeO2Sxu238M2Sum_01",   "hTeO2Sxu238M2Sum_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238M2Sum_1       = new TH1D("hTeO2Sxu238M2Sum_1",    "hTeO2Sxu238M2Sum_1",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238M2Sum_10      = new TH1D("hTeO2Sxu238M2Sum_10",   "hTeO2Sxu238M2Sum_10",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238M2Sum_100      = new TH1D("hTeO2Sxu238M2Sum_100",   "hTeO2Sxu238M2Sum_100",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M2Sum_100     = new TH1D("hTeO2Sxth232M2Sum_100",  "hTeO2Sxth232M2Sum_100",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2Sum_100     = new TH1D("hTeO2Sxpb210M2Sum_100",  "hTeO2Sxpb210M2Sum_100",  dNBins, dMinEnergy, dMaxEnergy);

  hTeO2th232onlyM2Sum      = new TH1D("hTeO2th232onlyM2Sum", "hTeO2th232onlyM2Sum",     dNBins, dMinEnergy, dMaxEnergy);
  hTeO2ra228pb208M2Sum     = new TH1D("hTeO2ra228pb208M2Sum", "hTeO2ra228pb208M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2u238th230M2Sum      = new TH1D("hTeO2u238th230M2Sum", "hTeO2u238th230M2Sum",     dNBins, dMinEnergy, dMaxEnergy);
  hTeO2th230onlyM2Sum      = new TH1D("hTeO2th230onlyM2Sum", "hTeO2th230onlyM2Sum",     dNBins, dMinEnergy, dMaxEnergy);
  hTeO2ra226pb210M2Sum     = new TH1D("hTeO2ra226pb210M2Sum", "hTeO2ra226pb210M2Sum",   dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Sxth232onlyM2Sum_001  = new TH1D("hTeO2Sxth232onlyM2Sum_001", "hTeO2Sxth232onlyM2Sum_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra228pb208M2Sum_001 = new TH1D("hTeO2Sxra228pb208M2Sum_001", "hTeO2Sxra228pb208M2Sum_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238th230M2Sum_001  = new TH1D("hTeO2Sxu238th230M2Sum_001", "hTeO2Sxu238th230M2Sum_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth230onlyM2Sum_001  = new TH1D("hTeO2Sxth230onlyM2Sum_001", "hTeO2Sxth230onlyM2Sum_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra226pb210M2Sum_001 = new TH1D("hTeO2Sxra226pb210M2Sum_001", "hTeO2Sxra226pb210M2Sum_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2Sum_0001     = new TH1D("hTeO2Sxpb210M2Sum_0001", "hTeO2Sxpb210M2Sum_0001",         dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Sxth232onlyM2Sum_01  = new TH1D("hTeO2Sxth232onlyM2Sum_01", "hTeO2Sxth232onlyM2Sum_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra228pb208M2Sum_01 = new TH1D("hTeO2Sxra228pb208M2Sum_01", "hTeO2Sxra228pb208M2Sum_01", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238th230M2Sum_01  = new TH1D("hTeO2Sxu238th230M2Sum_01", "hTeO2Sxu238th230M2Sum_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth230onlyM2Sum_01  = new TH1D("hTeO2Sxth230onlyM2Sum_01", "hTeO2Sxth230onlyM2Sum_01",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra226pb210M2Sum_01 = new TH1D("hTeO2Sxra226pb210M2Sum_01", "hTeO2Sxra226pb210M2Sum_01", dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Sxth232onlyM2Sum_0001  = new TH1D("hTeO2Sxth232onlyM2Sum_0001", "hTeO2Sxth232onlyM2Sum_0001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra228pb208M2Sum_0001 = new TH1D("hTeO2Sxra228pb208M2Sum_0001", "hTeO2Sxra228pb208M2Sum_0001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238th230M2Sum_0001  = new TH1D("hTeO2Sxu238th230M2Sum_0001", "hTeO2Sxu238th230M2Sum_0001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth230onlyM2Sum_0001  = new TH1D("hTeO2Sxth230onlyM2Sum_0001", "hTeO2Sxth230onlyM2Sum_0001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra226pb210M2Sum_0001 = new TH1D("hTeO2Sxra226pb210M2Sum_0001", "hTeO2Sxra226pb210M2Sum_0001", dNBins, dMinEnergy, dMaxEnergy);

/////////// CuBox + CuFrame M1 and M2
  hCuBox_CuFrameco60M1      = new TH1D("hCuBox_CuFrameco60M1", "hCuBox_CuFrameco60M1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramek40M1       = new TH1D("hCuBox_CuFramek40M1", "hCuBox_CuFramek40M1",      dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M1     = new TH1D("hCuBox_CuFrameth232M1", "hCuBox_CuFrameth232M1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M1      = new TH1D("hCuBox_CuFrameu238M1", "hCuBox_CuFrameu238M1",    dNBins, dMinEnergy, dMaxEnergy);

  hCuBox_CuFrameth232M1_10  = new TH1D("hCuBox_CuFrameth232M1_10", "hCuBox_CuFrameth232M1_10",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M1_10   = new TH1D("hCuBox_CuFrameu238M1_10", "hCuBox_CuFrameu238M1_10",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M1_10  = new TH1D("hCuBox_CuFramepb210M1_10", "hCuBox_CuFramepb210M1_10",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M1_1   = new TH1D("hCuBox_CuFramepb210M1_1", "hCuBox_CuFramepb210M1_1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M1_01  = new TH1D("hCuBox_CuFramepb210M1_01", "hCuBox_CuFramepb210M1_01",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M1_001 = new TH1D("hCuBox_CuFramepb210M1_001", "hCuBox_CuFramepb210M1_001",  dNBins, dMinEnergy, dMaxEnergy);

  hCuBox_CuFrameth232M1_1  = new TH1D("hCuBox_CuFrameth232M1_1", "hCuBox_CuFrameth232M1_1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M1_1   = new TH1D("hCuBox_CuFrameu238M1_1", "hCuBox_CuFrameu238M1_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M1_01  = new TH1D("hCuBox_CuFrameth232M1_01", "hCuBox_CuFrameth232M1_01",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M1_01   = new TH1D("hCuBox_CuFrameu238M1_01", "hCuBox_CuFrameu238M1_01",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M1_001  = new TH1D("hCuBox_CuFrameth232M1_001", "hCuBox_CuFrameth232M1_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M1_001   = new TH1D("hCuBox_CuFrameu238M1_001", "hCuBox_CuFrameu238M1_001",    dNBins, dMinEnergy, dMaxEnergy);

  hCuBox_CuFramepb210M1_100  = new TH1D("hCuBox_CuFramepb210M1_100", "hCuBox_CuFramepb210M1_100",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M1_100  = new TH1D("hCuBox_CuFrameth232M1_100", "hCuBox_CuFrameth232M1_100",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M1_100   = new TH1D("hCuBox_CuFrameu238M1_100", "hCuBox_CuFrameu238M1_100",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M1_50  = new TH1D("hCuBox_CuFramepb210M1_50", "hCuBox_CuFramepb210M1_50",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M1_50  = new TH1D("hCuBox_CuFrameth232M1_50", "hCuBox_CuFrameth232M1_50",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M1_50   = new TH1D("hCuBox_CuFrameu238M1_50", "hCuBox_CuFrameu238M1_50",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M1_5  = new TH1D("hCuBox_CuFramepb210M1_5", "hCuBox_CuFramepb210M1_5",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M1_5  = new TH1D("hCuBox_CuFrameth232M1_5", "hCuBox_CuFrameth232M1_5",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M1_5   = new TH1D("hCuBox_CuFrameu238M1_5", "hCuBox_CuFrameu238M1_5",    dNBins, dMinEnergy, dMaxEnergy);

  hCuBox_CuFrameco60M2      = new TH1D("hCuBox_CuFrameco60M2", "hCuBox_CuFrameco60M2",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramek40M2       = new TH1D("hCuBox_CuFramek40M2", "hCuBox_CuFramek40M2",      dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M2     = new TH1D("hCuBox_CuFrameth232M2", "hCuBox_CuFrameth232M2",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2      = new TH1D("hCuBox_CuFrameu238M2", "hCuBox_CuFrameu238M2",    dNBins, dMinEnergy, dMaxEnergy);

  hCuBox_CuFrameth232M2_10  = new TH1D("hCuBox_CuFrameth232M2_10", "hCuBox_CuFrameth232M2_10",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2_10   = new TH1D("hCuBox_CuFrameu238M2_10", "hCuBox_CuFrameu238M2_10",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M2_10  = new TH1D("hCuBox_CuFramepb210M2_10", "hCuBox_CuFramepb210M2_10",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M2_1   = new TH1D("hCuBox_CuFramepb210M2_1", "hCuBox_CuFramepb210M2_1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M2_01  = new TH1D("hCuBox_CuFramepb210M2_01", "hCuBox_CuFramepb210M2_01",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M2_001 = new TH1D("hCuBox_CuFramepb210M2_001", "hCuBox_CuFramepb210M2_001",  dNBins, dMinEnergy, dMaxEnergy);

  hCuBox_CuFrameth232M2_1  = new TH1D("hCuBox_CuFrameth232M2_1", "hCuBox_CuFrameth232M2_1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2_1   = new TH1D("hCuBox_CuFrameu238M2_1", "hCuBox_CuFrameu238M2_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M2_01  = new TH1D("hCuBox_CuFrameth232M2_01", "hCuBox_CuFrameth232M2_01",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2_01   = new TH1D("hCuBox_CuFrameu238M2_01", "hCuBox_CuFrameu238M2_01",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M2_001  = new TH1D("hCuBox_CuFrameth232M2_001", "hCuBox_CuFrameth232M2_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2_001   = new TH1D("hCuBox_CuFrameu238M2_001", "hCuBox_CuFrameu238M2_001",    dNBins, dMinEnergy, dMaxEnergy);

  hCuBox_CuFramepb210M2_100  = new TH1D("hCuBox_CuFramepb210M2_100", "hCuBox_CuFramepb210M2_100",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M2_100  = new TH1D("hCuBox_CuFrameth232M2_100", "hCuBox_CuFrameth232M2_100",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2_100   = new TH1D("hCuBox_CuFrameu238M2_100", "hCuBox_CuFrameu238M2_100",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M2_50  = new TH1D("hCuBox_CuFramepb210M2_50", "hCuBox_CuFramepb210M2_50",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M2_50  = new TH1D("hCuBox_CuFrameth232M2_50", "hCuBox_CuFrameth232M2_50",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2_50   = new TH1D("hCuBox_CuFrameu238M2_50", "hCuBox_CuFrameu238M2_50",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M2_5  = new TH1D("hCuBox_CuFramepb210M2_5", "hCuBox_CuFramepb210M2_5",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M2_5  = new TH1D("hCuBox_CuFrameth232M2_5", "hCuBox_CuFrameth232M2_5",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2_5   = new TH1D("hCuBox_CuFrameu238M2_5", "hCuBox_CuFrameu238M2_5",    dNBins, dMinEnergy, dMaxEnergy);

  hCuBox_CuFrameco60M2Sum   = new TH1D("hCuBox_CuFrameco60M2Sum", "hCuBox_CuFrameco60M2Sum",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramek40M2Sum    = new TH1D("hCuBox_CuFramek40M2Sum", "hCuBox_CuFramek40M2Sum",      dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M2Sum  = new TH1D("hCuBox_CuFrameth232M2Sum", "hCuBox_CuFrameth232M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2Sum   = new TH1D("hCuBox_CuFrameu238M2Sum", "hCuBox_CuFrameu238M2Sum",    dNBins, dMinEnergy, dMaxEnergy);

  hCuBox_CuFrameth232M2Sum_10   = new TH1D("hCuBox_CuFrameth232M2Sum_10", "hCuBox_CuFrameth232M2Sum_10",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2Sum_10    = new TH1D("hCuBox_CuFrameu238M2Sum_10", "hCuBox_CuFrameu238M2Sum_10",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M2Sum_10   = new TH1D("hCuBox_CuFramepb210M2Sum_10", "hCuBox_CuFramepb210M2Sum_10",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M2Sum_1   = new TH1D("hCuBox_CuFramepb210M2Sum_1", "hCuBox_CuFramepb210M2Sum_1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M2Sum_01  = new TH1D("hCuBox_CuFramepb210M2Sum_01", "hCuBox_CuFramepb210M2Sum_01",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M2Sum_001 = new TH1D("hCuBox_CuFramepb210M2Sum_001", "hCuBox_CuFramepb210M2Sum_001",  dNBins, dMinEnergy, dMaxEnergy);

  hCuBox_CuFrameth232M2Sum_1  = new TH1D("hCuBox_CuFrameth232M2Sum_1", "hCuBox_CuFrameth232M2Sum_1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2Sum_1   = new TH1D("hCuBox_CuFrameu238M2Sum_1", "hCuBox_CuFrameu238M2Sum_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M2Sum_01  = new TH1D("hCuBox_CuFrameth232M2Sum_01", "hCuBox_CuFrameth232M2Sum_01",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2Sum_01   = new TH1D("hCuBox_CuFrameu238M2Sum_01", "hCuBox_CuFrameu238M2Sum_01",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M2Sum_001  = new TH1D("hCuBox_CuFrameth232M2Sum_001", "hCuBox_CuFrameth232M2Sum_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2Sum_001   = new TH1D("hCuBox_CuFrameu238M2Sum_001", "hCuBox_CuFrameu238M2Sum_001",    dNBins, dMinEnergy, dMaxEnergy);

  hCuBox_CuFramepb210M2Sum_100  = new TH1D("hCuBox_CuFramepb210M2Sum_100", "hCuBox_CuFramepb210M2Sum_100",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M2Sum_100  = new TH1D("hCuBox_CuFrameth232M2Sum_100", "hCuBox_CuFrameth232M2Sum_100",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2Sum_100   = new TH1D("hCuBox_CuFrameu238M2Sum_100", "hCuBox_CuFrameu238M2Sum_100",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M2Sum_50  = new TH1D("hCuBox_CuFramepb210M2Sum_50", "hCuBox_CuFramepb210M2Sum_50",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M2Sum_50  = new TH1D("hCuBox_CuFrameth232M2Sum_50", "hCuBox_CuFrameth232M2Sum_50",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2Sum_50   = new TH1D("hCuBox_CuFrameu238M2Sum_50", "hCuBox_CuFrameu238M2Sum_50",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFramepb210M2Sum_5  = new TH1D("hCuBox_CuFramepb210M2Sum_5", "hCuBox_CuFramepb210M2Sum_5",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameth232M2Sum_5  = new TH1D("hCuBox_CuFrameth232M2Sum_5", "hCuBox_CuFrameth232M2Sum_5",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_CuFrameu238M2Sum_5   = new TH1D("hCuBox_CuFrameu238M2Sum_5", "hCuBox_CuFrameu238M2Sum_5",    dNBins, dMinEnergy, dMaxEnergy);

/////////// 50mK M1 and M2
  h50mKcs137M1     = new TH1D("h50mKcs137M1",  "h50mKcs137M1",  dNBins, dMinEnergy, dMaxEnergy);
  h50mKcs137M2     = new TH1D("h50mKcs137M2",  "h50mKcs137M2",  dNBins, dMinEnergy, dMaxEnergy);
  h50mKcs137M2Sum     = new TH1D("h50mKcs137M2Sum",  "h50mKcs137M2Sum",  dNBins, dMinEnergy, dMaxEnergy);

/////////// Internal Shields M1 and M2
  hInternalco60M1   = new TH1D("hInternalco60M1", "hInternalco60M1",    dNBins, dMinEnergy, dMaxEnergy);
  hInternalk40M1    = new TH1D("hInternalk40M1", "hInternalk40M1",      dNBins, dMinEnergy, dMaxEnergy);
  hInternalth232M1  = new TH1D("hInternalth232M1", "hInternalth232M1",  dNBins, dMinEnergy, dMaxEnergy);
  hInternalu238M1   = new TH1D("hInternalu238M1", "hInternalu238M1",    dNBins, dMinEnergy, dMaxEnergy);

  hInternalco60M2   = new TH1D("hInternalco60M2", "hInternalco60M2",    dNBins, dMinEnergy, dMaxEnergy);
  hInternalk40M2    = new TH1D("hInternalk40M2", "hInternalk40M2",      dNBins, dMinEnergy, dMaxEnergy);
  hInternalth232M2  = new TH1D("hInternalth232M2", "hInternalth232M2",  dNBins, dMinEnergy, dMaxEnergy);
  hInternalu238M2   = new TH1D("hInternalu238M2", "hInternalu238M2",    dNBins, dMinEnergy, dMaxEnergy);

  hInternalco60M2Sum  = new TH1D("hInternalco60M2Sum", "hInternalco60M2Sum",    dNBins, dMinEnergy, dMaxEnergy);
  hInternalk40M2Sum   = new TH1D("hInternalk40M2Sum", "hInternalk40M2Sum",      dNBins, dMinEnergy, dMaxEnergy);
  hInternalth232M2Sum = new TH1D("hInternalth232M2Sum", "hInternalth232M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  hInternalu238M2Sum  = new TH1D("hInternalu238M2Sum", "hInternalu238M2Sum",    dNBins, dMinEnergy, dMaxEnergy);


////////// Roman Lead M1 and M2
  hPbRomco60M1      = new TH1D("hPbRomco60M1",   "hPbRomco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hPbRomcs137M1     = new TH1D("hPbRomcs137M1",  "hPbRomcs137M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomk40M1       = new TH1D("hPbRomk40M1",    "hPbRomk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hPbRompb210M1     = new TH1D("hPbRompb210M1",  "hPbRompb210M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomth232M1     = new TH1D("hPbRomth232M1",  "hPbRomth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomu238M1      = new TH1D("hPbRomu238M1",   "hPbRomu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hPbRomco60M2      = new TH1D("hPbRomco60M2",   "hPbRomco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hPbRomcs137M2     = new TH1D("hPbRomcs137M2",  "hPbRomcs137M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomk40M2       = new TH1D("hPbRomk40M2",    "hPbRomk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hPbRompb210M2     = new TH1D("hPbRompb210M2",  "hPbRompb210M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomth232M2     = new TH1D("hPbRomth232M2",  "hPbRomth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomu238M2      = new TH1D("hPbRomu238M2",   "hPbRomu238M2",   dNBins, dMinEnergy, dMaxEnergy);

  hPbRomco60M2Sum      = new TH1D("hPbRomco60M2Sum",   "hPbRomco60M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hPbRomcs137M2Sum     = new TH1D("hPbRomcs137M2Sum",  "hPbRomcs137M2Sum",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomk40M2Sum       = new TH1D("hPbRomk40M2Sum",    "hPbRomk40M2Sum",    dNBins, dMinEnergy, dMaxEnergy);
  hPbRompb210M2Sum     = new TH1D("hPbRompb210M2Sum",  "hPbRompb210M2Sum",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomth232M2Sum     = new TH1D("hPbRomth232M2Sum",  "hPbRomth232M2Sum",  dNBins, dMinEnergy, dMaxEnergy);  
  hPbRomu238M2Sum      = new TH1D("hPbRomu238M2Sum",   "hPbRomu238M2Sum",   dNBins, dMinEnergy, dMaxEnergy);

///////// OVC M1 and M2
  hOVCco60M1      = new TH1D("hOVCco60M1",   "hOVCco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hOVCk40M1       = new TH1D("hOVCk40M1",    "hOVCk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hOVCth232M1     = new TH1D("hOVCth232M1",  "hOVCth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hOVCu238M1      = new TH1D("hOVCu238M1",   "hOVCu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hOVCco60M2      = new TH1D("hOVCco60M2",   "hOVCco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hOVCk40M2       = new TH1D("hOVCk40M2",    "hOVCk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hOVCth232M2     = new TH1D("hOVCth232M2",  "hOVCth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hOVCu238M2      = new TH1D("hOVCu238M2",   "hOVCu238M2",   dNBins, dMinEnergy, dMaxEnergy);  

  hOVCco60M2Sum      = new TH1D("hOVCco60M2Sum",   "hOVCco60M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hOVCk40M2Sum       = new TH1D("hOVCk40M2Sum",    "hOVCk40M2Sum",    dNBins, dMinEnergy, dMaxEnergy);
  hOVCth232M2Sum     = new TH1D("hOVCth232M2Sum",  "hOVCth232M2Sum",  dNBins, dMinEnergy, dMaxEnergy);  
  hOVCu238M2Sum      = new TH1D("hOVCu238M2Sum",   "hOVCu238M2Sum",   dNBins, dMinEnergy, dMaxEnergy);  


////////// External Shields
  hExtPbbi210M1 = new TH1D("hExtPbbi210M1", "hExtPbbi210M1", dNBins, dMinEnergy, dMaxEnergy);

  hExtPbbi210M2 = new TH1D("hExtPbbi210M2", "hExtPbbi210M2", dNBins, dMinEnergy, dMaxEnergy);
  
  hExtPbbi210M2Sum = new TH1D("hExtPbbi210M2Sum", "hExtPbbi210M2Sum", dNBins, dMinEnergy, dMaxEnergy);


/////////// Fudge Factors
  hOVC804M1 = new TH1D("hOVC804M1", "hOVC804M1", dNBins, dMinEnergy, dMaxEnergy);
  hOVC835M1 = new TH1D("hOVC835M1", "hOVC835M1", dNBins, dMinEnergy, dMaxEnergy);
  hOVC1063M1 = new TH1D("hOVC1063M1", "hOVC1063M1", dNBins, dMinEnergy, dMaxEnergy);

  hOVC804M2 = new TH1D("hOVC804M2", "hOVC804M2", dNBins, dMinEnergy, dMaxEnergy);
  hOVC835M2 = new TH1D("hOVC835M2", "hOVC835M2", dNBins, dMinEnergy, dMaxEnergy);
  hOVC1063M2 = new TH1D("hOVC1063M2", "hOVC1063M2", dNBins, dMinEnergy, dMaxEnergy);

  hOVC804M2Sum = new TH1D("hOVC804M2Sum", "hOVC804M2Sum", dNBins, dMinEnergy, dMaxEnergy);
  hOVC835M2Sum = new TH1D("hOVC835M2Sum", "hOVC835M2Sum", dNBins, dMinEnergy, dMaxEnergy);
  hOVC1063M2Sum = new TH1D("hOVC1063M2Sum", "hOVC1063M2Sum", dNBins, dMinEnergy, dMaxEnergy);

// Mess with rebinning here 
  // fDataHistoM1->Rebin(1);
  // fDataHistoM2->Rebin(1);
  // fDataHistoM2Sum->Rebin(1);
  dBaseBinSize = dBinSize*1;

/////// Adaptive binning
 // Calculates adaptive binning vectors
  dAdaptiveVectorM1 = AdaptiveBinning(fDataHistoM1, dBinBase);
  dAdaptiveBinsM1 = dAdaptiveVectorM1.size() - 1;
  dAdaptiveArrayM1 = &dAdaptiveVectorM1[0];
  dAdaptiveVectorM2 = AdaptiveBinning(fDataHistoM2, dBinBase);
  dAdaptiveBinsM2 = dAdaptiveVectorM2.size() - 1;
  dAdaptiveArrayM2 = &dAdaptiveVectorM2[0];
  dAdaptiveVectorM2Sum = AdaptiveBinning(fDataHistoM2Sum, dBinBase);
  dAdaptiveBinsM2Sum = dAdaptiveVectorM2Sum.size() - 1;
  dAdaptiveArrayM2Sum = &dAdaptiveVectorM2Sum[0];

 
  // Adaptive binning data
  fAdapDataHistoM1   = new TH1D("fAdapDataHistoM1",   "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapDataHistoM2   = new TH1D("fAdapDataHistoM2",   "", dAdaptiveBinsM2, dAdaptiveArrayM2);
  fAdapDataHistoM2Sum   = new TH1D("fAdapDataHistoM2Sum",   "", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  
  hnewM1 = fDataHistoM1->Rebin(dAdaptiveBinsM1, "hnewM1", dAdaptiveArrayM1);
  hnewM2 = fDataHistoM2->Rebin(dAdaptiveBinsM2, "hnewM2", dAdaptiveArrayM2);
  hnewM2Sum = fDataHistoM2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewM2Sum", dAdaptiveArrayM2Sum);

  for(int i = 1; i <= dAdaptiveBinsM1; i++)
  {
    fAdapDataHistoM1->SetBinContent(i, hnewM1->GetBinContent(i)/hnewM1->GetBinWidth(i));
    fAdapDataHistoM1->SetBinError(i, TMath::Sqrt(hnewM1->GetBinContent(i))/hnewM1->GetBinWidth(i));
    // fAdapDataHistoM1->SetBinError(i, 0); // If I don't want errors for some reason
  }

  for(int i = 1; i <= dAdaptiveBinsM2; i++)
  {
    fAdapDataHistoM2->SetBinContent(i, hnewM2->GetBinContent(i)/hnewM2->GetBinWidth(i));
    fAdapDataHistoM2->SetBinError(i, TMath::Sqrt(hnewM2->GetBinContent(i))/hnewM2->GetBinWidth(i));
    // fAdapDataHistoM2->SetBinError(i, 0); // If I don't want errors for some reason
  }

  for(int i = 1; i <= dAdaptiveBinsM2Sum; i++)
  {
    fAdapDataHistoM2Sum->SetBinContent(i, hnewM2Sum->GetBinContent(i)/hnewM2Sum->GetBinWidth(i));
    fAdapDataHistoM2Sum->SetBinError(i, TMath::Sqrt(hnewM2Sum->GetBinContent(i))/hnewM2Sum->GetBinWidth(i));
    // fAdapDataHistoM2Sum->SetBinError(i, 0); // If I don't want errors for some reason
  }

  // dDataIntegral = fDataHistoM1->Integral(1, dNBins);
  // dDataIntegralM1 = fDataHistoM1->Integral(1, 10000/dBinSize);
  // dDataIntegralM2 = fDataHistoM2->Integral(1, 10000/dBinSize);
  // dDataIntegralM2Sum = fDataHistoM2Sum->Integral(1, 10000/dBinSize);
  dDataIntegralTot = qtree->GetEntries();

  dDataIntegralM1 = fAdapDataHistoM1->Integral("width");
  dDataIntegralM2 = fAdapDataHistoM2->Integral("width");
  dDataIntegralM2Sum = fAdapDataHistoM2Sum->Integral("width");

  // dDataIntegralM1 = dDataIntegralM1;
  // dDataIntegralM2 = dDataIntegralM2;

  dFitMinBinM1 = fAdapDataHistoM1->FindBin(dFitMin);
  dFitMinBinM2 = fAdapDataHistoM2->FindBin(dFitMin);
  dFitMinBinM2Sum = fAdapDataHistoM2Sum->FindBin(dFitMin);
  dFitMaxBinM1 = fAdapDataHistoM1->FindBin(dFitMax);
  dFitMaxBinM2 = fAdapDataHistoM2->FindBin(dFitMax);
  dFitMaxBinM2Sum = fAdapDataHistoM2Sum->FindBin(dFitMax);

  // Outputs on screen
  cout << "Fit M1 from bin: " << dFitMinBinM1 << " to " << dFitMaxBinM1 << endl;
  cout << "Fit M2 from bin: " << dFitMinBinM2 << " to " << dFitMaxBinM2 << endl;
  cout << "Fit M2Sum from bin: " << dFitMinBinM2Sum << " to " << dFitMaxBinM2Sum << endl;

  cout << "Total Events in background spectrum: " << dDataIntegralTot << endl; 
  cout << "Events in background spectrum (M1): " << dDataIntegralM1 << endl;
  cout << "Events in background spectrum (M2): " << dDataIntegralM2 << endl;
  cout << "Events in background spectrum (M2Sum): " << dDataIntegralM2Sum << endl;
  cout << "Livetime of background: " << dLivetimeYr << endl;


//////////////// Adaptive binned histograms
////////// Total Adaptive binning histograms
  fModelTotAdapM1      = new TH1D("fModelTotAdapM1",      "Total PDF M1", dAdaptiveBinsM1, dAdaptiveArrayM1);  
  fModelTotAdapM2      = new TH1D("fModelTotAdapM2",      "Total PDF M2", dAdaptiveBinsM2, dAdaptiveArrayM2);  
  fModelTotAdapM2Sum      = new TH1D("fModelTotAdapM2Sum",      "Total PDF M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  

/////////// Crystal M1 and M2
  hAdapTeO20nuM1       = new TH1D("hAdapTeO20nuM1",    "TeO2 Bulk 0nu M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO22nuM1       = new TH1D("hAdapTeO22nuM1",    "TeO2 Bulk 2nu M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2co60M1      = new TH1D("hAdapTeO2co60M1",   "TeO2 Bulk co60 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2k40M1       = new TH1D("hAdapTeO2k40M1",    "TeO2 Bulk k40 M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2pb210M1     = new TH1D("hAdapTeO2pb210M1",  "TeO2 Bulk pb210 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2po210M1     = new TH1D("hAdapTeO2po210M1",  "TeO2 Bulk po210 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2te125M1     = new TH1D("hAdapTeO2te125M1",  "TeO2 Bulk te125 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2th232M1     = new TH1D("hAdapTeO2th232M1",  "TeO2 Bulk th232 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapTeO2th228M1     = new TH1D("hAdapTeO2th228M1",  "TeO2 Bulk th228 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2ra226M1     = new TH1D("hAdapTeO2ra226M1",  "TeO2 Bulk ra226 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2rn222M1     = new TH1D("hAdapTeO2rn222M1",  "TeO2 Bulk rn222 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2u238M1      = new TH1D("hAdapTeO2u238M1",   "TeO2 Bulk u238 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2th230M1     = new TH1D("hAdapTeO2th230M1",  "TeO2 Bulk th230M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2u234M1      = new TH1D("hAdapTeO2u234M1",   "TeO2 Bulk u234 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2sb125M1     = new TH1D("hAdapTeO2sb125M1",  "TeO2 Bulk sb125 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapTeO2Spb210M1_01      = new TH1D("hAdapTeO2Spb210M1_01",   "TeO2 S pb210M1 0.1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Spo210M1_001     = new TH1D("hAdapTeO2Spo210M1_001",  "TeO2 S po210M1 0.01 #mum",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Spo210M1_01      = new TH1D("hAdapTeO2Spo210M1_01",   "TeO2 S po210M1 0.1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sth232M1_01      = new TH1D("hAdapTeO2Sth232M1_01",   "TeO2 S th232M1 0.1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Su238M1_01       = new TH1D("hAdapTeO2Su238M1_01",    "TeO2 S u238M1 0.1 #mum",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpb210M1_001    = new TH1D("hAdapTeO2Sxpb210M1_001", "TeO2 Sx pb210 M1 0.01 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpb210M1_01     = new TH1D("hAdapTeO2Sxpb210M1_01",  "TeO2 Sx pb210 M1 0.1 #mum",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpb210M1_1      = new TH1D("hAdapTeO2Sxpb210M1_1",   "TeO2 Sx pb210 M1 1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpb210M1_10     = new TH1D("hAdapTeO2Sxpb210M1_10",  "TeO2 Sx pb210 M1 10 #mum",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpo210M1_001    = new TH1D("hAdapTeO2Sxpo210M1_001", "TeO2 Sx po210 M1 0.01 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpo210M1_01     = new TH1D("hAdapTeO2Sxpo210M1_01",  "TeO2 Sx po210 M1 0.1 #mum",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpo210M1_1      = new TH1D("hAdapTeO2Sxpo210M1_1",   "TeO2 Sx po210 M1 1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxth232M1_001    = new TH1D("hAdapTeO2Sxth232M1_001", "TeO2 Sx th232 M1 0.01 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxth232M1_01     = new TH1D("hAdapTeO2Sxth232M1_01",  "TeO2 Sx th232 M1 0.1 #mum",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxth232M1_1      = new TH1D("hAdapTeO2Sxth232M1_1",   "TeO2 Sx th232 M1 1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxth232M1_10     = new TH1D("hAdapTeO2Sxth232M1_10",  "TeO2 Sx th232 M1 10 #mum",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxu238M1_001     = new TH1D("hAdapTeO2Sxu238M1_001",  "TeO2 Sx u238 M1 0.01 #mum",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxu238M1_01      = new TH1D("hAdapTeO2Sxu238M1_01",   "TeO2 Sx u238 M1 0.1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxu238M1_1       = new TH1D("hAdapTeO2Sxu238M1_1",    "TeO2 Sx u238 M1 1 #mum",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxu238M1_10      = new TH1D("hAdapTeO2Sxu238M1_10",   "TeO2 Sx u238 M1 10 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapTeO2Sxu238M1_100      = new TH1D("hAdapTeO2Sxu238M1_100",   "TeO2 Sx u238 M1 100 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxth232M1_100     = new TH1D("hAdapTeO2Sxth232M1_100",  "TeO2 Sx th232 M1 100 #mum",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpb210M1_100     = new TH1D("hAdapTeO2Sxpb210M1_100",  "TeO2 Sx pb210 M1 100 #mum",  dAdaptiveBinsM1, dAdaptiveArrayM1);


  hAdapTeO2th232onlyM1      = new TH1D("hAdapTeO2th232onlyM1",   "TeO2 Bulk th232 only M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2ra228pb208M1     = new TH1D("hAdapTeO2ra228pb208M1",  "TeO2 Bulk ra228-pb208 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2th230onlyM1      = new TH1D("hAdapTeO2th230onlyM1", "TeO2 Bulk th232 only M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2u238th230M1      = new TH1D("hAdapTeO2u238th230M1", "TeO2 Bulk u238-th230 M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2ra226pb210M1     = new TH1D("hAdapTeO2ra226pb210M1", "TeO2 Bulk ra226-pb210 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapTeO2Sxth232onlyM1_001  = new TH1D("hAdapTeO2Sxth232onlyM1_001", "TeO2 Sx th232 only M1 0.01 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxra228pb208M1_001 = new TH1D("hAdapTeO2Sxra228pb208M1_001", "TeO2 Sx ra228-pb208 M1 0.01 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxu238th230M1_001  = new TH1D("hAdapTeO2Sxu238th230M1_001", "TeO2 Sx u238-th230 M1 0.01 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxth230onlyM1_001  = new TH1D("hAdapTeO2Sxth230onlyM1_001", "TeO2 Sx th230-only M1 0.01 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxra226pb210M1_001 = new TH1D("hAdapTeO2Sxra226pb210M1_001", "TeO2 Sx ra226-pb210 M1 0.01 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpb210M1_0001     = new TH1D("hAdapTeO2Sxpb210M1_0001", "TeO2 Sx pb210 M1 0.001 #mum",         dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapTeO2Sxth232onlyM1_01  = new TH1D("hAdapTeO2Sxth232onlyM1_01", "TeO2 Sx th232 only M1 0.1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxra228pb208M1_01 = new TH1D("hAdapTeO2Sxra228pb208M1_01", "TeO2 Sx ra228-pb208 M1 0.1 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxu238th230M1_01  = new TH1D("hAdapTeO2Sxu238th230M1_01", "TeO2 Sx u238-th230 M1 0.1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxth230onlyM1_01  = new TH1D("hAdapTeO2Sxth230onlyM1_01", "TeO2 Sx th230-only M1 0.1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxra226pb210M1_01 = new TH1D("hAdapTeO2Sxra226pb210M1_01", "TeO2 Sx ra226-pb210 M1 0.1 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapTeO2Sxth232onlyM1_0001  = new TH1D("hAdapTeO2Sxth232onlyM1_0001", "TeO2 Sx th232 only M1 0.001 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxra228pb208M1_0001 = new TH1D("hAdapTeO2Sxra228pb208M1_0001", "TeO2 Sx ra228-pb208 M1 0.001 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxu238th230M1_0001  = new TH1D("hAdapTeO2Sxu238th230M1_0001", "TeO2 Sx u238-th230 M1 0.001 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxth230onlyM1_0001  = new TH1D("hAdapTeO2Sxth230onlyM1_0001", "TeO2 Sx th230-only M1 0.001 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxra226pb210M1_0001 = new TH1D("hAdapTeO2Sxra226pb210M1_0001", "TeO2 Sx ra226-pb210 M1 0.001 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapTeO20nuM2       = new TH1D("hAdapTeO20nuM2",    "TeO2 Bulk 0nu M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO22nuM2       = new TH1D("hAdapTeO22nuM2",    "TeO2 Bulk 2nu M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2co60M2      = new TH1D("hAdapTeO2co60M2",   "TeO2 Bulk co60 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2k40M2       = new TH1D("hAdapTeO2k40M2",    "TeO2 Bulk k40 M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2pb210M2     = new TH1D("hAdapTeO2pb210M2",  "TeO2 Bulk pb210 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2po210M2     = new TH1D("hAdapTeO2po210M2",  "TeO2 Bulk po210 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2te125M2     = new TH1D("hAdapTeO2te125M2",  "TeO2 Bulk te125 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2th232M2     = new TH1D("hAdapTeO2th232M2",  "TeO2 Bulk th232 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapTeO2th228M2     = new TH1D("hAdapTeO2th228M2",  "TeO2 Bulk th228 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2ra226M2     = new TH1D("hAdapTeO2ra226M2",  "TeO2 Bulk ra226 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2rn222M2     = new TH1D("hAdapTeO2rn222M2",  "TeO2 Bulk rn222 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2u238M2      = new TH1D("hAdapTeO2u238M2",   "TeO2 Bulk u238 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2th230M2     = new TH1D("hAdapTeO2th230M2",  "TeO2 Bulk th230 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2u234M2      = new TH1D("hAdapTeO2u234M2",   "TeO2 Bulk u234 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2sb125M2     = new TH1D("hAdapTeO2sb125M2",  "TeO2 Bulk sb125 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapTeO2Spb210M2_01      = new TH1D("hAdapTeO2Spb210M2_01",   "TeO2 S pb210 M2 0.1 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Spo210M2_001     = new TH1D("hAdapTeO2Spo210M2_001",  "TeO2 S po210 M2 0.01 #mum",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Spo210M2_01      = new TH1D("hAdapTeO2Spo210M2_01",   "TeO2 S po210 M2 0.1 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sth232M2_01      = new TH1D("hAdapTeO2Sth232M2_01",   "TeO2 S th232 M2 0.1 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Su238M2_01       = new TH1D("hAdapTeO2Su238M2_01",    "TeO2 S u238 M2 0.1 #mum",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpb210M2_001    = new TH1D("hAdapTeO2Sxpb210M2_001", "TeO2 Sx pb210 M2 0.01 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpb210M2_01     = new TH1D("hAdapTeO2Sxpb210M2_01",  "TeO2 Sx pb210 M2 0.1 #mum",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpb210M2_1      = new TH1D("hAdapTeO2Sxpb210M2_1",   "TeO2 Sx pb210 M2 1 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpb210M2_10     = new TH1D("hAdapTeO2Sxpb210M2_10",  "TeO2 Sx pb210 M2 10 #mum",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpo210M2_001    = new TH1D("hAdapTeO2Sxpo210M2_001", "TeO2 Sx po210 M2 0.01 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpo210M2_01     = new TH1D("hAdapTeO2Sxpo210M2_01",  "TeO2 Sx po210 M2 0.1 #mum",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpo210M2_1      = new TH1D("hAdapTeO2Sxpo210M2_1",   "TeO2 Sx po210 M2 1",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxth232M2_001    = new TH1D("hAdapTeO2Sxth232M2_001", "TeO2 Sx th232 M2 0.01 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxth232M2_01     = new TH1D("hAdapTeO2Sxth232M2_01",  "TeO2 Sx th232 M2 0.1 #mum",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxth232M2_1      = new TH1D("hAdapTeO2Sxth232M2_1",   "TeO2 Sx th232 M2 1",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxth232M2_10     = new TH1D("hAdapTeO2Sxth232M2_10",  "TeO2 Sx th232 M2 10 #mum",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxu238M2_001     = new TH1D("hAdapTeO2Sxu238M2_001",  "TeO2 Sx u238 M2 0.01 #mum",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxu238M2_01      = new TH1D("hAdapTeO2Sxu238M2_01",   "TeO2 Sx u238 M2 0.1 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxu238M2_1       = new TH1D("hAdapTeO2Sxu238M2_1",    "TeO2 Sx u238 M2 1 #mum",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxu238M2_10      = new TH1D("hAdapTeO2Sxu238M2_10",   "TeO2 Sx u238 M2 10 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapTeO2Sxu238M2_100      = new TH1D("hAdapTeO2Sxu238M2_100",   "TeO2 Sx u238 M2 100 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxth232M2_100     = new TH1D("hAdapTeO2Sxth232M2_100",  "TeO2 Sx th232 M2 100 #mum",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpb210M2_100     = new TH1D("hAdapTeO2Sxpb210M2_100",  "TeO2 Sx pb210 M2 100 #mum",  dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapTeO2th232onlyM2      = new TH1D("hAdapTeO2th232onlyM2",   "TeO2 Bulk th232 only M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2ra228pb208M2     = new TH1D("hAdapTeO2ra228pb208M2",  "TeO2 Bulk ra228-pb208 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2th230onlyM2      = new TH1D("hAdapTeO2th230onlyM2", "TeO2 Bulk th232 only M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2u238th230M2      = new TH1D("hAdapTeO2u238th230M2", "TeO2 Bulk u238-th230 M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2ra226pb210M2     = new TH1D("hAdapTeO2ra226pb210M2", "TeO2 Bulk ra226-pb210 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapTeO2Sxth232onlyM2_001  = new TH1D("hAdapTeO2Sxth232onlyM2_001", "TeO2 Sx th232 only M2 0.01 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxra228pb208M2_001 = new TH1D("hAdapTeO2Sxra228pb208M2_001", "TeO2 Sx ra228-pb208 M2 0.01 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxu238th230M2_001  = new TH1D("hAdapTeO2Sxu238th230M2_001", "TeO2 Sx u238-th230 M2 0.01 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxth230onlyM2_001  = new TH1D("hAdapTeO2Sxth230onlyM2_001", "TeO2 Sx th230 only M2 0.01 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxra226pb210M2_001 = new TH1D("hAdapTeO2Sxra226pb210M2_001", "TeO2 Sx ra226-pb210 M2 0.01 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpb210M2_0001     = new TH1D("hAdapTeO2Sxpb210M2_0001", "TeO2 Sx pb210 M2 0.001 #mum",         dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapTeO2Sxth232onlyM2_01  = new TH1D("hAdapTeO2Sxth232onlyM2_01", "TeO2 Sx th232 only M2 0.1 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxra228pb208M2_01 = new TH1D("hAdapTeO2Sxra228pb208M2_01", "TeO2 Sx ra228-pb208 M2 0.1 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxu238th230M2_01  = new TH1D("hAdapTeO2Sxu238th230M2_01", "TeO2 Sx u238-th230 M2 0.1 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxth230onlyM2_01  = new TH1D("hAdapTeO2Sxth230onlyM2_01", "TeO2 Sx th230-only M2 0.1 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxra226pb210M2_01 = new TH1D("hAdapTeO2Sxra226pb210M2_01", "TeO2 Sx ra226-pb210 M2 0.1 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapTeO2Sxth232onlyM2_0001  = new TH1D("hAdapTeO2Sxth232onlyM2_0001", "TeO2 Sx th232 only M2 0.001 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxra228pb208M2_0001 = new TH1D("hAdapTeO2Sxra228pb208M2_0001", "TeO2 Sx ra228-pb208 M2 0.001 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxu238th230M2_0001  = new TH1D("hAdapTeO2Sxu238th230M2_0001", "TeO2 Sx u238-th230 M2 0.001 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxth230onlyM2_0001  = new TH1D("hAdapTeO2Sxth230onlyM2_0001", "TeO2 Sx th230-only M2 0.001 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxra226pb210M2_0001 = new TH1D("hAdapTeO2Sxra226pb210M2_0001", "TeO2 Sx ra226-pb210 M2 0.001 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);


  hAdapTeO20nuM2Sum       = new TH1D("hAdapTeO20nuM2Sum",    "TeO2 Bulk 0nu M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO22nuM2Sum       = new TH1D("hAdapTeO22nuM2Sum",    "TeO2 Bulk 2nu M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2co60M2Sum      = new TH1D("hAdapTeO2co60M2Sum",   "TeO2 Bulk co60 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2k40M2Sum       = new TH1D("hAdapTeO2k40M2Sum",    "TeO2 Bulk k40 M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2pb210M2Sum     = new TH1D("hAdapTeO2pb210M2Sum",  "TeO2 Bulk pb210 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2po210M2Sum     = new TH1D("hAdapTeO2po210M2Sum",  "TeO2 Bulk po210 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2te125M2Sum     = new TH1D("hAdapTeO2te125M2Sum",  "TeO2 Bulk te125 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2th232M2Sum     = new TH1D("hAdapTeO2th232M2Sum",  "TeO2 Bulk th232 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapTeO2th228M2Sum     = new TH1D("hAdapTeO2th228M2Sum",  "TeO2 Bulk th228 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2ra226M2Sum     = new TH1D("hAdapTeO2ra226M2Sum",  "TeO2 Bulk ra226 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2rn222M2Sum     = new TH1D("hAdapTeO2rn222M2Sum",  "TeO2 Bulk rn222 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2u238M2Sum      = new TH1D("hAdapTeO2u238M2Sum",   "TeO2 Bulk u238 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2th230M2Sum     = new TH1D("hAdapTeO2th230M2Sum",  "TeO2 Bulk th230 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2u234M2Sum      = new TH1D("hAdapTeO2u234M2Sum",   "TeO2 Bulk u234 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2sb125M2Sum     = new TH1D("hAdapTeO2sb125M2Sum",  "TeO2 Bulk sb125 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapTeO2Spb210M2Sum_01      = new TH1D("hAdapTeO2Spb210M2Sum_01",   "TeO2 S pb210 M2Sum 0.1 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Spo210M2Sum_001     = new TH1D("hAdapTeO2Spo210M2Sum_001",  "TeO2 S po210 M2Sum 0.01 #mum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Spo210M2Sum_01      = new TH1D("hAdapTeO2Spo210M2Sum_01",   "TeO2 S po210 M2Sum 0.1 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sth232M2Sum_01      = new TH1D("hAdapTeO2Sth232M2Sum_01",   "TeO2 S th232 M2Sum 0.1 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Su238M2Sum_01       = new TH1D("hAdapTeO2Su238M2Sum_01",    "TeO2 S u238 M2Sum 0.1 #mum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpb210M2Sum_001    = new TH1D("hAdapTeO2Sxpb210M2Sum_001", "TeO2 Sx pb210 M2Sum 0.01 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpb210M2Sum_01     = new TH1D("hAdapTeO2Sxpb210M2Sum_01",  "TeO2 Sx pb210 M2Sum 0.1 #mum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpb210M2Sum_1      = new TH1D("hAdapTeO2Sxpb210M2Sum_1",   "TeO2 Sx pb210 M2Sum 1 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpb210M2Sum_10     = new TH1D("hAdapTeO2Sxpb210M2Sum_10",  "TeO2 Sx pb210 M2Sum 10 #mum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpo210M2Sum_001    = new TH1D("hAdapTeO2Sxpo210M2Sum_001", "TeO2 Sx po210 M2Sum 0.01 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpo210M2Sum_01     = new TH1D("hAdapTeO2Sxpo210M2Sum_01",  "TeO2 Sx po210 M2Sum 0.1 #mum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpo210M2Sum_1      = new TH1D("hAdapTeO2Sxpo210M2Sum_1",   "TeO2 Sx po210 M2Sum 1",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxth232M2Sum_001    = new TH1D("hAdapTeO2Sxth232M2Sum_001", "TeO2 Sx th232 M2Sum 0.01 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxth232M2Sum_01     = new TH1D("hAdapTeO2Sxth232M2Sum_01",  "TeO2 Sx th232 M2Sum 0.1 #mum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxth232M2Sum_1      = new TH1D("hAdapTeO2Sxth232M2Sum_1",   "TeO2 Sx th232 M2Sum 1",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxth232M2Sum_10     = new TH1D("hAdapTeO2Sxth232M2Sum_10",  "TeO2 Sx th232 M2Sum 10 #mum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxu238M2Sum_001     = new TH1D("hAdapTeO2Sxu238M2Sum_001",  "TeO2 Sx u238 M2Sum 0.01 #mum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxu238M2Sum_01      = new TH1D("hAdapTeO2Sxu238M2Sum_01",   "TeO2 Sx u238 M2Sum 0.1 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxu238M2Sum_1       = new TH1D("hAdapTeO2Sxu238M2Sum_1",    "TeO2 Sx u238 M2Sum 1 #mum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxu238M2Sum_10      = new TH1D("hAdapTeO2Sxu238M2Sum_10",   "TeO2 Sx u238 M2Sum 10 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapTeO2Sxu238M2Sum_100      = new TH1D("hAdapTeO2Sxu238M2Sum_100",   "TeO2 Sx u238 M2Sum 100 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxth232M2Sum_100     = new TH1D("hAdapTeO2Sxth232M2Sum_100",  "TeO2 Sx th232 M2Sum 100 #mum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpb210M2Sum_100     = new TH1D("hAdapTeO2Sxpb210M2Sum_100",  "TeO2 Sx pb210 M2Sum 100 #mum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapTeO2th232onlyM2Sum      = new TH1D("hAdapTeO2th232onlyM2Sum",   "TeO2 Bulk th232 only M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2ra228pb208M2Sum     = new TH1D("hAdapTeO2ra228pb208M2Sum",  "TeO2 Bulk ra228-pb208 M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2th230onlyM2Sum      = new TH1D("hAdapTeO2th230onlyM2Sum", "TeO2 Bulk th232 only M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2u238th230M2Sum      = new TH1D("hAdapTeO2u238th230M2Sum", "TeO2 Bulk u238-th230 M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2ra226pb210M2Sum     = new TH1D("hAdapTeO2ra226pb210M2Sum", "TeO2 Bulk ra226-pb210 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapTeO2Sxth232onlyM2Sum_001  = new TH1D("hAdapTeO2Sxth232onlyM2Sum_001", "TeO2 Sx th232 only M2Sum 0.01 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxra228pb208M2Sum_001 = new TH1D("hAdapTeO2Sxra228pb208M2Sum_001", "TeO2 Sx ra228-pb208 M2Sum 0.01 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxu238th230M2Sum_001  = new TH1D("hAdapTeO2Sxu238th230M2Sum_001", "TeO2 Sx u238-th230 M2Sum 0.01 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxth230onlyM2Sum_001  = new TH1D("hAdapTeO2Sxth230onlyM2Sum_001", "TeO2 Sx th230 only M2Sum 0.01 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxra226pb210M2Sum_001 = new TH1D("hAdapTeO2Sxra226pb210M2Sum_001", "TeO2 Sx ra226-pb210 M2Sum 0.01 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpb210M2Sum_0001     = new TH1D("hAdapTeO2Sxpb210M2Sum_0001", "TeO2 Sx pb210 M2Sum 0.001 #mum",         dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapTeO2Sxth232onlyM2Sum_01  = new TH1D("hAdapTeO2Sxth232onlyM2Sum_01", "TeO2 Sx th232 only M2Sum 0.1 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxra228pb208M2Sum_01 = new TH1D("hAdapTeO2Sxra228pb208M2Sum_01", "TeO2 Sx ra228-pb208 M2Sum 0.1 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxu238th230M2Sum_01  = new TH1D("hAdapTeO2Sxu238th230M2Sum_01", "TeO2 Sx u238-th230 M2Sum 0.1 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxth230onlyM2Sum_01  = new TH1D("hAdapTeO2Sxth230onlyM2Sum_01", "TeO2 Sx th230-only M2Sum 0.1 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxra226pb210M2Sum_01 = new TH1D("hAdapTeO2Sxra226pb210M2Sum_01", "TeO2 Sx ra226-pb210 M2Sum 0.1 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapTeO2Sxth232onlyM2Sum_0001  = new TH1D("hAdapTeO2Sxth232onlyM2Sum_0001", "TeO2 Sx th232 only M2Sum 0.001 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxra228pb208M2Sum_0001 = new TH1D("hAdapTeO2Sxra228pb208M2Sum_0001", "TeO2 Sx ra228-pb208 M2Sum 0.001 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxu238th230M2Sum_0001  = new TH1D("hAdapTeO2Sxu238th230M2Sum_0001", "TeO2 Sx u238-th230 M2Sum 0.001 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxth230onlyM2Sum_0001  = new TH1D("hAdapTeO2Sxth230onlyM2Sum_0001", "TeO2 Sx th230-only M2Sum 0.001 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxra226pb210M2Sum_0001 = new TH1D("hAdapTeO2Sxra226pb210M2Sum_0001", "TeO2 Sx ra226-pb210 M2Sum 0.001 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

//////////// CuBox + CuFrame M1 and M2

  hAdapCuBox_CuFrameco60M1 = new TH1D("hAdapCuBox_CuFrameco60M1", "CuBox+CuFrame Bulk co60 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramek40M1 = new TH1D("hAdapCuBox_CuFramek40M1", "CuBox+CuFrame Bulk k40 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameth232M1 = new TH1D("hAdapCuBox_CuFrameth232M1", "CuBox+CuFrame Bulk th232 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameu238M1 = new TH1D("hAdapCuBox_CuFrameu238M1", "CuBox+CuFrame Bulk u238 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuBox_CuFrameth232M1_10 = new TH1D("hAdapCuBox_CuFrameth232M1_10", "CuBox+CuFrame Sx th232 M1 10 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameu238M1_10 = new TH1D("hAdapCuBox_CuFrameu238M1_10", "CuBox+CuFrame Sx u238 M1 10 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramepb210M1_10 = new TH1D("hAdapCuBox_CuFramepb210M1_10", "CuBox+CuFrame Sx pb210 M1 10 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramepb210M1_1 = new TH1D("hAdapCuBox_CuFramepb210M1_1", "CuBox+CuFrame Sx pb210 M1 1 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramepb210M1_01 = new TH1D("hAdapCuBox_CuFramepb210M1_01", "CuBox+CuFrame Sx pb210 M1 0.1 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramepb210M1_001 = new TH1D("hAdapCuBox_CuFramepb210M1_001", "CuBox+CuFrame Sx pb210 M1 0.01 #mum", dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuBox_CuFrameth232M1_1  = new TH1D("hAdapCuBox_CuFrameth232M1_1", "hAdapCuBox_CuFrameth232M1_1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameu238M1_1   = new TH1D("hAdapCuBox_CuFrameu238M1_1", "hAdapCuBox_CuFrameu238M1_1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameth232M1_01  = new TH1D("hAdapCuBox_CuFrameth232M1_01", "hAdapCuBox_CuFrameth232M1_01",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameu238M1_01   = new TH1D("hAdapCuBox_CuFrameu238M1_01", "hAdapCuBox_CuFrameu238M1_01",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameth232M1_001  = new TH1D("hAdapCuBox_CuFrameth232M1_001", "hAdapCuBox_CuFrameth232M1_001",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameu238M1_001   = new TH1D("hAdapCuBox_CuFrameu238M1_001", "hAdapCuBox_CuFrameu238M1_001",    dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuBox_CuFramepb210M1_100  = new TH1D("hAdapCuBox_CuFramepb210M1_100", "hAdapCuBox_CuFramepb210M1_100",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameth232M1_100  = new TH1D("hAdapCuBox_CuFrameth232M1_100", "hAdapCuBox_CuFrameth232M1_100",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameu238M1_100   = new TH1D("hAdapCuBox_CuFrameu238M1_100", "hAdapCuBox_CuFrameu238M1_100",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramepb210M1_50  = new TH1D("hAdapCuBox_CuFramepb210M1_50", "hAdapCuBox_CuFramepb210M1_50",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameth232M1_50  = new TH1D("hAdapCuBox_CuFrameth232M1_50", "hAdapCuBox_CuFrameth232M1_50",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameu238M1_50   = new TH1D("hAdapCuBox_CuFrameu238M1_50", "hAdapCuBox_CuFrameu238M1_50",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramepb210M1_5  = new TH1D("hAdapCuBox_CuFramepb210M1_5", "hAdapCuBox_CuFramepb210M1_5",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameth232M1_5  = new TH1D("hAdapCuBox_CuFrameth232M1_5", "hAdapCuBox_CuFrameth232M1_5",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameu238M1_5   = new TH1D("hAdapCuBox_CuFrameu238M1_5", "hAdapCuBox_CuFrameu238M1_5",    dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuBox_CuFrameco60M2 = new TH1D("hAdapCuBox_CuFrameco60M2", "CuBox+CuFrame Bulk co60 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFramek40M2 = new TH1D("hAdapCuBox_CuFramek40M2", "CuBox+CuFrame Bulk k40 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameth232M2 = new TH1D("hAdapCuBox_CuFrameth232M2", "CuBox+CuFrame Bulk th232 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameu238M2 = new TH1D("hAdapCuBox_CuFrameu238M2", "CuBox+CuFrame Bulk u238 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapCuBox_CuFrameth232M2_10 = new TH1D("hAdapCuBox_CuFrameth232M2_10", "CuBox+CuFrame Sx th232 M2 10 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameu238M2_10 = new TH1D("hAdapCuBox_CuFrameu238M2_10", "CuBox+CuFrame Sx u238 M2 10 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFramepb210M2_10 = new TH1D("hAdapCuBox_CuFramepb210M2_10", "CuBox+CuFrame Sx pb210 M2 10 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFramepb210M2_1 = new TH1D("hAdapCuBox_CuFramepb210M2_1", "CuBox+CuFrame Sx pb210 M2 1 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFramepb210M2_01 = new TH1D("hAdapCuBox_CuFramepb210M2_01", "CuBox+CuFrame Sx pb210 M2 0.1 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFramepb210M2_001 = new TH1D("hAdapCuBox_CuFramepb210M2_001", "CuBox+CuFrame Sx pb210 M2 0.01 #mum", dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapCuBox_CuFrameth232M2_1  = new TH1D("hAdapCuBox_CuFrameth232M2_1", "hAdapCuBox_CuFrameth232M2_1",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameu238M2_1   = new TH1D("hAdapCuBox_CuFrameu238M2_1", "hAdapCuBox_CuFrameu238M2_1",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameth232M2_01  = new TH1D("hAdapCuBox_CuFrameth232M2_01", "hAdapCuBox_CuFrameth232M2_01",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameu238M2_01   = new TH1D("hAdapCuBox_CuFrameu238M2_01", "hAdapCuBox_CuFrameu238M2_01",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameth232M2_001  = new TH1D("hAdapCuBox_CuFrameth232M2_001", "hAdapCuBox_CuFrameth232M2_001",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameu238M2_001   = new TH1D("hAdapCuBox_CuFrameu238M2_001", "hAdapCuBox_CuFrameu238M2_001",    dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapCuBox_CuFramepb210M2_100  = new TH1D("hAdapCuBox_CuFramepb210M2_100", "hAdapCuBox_CuFramepb210M2_100",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameth232M2_100  = new TH1D("hAdapCuBox_CuFrameth232M2_100", "hAdapCuBox_CuFrameth232M2_100",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameu238M2_100   = new TH1D("hAdapCuBox_CuFrameu238M2_100", "hAdapCuBox_CuFrameu238M2_100",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFramepb210M2_50  = new TH1D("hAdapCuBox_CuFramepb210M2_50", "hAdapCuBox_CuFramepb210M2_50",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameth232M2_50  = new TH1D("hAdapCuBox_CuFrameth232M2_50", "hAdapCuBox_CuFrameth232M2_50",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameu238M2_50   = new TH1D("hAdapCuBox_CuFrameu238M2_50", "hAdapCuBox_CuFrameu238M2_50",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFramepb210M2_5  = new TH1D("hAdapCuBox_CuFramepb210M2_5", "hAdapCuBox_CuFramepb210M2_5",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameth232M2_5  = new TH1D("hAdapCuBox_CuFrameth232M2_5", "hAdapCuBox_CuFrameth232M2_5",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameu238M2_5   = new TH1D("hAdapCuBox_CuFrameu238M2_5", "hAdapCuBox_CuFrameu238M2_5",    dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapCuBox_CuFrameco60M2Sum = new TH1D("hAdapCuBox_CuFrameco60M2Sum", "CuBox+CuFrame Bulk co60 M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFramek40M2Sum = new TH1D("hAdapCuBox_CuFramek40M2Sum", "CuBox+CuFrame Bulk k40 M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameth232M2Sum = new TH1D("hAdapCuBox_CuFrameth232M2Sum", "CuBox+CuFrame Bulk th232 M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameu238M2Sum = new TH1D("hAdapCuBox_CuFrameu238M2Sum", "CuBox+CuFrame Bulk u238 M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapCuBox_CuFrameth232M2Sum_10 = new TH1D("hAdapCuBox_CuFrameth232M2Sum_10", "CuBox+CuFrame Sx th232 M2Sum 10 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameu238M2Sum_10 = new TH1D("hAdapCuBox_CuFrameu238M2Sum_10", "CuBox+CuFrame Sx u238 M2Sum 10 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFramepb210M2Sum_10 = new TH1D("hAdapCuBox_CuFramepb210M2Sum_10", "CuBox+CuFrame Sx pb210 M2Sum 10 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFramepb210M2Sum_1 = new TH1D("hAdapCuBox_CuFramepb210M2Sum_1", "CuBox+CuFrame Sx pb210 M2Sum 1 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFramepb210M2Sum_01 = new TH1D("hAdapCuBox_CuFramepb210M2Sum_01", "CuBox+CuFrame Sx pb210 M2Sum 0.1 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFramepb210M2Sum_001 = new TH1D("hAdapCuBox_CuFramepb210M2Sum_001", "CuBox+CuFrame Sx pb210 M2Sum 0.01 #mum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapCuBox_CuFrameth232M2Sum_1  = new TH1D("hAdapCuBox_CuFrameth232M2Sum_1", "hAdapCuBox_CuFrameth232M2Sum_1",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameu238M2Sum_1   = new TH1D("hAdapCuBox_CuFrameu238M2Sum_1", "hAdapCuBox_CuFrameu238M2Sum_1",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameth232M2Sum_01  = new TH1D("hAdapCuBox_CuFrameth232M2Sum_01", "hAdapCuBox_CuFrameth232M2Sum_01",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameu238M2Sum_01   = new TH1D("hAdapCuBox_CuFrameu238M2Sum_01", "hAdapCuBox_CuFrameu238M2Sum_01",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameth232M2Sum_001  = new TH1D("hAdapCuBox_CuFrameth232M2Sum_001", "hAdapCuBox_CuFrameth232M2Sum_001",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameu238M2Sum_001   = new TH1D("hAdapCuBox_CuFrameu238M2Sum_001", "hAdapCuBox_CuFrameu238M2Sum_001",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapCuBox_CuFramepb210M2Sum_100  = new TH1D("hAdapCuBox_CuFramepb210M2Sum_100", "hAdapCuBox_CuFramepb210M2Sum_100",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameth232M2Sum_100  = new TH1D("hAdapCuBox_CuFrameth232M2Sum_100", "hAdapCuBox_CuFrameth232M2Sum_100",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameu238M2Sum_100   = new TH1D("hAdapCuBox_CuFrameu238M2Sum_100", "hAdapCuBox_CuFrameu238M2Sum_100",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFramepb210M2Sum_50  = new TH1D("hAdapCuBox_CuFramepb210M2Sum_50", "hAdapCuBox_CuFramepb210M2Sum_50",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameth232M2Sum_50  = new TH1D("hAdapCuBox_CuFrameth232M2Sum_50", "hAdapCuBox_CuFrameth232M2Sum_50",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameu238M2Sum_50   = new TH1D("hAdapCuBox_CuFrameu238M2Sum_50", "hAdapCuBox_CuFrameu238M2Sum_50",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFramepb210M2Sum_5  = new TH1D("hAdapCuBox_CuFramepb210M2Sum_5", "hAdapCuBox_CuFramepb210M2Sum_5",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameth232M2Sum_5  = new TH1D("hAdapCuBox_CuFrameth232M2Sum_5", "hAdapCuBox_CuFrameth232M2Sum_5",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameu238M2Sum_5   = new TH1D("hAdapCuBox_CuFrameu238M2Sum_5", "hAdapCuBox_CuFrameu238M2Sum_5",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

////////// 50mK M1 and M2
  hAdap50mKcs137M1     = new TH1D("hAdap50mKcs137M1",  "50mK Bulk cs137 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKcs137M2     = new TH1D("hAdap50mKcs137M2",  "50mK Bulk cs137 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap50mKcs137M2Sum     = new TH1D("hAdap50mKcs137M2Sum",  "50mK Bulk cs137 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

////////// Roman Lead M1 and M2
  hAdapPbRombi207M1     = new TH1D("hAdapPbRombi207M1",  "Roman Lead Bulk bi207 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapPbRomco60M1      = new TH1D("hAdapPbRomco60M1",   "Roman Lead Bulk co60 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapPbRomcs137M1     = new TH1D("hAdapPbRomcs137M1",  "Roman Lead Bulk cs137 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapPbRomk40M1       = new TH1D("hAdapPbRomk40M1",    "Roman Lead Bulk k40 M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapPbRompb210M1     = new TH1D("hAdapPbRompb210M1",  "Roman Lead Bulk pb210 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapPbRomth232M1     = new TH1D("hAdapPbRomth232M1",  "Roman Lead Bulk th232 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapPbRomu238M1      = new TH1D("hAdapPbRomu238M1",   "Roman Lead Bulk u238 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapPbRombi207M2     = new TH1D("hAdapPbRombi207M2",  "Roman Lead Bulk bi207 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapPbRomco60M2      = new TH1D("hAdapPbRomco60M2",   "Roman Lead Bulk co60 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapPbRomcs137M2     = new TH1D("hAdapPbRomcs137M2",  "Roman Lead Bulk cs137 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapPbRomk40M2       = new TH1D("hAdapPbRomk40M2",    "Roman Lead Bulk k40 M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapPbRompb210M2     = new TH1D("hAdapPbRompb210M2",  "Roman Lead Bulk pb210 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapPbRomth232M2     = new TH1D("hAdapPbRomth232M2",  "Roman Lead Bulk th232 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapPbRomu238M2      = new TH1D("hAdapPbRomu238M2",   "Roman Lead Bulk u238 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapPbRombi207M2Sum     = new TH1D("hAdapPbRombi207M2Sum",  "Roman Lead Bulk bi207 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapPbRomco60M2Sum      = new TH1D("hAdapPbRomco60M2Sum",   "Roman Lead Bulk co60 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapPbRomcs137M2Sum     = new TH1D("hAdapPbRomcs137M2Sum",  "Roman Lead Bulk cs137 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapPbRomk40M2Sum       = new TH1D("hAdapPbRomk40M2Sum",    "Roman Lead Bulk k40 M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapPbRompb210M2Sum     = new TH1D("hAdapPbRompb210M2Sum",  "Roman Lead Bulk pb210 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapPbRomth232M2Sum     = new TH1D("hAdapPbRomth232M2Sum",  "Roman Lead Bulk th232 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapPbRomu238M2Sum      = new TH1D("hAdapPbRomu238M2Sum",   "Roman Lead Bulk u238 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);


///////// Internal Shields M1 and M2
  hAdapInternalco60M1 = new TH1D("hAdapInternalco60M1", "Internal Shields Bulk co60 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapInternalk40M1 = new TH1D("hAdapInternalk40M1", "Internal Shields Bulk k40 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapInternalth232M1 = new TH1D("hAdapInternalth232M1", "Internal Shields Bulk th232 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapInternalu238M1 = new TH1D("hAdapInternalu238M1", "Internal Shields Bulk u238 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapInternalco60M2 = new TH1D("hAdapInternalco60M2", "Internal Shields Bulk co60 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapInternalk40M2 = new TH1D("hAdapInternalk40M2", "Internal Shields Bulk k40 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapInternalth232M2 = new TH1D("hAdapInternalth232M2", "Internal Shields Bulk th232 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapInternalu238M2 = new TH1D("hAdapInternalu238M2", "Internal Shields Bulk u238 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapInternalco60M2Sum = new TH1D("hAdapInternalco60M2Sum", "Internal Shields Bulk co60 M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapInternalk40M2Sum = new TH1D("hAdapInternalk40M2Sum", "Internal Shields Bulk k40 M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapInternalth232M2Sum = new TH1D("hAdapInternalth232M2Sum", "Internal Shields Bulk th232 M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapInternalu238M2Sum = new TH1D("hAdapInternalu238M2Sum", "Internal Shields Bulk u238 M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  

////////////// OVC M1 and M2
  hAdapOVCco60M1      = new TH1D("hAdapOVCco60M1",   "OVC Bulk co60 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapOVCk40M1       = new TH1D("hAdapOVCk40M1",    "OVC Bulk k40 M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapOVCth232M1     = new TH1D("hAdapOVCth232M1",  "OVC Bulk th232 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapOVCu238M1      = new TH1D("hAdapOVCu238M1",   "OVC Bulk u238 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapOVCco60M2      = new TH1D("hAdapOVCco60M2",   "OVC Bulk co60 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapOVCk40M2       = new TH1D("hAdapOVCk40M2",    "OVC Bulk k40 M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapOVCth232M2     = new TH1D("hAdapOVCth232M2",  "OVC Bulk th232 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapOVCu238M2      = new TH1D("hAdapOVCu238M2",   "OVC Bulk u238 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  

  hAdapOVCco60M2Sum      = new TH1D("hAdapOVCco60M2Sum",   "OVC Bulk co60 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapOVCk40M2Sum       = new TH1D("hAdapOVCk40M2Sum",    "OVC Bulk k40 M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapOVCth232M2Sum     = new TH1D("hAdapOVCth232M2Sum",  "OVC Bulk th232 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapOVCu238M2Sum      = new TH1D("hAdapOVCu238M2Sum",   "OVC Bulk u238 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  

/////////// External Sources M1 and M2
  hAdapExtPbbi210M1 = new TH1D("hAdapExtPbbi210M1", "External Lead Bulk bi210 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapExtPbbi210M2 = new TH1D("hAdapExtPbbi210M2", "External Lead Bulk bi210 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapExtPbbi210M2Sum = new TH1D("hAdapExtPbbi210M2Sum", "External Lead Bulk bi210 M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

/////////// Fudge Factors
  hAdapOVC804M1 = new TH1D("hAdapOVC804M1", "hAdapOVC804M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapOVC835M1 = new TH1D("hAdapOVC835M1", "hAdapOVC835M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapOVC1063M1 = new TH1D("hAdapOVC1063M1", "hAdapOVC1063M1", dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapOVC804M2 = new TH1D("hAdapOVC804M2", "hAdapOVC804M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapOVC835M2 = new TH1D("hAdapOVC835M2", "hAdapOVC835M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapOVC1063M2 = new TH1D("hAdapOVC1063M2", "hAdapOVC1063M2", dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapOVC804M2Sum = new TH1D("hAdapOVC804M2Sum", "hAdapOVC804M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapOVC835M2Sum = new TH1D("hAdapOVC835M2Sum", "hAdapOVC835M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapOVC1063M2Sum = new TH1D("hAdapOVC1063M2Sum", "hAdapOVC1063M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hEnergyScaleDummyM1 = new TH1D("hEnergyScaleDummyM1",   "Energy Scale M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hEnergyScaleDummyM2 = new TH1D("hEnergyScaleDummyM2",   "Energy Scale M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hEnergyScaleDummyM2Sum = new TH1D("hEnergyScaleDummyM2Sum",   "Energy Scale M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hChiSquaredProgressM1 = new TH1D("hChiSquaredProgressM1", "Chi Squared Progress", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hChiSquaredProgressM2 = new TH1D("hChiSquaredProgressM2", "Chi Squared Progress", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hChiSquaredProgressM2Sum = new TH1D("hChiSquaredProgressM2Sum", "Chi Squared Progress", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);


  // Loads all of the PDFs from file
  Initialize();


  // Add 2nbb events to background
/*
  SanityCheck();
  // Set Error Bars to 0 for MC histograms so I don't go crazy
  for(int i = 1; i <= dAdaptiveBinsM1; i++)
  {
    fAdapDataHistoM1->SetBinError(i, 0);
  }

  for(int i = 1; i <= dAdaptiveBinsM2; i++)
  {
    fAdapDataHistoM2->SetBinError(i, 0);
  }

  for(int i = 1; i <= dAdaptiveBinsM2Sum; i++)
  {
    fAdapDataHistoM2Sum->SetBinError(i, 0);
  }
  // Recalculate integrals
  dDataIntegralM1 = fAdapDataHistoM1->Integral("width");
  dDataIntegralM2 = fAdapDataHistoM2->Integral("width");
  dDataIntegralM2Sum = fAdapDataHistoM2Sum->Integral("width");

  cout << "Recalculated integrals: " << endl;
  cout << "Total Events in background spectrum: " << dDataIntegralTot << endl; 
  cout << "Events in background spectrum (M1): " << dDataIntegralM1 << endl;
  cout << "Events in background spectrum (M2): " << dDataIntegralM2 << endl;
  cout << "Events in background spectrum (M2Sum): " << dDataIntegralM2Sum << endl;
  cout << "Livetime of background: " << dLivetimeYr << endl;
*/


  dBestChiSq = 0; // Chi-Squared from best fit (for ProfileNLL calculation)
  // Do the fit now if no other tests are needed 
  nLoop = 0;

  GenerateParameters();
  DoTheFitAdaptive();

  // Try it out
  // ProfileNLL(0,0);
  // ProfileNLL2D(0,0,0);
  // Number of Toy fits
  if(bToyData)ToyFit(1);

}

// Probably needs updating  
TBackgroundModel::~TBackgroundModel()
{
  delete fDataHistoTot;
  delete fDataHistoM1;
  delete fDataHistoM2;
  delete fDataHistoM2Sum;

  delete fAdapDataHistoM1;
  delete fAdapDataHistoM2;
  delete fAdapDataHistoM2Sum;


  // Updated 01-20-2015
  // Total PDFs M1
  delete fModelTotM1;
  delete fModelTotAdapM1;

  // Total PDFs M2
  delete fModelTotM2;
  delete fModelTotAdapM2;

  // Total PDFs M2Sum
  delete fModelTotM2Sum;
  delete fModelTotAdapM2Sum;

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

  delete hTeO2th232onlyM1;
  delete hTeO2ra228pb208M1;
  delete hTeO2th230onlyM1;

  delete hTeO2Sxth232onlyM1_001;
  delete hTeO2Sxra228pb208M1_001;
  delete hTeO2Sxu238th230M1_001;
  delete hTeO2Sxth230onlyM1_001;
  delete hTeO2Sxra226pb210M1_001;
  delete hTeO2Sxpb210M1_0001;

  delete hTeO2Spb210M1_01;
  delete hTeO2Spo210M1_001;
  delete hTeO2Spo210M1_01;
  delete hTeO2Sth232M1_01;
  delete hTeO2Su238M1_01;
  delete hTeO2Sxpb210M1_001;
  delete hTeO2Sxpb210M1_01;
  delete hTeO2Sxpb210M1_1;
  delete hTeO2Sxpb210M1_10;
  delete hTeO2Sxpo210M1_001;
  delete hTeO2Sxpo210M1_01;
  delete hTeO2Sxpo210M1_1;
  delete hTeO2Sxth232M1_001;
  delete hTeO2Sxth232M1_01;
  delete hTeO2Sxth232M1_1;
  delete hTeO2Sxth232M1_10;
  delete hTeO2Sxu238M1_001;
  delete hTeO2Sxu238M1_01;
  delete hTeO2Sxu238M1_1;
  delete hTeO2Sxu238M1_10;

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

  delete hTeO2th232onlyM2;
  delete hTeO2ra228pb208M2;
  delete hTeO2th230onlyM2;

  delete hTeO2Sxth232onlyM2_001;
  delete hTeO2Sxra228pb208M2_001;
  delete hTeO2Sxu238th230M2_001;
  delete hTeO2Sxth230onlyM2_001;
  delete hTeO2Sxra226pb210M2_001;
  delete hTeO2Sxpb210M2_0001;


  delete hTeO2Spb210M2_01;
  delete hTeO2Spo210M2_001;
  delete hTeO2Spo210M2_01;
  delete hTeO2Sth232M2_01;
  delete hTeO2Su238M2_01;
  delete hTeO2Sxpb210M2_001;
  delete hTeO2Sxpb210M2_01;
  delete hTeO2Sxpb210M2_1;
  delete hTeO2Sxpb210M2_10;
  delete hTeO2Sxpo210M2_001;
  delete hTeO2Sxpo210M2_01;
  delete hTeO2Sxpo210M2_1;
  delete hTeO2Sxth232M2_001;
  delete hTeO2Sxth232M2_01;
  delete hTeO2Sxth232M2_1;
  delete hTeO2Sxth232M2_10;
  delete hTeO2Sxu238M2_001;
  delete hTeO2Sxu238M2_01;
  delete hTeO2Sxu238M2_1;
  delete hTeO2Sxu238M2_10;

  delete hTeO20nuM2Sum;
  delete hTeO22nuM2Sum;
  delete hTeO2co60M2Sum;
  delete hTeO2k40M2Sum;
  delete hTeO2pb210M2Sum;
  delete hTeO2po210M2Sum;
  delete hTeO2te125M2Sum;
  delete hTeO2th232M2Sum;
  delete hTeO2th228M2Sum;
  delete hTeO2ra226M2Sum;
  delete hTeO2rn222M2Sum;
  delete hTeO2u238M2Sum;
  delete hTeO2th230M2Sum;
  delete hTeO2u234M2Sum;

  delete hTeO2th232onlyM2Sum;
  delete hTeO2ra228pb208M2Sum;
  delete hTeO2th230onlyM2Sum;

  delete hTeO2Sxth232onlyM2Sum_001;
  delete hTeO2Sxra228pb208M2Sum_001;
  delete hTeO2Sxu238th230M2Sum_001;
  delete hTeO2Sxth230onlyM2Sum_001;
  delete hTeO2Sxra226pb210M2Sum_001;
  delete hTeO2Sxpb210M2Sum_0001;

  delete hTeO2Spb210M2Sum_01;
  delete hTeO2Spo210M2Sum_001;
  delete hTeO2Spo210M2Sum_01;
  delete hTeO2Sth232M2Sum_01;
  delete hTeO2Su238M2Sum_01;
  delete hTeO2Sxpb210M2Sum_001;
  delete hTeO2Sxpb210M2Sum_01;
  delete hTeO2Sxpb210M2Sum_1;
  delete hTeO2Sxpb210M2Sum_10;
  delete hTeO2Sxpo210M2Sum_001;
  delete hTeO2Sxpo210M2Sum_01;
  delete hTeO2Sxpo210M2Sum_1;
  delete hTeO2Sxth232M2Sum_001;
  delete hTeO2Sxth232M2Sum_01;
  delete hTeO2Sxth232M2Sum_1;
  delete hTeO2Sxth232M2Sum_10;
  delete hTeO2Sxu238M2Sum_001;
  delete hTeO2Sxu238M2Sum_01;
  delete hTeO2Sxu238M2Sum_1;
  delete hTeO2Sxu238M2Sum_10;

  delete hCuBox_CuFrameco60M1;
  delete hCuBox_CuFramek40M1;
  delete hCuBox_CuFrameth232M1;
  delete hCuBox_CuFrameu238M1;

  delete hCuBox_CuFrameth232M1_10;
  delete hCuBox_CuFrameu238M1_10;
  delete hCuBox_CuFramepb210M1_10;
  delete hCuBox_CuFramepb210M1_01;

  delete hCuBox_CuFrameco60M2;
  delete hCuBox_CuFramek40M2;
  delete hCuBox_CuFrameth232M2;
  delete hCuBox_CuFrameu238M2;

  delete hCuBox_CuFrameth232M2_10;
  delete hCuBox_CuFrameu238M2_10;
  delete hCuBox_CuFramepb210M2_10;
  delete hCuBox_CuFramepb210M2_01;

  delete hCuBox_CuFrameco60M2Sum;
  delete hCuBox_CuFramek40M2Sum;
  delete hCuBox_CuFrameth232M2Sum;
  delete hCuBox_CuFrameu238M2Sum;

  delete hCuBox_CuFrameth232M2Sum_10;
  delete hCuBox_CuFrameu238M2Sum_10;
  delete hCuBox_CuFramepb210M2Sum_10;
  delete hCuBox_CuFramepb210M2Sum_01;

  delete h50mKcs137M1; 
  delete h50mKcs137M2;
  delete h50mKcs137M2Sum;

  delete hInternalco60M1;
  delete hInternalk40M1;
  delete hInternalth232M1;
  delete hInternalu238M1;

  delete hInternalco60M2;
  delete hInternalk40M2;
  delete hInternalth232M2;
  delete hInternalu238M2;

  delete hInternalco60M2Sum;
  delete hInternalk40M2Sum;
  delete hInternalth232M2Sum;
  delete hInternalu238M2Sum;  

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

  delete hPbRombi207M2Sum;
  delete hPbRomco60M2Sum;
  delete hPbRomcs137M2Sum;
  delete hPbRomk40M2Sum;
  delete hPbRompb210M2Sum;
  delete hPbRomth232M2Sum;
  delete hPbRomu238M2Sum; 

  delete hExtPbbi210M1;
  delete hExtPbbi210M2;
  delete hExtPbbi210M2Sum;

  delete hOVCco60M1;
  delete hOVCk40M1;
  delete hOVCth232M1;
  delete hOVCu238M1;    

  delete hOVCco60M2;
  delete hOVCk40M2;
  delete hOVCth232M2;
  delete hOVCu238M2;  

  delete hOVCco60M2Sum;
  delete hOVCk40M2Sum;
  delete hOVCth232M2Sum;
  delete hOVCu238M2Sum; 

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

  delete hAdapTeO2th232onlyM1;
  delete hAdapTeO2ra228pb208M1;
  delete hAdapTeO2th230onlyM1;

  delete hAdapTeO2Sxth232onlyM1_001;
  delete hAdapTeO2Sxra228pb208M1_001;
  delete hAdapTeO2Sxu238th230M1_001;
  delete hAdapTeO2Sxth230onlyM1_001;
  delete hAdapTeO2Sxra226pb210M1_001;
  delete hAdapTeO2Sxpb210M1_0001;

  delete hAdapTeO2Spb210M1_01;
  delete hAdapTeO2Spo210M1_001;
  delete hAdapTeO2Spo210M1_01;
  delete hAdapTeO2Sth232M1_01;
  delete hAdapTeO2Su238M1_01;
  delete hAdapTeO2Sxpb210M1_001;
  delete hAdapTeO2Sxpb210M1_01;
  delete hAdapTeO2Sxpb210M1_1;
  delete hAdapTeO2Sxpb210M1_10;
  delete hAdapTeO2Sxpo210M1_001;
  delete hAdapTeO2Sxpo210M1_01;
  delete hAdapTeO2Sxpo210M1_1;
  delete hAdapTeO2Sxth232M1_001;
  delete hAdapTeO2Sxth232M1_01;
  delete hAdapTeO2Sxth232M1_1;
  delete hAdapTeO2Sxth232M1_10;
  delete hAdapTeO2Sxu238M1_001;
  delete hAdapTeO2Sxu238M1_01;
  delete hAdapTeO2Sxu238M1_1;
  delete hAdapTeO2Sxu238M1_10;

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

  delete hAdapTeO2th232onlyM2;
  delete hAdapTeO2ra228pb208M2;
  delete hAdapTeO2th230onlyM2;

  delete hAdapTeO2Sxth232onlyM2_001;
  delete hAdapTeO2Sxra228pb208M2_001;
  delete hAdapTeO2Sxu238th230M2_001;
  delete hAdapTeO2Sxth230onlyM2_001;
  delete hAdapTeO2Sxra226pb210M2_001;
  delete hAdapTeO2Sxpb210M2_0001;

  delete hAdapTeO2Spb210M2_01;
  delete hAdapTeO2Spo210M2_001;
  delete hAdapTeO2Spo210M2_01;
  delete hAdapTeO2Sth232M2_01;
  delete hAdapTeO2Su238M2_01;
  delete hAdapTeO2Sxpb210M2_001;
  delete hAdapTeO2Sxpb210M2_01;
  delete hAdapTeO2Sxpb210M2_1;
  delete hAdapTeO2Sxpb210M2_10;
  delete hAdapTeO2Sxpo210M2_001;
  delete hAdapTeO2Sxpo210M2_01;
  delete hAdapTeO2Sxpo210M2_1;
  delete hAdapTeO2Sxth232M2_001;
  delete hAdapTeO2Sxth232M2_01;
  delete hAdapTeO2Sxth232M2_1;
  delete hAdapTeO2Sxth232M2_10;
  delete hAdapTeO2Sxu238M2_001;
  delete hAdapTeO2Sxu238M2_01;
  delete hAdapTeO2Sxu238M2_1;
  delete hAdapTeO2Sxu238M2_10;

  delete hAdapTeO20nuM2Sum;
  delete hAdapTeO22nuM2Sum;
  delete hAdapTeO2co60M2Sum;
  delete hAdapTeO2k40M2Sum;
  delete hAdapTeO2pb210M2Sum;
  delete hAdapTeO2po210M2Sum;
  delete hAdapTeO2te125M2Sum;
  delete hAdapTeO2th232M2Sum;
  delete hAdapTeO2th228M2Sum;
  delete hAdapTeO2ra226M2Sum;
  delete hAdapTeO2rn222M2Sum;
  delete hAdapTeO2u238M2Sum;
  delete hAdapTeO2th230M2Sum;
  delete hAdapTeO2u234M2Sum;

  delete hAdapTeO2th232onlyM2Sum;
  delete hAdapTeO2ra228pb208M2Sum;
  delete hAdapTeO2th230onlyM2Sum;

  delete hAdapTeO2Sxth232onlyM2Sum_001;
  delete hAdapTeO2Sxra228pb208M2Sum_001;
  delete hAdapTeO2Sxu238th230M2Sum_001;
  delete hAdapTeO2Sxth230onlyM2Sum_001;
  delete hAdapTeO2Sxra226pb210M2Sum_001;
  delete hAdapTeO2Sxpb210M2Sum_0001;

  delete hAdapTeO2Spb210M2Sum_01;
  delete hAdapTeO2Spo210M2Sum_001;
  delete hAdapTeO2Spo210M2Sum_01;
  delete hAdapTeO2Sth232M2Sum_01;
  delete hAdapTeO2Su238M2Sum_01;
  delete hAdapTeO2Sxpb210M2Sum_001;
  delete hAdapTeO2Sxpb210M2Sum_01;
  delete hAdapTeO2Sxpb210M2Sum_1;
  delete hAdapTeO2Sxpb210M2Sum_10;
  delete hAdapTeO2Sxpo210M2Sum_001;
  delete hAdapTeO2Sxpo210M2Sum_01;
  delete hAdapTeO2Sxpo210M2Sum_1;
  delete hAdapTeO2Sxth232M2Sum_001;
  delete hAdapTeO2Sxth232M2Sum_01;
  delete hAdapTeO2Sxth232M2Sum_1;
  delete hAdapTeO2Sxth232M2Sum_10;
  delete hAdapTeO2Sxu238M2Sum_001;
  delete hAdapTeO2Sxu238M2Sum_01;
  delete hAdapTeO2Sxu238M2Sum_1;
  delete hAdapTeO2Sxu238M2Sum_10;

  delete hAdapCuBox_CuFrameco60M1;
  delete hAdapCuBox_CuFramek40M1;
  delete hAdapCuBox_CuFrameth232M1;
  delete hAdapCuBox_CuFrameu238M1;

  delete hAdapCuBox_CuFrameth232M1_10;
  delete hAdapCuBox_CuFrameu238M1_10;
  delete hAdapCuBox_CuFramepb210M1_10;
  delete hAdapCuBox_CuFramepb210M1_01;

  delete hAdapCuBox_CuFrameco60M2;
  delete hAdapCuBox_CuFramek40M2;
  delete hAdapCuBox_CuFrameth232M2;
  delete hAdapCuBox_CuFrameu238M2;

  delete hAdapCuBox_CuFrameth232M2_10;
  delete hAdapCuBox_CuFrameu238M2_10;
  delete hAdapCuBox_CuFramepb210M2_10;
  delete hAdapCuBox_CuFramepb210M2_01;

  delete hAdapCuBox_CuFrameco60M2Sum;
  delete hAdapCuBox_CuFramek40M2Sum;
  delete hAdapCuBox_CuFrameth232M2Sum;
  delete hAdapCuBox_CuFrameu238M2Sum;

  delete hAdapCuBox_CuFrameth232M2Sum_10;
  delete hAdapCuBox_CuFrameu238M2Sum_10;
  delete hAdapCuBox_CuFramepb210M2Sum_10;
  delete hAdapCuBox_CuFramepb210M2Sum_01;

  delete hAdap50mKcs137M1;
  delete hAdap50mKcs137M2;
  delete hAdap50mKcs137M2Sum;

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

  delete hAdapPbRombi207M2Sum;
  delete hAdapPbRomco60M2Sum;
  delete hAdapPbRomcs137M2Sum;
  delete hAdapPbRomk40M2Sum;
  delete hAdapPbRompb210M2Sum;
  delete hAdapPbRomth232M2Sum;
  delete hAdapPbRomu238M2Sum; 

  delete hAdapInternalco60M1;
  delete hAdapInternalk40M1;
  delete hAdapInternalth232M1;
  delete hAdapInternalu238M1;

  delete hAdapInternalco60M2;
  delete hAdapInternalk40M2;
  delete hAdapInternalth232M2;
  delete hAdapInternalu238M2;

  delete hAdapInternalco60M2Sum;
  delete hAdapInternalk40M2Sum;
  delete hAdapInternalth232M2Sum;
  delete hAdapInternalu238M2Sum;  

  delete hAdapOVCco60M1;
  delete hAdapOVCk40M1;
  delete hAdapOVCth232M1;
  delete hAdapOVCu238M1;    

  delete hAdapOVCco60M2;
  delete hAdapOVCk40M2;
  delete hAdapOVCth232M2;
  delete hAdapOVCu238M2;  

  delete hAdapOVCco60M2Sum;
  delete hAdapOVCk40M2Sum;
  delete hAdapOVCth232M2Sum;
  delete hAdapOVCu238M2Sum;

  delete hAdapExtPbbi210M1;

  delete hAdapExtPbbi210M2;
  
  delete hAdapExtPbbi210M2Sum;
}

// Returns vector of bin low-edge for adaptive binning
vector<double> TBackgroundModel::AdaptiveBinning(TH1D *h1, int dBinBase)
{
  // dBinBase is the minimal number of counts per bin
  vector<double> dBinArrayThing;

  double dDummy = 0;
  double dDummyFill = 0;
  int j = 0;
  for(int i = h1->FindBin(0); i < h1->FindBin(50); i++)
  {
    dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(i));
  }
  for(int i = h1->FindBin(50); i < h1->GetNbinsX(); i++)
  {
    // Added per each peak
/*
    // K40
    if(i >= h1->FindBin(1445) && i < h1->FindBin(1470))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(1445)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(1470-1445)/dBaseBinSize;
    }    
    // Tl-208
    if(i >= h1->FindBin(2600) && i < h1->FindBin(2630))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(2600)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(2630-2600)/dBaseBinSize;
    }
*/

    // Pt peak 3150 - 3400
    if(i >= h1->FindBin(3050) && i < h1->FindBin(3400))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(3050)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(3400-3050)/dBaseBinSize;
    }
    // 4050 - 4200
    if(i >= h1->FindBin(4025) && i < h1->FindBin(4170))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4025)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(4170-4025)/dBaseBinSize;
    }
    // 4200 - 4350
    if(i >= h1->FindBin(4170) && i < h1->FindBin(4350))
    // if(i >= h1->FindBin(4170) && i < h1->FindBin(4400))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4170)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(4400-4170)/dBaseBinSize;
    }    

    // 4700 - 4850
    if(i >= h1->FindBin(4700) && i < h1->FindBin(4850))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4700)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(4850-4700)/dBaseBinSize;
    }

    // 4850 - 4950
    if(i >= h1->FindBin(4850) && i < h1->FindBin(5000))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4850)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5000-4850)/dBaseBinSize;
    }

    // 5200 - 5400
    if(i >= h1->FindBin(5000) && i < h1->FindBin(5100))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5000)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5100-5000)/dBaseBinSize;
    }

    // 5200 - 5400
    if(i >= h1->FindBin(5200) && i < h1->FindBin(5400))
    // if(i >= h1->FindBin(5170) && i < h1->FindBin(5400))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5200)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5400-5200)/dBaseBinSize;
    }
    // 5400 - 5650
    if(i >= h1->FindBin(5400) && i < h1->FindBin(5650))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5400)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5650-5400)/dBaseBinSize;
    }

    // 5650 - 5800
    if(i >= h1->FindBin(5650) && i < h1->FindBin(5780))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5650)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5780-5650)/dBaseBinSize;
    }

    // 5800 - 6050
    if(i >= h1->FindBin(5780) && i < h1->FindBin(6100))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5780)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(6100-5780)/dBaseBinSize;
    }    

    // 6050 - 6350
    if(i >= h1->FindBin(6100) && i < h1->FindBin(6450))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(6100)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(6450-6100)/dBaseBinSize;
    }    

    // 6700 - 6900
    if(i >= h1->FindBin(6700) && i < h1->FindBin(7000))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(6700)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(7000-6700)/dBaseBinSize;
    }    

    // 7500 - 8000
    if(i >= h1->FindBin(7000) && i < h1->FindBin(8000))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(7000)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(8000-7000)/dBaseBinSize;
    }    

    // 8000 to the end 10000
    if(i >= h1->FindBin(8000) && i < h1->FindBin(10000))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(8000)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(10000-8000)/dBaseBinSize;
    }    


    dDummy = h1->GetBinContent(i);
    dDummyFill += dDummy;


    if(dDummyFill >= dBinBase)
    {
      dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(i-j));
      dDummy = 0;
      dDummyFill = 0;
      j = 0;
    }
    else if(i == h1->GetNbinsX()-1) // for the very end if it doesn't reach 50 events (which it won't)
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
    cout << " Residuals: min bin >= max bin" << endl;
  }

  if(dMult == 1)
  {
  TH1D  *hOut       = new TH1D("hOutResidualM1", "Fit Residuals M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  }
  if(dMult == 2)
  {
  TH1D  *hOut       = new TH1D("hOutResidualM2", "Fit Residuals M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  }
  if(dMult == 3)
  {
  TH1D  *hOut       = new TH1D("hOutResidualM2Sum", "Fit Residuals M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  }

  // Clone histograms for rebinning
  TH1D  *hCloneBkg    = (TH1D*)h1->Clone("hCloneBkg");
  TH1D  *hCloneMC   = (TH1D*)h2->Clone("hCloneMC");

  // Variables used in Residual calculations
  double dResidualX = 0, dResidualY = 0, dResidualXErr = 0, dResidualYErr = 0;

  double datam1_i = 0, modelm1_i = 0;

  // Residual plot and distribution
  for(int j = binMin ; j < binMax ; j++)
  {

    if( hCloneBkg->GetBinCenter(j) >= 3050 && hCloneBkg->GetBinCenter(j) <= 3400)continue;    

    dResidualX    = hCloneBkg->GetBinCenter(j);

    datam1_i = h1->GetBinContent(j)*h1->GetBinWidth(j);
    modelm1_i = h2->GetBinContent(j)*h1->GetBinWidth(j);
    
    // Re-multiply bin content by bin width (for # of counts)
    dResidualY = (datam1_i - modelm1_i)/TMath::Sqrt(datam1_i);

    hOut->SetBinContent(j, dResidualY);
    hOut->SetBinError(j, 1.0);
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

  TLegend *leg = new TLegend(0.75,0.75,0.97,0.97);
  leg->AddEntry(fDataHistoTot, "Total", "l");
  leg->AddEntry(fDataHistoM1, "M1", "l");
  leg->AddEntry(fDataHistoM2, "M2", "l");
  leg->Draw();

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1800);
  c1->Divide(1, 2);
  c1->cd(1);
  gPad->SetLogy();
  fAdapDataHistoM1->SetTitle("M1 Spectrum");
  fAdapDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM1->GetYaxis()->SetTitle("Counts/bin");  
  fAdapDataHistoM1->Draw("");


  c1->cd(2); 
  gPad->SetLogy();
  fAdapDataHistoM2->SetTitle("M2 Spectrum");
  fAdapDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM2->GetYaxis()->SetTitle("Counts/bin");   
  fAdapDataHistoM2->Draw("");

}

double TBackgroundModel::GetChiSquareAdaptive()
{
  double chiSquare = 0.;
  double datam1_i, errm1_i;
  double datam2_i, errm2_i;
  double datam2sum_i, errm2sum_i;
  double modelm1_i, modelm2_i, modelm2sum_i;


  for(int i = dFitMinBinM1; i < dFitMaxBinM1; i++)
  {
    // Dividing by base bin size in chi-squared because the weight is width/base bin size when filling
    if( fAdapDataHistoM1->GetBinCenter(i) >= 3050 && fAdapDataHistoM1->GetBinCenter(i) <= 3400)continue;

    datam1_i = fAdapDataHistoM1->GetBinContent(i)*fAdapDataHistoM1->GetBinWidth(i);
    modelm1_i = fModelTotAdapM1->GetBinContent(i)*fAdapDataHistoM1->GetBinWidth(i);
    
    if(modelm1_i != 0 && datam1_i != 0)
    // if(modelm1_i != 0)
    {
      // Log-likelihood
      chiSquare += 2 * (modelm1_i - datam1_i + datam1_i * TMath::Log(datam1_i/modelm1_i));

      // Pearson chi-squared
      // chiSquare += (datam1_i - modelm1_i)*(datam1_i - modelm1_i)/(datam1_i);
      // Neyman chi-squared
      // chiSquare += (datam1_i - modelm1_i)*(datam1_i - modelm1_i)/modelm1_i;
    }
  }


  for(int i = dFitMinBinM2; i < dFitMaxBinM2; i++)
  {
    if( fAdapDataHistoM2->GetBinCenter(i) >= 3050 && fAdapDataHistoM2->GetBinCenter(i) <= 3400)continue;

    datam2_i = fAdapDataHistoM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i);
    modelm2_i = fModelTotAdapM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i);
    
    if(modelm2_i != 0 && datam2_i != 0)
    // if(modelm2_i != 0)      
    {
      // Log-likelihood
      chiSquare += 2 * (modelm2_i - datam2_i + datam2_i * TMath::Log(datam2_i/modelm2_i));

      // Pearson chi-squared
      // chiSquare += (datam2_i - modelm2_i)*(datam2_i - modelm2_i)/(datam2_i);
      
      // Neyman chi-squared
      // chiSquare += (datam2_i - modelm2_i)*(datam2_i - modelm2_i)/modelm2_i;

    }
  }


/*
  for(int i = dFitMinBinM2Sum; i < dFitMaxBinM2Sum; i++)
  {
    if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 3150 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 3400)continue;
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 800 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 808)continue;
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 1060 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 1068)continue;

    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 506 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 515)continue;
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 579 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 589)continue;
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 605 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 615)continue;
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 906 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 917)continue;
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 1450 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 1475)continue;  
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 1755 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 1780)continue;  
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 2090 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 2130)continue;  
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 2200 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 2220)continue;  
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 2440 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 2450)continue;  
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 2600 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 2630)continue;  

    datam2sum_i = fAdapDataHistoM2Sum->GetBinContent(i)*fAdapDataHistoM2Sum->GetBinWidth(i);
    modelm2sum_i = fModelTotAdapM2Sum->GetBinContent(i)*fAdapDataHistoM2Sum->GetBinWidth(i);
    
    if(modelm2sum_i != 0 && datam2sum_i != 0)
    {
      // Log-likelihood
      chiSquare += 2 * (modelm2sum_i - datam2sum_i + datam2sum_i * TMath::Log(datam2sum_i/modelm2sum_i));

      // Pearson chi-squared
      // chiSquare += (datam2sum_i - modelm2sum_i)*(datam2sum_i - modelm2sum_i)/(datam2sum_i);
      
      // Neyman chi-squared
      // chiSquare += (datam2sum_i - modelm2sum_i)*(datam2sum_i - modelm2sum_i)/modelm2sum_i;

    }
  }
*/

  return chiSquare;
}

void TBackgroundModel::Initialize()
{	

  // Loads PDFs from file
  cout << "Loading PDF Histograms from file" << endl;
  cout << "Directory " << Form("%s/", dMCDir.c_str()) << endl;

  fBulkInner = new TFile(Form("%s/OldProd/MCProduction_BulkInner_1keV.root", dMCDir.c_str()));
  fBulkInnerOld = new TFile(Form("%s/OldProd/MCProduction_BulkInner_1keV.root", dMCDir.c_str()));
  fBulkOuter = new TFile(Form("%s/OldProd/MCProduction_BulkOuter_1keV.root", dMCDir.c_str()));
  fBulkOuterOld = new TFile(Form("%s/OldProd/MCProduction_BulkOuter_1keV.root", dMCDir.c_str()));
  fSurfaceCrystal = new TFile(Form("%s/OldProd/MCProduction_SurfaceCrystal_1keV_new.root", dMCDir.c_str())); // Has more PDFs
  fSurfaceCrystalOld = new TFile(Form("%s/OldProd/MCProduction_SurfaceCrystal_1keV_new.root", dMCDir.c_str())); // Has more PDFs
  fSurfaceOther = new TFile(Form("%s/OldProd/MCProduction_SurfaceOther_1keV.root", dMCDir.c_str()));
  fSurfaceOtherOld = new TFile(Form("%s/OldProd/MCProduction_SurfaceOther_1keV.root", dMCDir.c_str()));

  fBulk = new TFile(Form("%s/MCProduction_Bulk_1keV.root", dMCDir.c_str()));
  fSurface = new TFile(Form("%s/MCProduction_Surface_1keV.root", dMCDir.c_str()));



  fBulkInnerM2Sum = new TFile(Form("%s/OldProd/MCProduction_BulkInnerM2Sum_1keV.root", dMCDir.c_str()));
  fBulkOuterM2Sum = new TFile(Form("%s/OldProd/MCProduction_BulkOuterM2Sum_1keV.root", dMCDir.c_str()));

  fFudge = new TFile(Form("%s/MCProduction_FudgeFactor_1keV.root", dMCDir.c_str()));


///////////// Bulk Histograms
/////// Crystal M1 and M2
  hTeO20nuM1     = (TH1D*)fBulk->Get("hTeO20nuM1");
  hTeO22nuM1     = (TH1D*)fBulk->Get("hTeO22nuM1");
  hTeO2co60M1    = (TH1D*)fBulk->Get("hTeO2co60M1");
  hTeO2k40M1     = (TH1D*)fBulk->Get("hTeO2k40M1");
  hTeO2pb210M1   = (TH1D*)fBulk->Get("hTeO2pb210M1");
  hTeO2po210M1   = (TH1D*)fBulk->Get("hTeO2po210M1");
  hTeO2te125M1   = (TH1D*)fBulk->Get("hTeO2te125M1");
  hTeO2th232M1   = (TH1D*)fBulk->Get("hTeO2th232M1");
  // hTeO2th228M1   = (TH1D*)fBulk->Get("hTeO2th228M1");
  // hTeO2ra226M1   = (TH1D*)fBulk->Get("hTeO2ra226M1");
  // hTeO2rn222M1   = (TH1D*)fBulk->Get("hTeO2rn222M1");
  hTeO2u238M1    = (TH1D*)fBulk->Get("hTeO2u238M1");
  // hTeO2th230M1   = (TH1D*)fBulk->Get("hTeO2th230M1");
  // hTeO2u234M1    = (TH1D*)fBulk->Get("hTeO2u234M1");
  hTeO2sb125M1   = (TH1D*)fBulk->Get("hTeO2sb125M1");

  hTeO2th232onlyM1 = (TH1D*)fBulk->Get("hTeO2th232onlyM1");
  hTeO2ra228pb208M1 = (TH1D*)fBulk->Get("hTeO2ra228pb208M1");
  hTeO2th230onlyM1 = (TH1D*)fBulk->Get("hTeO2th230onlyM1");
  hTeO2u238th230M1 = (TH1D*)fBulk->Get("hTeO2u238th230M1");
  hTeO2ra226pb210M1 = (TH1D*)fBulk->Get("hTeO2ra226pb210M1");

  hTeO20nuM2     = (TH1D*)fBulk->Get("hTeO20nuM2");
  hTeO22nuM2     = (TH1D*)fBulk->Get("hTeO22nuM2");
  hTeO2co60M2    = (TH1D*)fBulk->Get("hTeO2co60M2");
  hTeO2k40M2     = (TH1D*)fBulk->Get("hTeO2k40M2");
  hTeO2pb210M2   = (TH1D*)fBulk->Get("hTeO2pb210M2");
  hTeO2po210M2   = (TH1D*)fBulk->Get("hTeO2po210M2");
  hTeO2te125M2   = (TH1D*)fBulk->Get("hTeO2te125M2");
  hTeO2th232M2   = (TH1D*)fBulk->Get("hTeO2th232M2");
  // hTeO2th228M2   = (TH1D*)fBulk->Get("hTeO2th228M2");
  // hTeO2ra226M2   = (TH1D*)fBulk->Get("hTeO2ra226M2");
  // hTeO2rn222M2   = (TH1D*)fBulk->Get("hTeO2rn222M2");
  hTeO2u238M2    = (TH1D*)fBulk->Get("hTeO2u238M2");
  // hTeO2th230M2   = (TH1D*)fBulk->Get("hTeO2th230M2");
  // hTeO2u234M2    = (TH1D*)fBulk->Get("hTeO2u234M2");
  hTeO2sb125M2   = (TH1D*)fBulk->Get("hTeO2sb125M2");
  hTeO2u238th230M2 = (TH1D*)fBulk->Get("hTeO2u238th230M2");

  hTeO2th232onlyM2 = (TH1D*)fBulk->Get("hTeO2th232onlyM2");
  hTeO2ra228pb208M2 = (TH1D*)fBulk->Get("hTeO2ra228pb208M2");
  hTeO2th230onlyM2 = (TH1D*)fBulk->Get("hTeO2th230onlyM2");
  hTeO2u238th230M2 = (TH1D*)fBulk->Get("hTeO2u238th230M2");
  hTeO2ra226pb210M2 = (TH1D*)fBulk->Get("hTeO2ra226pb210M2");

  hTeO20nuM2Sum     = (TH1D*)fBulkInnerM2Sum->Get("hTeO20nuM2Sum");
  hTeO22nuM2Sum     = (TH1D*)fBulkInnerM2Sum->Get("hTeO22nuM2Sum");
  hTeO2co60M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hTeO2co60M2Sum");
  hTeO2k40M2Sum     = (TH1D*)fBulkInnerM2Sum->Get("hTeO2k40M2Sum");
  hTeO2pb210M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hTeO2pb210M2Sum");
  hTeO2po210M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hTeO2po210M2Sum");
  hTeO2te125M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hTeO2te125M2Sum");
  hTeO2th232M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hTeO2th232M2Sum");
  // hTeO2th228M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hTeO2th228M2Sum");
  // hTeO2ra226M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hTeO2ra226M2Sum");
  // hTeO2rn222M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hTeO2rn222M2Sum");
  hTeO2u238M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hTeO2u238M2Sum");
  // hTeO2th230M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hTeO2th230M2Sum");
  // hTeO2u234M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hTeO2u234M2Sum");

  hTeO2th232onlyM2Sum = (TH1D*)fBulkInnerM2Sum->Get("hTeO2th232onlyM2Sum");
  hTeO2ra228pb208M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hTeO2ra228pb208M2Sum");
  hTeO2th230onlyM2Sum = (TH1D*)fBulkInnerM2Sum->Get("hTeO2th230onlyM2Sum");
  hTeO2u238th230M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hTeO2u238th230M2Sum");
  hTeO2ra226pb210M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hTeO2ra226pb210M2Sum");

///////// CuBox + CuFrame M1 and M2
  hCuBox_CuFrameco60M1 = (TH1D*)fBulk->Get("hCuBox_CuFrameco60M1");
  hCuBox_CuFramek40M1 = (TH1D*)fBulk->Get("hCuBox_CuFramek40M1");
  hCuBox_CuFrameth232M1 = (TH1D*)fBulk->Get("hCuBox_CuFrameth232M1");
  hCuBox_CuFrameu238M1 = (TH1D*)fBulk->Get("hCuBox_CuFrameu238M1");

  hCuBox_CuFrameco60M2 = (TH1D*)fBulk->Get("hCuBox_CuFrameco60M2");
  hCuBox_CuFramek40M2 = (TH1D*)fBulk->Get("hCuBox_CuFramek40M2");
  hCuBox_CuFrameth232M2 = (TH1D*)fBulk->Get("hCuBox_CuFrameth232M2");
  hCuBox_CuFrameu238M2 = (TH1D*)fBulk->Get("hCuBox_CuFrameu238M2");

  hCuBox_CuFrameco60M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hCuBox_CuFrameco60M2Sum");
  hCuBox_CuFramek40M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hCuBox_CuFramek40M2Sum");
  hCuBox_CuFrameth232M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hCuBox_CuFrameth232M2Sum");
  hCuBox_CuFrameu238M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hCuBox_CuFrameu238M2Sum");


///////// 50mK M1 and M2
  h50mKcs137M1   = (TH1D*)fBulkOuterOld->Get("h50mKcs137M1");
  h50mKcs137M2   = (TH1D*)fBulkOuterOld->Get("h50mKcs137M2");
  h50mKcs137M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("h50mKcs137M2Sum");

//////// Roman Lead M1 and M2
  hPbRomk40M1     = (TH1D*)fBulk->Get("hPbRomk40M1");
  hPbRomk40M2     = (TH1D*)fBulk->Get("hPbRomk40M2");
  hPbRomk40M2Sum     = (TH1D*)fBulkOuterM2Sum->Get("hPbRomk40M2Sum");


  hPbRomco60M1    = (TH1D*)fBulkOuterOld->Get("hPbRomco60M1");
  hPbRomcs137M1   = (TH1D*)fBulkOuterOld->Get("hPbRomcs137M1");
  hPbRompb210M1   = (TH1D*)fBulkOuterOld->Get("hPbRompb210M1");
  hPbRomth232M1   = (TH1D*)fBulkOuterOld->Get("hPbRomth232M1");
  hPbRomu238M1    = (TH1D*)fBulkOuterOld->Get("hPbRomu238M1");
   
  hPbRomco60M2    = (TH1D*)fBulkOuterOld->Get("hPbRomco60M2");
  hPbRomcs137M2   = (TH1D*)fBulkOuterOld->Get("hPbRomcs137M2");
  hPbRompb210M2   = (TH1D*)fBulkOuterOld->Get("hPbRompb210M2");
  hPbRomth232M2   = (TH1D*)fBulkOuterOld->Get("hPbRomth232M2");
  hPbRomu238M2    = (TH1D*)fBulkOuterOld->Get("hPbRomu238M2");

  hPbRomco60M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hPbRomco60M2Sum");
  hPbRomcs137M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("hPbRomcs137M2Sum");
  hPbRompb210M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("hPbRompb210M2Sum");
  hPbRomth232M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("hPbRomth232M2Sum");
  hPbRomu238M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hPbRomu238M2Sum");

////// Internal Shields M1 and M2
  // Old Production
  hInternalco60M1 = (TH1D*)fBulk->Get("hInternalco60M1");
  hInternalk40M1 = (TH1D*)fBulk->Get("hInternalk40M1");
  hInternalth232M1 = (TH1D*)fBulk->Get("hInternalth232M1");
  hInternalu238M1 = (TH1D*)fBulk->Get("hInternalu238M1");

  hInternalco60M2 = (TH1D*)fBulk->Get("hInternalco60M2");
  hInternalk40M2 = (TH1D*)fBulk->Get("hInternalk40M2");
  hInternalth232M2 = (TH1D*)fBulk->Get("hInternalth232M2");
  hInternalu238M2 = (TH1D*)fBulk->Get("hInternalu238M2");

  hInternalco60M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hInternalco60M2Sum");
  hInternalk40M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hInternalk40M2Sum");
  hInternalth232M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hInternalth232M2Sum");
  hInternalu238M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hInternalu238M2Sum");

/////// OVC M1 and M2
  hOVCco60M1    = (TH1D*)fBulk->Get("hOVCco60M1");
  hOVCk40M1     = (TH1D*)fBulk->Get("hOVCk40M1");
  hOVCth232M1   = (TH1D*)fBulk->Get("hOVCth232M1");
  hOVCu238M1    = (TH1D*)fBulk->Get("hOVCu238M1");
   
  hOVCco60M2    = (TH1D*)fBulk->Get("hOVCco60M2");
  hOVCk40M2     = (TH1D*)fBulk->Get("hOVCk40M2");
  hOVCth232M2   = (TH1D*)fBulk->Get("hOVCth232M2");
  hOVCu238M2    = (TH1D*)fBulk->Get("hOVCu238M2");

  hOVCco60M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hOVCco60M2Sum");
  hOVCk40M2Sum     = (TH1D*)fBulkOuterM2Sum->Get("hOVCk40M2Sum");
  hOVCth232M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("hOVCth232M2Sum");
  hOVCu238M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hOVCu238M2Sum");

/////// External Sources M1 and M2
  hExtPbbi210M1 = (TH1D*)fBulk->Get("hExtPbbi210M1");
 
  hExtPbbi210M2 = (TH1D*)fBulk->Get("hExtPbbi210M2");

  // hExtPbbi210M2Sum = (TH1D*)fBulkOuterM2Sum->Get("hExtPbbi210M2Sum");

////////// Fudge Factors
  hOVC804M1 = (TH1D*)fFudge->Get("hOVC804M1");
  hOVC835M1 = (TH1D*)fFudge->Get("hOVC835M1");
  hOVC1063M1 = (TH1D*)fFudge->Get("hOVC1063M1");

  hOVC804M2 = (TH1D*)fFudge->Get("hOVC804M2");
  hOVC835M2 = (TH1D*)fFudge->Get("hOVC835M2");
  hOVC1063M2 = (TH1D*)fFudge->Get("hOVC1063M2");
  
  hOVC804M2Sum = (TH1D*)fFudge->Get("hOVC804M2Sum");
  hOVC835M2Sum = (TH1D*)fFudge->Get("hOVC835M2Sum");
  hOVC1063M2Sum = (TH1D*)fFudge->Get("hOVC1063M2Sum");
  
//////////// Surface PDFs
///// Crystal M1 and M2
  hTeO2Sxpb210M1_001  = (TH1D*)fSurface->Get("hTeO2Sxpb210M1_001");
  hTeO2Sxpb210M1_01   = (TH1D*)fSurface->Get("hTeO2Sxpb210M1_01");
  hTeO2Sxpb210M1_1    = (TH1D*)fSurface->Get("hTeO2Sxpb210M1_1");
  hTeO2Sxpb210M1_10   = (TH1D*)fSurface->Get("hTeO2Sxpb210M1_10");
  hTeO2Sxpo210M1_001  = (TH1D*)fSurface->Get("hTeO2Sxpo210M1_001");
  hTeO2Sxpo210M1_01   = (TH1D*)fSurface->Get("hTeO2Sxpo210M1_01");
  hTeO2Sxpo210M1_1    = (TH1D*)fSurface->Get("hTeO2Sxpo210M1_1");
  hTeO2Sxth232M1_001  = (TH1D*)fSurface->Get("hTeO2Sxth232M1_001");
  hTeO2Sxth232M1_01   = (TH1D*)fSurface->Get("hTeO2Sxth232M1_01");
  hTeO2Sxth232M1_1    = (TH1D*)fSurface->Get("hTeO2Sxth232M1_1");
  hTeO2Sxth232M1_10   = (TH1D*)fSurface->Get("hTeO2Sxth232M1_10");
  hTeO2Sxu238M1_001   = (TH1D*)fSurface->Get("hTeO2Sxu238M1_001");
  hTeO2Sxu238M1_01    = (TH1D*)fSurface->Get("hTeO2Sxu238M1_01");
  hTeO2Sxu238M1_1     = (TH1D*)fSurface->Get("hTeO2Sxu238M1_1");
  hTeO2Sxu238M1_10    = (TH1D*)fSurface->Get("hTeO2Sxu238M1_10");

  hTeO2Sxu238M1_100    = (TH1D*)fSurface->Get("hTeO2Sxu238M1_100");
  hTeO2Sxth232M1_100   = (TH1D*)fSurface->Get("hTeO2Sxth232M1_100");
  hTeO2Sxpb210M1_100   = (TH1D*)fSurface->Get("hTeO2Sxpb210M1_100");

  hTeO2Sxth232onlyM1_001 = (TH1D*)fSurface->Get("hTeO2Sxth232onlyM1_001");
  hTeO2Sxra228pb208M1_001 = (TH1D*)fSurface->Get("hTeO2Sxra228pb208M1_001");
  hTeO2Sxu238th230M1_001 = (TH1D*)fSurface->Get("hTeO2Sxu238th230M1_001");
  hTeO2Sxth230onlyM1_001 = (TH1D*)fSurface->Get("hTeO2Sxth230onlyM1_001");
  hTeO2Sxra226pb210M1_001 = (TH1D*)fSurface->Get("hTeO2Sxra226pb210M1_001");
  hTeO2Sxpb210M1_0001 = (TH1D*)fSurface->Get("hTeO2Sxpb210M1_0001");

  hTeO2Sxth232onlyM1_01 = (TH1D*)fSurface->Get("hTeO2Sxth232onlyM1_01");
  hTeO2Sxra228pb208M1_01 = (TH1D*)fSurface->Get("hTeO2Sxra228pb208M1_01");
  hTeO2Sxu238th230M1_01 = (TH1D*)fSurface->Get("hTeO2Sxu238th230M1_01");
  hTeO2Sxth230onlyM1_01 = (TH1D*)fSurface->Get("hTeO2Sxth230onlyM1_01");
  hTeO2Sxra226pb210M1_01 = (TH1D*)fSurface->Get("hTeO2Sxra226pb210M1_01");

  hTeO2Sxth232onlyM1_0001 = (TH1D*)fSurface->Get("hTeO2Sxth232onlyM1_0001");
  hTeO2Sxra228pb208M1_0001 = (TH1D*)fSurface->Get("hTeO2Sxra228pb208M1_0001");
  hTeO2Sxu238th230M1_0001 = (TH1D*)fSurface->Get("hTeO2Sxu238th230M1_0001");
  hTeO2Sxth230onlyM1_0001 = (TH1D*)fSurface->Get("hTeO2Sxth230onlyM1_0001");
  hTeO2Sxra226pb210M1_0001 = (TH1D*)fSurface->Get("hTeO2Sxra226pb210M1_0001");

  hTeO2Sxpb210M2_001  = (TH1D*)fSurface->Get("hTeO2Sxpb210M2_001");
  hTeO2Sxpb210M2_01   = (TH1D*)fSurface->Get("hTeO2Sxpb210M2_01");
  hTeO2Sxpb210M2_1    = (TH1D*)fSurface->Get("hTeO2Sxpb210M2_1");
  hTeO2Sxpb210M2_10   = (TH1D*)fSurface->Get("hTeO2Sxpb210M2_10");
  hTeO2Sxpo210M2_001  = (TH1D*)fSurface->Get("hTeO2Sxpo210M2_001");
  hTeO2Sxpo210M2_01   = (TH1D*)fSurface->Get("hTeO2Sxpo210M2_01");
  hTeO2Sxpo210M2_1    = (TH1D*)fSurface->Get("hTeO2Sxpo210M2_1");
  hTeO2Sxth232M2_001  = (TH1D*)fSurface->Get("hTeO2Sxth232M2_001");
  hTeO2Sxth232M2_01   = (TH1D*)fSurface->Get("hTeO2Sxth232M2_01");
  hTeO2Sxth232M2_1    = (TH1D*)fSurface->Get("hTeO2Sxth232M2_1");
  hTeO2Sxth232M2_10   = (TH1D*)fSurface->Get("hTeO2Sxth232M2_10");
  hTeO2Sxu238M2_001   = (TH1D*)fSurface->Get("hTeO2Sxu238M2_001");
  hTeO2Sxu238M2_01    = (TH1D*)fSurface->Get("hTeO2Sxu238M2_01");
  hTeO2Sxu238M2_1     = (TH1D*)fSurface->Get("hTeO2Sxu238M2_1");
  hTeO2Sxu238M2_10    = (TH1D*)fSurface->Get("hTeO2Sxu238M2_10");

  hTeO2Sxu238M2_100    = (TH1D*)fSurface->Get("hTeO2Sxu238M2_100");
  hTeO2Sxth232M2_100   = (TH1D*)fSurface->Get("hTeO2Sxth232M2_100");
  hTeO2Sxpb210M2_100   = (TH1D*)fSurface->Get("hTeO2Sxpb210M2_100");

  hTeO2Sxth232onlyM2_001 = (TH1D*)fSurface->Get("hTeO2Sxth232onlyM2_001");
  hTeO2Sxra228pb208M2_001 = (TH1D*)fSurface->Get("hTeO2Sxra228pb208M2_001");
  hTeO2Sxu238th230M2_001 = (TH1D*)fSurface->Get("hTeO2Sxu238th230M2_001");
  hTeO2Sxth230onlyM2_001 = (TH1D*)fSurface->Get("hTeO2Sxth230onlyM2_001");
  hTeO2Sxra226pb210M2_001 = (TH1D*)fSurface->Get("hTeO2Sxra226pb210M2_001");
  hTeO2Sxpb210M2_0001 = (TH1D*)fSurface->Get("hTeO2Sxpb210M2_0001");

  hTeO2Sxth232onlyM2_01 = (TH1D*)fSurface->Get("hTeO2Sxth232onlyM2_01");
  hTeO2Sxra228pb208M2_01 = (TH1D*)fSurface->Get("hTeO2Sxra228pb208M2_01");
  hTeO2Sxu238th230M2_01 = (TH1D*)fSurface->Get("hTeO2Sxu238th230M2_01");
  hTeO2Sxth230onlyM2_01 = (TH1D*)fSurface->Get("hTeO2Sxth230onlyM2_01");
  hTeO2Sxra226pb210M2_01 = (TH1D*)fSurface->Get("hTeO2Sxra226pb210M2_01");

  hTeO2Sxth232onlyM2_0001 = (TH1D*)fSurface->Get("hTeO2Sxth232onlyM2_0001");
  hTeO2Sxra228pb208M2_0001 = (TH1D*)fSurface->Get("hTeO2Sxra228pb208M2_0001");
  hTeO2Sxu238th230M2_0001 = (TH1D*)fSurface->Get("hTeO2Sxu238th230M2_0001");
  hTeO2Sxth230onlyM2_0001 = (TH1D*)fSurface->Get("hTeO2Sxth230onlyM2_0001");
  hTeO2Sxra226pb210M2_0001 = (TH1D*)fSurface->Get("hTeO2Sxra226pb210M2_0001");

  hTeO2Sxpb210M2Sum_001  = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M2Sum_001");
  hTeO2Sxpb210M2Sum_01   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M2Sum_01");
  hTeO2Sxpb210M2Sum_1    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M2Sum_1");
  hTeO2Sxpb210M2Sum_10   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M2Sum_10");
  hTeO2Sxpo210M2Sum_001  = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpo210M2Sum_001");
  hTeO2Sxpo210M2Sum_01   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpo210M2Sum_01");
  hTeO2Sxpo210M2Sum_1    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpo210M2Sum_1");
  hTeO2Sxth232M2Sum_001  = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M2Sum_001");
  hTeO2Sxth232M2Sum_01   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M2Sum_01");
  hTeO2Sxth232M2Sum_1    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M2Sum_1");
  hTeO2Sxth232M2Sum_10   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M2Sum_10");
  hTeO2Sxu238M2Sum_001   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M2Sum_001");
  hTeO2Sxu238M2Sum_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M2Sum_01");
  hTeO2Sxu238M2Sum_1     = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M2Sum_1");
  hTeO2Sxu238M2Sum_10    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M2Sum_10");

  hTeO2Sxu238M2Sum_100    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M2Sum_100");
  hTeO2Sxth232M2Sum_100   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M2Sum_100");
  hTeO2Sxpb210M2Sum_100   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M2Sum_100");

  hTeO2Sxth232onlyM2Sum_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232onlyM2Sum_001");
  hTeO2Sxra228pb208M2Sum_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra228pb208M2Sum_001");
  hTeO2Sxu238th230M2Sum_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238th230M2Sum_001");
  hTeO2Sxth230onlyM2Sum_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth230onlyM2Sum_001");
  hTeO2Sxra226pb210M2Sum_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra226pb210M2Sum_001");
  hTeO2Sxpb210M2Sum_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M2Sum_0001");

  hTeO2Sxth232onlyM2Sum_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232onlyM2Sum_01");
  hTeO2Sxra228pb208M2Sum_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra228pb208M2Sum_01");
  hTeO2Sxu238th230M2Sum_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238th230M2Sum_01");
  hTeO2Sxth230onlyM2Sum_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth230onlyM2Sum_01");
  hTeO2Sxra226pb210M2Sum_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra226pb210M2Sum_01");

  hTeO2Sxth232onlyM2Sum_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232onlyM2Sum_0001");
  hTeO2Sxra228pb208M2Sum_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra228pb208M2Sum_0001");
  hTeO2Sxu238th230M2Sum_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238th230M2Sum_0001");
  hTeO2Sxth230onlyM2Sum_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth230onlyM2Sum_0001");
  hTeO2Sxra226pb210M2Sum_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra226pb210M2Sum_0001");

/////// CuBox + CuFrame

  hCuBox_CuFrameth232M1_10 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M1_10");
  hCuBox_CuFrameu238M1_10 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M1_10");
  hCuBox_CuFramepb210M1_10 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M1_10");
  hCuBox_CuFramepb210M1_1 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M1_1");
  hCuBox_CuFramepb210M1_01 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M1_01");
  hCuBox_CuFramepb210M1_001 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M1_001");

  hCuBox_CuFrameth232M1_1 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M1_1");
  hCuBox_CuFrameu238M1_1 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M1_1");
  hCuBox_CuFrameth232M1_01 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M1_01");
  hCuBox_CuFrameu238M1_01 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M1_01");
  hCuBox_CuFrameth232M1_001 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M1_001");
  hCuBox_CuFrameu238M1_001 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M1_001");  


  hCuBox_CuFrameth232M1_100 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M1_100");
  hCuBox_CuFrameu238M1_100 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M1_100");
  hCuBox_CuFramepb210M1_100 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M1_100");
  hCuBox_CuFrameth232M1_50 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M1_50");
  hCuBox_CuFrameu238M1_50 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M1_50");
  hCuBox_CuFramepb210M1_50 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M1_50");
  hCuBox_CuFrameth232M1_5 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M1_5");
  hCuBox_CuFrameu238M1_5 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M1_5");
  hCuBox_CuFramepb210M1_5 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M1_5");

  hCuBox_CuFrameth232M2_10 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M2_10");
  hCuBox_CuFrameu238M2_10 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M2_10");
  hCuBox_CuFramepb210M2_10 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M2_10");
  hCuBox_CuFramepb210M2_1 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M2_1");
  hCuBox_CuFramepb210M2_01 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M2_01");
  hCuBox_CuFramepb210M2_001 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M2_001");

  hCuBox_CuFrameth232M2_1 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M2_1");
  hCuBox_CuFrameu238M2_1 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M2_1");
  hCuBox_CuFrameth232M2_01 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M2_01");
  hCuBox_CuFrameu238M2_01 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M2_01");
  hCuBox_CuFrameth232M2_001 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M2_001");
  hCuBox_CuFrameu238M2_001 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M2_001"); 

  hCuBox_CuFrameth232M2_100 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M2_100");
  hCuBox_CuFrameu238M2_100 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M2_100");
  hCuBox_CuFramepb210M2_100 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M2_100");
  hCuBox_CuFrameth232M2_50 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M2_50");
  hCuBox_CuFrameu238M2_50 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M2_50");
  hCuBox_CuFramepb210M2_50 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M2_50");
  hCuBox_CuFrameth232M2_5 = (TH1D*)fSurface->Get("hCuBox_CuFrameth232M2_5");
  hCuBox_CuFrameu238M2_5 = (TH1D*)fSurface->Get("hCuBox_CuFrameu238M2_5");
  hCuBox_CuFramepb210M2_5 = (TH1D*)fSurface->Get("hCuBox_CuFramepb210M2_5");

  hCuBox_CuFrameth232M2Sum_10 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2Sum_10");
  hCuBox_CuFrameu238M2Sum_10 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameu238M2Sum_10");
  hCuBox_CuFramepb210M2Sum_10 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M2Sum_10");
  hCuBox_CuFramepb210M2Sum_1 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M2Sum_1");
  hCuBox_CuFramepb210M2Sum_01 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M2Sum_01");
  hCuBox_CuFramepb210M2Sum_001 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M2Sum_001"); 

  hCuBox_CuFrameth232M2Sum_1 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2Sum_1");
  hCuBox_CuFrameu238M2Sum_1 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M2Sum_1");
  hCuBox_CuFrameth232M2Sum_01 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2Sum_01");
  hCuBox_CuFrameu238M2Sum_01 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M2Sum_01");
  hCuBox_CuFrameth232M2Sum_001 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2Sum_001");
  hCuBox_CuFrameu238M2Sum_001 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M2Sum_001"); 

  hCuBox_CuFrameth232M2Sum_100 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2Sum_100");
  hCuBox_CuFrameu238M2Sum_100 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M2Sum_100");
  hCuBox_CuFramepb210M2Sum_100 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M2Sum_100");
  hCuBox_CuFrameth232M2Sum_50 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2Sum_50");
  hCuBox_CuFrameu238M2Sum_50 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M2Sum_50");
  hCuBox_CuFramepb210M2Sum_50 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M2Sum_50");
  hCuBox_CuFrameth232M2Sum_5 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2Sum_5");
  hCuBox_CuFrameu238M2Sum_5 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M2Sum_5");
  hCuBox_CuFramepb210M2Sum_5 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M2Sum_5");

///////////// Get adaptive binned histograms
//////// Crystal M1 and M2
  hnewTeO20nuM1 = hTeO20nuM1->Rebin(dAdaptiveBinsM1, "hnewTeO20nuM1", dAdaptiveArrayM1);
  hnewTeO22nuM1 = hTeO22nuM1->Rebin(dAdaptiveBinsM1, "hnewTeO22nuM1", dAdaptiveArrayM1);
  hnewTeO2co60M1 = hTeO2co60M1->Rebin(dAdaptiveBinsM1, "hnewTeO2co60M1", dAdaptiveArrayM1);
  hnewTeO2k40M1 = hTeO2k40M1->Rebin(dAdaptiveBinsM1, "hnewTeO2k40M1", dAdaptiveArrayM1);
  hnewTeO2pb210M1 = hTeO2pb210M1->Rebin(dAdaptiveBinsM1, "hnewTeO2pb210M1", dAdaptiveArrayM1);
  hnewTeO2po210M1 = hTeO2po210M1->Rebin(dAdaptiveBinsM1, "hnewTeO2po210M1", dAdaptiveArrayM1);
  hnewTeO2te125M1 = hTeO2te125M1->Rebin(dAdaptiveBinsM1, "hnewTeO2te125M1", dAdaptiveArrayM1);
  hnewTeO2th232M1 = hTeO2th232M1->Rebin(dAdaptiveBinsM1, "hnewTeO2th232M1", dAdaptiveArrayM1);
  hnewTeO2u238M1 = hTeO2u238M1->Rebin(dAdaptiveBinsM1, "hnewTeO2u238M1", dAdaptiveArrayM1);
  hnewTeO2sb125M1 = hTeO2sb125M1->Rebin(dAdaptiveBinsM1, "hnewTeO2sb125M1", dAdaptiveArrayM1);

  hnewTeO2th232onlyM1 = hTeO2th232onlyM1->Rebin(dAdaptiveBinsM1, "hnewTeO2th232onlyM1", dAdaptiveArrayM1);
  hnewTeO2ra228pb208M1 = hTeO2ra228pb208M1->Rebin(dAdaptiveBinsM1, "hnewTeO2ra228pb208M1", dAdaptiveArrayM1);
  hnewTeO2th230onlyM1 = hTeO2th230onlyM1->Rebin(dAdaptiveBinsM1, "hnewTeO2th230onlyM1", dAdaptiveArrayM1);
  hnewTeO2u238th230M1 = hTeO2u238th230M1->Rebin(dAdaptiveBinsM1, "hnewTeO2u238th230M1", dAdaptiveArrayM1);
  hnewTeO2ra226pb210M1 = hTeO2ra226pb210M1->Rebin(dAdaptiveBinsM1, "hnewTeO2ra226pb210M1", dAdaptiveArrayM1);

  hnewTeO2Sxpb210M1_001 = hTeO2Sxpb210M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpb210M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxpb210M1_01 = hTeO2Sxpb210M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpb210M1_01", dAdaptiveArrayM1);
  hnewTeO2Sxpb210M1_1 = hTeO2Sxpb210M1_1->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpb210M1_1", dAdaptiveArrayM1);
  hnewTeO2Sxpb210M1_10 = hTeO2Sxpb210M1_10->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpb210M1_10", dAdaptiveArrayM1);
  hnewTeO2Sxpo210M1_001 = hTeO2Sxpo210M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpo210M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxpo210M1_01 = hTeO2Sxpo210M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpo210M1_01", dAdaptiveArrayM1);
  hnewTeO2Sxpo210M1_1 = hTeO2Sxpo210M1_1->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpo210M1_1", dAdaptiveArrayM1);
  hnewTeO2Sxth232M1_001 = hTeO2Sxth232M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxth232M1_01 = hTeO2Sxth232M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232M1_01", dAdaptiveArrayM1);
  hnewTeO2Sxth232M1_1 = hTeO2Sxth232M1_1->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232M1_1", dAdaptiveArrayM1);
  hnewTeO2Sxth232M1_10 = hTeO2Sxth232M1_10->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232M1_10", dAdaptiveArrayM1);
  hnewTeO2Sxu238M1_001 = hTeO2Sxu238M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxu238M1_01 = hTeO2Sxu238M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238M1_01", dAdaptiveArrayM1);
  hnewTeO2Sxu238M1_1  = hTeO2Sxu238M1_1->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238M1_1", dAdaptiveArrayM1);
  hnewTeO2Sxu238M1_10 = hTeO2Sxu238M1_10->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238M1_10", dAdaptiveArrayM1);

  hnewTeO2Sxu238M1_100 = hTeO2Sxu238M1_100->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238M1_100", dAdaptiveArrayM1);
  hnewTeO2Sxth232M1_100 = hTeO2Sxth232M1_100->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232M1_100", dAdaptiveArrayM1);
  hnewTeO2Sxpb210M1_100 = hTeO2Sxpb210M1_100->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpb210M1_100", dAdaptiveArrayM1);


  hnewTeO2th232onlyM1 = hTeO2th232onlyM1->Rebin(dAdaptiveBinsM1, "hnewTeO2th232onlyM1", dAdaptiveArrayM1);
  hnewTeO2ra228pb208M1 = hTeO2ra228pb208M1->Rebin(dAdaptiveBinsM1, "hnewTeO2ra228pb208M1", dAdaptiveArrayM1);
  hnewTeO2th230onlyM1 = hTeO2th230onlyM1->Rebin(dAdaptiveBinsM1, "hnewTeO2th230onlyM1", dAdaptiveArrayM1);

  hnewTeO2Sxth232onlyM1_001 = hTeO2Sxth232onlyM1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232onlyM1_001", dAdaptiveArrayM1);
  hnewTeO2Sxra228pb208M1_001 = hTeO2Sxra228pb208M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxra228pb208M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxu238th230M1_001 = hTeO2Sxu238th230M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238th230M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxth230onlyM1_001 = hTeO2Sxth230onlyM1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth230onlyM1_001", dAdaptiveArrayM1);
  hnewTeO2Sxra226pb210M1_001 = hTeO2Sxra226pb210M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxra226pb210M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxpb210M1_0001 = hTeO2Sxpb210M1_0001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpb210M1_0001", dAdaptiveArrayM1);

  hnewTeO2Sxth232onlyM1_01 = hTeO2Sxth232onlyM1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232onlyM1_01", dAdaptiveArrayM1);
  hnewTeO2Sxra228pb208M1_01 = hTeO2Sxra228pb208M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxra228pb208M1_01", dAdaptiveArrayM1);
  hnewTeO2Sxu238th230M1_01 = hTeO2Sxu238th230M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238th230M1_01", dAdaptiveArrayM1);
  hnewTeO2Sxth230onlyM1_01 = hTeO2Sxth230onlyM1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth230onlyM1_01", dAdaptiveArrayM1);
  hnewTeO2Sxra226pb210M1_01 = hTeO2Sxra226pb210M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxra226pb210M1_01", dAdaptiveArrayM1);

  hnewTeO2Sxth232onlyM1_0001 = hTeO2Sxth232onlyM1_0001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232onlyM1_0001", dAdaptiveArrayM1);
  hnewTeO2Sxra228pb208M1_0001 = hTeO2Sxra228pb208M1_0001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxra228pb208M1_0001", dAdaptiveArrayM1);
  hnewTeO2Sxu238th230M1_0001 = hTeO2Sxu238th230M1_0001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238th230M1_0001", dAdaptiveArrayM1);
  hnewTeO2Sxth230onlyM1_0001 = hTeO2Sxth230onlyM1_0001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth230onlyM1_0001", dAdaptiveArrayM1);
  hnewTeO2Sxra226pb210M1_0001 = hTeO2Sxra226pb210M1_0001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxra226pb210M1_0001", dAdaptiveArrayM1);  

  hnewTeO20nuM2 = hTeO20nuM2->Rebin(dAdaptiveBinsM2, "hnewTeO20nuM2", dAdaptiveArrayM2);
  hnewTeO22nuM2 = hTeO22nuM2->Rebin(dAdaptiveBinsM2, "hnewTeO22nuM2", dAdaptiveArrayM2);
  hnewTeO2co60M2 = hTeO2co60M2->Rebin(dAdaptiveBinsM2, "hnewTeO2co60M2", dAdaptiveArrayM2);
  hnewTeO2k40M2 = hTeO2k40M2->Rebin(dAdaptiveBinsM2, "hnewTeO2k40M2", dAdaptiveArrayM2);
  hnewTeO2pb210M2 = hTeO2pb210M2->Rebin(dAdaptiveBinsM2, "hnewTeO2pb210M2", dAdaptiveArrayM2);
  hnewTeO2po210M2 = hTeO2po210M2->Rebin(dAdaptiveBinsM2, "hnewTeO2po210M2", dAdaptiveArrayM2);
  hnewTeO2te125M2 = hTeO2te125M2->Rebin(dAdaptiveBinsM2, "hnewTeO2te125M2", dAdaptiveArrayM2);
  hnewTeO2th232M2 = hTeO2th232M2->Rebin(dAdaptiveBinsM2, "hnewTeO2th232M2", dAdaptiveArrayM2);
  hnewTeO2u238M2 = hTeO2u238M2->Rebin(dAdaptiveBinsM2, "hnewTeO2u238M2", dAdaptiveArrayM2);
  hnewTeO2sb125M2 = hTeO2sb125M2->Rebin(dAdaptiveBinsM2, "hnewTeO2sb125M2", dAdaptiveArrayM2);

  hnewTeO2th232onlyM2 = hTeO2th232onlyM2->Rebin(dAdaptiveBinsM2, "hnewTeO2th232onlyM2", dAdaptiveArrayM2);
  hnewTeO2ra228pb208M2 = hTeO2ra228pb208M2->Rebin(dAdaptiveBinsM2, "hnewTeO2ra228pb208M2", dAdaptiveArrayM2);
  hnewTeO2th230onlyM2 = hTeO2th230onlyM2->Rebin(dAdaptiveBinsM2, "hnewTeO2th230onlyM2", dAdaptiveArrayM2);
  hnewTeO2u238th230M2 = hTeO2u238th230M2->Rebin(dAdaptiveBinsM2, "hnewTeO2u238th230M2", dAdaptiveArrayM2);
  hnewTeO2ra226pb210M2 = hTeO2ra226pb210M2->Rebin(dAdaptiveBinsM2, "hnewTeO2ra226pb210M2", dAdaptiveArrayM2);

  hnewTeO2Sxpb210M2_001 = hTeO2Sxpb210M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpb210M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxpb210M2_01 = hTeO2Sxpb210M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpb210M2_01", dAdaptiveArrayM2);
  hnewTeO2Sxpb210M2_1 = hTeO2Sxpb210M2_1->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpb210M2_1", dAdaptiveArrayM2);
  hnewTeO2Sxpb210M2_10 = hTeO2Sxpb210M2_10->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpb210M2_10", dAdaptiveArrayM2);
  hnewTeO2Sxpo210M2_001 = hTeO2Sxpo210M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpo210M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxpo210M2_01 = hTeO2Sxpo210M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpo210M2_01", dAdaptiveArrayM2);
  hnewTeO2Sxpo210M2_1 = hTeO2Sxpo210M2_1->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpo210M2_1", dAdaptiveArrayM2);
  hnewTeO2Sxth232M2_001 = hTeO2Sxth232M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxth232M2_01 = hTeO2Sxth232M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232M2_01", dAdaptiveArrayM2);
  hnewTeO2Sxth232M2_1 = hTeO2Sxth232M2_1->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232M2_1", dAdaptiveArrayM2);
  hnewTeO2Sxth232M2_10 = hTeO2Sxth232M2_10->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232M2_10", dAdaptiveArrayM2);
  hnewTeO2Sxu238M2_001 = hTeO2Sxu238M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxu238M2_01 = hTeO2Sxu238M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238M2_01", dAdaptiveArrayM2);
  hnewTeO2Sxu238M2_1 = hTeO2Sxu238M2_1->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238M2_1", dAdaptiveArrayM2);
  hnewTeO2Sxu238M2_10 = hTeO2Sxu238M2_10->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238M2_10", dAdaptiveArrayM2);

  hnewTeO2Sxu238M2_100 = hTeO2Sxu238M2_100->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238M2_100", dAdaptiveArrayM2);
  hnewTeO2Sxth232M2_100 = hTeO2Sxth232M2_100->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232M2_100", dAdaptiveArrayM2);
  hnewTeO2Sxpb210M2_100 = hTeO2Sxpb210M2_100->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpb210M2_100", dAdaptiveArrayM2);

  hnewTeO2th232onlyM2 = hTeO2th232onlyM2->Rebin(dAdaptiveBinsM2, "hnewTeO2th232onlyM2", dAdaptiveArrayM2);
  hnewTeO2ra228pb208M2 = hTeO2ra228pb208M2->Rebin(dAdaptiveBinsM2, "hnewTeO2ra228pb208M2", dAdaptiveArrayM2);
  hnewTeO2th230onlyM2 = hTeO2th230onlyM2->Rebin(dAdaptiveBinsM2, "hnewTeO2th230onlyM2", dAdaptiveArrayM2);

  hnewTeO2Sxth232onlyM2_001 = hTeO2Sxth232onlyM2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232onlyM2_001", dAdaptiveArrayM2);
  hnewTeO2Sxra228pb208M2_001 = hTeO2Sxra228pb208M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxra228pb208M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxu238th230M2_001 = hTeO2Sxu238th230M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238th230M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxth230onlyM2_001 = hTeO2Sxth230onlyM2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth230onlyM2_001", dAdaptiveArrayM2);
  hnewTeO2Sxra226pb210M2_001 = hTeO2Sxra226pb210M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxra226pb210M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxpb210M2_0001 = hTeO2Sxpb210M2_0001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpb210M2_0001", dAdaptiveArrayM2);

  hnewTeO2Sxth232onlyM2_01 = hTeO2Sxth232onlyM2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232onlyM2_01", dAdaptiveArrayM2);
  hnewTeO2Sxra228pb208M2_01 = hTeO2Sxra228pb208M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxra228pb208M2_01", dAdaptiveArrayM2);
  hnewTeO2Sxu238th230M2_01 = hTeO2Sxu238th230M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238th230M2_01", dAdaptiveArrayM2);
  hnewTeO2Sxth230onlyM2_01 = hTeO2Sxth230onlyM2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth230onlyM2_01", dAdaptiveArrayM2);
  hnewTeO2Sxra226pb210M2_01 = hTeO2Sxra226pb210M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxra226pb210M2_01", dAdaptiveArrayM2);

  hnewTeO2Sxth232onlyM2_0001 = hTeO2Sxth232onlyM2_0001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232onlyM2_0001", dAdaptiveArrayM2);
  hnewTeO2Sxra228pb208M2_0001 = hTeO2Sxra228pb208M2_0001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxra228pb208M2_0001", dAdaptiveArrayM2);
  hnewTeO2Sxu238th230M2_0001 = hTeO2Sxu238th230M2_0001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238th230M2_0001", dAdaptiveArrayM2);
  hnewTeO2Sxth230onlyM2_0001 = hTeO2Sxth230onlyM2_0001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth230onlyM2_0001", dAdaptiveArrayM2);
  hnewTeO2Sxra226pb210M2_0001 = hTeO2Sxra226pb210M2_0001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxra226pb210M2_0001", dAdaptiveArrayM2);  

  hnewTeO20nuM2Sum = hTeO20nuM2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO20nuM2Sum", dAdaptiveArrayM2Sum);
  hnewTeO22nuM2Sum = hTeO22nuM2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO22nuM2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2co60M2Sum = hTeO2co60M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2co60M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2k40M2Sum = hTeO2k40M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2k40M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2pb210M2Sum = hTeO2pb210M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2pb210M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2po210M2Sum = hTeO2po210M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2po210M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2te125M2Sum = hTeO2te125M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2te125M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2th232M2Sum = hTeO2th232M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2th232M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2u238M2Sum = hTeO2u238M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2u238M2Sum", dAdaptiveArrayM2Sum);

  hnewTeO2th232onlyM2Sum = hTeO2th232onlyM2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2th232onlyM2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2ra228pb208M2Sum = hTeO2ra228pb208M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2ra228pb208M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2th230onlyM2Sum = hTeO2th230onlyM2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2th230onlyM2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2u238th230M2Sum = hTeO2u238th230M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2u238th230M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2ra226pb210M2Sum = hTeO2ra226pb210M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2ra226pb210M2Sum", dAdaptiveArrayM2Sum);

  hnewTeO2Sxpb210M2Sum_001 = hTeO2Sxpb210M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxpb210M2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxpb210M2Sum_01 = hTeO2Sxpb210M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxpb210M2Sum_01", dAdaptiveArrayM2Sum);
  hnewTeO2Sxpb210M2Sum_1 = hTeO2Sxpb210M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxpb210M2Sum_1", dAdaptiveArrayM2Sum);
  hnewTeO2Sxpb210M2Sum_10 = hTeO2Sxpb210M2Sum_10->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxpb210M2Sum_10", dAdaptiveArrayM2Sum);
  hnewTeO2Sxpo210M2Sum_001 = hTeO2Sxpo210M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxpo210M2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxpo210M2Sum_01 = hTeO2Sxpo210M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxpo210M2Sum_01", dAdaptiveArrayM2Sum);
  hnewTeO2Sxpo210M2Sum_1 = hTeO2Sxpo210M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxpo210M2Sum_1", dAdaptiveArrayM2Sum);
  hnewTeO2Sxth232M2Sum_001 = hTeO2Sxth232M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxth232M2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxth232M2Sum_01 = hTeO2Sxth232M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxth232M2Sum_01", dAdaptiveArrayM2Sum);
  hnewTeO2Sxth232M2Sum_1 = hTeO2Sxth232M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxth232M2Sum_1", dAdaptiveArrayM2Sum);
  hnewTeO2Sxth232M2Sum_10 = hTeO2Sxth232M2Sum_10->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxth232M2Sum_10", dAdaptiveArrayM2Sum);
  hnewTeO2Sxu238M2Sum_001 = hTeO2Sxu238M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxu238M2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxu238M2Sum_01 = hTeO2Sxu238M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxu238M2Sum_01", dAdaptiveArrayM2Sum);
  hnewTeO2Sxu238M2Sum_1 = hTeO2Sxu238M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxu238M2Sum_1", dAdaptiveArrayM2Sum);
  hnewTeO2Sxu238M2Sum_10 = hTeO2Sxu238M2Sum_10->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxu238M2Sum_10", dAdaptiveArrayM2Sum);

  hnewTeO2Sxu238M2Sum_100 = hTeO2Sxu238M2Sum_100->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxu238M2Sum_100", dAdaptiveArrayM2Sum);
  hnewTeO2Sxth232M2Sum_100 = hTeO2Sxth232M2Sum_100->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxth232M2Sum_100", dAdaptiveArrayM2Sum);
  hnewTeO2Sxpb210M2Sum_100 = hTeO2Sxpb210M2Sum_100->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxpb210M2Sum_100", dAdaptiveArrayM2Sum);

  hnewTeO2th232onlyM2Sum = hTeO2th232onlyM2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2th232onlyM2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2ra228pb208M2Sum = hTeO2ra228pb208M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2ra228pb208M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2th230onlyM2Sum = hTeO2th230onlyM2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2th230onlyM2Sum", dAdaptiveArrayM2Sum);

  hnewTeO2Sxth232onlyM2Sum_001 = hTeO2Sxth232onlyM2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxth232onlyM2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxra228pb208M2Sum_001 = hTeO2Sxra228pb208M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxra228pb208M2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxu238th230M2Sum_001 = hTeO2Sxu238th230M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxu238th230M2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxth230onlyM2Sum_001 = hTeO2Sxth230onlyM2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxth230onlyM2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxra226pb210M2Sum_001 = hTeO2Sxra226pb210M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxra226pb210M2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxpb210M2Sum_0001 = hTeO2Sxpb210M2Sum_0001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxpb210M2Sum_0001", dAdaptiveArrayM2Sum);

  hnewTeO2Sxth232onlyM2Sum_01 = hTeO2Sxth232onlyM2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxth232onlyM2Sum_01", dAdaptiveArrayM2Sum);
  hnewTeO2Sxra228pb208M2Sum_01 = hTeO2Sxra228pb208M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxra228pb208M2Sum_01", dAdaptiveArrayM2Sum);
  hnewTeO2Sxu238th230M2Sum_01 = hTeO2Sxu238th230M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxu238th230M2Sum_01", dAdaptiveArrayM2Sum);
  hnewTeO2Sxth230onlyM2Sum_01 = hTeO2Sxth230onlyM2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxth230onlyM2Sum_01", dAdaptiveArrayM2Sum);
  hnewTeO2Sxra226pb210M2Sum_01 = hTeO2Sxra226pb210M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxra226pb210M2Sum_01", dAdaptiveArrayM2Sum);

  hnewTeO2Sxth232onlyM2Sum_0001 = hTeO2Sxth232onlyM2Sum_0001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxth232onlyM2Sum_0001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxra228pb208M2Sum_0001 = hTeO2Sxra228pb208M2Sum_0001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxra228pb208M2Sum_0001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxu238th230M2Sum_0001 = hTeO2Sxu238th230M2Sum_0001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxu238th230M2Sum_0001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxth230onlyM2Sum_0001 = hTeO2Sxth230onlyM2Sum_0001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxth230onlyM2Sum_0001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxra226pb210M2Sum_0001 = hTeO2Sxra226pb210M2Sum_0001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxra226pb210M2Sum_0001", dAdaptiveArrayM2Sum);  

///////// CuBox + CuFrame M1 and M2
  hnewCuBox_CuFrameco60M1 = hCuBox_CuFrameco60M1->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameco60M1", dAdaptiveArrayM1);
  hnewCuBox_CuFramek40M1 = hCuBox_CuFramek40M1->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFramek40M1", dAdaptiveArrayM1);
  hnewCuBox_CuFrameth232M1 = hCuBox_CuFrameth232M1->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameth232M1", dAdaptiveArrayM1);
  hnewCuBox_CuFrameu238M1 = hCuBox_CuFrameu238M1->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameu238M1", dAdaptiveArrayM1);

  hnewCuBox_CuFrameth232M1_10 = hCuBox_CuFrameth232M1_10->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameth232M1_10", dAdaptiveArrayM1);
  hnewCuBox_CuFrameu238M1_10 = hCuBox_CuFrameu238M1_10->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameu238M1_10", dAdaptiveArrayM1);
  hnewCuBox_CuFramepb210M1_10 = hCuBox_CuFramepb210M1_10->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFramepb210M1_10", dAdaptiveArrayM1);
  hnewCuBox_CuFramepb210M1_1 = hCuBox_CuFramepb210M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFramepb210M1_1", dAdaptiveArrayM1);
  hnewCuBox_CuFramepb210M1_01 = hCuBox_CuFramepb210M1_01->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFramepb210M1_01", dAdaptiveArrayM1);
  hnewCuBox_CuFramepb210M1_001 = hCuBox_CuFramepb210M1_001->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFramepb210M1_001", dAdaptiveArrayM1);

  hnewCuBox_CuFrameth232M1_1 = hCuBox_CuFrameth232M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameth232M1_1", dAdaptiveArrayM1);
  hnewCuBox_CuFrameu238M1_1 = hCuBox_CuFrameu238M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameu238M1_1", dAdaptiveArrayM1);
  hnewCuBox_CuFrameth232M1_01 = hCuBox_CuFrameth232M1_01->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameth232M1_01", dAdaptiveArrayM1);
  hnewCuBox_CuFrameu238M1_01 = hCuBox_CuFrameu238M1_01->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameu238M1_01", dAdaptiveArrayM1);
  hnewCuBox_CuFrameth232M1_001 = hCuBox_CuFrameth232M1_001->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameth232M1_001", dAdaptiveArrayM1);
  hnewCuBox_CuFrameu238M1_001 = hCuBox_CuFrameu238M1_001->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameu238M1_001", dAdaptiveArrayM1);    

  hnewCuBox_CuFrameth232M1_100 = hCuBox_CuFrameth232M1_100->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameth232M1_100", dAdaptiveArrayM1);
  hnewCuBox_CuFrameu238M1_100 = hCuBox_CuFrameu238M1_100->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameu238M1_100", dAdaptiveArrayM1);
  hnewCuBox_CuFramepb210M1_100 = hCuBox_CuFramepb210M1_100->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFramepb210M1_100", dAdaptiveArrayM1);
  hnewCuBox_CuFrameth232M1_50 = hCuBox_CuFrameth232M1_50->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameth232M1_50", dAdaptiveArrayM1);
  hnewCuBox_CuFrameu238M1_50 = hCuBox_CuFrameu238M1_50->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameu238M1_50", dAdaptiveArrayM1);
  hnewCuBox_CuFramepb210M1_50 = hCuBox_CuFramepb210M1_50->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFramepb210M1_50", dAdaptiveArrayM1);
  hnewCuBox_CuFrameth232M1_5 = hCuBox_CuFrameth232M1_5->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameth232M1_5", dAdaptiveArrayM1);
  hnewCuBox_CuFrameu238M1_5 = hCuBox_CuFrameu238M1_5->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameu238M1_5", dAdaptiveArrayM1);
  hnewCuBox_CuFramepb210M1_5 = hCuBox_CuFramepb210M1_5->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFramepb210M1_5", dAdaptiveArrayM1);

  hnewCuBox_CuFrameco60M2 = hCuBox_CuFrameco60M2->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameco60M2", dAdaptiveArrayM2);
  hnewCuBox_CuFramek40M2 = hCuBox_CuFramek40M2->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFramek40M2", dAdaptiveArrayM2);
  hnewCuBox_CuFrameth232M2 = hCuBox_CuFrameth232M2->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameth232M2", dAdaptiveArrayM2);
  hnewCuBox_CuFrameu238M2 = hCuBox_CuFrameu238M2->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameu238M2", dAdaptiveArrayM2);

  hnewCuBox_CuFrameth232M2_10 = hCuBox_CuFrameth232M2_10->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameth232M2_10", dAdaptiveArrayM2);
  hnewCuBox_CuFrameu238M2_10 = hCuBox_CuFrameu238M2_10->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameu238M2_10", dAdaptiveArrayM2);
  hnewCuBox_CuFramepb210M2_10 = hCuBox_CuFramepb210M2_10->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFramepb210M2_10", dAdaptiveArrayM2);
  hnewCuBox_CuFramepb210M2_1 = hCuBox_CuFramepb210M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFramepb210M2_1", dAdaptiveArrayM2);
  hnewCuBox_CuFramepb210M2_01 = hCuBox_CuFramepb210M2_01->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFramepb210M2_01", dAdaptiveArrayM2);
  hnewCuBox_CuFramepb210M2_001 = hCuBox_CuFramepb210M2_001->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFramepb210M2_001", dAdaptiveArrayM2);

  hnewCuBox_CuFrameth232M2_1 = hCuBox_CuFrameth232M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameth232M2_1", dAdaptiveArrayM2);
  hnewCuBox_CuFrameu238M2_1 = hCuBox_CuFrameu238M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameu238M2_1", dAdaptiveArrayM2);
  hnewCuBox_CuFrameth232M2_01 = hCuBox_CuFrameth232M2_01->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameth232M2_01", dAdaptiveArrayM2);
  hnewCuBox_CuFrameu238M2_01 = hCuBox_CuFrameu238M2_01->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameu238M2_01", dAdaptiveArrayM2);
  hnewCuBox_CuFrameth232M2_001 = hCuBox_CuFrameth232M2_001->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameth232M2_001", dAdaptiveArrayM2);
  hnewCuBox_CuFrameu238M2_001 = hCuBox_CuFrameu238M2_001->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameu238M2_001", dAdaptiveArrayM2);  

  hnewCuBox_CuFrameth232M2_100 = hCuBox_CuFrameth232M2_100->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameth232M2_100", dAdaptiveArrayM2);
  hnewCuBox_CuFrameu238M2_100 = hCuBox_CuFrameu238M2_100->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameu238M2_100", dAdaptiveArrayM2);
  hnewCuBox_CuFramepb210M2_100 = hCuBox_CuFramepb210M2_100->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFramepb210M2_100", dAdaptiveArrayM2);
  hnewCuBox_CuFrameth232M2_50 = hCuBox_CuFrameth232M2_50->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameth232M2_50", dAdaptiveArrayM2);
  hnewCuBox_CuFrameu238M2_50 = hCuBox_CuFrameu238M2_50->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameu238M2_50", dAdaptiveArrayM2);
  hnewCuBox_CuFramepb210M2_50 = hCuBox_CuFramepb210M2_50->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFramepb210M2_50", dAdaptiveArrayM2);
  hnewCuBox_CuFrameth232M2_5 = hCuBox_CuFrameth232M2_5->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameth232M2_5", dAdaptiveArrayM2);
  hnewCuBox_CuFrameu238M2_5 = hCuBox_CuFrameu238M2_5->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFrameu238M2_5", dAdaptiveArrayM2);
  hnewCuBox_CuFramepb210M2_5 = hCuBox_CuFramepb210M2_5->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFramepb210M2_5", dAdaptiveArrayM2);

  hnewCuBox_CuFrameco60M2Sum = hCuBox_CuFrameco60M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameco60M2Sum", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFramek40M2Sum = hCuBox_CuFramek40M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFramek40M2Sum", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFrameth232M2Sum = hCuBox_CuFrameth232M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameth232M2Sum", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFrameu238M2Sum = hCuBox_CuFrameu238M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameu238M2Sum", dAdaptiveArrayM2Sum);

  hnewCuBox_CuFrameth232M2Sum_10 = hCuBox_CuFrameth232M2Sum_10->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameth232M2Sum_10", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFrameu238M2Sum_10 = hCuBox_CuFrameu238M2Sum_10->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameu238M2Sum_10", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFramepb210M2Sum_10 = hCuBox_CuFramepb210M2Sum_10->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFramepb210M2Sum_10", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFramepb210M2Sum_1 = hCuBox_CuFramepb210M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFramepb210M2Sum_1", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFramepb210M2Sum_01 = hCuBox_CuFramepb210M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFramepb210M2Sum_01", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFramepb210M2Sum_001 = hCuBox_CuFramepb210M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFramepb210M2Sum_001", dAdaptiveArrayM2Sum);

  hnewCuBox_CuFrameth232M2Sum_1 = hCuBox_CuFrameth232M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameth232M2Sum_1", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFrameu238M2Sum_1 = hCuBox_CuFrameu238M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameu238M2Sum_1", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFrameth232M2Sum_01 = hCuBox_CuFrameth232M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameth232M2Sum_01", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFrameu238M2Sum_01 = hCuBox_CuFrameu238M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameu238M2Sum_01", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFrameth232M2Sum_001 = hCuBox_CuFrameth232M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameth232M2Sum_001", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFrameu238M2Sum_001 = hCuBox_CuFrameu238M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameu238M2Sum_001", dAdaptiveArrayM2Sum);  

  hnewCuBox_CuFrameth232M2Sum_100 = hCuBox_CuFrameth232M2Sum_100->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameth232M2Sum_100", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFrameu238M2Sum_100 = hCuBox_CuFrameu238M2Sum_100->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameu238M2Sum_100", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFramepb210M2Sum_100 = hCuBox_CuFramepb210M2Sum_100->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFramepb210M2Sum_100", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFrameth232M2Sum_50 = hCuBox_CuFrameth232M2Sum_50->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameth232M2Sum_50", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFrameu238M2Sum_50 = hCuBox_CuFrameu238M2Sum_50->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameu238M2Sum_50", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFramepb210M2Sum_50 = hCuBox_CuFramepb210M2Sum_50->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFramepb210M2Sum_50", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFrameth232M2Sum_5 = hCuBox_CuFrameth232M2Sum_5->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameth232M2Sum_5", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFrameu238M2Sum_5 = hCuBox_CuFrameu238M2Sum_5->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFrameu238M2Sum_5", dAdaptiveArrayM2Sum);
  hnewCuBox_CuFramepb210M2Sum_5 = hCuBox_CuFramepb210M2Sum_5->Rebin(dAdaptiveBinsM2Sum, "hnewCuBox_CuFramepb210M2Sum_5", dAdaptiveArrayM2Sum);

////////// 50mK M1 and M2
  hnew50mKcs137M1 = h50mKcs137M1->Rebin(dAdaptiveBinsM1, "hnew50mKcs137M1", dAdaptiveArrayM1);
  hnew50mKcs137M2 = h50mKcs137M2->Rebin(dAdaptiveBinsM2, "hnew50mKcs137M2", dAdaptiveArrayM2);
  hnew50mKcs137M2Sum = h50mKcs137M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnew50mKcs137M2Sum", dAdaptiveArrayM2Sum);
  
///////// Roman Lead M1 and M2
  hnewPbRomco60M1 = hPbRomco60M1->Rebin(dAdaptiveBinsM1, "hnewPbRomco60M1", dAdaptiveArrayM1);
  hnewPbRomcs137M1 = hPbRomcs137M1->Rebin(dAdaptiveBinsM1, "hnewPbRomcs137M1", dAdaptiveArrayM1);
  hnewPbRomk40M1 = hPbRomk40M1->Rebin(dAdaptiveBinsM1, "hnewPbRomk40M1", dAdaptiveArrayM1);
  hnewPbRompb210M1 = hPbRompb210M1->Rebin(dAdaptiveBinsM1, "hnewPbRompb210M1", dAdaptiveArrayM1);
  hnewPbRomth232M1 = hPbRomth232M1->Rebin(dAdaptiveBinsM1, "hnewPbRomth232M1", dAdaptiveArrayM1);
  hnewPbRomu238M1 = hPbRomu238M1->Rebin(dAdaptiveBinsM1, "hnewPbRomu238M1", dAdaptiveArrayM1);

  hnewPbRomco60M2 = hPbRomco60M2->Rebin(dAdaptiveBinsM2, "hnewPbRomco60M2", dAdaptiveArrayM2);
  hnewPbRomcs137M2 = hPbRomcs137M2->Rebin(dAdaptiveBinsM2, "hnewPbRomcs137M2", dAdaptiveArrayM2);
  hnewPbRomk40M2 = hPbRomk40M2->Rebin(dAdaptiveBinsM2, "hnewPbRomk40M2", dAdaptiveArrayM2);
  hnewPbRompb210M2 = hPbRompb210M2->Rebin(dAdaptiveBinsM2, "hnewPbRompb210M2", dAdaptiveArrayM2);
  hnewPbRomth232M2 = hPbRomth232M2->Rebin(dAdaptiveBinsM2, "hnewPbRomth232M2", dAdaptiveArrayM2);
  hnewPbRomu238M2 = hPbRomu238M2->Rebin(dAdaptiveBinsM2, "hnewPbRomu238M2", dAdaptiveArrayM2);

  hnewPbRomco60M2Sum = hPbRomco60M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewPbRomco60M2Sum", dAdaptiveArrayM2Sum);
  hnewPbRomcs137M2Sum = hPbRomcs137M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewPbRomcs137M2Sum", dAdaptiveArrayM2Sum);
  hnewPbRomk40M2Sum = hPbRomk40M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewPbRomk40M2Sum", dAdaptiveArrayM2Sum);
  hnewPbRompb210M2Sum = hPbRompb210M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewPbRompb210M2Sum", dAdaptiveArrayM2Sum);
  hnewPbRomth232M2Sum = hPbRomth232M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewPbRomth232M2Sum", dAdaptiveArrayM2Sum);
  hnewPbRomu238M2Sum = hPbRomu238M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewPbRomu238M2Sum", dAdaptiveArrayM2Sum);

////////// Internal Shields
  hnewInternalco60M1 = hInternalco60M1->Rebin(dAdaptiveBinsM1, "hnewInternalco60M1", dAdaptiveArrayM1);
  hnewInternalk40M1 = hInternalk40M1->Rebin(dAdaptiveBinsM1, "hnewInternalk40M1", dAdaptiveArrayM1);
  hnewInternalth232M1 = hInternalth232M1->Rebin(dAdaptiveBinsM1, "hnewInternalth232M1", dAdaptiveArrayM1);
  hnewInternalu238M1 = hInternalu238M1->Rebin(dAdaptiveBinsM1, "hnewInternalu238M1", dAdaptiveArrayM1);

  hnewInternalco60M2 = hInternalco60M2->Rebin(dAdaptiveBinsM2, "hnewInternalco60M2", dAdaptiveArrayM2);
  hnewInternalk40M2 = hInternalk40M2->Rebin(dAdaptiveBinsM2, "hnewInternalk40M2", dAdaptiveArrayM2);
  hnewInternalth232M2 = hInternalth232M2->Rebin(dAdaptiveBinsM2, "hnewInternalth232M2", dAdaptiveArrayM2);
  hnewInternalu238M2 = hInternalu238M2->Rebin(dAdaptiveBinsM2, "hnewInternalu238M2", dAdaptiveArrayM2);

  hnewInternalco60M2Sum = hInternalco60M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewInternalco60M2Sum", dAdaptiveArrayM2Sum);
  hnewInternalk40M2Sum = hInternalk40M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewInternalk40M2Sum", dAdaptiveArrayM2Sum);
  hnewInternalth232M2Sum = hInternalth232M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewInternalth232M2Sum", dAdaptiveArrayM2Sum);
  hnewInternalu238M2Sum = hInternalu238M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewInternalu238M2Sum", dAdaptiveArrayM2Sum);

////////// OVC M1 and M2
  hnewOVCco60M1 = hOVCco60M1->Rebin(dAdaptiveBinsM1, "hnewOVCco60M1", dAdaptiveArrayM1);
  hnewOVCk40M1 = hOVCk40M1->Rebin(dAdaptiveBinsM1, "hnewOVCk40M1", dAdaptiveArrayM1);
  hnewOVCth232M1 = hOVCth232M1->Rebin(dAdaptiveBinsM1, "hnewOVCth232M1", dAdaptiveArrayM1);
  hnewOVCu238M1 = hOVCu238M1->Rebin(dAdaptiveBinsM1, "hnewOVCu238M1", dAdaptiveArrayM1);

  hnewOVCco60M2 = hOVCco60M2->Rebin(dAdaptiveBinsM2, "hnewOVCco60M2", dAdaptiveArrayM2);
  hnewOVCk40M2 = hOVCk40M2->Rebin(dAdaptiveBinsM2, "hnewOVCk40M2", dAdaptiveArrayM2);
  hnewOVCth232M2 = hOVCth232M2->Rebin(dAdaptiveBinsM2, "hnewOVCth232M2", dAdaptiveArrayM2);
  hnewOVCu238M2 = hOVCu238M2->Rebin(dAdaptiveBinsM2, "hnewOVCu238M2", dAdaptiveArrayM2);

  hnewOVCco60M2Sum = hOVCco60M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewOVCco60M2Sum", dAdaptiveArrayM2Sum);
  hnewOVCk40M2Sum = hOVCk40M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewOVCk40M2Sum", dAdaptiveArrayM2Sum);
  hnewOVCth232M2Sum = hOVCth232M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewOVCth232M2Sum", dAdaptiveArrayM2Sum);
  hnewOVCu238M2Sum = hOVCu238M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewOVCu238M2Sum", dAdaptiveArrayM2Sum);

  hnewExtPbbi210M1 = hExtPbbi210M1->Rebin(dAdaptiveBinsM1, "hnewExtPbbi210M1", dAdaptiveArrayM1);
  hnewExtPbbi210M2 = hExtPbbi210M2->Rebin(dAdaptiveBinsM2, "hnewExtPbbi210M2", dAdaptiveArrayM2);
  hnewExtPbbi210M2Sum = hExtPbbi210M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewExtPbbi210M2Sum", dAdaptiveArrayM2Sum);


  ///////// Fudge Factors
  hnewOVC804M1 = hOVC804M1->Rebin(dAdaptiveBinsM1, "hnewOVC804M1", dAdaptiveArrayM1);
  hnewOVC835M1 = hOVC835M1->Rebin(dAdaptiveBinsM1, "hnewOVC835M1", dAdaptiveArrayM1);
  hnewOVC1063M1 = hOVC1063M1->Rebin(dAdaptiveBinsM1, "hnewOVC1063M1", dAdaptiveArrayM1);

  hnewOVC804M2 = hOVC804M2->Rebin(dAdaptiveBinsM2, "hnewOVC804M2", dAdaptiveArrayM2);
  hnewOVC835M2 = hOVC835M2->Rebin(dAdaptiveBinsM2, "hnewOVC835M2", dAdaptiveArrayM2);
  hnewOVC1063M2 = hOVC1063M2->Rebin(dAdaptiveBinsM2, "hnewOVC1063M2", dAdaptiveArrayM2);

  hnewOVC804M2Sum = hOVC804M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewOVC804M2Sum", dAdaptiveArrayM2Sum);
  hnewOVC835M2Sum = hOVC835M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewOVC835M2Sum", dAdaptiveArrayM2Sum);
  hnewOVC1063M2Sum = hOVC1063M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewOVC1063M2Sum", dAdaptiveArrayM2Sum);

  // Fill adaptive binning histograms
  for(int i = 1; i <= dAdaptiveBinsM1; i++)
  {
    hAdapTeO20nuM1->SetBinContent(i, hnewTeO20nuM1->GetBinContent(i)/hnewTeO20nuM1->GetBinWidth(i));
    hAdapTeO22nuM1->SetBinContent(i, hnewTeO22nuM1->GetBinContent(i)/hnewTeO22nuM1->GetBinWidth(i));
    hAdapTeO2co60M1->SetBinContent(i, hnewTeO2co60M1->GetBinContent(i)/hnewTeO2co60M1->GetBinWidth(i));
    hAdapTeO2k40M1->SetBinContent(i, hnewTeO2k40M1->GetBinContent(i)/hnewTeO2k40M1->GetBinWidth(i));
    hAdapTeO2pb210M1->SetBinContent(i, hnewTeO2pb210M1->GetBinContent(i)/hnewTeO2pb210M1->GetBinWidth(i));
    hAdapTeO2po210M1->SetBinContent(i, hnewTeO2po210M1->GetBinContent(i)/hnewTeO2po210M1->GetBinWidth(i));
    hAdapTeO2te125M1->SetBinContent(i, hnewTeO2te125M1->GetBinContent(i)/hnewTeO2te125M1->GetBinWidth(i));
    hAdapTeO2th232M1->SetBinContent(i, hnewTeO2th232M1->GetBinContent(i)/hnewTeO2th232M1->GetBinWidth(i));
    // hAdapTeO2th228M1->SetBinContent(i, hnewTeO2th228M1->GetBinContent(i)/hnewTeO2th228M1->GetBinWidth(i));
    // hAdapTeO2ra226M1->SetBinContent(i, hnewTeO2ra226M1->GetBinContent(i)/hnewTeO2ra226M1->GetBinWidth(i));
    // hAdapTeO2rn222M1->SetBinContent(i, hnewTeO2rn222M1->GetBinContent(i)/hnewTeO2rn222M1->GetBinWidth(i));
    hAdapTeO2u238M1->SetBinContent(i, hnewTeO2u238M1->GetBinContent(i)/hnewTeO2u238M1->GetBinWidth(i));
    // hAdapTeO2th230M1->SetBinContent(i, hnewTeO2th230M1->GetBinContent(i)/hnewTeO2th230M1->GetBinWidth(i));
    // hAdapTeO2u234M1->SetBinContent(i, hnewTeO2u234M1->GetBinContent(i)/hnewTeO2u234M1->GetBinWidth(i));
    hAdapTeO2sb125M1->SetBinContent(i, hnewTeO2sb125M1->GetBinContent(i)/hnewTeO2sb125M1->GetBinWidth(i));

    hAdapTeO2th232onlyM1->SetBinContent(i, hnewTeO2th232onlyM1->GetBinContent(i)/hnewTeO2th232onlyM1->GetBinWidth(i));
    hAdapTeO2ra228pb208M1->SetBinContent(i, hnewTeO2ra228pb208M1->GetBinContent(i)/hnewTeO2ra228pb208M1->GetBinWidth(i));
    hAdapTeO2th230onlyM1->SetBinContent(i, hnewTeO2th230onlyM1->GetBinContent(i)/hnewTeO2th230onlyM1->GetBinWidth(i));
    hAdapTeO2u238th230M1->SetBinContent(i, hnewTeO2u238th230M1->GetBinContent(i)/hnewTeO2u238th230M1->GetBinWidth(i));
    hAdapTeO2ra226pb210M1->SetBinContent(i, hnewTeO2ra226pb210M1->GetBinContent(i)/hnewTeO2ra226pb210M1->GetBinWidth(i));

    hAdapTeO2Sxth232onlyM1_001->SetBinContent(i, hnewTeO2Sxth232onlyM1_001->GetBinContent(i)/hnewTeO2Sxth232onlyM1_001->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M1_001->SetBinContent(i, hnewTeO2Sxra228pb208M1_001->GetBinContent(i)/hnewTeO2Sxra228pb208M1_001->GetBinWidth(i));
    hAdapTeO2Sxu238th230M1_001->SetBinContent(i, hnewTeO2Sxu238th230M1_001->GetBinContent(i)/hnewTeO2Sxu238th230M1_001->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM1_001->SetBinContent(i, hnewTeO2Sxth230onlyM1_001->GetBinContent(i)/hnewTeO2Sxth230onlyM1_001->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M1_001->SetBinContent(i, hnewTeO2Sxra226pb210M1_001->GetBinContent(i)/hnewTeO2Sxra226pb210M1_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M1_0001->SetBinContent(i, hnewTeO2Sxpb210M1_0001->GetBinContent(i)/hnewTeO2Sxpb210M1_0001->GetBinWidth(i));

    hAdapTeO2Sxth232onlyM1_01->SetBinContent(i, hnewTeO2Sxth232onlyM1_01->GetBinContent(i)/hnewTeO2Sxth232onlyM1_01->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M1_01->SetBinContent(i, hnewTeO2Sxra228pb208M1_01->GetBinContent(i)/hnewTeO2Sxra228pb208M1_01->GetBinWidth(i));
    hAdapTeO2Sxu238th230M1_01->SetBinContent(i, hnewTeO2Sxu238th230M1_01->GetBinContent(i)/hnewTeO2Sxu238th230M1_01->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM1_01->SetBinContent(i, hnewTeO2Sxth230onlyM1_01->GetBinContent(i)/hnewTeO2Sxth230onlyM1_01->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M1_01->SetBinContent(i, hnewTeO2Sxra226pb210M1_01->GetBinContent(i)/hnewTeO2Sxra226pb210M1_01->GetBinWidth(i));

    hAdapTeO2Sxth232onlyM1_0001->SetBinContent(i, hnewTeO2Sxth232onlyM1_0001->GetBinContent(i)/hnewTeO2Sxth232onlyM1_0001->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M1_0001->SetBinContent(i, hnewTeO2Sxra228pb208M1_0001->GetBinContent(i)/hnewTeO2Sxra228pb208M1_0001->GetBinWidth(i));
    hAdapTeO2Sxu238th230M1_0001->SetBinContent(i, hnewTeO2Sxu238th230M1_0001->GetBinContent(i)/hnewTeO2Sxu238th230M1_0001->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM1_0001->SetBinContent(i, hnewTeO2Sxth230onlyM1_0001->GetBinContent(i)/hnewTeO2Sxth230onlyM1_0001->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M1_0001->SetBinContent(i, hnewTeO2Sxra226pb210M1_0001->GetBinContent(i)/hnewTeO2Sxra226pb210M1_0001->GetBinWidth(i));

    hAdapTeO2Sxpb210M1_001->SetBinContent(i, hnewTeO2Sxpb210M1_001->GetBinContent(i)/hnewTeO2Sxpb210M1_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M1_01->SetBinContent(i, hnewTeO2Sxpb210M1_01->GetBinContent(i)/hnewTeO2Sxpb210M1_01->GetBinWidth(i));
    hAdapTeO2Sxpb210M1_1->SetBinContent(i, hnewTeO2Sxpb210M1_1->GetBinContent(i)/hnewTeO2Sxpb210M1_1->GetBinWidth(i));
    hAdapTeO2Sxpb210M1_10->SetBinContent(i, hnewTeO2Sxpb210M1_10->GetBinContent(i)/hnewTeO2Sxpb210M1_10->GetBinWidth(i));
    hAdapTeO2Sxpo210M1_001->SetBinContent(i, hnewTeO2Sxpo210M1_001->GetBinContent(i)/hnewTeO2Sxpo210M1_001->GetBinWidth(i));
    hAdapTeO2Sxpo210M1_01->SetBinContent(i, hnewTeO2Sxpo210M1_01->GetBinContent(i)/hnewTeO2Sxpo210M1_01->GetBinWidth(i));
    hAdapTeO2Sxpo210M1_1->SetBinContent(i, hnewTeO2Sxpo210M1_1->GetBinContent(i)/hnewTeO2Sxpo210M1_1->GetBinWidth(i));
    hAdapTeO2Sxth232M1_001->SetBinContent(i, hnewTeO2Sxth232M1_001->GetBinContent(i)/hnewTeO2Sxth232M1_001->GetBinWidth(i));
    hAdapTeO2Sxth232M1_01->SetBinContent(i, hnewTeO2Sxth232M1_01->GetBinContent(i)/hnewTeO2Sxth232M1_01->GetBinWidth(i));
    hAdapTeO2Sxth232M1_1->SetBinContent(i, hnewTeO2Sxth232M1_1->GetBinContent(i)/hnewTeO2Sxth232M1_1->GetBinWidth(i));
    hAdapTeO2Sxth232M1_10->SetBinContent(i, hnewTeO2Sxth232M1_10->GetBinContent(i)/hnewTeO2Sxth232M1_10->GetBinWidth(i));
    hAdapTeO2Sxu238M1_001->SetBinContent(i, hnewTeO2Sxu238M1_001->GetBinContent(i)/hnewTeO2Sxu238M1_001->GetBinWidth(i));
    hAdapTeO2Sxu238M1_01->SetBinContent(i, hnewTeO2Sxu238M1_01->GetBinContent(i)/hnewTeO2Sxu238M1_01->GetBinWidth(i));
    hAdapTeO2Sxu238M1_1->SetBinContent(i, hnewTeO2Sxu238M1_1->GetBinContent(i)/hnewTeO2Sxu238M1_1->GetBinWidth(i));
    hAdapTeO2Sxu238M1_10->SetBinContent(i, hnewTeO2Sxu238M1_10->GetBinContent(i)/hnewTeO2Sxu238M1_10->GetBinWidth(i));

    hAdapTeO2Sxu238M1_100->SetBinContent(i, hnewTeO2Sxu238M1_100->GetBinContent(i)/hnewTeO2Sxu238M1_100->GetBinWidth(i));
    hAdapTeO2Sxth232M1_100->SetBinContent(i, hnewTeO2Sxth232M1_100->GetBinContent(i)/hnewTeO2Sxth232M1_100->GetBinWidth(i));
    hAdapTeO2Sxpb210M1_100->SetBinContent(i, hnewTeO2Sxpb210M1_100->GetBinContent(i)/hnewTeO2Sxpb210M1_100->GetBinWidth(i));

    hAdapCuBox_CuFrameco60M1->SetBinContent(i, hnewCuBox_CuFrameco60M1->GetBinContent(i)/hnewCuBox_CuFrameco60M1->GetBinWidth(i));
    hAdapCuBox_CuFramek40M1->SetBinContent(i, hnewCuBox_CuFramek40M1->GetBinContent(i)/hnewCuBox_CuFramek40M1->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M1->SetBinContent(i, hnewCuBox_CuFrameth232M1->GetBinContent(i)/hnewCuBox_CuFrameth232M1->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M1->SetBinContent(i, hnewCuBox_CuFrameu238M1->GetBinContent(i)/hnewCuBox_CuFrameu238M1->GetBinWidth(i));

    hAdapCuBox_CuFrameth232M1_10->SetBinContent(i, hnewCuBox_CuFrameth232M1_10->GetBinContent(i)/hnewCuBox_CuFrameth232M1_10->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M1_10->SetBinContent(i, hnewCuBox_CuFrameu238M1_10->GetBinContent(i)/hnewCuBox_CuFrameu238M1_10->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M1_10->SetBinContent(i, hnewCuBox_CuFramepb210M1_10->GetBinContent(i)/hnewCuBox_CuFramepb210M1_10->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M1_1->SetBinContent(i, hnewCuBox_CuFramepb210M1_1->GetBinContent(i)/hnewCuBox_CuFramepb210M1_1->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M1_01->SetBinContent(i, hnewCuBox_CuFramepb210M1_01->GetBinContent(i)/hnewCuBox_CuFramepb210M1_01->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M1_001->SetBinContent(i, hnewCuBox_CuFramepb210M1_001->GetBinContent(i)/hnewCuBox_CuFramepb210M1_001->GetBinWidth(i));

    hAdapCuBox_CuFrameth232M1_1->SetBinContent(i, hnewCuBox_CuFrameth232M1_1->GetBinContent(i)/hnewCuBox_CuFrameth232M1_1->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M1_1->SetBinContent(i, hnewCuBox_CuFrameu238M1_1->GetBinContent(i)/hnewCuBox_CuFrameu238M1_1->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M1_01->SetBinContent(i, hnewCuBox_CuFrameth232M1_01->GetBinContent(i)/hnewCuBox_CuFrameth232M1_01->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M1_01->SetBinContent(i, hnewCuBox_CuFrameu238M1_01->GetBinContent(i)/hnewCuBox_CuFrameu238M1_01->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M1_001->SetBinContent(i, hnewCuBox_CuFrameth232M1_001->GetBinContent(i)/hnewCuBox_CuFrameth232M1_001->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M1_001->SetBinContent(i, hnewCuBox_CuFrameu238M1_001->GetBinContent(i)/hnewCuBox_CuFrameu238M1_001->GetBinWidth(i));        

    hAdapCuBox_CuFrameth232M1_100->SetBinContent(i, hnewCuBox_CuFrameth232M1_100->GetBinContent(i)/hnewCuBox_CuFrameth232M1_100->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M1_100->SetBinContent(i, hnewCuBox_CuFrameu238M1_100->GetBinContent(i)/hnewCuBox_CuFrameu238M1_100->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M1_100->SetBinContent(i, hnewCuBox_CuFramepb210M1_100->GetBinContent(i)/hnewCuBox_CuFramepb210M1_100->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M1_50->SetBinContent(i, hnewCuBox_CuFrameth232M1_50->GetBinContent(i)/hnewCuBox_CuFrameth232M1_50->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M1_50->SetBinContent(i, hnewCuBox_CuFrameu238M1_50->GetBinContent(i)/hnewCuBox_CuFrameu238M1_50->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M1_50->SetBinContent(i, hnewCuBox_CuFramepb210M1_50->GetBinContent(i)/hnewCuBox_CuFramepb210M1_50->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M1_5->SetBinContent(i, hnewCuBox_CuFrameth232M1_5->GetBinContent(i)/hnewCuBox_CuFrameth232M1_5->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M1_5->SetBinContent(i, hnewCuBox_CuFrameu238M1_5->GetBinContent(i)/hnewCuBox_CuFrameu238M1_5->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M1_5->SetBinContent(i, hnewCuBox_CuFramepb210M1_5->GetBinContent(i)/hnewCuBox_CuFramepb210M1_5->GetBinWidth(i));

    hAdap50mKcs137M1->SetBinContent(i, hnew50mKcs137M1->GetBinContent(i)/hnew50mKcs137M1->GetBinWidth(i));

    hAdapPbRomco60M1->SetBinContent(i, hnewPbRomco60M1->GetBinContent(i)/hnewPbRomco60M1->GetBinWidth(i));
    hAdapPbRomcs137M1->SetBinContent(i, hnewPbRomcs137M1->GetBinContent(i)/hnewPbRomcs137M1->GetBinWidth(i));
    hAdapPbRomk40M1->SetBinContent(i, hnewPbRomk40M1->GetBinContent(i)/hnewPbRomk40M1->GetBinWidth(i));
    hAdapPbRompb210M1->SetBinContent(i, hnewPbRompb210M1->GetBinContent(i)/hnewPbRompb210M1->GetBinWidth(i));
    hAdapPbRomth232M1->SetBinContent(i, hnewPbRomth232M1->GetBinContent(i)/hnewPbRomth232M1->GetBinWidth(i));
    hAdapPbRomu238M1->SetBinContent(i, hnewPbRomu238M1->GetBinContent(i)/hnewPbRomu238M1->GetBinWidth(i));

    hAdapInternalco60M1->SetBinContent(i, hnewInternalco60M1->GetBinContent(i)/hnewInternalco60M1->GetBinWidth(i));
    hAdapInternalk40M1->SetBinContent(i, hnewInternalk40M1->GetBinContent(i)/hnewInternalk40M1->GetBinWidth(i));
    hAdapInternalth232M1->SetBinContent(i, hnewInternalth232M1->GetBinContent(i)/hnewInternalth232M1->GetBinWidth(i));
    hAdapInternalu238M1->SetBinContent(i, hnewInternalu238M1->GetBinContent(i)/hnewInternalu238M1->GetBinWidth(i));

    hAdapOVCco60M1->SetBinContent(i, hnewOVCco60M1->GetBinContent(i)/hnewOVCco60M1->GetBinWidth(i));
    hAdapOVCk40M1->SetBinContent(i, hnewOVCk40M1->GetBinContent(i)/hnewOVCk40M1->GetBinWidth(i));
    hAdapOVCth232M1->SetBinContent(i, hnewOVCth232M1->GetBinContent(i)/hnewOVCth232M1->GetBinWidth(i));
    hAdapOVCu238M1->SetBinContent(i, hnewOVCu238M1->GetBinContent(i)/hnewOVCu238M1->GetBinWidth(i));

    hAdapExtPbbi210M1->SetBinContent(i, hnewExtPbbi210M1->GetBinContent(i)/hnewExtPbbi210M1->GetBinWidth(i));

    // Fudge
    hAdapOVC804M1->SetBinContent(i, hnewOVC804M1->GetBinContent(i)/hnewOVC804M1->GetBinWidth(i));
    hAdapOVC835M1->SetBinContent(i, hnewOVC835M1->GetBinContent(i)/hnewOVC835M1->GetBinWidth(i));
    hAdapOVC1063M1->SetBinContent(i, hnewOVC1063M1->GetBinContent(i)/hnewOVC1063M1->GetBinWidth(i));
  }

  for(int i = 1; i <= dAdaptiveBinsM2; i++)
  {
    hAdapTeO20nuM2->SetBinContent(i, hnewTeO20nuM2->GetBinContent(i)/hnewTeO20nuM2->GetBinWidth(i));
    hAdapTeO22nuM2->SetBinContent(i, hnewTeO22nuM2->GetBinContent(i)/hnewTeO22nuM2->GetBinWidth(i));
    hAdapTeO2co60M2->SetBinContent(i, hnewTeO2co60M2->GetBinContent(i)/hnewTeO2co60M2->GetBinWidth(i));
    hAdapTeO2k40M2->SetBinContent(i, hnewTeO2k40M2->GetBinContent(i)/hnewTeO2k40M2->GetBinWidth(i));
    hAdapTeO2pb210M2->SetBinContent(i, hnewTeO2pb210M2->GetBinContent(i)/hnewTeO2pb210M2->GetBinWidth(i));
    hAdapTeO2po210M2->SetBinContent(i, hnewTeO2po210M2->GetBinContent(i)/hnewTeO2po210M2->GetBinWidth(i));
    hAdapTeO2te125M2->SetBinContent(i, hnewTeO2te125M2->GetBinContent(i)/hnewTeO2te125M2->GetBinWidth(i));
    hAdapTeO2th232M2->SetBinContent(i, hnewTeO2th232M2->GetBinContent(i)/hnewTeO2th232M2->GetBinWidth(i));
    // hAdapTeO2th228M2->SetBinContent(i, hnewTeO2th228M2->GetBinContent(i)/hnewTeO2th228M2->GetBinWidth(i));
    // hAdapTeO2ra226M2->SetBinContent(i, hnewTeO2ra226M2->GetBinContent(i)/hnewTeO2ra226M2->GetBinWidth(i));
    // hAdapTeO2rn222M2->SetBinContent(i, hnewTeO2rn222M2->GetBinContent(i)/hnewTeO2rn222M2->GetBinWidth(i));
    hAdapTeO2u238M2->SetBinContent(i, hnewTeO2u238M2->GetBinContent(i)/hnewTeO2u238M2->GetBinWidth(i));
    // hAdapTeO2th230M2->SetBinContent(i, hnewTeO2th230M2->GetBinContent(i)/hnewTeO2th230M2->GetBinWidth(i));
    // hAdapTeO2u234M2->SetBinContent(i, hnewTeO2u234M2->GetBinContent(i)/hnewTeO2u234M2->GetBinWidth(i));
    hAdapTeO2sb125M2->SetBinContent(i, hnewTeO2sb125M2->GetBinContent(i)/hnewTeO2sb125M2->GetBinWidth(i));

    hAdapTeO2Sxpb210M2_001->SetBinContent(i, hnewTeO2Sxpb210M2_001->GetBinContent(i)/hnewTeO2Sxpb210M2_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_01->SetBinContent(i, hnewTeO2Sxpb210M2_01->GetBinContent(i)/hnewTeO2Sxpb210M2_01->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_1->SetBinContent(i, hnewTeO2Sxpb210M2_1->GetBinContent(i)/hnewTeO2Sxpb210M2_1->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_10->SetBinContent(i, hnewTeO2Sxpb210M2_10->GetBinContent(i)/hnewTeO2Sxpb210M2_10->GetBinWidth(i));
    hAdapTeO2Sxpo210M2_001->SetBinContent(i, hnewTeO2Sxpo210M2_001->GetBinContent(i)/hnewTeO2Sxpo210M2_001->GetBinWidth(i));
    hAdapTeO2Sxpo210M2_01->SetBinContent(i, hnewTeO2Sxpo210M2_01->GetBinContent(i)/hnewTeO2Sxpo210M2_01->GetBinWidth(i));
    hAdapTeO2Sxpo210M2_1->SetBinContent(i, hnewTeO2Sxpo210M2_1->GetBinContent(i)/hnewTeO2Sxpo210M2_1->GetBinWidth(i));
    hAdapTeO2Sxth232M2_001->SetBinContent(i, hnewTeO2Sxth232M2_001->GetBinContent(i)/hnewTeO2Sxth232M2_001->GetBinWidth(i));
    hAdapTeO2Sxth232M2_01->SetBinContent(i, hnewTeO2Sxth232M2_01->GetBinContent(i)/hnewTeO2Sxth232M2_01->GetBinWidth(i));
    hAdapTeO2Sxth232M2_1->SetBinContent(i, hnewTeO2Sxth232M2_1->GetBinContent(i)/hnewTeO2Sxth232M2_1->GetBinWidth(i));
    hAdapTeO2Sxth232M2_10->SetBinContent(i, hnewTeO2Sxth232M2_10->GetBinContent(i)/hnewTeO2Sxth232M2_10->GetBinWidth(i));
    hAdapTeO2Sxu238M2_001->SetBinContent(i, hnewTeO2Sxu238M2_001->GetBinContent(i)/hnewTeO2Sxu238M2_001->GetBinWidth(i));
    hAdapTeO2Sxu238M2_01->SetBinContent(i, hnewTeO2Sxu238M2_01->GetBinContent(i)/hnewTeO2Sxu238M2_01->GetBinWidth(i));
    hAdapTeO2Sxu238M2_1->SetBinContent(i, hnewTeO2Sxu238M2_1->GetBinContent(i)/hnewTeO2Sxu238M2_1->GetBinWidth(i));
    hAdapTeO2Sxu238M2_10->SetBinContent(i, hnewTeO2Sxu238M2_10->GetBinContent(i)/hnewTeO2Sxu238M2_10->GetBinWidth(i));

    hAdapTeO2th232onlyM2->SetBinContent(i, hnewTeO2th232onlyM2->GetBinContent(i)/hnewTeO2th232onlyM2->GetBinWidth(i));
    hAdapTeO2ra228pb208M2->SetBinContent(i, hnewTeO2ra228pb208M2->GetBinContent(i)/hnewTeO2ra228pb208M2->GetBinWidth(i));
    hAdapTeO2th230onlyM2->SetBinContent(i, hnewTeO2th230onlyM2->GetBinContent(i)/hnewTeO2th230onlyM2->GetBinWidth(i));
    hAdapTeO2u238th230M2->SetBinContent(i, hnewTeO2u238th230M2->GetBinContent(i)/hnewTeO2u238th230M2->GetBinWidth(i));
    hAdapTeO2ra226pb210M2->SetBinContent(i, hnewTeO2ra226pb210M2->GetBinContent(i)/hnewTeO2ra226pb210M2->GetBinWidth(i));

    hAdapTeO2Sxu238M2_100->SetBinContent(i, hnewTeO2Sxu238M2_100->GetBinContent(i)/hnewTeO2Sxu238M2_100->GetBinWidth(i));
    hAdapTeO2Sxth232M2_100->SetBinContent(i, hnewTeO2Sxth232M2_100->GetBinContent(i)/hnewTeO2Sxth232M2_100->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_100->SetBinContent(i, hnewTeO2Sxpb210M2_100->GetBinContent(i)/hnewTeO2Sxpb210M2_100->GetBinWidth(i));

    hAdapTeO2Sxth232onlyM2_001->SetBinContent(i, hnewTeO2Sxth232onlyM2_001->GetBinContent(i)/hnewTeO2Sxth232onlyM2_001->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M2_001->SetBinContent(i, hnewTeO2Sxra228pb208M2_001->GetBinContent(i)/hnewTeO2Sxra228pb208M2_001->GetBinWidth(i));
    hAdapTeO2Sxu238th230M2_001->SetBinContent(i, hnewTeO2Sxu238th230M2_001->GetBinContent(i)/hnewTeO2Sxu238th230M2_001->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM2_001->SetBinContent(i, hnewTeO2Sxth230onlyM2_001->GetBinContent(i)/hnewTeO2Sxth230onlyM2_001->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M2_001->SetBinContent(i, hnewTeO2Sxra226pb210M2_001->GetBinContent(i)/hnewTeO2Sxra226pb210M2_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_0001->SetBinContent(i, hnewTeO2Sxpb210M2_0001->GetBinContent(i)/hnewTeO2Sxpb210M2_0001->GetBinWidth(i));

    hAdapTeO2Sxth232onlyM2_01->SetBinContent(i, hnewTeO2Sxth232onlyM2_01->GetBinContent(i)/hnewTeO2Sxth232onlyM2_01->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M2_01->SetBinContent(i, hnewTeO2Sxra228pb208M2_01->GetBinContent(i)/hnewTeO2Sxra228pb208M2_01->GetBinWidth(i));
    hAdapTeO2Sxu238th230M2_01->SetBinContent(i, hnewTeO2Sxu238th230M2_01->GetBinContent(i)/hnewTeO2Sxu238th230M2_01->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM2_01->SetBinContent(i, hnewTeO2Sxth230onlyM2_01->GetBinContent(i)/hnewTeO2Sxth230onlyM2_01->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M2_01->SetBinContent(i, hnewTeO2Sxra226pb210M2_01->GetBinContent(i)/hnewTeO2Sxra226pb210M2_01->GetBinWidth(i));

    hAdapTeO2Sxth232onlyM2_0001->SetBinContent(i, hnewTeO2Sxth232onlyM2_0001->GetBinContent(i)/hnewTeO2Sxth232onlyM2_0001->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M2_0001->SetBinContent(i, hnewTeO2Sxra228pb208M2_0001->GetBinContent(i)/hnewTeO2Sxra228pb208M2_0001->GetBinWidth(i));
    hAdapTeO2Sxu238th230M2_0001->SetBinContent(i, hnewTeO2Sxu238th230M2_0001->GetBinContent(i)/hnewTeO2Sxu238th230M2_0001->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM2_0001->SetBinContent(i, hnewTeO2Sxth230onlyM2_0001->GetBinContent(i)/hnewTeO2Sxth230onlyM2_0001->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M2_0001->SetBinContent(i, hnewTeO2Sxra226pb210M2_0001->GetBinContent(i)/hnewTeO2Sxra226pb210M2_0001->GetBinWidth(i));

    hAdapCuBox_CuFrameco60M2->SetBinContent(i, hnewCuBox_CuFrameco60M2->GetBinContent(i)/hnewCuBox_CuFrameco60M2->GetBinWidth(i));
    hAdapCuBox_CuFramek40M2->SetBinContent(i, hnewCuBox_CuFramek40M2->GetBinContent(i)/hnewCuBox_CuFramek40M2->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M2->SetBinContent(i, hnewCuBox_CuFrameth232M2->GetBinContent(i)/hnewCuBox_CuFrameth232M2->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2->SetBinContent(i, hnewCuBox_CuFrameu238M2->GetBinContent(i)/hnewCuBox_CuFrameu238M2->GetBinWidth(i));

    hAdapCuBox_CuFrameth232M2_10->SetBinContent(i, hnewCuBox_CuFrameth232M2_10->GetBinContent(i)/hnewCuBox_CuFrameth232M2_10->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2_10->SetBinContent(i, hnewCuBox_CuFrameu238M2_10->GetBinContent(i)/hnewCuBox_CuFrameu238M2_10->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2_10->SetBinContent(i, hnewCuBox_CuFramepb210M2_10->GetBinContent(i)/hnewCuBox_CuFramepb210M2_10->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2_1->SetBinContent(i, hnewCuBox_CuFramepb210M2_1->GetBinContent(i)/hnewCuBox_CuFramepb210M2_1->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2_01->SetBinContent(i, hnewCuBox_CuFramepb210M2_01->GetBinContent(i)/hnewCuBox_CuFramepb210M2_01->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2_001->SetBinContent(i, hnewCuBox_CuFramepb210M2_001->GetBinContent(i)/hnewCuBox_CuFramepb210M2_001->GetBinWidth(i));

    hAdapCuBox_CuFrameth232M2_1->SetBinContent(i, hnewCuBox_CuFrameth232M2_1->GetBinContent(i)/hnewCuBox_CuFrameth232M2_1->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2_1->SetBinContent(i, hnewCuBox_CuFrameu238M2_1->GetBinContent(i)/hnewCuBox_CuFrameu238M2_1->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M2_01->SetBinContent(i, hnewCuBox_CuFrameth232M2_01->GetBinContent(i)/hnewCuBox_CuFrameth232M2_01->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2_01->SetBinContent(i, hnewCuBox_CuFrameu238M2_01->GetBinContent(i)/hnewCuBox_CuFrameu238M2_01->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M2_001->SetBinContent(i, hnewCuBox_CuFrameth232M2_001->GetBinContent(i)/hnewCuBox_CuFrameth232M2_001->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2_001->SetBinContent(i, hnewCuBox_CuFrameu238M2_001->GetBinContent(i)/hnewCuBox_CuFrameu238M2_001->GetBinWidth(i));  

    hAdapCuBox_CuFrameth232M2_100->SetBinContent(i, hnewCuBox_CuFrameth232M2_100->GetBinContent(i)/hnewCuBox_CuFrameth232M2_100->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2_100->SetBinContent(i, hnewCuBox_CuFrameu238M2_100->GetBinContent(i)/hnewCuBox_CuFrameu238M2_100->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2_100->SetBinContent(i, hnewCuBox_CuFramepb210M2_100->GetBinContent(i)/hnewCuBox_CuFramepb210M2_100->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M2_50->SetBinContent(i, hnewCuBox_CuFrameth232M2_50->GetBinContent(i)/hnewCuBox_CuFrameth232M2_50->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2_50->SetBinContent(i, hnewCuBox_CuFrameu238M2_50->GetBinContent(i)/hnewCuBox_CuFrameu238M2_50->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2_50->SetBinContent(i, hnewCuBox_CuFramepb210M2_50->GetBinContent(i)/hnewCuBox_CuFramepb210M2_50->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M2_5->SetBinContent(i, hnewCuBox_CuFrameth232M2_5->GetBinContent(i)/hnewCuBox_CuFrameth232M2_5->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2_5->SetBinContent(i, hnewCuBox_CuFrameu238M2_5->GetBinContent(i)/hnewCuBox_CuFrameu238M2_5->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2_5->SetBinContent(i, hnewCuBox_CuFramepb210M2_5->GetBinContent(i)/hnewCuBox_CuFramepb210M2_5->GetBinWidth(i));

    hAdap50mKcs137M2->SetBinContent(i, hnew50mKcs137M2->GetBinContent(i)/hnew50mKcs137M2->GetBinWidth(i));

    hAdapPbRomco60M2->SetBinContent(i, hnewPbRomco60M2->GetBinContent(i)/hnewPbRomco60M2->GetBinWidth(i));
    hAdapPbRomcs137M2->SetBinContent(i, hnewPbRomcs137M2->GetBinContent(i)/hnewPbRomcs137M2->GetBinWidth(i));
    hAdapPbRomk40M2->SetBinContent(i, hnewPbRomk40M2->GetBinContent(i)/hnewPbRomk40M2->GetBinWidth(i));
    hAdapPbRompb210M2->SetBinContent(i, hnewPbRompb210M2->GetBinContent(i)/hnewPbRompb210M2->GetBinWidth(i));
    hAdapPbRomth232M2->SetBinContent(i, hnewPbRomth232M2->GetBinContent(i)/hnewPbRomth232M2->GetBinWidth(i));
    hAdapPbRomu238M2->SetBinContent(i, hnewPbRomu238M2->GetBinContent(i)/hnewPbRomu238M2->GetBinWidth(i));

    hAdapInternalco60M2->SetBinContent(i, hnewInternalco60M2->GetBinContent(i)/hnewInternalco60M2->GetBinWidth(i));
    hAdapInternalk40M2->SetBinContent(i, hnewInternalk40M2->GetBinContent(i)/hnewInternalk40M2->GetBinWidth(i));
    hAdapInternalth232M2->SetBinContent(i, hnewInternalth232M2->GetBinContent(i)/hnewInternalth232M2->GetBinWidth(i));
    hAdapInternalu238M2->SetBinContent(i, hnewInternalu238M2->GetBinContent(i)/hnewInternalu238M2->GetBinWidth(i));

    hAdapOVCco60M2->SetBinContent(i, hnewOVCco60M2->GetBinContent(i)/hnewOVCco60M2->GetBinWidth(i));
    hAdapOVCk40M2->SetBinContent(i, hnewOVCk40M2->GetBinContent(i)/hnewOVCk40M2->GetBinWidth(i));
    hAdapOVCth232M2->SetBinContent(i, hnewOVCth232M2->GetBinContent(i)/hnewOVCth232M2->GetBinWidth(i));
    hAdapOVCu238M2->SetBinContent(i, hnewOVCu238M2->GetBinContent(i)/hnewOVCu238M2->GetBinWidth(i));

    hAdapExtPbbi210M2->SetBinContent(i, hnewExtPbbi210M2->GetBinContent(i)/hnewExtPbbi210M2->GetBinWidth(i));

    // Fudge
    hAdapOVC804M2->SetBinContent(i, hnewOVC804M2->GetBinContent(i)/hnewOVC804M2->GetBinWidth(i));
    hAdapOVC835M2->SetBinContent(i, hnewOVC835M2->GetBinContent(i)/hnewOVC835M2->GetBinWidth(i));
    hAdapOVC1063M2->SetBinContent(i, hnewOVC1063M2->GetBinContent(i)/hnewOVC1063M2->GetBinWidth(i));
  }


  for(int i = 1; i <= dAdaptiveBinsM2Sum; i++)
  {
    hAdapTeO20nuM2Sum->SetBinContent(i, hnewTeO20nuM2Sum->GetBinContent(i)/hnewTeO20nuM2Sum->GetBinWidth(i));
    hAdapTeO22nuM2Sum->SetBinContent(i, hnewTeO22nuM2Sum->GetBinContent(i)/hnewTeO22nuM2Sum->GetBinWidth(i));
    hAdapTeO2co60M2Sum->SetBinContent(i, hnewTeO2co60M2Sum->GetBinContent(i)/hnewTeO2co60M2Sum->GetBinWidth(i));
    hAdapTeO2k40M2Sum->SetBinContent(i, hnewTeO2k40M2Sum->GetBinContent(i)/hnewTeO2k40M2Sum->GetBinWidth(i));
    hAdapTeO2pb210M2Sum->SetBinContent(i, hnewTeO2pb210M2Sum->GetBinContent(i)/hnewTeO2pb210M2Sum->GetBinWidth(i));
    hAdapTeO2po210M2Sum->SetBinContent(i, hnewTeO2po210M2Sum->GetBinContent(i)/hnewTeO2po210M2Sum->GetBinWidth(i));
    hAdapTeO2te125M2Sum->SetBinContent(i, hnewTeO2te125M2Sum->GetBinContent(i)/hnewTeO2te125M2Sum->GetBinWidth(i));
    hAdapTeO2th232M2Sum->SetBinContent(i, hnewTeO2th232M2Sum->GetBinContent(i)/hnewTeO2th232M2Sum->GetBinWidth(i));
    // hAdapTeO2th228M2Sum->SetBinContent(i, hnewTeO2th228M2Sum->GetBinContent(i)/hnewTeO2th228M2Sum->GetBinWidth(i));
    // hAdapTeO2ra226M2Sum->SetBinContent(i, hnewTeO2ra226M2Sum->GetBinContent(i)/hnewTeO2ra226M2Sum->GetBinWidth(i));
    // hAdapTeO2rn222M2Sum->SetBinContent(i, hnewTeO2rn222M2Sum->GetBinContent(i)/hnewTeO2rn222M2Sum->GetBinWidth(i));
    hAdapTeO2u238M2Sum->SetBinContent(i, hnewTeO2u238M2Sum->GetBinContent(i)/hnewTeO2u238M2Sum->GetBinWidth(i));
    // hAdapTeO2th230M2Sum->SetBinContent(i, hnewTeO2th230M2Sum->GetBinContent(i)/hnewTeO2th230M2Sum->GetBinWidth(i));
    // hAdapTeO2u234M2Sum->SetBinContent(i, hnewTeO2u234M2Sum->GetBinContent(i)/hnewTeO2u234M2Sum->GetBinWidth(i));

    hAdapTeO2Sxpb210M2Sum_001->SetBinContent(i, hnewTeO2Sxpb210M2Sum_001->GetBinContent(i)/hnewTeO2Sxpb210M2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M2Sum_01->SetBinContent(i, hnewTeO2Sxpb210M2Sum_01->GetBinContent(i)/hnewTeO2Sxpb210M2Sum_01->GetBinWidth(i));
    hAdapTeO2Sxpb210M2Sum_1->SetBinContent(i, hnewTeO2Sxpb210M2Sum_1->GetBinContent(i)/hnewTeO2Sxpb210M2Sum_1->GetBinWidth(i));
    hAdapTeO2Sxpb210M2Sum_10->SetBinContent(i, hnewTeO2Sxpb210M2Sum_10->GetBinContent(i)/hnewTeO2Sxpb210M2Sum_10->GetBinWidth(i));
    hAdapTeO2Sxpo210M2Sum_001->SetBinContent(i, hnewTeO2Sxpo210M2Sum_001->GetBinContent(i)/hnewTeO2Sxpo210M2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxpo210M2Sum_01->SetBinContent(i, hnewTeO2Sxpo210M2Sum_01->GetBinContent(i)/hnewTeO2Sxpo210M2Sum_01->GetBinWidth(i));
    hAdapTeO2Sxpo210M2Sum_1->SetBinContent(i, hnewTeO2Sxpo210M2Sum_1->GetBinContent(i)/hnewTeO2Sxpo210M2Sum_1->GetBinWidth(i));
    hAdapTeO2Sxth232M2Sum_001->SetBinContent(i, hnewTeO2Sxth232M2Sum_001->GetBinContent(i)/hnewTeO2Sxth232M2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxth232M2Sum_01->SetBinContent(i, hnewTeO2Sxth232M2Sum_01->GetBinContent(i)/hnewTeO2Sxth232M2Sum_01->GetBinWidth(i));
    hAdapTeO2Sxth232M2Sum_1->SetBinContent(i, hnewTeO2Sxth232M2Sum_1->GetBinContent(i)/hnewTeO2Sxth232M2Sum_1->GetBinWidth(i));
    hAdapTeO2Sxth232M2Sum_10->SetBinContent(i, hnewTeO2Sxth232M2Sum_10->GetBinContent(i)/hnewTeO2Sxth232M2Sum_10->GetBinWidth(i));
    hAdapTeO2Sxu238M2Sum_001->SetBinContent(i, hnewTeO2Sxu238M2Sum_001->GetBinContent(i)/hnewTeO2Sxu238M2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxu238M2Sum_01->SetBinContent(i, hnewTeO2Sxu238M2Sum_01->GetBinContent(i)/hnewTeO2Sxu238M2Sum_01->GetBinWidth(i));
    hAdapTeO2Sxu238M2Sum_1->SetBinContent(i, hnewTeO2Sxu238M2Sum_1->GetBinContent(i)/hnewTeO2Sxu238M2Sum_1->GetBinWidth(i));
    hAdapTeO2Sxu238M2Sum_10->SetBinContent(i, hnewTeO2Sxu238M2Sum_10->GetBinContent(i)/hnewTeO2Sxu238M2Sum_10->GetBinWidth(i));

    hAdapTeO2Sxu238M2Sum_100->SetBinContent(i, hnewTeO2Sxu238M2Sum_100->GetBinContent(i)/hnewTeO2Sxu238M2Sum_100->GetBinWidth(i));
    hAdapTeO2Sxth232M2Sum_100->SetBinContent(i, hnewTeO2Sxth232M2Sum_100->GetBinContent(i)/hnewTeO2Sxth232M2Sum_100->GetBinWidth(i));
    hAdapTeO2Sxpb210M2Sum_100->SetBinContent(i, hnewTeO2Sxpb210M2Sum_100->GetBinContent(i)/hnewTeO2Sxpb210M2Sum_100->GetBinWidth(i));

    hAdapTeO2th232onlyM2Sum->SetBinContent(i, hnewTeO2th232onlyM2Sum->GetBinContent(i)/hnewTeO2th232onlyM2Sum->GetBinWidth(i));
    hAdapTeO2ra228pb208M2Sum->SetBinContent(i, hnewTeO2ra228pb208M2Sum->GetBinContent(i)/hnewTeO2ra228pb208M2Sum->GetBinWidth(i));
    hAdapTeO2th230onlyM2Sum->SetBinContent(i, hnewTeO2th230onlyM2Sum->GetBinContent(i)/hnewTeO2th230onlyM2Sum->GetBinWidth(i));
    hAdapTeO2u238th230M2Sum->SetBinContent(i, hnewTeO2u238th230M2Sum->GetBinContent(i)/hnewTeO2u238th230M2Sum->GetBinWidth(i));
    hAdapTeO2ra226pb210M2Sum->SetBinContent(i, hnewTeO2ra226pb210M2Sum->GetBinContent(i)/hnewTeO2ra226pb210M2Sum->GetBinWidth(i));

    hAdapTeO2Sxth232onlyM2Sum_001->SetBinContent(i, hnewTeO2Sxth232onlyM2Sum_001->GetBinContent(i)/hnewTeO2Sxth232onlyM2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M2Sum_001->SetBinContent(i, hnewTeO2Sxra228pb208M2Sum_001->GetBinContent(i)/hnewTeO2Sxra228pb208M2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxu238th230M2Sum_001->SetBinContent(i, hnewTeO2Sxu238th230M2Sum_001->GetBinContent(i)/hnewTeO2Sxu238th230M2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM2Sum_001->SetBinContent(i, hnewTeO2Sxth230onlyM2Sum_001->GetBinContent(i)/hnewTeO2Sxth230onlyM2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M2Sum_001->SetBinContent(i, hnewTeO2Sxra226pb210M2Sum_001->GetBinContent(i)/hnewTeO2Sxra226pb210M2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M2Sum_0001->SetBinContent(i, hnewTeO2Sxpb210M2Sum_0001->GetBinContent(i)/hnewTeO2Sxpb210M2Sum_0001->GetBinWidth(i));

    hAdapTeO2Sxth232onlyM2Sum_01->SetBinContent(i, hnewTeO2Sxth232onlyM2Sum_01->GetBinContent(i)/hnewTeO2Sxth232onlyM2Sum_01->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M2Sum_01->SetBinContent(i, hnewTeO2Sxra228pb208M2Sum_01->GetBinContent(i)/hnewTeO2Sxra228pb208M2Sum_01->GetBinWidth(i));
    hAdapTeO2Sxu238th230M2Sum_01->SetBinContent(i, hnewTeO2Sxu238th230M2Sum_01->GetBinContent(i)/hnewTeO2Sxu238th230M2Sum_01->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM2Sum_01->SetBinContent(i, hnewTeO2Sxth230onlyM2Sum_01->GetBinContent(i)/hnewTeO2Sxth230onlyM2Sum_01->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M2Sum_01->SetBinContent(i, hnewTeO2Sxra226pb210M2Sum_01->GetBinContent(i)/hnewTeO2Sxra226pb210M2Sum_01->GetBinWidth(i));

    hAdapTeO2Sxth232onlyM2Sum_0001->SetBinContent(i, hnewTeO2Sxth232onlyM2Sum_0001->GetBinContent(i)/hnewTeO2Sxth232onlyM2Sum_0001->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M2Sum_0001->SetBinContent(i, hnewTeO2Sxra228pb208M2Sum_0001->GetBinContent(i)/hnewTeO2Sxra228pb208M2Sum_0001->GetBinWidth(i));
    hAdapTeO2Sxu238th230M2Sum_0001->SetBinContent(i, hnewTeO2Sxu238th230M2Sum_0001->GetBinContent(i)/hnewTeO2Sxu238th230M2Sum_0001->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM2Sum_0001->SetBinContent(i, hnewTeO2Sxth230onlyM2Sum_0001->GetBinContent(i)/hnewTeO2Sxth230onlyM2Sum_0001->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M2Sum_0001->SetBinContent(i, hnewTeO2Sxra226pb210M2Sum_0001->GetBinContent(i)/hnewTeO2Sxra226pb210M2Sum_0001->GetBinWidth(i));

    hAdapCuBox_CuFrameco60M2Sum->SetBinContent(i, hnewCuBox_CuFrameco60M2Sum->GetBinContent(i)/hnewCuBox_CuFrameco60M2Sum->GetBinWidth(i));
    hAdapCuBox_CuFramek40M2Sum->SetBinContent(i, hnewCuBox_CuFramek40M2Sum->GetBinContent(i)/hnewCuBox_CuFramek40M2Sum->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M2Sum->SetBinContent(i, hnewCuBox_CuFrameth232M2Sum->GetBinContent(i)/hnewCuBox_CuFrameth232M2Sum->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2Sum->SetBinContent(i, hnewCuBox_CuFrameu238M2Sum->GetBinContent(i)/hnewCuBox_CuFrameu238M2Sum->GetBinWidth(i));

    hAdapCuBox_CuFrameth232M2Sum_10->SetBinContent(i, hnewCuBox_CuFrameth232M2Sum_10->GetBinContent(i)/hnewCuBox_CuFrameth232M2Sum_10->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2Sum_10->SetBinContent(i, hnewCuBox_CuFrameu238M2Sum_10->GetBinContent(i)/hnewCuBox_CuFrameu238M2Sum_10->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2Sum_10->SetBinContent(i, hnewCuBox_CuFramepb210M2Sum_10->GetBinContent(i)/hnewCuBox_CuFramepb210M2Sum_10->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2Sum_1->SetBinContent(i, hnewCuBox_CuFramepb210M2Sum_1->GetBinContent(i)/hnewCuBox_CuFramepb210M2Sum_1->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2Sum_01->SetBinContent(i, hnewCuBox_CuFramepb210M2Sum_01->GetBinContent(i)/hnewCuBox_CuFramepb210M2Sum_01->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2Sum_001->SetBinContent(i, hnewCuBox_CuFramepb210M2Sum_001->GetBinContent(i)/hnewCuBox_CuFramepb210M2Sum_001->GetBinWidth(i));

    hAdapCuBox_CuFrameth232M2Sum_1->SetBinContent(i, hnewCuBox_CuFrameth232M2Sum_1->GetBinContent(i)/hnewCuBox_CuFrameth232M2Sum_1->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2Sum_1->SetBinContent(i, hnewCuBox_CuFrameu238M2Sum_1->GetBinContent(i)/hnewCuBox_CuFrameu238M2Sum_1->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M2Sum_01->SetBinContent(i, hnewCuBox_CuFrameth232M2Sum_01->GetBinContent(i)/hnewCuBox_CuFrameth232M2Sum_01->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2Sum_01->SetBinContent(i, hnewCuBox_CuFrameu238M2Sum_01->GetBinContent(i)/hnewCuBox_CuFrameu238M2Sum_01->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M2Sum_001->SetBinContent(i, hnewCuBox_CuFrameth232M2Sum_001->GetBinContent(i)/hnewCuBox_CuFrameth232M2Sum_001->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2Sum_001->SetBinContent(i, hnewCuBox_CuFrameu238M2Sum_001->GetBinContent(i)/hnewCuBox_CuFrameu238M2Sum_001->GetBinWidth(i));  

    hAdapCuBox_CuFrameth232M2Sum_100->SetBinContent(i, hnewCuBox_CuFrameth232M2Sum_100->GetBinContent(i)/hnewCuBox_CuFrameth232M2Sum_100->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2Sum_100->SetBinContent(i, hnewCuBox_CuFrameu238M2Sum_100->GetBinContent(i)/hnewCuBox_CuFrameu238M2Sum_100->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2Sum_100->SetBinContent(i, hnewCuBox_CuFramepb210M2Sum_100->GetBinContent(i)/hnewCuBox_CuFramepb210M2Sum_100->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M2Sum_50->SetBinContent(i, hnewCuBox_CuFrameth232M2Sum_50->GetBinContent(i)/hnewCuBox_CuFrameth232M2Sum_50->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2Sum_50->SetBinContent(i, hnewCuBox_CuFrameu238M2Sum_50->GetBinContent(i)/hnewCuBox_CuFrameu238M2Sum_50->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2Sum_50->SetBinContent(i, hnewCuBox_CuFramepb210M2Sum_50->GetBinContent(i)/hnewCuBox_CuFramepb210M2Sum_50->GetBinWidth(i));
    hAdapCuBox_CuFrameth232M2Sum_5->SetBinContent(i, hnewCuBox_CuFrameth232M2Sum_5->GetBinContent(i)/hnewCuBox_CuFrameth232M2Sum_5->GetBinWidth(i));
    hAdapCuBox_CuFrameu238M2Sum_5->SetBinContent(i, hnewCuBox_CuFrameu238M2Sum_5->GetBinContent(i)/hnewCuBox_CuFrameu238M2Sum_5->GetBinWidth(i));
    hAdapCuBox_CuFramepb210M2Sum_5->SetBinContent(i, hnewCuBox_CuFramepb210M2Sum_5->GetBinContent(i)/hnewCuBox_CuFramepb210M2Sum_5->GetBinWidth(i));


    hAdap50mKcs137M2Sum->SetBinContent(i, hnew50mKcs137M2Sum->GetBinContent(i)/hnew50mKcs137M2Sum->GetBinWidth(i));

    hAdapPbRomco60M2Sum->SetBinContent(i, hnewPbRomco60M2Sum->GetBinContent(i)/hnewPbRomco60M2Sum->GetBinWidth(i));
    hAdapPbRomcs137M2Sum->SetBinContent(i, hnewPbRomcs137M2Sum->GetBinContent(i)/hnewPbRomcs137M2Sum->GetBinWidth(i));
    hAdapPbRomk40M2Sum->SetBinContent(i, hnewPbRomk40M2Sum->GetBinContent(i)/hnewPbRomk40M2Sum->GetBinWidth(i));
    hAdapPbRompb210M2Sum->SetBinContent(i, hnewPbRompb210M2Sum->GetBinContent(i)/hnewPbRompb210M2Sum->GetBinWidth(i));
    hAdapPbRomth232M2Sum->SetBinContent(i, hnewPbRomth232M2Sum->GetBinContent(i)/hnewPbRomth232M2Sum->GetBinWidth(i));
    hAdapPbRomu238M2Sum->SetBinContent(i, hnewPbRomu238M2Sum->GetBinContent(i)/hnewPbRomu238M2Sum->GetBinWidth(i));

    hAdapInternalco60M2Sum->SetBinContent(i, hnewInternalco60M2Sum->GetBinContent(i)/hnewInternalco60M2Sum->GetBinWidth(i));
    hAdapInternalk40M2Sum->SetBinContent(i, hnewInternalk40M2Sum->GetBinContent(i)/hnewInternalk40M2Sum->GetBinWidth(i));
    hAdapInternalth232M2Sum->SetBinContent(i, hnewInternalth232M2Sum->GetBinContent(i)/hnewInternalth232M2Sum->GetBinWidth(i));
    hAdapInternalu238M2Sum->SetBinContent(i, hnewInternalu238M2Sum->GetBinContent(i)/hnewInternalu238M2Sum->GetBinWidth(i));

    hAdapOVCco60M2Sum->SetBinContent(i, hnewOVCco60M2Sum->GetBinContent(i)/hnewOVCco60M2Sum->GetBinWidth(i));
    hAdapOVCk40M2Sum->SetBinContent(i, hnewOVCk40M2Sum->GetBinContent(i)/hnewOVCk40M2Sum->GetBinWidth(i));
    hAdapOVCth232M2Sum->SetBinContent(i, hnewOVCth232M2Sum->GetBinContent(i)/hnewOVCth232M2Sum->GetBinWidth(i));
    hAdapOVCu238M2Sum->SetBinContent(i, hnewOVCu238M2Sum->GetBinContent(i)/hnewOVCu238M2Sum->GetBinWidth(i));

    hAdapExtPbbi210M2Sum->SetBinContent(i, hnewExtPbbi210M2Sum->GetBinContent(i)/hnewExtPbbi210M2Sum->GetBinWidth(i));

    // Fudge
    hAdapOVC804M2Sum->SetBinContent(i, hnewOVC804M2Sum->GetBinContent(i)/hnewOVC804M2Sum->GetBinWidth(i));
    hAdapOVC835M2Sum->SetBinContent(i, hnewOVC835M2Sum->GetBinContent(i)/hnewOVC835M2Sum->GetBinWidth(i));
    hAdapOVC1063M2Sum->SetBinContent(i, hnewOVC1063M2Sum->GetBinContent(i)/hnewOVC1063M2Sum->GetBinWidth(i));
  }

  SetParEfficiency();
}

// Loads the background data
void TBackgroundModel::LoadData()
{
  switch(dDataSet)
  { 
  case 1:
    // qtree->Add(Form("%s/Unblinded/ReducedB-ds2049.root", dDataDir.c_str()));   
    // qtree->Add(Form("%s/Unblinded/ReducedB-ds2061.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2064.root", dDataDir.c_str()));   
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2067.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2070.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2073.root", dDataDir.c_str())); 
    // qtree->Add(Form("%s/Unblinded/ReducedB-ds2076.root", dDataDir.c_str())); 
    // dLivetime = 6042498; // DR 1 
    // dLivetime = 1082802+867904+429707+1288128+1146692+291436; // DR 1 without DS 2049
    // dLivetime = 5106669;
    dLivetime = 867904+429707+1288128+1146692; // Removing runs with high radon (2061 and 2076)
    cout << "Using Data Release 1" << endl;
  break;

  case 2:
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2079.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2085.root", dDataDir.c_str())); 
    // qtree->Add(Form("%s/Unblinded/ReducedB-ds2088.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2091.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2097.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2100.root", dDataDir.c_str())); 
    // dLivetime = 9387524; // DR 2
    dLivetime = 753416+2006480+1675946+1558268+1835197;// Removing high radon (2088)
    cout << "Using Data Release 2" << endl;
  break;

  case 3:
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2103.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2109.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2118.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2124.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2130.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2133.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2139.root", dDataDir.c_str()));
    // dLivetime = 1967358+2374242+1301016+2005292+1511379+901597+1639017; // DR 3
    dLivetime = 11699901;
    cout << "Using Data Release 3" << endl;
  break;

  // Combine DR 2 and 3
  case 4:
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2079.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2085.root", dDataDir.c_str())); 
    // qtree->Add(Form("%s/Unblinded/ReducedB-ds2088.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2091.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2097.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2100.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2103.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2109.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2118.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2124.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2130.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2133.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2139.root", dDataDir.c_str()));
    // dLivetime = 1967358+2374242+1301016+2005292+1511379+901597+1639017; // DR 3
    // dLivetime = 9387524+11699901; 
    dLivetime = 753416+2006480+1675946+1558268+1835197+11699901; // Removing high Radon runs (2088)
    cout << "Using Data Release 2+3" << endl;
  break;

  case -1:
    cout << "Using Toy data" << endl;
    bToyData = true;
    qtree->Add(Form("%s/Unblinded/ReducedB-ds*.root", dDataDir.c_str())); 
    // dLivetime = 23077930+1511379+901597; // seconds of livetime (DR1 to DR3)
    dLivetime = 11699901+9387524+5106669;
  break;

  default:   
    // qtree->Add(Form("%s/Unblinded/ReducedB-ds2061.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2064.root", dDataDir.c_str()));   
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2067.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2070.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2073.root", dDataDir.c_str())); 
    // qtree->Add(Form("%s/Unblinded/ReducedB-ds2076.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2079.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2085.root", dDataDir.c_str())); 
    // qtree->Add(Form("%s/Unblinded/ReducedB-ds2088.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2091.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2097.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2100.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2103.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2109.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2118.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2124.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2130.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2133.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2139.root", dDataDir.c_str()));
    // dLivetime = 23077930+1511379+901597; // seconds of livetime (DR1 to DR3)
    // dLivetime = 11699901+9387524+5106669;
    dLivetime = 867904+429707+1288128+1146692+753416+2006480+1675946+1558268+1835197+11699901; // Removing high radon
    cout << "Using Total Dataset" << endl;
  }

  qtree->Project("fDataHistoTot", "Energy", base_cut);
  qtree->Project("fDataHistoM1",  "Energy", base_cut && "Multiplicity_Sync == 1");
  qtree->Project("fDataHistoM2",  "Energy", base_cut && "Multiplicity_Sync == 2");
  qtree->Project("fDataHistoM2Sum",  "TotalEnergy_Sync", base_cut && "Multiplicity_Sync == 2");

  dLivetimeYr = dLivetime*dSecToYears;  
}

// Prints parameters, needs update 11-06-2014
void TBackgroundModel::PrintParameters()
{
  for(int i = 0; i < TBackgroundModel::dNParam; i++)
  {
    cout << i << " " << minuit->fCpnam[i] << ": " << fParameters[i] << " +/- " << fParError[i] << endl;
  }
}

void TBackgroundModel::PrintParActivity()
{ 
  // This only gets the number of counts corrected for detector efficiency
  cout << "Activity in Counts/Year (M1)" << endl;
  for(int i = 0; i < TBackgroundModel::dNParam; i++)
  {
  // cout << i << " " << minuit->fCpnam[i] << " activity: " << fParActivityM1[i] << endl;
  // cout << fParActivityM1[i] << "\t" <<  dDataIntegralM1*fParError[i] << endl;
    if(!bFixedArray[i])
    {
      cout << i << " " << fParActivityM1[i] << " " << fParActivityM1[i]*fParError[i]/fParameters[i] << " " << fParActivityM2[i] <<  " " << fParActivityM2[i]*fParError[i]/fParameters[i]  << endl;
    }
  }
}

void TBackgroundModel::GenerateParameters()
{
  // Initialization (Name, Index, Initial Value, Min Limit, Max Limit, pointer to histograms.. )
  // TeO2 Bulk
  BkgPar[0] = new TBkgModelParameter( "TeO2 2nu", 0, 0, 1E-7, 0, 1.0, hAdapTeO22nuM1, hAdapTeO22nuM2 );

  // Total data Alpha
  BkgPar[1] = new TBkgModelParameter( "TeO2 k40", 1, 0., 1E-7, 0, 1.0, hAdapTeO2k40M1, hAdapTeO2k40M2 ); 
  BkgPar[2] = new TBkgModelParameter( "TeO2 co60", 2, 0., 1E-7, 0, 1.0, hAdapTeO2co60M1 , hAdapTeO2co60M2 );  
  BkgPar[3] = new TBkgModelParameter( "TeO2 sb125", 3, 0., 1E-7, 0, 1.0, hAdapTeO2sb125M1 , hAdapTeO2sb125M2 );

  // TeO2 Surface
  BkgPar[4] = new TBkgModelParameter( "TeO2 Sx th232 1", 4, 2.85569e-03, 1E-7, 0, 1.0, hAdapTeO2Sxth232M1_1 , hAdapTeO2Sxth232M2_1 );
  BkgPar[5] = new TBkgModelParameter( "TeO2 Sx th232 only 0.01", 5, 1.39386e-03, 1E-7, 0, 1.0, hAdapTeO2Sxth232M1_001 , hAdapTeO2Sxth232M2_001 );
  BkgPar[6] = new TBkgModelParameter( "TeO2 Sx ra228-pb208 0.01", 6, 2.15768e-03, 1E-7, 0, 1.0, hAdapTeO2Sxra228pb208M1_001 , hAdapTeO2Sxra228pb208M2_001 );
  BkgPar[7] = new TBkgModelParameter( "TeO2 Sx u238-th230 0.01", 7, 2.05704e-03, 1E-7, 0, 1.0, hAdapTeO2Sxu238th230M1_001 , hAdapTeO2Sxu238th230M2_001 );
  BkgPar[8] = new TBkgModelParameter( "TeO2 Sx th230 only 0.01", 8, 1.42958e-03, 1E-7, 0, 1.0, hAdapTeO2Sxth230onlyM1_001 , hAdapTeO2Sxth230onlyM2_001 );
  BkgPar[9] = new TBkgModelParameter( "TeO2 Sx ra226-pb210 0.01", 9, 3.08708e-03, 1E-7, 0, 1.0, hAdapTeO2Sxra226pb210M1_001 , hAdapTeO2Sxra226pb210M2_001 );
  BkgPar[10] = new TBkgModelParameter( "TeO2 Sx pb210 0.1", 10, 4.86040e-02, 1E-7, 0, 1.0, hAdapTeO2Sxpb210M1_01 , hAdapTeO2Sxpb210M2_01 ); // Works better than 0.01

  // Cu Holder Surface
  BkgPar[11] = new TBkgModelParameter( "CuBox+CuFrame Sx u238 100", 11, 1.00599e-02, 1E-7, 0, 1.0, hAdapCuBox_CuFrameu238M1_100 , hAdapCuBox_CuFrameu238M2_100 );
  BkgPar[12] = new TBkgModelParameter( "CuBox+CuFrame Sx th232 100", 12, 4.18960e-02, 1E-7, 0, 1.0, hAdapCuBox_CuFrameth232M1_100 , hAdapCuBox_CuFrameth232M2_100 );
  BkgPar[13] = new TBkgModelParameter( "CuBox+CuFrame Sx u238 10", 13, 6.28868e-03, 1E-7, 0, 1.0, hAdapCuBox_CuFrameu238M1_10 , hAdapCuBox_CuFrameu238M2_10 );
  BkgPar[14] = new TBkgModelParameter( "CuBox+CuFrame Sx pb210 1", 14, 2.31940e-03, 1E-7, 0, 1.0, hAdapCuBox_CuFramepb210M1_1 , hAdapCuBox_CuFramepb210M2_1 );
  BkgPar[15] = new TBkgModelParameter( "CuBox+CuFrame Sx pb210 0.1", 15, 3.50695e-03, 1E-7, 0, 1.0, hAdapCuBox_CuFramepb210M1_01 , hAdapCuBox_CuFramepb210M2_01 );
  BkgPar[16] = new TBkgModelParameter( "CuBox+CuFrame Sx pb210 0.01", 16, 2.25633e-02, 1E-7, 0, 1.0, hAdapCuBox_CuFramepb210M1_001 , hAdapCuBox_CuFramepb210M2_001 );
  BkgPar[17] = new TBkgModelParameter( "CuBox+CuFrame Sx th232 0.01", 17, 8.90411e-05, 1E-7, 0, 1.0, hAdapCuBox_CuFrameth232M1_001 , hAdapCuBox_CuFrameth232M2_001 );

  // Bulk Contamination
  BkgPar[18] = new TBkgModelParameter( "PbRom k40", 18, 0., 1E-7, 0, 1.0, hAdapPbRomk40M1 , hAdapPbRomk40M2 );
  BkgPar[19] = new TBkgModelParameter( "50mK cs137", 19, 1.90062e-03, 1E-7, 0, 1.0, hAdap50mKcs137M1 , hAdap50mKcs137M2 );

  BkgPar[20] = new TBkgModelParameter( "OVC u238", 20, 1.41680e-01, 1E-7, 0, 1.0, hAdapOVCu238M1 , hAdapOVCu238M2 );
  BkgPar[21] = new TBkgModelParameter( "OVC th232", 21, 1.00966e-01, 1E-7, 0, 1.0, hAdapOVCth232M1 , hAdapOVCth232M2 );
  BkgPar[22] = new TBkgModelParameter( "OVC k40", 22, 7.11104e-02, 1E-7, 0, 1.0, hAdapOVCk40M1 , hAdapOVCk40M2 );
  BkgPar[23] = new TBkgModelParameter( "OVC co60", 23, 3.30743e-02, 1E-7, 0, 1.0, hAdapOVCco60M1 , hAdapOVCco60M2 );
  BkgPar[24] = new TBkgModelParameter( "External Lead bi210", 24, 3.35197e-07, 1E-7, 0, 1.0, hAdapExtPbbi210M1 , hAdapExtPbbi210M2 );
  BkgPar[25] = new TBkgModelParameter( "OVC 804", 25, 3.77822e-03, 1E-7, 0, 1.0, hAdapOVC804M1 , hAdapOVC804M2 );
  BkgPar[26] = new TBkgModelParameter( "OVC 835", 26, 2.69330e-03, 1E-7, 0, 1.0, hAdapOVC835M1 , hAdapOVC835M2 );
  BkgPar[27] = new TBkgModelParameter( "OVC 1063", 27, 5.20414e-03, 1E-7, 0, 1.0, hAdapOVC1063M1 , hAdapOVC1063M2 ); 


  bFixedArray[1] = true;

  // bFixedArray[7] = true;
  // bFixedArray[8] = true;
  // bFixedArray[9] = true;
  
  for(int i = 4; i < 18; i ++)
  {
    bFixedArray[i] = true;
  }


  // for(int i = 37; i < 59; i ++)
  // {
  //   bFixedArray[i] = true;
  // }


}


// Resets all parameters to 0
// This resets the combined histograms (the model)
void TBackgroundModel::ResetParameters()
{
  dNumCalls = 0;
  dChiSquare = 0;
  
  fModelTotM1->Reset();
  fModelTotAdapM1->Reset();

  // Total PDFs M2
  fModelTotM2->Reset();
  fModelTotAdapM2->Reset();

  // Total PDFs M2Sum
  fModelTotM2Sum->Reset();
  fModelTotAdapM2Sum->Reset();

  for(int i = 0; i < TBackgroundModel::dNParam; i++)
  {
    fParameters[i] = 0;
    fParError[i] = 0;
  }
}

// Set Parameters in Model
void TBackgroundModel::SetParameters(int index, double value)
{
	// Change the index max depending on model
	if(index > dNParam) cout << "Index too large" << endl;
	else fParameters[index] = value;
}

void TBackgroundModel::SetParEfficiency()
{
   fParEfficiencyM1[0] = 0.9292826; // TeO2 2nu
   fParEfficiencyM1[1] = 0.5614488; // CuBox+Frame co60
   fParEfficiencyM1[2] = 0.941975; // TeO2 th232 only
   fParEfficiencyM1[3] = 0.9419295; // TeO2 th230 only
   fParEfficiencyM1[4] = 0.4158284; // TeO2 Sx th232 0.01
   fParEfficiencyM1[5] = 0.3802244; // TeO2 Sx ra228 to pb208 0.01 
   fParEfficiencyM1[6] = 0.7457738; // TeO2 Sx u238 to th230 0.01  
   fParEfficiencyM1[7] = 0.7379278; // TeO2 Sx th230 only 0.01
   fParEfficiencyM1[8] = 0.5154646; // TeO2 Sx ra226 to pb210 0.01
   fParEfficiencyM1[9] = 0.4958653 ; // TeO2 Sx pb210 1 ==> necessary for bin below Po210 peak in M2
   fParEfficiencyM1[10] = 0.4457218; // TeO2 Sx pb210 0.01 ==> completely necessary for M2 spectrum
   fParEfficiencyM1[11] = 0.1013459; // CuBox+CuFrame Sx th232 10
   fParEfficiencyM1[12] = 0.0881936; // CuBox+CuFrame Sx u238 10
   fParEfficiencyM1[13] = 0.0725138; // CuBox+CuFrame Sx pb210 0.1 ==> useful for below Po210 peak in M1 but doesn't seem absolutely necessary
   fParEfficiencyM1[14] = 0.0742276; // CuBox+CuFrame Sx pb210 0.01 => necessary for below Po210 peak in M1
   fParEfficiencyM1[15] = 0.00484976; // PbRom k40
   fParEfficiencyM1[16] = 0.000465822; // OVC th232
   fParEfficiencyM1[17] = 0.000293556; // OVC u238
   fParEfficiencyM1[18] = 0.01430904; // OVC co60
   fParEfficiencyM1[19] = 0.00043995; // OVC k40
   fParEfficiencyM1[20] = 0.000191901; // External Lead bi210 
   fParEfficiencyM1[21] = 0.0510853; // CuBox+Frame th232 
   fParEfficiencyM1[22] = 0.0354736; // CuBox+Frame u238
   fParEfficiencyM1[23] = 0.0237369; // PbRom cs137

   fParMass[0] = 39000/1000.; // TeO2 2nu
   fParMass[1] = (2610.04+6929.71)/1000.; // CuBox+Frame co60
   fParMass[2] = 39000/1000.; // TeO2 th232 only
   fParMass[3] = 39000/1000.; // TeO2 th230 only
   fParMass[4] = 39000/1000.; // TeO2 Sx th232 0.01
   fParMass[5] = 39000/1000.; // TeO2 Sx ra228 to pb208 0.01 
   fParMass[6] = 39000/1000.; // TeO2 Sx u238 to th230 0.01  
   fParMass[7] = 39000/1000.; // TeO2 Sx th230 only 0.01
   fParMass[8] = 39000/1000.; // TeO2 Sx ra226 to pb210 0.01
   fParMass[9] = 39000/1000. ; // TeO2 Sx pb210 1 ==> necessary for bin below Po210 peak in M2
   fParMass[10] = 39000/1000.; // TeO2 Sx pb210 0.01 ==> completely necessary for M2 spectrum
   fParMass[11] = (2610.04+6929.71)/1000.; // CuBox+CuFrame Sx th232 10
   fParMass[12] = (2610.04+6929.71)/1000.; // CuBox+CuFrame Sx u238 10
   fParMass[13] = (2610.04+6929.71)/1000.; // CuBox+CuFrame Sx pb210 0.1 ==> useful for below Po210 peak in M1 but doesn't seem absolutely necessary
   fParMass[14] = (2610.04+6929.71)/1000.; // CuBox+CuFrame Sx pb210 0.01 => necessary for below Po210 peak in M1
   fParMass[15] = 202294.46/1000.; // PbRom k40
   fParMass[16] = 180704.38/1000.; // OVC th232
   fParMass[17] = 180704.38/1000.; // OVC u238
   fParMass[18] = 180704.38/1000.; // OVC co60
   fParMass[19] = 180704.38/1000.; // OVC k40
   fParMass[20] = 24652026/1000.; // External Lead bi210 
   fParMass[21] = (2610.04+6929.71)/1000.; // CuBox+Frame th232 
   fParMass[22] = (2610.04+6929.71)/1000.; // CuBox+Frame u238
   fParMass[23] = 202294.46/1000.; // PbRom cs137

   // fParSurfaceArea[0] = 7800.; // TeO2 0nu
   fParSurfaceArea[0] = 7800.; // TeO2 2nu
   fParSurfaceArea[1] = 2314.02+9467.18; // CuBox+Frame co60
   fParSurfaceArea[2] = 7800.; // TeO2 th232 only
   fParSurfaceArea[3] = 7800.; // TeO2 th230 only
   fParSurfaceArea[4] = 7800.; // TeO2 Sx th232 0.01
   fParSurfaceArea[5] = 7800.; // TeO2 Sx ra228 to pb208 0.01 
   fParSurfaceArea[6] = 7800.; // TeO2 Sx u238 to th230 0.01  
   fParSurfaceArea[7] = 7800.; // TeO2 Sx th230 only 0.01
   fParSurfaceArea[8] = 7800.; // TeO2 Sx ra226 to pb210 0.01
   fParSurfaceArea[9] = 7800. ; // TeO2 Sx pb210 1 ==> necessary for bin below Po210 peak in M2
   fParSurfaceArea[10] = 7800.; // TeO2 Sx pb210 0.01 ==> completely necessary for M2 spectrum
   fParSurfaceArea[11] = 2314.02+9467.18; // CuBox+CuFrame Sx th232 10
   fParSurfaceArea[12] = 2314.02+9467.18; // CuBox+CuFrame Sx u238 10
   fParSurfaceArea[13] = 2314.02+9467.18; // CuBox+CuFrame Sx pb210 0.1 ==> useful for below Po210 peak in M1 but doesn't seem absolutely necessary
   fParSurfaceArea[14] = 2314.02+9467.18; // CuBox+CuFrame Sx pb210 0.01 => necessary for below Po210 peak in M1
   fParSurfaceArea[15] = 20898.8; // PbRom k40
   fParSurfaceArea[16] = 87370.2; // OVC th232
   fParSurfaceArea[17] = 87370.2; // OVC u238
   fParSurfaceArea[18] = 87370.2; // OVC co60
   fParSurfaceArea[19] = 87370.2; // OVC k40
   fParSurfaceArea[20] = 2.38E+005; // External Lead bi210 
   fParSurfaceArea[21] = 2314.02+9467.18; // CuBox+Frame th232 
   fParSurfaceArea[22] = 2314.02+9467.18; // CuBox+Frame u238
   fParSurfaceArea[23] = 20898.8; // PbRom cs137
}

// Creates/updates the background model
void TBackgroundModel::UpdateModelAdaptive()
{
  if(fModelTotAdapM1 == NULL || fModelTotAdapM2 == NULL || fModelTotAdapM2Sum == NULL) 
  {
    cout << "Model Histogram Not Created" << endl;
    return;
  }
  // Reset all bins in model histograms in every loop
  fModelTotAdapM1->Reset();
  fModelTotAdapM2->Reset();
  // fModelTotAdapM2Sum->Reset();

  dNumCalls++;
  if(dNumCalls%1000==0)
  {
    cout << "Call #: "<< dNumCalls << endl;
  }
  // Create model
  for(int i = 0; i < dNParam; i++)
  {
    fModelTotAdapM1->Add( BkgPar[i]->GetHistM1(),              dDataIntegralM1*fParameters[i]);
    fModelTotAdapM2->Add( BkgPar[i]->GetHistM2(),              dDataIntegralM1*fParameters[i]);
    // fModelTotAdapM2Sum->Add( BkgPar[i]->GetHistM2Sum(),        dDataIntegralM1*fParameters[i]);
  }
}
  
// This method sets up minuit and does the fit
// fVariableValue is for contours
bool TBackgroundModel::DoTheFitAdaptive()
{ 
  // Reset initial parameter/error values
  ResetParameters();

  gStyle->SetOptStat(0);
  // gStyle->SetOptTitle(0);
  gStyle->SetOptFit();

  // Reduce Minuit Output
  minuit->SetPrintLevel(0); // Print level -1 (Quiet), 0 (Normal), 1 (Verbose)
  minuit->Command("SET STRategy 2"); // Sets strategy of fit
  minuit->SetMaxIterations(10000);
  minuit->SetObjectFit(this); //see the external FCN  above

  // Create and fix parameters
  for(int i = 0; i < dNParam; i++)
  {
    minuit->DefineParameter(BkgPar[i]->GetParIndex(), BkgPar[i]->GetParName(), BkgPar[i]->GetParInitial(), BkgPar[i]->GetParInitErr(), BkgPar[i]->GetParMin(), BkgPar[i]->GetParMax());
    if(bFixedArray[i]) 
    {
      minuit->FixParameter(i);
    // For debugging purposes only
    // cout << "Parameter " << fParameterName[i] << " is fixed" << endl; 
    }
  }

   // Number of free Parameters (for Chi-squared/NDF calculation only)
   dNumFreeParameters = minuit->GetNumPars() - minuit->GetNumFixedPars();
   dNDF = (dFitMaxBinM1+dFitMaxBinM2-dFitMinBinM1-dFitMinBinM2-dNumFreeParameters); // M1 and M2
   // Tell minuit what external function to use 
   minuit->SetFCN(myExternal_FCNAdap);
   int status = minuit->Command("MINImize 500000 0.1"); // Command that actually does the minimization

  // Get final parameters from fit
  for(int i = 0; i < dNParam; i++)
  {
    minuit->GetParameter(i, fParameters[i], fParError[i]);
  }

  // Update model with final parameters
  UpdateModelAdaptive();
  
  dChiSquare = GetChiSquareAdaptive();

  // cout << "Total number of calls = " << dNumCalls << "\t" << "ChiSq/NDF = " << dChiSquare/dNDF << endl; // for M1 and M2
  // cout << "ChiSq = " << dChiSquare << "\t" << "NDF = " << dNDF << endl;
  // cout << "Probability = " << TMath::Prob(dChiSquare, dNDF ) << endl;

  // Integrals
  for(int i = 0; i < dNParam; i++)
  {
    BkgPar[i]->GetHistM1()->Scale( dDataIntegralM1*fParameters[i]);
    BkgPar[i]->GetHistM2()->Scale( dDataIntegralM2*fParameters[i]);
    fParActivityM1[i] = BkgPar[i]->GetHistM1()->Integral("width")/dLivetimeYr;
    fParActivityM2[i] = BkgPar[i]->GetHistM2()->Integral("width")/dLivetimeYr;    
  }

  // ///// Draw Data M1
  fAdapDataHistoM1->SetLineColor(kBlack);
  fAdapDataHistoM1->SetLineWidth(2);
  // fAdapDataHistoM1->SetMaximum(9000);
  fAdapDataHistoM1->SetMinimum(0.02);  
  // hAdapTeO2k40M1->SetLineColor(3);
  // hAdapOVCk40M1->SetLineColor(4);
  // hAdapPbRomk40M1->SetLineColor(6);
  // hAdapInternalk40M1->SetLineColor(7);

  // Draw out the Alpha
  // TeO2 Surface
  // hAdapTeO2Sxu238M1_10->SetLineColor(5);
  // hAdapTeO2Sxu238M1_1->SetLineColor(5);

  // hAdapTeO2Sxra228pb208M1_001->SetLineColor(6);
  // hAdapTeO2Sxu238th230M1_001->SetLineColor(6);
  // hAdapTeO2Sxth230onlyM1_001->SetLineColor(6);
  // hAdapTeO2Sxra226pb210M1_001->SetLineColor(6);

  // hAdapTeO2Sxth232M1_10->SetLineColor(7);
  // hAdapTeO2Sxth232M1_1->SetLineColor(7);
  // hAdapTeO2Sxth232M1_001->SetLineColor(7);

  // hAdapTeO2po210M1->SetLineColor(3);
  hAdapTeO2Sxpb210M1_10->SetLineColor(3);
  hAdapTeO2Sxpb210M1_1->SetLineColor(4);
  hAdapTeO2Sxpb210M1_01->SetLineColor(5);

  // Cu Holder Surface
  // hAdapCuBox_CuFrameu238M1_100->SetLineColor(8);
  // hAdapCuBox_CuFrameu238M1_10->SetLineColor(8);
  // hAdapCuBox_CuFrameu238M1_1->SetLineColor(8);
  // hAdapCuBox_CuFrameu238M1_01->SetLineColor(8);
  // hAdapCuBox_CuFrameu238M1_001->SetLineColor(8);

  hAdapCuBox_CuFramepb210M1_100->SetLineColor(6);
  hAdapCuBox_CuFramepb210M1_10->SetLineColor(7);
  hAdapCuBox_CuFramepb210M1_1->SetLineColor(8);
  hAdapCuBox_CuFramepb210M1_01->SetLineColor(9);
  hAdapCuBox_CuFramepb210M1_001->SetLineColor(40);

  // hAdapCuBox_CuFrameth232M1_100->SetLineColor(46);  
  // hAdapCuBox_CuFrameth232M1_10->SetLineColor(46);
  // hAdapCuBox_CuFrameth232M1_1->SetLineColor(46);
  // hAdapCuBox_CuFrameth232M1_01->SetLineColor(46);
  // hAdapCuBox_CuFrameth232M1_001->SetLineColor(46);

  fModelTotAdapM1->SetLineColor(2);
  fModelTotAdapM1->SetLineWidth(1);


  TLegend *legfit1 = new TLegend(0.75,0.58,0.95,0.92);
  legfit1->SetFillColor(0);
  legfit1->SetTextSize(0.02);
  legfit1->AddEntry(fModelTotAdapM1, "Total PDF", "l");
  legfit1->AddEntry(hAdapTeO2Sxpb210M1_10, "TeO2 Pb210 Sx 10 #mum", "l" );
  legfit1->AddEntry(hAdapTeO2Sxpb210M1_1, "TeO2 Pb210 Sx 1 #mum", "l" );
  legfit1->AddEntry(hAdapTeO2Sxpb210M1_01, "TeO2 Pb210 Sx 0.1 #mum", "l" );
  legfit1->AddEntry(hAdapCuBox_CuFramepb210M1_100, "Copper Pb210 Sx 100 #mum", "l" );
  legfit1->AddEntry(hAdapCuBox_CuFramepb210M1_10, "Copper Pb210 Sx 10 #mum", "l" );
  legfit1->AddEntry(hAdapCuBox_CuFramepb210M1_1, "Copper Pb210 Sx 1 #mum", "l" );
  legfit1->AddEntry(hAdapCuBox_CuFramepb210M1_01, "Copper Pb210 Sx 0.1 #mum", "l" );
  legfit1->AddEntry(hAdapCuBox_CuFramepb210M1_001, "Copper Pb210 Sx 0.01 #mum", "l" );


  int nHisto = 2;  
  double width1 = 0.02;
  double width2 = 0.98;
  double canBotMargin = 0.02;
  double canTopMargin = 0.02;
  double padHeight = (1.-canTopMargin-canBotMargin)/nHisto;


  TCanvas *cadap1 = new TCanvas("cadap1", "cadap1", 1200, 800);
  cadap1->SetLogy();
  
  TPad* p1m1 = new TPad("p1m1","p1m1",width1,canBotMargin,width2,canBotMargin+padHeight/1.5,0,0);
  p1m1->SetTopMargin(0.);
  p1m1->SetBottomMargin(0.175);
  p1m1->SetLeftMargin(0.1);
  p1m1->SetRightMargin(0.05);
  p1m1->SetFillColor(0);
  p1m1->SetBorderMode(0);
  p1m1->SetBorderSize(0);
  p1m1->Draw();
  
  // p2 is on top!
  TPad* p2m1 = new TPad("p2m1","p2m1",width1,canBotMargin+padHeight/1.5,width2,canBotMargin+2*padHeight,0,0);
  p2m1->SetBottomMargin(0.);
  p2m1->SetTopMargin(0.08);
  p2m1->SetLeftMargin(0.1);
  p2m1->SetRightMargin(0.05);
  p2m1->SetFillColor(0);
  p2m1->SetBorderMode(0);
  p2m1->SetBorderSize(0);
  p2m1->Draw();


  p1m1->cd();
  TH1D *hRatioM1_1 = (TH1D*)fAdapDataHistoM1->Clone("hRatioM1_1");
  TH1D *hRatioM1_2 = (TH1D*)fAdapDataHistoM1->Clone("hRatioM1_2");
  TH1D *hRatioM1_3 = (TH1D*)fAdapDataHistoM1->Clone("hRatioM1_3");

  // hRatioM1_1->Add(fModelTotAdapM1, -1);
  hRatioM1_1->Divide(fModelTotAdapM1);
  hRatioM1_2->Divide(fModelTotAdapM1);
  hRatioM1_3->Divide(fModelTotAdapM1);
  hRatioM1_3->SetMaximum(2.9);
  hRatioM1_3->SetMinimum(-0.9);
  for(int i = 1; i <= hRatioM1_1->GetNbinsX(); i++)
  {
    hRatioM1_1->SetBinError(i, fAdapDataHistoM1->GetBinError(i)/fModelTotAdapM1->GetBinContent(i) );
    hRatioM1_2->SetBinError(i, 2*fAdapDataHistoM1->GetBinError(i)/fModelTotAdapM1->GetBinContent(i) );
    hRatioM1_3->SetBinError(i, 3*fAdapDataHistoM1->GetBinError(i)/fModelTotAdapM1->GetBinContent(i) );

  }
  TLine *LineM1 = new TLine();
  hRatioM1_3->GetXaxis()->SetRange(fAdapDataHistoM1->FindBin(dFitMin), fAdapDataHistoM1->FindBin(dFitMax-1));  
  // hRatioM1_3->GetXaxis()->SetRange(fAdapDataHistoM1->FindBin(1440), fAdapDataHistoM1->FindBin(1480));
  hRatioM1_3->SetMarkerStyle(6);
  hRatioM1_3->SetTitle("");
  hRatioM1_3->GetXaxis()->SetTitle("Energy (keV)");
  hRatioM1_3->GetYaxis()->SetTitle("Ratio (Data/MC)");
  hRatioM1_3->GetXaxis()->SetLabelSize(0.07);
  hRatioM1_3->GetYaxis()->SetLabelSize(0.07);
  hRatioM1_3->GetXaxis()->SetTitleSize(0.07);
  hRatioM1_3->GetYaxis()->SetTitleSize(0.07);
  hRatioM1_3->GetYaxis()->SetTitleOffset(0.45);
  hRatioM1_1->SetFillColor(kMagenta-9);
  hRatioM1_2->SetFillColor(kGreen-8);
  hRatioM1_3->SetFillColor(kCyan+3);
  hRatioM1_3->Draw("pe2");
  hRatioM1_2->Draw("SAME e2");
  hRatioM1_1->Draw("SAME e2");
  // LineM1->DrawLine(dFitMin, 1, dFitMax-1, 1);
  LineM1->DrawLineNDC(0.1, 0.59, 0.95, 0.59);

  p2m1->cd();
  p2m1->SetLogy();
  
  fAdapDataHistoM1->GetXaxis()->SetRange(fAdapDataHistoM1->FindBin(dFitMin), fAdapDataHistoM1->FindBin(dFitMax-1));
  // fAdapDataHistoM1->GetXaxis()->SetRange(fAdapDataHistoM1->FindBin(1), fAdapDataHistoM1->FindBin(dFitMax-1));
  // fAdapDataHistoM1->GetXaxis()->SetRange(fAdapDataHistoM1->FindBin(1440), fAdapDataHistoM1->FindBin(1480));  
  fAdapDataHistoM1->SetTitle("Total Model (M1)");
  // fAdapDataHistoM1->SetTitleOffset(1.5);
  // fAdapDataHistoM1->SetTitleSize(0.005);
  fAdapDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM1->GetYaxis()->SetTitle("Counts/Bin");
  fAdapDataHistoM1->Draw("E");
  fModelTotAdapM1->Draw("SAME");


  //// Draw additional histograms here
  // hAdapTeO2k40M1->Draw("SAME");
  // hAdapOVCk40M1->Draw("SAME");
  // hAdapInternalk40M1->Draw("SAME");
  // hAdapPbRomk40M1->Draw("SAME");  
  // legfit1->Draw();
  // hAdapCuFramek40M1->Draw("SAME");
  // TeO2 Surface
  // hAdapTeO2Sxu238M1_10->Draw("SAME");
  // hAdapTeO2Sxu238M1_1->Draw("SAME");

  // hAdapTeO2Sxra228pb208M1_001->Draw("SAME");
  // hAdapTeO2Sxu238th230M1_001->Draw("SAME");
  // hAdapTeO2Sxth230onlyM1_001->Draw("SAME");
  // hAdapTeO2Sxra226pb210M1_001->Draw("SAME");

  // hAdapTeO2Sxth232M1_10->Draw("SAME");
  // hAdapTeO2Sxth232M1_1->Draw("SAME");
  // hAdapTeO2Sxth232M1_001->Draw("SAME");


  // hAdapTeO2po210M1->Draw("SAME");
  // hAdapTeO2Sxpb210M1_10->Draw("SAME");
  // hAdapTeO2Sxpb210M1_1->Draw("SAME");
  // hAdapTeO2Sxpb210M1_01->Draw("SAME");

  // Cu Holder Surface
  // hAdapCuBox_CuFrameu238M1_100->Draw("SAME");
  // hAdapCuBox_CuFrameu238M1_10->Draw("SAME");
  // hAdapCuBox_CuFrameu238M1_1->Draw("SAME");
  // hAdapCuBox_CuFrameu238M1_01->Draw("SAME");
  // hAdapCuBox_CuFrameu238M1_001->Draw("SAME");

  // hAdapCuBox_CuFramepb210M1_100->Draw("SAME");
  // hAdapCuBox_CuFramepb210M1_10->Draw("SAME");
  // hAdapCuBox_CuFramepb210M1_1->Draw("SAME");
  // hAdapCuBox_CuFramepb210M1_01->Draw("SAME");
  // hAdapCuBox_CuFramepb210M1_001->Draw("SAME");

  // hAdapCuBox_CuFrameth232M1_100->Draw("SAME");
  // hAdapCuBox_CuFrameth232M1_10->Draw("SAME");
  // hAdapCuBox_CuFrameth232M1_1->Draw("SAME");
  // hAdapCuBox_CuFrameth232M1_01->Draw("SAME");
  // hAdapCuBox_CuFrameth232M1_001->Draw("SAME");

  ///// Draw Data M2
  fAdapDataHistoM2->SetLineColor(kBlack);
  fAdapDataHistoM2->SetLineWidth(2);
  // fAdapDataHistoM2->SetMinimum(0.002);
  // fAdapDataHistoM2->SetMaximum(9000);
  // TeO2 Surface
  // hAdapTeO2Sxu238M2_10->SetLineColor(5);
  // hAdapTeO2Sxu238M2_1->SetLineColor(5);

  // hAdapTeO2Sxra228pb208M2_001->SetLineColor(6);
  // hAdapTeO2Sxu238th230M2_001->SetLineColor(6);
  // hAdapTeO2Sxth230onlyM2_001->SetLineColor(6);
  // hAdapTeO2Sxra226pb210M2_001->SetLineColor(6);

  // hAdapTeO2Sxth232M2_10->SetLineColor(7);
  // hAdapTeO2Sxth232M2_1->SetLineColor(7);
  // hAdapTeO2Sxth232M2_001->SetLineColor(7);


  // hAdapTeO2po210M2->SetLineColor(3);
  hAdapTeO2Sxpb210M2_10->SetLineColor(3);
  hAdapTeO2Sxpb210M2_1->SetLineColor(4);
  hAdapTeO2Sxpb210M2_01->SetLineColor(5);

  // Cu Holder Surface
  // hAdapCuBox_CuFrameu238M2_100->SetLineColor(8);
  // hAdapCuBox_CuFrameu238M2_10->SetLineColor(8);
  // hAdapCuBox_CuFrameu238M2_1->SetLineColor(8);
  // hAdapCuBox_CuFrameu238M2_01->SetLineColor(8);
  // hAdapCuBox_CuFrameu238M2_001->SetLineColor(8);

  hAdapCuBox_CuFramepb210M2_100->SetLineColor(6);
  hAdapCuBox_CuFramepb210M2_10->SetLineColor(7);
  hAdapCuBox_CuFramepb210M2_1->SetLineColor(8);
  hAdapCuBox_CuFramepb210M2_01->SetLineColor(9);
  hAdapCuBox_CuFramepb210M2_001->SetLineColor(40);

  // hAdapCuBox_CuFrameth232M2_100->SetLineColor(46);  
  // hAdapCuBox_CuFrameth232M2_10->SetLineColor(46);
  // hAdapCuBox_CuFrameth232M2_1->SetLineColor(46);
  // hAdapCuBox_CuFrameth232M2_01->SetLineColor(46);
  // hAdapCuBox_CuFrameth232M2_001->SetLineColor(46);

  fModelTotAdapM2->SetLineColor(2);
  fModelTotAdapM2->SetLineWidth(1);

  TLegend *legfit2 = new TLegend(0.75,0.6,0.95,0.92);
  legfit2->SetFillColor(0);
  legfit2->SetTextSize(0.02);
  legfit2->AddEntry(fModelTotAdapM2, "Total PDF", "l");
  legfit2->AddEntry(hAdapTeO2Sxpb210M2_10, "TeO2 Pb210 Sx 10 #mum", "l" );
  legfit2->AddEntry(hAdapTeO2Sxpb210M2_1, "TeO2 Pb210 Sx 1 #mum", "l" );
  legfit2->AddEntry(hAdapTeO2Sxpb210M2_01, "TeO2 Pb210 Sx 0.1 #mum", "l" );
  legfit2->AddEntry(hAdapCuBox_CuFramepb210M2_100, "Copper Pb210 Sx 100 #mum", "l" );
  legfit2->AddEntry(hAdapCuBox_CuFramepb210M2_10, "Copper Pb210 Sx 10 #mum", "l" );
  legfit2->AddEntry(hAdapCuBox_CuFramepb210M2_1, "Copper Pb210 Sx 1 #mum", "l" );
  legfit2->AddEntry(hAdapCuBox_CuFramepb210M2_01, "Copper Pb210 Sx 0.1 #mum", "l" );
  legfit2->AddEntry(hAdapCuBox_CuFramepb210M2_001, "Copper Pb210 Sx 0.01 #mum", "l" );

  TCanvas *cadap2 = new TCanvas("cadap2", "cadap2", 1200, 800);  
  TPad* p1m2 = new TPad("p1m2","p1m2",width1,canBotMargin,width2,canBotMargin+padHeight/1.5,0,0);
  p1m2->SetTopMargin(0.);
  p1m2->SetBottomMargin(0.175);
  p1m2->SetLeftMargin(0.1);
  p1m2->SetRightMargin(0.05);
  p1m2->SetFillColor(0);
  p1m2->SetBorderMode(0);
  p1m2->SetBorderSize(0);
  p1m2->Draw();
  
  // p2 is on top!
  TPad* p2m2 = new TPad("p2m2","p2m2",width1,canBotMargin+padHeight/1.5,width2,canBotMargin+2*padHeight,0,0);
  p2m2->SetBottomMargin(0.);
  p2m2->SetTopMargin(0.08);
  p2m2->SetLeftMargin(0.1);
  p2m2->SetRightMargin(0.05);
  p2m2->SetFillColor(0);
  p2m2->SetBorderMode(0);
  p2m2->SetBorderSize(0);
  p2m2->Draw();


  p1m2->cd();
  TH1D *hRatioM2_1 = (TH1D*)fAdapDataHistoM2->Clone("hRatioM2_1");
  TH1D *hRatioM2_2 = (TH1D*)fAdapDataHistoM2->Clone("hRatioM2_2");
  TH1D *hRatioM2_3 = (TH1D*)fAdapDataHistoM2->Clone("hRatioM2_3");
  hRatioM2_1->Divide(fModelTotAdapM2);
  hRatioM2_2->Divide(fModelTotAdapM2);
  hRatioM2_3->Divide(fModelTotAdapM2);
  hRatioM2_3->SetMaximum(2.9);
  hRatioM2_3->SetMinimum(-0.9);
  for(int i = 1; i <= hRatioM2_1->GetNbinsX(); i++)
  {
    hRatioM2_1->SetBinError(i, fAdapDataHistoM2->GetBinError(i)/fModelTotAdapM2->GetBinContent(i) );
    hRatioM2_2->SetBinError(i, 2*fAdapDataHistoM2->GetBinError(i)/fModelTotAdapM2->GetBinContent(i) );
    hRatioM2_3->SetBinError(i, 3*fAdapDataHistoM2->GetBinError(i)/fModelTotAdapM2->GetBinContent(i) );
  }
  // hRatioM2_3->GetXaxis()->SetRange(fAdapDataHistoM2->FindBin(1400), fAdapDataHistoM2->FindBin(1480));
  hRatioM2_3->GetXaxis()->SetRange(fAdapDataHistoM2->FindBin(dFitMin), fAdapDataHistoM2->FindBin(dFitMax-1));
  hRatioM2_3->SetMarkerStyle(6);
  hRatioM2_3->SetTitle("");
  hRatioM2_3->GetXaxis()->SetTitle("Energy (keV)");
  hRatioM2_3->GetYaxis()->SetTitle("Ratio (Data/MC)");  
  hRatioM2_3->GetXaxis()->SetLabelSize(0.07);
  hRatioM2_3->GetYaxis()->SetLabelSize(0.07);
  hRatioM2_3->GetXaxis()->SetTitleSize(0.07);
  hRatioM2_3->GetYaxis()->SetTitleSize(0.07);
  hRatioM2_3->GetYaxis()->SetTitleOffset(0.45);
  hRatioM2_1->SetFillColor(kMagenta-9);
  hRatioM2_2->SetFillColor(kGreen-8);
  hRatioM2_3->SetFillColor(kCyan+3);
  hRatioM2_3->Draw("pE2");
  hRatioM2_2->Draw("SAME e2");
  hRatioM2_1->Draw("SAME e2");
  // LineM1->DrawLine(dFitMin, 1, dFitMax-1, 1);
  LineM1->DrawLineNDC(0.1, 0.59, 0.95, 0.59);

  p2m2->cd();
  p2m2->SetLogy();
  // fAdapDataHistoM2->GetXaxis()->SetRange(fAdapDataHistoM2->FindBin(1440), fAdapDataHistoM2->FindBin(1480));
  fAdapDataHistoM2->GetXaxis()->SetRange(fAdapDataHistoM2->FindBin(dFitMin), fAdapDataHistoM2->FindBin(dFitMax-1));
  // fAdapDataHistoM2->GetXaxis()->SetRange(fAdapDataHistoM2->FindBin(1), fAdapDataHistoM2->FindBin(dFitMax-1));
  // fAdapDataHistoM2->SetTitleOffset(0.45);
  // fAdapDataHistoM2->SetTitleSize(0.01);
  fAdapDataHistoM2->SetTitle("Total Model (M2)");
  fAdapDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM2->GetYaxis()->SetTitle("Counts/Bin");
  fAdapDataHistoM2->Draw("E");
  fModelTotAdapM2->Draw("SAME");
  // legfit2->Draw();
  //// Draw additional histograms
  // TeO2 Surface
  // hAdapTeO2Sxu238M2_10->Draw("SAME");
  // hAdapTeO2Sxu238M2_1->Draw("SAME");

  // hAdapTeO2Sxra228pb208M2_001->Draw("SAME");
  // hAdapTeO2Sxu238th230M2_001->Draw("SAME");
  // hAdapTeO2Sxth230onlyM2_001->Draw("SAME");
  // hAdapTeO2Sxra226pb210M2_001->Draw("SAME");

  // hAdapTeO2Sxth232M2_10->Draw("SAME");
  // hAdapTeO2Sxth232M2_1->Draw("SAME");
  // hAdapTeO2Sxth232M2_001->Draw("SAME");


  // hAdapTeO2po210M2->Draw("SAME");
  // hAdapTeO2Sxpb210M2_10->Draw("SAME");
  // hAdapTeO2Sxpb210M2_1->Draw("SAME");
  // hAdapTeO2Sxpb210M2_01->Draw("SAME");

  // Cu Holder Surface
  // hAdapCuBox_CuFrameu238M2_100->Draw("SAME");
  // hAdapCuBox_CuFrameu238M2_10->Draw("SAME");
  // hAdapCuBox_CuFrameu238M2_1->Draw("SAME");
  // hAdapCuBox_CuFrameu238M2_01->Draw("SAME");
  // hAdapCuBox_CuFrameu238M2_001->Draw("SAME");

  // hAdapCuBox_CuFramepb210M2_100->Draw("SAME");
  // hAdapCuBox_CuFramepb210M2_10->Draw("SAME");
  // hAdapCuBox_CuFramepb210M2_1->Draw("SAME");
  // hAdapCuBox_CuFramepb210M2_01->Draw("SAME");
  // hAdapCuBox_CuFramepb210M2_001->Draw("SAME");

  // hAdapCuBox_CuFrameth232M2_100->Draw("SAME");
  // hAdapCuBox_CuFrameth232M2_10->Draw("SAME");
  // hAdapCuBox_CuFrameth232M2_1->Draw("SAME");
  // hAdapCuBox_CuFrameth232M2_01->Draw("SAME");
  // hAdapCuBox_CuFrameth232M2_001->Draw("SAME");

/*
  ///// Draw Data M2Sum
  fAdapDataHistoM2Sum->SetLineColor(kBlack);
  fAdapDataHistoM2Sum->SetLineWidth(2);
  // fAdapDataHistoM2Sum->SetMaximum(9000);
  // fAdapDataHistoM2Sum->GetXaxis()->SetRange(1, fAdapDataHistoM2Sum->FindBin(3000));
  
  fModelTotAdapM2Sum->SetLineColor(2);
  fModelTotAdapM2Sum->SetLineWidth(1);

  TCanvas *cadap2sum = new TCanvas("cadap2sum", "cadap2sum", 1200, 800);
  TPad* p1m2sum = new TPad("p1m2sum","p1m2sum",width1,canBotMargin,width2,canBotMargin+padHeight/1.5,0,0);
  p1m2sum->SetTopMargin(0.);
  p1m2sum->SetBottomMargin(0.175);
  p1m2sum->SetLeftMargin(0.1);
  p1m2sum->SetRightMargin(0.05);
  p1m2sum->SetFillColor(0);
  p1m2sum->SetBorderMode(0);
  p1m2sum->SetBorderSize(0);
  p1m2sum->Draw();
  
  // p2 is on top!
  TPad* p2m2sum = new TPad("p2m2sum","p2m2sum",width1,canBotMargin+padHeight/1.5,width2,canBotMargin+2*padHeight,0,0);
  p2m2sum->SetBottomMargin(0.);
  p2m2sum->SetTopMargin(0.08);
  p2m2sum->SetLeftMargin(0.1);
  p2m2sum->SetRightMargin(0.05);
  p2m2sum->SetFillColor(0);
  p2m2sum->SetBorderMode(0);
  p2m2sum->SetBorderSize(0);
  p2m2sum->Draw();

  p1m2sum->cd();
  TH1D *hRatioM2Sum_1 = (TH1D*)fAdapDataHistoM2Sum->Clone("hRatioM2Sum_1");
  TH1D *hRatioM2Sum_2 = (TH1D*)fAdapDataHistoM2Sum->Clone("hRatioM2Sum_2");
  TH1D *hRatioM2Sum_3 = (TH1D*)fAdapDataHistoM2Sum->Clone("hRatioM2Sum_3");
  hRatioM2Sum_1->Divide(fModelTotAdapM2Sum);
  hRatioM2Sum_2->Divide(fModelTotAdapM2Sum);
  hRatioM2Sum_3->Divide(fModelTotAdapM2Sum);
  hRatioM2Sum_3->SetMaximum(2.9);
  hRatioM2Sum_3->SetMinimum(-0.9);
  for(int i = 1; i <= hRatioM2Sum_1->GetNbinsX(); i++)
  {
    hRatioM2Sum_1->SetBinError(i, fAdapDataHistoM2Sum->GetBinError(i)/fModelTotAdapM2Sum->GetBinContent(i) );
    hRatioM2Sum_2->SetBinError(i, 2*fAdapDataHistoM2Sum->GetBinError(i)/fModelTotAdapM2Sum->GetBinContent(i) );
    hRatioM2Sum_3->SetBinError(i, 3*fAdapDataHistoM2Sum->GetBinError(i)/fModelTotAdapM2Sum->GetBinContent(i) );
  }
  hRatioM2Sum_3->GetXaxis()->SetRange(fAdapDataHistoM2Sum->FindBin(dFitMin), fAdapDataHistoM2Sum->FindBin(dFitMax-1));
  hRatioM2Sum_3->SetMarkerStyle(6);
  hRatioM2Sum_3->SetTitle("");
  hRatioM2Sum_3->GetXaxis()->SetTitle("Energy (keV)");
  hRatioM2Sum_3->GetYaxis()->SetTitle("Ratio (Data/MC)");  
  hRatioM2Sum_3->GetXaxis()->SetLabelSize(0.07);
  hRatioM2Sum_3->GetYaxis()->SetLabelSize(0.07);
  hRatioM2Sum_3->GetXaxis()->SetTitleSize(0.07);
  hRatioM2Sum_3->GetYaxis()->SetTitleSize(0.07);
  hRatioM2Sum_3->GetYaxis()->SetTitleOffset(0.45);
  hRatioM2Sum_1->SetFillColor(kMagenta-9);
  hRatioM2Sum_2->SetFillColor(kGreen-8);
  hRatioM2Sum_3->SetFillColor(kCyan+3);
  hRatioM2Sum_3->Draw("pE2");
  hRatioM2Sum_2->Draw("SAME e2");
  hRatioM2Sum_1->Draw("SAME e2");
  LineM1->DrawLine(dFitMin, 1, dFitMax-1, 1);

  p2m2sum->cd();
  p2m2sum->SetLogy();
  fAdapDataHistoM2Sum->GetXaxis()->SetRange(fAdapDataHistoM2Sum->FindBin(dFitMin), fAdapDataHistoM2Sum->FindBin(dFitMax-1));
  // fAdapDataHistoM2Sum->SetTitleOffset(0.45);
  // fAdapDataHistoM2Sum->SetTitleSize(0.01);
  fAdapDataHistoM2Sum->SetTitle("Total Model (M2Sum)");
  fAdapDataHistoM2Sum->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM2Sum->GetYaxis()->SetTitle("Counts/Bin");
  fAdapDataHistoM2Sum->Draw("E");
  fModelTotAdapM2Sum->Draw("SAME");

*/




  TCanvas *CanvasM1 = new TCanvas("CanvasM1", "CanvasM1", 1200, 1000);
  CanvasM1->SetLogy();
  fAdapDataHistoM1->Draw("E");
  fModelTotAdapM1->Draw("SAME");


  TCanvas *CanvasM2 = new TCanvas("CanvasM2", "CanvasM2", 1200, 1000);
  CanvasM2->SetLogy();
  fAdapDataHistoM2->Draw("E");
  fModelTotAdapM2->Draw("SAME");







  // Residuals
  TCanvas *cResidual1 = new TCanvas("cResidual1", "cResidual1", 1200, 800);

  hResidualGausM1 = new TH1D("hResidualGausM1", "M1", 100, -50, 50);
  hResidualDistM1 = CalculateResidualsAdaptive(fAdapDataHistoM1, fModelTotAdapM1, hResidualGausM1, dFitMinBinM1, dFitMaxBinM1, 1);
  hResidualDistM1->SetLineColor(kBlack);
  hResidualDistM1->SetName("Residuals");
  hResidualDistM1->SetTitle("Normalized Residuals (M1)");
  hResidualDistM1->SetMarkerStyle(25);
  hResidualDistM1->GetXaxis()->SetTitle("Energy (keV)");
  // hResidualDistM1->GetXaxis()->SetTitleSize(0.04);
  // hResidualDistM1->GetXaxis()->SetLabelSize(0.05);
  // hResidualDistM1->GetYaxis()->SetLabelSize(0.05);
  // hResidualDistM1->GetYaxis()->SetTitleSize(0.04);
  // hResidualDistM1->GetXaxis()->SetRange(1, fAdapDataHistoM2->FindBin(3000));
  hResidualDistM1->GetYaxis()->SetTitle("Residuals (#sigma)");
  hResidualDistM1->Draw();

  TCanvas *cResidual2 = new TCanvas("cResidual2", "cResidual2", 1200, 800);

  hResidualGausM2 = new TH1D("hResidualGausM2", "M2", 100, -50, 50);  
  hResidualDistM2 = CalculateResidualsAdaptive(fAdapDataHistoM2, fModelTotAdapM2, hResidualGausM2, dFitMinBinM2, dFitMaxBinM2, 2);
  hResidualDistM2->SetLineColor(kBlack);
  hResidualDistM2->SetName("Residuals");
  hResidualDistM2->SetTitle("Normalized Residuals (M2)");
  hResidualDistM2->SetMarkerStyle(25);
  hResidualDistM2->GetXaxis()->SetTitle("Energy (keV)");
  // hResidualDistM2->GetXaxis()->SetTitleSize(0.04);
  // hResidualDistM2->GetXaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetTitleSize(0.04); 
  // hResidualDistM2->GetXaxis()->SetRange(1, fAdapDataHistoM2->FindBin(3000));
  hResidualDistM2->GetYaxis()->SetTitle("Residuals (#sigma)");
  hResidualDistM2->Draw();


/*
  TCanvas *cResidual2Sum = new TCanvas("cResidual2Sum", "cResidual2Sum", 1200, 800);

  hResidualGausM2Sum = new TH1D("hResidualGausM2Sum", "M2Sum", 100, -50, 50);  
  hResidualDistM2Sum = CalculateResidualsAdaptive(fAdapDataHistoM2Sum, fModelTotAdapM2Sum, hResidualGausM2Sum, dFitMinBinM2Sum, dFitMaxBinM2Sum, 3);
  hResidualDistM2Sum->SetLineColor(kBlack);
  hResidualDistM2Sum->SetName("Residuals");
  hResidualDistM2Sum->SetTitle("Normalized Residuals (M2Sum)");
  hResidualDistM2Sum->SetMarkerStyle(25);
  hResidualDistM2Sum->GetXaxis()->SetTitle("Energy (keV)");
  // hResidualDistM2Sum->GetXaxis()->SetTitleSize(0.04);
  // hResidualDistM2Sum->GetXaxis()->SetLabelSize(0.05);
  // hResidualDistM2Sum->GetYaxis()->SetLabelSize(0.05);
  // hResidualDistM2Sum->GetYaxis()->SetTitleSize(0.04); 
  // hResidualDistM2Sum->GetXaxis()->SetRange(1, fAdapDataHistoM2Sum->FindBin(3000));
  hResidualDistM2Sum->GetYaxis()->SetTitle("Residuals (#sigma)");
  hResidualDistM2Sum->Draw();
*/

  TCanvas *cres1 = new TCanvas("cres1", "cres1", 1600, 600);
  cres1->Divide(2,1);
  cres1->cd(1);

  // TCanvas *cres1 = new TCanvas("cres1", "cres1", 1200, 800);
  hResidualGausM1->Fit("gaus");
  hResidualGausM1->Draw();
  cres1->cd(2);
  // TCanvas *cres2 = new TCanvas("cres2", "cres2", 1200, 800);
  hResidualGausM2->Fit("gaus");
  hResidualGausM2->Draw();


  for (int i = dFitMinBinM1; i < dFitMaxBinM1; i++)
  {
    dResidualRMSM1 += hResidualDistM1->GetBinContent(i)*hResidualDistM1->GetBinContent(i);
  }
  for (int i = dFitMinBinM2; i < dFitMaxBinM2; i++)
  {
    dResidualRMSM2 += hResidualDistM2->GetBinContent(i)*hResidualDistM2->GetBinContent(i);
  }
  // for (int i = dFitMinBinM2Sum; i < dFitMaxBinM2Sum; i++)
  // {
    // dResidualRMSM2Sum += hResidualDistM2Sum->GetBinContent(i)*hResidualDistM2Sum->GetBinContent(i);
  // }

  dResidualRMSTot = TMath::Sqrt( (dResidualRMSM1 + dResidualRMSM2)/ (dNDF + dNumFreeParameters) );

  dResidualRMSM1 = TMath::Sqrt(dResidualRMSM1/(dFitMaxBinM1-dFitMinBinM1));
  dResidualRMSM2 = TMath::Sqrt(dResidualRMSM2/(dFitMaxBinM2-dFitMinBinM2));
  // dResidualRMSM2Sum = TMath::Sqrt(dResidualRMSM2Sum/(dFitMaxBinM2Sum-dFitMinBinM2Sum));


  //
  double dROIRange = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2570))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2570)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2470)); 
  double d2nbbRange = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2000))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2000)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(500));
  double dGammaRange = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(3000))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(3000)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(100));
  double dAlphaRangeM1 = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(8000))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(8000)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(3000));
  double dAlphaRangeM2 = fAdapDataHistoM2->GetBinLowEdge(fAdapDataHistoM2->FindBin(8000))+fAdapDataHistoM2->GetBinWidth(fAdapDataHistoM2->FindBin(8000)) - fAdapDataHistoM2->GetBinLowEdge(fAdapDataHistoM2->FindBin(3000));

  double d2nbbData = fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double d2nbbDataErr = TMath::Sqrt(fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);
  double d2nbbModel = fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double d2nbbModelErr = TMath::Sqrt(fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);
  double d2nbbPDF = hAdapTeO22nuM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double d2nbbPDFErr = TMath::Sqrt(hAdapTeO22nuM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);
  double dROIData = fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr);
  double dROIDataErr = TMath::Sqrt(fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr);
  double dROIModel = fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr);
  double dROIModelErr = TMath::Sqrt(fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr);
  
  double dGammaModel = fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(100), fAdapDataHistoM1->FindBin(3000), "width" )/(dGammaRange*dLivetimeYr);
  double dGammaModelErr = TMath::Sqrt(fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(100), fAdapDataHistoM1->FindBin(3000), "width" ))/(dGammaRange*dLivetimeYr);
  double dAlphaModelM1 = fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(3000), fAdapDataHistoM1->FindBin(8000), "width" )/(dAlphaRangeM1*dLivetimeYr);
  double dAlphaModelErrM1 = TMath::Sqrt(fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(3000), fAdapDataHistoM1->FindBin(8000), "width" ))/(dAlphaRangeM1*dLivetimeYr);
  double dAlphaModelM2 = fModelTotAdapM2->Integral( fAdapDataHistoM2->FindBin(3000), fAdapDataHistoM2->FindBin(8000), "width" )/(dAlphaRangeM2*dLivetimeYr);
  double dAlphaModelErrM2 = TMath::Sqrt(fModelTotAdapM2->Integral( fAdapDataHistoM2->FindBin(3000), fAdapDataHistoM2->FindBin(8000), "width" ))/(dAlphaRangeM2*dLivetimeYr);

  // Output integrals of stuff for limits
  cout << "ROI range: " << fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2470)) << " " << fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2570))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2570)) << " keV" << endl; // 2470 to 2572
  cout << "Integral Data in ROI: " << fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt( fAdapDataHistoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << "Integral Total PDF in ROI: " << fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width") << " +/- " << sqrt( fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << endl;
  cout << endl;
  cout << "Integral Data in ROI (counts/keV/y): " << fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " +/- " << sqrt( fAdapDataHistoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << endl;
  cout << "Integral Total PDF in ROI (counts/keV/y): " << fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " +/- " << sqrt( fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << endl;
  // cout << "Integral Total Po-210 PDF in ROI (counts/keV): " << fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " +/- " << sqrt(fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << endl;  
  // cout << "Integral Total 0NDBD PDF in ROI (counts/keV/y): " << fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " +/- " << sqrt(fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << endl;
  // cout << "Number of 2nbb: " << fParameters[0]*dDataIntegralM1 << " +/- " << fParError[0]*dDataIntegralM1 << "\t 2nbb half life: " << (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1) << " +/- " << (fParError[0]/fParameters[0]) * (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1) << endl;
  cout << "Counts in 2nbb (M1 + M2): " << hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width") + hAdapTeO22nuM2->Integral(1, fAdapDataHistoM2->FindBin(2700) , "width")/2 << "\t Half-Life " << (0.69314718056)*(4.726e25 * dLivetimeYr)/(hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width") + hAdapTeO22nuM2->Integral(1, fAdapDataHistoM2->FindBin(2700) , "width")/2) << endl;
  cout << endl;
  cout << endl;
  cout << "Data in 2nbb region (c/keV/y): " << d2nbbData << " $\\pm$ " << d2nbbDataErr << endl;  
  cout << "Model in 2nbb region (c/keV/y): " << d2nbbModel << " $\\pm$ " << d2nbbModelErr << endl;
  cout << "Model of 2nbb PDF (c/keV/y): " << d2nbbPDF << " $\\pm$ " << d2nbbPDFErr << endl;
  cout << "Model of PDFs in 2nbb region without 2nbb (c/keV/y): " << d2nbbModel - d2nbbPDF << " $\\pm$ " << d2nbbModelErr << endl;
  cout << "Data in 0nbb region (c/keV/y): " << dROIData << " $\\pm$ " << dROIDataErr << endl;
  cout << "Model in 0nbb region (c/keV/y): " << dROIModel << " $\\pm$ " << dROIModelErr << endl;
  cout << "Model/Data in 2nbb region: " << d2nbbModel/d2nbbData << endl;
  cout << "Model/Data in 0nbb region: " << dROIModel/dROIData << endl;
  cout << "Model of PDFs in 2nbb region without 2nbb/Data : " << (d2nbbModel - d2nbbPDF)/d2nbbData << endl;
  cout << "Model in Gamma region (c/keV/y): " << dGammaModel << " $\\pm$ " << dGammaModelErr << endl;
  cout << "Model in Alpha region M1 (c/keV/y): " << dAlphaModelM1 << " $\\pm$ " << dAlphaModelErrM1 << endl;
  cout << "Model in Alpha region M2 (c/keV/y): " << dAlphaModelM2 << " $\\pm$ " << dAlphaModelErrM2 << endl;

  cout << endl;  
  cout << endl;  
  cout << "Residual RMS (Tot): " << dResidualRMSTot << endl;
  cout << "Residual RMS (M1): " << dResidualRMSM1 << "\t" << "Residual RMS (M2): " << dResidualRMSM2 << "\t Residual RMS (M2Sum): "  << dResidualRMSM2Sum << endl;

  dChiSquare = GetChiSquareAdaptive();

  cout << "Total number of calls = " << dNumCalls << "\t" << "ChiSq/NDF = " << dChiSquare/dNDF << endl; // for M1 and M2
  cout << "ChiSq = " << dChiSquare << "\t" << "NumFreeParamters = " << dNumFreeParameters << "\t" << "NDF = " << dNDF << endl;
  cout << "Probability = " << TMath::Prob(dChiSquare, dNDF) << endl;

/*
    2nbb calculation:
     - TeO2 molar mass: 159.6 g/mol
     - half life is 9.81 * 10^20 years
     - how many in Q0 data so far? 1/rate = half life/ln(2) -> rate = ln(2)/half life = 7.066*10^-22 decays/year (Laura's thesis)
     - Moles = 750g * 49 crystals * 0.3408 abundance/159.6 g/mol = 78.474 mol total
     - N_TeO2 = 78.474 * N_A = 4.726*10^25 nuclei of Te130
     - N_2nbb = N_TeO2 * rate * livetime = 1.551*10^4 events
     - half life = rate * ln(2) = ln(2) * N_TeO2 * livetime / N_2nbb
*/


  // Correlation Matrix section
  TMatrixT<double> mCorrMatrix;
  mCorrMatrix.ResizeTo(TBackgroundModel::dNParam, TBackgroundModel::dNParam);
  minuit->mnemat(mCorrMatrix.GetMatrixArray(), TBackgroundModel::dNParam);

  TMatrixT<double> mReducedMatrix;
  mReducedMatrix.ResizeTo(dNumFreeParameters, dNumFreeParameters);

  for(int i = mCorrMatrix.GetRowLwb(); i <= mCorrMatrix.GetRowUpb(); i++)
    for(int j = mCorrMatrix.GetColLwb(); j <= mCorrMatrix.GetColUpb(); j++)
      mCorrMatrix(i,j) = mCorrMatrix(i,j)/(fParError[i]*fParError[j]);

  for(int i = mReducedMatrix.GetRowLwb(); i <= mReducedMatrix.GetRowUpb(); i++)
    for(int j = mReducedMatrix.GetColLwb(); j <= mReducedMatrix.GetColUpb(); j++)
      mReducedMatrix(i,j) = mCorrMatrix(i,j);


  // for(int i = mCorrMatrix.GetRowLwb(); i <= mCorrMatrix.GetRowUpb(); i++)
    // for(int j = mCorrMatrix.GetColLwb(); j <= mCorrMatrix.GetColUpb(); j++)
      // mCorrMatrixInverse(i,j) = mCorrMatrix(TBackgroundModel::dNParam-i-1, j); 

  TCanvas *cMatrix = new TCanvas("cMatrix", "cMatrix", 1800, 1000);
  TPad* pM1 = new TPad("pM1","pM1",width1,canBotMargin,width2*0.75,canBotMargin+2*padHeight,0,0);
  pM1->SetTopMargin(0.05);
  pM1->SetBottomMargin(0.05);
  pM1->SetLeftMargin(0.05);
  pM1->SetRightMargin(0.15);
  pM1->SetFillColor(0);
  pM1->SetBorderMode(0);
  pM1->SetBorderSize(0);
  pM1->Draw();

  TPad* pM2 = new TPad("pM2","pM2",width2*0.75,canBotMargin,width2,canBotMargin+2*padHeight,0,0);
  pM2->SetTopMargin(0.05);
  pM2->SetBottomMargin(0.05);
  pM2->SetLeftMargin(0.1);
  pM2->SetRightMargin(0.05);
  pM2->SetFillColor(0);
  pM2->SetBorderMode(0);
  pM2->SetBorderSize(0);
  pM2->Draw();

  pM1->cd();
  pM1->SetGrid();
  pM1->SetFillStyle(4000);
  mCorrMatrix.Draw("colz");


  // If want inverted y-axis for matrix
  // pM1->SetGrid();
  // TH2D *h2Dummy = new TH2D("h2Dummy","", TBackgroundModel::dNParam, 0, TBackgroundModel::dNParam, TBackgroundModel::dNParam, 0, TBackgroundModel::dNParam);
  // h2Dummy->Draw();
  // mCorrMatrix.Draw("colz SAME");

  // TH2C *hgrid = new TH2C("hgrid","",dNParam,0.,dNParam,dNParam,0.,dNParam);
  // hgrid->Draw("SAME");
  // hgrid->GetXaxis()->SetNdivisions(dNParam);
  // hgrid->GetYaxis()->SetNdivisions(dNParam);

  // mCorrMatrixInverse.Draw("colzSAME");
  // Setting axis (inverted axis for matrix)
       // h2Dummy->GetYaxis()->SetLabelOffset(999);
       // // Redraw the new axis
       // gPad->Update();
       // TGaxis *newaxis = new TGaxis(gPad->GetUxmin(),
       //                              gPad->GetUymax(),
       //                              gPad->GetUxmin()-0.0001,
       //                              gPad->GetUymin(),
       //                              h2Dummy->GetYaxis()->GetXmin(),
       //                              h2Dummy->GetYaxis()->GetXmax(),
       //                              505,"+");
       // newaxis->SetLabelOffset(-0.03);
       // newaxis->Draw();
  
       // h2Dummy->GetXaxis()->SetLabelOffset(999);
       // // Redraw the new axis
       // gPad->Update();
       // TGaxis *newaxis2 = new TGaxis(gPad->GetUxmin(),
       //                              gPad->GetUymin(),
       //                              gPad->GetUxmax(),
       //                              gPad->GetUymin(),
       //                              h2Dummy->GetXaxis()->GetXmin(),
       //                              h2Dummy->GetXaxis()->GetXmax(),
       //                              505,"+");
       // newaxis2->SetLabelOffset(0.01);
       // newaxis2->Draw();

  pM2->cd();
  TPaveText *pPave = new TPaveText(0,0,1,1, "NB"); // Text for matrix
  pPave->SetTextSize(0.04);
  pPave->SetFillColor(0);
  pPave->SetBorderSize(0);
  for(int i=0; i < TBackgroundModel::dNParam; i++)
  {
    pPave->AddText(Form("%d: %s", i, minuit->fCpnam[i].Data() ) );
  }
  pPave->Draw();


  double dProgressM1 = 0;
  double dProgressM2 = 0;
  double dProgressM2Sum = 0;
  double dataM1_i = 0, modelM1_i = 0;
  double dataM2_i = 0, modelM2_i = 0;
  double dataM2Sum_i = 0, modelM2Sum_i = 0;

  for(int i = dFitMinBinM1; i < dFitMaxBinM1; i++)
  {
    if( fAdapDataHistoM1->GetBinCenter(i) >= 3050 && fAdapDataHistoM1->GetBinCenter(i) <= 3400)continue;

    // Skipping unknown peaks
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 800 && fAdapDataHistoM1->GetBinCenter(i) <= 808)continue;
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 1060 && fAdapDataHistoM1->GetBinCenter(i) <= 1068)continue; 

    // if( fAdapDataHistoM1->GetBinCenter(i) >= 506 && fAdapDataHistoM1->GetBinCenter(i) <= 515)continue;
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 579 && fAdapDataHistoM1->GetBinCenter(i) <= 589)continue;
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 605 && fAdapDataHistoM1->GetBinCenter(i) <= 615)continue;
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 906 && fAdapDataHistoM1->GetBinCenter(i) <= 917)continue;
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 1445 && fAdapDataHistoM1->GetBinCenter(i) <= 1475)continue;  
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 1755 && fAdapDataHistoM1->GetBinCenter(i) <= 1780)continue;  
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 2090 && fAdapDataHistoM1->GetBinCenter(i) <= 2130)continue;  
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 2200 && fAdapDataHistoM1->GetBinCenter(i) <= 2220)continue;  
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 2440 && fAdapDataHistoM1->GetBinCenter(i) <= 2450)continue;  
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 2600 && fAdapDataHistoM1->GetBinCenter(i) <= 2630)continue;  

    dataM1_i = fAdapDataHistoM1->GetBinContent(i)*fAdapDataHistoM1->GetBinWidth(i);
    modelM1_i = fModelTotAdapM1->GetBinContent(i)*fAdapDataHistoM1->GetBinWidth(i);
    dProgressM1 += 2 * (modelM1_i - dataM1_i + dataM1_i * TMath::Log(dataM1_i/modelM1_i));
    hChiSquaredProgressM1->SetBinContent(i, dProgressM1);
  }
  for(int i = dFitMinBinM2; i < dFitMaxBinM2; i++)
  {
    if( fAdapDataHistoM2->GetBinCenter(i) >= 3050 && fAdapDataHistoM2->GetBinCenter(i) <= 3400)continue;

    // Skipping unknown peaks
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 800 && fAdapDataHistoM2->GetBinCenter(i) <= 808)continue;
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 1060 && fAdapDataHistoM2->GetBinCenter(i) <= 1068)continue; 

    // if( fAdapDataHistoM2->GetBinCenter(i) >= 506 && fAdapDataHistoM2->GetBinCenter(i) <= 515)continue;
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 579 && fAdapDataHistoM2->GetBinCenter(i) <= 589)continue;
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 605 && fAdapDataHistoM2->GetBinCenter(i) <= 615)continue;
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 906 && fAdapDataHistoM2->GetBinCenter(i) <= 917)continue;
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 1450 && fAdapDataHistoM2->GetBinCenter(i) <= 1475)continue;  
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 1755 && fAdapDataHistoM2->GetBinCenter(i) <= 1780)continue;  
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 2090 && fAdapDataHistoM2->GetBinCenter(i) <= 2130)continue;  
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 2200 && fAdapDataHistoM2->GetBinCenter(i) <= 2220)continue;  
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 2440 && fAdapDataHistoM2->GetBinCenter(i) <= 2450)continue;  
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 2600 && fAdapDataHistoM2->GetBinCenter(i) <= 2630)continue;  

    dataM2_i = fAdapDataHistoM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i);
    modelM2_i = fModelTotAdapM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i);
    dProgressM2 += 2 * (modelM2_i - dataM2_i + dataM2_i * TMath::Log(dataM2_i/modelM2_i));
    hChiSquaredProgressM2->SetBinContent(i, dProgressM2);
  }
  for(int i = dFitMinBinM2Sum; i < dFitMaxBinM2Sum; i++)
  {
    if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 3050 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 3400)continue;

    // Skipping unknown peaks
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 800 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 808)continue;
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 1060 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 1068)continue; 

    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 506 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 515)continue;
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 579 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 589)continue;
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 605 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 615)continue;
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 906 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 917)continue;
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 1450 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 1475)continue;  
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 1755 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 1780)continue;  
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 2090 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 2130)continue;  
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 2200 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 2220)continue;  
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 2440 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 2450)continue;  
    // if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 2600 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 2630)continue;  

    dataM2Sum_i = fAdapDataHistoM2Sum->GetBinContent(i)*fAdapDataHistoM2Sum->GetBinWidth(i);
    modelM2Sum_i = fModelTotAdapM2Sum->GetBinContent(i)*fAdapDataHistoM2Sum->GetBinWidth(i);
    dProgressM2Sum += 2 * (modelM2Sum_i - dataM2Sum_i + dataM2Sum_i * TMath::Log(dataM2Sum_i/modelM2Sum_i));
    hChiSquaredProgressM2Sum->SetBinContent(i, dProgressM2Sum);
  }

  TCanvas *cProgress = new TCanvas("cProgress", "cProgress", 1600, 1600);
  cProgress->Divide(1, 2);
  cProgress->cd(1);
  hChiSquaredProgressM1->SetTitle("Increase in #chi^{2} (M1)");
  hChiSquaredProgressM1->GetXaxis()->SetTitle("Energy (keV)");
  hChiSquaredProgressM1->GetYaxis()->SetTitle("#chi^{2}");
  hChiSquaredProgressM1->Draw();

  cProgress->cd(2);
  hChiSquaredProgressM2->SetTitle("Increase in #chi^{2} (M2)");  
  hChiSquaredProgressM2->GetXaxis()->SetTitle("Energy (keV)");
  hChiSquaredProgressM2->GetYaxis()->SetTitle("#chi^{2}");  
  hChiSquaredProgressM2->Draw();

  // cProgress->cd(3);
  // hChiSquaredProgressM2Sum->SetTitle("Increase in #chi^{2} (M2Sum)");  
  // hChiSquaredProgressM2Sum->GetXaxis()->SetTitle("Energy (keV)");
  // hChiSquaredProgressM2Sum->GetYaxis()->SetTitle("#chi^{2}");  
  // hChiSquaredProgressM2Sum->Draw();


  
  if(bSave)
  {
  // Save matrix to file
  // ofstream OutMatrix;
  // OutMatrix.open(Form("%s/FitResults/Test/OutMatrix_%d_%d.txt", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime()));
  // for(int j = mCorrMatrix.GetColLwb(); j <= mCorrMatrix.GetColUpb(); j++)
  // {
  //   for(int i = mCorrMatrix.GetRowLwb(); i <= mCorrMatrix.GetRowUpb(); i++)
  //   {  
  //     OutMatrix << mCorrMatrix(i,j) << "\t";
  //   }
  //   OutMatrix << endl;
  //   OutMatrix << endl;
  // }
  // OutMatrix.close();

  // // Saving plots
    // cadap1->SaveAs(Form("%s/Results/FitM1_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), dBinBase));
    // cadap2->SaveAs(Form("%s/Results/FitM2_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), dBinBase));
    // cResidual1->SaveAs(Form("%s/Results/FitM1Residual_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), dBinBase));
    // cResidual2->SaveAs(Form("%s/Results/FitM2Residual_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), dBinBase));
    // cres1->SaveAs(Form("%s/Results/FitResidualDist_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), dBinBase));
    // cProgress->SaveAs(Form("%s/Results/ChiSquareProgress_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), dBinBase));
    
    // cMatrix->SaveAs(Form("%s/Results/FitCovMatrix_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), dBinBase));
 
    LatexResultTable(0);

    // Save histograms to file
    TFile *fSaveResult = new TFile(Form("%s/Results/FitResult_%d_%d.root", dSaveDir.c_str(), tTime->GetDate(), dBinBase), "RECREATE");
    fSaveResult->Add(fDataHistoM1);
    fSaveResult->Add(fDataHistoM2);
    fSaveResult->Add(fAdapDataHistoM1);
    fSaveResult->Add(fAdapDataHistoM2);
    fSaveResult->Add(fModelTotAdapM1);
    fSaveResult->Add(fModelTotAdapM2);

    for(int i = 0; i < dNParam; i++)
    {
      fSaveResult->Add( BkgPar[i]->GetHistM1() );
      fSaveResult->Add( BkgPar[i]->GetHistM2() );
    
    }
    // fSaveResult->Add(&mCorrMatrix);

    fSaveResult->Write(); 

  // cadap1->SaveAs(Form("%s/FitResults/CovMatrix/FitM1_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cadap2->SaveAs(Form("%s/FitResults/CovMatrix/FitM2_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cResidual1->SaveAs(Form("%s/FitResults/CovMatrix/FitM1Residual_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cResidual2->SaveAs(Form("%s/FitResults/CovMatrix/FitM2Residual_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cres1->SaveAs(Form("%s/FitResults/CovMatrix/FitResidualDist_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cMatrix->SaveAs(Form("%s/FitResults/CovMatrix/FitCovMatrix_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));

    // fSaveResult->Close();
  } // end bSave



  return true;
}

// Txt file with useful stuff
void TBackgroundModel::LatexResultTable(double fValue)
{

  OutFile.open(Form("%s/Results/FitOutputTable_%d_%d.txt", dSaveDir.c_str(), tTime->GetDate(), dBinBase));
  OutFile << "Name -- Normalization -- Normalization Err -- Percent Err (for Norm)" << endl;
  for(int i = 0; i < TBackgroundModel::dNParam; i++)
  {
    OutFile << minuit->fCpnam[i] << "\t" << fParameters[i] << "\t" << fParError[i] << "\t" << fParError[i]/fParameters[i] << endl;
  }
  OutFile.close();



  // Outputs table of best fit values
  // OutFile << "// Name - BestFit - BestFit Err - Integral - Integral Err" << endl;
  // for(int i = 0; i < TBackgroundModel::dNParam; i++)
  // {
  //   OutFile << minuit->fCpnam[i] << Form(" & %.4e$\\pm$%.4e \t\t %.4e$\\pm$%.4e \\\\ ", fParameters[i], fParError[i], dDataIntegralM1*fParameters[i], dDataIntegralM1*fParError[i] ) << endl;
  // }
  // // Outputs table of activities of each cryostat element
  // OutFile <<  "// Name - Activity(events) - Activity (Bq/kg) - Activity (Bq/cm2)" << endl;
  // for(int i = 0; i < TBackgroundModel::dNParam; i++)
  // {
  //   OutFile << minuit->fCpnam[i] << Form(" & %.4e \\pm %.4e", fParActivityM1[i], fParActivityErr[i]) << "\t \t" <<  Form("& %.4e \\pm %.4e", fParActivityM1[i]/fParMass[i]/dLivetime, fParActivityM1Err[i]/fParMass[i]/dLivetime) << "\t \t" << Form("& %.4e \\pm %.4e \\\\", fParActivityM1[i]/fParSurfaceArea[i]/dLivetime, fParActivityM1Err[i]/fParSurfaceArea[i]/dLivetime) << endl;
  // }

}

// Adds a random percentage of events of 2nbb into spectrum
void TBackgroundModel::SanityCheck()
{
  double dM1 = fDataHistoM1->Integral(1, 10000/dBaseBinSize, "width");
  double dM2 = fDataHistoM2->Integral(1, 10000/dBaseBinSize, "width");
  // // Sanity check, adding to background a set amount of 2nbb events, see if code can reconstruct it properly
  // fDataHistoM1->Add(hTeO22nuM1, 0.1*dM1);
  // fDataHistoM2->Add(hTeO22nuM2, 0.1*dM2);

  fAdapDataHistoM1->Add(hAdapTeO22nuM1, 0.25*dM1);
  fAdapDataHistoM2->Add(hAdapTeO22nuM2, 0.25*dM2);
  // fDataHistoM2Sum->Add(hTeO22nuM2Sum, 1.0*dDataIntegralM2Sum);

  cout << "Adding " << (0.25*dM1)*hAdapTeO22nuM1->Integral("width") + (0.25*dM2)*hAdapTeO22nuM2->Integral("width") << " 2nbb events to spectrum" << endl;

}

// Only run in batch mode and make sure to have the 2nbb normalization FIXED
// Probably best to run Minuit in quiet mode as well
void TBackgroundModel::ProfileNLL(double fBestFitInit, double fBestFitChiSq)
{
  // Do the fit normally once first
  DoTheFitAdaptive(0,0);

  // Fix 2nbb value now
  minuit->FixParameter(0);
  
  dBestChiSq = dChiSquare; // Chi-Squared from best fit (for ProfileNLL calculation)
  // Do the fit now if no other tests are needed 
  nLoop = 0;
  for(int i = -25; i < 30; i++)
  {
    fInitValues.push_back(fParameters[0] + fParameters[0]/100*i);
  }



  OutPNLL.open(Form("%s/FitResults/ProfileNLL/ProfileNLL_%d_DR%d.C", dSaveDir.c_str(), tTime->GetDate(), dDataSet ));
  OutPNLL << "{" << endl;
  OutPNLL << "vector<double> dX;" << endl;
  OutPNLL << "vector<double> dT;" << endl;

  for(std::vector<double>::const_iterator iter = fInitValues.begin(); iter!=fInitValues.end(); iter++)
  {
    // cout << "Loop: " << nLoop << endl;
    DoTheFitAdaptive(*iter, 0);
    // LatexResultTable(*iter);
    cout << "delta ChiSq = " << dChiSquare - dBestChiSq << endl; // Needs to be entered, otherwise just 0
    // OutPNLL << Form("dX.push_back(%f); dT.push_back(%f);", (dChiSquare-dBestChiSq)/2., (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1) ) << endl;
    OutPNLL << Form("dX.push_back(%f); dT.push_back(%f);", (dChiSquare-dBestChiSq)/2., (0.69314718056)*(4.726e25 * dLivetimeYr)/(fModelTotAdap2NDBDM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width") + fModelTotAdap2NDBDM2->Integral(1, fAdapDataHistoM2->FindBin(2700) , "width")/2) ) << endl;

    nLoop++; // This is purely for file names and to keep track of number of loops
  }

  OutPNLL << "int n = dX.size();" << endl;
  OutPNLL << "double *y = &dX[0];" << endl;
  OutPNLL << "double *x = &dT[0];" << endl;
  OutPNLL << "TCanvas *cNLL = new TCanvas(\"cNLL\", \"cNLL\", 1200, 800);" << endl;
  OutPNLL << "TGraph *g1 = new TGraph(n, x, y);" << endl;
  OutPNLL << "g1->SetLineColor(kBlue);" << endl;
  OutPNLL << "g1->SetLineWidth(2);" << endl;
  OutPNLL << "g1->SetTitle(\"2#nu#beta#beta Profile Negative Log-Likelihood\");" << endl;
  OutPNLL << "g1->GetYaxis()->SetTitle(\"#Delta#chi^{2}/2\");" << endl;
  OutPNLL << "g1->GetXaxis()->SetTitle(\"t_{1/2} (y)\");" << endl;
  OutPNLL << "g1->Draw(\"AC\");" << endl;
  OutPNLL << "}" << endl;

  OutPNLL.close();
}

void TBackgroundModel::ProfileNLL2D(double fBestFitInit, double fBestFitInit2, double fBestFitChiSq)
{
  // dBestChiSq = fBestFitChiSq; // Chi-Squared from best fit (for ProfileNLL calculation)
  // Do the fit now if no other tests are needed 
  
  // Do the fit normally once first
  DoTheFitAdaptive(0,8.48360e-02);
  dBestChiSq = dChiSquare;

  // Fix 2nbb value now
  minuit->FixParameter(0); 
  // Fix for OVC K40  
  minuit->FixParameter(19);

  nLoop = 0;
  for(int i = -5; i < 5; i++)
  {
    fInitValues.push_back(fParameters[0] + fParameters[0]/75*i);
  }
  for(int j = -5; j < 5; j++)
  {
    fInitValues2.push_back(fParameters[19] + fParameters[19]/75*j);
  }


  OutPNLL.open(Form("%s/FitResults/ProfileNLL/ProfileNLL2D_%d_DR%d.C", dSaveDir.c_str(), tTime->GetDate(), dDataSet ));
  OutPNLL << "{" << endl;
  OutPNLL << "vector<double> dX;" << endl;
  OutPNLL << "vector<double> dY;" << endl;
  OutPNLL << "vector<double> dT;" << endl;

  for(std::vector<double>::const_iterator iter = fInitValues.begin(); iter != fInitValues.end(); iter++)
  {
    for(std::vector<double>::const_iterator iter2 = fInitValues2.begin(); iter2 != fInitValues2.end(); iter2++)
    {
    cout << "Step: " << *iter << " " << *iter2 << endl;    
    DoTheFitAdaptive(*iter, *iter2);
    cout << "delta ChiSq = " << dChiSquare - dBestChiSq << endl; // Needs to be entered, otherwise just 0
    OutPNLL << Form("dX.push_back(%f); dY.push_back(%.10f); dT.push_back(%f);", dChiSquare-dBestChiSq, *iter2, (0.69314718056)*(4.726e25 * dLivetimeYr)/(fModelTotAdap2NDBDM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width") + fModelTotAdap2NDBDM2->Integral(1, fAdapDataHistoM2->FindBin(2700) , "width")/2) ) << endl;

    nLoop++; // This is purely for file names and to keep track of number of loops
    }
  }

  OutPNLL << "int n = dX.size();" << endl;
  OutPNLL << "double *z = &dX[0];" << endl;
  OutPNLL << "double *y = &dY[0];" << endl;
  OutPNLL << "double *x = &dT[0];" << endl;
  OutPNLL << "TCanvas *cNLL = new TCanvas(\"cNLL\", \"cNLL\", 1200, 800);" << endl;
  OutPNLL << "TGraph2D *g1 = new TGraph2D(n, x, y, z);" << endl;
  OutPNLL << "g1->SetLineColor(kBlue);" << endl;
  OutPNLL << "g1->SetLineWidth(2);" << endl;
  OutPNLL << "g1->SetTitle(\"2#nu#beta#beta 2D Profile Negative Log-Likelihood\");" << endl;
  OutPNLL << "g1->GetZaxis()->SetTitle(\"#Delta#chi^{2}\");" << endl;
  OutPNLL << "g1->GetYaxis()->SetTitle(\"OVC K-40\");" << endl;
  OutPNLL << "g1->GetXaxis()->SetTitle(\"t_{1/2} (y)\");" << endl;
  // OutPNLL << "g1->Draw(\"surf1\");" << endl;
  OutPNLL << endl;
  OutPNLL << endl;

  OutPNLL << "TH2D *h1;" << endl;
  OutPNLL << "h1 = g1->GetHistogram();" << endl;
  OutPNLL << "double levels[] = {1, 4, 9};" << endl;
  OutPNLL << "h1->SetContour(3, levels);" << endl;
  OutPNLL << "int colors[] = {kGreen, kYellow, kRed};" << endl;
  OutPNLL << "gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);" << endl;
  OutPNLL << "h1->SetLineWidth(2);" << endl;
  OutPNLL << "h1->Draw(\"cont1\");" << endl;
  OutPNLL << "}" << endl;
  OutPNLL << endl;

  OutPNLL.close();
}


void TBackgroundModel::ToyFit(int fNumFits)
{
    OutToy.open(Form("%s/FitResults/ToyMC/Toy_%d.C", dSaveDir.c_str(), tTime->GetDate() ));

    double dInitial2nbbRate = (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1);
    double dInitial2nbbRateErr = fParError[0]/fParameters[0]*(0.69314718056)*(4.726e25*dLivetimeYr)/(fParameters[0]*dDataIntegralM1);

    cout << "Initial 2nbb Rate and Error: " << dInitial2nbbRate << " +/- " <<  dInitial2nbbRateErr << endl;

    OutToy << "{" << endl;
    OutToy << endl;

    OutToy << "hPullDist = new TH1D(\"hPullDist\", \"Pull Distribution\", 20, -5, 5);" << endl;
    OutToy << "hToy2nbbRate = new TH1D(\"hToy2nbbRate\", \"Toy Monte Carlo half-life fit values\", 100, 6.7e+20, 6.9e+20);" << endl;
    OutToy << "hToy2nbbError = new TH1D(\"hToy2nbbError\", \"Toy Monte Carlo half-life error values\", 50, 1.e+19, 1.2e+19);" << endl;
    OutToy << "hChiSquared = new TH1D(\"hChiSquared\", \"Distribution of Chi-Squared values\", 60, 0, 20);" << endl;

    // cout << "Number of Loops " << fNumFits << endl;

    double dToy2nbbRate;
    double dToy2nbbRateErr;

    // Number of toy fits
    for(int i = 1; i <= 500; i++)
    {
      cout << "Toy: " << i << endl;
      fAdapDataHistoM1->Reset();
      fAdapDataHistoM2->Reset();
      
      fToyData = new TFile(Form("%s/Toy/NoFudge/Toy_p%d.root", dMCDir.c_str(), i));
      fAdapDataHistoM1 = (TH1D*)fToyData->Get("fAdapDataHistoM1");
      fAdapDataHistoM2 = (TH1D*)fToyData->Get("fAdapDataHistoM2");
    
    // fAdapDataHistoM1->Scale(1.41872984571782959e+05/fAdapDataHistoM1->Integral("width"));
    // fAdapDataHistoM2->Scale(2.66538587693367845e+04/fAdapDataHistoM2->Integral("width"));
      // fAdapDataHistoM1->Scale(274875./fAdapDataHistoM1->Integral());
      // fAdapDataHistoM2->Scale(97635.6/fAdapDataHistoM2->Integral());

      for(int j = 1; j <= dAdaptiveBinsM1; j++)
      {
        fAdapDataHistoM1->SetBinError(j, TMath::Sqrt(fAdapDataHistoM1->GetBinContent(j))/fAdapDataHistoM1->GetBinWidth(j));
      }
      for(int k = 1; k <= dAdaptiveBinsM2; k++)
      {
        fAdapDataHistoM2->SetBinError(k, TMath::Sqrt(fAdapDataHistoM2->GetBinContent(k))/fAdapDataHistoM2->GetBinWidth(k));
      }

    // M1 - 450231
    // M2 - 135452
    // M1 - Adap 450162
    // M2 - Adap 135379
    // Scale histograms to have same integral as before
    // dDataIntegralM1 = 450162;
    // dDataIntegralM2 = 135379;
    // dDataIntegralM1 = fAdapDataHistoM1->Integral("width");
    // dDataIntegralM2 = fAdapDataHistoM2->Integral("width");
      DoTheFitAdaptive(0,0);
      // Probably need to add a way to input the best-fit half-life
      // Need to manually input best-fit value right now

      dToy2nbbRate = (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1);
      dToy2nbbRateErr = fParError[0]/fParameters[0]*(0.69314718056)*(4.726e25*dLivetimeYr)/(fParameters[0]*dDataIntegralM1);

      OutToy << Form("hToy2nbbRate->Fill(%.5e);", dToy2nbbRate ) << endl;
      OutToy << Form("hToy2nbbError->Fill(%.5e);", dToy2nbbRateErr ) << endl;
      OutToy << Form("hPullDist->Fill(%5e);", (dToy2nbbRate - dInitial2nbbRate)/(dToy2nbbRateErr) ) << endl;
      OutToy << Form("hChiSquared->Fill(%f);", dChiSquare) << endl;
    }
    OutToy << endl;
    OutToy << endl;
    OutToy << "TCanvas *c1 = new TCanvas(\"c1\", \"c1\", 1600, 1200);" << endl;
    OutToy << "c1->Divide(2,2);" << endl;
    OutToy << "c1->cd(1);" << endl;
    OutToy << "hPullDist->GetXaxis()->SetTitle(\"#Deltat_{1/2}/#sigma_{t_{1/2}}\");" << endl;
    OutToy << "hPullDist->Draw();" << endl;

    OutToy << "c1->cd(2);" << endl;    
    OutToy << "hToy2nbbRate->GetXaxis()->SetTitle(\"t_{1/2}\");" << endl;
    OutToy << "hToy2nbbRate->Draw();" << endl;

    OutToy << "c1->cd(3);" << endl;    
    OutToy << "hToy2nbbError->GetXaxis()->SetTitle(\"#sigma_{t_{1/2}}\");" << endl;
    OutToy << "hToy2nbbError->Draw();" << endl;

    OutToy << "c1->cd(4);" << endl;    
    OutToy << "hChiSquared->GetXaxis()->SetTitle(\"#chi^{2}\");" << endl;
    OutToy << "hChiSquared->Draw();" << endl;

    OutToy << endl;
    OutToy << endl;
    OutToy << "}" << endl;
    OutToy << endl;
    OutToy.close();
}

TH1D *TBackgroundModel::Kernal(TH1D *hMC, TH1D *hSMC)
{
  hSMC->Reset();
  cout << "Smearing: " << hMC->GetName() << endl;

  double dResolution = 1;

  double dNorm;
  double dSmearedValue;

  gaus = new TF1("gaus","gaus(0)", dMinEnergy, dMaxEnergy);


  for(int i = 1; i <= hMC->GetNbinsX(); i++)
  {
    for(int j = 1; j <= hMC->GetNbinsX(); j++)
    {
      // Area of gaussian
      dNorm = dBaseBinSize*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*dResolution);

      // Set parameters of gaussian
      gaus->SetParameters(dNorm, hMC->GetBinCenter(j), dResolution);

      // Smeared contribution from 2nd derivative of gaussian centered at bin j for bin i 
      dSmearedValue = gaus->Derivative2( hSMC->GetBinCenter(i) );

      // Fill bin i with contribution from gaussian centered at bin j
      // cout << "Filling Bin: " << i << " with " << dSmearedValue << endl;
      hSMC->Fill(hSMC->GetBinCenter(i), dSmearedValue);
    }
  }
  
  return hSMC;
}

