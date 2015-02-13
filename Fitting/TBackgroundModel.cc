#include "TMinuit.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TBackgroundModel.hh"
#include "TApplication.h"
#include "TRandom3.h"
#include "TPaveText.h"
#include "TAxis.h"
#include "TLine.h"
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
  for(int i = 0; i < 33; i++ )
  {
    Obj->SetParameters(i, x[i]);
  }
  // Implement a method in your class that calculates the quantity you want to minimize, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimize fval
    Obj->UpdateModelAdaptive();
    fval = Obj->GetChiSquareAdaptive();
}

TBackgroundModel::TBackgroundModel()
{

}

TBackgroundModel::TBackgroundModel(double fFitMin, double fFitMax, int fBinBase, int fDataSet)
{

  tTime = new TDatime();
  dNParam = 33; // number of fitting parameters
  dNumCalls = 0;
  dSecToYears = 1./(60*60*24*365);

  // Data directories depending on QCC/local
  // dDataDir =  "/Users/brian/macros/Simulations/Production/";
  dDataDir =  "/Users/brian/macros/CUOREZ/Bkg";
   // dDataDir = "/cuore/user/zhubrian/CUORE0/scratch/Sync";
  dMCDir = "/Users/brian/macros/Simulations/Production/OldProd";
   // dMCDir = "/cuore/user/zhubrian/MC/scratch/OldProd";
  dSaveDir = "/Users/brian/Dropbox/code/Fitting";
   // dSaveDir = "/cuore/user/zhubrian/code/Fitting";
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

  fEfficiency = new TF1("fEfficiency", "[0]+[1]/(1+[2]*exp(-[3]*x)) + [4]/(1+[5]*exp(-[6]*x))", dMinEnergy, dMaxEnergy);
  fEfficiency->SetParameters(-4.71e-2, 1.12e-1, 2.29, -8.81e-5, 9.68e-1, 2.09, 1.58e-2);

  hEfficiency = new TH1D("hEfficiency", "", dNBins, dMinEnergy, dMaxEnergy);

  for(int i = 1; i <= hEfficiency->GetNbinsX(); i++)
  {
    hEfficiency->SetBinContent(i, fEfficiency->Eval(hEfficiency->GetBinCenter(i)));
  }

  fDataHistoM1->Divide( hEfficiency );
  fDataHistoM2->Divide( hEfficiency );


  // Total model histograms M1
  // M2Sum histograms created but not used (in case I need them someday...)
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

  fModelTotExtM1   = new TH1D("fModelTotExtM1",   "Total External",  dNBins, dMinEnergy, dMaxEnergy);

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

  fModelTotExtM2   = new TH1D("fModelTotExtM2",   "Total External",  dNBins, dMinEnergy, dMaxEnergy);


  // Total model histograms M2Sum
  fModelTotM2Sum      = new TH1D("fModelTotM2Sum",      "Frame",        dNBins, dMinEnergy, dMaxEnergy);  
  fModelTotthM2Sum    = new TH1D("fModelTotthM2Sum",    "Total th232",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotuM2Sum     = new TH1D("fModelTotuM2Sum",     "Total u238",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotkM2Sum     = new TH1D("fModelTotkM2Sum",     "Total k40",    dNBins, dMinEnergy, dMaxEnergy);
  fModelTotcoM2Sum    = new TH1D("fModelTotcoM2Sum",    "Total co60",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotmnM2Sum    = new TH1D("fModelTotmnM2Sum",    "Total mn54",   dNBins, dMinEnergy, dMaxEnergy);

  fModelTotNDBDM2Sum  = new TH1D("fModelTotNDBDM2Sum",  "Total NDBD",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTot2NDBDM2Sum = new TH1D("fModelTot2NDBDM2Sum", "Total 2NDBD",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotbiM2Sum    = new TH1D("fModelTotbiM2Sum",    "Total bi207",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotbi2M2Sum   = new TH1D("fModelTotbi2M2Sum",   "Total bi210",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotptM2Sum    = new TH1D("fModelTotptM2Sum",    "Total pt190",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotpbM2Sum    = new TH1D("fModelTotpbM2Sum",    "Total pb210",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotcsM2Sum    = new TH1D("fModelTotcsM2Sum",    "Total cs137",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotco2M2Sum   = new TH1D("fModelTotco2M2Sum",   "Total co58",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotteo2M2Sum  = new TH1D("fModelTotteo2M2Sum",  "Total TeO2",   dNBins, dMinEnergy, dMaxEnergy);

  fModelTotSthM2Sum   = new TH1D("fModelTotSthM2Sum",   "Total S th232",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotSuM2Sum    = new TH1D("fModelTotSuM2Sum",    "Total S u238",   dNBins, dMinEnergy, dMaxEnergy);
  fModelTotSpbM2Sum   = new TH1D("fModelTotSpbM2Sum",   "Total S pb210",  dNBins, dMinEnergy, dMaxEnergy);
  fModelTotSpoM2Sum   = new TH1D("fModelTotSpoM2Sum",   "Total S po210",  dNBins, dMinEnergy, dMaxEnergy);

  fModelTotExtM2Sum   = new TH1D("fModelTotExtM2Sum",   "Total External",  dNBins, dMinEnergy, dMaxEnergy);

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

  hTeO2th232onlyM1      = new TH1D("hTeO2th232onlyM1", "hTeO2th232onlyM1",     dNBins, dMinEnergy, dMaxEnergy);
  hTeO2ra228pb208M1     = new TH1D("hTeO2ra228pb208M1", "hTeO2ra228pb208M1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2th230onlyM1      = new TH1D("hTeO2th230onlyM1", "hTeO2th230onlyM1",     dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Sxth232onlyM1_001  = new TH1D("hTeO2Sxth232onlyM1_001", "hTeO2Sxth232onlyM1_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra228pb208M1_001 = new TH1D("hTeO2Sxra228pb208M1_001", "hTeO2Sxra228pb208M1_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238th230M1_001  = new TH1D("hTeO2Sxu238th230M1_001", "hTeO2Sxu238th230M1_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth230onlyM1_001  = new TH1D("hTeO2Sxth230onlyM1_001", "hTeO2Sxth230onlyM1_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra226pb210M1_001 = new TH1D("hTeO2Sxra226pb210M1_001", "hTeO2Sxra226pb210M1_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M1_0001     = new TH1D("hTeO2Sxpb210M1_0001", "hTeO2Sxpb210M1_0001",         dNBins, dMinEnergy, dMaxEnergy);


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

  hTeO2th232onlyM2      = new TH1D("hTeO2th232onlyM2", "hTeO2th232onlyM2",     dNBins, dMinEnergy, dMaxEnergy);
  hTeO2ra228pb208M2     = new TH1D("hTeO2ra228pb208M2", "hTeO2ra228pb208M2",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2th230onlyM2      = new TH1D("hTeO2th230onlyM2", "hTeO2th230onlyM2",     dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Sxth232onlyM2_001  = new TH1D("hTeO2Sxth232onlyM2_001", "hTeO2Sxth232onlyM2_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra228pb208M2_001 = new TH1D("hTeO2Sxra228pb208M2_001", "hTeO2Sxra228pb208M2_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238th230M2_001  = new TH1D("hTeO2Sxu238th230M2_001", "hTeO2Sxu238th230M2_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth230onlyM2_001  = new TH1D("hTeO2Sxth230onlyM2_001", "hTeO2Sxth230onlyM2_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra226pb210M2_001 = new TH1D("hTeO2Sxra226pb210M2_001", "hTeO2Sxra226pb210M2_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2_0001     = new TH1D("hTeO2Sxpb210M2_0001", "hTeO2Sxpb210M2_0001",         dNBins, dMinEnergy, dMaxEnergy);


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

  hTeO2th232onlyM2Sum      = new TH1D("hTeO2th232onlyM2Sum", "hTeO2th232onlyM2Sum",     dNBins, dMinEnergy, dMaxEnergy);
  hTeO2ra228pb208M2Sum     = new TH1D("hTeO2ra228pb208M2Sum", "hTeO2ra228pb208M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2th230onlyM2Sum      = new TH1D("hTeO2th230onlyM2Sum", "hTeO2th230onlyM2Sum",     dNBins, dMinEnergy, dMaxEnergy);

  hTeO2Sxth232onlyM2Sum_001  = new TH1D("hTeO2Sxth232onlyM2Sum_001", "hTeO2Sxth232onlyM2Sum_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra228pb208M2Sum_001 = new TH1D("hTeO2Sxra228pb208M2Sum_001", "hTeO2Sxra228pb208M2Sum_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxu238th230M2Sum_001  = new TH1D("hTeO2Sxu238th230M2Sum_001", "hTeO2Sxu238th230M2Sum_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth230onlyM2Sum_001  = new TH1D("hTeO2Sxth230onlyM2Sum_001", "hTeO2Sxth230onlyM2Sum_001",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxra226pb210M2Sum_001 = new TH1D("hTeO2Sxra226pb210M2Sum_001", "hTeO2Sxra226pb210M2Sum_001", dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M2Sum_0001     = new TH1D("hTeO2Sxpb210M2Sum_0001", "hTeO2Sxpb210M2Sum_0001",         dNBins, dMinEnergy, dMaxEnergy);

////////////// Frame M1 and M2
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

  hCuFrameco58M2Sum      = new TH1D("hCuFrameco58M2Sum",   "hCuFrameco58M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameco60M2Sum      = new TH1D("hCuFrameco60M2Sum",   "hCuFrameco60M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFramecs137M2Sum     = new TH1D("hCuFramecs137M2Sum",  "hCuFramecs137M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFramek40M2Sum       = new TH1D("hCuFramek40M2Sum",    "hCuFramek40M2Sum",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFramemn54M2Sum      = new TH1D("hCuFramemn54M2Sum",   "hCuFramemn54M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFramepb210M2Sum     = new TH1D("hCuFramepb210M2Sum",  "hCuFramepb210M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameth232M2Sum     = new TH1D("hCuFrameth232M2Sum",  "hCuFrameth232M2Sum",  dNBins, dMinEnergy, dMaxEnergy);  
  hCuFrameu238M2Sum      = new TH1D("hCuFrameu238M2Sum",   "hCuFrameu238M2Sum",   dNBins, dMinEnergy, dMaxEnergy);

  hCuFrameSth232M2Sum_1     = new TH1D("hCuFrameSth232M2Sum_1",    "hCuFrameSth232M2Sum_1",     dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSu238M2Sum_1      = new TH1D("hCuFrameSu238M2Sum_1",     "hCuFrameSu238M2Sum_1",      dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxpb210M2Sum_001  = new TH1D("hCuFrameSxpb210M2Sum_001", "hCuFrameSxpb210M2Sum_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxpb210M2Sum_01   = new TH1D("hCuFrameSxpb210M2Sum_01",  "hCuFrameSxpb210M2Sum_01",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxpb210M2Sum_1    = new TH1D("hCuFrameSxpb210M2Sum_1",   "hCuFrameSxpb210M2Sum_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxpb210M2Sum_10   = new TH1D("hCuFrameSxpb210M2Sum_10",  "hCuFrameSxpb210M2Sum_10",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxth232M2Sum_001  = new TH1D("hCuFrameSxth232M2Sum_001", "hCuFrameSxth232M2Sum_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxth232M2Sum_01   = new TH1D("hCuFrameSxth232M2Sum_01",  "hCuFrameSxth232M2Sum_01",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxth232M2Sum_1    = new TH1D("hCuFrameSxth232M2Sum_1",   "hCuFrameSxth232M2Sum_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxth232M2Sum_10   = new TH1D("hCuFrameSxth232M2Sum_10",  "hCuFrameSxth232M2Sum_10",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxu238M2Sum_001   = new TH1D("hCuFrameSxu238M2Sum_001",  "hCuFrameSxu238M2Sum_001",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxu238M2Sum_01    = new TH1D("hCuFrameSxu238M2Sum_01",   "hCuFrameSxu238M2Sum_01",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxu238M2Sum_1     = new TH1D("hCuFrameSxu238M2Sum_1",    "hCuFrameSxu238M2Sum_1",     dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameSxu238M2Sum_10    = new TH1D("hCuFrameSxu238M2Sum_10",   "hCuFrameSxu238M2Sum_10",    dNBins, dMinEnergy, dMaxEnergy);

/////////// CuBox (TShield) M1 and M2
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

  hCuBoxco58M2Sum      = new TH1D("hCuBoxco58M2Sum",   "hCuBoxco58M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxco60M2Sum      = new TH1D("hCuBoxco60M2Sum",   "hCuBoxco60M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxcs137M2Sum     = new TH1D("hCuBoxcs137M2Sum",  "hCuBoxcs137M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxk40M2Sum       = new TH1D("hCuBoxk40M2Sum",    "hCuBoxk40M2Sum",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxmn54M2Sum      = new TH1D("hCuBoxmn54M2Sum",   "hCuBoxmn54M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxpb210M2Sum     = new TH1D("hCuBoxpb210M2Sum",  "hCuBoxpb210M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxth232M2Sum     = new TH1D("hCuBoxth232M2Sum",  "hCuBoxth232M2Sum",  dNBins, dMinEnergy, dMaxEnergy);  
  hCuBoxu238M2Sum      = new TH1D("hCuBoxu238M2Sum",   "hCuBoxu238M2Sum",   dNBins, dMinEnergy, dMaxEnergy);

  hCuBoxSth232M2Sum_1     = new TH1D("hCuBoxSth232M2Sum_1",    "hCuBoxSth232M2Sum_1",     dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSu238M2Sum_1      = new TH1D("hCuBoxSu238M2Sum_1",     "hCuBoxSu238M2Sum_1",      dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxpb210M2Sum_001  = new TH1D("hCuBoxSxpb210M2Sum_001", "hCuBoxSxpb210M2Sum_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxpb210M2Sum_01   = new TH1D("hCuBoxSxpb210M2Sum_01",  "hCuBoxSxpb210M2Sum_01",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxpb210M2Sum_1    = new TH1D("hCuBoxSxpb210M2Sum_1",   "hCuBoxSxpb210M2Sum_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxpb210M2Sum_10   = new TH1D("hCuBoxSxpb210M2Sum_10",  "hCuBoxSxpb210M2Sum_10",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxth232M2Sum_001  = new TH1D("hCuBoxSxth232M2Sum_001", "hCuBoxSxth232M2Sum_001",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxth232M2Sum_01   = new TH1D("hCuBoxSxth232M2Sum_01",  "hCuBoxSxth232M2Sum_01",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxth232M2Sum_1    = new TH1D("hCuBoxSxth232M2Sum_1",   "hCuBoxSxth232M2Sum_1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxth232M2Sum_10   = new TH1D("hCuBoxSxth232M2Sum_10",  "hCuBoxSxth232M2Sum_10",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxu238M2Sum_001   = new TH1D("hCuBoxSxu238M2Sum_001",  "hCuBoxSxu238M2Sum_001",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxu238M2Sum_01    = new TH1D("hCuBoxSxu238M2Sum_01",   "hCuBoxSxu238M2Sum_01",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxu238M2Sum_1     = new TH1D("hCuBoxSxu238M2Sum_1",    "hCuBoxSxu238M2Sum_1",     dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxSxu238M2Sum_10    = new TH1D("hCuBoxSxu238M2Sum_10",   "hCuBoxSxu238M2Sum_10",    dNBins, dMinEnergy, dMaxEnergy);

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

/////////// 50mK M1 and M2
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

  h50mKco58M2Sum      = new TH1D("h50mKco58M2Sum",   "h50mKco58M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  h50mKco60M2Sum      = new TH1D("h50mKco60M2Sum",   "h50mKco60M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  h50mKcs137M2Sum     = new TH1D("h50mKcs137M2Sum",  "h50mKcs137M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  h50mKk40M2Sum       = new TH1D("h50mKk40M2Sum",    "h50mKk40M2Sum",    dNBins, dMinEnergy, dMaxEnergy);
  h50mKmn54M2Sum      = new TH1D("h50mKmn54M2Sum",   "h50mKmn54M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  h50mKpb210M2Sum     = new TH1D("h50mKpb210M2Sum",  "h50mKpb210M2Sum",  dNBins, dMinEnergy, dMaxEnergy);
  h50mKth232M2Sum     = new TH1D("h50mKth232M2Sum",  "h50mKth232M2Sum",  dNBins, dMinEnergy, dMaxEnergy);  
  h50mKu238M2Sum      = new TH1D("h50mKu238M2Sum",   "h50mKu238M2Sum",   dNBins, dMinEnergy, dMaxEnergy);

/////////// 600mK M1 and M2
  h600mKco60M1      = new TH1D("h600mKco60M1",   "h600mKco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  h600mKk40M1       = new TH1D("h600mKk40M1",    "h600mKk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  h600mKth232M1     = new TH1D("h600mKth232M1",  "h600mKth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  h600mKu238M1      = new TH1D("h600mKu238M1",   "h600mKu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  h600mKco60M2      = new TH1D("h600mKco60M2",   "h600mKco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  h600mKk40M2       = new TH1D("h600mKk40M2",    "h600mKk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  h600mKth232M2     = new TH1D("h600mKth232M2",  "h600mKth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  h600mKu238M2      = new TH1D("h600mKu238M2",   "h600mKu238M2",   dNBins, dMinEnergy, dMaxEnergy);  

  h600mKco60M2Sum      = new TH1D("h600mKco60M2Sum",   "h600mKco60M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  h600mKk40M2Sum       = new TH1D("h600mKk40M2Sum",    "h600mKk40M2Sum",    dNBins, dMinEnergy, dMaxEnergy);
  h600mKth232M2Sum     = new TH1D("h600mKth232M2Sum",  "h600mKth232M2Sum",  dNBins, dMinEnergy, dMaxEnergy);  
  h600mKu238M2Sum      = new TH1D("h600mKu238M2Sum",   "h600mKu238M2Sum",   dNBins, dMinEnergy, dMaxEnergy); 

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

///////////// Main bath M1 and M2
  hMBco60M1      = new TH1D("hMBco60M1",   "hMBco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hMBk40M1       = new TH1D("hMBk40M1",    "hMBk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hMBth232M1     = new TH1D("hMBth232M1",  "hMBth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hMBu238M1      = new TH1D("hMBu238M1",   "hMBu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hMBco60M2      = new TH1D("hMBco60M2",   "hMBco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hMBk40M2       = new TH1D("hMBk40M2",    "hMBk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hMBth232M2     = new TH1D("hMBth232M2",  "hMBth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hMBu238M2      = new TH1D("hMBu238M2",   "hMBu238M2",   dNBins, dMinEnergy, dMaxEnergy);  

  hMBco60M2Sum      = new TH1D("hMBco60M2Sum",   "hMBco60M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hMBk40M2Sum       = new TH1D("hMBk40M2Sum",    "hMBk40M2Sum",    dNBins, dMinEnergy, dMaxEnergy);
  hMBth232M2Sum     = new TH1D("hMBth232M2Sum",  "hMBth232M2Sum",  dNBins, dMinEnergy, dMaxEnergy);  
  hMBu238M2Sum      = new TH1D("hMBu238M2Sum",   "hMBu238M2Sum",   dNBins, dMinEnergy, dMaxEnergy);  

///////// Super Insulation M1 and M2
  hSIk40M1       = new TH1D("hSIk40M1",    "hSIk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hSIth232M1     = new TH1D("hSIth232M1",  "hSIth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hSIu238M1      = new TH1D("hSIu238M1",   "hSIu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hSIk40M2       = new TH1D("hSIk40M2",    "hSIk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hSIth232M2     = new TH1D("hSIth232M2",  "hSIth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hSIu238M2      = new TH1D("hSIu238M2",   "hSIu238M2",   dNBins, dMinEnergy, dMaxEnergy);


///////// IVC M1 and M2
  hIVCco60M1      = new TH1D("hIVCco60M1",   "hIVCco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hIVCk40M1       = new TH1D("hIVCk40M1",    "hIVCk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hIVCth232M1     = new TH1D("hIVCth232M1",  "hIVCth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hIVCu238M1      = new TH1D("hIVCu238M1",   "hIVCu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hIVCco60M2Sum      = new TH1D("hIVCco60M2Sum",   "hIVCco60M2Sum",   dNBins, dMinEnergy, dMaxEnergy);
  hIVCk40M2Sum       = new TH1D("hIVCk40M2Sum",    "hIVCk40M2Sum",    dNBins, dMinEnergy, dMaxEnergy);
  hIVCth232M2Sum     = new TH1D("hIVCth232M2Sum",  "hIVCth232M2Sum",  dNBins, dMinEnergy, dMaxEnergy);  
  hIVCu238M2Sum      = new TH1D("hIVCu238M2Sum",   "hIVCu238M2Sum",   dNBins, dMinEnergy, dMaxEnergy);  

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



// Mess with rebinning here 
  // fDataHistoM1->Rebin(20);
  // fDataHistoM2->Rebin(20);
  // fDataHistoM2Sum->Rebin(20);
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

  dFitMinBinM1 = fAdapDataHistoM1->FindBin(dFitMin);
  dFitMinBinM2 = fAdapDataHistoM2->FindBin(dFitMin);
  dFitMinBinM2Sum = fAdapDataHistoM2Sum->FindBin(dFitMin);
  dFitMaxBinM1 = fAdapDataHistoM1->FindBin(dFitMax);
  dFitMaxBinM2 = fAdapDataHistoM2->FindBin(dFitMax);
  dFitMaxBinM2Sum = fAdapDataHistoM2Sum->FindBin(dFitMax);

  cout << "Fit M1 from bin: " << dFitMinBinM1 << " to " << dFitMaxBinM1 << endl;
  cout << "Fit M2 from bin: " << dFitMinBinM2 << " to " << dFitMaxBinM2 << endl;
  cout << "Fit M2Sum from bin: " << dFitMinBinM2Sum << " to " << dFitMaxBinM2Sum << endl;


//////////////// Adaptive binned histograms
////////// Total Adaptive binning histograms M1
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

  fModelTotAdapExtM1   = new TH1D("fModelTotAdapExtM1", "Total External", dAdaptiveBinsM1, dAdaptiveArrayM1);

/////////// Total Adaptive binning histograms M2
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

  fModelTotAdapExtM2   = new TH1D("fModelTotAdapExtM2", "Total External", dAdaptiveBinsM2, dAdaptiveArrayM2);

/////////// Total Adaptive binning histograms M2Sum
  fModelTotAdapM2Sum      = new TH1D("fModelTotAdapM2Sum",      "Total PDF M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  fModelTotAdapthM2Sum    = new TH1D("fModelTotAdapthM2Sum",    "Total th232",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapuM2Sum     = new TH1D("fModelTotAdapuM2Sum",     "Total u238",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapkM2Sum     = new TH1D("fModelTotAdapkM2Sum",     "Total k40",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapcoM2Sum    = new TH1D("fModelTotAdapcoM2Sum",    "Total co60",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapmnM2Sum    = new TH1D("fModelTotAdapmnM2Sum",    "Total mn54",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  fModelTotAdapNDBDM2Sum  = new TH1D("fModelTotAdapNDBDM2Sum",  "Total NDBD",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdap2NDBDM2Sum = new TH1D("fModelTotAdap2NDBDM2Sum", "Total 2NDBD",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapbiM2Sum    = new TH1D("fModelTotAdapbiM2Sum",    "Total bi207",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapbi2M2Sum   = new TH1D("fModelTotAdapbi2M2Sum",   "Total bi210",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapptM2Sum    = new TH1D("fModelTotAdapotM2Sum",    "Total pt190",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdappbM2Sum    = new TH1D("fModelTotAdappbM2Sum",    "Total pb210",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapcsM2Sum    = new TH1D("fModelTotAdapcsM2Sum",    "Total cs137",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapco2M2Sum   = new TH1D("fModelTotAdapco2M2Sum",   "Total co58",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapteo2M2Sum  = new TH1D("fModelTotAdapteo2M2Sum",  "Total TeO2",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  fModelTotAdapSthM2Sum   = new TH1D("fModelTotAdapSthM2Sum",   "Total S th232",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapSuM2Sum    = new TH1D("fModelTotAdapSuM2Sum",    "Total S u238",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapSpbM2Sum   = new TH1D("fModelTotAdapSpbM2Sum",   "Total S pb210",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapSpoM2Sum   = new TH1D("fModelTotAdapSpoM2Sum",   "Total S po210",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  fModelTotAdapExtM2Sum   = new TH1D("fModelTotAdapExtM2Sum", "Total External", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);


/////////// Crystal M1 and M2
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

  hAdapTeO2th232onlyM1      = new TH1D("hAdapTeO2th232onlyM1",   "hAdapTeO2th232onlyM1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2ra228pb208M1     = new TH1D("hAdapTeO2ra228pb208M1",  "hAdapTeO2ra228pb208M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2th230onlyM1      = new TH1D("hAdapTeO2th230onlyM1",   "hAdapTeO2th230onlyM1",  dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapTeO2Sxth232onlyM1_001  = new TH1D("hAdapTeO2Sxth232onlyM1_001", "hAdapTeO2Sxth232onlyM1_001",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxra228pb208M1_001 = new TH1D("hAdapTeO2Sxra228pb208M1_001", "hAdapTeO2Sxra228pb208M1_001", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxu238th230M1_001  = new TH1D("hAdapTeO2Sxu238th230M1_001", "hAdapTeO2Sxu238th230M1_001",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxth230onlyM1_001  = new TH1D("hAdapTeO2Sxth230onlyM1_001", "hAdapTeO2Sxth230onlyM1_001",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxra226pb210M1_001 = new TH1D("hAdapTeO2Sxra226pb210M1_001", "hAdapTeO2Sxra226pb210M1_001", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2Sxpb210M1_0001     = new TH1D("hAdapTeO2Sxpb210M1_0001", "hAdapTeO2Sxpb210M1_0001",         dAdaptiveBinsM1, dAdaptiveArrayM1);

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

  hAdapTeO2th232onlyM2      = new TH1D("hAdapTeO2th232onlyM2",   "hAdapTeO2th232onlyM2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2ra228pb208M2     = new TH1D("hAdapTeO2ra228pb208M2",  "hAdapTeO2ra228pb208M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2th230onlyM2      = new TH1D("hAdapTeO2th230onlyM2",   "hAdapTeO2th230onlyM2",  dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapTeO2Sxth232onlyM2_001  = new TH1D("hAdapTeO2Sxth232onlyM2_001", "hAdapTeO2Sxth232onlyM2_001",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxra228pb208M2_001 = new TH1D("hAdapTeO2Sxra228pb208M2_001", "hAdapTeO2Sxra228pb208M2_001", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxu238th230M2_001  = new TH1D("hAdapTeO2Sxu238th230M2_001", "hAdapTeO2Sxu238th230M2_001",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxth230onlyM2_001  = new TH1D("hAdapTeO2Sxth230onlyM2_001", "hAdapTeO2Sxth230onlyM2_001",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxra226pb210M2_001 = new TH1D("hAdapTeO2Sxra226pb210M2_001", "hAdapTeO2Sxra226pb210M2_001", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapTeO2Sxpb210M2_0001     = new TH1D("hAdapTeO2Sxpb210M2_0001", "hAdapTeO2Sxpb210M2_0001",         dAdaptiveBinsM2, dAdaptiveArrayM2);


  hAdapTeO20nuM2Sum       = new TH1D("hAdapTeO20nuM2Sum",    "hAdapTeO20nuM2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO22nuM2Sum       = new TH1D("hAdapTeO22nuM2Sum",    "hAdapTeO22nuM2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2co60M2Sum      = new TH1D("hAdapTeO2co60M2Sum",   "hAdapTeO2co60M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2k40M2Sum       = new TH1D("hAdapTeO2k40M2Sum",    "hAdapTeO2k40M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2pb210M2Sum     = new TH1D("hAdapTeO2pb210M2Sum",  "hAdapTeO2pb210M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2po210M2Sum     = new TH1D("hAdapTeO2po210M2Sum",  "hAdapTeO2po210M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2te125M2Sum     = new TH1D("hAdapTeO2te125M2Sum",  "hAdapTeO2te125M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2th232M2Sum     = new TH1D("hAdapTeO2th232M2Sum",  "hAdapTeO2th232M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapTeO2th228M2Sum     = new TH1D("hAdapTeO2th228M2Sum",  "hAdapTeO2th228M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2ra226M2Sum     = new TH1D("hAdapTeO2ra226M2Sum",  "hAdapTeO2ra226M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2rn222M2Sum     = new TH1D("hAdapTeO2rn222M2Sum",  "hAdapTeO2rn222M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2u238M2Sum      = new TH1D("hAdapTeO2u238M2Sum",   "hAdapTeO2u238M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2th230M2Sum     = new TH1D("hAdapTeO2th230M2Sum",  "hAdapTeO2th230M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2u234M2Sum      = new TH1D("hAdapTeO2u234M2Sum",   "hAdapTeO2u234M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapTeO2Spb210M2Sum_01      = new TH1D("hAdapTeO2Spb210M2Sum_01",   "hAdapTeO2Spb210M2Sum_01",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Spo210M2Sum_001     = new TH1D("hAdapTeO2Spo210M2Sum_001",  "hAdapTeO2Spo210M2Sum_001",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Spo210M2Sum_01      = new TH1D("hAdapTeO2Spo210M2Sum_01",   "hAdapTeO2Spo210M2Sum_01",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sth232M2Sum_01      = new TH1D("hAdapTeO2Sth232M2Sum_01",   "hAdapTeO2Sth232M2Sum_01",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Su238M2Sum_01       = new TH1D("hAdapTeO2Su238M2Sum_01",    "hAdapTeO2Su238M2Sum_01",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpb210M2Sum_001    = new TH1D("hAdapTeO2Sxpb210M2Sum_001", "hAdapTeO2Sxpb210M2Sum_001", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpb210M2Sum_01     = new TH1D("hAdapTeO2Sxpb210M2Sum_01",  "hAdapTeO2Sxpb210M2Sum_01",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpb210M2Sum_1      = new TH1D("hAdapTeO2Sxpb210M2Sum_1",   "hAdapTeO2Sxpb210M2Sum_1",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpb210M2Sum_10     = new TH1D("hAdapTeO2Sxpb210M2Sum_10",  "hAdapTeO2Sxpb210M2Sum_10",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpo210M2Sum_001    = new TH1D("hAdapTeO2Sxpo210M2Sum_001", "hAdapTeO2Sxpo210M2Sum_001", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpo210M2Sum_01     = new TH1D("hAdapTeO2Sxpo210M2Sum_01",  "hAdapTeO2Sxpo210M2Sum_01",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpo210M2Sum_1      = new TH1D("hAdapTeO2Sxpo210M2Sum_1",   "hAdapTeO2Sxpo210M2Sum_1",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxth232M2Sum_001    = new TH1D("hAdapTeO2Sxth232M2Sum_001", "hAdapTeO2Sxth232M2Sum_001", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxth232M2Sum_01     = new TH1D("hAdapTeO2Sxth232M2Sum_01",  "hAdapTeO2Sxth232M2Sum_01",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxth232M2Sum_1      = new TH1D("hAdapTeO2Sxth232M2Sum_1",   "hAdapTeO2Sxth232M2Sum_1",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxth232M2Sum_10     = new TH1D("hAdapTeO2Sxth232M2Sum_10",  "hAdapTeO2Sxth232M2Sum_10",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxu238M2Sum_001     = new TH1D("hAdapTeO2Sxu238M2Sum_001",  "hAdapTeO2Sxu238M2Sum_001",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxu238M2Sum_01      = new TH1D("hAdapTeO2Sxu238M2Sum_01",   "hAdapTeO2Sxu238M2Sum_01",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxu238M2Sum_1       = new TH1D("hAdapTeO2Sxu238M2Sum_1",    "hAdapTeO2Sxu238M2Sum_1",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxu238M2Sum_10      = new TH1D("hAdapTeO2Sxu238M2Sum_10",   "hAdapTeO2Sxu238M2Sum_10",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapTeO2th232onlyM2Sum      = new TH1D("hAdapTeO2th232onlyM2Sum",   "hAdapTeO2th232onlyM2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2ra228pb208M2Sum     = new TH1D("hAdapTeO2ra228pb208M2Sum",  "hAdapTeO2ra228pb208M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2th230onlyM2Sum      = new TH1D("hAdapTeO2th230onlyM2Sum",   "hAdapTeO2th230onlyM2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapTeO2Sxth232onlyM2Sum_001  = new TH1D("hAdapTeO2Sxth232onlyM2Sum_001", "hAdapTeO2Sxth232onlyM2Sum_001",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxra228pb208M2Sum_001 = new TH1D("hAdapTeO2Sxra228pb208M2Sum_001", "hAdapTeO2Sxra228pb208M2Sum_001", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxu238th230M2Sum_001  = new TH1D("hAdapTeO2Sxu238th230M2Sum_001", "hAdapTeO2Sxu238th230M2Sum_001",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxth230onlyM2Sum_001  = new TH1D("hAdapTeO2Sxth230onlyM2Sum_001", "hAdapTeO2Sxth230onlyM2Sum_001",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxra226pb210M2Sum_001 = new TH1D("hAdapTeO2Sxra226pb210M2Sum_001", "hAdapTeO2Sxra226pb210M2Sum_001", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapTeO2Sxpb210M2Sum_0001     = new TH1D("hAdapTeO2Sxpb210M2Sum_0001", "hAdapTeO2Sxpb210M2Sum_0001",         dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);



////////// Frame M1 and M2
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

  hAdapCuFrameco58M2Sum      = new TH1D("hAdapCuFrameco58M2Sum",   "hAdapCuFrameco58M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameco60M2Sum      = new TH1D("hAdapCuFrameco60M2Sum",   "hAdapCuFrameco60M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFramecs137M2Sum     = new TH1D("hAdapCuFramecs137M2Sum",  "hAdapCuFramecs137M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFramek40M2Sum       = new TH1D("hAdapCuFramek40M2Sum",    "hAdapCuFramek40M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFramemn54M2Sum      = new TH1D("hAdapCuFramemn54M2Sum",   "hAdapCuFramemn54M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFramepb210M2Sum     = new TH1D("hAdapCuFramepb210M2Sum",  "hAdapCuFramepb210M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameth232M2Sum     = new TH1D("hAdapCuFrameth232M2Sum",  "hAdapCuFrameth232M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapCuFrameu238M2Sum      = new TH1D("hAdapCuFrameu238M2Sum",   "hAdapCuFrameu238M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapCuFrameSth232M2Sum_1     = new TH1D("hAdapCuFrameSth232M2Sum_1",    "hAdapCuFrameSth232M2Sum_1",     dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSu238M2Sum_1      = new TH1D("hAdapCuFrameSu238M2Sum_1",     "hAdapCuFrameSu238M2Sum_1",      dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxpb210M2Sum_001  = new TH1D("hAdapCuFrameSxpb210M2Sum_001", "hAdapCuFrameSxpb210M2Sum_001",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxpb210M2Sum_01   = new TH1D("hAdapCuFrameSxpb210M2Sum_01",  "hAdapCuFrameSxpb210M2Sum_01",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxpb210M2Sum_1    = new TH1D("hAdapCuFrameSxpb210M2Sum_1",   "hAdapCuFrameSxpb210M2Sum_1",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxpb210M2Sum_10   = new TH1D("hAdapCuFrameSxpb210M2Sum_10",  "hAdapCuFrameSxpb210M2Sum_10",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxth232M2Sum_001  = new TH1D("hAdapCuFrameSxth232M2Sum_001", "hAdapCuFrameSxth232M2Sum_001",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxth232M2Sum_01   = new TH1D("hAdapCuFrameSxth232M2Sum_01",  "hAdapCuFrameSxth232M2Sum_01",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxth232M2Sum_1    = new TH1D("hAdapCuFrameSxth232M2Sum_1",   "hAdapCuFrameSxth232M2Sum_1",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxth232M2Sum_10   = new TH1D("hAdapCuFrameSxth232M2Sum_10",  "hAdapCuFrameSxth232M2Sum_10",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxu238M2Sum_001   = new TH1D("hAdapCuFrameSxu238M2Sum_001",  "hAdapCuFrameSxu238M2Sum_001",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxu238M2Sum_01    = new TH1D("hAdapCuFrameSxu238M2Sum_01",   "hAdapCuFrameSxu238M2Sum_01",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxu238M2Sum_1     = new TH1D("hAdapCuFrameSxu238M2Sum_1",    "hAdapCuFrameSxu238M2Sum_1",     dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxu238M2Sum_10    = new TH1D("hAdapCuFrameSxu238M2Sum_10",   "hAdapCuFrameSxu238M2Sum_10",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

////////// CuBox (TShield) M1 and M2
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

  hAdapCuBoxco58M2Sum      = new TH1D("hAdapCuBoxco58M2Sum",   "hAdapCuBoxco58M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxco60M2Sum      = new TH1D("hAdapCuBoxco60M2Sum",   "hAdapCuBoxco60M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxcs137M2Sum     = new TH1D("hAdapCuBoxcs137M2Sum",  "hAdapCuBoxcs137M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxk40M2Sum       = new TH1D("hAdapCuBoxk40M2Sum",    "hAdapCuBoxk40M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxmn54M2Sum      = new TH1D("hAdapCuBoxmn54M2Sum",   "hAdapCuBoxmn54M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxpb210M2Sum     = new TH1D("hAdapCuBoxpb210M2Sum",  "hAdapCuBoxpb210M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxth232M2Sum     = new TH1D("hAdapCuBoxth232M2Sum",  "hAdapCuBoxth232M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapCuBoxu238M2Sum      = new TH1D("hAdapCuBoxu238M2Sum",   "hAdapCuBoxu238M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapCuBoxSth232M2Sum_1     = new TH1D("hAdapCuBoxSth232M2Sum_1",    "hAdapCuBoxSth232M2Sum_1",     dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSu238M2Sum_1      = new TH1D("hAdapCuBoxSu238M2Sum_1",     "hAdapCuBoxSu238M2Sum_1",      dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxpb210M2Sum_001  = new TH1D("hAdapCuBoxSxpb210M2Sum_001", "hAdapCuBoxSxpb210M2Sum_001",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxpb210M2Sum_01   = new TH1D("hAdapCuBoxSxpb210M2Sum_01",  "hAdapCuBoxSxpb210M2Sum_01",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxpb210M2Sum_1    = new TH1D("hAdapCuBoxSxpb210M2Sum_1",   "hAdapCuBoxSxpb210M2Sum_1",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxpb210M2Sum_10   = new TH1D("hAdapCuBoxSxpb210M2Sum_10",  "hAdapCuBoxSxpb210M2Sum_10",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxth232M2Sum_001  = new TH1D("hAdapCuBoxSxth232M2Sum_001", "hAdapCuBoxSxth232M2Sum_001",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxth232M2Sum_01   = new TH1D("hAdapCuBoxSxth232M2Sum_01",  "hAdapCuBoxSxth232M2Sum_01",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxth232M2Sum_1    = new TH1D("hAdapCuBoxSxth232M2Sum_1",   "hAdapCuBoxSxth232M2Sum_1",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxth232M2Sum_10   = new TH1D("hAdapCuBoxSxth232M2Sum_10",  "hAdapCuBoxSxth232M2Sum_10",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxu238M2Sum_001   = new TH1D("hAdapCuBoxSxu238M2Sum_001",  "hAdapCuBoxSxu238M2Sum_001",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxu238M2Sum_01    = new TH1D("hAdapCuBoxSxu238M2Sum_01",   "hAdapCuBoxSxu238M2Sum_01",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxu238M2Sum_1     = new TH1D("hAdapCuBoxSxu238M2Sum_1",    "hAdapCuBoxSxu238M2Sum_1",     dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxu238M2Sum_10    = new TH1D("hAdapCuBoxSxu238M2Sum_10",   "hAdapCuBoxSxu238M2Sum_10",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

//////////// CuBox + CuFrame M1 and M2

  hAdapCuBox_CuFrameco60M1 = new TH1D("hAdapCuBox_CuFrameco60M1", "hAdapCuBox_CuFrameco60M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramek40M1 = new TH1D("hAdapCuBox_CuFramek40M1", "hAdapCuBox_CuFramek40M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameth232M1 = new TH1D("hAdapCuBox_CuFrameth232M1", "hAdapCuBox_CuFrameth232M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameu238M1 = new TH1D("hAdapCuBox_CuFrameu238M1", "hAdapCuBox_CuFrameu238M1", dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuBox_CuFrameth232M1_10 = new TH1D("hAdapCuBox_CuFrameth232M1_10", "hAdapCuBox_CuFrameth232M1_10", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameu238M1_10 = new TH1D("hAdapCuBox_CuFrameu238M1_10", "hAdapCuBox_CuFrameu238M1_10", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramepb210M1_10 = new TH1D("hAdapCuBox_CuFramepb210M1_10", "hAdapCuBox_CuFramepb210M1_10", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramepb210M1_1 = new TH1D("hAdapCuBox_CuFramepb210M1_1", "hAdapCuBox_CuFramepb210M1_1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramepb210M1_01 = new TH1D("hAdapCuBox_CuFramepb210M1_01", "hAdapCuBox_CuFramepb210M1_01", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramepb210M1_001 = new TH1D("hAdapCuBox_CuFramepb210M1_001", "hAdapCuBox_CuFramepb210M1_001", dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuBox_CuFrameco60M2 = new TH1D("hAdapCuBox_CuFrameco60M2", "hAdapCuBox_CuFrameco60M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFramek40M2 = new TH1D("hAdapCuBox_CuFramek40M2", "hAdapCuBox_CuFramek40M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameth232M2 = new TH1D("hAdapCuBox_CuFrameth232M2", "hAdapCuBox_CuFrameth232M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameu238M2 = new TH1D("hAdapCuBox_CuFrameu238M2", "hAdapCuBox_CuFrameu238M2", dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapCuBox_CuFrameth232M2_10 = new TH1D("hAdapCuBox_CuFrameth232M2_10", "hAdapCuBox_CuFrameth232M2_10", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFrameu238M2_10 = new TH1D("hAdapCuBox_CuFrameu238M2_10", "hAdapCuBox_CuFrameu238M2_10", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFramepb210M2_10 = new TH1D("hAdapCuBox_CuFramepb210M2_10", "hAdapCuBox_CuFramepb210M2_10", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFramepb210M2_1 = new TH1D("hAdapCuBox_CuFramepb210M2_1", "hAdapCuBox_CuFramepb210M2_1", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFramepb210M2_01 = new TH1D("hAdapCuBox_CuFramepb210M2_01", "hAdapCuBox_CuFramepb210M2_01", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFramepb210M2_001 = new TH1D("hAdapCuBox_CuFramepb210M2_001", "hAdapCuBox_CuFramepb210M2_001", dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapCuBox_CuFrameco60M2Sum = new TH1D("hAdapCuBox_CuFrameco60M2Sum", "hAdapCuBox_CuFrameco60M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFramek40M2Sum = new TH1D("hAdapCuBox_CuFramek40M2Sum", "hAdapCuBox_CuFramek40M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameth232M2Sum = new TH1D("hAdapCuBox_CuFrameth232M2Sum", "hAdapCuBox_CuFrameth232M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameu238M2Sum = new TH1D("hAdapCuBox_CuFrameu238M2Sum", "hAdapCuBox_CuFrameu238M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapCuBox_CuFrameth232M2Sum_10 = new TH1D("hAdapCuBox_CuFrameth232M2Sum_10", "hAdapCuBox_CuFrameth232M2Sum_10", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFrameu238M2Sum_10 = new TH1D("hAdapCuBox_CuFrameu238M2Sum_10", "hAdapCuBox_CuFrameu238M2Sum_10", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFramepb210M2Sum_10 = new TH1D("hAdapCuBox_CuFramepb210M2Sum_10", "hAdapCuBox_CuFramepb210M2Sum_10", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFramepb210M2Sum_1 = new TH1D("hAdapCuBox_CuFramepb210M2Sum_1", "hAdapCuBox_CuFramepb210M2Sum_1", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFramepb210M2Sum_01 = new TH1D("hAdapCuBox_CuFramepb210M2Sum_01", "hAdapCuBox_CuFramepb210M2Sum_01", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBox_CuFramepb210M2Sum_001 = new TH1D("hAdapCuBox_CuFramepb210M2Sum_001", "hAdapCuBox_CuFramepb210M2Sum_001", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);


////////// 50mK M1 and M2
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

  hAdap50mKco58M2Sum      = new TH1D("hAdap50mKco58M2Sum",   "hAdap50mKco58M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap50mKco60M2Sum      = new TH1D("hAdap50mKco60M2Sum",   "hAdap50mKco60M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap50mKcs137M2Sum     = new TH1D("hAdap50mKcs137M2Sum",  "hAdap50mKcs137M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap50mKk40M2Sum       = new TH1D("hAdap50mKk40M2Sum",    "hAdap50mKk40M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap50mKmn54M2Sum      = new TH1D("hAdap50mKmn54M2Sum",   "hAdap50mKmn54M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap50mKpb210M2Sum     = new TH1D("hAdap50mKpb210M2Sum",  "hAdap50mKpb210M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap50mKth232M2Sum     = new TH1D("hAdap50mKth232M2Sum",  "hAdap50mKth232M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdap50mKu238M2Sum      = new TH1D("hAdap50mKu238M2Sum",   "hAdap50mKu238M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

////////// 600mK M1 and M2
  hAdap600mKco60M1      = new TH1D("hAdap600mKco60M1",   "hAdap600mKco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap600mKk40M1       = new TH1D("hAdap600mKk40M1",    "hAdap600mKk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap600mKth232M1     = new TH1D("hAdap600mKth232M1",  "hAdap600mKth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdap600mKu238M1      = new TH1D("hAdap600mKu238M1",   "hAdap600mKu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdap600mKco60M2      = new TH1D("hAdap600mKco60M2",   "hAdap600mKco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap600mKk40M2       = new TH1D("hAdap600mKk40M2",    "hAdap600mKk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap600mKth232M2     = new TH1D("hAdap600mKth232M2",  "hAdap600mKth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdap600mKu238M2      = new TH1D("hAdap600mKu238M2",   "hAdap600mKu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  

  hAdap600mKco60M2Sum      = new TH1D("hAdap600mKco60M2Sum",   "hAdap600mKco60M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap600mKk40M2Sum       = new TH1D("hAdap600mKk40M2Sum",    "hAdap600mKk40M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap600mKth232M2Sum     = new TH1D("hAdap600mKth232M2Sum",  "hAdap600mKth232M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdap600mKu238M2Sum      = new TH1D("hAdap600mKu238M2Sum",   "hAdap600mKu238M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

////////// Roman Lead M1 and M2
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

  hAdapPbRombi207M2Sum     = new TH1D("hAdapPbRombi207M2Sum",  "hAdapPbRombi207M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapPbRomco60M2Sum      = new TH1D("hAdapPbRomco60M2Sum",   "hAdapPbRomco60M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapPbRomcs137M2Sum     = new TH1D("hAdapPbRomcs137M2Sum",  "hAdapPbRomcs137M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapPbRomk40M2Sum       = new TH1D("hAdapPbRomk40M2Sum",    "hAdapPbRomk40M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapPbRompb210M2Sum     = new TH1D("hAdapPbRompb210M2Sum",  "hAdapPbRompb210M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapPbRomth232M2Sum     = new TH1D("hAdapPbRomth232M2Sum",  "hAdapPbRomth232M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapPbRomu238M2Sum      = new TH1D("hAdapPbRomu238M2Sum",   "hAdapPbRomu238M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);


///////// Internal Shields M1 and M2
  hAdapInternalco60M1 = new TH1D("hAdapInternalco60M1", "hAdapInternalco60M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapInternalk40M1 = new TH1D("hAdapInternalk40M1", "hAdapInternalk40M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapInternalth232M1 = new TH1D("hAdapInternalth232M1", "hAdapInternalth232M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapInternalu238M1 = new TH1D("hAdapInternalu238M1", "hAdapInternalu238M1", dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapInternalco60M2 = new TH1D("hAdapInternalco60M2", "hAdapInternalco60M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapInternalk40M2 = new TH1D("hAdapInternalk40M2", "hAdapInternalk40M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapInternalth232M2 = new TH1D("hAdapInternalth232M2", "hAdapInternalth232M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapInternalu238M2 = new TH1D("hAdapInternalu238M2", "hAdapInternalu238M2", dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapInternalco60M2Sum = new TH1D("hAdapInternalco60M2Sum", "hAdapInternalco60M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapInternalk40M2Sum = new TH1D("hAdapInternalk40M2Sum", "hAdapInternalk40M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapInternalth232M2Sum = new TH1D("hAdapInternalth232M2Sum", "hAdapInternalth232M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapInternalu238M2Sum = new TH1D("hAdapInternalu238M2Sum", "hAdapInternalu238M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  

////////////// Main bath M1 and M2
  hAdapMBco60M1      = new TH1D("hAdapMBco60M1",   "hAdapMBco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapMBk40M1       = new TH1D("hAdapMBk40M1",    "hAdapMBk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapMBth232M1     = new TH1D("hAdapMBth232M1",  "hAdapMBth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapMBu238M1      = new TH1D("hAdapMBu238M1",   "hAdapMBu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapMBco60M2      = new TH1D("hAdapMBco60M2",   "hAdapMBco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapMBk40M2       = new TH1D("hAdapMBk40M2",    "hAdapMBk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapMBth232M2     = new TH1D("hAdapMBth232M2",  "hAdapMBth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapMBu238M2      = new TH1D("hAdapMBu238M2",   "hAdapMBu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  

  hAdapMBco60M2Sum      = new TH1D("hAdapMBco60M2Sum",   "hAdapMBco60M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapMBk40M2Sum       = new TH1D("hAdapMBk40M2Sum",    "hAdapMBk40M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapMBth232M2Sum     = new TH1D("hAdapMBth232M2Sum",  "hAdapMBth232M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapMBu238M2Sum      = new TH1D("hAdapMBu238M2Sum",   "hAdapMBu238M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  

/////////// Super Insulation M1 and M2
  hAdapSIk40M1       = new TH1D("hAdapSIk40M1",    "hAdapSIk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1); 
  hAdapSIth232M1     = new TH1D("hAdapSIth232M1",  "hAdapSIth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapSIu238M1      = new TH1D("hAdapSIu238M1",   "hAdapSIu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapSIk40M2       = new TH1D("hAdapSIk40M2",    "hAdapSIk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapSIth232M2     = new TH1D("hAdapSIth232M2",  "hAdapSIth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapSIu238M2      = new TH1D("hAdapSIu238M2",   "hAdapSIu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapSIk40M2Sum       = new TH1D("hAdapSIk40M2Sum",    "hAdapSIk40M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapSIth232M2Sum     = new TH1D("hAdapSIth232M2Sum",  "hAdapSIth232M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapSIu238M2Sum      = new TH1D("hAdapSIu238M2Sum",   "hAdapSIu238M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

//////////// IVC M1 and M2
  hAdapIVCco60M1      = new TH1D("hAdapIVCco60M1",   "hAdapIVCco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapIVCk40M1       = new TH1D("hAdapIVCk40M1",    "hAdapIVCk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapIVCth232M1     = new TH1D("hAdapIVCth232M1",  "hAdapIVCth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapIVCu238M1      = new TH1D("hAdapIVCu238M1",   "hAdapIVCu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapIVCco60M2      = new TH1D("hAdapIVCco60M2",   "hAdapIVCco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapIVCk40M2       = new TH1D("hAdapIVCk40M2",    "hAdapIVCk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapIVCth232M2     = new TH1D("hAdapIVCth232M2",  "hAdapIVCth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapIVCu238M2      = new TH1D("hAdapIVCu238M2",   "hAdapIVCu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  

  hAdapIVCco60M2Sum      = new TH1D("hAdapIVCco60M2Sum",   "hAdapIVCco60M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapIVCk40M2Sum       = new TH1D("hAdapIVCk40M2Sum",    "hAdapIVCk40M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapIVCth232M2Sum     = new TH1D("hAdapIVCth232M2Sum",  "hAdapIVCth232M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapIVCu238M2Sum      = new TH1D("hAdapIVCu238M2Sum",   "hAdapIVCu238M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum); 

////////////// OVC M1 and M2
  hAdapOVCco60M1      = new TH1D("hAdapOVCco60M1",   "hAdapOVCco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapOVCk40M1       = new TH1D("hAdapOVCk40M1",    "hAdapOVCk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapOVCth232M1     = new TH1D("hAdapOVCth232M1",  "hAdapOVCth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapOVCu238M1      = new TH1D("hAdapOVCu238M1",   "hAdapOVCu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapOVCco60M2      = new TH1D("hAdapOVCco60M2",   "hAdapOVCco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapOVCk40M2       = new TH1D("hAdapOVCk40M2",    "hAdapOVCk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapOVCth232M2     = new TH1D("hAdapOVCth232M2",  "hAdapOVCth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapOVCu238M2      = new TH1D("hAdapOVCu238M2",   "hAdapOVCu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  

  hAdapOVCco60M2Sum      = new TH1D("hAdapOVCco60M2Sum",   "hAdapOVCco60M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapOVCk40M2Sum       = new TH1D("hAdapOVCk40M2Sum",    "hAdapOVCk40M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapOVCth232M2Sum     = new TH1D("hAdapOVCth232M2Sum",  "hAdapOVCth232M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapOVCu238M2Sum      = new TH1D("hAdapOVCu238M2Sum",   "hAdapOVCu238M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  

/////////// External Sources M1 and M2
  hAdapExtPbbi210M1 = new TH1D("hAdapExtPbbi210M1", "hAdapExtPbbi210M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapExtPbbi210M2 = new TH1D("hAdapExtPbbi210M2", "hAdapExtPbbi210M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapExtPbbi210M2Sum = new TH1D("hAdapExtPbbi210M2Sum", "hAdapExtPbbi210M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);


  hEnergyScaleDummyM1 = new TH1D("hEnergyScaleDummyM1",   "Energy Scale M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hEnergyScaleDummyM2 = new TH1D("hEnergyScaleDummyM2",   "Energy Scale M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hEnergyScaleDummyM2Sum = new TH1D("hEnergyScaleDummyM2Sum",   "Energy Scale M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);


  // Loads all of the PDFs from file
  Initialize();

/*
  // Add 2nbb events to background
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
*/

  dDataIntegral = fDataHistoM1->Integral(1, dNBins);
  dDataIntegralM1 = fDataHistoM1->Integral(1, 10000/dBinSize);
  dDataIntegralM2 = fDataHistoM2->Integral(1, 10000/dBinSize);
  dDataIntegralM2Sum = fDataHistoM2Sum->Integral(1, 10000/dBinSize);
  dDataIntegralTot = qtree->GetEntries();

  cout << "Total Events in background spectrum: " << dDataIntegralTot << endl; 
  cout << "Events in background spectrum (M1): " << dDataIntegralM1 << endl;
  cout << "Events in background spectrum (M2): " << dDataIntegralM2 << endl;
  cout << "Events in background spectrum (M2Sum): " << dDataIntegralM2Sum << endl;
  cout << "Livetime of background: " << dLivetimeYr << endl;

  dBestChiSq = 0; // Chi-Squared from best fit (for ProfileNLL calculation)
  // Do the fit now if no other tests are needed 
  nLoop = 0;
  DoTheFitAdaptive(0);
  // LatexResultTable(0); 
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
  delete fModelTotthM1;
  delete fModelTotuM1;
  delete fModelTotkM1;
  delete fModelTotcoM1;
  delete fModelTotmnM1;
  delete fModelTotNDBDM1;
  delete fModelTot2NDBDM1;
  delete fModelTotbiM1;
  delete fModelTotbi2M1;  
  delete fModelTotptM1;
  delete fModelTotpbM1;
  delete fModelTotcsM1;
  delete fModelTotco2M1;
  delete fModelTotteo2M1;

  delete fModelTotSthM1;
  delete fModelTotSuM1;
  delete fModelTotSpoM1;
  delete fModelTotSpbM1;
  delete fModelTotExtM1;

  delete fModelTotAdapM1;
  delete fModelTotAdapthM1;
  delete fModelTotAdapuM1;
  delete fModelTotAdapkM1;
  delete fModelTotAdapcoM1;
  delete fModelTotAdapmnM1;
  delete fModelTotAdapNDBDM1;
  delete fModelTotAdap2NDBDM1;
  delete fModelTotAdapbiM1;
  delete fModelTotAdapbi2M1;
  delete fModelTotAdapptM1;
  delete fModelTotAdappbM1;
  delete fModelTotAdapcsM1;
  delete fModelTotAdapco2M1;
  delete fModelTotAdapteo2M1;

  delete fModelTotAdapSthM1;
  delete fModelTotAdapSuM1;
  delete fModelTotAdapSpoM1;
  delete fModelTotAdapSpbM1;
  delete fModelTotAdapExtM1;


  delete fModelTotAdapAlphaM1;
  delete fModelTotAdapAlphaHighM1;
  delete fModelTotAdapAlphaLowM1;

  // Total PDFs M2
  delete fModelTotM2;
  delete fModelTotthM2;
  delete fModelTotuM2;
  delete fModelTotkM2;
  delete fModelTotcoM2;
  delete fModelTotmnM2;
  delete fModelTotNDBDM2;
  delete fModelTot2NDBDM2;
  delete fModelTotbiM2;
  delete fModelTotbi2M2;
  delete fModelTotptM2;
  delete fModelTotpbM2;
  delete fModelTotcsM2;
  delete fModelTotco2M2;
  delete fModelTotteo2M2;

  delete fModelTotSthM2;
  delete fModelTotSuM2;
  delete fModelTotSpoM2;
  delete fModelTotSpbM2;
  delete fModelTotExtM2;

  delete fModelTotAdapM2;
  delete fModelTotAdapthM2;
  delete fModelTotAdapuM2;
  delete fModelTotAdapkM2;
  delete fModelTotAdapcoM2;
  delete fModelTotAdapmnM2;
  delete fModelTotAdapNDBDM2;
  delete fModelTotAdap2NDBDM2;
  delete fModelTotAdapbiM2;
  delete fModelTotAdapbi2M2;
  delete fModelTotAdapptM2;
  delete fModelTotAdappbM2;
  delete fModelTotAdapcsM2;
  delete fModelTotAdapco2M2;
  delete fModelTotAdapteo2M2;

  delete fModelTotAdapSthM2;
  delete fModelTotAdapSuM2;
  delete fModelTotAdapSpoM2;
  delete fModelTotAdapSpbM2;
  delete fModelTotAdapExtM2;

  // Total PDFs M2Sum
  delete fModelTotM2Sum;
  delete fModelTotthM2Sum;
  delete fModelTotuM2Sum;
  delete fModelTotkM2Sum;
  delete fModelTotcoM2Sum;
  delete fModelTotmnM2Sum;
  delete fModelTotNDBDM2Sum;
  delete fModelTot2NDBDM2Sum;
  delete fModelTotbiM2Sum;
  delete fModelTotbi2M2Sum;
  delete fModelTotptM2Sum;
  delete fModelTotpbM2Sum;
  delete fModelTotcsM2Sum;
  delete fModelTotco2M2Sum;
  delete fModelTotteo2M2Sum;

  delete fModelTotSthM2Sum;
  delete fModelTotSuM2Sum;
  delete fModelTotSpoM2Sum;
  delete fModelTotSpbM2Sum;
  delete fModelTotExtM2Sum;

  delete fModelTotAdapM2Sum;
  delete fModelTotAdapthM2Sum;
  delete fModelTotAdapuM2Sum;
  delete fModelTotAdapkM2Sum;
  delete fModelTotAdapcoM2Sum;
  delete fModelTotAdapmnM2Sum;
  delete fModelTotAdapNDBDM2Sum;
  delete fModelTotAdap2NDBDM2Sum;
  delete fModelTotAdapbiM2Sum;
  delete fModelTotAdapbi2M2Sum;
  delete fModelTotAdapptM2Sum;
  delete fModelTotAdappbM2Sum;
  delete fModelTotAdapcsM2Sum;
  delete fModelTotAdapco2M2Sum;
  delete fModelTotAdapteo2M2Sum;

  delete fModelTotAdapSthM2Sum;
  delete fModelTotAdapSuM2Sum;
  delete fModelTotAdapSpoM2Sum;
  delete fModelTotAdapSpbM2Sum;
  delete fModelTotAdapExtM2Sum;

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

  delete hCuFrameco58M1;
  delete hCuFrameco60M1;
  delete hCuFramecs137M1;
  delete hCuFramek40M1;
  delete hCuFramemn54M1;
  delete hCuFramepb210M1;
  delete hCuFrameth232M1;
  delete hCuFrameu238M1;

  delete hCuFrameSth232M1_1;
  delete hCuFrameSu238M1_1;
  delete hCuFrameSxpb210M1_001;
  delete hCuFrameSxpb210M1_01;
  delete hCuFrameSxpb210M1_1;
  delete hCuFrameSxpb210M1_10;
  delete hCuFrameSxth232M1_001;
  delete hCuFrameSxth232M1_01;
  delete hCuFrameSxth232M1_1;
  delete hCuFrameSxth232M1_10;
  delete hCuFrameSxu238M1_001;
  delete hCuFrameSxu238M1_01;
  delete hCuFrameSxu238M1_1;
  delete hCuFrameSxu238M1_10;

  delete hCuFrameco58M2;
  delete hCuFrameco60M2;
  delete hCuFramecs137M2;
  delete hCuFramek40M2;
  delete hCuFramemn54M2;
  delete hCuFramepb210M2;
  delete hCuFrameth232M2;
  delete hCuFrameu238M2;

  delete hCuFrameSth232M2_1;
  delete hCuFrameSu238M2_1;
  delete hCuFrameSxpb210M2_001;
  delete hCuFrameSxpb210M2_01;
  delete hCuFrameSxpb210M2_1;
  delete hCuFrameSxpb210M2_10;
  delete hCuFrameSxth232M2_001;
  delete hCuFrameSxth232M2_01;
  delete hCuFrameSxth232M2_1;
  delete hCuFrameSxth232M2_10;
  delete hCuFrameSxu238M2_001;
  delete hCuFrameSxu238M2_01;
  delete hCuFrameSxu238M2_1;
  delete hCuFrameSxu238M2_10;

  delete hCuFrameco58M2Sum;
  delete hCuFrameco60M2Sum;
  delete hCuFramecs137M2Sum;
  delete hCuFramek40M2Sum;
  delete hCuFramemn54M2Sum;
  delete hCuFramepb210M2Sum;
  delete hCuFrameth232M2Sum;
  delete hCuFrameu238M2Sum;

  delete hCuFrameSth232M2Sum_1;
  delete hCuFrameSu238M2Sum_1;
  delete hCuFrameSxpb210M2Sum_001;
  delete hCuFrameSxpb210M2Sum_01;
  delete hCuFrameSxpb210M2Sum_1;
  delete hCuFrameSxpb210M2Sum_10;
  delete hCuFrameSxth232M2Sum_001;
  delete hCuFrameSxth232M2Sum_01;
  delete hCuFrameSxth232M2Sum_1;
  delete hCuFrameSxth232M2Sum_10;
  delete hCuFrameSxu238M2Sum_001;
  delete hCuFrameSxu238M2Sum_01;
  delete hCuFrameSxu238M2Sum_1;
  delete hCuFrameSxu238M2Sum_10;

  delete hCuBoxco58M1;
  delete hCuBoxco60M1;
  delete hCuBoxcs137M1;
  delete hCuBoxk40M1;
  delete hCuBoxmn54M1;
  delete hCuBoxpb210M1;
  delete hCuBoxth232M1;
  delete hCuBoxu238M1;  

  delete hCuBoxSth232M1_1;
  delete hCuBoxSu238M1_1;
  delete hCuBoxSxpb210M1_001;
  delete hCuBoxSxpb210M1_01;
  delete hCuBoxSxpb210M1_1;
  delete hCuBoxSxpb210M1_10;
  delete hCuBoxSxth232M1_001;
  delete hCuBoxSxth232M1_01;
  delete hCuBoxSxth232M1_1;
  delete hCuBoxSxth232M1_10;
  delete hCuBoxSxu238M1_001;
  delete hCuBoxSxu238M1_01;
  delete hCuBoxSxu238M1_1;
  delete hCuBoxSxu238M1_10;

  delete hCuBoxco58M2;
  delete hCuBoxco60M2;
  delete hCuBoxcs137M2;
  delete hCuBoxk40M2;
  delete hCuBoxmn54M2;
  delete hCuBoxpb210M2;
  delete hCuBoxth232M2;
  delete hCuBoxu238M2;  

  delete hCuBoxSth232M2_1;
  delete hCuBoxSu238M2_1;
  delete hCuBoxSxpb210M2_001;
  delete hCuBoxSxpb210M2_01;
  delete hCuBoxSxpb210M2_1;
  delete hCuBoxSxpb210M2_10;
  delete hCuBoxSxth232M2_001;
  delete hCuBoxSxth232M2_01;
  delete hCuBoxSxth232M2_1;
  delete hCuBoxSxth232M2_10;
  delete hCuBoxSxu238M2_001;
  delete hCuBoxSxu238M2_01;
  delete hCuBoxSxu238M2_1;
  delete hCuBoxSxu238M2_10;

  delete hCuBoxco58M2Sum;
  delete hCuBoxco60M2Sum;
  delete hCuBoxcs137M2Sum;
  delete hCuBoxk40M2Sum;
  delete hCuBoxmn54M2Sum;
  delete hCuBoxpb210M2Sum;
  delete hCuBoxth232M2Sum;
  delete hCuBoxu238M2Sum; 

  delete hCuBoxSth232M2Sum_1;
  delete hCuBoxSu238M2Sum_1;
  delete hCuBoxSxpb210M2Sum_001;
  delete hCuBoxSxpb210M2Sum_01;
  delete hCuBoxSxpb210M2Sum_1;
  delete hCuBoxSxpb210M2Sum_10;
  delete hCuBoxSxth232M2Sum_001;
  delete hCuBoxSxth232M2Sum_01;
  delete hCuBoxSxth232M2Sum_1;
  delete hCuBoxSxth232M2Sum_10;
  delete hCuBoxSxu238M2Sum_001;
  delete hCuBoxSxu238M2Sum_01;
  delete hCuBoxSxu238M2Sum_1;
  delete hCuBoxSxu238M2Sum_10;

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

  delete h50mKco58M2Sum;
  delete h50mKco60M2Sum;
  delete h50mKcs137M2Sum;
  delete h50mKk40M2Sum;
  delete h50mKmn54M2Sum;
  delete h50mKpb210M2Sum;
  delete h50mKth232M2Sum;
  delete h50mKu238M2Sum;

  delete h600mKco60M1;
  delete h600mKk40M1;
  delete h600mKth232M1;
  delete h600mKu238M1;    

  delete h600mKco60M2;
  delete h600mKk40M2;
  delete h600mKth232M2;
  delete h600mKu238M2;  

  delete h600mKco60M2Sum;
  delete h600mKk40M2Sum;
  delete h600mKth232M2Sum;
  delete h600mKu238M2Sum; 

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

  delete hMBco60M1;
  delete hMBk40M1;
  delete hMBth232M1;
  delete hMBu238M1;   

  delete hMBco60M2;
  delete hMBk40M2;
  delete hMBth232M2;
  delete hMBu238M2; 

  delete hMBco60M2Sum;
  delete hMBk40M2Sum;
  delete hMBth232M2Sum;
  delete hMBu238M2Sum;  

  delete hSIk40M1;
  delete hSIth232M1;
  delete hSIu238M1;

  delete hSIk40M2;
  delete hSIth232M2;
  delete hSIu238M2;

  delete hSIk40M2Sum;
  delete hSIth232M2Sum;
  delete hSIu238M2Sum;

  delete hExtPbbi210M1;

  delete hExtPbbi210M2;
  
  delete hExtPbbi210M2Sum;

  delete hIVCco60M1;
  delete hIVCk40M1;
  delete hIVCth232M1;
  delete hIVCu238M1;    

  delete hIVCco60M2;
  delete hIVCk40M2;
  delete hIVCth232M2;
  delete hIVCu238M2;  

  delete hIVCco60M2Sum;
  delete hIVCk40M2Sum;
  delete hIVCth232M2Sum;
  delete hIVCu238M2Sum; 

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

  delete hAdapCuFrameco58M1;
  delete hAdapCuFrameco60M1;
  delete hAdapCuFramecs137M1;
  delete hAdapCuFramek40M1;
  delete hAdapCuFramemn54M1;
  delete hAdapCuFramepb210M1;
  delete hAdapCuFrameth232M1;
  delete hAdapCuFrameu238M1;

  delete hAdapCuFrameSth232M1_1;
  delete hAdapCuFrameSu238M1_1;
  delete hAdapCuFrameSxpb210M1_001;
  delete hAdapCuFrameSxpb210M1_01;
  delete hAdapCuFrameSxpb210M1_1;
  delete hAdapCuFrameSxpb210M1_10;
  delete hAdapCuFrameSxth232M1_001;
  delete hAdapCuFrameSxth232M1_01;
  delete hAdapCuFrameSxth232M1_1;
  delete hAdapCuFrameSxth232M1_10;
  delete hAdapCuFrameSxu238M1_001;
  delete hAdapCuFrameSxu238M1_01;
  delete hAdapCuFrameSxu238M1_1;
  delete hAdapCuFrameSxu238M1_10;

  delete hAdapCuFrameco58M2;
  delete hAdapCuFrameco60M2;
  delete hAdapCuFramecs137M2;
  delete hAdapCuFramek40M2;
  delete hAdapCuFramemn54M2;
  delete hAdapCuFramepb210M2;
  delete hAdapCuFrameth232M2;
  delete hAdapCuFrameu238M2;

  delete hAdapCuFrameSth232M2_1;
  delete hAdapCuFrameSu238M2_1;
  delete hAdapCuFrameSxpb210M2_001;
  delete hAdapCuFrameSxpb210M2_01;
  delete hAdapCuFrameSxpb210M2_1;
  delete hAdapCuFrameSxpb210M2_10;
  delete hAdapCuFrameSxth232M2_001;
  delete hAdapCuFrameSxth232M2_01;
  delete hAdapCuFrameSxth232M2_1;
  delete hAdapCuFrameSxth232M2_10;
  delete hAdapCuFrameSxu238M2_001;
  delete hAdapCuFrameSxu238M2_01;
  delete hAdapCuFrameSxu238M2_1;
  delete hAdapCuFrameSxu238M2_10;

  delete hAdapCuFrameco58M2Sum;
  delete hAdapCuFrameco60M2Sum;
  delete hAdapCuFramecs137M2Sum;
  delete hAdapCuFramek40M2Sum;
  delete hAdapCuFramemn54M2Sum;
  delete hAdapCuFramepb210M2Sum;
  delete hAdapCuFrameth232M2Sum;
  delete hAdapCuFrameu238M2Sum;

  delete hAdapCuFrameSth232M2Sum_1;
  delete hAdapCuFrameSu238M2Sum_1;
  delete hAdapCuFrameSxpb210M2Sum_001;
  delete hAdapCuFrameSxpb210M2Sum_01;
  delete hAdapCuFrameSxpb210M2Sum_1;
  delete hAdapCuFrameSxpb210M2Sum_10;
  delete hAdapCuFrameSxth232M2Sum_001;
  delete hAdapCuFrameSxth232M2Sum_01;
  delete hAdapCuFrameSxth232M2Sum_1;
  delete hAdapCuFrameSxth232M2Sum_10;
  delete hAdapCuFrameSxu238M2Sum_001;
  delete hAdapCuFrameSxu238M2Sum_01;
  delete hAdapCuFrameSxu238M2Sum_1;
  delete hAdapCuFrameSxu238M2Sum_10;

  delete hAdapCuBoxco58M1;
  delete hAdapCuBoxco60M1;
  delete hAdapCuBoxcs137M1;
  delete hAdapCuBoxk40M1;
  delete hAdapCuBoxmn54M1;
  delete hAdapCuBoxpb210M1;
  delete hAdapCuBoxth232M1;
  delete hAdapCuBoxu238M1;  

  delete hAdapCuBoxSth232M1_1;
  delete hAdapCuBoxSu238M1_1;
  delete hAdapCuBoxSxpb210M1_001;
  delete hAdapCuBoxSxpb210M1_01;
  delete hAdapCuBoxSxpb210M1_1;
  delete hAdapCuBoxSxpb210M1_10;
  delete hAdapCuBoxSxth232M1_001;
  delete hAdapCuBoxSxth232M1_01;
  delete hAdapCuBoxSxth232M1_1;
  delete hAdapCuBoxSxth232M1_10;
  delete hAdapCuBoxSxu238M1_001;
  delete hAdapCuBoxSxu238M1_01;
  delete hAdapCuBoxSxu238M1_1;
  delete hAdapCuBoxSxu238M1_10;

  delete hAdapCuBoxco58M2;
  delete hAdapCuBoxco60M2;
  delete hAdapCuBoxcs137M2;
  delete hAdapCuBoxk40M2;
  delete hAdapCuBoxmn54M2;
  delete hAdapCuBoxpb210M2;
  delete hAdapCuBoxth232M2;
  delete hAdapCuBoxu238M2;  

  delete hAdapCuBoxSth232M2_1;
  delete hAdapCuBoxSu238M2_1;
  delete hAdapCuBoxSxpb210M2_001;
  delete hAdapCuBoxSxpb210M2_01;
  delete hAdapCuBoxSxpb210M2_1;
  delete hAdapCuBoxSxpb210M2_10;
  delete hAdapCuBoxSxth232M2_001;
  delete hAdapCuBoxSxth232M2_01;
  delete hAdapCuBoxSxth232M2_1;
  delete hAdapCuBoxSxth232M2_10;
  delete hAdapCuBoxSxu238M2_001;
  delete hAdapCuBoxSxu238M2_01;
  delete hAdapCuBoxSxu238M2_1;
  delete hAdapCuBoxSxu238M2_10;

  delete hAdapCuBoxco58M2Sum;
  delete hAdapCuBoxco60M2Sum;
  delete hAdapCuBoxcs137M2Sum;
  delete hAdapCuBoxk40M2Sum;
  delete hAdapCuBoxmn54M2Sum;
  delete hAdapCuBoxpb210M2Sum;
  delete hAdapCuBoxth232M2Sum;
  delete hAdapCuBoxu238M2Sum; 

  delete hAdapCuBoxSth232M2Sum_1;
  delete hAdapCuBoxSu238M2Sum_1;
  delete hAdapCuBoxSxpb210M2Sum_001;
  delete hAdapCuBoxSxpb210M2Sum_01;
  delete hAdapCuBoxSxpb210M2Sum_1;
  delete hAdapCuBoxSxpb210M2Sum_10;
  delete hAdapCuBoxSxth232M2Sum_001;
  delete hAdapCuBoxSxth232M2Sum_01;
  delete hAdapCuBoxSxth232M2Sum_1;
  delete hAdapCuBoxSxth232M2Sum_10;
  delete hAdapCuBoxSxu238M2Sum_001;
  delete hAdapCuBoxSxu238M2Sum_01;
  delete hAdapCuBoxSxu238M2Sum_1;
  delete hAdapCuBoxSxu238M2Sum_10;

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

  delete hAdap50mKco58M2Sum;
  delete hAdap50mKco60M2Sum;
  delete hAdap50mKcs137M2Sum;
  delete hAdap50mKk40M2Sum;
  delete hAdap50mKmn54M2Sum;
  delete hAdap50mKpb210M2Sum;
  delete hAdap50mKth232M2Sum;
  delete hAdap50mKu238M2Sum;  

  delete hAdap600mKco60M1;
  delete hAdap600mKk40M1;
  delete hAdap600mKth232M1;
  delete hAdap600mKu238M1;    

  delete hAdap600mKco60M2;
  delete hAdap600mKk40M2;
  delete hAdap600mKth232M2;
  delete hAdap600mKu238M2;  

  delete hAdap600mKco60M2Sum;
  delete hAdap600mKk40M2Sum;
  delete hAdap600mKth232M2Sum;
  delete hAdap600mKu238M2Sum; 

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

  delete hAdapMBco60M1;
  delete hAdapMBk40M1;
  delete hAdapMBth232M1;
  delete hAdapMBu238M1;   

  delete hAdapMBco60M2;
  delete hAdapMBk40M2;
  delete hAdapMBth232M2;
  delete hAdapMBu238M2; 

  delete hAdapMBco60M2Sum;
  delete hAdapMBk40M2Sum;
  delete hAdapMBth232M2Sum;
  delete hAdapMBu238M2Sum;  

  delete hAdapSIk40M1;
  delete hAdapSIth232M1;
  delete hAdapSIu238M1;

  delete hAdapSIk40M2;
  delete hAdapSIth232M2;
  delete hAdapSIu238M2;

  delete hAdapSIk40M2Sum;
  delete hAdapSIth232M2Sum;
  delete hAdapSIu238M2Sum;

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

  delete hAdapIVCco60M1;
  delete hAdapIVCk40M1;
  delete hAdapIVCth232M1;
  delete hAdapIVCu238M1;    

  delete hAdapIVCco60M2;
  delete hAdapIVCk40M2;
  delete hAdapIVCth232M2;
  delete hAdapIVCu238M2;  

  delete hAdapIVCco60M2Sum;
  delete hAdapIVCk40M2Sum;
  delete hAdapIVCth232M2Sum;
  delete hAdapIVCu238M2Sum; 

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
  // 25 since start at 50 keV with 2 keV bin width
  for(int i = h1->FindBin(0); i < h1->FindBin(50); i++)
  {
    dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(i));
  }
  for(int i = h1->FindBin(50); i < h1->GetNbinsX(); i++)
  {
    // Added per each peak

    // K40
    if(i >= h1->FindBin(1440) && i < h1->FindBin(1480))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(1440)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(1480-1440)/dBaseBinSize;
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


    // Pt peak 3150 - 3400
    if(i >= h1->FindBin(3150) && i < h1->FindBin(3400))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(3150)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(3400-3150)/dBaseBinSize;
    }
    // 4050 - 4150
    if(i >= h1->FindBin(4050) && i < h1->FindBin(4150))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4050)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(4150-4050)/dBaseBinSize;
    }
    // 4200 - 4350
    if(i >= h1->FindBin(4200) && i < h1->FindBin(4350))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4200)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(4350-4200)/dBaseBinSize;
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
    if(i >= h1->FindBin(4850) && i < h1->FindBin(4950))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4850)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(4950-4850)/dBaseBinSize;
    }

    // 5200 - 5400
    if(i >= h1->FindBin(5200) && i < h1->FindBin(5400))
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

    // 5800 - 6050
    if(i >= h1->FindBin(5800) && i < h1->FindBin(6050))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5800)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(6050-5800)/dBaseBinSize;
    }    

    // 6050 - 6350
    if(i >= h1->FindBin(6050) && i < h1->FindBin(6350))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(6050)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(6350-6050)/dBaseBinSize;
    }    

    // 6700 - 6900
    if(i >= h1->FindBin(6700) && i < h1->FindBin(6900))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(6700)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(6900-6700)/dBaseBinSize;
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

  // Residual plot and distribution
  for(int j = binMin ; j < binMax ; j++)
  {

    if( hCloneBkg->GetBinCenter(j) >= 3150 && hCloneBkg->GetBinCenter(j) <= 3400)continue;    
    // if( hCloneBkg->GetBinCenter(j) >= 658 && hCloneBkg->GetBinCenter(j) <= 664)continue;
    // if( hCloneBkg->GetBinCenter(j) >= 800 && hCloneBkg->GetBinCenter(j) <= 808)continue;
    // if( hCloneBkg->GetBinCenter(j) >= 1060 && hCloneBkg->GetBinCenter(j) <= 1068)continue;

    dResidualX    = hCloneBkg->GetBinCenter(j);
    
    // Re-multiply bin content by bin width (for # of counts)
    
    dResidualY    = (hCloneBkg->GetBinContent(j) - hCloneMC->GetBinContent(j))*hCloneBkg->GetBinWidth(j) /
              (TMath::Sqrt((hCloneBkg->GetBinContent(j))*hCloneBkg->GetBinWidth(j)) ); // Sqrt(MC + data) = sigma for poisson distribution
    
    // dResidualY    = (hCloneBkg->GetBinContent(j) - hCloneMC->GetBinContent(j))*hCloneBkg->GetBinWidth(j) /
              // (TMath::Sqrt(hCloneBkg->GetBinContent(i))/hCloneBkg->GetBinWidth(i) ); // Sqrt(MC + data) = sigma for poisson distribution
    
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

/*
  TCanvas *cBkg1 = new TCanvas("cBkg1", "cBkg1", 1200, 1200);
  cBkg1->SetLogy();

  int dRebin = 1;

  fDataHistoTot->SetLineColor(kBlack);
  fDataHistoM1->SetLineColor(kBlue);
  fDataHistoM2->SetLineColor(kRed);
  fDataHistoM2Sum->SetLineColor(kGreen+1);

  // fDataHistoTot->Rebin(dRebin);
  // fDataHistoM1->Rebin(dRebin);
  // fDataHistoM2->Rebin(dRebin);
  // fDataHistoM2Sum->Rebin(dRebin);

  fDataHistoTot->GetXaxis()->SetTitle("Energy (keV)");
  fDataHistoTot->GetYaxis()->SetTitle(Form("Counts/(%d keV)/yr", dBinSize*dRebin));
  fDataHistoTot->Draw("E0");
  fDataHistoM1->Draw("E0SAME");
  fDataHistoM2->Draw("E0SAME");
  fDataHistoM2Sum->Draw("E0SAME");
*/

  TLegend *leg = new TLegend(0.75,0.75,0.97,0.97);
  leg->AddEntry(fDataHistoTot, "Total", "l");
  leg->AddEntry(fDataHistoM1, "M1", "l");
  leg->AddEntry(fDataHistoM2, "M2", "l");
  leg->Draw();

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1800);
  c1->Divide(1, 2);
  c1->cd(1);
  gPad->SetLogy();
  fAdapDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM1->GetYaxis()->SetTitle("Counts/bin");  
  fAdapDataHistoM1->Draw("");


  c1->cd(2); 
  gPad->SetLogy();
  fAdapDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM2->GetYaxis()->SetTitle("Counts/bin");   
  fAdapDataHistoM2->Draw("");

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
  double datam2sum_i, errm2sum_i;
  double modelm1_i, modelm2_i, modelm2sum_i;

  for(int i = dFitMinBinM1 ; i <= dFitMaxBinM1; i++)
  {

    // Dividing by base bin size in chi-squared because the weight is width/base bin size when filling
    if( fAdapDataHistoM1->GetBinCenter(i) >= 3150 && fAdapDataHistoM1->GetBinCenter(i) <= 3400)continue;
    // Skipping unknown peaks
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 658 && fAdapDataHistoM1->GetBinCenter(i) <= 664)continue;
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 800 && fAdapDataHistoM1->GetBinCenter(i) <= 808)continue;
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 1060 && fAdapDataHistoM1->GetBinCenter(i) <= 1068)continue;  

    datam1_i = fAdapDataHistoM1->GetBinContent(i)*fAdapDataHistoM1->GetBinWidth(i);
    modelm1_i = fModelTotAdapM1->GetBinContent(i)*fAdapDataHistoM1->GetBinWidth(i);

    if(modelm1_i != 0 && datam1_i != 0)
    {
      chiSquare += 2 * (modelm1_i - datam1_i + datam1_i * TMath::Log(datam1_i/modelm1_i));
    }
  }



  for(int i = dFitMinBinM2; i <= dFitMaxBinM2; i++)
  {
    if( fAdapDataHistoM2->GetBinCenter(i) >= 3150 && fAdapDataHistoM2->GetBinCenter(i) <= 3400)continue;
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 658 && fAdapDataHistoM2->GetBinCenter(i) <= 664)continue;
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 800 && fAdapDataHistoM2->GetBinCenter(i) <= 808)continue;
    // if( fAdapDataHistoM2->GetBinCenter(i) >= 1060 && fAdapDataHistoM2->GetBinCenter(i) <= 1068)continue;     

    datam2_i = fAdapDataHistoM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i);
    modelm2_i = fModelTotAdapM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i);

    if(modelm2_i != 0 && datam2_i != 0)
    {
      chiSquare += 2 * (modelm2_i - datam2_i + datam2_i * TMath::Log(datam2_i/modelm2_i));
    }

  }

  /*
  // Penalty terms for alpha region
  // Th232 term (teo2 th232 only, th232 0.01, th232 only 0.01, ra228-pb208; cubox+cuframe th232 10)
  double dP1 = fParameters[5];
  double dC1 = 0.0002930699;
  double dCE1 = TMath::Sqrt(0.0000567116*0.0000567116);
  chiSquare += (dP1 - dC1)(dP1 - dC1)/(dCE1*dCE1);

  double dP2 = fParameters[7];
  double dC2 = 0.0016913862;
  double dCE2 = TMath::Sqrt(0.0004400189*0.0004400189);
  chiSquare += (dP2 - dC2)(dP2 - dC2)/(dCE2*dCE2);

  double dP3 = fParameters[9];
  double dC3 = 0.0025633479;
  double dCE3 = TMath::Sqrt(0.0004222949*0.0004222949);
  chiSquare += (dP3 - dC3)(dP3 - dC3)/(dCE3*dCE3);

  double dParaCuTh232 = fParameters[15];
  double dConsCuTh232 = 0.0143128538;
  double dConsErrCuTh232 = 0.0025322971;
  chiSquare += (dParaCuTh232 - dConsCuTh232)(dParaCuTh232 - dConsCuTh232)/(dConsErrCuTh232*dConsErrCuTh232);

  // U238 term (teo2 th230 only, u238-th230 0.01, teO2 th232 only 0.01, ra226-pb210 0.01, cubox+cuframe u238 10)
  double dP5 = fParameters[6];
  double dC5 = 0.0002989930;
  double dCE5 = TMath::Sqrt( 0.0001164103*0.0001164103 );
  chiSquare += (dP5 - dC5)(dP5 - dC5)/(dCE5*dCE5);

  double dP6 = fParameters[10];
  double dC6 = 0.0017130394;
  double dCE6 = TMath::Sqrt( 0.0001216715*0.0001216715 );
  chiSquare += (dP6 - dC6)(dP6 - dC6)/(dCE6*dCE6);

  double dP7 = fParameters[11];
  double dC7 = 0.0007863522;
  double dCE7 = TMath::Sqrt( 0.0001648186*0.0001648186 );
  chiSquare += (dP7 - dC7)(dP7 - dC7)/(dCE7*dCE7);

  double dP8 = fParameters[12];
  double dC8 = 0.0030621684;
  double dCE8 = TMath::Sqrt( 0.0001819745*0.0001819745 );
  chiSquare += (dP8 - dC8)(dP8 - dC8)/(dCE8*dCE8);

  double dParaCuU238 = fParameters[16];
  double dConsCuU238 = 0.0012995965;
  double dConsErrCuU238 = 0.0030367388;
  chiSquare += (dParaCuU238 - dConsCuU238)(dParaCuU238 - dConsCuU238)/(dConsErrCuU238*dConsErrCuU238);


  // TeO2 Po210 term (teo2 pb210 1, pb210 0.01; cubox+cuframe pb210 10, pb210 0.1, pb210 0.01)
  double dParaTeO21Pb210 = fParameters[13];
  double dConsTeO21Pb210 = 0.0054316786;
  double dConsErrTeO21Pb210 = 0.0004570424;
  chiSquare += (dParaTeO21Pb210 - dConsTeO21Pb210)(dParaTeO21Pb210 - dConsTeO21Pb210)/(dConsErrTeO21Pb210*dConsErrTeO21Pb210);

  double dParaTeO22Pb210 = fParameters[14];
  double dConsTeO22Pb210 = 0.0420244871;
  double dConsErrTeO22Pb210 = 0.0005936196;
  chiSquare += (dParaTeO22Pb210 - dConsTeO22Pb210)(dParaTeO22Pb210 - dConsTeO22Pb210)/(dConsErrTeO22Pb210*dConsErrTeO22Pb210);

  // CuBox + CuFrame Po210 term
  double dP9 = fParameters[17];
  double dC9 = 0.0043300304;
  double dCE9 = TMath::Sqrt(0.0019620557*0.0019620557);
  chiSquare += (dP9 - dC9)(dP9 - dC9)/(dCE9*dCE9);

  double dP10 = fParameters[18];
  double dC10 = 0.0044207691;
  double dCE10 = TMath::Sqrt(0.0008941607*0.0008941607);
  chiSquare += (dP10 - dC10)(dP10 - dC10)/(dCE10*dCE10);

  double dP11 = fParameters[19];
  double dC11 = 0.0181033015;
  double dCE11 = TMath::Sqrt(0.0008256090*0.0008256090);
  chiSquare += (dP11 - dC11)(dP11 - dC11)/(dCE11*dCE11);
  */

  return chiSquare;
}

void TBackgroundModel::Initialize()
{	

  // Loads PDFs from file
  cout << "Loading PDF Histograms from file" << endl;
  cout << "Directory " << Form("%s/MCProduction", dMCDir.c_str()) << endl;

  fBulkInner = new TFile(Form("%s/MCProduction_BulkInner_1keV.root", dMCDir.c_str()));
  fBulkInnerM2Sum = new TFile(Form("%s/MCProduction_BulkInnerM2Sum_1keV.root", dMCDir.c_str()));

  fBulkOuter = new TFile(Form("%s/MCProduction_BulkOuter_1keV.root", dMCDir.c_str()));
  fBulkOuterM2Sum = new TFile(Form("%s/MCProduction_BulkOuterM2Sum_1keV.root", dMCDir.c_str()));

  fSurfaceCrystal = new TFile(Form("%s/MCProduction_SurfaceCrystal_1keV.root", dMCDir.c_str()));
  fSurfaceOther = new TFile(Form("%s/MCProduction_SurfaceOther_1keV.root", dMCDir.c_str()));

///////////// Bulk Histograms
/////// Crystal M1 and M2
  hTeO20nuM1     = (TH1D*)fBulkInner->Get("hTeO20nuM1");
  hTeO22nuM1     = (TH1D*)fBulkInner->Get("hTeO22nuM1");
  hTeO2co60M1    = (TH1D*)fBulkInner->Get("hTeO2co60M1");
  hTeO2k40M1     = (TH1D*)fBulkInner->Get("hTeO2k40M1");
  hTeO2pb210M1   = (TH1D*)fBulkInner->Get("hTeO2pb210M1");
  hTeO2po210M1   = (TH1D*)fBulkInner->Get("hTeO2po210M1");
  hTeO2te125M1   = (TH1D*)fBulkInner->Get("hTeO2te125M1");
  hTeO2th232M1   = (TH1D*)fBulkInner->Get("hTeO2th232M1");
  // hTeO2th228M1   = (TH1D*)fBulkInner->Get("hTeO2th228M1");
  // hTeO2ra226M1   = (TH1D*)fBulkInner->Get("hTeO2ra226M1");
  // hTeO2rn222M1   = (TH1D*)fBulkInner->Get("hTeO2rn222M1");
  hTeO2u238M1    = (TH1D*)fBulkInner->Get("hTeO2u238M1");
  // hTeO2th230M1   = (TH1D*)fBulkInner->Get("hTeO2th230M1");
  // hTeO2u234M1    = (TH1D*)fBulkInner->Get("hTeO2u234M1");

  hTeO2th232onlyM1 = (TH1D*)fBulkInner->Get("hTeO2th232onlyM1");
  hTeO2ra228pb208M1 = (TH1D*)fBulkInner->Get("hTeO2ra228pb208M1");
  hTeO2th230onlyM1 = (TH1D*)fBulkInner->Get("hTeO2th230onlyM1");

  hTeO20nuM2     = (TH1D*)fBulkInner->Get("hTeO20nuM2");
  hTeO22nuM2     = (TH1D*)fBulkInner->Get("hTeO22nuM2");
  hTeO2co60M2    = (TH1D*)fBulkInner->Get("hTeO2co60M2");
  hTeO2k40M2     = (TH1D*)fBulkInner->Get("hTeO2k40M2");
  hTeO2pb210M2   = (TH1D*)fBulkInner->Get("hTeO2pb210M2");
  hTeO2po210M2   = (TH1D*)fBulkInner->Get("hTeO2po210M2");
  hTeO2te125M2   = (TH1D*)fBulkInner->Get("hTeO2te125M2");
  hTeO2th232M2   = (TH1D*)fBulkInner->Get("hTeO2th232M2");
  // hTeO2th228M2   = (TH1D*)fBulkInner->Get("hTeO2th228M2");
  // hTeO2ra226M2   = (TH1D*)fBulkInner->Get("hTeO2ra226M2");
  // hTeO2rn222M2   = (TH1D*)fBulkInner->Get("hTeO2rn222M2");
  hTeO2u238M2    = (TH1D*)fBulkInner->Get("hTeO2u238M2");
  // hTeO2th230M2   = (TH1D*)fBulkInner->Get("hTeO2th230M2");
  // hTeO2u234M2    = (TH1D*)fBulkInner->Get("hTeO2u234M2");

  hTeO2th232onlyM2 = (TH1D*)fBulkInner->Get("hTeO2th232onlyM2");
  hTeO2ra228pb208M2 = (TH1D*)fBulkInner->Get("hTeO2ra228pb208M2");
  hTeO2th230onlyM2 = (TH1D*)fBulkInner->Get("hTeO2th230onlyM2");

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


///////// Frame M1 and M2
  hCuFrameco58M1    = (TH1D*)fBulkInner->Get("hCuFrameco58M1");
  hCuFrameco60M1    = (TH1D*)fBulkInner->Get("hCuFrameco60M1");
  hCuFramecs137M1   = (TH1D*)fBulkInner->Get("hCuFramecs137M1");
  hCuFramek40M1     = (TH1D*)fBulkInner->Get("hCuFramek40M1");
  hCuFramemn54M1    = (TH1D*)fBulkInner->Get("hCuFramemn54M1");
  hCuFramepb210M1   = (TH1D*)fBulkInner->Get("hCuFramepb210M1");
  hCuFrameth232M1   = (TH1D*)fBulkInner->Get("hCuFrameth232M1");
  hCuFrameu238M1    = (TH1D*)fBulkInner->Get("hCuFrameu238M1");
   
  hCuFrameco58M2    = (TH1D*)fBulkInner->Get("hCuFrameco58M2");
  hCuFrameco60M2    = (TH1D*)fBulkInner->Get("hCuFrameco60M2");
  hCuFramecs137M2   = (TH1D*)fBulkInner->Get("hCuFramecs137M2");
  hCuFramek40M2     = (TH1D*)fBulkInner->Get("hCuFramek40M2");
  hCuFramemn54M2    = (TH1D*)fBulkInner->Get("hCuFramemn54M2");
  hCuFramepb210M2   = (TH1D*)fBulkInner->Get("hCuFramepb210M2");
  hCuFrameth232M2   = (TH1D*)fBulkInner->Get("hCuFrameth232M2");
  hCuFrameu238M2    = (TH1D*)fBulkInner->Get("hCuFrameu238M2");

  hCuFrameco58M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuFrameco58M2Sum");
  hCuFrameco60M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuFrameco60M2Sum");
  hCuFramecs137M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hCuFramecs137M2Sum");
  hCuFramek40M2Sum     = (TH1D*)fBulkInnerM2Sum->Get("hCuFramek40M2Sum");
  hCuFramemn54M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuFramemn54M2Sum");
  hCuFramepb210M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hCuFramepb210M2Sum");
  hCuFrameth232M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hCuFrameth232M2Sum");
  hCuFrameu238M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuFrameu238M2Sum");

///////// CuBox (TShield) M1 and M2
  hCuBoxco58M1    = (TH1D*)fBulkInner->Get("hCuBoxco58M1");
  hCuBoxco60M1    = (TH1D*)fBulkInner->Get("hCuBoxco60M1");
  hCuBoxcs137M1   = (TH1D*)fBulkInner->Get("hCuBoxcs137M1");
  hCuBoxk40M1     = (TH1D*)fBulkInner->Get("hCuBoxk40M1");
  hCuBoxmn54M1    = (TH1D*)fBulkInner->Get("hCuBoxmn54M1");
  hCuBoxpb210M1   = (TH1D*)fBulkInner->Get("hCuBoxpb210M1");
  hCuBoxth232M1   = (TH1D*)fBulkInner->Get("hCuBoxth232M1");
  hCuBoxu238M1    = (TH1D*)fBulkInner->Get("hCuBoxu238M1");
   
  hCuBoxco58M2    = (TH1D*)fBulkInner->Get("hCuBoxco58M2");
  hCuBoxco60M2    = (TH1D*)fBulkInner->Get("hCuBoxco60M2");
  hCuBoxcs137M2   = (TH1D*)fBulkInner->Get("hCuBoxcs137M2");
  hCuBoxk40M2     = (TH1D*)fBulkInner->Get("hCuBoxk40M2");
  hCuBoxmn54M2    = (TH1D*)fBulkInner->Get("hCuBoxmn54M2");
  hCuBoxpb210M2   = (TH1D*)fBulkInner->Get("hCuBoxpb210M2");
  hCuBoxth232M2   = (TH1D*)fBulkInner->Get("hCuBoxth232M2");
  hCuBoxu238M2    = (TH1D*)fBulkInner->Get("hCuBoxu238M2");

  hCuBoxco58M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxco58M2Sum");
  hCuBoxco60M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxco60M2Sum");
  hCuBoxcs137M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxcs137M2Sum");
  hCuBoxk40M2Sum     = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxk40M2Sum");
  hCuBoxmn54M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxmn54M2Sum");
  hCuBoxpb210M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxpb210M2Sum");
  hCuBoxth232M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxth232M2Sum");
  hCuBoxu238M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxu238M2Sum");

///////// CuBox + CuFrame M1 and M2
  hCuBox_CuFrameco60M1 = (TH1D*)fBulkInner->Get("hCuBox_CuFrameco60M1");
  hCuBox_CuFramek40M1 = (TH1D*)fBulkInner->Get("hCuBox_CuFramek40M1");
  hCuBox_CuFrameth232M1 = (TH1D*)fBulkInner->Get("hCuBox_CuFrameth232M1");
  hCuBox_CuFrameu238M1 = (TH1D*)fBulkInner->Get("hCuBox_CuFrameu238M1");

  hCuBox_CuFrameco60M2 = (TH1D*)fBulkInner->Get("hCuBox_CuFrameco60M2");
  hCuBox_CuFramek40M2 = (TH1D*)fBulkInner->Get("hCuBox_CuFramek40M2");
  hCuBox_CuFrameth232M2 = (TH1D*)fBulkInner->Get("hCuBox_CuFrameth232M2");
  hCuBox_CuFrameu238M2 = (TH1D*)fBulkInner->Get("hCuBox_CuFrameu238M2");

  hCuBox_CuFrameco60M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hCuBox_CuFrameco60M2Sum");
  hCuBox_CuFramek40M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hCuBox_CuFramek40M2Sum");
  hCuBox_CuFrameth232M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hCuBox_CuFrameth232M2Sum");
  hCuBox_CuFrameu238M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hCuBox_CuFrameu238M2Sum");

/*
///////// 50mK M1 and M2
  h50mKco58M1    = (TH1D*)fBulkOuter->Get("h50mKco58M1");
  h50mKco60M1    = (TH1D*)fBulkOuter->Get("h50mKco60M1");
  h50mKcs137M1   = (TH1D*)fBulkOuter->Get("h50mKcs137M1");
  h50mKk40M1     = (TH1D*)fBulkOuter->Get("h50mKk40M1");
  h50mKmn54M1    = (TH1D*)fBulkOuter->Get("h50mKmn54M1");
  h50mKpb210M1   = (TH1D*)fBulkOuter->Get("h50mKpb210M1");
  h50mKth232M1   = (TH1D*)fBulkOuter->Get("h50mKth232M1");
  h50mKu238M1    = (TH1D*)fBulkOuter->Get("h50mKu238M1");
   
  h50mKco58M2    = (TH1D*)fBulkOuter->Get("h50mKco58M2");
  h50mKco60M2    = (TH1D*)fBulkOuter->Get("h50mKco60M2");
  h50mKcs137M2   = (TH1D*)fBulkOuter->Get("h50mKcs137M2");
  h50mKk40M2     = (TH1D*)fBulkOuter->Get("h50mKk40M2");
  h50mKmn54M2    = (TH1D*)fBulkOuter->Get("h50mKmn54M2");
  h50mKpb210M2   = (TH1D*)fBulkOuter->Get("h50mKpb210M2");
  h50mKth232M2   = (TH1D*)fBulkOuter->Get("h50mKth232M2");
  h50mKu238M2    = (TH1D*)fBulkOuter->Get("h50mKu238M2");

  h50mKco58M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("h50mKco58M2Sum");
  h50mKco60M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("h50mKco60M2Sum");
  h50mKcs137M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("h50mKcs137M2Sum");
  h50mKk40M2Sum     = (TH1D*)fBulkOuterM2Sum->Get("h50mKk40M2Sum");
  h50mKmn54M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("h50mKmn54M2Sum");
  h50mKpb210M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("h50mKpb210M2Sum");
  h50mKth232M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("h50mKth232M2Sum");
  h50mKu238M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("h50mKu238M2Sum");

///////// 600mK M1 and M2
  h600mKco60M1    = (TH1D*)fBulkOuter->Get("h600mKco60M1");
  h600mKk40M1     = (TH1D*)fBulkOuter->Get("h600mKk40M1");
  h600mKth232M1   = (TH1D*)fBulkOuter->Get("h600mKth232M1");
  h600mKu238M1    = (TH1D*)fBulkOuter->Get("h600mKu238M1");
   
  h600mKco60M2    = (TH1D*)fBulkOuter->Get("h600mKco60M2");
  h600mKk40M2     = (TH1D*)fBulkOuter->Get("h600mKk40M2");
  h600mKth232M2   = (TH1D*)fBulkOuter->Get("h600mKth232M2");
  h600mKu238M2    = (TH1D*)fBulkOuter->Get("h600mKu238M2");

  h600mKco60M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("h600mKco60M2Sum");
  h600mKk40M2Sum     = (TH1D*)fBulkOuterM2Sum->Get("h600mKk40M2Sum");
  h600mKth232M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("h600mKth232M2Sum");
  h600mKu238M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("h600mKu238M2Sum");
*/
//////// Roman Lead M1 and M2
  // hPbRombi207M1   = (TH1D*)fBulkOuter->Get("hPbRombi207M1");
  hPbRomco60M1    = (TH1D*)fBulkOuter->Get("hPbRomco60M1");
  hPbRomcs137M1   = (TH1D*)fBulkOuter->Get("hPbRomcs137M1");
  hPbRomk40M1     = (TH1D*)fBulkOuter->Get("hPbRomk40M1");
  hPbRompb210M1   = (TH1D*)fBulkOuter->Get("hPbRompb210M1");
  hPbRomth232M1   = (TH1D*)fBulkOuter->Get("hPbRomth232M1");
  hPbRomu238M1    = (TH1D*)fBulkOuter->Get("hPbRomu238M1");
   
  // hPbRombi207M2   = (TH1D*)fBulkOuter->Get("hPbRombi207M2");
  hPbRomco60M2    = (TH1D*)fBulkOuter->Get("hPbRomco60M2");
  hPbRomcs137M2   = (TH1D*)fBulkOuter->Get("hPbRomcs137M2");
  hPbRomk40M2     = (TH1D*)fBulkOuter->Get("hPbRomk40M2");
  hPbRompb210M2   = (TH1D*)fBulkOuter->Get("hPbRompb210M2");
  hPbRomth232M2   = (TH1D*)fBulkOuter->Get("hPbRomth232M2");
  hPbRomu238M2    = (TH1D*)fBulkOuter->Get("hPbRomu238M2");

  // hPbRombi207M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("hPbRombi207M2Sum");
  hPbRomco60M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hPbRomco60M2Sum");
  hPbRomcs137M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("hPbRomcs137M2Sum");
  hPbRomk40M2Sum     = (TH1D*)fBulkOuterM2Sum->Get("hPbRomk40M2Sum");
  hPbRompb210M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("hPbRompb210M2Sum");
  hPbRomth232M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("hPbRomth232M2Sum");
  hPbRomu238M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hPbRomu238M2Sum");

/*
/////// Main Bath M1 and M2
  hMBco60M1    = (TH1D*)fBulkOuter->Get("hMBco60M1");
  hMBk40M1     = (TH1D*)fBulkOuter->Get("hMBk40M1");
  hMBth232M1   = (TH1D*)fBulkOuter->Get("hMBth232M1");
  hMBu238M1    = (TH1D*)fBulkOuter->Get("hMBu238M1");
   
  hMBco60M2    = (TH1D*)fBulkOuter->Get("hMBco60M2");
  hMBk40M2     = (TH1D*)fBulkOuter->Get("hMBk40M2");
  hMBth232M2   = (TH1D*)fBulkOuter->Get("hMBth232M2");
  hMBu238M2    = (TH1D*)fBulkOuter->Get("hMBu238M2");

  hMBco60M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hMBco60M2Sum");
  hMBk40M2Sum     = (TH1D*)fBulkOuterM2Sum->Get("hMBk40M2Sum");
  hMBth232M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("hMBth232M2Sum");
  hMBu238M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hMBu238M2Sum");
*/
////// Internal Shields M1 and M2
  hInternalco60M1 = (TH1D*)fBulkInner->Get("hInternalco60M1");
  hInternalk40M1 = (TH1D*)fBulkInner->Get("hInternalk40M1");
  hInternalth232M1 = (TH1D*)fBulkInner->Get("hInternalth232M1");
  hInternalu238M1 = (TH1D*)fBulkInner->Get("hInternalu238M1");

  hInternalco60M2 = (TH1D*)fBulkInner->Get("hInternalco60M2");
  hInternalk40M2 = (TH1D*)fBulkInner->Get("hInternalk40M2");
  hInternalth232M2 = (TH1D*)fBulkInner->Get("hInternalth232M2");
  hInternalu238M2 = (TH1D*)fBulkInner->Get("hInternalu238M2");

  hInternalco60M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hInternalco60M2Sum");
  hInternalk40M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hInternalk40M2Sum");
  hInternalth232M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hInternalth232M2Sum");
  hInternalu238M2Sum = (TH1D*)fBulkInnerM2Sum->Get("hInternalu238M2Sum");

/*
/////// Super Insulation M1 and M2
  hSIk40M1     = (TH1D*)fBulkOuter->Get("hSIk40M1");
  hSIth232M1   = (TH1D*)fBulkOuter->Get("hSIth232M1");
  hSIu238M1    = (TH1D*)fBulkOuter->Get("hSIu238M1");  

  hSIk40M2     = (TH1D*)fBulkOuter->Get("hSIk40M2");
  hSIth232M2   = (TH1D*)fBulkOuter->Get("hSIth232M2");
  hSIu238M2    = (TH1D*)fBulkOuter->Get("hSIu238M2"); 

  hSIk40M2Sum     = (TH1D*)fBulkOuterM2Sum->Get("hSIk40M2Sum");
  hSIth232M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("hSIth232M2Sum");
  hSIu238M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hSIu238M2Sum"); 

/////// IVC M1 and M2
  hIVCco60M1    = (TH1D*)fBulkOuter->Get("hIVCco60M1");
  hIVCk40M1     = (TH1D*)fBulkOuter->Get("hIVCk40M1");
  hIVCth232M1   = (TH1D*)fBulkOuter->Get("hIVCth232M1");
  hIVCu238M1    = (TH1D*)fBulkOuter->Get("hIVCu238M1");
   
  hIVCco60M2    = (TH1D*)fBulkOuter->Get("hIVCco60M2");
  hIVCk40M2     = (TH1D*)fBulkOuter->Get("hIVCk40M2");
  hIVCth232M2   = (TH1D*)fBulkOuter->Get("hIVCth232M2");
  hIVCu238M2    = (TH1D*)fBulkOuter->Get("hIVCu238M2");

  hIVCco60M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hIVCco60M2Sum");
  hIVCk40M2Sum     = (TH1D*)fBulkOuterM2Sum->Get("hIVCk40M2Sum");
  hIVCth232M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("hIVCth232M2Sum");
  hIVCu238M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hIVCu238M2Sum");
*/
/////// OVC M1 and M2
  hOVCco60M1    = (TH1D*)fBulkOuter->Get("hOVCco60M1");
  hOVCk40M1     = (TH1D*)fBulkOuter->Get("hOVCk40M1");
  hOVCth232M1   = (TH1D*)fBulkOuter->Get("hOVCth232M1");
  hOVCu238M1    = (TH1D*)fBulkOuter->Get("hOVCu238M1");
   
  hOVCco60M2    = (TH1D*)fBulkOuter->Get("hOVCco60M2");
  hOVCk40M2     = (TH1D*)fBulkOuter->Get("hOVCk40M2");
  hOVCth232M2   = (TH1D*)fBulkOuter->Get("hOVCth232M2");
  hOVCu238M2    = (TH1D*)fBulkOuter->Get("hOVCu238M2");

  hOVCco60M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hOVCco60M2Sum");
  hOVCk40M2Sum     = (TH1D*)fBulkOuterM2Sum->Get("hOVCk40M2Sum");
  hOVCth232M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("hOVCth232M2Sum");
  hOVCu238M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hOVCu238M2Sum");

/////// External Sources M1 and M2
  hExtPbbi210M1 = (TH1D*)fBulkOuter->Get("hExtPbbi210M1");
 
  hExtPbbi210M2 = (TH1D*)fBulkOuter->Get("hExtPbbi210M2");

  hExtPbbi210M2Sum = (TH1D*)fBulkOuterM2Sum->Get("hExtPbbi210M2Sum");

//////////// Surface PDFs
///// Crystal M1 and M2
  hTeO2Spb210M1_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Spb210M1_01");
  hTeO2Spo210M1_001   = (TH1D*)fSurfaceCrystal->Get("hTeO2Spo210M1_001");
  hTeO2Spo210M1_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Spo210M1_01");
  hTeO2Sth232M1_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sth232M1_01");
  hTeO2Su238M1_01     = (TH1D*)fSurfaceCrystal->Get("hTeO2Su238M1_01");
  hTeO2Sxpb210M1_001  = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M1_001");
  hTeO2Sxpb210M1_01   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M1_01");
  hTeO2Sxpb210M1_1    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M1_1");
  hTeO2Sxpb210M1_10   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M1_10");
  hTeO2Sxpo210M1_001  = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpo210M1_001");
  hTeO2Sxpo210M1_01   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpo210M1_01");
  hTeO2Sxpo210M1_1    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpo210M1_1");
  hTeO2Sxth232M1_001  = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M1_001");
  hTeO2Sxth232M1_01   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M1_01");
  hTeO2Sxth232M1_1    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M1_1");
  hTeO2Sxth232M1_10   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M1_10");
  hTeO2Sxu238M1_001   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M1_001");
  hTeO2Sxu238M1_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M1_01");
  hTeO2Sxu238M1_1     = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M1_1");
  hTeO2Sxu238M1_10    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M1_10");

  hTeO2Sxth232onlyM1_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232onlyM1_001");
  hTeO2Sxra228pb208M1_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra228pb208M1_001");
  hTeO2Sxu238th230M1_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238th230M1_001");
  hTeO2Sxth230onlyM1_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth230onlyM1_001");
  hTeO2Sxra226pb210M1_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra226pb210M1_001");
  hTeO2Sxpb210M1_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M1_0001");

  hTeO2Spb210M2_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Spb210M2_01");
  hTeO2Spo210M2_001   = (TH1D*)fSurfaceCrystal->Get("hTeO2Spo210M2_001");
  hTeO2Spo210M2_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Spo210M2_01");
  hTeO2Sth232M2_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sth232M2_01");
  hTeO2Su238M2_01     = (TH1D*)fSurfaceCrystal->Get("hTeO2Su238M2_01");
  hTeO2Sxpb210M2_001  = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M2_001");
  hTeO2Sxpb210M2_01   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M2_01");
  hTeO2Sxpb210M2_1    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M2_1");
  hTeO2Sxpb210M2_10   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M2_10");
  hTeO2Sxpo210M2_001  = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpo210M2_001");
  hTeO2Sxpo210M2_01   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpo210M2_01");
  hTeO2Sxpo210M2_1    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpo210M2_1");
  hTeO2Sxth232M2_001  = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M2_001");
  hTeO2Sxth232M2_01   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M2_01");
  hTeO2Sxth232M2_1    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M2_1");
  hTeO2Sxth232M2_10   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M2_10");
  hTeO2Sxu238M2_001   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M2_001");
  hTeO2Sxu238M2_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M2_01");
  hTeO2Sxu238M2_1     = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M2_1");
  hTeO2Sxu238M2_10    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M2_10");

  hTeO2Sxth232onlyM2_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232onlyM2_001");
  hTeO2Sxra228pb208M2_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra228pb208M2_001");
  hTeO2Sxu238th230M2_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238th230M2_001");
  hTeO2Sxth230onlyM2_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth230onlyM2_001");
  hTeO2Sxra226pb210M2_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra226pb210M2_001");
  hTeO2Sxpb210M2_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M2_0001");

  hTeO2Spb210M2Sum_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Spb210M2Sum_01");
  hTeO2Spo210M2Sum_001   = (TH1D*)fSurfaceCrystal->Get("hTeO2Spo210M2Sum_001");
  hTeO2Spo210M2Sum_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Spo210M2Sum_01");
  hTeO2Sth232M2Sum_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sth232M2Sum_01");
  hTeO2Su238M2Sum_01     = (TH1D*)fSurfaceCrystal->Get("hTeO2Su238M2Sum_01");
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

  hTeO2Sxth232onlyM2Sum_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232onlyM2Sum_001");
  hTeO2Sxra228pb208M2Sum_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra228pb208M2Sum_001");
  hTeO2Sxu238th230M2Sum_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238th230M2Sum_001");
  hTeO2Sxth230onlyM2Sum_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth230onlyM2Sum_001");
  hTeO2Sxra226pb210M2Sum_001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra226pb210M2Sum_001");
  hTeO2Sxpb210M2Sum_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M2Sum_0001");


//////// Frame M1 and M2
  hCuFrameSth232M1_1    = (TH1D*)fSurfaceOther->Get("hCuFrameSth232M1_1");
  hCuFrameSu238M1_1     = (TH1D*)fSurfaceOther->Get("hCuFrameSu238M1_1");
  hCuFrameSxpb210M1_001 = (TH1D*)fSurfaceOther->Get("hCuFrameSxpb210M1_001");
  hCuFrameSxpb210M1_01  = (TH1D*)fSurfaceOther->Get("hCuFrameSxpb210M1_01");
  hCuFrameSxpb210M1_1   = (TH1D*)fSurfaceOther->Get("hCuFrameSxpb210M1_1");
  hCuFrameSxpb210M1_10  = (TH1D*)fSurfaceOther->Get("hCuFrameSxpb210M1_10");
  hCuFrameSxth232M1_001 = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M1_001");
  hCuFrameSxth232M1_01  = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M1_01");
  hCuFrameSxth232M1_1   = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M1_1");
  hCuFrameSxth232M1_10  = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M1_10");
  hCuFrameSxu238M1_001  = (TH1D*)fSurfaceOther->Get("hCuFrameSxu238M1_001");
  hCuFrameSxu238M1_01   = (TH1D*)fSurfaceOther->Get("hCuFrameSxu238M1_01");
  hCuFrameSxu238M1_1    = (TH1D*)fSurfaceOther->Get("hCuFrameSxu238M1_1");
  hCuFrameSxu238M1_10   = (TH1D*)fSurfaceOther->Get("hCuFrameSxu238M1_10");

  hCuFrameSth232M2_1    = (TH1D*)fSurfaceOther->Get("hCuFrameSth232M2_1");
  hCuFrameSu238M2_1     = (TH1D*)fSurfaceOther->Get("hCuFrameSu238M2_1");
  hCuFrameSxpb210M2_001 = (TH1D*)fSurfaceOther->Get("hCuFrameSxpb210M2_001");
  hCuFrameSxpb210M2_01  = (TH1D*)fSurfaceOther->Get("hCuFrameSxpb210M2_01");
  hCuFrameSxpb210M2_1   = (TH1D*)fSurfaceOther->Get("hCuFrameSxpb210M2_1");
  hCuFrameSxpb210M2_10  = (TH1D*)fSurfaceOther->Get("hCuFrameSxpb210M2_10");
  hCuFrameSxth232M2_001 = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M2_001");
  hCuFrameSxth232M2_01  = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M2_01");
  hCuFrameSxth232M2_1   = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M2_1");
  hCuFrameSxth232M2_10  = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M2_10");
  hCuFrameSxu238M2_001  = (TH1D*)fSurfaceOther->Get("hCuFrameSxu238M2_001");
  hCuFrameSxu238M2_01   = (TH1D*)fSurfaceOther->Get("hCuFrameSxu238M2_01");
  hCuFrameSxu238M2_1    = (TH1D*)fSurfaceOther->Get("hCuFrameSxu238M2_1");
  hCuFrameSxu238M2_10   = (TH1D*)fSurfaceOther->Get("hCuFrameSxu238M2_10");

  hCuFrameSth232M2Sum_1    = (TH1D*)fSurfaceOther->Get("hCuFrameSth232M2Sum_1");
  hCuFrameSu238M2Sum_1     = (TH1D*)fSurfaceOther->Get("hCuFrameSu238M2Sum_1");
  hCuFrameSxpb210M2Sum_001 = (TH1D*)fSurfaceOther->Get("hCuFrameSxpb210M2Sum_001");
  hCuFrameSxpb210M2Sum_01  = (TH1D*)fSurfaceOther->Get("hCuFrameSxpb210M2Sum_01");
  hCuFrameSxpb210M2Sum_1   = (TH1D*)fSurfaceOther->Get("hCuFrameSxpb210M2Sum_1");
  hCuFrameSxpb210M2Sum_10  = (TH1D*)fSurfaceOther->Get("hCuFrameSxpb210M2Sum_10");
  hCuFrameSxth232M2Sum_001 = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M2Sum_001");
  hCuFrameSxth232M2Sum_01  = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M2Sum_01");
  hCuFrameSxth232M2Sum_1   = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M2Sum_1");
  hCuFrameSxth232M2Sum_10  = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M2Sum_10");
  hCuFrameSxu238M2Sum_001  = (TH1D*)fSurfaceOther->Get("hCuFrameSxu238M2Sum_001");
  hCuFrameSxu238M2Sum_01   = (TH1D*)fSurfaceOther->Get("hCuFrameSxu238M2Sum_01");
  hCuFrameSxu238M2Sum_1    = (TH1D*)fSurfaceOther->Get("hCuFrameSxu238M2Sum_1");
  hCuFrameSxu238M2Sum_10   = (TH1D*)fSurfaceOther->Get("hCuFrameSxu238M2Sum_10");

/////// CuBox M1 and M2
  hCuBoxSth232M1_1    = (TH1D*)fSurfaceOther->Get("hCuBoxSth232M1_1");
  hCuBoxSu238M1_1     = (TH1D*)fSurfaceOther->Get("hCuBoxSu238M1_1");
  hCuBoxSxpb210M1_001 = (TH1D*)fSurfaceOther->Get("hCuBoxSxpb210M1_001");
  hCuBoxSxpb210M1_01  = (TH1D*)fSurfaceOther->Get("hCuBoxSxpb210M1_01");
  hCuBoxSxpb210M1_1   = (TH1D*)fSurfaceOther->Get("hCuBoxSxpb210M1_1");
  hCuBoxSxpb210M1_10  = (TH1D*)fSurfaceOther->Get("hCuBoxSxpb210M1_10");
  hCuBoxSxth232M1_001 = (TH1D*)fSurfaceOther->Get("hCuBoxSxth232M1_001");
  hCuBoxSxth232M1_01  = (TH1D*)fSurfaceOther->Get("hCuBoxSxth232M1_01");
  hCuBoxSxth232M1_1   = (TH1D*)fSurfaceOther->Get("hCuBoxSxth232M1_1");
  hCuBoxSxth232M1_10  = (TH1D*)fSurfaceOther->Get("hCuBoxSxth232M1_10");
  hCuBoxSxu238M1_001  = (TH1D*)fSurfaceOther->Get("hCuBoxSxu238M1_001");
  hCuBoxSxu238M1_01   = (TH1D*)fSurfaceOther->Get("hCuBoxSxu238M1_01");
  hCuBoxSxu238M1_1    = (TH1D*)fSurfaceOther->Get("hCuBoxSxu238M1_1");
  hCuBoxSxu238M1_10   = (TH1D*)fSurfaceOther->Get("hCuBoxSxu238M1_10");

  hCuBoxSth232M2_1    = (TH1D*)fSurfaceOther->Get("hCuBoxSth232M2_1");
  hCuBoxSu238M2_1     = (TH1D*)fSurfaceOther->Get("hCuBoxSu238M2_1");
  hCuBoxSxpb210M2_001 = (TH1D*)fSurfaceOther->Get("hCuBoxSxpb210M2_001");
  hCuBoxSxpb210M2_01  = (TH1D*)fSurfaceOther->Get("hCuBoxSxpb210M2_01");
  hCuBoxSxpb210M2_1   = (TH1D*)fSurfaceOther->Get("hCuBoxSxpb210M2_1");
  hCuBoxSxpb210M2_10  = (TH1D*)fSurfaceOther->Get("hCuBoxSxpb210M2_10");
  hCuBoxSxth232M2_001 = (TH1D*)fSurfaceOther->Get("hCuBoxSxth232M2_001");
  hCuBoxSxth232M2_01  = (TH1D*)fSurfaceOther->Get("hCuBoxSxth232M2_01");
  hCuBoxSxth232M2_1   = (TH1D*)fSurfaceOther->Get("hCuBoxSxth232M2_1");
  hCuBoxSxth232M2_10  = (TH1D*)fSurfaceOther->Get("hCuBoxSxth232M2_10");
  hCuBoxSxu238M2_001  = (TH1D*)fSurfaceOther->Get("hCuBoxSxu238M2_001");
  hCuBoxSxu238M2_01   = (TH1D*)fSurfaceOther->Get("hCuBoxSxu238M2_01");
  hCuBoxSxu238M2_1    = (TH1D*)fSurfaceOther->Get("hCuBoxSxu238M2_1");
  hCuBoxSxu238M2_10   = (TH1D*)fSurfaceOther->Get("hCuBoxSxu238M2_10");

  hCuBoxSth232M2Sum_1    = (TH1D*)fSurfaceOther->Get("hCuBoxSth232M2Sum_1");
  hCuBoxSu238M2Sum_1     = (TH1D*)fSurfaceOther->Get("hCuBoxSu238M2Sum_1");
  hCuBoxSxpb210M2Sum_001 = (TH1D*)fSurfaceOther->Get("hCuBoxSxpb210M2Sum_001");
  hCuBoxSxpb210M2Sum_01  = (TH1D*)fSurfaceOther->Get("hCuBoxSxpb210M2Sum_01");
  hCuBoxSxpb210M2Sum_1   = (TH1D*)fSurfaceOther->Get("hCuBoxSxpb210M2Sum_1");
  hCuBoxSxpb210M2Sum_10  = (TH1D*)fSurfaceOther->Get("hCuBoxSxpb210M2Sum_10");
  hCuBoxSxth232M2Sum_001 = (TH1D*)fSurfaceOther->Get("hCuBoxSxth232M2Sum_001");
  hCuBoxSxth232M2Sum_01  = (TH1D*)fSurfaceOther->Get("hCuBoxSxth232M2Sum_01");
  hCuBoxSxth232M2Sum_1   = (TH1D*)fSurfaceOther->Get("hCuBoxSxth232M2Sum_1");
  hCuBoxSxth232M2Sum_10  = (TH1D*)fSurfaceOther->Get("hCuBoxSxth232M2Sum_10");
  hCuBoxSxu238M2Sum_001  = (TH1D*)fSurfaceOther->Get("hCuBoxSxu238M2Sum_001");
  hCuBoxSxu238M2Sum_01   = (TH1D*)fSurfaceOther->Get("hCuBoxSxu238M2Sum_01");
  hCuBoxSxu238M2Sum_1    = (TH1D*)fSurfaceOther->Get("hCuBoxSxu238M2Sum_1");
  hCuBoxSxu238M2Sum_10   = (TH1D*)fSurfaceOther->Get("hCuBoxSxu238M2Sum_10");

/////// CuBox + CuFrame

  hCuBox_CuFrameth232M1_10 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M1_10");
  hCuBox_CuFrameu238M1_10 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameu238M1_10");
  hCuBox_CuFramepb210M1_10 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M1_10");
  hCuBox_CuFramepb210M1_1 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M1_1");
  hCuBox_CuFramepb210M1_01 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M1_01");
  hCuBox_CuFramepb210M1_001 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M1_001");

  hCuBox_CuFrameth232M2_10 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2_10");
  hCuBox_CuFrameu238M2_10 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameu238M2_10");
  hCuBox_CuFramepb210M2_10 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M2_10");
  hCuBox_CuFramepb210M2_1 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M2_1");
  hCuBox_CuFramepb210M2_01 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M2_01");
  hCuBox_CuFramepb210M2_001 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M2_001");

  hCuBox_CuFrameth232M2Sum_10 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2Sum_10");
  hCuBox_CuFrameu238M2Sum_10 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameu238M2Sum_10");
  hCuBox_CuFramepb210M2Sum_10 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M2Sum_10");
  hCuBox_CuFramepb210M2Sum_1 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M2Sum_1");
  hCuBox_CuFramepb210M2Sum_01 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M2Sum_01");
  hCuBox_CuFramepb210M2Sum_001 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFramepb210M2Sum_001"); 

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

  hnewTeO2Spb210M1_01 = hTeO2Spb210M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Spb210M1_01", dAdaptiveArrayM1);
  hnewTeO2Spo210M1_001 = hTeO2Spo210M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Spo210M1_001", dAdaptiveArrayM1);
  hnewTeO2Spo210M1_01 = hTeO2Spo210M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Spo210M1_01", dAdaptiveArrayM1);
  hnewTeO2Sth232M1_01 = hTeO2Sth232M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sth232M1_01", dAdaptiveArrayM1);
  hnewTeO2Su238M1_01 = hTeO2Su238M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Su238M1_01", dAdaptiveArrayM1);
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

  hnewTeO2th232onlyM1 = hTeO2th232onlyM1->Rebin(dAdaptiveBinsM1, "hnewTeO2th232onlyM1", dAdaptiveArrayM1);
  hnewTeO2ra228pb208M1 = hTeO2ra228pb208M1->Rebin(dAdaptiveBinsM1, "hnewTeO2ra228pb208M1", dAdaptiveArrayM1);
  hnewTeO2th230onlyM1 = hTeO2th230onlyM1->Rebin(dAdaptiveBinsM1, "hnewTeO2th230onlyM1", dAdaptiveArrayM1);

  hnewTeO2Sxth232onlyM1_001 = hTeO2Sxth232onlyM1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232onlyM1_001", dAdaptiveArrayM1);
  hnewTeO2Sxra228pb208M1_001 = hTeO2Sxra228pb208M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxra228pb208M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxu238th230M1_001 = hTeO2Sxu238th230M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238th230M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxth230onlyM1_001 = hTeO2Sxth230onlyM1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth230onlyM1_001", dAdaptiveArrayM1);
  hnewTeO2Sxra226pb210M1_001 = hTeO2Sxra226pb210M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxra226pb210M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxpb210M1_0001 = hTeO2Sxpb210M1_0001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpb210M1_0001", dAdaptiveArrayM1);


  hnewTeO20nuM2 = hTeO20nuM2->Rebin(dAdaptiveBinsM2, "hnewTeO20nuM2", dAdaptiveArrayM2);
  hnewTeO22nuM2 = hTeO22nuM2->Rebin(dAdaptiveBinsM2, "hnewTeO22nuM2", dAdaptiveArrayM2);
  hnewTeO2co60M2 = hTeO2co60M2->Rebin(dAdaptiveBinsM2, "hnewTeO2co60M2", dAdaptiveArrayM2);
  hnewTeO2k40M2 = hTeO2k40M2->Rebin(dAdaptiveBinsM2, "hnewTeO2k40M2", dAdaptiveArrayM2);
  hnewTeO2pb210M2 = hTeO2pb210M2->Rebin(dAdaptiveBinsM2, "hnewTeO2pb210M2", dAdaptiveArrayM2);
  hnewTeO2po210M2 = hTeO2po210M2->Rebin(dAdaptiveBinsM2, "hnewTeO2po210M2", dAdaptiveArrayM2);
  hnewTeO2te125M2 = hTeO2te125M2->Rebin(dAdaptiveBinsM2, "hnewTeO2te125M2", dAdaptiveArrayM2);
  hnewTeO2th232M2 = hTeO2th232M2->Rebin(dAdaptiveBinsM2, "hnewTeO2th232M2", dAdaptiveArrayM2);
  hnewTeO2u238M2 = hTeO2u238M2->Rebin(dAdaptiveBinsM2, "hnewTeO2u238M2", dAdaptiveArrayM2);

  hnewTeO2Spb210M2_01 = hTeO2Spb210M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Spb210M2_01", dAdaptiveArrayM2);
  hnewTeO2Spo210M2_001 = hTeO2Spo210M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Spo210M2_001", dAdaptiveArrayM2);
  hnewTeO2Spo210M2_01 = hTeO2Spo210M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Spo210M2_01", dAdaptiveArrayM2);
  hnewTeO2Sth232M2_01 = hTeO2Sth232M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sth232M2_01", dAdaptiveArrayM2);
  hnewTeO2Su238M2_01 = hTeO2Su238M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Su238M2_01", dAdaptiveArrayM2);
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

  hnewTeO2th232onlyM2 = hTeO2th232onlyM2->Rebin(dAdaptiveBinsM2, "hnewTeO2th232onlyM2", dAdaptiveArrayM2);
  hnewTeO2ra228pb208M2 = hTeO2ra228pb208M2->Rebin(dAdaptiveBinsM2, "hnewTeO2ra228pb208M2", dAdaptiveArrayM2);
  hnewTeO2th230onlyM2 = hTeO2th230onlyM2->Rebin(dAdaptiveBinsM2, "hnewTeO2th230onlyM2", dAdaptiveArrayM2);

  hnewTeO2Sxth232onlyM2_001 = hTeO2Sxth232onlyM2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232onlyM2_001", dAdaptiveArrayM2);
  hnewTeO2Sxra228pb208M2_001 = hTeO2Sxra228pb208M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxra228pb208M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxu238th230M2_001 = hTeO2Sxu238th230M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238th230M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxth230onlyM2_001 = hTeO2Sxth230onlyM2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth230onlyM2_001", dAdaptiveArrayM2);
  hnewTeO2Sxra226pb210M2_001 = hTeO2Sxra226pb210M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxra226pb210M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxpb210M2_0001 = hTeO2Sxpb210M2_0001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpb210M2_0001", dAdaptiveArrayM2);


  hnewTeO20nuM2Sum = hTeO20nuM2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO20nuM2Sum", dAdaptiveArrayM2Sum);
  hnewTeO22nuM2Sum = hTeO22nuM2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO22nuM2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2co60M2Sum = hTeO2co60M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2co60M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2k40M2Sum = hTeO2k40M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2k40M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2pb210M2Sum = hTeO2pb210M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2pb210M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2po210M2Sum = hTeO2po210M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2po210M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2te125M2Sum = hTeO2te125M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2te125M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2th232M2Sum = hTeO2th232M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2th232M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2u238M2Sum = hTeO2u238M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2u238M2Sum", dAdaptiveArrayM2Sum);

  hnewTeO2Spb210M2Sum_01 = hTeO2Spb210M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Spb210M2Sum_01", dAdaptiveArrayM2Sum);
  hnewTeO2Spo210M2Sum_001 = hTeO2Spo210M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Spo210M2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Spo210M2Sum_01 = hTeO2Spo210M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Spo210M2Sum_01", dAdaptiveArrayM2Sum);
  hnewTeO2Sth232M2Sum_01 = hTeO2Sth232M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sth232M2Sum_01", dAdaptiveArrayM2Sum);
  hnewTeO2Su238M2Sum_01 = hTeO2Su238M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Su238M2Sum_01", dAdaptiveArrayM2Sum);
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

  hnewTeO2th232onlyM2Sum = hTeO2th232onlyM2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2th232onlyM2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2ra228pb208M2Sum = hTeO2ra228pb208M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2ra228pb208M2Sum", dAdaptiveArrayM2Sum);
  hnewTeO2th230onlyM2Sum = hTeO2th230onlyM2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2th230onlyM2Sum", dAdaptiveArrayM2Sum);

  hnewTeO2Sxth232onlyM2Sum_001 = hTeO2Sxth232onlyM2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxth232onlyM2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxra228pb208M2Sum_001 = hTeO2Sxra228pb208M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxra228pb208M2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxu238th230M2Sum_001 = hTeO2Sxu238th230M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxu238th230M2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxth230onlyM2Sum_001 = hTeO2Sxth230onlyM2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxth230onlyM2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxra226pb210M2Sum_001 = hTeO2Sxra226pb210M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxra226pb210M2Sum_001", dAdaptiveArrayM2Sum);
  hnewTeO2Sxpb210M2Sum_0001 = hTeO2Sxpb210M2Sum_0001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sxpb210M2Sum_0001", dAdaptiveArrayM2Sum);


////////// Frame M1 and M2
  hnewCuFrameco58M1 = hCuFrameco58M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameco58M1", dAdaptiveArrayM1);
  hnewCuFrameco60M1 = hCuFrameco60M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameco60M1", dAdaptiveArrayM1);
  hnewCuFramecs137M1 = hCuFramecs137M1->Rebin(dAdaptiveBinsM1, "hnewCuFramecs137M1", dAdaptiveArrayM1);
  hnewCuFramek40M1 = hCuFramek40M1->Rebin(dAdaptiveBinsM1, "hnewCuFramek40M1", dAdaptiveArrayM1);
  hnewCuFramemn54M1 = hCuFramemn54M1->Rebin(dAdaptiveBinsM1, "hnewCuFramemn54M1", dAdaptiveArrayM1);
  hnewCuFramepb210M1 = hCuFramepb210M1->Rebin(dAdaptiveBinsM1, "hnewCuFramepb210M1", dAdaptiveArrayM1);
  hnewCuFrameth232M1 = hCuFrameth232M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameth232M1", dAdaptiveArrayM1);
  hnewCuFrameu238M1 = hCuFrameu238M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameu238M1", dAdaptiveArrayM1);

  hnewCuFrameSth232M1_1  = hCuFrameSth232M1_1->Rebin(dAdaptiveBinsM1, "hnewCuFrameSth232M1_1", dAdaptiveArrayM1);
  hnewCuFrameSu238M1_1 = hCuFrameSu238M1_1->Rebin(dAdaptiveBinsM1, "hnewCuFrameSu238M1_1", dAdaptiveArrayM1);
  hnewCuFrameSxpb210M1_001 = hCuFrameSxpb210M1_001->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxpb210M1_001", dAdaptiveArrayM1);
  hnewCuFrameSxpb210M1_01 = hCuFrameSxpb210M1_01->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxpb210M1_01", dAdaptiveArrayM1);
  hnewCuFrameSxpb210M1_1 = hCuFrameSxpb210M1_1->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxpb210M1_1", dAdaptiveArrayM1);
  hnewCuFrameSxpb210M1_10 = hCuFrameSxpb210M1_10->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxpb210M1_10", dAdaptiveArrayM1);
  hnewCuFrameSxth232M1_001 = hCuFrameSxth232M1_001->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxth232M1_001", dAdaptiveArrayM1);
  hnewCuFrameSxth232M1_01 = hCuFrameSxth232M1_01->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxth232M1_01", dAdaptiveArrayM1);
  hnewCuFrameSxth232M1_1 = hCuFrameSxth232M1_1->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxth232M1_1", dAdaptiveArrayM1);
  hnewCuFrameSxth232M1_10 = hCuFrameSxth232M1_10->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxth232M1_10", dAdaptiveArrayM1);
  hnewCuFrameSxu238M1_001 = hCuFrameSxu238M1_001->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxu238M1_001", dAdaptiveArrayM1);
  hnewCuFrameSxu238M1_01 = hCuFrameSxu238M1_01->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxu238M1_01", dAdaptiveArrayM1);
  hnewCuFrameSxu238M1_1 = hCuFrameSxu238M1_1->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxu238M1_1", dAdaptiveArrayM1);
  hnewCuFrameSxu238M1_10 = hCuFrameSxu238M1_10->Rebin(dAdaptiveBinsM1, "hnewCuFrameSxu238M1_10", dAdaptiveArrayM1);

  hnewCuFrameco58M2 = hCuFrameco58M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameco58M2", dAdaptiveArrayM2);
  hnewCuFrameco60M2 = hCuFrameco60M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameco60M2", dAdaptiveArrayM2);
  hnewCuFramecs137M2 = hCuFramecs137M2->Rebin(dAdaptiveBinsM2, "hnewCuFramecs137M2", dAdaptiveArrayM2);
  hnewCuFramek40M2 = hCuFramek40M2->Rebin(dAdaptiveBinsM2, "hnewCuFramek40M2", dAdaptiveArrayM2);
  hnewCuFramemn54M2 = hCuFramemn54M2->Rebin(dAdaptiveBinsM2, "hnewCuFramemn54M2", dAdaptiveArrayM2);
  hnewCuFramepb210M2 = hCuFramepb210M2->Rebin(dAdaptiveBinsM2, "hnewCuFramepb210M2", dAdaptiveArrayM2);
  hnewCuFrameth232M2 = hCuFrameth232M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameth232M2", dAdaptiveArrayM2);
  hnewCuFrameu238M2 = hCuFrameu238M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameu238M2", dAdaptiveArrayM2);

  hnewCuFrameSth232M2_1 = hCuFrameSth232M2_1->Rebin(dAdaptiveBinsM2, "hnewCuFrameSth232M2_1", dAdaptiveArrayM2);
  hnewCuFrameSu238M2_1 = hCuFrameSu238M2_1->Rebin(dAdaptiveBinsM2, "hnewCuFrameSu238M2_1", dAdaptiveArrayM2);
  hnewCuFrameSxpb210M2_001 = hCuFrameSxpb210M2_001->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxpb210M2_001", dAdaptiveArrayM2);
  hnewCuFrameSxpb210M2_01 = hCuFrameSxpb210M2_01->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxpb210M2_01", dAdaptiveArrayM2);
  hnewCuFrameSxpb210M2_1 = hCuFrameSxpb210M2_1->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxpb210M2_1", dAdaptiveArrayM2);
  hnewCuFrameSxpb210M2_10 = hCuFrameSxpb210M2_10->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxpb210M2_10", dAdaptiveArrayM2);
  hnewCuFrameSxth232M2_001 = hCuFrameSxth232M2_001->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxth232M2_001", dAdaptiveArrayM2);
  hnewCuFrameSxth232M2_01 = hCuFrameSxth232M2_01->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxth232M2_01", dAdaptiveArrayM2);
  hnewCuFrameSxth232M2_1 = hCuFrameSxth232M2_1->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxth232M2_1", dAdaptiveArrayM2);
  hnewCuFrameSxth232M2_10 = hCuFrameSxth232M2_10->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxth232M2_10", dAdaptiveArrayM2);
  hnewCuFrameSxu238M2_001 = hCuFrameSxu238M2_001->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxu238M2_001", dAdaptiveArrayM2);
  hnewCuFrameSxu238M2_01 = hCuFrameSxu238M2_01->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxu238M2_01", dAdaptiveArrayM2);
  hnewCuFrameSxu238M2_1 = hCuFrameSxu238M2_1->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxu238M2_1", dAdaptiveArrayM2);
  hnewCuFrameSxu238M2_10 = hCuFrameSxu238M2_10->Rebin(dAdaptiveBinsM2, "hnewCuFrameSxu238M2_10", dAdaptiveArrayM2);

  hnewCuFrameco58M2Sum = hCuFrameco58M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameco58M2Sum", dAdaptiveArrayM2Sum);
  hnewCuFrameco60M2Sum = hCuFrameco60M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameco60M2Sum", dAdaptiveArrayM2Sum);
  hnewCuFramecs137M2Sum = hCuFramecs137M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuFramecs137M2Sum", dAdaptiveArrayM2Sum);
  hnewCuFramek40M2Sum = hCuFramek40M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuFramek40M2Sum", dAdaptiveArrayM2Sum);
  hnewCuFramemn54M2Sum = hCuFramemn54M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuFramemn54M2Sum", dAdaptiveArrayM2Sum);
  hnewCuFramepb210M2Sum = hCuFramepb210M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuFramepb210M2Sum", dAdaptiveArrayM2Sum);
  hnewCuFrameth232M2Sum = hCuFrameth232M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameth232M2Sum", dAdaptiveArrayM2Sum);
  hnewCuFrameu238M2Sum = hCuFrameu238M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameu238M2Sum", dAdaptiveArrayM2Sum);

  hnewCuFrameSth232M2Sum_1 = hCuFrameSth232M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSth232M2Sum_1", dAdaptiveArrayM2Sum);
  hnewCuFrameSu238M2Sum_1 = hCuFrameSu238M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSu238M2Sum_1", dAdaptiveArrayM2Sum);
  hnewCuFrameSxpb210M2Sum_001 = hCuFrameSxpb210M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSxpb210M2Sum_001", dAdaptiveArrayM2Sum);
  hnewCuFrameSxpb210M2Sum_01 = hCuFrameSxpb210M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSxpb210M2Sum_01", dAdaptiveArrayM2Sum);
  hnewCuFrameSxpb210M2Sum_1 = hCuFrameSxpb210M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSxpb210M2Sum_1", dAdaptiveArrayM2Sum);
  hnewCuFrameSxpb210M2Sum_10 = hCuFrameSxpb210M2Sum_10->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSxpb210M2Sum_10", dAdaptiveArrayM2Sum);
  hnewCuFrameSxth232M2Sum_001 = hCuFrameSxth232M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSxth232M2Sum_001", dAdaptiveArrayM2Sum);
  hnewCuFrameSxth232M2Sum_01 = hCuFrameSxth232M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSxth232M2Sum_01", dAdaptiveArrayM2Sum);
  hnewCuFrameSxth232M2Sum_1 = hCuFrameSxth232M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSxth232M2Sum_1", dAdaptiveArrayM2Sum);
  hnewCuFrameSxth232M2Sum_10 = hCuFrameSxth232M2Sum_10->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSxth232M2Sum_10", dAdaptiveArrayM2Sum);
  hnewCuFrameSxu238M2Sum_001 = hCuFrameSxu238M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSxu238M2Sum_001", dAdaptiveArrayM2Sum);
  hnewCuFrameSxu238M2Sum_01 = hCuFrameSxu238M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSxu238M2Sum_01", dAdaptiveArrayM2Sum);
  hnewCuFrameSxu238M2Sum_1 = hCuFrameSxu238M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSxu238M2Sum_1", dAdaptiveArrayM2Sum);
  hnewCuFrameSxu238M2Sum_10 = hCuFrameSxu238M2Sum_10->Rebin(dAdaptiveBinsM2Sum, "hnewCuFrameSxu238M2Sum_10", dAdaptiveArrayM2Sum);

///////// CuBox (TShield) M1 and M2
  hnewCuBoxco58M1 = hCuBoxco58M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxco58M1", dAdaptiveArrayM1);
  hnewCuBoxco60M1 = hCuBoxco60M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxco60M1", dAdaptiveArrayM1);
  hnewCuBoxcs137M1 = hCuBoxcs137M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxcs137M1", dAdaptiveArrayM1);
  hnewCuBoxk40M1 = hCuBoxk40M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxk40M1", dAdaptiveArrayM1);
  hnewCuBoxmn54M1 = hCuBoxmn54M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxmn54M1", dAdaptiveArrayM1);
  hnewCuBoxpb210M1 = hCuBoxpb210M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxpb210M1", dAdaptiveArrayM1);
  hnewCuBoxth232M1 = hCuBoxth232M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxth232M1", dAdaptiveArrayM1);
  hnewCuBoxu238M1 = hCuBoxu238M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxu238M1", dAdaptiveArrayM1);

  hnewCuBoxSth232M1_1 = hCuBoxSth232M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBoxSth232M1_1", dAdaptiveArrayM1);
  hnewCuBoxSu238M1_1 = hCuBoxSu238M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBoxSu238M1_1", dAdaptiveArrayM1);
  hnewCuBoxSxpb210M1_001 = hCuBoxSxpb210M1_001->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxpb210M1_001", dAdaptiveArrayM1);
  hnewCuBoxSxpb210M1_01 = hCuBoxSxpb210M1_01->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxpb210M1_01", dAdaptiveArrayM1);
  hnewCuBoxSxpb210M1_1 = hCuBoxSxpb210M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxpb210M1_1", dAdaptiveArrayM1);
  hnewCuBoxSxpb210M1_10 = hCuBoxSxpb210M1_10->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxpb210M1_10", dAdaptiveArrayM1);
  hnewCuBoxSxth232M1_001 = hCuBoxSxth232M1_001->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxth232M1_001", dAdaptiveArrayM1);
  hnewCuBoxSxth232M1_01 = hCuBoxSxth232M1_01->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxth232M1_01", dAdaptiveArrayM1);
  hnewCuBoxSxth232M1_1 = hCuBoxSxth232M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxth232M1_1", dAdaptiveArrayM1);
  hnewCuBoxSxth232M1_10 = hCuBoxSxth232M1_10->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxth232M1_10", dAdaptiveArrayM1);
  hnewCuBoxSxu238M1_001 = hCuBoxSxu238M1_001->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxu238M1_001", dAdaptiveArrayM1);
  hnewCuBoxSxu238M1_01 = hCuBoxSxu238M1_01->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxu238M1_01", dAdaptiveArrayM1);
  hnewCuBoxSxu238M1_1 = hCuBoxSxu238M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxu238M1_1", dAdaptiveArrayM1);
  hnewCuBoxSxu238M1_10 = hCuBoxSxu238M1_10->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxu238M1_10", dAdaptiveArrayM1);

  hnewCuBoxco58M2 = hCuBoxco58M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxco58M2", dAdaptiveArrayM2);
  hnewCuBoxco60M2 = hCuBoxco60M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxco60M2", dAdaptiveArrayM2);
  hnewCuBoxcs137M2 = hCuBoxcs137M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxcs137M2", dAdaptiveArrayM2);
  hnewCuBoxk40M2 = hCuBoxk40M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxk40M2", dAdaptiveArrayM2);
  hnewCuBoxmn54M2 = hCuBoxmn54M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxmn54M2", dAdaptiveArrayM2);
  hnewCuBoxpb210M2 = hCuBoxpb210M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxpb210M2", dAdaptiveArrayM2);
  hnewCuBoxth232M2 = hCuBoxth232M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxth232M2", dAdaptiveArrayM2);
  hnewCuBoxu238M2 = hCuBoxu238M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxu238M2", dAdaptiveArrayM2);

  hnewCuBoxSth232M2_1 = hCuBoxSth232M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBoxSth232M2_1", dAdaptiveArrayM2);
  hnewCuBoxSu238M2_1 = hCuBoxSu238M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBoxSu238M2_1", dAdaptiveArrayM2);
  hnewCuBoxSxpb210M2_001 = hCuBoxSxpb210M2_001->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxpb210M2_001", dAdaptiveArrayM2);
  hnewCuBoxSxpb210M2_01 = hCuBoxSxpb210M2_01->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxpb210M2_01", dAdaptiveArrayM2);
  hnewCuBoxSxpb210M2_1 = hCuBoxSxpb210M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxpb210M2_1", dAdaptiveArrayM2);
  hnewCuBoxSxpb210M2_10 = hCuBoxSxpb210M2_10->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxpb210M2_10", dAdaptiveArrayM2);
  hnewCuBoxSxth232M2_001 = hCuBoxSxth232M2_001->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxth232M2_001", dAdaptiveArrayM2);
  hnewCuBoxSxth232M2_01 = hCuBoxSxth232M2_01->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxth232M2_01", dAdaptiveArrayM2);
  hnewCuBoxSxth232M2_1 = hCuBoxSxth232M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxth232M2_1", dAdaptiveArrayM2);
  hnewCuBoxSxth232M2_10 = hCuBoxSxth232M2_10->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxth232M2_10", dAdaptiveArrayM2);
  hnewCuBoxSxu238M2_001 = hCuBoxSxu238M2_001->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxu238M2_001", dAdaptiveArrayM2);
  hnewCuBoxSxu238M2_01 = hCuBoxSxu238M2_01->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxu238M2_01", dAdaptiveArrayM2);
  hnewCuBoxSxu238M2_1 = hCuBoxSxu238M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxu238M2_1", dAdaptiveArrayM2);
  hnewCuBoxSxu238M2_10 = hCuBoxSxu238M2_10->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxu238M2_10", dAdaptiveArrayM2);


  hnewCuBoxco58M2Sum = hCuBoxco58M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxco58M2Sum", dAdaptiveArrayM2Sum);
  hnewCuBoxco60M2Sum = hCuBoxco60M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxco60M2Sum", dAdaptiveArrayM2Sum);
  hnewCuBoxcs137M2Sum = hCuBoxcs137M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxcs137M2Sum", dAdaptiveArrayM2Sum);
  hnewCuBoxk40M2Sum = hCuBoxk40M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxk40M2Sum", dAdaptiveArrayM2Sum);
  hnewCuBoxmn54M2Sum = hCuBoxmn54M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxmn54M2Sum", dAdaptiveArrayM2Sum);
  hnewCuBoxpb210M2Sum = hCuBoxpb210M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxpb210M2Sum", dAdaptiveArrayM2Sum);
  hnewCuBoxth232M2Sum = hCuBoxth232M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxth232M2Sum", dAdaptiveArrayM2Sum);
  hnewCuBoxu238M2Sum = hCuBoxu238M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxu238M2Sum", dAdaptiveArrayM2Sum);

  hnewCuBoxSth232M2Sum_1 = hCuBoxSth232M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSth232M2Sum_1", dAdaptiveArrayM2Sum);
  hnewCuBoxSu238M2Sum_1 = hCuBoxSu238M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSu238M2Sum_1", dAdaptiveArrayM2Sum);
  hnewCuBoxSxpb210M2Sum_001 = hCuBoxSxpb210M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSxpb210M2Sum_001", dAdaptiveArrayM2Sum);
  hnewCuBoxSxpb210M2Sum_01 = hCuBoxSxpb210M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSxpb210M2Sum_01", dAdaptiveArrayM2Sum);
  hnewCuBoxSxpb210M2Sum_1 = hCuBoxSxpb210M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSxpb210M2Sum_1", dAdaptiveArrayM2Sum);
  hnewCuBoxSxpb210M2Sum_10 = hCuBoxSxpb210M2Sum_10->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSxpb210M2Sum_10", dAdaptiveArrayM2Sum);
  hnewCuBoxSxth232M2Sum_001 = hCuBoxSxth232M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSxth232M2Sum_001", dAdaptiveArrayM2Sum);
  hnewCuBoxSxth232M2Sum_01 = hCuBoxSxth232M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSxth232M2Sum_01", dAdaptiveArrayM2Sum);
  hnewCuBoxSxth232M2Sum_1 = hCuBoxSxth232M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSxth232M2Sum_1", dAdaptiveArrayM2Sum);
  hnewCuBoxSxth232M2Sum_10 = hCuBoxSxth232M2Sum_10->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSxth232M2Sum_10", dAdaptiveArrayM2Sum);
  hnewCuBoxSxu238M2Sum_001 = hCuBoxSxu238M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSxu238M2Sum_001", dAdaptiveArrayM2Sum);
  hnewCuBoxSxu238M2Sum_01 = hCuBoxSxu238M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSxu238M2Sum_01", dAdaptiveArrayM2Sum);
  hnewCuBoxSxu238M2Sum_1 = hCuBoxSxu238M2Sum_1->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSxu238M2Sum_1", dAdaptiveArrayM2Sum);
  hnewCuBoxSxu238M2Sum_10 = hCuBoxSxu238M2Sum_10->Rebin(dAdaptiveBinsM2Sum, "hnewCuBoxSxu238M2Sum_10", dAdaptiveArrayM2Sum);


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

/*
////////// 50mK M1 and M2
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

  h50mKco58M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnew50mKco58M2Sum", dAdaptiveArrayM2Sum);
  h50mKco60M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnew50mKco60M2Sum", dAdaptiveArrayM2Sum);
  h50mKcs137M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnew50mKcs137M2Sum", dAdaptiveArrayM2Sum);
  h50mKk40M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnew50mKk40M2Sum", dAdaptiveArrayM2Sum);
  h50mKmn54M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnew50mKmn54M2Sum", dAdaptiveArrayM2Sum);
  h50mKpb210M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnew50mKpb210M2Sum", dAdaptiveArrayM2Sum);
  h50mKth232M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnew50mKth232M2Sum", dAdaptiveArrayM2Sum);
  h50mKu238M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnew50mKu238M2Sum", dAdaptiveArrayM2Sum);
*/
//////// 600mK
  hnew600mKco60M1 = h600mKco60M1->Rebin(dAdaptiveBinsM1, "hnew600mKco60M1", dAdaptiveArrayM1);
  hnew600mKk40M1 = h600mKk40M1->Rebin(dAdaptiveBinsM1, "hnew600mKk40M1", dAdaptiveArrayM1);
  hnew600mKth232M1 = h600mKth232M1->Rebin(dAdaptiveBinsM1, "hnew600mKth232M1", dAdaptiveArrayM1);
  hnew600mKu238M1 = h600mKu238M1->Rebin(dAdaptiveBinsM1, "hnew600mKu238M1", dAdaptiveArrayM1);

  hnew600mKco60M2 = h600mKco60M2->Rebin(dAdaptiveBinsM2, "hnew600mKco60M2", dAdaptiveArrayM2);
  hnew600mKk40M2 = h600mKk40M2->Rebin(dAdaptiveBinsM2, "hnew600mKk40M2", dAdaptiveArrayM2);
  hnew600mKth232M2 = h600mKth232M2->Rebin(dAdaptiveBinsM2, "hnew600mKth232M2", dAdaptiveArrayM2);
  hnew600mKu238M2 = h600mKu238M2->Rebin(dAdaptiveBinsM2, "hnew600mKu238M2", dAdaptiveArrayM2);

  hnew600mKco60M2Sum = h600mKco60M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnew600mKco60M2Sum", dAdaptiveArrayM2Sum);
  hnew600mKk40M2Sum = h600mKk40M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnew600mKk40M2Sum", dAdaptiveArrayM2Sum);
  hnew600mKth232M2Sum = h600mKth232M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnew600mKth232M2Sum", dAdaptiveArrayM2Sum);
  hnew600mKu238M2Sum = h600mKu238M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnew600mKu238M2Sum", dAdaptiveArrayM2Sum);


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

/*
////////// Main Bath M1 and M2
  hMBco60M1->Rebin(dAdaptiveBinsM1, "hnewMBco60M1", dAdaptiveArrayM1);
  hMBk40M1->Rebin(dAdaptiveBinsM1, "hnewMBk40M1", dAdaptiveArrayM1);
  hMBth232M1->Rebin(dAdaptiveBinsM1, "hnewMBth232M1", dAdaptiveArrayM1);
  hMBu238M1->Rebin(dAdaptiveBinsM1, "hnewMBu238M1", dAdaptiveArrayM1);

  hMBco60M2->Rebin(dAdaptiveBinsM2, "hnewMBco60M2", dAdaptiveArrayM2);
  hMBk40M2->Rebin(dAdaptiveBinsM2, "hnewMBk40M2", dAdaptiveArrayM2);
  hMBth232M2->Rebin(dAdaptiveBinsM2, "hnewMBth232M2", dAdaptiveArrayM2);
  hMBu238M2->Rebin(dAdaptiveBinsM2, "hnewMBu238M2", dAdaptiveArrayM2);

  hMBco60M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewMBco60M2Sum", dAdaptiveArrayM2Sum);
  hMBk40M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewMBk40M2Sum", dAdaptiveArrayM2Sum);
  hMBth232M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewMBth232M2Sum", dAdaptiveArrayM2Sum);
  hMBu238M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewMBu238M2Sum", dAdaptiveArrayM2Sum);

///////// Super Insulation M1 and M2
  hSIk40M1->Rebin(dAdaptiveBinsM1, "hnewSIk40M1", dAdaptiveArrayM1);
  hSIth232M1->Rebin(dAdaptiveBinsM1, "hnewSIth232M1", dAdaptiveArrayM1);
  hSIu238M1->Rebin(dAdaptiveBinsM1, "hnewSIu238M1", dAdaptiveArrayM1);

  hSIk40M2->Rebin(dAdaptiveBinsM2, "hnewSIk40M2", dAdaptiveArrayM2);
  hSIth232M2->Rebin(dAdaptiveBinsM2, "hnewSIth232M2", dAdaptiveArrayM2);
  hSIu238M2->Rebin(dAdaptiveBinsM2, "hnewSIu238M2", dAdaptiveArrayM2);

  hSIk40M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewSIk40M2Sum", dAdaptiveArrayM2Sum);
  hSIth232M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewSIth232M2Sum", dAdaptiveArrayM2Sum);
  hSIu238M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewSIu238M2Sum", dAdaptiveArrayM2Sum);

///////// IVC M1 and M2
  hIVCco60M1->Rebin(dAdaptiveBinsM1, "hnewIVCco60M1", dAdaptiveArrayM1);
  hIVCk40M1->Rebin(dAdaptiveBinsM1, "hnewIVCk40M1", dAdaptiveArrayM1);
  hIVCth232M1->Rebin(dAdaptiveBinsM1, "hnewIVCth232M1", dAdaptiveArrayM1);
  hIVCu238M1->Rebin(dAdaptiveBinsM1, "hnewIVCu238M1", dAdaptiveArrayM1);

  hIVCco60M2->Rebin(dAdaptiveBinsM2, "hnewIVCco60M2", dAdaptiveArrayM2);
  hIVCk40M2->Rebin(dAdaptiveBinsM2, "hnewIVCk40M2", dAdaptiveArrayM2);
  hIVCth232M2->Rebin(dAdaptiveBinsM2, "hnewIVCth232M2", dAdaptiveArrayM2);
  hIVCu238M2->Rebin(dAdaptiveBinsM2, "hnewIVCu238M2", dAdaptiveArrayM2);

  hIVCco60M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewIVCco60M2Sum", dAdaptiveArrayM2Sum);
  hIVCk40M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewIVCk40M2Sum", dAdaptiveArrayM2Sum);
  hIVCth232M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewIVCth232M2Sum", dAdaptiveArrayM2Sum);
  hIVCu238M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewIVCu238M2Sum", dAdaptiveArrayM2Sum);
*/

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

    hAdapTeO2th232onlyM1->SetBinContent(i, hnewTeO2th232onlyM1->GetBinContent(i)/hnewTeO2th232onlyM1->GetBinWidth(i));
    hAdapTeO2ra228pb208M1->SetBinContent(i, hnewTeO2ra228pb208M1->GetBinContent(i)/hnewTeO2ra228pb208M1->GetBinWidth(i));
    hAdapTeO2th230onlyM1->SetBinContent(i, hnewTeO2th230onlyM1->GetBinContent(i)/hnewTeO2th230onlyM1->GetBinWidth(i));

    hAdapTeO2Sxth232onlyM1_001->SetBinContent(i, hnewTeO2Sxth232onlyM1_001->GetBinContent(i)/hnewTeO2Sxth232onlyM1_001->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M1_001->SetBinContent(i, hnewTeO2Sxra228pb208M1_001->GetBinContent(i)/hnewTeO2Sxra228pb208M1_001->GetBinWidth(i));
    hAdapTeO2Sxu238th230M1_001->SetBinContent(i, hnewTeO2Sxu238th230M1_001->GetBinContent(i)/hnewTeO2Sxu238th230M1_001->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM1_001->SetBinContent(i, hnewTeO2Sxth230onlyM1_001->GetBinContent(i)/hnewTeO2Sxth230onlyM1_001->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M1_001->SetBinContent(i, hnewTeO2Sxra226pb210M1_001->GetBinContent(i)/hnewTeO2Sxra226pb210M1_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M1_0001->SetBinContent(i, hnewTeO2Sxpb210M1_0001->GetBinContent(i)/hnewTeO2Sxpb210M1_0001->GetBinWidth(i));

    hAdapTeO2Spb210M1_01->SetBinContent(i, hnewTeO2Spb210M1_01->GetBinContent(i)/hnewTeO2Spb210M1_01->GetBinWidth(i));
    hAdapTeO2Spo210M1_001->SetBinContent(i, hnewTeO2Spo210M1_001->GetBinContent(i)/hnewTeO2Spo210M1_001->GetBinWidth(i));
    hAdapTeO2Spo210M1_01->SetBinContent(i, hnewTeO2Spo210M1_01->GetBinContent(i)/hnewTeO2Spo210M1_01->GetBinWidth(i));
    hAdapTeO2Sth232M1_01->SetBinContent(i, hnewTeO2Sth232M1_01->GetBinContent(i)/hnewTeO2Sth232M1_01->GetBinWidth(i));
    hAdapTeO2Su238M1_01->SetBinContent(i, hnewTeO2Su238M1_01->GetBinContent(i)/hnewTeO2Su238M1_01->GetBinWidth(i));
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


    hAdapCuFrameco58M1->SetBinContent(i, hnewCuFrameco58M1->GetBinContent(i)/hnewCuFrameco58M1->GetBinWidth(i));
    hAdapCuFrameco60M1->SetBinContent(i, hnewCuFrameco60M1->GetBinContent(i)/hnewCuFrameco60M1->GetBinWidth(i));
    hAdapCuFramecs137M1->SetBinContent(i, hnewCuFramecs137M1->GetBinContent(i)/hnewCuFramecs137M1->GetBinWidth(i));
    hAdapCuFramek40M1->SetBinContent(i, hnewCuFramek40M1->GetBinContent(i)/hnewCuFramek40M1->GetBinWidth(i));
    hAdapCuFramemn54M1->SetBinContent(i, hnewCuFramemn54M1->GetBinContent(i)/hnewCuFramemn54M1->GetBinWidth(i));
    hAdapCuFramepb210M1->SetBinContent(i, hnewCuFramepb210M1->GetBinContent(i)/hnewCuFramepb210M1->GetBinWidth(i));
    hAdapCuFrameth232M1->SetBinContent(i, hnewCuFrameth232M1->GetBinContent(i)/hnewCuFrameth232M1->GetBinWidth(i));
    hAdapCuFrameu238M1->SetBinContent(i, hnewCuFrameu238M1->GetBinContent(i)/hnewCuFrameu238M1->GetBinWidth(i));

    hAdapCuFrameSth232M1_1->SetBinContent(i, hnewCuFrameSth232M1_1->GetBinContent(i)/hnewCuFrameSth232M1_1->GetBinWidth(i));
    hAdapCuFrameSu238M1_1->SetBinContent(i, hnewCuFrameSu238M1_1->GetBinContent(i)/hnewCuFrameSu238M1_1->GetBinWidth(i));
    hAdapCuFrameSxpb210M1_001->SetBinContent(i, hnewCuFrameSxpb210M1_001->GetBinContent(i)/hnewCuFrameSxpb210M1_001->GetBinWidth(i));
    hAdapCuFrameSxpb210M1_01->SetBinContent(i, hnewCuFrameSxpb210M1_01->GetBinContent(i)/hnewCuFrameSxpb210M1_01->GetBinWidth(i));
    hAdapCuFrameSxpb210M1_1->SetBinContent(i, hnewCuFrameSxpb210M1_1->GetBinContent(i)/hnewCuFrameSxpb210M1_1->GetBinWidth(i));
    hAdapCuFrameSxpb210M1_10->SetBinContent(i, hnewCuFrameSxpb210M1_10->GetBinContent(i)/hnewCuFrameSxpb210M1_10->GetBinWidth(i));
    hAdapCuFrameSxth232M1_001->SetBinContent(i, hnewCuFrameSxth232M1_001->GetBinContent(i)/hnewCuFrameSxth232M1_001->GetBinWidth(i));
    hAdapCuFrameSxth232M1_01->SetBinContent(i, hnewCuFrameSxth232M1_01->GetBinContent(i)/hnewCuFrameSxth232M1_01->GetBinWidth(i));
    hAdapCuFrameSxth232M1_1->SetBinContent(i, hnewCuFrameSxth232M1_1->GetBinContent(i)/hnewCuFrameSxth232M1_1->GetBinWidth(i));
    hAdapCuFrameSxth232M1_10->SetBinContent(i, hnewCuFrameSxth232M1_10->GetBinContent(i)/hnewCuFrameSxth232M1_10->GetBinWidth(i));
    hAdapCuFrameSxu238M1_001->SetBinContent(i, hnewCuFrameSxu238M1_001->GetBinContent(i)/hnewCuFrameSxu238M1_001->GetBinWidth(i));
    hAdapCuFrameSxu238M1_01->SetBinContent(i, hnewCuFrameSxu238M1_01->GetBinContent(i)/hnewCuFrameSxu238M1_01->GetBinWidth(i));
    hAdapCuFrameSxu238M1_1->SetBinContent(i, hnewCuFrameSxu238M1_1->GetBinContent(i)/hnewCuFrameSxu238M1_1->GetBinWidth(i));
    hAdapCuFrameSxu238M1_10->SetBinContent(i, hnewCuFrameSxu238M1_10->GetBinContent(i)/hnewCuFrameSxu238M1_10->GetBinWidth(i));

    hAdapCuBoxco58M1->SetBinContent(i, hnewCuBoxco58M1->GetBinContent(i)/hnewCuBoxco58M1->GetBinWidth(i));
    hAdapCuBoxco60M1->SetBinContent(i, hnewCuBoxco60M1->GetBinContent(i)/hnewCuBoxco60M1->GetBinWidth(i));
    hAdapCuBoxcs137M1->SetBinContent(i, hnewCuBoxcs137M1->GetBinContent(i)/hnewCuBoxcs137M1->GetBinWidth(i));
    hAdapCuBoxk40M1->SetBinContent(i, hnewCuBoxk40M1->GetBinContent(i)/hnewCuBoxk40M1->GetBinWidth(i));
    hAdapCuBoxmn54M1->SetBinContent(i, hnewCuBoxmn54M1->GetBinContent(i)/hnewCuBoxmn54M1->GetBinWidth(i));
    hAdapCuBoxpb210M1->SetBinContent(i, hnewCuBoxpb210M1->GetBinContent(i)/hnewCuBoxpb210M1->GetBinWidth(i));
    hAdapCuBoxth232M1->SetBinContent(i, hnewCuBoxth232M1->GetBinContent(i)/hnewCuBoxth232M1->GetBinWidth(i));
    hAdapCuBoxu238M1->SetBinContent(i, hnewCuBoxu238M1->GetBinContent(i)/hnewCuBoxu238M1->GetBinWidth(i));

    hAdapCuBoxSth232M1_1->SetBinContent(i, hnewCuBoxSth232M1_1->GetBinContent(i)/hnewCuBoxSth232M1_1->GetBinWidth(i));
    hAdapCuBoxSu238M1_1->SetBinContent(i, hnewCuBoxSu238M1_1->GetBinContent(i)/hnewCuBoxSu238M1_1->GetBinWidth(i));
    hAdapCuBoxSxpb210M1_001->SetBinContent(i, hnewCuBoxSxpb210M1_001->GetBinContent(i)/hnewCuBoxSxpb210M1_001->GetBinWidth(i));
    hAdapCuBoxSxpb210M1_01->SetBinContent(i, hnewCuBoxSxpb210M1_01->GetBinContent(i)/hnewCuBoxSxpb210M1_01->GetBinWidth(i));
    hAdapCuBoxSxpb210M1_1->SetBinContent(i, hnewCuBoxSxpb210M1_1->GetBinContent(i)/hnewCuBoxSxpb210M1_1->GetBinWidth(i));
    hAdapCuBoxSxpb210M1_10->SetBinContent(i, hnewCuBoxSxpb210M1_10->GetBinContent(i)/hnewCuBoxSxpb210M1_10->GetBinWidth(i));
    hAdapCuBoxSxth232M1_001->SetBinContent(i, hnewCuBoxSxth232M1_001->GetBinContent(i)/hnewCuBoxSxth232M1_001->GetBinWidth(i));
    hAdapCuBoxSxth232M1_01->SetBinContent(i, hnewCuBoxSxth232M1_01->GetBinContent(i)/hnewCuBoxSxth232M1_01->GetBinWidth(i));
    hAdapCuBoxSxth232M1_1->SetBinContent(i, hnewCuBoxSxth232M1_1->GetBinContent(i)/hnewCuBoxSxth232M1_1->GetBinWidth(i));
    hAdapCuBoxSxth232M1_10->SetBinContent(i, hnewCuBoxSxth232M1_10->GetBinContent(i)/hnewCuBoxSxth232M1_10->GetBinWidth(i));
    hAdapCuBoxSxu238M1_001->SetBinContent(i, hnewCuBoxSxu238M1_001->GetBinContent(i)/hnewCuBoxSxu238M1_001->GetBinWidth(i));
    hAdapCuBoxSxu238M1_01->SetBinContent(i, hnewCuBoxSxu238M1_01->GetBinContent(i)/hnewCuBoxSxu238M1_01->GetBinWidth(i));
    hAdapCuBoxSxu238M1_1->SetBinContent(i, hnewCuBoxSxu238M1_1->GetBinContent(i)/hnewCuBoxSxu238M1_1->GetBinWidth(i));
    hAdapCuBoxSxu238M1_10->SetBinContent(i, hnewCuBoxSxu238M1_10->GetBinContent(i)/hnewCuBoxSxu238M1_10->GetBinWidth(i));

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

/*
    hAdap50mKco58M1->SetBinContent(i, hnew50mKco58M1->GetBinContent(i)/hnew50mKco58M1->GetBinWidth(i));
    hAdap50mKco60M1->SetBinContent(i, hnew50mKco60M1->GetBinContent(i)/hnew50mKco60M1->GetBinWidth(i));
    hAdap50mKcs137M1->SetBinContent(i, hnew50mKcs137M1->GetBinContent(i)/hnew50mKcs137M1->GetBinWidth(i));
    hAdap50mKk40M1->SetBinContent(i, hnew50mKk40M1->GetBinContent(i)/hnew50mKk40M1->GetBinWidth(i));
    hAdap50mKmn54M1->SetBinContent(i, hnew50mKmn54M1->GetBinContent(i)/hnew50mKmn54M1->GetBinWidth(i));
    hAdap50mKpb210M1->SetBinContent(i, hnew50mKpb210M1->GetBinContent(i)/hnew50mKpb210M1->GetBinWidth(i));
    hAdap50mKth232M1->SetBinContent(i, hnew50mKth232M1->GetBinContent(i)/hnew50mKth232M1->GetBinWidth(i));
    hAdap50mKu238M1->SetBinContent(i, hnew50mKu238M1->GetBinContent(i)/hnew50mKu238M1->GetBinWidth(i));

    hAdap600mKco60M1->SetBinContent(i, hnew600mKco60M1->GetBinContent(i)/hnew600mKco60M1->GetBinWidth(i));
    hAdap600mKk40M1->SetBinContent(i, hnew600mKk40M1->GetBinContent(i)/hnew600mKk40M1->GetBinWidth(i));
    hAdap600mKth232M1->SetBinContent(i, hnew600mKth232M1->GetBinContent(i)/hnew600mKth232M1->GetBinWidth(i));
    hAdap600mKu238M1->SetBinContent(i, hnew600mKu238M1->GetBinContent(i)/hnew600mKu238M1->GetBinWidth(i));
*/

    // hAdapPbRombi207M1->SetBinContent(i, hnewPbRombi207M1->GetBinContent(i)/hnewPbRombi207M1->GetBinWidth(i));
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

/*
    hAdapMBco60M1->SetBinContent(i, hnewMBco60M1->GetBinContent(i)/hnewMBco60M1->GetBinWidth(i));
    hAdapMBk40M1->SetBinContent(i, hnewMBk40M1->GetBinContent(i)/hnewMBk40M1->GetBinWidth(i));
    hAdapMBth232M1->SetBinContent(i, hnewMBth232M1->GetBinContent(i)/hnewMBth232M1->GetBinWidth(i));
    hAdapMBu238M1->SetBinContent(i, hnewMBu238M1->GetBinContent(i)/hnewMBu238M1->GetBinWidth(i));

    hAdapSIk40M1->SetBinContent(i, hnewSIk40M1->GetBinContent(i)/hnewSIk40M1->GetBinWidth(i));
    hAdapSIth232M1->SetBinContent(i, hnewSIth232M1->GetBinContent(i)/hnewSIth232M1->GetBinWidth(i));
    hAdapSIu238M1->SetBinContent(i, hnewSIu238M1->GetBinContent(i)/hnewSIu238M1->GetBinWidth(i));

    hAdapIVCco60M1->SetBinContent(i, hnewIVCco60M1->GetBinContent(i)/hnewIVCco60M1->GetBinWidth(i));
    hAdapIVCk40M1->SetBinContent(i, hnewIVCk40M1->GetBinContent(i)/hnewIVCk40M1->GetBinWidth(i));
    hAdapIVCth232M1->SetBinContent(i, hnewIVCth232M1->GetBinContent(i)/hnewIVCth232M1->GetBinWidth(i));
    hAdapIVCu238M1->SetBinContent(i, hnewIVCu238M1->GetBinContent(i)/hnewIVCu238M1->GetBinWidth(i));
*/
    hAdapOVCco60M1->SetBinContent(i, hnewOVCco60M1->GetBinContent(i)/hnewOVCco60M1->GetBinWidth(i));
    hAdapOVCk40M1->SetBinContent(i, hnewOVCk40M1->GetBinContent(i)/hnewOVCk40M1->GetBinWidth(i));
    hAdapOVCth232M1->SetBinContent(i, hnewOVCth232M1->GetBinContent(i)/hnewOVCth232M1->GetBinWidth(i));
    hAdapOVCu238M1->SetBinContent(i, hnewOVCu238M1->GetBinContent(i)/hnewOVCu238M1->GetBinWidth(i));

    hAdapExtPbbi210M1->SetBinContent(i, hnewExtPbbi210M1->GetBinContent(i)/hnewExtPbbi210M1->GetBinWidth(i));
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

    hAdapTeO2Spb210M2_01->SetBinContent(i, hnewTeO2Spb210M2_01->GetBinContent(i)/hnewTeO2Spb210M2_01->GetBinWidth(i));
    hAdapTeO2Spo210M2_001->SetBinContent(i, hnewTeO2Spo210M2_001->GetBinContent(i)/hnewTeO2Spo210M2_001->GetBinWidth(i));
    hAdapTeO2Spo210M2_01->SetBinContent(i, hnewTeO2Spo210M2_01->GetBinContent(i)/hnewTeO2Spo210M2_01->GetBinWidth(i));
    hAdapTeO2Sth232M2_01->SetBinContent(i, hnewTeO2Sth232M2_01->GetBinContent(i)/hnewTeO2Sth232M2_01->GetBinWidth(i));
    hAdapTeO2Su238M2_01->SetBinContent(i, hnewTeO2Su238M2_01->GetBinContent(i)/hnewTeO2Su238M2_01->GetBinWidth(i));
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

    hAdapTeO2Sxth232onlyM2_001->SetBinContent(i, hnewTeO2Sxth232onlyM2_001->GetBinContent(i)/hnewTeO2Sxth232onlyM2_001->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M2_001->SetBinContent(i, hnewTeO2Sxra228pb208M2_001->GetBinContent(i)/hnewTeO2Sxra228pb208M2_001->GetBinWidth(i));
    hAdapTeO2Sxu238th230M2_001->SetBinContent(i, hnewTeO2Sxu238th230M2_001->GetBinContent(i)/hnewTeO2Sxu238th230M2_001->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM2_001->SetBinContent(i, hnewTeO2Sxth230onlyM2_001->GetBinContent(i)/hnewTeO2Sxth230onlyM2_001->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M2_001->SetBinContent(i, hnewTeO2Sxra226pb210M2_001->GetBinContent(i)/hnewTeO2Sxra226pb210M2_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_0001->SetBinContent(i, hnewTeO2Sxpb210M2_0001->GetBinContent(i)/hnewTeO2Sxpb210M2_0001->GetBinWidth(i));


    hAdapCuFrameco58M2->SetBinContent(i, hnewCuFrameco58M2->GetBinContent(i)/hnewCuFrameco58M2->GetBinWidth(i));
    hAdapCuFrameco60M2->SetBinContent(i, hnewCuFrameco60M2->GetBinContent(i)/hnewCuFrameco60M2->GetBinWidth(i));
    hAdapCuFramecs137M2->SetBinContent(i, hnewCuFramecs137M2->GetBinContent(i)/hnewCuFramecs137M2->GetBinWidth(i));
    hAdapCuFramek40M2->SetBinContent(i, hnewCuFramek40M2->GetBinContent(i)/hnewCuFramek40M2->GetBinWidth(i));
    hAdapCuFramemn54M2->SetBinContent(i, hnewCuFramemn54M2->GetBinContent(i)/hnewCuFramemn54M2->GetBinWidth(i));
    hAdapCuFramepb210M2->SetBinContent(i, hnewCuFramepb210M2->GetBinContent(i)/hnewCuFramepb210M2->GetBinWidth(i));
    hAdapCuFrameth232M2->SetBinContent(i, hnewCuFrameth232M2->GetBinContent(i)/hnewCuFrameth232M2->GetBinWidth(i));
    hAdapCuFrameu238M2->SetBinContent(i, hnewCuFrameu238M2->GetBinContent(i)/hnewCuFrameu238M2->GetBinWidth(i));

    hAdapCuFrameSth232M2_1->SetBinContent(i, hnewCuFrameSth232M2_1->GetBinContent(i)/hnewCuFrameSth232M2_1->GetBinWidth(i));
    hAdapCuFrameSu238M2_1->SetBinContent(i, hnewCuFrameSu238M2_1->GetBinContent(i)/hnewCuFrameSu238M2_1->GetBinWidth(i));
    hAdapCuFrameSxpb210M2_001->SetBinContent(i, hnewCuFrameSxpb210M2_001->GetBinContent(i)/hnewCuFrameSxpb210M2_001->GetBinWidth(i));
    hAdapCuFrameSxpb210M2_01->SetBinContent(i, hnewCuFrameSxpb210M2_01->GetBinContent(i)/hnewCuFrameSxpb210M2_01->GetBinWidth(i));
    hAdapCuFrameSxpb210M2_1->SetBinContent(i, hnewCuFrameSxpb210M2_1->GetBinContent(i)/hnewCuFrameSxpb210M2_1->GetBinWidth(i));
    hAdapCuFrameSxpb210M2_10->SetBinContent(i, hnewCuFrameSxpb210M2_10->GetBinContent(i)/hnewCuFrameSxpb210M2_10->GetBinWidth(i));
    hAdapCuFrameSxth232M2_001->SetBinContent(i, hnewCuFrameSxth232M2_001->GetBinContent(i)/hnewCuFrameSxth232M2_001->GetBinWidth(i));
    hAdapCuFrameSxth232M2_01->SetBinContent(i, hnewCuFrameSxth232M2_01->GetBinContent(i)/hnewCuFrameSxth232M2_01->GetBinWidth(i));
    hAdapCuFrameSxth232M2_1->SetBinContent(i, hnewCuFrameSxth232M2_1->GetBinContent(i)/hnewCuFrameSxth232M2_1->GetBinWidth(i));
    hAdapCuFrameSxth232M2_10->SetBinContent(i, hnewCuFrameSxth232M2_10->GetBinContent(i)/hnewCuFrameSxth232M2_10->GetBinWidth(i));
    hAdapCuFrameSxu238M2_001->SetBinContent(i, hnewCuFrameSxu238M2_001->GetBinContent(i)/hnewCuFrameSxu238M2_001->GetBinWidth(i));
    hAdapCuFrameSxu238M2_01->SetBinContent(i, hnewCuFrameSxu238M2_01->GetBinContent(i)/hnewCuFrameSxu238M2_01->GetBinWidth(i));
    hAdapCuFrameSxu238M2_1->SetBinContent(i, hnewCuFrameSxu238M2_1->GetBinContent(i)/hnewCuFrameSxu238M2_1->GetBinWidth(i));
    hAdapCuFrameSxu238M2_10->SetBinContent(i, hnewCuFrameSxu238M2_10->GetBinContent(i)/hnewCuFrameSxu238M2_10->GetBinWidth(i));

    hAdapCuBoxco58M2->SetBinContent(i, hnewCuBoxco58M2->GetBinContent(i)/hnewCuBoxco58M2->GetBinWidth(i));
    hAdapCuBoxco60M2->SetBinContent(i, hnewCuBoxco60M2->GetBinContent(i)/hnewCuBoxco60M2->GetBinWidth(i));
    hAdapCuBoxcs137M2->SetBinContent(i, hnewCuBoxcs137M2->GetBinContent(i)/hnewCuBoxcs137M2->GetBinWidth(i));
    hAdapCuBoxk40M2->SetBinContent(i, hnewCuBoxk40M2->GetBinContent(i)/hnewCuBoxk40M2->GetBinWidth(i));
    hAdapCuBoxmn54M2->SetBinContent(i, hnewCuBoxmn54M2->GetBinContent(i)/hnewCuBoxmn54M2->GetBinWidth(i));
    hAdapCuBoxpb210M2->SetBinContent(i, hnewCuBoxpb210M2->GetBinContent(i)/hnewCuBoxpb210M2->GetBinWidth(i));
    hAdapCuBoxth232M2->SetBinContent(i, hnewCuBoxth232M2->GetBinContent(i)/hnewCuBoxth232M2->GetBinWidth(i));
    hAdapCuBoxu238M2->SetBinContent(i, hnewCuBoxu238M2->GetBinContent(i)/hnewCuBoxu238M2->GetBinWidth(i));

    hAdapCuBoxSth232M2_1->SetBinContent(i, hnewCuBoxSth232M2_1->GetBinContent(i)/hnewCuBoxSth232M2_1->GetBinWidth(i));
    hAdapCuBoxSu238M2_1->SetBinContent(i, hnewCuBoxSu238M2_1->GetBinContent(i)/hnewCuBoxSu238M2_1->GetBinWidth(i));
    hAdapCuBoxSxpb210M2_001->SetBinContent(i, hnewCuBoxSxpb210M2_001->GetBinContent(i)/hnewCuBoxSxpb210M2_001->GetBinWidth(i));
    hAdapCuBoxSxpb210M2_01->SetBinContent(i, hnewCuBoxSxpb210M2_01->GetBinContent(i)/hnewCuBoxSxpb210M2_01->GetBinWidth(i));
    hAdapCuBoxSxpb210M2_1->SetBinContent(i, hnewCuBoxSxpb210M2_1->GetBinContent(i)/hnewCuBoxSxpb210M2_1->GetBinWidth(i));
    hAdapCuBoxSxpb210M2_10->SetBinContent(i, hnewCuBoxSxpb210M2_10->GetBinContent(i)/hnewCuBoxSxpb210M2_10->GetBinWidth(i));
    hAdapCuBoxSxth232M2_001->SetBinContent(i, hnewCuBoxSxth232M2_001->GetBinContent(i)/hnewCuBoxSxth232M2_001->GetBinWidth(i));
    hAdapCuBoxSxth232M2_01->SetBinContent(i, hnewCuBoxSxth232M2_01->GetBinContent(i)/hnewCuBoxSxth232M2_01->GetBinWidth(i));
    hAdapCuBoxSxth232M2_1->SetBinContent(i, hnewCuBoxSxth232M2_1->GetBinContent(i)/hnewCuBoxSxth232M2_1->GetBinWidth(i));
    hAdapCuBoxSxth232M2_10->SetBinContent(i, hnewCuBoxSxth232M2_10->GetBinContent(i)/hnewCuBoxSxth232M2_10->GetBinWidth(i));
    hAdapCuBoxSxu238M2_001->SetBinContent(i, hnewCuBoxSxu238M2_001->GetBinContent(i)/hnewCuBoxSxu238M2_001->GetBinWidth(i));
    hAdapCuBoxSxu238M2_01->SetBinContent(i, hnewCuBoxSxu238M2_01->GetBinContent(i)/hnewCuBoxSxu238M2_01->GetBinWidth(i));
    hAdapCuBoxSxu238M2_1->SetBinContent(i, hnewCuBoxSxu238M2_1->GetBinContent(i)/hnewCuBoxSxu238M2_1->GetBinWidth(i));
    hAdapCuBoxSxu238M2_10->SetBinContent(i, hnewCuBoxSxu238M2_10->GetBinContent(i)/hnewCuBoxSxu238M2_10->GetBinWidth(i));


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
/*
    hAdap50mKco58M2->SetBinContent(i, hnew50mKco58M2->GetBinContent(i)/hnew50mKco58M2->GetBinWidth(i));
    hAdap50mKco60M2->SetBinContent(i, hnew50mKco60M2->GetBinContent(i)/hnew50mKco60M2->GetBinWidth(i));
    hAdap50mKcs137M2->SetBinContent(i, hnew50mKcs137M2->GetBinContent(i)/hnew50mKcs137M2->GetBinWidth(i));
    hAdap50mKk40M2->SetBinContent(i, hnew50mKk40M2->GetBinContent(i)/hnew50mKk40M2->GetBinWidth(i));
    hAdap50mKmn54M2->SetBinContent(i, hnew50mKmn54M2->GetBinContent(i)/hnew50mKmn54M2->GetBinWidth(i));
    hAdap50mKpb210M2->SetBinContent(i, hnew50mKpb210M2->GetBinContent(i)/hnew50mKpb210M2->GetBinWidth(i));
    hAdap50mKth232M2->SetBinContent(i, hnew50mKth232M2->GetBinContent(i)/hnew50mKth232M2->GetBinWidth(i));
    hAdap50mKu238M2->SetBinContent(i, hnew50mKu238M2->GetBinContent(i)/hnew50mKu238M2->GetBinWidth(i));

    hAdap600mKco60M2->SetBinContent(i, hnew600mKco60M2->GetBinContent(i)/hnew600mKco60M2->GetBinWidth(i));
    hAdap600mKk40M2->SetBinContent(i, hnew600mKk40M2->GetBinContent(i)/hnew600mKk40M2->GetBinWidth(i));
    hAdap600mKth232M2->SetBinContent(i, hnew600mKth232M2->GetBinContent(i)/hnew600mKth232M2->GetBinWidth(i));
    hAdap600mKu238M2->SetBinContent(i, hnew600mKu238M2->GetBinContent(i)/hnew600mKu238M2->GetBinWidth(i));
*/

    // hAdapPbRombi207M2->SetBinContent(i, hnewPbRombi207M2->GetBinContent(i)/hnewPbRombi207M2->GetBinWidth(i));
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

/*
    hAdapMBco60M2->SetBinContent(i, hnewMBco60M2->GetBinContent(i)/hnewMBco60M2->GetBinWidth(i));
    hAdapMBk40M2->SetBinContent(i, hnewMBk40M2->GetBinContent(i)/hnewMBk40M2->GetBinWidth(i));
    hAdapMBth232M2->SetBinContent(i, hnewMBth232M2->GetBinContent(i)/hnewMBth232M2->GetBinWidth(i));
    hAdapMBu238M2->SetBinContent(i, hnewMBu238M2->GetBinContent(i)/hnewMBu238M2->GetBinWidth(i));

    hAdapSIk40M2->SetBinContent(i, hnewSIk40M2->GetBinContent(i)/hnewSIk40M2->GetBinWidth(i));
    hAdapSIth232M2->SetBinContent(i, hnewSIth232M2->GetBinContent(i)/hnewSIth232M2->GetBinWidth(i));
    hAdapSIu238M2->SetBinContent(i, hnewSIu238M2->GetBinContent(i)/hnewSIu238M2->GetBinWidth(i));


    hAdapIVCco60M2->SetBinContent(i, hnewIVCco60M2->GetBinContent(i)/hnewIVCco60M2->GetBinWidth(i));
    hAdapIVCk40M2->SetBinContent(i, hnewIVCk40M2->GetBinContent(i)/hnewIVCk40M2->GetBinWidth(i));
    hAdapIVCth232M2->SetBinContent(i, hnewIVCth232M2->GetBinContent(i)/hnewIVCth232M2->GetBinWidth(i));
    hAdapIVCu238M2->SetBinContent(i, hnewIVCu238M2->GetBinContent(i)/hnewIVCu238M2->GetBinWidth(i));
*/

    hAdapOVCco60M2->SetBinContent(i, hnewOVCco60M2->GetBinContent(i)/hnewOVCco60M2->GetBinWidth(i));
    hAdapOVCk40M2->SetBinContent(i, hnewOVCk40M2->GetBinContent(i)/hnewOVCk40M2->GetBinWidth(i));
    hAdapOVCth232M2->SetBinContent(i, hnewOVCth232M2->GetBinContent(i)/hnewOVCth232M2->GetBinWidth(i));
    hAdapOVCu238M2->SetBinContent(i, hnewOVCu238M2->GetBinContent(i)/hnewOVCu238M2->GetBinWidth(i));

    hAdapExtPbbi210M2->SetBinContent(i, hnewExtPbbi210M2->GetBinContent(i)/hnewExtPbbi210M2->GetBinWidth(i));
  }

/*
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

    hAdapTeO2Spb210M2Sum_01->SetBinContent(i, hnewTeO2Spb210M2Sum_01->GetBinContent(i)/hnewTeO2Spb210M2Sum_01->GetBinWidth(i));
    hAdapTeO2Spo210M2Sum_001->SetBinContent(i, hnewTeO2Spo210M2Sum_001->GetBinContent(i)/hnewTeO2Spo210M2Sum_001->GetBinWidth(i));
    hAdapTeO2Spo210M2Sum_01->SetBinContent(i, hnewTeO2Spo210M2Sum_01->GetBinContent(i)/hnewTeO2Spo210M2Sum_01->GetBinWidth(i));
    hAdapTeO2Sth232M2Sum_01->SetBinContent(i, hnewTeO2Sth232M2Sum_01->GetBinContent(i)/hnewTeO2Sth232M2Sum_01->GetBinWidth(i));
    hAdapTeO2Su238M2Sum_01->SetBinContent(i, hnewTeO2Su238M2Sum_01->GetBinContent(i)/hnewTeO2Su238M2Sum_01->GetBinWidth(i));
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

    hAdapTeO2th232onlyM2Sum->SetBinContent(i, hnewTeO2th232onlyM2Sum->GetBinContent(i)/hnewTeO2th232onlyM2Sum->GetBinWidth(i));
    hAdapTeO2ra228pb208M2Sum->SetBinContent(i, hnewTeO2ra228pb208M2Sum->GetBinContent(i)/hnewTeO2ra228pb208M2Sum->GetBinWidth(i));
    hAdapTeO2th230onlyM2Sum->SetBinContent(i, hnewTeO2th230onlyM2Sum->GetBinContent(i)/hnewTeO2th230onlyM2Sum->GetBinWidth(i));

    hAdapTeO2Sxth232onlyM2Sum_001->SetBinContent(i, hnewTeO2Sxth232onlyM2Sum_001->GetBinContent(i)/hnewTeO2Sxth232onlyM2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M2Sum_001->SetBinContent(i, hnewTeO2Sxra228pb208M2Sum_001->GetBinContent(i)/hnewTeO2Sxra228pb208M2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxu238th230M2Sum_001->SetBinContent(i, hnewTeO2Sxu238th230M2Sum_001->GetBinContent(i)/hnewTeO2Sxu238th230M2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM2Sum_001->SetBinContent(i, hnewTeO2Sxth230onlyM2Sum_001->GetBinContent(i)/hnewTeO2Sxth230onlyM2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M2Sum_001->SetBinContent(i, hnewTeO2Sxra226pb210M2Sum_001->GetBinContent(i)/hnewTeO2Sxra226pb210M2Sum_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M2Sum_0001->SetBinContent(i, hnewTeO2Sxpb210M2Sum_0001->GetBinContent(i)/hnewTeO2Sxpb210M2Sum_0001->GetBinWidth(i));


    hAdapCuFrameco58M2Sum->SetBinContent(i, hnewCuFrameco58M2Sum->GetBinContent(i)/hnewCuFrameco58M2Sum->GetBinWidth(i));
    hAdapCuFrameco60M2Sum->SetBinContent(i, hnewCuFrameco60M2Sum->GetBinContent(i)/hnewCuFrameco60M2Sum->GetBinWidth(i));
    hAdapCuFramecs137M2Sum->SetBinContent(i, hnewCuFramecs137M2Sum->GetBinContent(i)/hnewCuFramecs137M2Sum->GetBinWidth(i));
    hAdapCuFramek40M2Sum->SetBinContent(i, hnewCuFramek40M2Sum->GetBinContent(i)/hnewCuFramek40M2Sum->GetBinWidth(i));
    hAdapCuFramemn54M2Sum->SetBinContent(i, hnewCuFramemn54M2Sum->GetBinContent(i)/hnewCuFramemn54M2Sum->GetBinWidth(i));
    hAdapCuFramepb210M2Sum->SetBinContent(i, hnewCuFramepb210M2Sum->GetBinContent(i)/hnewCuFramepb210M2Sum->GetBinWidth(i));
    hAdapCuFrameth232M2Sum->SetBinContent(i, hnewCuFrameth232M2Sum->GetBinContent(i)/hnewCuFrameth232M2Sum->GetBinWidth(i));
    hAdapCuFrameu238M2Sum->SetBinContent(i, hnewCuFrameu238M2Sum->GetBinContent(i)/hnewCuFrameu238M2Sum->GetBinWidth(i));

    hAdapCuFrameSth232M2Sum_1->SetBinContent(i, hnewCuFrameSth232M2Sum_1->GetBinContent(i)/hnewCuFrameSth232M2Sum_1->GetBinWidth(i));
    hAdapCuFrameSu238M2Sum_1->SetBinContent(i, hnewCuFrameSu238M2Sum_1->GetBinContent(i)/hnewCuFrameSu238M2Sum_1->GetBinWidth(i));
    hAdapCuFrameSxpb210M2Sum_001->SetBinContent(i, hnewCuFrameSxpb210M2Sum_001->GetBinContent(i)/hnewCuFrameSxpb210M2Sum_001->GetBinWidth(i));
    hAdapCuFrameSxpb210M2Sum_01->SetBinContent(i, hnewCuFrameSxpb210M2Sum_01->GetBinContent(i)/hnewCuFrameSxpb210M2Sum_01->GetBinWidth(i));
    hAdapCuFrameSxpb210M2Sum_1->SetBinContent(i, hnewCuFrameSxpb210M2Sum_1->GetBinContent(i)/hnewCuFrameSxpb210M2Sum_1->GetBinWidth(i));
    hAdapCuFrameSxpb210M2Sum_10->SetBinContent(i, hnewCuFrameSxpb210M2Sum_10->GetBinContent(i)/hnewCuFrameSxpb210M2Sum_10->GetBinWidth(i));
    hAdapCuFrameSxth232M2Sum_001->SetBinContent(i, hnewCuFrameSxth232M2Sum_001->GetBinContent(i)/hnewCuFrameSxth232M2Sum_001->GetBinWidth(i));
    hAdapCuFrameSxth232M2Sum_01->SetBinContent(i, hnewCuFrameSxth232M2Sum_01->GetBinContent(i)/hnewCuFrameSxth232M2Sum_01->GetBinWidth(i));
    hAdapCuFrameSxth232M2Sum_1->SetBinContent(i, hnewCuFrameSxth232M2Sum_1->GetBinContent(i)/hnewCuFrameSxth232M2Sum_1->GetBinWidth(i));
    hAdapCuFrameSxth232M2Sum_10->SetBinContent(i, hnewCuFrameSxth232M2Sum_10->GetBinContent(i)/hnewCuFrameSxth232M2Sum_10->GetBinWidth(i));
    hAdapCuFrameSxu238M2Sum_001->SetBinContent(i, hnewCuFrameSxu238M2Sum_001->GetBinContent(i)/hnewCuFrameSxu238M2Sum_001->GetBinWidth(i));
    hAdapCuFrameSxu238M2Sum_01->SetBinContent(i, hnewCuFrameSxu238M2Sum_01->GetBinContent(i)/hnewCuFrameSxu238M2Sum_01->GetBinWidth(i));
    hAdapCuFrameSxu238M2Sum_1->SetBinContent(i, hnewCuFrameSxu238M2Sum_1->GetBinContent(i)/hnewCuFrameSxu238M2Sum_1->GetBinWidth(i));
    hAdapCuFrameSxu238M2Sum_10->SetBinContent(i, hnewCuFrameSxu238M2Sum_10->GetBinContent(i)/hnewCuFrameSxu238M2Sum_10->GetBinWidth(i));

    hAdapCuBoxco58M2Sum->SetBinContent(i, hnewCuBoxco58M2Sum->GetBinContent(i)/hnewCuBoxco58M2Sum->GetBinWidth(i));
    hAdapCuBoxco60M2Sum->SetBinContent(i, hnewCuBoxco60M2Sum->GetBinContent(i)/hnewCuBoxco60M2Sum->GetBinWidth(i));
    hAdapCuBoxcs137M2Sum->SetBinContent(i, hnewCuBoxcs137M2Sum->GetBinContent(i)/hnewCuBoxcs137M2Sum->GetBinWidth(i));
    hAdapCuBoxk40M2Sum->SetBinContent(i, hnewCuBoxk40M2Sum->GetBinContent(i)/hnewCuBoxk40M2Sum->GetBinWidth(i));
    hAdapCuBoxmn54M2Sum->SetBinContent(i, hnewCuBoxmn54M2Sum->GetBinContent(i)/hnewCuBoxmn54M2Sum->GetBinWidth(i));
    hAdapCuBoxpb210M2Sum->SetBinContent(i, hnewCuBoxpb210M2Sum->GetBinContent(i)/hnewCuBoxpb210M2Sum->GetBinWidth(i));
    hAdapCuBoxth232M2Sum->SetBinContent(i, hnewCuBoxth232M2Sum->GetBinContent(i)/hnewCuBoxth232M2Sum->GetBinWidth(i));
    hAdapCuBoxu238M2Sum->SetBinContent(i, hnewCuBoxu238M2Sum->GetBinContent(i)/hnewCuBoxu238M2Sum->GetBinWidth(i));

    hAdapCuBoxSth232M2Sum_1->SetBinContent(i, hnewCuBoxSth232M2Sum_1->GetBinContent(i)/hnewCuBoxSth232M2Sum_1->GetBinWidth(i));
    hAdapCuBoxSu238M2Sum_1->SetBinContent(i, hnewCuBoxSu238M2Sum_1->GetBinContent(i)/hnewCuBoxSu238M2Sum_1->GetBinWidth(i));
    hAdapCuBoxSxpb210M2Sum_001->SetBinContent(i, hnewCuBoxSxpb210M2Sum_001->GetBinContent(i)/hnewCuBoxSxpb210M2Sum_001->GetBinWidth(i));
    hAdapCuBoxSxpb210M2Sum_01->SetBinContent(i, hnewCuBoxSxpb210M2Sum_01->GetBinContent(i)/hnewCuBoxSxpb210M2Sum_01->GetBinWidth(i));
    hAdapCuBoxSxpb210M2Sum_1->SetBinContent(i, hnewCuBoxSxpb210M2Sum_1->GetBinContent(i)/hnewCuBoxSxpb210M2Sum_1->GetBinWidth(i));
    hAdapCuBoxSxpb210M2Sum_10->SetBinContent(i, hnewCuBoxSxpb210M2Sum_10->GetBinContent(i)/hnewCuBoxSxpb210M2Sum_10->GetBinWidth(i));
    hAdapCuBoxSxth232M2Sum_001->SetBinContent(i, hnewCuBoxSxth232M2Sum_001->GetBinContent(i)/hnewCuBoxSxth232M2Sum_001->GetBinWidth(i));
    hAdapCuBoxSxth232M2Sum_01->SetBinContent(i, hnewCuBoxSxth232M2Sum_01->GetBinContent(i)/hnewCuBoxSxth232M2Sum_01->GetBinWidth(i));
    hAdapCuBoxSxth232M2Sum_1->SetBinContent(i, hnewCuBoxSxth232M2Sum_1->GetBinContent(i)/hnewCuBoxSxth232M2Sum_1->GetBinWidth(i));
    hAdapCuBoxSxth232M2Sum_10->SetBinContent(i, hnewCuBoxSxth232M2Sum_10->GetBinContent(i)/hnewCuBoxSxth232M2Sum_10->GetBinWidth(i));
    hAdapCuBoxSxu238M2Sum_001->SetBinContent(i, hnewCuBoxSxu238M2Sum_001->GetBinContent(i)/hnewCuBoxSxu238M2Sum_001->GetBinWidth(i));
    hAdapCuBoxSxu238M2Sum_01->SetBinContent(i, hnewCuBoxSxu238M2Sum_01->GetBinContent(i)/hnewCuBoxSxu238M2Sum_01->GetBinWidth(i));
    hAdapCuBoxSxu238M2Sum_1->SetBinContent(i, hnewCuBoxSxu238M2Sum_1->GetBinContent(i)/hnewCuBoxSxu238M2Sum_1->GetBinWidth(i));
    hAdapCuBoxSxu238M2Sum_10->SetBinContent(i, hnewCuBoxSxu238M2Sum_10->GetBinContent(i)/hnewCuBoxSxu238M2Sum_10->GetBinWidth(i));


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

    hAdap50mKco58M2Sum->SetBinContent(i, hnew50mKco58M2Sum->GetBinContent(i)/hnew50mKco58M2Sum->GetBinWidth(i));
    hAdap50mKco60M2Sum->SetBinContent(i, hnew50mKco60M2Sum->GetBinContent(i)/hnew50mKco60M2Sum->GetBinWidth(i));
    hAdap50mKcs137M2Sum->SetBinContent(i, hnew50mKcs137M2Sum->GetBinContent(i)/hnew50mKcs137M2Sum->GetBinWidth(i));
    hAdap50mKk40M2Sum->SetBinContent(i, hnew50mKk40M2Sum->GetBinContent(i)/hnew50mKk40M2Sum->GetBinWidth(i));
    hAdap50mKmn54M2Sum->SetBinContent(i, hnew50mKmn54M2Sum->GetBinContent(i)/hnew50mKmn54M2Sum->GetBinWidth(i));
    hAdap50mKpb210M2Sum->SetBinContent(i, hnew50mKpb210M2Sum->GetBinContent(i)/hnew50mKpb210M2Sum->GetBinWidth(i));
    hAdap50mKth232M2Sum->SetBinContent(i, hnew50mKth232M2Sum->GetBinContent(i)/hnew50mKth232M2Sum->GetBinWidth(i));
    hAdap50mKu238M2Sum->SetBinContent(i, hnew50mKu238M2Sum->GetBinContent(i)/hnew50mKu238M2Sum->GetBinWidth(i));

    hAdap600mKco60M2Sum->SetBinContent(i, hnew600mKco60M2Sum->GetBinContent(i)/hnew600mKco60M2Sum->GetBinWidth(i));
    hAdap600mKk40M2Sum->SetBinContent(i, hnew600mKk40M2Sum->GetBinContent(i)/hnew600mKk40M2Sum->GetBinWidth(i));
    hAdap600mKth232M2Sum->SetBinContent(i, hnew600mKth232M2Sum->GetBinContent(i)/hnew600mKth232M2Sum->GetBinWidth(i));
    hAdap600mKu238M2Sum->SetBinContent(i, hnew600mKu238M2Sum->GetBinContent(i)/hnew600mKu238M2Sum->GetBinWidth(i));


    // hAdapPbRombi207M2Sum->SetBinContent(i, hnewPbRombi207M2Sum->GetBinContent(i)/hnewPbRombi207M2Sum->GetBinWidth(i));
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


    hAdapMBco60M2Sum->SetBinContent(i, hnewMBco60M2Sum->GetBinContent(i)/hnewMBco60M2Sum->GetBinWidth(i));
    hAdapMBk40M2Sum->SetBinContent(i, hnewMBk40M2Sum->GetBinContent(i)/hnewMBk40M2Sum->GetBinWidth(i));
    hAdapMBth232M2Sum->SetBinContent(i, hnewMBth232M2Sum->GetBinContent(i)/hnewMBth232M2Sum->GetBinWidth(i));
    hAdapMBu238M2Sum->SetBinContent(i, hnewMBu238M2Sum->GetBinContent(i)/hnewMBu238M2Sum->GetBinWidth(i));

    hAdapSIk40M2Sum->SetBinContent(i, hnewSIk40M2Sum->GetBinContent(i)/hnewSIk40M2Sum->GetBinWidth(i));
    hAdapSIth232M2Sum->SetBinContent(i, hnewSIth232M2Sum->GetBinContent(i)/hnewSIth232M2Sum->GetBinWidth(i));
    hAdapSIu238M2Sum->SetBinContent(i, hnewSIu238M2Sum->GetBinContent(i)/hnewSIu238M2Sum->GetBinWidth(i));

    hAdapIVCco60M2Sum->SetBinContent(i, hnewIVCco60M2Sum->GetBinContent(i)/hnewIVCco60M2Sum->GetBinWidth(i));
    hAdapIVCk40M2Sum->SetBinContent(i, hnewIVCk40M2Sum->GetBinContent(i)/hnewIVCk40M2Sum->GetBinWidth(i));
    hAdapIVCth232M2Sum->SetBinContent(i, hnewIVCth232M2Sum->GetBinContent(i)/hnewIVCth232M2Sum->GetBinWidth(i));
    hAdapIVCu238M2Sum->SetBinContent(i, hnewIVCu238M2Sum->GetBinContent(i)/hnewIVCu238M2Sum->GetBinWidth(i));


    hAdapOVCco60M2Sum->SetBinContent(i, hnewOVCco60M2Sum->GetBinContent(i)/hnewOVCco60M2Sum->GetBinWidth(i));
    hAdapOVCk40M2Sum->SetBinContent(i, hnewOVCk40M2Sum->GetBinContent(i)/hnewOVCk40M2Sum->GetBinWidth(i));
    hAdapOVCth232M2Sum->SetBinContent(i, hnewOVCth232M2Sum->GetBinContent(i)/hnewOVCth232M2Sum->GetBinWidth(i));
    hAdapOVCu238M2Sum->SetBinContent(i, hnewOVCu238M2Sum->GetBinContent(i)/hnewOVCu238M2Sum->GetBinWidth(i));

    hAdapExtPbbi210M2Sum->SetBinContent(i, hnewExtPbbi210M2Sum->GetBinContent(i)/hnewExtPbbi210M2Sum->GetBinWidth(i));

  }
*/

  SetParEfficiency();
}

void TBackgroundModel::Initialize2()
{
  cout << "Loading PDF Histograms from file" << endl;  
  cout << "Directory " << Form("%s/MCProduction", dMCDir.c_str()) << endl;

  fBulkSmeared = new TFile(Form("%s/MCProduction_ReducedBulk.root", dMCDir.c_str()));
  fSurfaceSmeared = new TFile(Form("%s/MCProduction_ReducedSurface.root", dMCDir.c_str()));

  hTeO20nuM1     = (TH1D*)fBulkSmeared->Get("hTeO20nuM1");
  hTeO22nuM1     = (TH1D*)fBulkSmeared->Get("hTeO22nuM1");
  hTeO2co60M1    = (TH1D*)fBulkSmeared->Get("hTeO2co60M1");
  hTeO2k40M1     = (TH1D*)fBulkSmeared->Get("hTeO2k40M1");
  hTeO2po210M1   = (TH1D*)fBulkSmeared->Get("hTeO2po210M1");
  hTeO2th232onlyM1 = (TH1D*)fBulkSmeared->Get("hTeO2th232onlyM1");
  hTeO2ra228pb208M1 = (TH1D*)fBulkSmeared->Get("hTeO2ra228pb208M1");
  hTeO2th230onlyM1 = (TH1D*)fBulkSmeared->Get("hTeO2th230onlyM1");

  hTeO2Sxpb210M1_001  = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxpb210M1_001");
  hTeO2Sxpb210M1_1    = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxpb210M1_1");
  hTeO2Sxth232onlyM1_001 = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxth232onlyM1_001");
  hTeO2Sxra228pb208M1_001 = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxra228pb208M1_001");
  hTeO2Sxu238th230M1_001 = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxu238th230M1_001");
  hTeO2Sxth230onlyM1_001 = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxth230onlyM1_001");
  hTeO2Sxra226pb210M1_001 = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxra226pb210M1_001");
  hTeO2Sxth232M1_001 = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxth232M1_001");

  hCuBoxSxpb210M1_01  = (TH1D*)fSurfaceSmeared->Get("hCuBoxSxpb210M1_01");
  hCuBoxSxpb210M1_10  = (TH1D*)fSurfaceSmeared->Get("hCuBoxSxpb210M1_10");
  hCuBoxSxth232M1_10  = (TH1D*)fSurfaceSmeared->Get("hCuBoxSxth232M1_10");
  hCuBoxSxu238M1_10   = (TH1D*)fSurfaceSmeared->Get("hCuBoxSxu238M1_10");

  hCuFrameco60M1    = (TH1D*)fBulkSmeared->Get("hCuFrameco60M1");
  hCuFramek40M1     = (TH1D*)fBulkSmeared->Get("hCuFramek40M1");
  hCuFrameth232M1   = (TH1D*)fBulkSmeared->Get("hCuFrameth232M1");
  hCuFrameu238M1    = (TH1D*)fBulkSmeared->Get("hCuFrameu238M1");

  h600mKco60M1    = (TH1D*)fBulkSmeared->Get("h600mKco60M1");
  h600mKk40M1     = (TH1D*)fBulkSmeared->Get("h600mKk40M1");
  h600mKth232M1   = (TH1D*)fBulkSmeared->Get("h600mKth232M1");
  h600mKu238M1    = (TH1D*)fBulkSmeared->Get("h600mKu238M1");

  hPbRomco60M1    = (TH1D*)fBulkSmeared->Get("hPbRomco60M1");
  hPbRomk40M1     = (TH1D*)fBulkSmeared->Get("hPbRomk40M1");
  hPbRomth232M1   = (TH1D*)fBulkSmeared->Get("hPbRomth232M1");
  hPbRomu238M1    = (TH1D*)fBulkSmeared->Get("hPbRomu238M1");

  hOVCco60M1    = (TH1D*)fBulkSmeared->Get("hOVCco60M1");
  hOVCk40M1     = (TH1D*)fBulkSmeared->Get("hOVCk40M1");
  hOVCth232M1   = (TH1D*)fBulkSmeared->Get("hOVCth232M1");
  hOVCu238M1    = (TH1D*)fBulkSmeared->Get("hOVCu238M1");
  hExtPbbi210M1 = (TH1D*)fBulkSmeared->Get("hExtPbbi210M1");

  hTeO20nuM2     = (TH1D*)fBulkSmeared->Get("hTeO20nuM2");
  hTeO22nuM2     = (TH1D*)fBulkSmeared->Get("hTeO22nuM2");
  hTeO2co60M2    = (TH1D*)fBulkSmeared->Get("hTeO2co60M2");
  hTeO2k40M2     = (TH1D*)fBulkSmeared->Get("hTeO2k40M2");
  hTeO2po210M2   = (TH1D*)fBulkSmeared->Get("hTeO2po210M2");
  hTeO2th232onlyM2 = (TH1D*)fBulkSmeared->Get("hTeO2th232onlyM2");
  hTeO2ra228pb208M2 = (TH1D*)fBulkSmeared->Get("hTeO2ra228pb208M2");
  hTeO2th230onlyM2 = (TH1D*)fBulkSmeared->Get("hTeO2th230onlyM2");

  hTeO2Sxpb210M2_001  = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxpb210M2_001");
  hTeO2Sxpb210M2_1    = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxpb210M2_1");
  hTeO2Sxth232onlyM2_001 = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxth232onlyM2_001");
  hTeO2Sxra228pb208M2_001 = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxra228pb208M2_001");
  hTeO2Sxu238th230M2_001 = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxu238th230M2_001");
  hTeO2Sxth230onlyM2_001 = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxth230onlyM2_001");
  hTeO2Sxra226pb210M2_001 = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxra226pb210M2_001");
  hTeO2Sxth232M2_001 = (TH1D*)fSurfaceSmeared->Get("hTeO2Sxth232M2_001");

  hCuBoxSxpb210M2_01  = (TH1D*)fSurfaceSmeared->Get("hCuBoxSxpb210M2_01");
  hCuBoxSxpb210M2_10  = (TH1D*)fSurfaceSmeared->Get("hCuBoxSxpb210M2_10");
  hCuBoxSxth232M2_10  = (TH1D*)fSurfaceSmeared->Get("hCuBoxSxth232M2_10");
  hCuBoxSxu238M2_10   = (TH1D*)fSurfaceSmeared->Get("hCuBoxSxu238M2_10");

  hCuFrameco60M2    = (TH1D*)fBulkSmeared->Get("hCuFrameco60M2");
  hCuFramek40M2     = (TH1D*)fBulkSmeared->Get("hCuFramek40M2");
  hCuFrameth232M2   = (TH1D*)fBulkSmeared->Get("hCuFrameth232M2");
  hCuFrameu238M2    = (TH1D*)fBulkSmeared->Get("hCuFrameu238M2");

  h600mKco60M2    = (TH1D*)fBulkSmeared->Get("h600mKco60M2");
  h600mKk40M2     = (TH1D*)fBulkSmeared->Get("h600mKk40M2");
  h600mKth232M2   = (TH1D*)fBulkSmeared->Get("h600mKth232M2");
  h600mKu238M2    = (TH1D*)fBulkSmeared->Get("h600mKu238M2");

  hPbRomco60M2    = (TH1D*)fBulkSmeared->Get("hPbRomco60M2");
  hPbRomk40M2     = (TH1D*)fBulkSmeared->Get("hPbRomk40M2");
  hPbRomth232M2   = (TH1D*)fBulkSmeared->Get("hPbRomth232M2");
  hPbRomu238M2    = (TH1D*)fBulkSmeared->Get("hPbRomu238M2");

  hOVCco60M2    = (TH1D*)fBulkSmeared->Get("hOVCco60M2");
  hOVCk40M2     = (TH1D*)fBulkSmeared->Get("hOVCk40M2");
  hOVCth232M2   = (TH1D*)fBulkSmeared->Get("hOVCth232M2");
  hOVCu238M2    = (TH1D*)fBulkSmeared->Get("hOVCu238M2");

  hExtPbbi210M2 = (TH1D*)fBulkSmeared->Get("hExtPbbi210M2");


  hnewTeO20nuM1 = hTeO20nuM1->Rebin(dAdaptiveBinsM1, "hnewTeO20nuM1", dAdaptiveArrayM1);
  hnewTeO22nuM1 = hTeO22nuM1->Rebin(dAdaptiveBinsM1, "hnewTeO22nuM1", dAdaptiveArrayM1);
  hnewTeO2co60M1 = hTeO2co60M1->Rebin(dAdaptiveBinsM1, "hnewTeO2co60M1", dAdaptiveArrayM1);
  hnewTeO2k40M1 = hTeO2k40M1->Rebin(dAdaptiveBinsM1, "hnewTeO2k40M1", dAdaptiveArrayM1);
  hnewTeO2po210M1 = hTeO2po210M1->Rebin(dAdaptiveBinsM1, "hnewTeO2po210M1", dAdaptiveArrayM1);
  hnewTeO2th232onlyM1 = hTeO2th232onlyM1->Rebin(dAdaptiveBinsM1, "hnewTeO2th232onlyM1", dAdaptiveArrayM1);
  hnewTeO2ra228pb208M1 = hTeO2ra228pb208M1->Rebin(dAdaptiveBinsM1, "hnewTeO2ra228pb208M1", dAdaptiveArrayM1);
  hnewTeO2th230onlyM1 = hTeO2th230onlyM1->Rebin(dAdaptiveBinsM1, "hnewTeO2th230onlyM1", dAdaptiveArrayM1);

  hnewTeO2Sxpb210M1_001 = hTeO2Sxpb210M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpb210M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxpb210M1_1 = hTeO2Sxpb210M1_1->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxpb210M1_1", dAdaptiveArrayM1);
  hnewTeO2Sxth232onlyM1_001 = hTeO2Sxth232onlyM1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232onlyM1_001", dAdaptiveArrayM1);
  hnewTeO2Sxra228pb208M1_001 = hTeO2Sxra228pb208M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxra228pb208M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxu238th230M1_001 = hTeO2Sxu238th230M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxu238th230M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxth230onlyM1_001 = hTeO2Sxth230onlyM1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth230onlyM1_001", dAdaptiveArrayM1);
  hnewTeO2Sxra226pb210M1_001 = hTeO2Sxra226pb210M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxra226pb210M1_001", dAdaptiveArrayM1);
  hnewTeO2Sxth232M1_001 = hTeO2Sxth232M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Sxth232M1_001", dAdaptiveArrayM1);

  hnewCuBoxSxpb210M1_01 = hCuBoxSxpb210M1_01->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxpb210M1_01", dAdaptiveArrayM1);
  hnewCuBoxSxpb210M1_10 = hCuBoxSxpb210M1_10->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxpb210M1_10", dAdaptiveArrayM1);
  hnewCuBoxSxth232M1_10 = hCuBoxSxth232M1_10->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxth232M1_10", dAdaptiveArrayM1);
  hnewCuBoxSxu238M1_10 = hCuBoxSxu238M1_10->Rebin(dAdaptiveBinsM1, "hnewCuBoxSxu238M1_10", dAdaptiveArrayM1);

  hnewCuFrameco60M1 = hCuFrameco60M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameco60M1", dAdaptiveArrayM1);
  hnewCuFramek40M1 = hCuFramek40M1->Rebin(dAdaptiveBinsM1, "hnewCuFramek40M1", dAdaptiveArrayM1);
  hnewCuFrameth232M1 = hCuFrameth232M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameth232M1", dAdaptiveArrayM1);
  hnewCuFrameu238M1 = hCuFrameu238M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameu238M1", dAdaptiveArrayM1);

  hnew600mKco60M1 = h600mKco60M1->Rebin(dAdaptiveBinsM1, "hnew600mKco60M1", dAdaptiveArrayM1);
  hnew600mKk40M1 = h600mKk40M1->Rebin(dAdaptiveBinsM1, "hnew600mKk40M1", dAdaptiveArrayM1);
  hnew600mKth232M1 = h600mKth232M1->Rebin(dAdaptiveBinsM1, "hnew600mKth232M1", dAdaptiveArrayM1);
  hnew600mKu238M1 = h600mKu238M1->Rebin(dAdaptiveBinsM1, "hnew600mKu238M1", dAdaptiveArrayM1);

  hnewPbRomco60M1 = hPbRomco60M1->Rebin(dAdaptiveBinsM1, "hnewPbRomco60M1", dAdaptiveArrayM1);
  hnewPbRomk40M1 = hPbRomk40M1->Rebin(dAdaptiveBinsM1, "hnewPbRomk40M1", dAdaptiveArrayM1);
  hnewPbRomth232M1 = hPbRomth232M1->Rebin(dAdaptiveBinsM1, "hnewPbRomth232M1", dAdaptiveArrayM1);
  hnewPbRomu238M1 = hPbRomu238M1->Rebin(dAdaptiveBinsM1, "hnewPbRomu238M1", dAdaptiveArrayM1);

  hnewOVCco60M1 = hOVCco60M1->Rebin(dAdaptiveBinsM1, "hnewOVCco60M1", dAdaptiveArrayM1);
  hnewOVCk40M1 = hOVCk40M1->Rebin(dAdaptiveBinsM1, "hnewOVCk40M1", dAdaptiveArrayM1);
  hnewOVCth232M1 = hOVCth232M1->Rebin(dAdaptiveBinsM1, "hnewOVCth232M1", dAdaptiveArrayM1);
  hnewOVCu238M1 = hOVCu238M1->Rebin(dAdaptiveBinsM1, "hnewOVCu238M1", dAdaptiveArrayM1);

  hnewExtPbbi210M1 = hExtPbbi210M1->Rebin(dAdaptiveBinsM1, "hnewExtPbbi210M1", dAdaptiveArrayM1);


  hnewTeO20nuM2 = hTeO20nuM2->Rebin(dAdaptiveBinsM2, "hnewTeO20nuM2", dAdaptiveArrayM2);
  hnewTeO22nuM2 = hTeO22nuM2->Rebin(dAdaptiveBinsM2, "hnewTeO22nuM2", dAdaptiveArrayM2);
  hnewTeO2co60M2 = hTeO2co60M2->Rebin(dAdaptiveBinsM2, "hnewTeO2co60M2", dAdaptiveArrayM2);
  hnewTeO2k40M2 = hTeO2k40M2->Rebin(dAdaptiveBinsM2, "hnewTeO2k40M2", dAdaptiveArrayM2);
  hnewTeO2po210M2 = hTeO2po210M2->Rebin(dAdaptiveBinsM2, "hnewTeO2po210M2", dAdaptiveArrayM2);
  hnewTeO2th232onlyM2 = hTeO2th232onlyM2->Rebin(dAdaptiveBinsM2, "hnewTeO2th232onlyM2", dAdaptiveArrayM2);
  hnewTeO2ra228pb208M2 = hTeO2ra228pb208M2->Rebin(dAdaptiveBinsM2, "hnewTeO2ra228pb208M2", dAdaptiveArrayM2);
  hnewTeO2th230onlyM2 = hTeO2th230onlyM2->Rebin(dAdaptiveBinsM2, "hnewTeO2th230onlyM2", dAdaptiveArrayM2);

  hnewTeO2Sxpb210M2_001 = hTeO2Sxpb210M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpb210M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxpb210M2_1 = hTeO2Sxpb210M2_1->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxpb210M2_1", dAdaptiveArrayM2);
  hnewTeO2Sxth232onlyM2_001 = hTeO2Sxth232onlyM2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232onlyM2_001", dAdaptiveArrayM2);
  hnewTeO2Sxra228pb208M2_001 = hTeO2Sxra228pb208M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxra228pb208M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxu238th230M2_001 = hTeO2Sxu238th230M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxu238th230M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxth230onlyM2_001 = hTeO2Sxth230onlyM2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth230onlyM2_001", dAdaptiveArrayM2);
  hnewTeO2Sxra226pb210M2_001 = hTeO2Sxra226pb210M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxra226pb210M2_001", dAdaptiveArrayM2);
  hnewTeO2Sxth232M2_001 = hTeO2Sxth232M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Sxth232M2_001", dAdaptiveArrayM2);

  hnewCuBoxSxpb210M2_01 = hCuBoxSxpb210M2_01->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxpb210M2_01", dAdaptiveArrayM2);
  hnewCuBoxSxpb210M2_10 = hCuBoxSxpb210M2_10->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxpb210M2_10", dAdaptiveArrayM2);
  hnewCuBoxSxth232M2_10 = hCuBoxSxth232M2_10->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxth232M2_10", dAdaptiveArrayM2);
  hnewCuBoxSxu238M2_10 = hCuBoxSxu238M2_10->Rebin(dAdaptiveBinsM2, "hnewCuBoxSxu238M2_10", dAdaptiveArrayM2);

  hnewCuFrameco60M2 = hCuFrameco60M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameco60M2", dAdaptiveArrayM2);
  hnewCuFramek40M2 = hCuFramek40M2->Rebin(dAdaptiveBinsM2, "hnewCuFramek40M2", dAdaptiveArrayM2);
  hnewCuFrameth232M2 = hCuFrameth232M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameth232M2", dAdaptiveArrayM2);
  hnewCuFrameu238M2 = hCuFrameu238M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameu238M2", dAdaptiveArrayM2);

  hnew600mKco60M2 =h600mKco60M2->Rebin(dAdaptiveBinsM2, "hnew600mKco60M2", dAdaptiveArrayM2);
  hnew600mKk40M2 =h600mKk40M2->Rebin(dAdaptiveBinsM2, "hnew600mKk40M2", dAdaptiveArrayM2);
  hnew600mKth232M2 =h600mKth232M2->Rebin(dAdaptiveBinsM2, "hnew600mKth232M2", dAdaptiveArrayM2);
  hnew600mKu238M2 =h600mKu238M2->Rebin(dAdaptiveBinsM2, "hnew600mKu238M2", dAdaptiveArrayM2);

  hnewPbRomco60M2 = hPbRomco60M2->Rebin(dAdaptiveBinsM2, "hnewPbRomco60M2", dAdaptiveArrayM2);
  hnewPbRomk40M2 = hPbRomk40M2->Rebin(dAdaptiveBinsM2, "hnewPbRomk40M2", dAdaptiveArrayM2);
  hnewPbRomth232M2 = hPbRomth232M2->Rebin(dAdaptiveBinsM2, "hnewPbRomth232M2", dAdaptiveArrayM2);
  hnewPbRomu238M2 = hPbRomu238M2->Rebin(dAdaptiveBinsM2, "hnewPbRomu238M2", dAdaptiveArrayM2);

  hnewOVCco60M2 = hOVCco60M2->Rebin(dAdaptiveBinsM2, "hnewOVCco60M2", dAdaptiveArrayM2);
  hnewOVCk40M2 = hOVCk40M2->Rebin(dAdaptiveBinsM2, "hnewOVCk40M2", dAdaptiveArrayM2);
  hnewOVCth232M2 = hOVCth232M2->Rebin(dAdaptiveBinsM2, "hnewOVCth232M2", dAdaptiveArrayM2);
  hnewOVCu238M2 = hOVCu238M2->Rebin(dAdaptiveBinsM2, "hnewOVCu238M2", dAdaptiveArrayM2);

  hnewExtPbbi210M2 = hExtPbbi210M2->Rebin(dAdaptiveBinsM2, "hnewExtPbbi210M2", dAdaptiveArrayM2);


  for(int i = 1; i <= dAdaptiveBinsM1; i++)
  {
    hAdapTeO20nuM1->SetBinContent(i, dBinSize * hnewTeO20nuM1->GetBinContent(i)/hnewTeO20nuM1->GetBinWidth(i));
    hAdapTeO22nuM1->SetBinContent(i, dBinSize * hnewTeO22nuM1->GetBinContent(i)/hnewTeO22nuM1->GetBinWidth(i));
    hAdapTeO2co60M1->SetBinContent(i, dBinSize * hnewTeO2co60M1->GetBinContent(i)/hnewTeO2co60M1->GetBinWidth(i));
    hAdapTeO2k40M1->SetBinContent(i, dBinSize * hnewTeO2k40M1->GetBinContent(i)/hnewTeO2k40M1->GetBinWidth(i));
    hAdapTeO2po210M1->SetBinContent(i, dBinSize * hnewTeO2po210M1->GetBinContent(i)/hnewTeO2po210M1->GetBinWidth(i));
    hAdapTeO2th232onlyM1->SetBinContent(i, dBinSize * hnewTeO2th232onlyM1->GetBinContent(i)/hnewTeO2th232onlyM1->GetBinWidth(i));
    hAdapTeO2ra228pb208M1->SetBinContent(i, dBinSize * hnewTeO2ra228pb208M1->GetBinContent(i)/hnewTeO2ra228pb208M1->GetBinWidth(i));
    hAdapTeO2th230onlyM1->SetBinContent(i, dBinSize * hnewTeO2th230onlyM1->GetBinContent(i)/hnewTeO2th230onlyM1->GetBinWidth(i));

    hAdapTeO2Sxth232onlyM1_001->SetBinContent(i, dBinSize * hnewTeO2Sxth232onlyM1_001->GetBinContent(i)/hnewTeO2Sxth232onlyM1_001->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M1_001->SetBinContent(i, dBinSize * hnewTeO2Sxra228pb208M1_001->GetBinContent(i)/hnewTeO2Sxra228pb208M1_001->GetBinWidth(i));
    hAdapTeO2Sxu238th230M1_001->SetBinContent(i, dBinSize * hnewTeO2Sxu238th230M1_001->GetBinContent(i)/hnewTeO2Sxu238th230M1_001->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM1_001->SetBinContent(i, dBinSize * hnewTeO2Sxth230onlyM1_001->GetBinContent(i)/hnewTeO2Sxth230onlyM1_001->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M1_001->SetBinContent(i, dBinSize * hnewTeO2Sxra226pb210M1_001->GetBinContent(i)/hnewTeO2Sxra226pb210M1_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M1_001->SetBinContent(i, dBinSize * hnewTeO2Sxpb210M1_001->GetBinContent(i)/hnewTeO2Sxpb210M1_001->GetBinWidth(i));
    hAdapTeO2Sxth232M1_001->SetBinContent(i, dBinSize * hnewTeO2Sxth232M1_001->GetBinContent(i)/hnewTeO2Sxth232M1_001->GetBinWidth(i));

    hAdapTeO2Sxpb210M1_1->SetBinContent(i, dBinSize * hnewTeO2Sxpb210M1_1->GetBinContent(i)/hnewTeO2Sxpb210M1_1->GetBinWidth(i));
    hAdapCuBoxSxpb210M1_01->SetBinContent(i, dBinSize * hnewCuBoxSxpb210M1_01->GetBinContent(i)/hnewCuBoxSxpb210M1_01->GetBinWidth(i));
    hAdapCuBoxSxpb210M1_10->SetBinContent(i, dBinSize * hnewCuBoxSxpb210M1_10->GetBinContent(i)/hnewCuBoxSxpb210M1_10->GetBinWidth(i));
    hAdapCuBoxSxth232M1_10->SetBinContent(i, dBinSize * hnewCuBoxSxth232M1_10->GetBinContent(i)/hnewCuBoxSxth232M1_10->GetBinWidth(i));
    hAdapCuBoxSxu238M1_10->SetBinContent(i, dBinSize * hnewCuBoxSxu238M1_10->GetBinContent(i)/hnewCuBoxSxu238M1_10->GetBinWidth(i));


    hAdapCuFrameco60M1->SetBinContent(i, dBinSize * hnewCuFrameco60M1->GetBinContent(i)/hnewCuFrameco60M1->GetBinWidth(i));
    hAdapCuFramek40M1->SetBinContent(i, dBinSize * hnewCuFramek40M1->GetBinContent(i)/hnewCuFramek40M1->GetBinWidth(i));
    hAdapCuFrameth232M1->SetBinContent(i, dBinSize * hnewCuFrameth232M1->GetBinContent(i)/hnewCuFrameth232M1->GetBinWidth(i));
    hAdapCuFrameu238M1->SetBinContent(i, dBinSize * hnewCuFrameu238M1->GetBinContent(i)/hnewCuFrameu238M1->GetBinWidth(i));

    hAdap600mKco60M1->SetBinContent(i, dBinSize * hnew600mKco60M1->GetBinContent(i)/hnew600mKco60M1->GetBinWidth(i));
    hAdap600mKk40M1->SetBinContent(i, dBinSize * hnew600mKk40M1->GetBinContent(i)/hnew600mKk40M1->GetBinWidth(i));
    hAdap600mKth232M1->SetBinContent(i, dBinSize * hnew600mKth232M1->GetBinContent(i)/hnew600mKth232M1->GetBinWidth(i));
    hAdap600mKu238M1->SetBinContent(i, dBinSize * hnew600mKu238M1->GetBinContent(i)/hnew600mKu238M1->GetBinWidth(i));

    hAdapPbRomco60M1->SetBinContent(i, dBinSize * hnewPbRomco60M1->GetBinContent(i)/hnewPbRomco60M1->GetBinWidth(i));
    hAdapPbRomk40M1->SetBinContent(i, dBinSize * hnewPbRomk40M1->GetBinContent(i)/hnewPbRomk40M1->GetBinWidth(i));
    hAdapPbRomth232M1->SetBinContent(i, dBinSize * hnewPbRomth232M1->GetBinContent(i)/hnewPbRomth232M1->GetBinWidth(i));
    hAdapPbRomu238M1->SetBinContent(i, dBinSize * hnewPbRomu238M1->GetBinContent(i)/hnewPbRomu238M1->GetBinWidth(i));

    hAdapOVCco60M1->SetBinContent(i, dBinSize * hnewOVCco60M1->GetBinContent(i)/hnewOVCco60M1->GetBinWidth(i));
    hAdapOVCk40M1->SetBinContent(i, dBinSize * hnewOVCk40M1->GetBinContent(i)/hnewOVCk40M1->GetBinWidth(i));
    hAdapOVCth232M1->SetBinContent(i, dBinSize * hnewOVCth232M1->GetBinContent(i)/hnewOVCth232M1->GetBinWidth(i));
    hAdapOVCu238M1->SetBinContent(i, dBinSize * hnewOVCu238M1->GetBinContent(i)/hnewOVCu238M1->GetBinWidth(i));

    hAdapExtPbbi210M1->SetBinContent(i, dBinSize * hnewExtPbbi210M1->GetBinContent(i)/hnewExtPbbi210M1->GetBinWidth(i));


  }

  for(int i = 1; i <= dAdaptiveBinsM2; i++)
  {
    hAdapTeO20nuM2->SetBinContent(i, dBinSize * hnewTeO20nuM2->GetBinContent(i)/hnewTeO20nuM2->GetBinWidth(i));
    hAdapTeO22nuM2->SetBinContent(i, dBinSize * hnewTeO22nuM2->GetBinContent(i)/hnewTeO22nuM2->GetBinWidth(i));
    hAdapTeO2co60M2->SetBinContent(i, dBinSize * hnewTeO2co60M2->GetBinContent(i)/hnewTeO2co60M2->GetBinWidth(i));
    hAdapTeO2k40M2->SetBinContent(i, dBinSize * hnewTeO2k40M2->GetBinContent(i)/hnewTeO2k40M2->GetBinWidth(i));
    hAdapTeO2po210M2->SetBinContent(i, dBinSize * hnewTeO2po210M2->GetBinContent(i)/hnewTeO2po210M2->GetBinWidth(i));
    hAdapTeO2th232onlyM2->SetBinContent(i, dBinSize * hnewTeO2th232onlyM2->GetBinContent(i)/hnewTeO2th232onlyM2->GetBinWidth(i));
    hAdapTeO2ra228pb208M2->SetBinContent(i, dBinSize * hnewTeO2ra228pb208M2->GetBinContent(i)/hnewTeO2ra228pb208M2->GetBinWidth(i));
    hAdapTeO2th230onlyM2->SetBinContent(i, dBinSize * hnewTeO2th230onlyM2->GetBinContent(i)/hnewTeO2th230onlyM2->GetBinWidth(i));

    hAdapTeO2Sxth232onlyM2_001->SetBinContent(i, dBinSize * hnewTeO2Sxth232onlyM2_001->GetBinContent(i)/hnewTeO2Sxth232onlyM2_001->GetBinWidth(i));
    hAdapTeO2Sxra228pb208M2_001->SetBinContent(i, dBinSize * hnewTeO2Sxra228pb208M2_001->GetBinContent(i)/hnewTeO2Sxra228pb208M2_001->GetBinWidth(i));
    hAdapTeO2Sxu238th230M2_001->SetBinContent(i, dBinSize * hnewTeO2Sxu238th230M2_001->GetBinContent(i)/hnewTeO2Sxu238th230M2_001->GetBinWidth(i));
    hAdapTeO2Sxth230onlyM2_001->SetBinContent(i, dBinSize * hnewTeO2Sxth230onlyM2_001->GetBinContent(i)/hnewTeO2Sxth230onlyM2_001->GetBinWidth(i));
    hAdapTeO2Sxra226pb210M2_001->SetBinContent(i, dBinSize * hnewTeO2Sxra226pb210M2_001->GetBinContent(i)/hnewTeO2Sxra226pb210M2_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_001->SetBinContent(i, dBinSize * hnewTeO2Sxpb210M2_001->GetBinContent(i)/hnewTeO2Sxpb210M2_001->GetBinWidth(i));
    hAdapTeO2Sxth232M2_001->SetBinContent(i, dBinSize * hnewTeO2Sxth232M2_001->GetBinContent(i)/hnewTeO2Sxth232M2_001->GetBinWidth(i));

    hAdapTeO2Sxpb210M2_1->SetBinContent(i, dBinSize * hnewTeO2Sxpb210M2_1->GetBinContent(i)/hnewTeO2Sxpb210M2_1->GetBinWidth(i));
    hAdapCuBoxSxpb210M2_01->SetBinContent(i, dBinSize * hnewCuBoxSxpb210M2_01->GetBinContent(i)/hnewCuBoxSxpb210M2_01->GetBinWidth(i));
    hAdapCuBoxSxpb210M2_10->SetBinContent(i, dBinSize * hnewCuBoxSxpb210M2_10->GetBinContent(i)/hnewCuBoxSxpb210M2_10->GetBinWidth(i));
    hAdapCuBoxSxth232M2_10->SetBinContent(i, dBinSize * hnewCuBoxSxth232M2_10->GetBinContent(i)/hnewCuBoxSxth232M2_10->GetBinWidth(i));
    hAdapCuBoxSxu238M2_10->SetBinContent(i, dBinSize * hnewCuBoxSxu238M2_10->GetBinContent(i)/hnewCuBoxSxu238M2_10->GetBinWidth(i));

    hAdapCuFrameco60M2->SetBinContent(i, dBinSize * hnewCuFrameco60M2->GetBinContent(i)/hnewCuFrameco60M2->GetBinWidth(i));
    hAdapCuFramek40M2->SetBinContent(i, dBinSize * hnewCuFramek40M2->GetBinContent(i)/hnewCuFramek40M2->GetBinWidth(i));
    hAdapCuFrameth232M2->SetBinContent(i, dBinSize * hnewCuFrameth232M2->GetBinContent(i)/hnewCuFrameth232M2->GetBinWidth(i));
    hAdapCuFrameu238M2->SetBinContent(i, dBinSize * hnewCuFrameu238M2->GetBinContent(i)/hnewCuFrameu238M2->GetBinWidth(i));

    hAdap600mKco60M2->SetBinContent(i, dBinSize * hnew600mKco60M2->GetBinContent(i)/hnew600mKco60M2->GetBinWidth(i));
    hAdap600mKk40M2->SetBinContent(i, dBinSize * hnew600mKk40M2->GetBinContent(i)/hnew600mKk40M2->GetBinWidth(i));
    hAdap600mKth232M2->SetBinContent(i, dBinSize * hnew600mKth232M2->GetBinContent(i)/hnew600mKth232M2->GetBinWidth(i));
    hAdap600mKu238M2->SetBinContent(i, dBinSize * hnew600mKu238M2->GetBinContent(i)/hnew600mKu238M2->GetBinWidth(i));

    hAdapPbRomco60M2->SetBinContent(i, dBinSize * hnewPbRomco60M2->GetBinContent(i)/hnewPbRomco60M2->GetBinWidth(i));
    hAdapPbRomk40M2->SetBinContent(i, dBinSize * hnewPbRomk40M2->GetBinContent(i)/hnewPbRomk40M2->GetBinWidth(i));
    hAdapPbRomth232M2->SetBinContent(i, dBinSize * hnewPbRomth232M2->GetBinContent(i)/hnewPbRomth232M2->GetBinWidth(i));
    hAdapPbRomu238M2->SetBinContent(i, dBinSize * hnewPbRomu238M2->GetBinContent(i)/hnewPbRomu238M2->GetBinWidth(i));

    hAdapOVCco60M2->SetBinContent(i, dBinSize * hnewOVCco60M2->GetBinContent(i)/hnewOVCco60M2->GetBinWidth(i));
    hAdapOVCk40M2->SetBinContent(i, dBinSize * hnewOVCk40M2->GetBinContent(i)/hnewOVCk40M2->GetBinWidth(i));
    hAdapOVCth232M2->SetBinContent(i, dBinSize * hnewOVCth232M2->GetBinContent(i)/hnewOVCth232M2->GetBinWidth(i));
    hAdapOVCu238M2->SetBinContent(i, dBinSize * hnewOVCu238M2->GetBinContent(i)/hnewOVCu238M2->GetBinWidth(i));

    hAdapExtPbbi210M2->SetBinContent(i, dBinSize * hnewExtPbbi210M2->GetBinContent(i)/hnewExtPbbi210M2->GetBinWidth(i));


  }



  SetParEfficiency();
}

// Loads the background data
void TBackgroundModel::LoadData()
{
  switch(dDataSet)
  { 
  case 1:
    qtree->Add(Form("%s/ReducedB-ds2049.root", dDataDir.c_str()));   
    qtree->Add(Form("%s/ReducedB-ds2061.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/ReducedB-ds2064.root", dDataDir.c_str()));   
    qtree->Add(Form("%s/ReducedB-ds2067.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/ReducedB-ds2070.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/ReducedB-ds2073.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/ReducedB-ds2076.root", dDataDir.c_str())); 
    dLivetime = 6042498; // DR 1 
    cout << "Using Data Release 1" << endl;
  break;

  case 2:
    qtree->Add(Form("%s/ReducedB-ds2079.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/ReducedB-ds2085.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/ReducedB-ds2088.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/ReducedB-ds2091.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/ReducedB-ds2097.root", dDataDir.c_str()));
    qtree->Add(Form("%s/ReducedB-ds2100.root", dDataDir.c_str())); 
    dLivetime = 9387524; // DR 2
    cout << "Using Data Release 2" << endl;
  break;

  case 3:
    qtree->Add(Form("%s/ReducedB-ds2103.root", dDataDir.c_str()));
    qtree->Add(Form("%s/ReducedB-ds2109.root", dDataDir.c_str()));
    qtree->Add(Form("%s/ReducedB-ds2118.root", dDataDir.c_str()));
    qtree->Add(Form("%s/ReducedB-ds2124.root", dDataDir.c_str()));
    dLivetime = 7647908; // DR 3
    cout << "Using Data Release 3" << endl;
  break;

  default:   
    qtree->Add(Form("%s/ReducedB-ds*.root", dDataDir.c_str())); 
    dLivetime = 23077930; // seconds of livetime (DR1 to DR3)
    cout << "Using Total Dataset" << endl;

  }

  qtree->Project("fDataHistoTot", "Energy", base_cut);
  qtree->Project("fDataHistoM1",  "Energy", base_cut && "Multiplicity_Sync == 1");
  qtree->Project("fDataHistoM2",  "Energy", base_cut && "Multiplicity_Sync == 2");
  qtree->Project("fDataHistoM2Sum",  "TotalEnergy_Sync", base_cut && "Multiplicity_Sync == 2");

  dLivetimeYr = dLivetime*dSecToYears;  
}

// Prints parameters, needs update 11-06-2014
// void TBackgroundModel::PrintParameters()
// {
//   for(int i = 0; i < TBackgroundModel::dNParam; i++)
//   {
//     cout << minuit->fCpnam[i] << ": " << fParameters[i] << " +/- " << fParError[i] << endl;
//   }
// }

// void TBackgroundModel::PrintParActivity()
// { 
//   // This only gets the number of counts corrected for detector efficiency
//   for(int i = 0; i < TBackgroundModel::dNParam; i++)
//   {
//   cout << minuit->fCpnam[i] << " activity: " << fParActivity[i] << endl;
//   }
// 
// }

// Resets all parameters to 0
void TBackgroundModel::ResetParameters()
{
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
  fParEfficiencyM1[0] = 0.8897514;
  fParEfficiencyM1[1] = 0.9292826;
  fParEfficiencyM1[2] = 0.4424162;
  fParEfficiencyM1[3] = 0.8584163;
  fParEfficiencyM1[4] = 0.9418009;
  fParEfficiencyM1[5] = 0.941975;
  fParEfficiencyM1[6] = 0.4975102;
  fParEfficiencyM1[7] = 0.9419295;
  fParEfficiencyM1[8] = 0.7376795;
  fParEfficiencyM1[9] = 0.3802244;
  fParEfficiencyM1[10] = 0.7457738;
  fParEfficiencyM1[11] = 0.7379278;
  fParEfficiencyM1[12] = 0.3802244;
  fParEfficiencyM1[13] = 0.4457218;
  fParEfficiencyM1[14] = 0.0510853;
  fParEfficiencyM1[15] = 0.0354736;
  fParEfficiencyM1[16] = 0.2807244;
  fParEfficiencyM1[17] = 0.0369882;
  fParEfficiencyM1[18] = 0.1013459;
  fParEfficiencyM1[19] = 0.0881936;
  fParEfficiencyM1[20] = 0.046961;
  fParEfficiencyM1[21] = 0.0725138;
  fParEfficiencyM1[22] = 0.006387241;
  fParEfficiencyM1[23] = 0.012149995;
  fParEfficiencyM1[24] = 0.05125624;
  fParEfficiencyM1[25] = 0.00937862;
  fParEfficiencyM1[26] = 0.000191901;
  fParEfficiencyM1[27] = 0.07314378;
  fParEfficiencyM1[28] = 0.00484976;
  fParEfficiencyM1[29] = 0.006349752;
  fParEfficiencyM1[30] = 0.004160072;
  fParEfficiencyM1[31] = 0.00715452;
  fParEfficiencyM1[32] = 0.00043995;
  fParEfficiencyM1[33] = 0.000465822;
  fParEfficiencyM1[34] = 0.000293556;
  fParEfficiencyM1[35] = 0.4958653;
  fParEfficiencyM1[36] = 0.0742276;
}

// Creates/updates the background model
void TBackgroundModel::UpdateModelAdaptive()
{
  if(fModelTotAdapM1 == NULL || fModelTotAdapM2 == NULL) 
  {
    cout << "Model Histogram Not Created" << endl;
    return;
  }

  // Reset all bins in model histograms in every loop
  fModelTotAdapM1->Reset();
  fModelTotAdapM2->Reset();


  dNumCalls++;
  if(dNumCalls%1000==0)
  {
    cout << "Call #: "<< dNumCalls << endl;
  }

/////// Create model
/////// M1
  fModelTotAdapM1->Add( hAdapTeO20nuM1,              dDataIntegralM1*fParameters[0]);
  fModelTotAdapM1->Add( hAdapTeO22nuM1,              dDataIntegralM1*fParameters[1]);
  fModelTotAdapM1->Add( hAdapTeO2co60M1,             dDataIntegralM1*fParameters[2]);
  fModelTotAdapM1->Add( hAdapTeO2k40M1,              dDataIntegralM1*fParameters[3]);
  fModelTotAdapM1->Add( hAdapTeO2po210M1,            dDataIntegralM1*fParameters[4]);
  fModelTotAdapM1->Add( hAdapTeO2th232onlyM1,        dDataIntegralM1*fParameters[5]);
  fModelTotAdapM1->Add( hAdapTeO2th230onlyM1,        dDataIntegralM1*fParameters[6]);
  fModelTotAdapM1->Add( hAdapTeO2Sxth232M1_001,      dDataIntegralM1*fParameters[7]);
  fModelTotAdapM1->Add( hAdapTeO2Sxth232onlyM1_001,  dDataIntegralM1*fParameters[8]);
  fModelTotAdapM1->Add( hAdapTeO2Sxra228pb208M1_001, dDataIntegralM1*fParameters[9]);
  fModelTotAdapM1->Add( hAdapTeO2Sxu238th230M1_001,  dDataIntegralM1*fParameters[10]);
  fModelTotAdapM1->Add( hAdapTeO2Sxth230onlyM1_001,  dDataIntegralM1*fParameters[11]);
  fModelTotAdapM1->Add( hAdapTeO2Sxra226pb210M1_001, dDataIntegralM1*fParameters[12]);
  fModelTotAdapM1->Add( hAdapTeO2Sxpb210M1_1,        dDataIntegralM1*fParameters[13]);
  fModelTotAdapM1->Add( hAdapTeO2Sxpb210M1_001,      dDataIntegralM1*fParameters[14]);

  fModelTotAdapM1->Add( hAdapCuBox_CuFrameth232M1_10,   dDataIntegralM1*fParameters[15]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFrameu238M1_10,    dDataIntegralM1*fParameters[16]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFramepb210M1_10,   dDataIntegralM1*fParameters[17]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFramepb210M1_01,   dDataIntegralM1*fParameters[18]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFramepb210M1_001,  dDataIntegralM1*fParameters[19]);

  fModelTotAdapM1->Add( hAdapPbRomth232M1,     dDataIntegralM1*fParameters[20]);
  fModelTotAdapM1->Add( hAdapPbRomu238M1,      dDataIntegralM1*fParameters[21]);
  fModelTotAdapM1->Add( hAdapPbRomco60M1,      dDataIntegralM1*fParameters[22]);
  fModelTotAdapM1->Add( hAdapPbRomk40M1,       dDataIntegralM1*fParameters[23]);

  fModelTotAdapM1->Add( hAdapOVCth232M1,     dDataIntegralM1*fParameters[24]);
  fModelTotAdapM1->Add( hAdapOVCu238M1,      dDataIntegralM1*fParameters[25]);
  fModelTotAdapM1->Add( hAdapOVCco60M1,      dDataIntegralM1*fParameters[26]);
  fModelTotAdapM1->Add( hAdapOVCk40M1,       dDataIntegralM1*fParameters[27]);
  fModelTotAdapM1->Add( hAdapExtPbbi210M1,   dDataIntegralM1*fParameters[28]);

  fModelTotAdapM1->Add( hAdapInternalth232M1,     dDataIntegralM1*fParameters[29]);
  fModelTotAdapM1->Add( hAdapInternalu238M1,      dDataIntegralM1*fParameters[30]);
  fModelTotAdapM1->Add( hAdapInternalco60M1,      dDataIntegralM1*fParameters[31]);
  fModelTotAdapM1->Add( hAdapInternalk40M1,       dDataIntegralM1*fParameters[32]);

/////// M2
  fModelTotAdapM2->Add( hAdapTeO20nuM2,              dDataIntegralM1*fParameters[0]);
  fModelTotAdapM2->Add( hAdapTeO22nuM2,              dDataIntegralM1*fParameters[1]);
  fModelTotAdapM2->Add( hAdapTeO2co60M2,             dDataIntegralM1*fParameters[2]);
  fModelTotAdapM2->Add( hAdapTeO2k40M2,              dDataIntegralM1*fParameters[3]);
  fModelTotAdapM2->Add( hAdapTeO2po210M2,            dDataIntegralM1*fParameters[4]);
  fModelTotAdapM2->Add( hAdapTeO2th232onlyM2,        dDataIntegralM1*fParameters[5]);
  fModelTotAdapM2->Add( hAdapTeO2th230onlyM2,        dDataIntegralM1*fParameters[6]);
  fModelTotAdapM2->Add( hAdapTeO2Sxth232M2_001,      dDataIntegralM1*fParameters[7]);
  fModelTotAdapM2->Add( hAdapTeO2Sxth232onlyM2_001,  dDataIntegralM1*fParameters[8]);
  fModelTotAdapM2->Add( hAdapTeO2Sxra228pb208M2_001, dDataIntegralM1*fParameters[9]);
  fModelTotAdapM2->Add( hAdapTeO2Sxu238th230M2_001,  dDataIntegralM1*fParameters[10]);
  fModelTotAdapM2->Add( hAdapTeO2Sxth230onlyM2_001,  dDataIntegralM1*fParameters[11]);
  fModelTotAdapM2->Add( hAdapTeO2Sxra226pb210M2_001, dDataIntegralM1*fParameters[12]);
  fModelTotAdapM2->Add( hAdapTeO2Sxpb210M2_1,        dDataIntegralM1*fParameters[13]);
  fModelTotAdapM2->Add( hAdapTeO2Sxpb210M2_001,      dDataIntegralM1*fParameters[14]);

  fModelTotAdapM2->Add( hAdapCuBox_CuFrameth232M2_10,   dDataIntegralM1*fParameters[15]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFrameu238M2_10,    dDataIntegralM1*fParameters[16]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFramepb210M2_10,   dDataIntegralM1*fParameters[17]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFramepb210M2_01,   dDataIntegralM1*fParameters[18]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFramepb210M2_001,  dDataIntegralM1*fParameters[19]);

  fModelTotAdapM2->Add( hAdapPbRomth232M2,     dDataIntegralM1*fParameters[20]);
  fModelTotAdapM2->Add( hAdapPbRomu238M2,      dDataIntegralM1*fParameters[21]);
  fModelTotAdapM2->Add( hAdapPbRomco60M2,      dDataIntegralM1*fParameters[22]);
  fModelTotAdapM2->Add( hAdapPbRomk40M2,       dDataIntegralM1*fParameters[23]);

  fModelTotAdapM2->Add( hAdapOVCth232M2,     dDataIntegralM1*fParameters[24]);
  fModelTotAdapM2->Add( hAdapOVCu238M2,      dDataIntegralM1*fParameters[25]);
  fModelTotAdapM2->Add( hAdapOVCco60M2,      dDataIntegralM1*fParameters[26]);
  fModelTotAdapM2->Add( hAdapOVCk40M2,       dDataIntegralM1*fParameters[27]);
  fModelTotAdapM2->Add( hAdapExtPbbi210M2,   dDataIntegralM1*fParameters[28]);

  fModelTotAdapM2->Add( hAdapInternalth232M2,     dDataIntegralM1*fParameters[29]);
  fModelTotAdapM2->Add( hAdapInternalu238M2,      dDataIntegralM1*fParameters[30]);
  fModelTotAdapM2->Add( hAdapInternalco60M2,      dDataIntegralM1*fParameters[31]);
  fModelTotAdapM2->Add( hAdapInternalk40M2,       dDataIntegralM1*fParameters[32]);

}
  
// This method sets up minuit and does the fit
bool TBackgroundModel::DoTheFitAdaptive(double f2nuValue)
{ 
  // Reset initial parameter/error values
  ResetParameters();

  gStyle->SetOptStat(0);
  gStyle->SetOptFit();

  // Reduce Minuit Output
  minuit->SetPrintLevel(0); // Print level -1 (Quiet), 0 (Normal), 1 (Verbose)
  minuit->Command("SET STRategy 2"); // Sets strategy of fit
  minuit->SetMaxIterations(10000);
  minuit->SetObjectFit(this); //see the external FCN  above

  // With Initial Values
  minuit->DefineParameter(0, "TeO2 0nu",  0., 1E-6, 0, 1.0);
  minuit->DefineParameter(1, "TeO2 2nu",  f2nuValue, 1E-6, 0, 1.0);
  minuit->DefineParameter(2, "TeO2 co60",  0., 1E-6, 0, 1.0);
  minuit->DefineParameter(3, "TeO2 k40",  0., 1E-6, 0, 1.0);
  minuit->DefineParameter(4, "TeO2 po210",  0., 1E-6, 0, 1.0);
  minuit->DefineParameter(5, "TeO2 th232 only", 0.0002930441, 1E-6, 0, 1.0);
  minuit->DefineParameter(6, "TeO2 th230 only", 0.0002989944, 1E-6, 0, 1.0);
  minuit->DefineParameter(7, "TeO2 Sx th232 0.01", 0.0016914715, 1E-6, 0, 1.0);
  minuit->DefineParameter(8, "TeO2 Sx th232 only 0.01", 0., 1E-6, 0, 1.0);
  minuit->DefineParameter(9, "TeO2 Sx ra228 to pb208 0.01", 0.0025632619, 1E-6, 0, 1.0);
  minuit->DefineParameter(10, "TeO2 Sx u238 to th230 0.01", 0.0017130443, 1E-6, 0, 1.0);
  minuit->DefineParameter(11, "TeO2 Sx th230 only 0.01", 0.0007863532, 1E-6, 0, 1.0);
  minuit->DefineParameter(12, "TeO2 Sx ra226 to pb210 0.01", 0.0030621785, 1E-6, 0, 1.0);
  minuit->DefineParameter(13, "TeO2 Sx pb210 1", 0.0054316348, 1E-6, 0, 1.0);
  minuit->DefineParameter(14, "TeO2 Sx pb210 0.01", 0.0420245450, 1E-6, 0, 1.0);
  minuit->DefineParameter(15, "CuBox + CuFrame Sx th232 10", 0.0143126175, 1E-6, 0, 1.0);
  minuit->DefineParameter(16, "CuBox + CuFrame Sx u238 10 ", 0.0012999801, 1E-6, 0, 1.0);
  minuit->DefineParameter(17, "CuBox + CuFrame Sx pb210 10", 0.0043297445, 1E-6, 0, 1.0);
  minuit->DefineParameter(18, "CuBox + CuFrame Sx pb210 0.1", 0.0044209233, 1E-6, 0, 1.0);
  minuit->DefineParameter(19, "CuBox + CuFrame Sx pb210 0.01", 0.0181031801, 1E-6, 0, 1.0);

  minuit->DefineParameter(20, "PbRom th232",  0., 1E-6, 0, 1.0);
  minuit->DefineParameter(21, "PbRom u238",  0., 1E-6, 0, 1.0);
  minuit->DefineParameter(22, "PbRom co60",  0., 1E-6, 0, 1.0);
  minuit->DefineParameter(23, "PbRom k40",  0., 1E-6, 0, 1.0);

  minuit->DefineParameter(24, "OVC th232",  0., 1E-6, 0, 1.0);
  minuit->DefineParameter(25, "OVC u238",  0., 1E-6, 0, 1.0);
  minuit->DefineParameter(26, "OVC co60",  0., 1E-6, 0, 1.0);    
  minuit->DefineParameter(27, "OVC k40",  0., 1E-6, 0, 1.0);
  minuit->DefineParameter(28, "External Lead bi210", 0., 1E-6, 0, 1.0);

  minuit->DefineParameter(29, "Internal Shields th232",  0., 1E-6, 0, 1.0);
  minuit->DefineParameter(30, "Internal Shields u238",  0., 1E-6, 0, 1.0);
  minuit->DefineParameter(31, "Internal Shields co60",  0., 1E-6, 0, 1.0);
  minuit->DefineParameter(32, "Internal Shields k40",  0., 1E-6, 0, 1.0);

//////////////////////////////////////

   // Uncommend to fix parameters here
   // minuit->FixParameter(0); // TeO2 0nu
   // minuit->FixParameter(1); // TeO2 2nu
   // minuit->FixParameter(2); // TeO2 co60
   // minuit->FixParameter(3); // TeO2 k40
   // minuit->FixParameter(4); // TeO2 po210
   // minuit->FixParameter(5); // TeO2 th232 only
   // minuit->FixParameter(6); // TeO2 th230 only
   // minuit->FixParameter(7); // TeO2 Sx th232 0.01  
   // minuit->FixParameter(8); // TeO2 Sx th232 only 0.01
   // minuit->FixParameter(9); // TeO2 Sx ra228 to pb208 0.01
   // minuit->FixParameter(10); // TeO2 Sx u238 to th230 0.01
   // minuit->FixParameter(11); // TeO2 Sx th230 only 0.01
   // minuit->FixParameter(12); // TeO2 Sx ra226 to pb210 0.01
   // minuit->FixParameter(13); // TeO2 Sx pb210 1  
   // minuit->FixParameter(14); // TeO2 Sx pb210 0.01

   // minuit->FixParameter(15); // CuBox+CuFrame Sx th232 10
   // minuit->FixParameter(16); // CuBox+CuFrame Sx u238 10
   // minuit->FixParameter(17); // CuBox+CuFrame Sx pb210 10
   // minuit->FixParameter(18); // CuBox+CuFrame Sx pb210 0.1
   // minuit->FixParameter(19); // CuBox+CuFrame Sx pb210 0.01

   // minuit->FixParameter(20); // PbRom th232
   // minuit->FixParameter(21); // PbRom u238
   // minuit->FixParameter(22); // PbRom co60
   // minuit->FixParameter(23); // PbRom k40

   // minuit->FixParameter(24); // OVC th232
   // minuit->FixParameter(25); // OVC u238
   // minuit->FixParameter(26); // OVC co60
   // minuit->FixParameter(27); // OVC k40
   // minuit->FixParameter(28); // External Lead bi210 

   // minuit->FixParameter(29); // Internal th232
   // minuit->FixParameter(30); // Internal u238
   // minuit->FixParameter(31); // Internal co60
   // minuit->FixParameter(32); // Internal k40

   // Number of free Parameters (for Chi-squared/NDF calculation only)
   dNumFreeParameters = minuit->GetNumPars() - minuit->GetNumFixedPars();
   // Tell minuit what external function to use 
   minuit->SetFCN(myExternal_FCNAdap);
   int status = minuit->Command("MINImize 500000 0.1"); // Command that actually does the minimization

  minuit->Command("SHOw COVariance");

  // Get final parameters from fit
  for(int i = 0; i < dNParam; i++)
  {
    minuit->GetParameter(i, fParameters[i], fParError[i]);
  }

  // Update model with final parameters
  UpdateModelAdaptive();
  
  dChiSquare = GetChiSquareAdaptive();

  cout << "Total number of calls = " << dNumCalls << "\t" << "ChiSq/NDF = " << dChiSquare/(dFitMaxBinM1+dFitMaxBinM2-dFitMinBinM1-dFitMinBinM2-dNumFreeParameters) << endl; // for M1 and M2
  cout << "ChiSq = " << dChiSquare << "\t" << "NDF = " << (dFitMaxBinM1+dFitMaxBinM2-dFitMinBinM1-dFitMinBinM2-dNumFreeParameters) << endl;
  cout << "Probability = " << TMath::Prob(dChiSquare, (dFitMaxBinM1+dFitMaxBinM2-dFitMinBinM1-dFitMinBinM2-dNumFreeParameters) ) << endl;

// M1
  fModelTotAdapNDBDM1->Add( hAdapTeO20nuM1,              dDataIntegralM1*fParameters[0]);
  fModelTotAdap2NDBDM1->Add( hAdapTeO22nuM1,              dDataIntegralM1*fParameters[1]);
  fModelTotAdapcoM1->Add( hAdapTeO2co60M1,             dDataIntegralM1*fParameters[2]);
  fModelTotAdapkM1->Add( hAdapTeO2k40M1,              dDataIntegralM1*fParameters[3]);
  fModelTotAdapSpoM1->Add( hAdapTeO2po210M1,            dDataIntegralM1*fParameters[4]);
  fModelTotAdapthM1->Add( hAdapTeO2th232onlyM1,        dDataIntegralM1*fParameters[5]);
  fModelTotAdapuM1->Add( hAdapTeO2th230onlyM1,        dDataIntegralM1*fParameters[6]);
  fModelTotAdapthM1->Add( hAdapTeO2Sxth232M1_001,   dDataIntegralM1*fParameters[7]);
  fModelTotAdapthM1->Add( hAdapTeO2Sxth232onlyM1_001,   dDataIntegralM1*fParameters[8]);
  fModelTotAdapthM1->Add( hAdapTeO2Sxra228pb208M1_001, dDataIntegralM1*fParameters[9]);
  fModelTotAdapuM1->Add( hAdapTeO2Sxu238th230M1_001,  dDataIntegralM1*fParameters[10]);
  fModelTotAdapuM1->Add( hAdapTeO2Sxth230onlyM1_001,  dDataIntegralM1*fParameters[11]);
  fModelTotAdapuM1->Add( hAdapTeO2Sxra226pb210M1_001, dDataIntegralM1*fParameters[12]);
  fModelTotAdapSpbM1->Add( hAdapTeO2Sxpb210M1_1,     dDataIntegralM1*fParameters[13]);
  fModelTotAdapSpbM1->Add( hAdapTeO2Sxpb210M1_001,     dDataIntegralM1*fParameters[14]);

  fModelTotAdapthM1->Add( hAdapCuBox_CuFrameth232M1_10,  dDataIntegralM1*fParameters[15]);
  fModelTotAdapuM1->Add( hAdapCuBox_CuFrameu238M1_10,   dDataIntegralM1*fParameters[16]);
  fModelTotAdapSpbM1->Add( hAdapCuBox_CuFramepb210M1_10,  dDataIntegralM1*fParameters[17]);
  fModelTotAdapSpbM1->Add( hAdapCuBox_CuFramepb210M1_01,  dDataIntegralM1*fParameters[18]);
  fModelTotAdapSpbM1->Add( hAdapCuBox_CuFramepb210M1_001,  dDataIntegralM1*fParameters[19]);

  fModelTotAdapthM1->Add( hAdapPbRomth232M1,    dDataIntegralM1*fParameters[20]);
  fModelTotAdapuM1->Add( hAdapPbRomu238M1,     dDataIntegralM1*fParameters[21]);
  fModelTotAdapcoM1->Add( hAdapPbRomco60M1,     dDataIntegralM1*fParameters[22]);
  fModelTotAdapkM1->Add( hAdapPbRomk40M1,      dDataIntegralM1*fParameters[23]);

  fModelTotAdapthM1->Add( hAdapOVCth232M1,    dDataIntegralM1*fParameters[24]);
  fModelTotAdapuM1->Add( hAdapOVCu238M1,     dDataIntegralM1*fParameters[25]);
  fModelTotAdapcoM1->Add( hAdapOVCco60M1,     dDataIntegralM1*fParameters[26]);
  fModelTotAdapkM1->Add( hAdapOVCk40M1,      dDataIntegralM1*fParameters[27]);
  fModelTotAdapbiM1->Add( hAdapExtPbbi210M1,        dDataIntegralM1*fParameters[28]);

  fModelTotAdapthM1->Add( hAdapInternalth232M1,    dDataIntegralM1*fParameters[29]);
  fModelTotAdapuM1->Add( hAdapInternalu238M1,     dDataIntegralM1*fParameters[30]);
  fModelTotAdapcoM1->Add( hAdapInternalco60M1,     dDataIntegralM1*fParameters[31]);
  fModelTotAdapkM1->Add( hAdapInternalk40M1,      dDataIntegralM1*fParameters[32]);


// M2
  fModelTotAdapNDBDM2->Add( hAdapTeO20nuM2,              dDataIntegralM2*fParameters[0]);
  fModelTotAdap2NDBDM2->Add( hAdapTeO22nuM2,              dDataIntegralM2*fParameters[1]);
  fModelTotAdapcoM2->Add( hAdapTeO2co60M2,             dDataIntegralM2*fParameters[2]);
  fModelTotAdapkM2->Add( hAdapTeO2k40M2,              dDataIntegralM2*fParameters[3]);
  fModelTotAdapSpoM2->Add( hAdapTeO2po210M2,            dDataIntegralM2*fParameters[4]);
  fModelTotAdapthM2->Add( hAdapTeO2th232onlyM2,        dDataIntegralM2*fParameters[5]);
  fModelTotAdapuM2->Add( hAdapTeO2th230onlyM2,        dDataIntegralM2*fParameters[6]);
  fModelTotAdapthM2->Add( hAdapTeO2Sxth232M2_001,   dDataIntegralM2*fParameters[7]);
  fModelTotAdapthM2->Add( hAdapTeO2Sxth232onlyM2_001,   dDataIntegralM2*fParameters[8]);
  fModelTotAdapthM2->Add( hAdapTeO2Sxra228pb208M2_001, dDataIntegralM2*fParameters[9]);
  fModelTotAdapuM2->Add( hAdapTeO2Sxu238th230M2_001,  dDataIntegralM2*fParameters[10]);
  fModelTotAdapuM2->Add( hAdapTeO2Sxth230onlyM2_001,  dDataIntegralM2*fParameters[11]);
  fModelTotAdapuM2->Add( hAdapTeO2Sxra226pb210M2_001, dDataIntegralM2*fParameters[12]);
  fModelTotAdapSpbM2->Add( hAdapTeO2Sxpb210M2_1,     dDataIntegralM2*fParameters[13]);
  fModelTotAdapSpbM2->Add( hAdapTeO2Sxpb210M2_001,     dDataIntegralM2*fParameters[14]);

  fModelTotAdapthM2->Add( hAdapCuBox_CuFrameth232M2_10,  dDataIntegralM2*fParameters[15]);
  fModelTotAdapuM2->Add( hAdapCuBox_CuFrameu238M2_10,   dDataIntegralM2*fParameters[16]);
  fModelTotAdapSpbM2->Add( hAdapCuBox_CuFramepb210M2_10,  dDataIntegralM2*fParameters[17]);
  fModelTotAdapSpbM2->Add( hAdapCuBox_CuFramepb210M2_01,  dDataIntegralM2*fParameters[18]);
  fModelTotAdapSpbM2->Add( hAdapCuBox_CuFramepb210M2_001,  dDataIntegralM2*fParameters[19]);

  fModelTotAdapthM2->Add( hAdapPbRomth232M2,    dDataIntegralM2*fParameters[20]);
  fModelTotAdapuM2->Add( hAdapPbRomu238M2,     dDataIntegralM2*fParameters[21]);
  fModelTotAdapcoM2->Add( hAdapPbRomco60M2,     dDataIntegralM2*fParameters[22]);
  fModelTotAdapkM2->Add( hAdapPbRomk40M2,      dDataIntegralM2*fParameters[23]);

  fModelTotAdapthM2->Add( hAdapOVCth232M2,    dDataIntegralM2*fParameters[24]);
  fModelTotAdapuM2->Add( hAdapOVCu238M2,     dDataIntegralM2*fParameters[25]);
  fModelTotAdapcoM2->Add( hAdapOVCco60M2,     dDataIntegralM2*fParameters[26]);
  fModelTotAdapkM2->Add( hAdapOVCk40M2,      dDataIntegralM2*fParameters[27]);
  fModelTotAdapbiM2->Add( hAdapExtPbbi210M2,        dDataIntegralM2*fParameters[28]);

  fModelTotAdapthM2->Add( hAdapInternalth232M2,    dDataIntegralM2*fParameters[29]);
  fModelTotAdapuM2->Add( hAdapInternalu238M2,     dDataIntegralM2*fParameters[30]);
  fModelTotAdapcoM2->Add( hAdapInternalco60M2,     dDataIntegralM2*fParameters[31]);
  fModelTotAdapkM2->Add( hAdapInternalk40M2,      dDataIntegralM2*fParameters[32]);

  ///// Draw Data M1
  fAdapDataHistoM1->SetLineColor(kBlack);
  fAdapDataHistoM1->SetLineWidth(2);
  fAdapDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM1->GetYaxis()->SetTitle("Counts/Bin");
  // fAdapDataHistoM1->SetMaximum(90000);
  // fAdapDataHistoM1->GetXaxis()->SetRange(1, fAdapDataHistoM1->FindBin(3000));


  fModelTotAdapM1->SetLineColor(2);
  fModelTotAdapM1->SetLineWidth(1);
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
  fModelTotAdapbiM1->SetLineColor(5);
  fModelTotAdapbiM1->SetLineStyle(2);
  // fModelTotAdapmnM1->SetLineColor(40);
  // fModelTotAdapmnM1->SetLineStyle(2);

  fModelTotAdapSpoM1->SetLineStyle(2);
  fModelTotAdapSpoM1->SetLineColor(44);

  fModelTotAdapSpbM1->SetLineStyle(2);
  fModelTotAdapSpbM1->SetLineColor(8);

  // fModelTotAdapteo2M1->SetLineStyle(2);
  // fModelTotAdapteo2M1->SetLineColor(41);

  // fModelTotAdapteo2M1->Draw("SAME");
  // fModelTotAdapPbM1->Draw("SAME");

  TLegend *legfit1 = new TLegend(0.7,0.7,0.95,0.95);
  legfit1->SetFillColor(0);
  legfit1->SetTextSize(0.02);
  legfit1->AddEntry(fModelTotAdapM1, "Total PDF", "l");
  legfit1->AddEntry(fModelTotAdapthM1, "Total th-232", "l");
  legfit1->AddEntry(fModelTotAdapuM1, "Total u-238", "l");
  legfit1->AddEntry(fModelTotAdapkM1, "Total k-40", "l");
  legfit1->AddEntry(fModelTotAdapcoM1, "Total co-60", "l");
  legfit1->AddEntry(fModelTotAdapNDBDM1, "NDBD", "l");
  legfit1->AddEntry(fModelTotAdap2NDBDM1, "2NDBD", "l");
  legfit1->AddEntry(fModelTotAdapSpoM1, "Total po-210", "l");
  legfit1->AddEntry(fModelTotAdapSpbM1, "Surface pb-210", "l");
  legfit1->AddEntry(fModelTotAdapbiM1, "Total External bi-210", "l");  


  int nHisto = 2;  
  double width1 = 0.02;
  double width2 = 0.98;
  double canBotMargin = 0.02;
  double canTopMargin = 0.02;
  double padHeight = (1.-canTopMargin-canBotMargin)/nHisto;


  TCanvas *cadap1 = new TCanvas("cadap1", "cadap1", 1200, 800);
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
  p2m1->SetTopMargin(0.05);
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
  hRatioM1_3->GetXaxis()->SetRange(fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(7999));
  hRatioM1_3->SetMarkerStyle(6);
  hRatioM1_3->GetXaxis()->SetTitle("Energy (keV)");
  hRatioM1_3->GetYaxis()->SetTitle("Ratio (Data/MC)");
  hRatioM1_3->GetXaxis()->SetLabelSize(0.07);
  hRatioM1_3->GetYaxis()->SetLabelSize(0.07);
  hRatioM1_3->GetXaxis()->SetTitleSize(0.07);
  hRatioM1_3->GetYaxis()->SetTitleSize(0.07);
  hRatioM1_3->GetYaxis()->SetTitleOffset(0.45);
  hRatioM1_1->SetFillColor(kGreen);
  hRatioM1_2->SetFillColor(kYellow);
  hRatioM1_3->SetFillColor(kRed);
  hRatioM1_3->Draw("pe2");
  hRatioM1_2->Draw("SAME e2");
  hRatioM1_1->Draw("SAME e2");
  LineM1->DrawLine(500, 1, 7999, 1);

  p2m1->cd();
  p2m1->SetLogy();
  fAdapDataHistoM1->GetXaxis()->SetRange(fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(7999));
  fAdapDataHistoM1->Draw("E");
  fModelTotAdapM1->Draw("SAME");
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
  legfit1->Draw();



  ///// Draw Data M2
  fAdapDataHistoM2->SetLineColor(kBlack);
  fAdapDataHistoM2->SetLineWidth(2);
  fAdapDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM2->GetYaxis()->SetTitle("Counts/Bin");
  // fAdapDataHistoM2->SetMaximum(9000);
  // fAdapDataHistoM2->GetXaxis()->SetRange(1, fAdapDataHistoM2->FindBin(3000));
  
  fModelTotAdapM2->SetLineColor(2);
  fModelTotAdapM2->SetLineWidth(1);
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
  fModelTotAdapSpoM2->SetLineColor(44);

  fModelTotAdapSpbM2->SetLineStyle(2);
  fModelTotAdapSpbM2->SetLineColor(8);

  // fModelTotAdapteo2M2->SetLineStyle(2);
  // fModelTotAdapteo2M2->SetLineColor(41);

  TLegend *legfit2 = new TLegend(0.7,0.7,0.95,0.95);
  legfit2->SetFillColor(0);
  legfit2->SetTextSize(0.02);
  legfit2->AddEntry(fModelTotAdapM2, "Total PDF", "l");
  legfit2->AddEntry(fModelTotAdapthM2, "Total th-232", "l");
  legfit2->AddEntry(fModelTotAdapuM2, "Total u-238", "l");
  legfit2->AddEntry(fModelTotAdapkM2, "Total k-40", "l");
  legfit2->AddEntry(fModelTotAdapcoM2, "Total co-60", "l");
  legfit2->AddEntry(fModelTotAdapNDBDM2, "NDBD", "l");
  legfit2->AddEntry(fModelTotAdap2NDBDM2, "2NDBD", "l");
  legfit2->AddEntry(fModelTotAdapSpoM2, "Total po-210", "l");
  legfit2->AddEntry(fModelTotAdapSpbM2, "Surface pb-210", "l");
  legfit2->AddEntry(fModelTotAdapbiM2, "Total External bi-210", "l");  


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
  p2m2->SetTopMargin(0.05);
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
  hRatioM2_3->GetXaxis()->SetRange(fAdapDataHistoM2->FindBin(500), fAdapDataHistoM2->FindBin(7999));
  hRatioM2_3->SetMarkerStyle(6);
  hRatioM2_3->GetXaxis()->SetTitle("Energy (keV)");
  hRatioM2_3->GetYaxis()->SetTitle("Ratio (Data/MC)");  
  hRatioM2_3->GetXaxis()->SetLabelSize(0.07);
  hRatioM2_3->GetYaxis()->SetLabelSize(0.07);
  hRatioM2_3->GetXaxis()->SetTitleSize(0.07);
  hRatioM2_3->GetYaxis()->SetTitleSize(0.07);
  hRatioM2_3->GetYaxis()->SetTitleOffset(0.45);
  hRatioM2_3->SetFillColor(kRed);
  hRatioM2_1->SetFillColor(kGreen);
  hRatioM2_2->SetFillColor(kYellow);
  hRatioM2_3->Draw("pE2");
  hRatioM2_2->Draw("SAME e2");
  hRatioM2_1->Draw("SAME e2");
  LineM1->DrawLine(500, 1, 7999, 1);



  p2m2->cd();
  p2m2->SetLogy();
  fAdapDataHistoM2->GetXaxis()->SetRange(fAdapDataHistoM2->FindBin(500), fAdapDataHistoM2->FindBin(7999));
  fAdapDataHistoM2->Draw("E");
  fModelTotAdapM2->Draw("SAME");
  fModelTotAdapthM2->Draw("SAME");
  fModelTotAdapuM2->Draw("SAME");
  fModelTotAdapkM2->Draw("SAME");
  fModelTotAdapcoM2->Draw("SAME");
  fModelTotAdapNDBDM2->Draw("SAME");
  fModelTotAdap2NDBDM2->Draw("SAME");
  fModelTotAdapbiM2->Draw("SAME");
  fModelTotAdapmnM2->Draw("SAME");
  fModelTotAdapSpbM2->Draw("SAME");
  fModelTotAdapSpoM2->Draw("SAME");
  fModelTotAdappbM2->Draw("SAME");
  legfit2->Draw();




  // Residuals
  TCanvas *cResidual1 = new TCanvas("cResidual1", "cResidual1", 1200, 800);
  hResidualGausM1 = new TH1D("hResidualGausM1", "M1", 100, -50, 50);
  hResidualDistM1 = CalculateResidualsAdaptive(fAdapDataHistoM1, fModelTotAdapM1, hResidualGausM1, dFitMinBinM1, dFitMaxBinM1, 1);
  hResidualDistM1->SetLineColor(kBlack);
  hResidualDistM1->SetName("Residuals");
  hResidualDistM1->SetTitle("");
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
  hResidualDistM2->SetTitle("");
  hResidualDistM2->SetMarkerStyle(25);
  hResidualDistM2->GetXaxis()->SetTitle("Energy (keV)");
  // hResidualDistM2->GetXaxis()->SetTitleSize(0.04);
  // hResidualDistM2->GetXaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetTitleSize(0.04); 
  // hResidualDistM2->GetXaxis()->SetRange(1, fAdapDataHistoM2->FindBin(3000));
  hResidualDistM2->GetYaxis()->SetTitle("Residuals (#sigma)");

  hResidualDistM2->Draw();

  TCanvas *cres1 = new TCanvas("cres1", "cres1", 1600, 600);
  cres1->Divide(2,1);
  cres1->cd(1);
  hResidualGausM1->Fit("gaus");
  hResidualGausM1->Draw();
  cres1->cd(2);
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

  dResidualRMSTot = TMath::Sqrt( (dResidualRMSM1 + dResidualRMSM2)/ (hResidualDistM1->GetNbinsX()+hResidualDistM2->GetNbinsX()) );


  dResidualRMSM1 = TMath::Sqrt(dResidualRMSM1/hResidualDistM1->GetNbinsX());
  dResidualRMSM2 = TMath::Sqrt(dResidualRMSM2/hResidualDistM2->GetNbinsX());


  double dROIRange = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2570))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2570)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2486)); 
  // Output integrals of stuff for limits
  cout << "ROI range: " << fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2486)) << " " << fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2570))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2570)) << " keV" << endl; // 2486 to 2572
  cout << "Integral Data in ROI: " << fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt( fAdapDataHistoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << "Integral Total PDF in ROI: " << fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width") << " +/- " << sqrt( fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << "Integral Total Th-232 PDF in ROI: " << fModelTotAdapthM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt( fModelTotAdapthM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << "Integral Total U-238 PDF in ROI: " << fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << "Integral Total Co PDF in ROI: " << fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << "Integral Total Pb-210 PDF in ROI: " << fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << "Integral Total Po-210 PDF in ROI: " << fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;  
  cout << "Integral Total 0NDBD PDF in ROI: " << fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << endl;
  cout << "Integral Data in ROI (counts/keV): " << fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt( fAdapDataHistoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  cout << "Integral Total PDF in ROI (counts/keV): " << fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt( fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  cout << "Integral Total Th-232 PDF in ROI (counts/keV): " << fModelTotAdapthM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt( fModelTotAdapthM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  cout << "Integral Total U-238 PDF in ROI (counts/keV): " << fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  cout << "Integral Total Co PDF in ROI (counts/keV): " << fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  cout << "Integral Total Pb-210 PDF in ROI (counts/keV): " << fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  cout << "Integral Total Po-210 PDF in ROI (counts/keV): " << fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;  
  cout << "Integral Total 0NDBD PDF in ROI (counts/keV): " << fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  cout << "Number of 2nbb: " << fParameters[1]*dDataIntegralM1 << " +/- " << fParError[1]*dDataIntegralM1 << "\t 2nbb half life: " << (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[1]*dDataIntegralM1) << " +/- " << (fParError[1]/fParameters[1]) * (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[1]*dDataIntegralM1) << endl;
  cout << "Counts in 2nbb (M1 + M2): " << fModelTotAdap2NDBDM1->Integral(1, fAdapDataHistoM1->FindBin(3000), "width") + fModelTotAdap2NDBDM2->Integral(1, fAdapDataHistoM2->FindBin(3000) , "width")/2 << "\t Half-Life " << (0.69314718056)*(4.726e25 * dLivetimeYr)/(fModelTotAdap2NDBDM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width") + fModelTotAdap2NDBDM2->Integral(1, fAdapDataHistoM2->FindBin(2700) , "width")/2) << endl;

  cout << "Residual RMS (Tot): " << dResidualRMSTot << endl;
  cout << "Residual RMS (M1): " << dResidualRMSM1 << "\t" << "Residual RMS (M2): " << dResidualRMSM2 << endl;

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


  // This only gets the number of counts corrected for detector efficiency
  for(int i = 0; i < TBackgroundModel::dNParam; i++)
  {
    fParActivity[i] = fParameters[i]*dDataIntegralM1/fParEfficiencyM1[i];
  }

  // Saving plots
  // cadap1->SaveAs(Form("%s/FitResults/EffCorr/FitM1_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cadap2->SaveAs(Form("%s/FitResults/EffCorr/FitM2_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cResidual1->SaveAs(Form("%s/FitResults/EffCorr/FitM1Residual_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cResidual2->SaveAs(Form("%s/FitResults/EffCorr/FitM2Residual_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cres1->SaveAs(Form("%s/FitResults/EffCorr/FitResidualDist_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
/*
  cadap1->SaveAs(Form("%s/FitResults/EffCorr/FitM1_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  cadap2->SaveAs(Form("%s/FitResults/EffCorr/FitM2_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  cResidual1->SaveAs(Form("%s/FitResults/EffCorr/FitM1Residual_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  cResidual2->SaveAs(Form("%s/FitResults/EffCorr/FitM2Residual_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  cres1->SaveAs(Form("%s/FitResults/EffCorr/FitResidualDist_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
*/

  // Kernal Convolution
  // TH1D *hKernalConvM1 = new TH1D("hKernalConvM1", "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  // Kernal(fModelTotAdapM1, hKernalConvM1);

  // TCanvas *cKernal1 = new TCanvas("cKernal1", "cKernal1", 1200, 1200);
  // hKernalConvM1->Draw();


  return true;

}


void TBackgroundModel::DrawMC()
{
// Draws all MC spectra, must Initialize first!

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TLegend *legspb1 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legspb2 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legspb2sum = new TLegend(0.65,0.7,0.95,0.95);

  TLegend *legpb1 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legpb2 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legpb2sum = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legpo1 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legpo2 = new TLegend(0.65,0.7,0.95,0.95);
  TLegend *legpo2sum = new TLegend(0.65,0.7,0.95,0.95);


  TCanvas *cSPb2101 = new TCanvas("cSPb2101", "cSPb2101", 1200, 800);
  cSPb2101->SetLogy();

  hCuBox_CuFramepb210M1_10->SetLineColor(1);
  hCuBox_CuFramepb210M1_1->SetLineColor(2);
  hCuBox_CuFramepb210M1_01->SetLineColor(3);
  hCuBox_CuFramepb210M1_001->SetLineColor(4);

  hCuBox_CuFramepb210M1_001->GetXaxis()->SetTitle("Energy (keV)");
  hCuBox_CuFramepb210M1_001->GetYaxis()->SetTitle("Probability");  
  hCuBox_CuFramepb210M1_001->DrawNormalized();
  hCuBox_CuFramepb210M1_01->DrawNormalized("SAME");
  hCuBox_CuFramepb210M1_1->DrawNormalized("SAME");
  hCuBox_CuFramepb210M1_10->DrawNormalized("SAME");


  legspb1->AddEntry(hCuBox_CuFramepb210M1_001, "CuBox+CuFrame Sx Pb210 0.01", "l");  
  legspb1->AddEntry(hCuBox_CuFramepb210M1_01, "CuBox+CuFrame Sx Pb210 0.1", "l");
  legspb1->AddEntry(hCuBox_CuFramepb210M1_1, "CuBox+CuFrame Sx Pb210 1", "l");
  legspb1->AddEntry(hCuBox_CuFramepb210M1_10, "CuBox+CuFrame Sx Pb210 10", "l");
  legspb1->Draw();


  TCanvas *cSPb2102 = new TCanvas("cSPb2102", "cSPb2102", 1200, 800);
  cSPb2102->SetLogy();

  hCuBox_CuFramepb210M2_10->SetLineColor(1);
  hCuBox_CuFramepb210M2_1->SetLineColor(2);
  hCuBox_CuFramepb210M2_01->SetLineColor(3);
  hCuBox_CuFramepb210M2_001->SetLineColor(4);

  hCuBox_CuFramepb210M2_001->GetXaxis()->SetTitle("Energy (keV)");
  hCuBox_CuFramepb210M2_001->GetYaxis()->SetTitle("Probability");  
  hCuBox_CuFramepb210M2_001->DrawNormalized();
  hCuBox_CuFramepb210M2_01->DrawNormalized("SAME");
  hCuBox_CuFramepb210M2_1->DrawNormalized("SAME");
  hCuBox_CuFramepb210M2_10->DrawNormalized("SAME");



  legspb2->AddEntry(hCuBox_CuFramepb210M2_001, "CuBox+CuFrame Sx Pb210 0.01", "l");  
  legspb2->AddEntry(hCuBox_CuFramepb210M2_01, "CuBox+CuFrame Sx Pb210 0.1", "l");
  legspb2->AddEntry(hCuBox_CuFramepb210M2_1, "CuBox+CuFrame Sx Pb210 1", "l");
  legspb2->AddEntry(hCuBox_CuFramepb210M2_10, "CuBox+CuFrame Sx Pb210 10", "l");
  legspb2->Draw();


  TCanvas *cSPb2102Sum = new TCanvas("cSPb2102Sum", "cSPb2102Sum", 1200, 800);
  cSPb2102Sum->SetLogy();

  hCuBox_CuFramepb210M2Sum_10->SetLineColor(1);
  hCuBox_CuFramepb210M2Sum_1->SetLineColor(2);
  hCuBox_CuFramepb210M2Sum_01->SetLineColor(3);
  hCuBox_CuFramepb210M2Sum_001->SetLineColor(4);

  hCuBox_CuFramepb210M2Sum_001->GetXaxis()->SetTitle("Energy (keV)");
  hCuBox_CuFramepb210M2Sum_001->GetYaxis()->SetTitle("Probability");  
  hCuBox_CuFramepb210M2Sum_001->DrawNormalized();
  hCuBox_CuFramepb210M2Sum_01->DrawNormalized("SAME");
  hCuBox_CuFramepb210M2Sum_1->DrawNormalized("SAME");
  hCuBox_CuFramepb210M2Sum_10->DrawNormalized("SAME");



  legspb2sum->AddEntry(hCuBox_CuFramepb210M2Sum_001, "CuBox+CuFrame Sx Pb210 0.01", "l");  
  legspb2sum->AddEntry(hCuBox_CuFramepb210M2Sum_01, "CuBox+CuFrame Sx Pb210 0.1", "l");
  legspb2sum->AddEntry(hCuBox_CuFramepb210M2Sum_1, "CuBox+CuFrame Sx Pb210 1", "l");
  legspb2sum->AddEntry(hCuBox_CuFramepb210M2Sum_10, "CuBox+CuFrame Sx Pb210 10", "l");
  legspb2sum->Draw();
/*
  TCanvas *cPb2101 = new TCanvas("cPb2101", "cPb2101", 1200, 800);
  cPb2101->SetLogy();

  hTeO2Spb210M1_01->SetLineColor(1);
  hTeO2Sxpb210M1_10->SetLineColor(2);
  hTeO2Sxpb210M1_1->SetLineColor(3);
  hTeO2Sxpb210M1_01->SetLineColor(4);
  hTeO2Sxpb210M1_001->SetLineColor(6);
  hTeO2Sxpb210M1_0001->SetLineColor(7);

  hTeO2Spb210M1_01->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2Spb210M1_01->GetYaxis()->SetTitle("Probability");  
  hTeO2Spb210M1_01->DrawNormalized();
  hTeO2Sxpb210M1_10->DrawNormalized("SAME");
  hTeO2Sxpb210M1_1->DrawNormalized("SAME");
  hTeO2Sxpb210M1_01->DrawNormalized("SAME");
  hTeO2Sxpb210M1_001->DrawNormalized("SAME");
  hTeO2Sxpb210M1_0001->DrawNormalized("SAME");

  legpb1->AddEntry(hTeO2Spb210M1_01, "TeO2 S Pb210 0.1", "l");  
  legpb1->AddEntry(hTeO2Sxpb210M1_10, "TeO2 Sx Pb210 10", "l");
  legpb1->AddEntry(hTeO2Sxpb210M1_1, "TeO2 Sx Pb210 1", "l");
  legpb1->AddEntry(hTeO2Sxpb210M1_01, "TeO2 Sx Pb210 0.1", "l");
  legpb1->AddEntry(hTeO2Sxpb210M1_001, "TeO2 Sx Pb210 0.01", "l");
  legpb1->AddEntry(hTeO2Sxpb210M1_0001, "TeO2 Sx Pb210 0.001", "l");
  legpb1->Draw();

  TCanvas *cPb2102 = new TCanvas("cPb2102", "cPb2102", 1200, 800);
  cPb2102->SetLogy();

  hTeO2Spb210M2_01->SetLineColor(1);
  hTeO2Sxpb210M2_10->SetLineColor(2);
  hTeO2Sxpb210M2_1->SetLineColor(3);
  hTeO2Sxpb210M2_01->SetLineColor(4);
  hTeO2Sxpb210M2_001->SetLineColor(6);
  hTeO2Sxpb210M2_0001->SetLineColor(7);

  hTeO2Spb210M2_01->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2Spb210M2_01->GetYaxis()->SetTitle("Probability");  
  hTeO2Spb210M2_01->DrawNormalized();
  hTeO2Sxpb210M2_10->DrawNormalized("SAME");
  hTeO2Sxpb210M2_1->DrawNormalized("SAME");
  hTeO2Sxpb210M2_01->DrawNormalized("SAME");
  hTeO2Sxpb210M2_001->DrawNormalized("SAME");
  hTeO2Sxpb210M2_0001->DrawNormalized("SAME");

  legpb2->AddEntry(hTeO2Spb210M2_01, "TeO2 S Pb210 0.1", "l");  
  legpb2->AddEntry(hTeO2Sxpb210M2_10, "TeO2 Sx Pb210 10", "l");
  legpb2->AddEntry(hTeO2Sxpb210M2_1, "TeO2 Sx Pb210 1", "l");
  legpb2->AddEntry(hTeO2Sxpb210M2_01, "TeO2 Sx Pb210 0.1", "l");
  legpb2->AddEntry(hTeO2Sxpb210M2_001, "TeO2 Sx Pb210 0.01", "l");
  legpb2->AddEntry(hTeO2Sxpb210M2_0001, "TeO2 Sx Pb210 0.001", "l");
  legpb2->Draw();

  TCanvas *cPb2102Sum = new TCanvas("cPb2102Sum", "cPb2102Sum", 1200, 800);
  cPb2102Sum->SetLogy();

  hTeO2Spb210M2Sum_01->SetLineColor(1);
  hTeO2Sxpb210M2Sum_10->SetLineColor(2);
  hTeO2Sxpb210M2Sum_1->SetLineColor(3);
  hTeO2Sxpb210M2Sum_01->SetLineColor(4);
  hTeO2Sxpb210M2Sum_001->SetLineColor(6);
  hTeO2Sxpb210M2Sum_0001->SetLineColor(7);

  hTeO2Spb210M2Sum_01->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2Spb210M2Sum_01->GetYaxis()->SetTitle("Probability");  
  hTeO2Spb210M2Sum_01->DrawNormalized();
  hTeO2Sxpb210M2Sum_10->DrawNormalized("SAME");
  hTeO2Sxpb210M2Sum_1->DrawNormalized("SAME");
  hTeO2Sxpb210M2Sum_01->DrawNormalized("SAME");
  hTeO2Sxpb210M2Sum_001->DrawNormalized("SAME");
  hTeO2Sxpb210M2Sum_0001->DrawNormalized("SAME");

  legpb2sum->AddEntry(hTeO2Spb210M2Sum_01, "TeO2 S Pb210 0.1", "l");  
  legpb2sum->AddEntry(hTeO2Sxpb210M2Sum_10, "TeO2 Sx Pb210 10", "l");
  legpb2sum->AddEntry(hTeO2Sxpb210M2Sum_1, "TeO2 Sx Pb210 1", "l");
  legpb2sum->AddEntry(hTeO2Sxpb210M2Sum_01, "TeO2 Sx Pb210 0.1", "l");
  legpb2sum->AddEntry(hTeO2Sxpb210M2Sum_001, "TeO2 Sx Pb210 0.01", "l");
  legpb2sum->AddEntry(hTeO2Sxpb210M2Sum_0001, "TeO2 Sx Pb210 0.001", "l");
  legpb2sum->Draw();


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



  TCanvas *cPo2102Sum = new TCanvas("cPo2102Sum", "cPo2102Sum", 1200, 800);
  cPo2102Sum->SetLogy();

  hTeO2Spo210M2Sum_01->SetLineColor(1);
  hTeO2Spo210M2Sum_001->SetLineColor(2);
  hTeO2Sxpo210M2Sum_1->SetLineColor(3);
  hTeO2Sxpo210M2Sum_01->SetLineColor(4);
  hTeO2Sxpo210M2Sum_001->SetLineColor(6);

  hTeO2Spo210M2Sum_01->GetXaxis()->SetTitle("Energy (keV)");
  hTeO2Spo210M2Sum_01->GetYaxis()->SetTitle("Probability");  
  hTeO2Spo210M2Sum_01->DrawNormalized();
  hTeO2Spo210M2Sum_001->DrawNormalized("SAME");
  hTeO2Sxpo210M2Sum_1->DrawNormalized("SAME");
  hTeO2Sxpo210M2Sum_01->DrawNormalized("SAME");
  hTeO2Sxpo210M2Sum_001->DrawNormalized("SAME");

  legpo2sum->AddEntry(hTeO2Spo210M2Sum_01, "TeO2 S po210 0.1", "l");  
  legpo2sum->AddEntry(hTeO2Spo210M2Sum_001, "TeO2 S po210 0.01", "l");  
  legpo2sum->AddEntry(hTeO2Sxpo210M2Sum_1, "TeO2 Sx po210 1", "l");
  legpo2sum->AddEntry(hTeO2Sxpo210M2Sum_01, "TeO2 Sx po210 0.1", "l");
  legpo2sum->AddEntry(hTeO2Sxpo210M2Sum_001, "TeO2 Sx po210 0.01", "l");
  legpo2sum->Draw();

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


// Txt file with useful stuff
void TBackgroundModel::LatexResultTable(double fValue)
{

  OutFile.open(Form("%s/FitResults/EffCorr/FitOutputTable_%d_%d_%d.txt", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop ));
  OutFile << "Initial Value of 2nbb: " << fValue << endl;
  OutFile << "Fit Range: " << dFitMin << " to " << dFitMax << endl;
  OutFile << "Base binning: " << dBinBase << endl;
  OutFile << "Events in background spectrum (M1): " << dDataIntegralM1 << endl;
  OutFile << "Events in background spectrum (M2): " << dDataIntegralM2 << endl;
  OutFile << "Events in background spectrum (M2Sum): " << dDataIntegralM2Sum << endl;
  OutFile << "Livetime of background: " << dLivetimeYr << endl;
  OutFile << "Total number of calls = " << dNumCalls << "\t" << "ChiSq/NDF = " << dChiSquare/(dFitMaxBinM1+dFitMaxBinM2-dFitMinBinM1-dFitMinBinM2-dNumFreeParameters) << endl; // for M1 and M2
  OutFile << "ChiSq = " << dChiSquare << "\t" << "NDF = " << (dFitMaxBinM1+dFitMaxBinM2-dFitMinBinM1-dFitMinBinM2-dNumFreeParameters) << endl;
  OutFile << "Probability = " << TMath::Prob(dChiSquare, (dFitMaxBinM1+dFitMaxBinM2-dFitMinBinM1-dFitMinBinM2-dNumFreeParameters) ) << endl;

  OutFile << endl;
  OutFile << endl;

  double dROIRange = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2570))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2570)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2486)); 
  // Output integrals of stuff for limits
  OutFile << "ROI range: " << fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2486)) << " " << fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2570))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2570)) << " keV" << endl; // 2486 to 2572
  OutFile << "Integral Data in ROI: " << fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt( fAdapDataHistoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total PDF in ROI: " << fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width") << " +/- " << sqrt( fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total Th-232 PDF in ROI: " << fModelTotAdapthM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt( fModelTotAdapthM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total U-238 PDF in ROI: " << fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total Co PDF in ROI: " << fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total Pb-210 PDF in ROI: " << fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total Po-210 PDF in ROI: " << fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;  
  // OutFile << "Integral Total 2NDBD PDF in ROI: " << fModelTotAdap2NDBDM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdap2NDBDM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total 0NDBD PDF in ROI: " << fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << endl;
  OutFile << "Integral Data in ROI (counts/keV): " << fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt( fAdapDataHistoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Integral Total PDF in ROI (counts/keV): " << fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt( fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Integral Total Th-232 PDF in ROI (counts/keV): " << fModelTotAdapthM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt( fModelTotAdapthM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Integral Total U-238 PDF in ROI (counts/keV): " << fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Integral Total Co PDF in ROI (counts/keV): " << fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Integral Total Pb-210 PDF in ROI (counts/keV): " << fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Integral Total Po-210 PDF in ROI (counts/keV): " << fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;  
  // OutFile << "Integral Total 2NDBD PDF in ROI (counts/keV): " << fModelTotAdap2NDBDM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdap2NDBDM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Integral Total 0NDBD PDF in ROI (counts/keV): " << fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2486),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Number of 2nbb: " << fParameters[1]*dDataIntegralM1 << " +/- " << fParError[1]*dDataIntegralM1 << "\t 2nbb half life: " << (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[1]*dDataIntegralM1) << " +/- " << fParError[1]/fParameters[1] * (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[1]*dDataIntegralM1) << endl;
  OutFile << "Counts in 2nbb (M1 + M2): " << fModelTotAdap2NDBDM1->Integral(1, fAdapDataHistoM1->FindBin(3000), "width") + fModelTotAdap2NDBDM2->Integral(1, fAdapDataHistoM2->FindBin(3000) , "width")/2 << "\t Half-Life " << (0.69314718056)*(4.726e25 * dLivetimeYr)/(fModelTotAdap2NDBDM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width") + fModelTotAdap2NDBDM2->Integral(1, fAdapDataHistoM2->FindBin(2700) , "width")/2) << endl;
  
  OutFile << "Residual RMS (Tot): " << dResidualRMSTot << endl;
  OutFile << "Residual RMS (M1): " << dResidualRMSM1 << "\t" << "Residual RMS (M2): " << dResidualRMSM2 << endl;
  OutFile << endl;
  OutFile << endl;

  // for(int i = 0; i < TBackgroundModel::dNParam; i++)
  // {
  //   OutFile << minuit->fCpnam[i] << Form(" & %.10f$\\pm$%.10f \\\\ ", fParameters[i], fParError[i] ) << endl;
  // }

  OutFile << endl;
  OutFile << endl;

  // for(int i = 0; i < TBackgroundModel::dNParam; i++)
  // {
  //   OutFile << minuit->fCpnam[i] << " activity: " << fParActivity[i] << endl;
  // }
  OutFile.close();


}

// Adds a random percentage of events of 2nbb into spectrum
void TBackgroundModel::SanityCheck()
{
  double dM1 = fDataHistoM1->Integral(1, 10000/dBinSize);
  double dM2 = fDataHistoM2->Integral(1, 10000/dBinSize);
  // Sanity check, adding to background a set amount of 2nbb events, see if code can reconstruct it properly
  fDataHistoM1->Add(hTeO22nuM1, 1.0*dM1);
  fDataHistoM2->Add(hTeO22nuM2, 1.0*dM2);

  fAdapDataHistoM1->Add(hAdapTeO22nuM1, 1.0*dM1);
  fAdapDataHistoM2->Add(hAdapTeO22nuM2, 1.0*dM2);
  // fDataHistoM2Sum->Add(hTeO22nuM2Sum, 0.5*dDataIntegralM2Sum);
}

// Only run in batch mode and make sure to have the 2nbb normalization FIXED
// Probably best to run Minuit in quiet mode as well
void TBackgroundModel::ProfileNLL(double fBestFitInit, double fBestFitChiSq)
{
  dBestChiSq = fBestFitChiSq; // Chi-Squared from best fit (for ProfileNLL calculation)
  // Do the fit now if no other tests are needed 
  nLoop = 0;
  for(int i = -15; i < 15; i++)
  {
    fInitValues.push_back(fBestFitInit + fBestFitInit/100*i);
  }


  OutPNLL.open(Form("%s/FitResults/ProfileNLL/ProfileNLL_%d_DR%d.C", dSaveDir.c_str(), tTime->GetDate(), dDataSet ));
  OutPNLL << "{" << endl;
  OutPNLL << "vector<double> dX;" << endl;
  OutPNLL << "vector<double> dT;" << endl;

  for(std::vector<double>::const_iterator iter = fInitValues.begin(); iter!=fInitValues.end(); iter++)
  {
    // cout << "Loop: " << nLoop << endl;
    DoTheFitAdaptive(*iter);
    // LatexResultTable(*iter);
    cout << "delta ChiSq = " << dChiSquare - dBestChiSq << endl; // Needs to be entered, otherwise just 0
    OutPNLL << Form("dX.push_back(%f); dT.push_back(%f);", dChiSquare-dBestChiSq, (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[1]*dDataIntegralM1) ) << endl;

    dNumCalls = 0; // Resets number of calls (for saving purposes)
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
  OutPNLL << "g1->GetYaxis()->SetTitle(\"#Delta#chi^{2}\");" << endl;
  OutPNLL << "g1->GetXaxis()->SetTitle(\"t_{1/2} (y)\");" << endl;
  OutPNLL << "g1->Draw(\"AC\");" << endl;
  OutPNLL << "}" << endl;

  OutPNLL.close();
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