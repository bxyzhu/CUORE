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
    Obj->SetParValue(i, x[i]);
  }
  // Implement a method in your class that calculates the quantity you want to minimize, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimize fval
    Obj->UpdateModel();
    fval = Obj->GetChiSquare();
}

TBackgroundModel::TBackgroundModel()
{}

TBackgroundModel::TBackgroundModel(double fFitMin, double fFitMax, int fBinBase, int fDataSet, bool fSave) : 
dNParam(50), dNumCalls(0), dMass(36.75), dMinEnergy(0.), dMaxEnergy(10000.), dBinSize(1), dChiSquare(0.), dBestChiSq(0.), dSecToYears(1./(60*60*24*365))
{

  bSave = fSave;
  tTime = new TDatime();

  // Data directories depending on QCC/local
  dDataDir =  "/Users/brian/macros/CUOREZ/Bkg";
   // dDataDir = "/cuore/user/zhubrian/CUORE0/scratch/";
  // dDataDir = "/home/zhubrian/data/CUORE0";

  dMCDir = "/Users/brian/macros/Simulations/Production";
   // dMCDir = "/cuore/user/zhubrian/MC/Bkg";
  // dMCDir = "/home/zhubrian/data/Simulations";
  
  dSaveDir = "/Users/brian/Dropbox/code/Fitting";
   // dSaveDir = "/cuore/user/zhubrian/";
  // dSaveDir = "/home/zhubrian/code/Fitting";

  // Base number of counts
  dBinBase = fBinBase;
  dDataSet = fDataSet;

  if(fFitMin >= fFitMax)
  {
    cout << "Fit Min >= Fit Max!" << endl;
  }
  // Fitting range
  dFitMin = fFitMin;
  dFitMax = fFitMax;

  dNBins = (dMaxEnergy - dMinEnergy)/ dBinSize;

  // Data Histograms
  fDataHistoTot     = new TH1D("fDataHistoTot",  "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM1      = new TH1D("fDataHistoM1",   "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM2      = new TH1D("fDataHistoM2",   "", dNBins, dMinEnergy, dMaxEnergy);

  // Data cuts 
  qtree = new TChain("qtree");

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

// Mess with rebinning here 
  // Rebinning with 2-4 bins seems to be consistent
  // fDataHistoM1->Rebin(1);
  // fDataHistoM2->Rebin(1);
  dBaseBinSize = dBinSize*1;

/////// Adaptive binning
 // Calculates adaptive binning vectors
  dAdaptiveVectorM1 = AdaptiveBinningM1(fDataHistoM1, dBinBase);
  dAdaptiveBinsM1 = dAdaptiveVectorM1.size() - 1;
  dAdaptiveArrayM1 = &dAdaptiveVectorM1[0];
  dAdaptiveVectorM2 = AdaptiveBinningM2(fDataHistoM2, dBinBase);
  dAdaptiveBinsM2 = dAdaptiveVectorM2.size() - 1;
  dAdaptiveArrayM2 = &dAdaptiveVectorM2[0];
 
  // Adaptive binning data
  fAdapDataHistoM1   = new TH1D("fAdapDataHistoM1",   "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  fAdapDataHistoM2   = new TH1D("fAdapDataHistoM2",   "", dAdaptiveBinsM2, dAdaptiveArrayM2);
  
  hnewM1 = fDataHistoM1->Rebin(dAdaptiveBinsM1, "hnewM1", dAdaptiveArrayM1);
  hnewM2 = fDataHistoM2->Rebin(dAdaptiveBinsM2, "hnewM2", dAdaptiveArrayM2);

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

  dDataIntegralTot = qtree->GetEntries();

  dDataIntegralM1 = fAdapDataHistoM1->Integral("width");
  dDataIntegralM2 = fAdapDataHistoM2->Integral("width");

  dFitMinBinM1 = fAdapDataHistoM1->FindBin(dFitMin);
  dFitMinBinM2 = fAdapDataHistoM2->FindBin(dFitMin);
  dFitMaxBinM1 = fAdapDataHistoM1->FindBin(dFitMax);
  dFitMaxBinM2 = fAdapDataHistoM2->FindBin(dFitMax);

  // Outputs on screen
  cout << "Fit M1 from bin: " << dFitMinBinM1 << " to " << dFitMaxBinM1 << endl;
  cout << "Fit M2 from bin: " << dFitMinBinM2 << " to " << dFitMaxBinM2 << endl;

  cout << "Total Events in background spectrum: " << dDataIntegralTot << endl; 
  cout << "Events in background spectrum (M1): " << dDataIntegralM1 << endl;
  cout << "Events in background spectrum (M2): " << dDataIntegralM2 << endl;
  cout << "Events in fitting region (M1): " << fAdapDataHistoM1->Integral(dFitMinBinM1, dFitMaxBinM1, "width") << endl;
  cout << "Events in fitting region (M2): " << fAdapDataHistoM2->Integral(dFitMinBinM2, dFitMaxBinM2, "width") << endl;  
  cout << "Livetime of background: " << dLivetimeYr << endl;

  // Create Minuit object with dNParam parameters
  minuit = new TMinuit(dNParam);

  // Creates model histograms and then load all of the PDFs from file
  CreateModelHistograms();
  Initialize();
  // Generate and initialize parameters
  GenerateParameters();
  SetParEfficiency();

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

  // Recalculate integrals
  dDataIntegralM1 = fAdapDataHistoM1->Integral("width");
  dDataIntegralM2 = fAdapDataHistoM2->Integral("width");

  cout << "Recalculated integrals: " << endl;
  cout << "Total Events in background spectrum: " << dDataIntegralTot << endl; 
  cout << "Events in background spectrum (M1): " << dDataIntegralM1 << endl;
  cout << "Events in background spectrum (M2): " << dDataIntegralM2 << endl;
  cout << "Livetime of background: " << dLivetimeYr << endl;
*/

}

// Definitely needs updating  
TBackgroundModel::~TBackgroundModel()
{
  delete fDataHistoTot;
  delete fDataHistoM1;
  delete fDataHistoM2;
  delete fAdapDataHistoM1;
  delete fAdapDataHistoM2;

  // Updated 01-20-2015
  // Total PDFs M1
  delete fModelTotM1;
  delete fModelTotAdapM1;

  // Total PDFs M2
  delete fModelTotM2;
  delete fModelTotAdapM2;
}

// Returns vector of bin low-edge values for adaptive binning
// Is this the best way to bin the histograms? I don't know but it's the first way I thought up and never bothered re-thinking it
vector<double> TBackgroundModel::AdaptiveBinningM1(TH1D *h1, int dBinBase)
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
    // Pt peak 3150 - 3400
    if(i >= h1->FindBin(3150) && i < h1->FindBin(3400))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(3150)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(3400-3150)/dBaseBinSize;
    }
    // 4050 - 4200
    if(i >= h1->FindBin(4050) && i < h1->FindBin(4170))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4050)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(4170-4050)/dBaseBinSize;
    }
    // 4200 - 4350
    if(i >= h1->FindBin(4170) && i < h1->FindBin(4350))
    // if(i >= h1->FindBin(4170) && i < h1->FindBin(4400))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4170)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(4350-4170)/dBaseBinSize;
    }    

    // 4680 - 4850
    if(i >= h1->FindBin(4680) && i < h1->FindBin(4850))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4680)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(4850-4680)/dBaseBinSize;
    }

    // 4850 - 4950
    if(i >= h1->FindBin(4850) && i < h1->FindBin(5000))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4850)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5000-4850)/dBaseBinSize;
    }

    // 5200 - 5400
    if(i >= h1->FindBin(5000) && i < h1->FindBin(5100))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5000)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5100-5000)/dBaseBinSize;
    }

    // 5200 - 5400
    if(i >= h1->FindBin(5250) && i < h1->FindBin(5400))
    // if(i >= h1->FindBin(5170) && i < h1->FindBin(5400))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5250)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5400-5250)/dBaseBinSize;
    }
    // 5400 - 5650
    if(i >= h1->FindBin(5400) && i < h1->FindBin(5500))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5400)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5500-5400)/dBaseBinSize;
    }

    // 5400 - 5650
    if(i >= h1->FindBin(5500) && i < h1->FindBin(5650))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5500)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5650-5500)/dBaseBinSize;
    }

    // 5650 - 5800
    if(i >= h1->FindBin(5650) && i < h1->FindBin(5780))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5650)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5780-5650)/dBaseBinSize;
    }

    // 5800 - 6050
    if(i >= h1->FindBin(5780) && i < h1->FindBin(5900))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5780)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5900-5780)/dBaseBinSize;
    }    

    if(i >= h1->FindBin(5900) && i < h1->FindBin(6100))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5900)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(6100-5900)/dBaseBinSize;
    }    

    // 6050 - 6350
    if(i >= h1->FindBin(6100) && i < h1->FindBin(6450))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(6100)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(6450-6100)/dBaseBinSize;
    }    

    // 6700 - 6900
    if(i >= h1->FindBin(6700) && i < h1->FindBin(7000))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(6700)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(7000-6700)/dBaseBinSize;
    }    

    // 7000 - 8000
    if(i >= h1->FindBin(7000) && i < h1->FindBin(8000))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(7000)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(8000-7000)/dBaseBinSize;
    }    

    // 8000 to the end 10000
    if(i >= h1->FindBin(8000) && i < h1->FindBin(10000))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(8000)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(10000-8000)/dBaseBinSize;
    }    

    dDummy = h1->GetBinContent(i);
    dDummyFill += dDummy;

    // Set bin low edge if dummy counts goes above minimum bin content
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

// Separating M1 binning from M2 binning for now
// Really the only difference is I don't have the Pt peak in one bin for M2
vector<double> TBackgroundModel::AdaptiveBinningM2(TH1D *h1, int dBinBase)
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
    // 4050 - 4200
    if(i >= h1->FindBin(4000) && i < h1->FindBin(4150))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4000)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(4150-4000)/dBaseBinSize;
    }
    // 4200 - 4350
    if(i >= h1->FindBin(4150) && i < h1->FindBin(4400))
    // if(i >= h1->FindBin(4170) && i < h1->FindBin(4400))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4150)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(4400-4150)/dBaseBinSize;
    }    

    // 4680 - 4850
    if(i >= h1->FindBin(4680) && i < h1->FindBin(4850))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4680)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(4850-4680)/dBaseBinSize;
    }

    // 4850 - 4950
    if(i >= h1->FindBin(4850) && i < h1->FindBin(5000))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4850)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5000-4850)/dBaseBinSize;
    }

    // 5200 - 5400
    if(i >= h1->FindBin(5000) && i < h1->FindBin(5100))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5000)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5100-5000)/dBaseBinSize;
    }

    // 5200 - 5400
    if(i >= h1->FindBin(5200) && i < h1->FindBin(5450))
    // if(i >= h1->FindBin(5170) && i < h1->FindBin(5400))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5200)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5450-5200)/dBaseBinSize;
    }
    // 5400 - 5650
    if(i >= h1->FindBin(5450) && i < h1->FindBin(5600))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5450)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5600-5450)/dBaseBinSize;
    }
    // 5650 - 5800
    if(i >= h1->FindBin(5650) && i < h1->FindBin(5780))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5650)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5780-5650)/dBaseBinSize;
    }

    // 5800 - 6050
    if(i >= h1->FindBin(5780) && i < h1->FindBin(5900))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5780)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(5900-5780)/dBaseBinSize;
    }    

    if(i >= h1->FindBin(5900) && i < h1->FindBin(6100))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(5900)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(6100-5900)/dBaseBinSize;
    }    

    // 6050 - 6350
    if(i >= h1->FindBin(6100) && i < h1->FindBin(6450))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(6100)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(6450-6100)/dBaseBinSize;
    }    

    // 6700 - 6900
    if(i >= h1->FindBin(6700) && i < h1->FindBin(7000))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(6700)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(7000-6700)/dBaseBinSize;
    }    

    // 7000 - 8000
    if(i >= h1->FindBin(7000) && i < h1->FindBin(8000))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(7000)));
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(8000-7000)/dBaseBinSize;
    }    

    // 8000 to the end 10000
    if(i >= h1->FindBin(8000) && i < h1->FindBin(10000))
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(8000)));
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

// Returns a histogram with the normalized residuals as well as fills a histogram with residual distribution
TH1D *TBackgroundModel::CalculateResiduals(TH1D *fData, TH1D *fModel, TH1D *fResidual, int binMin, int binMax, int dMult)
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

  // Stuff for residual calculations
  double dResidualX = 0, dResidualY = 0, dResidualXErr = 0, dResidualYErr = 0;
  double datam1_i = 0, modelm1_i = 0;

  // Residual plot and distribution
  for(int j = binMin ; j < binMax ; j++)
  {

    if( fData->GetBinCenter(j) >= 3150 && fData->GetBinCenter(j) <= 3400)continue;    

    dResidualX    = fData->GetBinCenter(j);

    // Multiply bin content by bin width (for # of counts)
    datam1_i = fData->GetBinContent(j)*fData->GetBinWidth(j);
    modelm1_i = fModel->GetBinContent(j)*fData->GetBinWidth(j);
    
    dResidualY = (datam1_i - modelm1_i)/TMath::Sqrt(datam1_i);

    hOut->SetBinContent(j, dResidualY);
    hOut->SetBinError(j, 0); // Just for visualization
    fResidual->Fill(dResidualY); // Fills residual distribution
  }

  return hOut;
}

double TBackgroundModel::GetChiSquare()
{
  double chiSquare = 0.;
  double datam1_i, errm1_i;
  double datam2_i, errm2_i;
  double modelm1_i, modelm2_i;


  for(int i = dFitMinBinM1; i < dFitMaxBinM1; i++)
  {
    // Dividing by base bin size in chi-squared because the weight is width/base bin size when filling
    if( fAdapDataHistoM1->GetBinCenter(i) >= 3150 && fAdapDataHistoM1->GetBinCenter(i) <= 3400)continue;

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
    if( fAdapDataHistoM2->GetBinCenter(i) >= 3150 && fAdapDataHistoM2->GetBinCenter(i) <= 3400)continue;

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
  // Add Gaussian and Poisson constraints on parameters here
  // Gaussian constraints work fine, however the Poisson limits don't seem to work (yet)
  // chiSquare += 0.5*(fParameters[28]*dDataIntegralM1/dLivetimeYr/fParEfficiencyM1[28]/fParMass[28]/(365*24*60*60) - 7.0E-4)*(fParameters[28]*dDataIntegralM1/dLivetimeYr/fParEfficiencyM1[28]/fParMass[28]/(365*24*60*60) - 7.0E-4)/2.0E-4;
  // chiSquare += 0.5*(fParameters[29]*dDataIntegralM1/dLivetimeYr/fParEfficiencyM1[29]/fParMass[28]/(365*24*60*60) - 5.0E-4)*(fParameters[29]*dDataIntegralM1/dLivetimeYr/fParEfficiencyM1[28]/fParMass[28]/(365*24*60*60) - 5.0E-4)/1.0E-4;

  //fParPrior[28] = 0.; // 7+/-2E-4 measured
  // fParPrior[29] = 0.; // 5+/-1E-5

  // for(int i = 0; i < dNParam; i++)
  // {
    // chiSquare += AddConstraint(26);
  // }


  return chiSquare;
}

void TBackgroundModel::Initialize()
{	
  // Loads PDFs from file
  cout << "Loading PDF Histograms from file" << endl;
  cout << "Directory " << Form("%s/", dMCDir.c_str()) << endl;

  fBulk = new TFile(Form("%s/MCProduction_Bulk_1keV.root", dMCDir.c_str()));
  fSurface = new TFile(Form("%s/MCProduction_Surface_1keV.root", dMCDir.c_str()));
  fBulk_CDR = new TFile(Form("%s/MCProduction_BulkCDR_1keV.root", dMCDir.c_str()));
  fBulk_CDRInternal = new TFile(Form("%s/MCProduction_BulkCDRInternal_1keV.root", dMCDir.c_str()));

///////////// Bulk Histograms
/////// Crystal M1 and M2
  hTeO20nuM1     = (TH1D*)fBulk->Get("hTeO20nuM1");
  hTeO22nuM1     = (TH1D*)fBulk->Get("hTeO22nuM1");
  hTeO22nuIKM1     = (TH1D*)fBulk->Get("hTeO22nuIKM1");
  hTeO2co60M1    = (TH1D*)fBulk->Get("hTeO2co60M1");
  hTeO2k40M1     = (TH1D*)fBulk->Get("hTeO2k40M1");
  hTeO2pb210M1   = (TH1D*)fBulk->Get("hTeO2pb210M1");
  hTeO2po210M1   = (TH1D*)fBulk->Get("hTeO2po210M1");
  hTeO2te125M1   = (TH1D*)fBulk->Get("hTeO2te125M1");
  hTeO2th232M1   = (TH1D*)fBulk->Get("hTeO2th232M1");
  hTeO2u238M1    = (TH1D*)fBulk->Get("hTeO2u238M1");
  hTeO2sb125M1   = (TH1D*)fBulk->Get("hTeO2sb125M1");

  hTeO2th232onlyM1 = (TH1D*)fBulk->Get("hTeO2th232onlyM1");
  hTeO2ra228pb208M1 = (TH1D*)fBulk->Get("hTeO2ra228pb208M1");
  hTeO2th230onlyM1 = (TH1D*)fBulk->Get("hTeO2th230onlyM1");
  hTeO2u238th230M1 = (TH1D*)fBulk->Get("hTeO2u238th230M1");
  hTeO2ra226pb210M1 = (TH1D*)fBulk->Get("hTeO2ra226pb210M1");

  hTeO20nuM2     = (TH1D*)fBulk->Get("hTeO20nuM2");
  hTeO22nuM2     = (TH1D*)fBulk->Get("hTeO22nuM2");
  hTeO22nuIKM2     = (TH1D*)fBulk->Get("hTeO22nuIKM2");
  hTeO2co60M2    = (TH1D*)fBulk->Get("hTeO2co60M2");
  hTeO2k40M2     = (TH1D*)fBulk->Get("hTeO2k40M2");
  hTeO2pb210M2   = (TH1D*)fBulk->Get("hTeO2pb210M2");
  hTeO2po210M2   = (TH1D*)fBulk->Get("hTeO2po210M2");
  hTeO2te125M2   = (TH1D*)fBulk->Get("hTeO2te125M2");
  hTeO2th232M2   = (TH1D*)fBulk->Get("hTeO2th232M2");
  hTeO2u238M2    = (TH1D*)fBulk->Get("hTeO2u238M2");
  hTeO2sb125M2   = (TH1D*)fBulk->Get("hTeO2sb125M2");

  hTeO2th232onlyM2 = (TH1D*)fBulk->Get("hTeO2th232onlyM2");
  hTeO2ra228pb208M2 = (TH1D*)fBulk->Get("hTeO2ra228pb208M2");
  hTeO2th230onlyM2 = (TH1D*)fBulk->Get("hTeO2th230onlyM2");
  hTeO2u238th230M2 = (TH1D*)fBulk->Get("hTeO2u238th230M2");
  hTeO2ra226pb210M2 = (TH1D*)fBulk->Get("hTeO2ra226pb210M2");

/*
///////// CuBox + CuFrame M1 and M2
  hCuBox_CuFrameco60M1 = (TH1D*)fBulk->Get("hCuBox_CuFrameco60M1");
  hCuBox_CuFramek40M1 = (TH1D*)fBulk->Get("hCuBox_CuFramek40M1");
  hCuBox_CuFrameth232M1 = (TH1D*)fBulk->Get("hCuBox_CuFrameth232M1");
  hCuBox_CuFrameu238M1 = (TH1D*)fBulk->Get("hCuBox_CuFrameu238M1");
  hCuBox_CuFramebi207M1 = (TH1D*)fBulk->Get("hCuBox_CuFramebi207M1");
  hCuBox_CuFramemn54M1 = (TH1D*)fBulk->Get("hCuBox_CuFramemn54M1");

  hCuBox_CuFrameco60M2 = (TH1D*)fBulk->Get("hCuBox_CuFrameco60M2");
  hCuBox_CuFramek40M2 = (TH1D*)fBulk->Get("hCuBox_CuFramek40M2");
  hCuBox_CuFrameth232M2 = (TH1D*)fBulk->Get("hCuBox_CuFrameth232M2");
  hCuBox_CuFrameu238M2 = (TH1D*)fBulk->Get("hCuBox_CuFrameu238M2");
  hCuBox_CuFramebi207M2 = (TH1D*)fBulk->Get("hCuBox_CuFramebi207M2");
  hCuBox_CuFramemn54M2 = (TH1D*)fBulk->Get("hCuBox_CuFramemn54M2");


///////// 50mK M1 and M2
  h50mKcs137M1   = (TH1D*)fBulk->Get("hInternalcs137M1");
  h50mKcs137M2   = (TH1D*)fBulk->Get("hInternalcs137M2");

//////// Roman Lead M1 and M2
  hPbRomk40M1     = (TH1D*)fBulk->Get("hPbRomk40M1");
  hPbRomk40M2     = (TH1D*)fBulk->Get("hPbRomk40M2");
  hPbRomco60M1    = (TH1D*)fBulk->Get("hPbRomco60M1");
  hPbRomth232M1   = (TH1D*)fBulk->Get("hPbRomth232M1");
  hPbRomu238M1    = (TH1D*)fBulk->Get("hPbRomu238M1");
  hPbRomco60M2    = (TH1D*)fBulk->Get("hPbRomco60M2");
  hPbRomth232M2   = (TH1D*)fBulk->Get("hPbRomth232M2");
  hPbRomu238M2    = (TH1D*)fBulk->Get("hPbRomu238M2");

////// Internal Shields M1 and M2
  hInternalco60M1 = (TH1D*)fBulk->Get("hInternalco60M1");
  hInternalk40M1 = (TH1D*)fBulk->Get("hInternalk40M1");
  hInternalth232M1 = (TH1D*)fBulk->Get("hInternalth232M1");
  hInternalu238M1 = (TH1D*)fBulk->Get("hInternalu238M1");

  hInternalco60M2 = (TH1D*)fBulk->Get("hInternalco60M2");
  hInternalk40M2 = (TH1D*)fBulk->Get("hInternalk40M2");
  hInternalth232M2 = (TH1D*)fBulk->Get("hInternalth232M2");
  hInternalu238M2 = (TH1D*)fBulk->Get("hInternalu238M2");

/////// OVC M1 and M2
  hOVCco60M1    = (TH1D*)fBulk->Get("hOVCco60M1");
  hOVCk40M1     = (TH1D*)fBulk->Get("hOVCk40M1");
  hOVCth232M1   = (TH1D*)fBulk->Get("hOVCth232M1");
  hOVCu238M1    = (TH1D*)fBulk->Get("hOVCu238M1");
  hOVCbi207M1   = (TH1D*)fBulk->Get("hOVCbi207M1");
   
  hOVCco60M2    = (TH1D*)fBulk->Get("hOVCco60M2");
  hOVCk40M2     = (TH1D*)fBulk->Get("hOVCk40M2");
  hOVCth232M2   = (TH1D*)fBulk->Get("hOVCth232M2");
  hOVCu238M2    = (TH1D*)fBulk->Get("hOVCu238M2");
  hOVCbi207M2   = (TH1D*)fBulk->Get("hOVCbi207M2");

/////// External Sources M1 and M2
  hExtPbbi210M1 = (TH1D*)fBulk->Get("hExtPbbi210M1");
  hExtPbk40M1 = (TH1D*)fBulk->Get("hExtPbk40M1");
  hExtPbth232M1 = (TH1D*)fBulk->Get("hExtPbth232M1");
  hExtPbu238M1 = (TH1D*)fBulk->Get("hExtPbu238M1");
  hExtPbpb210M1 = (TH1D*)fBulk->Get("hExtPbpb210M1");
 
  hExtPbbi210M2 = (TH1D*)fBulk->Get("hExtPbbi210M2");
  hExtPbk40M2 = (TH1D*)fBulk->Get("hExtPbk40M2");
  hExtPbth232M2 = (TH1D*)fBulk->Get("hExtPbth232M2");
  hExtPbu238M2 = (TH1D*)fBulk->Get("hExtPbu238M2");
  hExtPbpb210M2 = (TH1D*)fBulk->Get("hExtPbpb210M2");

  hExtMuonM1 = (TH1D*)fBulk->Get("hExtMuonM1");
  hCuBox_th232spotM1 = (TH1D*)fBulk->Get("hCuBox_th232spotM1");
  hCuBox_k40spotM1 = (TH1D*)fBulk->Get("hCuBox_k40spotM1");
  hBotExtPb_k40spotM1 = (TH1D*)fBulk->Get("hBotExtPb_k40spotM1");

  hExtMuonM2 = (TH1D*)fBulk->Get("hExtMuonM2");
  hCuBox_th232spotM2 = (TH1D*)fBulk->Get("hCuBox_th232spotM2");
  hCuBox_k40spotM2 = (TH1D*)fBulk->Get("hCuBox_k40spotM2");
  hBotExtPb_k40spotM2 = (TH1D*)fBulk->Get("hBotExtPb_k40spotM2");
*/

////////// CDR
///////// CuBox + CuFrame M1 and M2
  hCuBox_CuFrameco60M1 = (TH1D*)fBulk_CDR->Get("hCuBox_CuFrameco60M1");
  hCuBox_CuFramek40M1 = (TH1D*)fBulk_CDR->Get("hCuBox_CuFramek40M1");
  hCuBox_CuFrameth232M1 = (TH1D*)fBulk_CDR->Get("hCuBox_CuFrameth232M1");
  hCuBox_CuFrameu238M1 = (TH1D*)fBulk_CDR->Get("hCuBox_CuFrameu238M1");
  hCuBox_CuFramebi207M1 = (TH1D*)fBulk->Get("hCuBox_CuFramebi207M1");
  hCuBox_CuFramemn54M1 = (TH1D*)fBulk_CDR->Get("hCuBox_CuFramemn54M1");

  hCuBox_CuFrameco60M2 = (TH1D*)fBulk_CDR->Get("hCuBox_CuFrameco60M2");
  hCuBox_CuFramek40M2 = (TH1D*)fBulk_CDR->Get("hCuBox_CuFramek40M2");
  hCuBox_CuFrameth232M2 = (TH1D*)fBulk_CDR->Get("hCuBox_CuFrameth232M2");
  hCuBox_CuFrameu238M2 = (TH1D*)fBulk_CDR->Get("hCuBox_CuFrameu238M2");
  hCuBox_CuFramebi207M2 = (TH1D*)fBulk->Get("hCuBox_CuFramebi207M2");
  hCuBox_CuFramemn54M2 = (TH1D*)fBulk_CDR->Get("hCuBox_CuFramemn54M2");

///////// 50mK M1 and M2
  h50mKcs137M1   = (TH1D*)fBulk_CDRInternal->Get("hInternalcs137M1");
  h50mKcs137M2   = (TH1D*)fBulk_CDRInternal->Get("hInternalcs137M2");
////// Internal Shields M1 and M2
  hInternalco60M1 = (TH1D*)fBulk_CDRInternal->Get("hInternalco60M1");
  hInternalk40M1 = (TH1D*)fBulk_CDRInternal->Get("hInternalk40M1");
  hInternalth232M1 = (TH1D*)fBulk_CDRInternal->Get("hInternalth232M1");
  hInternalu238M1 = (TH1D*)fBulk_CDRInternal->Get("hInternalu238M1");

  hInternalco60M2 = (TH1D*)fBulk_CDRInternal->Get("hInternalco60M2");
  hInternalk40M2 = (TH1D*)fBulk_CDRInternal->Get("hInternalk40M2");
  hInternalth232M2 = (TH1D*)fBulk_CDRInternal->Get("hInternalth232M2");
  hInternalu238M2 = (TH1D*)fBulk_CDRInternal->Get("hInternalu238M2");

//////// Roman Lead M1 and M2
  hPbRomk40M1     = (TH1D*)fBulk_CDR->Get("hPbRomk40M1");
  hPbRomk40M2     = (TH1D*)fBulk_CDR->Get("hPbRomk40M2");
  hPbRomco60M1    = (TH1D*)fBulk_CDR->Get("hPbRomco60M1");
  hPbRomth232M1   = (TH1D*)fBulk_CDR->Get("hPbRomth232M1");
  hPbRomu238M1    = (TH1D*)fBulk_CDR->Get("hPbRomu238M1");
  hPbRomco60M2    = (TH1D*)fBulk_CDR->Get("hPbRomco60M2");
  hPbRomth232M2   = (TH1D*)fBulk_CDR->Get("hPbRomth232M2");
  hPbRomu238M2    = (TH1D*)fBulk_CDR->Get("hPbRomu238M2");

/////// OVC M1 and M2
  hOVCco60M1    = (TH1D*)fBulk_CDR->Get("hOVCco60M1");
  hOVCk40M1     = (TH1D*)fBulk_CDR->Get("hOVCk40M1");
  hOVCth232M1   = (TH1D*)fBulk_CDR->Get("hOVCth232M1");
  hOVCu238M1    = (TH1D*)fBulk_CDR->Get("hOVCu238M1");
  hOVCbi207M1   = (TH1D*)fBulk_CDR->Get("hOVCbi207M1");
   
  hOVCco60M2    = (TH1D*)fBulk_CDR->Get("hOVCco60M2");
  hOVCk40M2     = (TH1D*)fBulk_CDR->Get("hOVCk40M2");
  hOVCth232M2   = (TH1D*)fBulk_CDR->Get("hOVCth232M2");
  hOVCu238M2    = (TH1D*)fBulk_CDR->Get("hOVCu238M2");
  hOVCbi207M2   = (TH1D*)fBulk_CDR->Get("hOVCbi207M2");

/////// External Sources M1 and M2
  hExtPbbi210M1 = (TH1D*)fBulk_CDR->Get("hExtPbbi210M1");
  hExtPbk40M1 = (TH1D*)fBulk_CDR->Get("hExtPbk40M1");
  hExtPbth232M1 = (TH1D*)fBulk_CDR->Get("hExtPbth232M1");
  hExtPbu238M1 = (TH1D*)fBulk_CDR->Get("hExtPbu238M1");
  hExtPbpb210M1 = (TH1D*)fBulk_CDR->Get("hExtPbpb210M1");
 
  hExtPbbi210M2 = (TH1D*)fBulk_CDR->Get("hExtPbbi210M2");
  hExtPbk40M2 = (TH1D*)fBulk_CDR->Get("hExtPbk40M2");
  hExtPbth232M2 = (TH1D*)fBulk_CDR->Get("hExtPbth232M2");
  hExtPbu238M2 = (TH1D*)fBulk_CDR->Get("hExtPbu238M2");
  hExtPbpb210M2 = (TH1D*)fBulk_CDR->Get("hExtPbpb210M2");

  hExtMuonM1 = (TH1D*)fBulk->Get("hExtMuonM1");
  hCuBox_th232spotM1 = (TH1D*)fBulk->Get("hCuBox_th232spotM1");
  hCuBox_k40spotM1 = (TH1D*)fBulk->Get("hCuBox_k40spotM1");
  hBotExtPb_k40spotM1 = (TH1D*)fBulk_CDR->Get("hBotExtPb_k40spotM1");

  hExtMuonM2 = (TH1D*)fBulk->Get("hExtMuonM2");
  hCuBox_th232spotM2 = (TH1D*)fBulk->Get("hCuBox_th232spotM2");
  hCuBox_k40spotM2 = (TH1D*)fBulk->Get("hCuBox_k40spotM2");
  hBotExtPb_k40spotM2 = (TH1D*)fBulk_CDR->Get("hBotExtPb_k40spotM2");
  
//////////// Surface PDFs
///// Crystal M1 and M2
  hTeO2Sxpb210M1_001  = (TH1D*)fSurface->Get("hTeO2Sxpb210M1_001");
  hTeO2Sxpb210M1_01   = (TH1D*)fSurface->Get("hTeO2Sxpb210M1_01");
  hTeO2Sxpb210M1_1    = (TH1D*)fSurface->Get("hTeO2Sxpb210M1_1");
  hTeO2Sxpb210M1_10   = (TH1D*)fSurface->Get("hTeO2Sxpb210M1_10");
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

///////////// Get adaptive binned histograms
//////// Crystal M1 and M2
  hnewTeO20nuM1 = hTeO20nuM1->Rebin(dAdaptiveBinsM1, "hnewTeO20nuM1", dAdaptiveArrayM1);
  hnewTeO22nuM1 = hTeO22nuM1->Rebin(dAdaptiveBinsM1, "hnewTeO22nuM1", dAdaptiveArrayM1);
  hnewTeO22nuIKM1 = hTeO22nuIKM1->Rebin(dAdaptiveBinsM1, "hnewTeO22nuIKM1", dAdaptiveArrayM1);
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
  hnewTeO22nuIKM2 = hTeO22nuIKM2->Rebin(dAdaptiveBinsM2, "hnewTeO22nuIKM2", dAdaptiveArrayM2); 
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

///////// CuBox + CuFrame M1 and M2
  hnewCuBox_CuFrameco60M1 = hCuBox_CuFrameco60M1->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameco60M1", dAdaptiveArrayM1);
  hnewCuBox_CuFramek40M1 = hCuBox_CuFramek40M1->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFramek40M1", dAdaptiveArrayM1);
  hnewCuBox_CuFrameth232M1 = hCuBox_CuFrameth232M1->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameth232M1", dAdaptiveArrayM1);
  hnewCuBox_CuFrameu238M1 = hCuBox_CuFrameu238M1->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFrameu238M1", dAdaptiveArrayM1);
  hnewCuBox_CuFramemn54M1 = hCuBox_CuFramemn54M1->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFramemn54M1", dAdaptiveArrayM1);
  hnewCuBox_CuFramebi207M1 = hCuBox_CuFramebi207M1->Rebin(dAdaptiveBinsM1, "hnewCuBox_CuFramebi207M1", dAdaptiveArrayM1);

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
  hnewCuBox_CuFramemn54M2 = hCuBox_CuFramemn54M2->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFramemn54M2", dAdaptiveArrayM2);
  hnewCuBox_CuFramebi207M2 = hCuBox_CuFramebi207M2->Rebin(dAdaptiveBinsM2, "hnewCuBox_CuFramebi207M2", dAdaptiveArrayM2);

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

////////// 50mK M1 and M2
  hnew50mKcs137M1 = h50mKcs137M1->Rebin(dAdaptiveBinsM1, "hnew50mKcs137M1", dAdaptiveArrayM1);
  hnew50mKcs137M2 = h50mKcs137M2->Rebin(dAdaptiveBinsM2, "hnew50mKcs137M2", dAdaptiveArrayM2);
  
///////// Roman Lead M1 and M2
  hnewPbRomco60M1 = hPbRomco60M1->Rebin(dAdaptiveBinsM1, "hnewPbRomco60M1", dAdaptiveArrayM1);
  hnewPbRomk40M1 = hPbRomk40M1->Rebin(dAdaptiveBinsM1, "hnewPbRomk40M1", dAdaptiveArrayM1);
  hnewPbRomth232M1 = hPbRomth232M1->Rebin(dAdaptiveBinsM1, "hnewPbRomth232M1", dAdaptiveArrayM1);
  hnewPbRomu238M1 = hPbRomu238M1->Rebin(dAdaptiveBinsM1, "hnewPbRomu238M1", dAdaptiveArrayM1);

  hnewPbRomco60M2 = hPbRomco60M2->Rebin(dAdaptiveBinsM2, "hnewPbRomco60M2", dAdaptiveArrayM2);
  hnewPbRomk40M2 = hPbRomk40M2->Rebin(dAdaptiveBinsM2, "hnewPbRomk40M2", dAdaptiveArrayM2);
  hnewPbRomth232M2 = hPbRomth232M2->Rebin(dAdaptiveBinsM2, "hnewPbRomth232M2", dAdaptiveArrayM2);
  hnewPbRomu238M2 = hPbRomu238M2->Rebin(dAdaptiveBinsM2, "hnewPbRomu238M2", dAdaptiveArrayM2);

////////// Internal Shields
  hnewInternalco60M1 = hInternalco60M1->Rebin(dAdaptiveBinsM1, "hnewInternalco60M1", dAdaptiveArrayM1);
  hnewInternalk40M1 = hInternalk40M1->Rebin(dAdaptiveBinsM1, "hnewInternalk40M1", dAdaptiveArrayM1);
  hnewInternalth232M1 = hInternalth232M1->Rebin(dAdaptiveBinsM1, "hnewInternalth232M1", dAdaptiveArrayM1);
  hnewInternalu238M1 = hInternalu238M1->Rebin(dAdaptiveBinsM1, "hnewInternalu238M1", dAdaptiveArrayM1);

  hnewInternalco60M2 = hInternalco60M2->Rebin(dAdaptiveBinsM2, "hnewInternalco60M2", dAdaptiveArrayM2);
  hnewInternalk40M2 = hInternalk40M2->Rebin(dAdaptiveBinsM2, "hnewInternalk40M2", dAdaptiveArrayM2);
  hnewInternalth232M2 = hInternalth232M2->Rebin(dAdaptiveBinsM2, "hnewInternalth232M2", dAdaptiveArrayM2);
  hnewInternalu238M2 = hInternalu238M2->Rebin(dAdaptiveBinsM2, "hnewInternalu238M2", dAdaptiveArrayM2);

////////// OVC M1 and M2
  hnewOVCco60M1 = hOVCco60M1->Rebin(dAdaptiveBinsM1, "hnewOVCco60M1", dAdaptiveArrayM1);
  hnewOVCk40M1 = hOVCk40M1->Rebin(dAdaptiveBinsM1, "hnewOVCk40M1", dAdaptiveArrayM1);
  hnewOVCth232M1 = hOVCth232M1->Rebin(dAdaptiveBinsM1, "hnewOVCth232M1", dAdaptiveArrayM1);
  hnewOVCu238M1 = hOVCu238M1->Rebin(dAdaptiveBinsM1, "hnewOVCu238M1", dAdaptiveArrayM1);
  hnewOVCbi207M1 = hOVCbi207M1->Rebin(dAdaptiveBinsM1, "hnewOVCbi207M1", dAdaptiveArrayM1);

  hnewOVCco60M2 = hOVCco60M2->Rebin(dAdaptiveBinsM2, "hnewOVCco60M2", dAdaptiveArrayM2);
  hnewOVCk40M2 = hOVCk40M2->Rebin(dAdaptiveBinsM2, "hnewOVCk40M2", dAdaptiveArrayM2);
  hnewOVCth232M2 = hOVCth232M2->Rebin(dAdaptiveBinsM2, "hnewOVCth232M2", dAdaptiveArrayM2);
  hnewOVCu238M2 = hOVCu238M2->Rebin(dAdaptiveBinsM2, "hnewOVCu238M2", dAdaptiveArrayM2);
  hnewOVCbi207M2 = hOVCbi207M2->Rebin(dAdaptiveBinsM2, "hnewOVCbi207M2", dAdaptiveArrayM2);

  hnewExtPbbi210M1 = hExtPbbi210M1->Rebin(dAdaptiveBinsM1, "hnewExtPbbi210M1", dAdaptiveArrayM1);
  hnewExtPbk40M1 = hExtPbk40M1->Rebin(dAdaptiveBinsM1, "hnewExtPbk40M1", dAdaptiveArrayM1);
  hnewExtPbth232M1 = hExtPbth232M1->Rebin(dAdaptiveBinsM1, "hnewExtPbth232M1", dAdaptiveArrayM1);
  hnewExtPbu238M1 = hExtPbu238M1->Rebin(dAdaptiveBinsM1, "hnewExtPbu238M1", dAdaptiveArrayM1);
  hnewExtPbpb210M1 = hExtPbpb210M1->Rebin(dAdaptiveBinsM1, "hnewExtPbpb210M1", dAdaptiveArrayM1);

  hnewExtPbbi210M2 = hExtPbbi210M2->Rebin(dAdaptiveBinsM2, "hnewExtPbbi210M2", dAdaptiveArrayM2);
  hnewExtPbk40M2 = hExtPbk40M2->Rebin(dAdaptiveBinsM2, "hnewExtPbk40M2", dAdaptiveArrayM2);
  hnewExtPbth232M2 = hExtPbth232M2->Rebin(dAdaptiveBinsM2, "hnewExtPbth232M2", dAdaptiveArrayM2);
  hnewExtPbu238M2 = hExtPbu238M2->Rebin(dAdaptiveBinsM2, "hnewExtPbu238M2", dAdaptiveArrayM2);
  hnewExtPbpb210M2 = hExtPbpb210M2->Rebin(dAdaptiveBinsM2, "hnewExtPbpb210M2", dAdaptiveArrayM2);

  hnewExtMuonM1 = hExtMuonM1->Rebin(dAdaptiveBinsM1, "hnewExtMuonM1", dAdaptiveArrayM1);
  hnewCuBox_th232spotM1 = hCuBox_th232spotM1->Rebin(dAdaptiveBinsM1, "hnewCuBox_th232spotM1", dAdaptiveArrayM1);
  hnewCuBox_k40spotM1 = hCuBox_k40spotM1->Rebin(dAdaptiveBinsM1, "hnewCuBox_k40spotM1", dAdaptiveArrayM1);
  hnewBotExtPb_k40spotM1 = hBotExtPb_k40spotM1->Rebin(dAdaptiveBinsM1, "hnewBotExtPb_k40spotM1", dAdaptiveArrayM1);

  hnewExtMuonM2 = hExtMuonM2->Rebin(dAdaptiveBinsM2, "hnewExtMuonM2", dAdaptiveArrayM2);
  hnewCuBox_th232spotM2 = hCuBox_th232spotM2->Rebin(dAdaptiveBinsM2, "hnewCuBox_th232spotM2", dAdaptiveArrayM2);
  hnewCuBox_k40spotM2 = hCuBox_k40spotM2->Rebin(dAdaptiveBinsM2, "hnewCuBox_k40spotM2", dAdaptiveArrayM2);
  hnewBotExtPb_k40spotM2 = hBotExtPb_k40spotM2->Rebin(dAdaptiveBinsM2, "hnewBotExtPb_k40spotM2", dAdaptiveArrayM2);

  // Fill adaptive binned histograms
  for(int i = 1; i <= dAdaptiveBinsM1; i++)
  {
    hAdapTeO20nuM1->SetBinContent(i, hnewTeO20nuM1->GetBinContent(i)/hnewTeO20nuM1->GetBinWidth(i));
    hAdapTeO22nuM1->SetBinContent(i, hnewTeO22nuM1->GetBinContent(i)/hnewTeO22nuM1->GetBinWidth(i));
    hAdapTeO22nuIKM1->SetBinContent(i, hnewTeO22nuIKM1->GetBinContent(i)/hnewTeO22nuIKM1->GetBinWidth(i));
    hAdapTeO2co60M1->SetBinContent(i, hnewTeO2co60M1->GetBinContent(i)/hnewTeO2co60M1->GetBinWidth(i));
    hAdapTeO2k40M1->SetBinContent(i, hnewTeO2k40M1->GetBinContent(i)/hnewTeO2k40M1->GetBinWidth(i));
    hAdapTeO2pb210M1->SetBinContent(i, hnewTeO2pb210M1->GetBinContent(i)/hnewTeO2pb210M1->GetBinWidth(i));
    hAdapTeO2po210M1->SetBinContent(i, hnewTeO2po210M1->GetBinContent(i)/hnewTeO2po210M1->GetBinWidth(i));
    hAdapTeO2te125M1->SetBinContent(i, hnewTeO2te125M1->GetBinContent(i)/hnewTeO2te125M1->GetBinWidth(i));
    hAdapTeO2th232M1->SetBinContent(i, hnewTeO2th232M1->GetBinContent(i)/hnewTeO2th232M1->GetBinWidth(i));
    hAdapTeO2u238M1->SetBinContent(i, hnewTeO2u238M1->GetBinContent(i)/hnewTeO2u238M1->GetBinWidth(i));
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
    hAdapCuBox_CuFramemn54M1->SetBinContent(i, hnewCuBox_CuFramemn54M1->GetBinContent(i)/hnewCuBox_CuFramemn54M1->GetBinWidth(i));
    hAdapCuBox_CuFramebi207M1->SetBinContent(i, hnewCuBox_CuFramebi207M1->GetBinContent(i)/hnewCuBox_CuFramebi207M1->GetBinWidth(i));

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
    hAdapPbRomk40M1->SetBinContent(i, hnewPbRomk40M1->GetBinContent(i)/hnewPbRomk40M1->GetBinWidth(i));
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
    hAdapOVCbi207M1->SetBinContent(i, hnewOVCbi207M1->GetBinContent(i)/hnewOVCbi207M1->GetBinWidth(i));

    hAdapExtPbbi210M1->SetBinContent(i, hnewExtPbbi210M1->GetBinContent(i)/hnewExtPbbi210M1->GetBinWidth(i));
    hAdapExtPbk40M1->SetBinContent(i, hnewExtPbk40M1->GetBinContent(i)/hnewExtPbk40M1->GetBinWidth(i));
    hAdapExtPbth232M1->SetBinContent(i, hnewExtPbth232M1->GetBinContent(i)/hnewExtPbth232M1->GetBinWidth(i));
    hAdapExtPbu238M1->SetBinContent(i, hnewExtPbu238M1->GetBinContent(i)/hnewExtPbu238M1->GetBinWidth(i));
    hAdapExtPbpb210M1->SetBinContent(i, hnewExtPbpb210M1->GetBinContent(i)/hnewExtPbpb210M1->GetBinWidth(i));

    hAdapExtMuonM1->SetBinContent(i, hnewExtMuonM1->GetBinContent(i)/hnewExtMuonM1->GetBinWidth(i));
    hAdapCuBox_th232spotM1->SetBinContent(i, hnewCuBox_th232spotM1->GetBinContent(i)/hnewCuBox_th232spotM1->GetBinWidth(i));
    hAdapCuBox_k40spotM1->SetBinContent(i, hnewCuBox_k40spotM1->GetBinContent(i)/hnewCuBox_k40spotM1->GetBinWidth(i));
    hAdapBotExtPb_k40spotM1->SetBinContent(i, hnewBotExtPb_k40spotM1->GetBinContent(i)/hnewBotExtPb_k40spotM1->GetBinWidth(i));

  }

  for(int i = 1; i <= dAdaptiveBinsM2; i++)
  {
    hAdapTeO20nuM2->SetBinContent(i, hnewTeO20nuM2->GetBinContent(i)/hnewTeO20nuM2->GetBinWidth(i));
    hAdapTeO22nuM2->SetBinContent(i, hnewTeO22nuM2->GetBinContent(i)/hnewTeO22nuM2->GetBinWidth(i));
    hAdapTeO22nuIKM2->SetBinContent(i, hnewTeO22nuIKM2->GetBinContent(i)/hnewTeO22nuIKM2->GetBinWidth(i));
    hAdapTeO2co60M2->SetBinContent(i, hnewTeO2co60M2->GetBinContent(i)/hnewTeO2co60M2->GetBinWidth(i));
    hAdapTeO2k40M2->SetBinContent(i, hnewTeO2k40M2->GetBinContent(i)/hnewTeO2k40M2->GetBinWidth(i));
    hAdapTeO2pb210M2->SetBinContent(i, hnewTeO2pb210M2->GetBinContent(i)/hnewTeO2pb210M2->GetBinWidth(i));
    hAdapTeO2po210M2->SetBinContent(i, hnewTeO2po210M2->GetBinContent(i)/hnewTeO2po210M2->GetBinWidth(i));
    hAdapTeO2te125M2->SetBinContent(i, hnewTeO2te125M2->GetBinContent(i)/hnewTeO2te125M2->GetBinWidth(i));
    hAdapTeO2th232M2->SetBinContent(i, hnewTeO2th232M2->GetBinContent(i)/hnewTeO2th232M2->GetBinWidth(i));
    hAdapTeO2u238M2->SetBinContent(i, hnewTeO2u238M2->GetBinContent(i)/hnewTeO2u238M2->GetBinWidth(i));
    hAdapTeO2sb125M2->SetBinContent(i, hnewTeO2sb125M2->GetBinContent(i)/hnewTeO2sb125M2->GetBinWidth(i));

    hAdapTeO2Sxpb210M2_001->SetBinContent(i, hnewTeO2Sxpb210M2_001->GetBinContent(i)/hnewTeO2Sxpb210M2_001->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_01->SetBinContent(i, hnewTeO2Sxpb210M2_01->GetBinContent(i)/hnewTeO2Sxpb210M2_01->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_1->SetBinContent(i, hnewTeO2Sxpb210M2_1->GetBinContent(i)/hnewTeO2Sxpb210M2_1->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_10->SetBinContent(i, hnewTeO2Sxpb210M2_10->GetBinContent(i)/hnewTeO2Sxpb210M2_10->GetBinWidth(i));
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
    hAdapCuBox_CuFramemn54M2->SetBinContent(i, hnewCuBox_CuFramemn54M2->GetBinContent(i)/hnewCuBox_CuFramemn54M2->GetBinWidth(i));
    hAdapCuBox_CuFramebi207M2->SetBinContent(i, hnewCuBox_CuFramebi207M2->GetBinContent(i)/hnewCuBox_CuFramebi207M2->GetBinWidth(i));

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
    hAdapPbRomk40M2->SetBinContent(i, hnewPbRomk40M2->GetBinContent(i)/hnewPbRomk40M2->GetBinWidth(i));
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
    hAdapOVCbi207M2->SetBinContent(i, hnewOVCbi207M2->GetBinContent(i)/hnewOVCbi207M2->GetBinWidth(i));

    hAdapExtPbbi210M2->SetBinContent(i, hnewExtPbbi210M2->GetBinContent(i)/hnewExtPbbi210M2->GetBinWidth(i));
    hAdapExtPbk40M2->SetBinContent(i, hnewExtPbk40M2->GetBinContent(i)/hnewExtPbk40M2->GetBinWidth(i));
    hAdapExtPbth232M2->SetBinContent(i, hnewExtPbth232M2->GetBinContent(i)/hnewExtPbth232M2->GetBinWidth(i));
    hAdapExtPbu238M2->SetBinContent(i, hnewExtPbu238M2->GetBinContent(i)/hnewExtPbu238M2->GetBinWidth(i));
    hAdapExtPbpb210M2->SetBinContent(i, hnewExtPbpb210M2->GetBinContent(i)/hnewExtPbpb210M2->GetBinWidth(i));

    hAdapExtMuonM2->SetBinContent(i, hnewExtMuonM2->GetBinContent(i)/hnewExtMuonM2->GetBinWidth(i));
    hAdapCuBox_th232spotM2->SetBinContent(i, hnewCuBox_th232spotM2->GetBinContent(i)/hnewCuBox_th232spotM2->GetBinWidth(i));
    hAdapCuBox_k40spotM2->SetBinContent(i, hnewCuBox_k40spotM2->GetBinContent(i)/hnewCuBox_k40spotM2->GetBinWidth(i));
    hAdapBotExtPb_k40spotM2->SetBinContent(i, hnewBotExtPb_k40spotM2->GetBinContent(i)/hnewBotExtPb_k40spotM2->GetBinWidth(i));

  }
}

// Load background data
// Data can be loaded in 
void TBackgroundModel::LoadData()
{
  // Set cuts
  base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
  base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
  // base_cut = base_cut && "Channel != 1 && Channel != 10"; 

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

  switch(dDataSet)
  { 
  // Livetimes only set for default because I'm lazy!
  case 1:
    // qtree->Add(Form("%s/Unblinded/ReducedB-ds2049.root", dDataDir.c_str()));   
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2061.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2064.root", dDataDir.c_str()));   
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2067.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2070.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2073.root", dDataDir.c_str())); 
    cout << "Using Data Release 1" << endl;
  break;

  case 2:
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2079.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2085.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2088.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2091.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2097.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2100.root", dDataDir.c_str())); 
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
    cout << "Using Data Release 3" << endl;
  break;

  // Combine DR 2 and 3 (the good data releases)
  case 4:
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2079.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2085.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2088.root", dDataDir.c_str())); 
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
    cout << "Using Data Release 2+3" << endl;
  break;

  default:   
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2061.root", dDataDir.c_str())); // Use this or no?
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2064.root", dDataDir.c_str()));   
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2067.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2070.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2073.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2076.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2079.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2085.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2088.root", dDataDir.c_str())); 
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
    cout << "Using Total Dataset" << endl;
  }

  qtree->Project("fDataHistoTot", "Energy", base_cut);
  qtree->Project("fDataHistoM1",  "Energy", base_cut && "Multiplicity_Sync == 1");
  qtree->Project("fDataHistoM2",  "Energy", base_cut && "Multiplicity_Sync == 2");

  dLivetimeYr = 0.8738; 
  dExposure = 33.4229;

  // Not using 2061, 2076, and 2088 for Radon
  // dExposure = 29.8967;
  // dLivetimeYr = 0.781613;
}

// Prints parameters, needs update 11-06-2014
void TBackgroundModel::PrintParameters()
{
  for(int i = 0; i < dNParam; i++)
  {
    cout << i << " " << minuit->fCpnam[i] << ": " << fParameters[i] << " +/- " << fParError[i] << endl;
  }
}

void TBackgroundModel::PrintParActivity()
{ 
  if(!bSave)
  {
    for(int i = 0; i < dNParam; i++)
    {
      BkgPar[i]->GetHistM1(1)->Scale( dDataIntegralM1*fParameters[i]);
      BkgPar[i]->GetHistM2(1)->Scale( dDataIntegralM2*fParameters[i]);
      fParCountsM1[i] = BkgPar[i]->GetHistM1(1)->Integral("width")/dLivetimeYr;
      fParCountsM2[i] = BkgPar[i]->GetHistM2(1)->Integral("width")/dLivetimeYr;    
      // Calculate Activity in terms of Bq/Kg
      fParActivityM1[i] = fParCountsM1[i]/fParEfficiencyM1[i]/fParMass[i]/(365*24*60*60);
    }
  }
  cout << "Counts/Year (M1) --- Error --- Activity (Bq/Kg)" << endl;
  for(int i = 0; i < dNParam; i++)
  {
    if(!bFixedArray[i])
    {
      cout << i << " " << fParCountsM1[i] << " +/- " << fParCountsM1[i]*fParError[i]/fParameters[i] << "\t" << fParActivityM1[i] << endl;
    }
  }
}

// Initializes the normalization values, also links the histograms to a variable index
void TBackgroundModel::GenerateParameters()
{
  // Initialization (Name, Index, Initial Value, Min Limit, Max Limit, pointer to histograms.. )
  // TeO2 Bulk
  // BkgPar[0] = new TBkgModelParameter( "TeO2 2$\\nu\\beta\\beta$", 0, 0, 1E-7, 0, 1.0, hAdapTeO22nuM1, hAdapTeO22nuM2 );
  BkgPar[0] = new TBkgModelParameter( "TeO2 2$\\nu\\beta\\beta$", 0, 0., 1E-7, 0, 1.0, hAdapTeO22nuIKM1, hAdapTeO22nuIKM2, hTeO22nuIKM1, hTeO22nuIKM2 );
  BkgPar[1] = new TBkgModelParameter( "TeO2 K40", 1, 0., 1E-7, 0, 1.0, hAdapTeO2k40M1, hAdapTeO2k40M2, hTeO2k40M1, hTeO2k40M2 ); 
  BkgPar[2] = new TBkgModelParameter( "TeO2 Co60", 2, 0., 1E-7, 0, 1.0, hAdapTeO2co60M1 , hAdapTeO2co60M2, hTeO2co60M1, hTeO2co60M2 );  
  BkgPar[3] = new TBkgModelParameter( "TeO2 Sb125", 3, 0, 1E-7, 0, 1.0, hAdapTeO2sb125M1 , hAdapTeO2sb125M2, hTeO2sb125M1, hTeO2sb125M2 ); 
  BkgPar[4] = new TBkgModelParameter( "TeO2 Th232 only", 4, 0., 1E-7, 0, 1.0, hAdapTeO2th232onlyM1 , hAdapTeO2th232onlyM2, hTeO2th232onlyM1, hTeO2th232onlyM2 ); 
  BkgPar[5] = new TBkgModelParameter( "TeO2 Ra228-Pb208", 5, 5.95362e-05, 1E-7, 0, 1.0, hAdapTeO2ra228pb208M1 , hAdapTeO2ra228pb208M2, hTeO2ra228pb208M1, hTeO2ra228pb208M2 ); 
  BkgPar[6] = new TBkgModelParameter( "TeO2 U238-Th230 ", 6, 0., 1E-7, 0, 1.0, hAdapTeO2u238th230M1 , hAdapTeO2u238th230M2, hTeO2u238th230M1, hTeO2u238th230M2 ); 
  BkgPar[7] = new TBkgModelParameter( "TeO2 Th230 only", 7, 5.49412e-04, 1E-7, 0, 1.0, hAdapTeO2th230onlyM1 , hAdapTeO2th230onlyM2, hTeO2th230onlyM1, hTeO2th230onlyM2 ); 
  BkgPar[8] = new TBkgModelParameter( "TeO2 Ra226-Pb210", 8, 4.39777e-04, 1E-7, 0, 1.0, hAdapTeO2ra226pb210M1 , hAdapTeO2ra226pb210M2, hTeO2ra226pb210M1, hTeO2ra226pb210M2 ); 

  // TeO2 Surface
  BkgPar[9] = new TBkgModelParameter( "TeO2 Sx Th232 only 0.001 $\\mu$m", 9, 8.31726e-05, 1E-7, 0, 1.0, hAdapTeO2Sxth232onlyM1_0001 , hAdapTeO2Sxth232onlyM2_0001, hTeO2Sxth232onlyM1_0001, hTeO2Sxth232onlyM2_0001 );
  BkgPar[10] = new TBkgModelParameter( "TeO2 Sx Ra228-Pb208 0.001 $\\mu$m", 10, 2.47160e-03, 1E-7, 0, 1.0, hAdapTeO2Sxra228pb208M1_0001 , hAdapTeO2Sxra228pb208M2_0001, hTeO2Sxra228pb208M1_0001, hTeO2Sxra228pb208M2_0001 );
  BkgPar[11] = new TBkgModelParameter( "TeO2 Sx U238-Th230 0.001 $\\mu$m", 11, 1.83844e-03, 1E-7, 0, 1.0, hAdapTeO2Sxu238th230M1_0001 , hAdapTeO2Sxu238th230M2_0001, hTeO2Sxu238th230M1_0001, hTeO2Sxu238th230M2_0001 );
  BkgPar[12] = new TBkgModelParameter( "TeO2 Sx Th230 only 0.001 $\\mu$m", 12, 3.12111e-04, 1E-7, 0, 1.0, hAdapTeO2Sxth230onlyM1_0001 , hAdapTeO2Sxth230onlyM2_0001, hTeO2Sxth230onlyM1_0001, hTeO2Sxth230onlyM2_0001 );
  BkgPar[13] = new TBkgModelParameter( "TeO2 Sx Ra226-pb210 0.001 $\\mu$m", 13, 3.11301e-03, 1E-7, 0, 1.0, hAdapTeO2Sxra226pb210M1_0001 , hAdapTeO2Sxra226pb210M2_0001, hTeO2Sxra226pb210M1_0001, hTeO2Sxra226pb210M2_0001  );
  BkgPar[14] = new TBkgModelParameter( "TeO2 Sx Pb210 0.01 $\\mu$m", 14, 4.28515e-02, 1E-7, 0, 1.0, hAdapTeO2Sxpb210M1_001 , hAdapTeO2Sxpb210M2_001, hTeO2Sxpb210M1_001, hTeO2Sxpb210M2_001 ); // Gives the peak at low energy
  BkgPar[15] = new TBkgModelParameter( "TeO2 Sx Pb210 1 $\\mu$m", 15, 4.23544e-03, 1E-7, 0, 1.0, hAdapTeO2Sxpb210M1_1 , hAdapTeO2Sxpb210M2_1, hTeO2Sxpb210M1_1, hTeO2Sxpb210M2_1 ); 
  BkgPar[16] = new TBkgModelParameter( "TeO2 Sx Th232 1", 16, 1.65735e-03, 1E-7, 0, 1.0, hAdapTeO2Sxth232M1_1 , hAdapTeO2Sxth232M2_1, hTeO2Sxth232M1_1, hTeO2Sxth232M2_1 );
  BkgPar[17] = new TBkgModelParameter( "TeO2 Sx Th232 10 $\\mu$m", 17, 1.29262e-03, 1E-7, 0, 1.0, hAdapTeO2Sxth232M1_10 , hAdapTeO2Sxth232M2_10, hTeO2Sxth232M1_10, hTeO2Sxth232M2_10 );   
  BkgPar[18] = new TBkgModelParameter( "Copper Holder Sx U238 100 $\\mu$m", 18, 4.70819e-03, 1E-7, 0, 1.0, hAdapCuBox_CuFrameu238M1_100 , hAdapCuBox_CuFrameu238M2_100, hCuBox_CuFrameu238M1_100, hCuBox_CuFrameu238M2_100 );
  BkgPar[19] = new TBkgModelParameter( "Copper Holder Sx Th232 100 $\\mu$m", 19, 2.18036e-02, 1E-7, 0, 1.0, hAdapCuBox_CuFrameth232M1_100 , hAdapCuBox_CuFrameth232M2_100, hCuBox_CuFrameth232M1_100, hCuBox_CuFrameth232M2_100 );
  BkgPar[20] = new TBkgModelParameter( "Copper Holder Sx U238 10 $\\mu$m", 20, 1.09059e-02, 1E-7, 0, 1.0, hAdapCuBox_CuFrameu238M1_10 , hAdapCuBox_CuFrameu238M2_10, hCuBox_CuFrameu238M1_10, hCuBox_CuFrameu238M2_10 );
  BkgPar[21] = new TBkgModelParameter( "Copper Holder Sx Pb210 1 $\\mu$m", 21, 1.62080e-03, 1E-7, 0, 1.0, hAdapCuBox_CuFramepb210M1_1 , hAdapCuBox_CuFramepb210M2_1, hCuBox_CuFramepb210M1_1, hCuBox_CuFramepb210M2_1 );
  BkgPar[22] = new TBkgModelParameter( "Copper Holder Sx Pb210 0.1 $\\mu$m", 22, 5.11887e-03, 1E-7, 0, 1.0, hAdapCuBox_CuFramepb210M1_01 , hAdapCuBox_CuFramepb210M2_01, hCuBox_CuFramepb210M1_01, hCuBox_CuFramepb210M2_01 );
  BkgPar[23] = new TBkgModelParameter( "Copper Holder Sx Pb210 0.01 $\\mu$m", 23, 1.93437e-02, 1E-7, 0, 1.0, hAdapCuBox_CuFramepb210M1_001 , hAdapCuBox_CuFramepb210M2_001, hCuBox_CuFramepb210M1_001, hCuBox_CuFramepb210M2_001 );
  BkgPar[24] = new TBkgModelParameter( "Copper Holder Sx Th232 0.01 $\\mu$m", 24, 4.67469e-04, 1E-7, 0, 1.0, hAdapCuBox_CuFrameth232M1_001 , hAdapCuBox_CuFrameth232M2_001, hCuBox_CuFrameth232M1_001, hCuBox_CuFrameth232M2_001 );

  // Bulk Contamination
  BkgPar[25] = new TBkgModelParameter( "Copper Holder Mn54", 25, 1.20782e-03, 1E-7, 0, 1.0, hAdapCuBox_CuFramemn54M1 , hAdapCuBox_CuFramemn54M2, hCuBox_CuFramemn54M1, hCuBox_CuFramemn54M2 ); // Mn54 (Cosmogenic activation of copper)
  BkgPar[26] = new TBkgModelParameter( "Copper Holder U238", 26, 0., 1E-7, 0, 1.0, hAdapCuBox_CuFrameu238M1 , hAdapCuBox_CuFrameu238M2, hCuBox_CuFrameu238M1, hCuBox_CuFrameu238M2 );
  BkgPar[27] = new TBkgModelParameter( "Copper Holder Th232", 27, 0., 1E-7, 0., 1.0, hAdapCuBox_CuFrameth232M1 , hAdapCuBox_CuFrameth232M2, hCuBox_CuFrameth232M1, hCuBox_CuFrameth232M2 );
  BkgPar[28] = new TBkgModelParameter( "Copper Holder K40", 28, 0., 1E-7, 0, 1.0, hAdapCuBox_CuFramek40M1 , hAdapCuBox_CuFramek40M2, hCuBox_CuFramek40M1, hCuBox_CuFramek40M2 );
  BkgPar[29] = new TBkgModelParameter( "Copper Holder Co60", 29, 0., 1E-7, 0, 1.0, hAdapCuBox_CuFrameco60M1 , hAdapCuBox_CuFrameco60M2, hCuBox_CuFrameco60M1, hCuBox_CuFrameco60M2 );
  BkgPar[30] = new TBkgModelParameter( "Internal Shields Cs137", 30, 1.88141e-03, 1E-7, 0, 1.0, hAdap50mKcs137M1 , hAdap50mKcs137M2, h50mKcs137M1, h50mKcs137M2 );
  BkgPar[31] = new TBkgModelParameter( "Internal Shields U238", 31, 0., 1E-7, 0, 1.0, hAdapInternalu238M1 , hAdapInternalu238M2, hInternalu238M1, hInternalu238M2 );
  BkgPar[32] = new TBkgModelParameter( "Internal Shields Th232", 32, 0., 1E-7, 0, 1.0, hAdapInternalth232M1 , hAdapInternalth232M2, hInternalth232M1, hInternalth232M2 );
  BkgPar[33] = new TBkgModelParameter( "Internal Shields K40", 33, 0., 1E-7, 0, 1.0, hAdapInternalk40M1 , hAdapInternalk40M2, hInternalk40M1, hInternalk40M2 );
  BkgPar[34] = new TBkgModelParameter( "Internal Shields Co60", 34, 0., 1E-7, 0, 1.0, hAdapInternalco60M1 , hAdapInternalco60M2, hInternalco60M1, hInternalco60M2 );
  BkgPar[35] = new TBkgModelParameter( "Roman Lead U238", 35, 0., 1E-7, 0, 1.0, hAdapPbRomu238M1 , hAdapPbRomu238M2, hPbRomu238M1, hPbRomu238M2 );
  BkgPar[36] = new TBkgModelParameter( "Roman Lead Th232", 36, 0., 1E-7, 0, 1.0, hAdapPbRomth232M1 , hAdapPbRomth232M2, hPbRomth232M1, hPbRomth232M2 );
  BkgPar[37] = new TBkgModelParameter( "Roman Lead K40", 37, 0., 1E-7, 0, 1.0, hAdapPbRomk40M1 , hAdapPbRomk40M2, hPbRomk40M1, hPbRomk40M2 );
  BkgPar[38] = new TBkgModelParameter( "Roman Lead Co60", 38, 0., 1E-7, 0, 1.0, hAdapPbRomco60M1 , hAdapPbRomco60M2, hPbRomco60M1, hPbRomco60M2 );
  BkgPar[39] = new TBkgModelParameter( "OVC U238", 39, 7.19053e-02, 1E-7, 0, 1.0, hAdapOVCu238M1 , hAdapOVCu238M2, hOVCu238M1, hOVCu238M2 );
  BkgPar[40] = new TBkgModelParameter( "OVC Th232", 40, 5.61000e-02, 1E-7, 0, 1.0, hAdapOVCth232M1 , hAdapOVCth232M2, hOVCth232M1, hOVCth232M2 );
  BkgPar[41] = new TBkgModelParameter( "OVC K40", 41, 5.41704e-02, 1E-7, 0, 1.0, hAdapOVCk40M1 , hAdapOVCk40M2, hOVCk40M1, hOVCk40M2 );
  BkgPar[42] = new TBkgModelParameter( "OVC Co60", 42, 2.49256e-02, 1E-7, 0, 1.0, hAdapOVCco60M1 , hAdapOVCco60M2, hOVCco60M1, hOVCco60M2 );
  BkgPar[43] = new TBkgModelParameter( "OVC Bi207", 43, 6.02468e-03, 1E-7, 0, 1.0, hAdapOVCbi207M1 , hAdapOVCbi207M2, hOVCbi207M1, hOVCbi207M2 );
  BkgPar[44] = new TBkgModelParameter( "External Lead Bi210", 44, 1.03086e-01, 1E-7, 0, 1.0, hAdapExtPbbi210M1 , hAdapExtPbbi210M2, hExtPbbi210M1, hExtPbbi210M2 );
  BkgPar[45] = new TBkgModelParameter( "External Lead K40", 45, 0, 1E-7, 0, 1.0, hAdapExtPbk40M1 , hAdapExtPbk40M2, hExtPbk40M1, hExtPbk40M2 );
  BkgPar[46] = new TBkgModelParameter( "External Lead Th232", 46, 2.61106e-02, 1E-7, 0, 1.0, hAdapExtPbth232M1 , hAdapExtPbth232M2, hExtPbth232M1, hExtPbth232M2 );
  BkgPar[47] = new TBkgModelParameter( "External Lead U238", 47, 3.93821e-02, 1E-7, 0, 1.0, hAdapExtPbu238M1 , hAdapExtPbu238M2, hExtPbu238M1, hExtPbu238M2 );
  BkgPar[48] = new TBkgModelParameter( "External Lead Pb210", 48, 3.96617e-03, 1E-7, 0, 1.0, hAdapExtPbpb210M1 , hAdapExtPbpb210M2, hExtPbpb210M1, hExtPbpb210M2 );
  BkgPar[49] = new TBkgModelParameter( "Bottom External Lead K40", 49, 3.69245e-02, 1E-7, 0, 1.0, hAdapBotExtPb_k40spotM1 , hAdapBotExtPb_k40spotM2, hBotExtPb_k40spotM1, hBotExtPb_k40spotM2 ); //
  

  //////////////
  // Don't fix parameter here when doing Profile Likelihood (except for other parameters that need fixing)
  // The Profile Likelihood function automatically fixes the parameter!
  // This section is just for fixing parameters for debugging only
  //////////////
  // bFixedArray[0] = true;
  // Fixing alphas
  dStartAlpha = 4; // These are just for loops
  dEndAlpha = 25;
  // for(int i = dStartAlpha; i < dEndAlpha; i ++)
  // {
  //   bFixedArray[i] = true;
  // }
}

// Resets all parameters to 0
// This method only resets the combined histograms, not the individual PDFs!
void TBackgroundModel::ResetParameters()
{
  dNumCalls = 0;
  dChiSquare = 0;
  
  fModelTotM1->Reset();
  fModelTotAdapM1->Reset();

  // Total PDFs M2
  fModelTotM2->Reset();
  fModelTotAdapM2->Reset();

  for(int i = 0; i < dNParam; i++)
  {
    fParameters[i] = 0;
    fParError[i] = 0;
  }
}

// Set Parameters in Model
void TBackgroundModel::SetParValue(int index, double value)
{
	// Change the index max depending on model
	if(index > dNParam) cout << "Index too large" << endl;
	else fParameters[index] = value;
}

// Creates/updates the background model
void TBackgroundModel::UpdateModel()
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
  if(dNumCalls%5000==0)
  {
    cout << "Call #: "<< dNumCalls << endl;
  }
  // Create model
  for(int i = 0; i < dNParam; i++)
  {
    // if(fParameters[i] <= 0 )continue;
    fModelTotAdapM1->Add( BkgPar[i]->GetHistM1(1),              dDataIntegralM1*fParameters[i]);
    fModelTotAdapM2->Add( BkgPar[i]->GetHistM2(1),              dDataIntegralM1*fParameters[i]);
  }
  // Add in Bulk Po-210 by hand
  // 8.5885e03 counts/year
  // double dNormPo210 = (8.5885e3)/dDataIntegralM1;
  // 1.77E3 counts assigned to Po-210 component
  fModelTotAdapM1->Add( hAdapTeO2po210M1, (1.77E+3)/dDataIntegralM1 );
  fModelTotAdapM2->Add( hAdapTeO2po210M2, (1.77E+3)/dDataIntegralM1 );

  // Adding Muon distribution, normalized to livetime of exposure
  fModelTotAdapM1->Add( hAdapExtMuonM1, (1/(7.33063e+10/(60*60*24*365))*dLivetimeYr) );
  fModelTotAdapM2->Add( hAdapExtMuonM2, (1/(7.33063e+10/(60*60*24*365))*dLivetimeYr) );

}
  
// This method sets up minuit and does the fit
bool TBackgroundModel::DoTheFit()
{ 
  // Reset initial parameter/error values
  ResetParameters();

  gStyle->SetOptStat(0);
  // gStyle->SetOptTitle(0);
  gStyle->SetOptFit();

  // Reduce Minuit Output
  minuit->SetPrintLevel(0); // Print level -1 (Quiet), 0 (Normal), 1 (Verbose)
  // minuit->Command("SET STRategy 1"); // Sets strategy of fit: 0, 1, or 2; 2 is the slowest but best
  minuit->Command("SET STRategy 2"); // Sets strategy of fit: 0, 1, or 2; 2 is the slowest but best  
  minuit->SetMaxIterations(100000);
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

   int status = minuit->Command("MINImize 100000 0.1"); // Command that actually does the minimization
   // int status = minuit->Command("MINImize 10000000 1"); // Command that actually does the minimization

  // Get final parameters from fit
  for(int i = 0; i < dNParam; i++)
  {
    minuit->GetParameter(i, fParameters[i], fParError[i]);
  }

  // Update model with final parameters
  UpdateModel();
  
  dChiSquare = GetChiSquare();

  // ///// Draw Data M1
  fAdapDataHistoM1->SetLineColor(kBlack);
  fAdapDataHistoM1->SetLineWidth(2);
  // fAdapDataHistoM1->SetMaximum(9000);
  fAdapDataHistoM1->SetMinimum(0.02);  

  fModelTotAdapM1->SetLineColor(2);
  fModelTotAdapM1->SetLineWidth(1);

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
  LineM1->DrawLineNDC(0.1, 0.59, 0.95, 0.59);

  TLegend *legfit1 = new TLegend(0.85,0.75,0.95,1);
  legfit1->SetFillColor(0);
  // legfit1->SetTextSize(0.02);
  legfit1->AddEntry(hRatioM1_1, "1 #sigma", "f");
  legfit1->AddEntry(hRatioM1_2, "2 #sigma", "f");
  legfit1->AddEntry(hRatioM1_3, "3 #sigma", "f");
  legfit1->Draw();

  p2m1->cd();
  p2m1->SetLogy();
  
  fAdapDataHistoM1->GetXaxis()->SetRange(fAdapDataHistoM1->FindBin(dFitMin), fAdapDataHistoM1->FindBin(dFitMax-1));
  fAdapDataHistoM1->SetTitle("Total Model (M1)");
  // fAdapDataHistoM1->SetTitleOffset(1.5);
  // fAdapDataHistoM1->SetTitleSize(0.005);
  fAdapDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM1->GetYaxis()->SetTitle("Counts/Bin");
  fAdapDataHistoM1->Draw("E");
  fModelTotAdapM1->Draw("SAME");

  TLegend *leg1 = new TLegend(0.78,0.75,0.95,0.92);
  leg1->SetFillColor(0);
  // leg1->SetTextSize(0.02);
  leg1->AddEntry(fAdapDataHistoM1, "Data", "l");
  leg1->AddEntry(fModelTotAdapM1, "Model", "l");
  leg1->Draw();


  ///// Draw Data M2
  fAdapDataHistoM2->SetLineColor(kBlack);
  fAdapDataHistoM2->SetLineWidth(2);
  // fAdapDataHistoM2->SetMinimum(0.002);
  // fAdapDataHistoM2->SetMaximum(9000);

  fModelTotAdapM2->SetLineColor(2);
  fModelTotAdapM2->SetLineWidth(1);

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
  LineM1->DrawLineNDC(0.1, 0.59, 0.95, 0.59);

  TLegend *legfit2 = new TLegend(0.85,0.75,0.95,1);
  legfit2->SetFillColor(0);
  // legfit2->SetTextSize(0.02);
  legfit2->AddEntry(hRatioM2_1, "1 #sigma", "f");
  legfit2->AddEntry(hRatioM2_2, "2 #sigma", "f");
  legfit2->AddEntry(hRatioM2_3, "3 #sigma", "f");
  legfit2->Draw();


  p2m2->cd();
  p2m2->SetLogy();
  fAdapDataHistoM2->GetXaxis()->SetRange(fAdapDataHistoM2->FindBin(dFitMin), fAdapDataHistoM2->FindBin(dFitMax-1));
  // fAdapDataHistoM2->SetTitleOffset(0.45);
  // fAdapDataHistoM2->SetTitleSize(0.01);
  fAdapDataHistoM2->SetTitle("Total Model (M2)");
  fAdapDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM2->GetYaxis()->SetTitle("Counts/Bin");
  fAdapDataHistoM2->Draw("E");
  fModelTotAdapM2->Draw("SAME");
  
  TLegend *leg2 = new TLegend(0.78,0.75,0.95,0.92);
  leg2->SetFillColor(0);
  // leg2->SetTextSize(0.02);
  leg2->AddEntry(fAdapDataHistoM2, "Data", "l");
  leg2->AddEntry(fModelTotAdapM2, "Model", "l");
  leg2->Draw();

  // Residuals
  TCanvas *cResidual1 = new TCanvas("cResidual1", "cResidual1", 1200, 800);
  hResidualGausM1 = new TH1D("hResidualGausM1", "M1", 100, -50, 50);
  hResidualDistM1 = CalculateResiduals(fAdapDataHistoM1, fModelTotAdapM1, hResidualGausM1, dFitMinBinM1, dFitMaxBinM1, 1);
  hResidualDistM1->SetLineColor(kBlack);
  hResidualDistM1->SetName("Residuals");
  hResidualDistM1->SetTitle("Normalized Residuals (M1)");
  hResidualDistM1->SetMarkerStyle(25);
  hResidualDistM1->GetXaxis()->SetTitle("Energy (keV)");
  // hResidualDistM1->GetXaxis()->SetTitleSize(0.04);
  // hResidualDistM1->GetXaxis()->SetLabelSize(0.05);
  // hResidualDistM1->GetYaxis()->SetLabelSize(0.05);
  // hResidualDistM1->GetYaxis()->SetTitleSize(0.04);
  hResidualDistM1->GetXaxis()->SetRange(fAdapDataHistoM1->FindBin(dFitMin), fAdapDataHistoM1->FindBin(dFitMax-1));
  hResidualDistM1->GetYaxis()->SetTitle("Residuals (#sigma)");
  hResidualDistM1->Draw("P");

  TCanvas *cResidual2 = new TCanvas("cResidual2", "cResidual2", 1200, 800);
  hResidualGausM2 = new TH1D("hResidualGausM2", "M2", 100, -50, 50);  
  hResidualDistM2 = CalculateResiduals(fAdapDataHistoM2, fModelTotAdapM2, hResidualGausM2, dFitMinBinM2, dFitMaxBinM2, 2);
  hResidualDistM2->SetLineColor(kBlack);
  hResidualDistM2->SetName("Residuals");
  hResidualDistM2->SetTitle("Normalized Residuals (M2)");
  hResidualDistM2->SetMarkerStyle(25);
  hResidualDistM2->GetXaxis()->SetTitle("Energy (keV)");
  // hResidualDistM2->GetXaxis()->SetTitleSize(0.04);
  // hResidualDistM2->GetXaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetLabelSize(0.05);
  // hResidualDistM2->GetYaxis()->SetTitleSize(0.04); 
  hResidualDistM2->GetXaxis()->SetRange(fAdapDataHistoM2->FindBin(dFitMin), fAdapDataHistoM2->FindBin(dFitMax-1));
  hResidualDistM2->GetYaxis()->SetTitle("Residuals (#sigma)");
  hResidualDistM2->Draw("P");

  // Draw residual distribution and fit with Gaussian
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

  dResidualRMSTot = TMath::Sqrt( (dResidualRMSM1 + dResidualRMSM2)/ (dNDF + dNumFreeParameters) );

  dResidualRMSM1 = TMath::Sqrt(dResidualRMSM1/(dFitMaxBinM1-dFitMinBinM1));
  dResidualRMSM2 = TMath::Sqrt(dResidualRMSM2/(dFitMaxBinM2-dFitMinBinM2));

  // Some basic checks on rates
  double dROIRange = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2570))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2570)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2470)); 
  double d2nbbRange = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2000))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2000)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(500));
  double dGammaRange = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(3000))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(3000)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(100));
  double dAlphaRangeM1 = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(8000))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(8000)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(3000));
  double dAlphaRangeM2 = fAdapDataHistoM2->GetBinLowEdge(fAdapDataHistoM2->FindBin(8000))+fAdapDataHistoM2->GetBinWidth(fAdapDataHistoM2->FindBin(8000)) - fAdapDataHistoM2->GetBinLowEdge(fAdapDataHistoM2->FindBin(3000));

  double d2nbbData = fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double d2nbbDataErr = TMath::Sqrt(fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);
  double d2nbbModel = fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double d2nbbModelErr = TMath::Sqrt(fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);
  double d2nbbPDF = dDataIntegralM1*fParameters[0]*hAdapTeO22nuM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
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
  
  if(!bFixedArray[0])
  {
    double d2nbbHL = fParEfficiencyM1[0]*(0.69314718056)*(4.93124196176940785e+25 * dLivetimeYr)/(dDataIntegralM1*fParameters[0]*hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width"));
    cout << "Counts in 2nbb (M1): " << dDataIntegralM1*fParameters[0]*hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width") << "\t Half-Life " << d2nbbHL << " +/- " << d2nbbHL*fParError[0]/fParameters[0] << endl;
  }
  
  cout << endl;
  cout << endl;
  cout << "Data in 2nbb region (c/keV/yr): " << d2nbbData << " $\\pm$ " << d2nbbDataErr << endl;  
  cout << "Model in 2nbb region (c/keV/yr): " << d2nbbModel << " $\\pm$ " << d2nbbModelErr << endl;
  cout << "Model of 2nbb PDF (c/keV/yr): " << d2nbbPDF << " $\\pm$ " << d2nbbPDFErr << endl;
  cout << "Model of PDFs in 2nbb region without 2nbb (c/keV/y): " << d2nbbModel - d2nbbPDF << " $\\pm$ " << d2nbbModelErr << endl;
  cout << "Data in 0nbb region (c/keV/yr): " << dROIData << " $\\pm$ " << dROIDataErr << endl;
  cout << "Model in 0nbb region (c/keV/yr): " << dROIModel << " $\\pm$ " << dROIModelErr << endl;
  cout << "Model/Data in 2nbb region: " << d2nbbModel/d2nbbData << endl;
  cout << "Model/Data in 0nbb region: " << dROIModel/dROIData << endl;
  cout << "Model of PDFs in 2nbb region without 2nbb/Data : " << (d2nbbModel - d2nbbPDF)/d2nbbData << endl;
  cout << "Model in Gamma region (c/keV/yr): " << dGammaModel << " $\\pm$ " << dGammaModelErr << endl;
  cout << "Model in Alpha region M1 (c/keV/yr): " << dAlphaModelM1 << " $\\pm$ " << dAlphaModelErrM1 << endl;
  cout << "Model in Alpha region M2 (c/keV/yr): " << dAlphaModelM2 << " $\\pm$ " << dAlphaModelErrM2 << endl;
  cout << endl;  
  cout << endl;    
  cout << "Residual RMS (Tot): " << dResidualRMSTot << endl;
  cout << "Residual RMS (M1): " << dResidualRMSM1 << "\t" << "Residual RMS (M2): " << dResidualRMSM2 << endl;

  dChiSquare = GetChiSquare();
  // cout << "Additional ChiSq: " << AddConstraint(26) << endl;
  cout << "Total number of calls = " << dNumCalls << "\t" << "ChiSq/NDF = " << dChiSquare/dNDF << endl; // for M1 and M2
  cout << "ChiSq = " << dChiSquare << "\t" << "NumFreeParamters = " << dNumFreeParameters << "\t" << "NDF = " << dNDF << endl;
  cout << "Probability = " << TMath::Prob(dChiSquare, dNDF) << endl;
  // cout << "Fit Status = " << status << endl;

/*
    2nbb calculation:
     - N_A = 6.0221409e+23
     - TeO2 molar mass: 159.6 g/mol
     - half life is 9.81 * 10^20 years
     - how many in Q0 data so far? 1/rate = half life/ln(2) -> rate = ln(2)/half life = 7.066*10^-22 decays/year (Laura's thesis)
     - Moles = 750g * 49 crystals * 0.3408 abundance/159.6 g/mol = 78.474 mol total
     - N_TeO2 = 78.474 * N_A = 4.726*10^25 nuclei of Te130
     - N_2nbb = N_TeO2 * rate * livetime = 1.551*10^4 events

     - Moles = 750g * 51 crystals * 0.3408 abundance/159.6 g/mol = 81.677 mol total
     - N_TeO2 = 81.677 * N_A = 4.9187e+25 nuclei of Te130
      -- 0.34167 abundance ==> N_TeO2 = 4.93124196176940785e+25
     - half life = rate * ln(2) = ln(2) * N_TeO2 * livetime / N_2nbb   
*/

  // Correlation Matrix section
  TMatrixT<double> mCorrMatrix;
  mCorrMatrix.ResizeTo(dNParam, dNParam);
  minuit->mnemat(mCorrMatrix.GetMatrixArray(), dNParam);

  for(int i = mCorrMatrix.GetRowLwb(); i <= mCorrMatrix.GetRowUpb(); i++)
    for(int j = mCorrMatrix.GetColLwb(); j <= mCorrMatrix.GetColUpb(); j++)
      mCorrMatrix(i,j) = mCorrMatrix(i,j)/(fParError[i]*fParError[j]);

  // Just the Matrix by itself
  gStyle->SetPalette(55);
  TCanvas *cMatrix2 = new TCanvas("cMatrix2", "cMatrix2", 1200, 1200);
  TPad* pM3 = new TPad("pM3","pM3",width1,canBotMargin,width2,canBotMargin+2*padHeight,0,0);
  pM3->SetTopMargin(0.05);
  pM3->SetBottomMargin(0.05);
  pM3->SetLeftMargin(0.05);
  pM3->SetRightMargin(0.15);
  pM3->SetFillColor(0);
  pM3->SetBorderMode(0);
  pM3->SetBorderSize(0);
  pM3->Draw();
  pM3->cd();
  pM3->SetGrid();
  pM3->SetTicks();
  pM3->SetFillStyle(4000);
  mCorrMatrix.Draw("colz");


  double dProgressM1 = 0;
  double dProgressM2 = 0;
  double dataM1_i = 0, modelM1_i = 0;
  double dataM2_i = 0, modelM2_i = 0;

  for(int i = dFitMinBinM1; i < dFitMaxBinM1; i++)
  {
    if( fAdapDataHistoM1->GetBinCenter(i) >= 3150 && fAdapDataHistoM1->GetBinCenter(i) <= 3400)continue;

    dataM1_i = fAdapDataHistoM1->GetBinContent(i)*fAdapDataHistoM1->GetBinWidth(i);
    modelM1_i = fModelTotAdapM1->GetBinContent(i)*fAdapDataHistoM1->GetBinWidth(i);
    dProgressM1 += 2 * (modelM1_i - dataM1_i + dataM1_i * TMath::Log(dataM1_i/modelM1_i));
    hChiSquaredProgressM1->SetBinContent(i, dProgressM1);
  }
  for(int i = dFitMinBinM2; i < dFitMaxBinM2; i++)
  {
    if( fAdapDataHistoM2->GetBinCenter(i) >= 3150 && fAdapDataHistoM2->GetBinCenter(i) <= 3400)continue;

    dataM2_i = fAdapDataHistoM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i);
    modelM2_i = fModelTotAdapM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i);
    dProgressM2 += 2 * (modelM2_i - dataM2_i + dataM2_i * TMath::Log(dataM2_i/modelM2_i));
    hChiSquaredProgressM2->SetBinContent(i, dProgressM2);
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

  // fParArray = new TVectorD(dNParam, fParameters);
  // fParArrayErr = new TVectorD(dNParam, fParError);

  if(bSave)
  {
    // Integrals
    for(int i = 0; i < dNParam; i++)
    {
      if(bFixedArray[i])continue;
      BkgPar[i]->GetHistM1(1)->Scale( dDataIntegralM1*fParameters[i]);
      BkgPar[i]->GetHistM2(1)->Scale( dDataIntegralM2*fParameters[i]);
      BkgPar[i]->GetHistM1(2)->Scale( dDataIntegralM1*fParameters[i]);
      BkgPar[i]->GetHistM2(2)->Scale( dDataIntegralM2*fParameters[i]);      
      fParCountsM1[i] = BkgPar[i]->GetHistM1(1)->Integral("width")/dLivetimeYr;
      fParCountsM2[i] = BkgPar[i]->GetHistM2(1)->Integral("width")/dLivetimeYr;    
      // Calculate Activity in terms of Bq/Kg
      fParActivityM1[i] = fParCountsM1[i]/fParEfficiencyM1[i]/fParMass[i]/(365*24*60*60);
    }
    // Scale Po-210 
    hTeO2po210M1->Scale((1.77E+3));
    hTeO2po210M2->Scale((1.77E+3));
    hAdapTeO2po210M1->Scale( (1.77E+3) );
    hAdapTeO2po210M2->Scale( (1.77E+3) );


  // Saving plots
    cadap1->SaveAs(Form("%s/Final/FitM1_%d.pdf", dSaveDir.c_str(), tTime->GetDate() ));
    cadap2->SaveAs(Form("%s/Final/FitM2_%d.pdf", dSaveDir.c_str(), tTime->GetDate() ));
    cResidual1->SaveAs(Form("%s/Final/FitM1Residual_%d.pdf", dSaveDir.c_str(), tTime->GetDate() ));
    cResidual2->SaveAs(Form("%s/Final/FitM2Residual_%d.pdf", dSaveDir.c_str(), tTime->GetDate() ));
    cres1->SaveAs(Form("%s/Final/FitResidualDist_%d.pdf", dSaveDir.c_str(), tTime->GetDate() ));
    cProgress->SaveAs(Form("%s/Final/ChiSquareProgress_%d.pdf", dSaveDir.c_str(), tTime->GetDate() ));
    cMatrix2->SaveAs(Form("%s/Final/FitCovMatrix_%d.pdf", dSaveDir.c_str(), tTime->GetDate() ));


    // Save histograms to file
    TFile *fSaveResult = new TFile(Form("%s/Final/FitResult_%d.root", dSaveDir.c_str(), tTime->GetDate() ), "RECREATE");
    fSaveResult->cd();
    fDataHistoM1->Write();
    fDataHistoM2->Write();
    fAdapDataHistoM1->Write();
    fAdapDataHistoM2->Write();
    fModelTotAdapM1->Write();
    fModelTotAdapM2->Write();
    for(int i = 0; i < dNParam; i++)
    {
      BkgPar[i]->GetHistM1(1)->Write();
      BkgPar[i]->GetHistM2(1)->Write();
      BkgPar[i]->GetHistM1(2)->Write();
      BkgPar[i]->GetHistM2(2)->Write();    
    }
    hAdapTeO2po210M1->Write();
    hAdapTeO2po210M2->Write();
    hTeO2po210M1->Write();
    hTeO2po210M2->Write();
    mCorrMatrix.Write();

    // fParArray->Write();
    // fParArrayErr->Write();
 
    fSaveResult->Close();
  } // end bSave

  return true;
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

  cout << "Adding " << (0.25*dM1)*hAdapTeO22nuM1->Integral("width") + (0.25*dM2)*hAdapTeO22nuM2->Integral("width") << " 2nbb events to spectrum" << endl;

}

// Only run in batch mode and make sure to have the 2nbb normalization FIXED
// Probably best to run Minuit in quiet mode as well (to avoid the output)
// eg: minuit->SetPrintLevel(-1);
void TBackgroundModel::ProfileNLL(int fParFixed)
{
  // Do the fit normally once first to obtain best fit value
  DoTheFit();

  // Fix 2nbb value now
  bFixedArray[fParFixed] = true;
  
  // Chi-Squared from best fit (for ProfileNLL calculation)
  dBestChiSq = dChiSquare; 

  // Do the fit now if no other tests are needed 
  // for(int i = -50; i < 50; i++)
  for(int i = -5; i < 5; i++)
  {
    fInitValues.push_back(fParameters[fParFixed] + fParameters[fParFixed]/50*i );
    cout << "Input initial value: " << fParameters[fParFixed] + fParameters[fParFixed]/50*i << endl;
  }

    int dStatus, dIndex = 0; 
    double deltaChiSq;
    double Toy2nbbHL, Toy2nbbHLErr;
    double Toy2nbbIntegral;
    double Toy2nbbIntegralErr;
    double Toy2nbbHLErr2, dPull2, Toy2nbbIntegralErr2;

    TFile *ProfileFile = new TFile(Form("%s/Final/ProfileNLL/ProfileNLL_Par%d_%d.root", dSaveDir.c_str(), fParFixed, tTime->GetDate() ), "RECREATE");
    TTree *ProfileTree = new TTree("ProfileTree", "Tree with Profile NLL Results");

    ProfileTree->Branch("Index", &dIndex, "Index/I");
    ProfileTree->Branch("ChiSq", &dChiSquare, "ChiSq/D");
    ProfileTree->Branch("deltaChiSq", &deltaChiSq, "Pull/D");
    ProfileTree->Branch("dDataIntegralM1", &dDataIntegralM1, "dDataIntegralM1/D");
    ProfileTree->Branch("fParameters", &fParameters, "fParameters[50]/D");
    ProfileTree->Branch("fParError", &fParError, "fParError[50]/D");    
    ProfileTree->Branch("FitStatus", &dStatus, "FitStatus/I");
    ProfileTree->Branch("Toy2nbbHL", &Toy2nbbHL, "Toy2nbbHL/D");
    ProfileTree->Branch("Toy2nbbHLErr", &Toy2nbbHLErr, "Toy2nbbHLErr/D");
    ProfileTree->Branch("Toy2nbbIntegral", &Toy2nbbIntegral, "Toy2nbbIntegral/D");
    ProfileTree->Branch("Toy2nbbIntegralErr", &Toy2nbbIntegralErr, "Toy2nbbIntegralErr/D");

    ProfileTree->Branch("fAdapDataHistoM1", "TH1D", &fAdapDataHistoM1, 32000, 0);
    ProfileTree->Branch("fAdapDataHistoM2", "TH1D", &fAdapDataHistoM2, 32000, 0);
    ProfileTree->Branch("fModelTotAdapM1", "TH1D", &fModelTotAdapM1, 32000, 0);
    ProfileTree->Branch("fModelTotAdapM2", "TH1D", &fModelTotAdapM2, 32000, 0);

    TH1::AddDirectory(kFALSE);

    // Write 0th entry
    ProfileTree->Write();

  for(std::vector<double>::const_iterator iter = fInitValues.begin(); iter!=fInitValues.end(); iter++)
  {
    // Set new initial value and repeat fit
    cout << "Using Initial Value: " << *iter << endl;
    BkgPar[fParFixed]->SetInitValue(*iter);
    DoTheFit();
    dIndex++;
    
    deltaChiSq = dChiSquare - dBestChiSq;
    cout << "delta ChiSq = " << dChiSquare - dBestChiSq << endl; // Needs to be entered, otherwise just 0

    ProfileTree->Write();

  }

  ProfileFile->cd();
  ProfileTree->Write();
  ProfileFile->Close();

}

// Gets limit for fitted limited parameter
void TBackgroundModel::SetLimit(int fParFixed)
{
  // First do the fit to get parameters
  DoTheFit();
  dBestChiSq = dChiSquare; // Chi-Squared from best fit (for ProfileNLL calculation)

  double  dChiSquareDummy;
  double  dDummyIntegral;

  dDummyIntegral = BkgPar[fParFixed]->GetHistM1(1)->Integral("width");


  OutFile.open(Form("%s/Final/Limits/Limit_Par%d_%d.txt", dSaveDir.c_str(), fParFixed, tTime->GetDate() ));


  for(int i = 0; i < 300 ; i++ )
  {
    fInitValues.push_back(fParError[fParFixed]/13000*i );
    // cout << "Input initial value: " << fParameters[fParFixed] + fParError[fParFixed]/200*i << endl;
  }

  // OutFile << "ChiSquare ---  delta ChiSquare/2 ---  Parameter  --- Activity (Bq/kg)" << endl;

  for(std::vector<double>::const_iterator iter = fInitValues.begin(); iter!=fInitValues.end(); iter++)
  {
    UpdateModel(); 
    // Set new value for parameter
    fParameters[fParFixed] = *iter;

    // cout << "Integral before scaling: " << dDummyIntegral << endl;
    // cout << "Scaling: " << *iter*dDataIntegralM1 << "\t" << "Integral: " << dDummyIntegral << endl;
    // Re-calculate Chi-squared
    dChiSquareDummy = GetChiSquare();
    OutFile << dChiSquareDummy << "\t" << (dChiSquareDummy - dBestChiSq)/2 << "\t" << *iter << "\t" << dDummyIntegral*dDataIntegralM1*(*iter)/dLivetimeYr/fParEfficiencyM1[fParFixed]/fParMass[fParFixed]/(365*24*60*60) << endl; 
    cout << dChiSquareDummy << "\t" << (dChiSquareDummy - dBestChiSq)/2 << "\t" << *iter << "\t" << dDummyIntegral*dDataIntegralM1*(*iter)/dLivetimeYr/fParEfficiencyM1[fParFixed]/fParMass[fParFixed]/(365*24*60*60) << endl; 
  }
}

// 2D Contour plot, only in comparison with 2nbb!
// Try out: 38, 42, 44, 47, 50
void TBackgroundModel::ProfileNLL2D(int fParFixed)
{
  // dBestChiSq = fBestFitChiSq; // Chi-Squared from best fit (for ProfileNLL calculation)
  // Do the fit now if no other tests are needed 
  
  // Do the fit normally once first
  DoTheFit();

    int dStatus, dIndex; 
    double dPull = 0;
    double Toy2nbbHL, Toy2nbbHLErr;
    double Toy2nbbIntegral;
    double Toy2nbbIntegralErr;
    double Toy2nbbHLErr2, dPull2, Toy2nbbIntegralErr2;

  // Fix 2nbb value and TeO2 K-40 value
  bFixedArray[0] = true;
  bFixedArray[fParFixed] = true;

  dBestChiSq = dChiSquare;

  for(int i = -5; i < 5; i++)
  {
    fInitValues.push_back(fParameters[0] + fParameters[0]/20*i);
  }
  for(int j = -30; j < 30; j++)
  {
    fInitValues2.push_back(fParameters[fParFixed] + fParameters[fParFixed]/10*j);
  }

    TFile *ProfileFile = new TFile(Form("%s/Final/ProfileNLL/ProfileFile_Par%d_%d.root", dSaveDir.c_str(), fParFixed, tTime->GetDate() ), "RECREATE");
    
    TTree *ProfileTree = new TTree("ProfileTree", "Tree with Profile NLL Results");

    ProfileTree->Branch("Index", &dIndex, "Index/I");
    ProfileTree->Branch("ChiSq", &dChiSquare, "ChiSq/D");
    ProfileTree->Branch("Pull", &dPull, "Pull/D");
    ProfileTree->Branch("dDataIntegralM1", &dDataIntegralM1, "dDataIntegralM1/D");
    ProfileTree->Branch("fParameters", &fParameters, "fParameters[50]/D");
    ProfileTree->Branch("fParError", &fParError, "fParError[50]/D");    
    ProfileTree->Branch("FitStatus", &dStatus, "FitStatus/I");
    ProfileTree->Branch("Toy2nbbHL", &Toy2nbbHL, "Toy2nbbHL/D");
    ProfileTree->Branch("Toy2nbbHLErr", &Toy2nbbHLErr, "Toy2nbbHLErr/D");
    ProfileTree->Branch("Toy2nbbIntegral", &Toy2nbbIntegral, "Toy2nbbIntegral/D");
    ProfileTree->Branch("Toy2nbbIntegralErr", &Toy2nbbIntegralErr, "Toy2nbbIntegralErr/D");

    ProfileTree->Branch("fAdapDataHistoM1", "TH1D", &fAdapDataHistoM1, 32000, 0);
    ProfileTree->Branch("fAdapDataHistoM2", "TH1D", &fAdapDataHistoM2, 32000, 0);
    ProfileTree->Branch("fModelTotAdapM1", "TH1D", &fModelTotAdapM1, 32000, 0);
    ProfileTree->Branch("fModelTotAdapM2", "TH1D", &fModelTotAdapM2, 32000, 0);

    TH1::AddDirectory(kFALSE);


  for(std::vector<double>::const_iterator iter = fInitValues.begin(); iter != fInitValues.end(); iter++)
  {
    for(std::vector<double>::const_iterator iter2 = fInitValues2.begin(); iter2 != fInitValues2.end(); iter2++)
    {
    cout << "Step: " << *iter << " " << *iter2 << endl;    
    BkgPar[0]->SetInitValue(*iter);
    BkgPar[fParFixed]->SetInitValue(*iter2);

    DoTheFit();
    cout << "delta ChiSq = " << dChiSquare - dBestChiSq << endl; // Needs to be entered, otherwise just 0
    
    ProfileTree->Fill();
    }
  }
    ProfileFile->cd();
    ProfileTree->Write();
    ProfileFile->Close();
}

void TBackgroundModel::ToyFit(int fStart, int fStop)
{
    cout << "Using Toy Data" << endl;
    int dStatus, dIndex; 
    double dPull = 0;
    double Toy2nbbHL, Toy2nbbHLErr;
    double Toy2nbbIntegral;
    double Toy2nbbIntegralErr;
    double Toy2nbbHLErr2, dPull2, Toy2nbbIntegralErr2;

    TFile *ToyTreeFile = new TFile(Form("%s/Final/ToyMC/ToyFile_%d.root", dSaveDir.c_str(), tTime->GetDate() ), "RECREATE");
    TTree *ToyTree = new TTree("ToyTree", "Tree with Toy Fit Results");

    ToyTree->Branch("Index", &dIndex, "Index/I");
    ToyTree->Branch("ChiSq", &dChiSquare, "ChiSq/D");
    ToyTree->Branch("Pull", &dPull, "Pull/D");
    ToyTree->Branch("dDataIntegralM1", &dDataIntegralM1, "dDataIntegralM1/D");
    ToyTree->Branch("fParameters", &fParameters, "fParameters[50]/D");
    ToyTree->Branch("fParError", &fParError, "fParError[50]/D");    
    ToyTree->Branch("FitStatus", &dStatus, "FitStatus/I");
    ToyTree->Branch("Toy2nbbHL", &Toy2nbbHL, "Toy2nbbHL/D");
    ToyTree->Branch("Toy2nbbHLErr", &Toy2nbbHLErr, "Toy2nbbHLErr/D");
    ToyTree->Branch("Toy2nbbIntegral", &Toy2nbbIntegral, "Toy2nbbIntegral/D");
    ToyTree->Branch("Toy2nbbIntegralErr", &Toy2nbbIntegralErr, "Toy2nbbIntegralErr/D");

    ToyTree->Branch("Pull2", &dPull2, "Pull2/D");
    ToyTree->Branch("Toy2nbbHLErr2", &Toy2nbbHLErr2, "Toy2nbbHLErr2/D");
    ToyTree->Branch("Toy2nbbIntegralErr2", &Toy2nbbIntegralErr2, "Toy2nbbIntegralErr2/D");

    ToyTree->Branch("fAdapDataHistoM1", "TH1D", &fAdapDataHistoM1, 32000, 0);
    ToyTree->Branch("fAdapDataHistoM2", "TH1D", &fAdapDataHistoM2, 32000, 0);
    ToyTree->Branch("fModelTotAdapM1", "TH1D", &fModelTotAdapM1, 32000, 0);
    ToyTree->Branch("fModelTotAdapM2", "TH1D", &fModelTotAdapM2, 32000, 0);

    // Initial Fit to get the rate
    dStatus = DoTheFit();

    double dInitial2nbbRate, dInitial2nbbRateErr;
    double dInitial2nbb, dInitial2nbbErr;
    dInitial2nbb = fParameters[0];
    dInitial2nbbErr = fParError[0];
    dInitial2nbbRate = fParEfficiencyM1[0]*(0.69314718056)*(4.93124196176940785e+25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1*hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width"));
    dInitial2nbbRateErr = fParError[0]/fParameters[0]*(fParEfficiencyM1[0]*(0.69314718056)*(4.93124196176940785e+25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1*hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width")) );

    dIndex = 0;
    Toy2nbbHL = dInitial2nbbRate;
    Toy2nbbHLErr = dInitial2nbbRateErr;
    Toy2nbbIntegral = fParameters[0]*dDataIntegralM1*hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width");
    Toy2nbbIntegralErr = fParError[0]*dDataIntegralM1*hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width");

    ToyTree->Fill();

    cout << "Initial 2nbb Rate and Error: " << dInitial2nbbRate << " +/- " <<  dInitial2nbbRateErr << endl;

    TH1::AddDirectory(kFALSE);

    TFile *fToyDataTest;

    for(int i = fStart; i <= fStop; i++)
    {
      cout << "Toy: " << i << endl;
      dIndex = i;
      // Reset loaded data
      fAdapDataHistoM1->Reset();
      fAdapDataHistoM2->Reset();
      
      // Load Toy data
      fToyDataTest = new TFile(Form("%s/Toy/ToyData_p%d.root", dMCDir.c_str(), i));
       // fToyDataTest = new TFile(Form("%s/Toy/HighStatToyData_p%d.root", dMCDir.c_str(), i));
      // fToyDataTest = TFile::Open(Form("%s/Toy/HighStatToyData_p%d.root", dMCDir.c_str(), i));
      // TFile fToyDataTest(Form("%s/Toy/HighStatToyData_p%d.root", dMCDir.c_str(), i));

      fAdapDataHistoM1 = (TH1D*)fToyDataTest->Get("fAdapDataHistoM1");
      fAdapDataHistoM2 = (TH1D*)fToyDataTest->Get("fAdapDataHistoM2");
      
      // Scaling for high statistics Toy MC
      // fAdapDataHistoM1->Scale(1./5000);
      // fAdapDataHistoM2->Scale(1./5000);
      // dDataIntegralM1 = fAdapDataHistoM1->Integral("width");

      dStatus = DoTheFit();

      Toy2nbbHL = fParEfficiencyM1[0]*(0.69314718056)*(4.93124196176940785e+25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1*hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width"));
      Toy2nbbHLErr = fParError[0]/fParameters[0]*(fParEfficiencyM1[0]*(0.69314718056)*(4.93124196176940785e+25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1*hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width")));
      Toy2nbbIntegral = fParameters[0]*dDataIntegralM1*hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width");
      Toy2nbbIntegralErr = fParError[0]*dDataIntegralM1*hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width");

      // Alternate calculation, used only for the high statistics
      Toy2nbbHLErr2 = fParError[0]/fParameters[0]*(fParEfficiencyM1[0]*(0.69314718056)*(4.93124196176940785e+25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1*hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width")))/5000;
      Toy2nbbIntegralErr2 = fParError[0]*dDataIntegralM1*hAdapTeO22nuM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width")/5000;

      dPull = (Toy2nbbHL - dInitial2nbbRate)/(Toy2nbbHLErr);
      dPull2 = (fParameters[0] - dInitial2nbb)/(fParError[0]);

      // cout << "Toy Integral: " << dDataIntegralM1 << "\t Toy Livetime: " << dDummyLivetimeYr << endl;
      cout << "Toy Counts: " << Toy2nbbIntegral << " +/- " <<  Toy2nbbIntegralErr << endl;
      cout << "Toy Half-Life: " << Toy2nbbHL << " +/- " <<  Toy2nbbHLErr << endl;
      cout << "Pull1 : " << dPull << "\t Pull2:" << dPull2 << endl;

      ToyTree->Fill();

      // delete fToyDataTest;
      fToyDataTest->Close();
    }

    ToyTreeFile->cd();
    ToyTree->Write();
    ToyTreeFile->Close();
}

// Add-on to the Chi-Square, more terms
double TBackgroundModel::AddConstraint(int fParInput)
{
  double dAdditionalChiSq = 0;

  // fParameters[fParInput]*dDataIntegralM1/dLivetimeYr/fParEfficiencyM1[i]/fParMass[i]/(365*24*60*60)
  if(fParPrior[fParInput] != 0)
  {
  double dModel_Integral = fParameters[fParInput]*dDataIntegralM1;
  double dData_Integral = fParPrior[fParInput]*dLivetimeYr*fParEfficiencyM1[fParInput]*fParMass[fParInput]*(365*24*60*60);
  dAdditionalChiSq += 2 * (dModel_Integral - dData_Integral + dData_Integral * TMath::Log(dData_Integral/dModel_Integral));
  }

  return dAdditionalChiSq;
}

// Sets the constants for MonteCarlo efficiency as well as Source mass
// Used to convert counting rate into activity in Bq/Kg
void TBackgroundModel::SetParEfficiency()
{
  fParEfficiencyM1[0] = 9.6696e-01; // For IK
  // fParEfficiencyM1[0] = 9.6677e-01; // For Primakoff Rosen => Basically the same
  fParEfficiencyM1[1] = 0.892441;
  fParEfficiencyM1[2] = 0.437839;
  fParEfficiencyM1[3] = 0.924953;
  fParEfficiencyM1[4] = 0.980596;
  fParEfficiencyM1[5] = 4.61506;
  fParEfficiencyM1[6] = 3.35465;
  fParEfficiencyM1[7] = 0.980505;
  fParEfficiencyM1[8] = 4.35318;

  fParEfficiencyM1[9] = 0.7101714;
  fParEfficiencyM1[10] = 3.30111;
  fParEfficiencyM1[11] = 2.17821;
  fParEfficiencyM1[12] = 0.6965972;
  fParEfficiencyM1[13] = 2.82943;
  // fParEfficiencyM1[13] = 1.41013;
  fParEfficiencyM1[14] = 1.36804;
  fParEfficiencyM1[15] = 1.5236;
  // fParEfficiencyM1[15] = 1.81667;
  
  fParEfficiencyM1[16] = 4.33309;
  fParEfficiencyM1[17] = 4.68744;
  // fParEfficiencyM1[20] = 8.99411;
  fParEfficiencyM1[18] = 0.729299;
  fParEfficiencyM1[19] = 0.721855;
  fParEfficiencyM1[20] = 1.1924;
  fParEfficiencyM1[21] = 0.335822;
  fParEfficiencyM1[22] = 0.367055;
  fParEfficiencyM1[23] = 0.380997;
  fParEfficiencyM1[24] = 1.95302;

  fParEfficiencyM1[25] = 0.208363;
  fParEfficiencyM1[26] = 0.48936;
  fParEfficiencyM1[27] = 0.540099;
  fParEfficiencyM1[28] = 0.0389089;
  fParEfficiencyM1[29] = 0.293154;
  fParEfficiencyM1[30] = 0.074132;
  fParEfficiencyM1[31] = 0.22195;
  fParEfficiencyM1[32] = 0.258695;
  fParEfficiencyM1[33] = 0.0127853;
  fParEfficiencyM1[34] = 0.158158;
  fParEfficiencyM1[35] = 0.0558653;
  fParEfficiencyM1[36] = 0.0653826;
  fParEfficiencyM1[37] = 0.00499584;
  fParEfficiencyM1[38] = 0.0751074;
  fParEfficiencyM1[39] = 0.00391235;
  fParEfficiencyM1[40] = 0.0047853;
  fParEfficiencyM1[41] = 0.000452146;
  fParEfficiencyM1[42] = 0.00734568;
  fParEfficiencyM1[43] = 0.00335353;
  fParEfficiencyM1[44] = 0.000228305;
  fParEfficiencyM1[45] = 0.000219387;
  fParEfficiencyM1[46] = 0.000501857;
  fParEfficiencyM1[47] = 0.000398748;
  fParEfficiencyM1[48] = 0.000309365;
  fParEfficiencyM1[49] = 0.000338529;

  // In Kg
  fParMass[0] = 38.25;
  fParMass[1] = 38.25;
  fParMass[2] = 38.25;
  fParMass[3] = 38.25;
  fParMass[4] = 38.25;
  fParMass[5] = 38.25;
  fParMass[6] = 38.25;
  fParMass[7] = 38.25;
  fParMass[8] = 38.25;
  fParMass[9] = 38.25;
  fParMass[10] = 38.25;
  fParMass[11] = 38.25;
  fParMass[12] = 38.25;
  fParMass[13] = 38.25;
  // fParMass[14] = 38.25;
  fParMass[14] = 38.25;
  fParMass[15] = 38.25;
  // fParMass[17] = 38.25;
  fParMass[16] = 38.25;
  fParMass[17] = 38.25;
  // fParMass[20] = 38.25;
  fParMass[18] = (2610.04+6929.71)/1000;
  fParMass[19] = (2610.04+6929.71)/1000;
  fParMass[20] = (2610.04+6929.71)/1000;
  fParMass[21] = (2610.04+6929.71)/1000;
  fParMass[22] = (2610.04+6929.71)/1000;
  fParMass[23] = (2610.04+6929.71)/1000;
  fParMass[24] = (2610.04+6929.71)/1000;

  fParMass[25] = (2610.04+6929.71)/1000;
  fParMass[26] = (2610.04+6929.71)/1000;
  fParMass[27] = (2610.04+6929.71)/1000;
  fParMass[28] = (2610.04+6929.71)/1000;
  fParMass[29] = (2610.04+6929.71)/1000;
  fParMass[30] = 82.26006;
  fParMass[31] = 82.26006;
  fParMass[32] = 82.26006;
  fParMass[33] = 82.26006;
  fParMass[34] = 82.26006;
  fParMass[35] = 202294.46/1000;
  fParMass[36] = 202294.46/1000;
  fParMass[37] = 202294.46/1000;
  fParMass[38] = 202294.46/1000;
  fParMass[39] = 180704.38/1000;
  fParMass[40] = 180704.38/1000;
  fParMass[41] = 180704.38/1000;
  fParMass[42] = 180704.38/1000;
  fParMass[43] = 180704.38/1000;
  fParMass[44] = 24693984./1000;
  fParMass[45] = 24693984./1000;
  fParMass[46] = 24693984./1000;
  fParMass[47] = 24693984./1000;
  fParMass[48] = 24693984./1000;
  fParMass[49] = 24693984./1000;

  // Sets Priors (in Bq/Kg)
  fParPrior[0] = 0.;
  fParPrior[1] = 1.0E-3;
  fParPrior[2] = 3.0E-7;
  fParPrior[3] = 0.;
  fParPrior[4] = 0.;
  fParPrior[5] = 0.;
  fParPrior[6] = 0.;
  fParPrior[7] = 0.;
  fParPrior[8] = 0.;
  fParPrior[9] = 0.;
  fParPrior[10] = 0.;
  fParPrior[11] = 0.;
  fParPrior[12] = 0.;
  fParPrior[13] = 0.;
  fParPrior[14] = 0.;
  fParPrior[15] = 0.;
  fParPrior[16] = 0.;
  fParPrior[17] = 0.;
  fParPrior[18] = 0.;
  fParPrior[19] = 0.;
  fParPrior[20] = 0.;
  fParPrior[21] = 0.;
  fParPrior[22] = 0.;
  fParPrior[23] = 0.;
  fParPrior[24] = 0.;
  fParPrior[25] = 0.;
  fParPrior[26] = 6.5E-5;
  fParPrior[27] = 2.0E-6;
  fParPrior[28] = 0.; // 7+/-2E-4 measured
  fParPrior[29] = 0.; // 5+/-1E-5
  fParPrior[30] = 0.;
  fParPrior[31] = 6.7E-4;
  fParPrior[32] = 4.4E-4;
  fParPrior[33] = 1.8E-4;
  fParPrior[34] = 0.;
  fParPrior[35] = 4.5E-5;
  fParPrior[36] = 3.3E-5;
  fParPrior[37] = 2.3E-5;
  fParPrior[38] = 0.;
  fParPrior[39] = 0.;
  fParPrior[40] = 0.;
  fParPrior[41] = 0.;
  fParPrior[42] = 4.2E-5;
  fParPrior[43] = 0.;
  fParPrior[44] = 0.;
  fParPrior[45] = 9.0E-3;
  fParPrior[46] = 2.6E-4;
  fParPrior[47] = 4.1E-4;
  fParPrior[48] = 0.;
  fParPrior[49] = 9.0E-3;
}


void TBackgroundModel::CreateModelHistograms()
{
  // Total model histograms
  fModelTotM1      = new TH1D("fModelTotM1",      "Total Model",        dNBins, dMinEnergy, dMaxEnergy);  
  fModelTotM2      = new TH1D("fModelTotM2",      "Total Model",        dNBins, dMinEnergy, dMaxEnergy);  

  /////////////////// Model histograms
  //////////// Crystal M1 and M2
  hTeO20nuM1       = new TH1D("hTeO20nuM1",    "hTeO20nuM1",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO22nuM1       = new TH1D("hTeO22nuM1",    "hTeO22nuM1",    dNBins, dMinEnergy, dMaxEnergy);
  hTeO22nuIKM1       = new TH1D("hTeO22nuIKM1",    "hTeO22nuIKM1",    dNBins, dMinEnergy, dMaxEnergy);
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
  hTeO22nuIKM2       = new TH1D("hTeO22nuIKM2",    "hTeO22nuIKM2",    dNBins, dMinEnergy, dMaxEnergy);
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

/////////// 50mK M1 and M2
  h50mKcs137M1     = new TH1D("h50mKcs137M1",  "h50mKcs137M1",  dNBins, dMinEnergy, dMaxEnergy);
  h50mKcs137M2     = new TH1D("h50mKcs137M2",  "h50mKcs137M2",  dNBins, dMinEnergy, dMaxEnergy);

/////////// Internal Shields M1 and M2
  hInternalco60M1   = new TH1D("hInternalco60M1", "hInternalco60M1",    dNBins, dMinEnergy, dMaxEnergy);
  hInternalk40M1    = new TH1D("hInternalk40M1", "hInternalk40M1",      dNBins, dMinEnergy, dMaxEnergy);
  hInternalth232M1  = new TH1D("hInternalth232M1", "hInternalth232M1",  dNBins, dMinEnergy, dMaxEnergy);
  hInternalu238M1   = new TH1D("hInternalu238M1", "hInternalu238M1",    dNBins, dMinEnergy, dMaxEnergy);

  hInternalco60M2   = new TH1D("hInternalco60M2", "hInternalco60M2",    dNBins, dMinEnergy, dMaxEnergy);
  hInternalk40M2    = new TH1D("hInternalk40M2", "hInternalk40M2",      dNBins, dMinEnergy, dMaxEnergy);
  hInternalth232M2  = new TH1D("hInternalth232M2", "hInternalth232M2",  dNBins, dMinEnergy, dMaxEnergy);
  hInternalu238M2   = new TH1D("hInternalu238M2", "hInternalu238M2",    dNBins, dMinEnergy, dMaxEnergy);


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

///////// OVC M1 and M2
  hOVCco60M1      = new TH1D("hOVCco60M1",   "hOVCco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hOVCk40M1       = new TH1D("hOVCk40M1",    "hOVCk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hOVCth232M1     = new TH1D("hOVCth232M1",  "hOVCth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hOVCu238M1      = new TH1D("hOVCu238M1",   "hOVCu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hOVCco60M2      = new TH1D("hOVCco60M2",   "hOVCco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hOVCk40M2       = new TH1D("hOVCk40M2",    "hOVCk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hOVCth232M2     = new TH1D("hOVCth232M2",  "hOVCth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hOVCu238M2      = new TH1D("hOVCu238M2",   "hOVCu238M2",   dNBins, dMinEnergy, dMaxEnergy);  

////////// External Shields
  hExtPbbi210M1 = new TH1D("hExtPbbi210M1", "hExtPbbi210M1", dNBins, dMinEnergy, dMaxEnergy);
  hExtPbk40M1 = new TH1D("hExtPbk40M1", "hExtPbk40M1", dNBins, dMinEnergy, dMaxEnergy);
  hExtPbth232M1 = new TH1D("hExtPbth232M1", "hExtPbth232M1", dNBins, dMinEnergy, dMaxEnergy);
  hExtPbu238M1 = new TH1D("hExtPbu238M1", "hExtPbu238M1", dNBins, dMinEnergy, dMaxEnergy);
  hExtPbpb210M1 = new TH1D("hExtPbpb210M1", "hExtPbpb210M1", dNBins, dMinEnergy, dMaxEnergy);

  hExtPbbi210M2 = new TH1D("hExtPbbi210M2", "hExtPbbi210M2", dNBins, dMinEnergy, dMaxEnergy);
  hExtPbk40M2 = new TH1D("hExtPbk40M2", "hExtPbk40M2", dNBins, dMinEnergy, dMaxEnergy);
  hExtPbth232M2 = new TH1D("hExtPbth232M2", "hExtPbth232M2", dNBins, dMinEnergy, dMaxEnergy);
  hExtPbu238M2 = new TH1D("hExtPbu238M2", "hExtPbu238M2", dNBins, dMinEnergy, dMaxEnergy);
  hExtPbpb210M2 = new TH1D("hExtPbpb210M2", "hExtPbpb210M2", dNBins, dMinEnergy, dMaxEnergy);

  hExtMuonM1 = new TH1D("hExtMuonM1", "hExtMuonM1", dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_th232spotM1 = new TH1D("hCuBox_th232spotM1", "hCuBox_th232spotM1", dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_k40spotM1 = new TH1D("hCuBox_k40spotM1", "hCuBox_k40spotM1", dNBins, dMinEnergy, dMaxEnergy);
  hBotExtPb_k40spotM1 = new TH1D("hBotExtPb_k40spotM1", "hBotExtPb_k40spotM1", dNBins, dMinEnergy, dMaxEnergy);

  hExtMuonM2 = new TH1D("hExtMuonM2", "hExtMuonM2", dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_th232spotM2 = new TH1D("hCuBox_th232spotM2", "hCuBox_th232spotM2", dNBins, dMinEnergy, dMaxEnergy);
  hCuBox_k40spotM2 = new TH1D("hCuBox_k40spotM2", "hCuBox_k40spotM2", dNBins, dMinEnergy, dMaxEnergy);
  hBotExtPb_k40spotM2 = new TH1D("hBotExtPb_k40spotM2", "hBotExtPb_k40spotM2", dNBins, dMinEnergy, dMaxEnergy);

//////////////// Adaptive binned histograms
////////// Total Adaptive binning histograms
  fModelTotAdapM1      = new TH1D("fModelTotAdapM1",      "Total PDF M1", dAdaptiveBinsM1, dAdaptiveArrayM1);  
  fModelTotAdapM2      = new TH1D("fModelTotAdapM2",      "Total PDF M2", dAdaptiveBinsM2, dAdaptiveArrayM2);  

/////////// Crystal M1 and M2
  hAdapTeO20nuM1       = new TH1D("hAdapTeO20nuM1",    "TeO2 Bulk 0nu M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO22nuM1       = new TH1D("hAdapTeO22nuM1",    "TeO2 Bulk 2nu M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO22nuIKM1       = new TH1D("hAdapTeO22nuIKM1",    "TeO2 Bulk 2nu M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
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
  hAdapTeO22nuIKM2       = new TH1D("hAdapTeO22nuIKM2",    "TeO2 Bulk 2nu M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
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

//////////// CuBox + CuFrame M1 and M2

  hAdapCuBox_CuFrameco60M1 = new TH1D("hAdapCuBox_CuFrameco60M1", "CuBox+CuFrame Bulk co60 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramek40M1 = new TH1D("hAdapCuBox_CuFramek40M1", "CuBox+CuFrame Bulk k40 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameth232M1 = new TH1D("hAdapCuBox_CuFrameth232M1", "CuBox+CuFrame Bulk th232 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFrameu238M1 = new TH1D("hAdapCuBox_CuFrameu238M1", "CuBox+CuFrame Bulk u238 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramemn54M1 = new TH1D("hAdapCuBox_CuFramemn54M1", "CuBox+CuFrame Bulk mn54 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_CuFramebi207M1 = new TH1D("hAdapCuBox_CuFramebi207M1", "CuBox+CuFrame Bulk bi207 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);


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
  hAdapCuBox_CuFramemn54M2 = new TH1D("hAdapCuBox_CuFramemn54M2", "CuBox+CuFrame Bulk mn54 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_CuFramebi207M2 = new TH1D("hAdapCuBox_CuFramebi207M2", "CuBox+CuFrame Bulk bi207 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);


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

////////// 50mK M1 and M2
  hAdap50mKcs137M1     = new TH1D("hAdap50mKcs137M1",  "50mK Bulk cs137 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKcs137M2     = new TH1D("hAdap50mKcs137M2",  "50mK Bulk cs137 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);

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

///////// Internal Shields M1 and M2
  hAdapInternalco60M1 = new TH1D("hAdapInternalco60M1", "Internal Shields Bulk co60 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapInternalk40M1 = new TH1D("hAdapInternalk40M1", "Internal Shields Bulk k40 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapInternalth232M1 = new TH1D("hAdapInternalth232M1", "Internal Shields Bulk th232 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapInternalu238M1 = new TH1D("hAdapInternalu238M1", "Internal Shields Bulk u238 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapInternalco60M2 = new TH1D("hAdapInternalco60M2", "Internal Shields Bulk co60 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapInternalk40M2 = new TH1D("hAdapInternalk40M2", "Internal Shields Bulk k40 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapInternalth232M2 = new TH1D("hAdapInternalth232M2", "Internal Shields Bulk th232 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapInternalu238M2 = new TH1D("hAdapInternalu238M2", "Internal Shields Bulk u238 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);

////////////// OVC M1 and M2
  hAdapOVCco60M1      = new TH1D("hAdapOVCco60M1",   "OVC Bulk co60 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapOVCk40M1       = new TH1D("hAdapOVCk40M1",    "OVC Bulk k40 M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapOVCth232M1     = new TH1D("hAdapOVCth232M1",  "OVC Bulk th232 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapOVCu238M1      = new TH1D("hAdapOVCu238M1",   "OVC Bulk u238 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapOVCbi207M1     = new TH1D("hAdapOVCbi207M1",  "OVC Bulk bi207 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  

  hAdapOVCco60M2      = new TH1D("hAdapOVCco60M2",   "OVC Bulk co60 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapOVCk40M2       = new TH1D("hAdapOVCk40M2",    "OVC Bulk k40 M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapOVCth232M2     = new TH1D("hAdapOVCth232M2",  "OVC Bulk th232 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapOVCu238M2      = new TH1D("hAdapOVCu238M2",   "OVC Bulk u238 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapOVCbi207M2     = new TH1D("hAdapOVCbi207M2",  "OVC Bulk bi207 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  

/////////// External Sources M1 and M2
  hAdapExtPbbi210M1 = new TH1D("hAdapExtPbbi210M1", "External Lead Bulk bi210 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapExtPbk40M1 = new TH1D("hAdapExtPbk40M1", "External Lead Bulk k40 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapExtPbth232M1 = new TH1D("hAdapExtPbth232M1", "External Lead Bulk th232 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapExtPbu238M1 = new TH1D("hAdapExtPbu238M1", "External Lead Bulk u238 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapExtPbpb210M1 = new TH1D("hAdapExtPbpb210M1", "External Lead Bulk pb210 M1", dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapExtPbbi210M2 = new TH1D("hAdapExtPbbi210M2", "External Lead Bulk bi210 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapExtPbk40M2 = new TH1D("hAdapExtPbk40M2", "External Lead Bulk k40 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapExtPbth232M2 = new TH1D("hAdapExtPbth232M2", "External Lead Bulk th232 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapExtPbu238M2 = new TH1D("hAdapExtPbu238M2", "External Lead Bulk u238 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapExtPbpb210M2 = new TH1D("hAdapExtPbpb210M2", "External Lead Bulk pb210 M2", dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapExtMuonM1 = new TH1D("hAdapExtMuonM1", "hAdapExtMuonM1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_th232spotM1 = new TH1D("hAdapCuBox_th232spotM1", "hAdapCuBox_th232spotM1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBox_k40spotM1 = new TH1D("hAdapCuBox_k40spotM1", "hAdapCuBox_k40spotM1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapBotExtPb_k40spotM1 = new TH1D("hAdapBotExtPb_k40spotM1", "hAdapBotExtPb_k40spotM1", dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapExtMuonM2 = new TH1D("hAdapExtMuonM2", "hAdapExtMuonM2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_th232spotM2 = new TH1D("hAdapCuBox_th232spotM2", "hAdapCuBox_th232spotM2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBox_k40spotM2 = new TH1D("hAdapCuBox_k40spotM2", "hAdapCuBox_k40spotM2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapBotExtPb_k40spotM2 = new TH1D("hAdapBotExtPb_k40spotM2", "hAdapBotExtPb_k40spotM2", dAdaptiveBinsM2, dAdaptiveArrayM2);

  hChiSquaredProgressM1 = new TH1D("hChiSquaredProgressM1", "Chi Squared Progress", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hChiSquaredProgressM2 = new TH1D("hChiSquaredProgressM2", "Chi Squared Progress", dAdaptiveBinsM2, dAdaptiveArrayM2);
}

