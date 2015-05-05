#include "TBkgModelSource.hh"

using namespace std;
using std::cout;
using std::endl;
using std::map;
using std::vector;

ClassImp(TBkgModelSource)

TBkgModelSource::TBkgModelSource()
{
	
}

// Initialize datasets
TBkgModelSource::TBkgModelSource(double fFitMin, double fFitMax, int fBinBase, int fDataSet)
{

  // Data directories depending on QCC/local
  // dDataDir =  "/Users/brian/macros/Simulations/Production/";
  dDataDir =  "/Users/brian/macros/CUOREZ/Bkg";
   // dDataDir = "/cuore/user/zhubrian/CUORE0/scratch/v2.30";
  dMCDir = "/Users/brian/macros/Simulations/Production";
   // dMCDir = "/cuore/user/zhubrian/MC/Bkg";
  // dSaveDir = "/Users/brian/Dropbox/code/Fitting";
   // dSaveDir = "/cuore/user/zhubrian/code/Fitting";

  dSecToYears = 1./(60*60*24*365);

  // Bin size (keV) -- base binning is 1 keV
  dBinSize = 1; 
  // Histogram range - from 0 to 10 MeV
  dMinEnergy = 0.;
  dMaxEnergy = 10000.;
  dBinBase = fBinBase;
  dDataSet = fDataSet;

  // Fitting range
  dFitMin = fFitMin;
  dFitMax = fFitMax;

  dNBins = (dMaxEnergy - dMinEnergy)/ dBinSize;  

  // Create Data Histograms
  fDataHistoTot     = new TH1D("fDataHistoTot",  "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM1      = new TH1D("fDataHistoM1",   "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM2      = new TH1D("fDataHistoM2",   "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM2Sum   = new TH1D("fDataHistoM2Sum",   "", dNBins, dMinEnergy, dMaxEnergy);


  // Data cuts 
  qtree = new TChain("qtree");
  // base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
  // base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
  base_cut = base_cut && "Channel != 1 && Channel != 10"; 

  // New PSA cuts
  base_cut = base_cut && "NormBaselineSlope < 4.8 && NormBaselineSlope > -4";
  base_cut = base_cut && "NormRiseTime < 4.8 && NormRiseTime > -4";
  base_cut = base_cut && "NormDecayTime < 4.8 && NormDecayTime > -4";
  base_cut = base_cut && "NormDelay < 4.8 && NormDelay > -4";
  base_cut = base_cut && "NormTVL < 5.3 && NormTVL > -6";
  base_cut = base_cut && "NormTVR < 5.3 && NormTVR > -6";

  dBaseBinSize = dBinSize*1;


  TBkgModelSource::LoadData();

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

  dDataIntegralTot = qtree->GetEntries();

  dDataIntegralM1 = fAdapDataHistoM1->Integral("width");
  dDataIntegralM2 = fAdapDataHistoM2->Integral("width");
  dDataIntegralM2Sum = fAdapDataHistoM2Sum->Integral("width");

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

  // Create Model Histograms
  TBkgModelSource::CreateModelHistograms();

  TBkgModelSource::LoadSources();

}


TBkgModelSource::~TBkgModelSource()
{}


// Returns vector of bin low-edge for adaptive binning
vector<double> TBkgModelSource::AdaptiveBinning(TH1D *h1, int dBinBase)
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
    {
     dBinArrayThing.push_back(h1->GetXaxis()->GetBinLowEdge(h1->FindBin(4170)));
     // Reset everything
     j = 0;
     dDummyFill = 0;
     dDummy = 0;
     i = i+(4350-4170)/dBaseBinSize;
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

void TBkgModelSource::CreateModelHistograms()
{
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
  hTeO2Sxu238M1_100      = new TH1D("hTeO2Sxu238M1_100",   "hTeO2Sxu238M1_100",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxth232M1_100     = new TH1D("hTeO2Sxth232M1_100",  "hTeO2Sxth232M1_100",  dNBins, dMinEnergy, dMaxEnergy);
  hTeO2Sxpb210M1_100     = new TH1D("hTeO2Sxpb210M1_100",  "hTeO2Sxpb210M1_100",  dNBins, dMinEnergy, dMaxEnergy);


  hTeO2th232onlyM1      = new TH1D("hTeO2th232onlyM1", "hTeO2th232onlyM1",     dNBins, dMinEnergy, dMaxEnergy);
  hTeO2ra228pb208M1     = new TH1D("hTeO2ra228pb208M1", "hTeO2ra228pb208M1",   dNBins, dMinEnergy, dMaxEnergy);
  hTeO2th230onlyM1      = new TH1D("hTeO2th230onlyM1", "hTeO2th230onlyM1",     dNBins, dMinEnergy, dMaxEnergy);

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
  hTeO2th230onlyM2      = new TH1D("hTeO2th230onlyM2", "hTeO2th230onlyM2",     dNBins, dMinEnergy, dMaxEnergy);

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
  hTeO2th230onlyM2Sum      = new TH1D("hTeO2th230onlyM2Sum", "hTeO2th230onlyM2Sum",     dNBins, dMinEnergy, dMaxEnergy);

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

/////////// Fudge Factors
  hOVC804M1 = new TH1D("hOVC804M1", "hOVC804M1", dNBins, dMinEnergy, dMaxEnergy);
  hOVC1063M1 = new TH1D("hOVC1063M1", "hOVC1063M1", dNBins, dMinEnergy, dMaxEnergy);
  hPbRom804M1 = new TH1D("hPbRom804M1", "hPbRom804M1", dNBins, dMinEnergy, dMaxEnergy);
  hPbRom1063M1 = new TH1D("hPbRom1063M1", "hPbRom1063M1", dNBins, dMinEnergy, dMaxEnergy);

  hOVC804M2 = new TH1D("hOVC804M2", "hOVC804M2", dNBins, dMinEnergy, dMaxEnergy);
  hOVC1063M2 = new TH1D("hOVC1063M2", "hOVC1063M2", dNBins, dMinEnergy, dMaxEnergy);
  hPbRom804M2 = new TH1D("hPbRom804M2", "hPbRom804M2", dNBins, dMinEnergy, dMaxEnergy);
  hPbRom1063M2 = new TH1D("hPbRom1063M2", "hPbRom1063M2", dNBins, dMinEnergy, dMaxEnergy);

  hOVC804M2Sum = new TH1D("hOVC804M2Sum", "hOVC804M2Sum", dNBins, dMinEnergy, dMaxEnergy);
  hOVC1063M2Sum = new TH1D("hOVC1063M2Sum", "hOVC1063M2Sum", dNBins, dMinEnergy, dMaxEnergy);
  hPbRom804M2Sum = new TH1D("hPbRom804M2Sum", "hPbRom804M2Sum", dNBins, dMinEnergy, dMaxEnergy);
  hPbRom1063M2Sum = new TH1D("hPbRom1063M2Sum", "hPbRom1063M2Sum", dNBins, dMinEnergy, dMaxEnergy);


//////////////// Adaptive binned histograms
////////// Total Adaptive binning histograms M1
  fModelTotAdapM1      = new TH1D("fModelTotAdapM1",      "Total PDF M1", dAdaptiveBinsM1, dAdaptiveArrayM1);  
  fModelTotAdapthM1    = new TH1D("fModelTotAdapthM1",    "Total Bulk th232",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapuM1     = new TH1D("fModelTotAdapuM1",     "Total Bulk u238",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapkM1     = new TH1D("fModelTotAdapkM1",     "Total Bulk k40",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapcoM1    = new TH1D("fModelTotAdapcoM1",    "Total Bulk co60",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapmnM1    = new TH1D("fModelTotAdapmnM1",    "Total Bulk mn54",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  fModelTotAdapNDBDM1  = new TH1D("fModelTotAdapNDBDM1",  "Total Bulk NDBD",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdap2NDBDM1 = new TH1D("fModelTotAdap2NDBDM1", "Total Bulk 2NDBD",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapbiM1    = new TH1D("fModelTotAdapbiM1",    "Total Bulk bi207",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapbi2M1   = new TH1D("fModelTotAdapbi2M1",   "Total Bulk bi210",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapptM1    = new TH1D("fModelTotAdapptM1",    "Total Bulk pt190",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdappbM1    = new TH1D("fModelTotAdappbM1",    "Total Bulk pb210",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapcsM1    = new TH1D("fModelTotAdapcsM1",    "Total Bulk cs137",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapco2M1   = new TH1D("fModelTotAdapco2M1",   "Total Bulk co58",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapteo2M1  = new TH1D("fModelTotAdapteo2M1",  "Total Bulk TeO2",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  fModelTotAdapSthM1   = new TH1D("fModelTotAdapSthM1",   "Total Surface th232",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapSuM1    = new TH1D("fModelTotAdapSuM1",    "Total Surface u238",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapSpbM1   = new TH1D("fModelTotAdapSpbM1",   "Total Surface pb210",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelTotAdapSpoM1   = new TH1D("fModelTotAdapSpoM1",   "Total Surface po210",  dAdaptiveBinsM1, dAdaptiveArrayM1);

  fModelTotAdapExtM1   = new TH1D("fModelTotAdapExtM1", "Total External", dAdaptiveBinsM1, dAdaptiveArrayM1);

  fModelTotAdapFudgeM1 = new TH1D("fModelTotAdapFudgeM1", "Total Fudge Factors", dAdaptiveBinsM1, dAdaptiveArrayM1);

/////////// Total Adaptive binning histograms M2
  fModelTotAdapM2      = new TH1D("fModelTotAdapM2",      "Total PDF M2", dAdaptiveBinsM2, dAdaptiveArrayM2);  
  fModelTotAdapthM2    = new TH1D("fModelTotAdapthM2",    "Total Bulk th232",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapuM2     = new TH1D("fModelTotAdapuM2",     "Total Bulk u238",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapkM2     = new TH1D("fModelTotAdapkM2",     "Total Bulk k40",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapcoM2    = new TH1D("fModelTotAdapcoM2",    "Total Bulk co60",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapmnM2    = new TH1D("fModelTotAdapmnM2",    "Total Bulk mn54",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  fModelTotAdapNDBDM2  = new TH1D("fModelTotAdapNDBDM2",  "Total Bulk NDBD",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdap2NDBDM2 = new TH1D("fModelTotAdap2NDBDM2", "Total Bulk 2NDBD",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapbiM2    = new TH1D("fModelTotAdapbiM2",    "Total Bulk bi207",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapbi2M2   = new TH1D("fModelTotAdapbi2M2",   "Total Bulk bi210",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapptM2    = new TH1D("fModelTotAdapotM2",    "Total Bulk pt190",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdappbM2    = new TH1D("fModelTotAdappbM2",    "Total Bulk pb210",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapcsM2    = new TH1D("fModelTotAdapcsM2",    "Total Bulk cs137",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapco2M2   = new TH1D("fModelTotAdapco2M2",   "Total Bulk co58",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapteo2M2  = new TH1D("fModelTotAdapteo2M2",  "Total Bulk TeO2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  fModelTotAdapSthM2   = new TH1D("fModelTotAdapSthM2",   "Total Surface th232",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapSuM2    = new TH1D("fModelTotAdapSuM2",    "Total Surface u238",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapSpbM2   = new TH1D("fModelTotAdapSpbM2",   "Total Surface pb210",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  fModelTotAdapSpoM2   = new TH1D("fModelTotAdapSpoM2",   "Total Surface po210",  dAdaptiveBinsM2, dAdaptiveArrayM2);

  fModelTotAdapExtM2   = new TH1D("fModelTotAdapExtM2", "Total External", dAdaptiveBinsM2, dAdaptiveArrayM2);

  fModelTotAdapFudgeM2 = new TH1D("fModelTotAdapFudgeM2", "Total Fudge Factors", dAdaptiveBinsM2, dAdaptiveArrayM2);

/////////// Total Adaptive binning histograms M2Sum

  fModelTotAdapM2Sum      = new TH1D("fModelTotAdapM2Sum",      "Total PDF M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  fModelTotAdapthM2Sum    = new TH1D("fModelTotAdapthM2Sum",    "Total Bulk th232",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapuM2Sum     = new TH1D("fModelTotAdapuM2Sum",     "Total Bulk u238",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapkM2Sum     = new TH1D("fModelTotAdapkM2Sum",     "Total Bulk k40",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapcoM2Sum    = new TH1D("fModelTotAdapcoM2Sum",    "Total Bulk co60",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapmnM2Sum    = new TH1D("fModelTotAdapmnM2Sum",    "Total Bulk mn54",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  fModelTotAdapNDBDM2Sum  = new TH1D("fModelTotAdapNDBDM2Sum",  "Total Bulk NDBD",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdap2NDBDM2Sum = new TH1D("fModelTotAdap2NDBDM2Sum", "Total Bulk 2NDBD",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapbiM2Sum    = new TH1D("fModelTotAdapbiM2Sum",    "Total Bulk bi207",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapbi2M2Sum   = new TH1D("fModelTotAdapbi2M2Sum",   "Total Bulk bi210",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapptM2Sum    = new TH1D("fModelTotAdapotM2Sum",    "Total Bulk pt190",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdappbM2Sum    = new TH1D("fModelTotAdappbM2Sum",    "Total Bulk pb210",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapcsM2Sum    = new TH1D("fModelTotAdapcsM2Sum",    "Total Bulk cs137",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapco2M2Sum   = new TH1D("fModelTotAdapco2M2Sum",   "Total Bulk co58",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapteo2M2Sum  = new TH1D("fModelTotAdapteo2M2Sum",  "Total Bulk TeO2",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  fModelTotAdapSthM2Sum   = new TH1D("fModelTotAdapSthM2Sum",   "Total Surface th232",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapSuM2Sum    = new TH1D("fModelTotAdapSuM2Sum",    "Total Surface u238",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapSpbM2Sum   = new TH1D("fModelTotAdapSpbM2Sum",   "Total Surface pb210",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  fModelTotAdapSpoM2Sum   = new TH1D("fModelTotAdapSpoM2Sum",   "Total Surface po210",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  fModelTotAdapExtM2Sum   = new TH1D("fModelTotAdapExtM2Sum", "Total External", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  fModelTotAdapFudgeM2Sum = new TH1D("fModelTotAdapFudgeM2Sum", "Total Fudge Factors", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

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
  hAdapTeO2th230M1     = new TH1D("hAdapTeO2th230M1",  "TeO2 Bulk th230M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapTeO2u234M1      = new TH1D("hAdapTeO2u234M1",   "TeO2 Bulk u234 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

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
  hAdapTeO2th230onlyM1      = new TH1D("hAdapTeO2th230onlyM1",   "TeO2 Bulk th230 only M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);

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
  hAdapTeO2th230onlyM2      = new TH1D("hAdapTeO2th230onlyM2",   "TeO2 Bulk th230 only M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);

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
  hAdapTeO2th230onlyM2Sum      = new TH1D("hAdapTeO2th230onlyM2Sum",   "TeO2 Bulk th230 only M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

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

////////// Frame M1 and M2
  hAdapCuFrameco58M1      = new TH1D("hAdapCuFrameco58M1",   "CuFrame Bulk co58 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameco60M1      = new TH1D("hAdapCuFrameco60M1",   "CuFrame Bulk co60 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFramecs137M1     = new TH1D("hAdapCuFramecs137M1",  "CuFrame Bulk cs137 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFramek40M1       = new TH1D("hAdapCuFramek40M1",    "CuFrame Bulk k40 M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFramemn54M1      = new TH1D("hAdapCuFramemn54M1",   "CuFrame Bulk mn54 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFramepb210M1     = new TH1D("hAdapCuFramepb210M1",  "CuFrame Bulk pb210 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameth232M1     = new TH1D("hAdapCuFrameth232M1",  "CuFrame Bulk th232 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapCuFrameu238M1      = new TH1D("hAdapCuFrameu238M1",   "CuFrame Bulk u238 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuFrameSth232M1_1     = new TH1D("hAdapCuFrameSth232M1_1",    "CuFrame S th232 M1 1 #mum",     dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSu238M1_1      = new TH1D("hAdapCuFrameSu238M1_1",     "CuFrame S u238 M1 1 #mum",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxpb210M1_001  = new TH1D("hAdapCuFrameSxpb210M1_001", "CuFrame Sx pb210 M1 0.01 #mum",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxpb210M1_01   = new TH1D("hAdapCuFrameSxpb210M1_01",  "CuFrame Sx pb210 M1 0.1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxpb210M1_1    = new TH1D("hAdapCuFrameSxpb210M1_1",   "CuFrame Sx pb210 M1 1 #mum",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxpb210M1_10   = new TH1D("hAdapCuFrameSxpb210M1_10",  "CuFrame Sx pb210 M1 10 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxth232M1_001  = new TH1D("hAdapCuFrameSxth232M1_001", "CuFrame Sx th232 M1 0.01 #mum",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxth232M1_01   = new TH1D("hAdapCuFrameSxth232M1_01",  "CuFrame Sx th232 M1 0.1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxth232M1_1    = new TH1D("hAdapCuFrameSxth232M1_1",   "CuFrame Sx th232 M1 1 #mum",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxth232M1_10   = new TH1D("hAdapCuFrameSxth232M1_10",  "CuFrame Sx th232 M1 10 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxu238M1_001   = new TH1D("hAdapCuFrameSxu238M1_001",  "CuFrame Sx u238 M1 0.01 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxu238M1_01    = new TH1D("hAdapCuFrameSxu238M1_01",   "CuFrame Sx u238 M1 0.1 #mum",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxu238M1_1     = new TH1D("hAdapCuFrameSxu238M1_1",    "CuFrame Sx u238 M1 1 #mum",     dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameSxu238M1_10    = new TH1D("hAdapCuFrameSxu238M1_10",   "CuFrame Sx u238 M1 10 #mum",    dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuFrameco58M2      = new TH1D("hAdapCuFrameco58M2",   "CuFrame Bulk co58 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameco60M2      = new TH1D("hAdapCuFrameco60M2",   "CuFrame Bulk co60 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFramecs137M2     = new TH1D("hAdapCuFramecs137M2",  "CuFrame Bulk cs137 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFramek40M2       = new TH1D("hAdapCuFramek40M2",    "CuFrame Bulk k40 M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFramemn54M2      = new TH1D("hAdapCuFramemn54M2",   "CuFrame Bulk mn54 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFramepb210M2     = new TH1D("hAdapCuFramepb210M2",  "CuFrame Bulk pb210 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameth232M2     = new TH1D("hAdapCuFrameth232M2",  "CuFrame Bulk th232 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapCuFrameu238M2      = new TH1D("hAdapCuFrameu238M2",   "CuFrame Bulk u238 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapCuFrameSth232M2_1     = new TH1D("hAdapCuFrameSth232M2_1",    "CuFrame S th232 M2 1 #mum",     dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSu238M2_1      = new TH1D("hAdapCuFrameSu238M2_1",     "CuFrame S u238 M2 1 #mum",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxpb210M2_001  = new TH1D("hAdapCuFrameSxpb210M2_001", "CuFrame Sx pb210 M2 0.01 #mum",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxpb210M2_01   = new TH1D("hAdapCuFrameSxpb210M2_01",  "CuFrame Sx pb210 M2 0.1 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxpb210M2_1    = new TH1D("hAdapCuFrameSxpb210M2_1",   "CuFrame Sx pb210 M2 1 #mum",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxpb210M2_10   = new TH1D("hAdapCuFrameSxpb210M2_10",  "CuFrame Sx pb210 M2 10 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxth232M2_001  = new TH1D("hAdapCuFrameSxth232M2_001", "CuFrame Sx th232 M2 0.01 #mum",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxth232M2_01   = new TH1D("hAdapCuFrameSxth232M2_01",  "CuFrame Sx th232 M2 0.1 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxth232M2_1    = new TH1D("hAdapCuFrameSxth232M2_1",   "CuFrame Sx th232 M2 1 #mum",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxth232M2_10   = new TH1D("hAdapCuFrameSxth232M2_10",  "CuFrame Sx th232 M2 10 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxu238M2_001   = new TH1D("hAdapCuFrameSxu238M2_001",  "CuFrame Sx u238 M2 0.01 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxu238M2_01    = new TH1D("hAdapCuFrameSxu238M2_01",   "CuFrame Sx u238 M2 0.1 #mum",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxu238M2_1     = new TH1D("hAdapCuFrameSxu238M2_1",    "CuFrame Sx u238 M2 1 #mum",     dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameSxu238M2_10    = new TH1D("hAdapCuFrameSxu238M2_10",   "CuFrame Sx u238 M2 10 #mum",    dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapCuFrameco58M2Sum      = new TH1D("hAdapCuFrameco58M2Sum",   "CuFrame Bulk co58 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameco60M2Sum      = new TH1D("hAdapCuFrameco60M2Sum",   "CuFrame Bulk co60 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFramecs137M2Sum     = new TH1D("hAdapCuFramecs137M2Sum",  "CuFrame Bulk cs137 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFramek40M2Sum       = new TH1D("hAdapCuFramek40M2Sum",    "CuFrame Bulk k40 M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFramemn54M2Sum      = new TH1D("hAdapCuFramemn54M2Sum",   "CuFrame Bulk mn54 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFramepb210M2Sum     = new TH1D("hAdapCuFramepb210M2Sum",  "CuFrame Bulk pb210 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameth232M2Sum     = new TH1D("hAdapCuFrameth232M2Sum",  "CuFrame Bulk th232 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapCuFrameu238M2Sum      = new TH1D("hAdapCuFrameu238M2Sum",   "CuFrame Bulk u238 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapCuFrameSth232M2Sum_1     = new TH1D("hAdapCuFrameSth232M2Sum_1",    "CuFrame S th232 M2Sum 1 #mum",     dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSu238M2Sum_1      = new TH1D("hAdapCuFrameSu238M2Sum_1",     "CuFrame S u238 M2Sum 1 #mum",      dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxpb210M2Sum_001  = new TH1D("hAdapCuFrameSxpb210M2Sum_001", "CuFrame Sx pb210 M2Sum 0.01 #mum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxpb210M2Sum_01   = new TH1D("hAdapCuFrameSxpb210M2Sum_01",  "CuFrame Sx pb210 M2Sum 0.1 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxpb210M2Sum_1    = new TH1D("hAdapCuFrameSxpb210M2Sum_1",   "CuFrame Sx pb210 M2Sum 1 #mum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxpb210M2Sum_10   = new TH1D("hAdapCuFrameSxpb210M2Sum_10",  "CuFrame Sx pb210 M2Sum 10 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxth232M2Sum_001  = new TH1D("hAdapCuFrameSxth232M2Sum_001", "CuFrame Sx th232 M2Sum 0.01 #mum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxth232M2Sum_01   = new TH1D("hAdapCuFrameSxth232M2Sum_01",  "CuFrame Sx th232 M2Sum 0.1 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxth232M2Sum_1    = new TH1D("hAdapCuFrameSxth232M2Sum_1",   "CuFrame Sx th232 M2Sum 1 #mum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxth232M2Sum_10   = new TH1D("hAdapCuFrameSxth232M2Sum_10",  "CuFrame Sx th232 M2Sum 10 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxu238M2Sum_001   = new TH1D("hAdapCuFrameSxu238M2Sum_001",  "CuFrame Sx u238 M2Sum 0.01 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxu238M2Sum_01    = new TH1D("hAdapCuFrameSxu238M2Sum_01",   "CuFrame Sx u238 M2Sum 0.1 #mum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxu238M2Sum_1     = new TH1D("hAdapCuFrameSxu238M2Sum_1",    "CuFrame Sx u238 M2Sum 1 #mum",     dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuFrameSxu238M2Sum_10    = new TH1D("hAdapCuFrameSxu238M2Sum_10",   "CuFrame Sx u238 M2Sum 10 #mum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

////////// CuBox (TShield) M1 and M2
  hAdapCuBoxco58M1      = new TH1D("hAdapCuBoxco58M1",   "CuBox Bulk co58 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxco60M1      = new TH1D("hAdapCuBoxco60M1",   "CuBox Bulk co60 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxcs137M1     = new TH1D("hAdapCuBoxcs137M1",  "CuBox Bulk cs137 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxk40M1       = new TH1D("hAdapCuBoxk40M1",    "CuBox Bulk k40 M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxmn54M1      = new TH1D("hAdapCuBoxmn54M1",   "CuBox Bulk mn54 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxpb210M1     = new TH1D("hAdapCuBoxpb210M1",  "CuBox Bulk pb210 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxth232M1     = new TH1D("hAdapCuBoxth232M1",  "CuBox Bulk th232 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapCuBoxu238M1      = new TH1D("hAdapCuBoxu238M1",   "CuBox Bulk u238 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuBoxSth232M1_1     = new TH1D("hAdapCuBoxSth232M1_1",    "CuBox S th232 M1 1 #mum",     dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSu238M1_1      = new TH1D("hAdapCuBoxSu238M1_1",     "CuBox S u238 M1 1 #mum",      dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxpb210M1_001  = new TH1D("hAdapCuBoxSxpb210M1_001", "CuBox Sx pb210 M1 0.01 #mum",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxpb210M1_01   = new TH1D("hAdapCuBoxSxpb210M1_01",  "CuBox Sx pb210 M1 0.1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxpb210M1_1    = new TH1D("hAdapCuBoxSxpb210M1_1",   "CuBox Sx pb210 M1 1 #mum",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxpb210M1_10   = new TH1D("hAdapCuBoxSxpb210M1_10",  "CuBox Sx pb210 M1 10 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxth232M1_001  = new TH1D("hAdapCuBoxSxth232M1_001", "CuBox Sx th232 M1 0.01 #mum",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxth232M1_01   = new TH1D("hAdapCuBoxSxth232M1_01",  "CuBox Sx th232 M1 0.1 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxth232M1_1    = new TH1D("hAdapCuBoxSxth232M1_1",   "CuBox Sx th232 M1 1 #mum",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxth232M1_10   = new TH1D("hAdapCuBoxSxth232M1_10",  "CuBox Sx th232 M1 10 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxu238M1_001   = new TH1D("hAdapCuBoxSxu238M1_001",  "CuBox Sx u238 M1 0.01 #mum",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxu238M1_01    = new TH1D("hAdapCuBoxSxu238M1_01",   "CuBox Sx u238 M1 0.1 #mum",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxu238M1_1     = new TH1D("hAdapCuBoxSxu238M1_1",    "CuBox Sx u238 M1 1 #mum",     dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxSxu238M1_10    = new TH1D("hAdapCuBoxSxu238M1_10",   "CuBox Sx u238 M1 10 #mum",    dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuBoxco58M2      = new TH1D("hAdapCuBoxco58M2",   "CuBox Bulk co58 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxco60M2      = new TH1D("hAdapCuBoxco60M2",   "CuBox Bulk co60 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxcs137M2     = new TH1D("hAdapCuBoxcs137M2",  "CuBox Bulk cs137 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxk40M2       = new TH1D("hAdapCuBoxk40M2",    "CuBox Bulk k40 M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxmn54M2      = new TH1D("hAdapCuBoxmn54M2",   "CuBox Bulk mn54 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxpb210M2     = new TH1D("hAdapCuBoxpb210M2",  "CuBox Bulk pb210 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxth232M2     = new TH1D("hAdapCuBoxth232M2",  "CuBox Bulk th232 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapCuBoxu238M2      = new TH1D("hAdapCuBoxu238M2",   "CuBox Bulk u238 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapCuBoxSth232M2_1     = new TH1D("hAdapCuBoxSth232M2_1",    "CuBox S th232 M2 1 #mum",     dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSu238M2_1      = new TH1D("hAdapCuBoxSu238M2_1",     "CuBox S u238 M2 1 #mum",      dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxpb210M2_001  = new TH1D("hAdapCuBoxSxpb210M2_001", "CuBox Sx pb210 M2 0.01 #mum",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxpb210M2_01   = new TH1D("hAdapCuBoxSxpb210M2_01",  "CuBox Sx pb210 M2 0.1 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxpb210M2_1    = new TH1D("hAdapCuBoxSxpb210M2_1",   "CuBox Sx pb210 M2 1 #mum",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxpb210M2_10   = new TH1D("hAdapCuBoxSxpb210M2_10",  "CuBox Sx pb210 M2 10 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxth232M2_001  = new TH1D("hAdapCuBoxSxth232M2_001", "CuBox Sx th232 M2 0.01 #mum",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxth232M2_01   = new TH1D("hAdapCuBoxSxth232M2_01",  "CuBox Sx th232 M2 0.1 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxth232M2_1    = new TH1D("hAdapCuBoxSxth232M2_1",   "CuBox Sx th232 M2 1 #mum",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxth232M2_10   = new TH1D("hAdapCuBoxSxth232M2_10",  "CuBox Sx th232 M2 10 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxu238M2_001   = new TH1D("hAdapCuBoxSxu238M2_001",  "CuBox Sx u238 M2 0.01 #mum",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxu238M2_01    = new TH1D("hAdapCuBoxSxu238M2_01",   "CuBox Sx u238 M2 0.1 #mum",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxu238M2_1     = new TH1D("hAdapCuBoxSxu238M2_1",    "CuBox Sx u238 M2 1 #mum",     dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxSxu238M2_10    = new TH1D("hAdapCuBoxSxu238M2_10",   "CuBox Sx u238 M2 10 #mum",    dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapCuBoxco58M2Sum      = new TH1D("hAdapCuBoxco58M2Sum",   "CuBox Bulk co58 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxco60M2Sum      = new TH1D("hAdapCuBoxco60M2Sum",   "CuBox Bulk co60 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxcs137M2Sum     = new TH1D("hAdapCuBoxcs137M2Sum",  "CuBox Bulk cs137 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxk40M2Sum       = new TH1D("hAdapCuBoxk40M2Sum",    "CuBox Bulk k40 M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxmn54M2Sum      = new TH1D("hAdapCuBoxmn54M2Sum",   "CuBox Bulk mn54 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxpb210M2Sum     = new TH1D("hAdapCuBoxpb210M2Sum",  "CuBox Bulk pb210 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxth232M2Sum     = new TH1D("hAdapCuBoxth232M2Sum",  "CuBox Bulk th232 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapCuBoxu238M2Sum      = new TH1D("hAdapCuBoxu238M2Sum",   "CuBox Bulk u238 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

  hAdapCuBoxSth232M2Sum_1     = new TH1D("hAdapCuBoxSth232M2Sum_1",    "CuBox S th232 M2Sum 1 #mum",     dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSu238M2Sum_1      = new TH1D("hAdapCuBoxSu238M2Sum_1",     "CuBox S u238 M2Sum 1 #mum",      dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxpb210M2Sum_001  = new TH1D("hAdapCuBoxSxpb210M2Sum_001", "CuBox Sx pb210 M2Sum 0.01 #mum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxpb210M2Sum_01   = new TH1D("hAdapCuBoxSxpb210M2Sum_01",  "CuBox Sx pb210 M2Sum 0.1 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxpb210M2Sum_1    = new TH1D("hAdapCuBoxSxpb210M2Sum_1",   "CuBox Sx pb210 M2Sum 1 #mum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxpb210M2Sum_10   = new TH1D("hAdapCuBoxSxpb210M2Sum_10",  "CuBox Sx pb210 M2Sum 10 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxth232M2Sum_001  = new TH1D("hAdapCuBoxSxth232M2Sum_001", "CuBox Sx th232 M2Sum 0.01 #mum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxth232M2Sum_01   = new TH1D("hAdapCuBoxSxth232M2Sum_01",  "CuBox Sx th232 M2Sum 0.1 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxth232M2Sum_1    = new TH1D("hAdapCuBoxSxth232M2Sum_1",   "CuBox Sx th232 M2Sum 1 #mum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxth232M2Sum_10   = new TH1D("hAdapCuBoxSxth232M2Sum_10",  "CuBox Sx th232 M2Sum 10 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxu238M2Sum_001   = new TH1D("hAdapCuBoxSxu238M2Sum_001",  "CuBox Sx u238 M2Sum 0.01 #mum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxu238M2Sum_01    = new TH1D("hAdapCuBoxSxu238M2Sum_01",   "CuBox Sx u238 M2Sum 0.1 #mum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxu238M2Sum_1     = new TH1D("hAdapCuBoxSxu238M2Sum_1",    "CuBox Sx u238 M2Sum 1 #mum",     dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapCuBoxSxu238M2Sum_10    = new TH1D("hAdapCuBoxSxu238M2Sum_10",   "CuBox Sx u238 M2Sum 10 #mum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

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
  hAdap50mKco58M1      = new TH1D("hAdap50mKco58M1",   "50mK Bulk co58 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKco60M1      = new TH1D("hAdap50mKco60M1",   "50mK Bulk co60 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKcs137M1     = new TH1D("hAdap50mKcs137M1",  "50mK Bulk cs137 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKk40M1       = new TH1D("hAdap50mKk40M1",    "50mK Bulk k40 M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKmn54M1      = new TH1D("hAdap50mKmn54M1",   "50mK Bulk mn54 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKpb210M1     = new TH1D("hAdap50mKpb210M1",  "50mK Bulk pb210 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap50mKth232M1     = new TH1D("hAdap50mKth232M1",  "50mK Bulk th232 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdap50mKu238M1      = new TH1D("hAdap50mKu238M1",   "50mK Bulk u238 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdap50mKco58M2      = new TH1D("hAdap50mKco58M2",   "50mK Bulk co58 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap50mKco60M2      = new TH1D("hAdap50mKco60M2",   "50mK Bulk co60 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap50mKcs137M2     = new TH1D("hAdap50mKcs137M2",  "50mK Bulk cs137 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap50mKk40M2       = new TH1D("hAdap50mKk40M2",    "50mK Bulk k40 M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap50mKmn54M2      = new TH1D("hAdap50mKmn54M2",   "50mK Bulk mn54 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap50mKpb210M2     = new TH1D("hAdap50mKpb210M2",  "50mK Bulk pb210 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap50mKth232M2     = new TH1D("hAdap50mKth232M2",  "50mK Bulk th232 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdap50mKu238M2      = new TH1D("hAdap50mKu238M2",   "50mK Bulk u238 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdap50mKco58M2Sum      = new TH1D("hAdap50mKco58M2Sum",   "50mK Bulk co58 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap50mKco60M2Sum      = new TH1D("hAdap50mKco60M2Sum",   "50mK Bulk co60 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap50mKcs137M2Sum     = new TH1D("hAdap50mKcs137M2Sum",  "50mK Bulk cs137 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap50mKk40M2Sum       = new TH1D("hAdap50mKk40M2Sum",    "50mK Bulk k40 M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap50mKmn54M2Sum      = new TH1D("hAdap50mKmn54M2Sum",   "50mK Bulk mn54 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap50mKpb210M2Sum     = new TH1D("hAdap50mKpb210M2Sum",  "50mK Bulk pb210 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap50mKth232M2Sum     = new TH1D("hAdap50mKth232M2Sum",  "50mK Bulk th232 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdap50mKu238M2Sum      = new TH1D("hAdap50mKu238M2Sum",   "50mK Bulk u238 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

////////// 600mK M1 and M2
  hAdap600mKco60M1      = new TH1D("hAdap600mKco60M1",   "600mK Bulk co60 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap600mKk40M1       = new TH1D("hAdap600mKk40M1",    "600mK Bulk k40 M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap600mKth232M1     = new TH1D("hAdap600mKth232M1",  "600mK Bulk th232 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdap600mKu238M1      = new TH1D("hAdap600mKu238M1",   "600mK Bulk u238 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdap600mKco60M2      = new TH1D("hAdap600mKco60M2",   "600mK Bulk co60 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap600mKk40M2       = new TH1D("hAdap600mKk40M2",    "600mK Bulk k40 M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap600mKth232M2     = new TH1D("hAdap600mKth232M2",  "600mK Bulk th232 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdap600mKu238M2      = new TH1D("hAdap600mKu238M2",   "600mK Bulk u238 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  

  hAdap600mKco60M2Sum      = new TH1D("hAdap600mKco60M2Sum",   "600mK Bulk co60 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap600mKk40M2Sum       = new TH1D("hAdap600mKk40M2Sum",    "600mK Bulk k40 M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdap600mKth232M2Sum     = new TH1D("hAdap600mKth232M2Sum",  "600mK Bulk th232 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdap600mKu238M2Sum      = new TH1D("hAdap600mKu238M2Sum",   "600mK Bulk u238 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

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

////////////// Main bath M1 and M2
  hAdapMBco60M1      = new TH1D("hAdapMBco60M1",   "Main Bath Bulk co60 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapMBk40M1       = new TH1D("hAdapMBk40M1",    "Main Bath Bulk k40 M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapMBth232M1     = new TH1D("hAdapMBth232M1",  "Main Bath Bulk th232 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapMBu238M1      = new TH1D("hAdapMBu238M1",   "Main Bath Bulk u238 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapMBco60M2      = new TH1D("hAdapMBco60M2",   "Main Bath Bulk co60 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapMBk40M2       = new TH1D("hAdapMBk40M2",    "Main Bath Bulk k40 M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapMBth232M2     = new TH1D("hAdapMBth232M2",  "Main Bath Bulk th232 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapMBu238M2      = new TH1D("hAdapMBu238M2",   "Main Bath Bulk u238 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  

  hAdapMBco60M2Sum      = new TH1D("hAdapMBco60M2Sum",   "Main Bath Bulk co60 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapMBk40M2Sum       = new TH1D("hAdapMBk40M2Sum",    "Main Bath Bulk k40 M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapMBth232M2Sum     = new TH1D("hAdapMBth232M2Sum",  "Main Bath Bulk th232 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapMBu238M2Sum      = new TH1D("hAdapMBu238M2Sum",   "Main Bath Bulk u238 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  

/////////// Super Insulation M1 and M2
  hAdapSIk40M1       = new TH1D("hAdapSIk40M1",    "Super Insulation Bulk k40 M1",    dAdaptiveBinsM1, dAdaptiveArrayM1); 
  hAdapSIth232M1     = new TH1D("hAdapSIth232M1",  "Super Insulation Bulk th232 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapSIu238M1      = new TH1D("hAdapSIu238M1",   "Super Insulation Bulk u238 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapSIk40M2       = new TH1D("hAdapSIk40M2",    "Super Insulation Bulk k40 M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapSIth232M2     = new TH1D("hAdapSIth232M2",  "Super Insulation Bulk th232 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapSIu238M2      = new TH1D("hAdapSIu238M2",   "Super Insulation Bulk u238 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapSIk40M2Sum       = new TH1D("hAdapSIk40M2Sum",    "Super Insulation Bulk k40 M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapSIth232M2Sum     = new TH1D("hAdapSIth232M2Sum",  "Super Insulation Bulk th232 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapSIu238M2Sum      = new TH1D("hAdapSIu238M2Sum",   "Super Insulation Bulk u238 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);

//////////// IVC M1 and M2
  hAdapIVCco60M1      = new TH1D("hAdapIVCco60M1",   "IVC Bulk co60 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapIVCk40M1       = new TH1D("hAdapIVCk40M1",    "IVC Bulk k40 M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapIVCth232M1     = new TH1D("hAdapIVCth232M1",  "IVC Bulk th232 M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapIVCu238M1      = new TH1D("hAdapIVCu238M1",   "IVC Bulk u238 M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapIVCco60M2      = new TH1D("hAdapIVCco60M2",   "IVC Bulk co60 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapIVCk40M2       = new TH1D("hAdapIVCk40M2",    "IVC Bulk k40 M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapIVCth232M2     = new TH1D("hAdapIVCth232M2",  "IVC Bulk th232 M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapIVCu238M2      = new TH1D("hAdapIVCu238M2",   "IVC Bulk u238 M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  

  hAdapIVCco60M2Sum      = new TH1D("hAdapIVCco60M2Sum",   "IVC Bulk co60 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapIVCk40M2Sum       = new TH1D("hAdapIVCk40M2Sum",    "IVC Bulk k40 M2Sum",    dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapIVCth232M2Sum     = new TH1D("hAdapIVCth232M2Sum",  "IVC Bulk th232 M2Sum",  dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);  
  hAdapIVCu238M2Sum      = new TH1D("hAdapIVCu238M2Sum",   "IVC Bulk u238 M2Sum",   dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum); 

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
  hAdapOVC1063M1 = new TH1D("hAdapOVC1063M1", "hAdapOVC1063M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapPbRom804M1 = new TH1D("hAdapPbRom804M1", "hAdapPbRom804M1", dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapPbRom1063M1 = new TH1D("hAdapPbRom1063M1", "hAdapPbRom1063M1", dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapOVC804M2 = new TH1D("hAdapOVC804M2", "hAdapOVC804M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapOVC1063M2 = new TH1D("hAdapOVC1063M2", "hAdapOVC1063M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapPbRom804M2 = new TH1D("hAdapPbRom804M2", "hAdapPbRom804M2", dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapPbRom1063M2 = new TH1D("hAdapPbRom1063M2", "hAdapPbRom1063M2", dAdaptiveBinsM2, dAdaptiveArrayM2);

  hAdapOVC804M2Sum = new TH1D("hAdapOVC804M2Sum", "hAdapOVC804M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapOVC1063M2Sum = new TH1D("hAdapOVC1063M2Sum", "hAdapOVC1063M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapPbRom804M2Sum = new TH1D("hAdapPbRom804M2Sum", "hAdapPbRom804M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);
  hAdapPbRom1063M2Sum = new TH1D("hAdapPbRom1063M2Sum", "hAdapPbRom1063M2Sum", dAdaptiveBinsM2Sum, dAdaptiveArrayM2Sum);


}


void TBkgModelSource::LoadData()
{
  switch(dDataSet)
  { 
  case 1:
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2049.root", dDataDir.c_str()));   
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2061.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2064.root", dDataDir.c_str()));   
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2067.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2070.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2073.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2076.root", dDataDir.c_str())); 
    dLivetime = 6042498; // DR 1 
    cout << "Using Data Release 1" << endl;
  break;

  case 2:
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2079.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2085.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2088.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2091.root", dDataDir.c_str())); 
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2097.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2100.root", dDataDir.c_str())); 
    dLivetime = 9387524; // DR 2
    cout << "Using Data Release 2" << endl;
  break;

  case 3:
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2103.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2109.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2118.root", dDataDir.c_str()));
    qtree->Add(Form("%s/Unblinded/ReducedB-ds2124.root", dDataDir.c_str()));
    dLivetime = 7647908; // DR 3
    cout << "Using Data Release 3" << endl;
  break;

  // case -1:
  //   cout << "Using Toy data" << endl;
  //   bToyData = true;
  //   qtree->Add(Form("%s/Unblinded/ReducedB-ds*.root", dDataDir.c_str())); 
  //   dLivetime = 23077930+1511379+901597; // seconds of livetime (DR1 to DR3)
  // break;

  default:   
    qtree->Add(Form("%s/Unblinded/ReducedB-ds*.root", dDataDir.c_str())); 
    dLivetime = 23077930+1511379+901597; // seconds of livetime (DR1 to DR3)
    cout << "Using Total Dataset" << endl;
  }

  qtree->Project("fDataHistoTot", "RescaledEnergy", base_cut);
  qtree->Project("fDataHistoM1",  "RescaledEnergy", base_cut && "Multiplicity_Sync == 1");
  qtree->Project("fDataHistoM2",  "RescaledEnergy", base_cut && "Multiplicity_Sync == 2");
  qtree->Project("fDataHistoM2Sum",  "TotalEnergy_Sync", base_cut && "Multiplicity_Sync == 2");

  dLivetimeYr = dLivetime*dSecToYears;  
}

void TBkgModelSource::LoadSources()
{	
  // Loads PDFs from file
  cout << "Loading PDF Histograms from file" << endl;
  cout << "Directory " << Form("%s/", dMCDir.c_str()) << endl;

  fBulkInner = new TFile(Form("%s/OldProd/MCProduction_BulkInner_1keV.root", dMCDir.c_str()));
  fBulkInnerOld = new TFile(Form("%s/OldProd/MCProduction_BulkInner_1keV.root", dMCDir.c_str()));
  fBulkInnerM2Sum = new TFile(Form("%s/OldProd/MCProduction_BulkInnerM2Sum_1keV.root", dMCDir.c_str()));

  fBulkOuter = new TFile(Form("%s/OldProd/MCProduction_BulkOuter_1keV.root", dMCDir.c_str()));
  fBulkOuterOld = new TFile(Form("%s/OldProd/MCProduction_BulkOuter_1keV.root", dMCDir.c_str()));
  fBulkOuterM2Sum = new TFile(Form("%s/OldProd/MCProduction_BulkOuterM2Sum_1keV.root", dMCDir.c_str()));

  fSurfaceCrystal = new TFile(Form("%s/OldProd/MCProduction_SurfaceCrystal_1keV_new.root", dMCDir.c_str()));
  fSurfaceCrystalOld = new TFile(Form("%s/OldProd/MCProduction_SurfaceCrystal_1keV_new.root", dMCDir.c_str()));
  fSurfaceOther = new TFile(Form("%s/OldProd/MCProduction_SurfaceOther_1keV.root", dMCDir.c_str()));
  fSurfaceOtherOld = new TFile(Form("%s/OldProd/MCProduction_SurfaceOther_1keV.root", dMCDir.c_str()));

  fFudge = new TFile(Form("%s/MCProduction_FudgeFactor_1keV.root", dMCDir.c_str()));

///////////// Bulk Histograms
/////// Crystal M1 and M2
  hTeO20nuM1     = (TH1D*)fBulkInnerOld->Get("hTeO20nuM1");
  hTeO22nuM1     = (TH1D*)fBulkInnerOld->Get("hTeO22nuM1");
  hTeO2co60M1    = (TH1D*)fBulkInnerOld->Get("hTeO2co60M1");
  hTeO2k40M1     = (TH1D*)fBulkInnerOld->Get("hTeO2k40M1");
  hTeO2pb210M1   = (TH1D*)fBulkInnerOld->Get("hTeO2pb210M1");
  hTeO2po210M1   = (TH1D*)fBulkInnerOld->Get("hTeO2po210M1");
  hTeO2te125M1   = (TH1D*)fBulkInnerOld->Get("hTeO2te125M1");
  hTeO2th232M1   = (TH1D*)fBulkInnerOld->Get("hTeO2th232M1");
  // hTeO2th228M1   = (TH1D*)fBulkInner->Get("hTeO2th228M1");
  // hTeO2ra226M1   = (TH1D*)fBulkInner->Get("hTeO2ra226M1");
  // hTeO2rn222M1   = (TH1D*)fBulkInner->Get("hTeO2rn222M1");
  hTeO2u238M1    = (TH1D*)fBulkInnerOld->Get("hTeO2u238M1");
  // hTeO2th230M1   = (TH1D*)fBulkInner->Get("hTeO2th230M1");
  // hTeO2u234M1    = (TH1D*)fBulkInner->Get("hTeO2u234M1");

  hTeO2th232onlyM1 = (TH1D*)fBulkInnerOld->Get("hTeO2th232onlyM1");
  hTeO2ra228pb208M1 = (TH1D*)fBulkInnerOld->Get("hTeO2ra228pb208M1");
  hTeO2th230onlyM1 = (TH1D*)fBulkInnerOld->Get("hTeO2th230onlyM1");

  hTeO20nuM2     = (TH1D*)fBulkInnerOld->Get("hTeO20nuM2");
  hTeO22nuM2     = (TH1D*)fBulkInnerOld->Get("hTeO22nuM2");
  hTeO2co60M2    = (TH1D*)fBulkInnerOld->Get("hTeO2co60M2");
  hTeO2k40M2     = (TH1D*)fBulkInnerOld->Get("hTeO2k40M2");
  hTeO2pb210M2   = (TH1D*)fBulkInnerOld->Get("hTeO2pb210M2");
  hTeO2po210M2   = (TH1D*)fBulkInnerOld->Get("hTeO2po210M2");
  hTeO2te125M2   = (TH1D*)fBulkInnerOld->Get("hTeO2te125M2");
  hTeO2th232M2   = (TH1D*)fBulkInnerOld->Get("hTeO2th232M2");
  // hTeO2th228M2   = (TH1D*)fBulkInner->Get("hTeO2th228M2");
  // hTeO2ra226M2   = (TH1D*)fBulkInner->Get("hTeO2ra226M2");
  // hTeO2rn222M2   = (TH1D*)fBulkInner->Get("hTeO2rn222M2");
  hTeO2u238M2    = (TH1D*)fBulkInnerOld->Get("hTeO2u238M2");
  // hTeO2th230M2   = (TH1D*)fBulkInner->Get("hTeO2th230M2");
  // hTeO2u234M2    = (TH1D*)fBulkInner->Get("hTeO2u234M2");

  hTeO2th232onlyM2 = (TH1D*)fBulkInnerOld->Get("hTeO2th232onlyM2");
  hTeO2ra228pb208M2 = (TH1D*)fBulkInnerOld->Get("hTeO2ra228pb208M2");
  hTeO2th230onlyM2 = (TH1D*)fBulkInnerOld->Get("hTeO2th230onlyM2");

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
  hCuFrameco58M1    = (TH1D*)fBulkInnerOld->Get("hCuFrameco58M1");
  hCuFrameco60M1    = (TH1D*)fBulkInnerOld->Get("hCuFrameco60M1");
  hCuFramecs137M1   = (TH1D*)fBulkInnerOld->Get("hCuFramecs137M1");
  hCuFramek40M1     = (TH1D*)fBulkInnerOld->Get("hCuFramek40M1");
  hCuFramemn54M1    = (TH1D*)fBulkInnerOld->Get("hCuFramemn54M1");
  hCuFramepb210M1   = (TH1D*)fBulkInnerOld->Get("hCuFramepb210M1");
  hCuFrameth232M1   = (TH1D*)fBulkInnerOld->Get("hCuFrameth232M1");
  hCuFrameu238M1    = (TH1D*)fBulkInnerOld->Get("hCuFrameu238M1");
   
  hCuFrameco58M2    = (TH1D*)fBulkInnerOld->Get("hCuFrameco58M2");
  hCuFrameco60M2    = (TH1D*)fBulkInnerOld->Get("hCuFrameco60M2");
  hCuFramecs137M2   = (TH1D*)fBulkInnerOld->Get("hCuFramecs137M2");
  hCuFramek40M2     = (TH1D*)fBulkInnerOld->Get("hCuFramek40M2");
  hCuFramemn54M2    = (TH1D*)fBulkInnerOld->Get("hCuFramemn54M2");
  hCuFramepb210M2   = (TH1D*)fBulkInnerOld->Get("hCuFramepb210M2");
  hCuFrameth232M2   = (TH1D*)fBulkInnerOld->Get("hCuFrameth232M2");
  hCuFrameu238M2    = (TH1D*)fBulkInnerOld->Get("hCuFrameu238M2");

  hCuFrameco58M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuFrameco58M2Sum");
  hCuFrameco60M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuFrameco60M2Sum");
  hCuFramecs137M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hCuFramecs137M2Sum");
  hCuFramek40M2Sum     = (TH1D*)fBulkInnerM2Sum->Get("hCuFramek40M2Sum");
  hCuFramemn54M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuFramemn54M2Sum");
  hCuFramepb210M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hCuFramepb210M2Sum");
  hCuFrameth232M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hCuFrameth232M2Sum");
  hCuFrameu238M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuFrameu238M2Sum");

///////// CuBox (TShield) M1 and M2
  hCuBoxco58M1    = (TH1D*)fBulkInnerOld->Get("hCuBoxco58M1");
  hCuBoxco60M1    = (TH1D*)fBulkInnerOld->Get("hCuBoxco60M1");
  hCuBoxcs137M1   = (TH1D*)fBulkInnerOld->Get("hCuBoxcs137M1");
  hCuBoxk40M1     = (TH1D*)fBulkInnerOld->Get("hCuBoxk40M1");
  hCuBoxmn54M1    = (TH1D*)fBulkInnerOld->Get("hCuBoxmn54M1");
  hCuBoxpb210M1   = (TH1D*)fBulkInnerOld->Get("hCuBoxpb210M1");
  hCuBoxth232M1   = (TH1D*)fBulkInnerOld->Get("hCuBoxth232M1");
  hCuBoxu238M1    = (TH1D*)fBulkInnerOld->Get("hCuBoxu238M1");
   
  hCuBoxco58M2    = (TH1D*)fBulkInnerOld->Get("hCuBoxco58M2");
  hCuBoxco60M2    = (TH1D*)fBulkInnerOld->Get("hCuBoxco60M2");
  hCuBoxcs137M2   = (TH1D*)fBulkInnerOld->Get("hCuBoxcs137M2");
  hCuBoxk40M2     = (TH1D*)fBulkInnerOld->Get("hCuBoxk40M2");
  hCuBoxmn54M2    = (TH1D*)fBulkInnerOld->Get("hCuBoxmn54M2");
  hCuBoxpb210M2   = (TH1D*)fBulkInnerOld->Get("hCuBoxpb210M2");
  hCuBoxth232M2   = (TH1D*)fBulkInnerOld->Get("hCuBoxth232M2");
  hCuBoxu238M2    = (TH1D*)fBulkInnerOld->Get("hCuBoxu238M2");

  hCuBoxco58M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxco58M2Sum");
  hCuBoxco60M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxco60M2Sum");
  hCuBoxcs137M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxcs137M2Sum");
  hCuBoxk40M2Sum     = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxk40M2Sum");
  hCuBoxmn54M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxmn54M2Sum");
  hCuBoxpb210M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxpb210M2Sum");
  hCuBoxth232M2Sum   = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxth232M2Sum");
  hCuBoxu238M2Sum    = (TH1D*)fBulkInnerM2Sum->Get("hCuBoxu238M2Sum");

///////// CuBox + CuFrame M1 and M2
  hCuBox_CuFrameco60M1 = (TH1D*)fBulkInnerOld->Get("hCuBox_CuFrameco60M1");
  hCuBox_CuFramek40M1 = (TH1D*)fBulkInnerOld->Get("hCuBox_CuFramek40M1");
  hCuBox_CuFrameth232M1 = (TH1D*)fBulkInnerOld->Get("hCuBox_CuFrameth232M1");
  hCuBox_CuFrameu238M1 = (TH1D*)fBulkInnerOld->Get("hCuBox_CuFrameu238M1");

  hCuBox_CuFrameco60M2 = (TH1D*)fBulkInnerOld->Get("hCuBox_CuFrameco60M2");
  hCuBox_CuFramek40M2 = (TH1D*)fBulkInnerOld->Get("hCuBox_CuFramek40M2");
  hCuBox_CuFrameth232M2 = (TH1D*)fBulkInnerOld->Get("hCuBox_CuFrameth232M2");
  hCuBox_CuFrameu238M2 = (TH1D*)fBulkInnerOld->Get("hCuBox_CuFrameu238M2");

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
  hPbRomco60M1    = (TH1D*)fBulkOuterOld->Get("hPbRomco60M1");
  hPbRomcs137M1   = (TH1D*)fBulkOuterOld->Get("hPbRomcs137M1");
  hPbRomk40M1     = (TH1D*)fBulkOuterOld->Get("hPbRomk40M1");
  hPbRompb210M1   = (TH1D*)fBulkOuterOld->Get("hPbRompb210M1");
  hPbRomth232M1   = (TH1D*)fBulkOuterOld->Get("hPbRomth232M1");
  hPbRomu238M1    = (TH1D*)fBulkOuterOld->Get("hPbRomu238M1");
   
  // hPbRombi207M2   = (TH1D*)fBulkOuter->Get("hPbRombi207M2");
  hPbRomco60M2    = (TH1D*)fBulkOuterOld->Get("hPbRomco60M2");
  hPbRomcs137M2   = (TH1D*)fBulkOuterOld->Get("hPbRomcs137M2");
  hPbRomk40M2     = (TH1D*)fBulkOuterOld->Get("hPbRomk40M2");
  hPbRompb210M2   = (TH1D*)fBulkOuterOld->Get("hPbRompb210M2");
  hPbRomth232M2   = (TH1D*)fBulkOuterOld->Get("hPbRomth232M2");
  hPbRomu238M2    = (TH1D*)fBulkOuterOld->Get("hPbRomu238M2");

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
  // hInternalco60M1 = (TH1D*)fBulkInnerOld->Get("hInternalco60M1");
  // hInternalk40M1 = (TH1D*)fBulkOuter->Get("hInternalk40M1");
  // hInternalth232M1 = (TH1D*)fBulkOuter->Get("hInternalth232M1");
  // hInternalu238M1 = (TH1D*)fBulkInnerOld->Get("hInternalu238M1");

  // hInternalco60M2 = (TH1D*)fBulkInnerOld->Get("hInternalco60M2");
  // hInternalk40M2 = (TH1D*)fBulkOuter->Get("hInternalk40M2");
  // hInternalth232M2 = (TH1D*)fBulkOuter->Get("hInternalth232M2");
  // hInternalu238M2 = (TH1D*)fBulkInnerOld->Get("hInternalu238M2");

  // Old Production
  hInternalco60M1 = (TH1D*)fBulkInnerOld->Get("hInternalco60M1");
  hInternalk40M1 = (TH1D*)fBulkInnerOld->Get("hInternalk40M1");
  hInternalth232M1 = (TH1D*)fBulkInnerOld->Get("hInternalth232M1");
  hInternalu238M1 = (TH1D*)fBulkInnerOld->Get("hInternalu238M1");

  hInternalco60M2 = (TH1D*)fBulkInnerOld->Get("hInternalco60M2");
  hInternalk40M2 = (TH1D*)fBulkInnerOld->Get("hInternalk40M2");
  hInternalth232M2 = (TH1D*)fBulkInnerOld->Get("hInternalth232M2");
  hInternalu238M2 = (TH1D*)fBulkInnerOld->Get("hInternalu238M2");

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
  hOVCco60M1    = (TH1D*)fBulkOuterOld->Get("hOVCco60M1");
  hOVCk40M1     = (TH1D*)fBulkOuterOld->Get("hOVCk40M1");
  hOVCth232M1   = (TH1D*)fBulkOuterOld->Get("hOVCth232M1");
  hOVCu238M1    = (TH1D*)fBulkOuterOld->Get("hOVCu238M1");
   
  hOVCco60M2    = (TH1D*)fBulkOuterOld->Get("hOVCco60M2");
  hOVCk40M2     = (TH1D*)fBulkOuterOld->Get("hOVCk40M2");
  hOVCth232M2   = (TH1D*)fBulkOuterOld->Get("hOVCth232M2");
  hOVCu238M2    = (TH1D*)fBulkOuterOld->Get("hOVCu238M2");

  hOVCco60M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hOVCco60M2Sum");
  hOVCk40M2Sum     = (TH1D*)fBulkOuterM2Sum->Get("hOVCk40M2Sum");
  hOVCth232M2Sum   = (TH1D*)fBulkOuterM2Sum->Get("hOVCth232M2Sum");
  hOVCu238M2Sum    = (TH1D*)fBulkOuterM2Sum->Get("hOVCu238M2Sum");

/////// External Sources M1 and M2
  hExtPbbi210M1 = (TH1D*)fBulkOuterOld->Get("hExtPbbi210M1");
 
  hExtPbbi210M2 = (TH1D*)fBulkOuterOld->Get("hExtPbbi210M2");

  // hExtPbbi210M2Sum = (TH1D*)fBulkOuterM2Sum->Get("hExtPbbi210M2Sum");

////////// Fudge Factors
  hOVC804M1 = (TH1D*)fFudge->Get("hOVC804M1");
  hOVC1063M1 = (TH1D*)fFudge->Get("hOVC1063M1");
  hPbRom804M1 = (TH1D*)fFudge->Get("hPbRom804M1");
  hPbRom1063M1 = (TH1D*)fFudge->Get("hPbRom1063M1");

  hOVC804M2 = (TH1D*)fFudge->Get("hOVC804M2");
  hOVC1063M2 = (TH1D*)fFudge->Get("hOVC1063M2");
  hPbRom804M2 = (TH1D*)fFudge->Get("hPbRom804M2");
  hPbRom1063M2 = (TH1D*)fFudge->Get("hPbRom1063M2");

  hOVC804M2Sum = (TH1D*)fFudge->Get("hOVC804M2Sum");
  hOVC1063M2Sum = (TH1D*)fFudge->Get("hOVC1063M2Sum");
  hPbRom804M2Sum = (TH1D*)fFudge->Get("hPbRom804M2Sum");
  hPbRom1063M2Sum = (TH1D*)fFudge->Get("hPbRom1063M2Sum");

//////////// Surface PDFs
///// Crystal M1 and M2
  // hTeO2Spb210M1_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Spb210M1_01");
  // hTeO2Spo210M1_001   = (TH1D*)fSurfaceCrystal->Get("hTeO2Spo210M1_001");
  // hTeO2Spo210M1_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Spo210M1_01");
  // hTeO2Sth232M1_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sth232M1_01");
  // hTeO2Su238M1_01     = (TH1D*)fSurfaceCrystal->Get("hTeO2Su238M1_01");
  hTeO2Sxpb210M1_001  = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpb210M1_001");
  hTeO2Sxpb210M1_01   = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpb210M1_01");
  hTeO2Sxpb210M1_1    = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpb210M1_1");
  hTeO2Sxpb210M1_10   = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpb210M1_10");
  hTeO2Sxpo210M1_001  = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpo210M1_001");
  hTeO2Sxpo210M1_01   = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpo210M1_01");
  hTeO2Sxpo210M1_1    = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpo210M1_1");
  hTeO2Sxth232M1_001  = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxth232M1_001");
  hTeO2Sxth232M1_01   = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxth232M1_01");
  hTeO2Sxth232M1_1    = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxth232M1_1");
  hTeO2Sxth232M1_10   = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxth232M1_10");
  hTeO2Sxu238M1_001   = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxu238M1_001");
  hTeO2Sxu238M1_01    = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxu238M1_01");
  hTeO2Sxu238M1_1     = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxu238M1_1");
  hTeO2Sxu238M1_10    = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxu238M1_10");

  hTeO2Sxu238M1_100    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M1_100");
  hTeO2Sxth232M1_100   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M1_100");
  hTeO2Sxpb210M1_100   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M1_100");

  hTeO2Sxth232onlyM1_001 = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxth232onlyM1_001");
  hTeO2Sxra228pb208M1_001 = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxra228pb208M1_001");
  hTeO2Sxu238th230M1_001 = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxu238th230M1_001");
  hTeO2Sxth230onlyM1_001 = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxth230onlyM1_001");
  hTeO2Sxra226pb210M1_001 = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxra226pb210M1_001");
  hTeO2Sxpb210M1_0001 = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpb210M1_0001");

  hTeO2Sxth232onlyM1_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232onlyM1_01");
  hTeO2Sxra228pb208M1_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra228pb208M1_01");
  hTeO2Sxu238th230M1_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238th230M1_01");
  hTeO2Sxth230onlyM1_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth230onlyM1_01");
  hTeO2Sxra226pb210M1_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra226pb210M1_01");

  hTeO2Sxth232onlyM1_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232onlyM1_0001");
  hTeO2Sxra228pb208M1_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra228pb208M1_0001");
  hTeO2Sxu238th230M1_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238th230M1_0001");
  hTeO2Sxth230onlyM1_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth230onlyM1_0001");
  hTeO2Sxra226pb210M1_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra226pb210M1_0001");

  // hTeO2Spb210M2_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Spb210M2_01");
  // hTeO2Spo210M2_001   = (TH1D*)fSurfaceCrystal->Get("hTeO2Spo210M2_001");
  // hTeO2Spo210M2_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Spo210M2_01");
  // hTeO2Sth232M2_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sth232M2_01");
  // hTeO2Su238M2_01     = (TH1D*)fSurfaceCrystal->Get("hTeO2Su238M2_01");
  hTeO2Sxpb210M2_001  = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpb210M2_001");
  hTeO2Sxpb210M2_01   = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpb210M2_01");
  hTeO2Sxpb210M2_1    = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpb210M2_1");
  hTeO2Sxpb210M2_10   = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpb210M2_10");
  hTeO2Sxpo210M2_001  = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpo210M2_001");
  hTeO2Sxpo210M2_01   = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpo210M2_01");
  hTeO2Sxpo210M2_1    = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpo210M2_1");
  hTeO2Sxth232M2_001  = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxth232M2_001");
  hTeO2Sxth232M2_01   = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxth232M2_01");
  hTeO2Sxth232M2_1    = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxth232M2_1");
  hTeO2Sxth232M2_10   = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxth232M2_10");
  hTeO2Sxu238M2_001   = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxu238M2_001");
  hTeO2Sxu238M2_01    = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxu238M2_01");
  hTeO2Sxu238M2_1     = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxu238M2_1");
  hTeO2Sxu238M2_10    = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxu238M2_10");

  hTeO2Sxu238M2_100    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238M2_100");
  hTeO2Sxth232M2_100   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232M2_100");
  hTeO2Sxpb210M2_100   = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxpb210M2_100");

  hTeO2Sxth232onlyM2_001 = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxth232onlyM2_001");
  hTeO2Sxra228pb208M2_001 = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxra228pb208M2_001");
  hTeO2Sxu238th230M2_001 = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxu238th230M2_001");
  hTeO2Sxth230onlyM2_001 = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxth230onlyM2_001");
  hTeO2Sxra226pb210M2_001 = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxra226pb210M2_001");
  hTeO2Sxpb210M2_0001 = (TH1D*)fSurfaceCrystalOld->Get("hTeO2Sxpb210M2_0001");

  hTeO2Sxth232onlyM2_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232onlyM2_01");
  hTeO2Sxra228pb208M2_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra228pb208M2_01");
  hTeO2Sxu238th230M2_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238th230M2_01");
  hTeO2Sxth230onlyM2_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth230onlyM2_01");
  hTeO2Sxra226pb210M2_01 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra226pb210M2_01");

  hTeO2Sxth232onlyM2_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth232onlyM2_0001");
  hTeO2Sxra228pb208M2_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra228pb208M2_0001");
  hTeO2Sxu238th230M2_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxu238th230M2_0001");
  hTeO2Sxth230onlyM2_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxth230onlyM2_0001");
  hTeO2Sxra226pb210M2_0001 = (TH1D*)fSurfaceCrystal->Get("hTeO2Sxra226pb210M2_0001");

  // hTeO2Spb210M2Sum_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Spb210M2Sum_01");
  // hTeO2Spo210M2Sum_001   = (TH1D*)fSurfaceCrystal->Get("hTeO2Spo210M2Sum_001");
  // hTeO2Spo210M2Sum_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Spo210M2Sum_01");
  // hTeO2Sth232M2Sum_01    = (TH1D*)fSurfaceCrystal->Get("hTeO2Sth232M2Sum_01");
  // hTeO2Su238M2Sum_01     = (TH1D*)fSurfaceCrystal->Get("hTeO2Su238M2Sum_01");
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

//////// Frame M1 and M2
  // hCuFrameSth232M1_1    = (TH1D*)fSurfaceOther->Get("hCuFrameSth232M1_1");
  // hCuFrameSu238M1_1     = (TH1D*)fSurfaceOther->Get("hCuFrameSu238M1_1");
  hCuFrameSxpb210M1_001 = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxpb210M1_001");
  hCuFrameSxpb210M1_01  = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxpb210M1_01");
  hCuFrameSxpb210M1_1   = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxpb210M1_1");
  hCuFrameSxpb210M1_10  = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxpb210M1_10");
  hCuFrameSxth232M1_001 = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M1_001");
  hCuFrameSxth232M1_01  = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M1_01");
  hCuFrameSxth232M1_1   = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M1_1");
  hCuFrameSxth232M1_10  = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M1_10");
  hCuFrameSxu238M1_001  = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxu238M1_001");
  hCuFrameSxu238M1_01   = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxu238M1_01");
  hCuFrameSxu238M1_1    = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxu238M1_1");
  hCuFrameSxu238M1_10   = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxu238M1_10");

  // hCuFrameSth232M2_1    = (TH1D*)fSurfaceOther->Get("hCuFrameSth232M2_1");
  // hCuFrameSu238M2_1     = (TH1D*)fSurfaceOther->Get("hCuFrameSu238M2_1");
  hCuFrameSxpb210M2_001 = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxpb210M2_001");
  hCuFrameSxpb210M2_01  = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxpb210M2_01");
  hCuFrameSxpb210M2_1   = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxpb210M2_1");
  hCuFrameSxpb210M2_10  = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxpb210M2_10");
  hCuFrameSxth232M2_001 = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M2_001");
  hCuFrameSxth232M2_01  = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M2_01");
  hCuFrameSxth232M2_1   = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M2_1");
  hCuFrameSxth232M2_10  = (TH1D*)fSurfaceOther->Get("hCuFrameSxth232M2_10");
  hCuFrameSxu238M2_001  = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxu238M2_001");
  hCuFrameSxu238M2_01   = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxu238M2_01");
  hCuFrameSxu238M2_1    = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxu238M2_1");
  hCuFrameSxu238M2_10   = (TH1D*)fSurfaceOtherOld->Get("hCuFrameSxu238M2_10");

  // hCuFrameSth232M2Sum_1    = (TH1D*)fSurfaceOther->Get("hCuFrameSth232M2Sum_1");
  // hCuFrameSu238M2Sum_1     = (TH1D*)fSurfaceOther->Get("hCuFrameSu238M2Sum_1");
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
  // hCuBoxSth232M1_1    = (TH1D*)fSurfaceOther->Get("hCuBoxSth232M1_1");
  // hCuBoxSu238M1_1     = (TH1D*)fSurfaceOther->Get("hCuBoxSu238M1_1");
  hCuBoxSxpb210M1_001 = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxpb210M1_001");
  hCuBoxSxpb210M1_01  = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxpb210M1_01");
  hCuBoxSxpb210M1_1   = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxpb210M1_1");
  hCuBoxSxpb210M1_10  = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxpb210M1_10");
  hCuBoxSxth232M1_001 = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxth232M1_001");
  hCuBoxSxth232M1_01  = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxth232M1_01");
  hCuBoxSxth232M1_1   = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxth232M1_1");
  hCuBoxSxth232M1_10  = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxth232M1_10");
  hCuBoxSxu238M1_001  = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxu238M1_001");
  hCuBoxSxu238M1_01   = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxu238M1_01");
  hCuBoxSxu238M1_1    = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxu238M1_1");
  hCuBoxSxu238M1_10   = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxu238M1_10");

  // hCuBoxSth232M2_1    = (TH1D*)fSurfaceOther->Get("hCuBoxSth232M2_1");
  // hCuBoxSu238M2_1     = (TH1D*)fSurfaceOther->Get("hCuBoxSu238M2_1");
  hCuBoxSxpb210M2_001 = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxpb210M2_001");
  hCuBoxSxpb210M2_01  = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxpb210M2_01");
  hCuBoxSxpb210M2_1   = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxpb210M2_1");
  hCuBoxSxpb210M2_10  = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxpb210M2_10");
  hCuBoxSxth232M2_001 = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxth232M2_001");
  hCuBoxSxth232M2_01  = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxth232M2_01");
  hCuBoxSxth232M2_1   = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxth232M2_1");
  hCuBoxSxth232M2_10  = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxth232M2_10");
  hCuBoxSxu238M2_001  = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxu238M2_001");
  hCuBoxSxu238M2_01   = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxu238M2_01");
  hCuBoxSxu238M2_1    = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxu238M2_1");
  hCuBoxSxu238M2_10   = (TH1D*)fSurfaceOtherOld->Get("hCuBoxSxu238M2_10");

  // hCuBoxSth232M2Sum_1    = (TH1D*)fSurfaceOther->Get("hCuBoxSth232M2Sum_1");
  // hCuBoxSu238M2Sum_1     = (TH1D*)fSurfaceOther->Get("hCuBoxSu238M2Sum_1");
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
  hCuBox_CuFrameu238M1_10 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M1_10");
  hCuBox_CuFramepb210M1_10 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M1_10");
  hCuBox_CuFramepb210M1_1 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M1_1");
  hCuBox_CuFramepb210M1_01 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M1_01");
  hCuBox_CuFramepb210M1_001 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M1_001");

  hCuBox_CuFrameth232M1_1 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M1_1");
  hCuBox_CuFrameu238M1_1 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M1_1");
  hCuBox_CuFrameth232M1_01 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M1_01");
  hCuBox_CuFrameu238M1_01 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M1_01");
  hCuBox_CuFrameth232M1_001 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M1_001");
  hCuBox_CuFrameu238M1_001 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M1_001");  


  hCuBox_CuFrameth232M1_100 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M1_100");
  hCuBox_CuFrameu238M1_100 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M1_100");
  hCuBox_CuFramepb210M1_100 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M1_100");
  hCuBox_CuFrameth232M1_50 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M1_50");
  hCuBox_CuFrameu238M1_50 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M1_50");
  hCuBox_CuFramepb210M1_50 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M1_50");
  hCuBox_CuFrameth232M1_5 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M1_5");
  hCuBox_CuFrameu238M1_5 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M1_5");
  hCuBox_CuFramepb210M1_5 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M1_5");

  hCuBox_CuFrameth232M2_10 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2_10");
  hCuBox_CuFrameu238M2_10 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M2_10");
  hCuBox_CuFramepb210M2_10 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M2_10");
  hCuBox_CuFramepb210M2_1 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M2_1");
  hCuBox_CuFramepb210M2_01 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M2_01");
  hCuBox_CuFramepb210M2_001 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M2_001");

  hCuBox_CuFrameth232M2_1 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2_1");
  hCuBox_CuFrameu238M2_1 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M2_1");
  hCuBox_CuFrameth232M2_01 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2_01");
  hCuBox_CuFrameu238M2_01 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M2_01");
  hCuBox_CuFrameth232M2_001 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2_001");
  hCuBox_CuFrameu238M2_001 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M2_001"); 

  hCuBox_CuFrameth232M2_100 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2_100");
  hCuBox_CuFrameu238M2_100 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M2_100");
  hCuBox_CuFramepb210M2_100 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M2_100");
  hCuBox_CuFrameth232M2_50 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2_50");
  hCuBox_CuFrameu238M2_50 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M2_50");
  hCuBox_CuFramepb210M2_50 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M2_50");
  hCuBox_CuFrameth232M2_5 = (TH1D*)fSurfaceOther->Get("hCuBox_CuFrameth232M2_5");
  hCuBox_CuFrameu238M2_5 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFrameu238M2_5");
  hCuBox_CuFramepb210M2_5 = (TH1D*)fSurfaceOtherOld->Get("hCuBox_CuFramepb210M2_5");

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

  // hnewTeO2Spb210M1_01 = hTeO2Spb210M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Spb210M1_01", dAdaptiveArrayM1);
  // hnewTeO2Spo210M1_001 = hTeO2Spo210M1_001->Rebin(dAdaptiveBinsM1, "hnewTeO2Spo210M1_001", dAdaptiveArrayM1);
  // hnewTeO2Spo210M1_01 = hTeO2Spo210M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Spo210M1_01", dAdaptiveArrayM1);
  // hnewTeO2Sth232M1_01 = hTeO2Sth232M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Sth232M1_01", dAdaptiveArrayM1);
  // hnewTeO2Su238M1_01 = hTeO2Su238M1_01->Rebin(dAdaptiveBinsM1, "hnewTeO2Su238M1_01", dAdaptiveArrayM1);
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

  // hnewTeO2Spb210M2_01 = hTeO2Spb210M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Spb210M2_01", dAdaptiveArrayM2);
  // hnewTeO2Spo210M2_001 = hTeO2Spo210M2_001->Rebin(dAdaptiveBinsM2, "hnewTeO2Spo210M2_001", dAdaptiveArrayM2);
  // hnewTeO2Spo210M2_01 = hTeO2Spo210M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Spo210M2_01", dAdaptiveArrayM2);
  // hnewTeO2Sth232M2_01 = hTeO2Sth232M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Sth232M2_01", dAdaptiveArrayM2);
  // hnewTeO2Su238M2_01 = hTeO2Su238M2_01->Rebin(dAdaptiveBinsM2, "hnewTeO2Su238M2_01", dAdaptiveArrayM2);
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

  // hnewTeO2Spb210M2Sum_01 = hTeO2Spb210M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Spb210M2Sum_01", dAdaptiveArrayM2Sum);
  // hnewTeO2Spo210M2Sum_001 = hTeO2Spo210M2Sum_001->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Spo210M2Sum_001", dAdaptiveArrayM2Sum);
  // hnewTeO2Spo210M2Sum_01 = hTeO2Spo210M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Spo210M2Sum_01", dAdaptiveArrayM2Sum);
  // hnewTeO2Sth232M2Sum_01 = hTeO2Sth232M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Sth232M2Sum_01", dAdaptiveArrayM2Sum);
  // hnewTeO2Su238M2Sum_01 = hTeO2Su238M2Sum_01->Rebin(dAdaptiveBinsM2Sum, "hnewTeO2Su238M2Sum_01", dAdaptiveArrayM2Sum);
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

////////// Frame M1 and M2
  hnewCuFrameco58M1 = hCuFrameco58M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameco58M1", dAdaptiveArrayM1);
  hnewCuFrameco60M1 = hCuFrameco60M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameco60M1", dAdaptiveArrayM1);
  hnewCuFramecs137M1 = hCuFramecs137M1->Rebin(dAdaptiveBinsM1, "hnewCuFramecs137M1", dAdaptiveArrayM1);
  hnewCuFramek40M1 = hCuFramek40M1->Rebin(dAdaptiveBinsM1, "hnewCuFramek40M1", dAdaptiveArrayM1);
  hnewCuFramemn54M1 = hCuFramemn54M1->Rebin(dAdaptiveBinsM1, "hnewCuFramemn54M1", dAdaptiveArrayM1);
  hnewCuFramepb210M1 = hCuFramepb210M1->Rebin(dAdaptiveBinsM1, "hnewCuFramepb210M1", dAdaptiveArrayM1);
  hnewCuFrameth232M1 = hCuFrameth232M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameth232M1", dAdaptiveArrayM1);
  hnewCuFrameu238M1 = hCuFrameu238M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameu238M1", dAdaptiveArrayM1);

  // hnewCuFrameSth232M1_1  = hCuFrameSth232M1_1->Rebin(dAdaptiveBinsM1, "hnewCuFrameSth232M1_1", dAdaptiveArrayM1);
  // hnewCuFrameSu238M1_1 = hCuFrameSu238M1_1->Rebin(dAdaptiveBinsM1, "hnewCuFrameSu238M1_1", dAdaptiveArrayM1);
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

  // hnewCuFrameSth232M2_1 = hCuFrameSth232M2_1->Rebin(dAdaptiveBinsM2, "hnewCuFrameSth232M2_1", dAdaptiveArrayM2);
  // hnewCuFrameSu238M2_1 = hCuFrameSu238M2_1->Rebin(dAdaptiveBinsM2, "hnewCuFrameSu238M2_1", dAdaptiveArrayM2);
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

  // hnewCuBoxSth232M1_1 = hCuBoxSth232M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBoxSth232M1_1", dAdaptiveArrayM1);
  // hnewCuBoxSu238M1_1 = hCuBoxSu238M1_1->Rebin(dAdaptiveBinsM1, "hnewCuBoxSu238M1_1", dAdaptiveArrayM1);
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

  // hnewCuBoxSth232M2_1 = hCuBoxSth232M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBoxSth232M2_1", dAdaptiveArrayM2);
  // hnewCuBoxSu238M2_1 = hCuBoxSu238M2_1->Rebin(dAdaptiveBinsM2, "hnewCuBoxSu238M2_1", dAdaptiveArrayM2);
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

  ///////// Fudge Factors
  hnewOVC804M1 = hOVC804M1->Rebin(dAdaptiveBinsM1, "hnewOVC804M1", dAdaptiveArrayM1);
  hnewOVC1063M1 = hOVC1063M1->Rebin(dAdaptiveBinsM1, "hnewOVC1063M1", dAdaptiveArrayM1);
  hnewPbRom804M1 = hPbRom804M1->Rebin(dAdaptiveBinsM1, "hnewPbRom804M1", dAdaptiveArrayM1);
  hnewPbRom1063M1 = hPbRom1063M1->Rebin(dAdaptiveBinsM1, "hnewPbRom1063M1", dAdaptiveArrayM1);

  hnewOVC804M2 = hOVC804M2->Rebin(dAdaptiveBinsM2, "hnewOVC804M2", dAdaptiveArrayM2);
  hnewOVC1063M2 = hOVC1063M2->Rebin(dAdaptiveBinsM2, "hnewOVC1063M2", dAdaptiveArrayM2);
  hnewPbRom804M2 = hPbRom804M2->Rebin(dAdaptiveBinsM2, "hnewPbRom804M2", dAdaptiveArrayM2);
  hnewPbRom1063M2 = hPbRom1063M2->Rebin(dAdaptiveBinsM2, "hnewPbRom1063M2", dAdaptiveArrayM2);

  hnewOVC804M2Sum = hOVC804M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewOVC804M2Sum", dAdaptiveArrayM2Sum);
  hnewOVC1063M2Sum = hOVC1063M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewOVC1063M2Sum", dAdaptiveArrayM2Sum);
  hnewPbRom804M2Sum = hPbRom804M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewPbRom804M2Sum", dAdaptiveArrayM2Sum);
  hnewPbRom1063M2Sum = hPbRom1063M2Sum->Rebin(dAdaptiveBinsM2Sum, "hnewPbRom1063M2Sum", dAdaptiveArrayM2Sum);

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

    // hAdapTeO2Spb210M1_01->SetBinContent(i, hnewTeO2Spb210M1_01->GetBinContent(i)/hnewTeO2Spb210M1_01->GetBinWidth(i));
    // hAdapTeO2Spo210M1_001->SetBinContent(i, hnewTeO2Spo210M1_001->GetBinContent(i)/hnewTeO2Spo210M1_001->GetBinWidth(i));
    // hAdapTeO2Spo210M1_01->SetBinContent(i, hnewTeO2Spo210M1_01->GetBinContent(i)/hnewTeO2Spo210M1_01->GetBinWidth(i));
    // hAdapTeO2Sth232M1_01->SetBinContent(i, hnewTeO2Sth232M1_01->GetBinContent(i)/hnewTeO2Sth232M1_01->GetBinWidth(i));
    // hAdapTeO2Su238M1_01->SetBinContent(i, hnewTeO2Su238M1_01->GetBinContent(i)/hnewTeO2Su238M1_01->GetBinWidth(i));
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


    hAdapCuFrameco58M1->SetBinContent(i, hnewCuFrameco58M1->GetBinContent(i)/hnewCuFrameco58M1->GetBinWidth(i));
    hAdapCuFrameco60M1->SetBinContent(i, hnewCuFrameco60M1->GetBinContent(i)/hnewCuFrameco60M1->GetBinWidth(i));
    hAdapCuFramecs137M1->SetBinContent(i, hnewCuFramecs137M1->GetBinContent(i)/hnewCuFramecs137M1->GetBinWidth(i));
    hAdapCuFramek40M1->SetBinContent(i, hnewCuFramek40M1->GetBinContent(i)/hnewCuFramek40M1->GetBinWidth(i));
    hAdapCuFramemn54M1->SetBinContent(i, hnewCuFramemn54M1->GetBinContent(i)/hnewCuFramemn54M1->GetBinWidth(i));
    hAdapCuFramepb210M1->SetBinContent(i, hnewCuFramepb210M1->GetBinContent(i)/hnewCuFramepb210M1->GetBinWidth(i));
    hAdapCuFrameth232M1->SetBinContent(i, hnewCuFrameth232M1->GetBinContent(i)/hnewCuFrameth232M1->GetBinWidth(i));
    hAdapCuFrameu238M1->SetBinContent(i, hnewCuFrameu238M1->GetBinContent(i)/hnewCuFrameu238M1->GetBinWidth(i));

    // hAdapCuFrameSth232M1_1->SetBinContent(i, hnewCuFrameSth232M1_1->GetBinContent(i)/hnewCuFrameSth232M1_1->GetBinWidth(i));
    // hAdapCuFrameSu238M1_1->SetBinContent(i, hnewCuFrameSu238M1_1->GetBinContent(i)/hnewCuFrameSu238M1_1->GetBinWidth(i));
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

    // hAdapCuBoxSth232M1_1->SetBinContent(i, hnewCuBoxSth232M1_1->GetBinContent(i)/hnewCuBoxSth232M1_1->GetBinWidth(i));
    // hAdapCuBoxSu238M1_1->SetBinContent(i, hnewCuBoxSu238M1_1->GetBinContent(i)/hnewCuBoxSu238M1_1->GetBinWidth(i));
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

    hAdapOVC804M1->SetBinContent(i, hnewOVC804M1->GetBinContent(i)/hnewOVC804M1->GetBinWidth(i));
    hAdapOVC1063M1->SetBinContent(i, hnewOVC1063M1->GetBinContent(i)/hnewOVC1063M1->GetBinWidth(i));
    hAdapPbRom804M1->SetBinContent(i, hnewPbRom804M1->GetBinContent(i)/hnewPbRom804M1->GetBinWidth(i));
    hAdapPbRom1063M1->SetBinContent(i, hnewPbRom1063M1->GetBinContent(i)/hnewPbRom1063M1->GetBinWidth(i));

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

    // hAdapTeO2Spb210M2_01->SetBinContent(i, hnewTeO2Spb210M2_01->GetBinContent(i)/hnewTeO2Spb210M2_01->GetBinWidth(i));
    // hAdapTeO2Spo210M2_001->SetBinContent(i, hnewTeO2Spo210M2_001->GetBinContent(i)/hnewTeO2Spo210M2_001->GetBinWidth(i));
    // hAdapTeO2Spo210M2_01->SetBinContent(i, hnewTeO2Spo210M2_01->GetBinContent(i)/hnewTeO2Spo210M2_01->GetBinWidth(i));
    // hAdapTeO2Sth232M2_01->SetBinContent(i, hnewTeO2Sth232M2_01->GetBinContent(i)/hnewTeO2Sth232M2_01->GetBinWidth(i));
    // hAdapTeO2Su238M2_01->SetBinContent(i, hnewTeO2Su238M2_01->GetBinContent(i)/hnewTeO2Su238M2_01->GetBinWidth(i));
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

    hAdapTeO2Sxu238M2_100->SetBinContent(i, hnewTeO2Sxu238M2_100->GetBinContent(i)/hnewTeO2Sxu238M2_100->GetBinWidth(i));
    hAdapTeO2Sxth232M2_100->SetBinContent(i, hnewTeO2Sxth232M2_100->GetBinContent(i)/hnewTeO2Sxth232M2_100->GetBinWidth(i));
    hAdapTeO2Sxpb210M2_100->SetBinContent(i, hnewTeO2Sxpb210M2_100->GetBinContent(i)/hnewTeO2Sxpb210M2_100->GetBinWidth(i));

    hAdapTeO2th232onlyM2->SetBinContent(i, hnewTeO2th232onlyM2->GetBinContent(i)/hnewTeO2th232onlyM2->GetBinWidth(i));
    hAdapTeO2ra228pb208M2->SetBinContent(i, hnewTeO2ra228pb208M2->GetBinContent(i)/hnewTeO2ra228pb208M2->GetBinWidth(i));
    hAdapTeO2th230onlyM2->SetBinContent(i, hnewTeO2th230onlyM2->GetBinContent(i)/hnewTeO2th230onlyM2->GetBinWidth(i));

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

    hAdapCuFrameco58M2->SetBinContent(i, hnewCuFrameco58M2->GetBinContent(i)/hnewCuFrameco58M2->GetBinWidth(i));
    hAdapCuFrameco60M2->SetBinContent(i, hnewCuFrameco60M2->GetBinContent(i)/hnewCuFrameco60M2->GetBinWidth(i));
    hAdapCuFramecs137M2->SetBinContent(i, hnewCuFramecs137M2->GetBinContent(i)/hnewCuFramecs137M2->GetBinWidth(i));
    hAdapCuFramek40M2->SetBinContent(i, hnewCuFramek40M2->GetBinContent(i)/hnewCuFramek40M2->GetBinWidth(i));
    hAdapCuFramemn54M2->SetBinContent(i, hnewCuFramemn54M2->GetBinContent(i)/hnewCuFramemn54M2->GetBinWidth(i));
    hAdapCuFramepb210M2->SetBinContent(i, hnewCuFramepb210M2->GetBinContent(i)/hnewCuFramepb210M2->GetBinWidth(i));
    hAdapCuFrameth232M2->SetBinContent(i, hnewCuFrameth232M2->GetBinContent(i)/hnewCuFrameth232M2->GetBinWidth(i));
    hAdapCuFrameu238M2->SetBinContent(i, hnewCuFrameu238M2->GetBinContent(i)/hnewCuFrameu238M2->GetBinWidth(i));

    // hAdapCuFrameSth232M2_1->SetBinContent(i, hnewCuFrameSth232M2_1->GetBinContent(i)/hnewCuFrameSth232M2_1->GetBinWidth(i));
    // hAdapCuFrameSu238M2_1->SetBinContent(i, hnewCuFrameSu238M2_1->GetBinContent(i)/hnewCuFrameSu238M2_1->GetBinWidth(i));
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

    // hAdapCuBoxSth232M2_1->SetBinContent(i, hnewCuBoxSth232M2_1->GetBinContent(i)/hnewCuBoxSth232M2_1->GetBinWidth(i));
    // hAdapCuBoxSu238M2_1->SetBinContent(i, hnewCuBoxSu238M2_1->GetBinContent(i)/hnewCuBoxSu238M2_1->GetBinWidth(i));
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

    hAdapOVC804M2->SetBinContent(i, hnewOVC804M2->GetBinContent(i)/hnewOVC804M2->GetBinWidth(i));
    hAdapOVC1063M2->SetBinContent(i, hnewOVC1063M2->GetBinContent(i)/hnewOVC1063M2->GetBinWidth(i));
    hAdapPbRom804M2->SetBinContent(i, hnewPbRom804M2->GetBinContent(i)/hnewPbRom804M2->GetBinWidth(i));
    hAdapPbRom1063M2->SetBinContent(i, hnewPbRom1063M2->GetBinContent(i)/hnewPbRom1063M2->GetBinWidth(i));    
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

    // hAdapTeO2Spb210M2Sum_01->SetBinContent(i, hnewTeO2Spb210M2Sum_01->GetBinContent(i)/hnewTeO2Spb210M2Sum_01->GetBinWidth(i));
    // hAdapTeO2Spo210M2Sum_001->SetBinContent(i, hnewTeO2Spo210M2Sum_001->GetBinContent(i)/hnewTeO2Spo210M2Sum_001->GetBinWidth(i));
    // hAdapTeO2Spo210M2Sum_01->SetBinContent(i, hnewTeO2Spo210M2Sum_01->GetBinContent(i)/hnewTeO2Spo210M2Sum_01->GetBinWidth(i));
    // hAdapTeO2Sth232M2Sum_01->SetBinContent(i, hnewTeO2Sth232M2Sum_01->GetBinContent(i)/hnewTeO2Sth232M2Sum_01->GetBinWidth(i));
    // hAdapTeO2Su238M2Sum_01->SetBinContent(i, hnewTeO2Su238M2Sum_01->GetBinContent(i)/hnewTeO2Su238M2Sum_01->GetBinWidth(i));
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

    hAdapCuFrameco58M2Sum->SetBinContent(i, hnewCuFrameco58M2Sum->GetBinContent(i)/hnewCuFrameco58M2Sum->GetBinWidth(i));
    hAdapCuFrameco60M2Sum->SetBinContent(i, hnewCuFrameco60M2Sum->GetBinContent(i)/hnewCuFrameco60M2Sum->GetBinWidth(i));
    hAdapCuFramecs137M2Sum->SetBinContent(i, hnewCuFramecs137M2Sum->GetBinContent(i)/hnewCuFramecs137M2Sum->GetBinWidth(i));
    hAdapCuFramek40M2Sum->SetBinContent(i, hnewCuFramek40M2Sum->GetBinContent(i)/hnewCuFramek40M2Sum->GetBinWidth(i));
    hAdapCuFramemn54M2Sum->SetBinContent(i, hnewCuFramemn54M2Sum->GetBinContent(i)/hnewCuFramemn54M2Sum->GetBinWidth(i));
    hAdapCuFramepb210M2Sum->SetBinContent(i, hnewCuFramepb210M2Sum->GetBinContent(i)/hnewCuFramepb210M2Sum->GetBinWidth(i));
    hAdapCuFrameth232M2Sum->SetBinContent(i, hnewCuFrameth232M2Sum->GetBinContent(i)/hnewCuFrameth232M2Sum->GetBinWidth(i));
    hAdapCuFrameu238M2Sum->SetBinContent(i, hnewCuFrameu238M2Sum->GetBinContent(i)/hnewCuFrameu238M2Sum->GetBinWidth(i));

    // hAdapCuFrameSth232M2Sum_1->SetBinContent(i, hnewCuFrameSth232M2Sum_1->GetBinContent(i)/hnewCuFrameSth232M2Sum_1->GetBinWidth(i));
    // hAdapCuFrameSu238M2Sum_1->SetBinContent(i, hnewCuFrameSu238M2Sum_1->GetBinContent(i)/hnewCuFrameSu238M2Sum_1->GetBinWidth(i));
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

    // hAdapCuBoxSth232M2Sum_1->SetBinContent(i, hnewCuBoxSth232M2Sum_1->GetBinContent(i)/hnewCuBoxSth232M2Sum_1->GetBinWidth(i));
    // hAdapCuBoxSu238M2Sum_1->SetBinContent(i, hnewCuBoxSu238M2Sum_1->GetBinContent(i)/hnewCuBoxSu238M2Sum_1->GetBinWidth(i));
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

/*
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
*/

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

/*
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
*/

    hAdapOVCco60M2Sum->SetBinContent(i, hnewOVCco60M2Sum->GetBinContent(i)/hnewOVCco60M2Sum->GetBinWidth(i));
    hAdapOVCk40M2Sum->SetBinContent(i, hnewOVCk40M2Sum->GetBinContent(i)/hnewOVCk40M2Sum->GetBinWidth(i));
    hAdapOVCth232M2Sum->SetBinContent(i, hnewOVCth232M2Sum->GetBinContent(i)/hnewOVCth232M2Sum->GetBinWidth(i));
    hAdapOVCu238M2Sum->SetBinContent(i, hnewOVCu238M2Sum->GetBinContent(i)/hnewOVCu238M2Sum->GetBinWidth(i));

    hAdapExtPbbi210M2Sum->SetBinContent(i, hnewExtPbbi210M2Sum->GetBinContent(i)/hnewExtPbbi210M2Sum->GetBinWidth(i));

    hAdapOVC804M2Sum->SetBinContent(i, hnewOVC804M2Sum->GetBinContent(i)/hnewOVC804M2Sum->GetBinWidth(i));
    hAdapOVC1063M2Sum->SetBinContent(i, hnewOVC1063M2Sum->GetBinContent(i)/hnewOVC1063M2Sum->GetBinWidth(i));
    hAdapPbRom804M2Sum->SetBinContent(i, hnewPbRom804M2Sum->GetBinContent(i)/hnewPbRom804M2Sum->GetBinWidth(i));
    hAdapPbRom1063M2Sum->SetBinContent(i, hnewPbRom1063M2Sum->GetBinContent(i)/hnewPbRom1063M2Sum->GetBinWidth(i));   

  }

}
