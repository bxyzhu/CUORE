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

// In an effort to make things more legible...


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

  // Implement a method in your class that calculates the quantity you want to minimize, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimise fval
  if(bAdaptiveBinning)
  {
    Obj->UpdateModelAdaptive();
    fval = Obj->GetChiSquareAdaptive();
  }
  else
  {
  	Obj->UpdateModel();
    fval = Obj->GetChiSquare();
  }

}


TBackgroundModel::TBackgroundModel(double fFitMin, double fFitMax)
{

  dNumCalls = 0;
  dSecToYears = 1./(60*60*24*365);

  dDataDir =  "/Users/brian/macros/Simulations/Bkg/Unsmeared/";
  dDataIntegral = 0;
  bAdaptiveBinning = false; // Start off false, can be turned on with DoTheFitAdaptive()

  // Bin size (keV) -- base binning is 2 keV
  dBinSize = 2; 
  // Histogram range - from 0 to 10 MeV
  dMinEnergy = 0.;
  dMaxEnergy = 10000.;

  if(fFitMin >= fFitMax)
  {
    cout << "Fit Min >= Fit Max!" << endl;
  }

  // Fitting range
  dFitMin = fFitMin;
  dFitMax = fFitMax;

  dNBins = (dMaxEnergy - dMinEnergy)/ dBinSize;

  // Data Histograms
  fDataHistoTot  = new TH1D("fDataHistoTot",  "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM1   = new TH1D("fDataHistoM1",   "", dNBins, dMinEnergy, dMaxEnergy);
  fDataHistoM2   = new TH1D("fDataHistoM2",   "", dNBins, dMinEnergy, dMaxEnergy);

  qtree = new TChain("qredtree");
  // Data cuts 
  // qtree = new TChain("qtree");
  // base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
  // base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
  // base_cut = base_cut && "abs(BaselineSlope)<0.1";
  // base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";

  // Load data here
  LoadData();

  cout << "Loaded Data" << endl;

  // Scaling by livetime, should change this part when using Toy data
  dDataIntegral = fDataHistoM1->Integral(1, dNBins);
  int dDataIntegralTot = qtree->GetEntries();

  cout << "Total Events in background spectrum: " << dDataIntegralTot << endl; 
  cout << "Events in background spectrum (M1): " << fDataHistoM1->Integral(1, 3000/dBinSize) << endl;
  cout << "Events in background spectrum (M2): " << fDataHistoM2->Integral(1, 3000/dBinSize) << endl;

  // Scale by Live-time (ds 2061 - 2100) 14647393.0 seconds
  fDataHistoM1->Scale(1/((936398+14647393.0) * dSecToYears));
  fDataHistoM2->Scale(1/((936398+14647393.0) * dSecToYears));  

  cout << "Normalized Data using Livetime of: " << (936398+14647393.0) * dSecToYears << " years" <<endl;




  // Calculates adaptive binning vectors
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


  // Frame M1 and M2
  hCuFrameco58M1      = new TH1D("hCuFrameco58M1",   "hCuFrameco58M1",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameco60M1      = new TH1D("hCuFrameco60M1",   "hCuFrameco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFramecs137M1     = new TH1D("hCuFramecs137M1",  "hCuFramecs137M1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFramek40M1       = new TH1D("hCuFramek40M1",    "hCuFramek40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFramemn54M1      = new TH1D("hCuFramemn54M1",   "hCuFramemn54M1",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFramepb210M1     = new TH1D("hCuFramepb210M1",  "hCuFramepb210M1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameth232M1     = new TH1D("hCuFrameth232M1",  "hCuFrameth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hCuFrameu238M1      = new TH1D("hCuFrameu238M1",   "hCuFrameu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hCuFrameco58M2      = new TH1D("hCuFrameco58M2",   "hCuFrameco58M2",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameco60M2      = new TH1D("hCuFrameco60M2",   "hCuFrameco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFramecs137M2     = new TH1D("hCuFramecs137M2",  "hCuFramecs137M2",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFramek40M2       = new TH1D("hCuFramek40M2",    "hCuFramek40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hCuFramemn54M2      = new TH1D("hCuFramemn54M2",   "hCuFramemn54M2",   dNBins, dMinEnergy, dMaxEnergy);
  hCuFramepb210M2     = new TH1D("hCuFramepb210M2",  "hCuFramepb210M2",  dNBins, dMinEnergy, dMaxEnergy);
  hCuFrameth232M2     = new TH1D("hCuFrameth232M2",  "hCuFrameth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hCuFrameu238M2      = new TH1D("hCuFrameu238M2",   "hCuFrameu238M2",   dNBins, dMinEnergy, dMaxEnergy);

  // CuBox (TShield) M1 and M2
  hCuBoxco58M1      = new TH1D("hCuBoxco58M1",   "hCuBoxco58M1",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxco60M1      = new TH1D("hCuBoxco60M1",   "hCuBoxco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxcs137M1     = new TH1D("hCuBoxcs137M1",  "hCuBoxcs137M1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxk40M1       = new TH1D("hCuBoxk40M1",    "hCuBoxk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxmn54M1      = new TH1D("hCuBoxmn54M1",   "hCuBoxmn54M1",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxpb210M1     = new TH1D("hCuBoxpb210M1",  "hCuBoxpb210M1",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxth232M1     = new TH1D("hCuBoxth232M1",  "hCuBoxth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hCuBoxu238M1      = new TH1D("hCuBoxu238M1",   "hCuBoxu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hCuBoxco58M2      = new TH1D("hCuBoxco58M2",   "hCuBoxco58M2",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxco60M2      = new TH1D("hCuBoxco60M2",   "hCuBoxco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxcs137M2     = new TH1D("hCuBoxcs137M2",  "hCuBoxcs137M2",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxk40M2       = new TH1D("hCuBoxk40M2",    "hCuBoxk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxmn54M2      = new TH1D("hCuBoxmn54M2",   "hCuBoxmn54M2",   dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxpb210M2     = new TH1D("hCuBoxpb210M2",  "hCuBoxpb210M2",  dNBins, dMinEnergy, dMaxEnergy);
  hCuBoxth232M2     = new TH1D("hCuBoxth232M2",  "hCuBoxth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hCuBoxu238M2      = new TH1D("hCuBoxu238M2",   "hCuBoxu238M2",   dNBins, dMinEnergy, dMaxEnergy);

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
  h600mKbi207M1     = new TH1D("h600mKbi207M1",  "h600mKbi207M1",  dNBins, dMinEnergy, dMaxEnergy);  
  h600mKco60M1      = new TH1D("h600mKco60M1",   "h600mKco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  h600mKcs137M1     = new TH1D("h600mKcs137M1",  "h600mKcs137M1",  dNBins, dMinEnergy, dMaxEnergy);  
  h600mKk40M1       = new TH1D("h600mKk40M1",    "h600mKk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  h600mKpb210M1     = new TH1D("h600mKpb210M1",  "h600mKpb210M1",  dNBins, dMinEnergy, dMaxEnergy);  
  h600mKth232M1     = new TH1D("h600mKth232M1",  "h600mKth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  h600mKu238M1      = new TH1D("h600mKu238M1",   "h600mKu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  h600mKbi207M2     = new TH1D("h600mKbi207M2",  "h600mKbi207M2",  dNBins, dMinEnergy, dMaxEnergy);  
  h600mKco60M2      = new TH1D("h600mKco60M2",   "h600mKco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  h600mKcs137M2     = new TH1D("h600mKcs137M2",  "h600mKcs137M2",  dNBins, dMinEnergy, dMaxEnergy);  
  h600mKk40M2       = new TH1D("h600mKk40M2",    "h600mKk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  h600mKpb210M2     = new TH1D("h600mKpb210M2",  "h600mKpb210M2",  dNBins, dMinEnergy, dMaxEnergy);  
  h600mKth232M2     = new TH1D("h600mKth232M2",  "h600mKth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  h600mKu238M2      = new TH1D("h600mKu238M2",   "h600mKu238M2",   dNBins, dMinEnergy, dMaxEnergy);

  // Main bath M1 and M2
  hMBco60M1      = new TH1D("hMBco60M1",   "hMBco60M1",   dNBins, dMinEnergy, dMaxEnergy);
  hMBk40M1       = new TH1D("hMBk40M1",    "hMBk40M1",    dNBins, dMinEnergy, dMaxEnergy);
  hMBth232M1     = new TH1D("hMBth232M1",  "hMBth232M1",  dNBins, dMinEnergy, dMaxEnergy);  
  hMBu238M1      = new TH1D("hMBu238M1",   "hMBu238M1",   dNBins, dMinEnergy, dMaxEnergy);

  hMBco60M2      = new TH1D("hMBco60M2",   "hMBco60M2",   dNBins, dMinEnergy, dMaxEnergy);
  hMBk40M2       = new TH1D("hMBk40M2",    "hMBk40M2",    dNBins, dMinEnergy, dMaxEnergy);
  hMBth232M2     = new TH1D("hMBth232M2",  "hMBth232M2",  dNBins, dMinEnergy, dMaxEnergy);  
  hMBu238M2      = new TH1D("hMBu238M2",   "hMBu238M2",   dNBins, dMinEnergy, dMaxEnergy);  

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

  
//////////////// Adaptive binned histograms
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


  // Frame M1 and M2
  hAdapCuFrameco58M1      = new TH1D("hAdapCuFrameco58M1",   "hAdapCuFrameco58M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameco60M1      = new TH1D("hAdapCuFrameco60M1",   "hAdapCuFrameco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFramecs137M1     = new TH1D("hAdapCuFramecs137M1",  "hAdapCuFramecs137M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFramek40M1       = new TH1D("hAdapCuFramek40M1",    "hAdapCuFramek40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFramemn54M1      = new TH1D("hAdapCuFramemn54M1",   "hAdapCuFramemn54M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFramepb210M1     = new TH1D("hAdapCuFramepb210M1",  "hAdapCuFramepb210M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuFrameth232M1     = new TH1D("hAdapCuFrameth232M1",  "hAdapCuFrameth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapCuFrameu238M1      = new TH1D("hAdapCuFrameu238M1",   "hAdapCuFrameu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuFrameco58M2      = new TH1D("hAdapCuFrameco58M2",   "hAdapCuFrameco58M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameco60M2      = new TH1D("hAdapCuFrameco60M2",   "hAdapCuFrameco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFramecs137M2     = new TH1D("hAdapCuFramecs137M2",  "hAdapCuFramecs137M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFramek40M2       = new TH1D("hAdapCuFramek40M2",    "hAdapCuFramek40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFramemn54M2      = new TH1D("hAdapCuFramemn54M2",   "hAdapCuFramemn54M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFramepb210M2     = new TH1D("hAdapCuFramepb210M2",  "hAdapCuFramepb210M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuFrameth232M2     = new TH1D("hAdapCuFrameth232M2",  "hAdapCuFrameth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapCuFrameu238M2      = new TH1D("hAdapCuFrameu238M2",   "hAdapCuFrameu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  // CuBox (TShield) M1 and M2
  hAdapCuBoxco58M1      = new TH1D("hAdapCuBoxco58M1",   "hAdapCuBoxco58M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxco60M1      = new TH1D("hAdapCuBoxco60M1",   "hAdapCuBoxco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxcs137M1     = new TH1D("hAdapCuBoxcs137M1",  "hAdapCuBoxcs137M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxk40M1       = new TH1D("hAdapCuBoxk40M1",    "hAdapCuBoxk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxmn54M1      = new TH1D("hAdapCuBoxmn54M1",   "hAdapCuBoxmn54M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxpb210M1     = new TH1D("hAdapCuBoxpb210M1",  "hAdapCuBoxpb210M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapCuBoxth232M1     = new TH1D("hAdapCuBoxth232M1",  "hAdapCuBoxth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapCuBoxu238M1      = new TH1D("hAdapCuBoxu238M1",   "hAdapCuBoxu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapCuBoxco58M2      = new TH1D("hAdapCuBoxco58M2",   "hAdapCuBoxco58M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxco60M2      = new TH1D("hAdapCuBoxco60M2",   "hAdapCuBoxco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxcs137M2     = new TH1D("hAdapCuBoxcs137M2",  "hAdapCuBoxcs137M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxk40M2       = new TH1D("hAdapCuBoxk40M2",    "hAdapCuBoxk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxmn54M2      = new TH1D("hAdapCuBoxmn54M2",   "hAdapCuBoxmn54M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxpb210M2     = new TH1D("hAdapCuBoxpb210M2",  "hAdapCuBoxpb210M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapCuBoxth232M2     = new TH1D("hAdapCuBoxth232M2",  "hAdapCuBoxth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapCuBoxu238M2      = new TH1D("hAdapCuBoxu238M2",   "hAdapCuBoxu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

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
  hAdap600mKbi207M1     = new TH1D("hAdap600mKbi207M1",  "hAdap600mKbi207M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdap600mKco60M1      = new TH1D("hAdap600mKco60M1",   "hAdap600mKco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap600mKcs137M1     = new TH1D("hAdap600mKcs137M1",  "hAdap600mKcs137M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdap600mKk40M1       = new TH1D("hAdap600mKk40M1",    "hAdap600mKk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdap600mKpb210M1     = new TH1D("hAdap600mKpb210M1",  "hAdap600mKpb210M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdap600mKth232M1     = new TH1D("hAdap600mKth232M1",  "hAdap600mKth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdap600mKu238M1      = new TH1D("hAdap600mKu238M1",   "hAdap600mKu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdap600mKbi207M2     = new TH1D("hAdap600mKbi207M2",  "hAdap600mKbi207M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdap600mKco60M2      = new TH1D("hAdap600mKco60M2",   "hAdap600mKco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap600mKcs137M2     = new TH1D("hAdap600mKcs137M2",  "hAdap600mKcs137M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdap600mKk40M2       = new TH1D("hAdap600mKk40M2",    "hAdap600mKk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdap600mKpb210M2     = new TH1D("hAdap600mKpb210M2",  "hAdap600mKpb210M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdap600mKth232M2     = new TH1D("hAdap600mKth232M2",  "hAdap600mKth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdap600mKu238M2      = new TH1D("hAdap600mKu238M2",   "hAdap600mKu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);

  // Main bath M1 and M2
  hAdapMBco60M1      = new TH1D("hAdapMBco60M1",   "hAdapMBco60M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapMBk40M1       = new TH1D("hAdapMBk40M1",    "hAdapMBk40M1",    dAdaptiveBinsM1, dAdaptiveArrayM1);
  hAdapMBth232M1     = new TH1D("hAdapMBth232M1",  "hAdapMBth232M1",  dAdaptiveBinsM1, dAdaptiveArrayM1);  
  hAdapMBu238M1      = new TH1D("hAdapMBu238M1",   "hAdapMBu238M1",   dAdaptiveBinsM1, dAdaptiveArrayM1);

  hAdapMBco60M2      = new TH1D("hAdapMBco60M2",   "hAdapMBco60M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapMBk40M2       = new TH1D("hAdapMBk40M2",    "hAdapMBk40M2",    dAdaptiveBinsM2, dAdaptiveArrayM2);
  hAdapMBth232M2     = new TH1D("hAdapMBth232M2",  "hAdapMBth232M2",  dAdaptiveBinsM2, dAdaptiveArrayM2);  
  hAdapMBu238M2      = new TH1D("hAdapMBu238M2",   "hAdapMBu238M2",   dAdaptiveBinsM2, dAdaptiveArrayM2);  

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
  // Correction for accidental coincidence
  // fFileCorrection = new TFile(Form("MCCorrection-%dkeV.root", dBinSize));
  fCorrectionM2 = (TH1D*)fFileCorrection->Get("fModelTotM1");

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

// Returns vector of bin low-edge for adaptive binning
vector<double> TBackgroundModel::AdaptiveBinning(TH1D *h1)
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
	gStyle->SetOptStat(0);
   // This method actually sets up minuit and does the fit
   // TMinuit minuit(14); //initialize minuit, n is the max number of parameters
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
  // fModelTotRaM1->Add(fModelOVCRaM1,     fParameters[3]);

  // fModelTotKM1->Add(fModelFrameKM1,     fParameters[4]);
  fModelTotKM1->Add(fModelTShieldKM1,   fParameters[4]);
  fModelTotKM1->Add(fModel50mKKM1,      fParameters[4]);
  fModelTotKM1->Add(fModel600mKKM1,     fParameters[4]);
  fModelTotKM1->Add(fModelIVCKM1,       fParameters[5]);
  // fModelTotKM1->Add(fModelOVCKM1,       fParameters[5]);

  // fModelTotCoM1->Add(fModelFrameCoM1,   fParameters[6]);
  fModelTotCoM1->Add(fModelTShieldCoM1, fParameters[6]);
  fModelTotCoM1->Add(fModel50mKCoM1,    fParameters[6]);
  fModelTotCoM1->Add(fModel600mKCoM1,   fParameters[6]);
  fModelTotCoM1->Add(fModelIVCCoM1,     fParameters[7]);
  // fModelTotCoM1->Add(fModelOVCCoM1,     fParameters[7]);

  fModelTotMnM1->Add(fModelTShieldMnM1, fParameters[11]);
  fModelTotMnM1->Add(fModelIVCMnM1,     fParameters[12]);

  fModelTotNDBDM1->Add(fModelNDBDM1,    fParameters[9]);
  fModelTot2NDBDM1->Add(fModel2NDBDM1,  fParameters[8]);
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
  // fModelTotRaM2->Add(fModelOVCRaM2,     fParameters[3]);

  // fModelTotKM2->Add(fModelFrameKM2,     fParameters[4]);
  fModelTotKM2->Add(fModelTShieldKM2,   fParameters[4]);
  fModelTotKM2->Add(fModel50mKKM2,      fParameters[4]);
  fModelTotKM2->Add(fModel600mKKM2,     fParameters[4]);
  fModelTotKM2->Add(fModelIVCKM2,       fParameters[5]);
  // fModelTotKM2->Add(fModelOVCKM2,       fParameters[5]);

  // fModelTotCoM2->Add(fModelFrameCoM2,   fParameters[6]);
  fModelTotCoM2->Add(fModelTShieldCoM2, fParameters[6]);
  fModelTotCoM2->Add(fModel50mKCoM2,    fParameters[6]);
  fModelTotCoM2->Add(fModel600mKCoM2,   fParameters[6]);
  fModelTotCoM2->Add(fModelIVCCoM2,     fParameters[7]);
  // fModelTotCoM2->Add(fModelOVCCoM2,     fParameters[7]);

  fModelTotMnM2->Add(fModelTShieldMnM2, fParameters[11]);
  fModelTotMnM2->Add(fModelIVCMnM2,     fParameters[12]);

  fModelTotNDBDM2->Add(fModelNDBDM2,    fParameters[9]);
  fModelTot2NDBDM2->Add(fModel2NDBDM2,  fParameters[8]);
  fModelTotBiM2->Add(fModelBiM2,      fParameters[10]);
*/

  ///////////////////////////////////////////
  //// Many Parameters
  ///////////////////////////////////////////
  /// Use only after previous step converges!
  // 

  // M1 Parameters
  fModelTotThM1->Add(fModelFrameThM1,   fParameters[0]);
  fModelTotThM1->Add(fModelTShieldThM1, fParameters[1]);
  fModelTotThM1->Add(fModel50mKThM1,    fParameters[13]);
  fModelTotThM1->Add(fModel600mKThM1,   fParameters[14]);
  fModelTotThM1->Add(fModelIVCThM1,     fParameters[15]);
  fModelTotThM1->Add(fModelOVCThM1,     fParameters[16]);

  fModelTotRaM1->Add(fModelFrameRaM1,   fParameters[2]);
  fModelTotRaM1->Add(fModelTShieldRaM1, fParameters[3]);
  fModelTotRaM1->Add(fModel50mKRaM1,    fParameters[17]);
  fModelTotRaM1->Add(fModel600mKRaM1,   fParameters[18]);
  fModelTotRaM1->Add(fModelIVCRaM1,     fParameters[19]);
  fModelTotRaM1->Add(fModelOVCRaM1,     fParameters[20]);

  fModelTotKM1->Add(fModelFrameKM1,     fParameters[4]);
  fModelTotKM1->Add(fModelTShieldKM1,   fParameters[4]);
  fModelTotKM1->Add(fModel50mKKM1,      fParameters[4]);
  fModelTotKM1->Add(fModel600mKKM1,     fParameters[4]);
  fModelTotKM1->Add(fModelIVCKM1,       fParameters[5]);
  fModelTotKM1->Add(fModelOVCKM1,       fParameters[5]);

  fModelTotCoM1->Add(fModelFrameCoM1,   fParameters[6]);
  fModelTotCoM1->Add(fModelTShieldCoM1, fParameters[7]);
  fModelTotCoM1->Add(fModel50mKCoM1,    fParameters[21]);
  fModelTotCoM1->Add(fModel600mKCoM1,   fParameters[22]);
  fModelTotCoM1->Add(fModelIVCCoM1,     fParameters[23]);
  fModelTotCoM1->Add(fModelOVCCoM1,     fParameters[24]);

  fModelTotMnM1->Add(fModelTShieldMnM1, fParameters[11]);
  fModelTotMnM1->Add(fModelIVCMnM1,     fParameters[12]);


  fModelTotNDBDM1->Add(fModelNDBDM1,    fParameters[9]);
  fModelTot2NDBDM1->Add(fModel2NDBDM1,  fParameters[8]);
  fModelTotBiM1->Add(fModelBiM1,        fParameters[10]);

  fModelTotPbM1->Add(fModelCrystalBi2M1, fParameters[25]);


  // M2 Parameters
  fModelTotThM2->Add(fModelFrameThM2,   fParameters[0]);
  fModelTotThM2->Add(fModelTShieldThM2, fParameters[1]);
  fModelTotThM2->Add(fModel50mKThM2,    fParameters[13]);
  fModelTotThM2->Add(fModel600mKThM2,   fParameters[14]);
  fModelTotThM2->Add(fModelIVCThM2,     fParameters[15]);
  fModelTotThM2->Add(fModelOVCThM2,     fParameters[16]);

  fModelTotRaM2->Add(fModelFrameRaM2,   fParameters[2]);
  fModelTotRaM2->Add(fModelTShieldRaM2, fParameters[3]);
  fModelTotRaM2->Add(fModel50mKRaM2,    fParameters[17]);
  fModelTotRaM2->Add(fModel600mKRaM2,   fParameters[18]);
  fModelTotRaM2->Add(fModelIVCRaM2,     fParameters[19]);
  fModelTotRaM2->Add(fModelOVCRaM2,     fParameters[20]);

  fModelTotKM2->Add(fModelFrameKM2,     fParameters[4]);
  fModelTotKM2->Add(fModelTShieldKM2,   fParameters[4]);
  fModelTotKM2->Add(fModel50mKKM2,      fParameters[4]);
  fModelTotKM2->Add(fModel600mKKM2,     fParameters[4]);
  fModelTotKM2->Add(fModelIVCKM2,       fParameters[5]);
  fModelTotKM2->Add(fModelOVCKM2,       fParameters[5]);

  fModelTotCoM2->Add(fModelFrameCoM2,   fParameters[6]);
  fModelTotCoM2->Add(fModelTShieldCoM2, fParameters[7]);
  fModelTotCoM2->Add(fModel50mKCoM2,    fParameters[21]);
  fModelTotCoM2->Add(fModel600mKCoM2,   fParameters[22]);
  fModelTotCoM2->Add(fModelIVCCoM2,     fParameters[23]);
  fModelTotCoM2->Add(fModelOVCCoM2,     fParameters[24]);

  fModelTotMnM2->Add(fModelTShieldMnM2, fParameters[11]);
  fModelTotMnM2->Add(fModelIVCMnM2,     fParameters[12]);

  fModelTotNDBDM2->Add(fModelNDBDM2,    fParameters[9]);
  fModelTot2NDBDM2->Add(fModel2NDBDM2,  fParameters[8]);
  fModelTotBiM2->Add(fModelBiM2,      fParameters[10]);

  fModelTotPbM2->Add( fModelCrystalBi2M2, fParameters[25]);

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

// Calculates ChiSquare... model parameters not set here!
double TBackgroundModel::GetChiSquare()
{
	double chiSquare = 0.;
	double datam1_i, errm1_i;
  double datam2_i, errm2_i;
  double modelm1_i, modelm2_i;	

	for(int i = dFitMin/dBinSize+1; i <= dFitMax/dBinSize; i++)
	{

		datam1_i = fDataHistoM1->GetBinContent(i); // For real data
    datam2_i = fDataHistoM2->GetBinContent(i); // For real data

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


void TBackgroundModel::Initialize()
{	

  // Loads PDFs from file
  cout << "Loading PDF Histograms from file" << endl;
  // fFile = new TFile(Form("MCHist-%dkeV.root", dBinSize));
  // fFile = new TFile("MCAdaptiveGamma.root");   
  // fFile = new TFile(Form("MCData-%dkeV.root", dBinSize));
  fFile = new TFile("/Users/brian/macros/Simulations/Production/MCProduction_Bulk_1104.root"); 


  fModelFrameThM1   = (TH1D*)fFile->Get("hCuFrameth232M1");
  fModelTShieldThM1 = (TH1D*)fFile->Get("hCuBoxth232M1");
  fModel50mKThM1    = (TH1D*)fFile->Get("h50mKth232M1");
  fModel600mKThM1   = (TH1D*)fFile->Get("h600mKth232M1");
  fModelIVCThM1     = (TH1D*)fFile->Get("hIVCth232M1");
  fModelOVCThM1     = (TH1D*)fFile->Get("hOVCth232M1");
  fModelMBThM1      = (TH1D*)fFile->Get("hMBth232M1");



  fModelFrameRaM1   = (TH1D*)fFile->Get("hCuFrameu238M1");
  fModelTShieldRaM1 = (TH1D*)fFile->Get("hCuBoxu238M1");
  fModel50mKRaM1    = (TH1D*)fFile->Get("hModel50mKRaM1");
  fModel600mKRaM1   = (TH1D*)fFile->Get("hModel600mKRaM1");
  fModelIVCRaM1     = (TH1D*)fFile->Get("hModelIVCRaM1");
  fModelOVCRaM1     = (TH1D*)fFile->Get("hModelOVCRaM1");

  fModelFrameCoM1   = (TH1D*)fFile->Get("hModelFrameCoM1");
  fModelTShieldCoM1 = (TH1D*)fFile->Get("hModelTShieldCoM1");
  fModel50mKCoM1    = (TH1D*)fFile->Get("hModel50mKCoM1");
  fModel600mKCoM1   = (TH1D*)fFile->Get("hModel600mKCoM1");
  fModelIVCCoM1     = (TH1D*)fFile->Get("hModelIVCCoM1");
  fModelOVCCoM1     = (TH1D*)fFile->Get("hModelOVCCoM1");    

  fModelFrameKM1    = (TH1D*)fFile->Get("hModelFrameKM1");
  fModelTShieldKM1  = (TH1D*)fFile->Get("hModelTShieldKM1");
  fModel50mKKM1     = (TH1D*)fFile->Get("hModel50mKKM1");
  fModel600mKKM1    = (TH1D*)fFile->Get("hModel600mKKM1");
  fModelIVCKM1      = (TH1D*)fFile->Get("hModelIVCKM1");
  fModelOVCKM1      = (TH1D*)fFile->Get("hModelOVCKM1");    

  fModelTShieldMnM1 = (TH1D*)fFile->Get("hModelTShieldMnM1");
  fModelIVCMnM1     = (TH1D*)fFile->Get("hModelIVCMnM1");

  fModelCrystalBi2M1 = (TH1D*)fFile->Get("hModelCrystalBi2M1");
  fModelFrameBi2M1 = (TH1D*)fFile->Get("hModelFrameBi2M1");


  fModel2NDBDM1     = (TH1D*)fFile->Get("hModel2NDBDM1");
  fModelNDBDM1      = (TH1D*)fFile->Get("hModelNDBDM1");
  fModelBiM1        = (TH1D*)fFile->Get("hModelBiM1");

  fModelCrystalPtM1 = (TH1D*)fFile->Get("hModelCrystalPtM1");
  fModelCrystalPbBM1 = (TH1D*)fFile->Get("hModelCrystalPbBM1");
  fModelCrystalPbS01M1 = (TH1D*)fFile->Get("hModelCrystalPbS01M1");
  fModelCrystalPbS1M1 = (TH1D*)fFile->Get("hModelCrystalPbS1M1");
  fModelCrystalPbS10M1 = (TH1D*)fFile->Get("hModelCrystalPbS10M1");
  fModelCrystalPbS100M1 = (TH1D*)fFile->Get("hModelCrystalPbS100M1");
  fModelFramePbBM1 = (TH1D*)fFile->Get("hModelFramePbBM1");
  fModelFramePbS01M1 = (TH1D*)fFile->Get("hModelFramePbS01M1");
  fModelFramePbS1M1 = (TH1D*)fFile->Get("hModelFramePbS1M1");
  fModelFramePbS10M1 = (TH1D*)fFile->Get("hModelFramePbS10M1");
  fModelFramePbS100M1 = (TH1D*)fFile->Get("hModelFramePbS100M1");

  fModelFrameThM2   = (TH1D*)fFile->Get("hModelFrameThM2");
  fModelTShieldThM2 = (TH1D*)fFile->Get("hModelTShieldThM2");
  fModel50mKThM2    = (TH1D*)fFile->Get("hModel50mKThM2");
  fModel600mKThM2   = (TH1D*)fFile->Get("hModel600mKThM2");
  fModelIVCThM2     = (TH1D*)fFile->Get("hModelIVCThM2");
  fModelOVCThM2     = (TH1D*)fFile->Get("hModelOVCThM2");

  fModelFrameRaM2   = (TH1D*)fFile->Get("hModelFrameRaM2");
  fModelTShieldRaM2 = (TH1D*)fFile->Get("hModelTShieldRaM2");
  fModel50mKRaM2    = (TH1D*)fFile->Get("hModel50mKRaM2");
  fModel600mKRaM2   = (TH1D*)fFile->Get("hModel600mKRaM2");
  fModelIVCRaM2     = (TH1D*)fFile->Get("hModelIVCRaM2");
  fModelOVCRaM2     = (TH1D*)fFile->Get("hModelOVCRaM2");

  fModelFrameCoM2   = (TH1D*)fFile->Get("hModelFrameCoM2");
  fModelTShieldCoM2 = (TH1D*)fFile->Get("hModelTShieldCoM2");
  fModel50mKCoM2    = (TH1D*)fFile->Get("hModel50mKCoM2");
  fModel600mKCoM2   = (TH1D*)fFile->Get("hModel600mKCoM2");
  fModelIVCCoM2     = (TH1D*)fFile->Get("hModelIVCCoM2");
  fModelOVCCoM2     = (TH1D*)fFile->Get("hModelOVCCoM2");    

  fModelFrameKM2    = (TH1D*)fFile->Get("hModelFrameKM2");
  fModelTShieldKM2  = (TH1D*)fFile->Get("hModelTShieldKM2");
  fModel50mKKM2     = (TH1D*)fFile->Get("hModel50mKKM2");
  fModel600mKKM2    = (TH1D*)fFile->Get("hModel600mKKM2");
  fModelIVCKM2      = (TH1D*)fFile->Get("hModelIVCKM2");
  fModelOVCKM2      = (TH1D*)fFile->Get("hModelOVCKM2");    

  fModelTShieldMnM2 = (TH1D*)fFile->Get("hModelTShieldMnM2");
  fModelIVCMnM2     = (TH1D*)fFile->Get("hModelIVCMnM2");

  fModelCrystalBi2M2 = (TH1D*)fFile->Get("hModelCrystalBi2M2");
  fModelFrameBi2M2 = (TH1D*)fFile->Get("hModelFrameBi2M2");

  fModel2NDBDM2     = (TH1D*)fFile->Get("hModel2NDBDM2");
  fModelNDBDM2      = (TH1D*)fFile->Get("hModelNDBDM2");
  fModelBiM2        = (TH1D*)fFile->Get("hModelBiM2");    
 
  fModelCrystalPtM2 = (TH1D*)fFile->Get("hModelCrystalPtM2");
  fModelCrystalPbBM2 = (TH1D*)fFile->Get("hModelCrystalPbBM2");
  fModelCrystalPbS01M2 = (TH1D*)fFile->Get("hModelCrystalPbS01M2");
  fModelCrystalPbS1M2 = (TH1D*)fFile->Get("hModelCrystalPbS1M2");
  fModelCrystalPbS10M2 = (TH1D*)fFile->Get("hModelCrystalPbS10M2");
  fModelCrystalPbS100M2 = (TH1D*)fFile->Get("hModelCrystalPbS100M2");
  fModelFramePbBM2 = (TH1D*)fFile->Get("hModelFramePbBM2");
  fModelFramePbS01M2 = (TH1D*)fFile->Get("hModelFramePbS01M2"); 
  fModelFramePbS1M2 = (TH1D*)fFile->Get("hModelFramePbS1M2"); 
  fModelFramePbS10M2 = (TH1D*)fFile->Get("hModelFramePbS10M2"); 
  fModelFramePbS100M2 = (TH1D*)fFile->Get("hModelFramePbS100M2"); 

  // Adaptively binned
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

  // Currently using Jon's reduced file
  qtree->Add("/Users/brian/macros/CUOREZ/Bkg/Q0_DR2_BackgroundSignalData.root"); 
  qtree->Project("fDataHistoTot", "Energy", base_cut);
  qtree->Project("fDataHistoM1",  "Energy", base_cut && "Multiplicity == 1");
  qtree->Project("fDataHistoM2",  "Energy", base_cut && "Multiplicity == 2");

  cout << "Loaded Data" << endl;
}


// Normalize single histogram
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

// Normalizes pair of histograms (for normalizing M2 with the same value as M1)
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
  // fModelTotM1->Add( fModelOVCRaM1,      fParameters[3]);

  // fModelTotM1->Add( fModelFrameKM1,    fParameters[4]);
  fModelTotM1->Add( fModelTShieldKM1,  fParameters[4]);
  fModelTotM1->Add( fModel50mKKM1,     fParameters[4]);
  fModelTotM1->Add( fModel600mKKM1,    fParameters[4]);
  fModelTotM1->Add( fModelIVCKM1,      fParameters[5]);
  // fModelTotM1->Add( fModelOVCKM1,      fParameters[5]); 

  // fModelTotM1->Add( fModelFrameCoM1,    fParameters[6]);
  fModelTotM1->Add( fModelTShieldCoM1,  fParameters[6]);
  fModelTotM1->Add( fModel50mKCoM1,     fParameters[6]);
  fModelTotM1->Add( fModel600mKCoM1,    fParameters[6]);
  fModelTotM1->Add( fModelIVCCoM1,      fParameters[7]);
  // fModelTotM1->Add( fModelOVCCoM1,      fParameters[7]);  

  fModelTotM1->Add( fModelTShieldMnM1,  fParameters[11]);
  fModelTotM1->Add( fModelIVCMnM1,      fParameters[12]);

  fModelTotM1->Add( fModelNDBDM1,      fParameters[9]);  
  fModelTotM1->Add( fModel2NDBDM1,     fParameters[8]);  
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
  // fModelTotM2->Add( fModelOVCRaM2,      fParameters[3]);

  // fModelTotM2->Add( fModelFrameKM2,    fParameters[4]);
  fModelTotM2->Add( fModelTShieldKM2,  fParameters[4]);
  fModelTotM2->Add( fModel50mKKM2,     fParameters[4]);
  fModelTotM2->Add( fModel600mKKM2,    fParameters[4]);
  fModelTotM2->Add( fModelIVCKM2,      fParameters[5]);
  // fModelTotM2->Add( fModelOVCKM2,      fParameters[5]); 

  // fModelTotM2->Add( fModelFrameCoM2,    fParameters[6]);
  fModelTotM2->Add( fModelTShieldCoM2,  fParameters[6]);
  fModelTotM2->Add( fModel50mKCoM2,     fParameters[6]);
  fModelTotM2->Add( fModel600mKCoM2,    fParameters[6]);
  fModelTotM2->Add( fModelIVCCoM2,      fParameters[7]);
  // fModelTotM2->Add( fModelOVCCoM2,      fParameters[7]);  

  fModelTotM2->Add( fModelTShieldMnM2,  fParameters[11]);
  fModelTotM2->Add( fModelIVCMnM2,      fParameters[12]);

  fModelTotM2->Add( fModelNDBDM2,      fParameters[9]);  
  fModelTotM2->Add( fModel2NDBDM2,     fParameters[8]);  
  fModelTotM2->Add( fModelBiM2,        fParameters[10]);  
*/


  /////////////////////////////////////
  //// Many parameters
  ////////////////////////////////////

  // M1
  fModelTotM1->Add( fModelFrameThM1,    fParameters[0]);
  fModelTotM1->Add( fModelTShieldThM1,  fParameters[1]);  
  fModelTotM1->Add( fModel50mKThM1,     fParameters[13]);
  fModelTotM1->Add( fModel600mKThM1,    fParameters[14]);
  fModelTotM1->Add( fModelIVCThM1,      fParameters[15]);
  fModelTotM1->Add( fModelOVCThM1,      fParameters[16]);

  fModelTotM1->Add( fModelFrameRaM1,    fParameters[2]);
  fModelTotM1->Add( fModelTShieldRaM1,  fParameters[3]);  
  fModelTotM1->Add( fModel50mKRaM1,     fParameters[17]);
  fModelTotM1->Add( fModel600mKRaM1,    fParameters[18]);
  fModelTotM1->Add( fModelIVCRaM1,      fParameters[19]);
  fModelTotM1->Add( fModelOVCRaM1,      fParameters[20]);

  fModelTotM1->Add( fModelFrameKM1,     fParameters[4]);
  fModelTotM1->Add( fModelTShieldKM1,   fParameters[4]);
  fModelTotM1->Add( fModel50mKKM1,      fParameters[4]);
  fModelTotM1->Add( fModel600mKKM1,     fParameters[4]);
  fModelTotM1->Add( fModelIVCKM1,       fParameters[5]);
  fModelTotM1->Add( fModelOVCKM1,       fParameters[5]); 

  fModelTotM1->Add( fModelFrameCoM1,    fParameters[6]);
  fModelTotM1->Add( fModelTShieldCoM1,  fParameters[7]);
  fModelTotM1->Add( fModel50mKCoM1,     fParameters[21]);
  fModelTotM1->Add( fModel600mKCoM1,    fParameters[22]);
  fModelTotM1->Add( fModelIVCCoM1,      fParameters[23]);
  fModelTotM1->Add( fModelOVCCoM1,      fParameters[24]);  

  fModelTotM1->Add( fModelTShieldMnM1,  fParameters[11]);
  fModelTotM1->Add( fModelIVCMnM1,      fParameters[12]);

  fModelTotM1->Add( fModelNDBDM1,      fParameters[9]);  
  fModelTotM1->Add( fModel2NDBDM1,     fParameters[8]);  
  fModelTotM1->Add( fModelBiM1,        fParameters[10]);  

  fModelTotM1->Add( fModelCrystalBi2M1, fParameters[25]);  

  // M2
  fModelTotM2->Add( fModelFrameThM2,    fParameters[0]);
  fModelTotM2->Add( fModelTShieldThM2,  fParameters[1]);  
  fModelTotM2->Add( fModel50mKThM2,     fParameters[13]);
  fModelTotM2->Add( fModel600mKThM2,    fParameters[14]);
  fModelTotM2->Add( fModelIVCThM2,      fParameters[15]);
  fModelTotM2->Add( fModelOVCThM2,      fParameters[16]);

  fModelTotM2->Add( fModelFrameRaM2,    fParameters[2]);
  fModelTotM2->Add( fModelTShieldRaM2,  fParameters[3]);  
  fModelTotM2->Add( fModel50mKRaM2,     fParameters[17]);
  fModelTotM2->Add( fModel600mKRaM2,    fParameters[18]);
  fModelTotM2->Add( fModelIVCRaM2,      fParameters[19]);
  fModelTotM2->Add( fModelOVCRaM2,      fParameters[20]);

  fModelTotM2->Add( fModelFrameKM2,     fParameters[4]);
  fModelTotM2->Add( fModelTShieldKM2,   fParameters[4]);
  fModelTotM2->Add( fModel50mKKM2,      fParameters[4]);
  fModelTotM2->Add( fModel600mKKM2,     fParameters[4]);
  fModelTotM2->Add( fModelIVCKM2,       fParameters[5]);
  fModelTotM2->Add( fModelOVCKM2,       fParameters[5]); 

  fModelTotM2->Add( fModelFrameCoM2,    fParameters[6]);
  fModelTotM2->Add( fModelTShieldCoM2,  fParameters[7]);
  fModelTotM2->Add( fModel50mKCoM2,     fParameters[21]);
  fModelTotM2->Add( fModel600mKCoM2,    fParameters[22]);
  fModelTotM2->Add( fModelIVCCoM2,      fParameters[23]);
  fModelTotM2->Add( fModelOVCCoM2,      fParameters[24]);  

  fModelTotM2->Add( fModelTShieldMnM2,  fParameters[11]);
  fModelTotM2->Add( fModelIVCMnM2,      fParameters[12]);

  fModelTotM2->Add( fModelNDBDM2,      fParameters[9]);  
  fModelTotM2->Add( fModel2NDBDM2,     fParameters[8]);  
  fModelTotM2->Add( fModelBiM2,        fParameters[10]);  

  fModelTotM1->Add( fModelCrystalBi2M2 , fParameters[25]);  


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
  bAdaptiveBinning = true;
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

/*
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
*/

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


