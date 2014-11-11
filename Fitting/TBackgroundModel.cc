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

// In an effort to make things more legible I'm trying to clean up the code...

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
  Obj->SetParameters(26, x[26]);
  Obj->SetParameters(27, x[27]);
  Obj->SetParameters(28, x[28]);
  Obj->SetParameters(29, x[29]);
  Obj->SetParameters(30, x[30]);
  Obj->SetParameters(31, x[31]);
  Obj->SetParameters(32, x[32]);
  Obj->SetParameters(33, x[33]);
  Obj->SetParameters(34, x[34]);
  Obj->SetParameters(35, x[35]);
  Obj->SetParameters(36, x[36]);
  Obj->SetParameters(37, x[37]);
  Obj->SetParameters(38, x[38]);
  Obj->SetParameters(39, x[39]);
  Obj->SetParameters(40, x[40]);
  Obj->SetParameters(41, x[41]);
  Obj->SetParameters(42, x[42]);
  Obj->SetParameters(43, x[43]);
  Obj->SetParameters(44, x[44]);
  Obj->SetParameters(45, x[45]);
  Obj->SetParameters(46, x[46]);
  Obj->SetParameters(47, x[47]);
  Obj->SetParameters(48, x[48]);
  Obj->SetParameters(49, x[49]);
  Obj->SetParameters(50, x[50]);
  Obj->SetParameters(51, x[51]);
  Obj->SetParameters(52, x[52]);
  Obj->SetParameters(53, x[53]);
  Obj->SetParameters(54, x[54]);
  Obj->SetParameters(55, x[55]);
  Obj->SetParameters(56, x[56]);
  Obj->SetParameters(57, x[57]);
  Obj->SetParameters(58, x[58]);
  Obj->SetParameters(59, x[59]);
  Obj->SetParameters(60, x[60]); 


  // Implement a method in your class that calculates the quantity you want to minimize, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimise fval
  // if(bAdaptiveBinning)
  // {
    Obj->UpdateModelAdaptive();
    fval = Obj->GetChiSquareAdaptive();
  // }
  // else
  // {
  	// Obj->UpdateModel();
    // fval = Obj->GetChiSquare();
  // }

}


TBackgroundModel::TBackgroundModel(double fFitMin, double fFitMax)
{

  dNumCalls = 0;
  dSecToYears = 1./(60*60*24*365);

  dDataDir =  "/Users/brian/macros/Simulations/Production/";
  // dDataDir =  "/Users/brian/macros/CUOREZ/Bkg/";
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

  // Data variables change when switching between Toy/Real data
  // qtree = new TChain("qredtree");
  qtree = new TChain("CombiTree");
  // Data cuts 
  // qtree = new TChain("qtree");
  // base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
  // base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
  // base_cut = base_cut && "abs(BaselineSlope)<0.1";
  // base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";

  // Load data here
  LoadData();

/*
  // Scaling by livetime, don't use with Toy data
  dDataIntegral = fDataHistoM1->Integral(1, dNBins);
  int dDataIntegralTot = qtree->GetEntries();

  cout << "Total Events in background spectrum: " << dDataIntegralTot << endl; 
  cout << "Events in background spectrum (M1): " << fDataHistoM1->Integral(1, 3000/dBinSize) << endl;
  cout << "Events in background spectrum (M2): " << fDataHistoM2->Integral(1, 3000/dBinSize) << endl;

  // Scale by Live-time (ds 2061 - 2100) 14647393.0 seconds
  fDataHistoM1->Scale(1/((936398+14647393.0) * dSecToYears));
  fDataHistoM2->Scale(1/((936398+14647393.0) * dSecToYears));  

  cout << "Normalized Data using Livetime of: " << (936398+14647393.0) * dSecToYears << " years" << endl;
*/

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


  // Loads all of the PDFs from file
  Initialize();
}
  
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
  delete fModelTotThM1;
  delete fModelTotRaM1;
  delete fModelTotKM1;
  delete fModelTotCoM1;
  delete fModelTotMnM1;
  delete fModelTotNDBDM1;
  delete fModelTot2NDBDM1;
  delete fModelTotBiM1;
  delete fModelTotPtM1;
  delete fModelTotPbM1;

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


  // Residuals
  TGraph      *gResidualM1;
  TGraph      *gResidualM2;

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

    if(dDummyFill >= 10)
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
   TMinuit minuit(61); // for more parameters

   // Reduce Minuit Output
   minuit.SetPrintLevel(1);
//   minuit.Command("SET MINImize 1000 0.001");
   minuit.Command("SET STRategy 2");
  //minuit.Command("SET IMProve 1000 ");
   minuit.SetMaxIterations(1000);
   minuit.SetObjectFit(this); //see the external FCN  above
   
   //define the parameters and set the ranges and initial guesses see ROOTs TMinuit documentation
   ////////////////////////////////////////////////
   // Using more parameters
   ////////////////////////////////////////////////
   minuit.DefineParameter(0, "hTeO20nuM1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(1, "hTeO22nuM1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(2, "hTeO2co60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(3, "hTeO2k40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(4, "hTeO2pb210M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(5, "hTeO2po210M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(6, "hTeO2te125M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(7, "hTeO2th232M1",  10., 10.0, 0., 1000000);
   minuit.DefineParameter(8, "hTeO2th228M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(9, "hTeO2ra226M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(10, "hTeO2rn222M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(11, "hTeO2u238M1",  10., 10.0, 0., 1000000);
   minuit.DefineParameter(12, "hTeO2th230M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(13, "hTeO2u234M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(14, "hCuFrameco58M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(15, "hCuFrameco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(16, "hCuFramecs137M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(17, "hCuFramek40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(18, "hCuFramemn54M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(19, "hCuFramepb210M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(20, "hCuFrameth232M1",  10., 10.0, 0., 1000000);
   minuit.DefineParameter(21, "hCuFrameu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(22, "hCuBoxco58M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(23, "hCuBoxco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(24, "hCuBoxcs137M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(25, "hCuBoxk40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(26, "hCuBoxmn54M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(27, "hCuBoxpb210M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(28, "hCuBoxth232M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(29, "hCuBoxu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(30, "h50mKco58M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(31, "h50mKco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(32, "h50mKcs137M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(33, "h50mKk40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(34, "h50mKmn54M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(35, "h50mKpb210M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(36, "h50mKth232M1",  10., 10.0, 0., 1000000);
   minuit.DefineParameter(37, "h50mKu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(38, "h600mKco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(39, "h600mKk40M1",  0., 10.0, 0., 1000000); 
   minuit.DefineParameter(40, "h600mKth232M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(41, "h600mKu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(42, "hPbRombi207M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(43, "hPbRomco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(44, "hPbRomcs137M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(45, "hPbRomk40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(46, "hPbRompb210M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(47, "hPbRomth232M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(48, "hPbRomu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(49, "hMBco60M1",  0., 10.0, 0., 1000000); 
   minuit.DefineParameter(50, "hMBk40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(51, "hMBth232M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(52, "hMBu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(53, "hIVCco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(54, "hIVCk40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(55, "hIVCth232M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(56, "hIVCu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(57, "hOVCco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(58, "hOVCk40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(59, "hOVCth232M1",  0., 10.0, 0., 1000000);    
   minuit.DefineParameter(60, "hOVCu238M1",  0., 10.0, 0., 1000000);

   // Fix parameters here
   minuit.FixParameter(0); // TeO2 0nu
   minuit.FixParameter(1); // TeO2 2nu
   minuit.FixParameter(2); // TeO2 co60
   minuit.FixParameter(3); // TeO2 k40
   minuit.FixParameter(4); // TeO2 pb210
   minuit.FixParameter(5); // TeO2 po210
   minuit.FixParameter(6); // TeO2 te125m
   // minuit.FixParameter(7); // TeO2 th232
   minuit.FixParameter(8); // TeO2 th232-th228
   minuit.FixParameter(9); // TeO2 u238-ra226
   minuit.FixParameter(10); // TeO2 u238-rn222
   // minuit.FixParameter(11); // TeO2 u238
   minuit.FixParameter(12); // TeO2 u238-th230
   minuit.FixParameter(13); // TeO2 u238-u234
   minuit.FixParameter(14); // Frame co58
   minuit.FixParameter(15); // Frame co60
   minuit.FixParameter(16); // Frame cs137
   minuit.FixParameter(17); // Frame k40
   minuit.FixParameter(18); // Frame mn54
   minuit.FixParameter(19); // Frame pb210
   // minuit.FixParameter(20); // Frame th232
   minuit.FixParameter(21); // Frame u238
   minuit.FixParameter(22); // CuBox co58
   minuit.FixParameter(23); // CuBox co60
   minuit.FixParameter(24); // CuBox cs137
   minuit.FixParameter(25); // CuBox k40
   minuit.FixParameter(26); // CuBox mn54
   minuit.FixParameter(27); // CuBox pb210
   minuit.FixParameter(28); // CuBox th232
   minuit.FixParameter(29); // CuBox u238
   minuit.FixParameter(30); // 50mK co58
   minuit.FixParameter(31); // 50mK co60
   minuit.FixParameter(32); // 50mK cs137
   minuit.FixParameter(33); // 50mK k40
   minuit.FixParameter(34); // 50mK mn54
   minuit.FixParameter(35); // 50mK pb210
   // minuit.FixParameter(36); // 50mK th232
   minuit.FixParameter(37); // 50mK u238
   minuit.FixParameter(38); // 600mK co60
   minuit.FixParameter(39); // 600mK k40
   minuit.FixParameter(40); // 600mK th232
   minuit.FixParameter(41); // 600mK u238
   minuit.FixParameter(42); // RLead bi207
   minuit.FixParameter(43); // RLead co60
   minuit.FixParameter(44); // RLead cs137
   minuit.FixParameter(45); // RLead k40
   minuit.FixParameter(46); // RLead pb210
   minuit.FixParameter(47); // RLead th232
   minuit.FixParameter(48); // RLead u238
   minuit.FixParameter(49); // MB co60
   minuit.FixParameter(50); // MB k40
   minuit.FixParameter(51); // MB th232
   minuit.FixParameter(52); // MB u238
   minuit.FixParameter(53); // IVC co60
   minuit.FixParameter(54); // IVC k40
   minuit.FixParameter(55); // IVC th232
   minuit.FixParameter(56); // IVC u238
   minuit.FixParameter(57); // OVC co60
   minuit.FixParameter(58); // OVC k40
   minuit.FixParameter(59); // OVC th232
   minuit.FixParameter(60); // OVC u238
   // Number of Parameters (for Chi-squared/NDF calculation)
   int dNumParameters = 4;

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
  minuit.GetParameter(26,  fParameters[26],   fParError[26]);
  minuit.GetParameter(27,  fParameters[27],   fParError[27]);
  minuit.GetParameter(28,  fParameters[28],   fParError[28]);
  minuit.GetParameter(29,  fParameters[29],   fParError[29]);
  minuit.GetParameter(30,  fParameters[30],   fParError[30]);
  minuit.GetParameter(31,  fParameters[31],   fParError[31]);
  minuit.GetParameter(32,  fParameters[32],   fParError[32]);
  minuit.GetParameter(33,  fParameters[33],   fParError[33]);
  minuit.GetParameter(34,  fParameters[34],   fParError[34]);
  minuit.GetParameter(35,  fParameters[35],   fParError[35]);
  minuit.GetParameter(36,  fParameters[36],   fParError[36]);
  minuit.GetParameter(37,  fParameters[37],   fParError[37]);
  minuit.GetParameter(38,  fParameters[38],   fParError[38]);
  minuit.GetParameter(39,  fParameters[39],   fParError[39]);
  minuit.GetParameter(40,  fParameters[40],   fParError[40]);
  minuit.GetParameter(41,  fParameters[41],   fParError[41]);
  minuit.GetParameter(42,  fParameters[42],   fParError[42]);
  minuit.GetParameter(43,  fParameters[43],   fParError[43]);
  minuit.GetParameter(44,  fParameters[44],   fParError[44]);
  minuit.GetParameter(45,  fParameters[45],   fParError[45]);
  minuit.GetParameter(46,  fParameters[46],   fParError[46]);
  minuit.GetParameter(47,  fParameters[47],   fParError[47]);
  minuit.GetParameter(48,  fParameters[48],   fParError[48]);
  minuit.GetParameter(49,  fParameters[49],   fParError[49]);
  minuit.GetParameter(50,  fParameters[50],   fParError[50]);
  minuit.GetParameter(51,  fParameters[51],   fParError[51]);
  minuit.GetParameter(52,  fParameters[52],   fParError[52]);
  minuit.GetParameter(53,  fParameters[53],   fParError[53]);
  minuit.GetParameter(54,  fParameters[54],   fParError[54]);
  minuit.GetParameter(55,  fParameters[55],   fParError[55]);
  minuit.GetParameter(56,  fParameters[56],   fParError[56]);
  minuit.GetParameter(57,  fParameters[57],   fParError[57]);
  minuit.GetParameter(58,  fParameters[58],   fParError[58]);
  minuit.GetParameter(59,  fParameters[59],   fParError[59]);
  minuit.GetParameter(60,  fParameters[60],   fParError[60]);
	UpdateModel();
	
	cout << "At the end; ChiSq/NDF = " << GetChiSquare()/(2*(dFitMax-dFitMin)/dBinSize - dNumParameters) << endl; // for M1 and M2
  // cout << "At the end; ChiSq/NDF = " << GetChiSquare()/((dFitMax-dFitMin)/dBinSize - dNumParameters) << endl;  // for M1
  cout << "Total number of calls = " << dNumCalls << endl;


  ///////////////////////////////////////////
  //// Many Parameters
  ///////////////////////////////////////////
  /// Use only after previous step converges!
  // 
  // M1 Parameters
  fModelTotthM1->Add( hTeO2th232M1,     fParameters[7]  );
  fModelTotthM1->Add( hCuFrameth232M1,  fParameters[20] );
  fModelTotthM1->Add( hCuBoxth232M1,    fParameters[28] );
  fModelTotthM1->Add( h50mKth232M1,     fParameters[36] );
  fModelTotthM1->Add( h600mKth232M1,    fParameters[40] );
  fModelTotthM1->Add( hPbRomth232M1,    fParameters[47] );
  fModelTotthM1->Add( hMBth232M1,       fParameters[51] );
  fModelTotthM1->Add( hIVCth232M1,      fParameters[55] );
  fModelTotthM1->Add( hOVCth232M1,      fParameters[59] );

  fModelTotuM1->Add( hTeO2u238M1,       fParameters[11] );
  fModelTotuM1->Add( hCuFrameu238M1,    fParameters[21] );
  fModelTotuM1->Add( hCuBoxu238M1,      fParameters[29] );
  fModelTotuM1->Add( h50mKu238M1,       fParameters[37] );
  fModelTotuM1->Add( h600mKu238M1,      fParameters[41] );
  fModelTotuM1->Add( hPbRomu238M1,      fParameters[48] );
  fModelTotuM1->Add( hMBu238M1,         fParameters[52] );
  fModelTotuM1->Add( hIVCu238M1,        fParameters[56] );
  fModelTotuM1->Add( hOVCu238M1,        fParameters[60] );

  fModelTotkM1->Add( hTeO2k40M1,        fParameters[3]  );
  fModelTotkM1->Add( hCuFramek40M1,     fParameters[17] );
  fModelTotkM1->Add( hCuBoxk40M1,       fParameters[25] );
  fModelTotkM1->Add( h50mKk40M1,        fParameters[33] );
  fModelTotkM1->Add( h600mKk40M1,       fParameters[39] );
  fModelTotkM1->Add( hPbRomk40M1,       fParameters[45] );
  fModelTotkM1->Add( hMBk40M1,          fParameters[50] );
  fModelTotkM1->Add( hIVCk40M1,         fParameters[54] );
  fModelTotkM1->Add( hOVCk40M1,         fParameters[58] ); 

  fModelTotcoM1->Add( hTeO2co60M1,      fParameters[2]  );
  fModelTotcoM1->Add( hCuFrameco60M1,   fParameters[15] );
  fModelTotcoM1->Add( hCuBoxco60M1,     fParameters[23] );
  fModelTotcoM1->Add( h50mKco60M1,      fParameters[31] );
  fModelTotcoM1->Add( h600mKco60M1,     fParameters[38] );
  fModelTotcoM1->Add( hPbRomco60M1,     fParameters[43] );
  fModelTotcoM1->Add( hMBco60M1,        fParameters[49] );
  fModelTotcoM1->Add( hIVCco60M1,       fParameters[53] );
  fModelTotcoM1->Add( hOVCco60M1,       fParameters[57] );

  fModelTotpbM1->Add( hTeO2pb210M1,     fParameters[4]  );
  fModelTotpbM1->Add( hCuFramepb210M1,  fParameters[15] );
  fModelTotpbM1->Add( hCuBoxpb210M1,    fParameters[27] );
  fModelTotpbM1->Add( h50mKpb210M1,     fParameters[35] );
  fModelTotpbM1->Add( hPbRompb210M1,    fParameters[46] );

  fModelTotcsM1->Add( hCuFramecs137M1,  fParameters[16] );
  fModelTotcsM1->Add( hCuBoxcs137M1,    fParameters[24] );
  fModelTotcsM1->Add( h50mKcs137M1,     fParameters[32] );
  fModelTotcsM1->Add( hPbRomcs137M1,    fParameters[44] );

  fModelTotcoM1->Add( hCuFrameco58M1,   fParameters[14] );
  fModelTotcoM1->Add( hCuBoxco58M1,     fParameters[22] );
  fModelTotcoM1->Add( h50mKco58M1,      fParameters[30] );

  fModelTotteo2M1->Add( hTeO2po210M1,   fParameters[5]  );
  fModelTotteo2M1->Add( hTeO2te125M1,   fParameters[6]  );
  fModelTotteo2M1->Add( hTeO2th228M1,   fParameters[8]  );
  fModelTotteo2M1->Add( hTeO2ra226M1,   fParameters[9]  );
  fModelTotteo2M1->Add( hTeO2rn222M1,   fParameters[10] );
  fModelTotteo2M1->Add( hTeO2th230M1,   fParameters[12] );
  fModelTotteo2M1->Add( hTeO2u234M1,    fParameters[13] );

  fModelTotmnM1->Add( hCuFramemn54M1,   fParameters[18] );
  fModelTotmnM1->Add( hCuBoxmn54M1,     fParameters[26] );
  fModelTotmnM1->Add( h50mKmn54M1,      fParameters[34] );

  fModelTotbiM1->Add( hPbRombi207M1,    fParameters[42] );
  fModelTotNDBDM1->Add( hTeO20nuM1,     fParameters[0]  );
  fModelTot2NDBDM1->Add( hTeO22nuM1,    fParameters[1]  );


// M2
  fModelTotthM2->Add( hTeO2th232M2,     fParameters[7]  );
  fModelTotthM2->Add( hCuFrameth232M2,  fParameters[20] );
  fModelTotthM2->Add( hCuBoxth232M2,    fParameters[28] );
  fModelTotthM2->Add( h50mKth232M2,     fParameters[36] );
  fModelTotthM2->Add( h600mKth232M2,    fParameters[40] );
  fModelTotthM2->Add( hPbRomth232M2,    fParameters[47] );
  fModelTotthM2->Add( hMBth232M2,       fParameters[51] );
  fModelTotthM2->Add( hIVCth232M2,      fParameters[55] );
  fModelTotthM2->Add( hOVCth232M2,      fParameters[59] );

  fModelTotuM2->Add( hTeO2u238M2,       fParameters[11] );
  fModelTotuM2->Add( hCuFrameu238M2,    fParameters[21] );
  fModelTotuM2->Add( hCuBoxu238M2,      fParameters[29] );
  fModelTotuM2->Add( h50mKu238M2,       fParameters[37] );
  fModelTotuM2->Add( h600mKu238M2,      fParameters[41] );
  fModelTotuM2->Add( hPbRomu238M2,      fParameters[48] );
  fModelTotuM2->Add( hMBu238M2,         fParameters[52] );
  fModelTotuM2->Add( hIVCu238M2,        fParameters[56] );
  fModelTotuM2->Add( hOVCu238M2,        fParameters[60] );

  fModelTotkM2->Add( hTeO2k40M2,        fParameters[3]  );
  fModelTotkM2->Add( hCuFramek40M2,     fParameters[17] );
  fModelTotkM2->Add( hCuBoxk40M2,       fParameters[25] );
  fModelTotkM2->Add( h50mKk40M2,        fParameters[33] );
  fModelTotkM2->Add( h600mKk40M2,       fParameters[39] );
  fModelTotkM2->Add( hPbRomk40M2,       fParameters[45] );
  fModelTotkM2->Add( hMBk40M2,          fParameters[50] );
  fModelTotkM2->Add( hIVCk40M2,         fParameters[54] );
  fModelTotkM2->Add( hOVCk40M2,         fParameters[58] ); 

  fModelTotcoM2->Add( hTeO2co60M2,      fParameters[2]  );
  fModelTotcoM2->Add( hCuFrameco60M2,   fParameters[15] );
  fModelTotcoM2->Add( hCuBoxco60M2,     fParameters[23] );
  fModelTotcoM2->Add( h50mKco60M2,      fParameters[31] );
  fModelTotcoM2->Add( h600mKco60M2,     fParameters[38] );
  fModelTotcoM2->Add( hPbRomco60M2,     fParameters[43] );
  fModelTotcoM2->Add( hMBco60M2,        fParameters[49] );
  fModelTotcoM2->Add( hIVCco60M2,       fParameters[53] );
  fModelTotcoM2->Add( hOVCco60M2,       fParameters[57] );

  fModelTotpbM2->Add( hTeO2pb210M2,     fParameters[4]  );
  fModelTotpbM2->Add( hCuFramepb210M2,  fParameters[15] );
  fModelTotpbM2->Add( hCuBoxpb210M2,    fParameters[27] );
  fModelTotpbM2->Add( h50mKpb210M2,     fParameters[35] );
  fModelTotpbM2->Add( hPbRompb210M2,    fParameters[46] );

  fModelTotcsM2->Add( hCuFramecs137M2,  fParameters[16] );
  fModelTotcsM2->Add( hCuBoxcs137M2,    fParameters[24] );
  fModelTotcsM2->Add( h50mKcs137M2,     fParameters[32] );
  fModelTotcsM2->Add( hPbRomcs137M2,    fParameters[44] );

  fModelTotcoM2->Add( hCuFrameco58M2,   fParameters[14] );
  fModelTotcoM2->Add( hCuBoxco58M2,     fParameters[22] );
  fModelTotcoM2->Add( h50mKco58M2,      fParameters[30] );

  fModelTotteo2M2->Add( hTeO2po210M2,   fParameters[5]  );
  fModelTotteo2M2->Add( hTeO2te125M2,   fParameters[6]  );
  fModelTotteo2M2->Add( hTeO2th228M2,   fParameters[8]  );
  fModelTotteo2M2->Add( hTeO2ra226M2,   fParameters[9]  );
  fModelTotteo2M2->Add( hTeO2rn222M2,   fParameters[10] );
  fModelTotteo2M2->Add( hTeO2th230M2,   fParameters[12] );
  fModelTotteo2M2->Add( hTeO2u234M2,    fParameters[13] );

  fModelTotmnM2->Add( hCuFramemn54M2,   fParameters[18] );
  fModelTotmnM2->Add( hCuBoxmn54M2,     fParameters[26] );
  fModelTotmnM2->Add( h50mKmn54M2,      fParameters[34] );

  fModelTotbiM2->Add( hPbRombi207M2,    fParameters[42] );
  fModelTotNDBDM2->Add( hTeO20nuM2,     fParameters[0]  );
  fModelTot2NDBDM2->Add( hTeO22nuM2,    fParameters[1]  );

  ////////// Only for testing
  // Correction for M2 spectra, it's the M1 spectra but scaled down by N_M1*1-Exp(2*R*T)
  // 2 because of double counting in M2 spectrum...
  // fTotCorrection->Add(fCorrectionM2, 180197*(1-TMath::Exp(-2*0.05*0.1)));




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

  fModelTotthM1->SetLineColor(3);
  fModelTotthM1->SetLineStyle(2);
  fModelTotuM1->SetLineColor(4);
  fModelTotuM1->SetLineStyle(2);
  fModelTotkM1->SetLineColor(6);
  fModelTotkM1->SetLineStyle(2);
  fModelTotcoM1->SetLineColor(7);
  fModelTotcoM1->SetLineStyle(2);
  fModelTotNDBDM1->SetLineColor(42);
  fModelTotNDBDM1->SetLineStyle(2);
  fModelTot2NDBDM1->SetLineColor(46);
  fModelTot2NDBDM1->SetLineStyle(2);
  fModelTotbiM1->SetLineColor(5);
  fModelTotbiM1->SetLineStyle(2);
  fModelTotmnM1->SetLineColor(40);
  fModelTotmnM1->SetLineStyle(20);


  fModelTotpbM1->SetLineStyle(2);
  fModelTotpbM1->SetLineColor(38);

  fModelTotthM1->Draw("SAME");
  fModelTotuM1->Draw("SAME");
  fModelTotkM1->Draw("SAME");
  fModelTotcoM1->Draw("SAME");
  fModelTotNDBDM1->Draw("SAME");
  fModelTot2NDBDM1->Draw("SAME");
  fModelTotbiM1->Draw("SAME");
  fModelTotmnM1->Draw("SAME");

  fModelTotpbM1->Draw("SAME");


/*
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
*/


  TLegend *legfit1 = new TLegend(0.8,0.8,0.97,0.97);
  legfit1->AddEntry(fModelTotM1, "Total PDF", "l");
  legfit1->AddEntry(fModelTotthM1, "Total th-232", "l");
  legfit1->AddEntry(fModelTotuM1, "Total u-238", "l");
  legfit1->AddEntry(fModelTotkM1, "Total k-40", "l");
  legfit1->AddEntry(fModelTotcoM1, "Total co-60", "l");
  legfit1->AddEntry(fModelTotNDBDM1, "NDBD", "l");
  legfit1->AddEntry(fModelTot2NDBDM1, "2NDBD", "l");
  legfit1->AddEntry(fModelTotbiM1, "bi-207", "l");
  legfit1->AddEntry(fModelTotmnM1, "mn-54", "l");
  legfit1->AddEntry(fModelTotpbM1 , "pb-210", "l");    
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

  fModelTotthM2->SetLineColor(3);
  fModelTotthM2->SetLineStyle(2);
  fModelTotuM2->SetLineColor(4);
  fModelTotuM2->SetLineStyle(2);
  fModelTotkM2->SetLineColor(6);
  fModelTotkM2->SetLineStyle(2);
  fModelTotcoM2->SetLineColor(7);
  fModelTotcoM2->SetLineStyle(2);
  fModelTotNDBDM2->SetLineColor(42);
  fModelTotNDBDM2->SetLineStyle(2);
  fModelTot2NDBDM2->SetLineColor(46);
  fModelTot2NDBDM2->SetLineStyle(2);
  fModelTotbiM2->SetLineColor(5);
  fModelTotbiM2->SetLineStyle(2);
  fModelTotmnM2->SetLineColor(40);
  fModelTotmnM2->SetLineStyle(2);

  fModelTotpbM2->SetLineStyle(2);
  fModelTotpbM2->SetLineColor(38);
  // fTotCorrection->SetLineStyle(2);
  // fTotCorrection->SetLineColor(kBlue+2);

  fModelTotthM2->Draw("SAME");
  fModelTotuM2->Draw("SAME");
  fModelTotkM2->Draw("SAME");
  fModelTotcoM2->Draw("SAME");
  fModelTotNDBDM2->Draw("SAME");
  fModelTot2NDBDM2->Draw("SAME");
  fModelTotbiM2->Draw("SAME");    
  fModelTotmnM2->Draw("SAME");

  fModelTotpbM2->Draw("SAME");
  // fTotCorrection->Draw("SAME");    

/*
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
*/



  TLegend *legfit2 = new TLegend(0.8,0.8,0.97,0.97);
  legfit2->AddEntry(fModelTotM2, "Total PDF", "l");
  legfit2->AddEntry(fModelTotthM2, "Total th-232", "l");
  legfit2->AddEntry(fModelTotuM2, "Total u-238", "l");
  legfit2->AddEntry(fModelTotkM2, "Total k-40", "l");
  legfit2->AddEntry(fModelTotcoM2, "Total co-60", "l");
  legfit2->AddEntry(fModelTotNDBDM2, "NDBD", "l");
  legfit2->AddEntry(fModelTot2NDBDM2, "2NDBD", "l");
  legfit2->AddEntry(fModelTotbiM2, "bi-207", "l");
  legfit2->AddEntry(fModelTotmnM2, "mn-54", "l");
  legfit2->AddEntry(fModelTotpbM2 , "bi-210", "l");
  // legfit2->AddEntry(fTotCorrection, "Accidental coincidence correction", "l");
  legfit2->Draw();



/*
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
  cout << "Integral Total Th PDF in ROI: " << fModelTotthM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTotthM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;
  cout << "Integral Total Ra PDF in ROI: " << fModelTotuM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTotuM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;
  cout << "Integral Total Co PDF in ROI: " << fModelTotcoM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTotcoM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;
  cout << "Integral Total K PDF in ROI: " << fModelTotkM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTotkM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;
  cout << "Integral Total Bi-207 PDF in ROI: " << fModelTotbiM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTotbiM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;  
  cout << "Integral Total 2NDBD PDF in ROI: " << fModelTot2NDBDM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTot2NDBDM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;
  cout << "Integral Total 0NDBD PDF in ROI: " << fModelTotNDBDM1->Integral(2470/dBinSize, 2570/dBinSize) << " +/- " << sqrt(fModelTotNDBDM1->Integral(2470/dBinSize, 2570/dBinSize)) << endl;

  cout << "M2/(M1+M2) = " << (double)fModelTotM2->Integral(300/dBinSize, 3000/dBinSize)/(fModelTotM1->Integral(300/dBinSize, 3000/dBinSize)+fModelTotM2->Integral(300/dBinSize, 3000/dBinSize)) << endl;
*/

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
  fFile = new TFile("/Users/brian/macros/Simulations/Production/MCProduction_Bulk_nonnormalized.root"); 
  // fFile = new TFile("/Users/brian/macros/Simulations/Production/MCProduction_Bulk.root"); 

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


  // Frame M1 and M2
  hCuFrameco58M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameco58M1", dAdaptiveArrayM1);
  hCuFrameco60M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameco60M1", dAdaptiveArrayM1);
  hCuFramecs137M1->Rebin(dAdaptiveBinsM1, "hnewCuFramecs137M1", dAdaptiveArrayM1);
  hCuFramek40M1->Rebin(dAdaptiveBinsM1, "hnewCuFramek40M1", dAdaptiveArrayM1);
  hCuFramemn54M1->Rebin(dAdaptiveBinsM1, "hnewCuFramemn54M1", dAdaptiveArrayM1);
  hCuFramepb210M1->Rebin(dAdaptiveBinsM1, "hnewCuFramepb210M1", dAdaptiveArrayM1);
  hCuFrameth232M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameth232M1", dAdaptiveArrayM1);
  hCuFrameu238M1->Rebin(dAdaptiveBinsM1, "hnewCuFrameu238M1", dAdaptiveArrayM1);

  hCuFrameco58M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameco58M2", dAdaptiveArrayM2);
  hCuFrameco60M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameco60M2", dAdaptiveArrayM2);
  hCuFramecs137M2->Rebin(dAdaptiveBinsM2, "hnewCuFramecs137M2", dAdaptiveArrayM2);
  hCuFramek40M2->Rebin(dAdaptiveBinsM2, "hnewCuFramek40M2", dAdaptiveArrayM2);
  hCuFramemn54M2->Rebin(dAdaptiveBinsM2, "hnewCuFramemn54M2", dAdaptiveArrayM2);
  hCuFramepb210M2->Rebin(dAdaptiveBinsM2, "hnewCuFramepb210M2", dAdaptiveArrayM2);
  hCuFrameth232M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameth232M2", dAdaptiveArrayM2);
  hCuFrameu238M2->Rebin(dAdaptiveBinsM2, "hnewCuFrameu238M2", dAdaptiveArrayM2);

  // CuBox (TShield) M1 and M2
  hCuBoxco58M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxco58M1", dAdaptiveArrayM1);
  hCuBoxco60M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxco60M1", dAdaptiveArrayM1);
  hCuBoxcs137M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxcs137M1", dAdaptiveArrayM1);
  hCuBoxk40M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxk40M1", dAdaptiveArrayM1);
  hCuBoxmn54M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxmn54M1", dAdaptiveArrayM1);
  hCuBoxpb210M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxpb210M1", dAdaptiveArrayM1);
  hCuBoxth232M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxth232M1", dAdaptiveArrayM1);
  hCuBoxu238M1->Rebin(dAdaptiveBinsM1, "hnewCuBoxu238M1", dAdaptiveArrayM1);

  hCuBoxco58M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxco58M2", dAdaptiveArrayM2);
  hCuBoxco60M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxco60M2", dAdaptiveArrayM2);
  hCuBoxcs137M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxcs137M2", dAdaptiveArrayM2);
  hCuBoxk40M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxk40M2", dAdaptiveArrayM2);
  hCuBoxmn54M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxmn54M2", dAdaptiveArrayM2);
  hCuBoxpb210M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxpb210M2", dAdaptiveArrayM2);
  hCuBoxth232M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxth232M2", dAdaptiveArrayM2);
  hCuBoxu238M2->Rebin(dAdaptiveBinsM2, "hnewCuBoxu238M2", dAdaptiveArrayM2);

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

    hAdapCuFrameco58M1->SetBinContent(i, dBinSize * hnewCuFrameco58M1->GetBinContent(i)/hnewCuFrameco58M1->GetBinWidth(i));
    hAdapCuFrameco60M1->SetBinContent(i, dBinSize * hnewCuFrameco60M1->GetBinContent(i)/hnewCuFrameco60M1->GetBinWidth(i));
    hAdapCuFramecs137M1->SetBinContent(i, dBinSize * hnewCuFramecs137M1->GetBinContent(i)/hnewCuFramecs137M1->GetBinWidth(i));
    hAdapCuFramek40M1->SetBinContent(i, dBinSize * hnewCuFramek40M1->GetBinContent(i)/hnewCuFramek40M1->GetBinWidth(i));
    hAdapCuFramemn54M1->SetBinContent(i, dBinSize * hnewCuFramemn54M1->GetBinContent(i)/hnewCuFramemn54M1->GetBinWidth(i));
    hAdapCuFramepb210M1->SetBinContent(i, dBinSize * hnewCuFramepb210M1->GetBinContent(i)/hnewCuFramepb210M1->GetBinWidth(i));
    hAdapCuFrameth232M1->SetBinContent(i, dBinSize * hnewCuFrameth232M1->GetBinContent(i)/hnewCuFrameth232M1->GetBinWidth(i));
    hAdapCuFrameu238M1->SetBinContent(i, dBinSize * hnewCuFrameu238M1->GetBinContent(i)/hnewCuFrameu238M1->GetBinWidth(i));

    hAdapCuBoxco58M1->SetBinContent(i, dBinSize * hnewCuBoxco58M1->GetBinContent(i)/hnewCuBoxco58M1->GetBinWidth(i));
    hAdapCuBoxco60M1->SetBinContent(i, dBinSize * hnewCuBoxco60M1->GetBinContent(i)/hnewCuBoxco60M1->GetBinWidth(i));
    hAdapCuBoxcs137M1->SetBinContent(i, dBinSize * hnewCuBoxcs137M1->GetBinContent(i)/hnewCuBoxcs137M1->GetBinWidth(i));
    hAdapCuBoxk40M1->SetBinContent(i, dBinSize * hnewCuBoxk40M1->GetBinContent(i)/hnewCuBoxk40M1->GetBinWidth(i));
    hAdapCuBoxmn54M1->SetBinContent(i, dBinSize * hnewCuBoxmn54M1->GetBinContent(i)/hnewCuBoxmn54M1->GetBinWidth(i));
    hAdapCuBoxpb210M1->SetBinContent(i, dBinSize * hnewCuBoxpb210M1->GetBinContent(i)/hnewCuBoxpb210M1->GetBinWidth(i));
    hAdapCuBoxth232M1->SetBinContent(i, dBinSize * hnewCuBoxth232M1->GetBinContent(i)/hnewCuBoxth232M1->GetBinWidth(i));
    hAdapCuBoxu238M1->SetBinContent(i, dBinSize * hnewCuBoxu238M1->GetBinContent(i)/hnewCuBoxu238M1->GetBinWidth(i));

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

    hAdapCuFrameco58M2->SetBinContent(i, dBinSize * hnewCuFrameco58M2->GetBinContent(i)/hnewCuFrameco58M2->GetBinWidth(i));
    hAdapCuFrameco60M2->SetBinContent(i, dBinSize * hnewCuFrameco60M2->GetBinContent(i)/hnewCuFrameco60M2->GetBinWidth(i));
    hAdapCuFramecs137M2->SetBinContent(i, dBinSize * hnewCuFramecs137M2->GetBinContent(i)/hnewCuFramecs137M2->GetBinWidth(i));
    hAdapCuFramek40M2->SetBinContent(i, dBinSize * hnewCuFramek40M2->GetBinContent(i)/hnewCuFramek40M2->GetBinWidth(i));
    hAdapCuFramemn54M2->SetBinContent(i, dBinSize * hnewCuFramemn54M2->GetBinContent(i)/hnewCuFramemn54M2->GetBinWidth(i));
    hAdapCuFramepb210M2->SetBinContent(i, dBinSize * hnewCuFramepb210M2->GetBinContent(i)/hnewCuFramepb210M2->GetBinWidth(i));
    hAdapCuFrameth232M2->SetBinContent(i, dBinSize * hnewCuFrameth232M2->GetBinContent(i)/hnewCuFrameth232M2->GetBinWidth(i));
    hAdapCuFrameu238M2->SetBinContent(i, dBinSize * hnewCuFrameu238M2->GetBinContent(i)/hnewCuFrameu238M2->GetBinWidth(i));

    hAdapCuBoxco58M2->SetBinContent(i, dBinSize * hnewCuBoxco58M2->GetBinContent(i)/hnewCuBoxco58M2->GetBinWidth(i));
    hAdapCuBoxco60M2->SetBinContent(i, dBinSize * hnewCuBoxco60M2->GetBinContent(i)/hnewCuBoxco60M2->GetBinWidth(i));
    hAdapCuBoxcs137M2->SetBinContent(i, dBinSize * hnewCuBoxcs137M2->GetBinContent(i)/hnewCuBoxcs137M2->GetBinWidth(i));
    hAdapCuBoxk40M2->SetBinContent(i, dBinSize * hnewCuBoxk40M2->GetBinContent(i)/hnewCuBoxk40M2->GetBinWidth(i));
    hAdapCuBoxmn54M2->SetBinContent(i, dBinSize * hnewCuBoxmn54M2->GetBinContent(i)/hnewCuBoxmn54M2->GetBinWidth(i));
    hAdapCuBoxpb210M2->SetBinContent(i, dBinSize * hnewCuBoxpb210M2->GetBinContent(i)/hnewCuBoxpb210M2->GetBinWidth(i));
    hAdapCuBoxth232M2->SetBinContent(i, dBinSize * hnewCuBoxth232M2->GetBinContent(i)/hnewCuBoxth232M2->GetBinWidth(i));
    hAdapCuBoxu238M2->SetBinContent(i, dBinSize * hnewCuBoxu238M2->GetBinContent(i)/hnewCuBoxu238M2->GetBinWidth(i));

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
  qtree->Add("/Users/brian/macros/Simulations/Toy/combi1/combi1.root"); 
  qtree->Project("fDataHistoTot", "Ener2");
  qtree->Project("fDataHistoM1",  "Ener2", "Multiplicity == 1");
  qtree->Project("fDataHistoM2",  "Ener2", "Multiplicity == 2");

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


// Prints parameters, needs update 11-06-2014
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
	if(index > 61) cout << "Index too large" << endl;
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
  dNumCalls++;
  if(dNumCalls%50==0)
  {
    cout << "Call #: "<< dNumCalls << endl;
  }

  /////////////////////////////////////
  //// Many parameters
  ////////////////////////////////////
  // M1
  fModelTotM1->Add( hTeO20nuM1,      fParameters[0]);
  fModelTotM1->Add( hTeO22nuM1,      fParameters[1]);
  fModelTotM1->Add( hTeO2co60M1,     fParameters[2]);
  fModelTotM1->Add( hTeO2k40M1,      fParameters[3]);
  fModelTotM1->Add( hTeO2pb210M1,    fParameters[4]);
  fModelTotM1->Add( hTeO2po210M1,    fParameters[5]);
  fModelTotM1->Add( hTeO2te125M1,    fParameters[6]);
  fModelTotM1->Add( hTeO2th232M1,    fParameters[7]);
  fModelTotM1->Add( hTeO2th228M1,    fParameters[8]);
  fModelTotM1->Add( hTeO2ra226M1,    fParameters[9]);
  fModelTotM1->Add( hTeO2rn222M1,    fParameters[10]);
  fModelTotM1->Add( hTeO2u238M1,     fParameters[11]);
  fModelTotM1->Add( hTeO2th230M1,    fParameters[12]);
  fModelTotM1->Add( hTeO2u234M1,     fParameters[13]);

  fModelTotM1->Add( hCuFrameco58M1,  fParameters[14]);
  fModelTotM1->Add( hCuFrameco60M1,  fParameters[15]);
  fModelTotM1->Add( hCuFramecs137M1, fParameters[16]);
  fModelTotM1->Add( hCuFramek40M1,   fParameters[17]);
  fModelTotM1->Add( hCuFramemn54M1,  fParameters[18]);
  fModelTotM1->Add( hCuFramepb210M1, fParameters[19]);
  fModelTotM1->Add( hCuFrameth232M1, fParameters[20]);
  fModelTotM1->Add( hCuFrameu238M1,  fParameters[21]);

  fModelTotM1->Add( hCuBoxco58M1,    fParameters[22]);
  fModelTotM1->Add( hCuBoxco60M1,    fParameters[23]);
  fModelTotM1->Add( hCuBoxcs137M1,   fParameters[24]);
  fModelTotM1->Add( hCuBoxk40M1,     fParameters[25]);
  fModelTotM1->Add( hCuBoxmn54M1,    fParameters[26]);
  fModelTotM1->Add( hCuBoxpb210M1,   fParameters[27]);
  fModelTotM1->Add( hCuBoxth232M1,   fParameters[28]);
  fModelTotM1->Add( hCuBoxu238M1,    fParameters[29]);

  fModelTotM1->Add( h50mKco58M1,     fParameters[30]);
  fModelTotM1->Add( h50mKco60M1,     fParameters[31]);
  fModelTotM1->Add( h50mKcs137M1,    fParameters[32]);
  fModelTotM1->Add( h50mKk40M1,      fParameters[33]);
  fModelTotM1->Add( h50mKmn54M1,     fParameters[34]);
  fModelTotM1->Add( h50mKpb210M1,    fParameters[35]);
  fModelTotM1->Add( h50mKth232M1,    fParameters[36]);
  fModelTotM1->Add( h50mKu238M1,     fParameters[37]);

  fModelTotM1->Add( h600mKco60M1,    fParameters[38]);
  fModelTotM1->Add( h600mKk40M1,     fParameters[39]);
  fModelTotM1->Add( h600mKth232M1,   fParameters[40]);
  fModelTotM1->Add( h600mKu238M1,    fParameters[41]);

  fModelTotM1->Add( hPbRombi207M1,   fParameters[42]);
  fModelTotM1->Add( hPbRomco60M1,    fParameters[43]);
  fModelTotM1->Add( hPbRomcs137M1,   fParameters[44]);
  fModelTotM1->Add( hPbRomk40M1,     fParameters[45]);
  fModelTotM1->Add( hPbRompb210M1,   fParameters[46]);
  fModelTotM1->Add( hPbRomth232M1,   fParameters[47]);
  fModelTotM1->Add( hPbRomu238M1,    fParameters[48]);

  fModelTotM1->Add( hMBco60M1,       fParameters[49]);
  fModelTotM1->Add( hMBk40M1,        fParameters[50]);
  fModelTotM1->Add( hMBth232M1,      fParameters[51]);
  fModelTotM1->Add( hMBu238M1,       fParameters[52]);

  fModelTotM1->Add( hIVCco60M1,      fParameters[53]);
  fModelTotM1->Add( hIVCk40M1,       fParameters[54]);
  fModelTotM1->Add( hIVCth232M1,     fParameters[55]);
  fModelTotM1->Add( hIVCu238M1,      fParameters[56]);

  fModelTotM1->Add( hOVCco60M1,      fParameters[57]);
  fModelTotM1->Add( hOVCk40M1,       fParameters[58]);
  fModelTotM1->Add( hOVCth232M1,     fParameters[59]);
  fModelTotM1->Add( hOVCu238M1,      fParameters[60]);

// M2
  fModelTotM2->Add( hTeO20nuM2,      fParameters[0]);
  fModelTotM2->Add( hTeO22nuM2,      fParameters[1]);
  fModelTotM2->Add( hTeO2co60M2,     fParameters[2]);
  fModelTotM2->Add( hTeO2k40M2,      fParameters[3]);
  fModelTotM2->Add( hTeO2pb210M2,    fParameters[4]);
  fModelTotM2->Add( hTeO2po210M2,    fParameters[5]);
  fModelTotM2->Add( hTeO2te125M2,    fParameters[6]);
  fModelTotM2->Add( hTeO2th232M2,    fParameters[7]);
  fModelTotM2->Add( hTeO2th228M2,    fParameters[8]);
  fModelTotM2->Add( hTeO2ra226M2,    fParameters[9]);
  fModelTotM2->Add( hTeO2rn222M2,    fParameters[10]);
  fModelTotM2->Add( hTeO2u238M2,     fParameters[11]);
  fModelTotM2->Add( hTeO2th230M2,    fParameters[12]);
  fModelTotM2->Add( hTeO2u234M2,     fParameters[13]);

  fModelTotM2->Add( hCuFrameco58M2,  fParameters[14]);
  fModelTotM2->Add( hCuFrameco60M2,  fParameters[15]);
  fModelTotM2->Add( hCuFramecs137M2, fParameters[16]);
  fModelTotM2->Add( hCuFramek40M2,   fParameters[17]);
  fModelTotM2->Add( hCuFramemn54M2,  fParameters[18]);
  fModelTotM2->Add( hCuFramepb210M2, fParameters[19]);
  fModelTotM2->Add( hCuFrameth232M2, fParameters[20]);
  fModelTotM2->Add( hCuFrameu238M2,  fParameters[21]);

  fModelTotM2->Add( hCuBoxco58M2,    fParameters[22]);
  fModelTotM2->Add( hCuBoxco60M2,    fParameters[23]);
  fModelTotM2->Add( hCuBoxcs137M2,   fParameters[24]);
  fModelTotM2->Add( hCuBoxk40M2,     fParameters[25]);
  fModelTotM2->Add( hCuBoxmn54M2,    fParameters[26]);
  fModelTotM2->Add( hCuBoxpb210M2,   fParameters[27]);
  fModelTotM2->Add( hCuBoxth232M2,   fParameters[28]);
  fModelTotM2->Add( hCuBoxu238M2,    fParameters[29]);

  fModelTotM2->Add( h50mKco58M2,     fParameters[30]);
  fModelTotM2->Add( h50mKco60M2,     fParameters[31]);
  fModelTotM2->Add( h50mKcs137M2,    fParameters[32]);
  fModelTotM2->Add( h50mKk40M2,      fParameters[33]);
  fModelTotM2->Add( h50mKmn54M2,     fParameters[34]);
  fModelTotM2->Add( h50mKpb210M2,    fParameters[35]);
  fModelTotM2->Add( h50mKth232M2,    fParameters[36]);
  fModelTotM2->Add( h50mKu238M2,     fParameters[37]);

  fModelTotM2->Add( h600mKco60M2,    fParameters[38]);
  fModelTotM2->Add( h600mKk40M2,     fParameters[39]);
  fModelTotM2->Add( h600mKth232M2,   fParameters[40]);
  fModelTotM2->Add( h600mKu238M2,    fParameters[41]);

  fModelTotM2->Add( hPbRombi207M2,   fParameters[42]);
  fModelTotM2->Add( hPbRomco60M2,    fParameters[43]);
  fModelTotM2->Add( hPbRomcs137M2,   fParameters[44]);
  fModelTotM2->Add( hPbRomk40M2,     fParameters[45]);
  fModelTotM2->Add( hPbRompb210M2,   fParameters[46]);
  fModelTotM2->Add( hPbRomth232M2,   fParameters[47]);
  fModelTotM2->Add( hPbRomu238M2,    fParameters[48]);

  fModelTotM2->Add( hMBco60M2,       fParameters[49]);
  fModelTotM2->Add( hMBk40M2,        fParameters[50]);
  fModelTotM2->Add( hMBth232M2,      fParameters[51]);
  fModelTotM2->Add( hMBu238M2,       fParameters[52]);

  fModelTotM2->Add( hIVCco60M2,      fParameters[53]);
  fModelTotM2->Add( hIVCk40M2,       fParameters[54]);
  fModelTotM2->Add( hIVCth232M2,     fParameters[55]);
  fModelTotM2->Add( hIVCu238M2,      fParameters[56]);

  fModelTotM2->Add( hOVCco60M2,      fParameters[57]);
  fModelTotM2->Add( hOVCk40M2,       fParameters[58]);
  fModelTotM2->Add( hOVCth232M2,     fParameters[59]);
  fModelTotM2->Add( hOVCu238M2,      fParameters[60]);  

}

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
  fModelTotAdapM1->Add( hAdapTeO20nuM1,      fParameters[0]);
  fModelTotAdapM1->Add( hAdapTeO22nuM1,      fParameters[1]);
  fModelTotAdapM1->Add( hAdapTeO2co60M1,     fParameters[2]);
  fModelTotAdapM1->Add( hAdapTeO2k40M1,      fParameters[3]);
  fModelTotAdapM1->Add( hAdapTeO2pb210M1,    fParameters[4]);
  fModelTotAdapM1->Add( hAdapTeO2po210M1,    fParameters[5]);
  fModelTotAdapM1->Add( hAdapTeO2te125M1,    fParameters[6]);
  fModelTotAdapM1->Add( hAdapTeO2th232M1,    fParameters[7]);
  fModelTotAdapM1->Add( hAdapTeO2th228M1,    fParameters[8]);
  fModelTotAdapM1->Add( hAdapTeO2ra226M1,    fParameters[9]);
  fModelTotAdapM1->Add( hAdapTeO2rn222M1,    fParameters[10]);
  fModelTotAdapM1->Add( hAdapTeO2u238M1,     fParameters[11]);
  fModelTotAdapM1->Add( hAdapTeO2th230M1,    fParameters[12]);
  fModelTotAdapM1->Add( hAdapTeO2u234M1,     fParameters[13]);

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

// M2
  fModelTotAdapM2->Add( hAdapTeO20nuM2,      fParameters[0]);
  fModelTotAdapM2->Add( hAdapTeO22nuM2,      fParameters[1]);
  fModelTotAdapM2->Add( hAdapTeO2co60M2,     fParameters[2]);
  fModelTotAdapM2->Add( hAdapTeO2k40M2,      fParameters[3]);
  fModelTotAdapM2->Add( hAdapTeO2pb210M2,    fParameters[4]);
  fModelTotAdapM2->Add( hAdapTeO2po210M2,    fParameters[5]);
  fModelTotAdapM2->Add( hAdapTeO2te125M2,    fParameters[6]);
  fModelTotAdapM2->Add( hAdapTeO2th232M2,    fParameters[7]);
  fModelTotAdapM2->Add( hAdapTeO2th228M2,    fParameters[8]);
  fModelTotAdapM2->Add( hAdapTeO2ra226M2,    fParameters[9]);
  fModelTotAdapM2->Add( hAdapTeO2rn222M2,    fParameters[10]);
  fModelTotAdapM2->Add( hAdapTeO2u238M2,     fParameters[11]);
  fModelTotAdapM2->Add( hAdapTeO2th230M2,    fParameters[12]);
  fModelTotAdapM2->Add( hAdapTeO2u234M2,     fParameters[13]);

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


   TMinuit minuit(61);

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
   minuit.DefineParameter(0, "hTeO20nuM1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(1, "hTeO22nuM1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(2, "hTeO2co60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(3, "hTeO2k40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(4, "hTeO2pb210M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(5, "hTeO2po210M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(6, "hTeO2te125M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(7, "hTeO2th232M1",  10., 10.0, 0., 1000000);
   minuit.DefineParameter(8, "hTeO2th228M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(9, "hTeO2ra226M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(10, "hTeO2rn222M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(11, "hTeO2u238M1",  10., 10.0, 0., 1000000);
   minuit.DefineParameter(12, "hTeO2th230M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(13, "hTeO2u234M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(14, "hCuFrameco58M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(15, "hCuFrameco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(16, "hCuFramecs137M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(17, "hCuFramek40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(18, "hCuFramemn54M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(19, "hCuFramepb210M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(20, "hCuFrameth232M1",  10., 10.0, 0., 1000000);
   minuit.DefineParameter(21, "hCuFrameu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(22, "hCuBoxco58M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(23, "hCuBoxco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(24, "hCuBoxcs137M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(25, "hCuBoxk40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(26, "hCuBoxmn54M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(27, "hCuBoxpb210M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(28, "hCuBoxth232M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(29, "hCuBoxu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(30, "h50mKco58M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(31, "h50mKco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(32, "h50mKcs137M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(33, "h50mKk40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(34, "h50mKmn54M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(35, "h50mKpb210M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(36, "h50mKth232M1",  10., 10.0, 0., 1000000);
   minuit.DefineParameter(37, "h50mKu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(38, "h600mKco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(39, "h600mKk40M1",  0., 10.0, 0., 1000000); 
   minuit.DefineParameter(40, "h600mKth232M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(41, "h600mKu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(42, "hPbRombi207M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(43, "hPbRomco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(44, "hPbRomcs137M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(45, "hPbRomk40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(46, "hPbRompb210M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(47, "hPbRomth232M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(48, "hPbRomu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(49, "hMBco60M1",  0., 10.0, 0., 1000000); 
   minuit.DefineParameter(50, "hMBk40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(51, "hMBth232M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(52, "hMBu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(53, "hIVCco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(54, "hIVCk40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(55, "hIVCth232M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(56, "hIVCu238M1",  0., 10.0, 0., 1000000);

   minuit.DefineParameter(57, "hOVCco60M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(58, "hOVCk40M1",  0., 10.0, 0., 1000000);
   minuit.DefineParameter(59, "hOVCth232M1",  0., 10.0, 0., 1000000);    
   minuit.DefineParameter(60, "hOVCu238M1",  0., 10.0, 0., 1000000);

   // Fix parameters here
   minuit.FixParameter(0); // TeO2 0nu
   minuit.FixParameter(1); // TeO2 2nu
   minuit.FixParameter(2); // TeO2 co60
   minuit.FixParameter(3); // TeO2 k40
   minuit.FixParameter(4); // TeO2 pb210
   minuit.FixParameter(5); // TeO2 po210
   minuit.FixParameter(6); // TeO2 te125m
   // minuit.FixParameter(7); // TeO2 th232
   minuit.FixParameter(8); // TeO2 th232-th228
   minuit.FixParameter(9); // TeO2 u238-ra226
   minuit.FixParameter(10); // TeO2 u238-rn222
   // minuit.FixParameter(11); // TeO2 u238
   minuit.FixParameter(12); // TeO2 u238-th230
   minuit.FixParameter(13); // TeO2 u238-u234
   minuit.FixParameter(14); // Frame co58
   minuit.FixParameter(15); // Frame co60
   minuit.FixParameter(16); // Frame cs137
   minuit.FixParameter(17); // Frame k40
   minuit.FixParameter(18); // Frame mn54
   minuit.FixParameter(19); // Frame pb210
   // minuit.FixParameter(20); // Frame th232
   minuit.FixParameter(21); // Frame u238
   minuit.FixParameter(22); // CuBox co58
   minuit.FixParameter(23); // CuBox co60
   minuit.FixParameter(24); // CuBox cs137
   minuit.FixParameter(25); // CuBox k40
   minuit.FixParameter(26); // CuBox mn54
   minuit.FixParameter(27); // CuBox pb210
   minuit.FixParameter(28); // CuBox th232
   minuit.FixParameter(29); // CuBox u238
   minuit.FixParameter(30); // 50mK co58
   minuit.FixParameter(31); // 50mK co60
   minuit.FixParameter(32); // 50mK cs137
   minuit.FixParameter(33); // 50mK k40
   minuit.FixParameter(34); // 50mK mn54
   minuit.FixParameter(35); // 50mK pb210
   // minuit.FixParameter(36); // 50mK th232
   minuit.FixParameter(37); // 50mK u238
   minuit.FixParameter(38); // 600mK co60
   minuit.FixParameter(39); // 600mK k40
   minuit.FixParameter(40); // 600mK th232
   minuit.FixParameter(41); // 600mK u238
   minuit.FixParameter(42); // RLead bi207
   minuit.FixParameter(43); // RLead co60
   minuit.FixParameter(44); // RLead cs137
   minuit.FixParameter(45); // RLead k40
   minuit.FixParameter(46); // RLead pb210
   minuit.FixParameter(47); // RLead th232
   minuit.FixParameter(48); // RLead u238
   minuit.FixParameter(49); // MB co60
   minuit.FixParameter(50); // MB k40
   minuit.FixParameter(51); // MB th232
   minuit.FixParameter(52); // MB u238
   minuit.FixParameter(53); // IVC co60
   minuit.FixParameter(54); // IVC k40
   minuit.FixParameter(55); // IVC th232
   minuit.FixParameter(56); // IVC u238
   minuit.FixParameter(57); // OVC co60
   minuit.FixParameter(58); // OVC k40
   minuit.FixParameter(59); // OVC th232
   minuit.FixParameter(60); // OVC u238
   // Number of Parameters (for Chi-squared/NDF calculation)
   int dNumParameters = 4;

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
  minuit.GetParameter(26,  fParameters[26],   fParError[26]);
  minuit.GetParameter(27,  fParameters[27],   fParError[27]);
  minuit.GetParameter(28,  fParameters[28],   fParError[28]);
  minuit.GetParameter(29,  fParameters[29],   fParError[29]);
  minuit.GetParameter(30,  fParameters[30],   fParError[30]);
  minuit.GetParameter(31,  fParameters[31],   fParError[31]);
  minuit.GetParameter(32,  fParameters[32],   fParError[32]);
  minuit.GetParameter(33,  fParameters[33],   fParError[33]);
  minuit.GetParameter(34,  fParameters[34],   fParError[34]);
  minuit.GetParameter(35,  fParameters[35],   fParError[35]);
  minuit.GetParameter(36,  fParameters[36],   fParError[36]);
  minuit.GetParameter(37,  fParameters[37],   fParError[37]);
  minuit.GetParameter(38,  fParameters[38],   fParError[38]);
  minuit.GetParameter(39,  fParameters[39],   fParError[39]);
  minuit.GetParameter(40,  fParameters[40],   fParError[40]);
  minuit.GetParameter(41,  fParameters[41],   fParError[41]);
  minuit.GetParameter(42,  fParameters[42],   fParError[42]);
  minuit.GetParameter(43,  fParameters[43],   fParError[43]);
  minuit.GetParameter(44,  fParameters[44],   fParError[44]);
  minuit.GetParameter(45,  fParameters[45],   fParError[45]);
  minuit.GetParameter(46,  fParameters[46],   fParError[46]);
  minuit.GetParameter(47,  fParameters[47],   fParError[47]);
  minuit.GetParameter(48,  fParameters[48],   fParError[48]);
  minuit.GetParameter(49,  fParameters[49],   fParError[49]);
  minuit.GetParameter(50,  fParameters[50],   fParError[50]);
  minuit.GetParameter(51,  fParameters[51],   fParError[51]);
  minuit.GetParameter(52,  fParameters[52],   fParError[52]);
  minuit.GetParameter(53,  fParameters[53],   fParError[53]);
  minuit.GetParameter(54,  fParameters[54],   fParError[54]);
  minuit.GetParameter(55,  fParameters[55],   fParError[55]);
  minuit.GetParameter(56,  fParameters[56],   fParError[56]);
  minuit.GetParameter(57,  fParameters[57],   fParError[57]);
  minuit.GetParameter(58,  fParameters[58],   fParError[58]);
  minuit.GetParameter(59,  fParameters[59],   fParError[59]);
  minuit.GetParameter(60,  fParameters[60],   fParError[60]);
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

  fModelTotAdapuM1->Add( hAdapTeO2u238M1,       fParameters[11] );
  fModelTotAdapuM1->Add( hAdapCuFrameu238M1,    fParameters[21] );
  fModelTotAdapuM1->Add( hAdapCuBoxu238M1,      fParameters[29] );
  fModelTotAdapuM1->Add( hAdap50mKu238M1,       fParameters[37] );
  fModelTotAdapuM1->Add( hAdap600mKu238M1,      fParameters[41] );
  fModelTotAdapuM1->Add( hAdapPbRomu238M1,      fParameters[48] );
  fModelTotAdapuM1->Add( hAdapMBu238M1,         fParameters[52] );
  fModelTotAdapuM1->Add( hAdapIVCu238M1,        fParameters[56] );
  fModelTotAdapuM1->Add( hAdapOVCu238M1,        fParameters[60] );

  fModelTotAdapkM1->Add( hAdapTeO2k40M1,        fParameters[3]  );
  fModelTotAdapkM1->Add( hAdapCuFramek40M1,     fParameters[17] );
  fModelTotAdapkM1->Add( hAdapCuBoxk40M1,       fParameters[25] );
  fModelTotAdapkM1->Add( hAdap50mKk40M1,        fParameters[33] );
  fModelTotAdapkM1->Add( hAdap600mKk40M1,       fParameters[39] );
  fModelTotAdapkM1->Add( hAdapPbRomk40M1,       fParameters[45] );
  fModelTotAdapkM1->Add( hAdapMBk40M1,          fParameters[50] );
  fModelTotAdapkM1->Add( hAdapIVCk40M1,         fParameters[54] );
  fModelTotAdapkM1->Add( hAdapOVCk40M1,         fParameters[58] ); 

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
  fModelTotAdappbM1->Add( hAdapCuFramepb210M1,  fParameters[15] );
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

  fModelTotAdapteo2M1->Add( hAdapTeO2po210M1,   fParameters[5] );
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

  fModelTotAdapuM2->Add( hAdapTeO2u238M2,       fParameters[11] );
  fModelTotAdapuM2->Add( hAdapCuFrameu238M2,    fParameters[21] );
  fModelTotAdapuM2->Add( hAdapCuBoxu238M2,      fParameters[29] );
  fModelTotAdapuM2->Add( hAdap50mKu238M2,       fParameters[37] );
  fModelTotAdapuM2->Add( hAdap600mKu238M2,      fParameters[41] );
  fModelTotAdapuM2->Add( hAdapPbRomu238M2,      fParameters[48] );
  fModelTotAdapuM2->Add( hAdapMBu238M2,         fParameters[52] );
  fModelTotAdapuM2->Add( hAdapIVCu238M2,        fParameters[56] );
  fModelTotAdapuM2->Add( hAdapOVCu238M2,        fParameters[60] );

  fModelTotAdapkM2->Add( hAdapTeO2k40M2,        fParameters[3]  );
  fModelTotAdapkM2->Add( hAdapCuFramek40M2,     fParameters[17] );
  fModelTotAdapkM2->Add( hAdapCuBoxk40M2,       fParameters[25] );
  fModelTotAdapkM2->Add( hAdap50mKk40M2,        fParameters[33] );
  fModelTotAdapkM2->Add( hAdap600mKk40M2,       fParameters[39] );
  fModelTotAdapkM2->Add( hAdapPbRomk40M2,       fParameters[45] );
  fModelTotAdapkM2->Add( hAdapMBk40M2,          fParameters[50] );
  fModelTotAdapkM2->Add( hAdapIVCk40M2,         fParameters[54] );
  fModelTotAdapkM2->Add( hAdapOVCk40M2,         fParameters[58] ); 

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
  fModelTotAdappbM2->Add( hAdapCuFramepb210M2,  fParameters[15] );
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

  fModelTotAdapteo2M2->Add( hAdapTeO2po210M2,   fParameters[5] );
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
  fModelTotAdapmnM1->SetLineColor(40);
  fModelTotAdapmnM1->SetLineStyle(2);


  fModelTotAdappbM1->SetLineStyle(2);
  fModelTotAdappbM1->SetLineColor(38);

  fModelTotAdapthM1->Draw("SAME");
  fModelTotAdapuM1->Draw("SAME");
  fModelTotAdapkM1->Draw("SAME");
  fModelTotAdapcoM1->Draw("SAME");
  fModelTotAdapNDBDM1->Draw("SAME");
  fModelTotAdap2NDBDM1->Draw("SAME");
  fModelTotAdapbiM1->Draw("SAME");
  fModelTotAdapmnM1->Draw("SAME");

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
  legfit1->AddEntry(fModelTotAdapbiM1, "bi-207", "l");
  legfit1->AddEntry(fModelTotAdapmnM1, "mn-54", "l");
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

  fModelTotAdapthM2->Draw("SAME");
  fModelTotAdapuM2->Draw("SAME");
  fModelTotAdapkM2->Draw("SAME");
  fModelTotAdapcoM2->Draw("SAME");
  fModelTotAdapNDBDM2->Draw("SAME");
  fModelTotAdap2NDBDM2->Draw("SAME");
  fModelTotAdapbiM2->Draw("SAME");    
  fModelTotAdapmnM2->Draw("SAME"); 

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
  legfit2->AddEntry(fModelTotAdapbiM2, "bi-207", "l");
  legfit2->AddEntry(fModelTotAdapmnM2, "mn-54", "l");

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


