#include "TMinuit.h"
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
  Obj->SetParameters(61, x[61]);
  Obj->SetParameters(62, x[62]);
  Obj->SetParameters(63, x[63]); 

  // Surface
  /*
  Obj->SetParameters(64, x[64]);
  Obj->SetParameters(65, x[65]);
  Obj->SetParameters(66, x[66]);
  Obj->SetParameters(67, x[67]);
  Obj->SetParameters(68, x[68]);
  Obj->SetParameters(69, x[69]);  
  Obj->SetParameters(70, x[70]);
  Obj->SetParameters(71, x[71]);
  Obj->SetParameters(72, x[72]);
  Obj->SetParameters(73, x[73]);
  Obj->SetParameters(74, x[74]);
  Obj->SetParameters(75, x[75]);
  Obj->SetParameters(76, x[76]);
  Obj->SetParameters(77, x[77]);
  Obj->SetParameters(78, x[78]);
  Obj->SetParameters(79, x[79]);
  Obj->SetParameters(80, x[80]);
  Obj->SetParameters(81, x[81]);
  Obj->SetParameters(82, x[82]);
  Obj->SetParameters(83, x[83]);
  Obj->SetParameters(84, x[84]);
  Obj->SetParameters(85, x[85]);
  Obj->SetParameters(86, x[86]);
  Obj->SetParameters(87, x[87]);
  Obj->SetParameters(88, x[88]);
  Obj->SetParameters(89, x[89]);
  Obj->SetParameters(90, x[90]);
  Obj->SetParameters(91, x[91]);
  Obj->SetParameters(92, x[92]);
  Obj->SetParameters(93, x[93]);
  Obj->SetParameters(94, x[94]);
  Obj->SetParameters(95, x[95]);
  Obj->SetParameters(96, x[96]);
  Obj->SetParameters(97, x[97]);
  Obj->SetParameters(98, x[98]);
  Obj->SetParameters(99, x[99]);
  Obj->SetParameters(100, x[100]);
  Obj->SetParameters(101, x[101]);
  Obj->SetParameters(102, x[102]);
  Obj->SetParameters(103, x[103]);
  Obj->SetParameters(104, x[104]);
  Obj->SetParameters(105, x[105]);
  Obj->SetParameters(106, x[106]);
  Obj->SetParameters(107, x[107]);
  Obj->SetParameters(108, x[108]);
  Obj->SetParameters(109, x[109]);
  Obj->SetParameters(110, x[110]);
  Obj->SetParameters(111, x[111]);
  */

  // Implement a method in your class that calculates the quantity you want to minimize, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimise fval
    Obj->UpdateModel();
    fval = Obj->GetChiSquare();
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
  qtree = new TChain("qredtree");
  // qtree = new TChain("CombiTree");
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

  TH1D  *hOut       = new TH1D("hOut", "Fit Residuals", dAdap);


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
   TMinuit minuit(112); // for more parameters

   // Reduce Minuit Output
   minuit.SetPrintLevel(1);
   // minuit.Command("SET MINImize 10000 0.001");
   minuit.Command("SET STRategy 2");
  //minuit.Command("SET IMProve 1000 ");
   minuit.SetMaxIterations(10000);
   // minuit.Command("SET MIGrad 10000")
   minuit.SetObjectFit(this); //see the external FCN  above
   
   //define the parameters and set the ranges and initial guesses see ROOTs TMinuit documentation
   ////////////////////////////////////////////////
   // Using more parameters
   ////////////////////////////////////////////////
   minuit.DefineParameter(0, "TeO2 0nu",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(1, "TeO2 2nu",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(2, "TeO2 co60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(3, "TeO2 k40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(4, "TeO2 pb210",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(5, "TeO2 po210",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(6, "TeO2 te125",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(7, "TeO2 th232",  0.1, 0.1, 0., 1000000);
   minuit.DefineParameter(8, "TeO2 th228",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(9, "TeO2 ra226",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(10, "TeO2 rn222",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(11, "TeO2 u238",  0.1, 0.1, 0., 1000000);
   minuit.DefineParameter(12, "TeO2 th230",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(13, "TeO2 u234",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(14, "CuFrame co58",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(15, "CuFrame co60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(16, "CuFrame cs137",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(17, "CuFrame k40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(18, "CuFrame mn54",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(19, "CuFrame pb210",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(20, "CuFrame th232",  0.1, 0.1, 0., 1000000);
   minuit.DefineParameter(21, "CuFrame u238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(22, "CuBoxco58",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(23, "CuBoxco60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(24, "CuBoxcs137",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(25, "CuBoxk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(26, "CuBoxmn54",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(27, "CuBoxpb210",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(28, "CuBoxth232",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(29, "CuBoxu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(30, "50mKco58",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(31, "50mKco60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(32, "50mKcs137",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(33, "50mKk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(34, "50mKmn54",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(35, "50mKpb210",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(36, "50mKth232",  0.1, 0.1, 0., 1000000);
   minuit.DefineParameter(37, "50mKu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(38, "600mKco60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(39, "600mKk40",  0., 0.1, 0., 1000000); 
   minuit.DefineParameter(40, "600mKth232",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(41, "600mKu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(42, "PbRombi207",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(43, "PbRomco60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(44, "PbRomcs137",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(45, "PbRomk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(46, "PbRompb210",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(47, "PbRomth232",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(48, "PbRomu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(49, "MBco60",  0., 0.1, 0., 1000000); 
   minuit.DefineParameter(50, "MBk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(51, "MBth232",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(52, "MBu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(53, "IVCco60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(54, "IVCk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(55, "IVCth232",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(56, "IVCu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(57, "OVCco60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(58, "OVCk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(59, "OVCth232",  0., 0.1, 0., 1000000);    
   minuit.DefineParameter(60, "OVCu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(61, "SIk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(62, "SIth232",  0., 0.1, 0., 1000000);    
   minuit.DefineParameter(63, "SIu238",  0., 0.1, 0., 1000000);

/*
   minuit.DefineParameter(64, "TeO2Spb210_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(65, "TeO2Spo210_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(66, "TeO2Spo210_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(67, "TeO2Sth232_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(68, "TeO2Su238_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(69, "TeO2Sxpb210_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(70, "TeO2Sxpb210_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(71, "TeO2Sxpb210_1",  0.03, 0.1, 0., 1000000);
   minuit.DefineParameter(72, "TeO2Sxpb210_10",  0., 0.1, 0., 1000000);    
   minuit.DefineParameter(73, "TeO2Sxpo210_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(74, "TeO2Sxpo210_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(75, "TeO2Sxpo210_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(76, "TeO2Sxth232_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(77, "TeO2Sxth232_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(78, "TeO2Sxth232_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(79, "TeO2Sxth232_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(80, "TeO2Sxu238_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(81, "TeO2Sxu238_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(82, "TeO2Sxu238_1",  0., 0.1, 0., 1000000);   
   minuit.DefineParameter(83, "TeO2Sxu238_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(84, "CuFrameSth232_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(85, "CuFrameSu238_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(86, "CuFrameSxpb210_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(87, "CuFrameSxpb210_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(88, "CuFrameSxpb210_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(89, "CuFrameSxpb210_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(90, "CuFrameSxth232_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(91, "CuFrameSxth232_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(92, "CuFrameSxth232_1",  0., 0.1, 0., 1000000);   
   minuit.DefineParameter(93, "CuFrameSxth232_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(94, "CuFrameSxu238_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(95, "CuFrameSxu238_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(96, "CuFrameSxu238_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(97, "CuFrameSxu238_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(98, "CuBoxSth232_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(99, "CuBoxSu238_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(100, "CuBoxSxpb210_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(101, "CuBoxSxpb210_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(102, "CuBoxSxpb210_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(103, "CuBoxSxpb210_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(104, "CuBoxSxth232_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(105, "CuBoxSxth232_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(106, "CuBoxSxth232_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(107, "CuBoxSxth232_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(108, "CuBoxSxu238_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(109, "CuBoxSxu238_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(110, "CuBoxSxu238_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(111, "CuBoxSxu238_10",  0., 0.1, 0., 1000000);
*/

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
   minuit.FixParameter(61); // SI k40
   minuit.FixParameter(62); // SI th232
   minuit.FixParameter(63); // SI u238   
   /*
   minuit.FixParameter(64); // TeO2 S pb210 01
   minuit.FixParameter(65); // TeO2 S po210 001
   minuit.FixParameter(66); // TeO2 S po210 01
   minuit.FixParameter(67); // TeO2 S th232 01
   minuit.FixParameter(68); // TeO2 S u238 01
   minuit.FixParameter(69); // TeO2 Sx pb210 001
   minuit.FixParameter(70); // TeO2 Sx pb210 01
   // minuit.FixParameter(71); // TeO2 Sx pb210 1
   minuit.FixParameter(72); // TeO2 S pb210 10
   minuit.FixParameter(73); // TeO2 Sx po210 001
   minuit.FixParameter(74); // TeO2 Sx po210 01
   minuit.FixParameter(75); // TeO2 Sx po210 1
   minuit.FixParameter(76); // TeO2 Sx th232 001
   minuit.FixParameter(77); // TeO2 Sx th232 01
   minuit.FixParameter(78); // TeO2 Sx th232 1
   minuit.FixParameter(79); // TeO2 Sx th232 10
   minuit.FixParameter(80); // TeO2 Sx u238 001
   minuit.FixParameter(81); // TeO2 Sx u238 01
   minuit.FixParameter(82); // TeO2 Sx u238 1
   minuit.FixParameter(83); // TeO2 Sx u238 10
   minuit.FixParameter(84); // Frame S th232 1
   minuit.FixParameter(85); // Frame S u238 1
   minuit.FixParameter(86); // Frame Sx pb210 001
   minuit.FixParameter(87); // Frame Sx pb210 01
   minuit.FixParameter(88); // Frame Sx pb210 1
   minuit.FixParameter(89); // Frame Sx pb210 10
   minuit.FixParameter(90); // Frame Sx th232 001
   minuit.FixParameter(91); // Frame Sx th232 01
   minuit.FixParameter(92); // Frame Sx th232 1
   minuit.FixParameter(93); // Frame Sx th232 10
   minuit.FixParameter(94); // Frame Sx u238 001
   minuit.FixParameter(95); // Frame Sx u238 01
   minuit.FixParameter(96); // Frame Sx u238 1
   minuit.FixParameter(97); // Frame Sx u238 10
   minuit.FixParameter(98); // CuBox S th232 1
   minuit.FixParameter(99); // CuBox S u238 1
   minuit.FixParameter(100); // CuBox Sx pb210 001
   minuit.FixParameter(101); // CuBox Sx pb210 01
   minuit.FixParameter(102); // CuBox Sx pb210 1
   minuit.FixParameter(103); // CuBox Sx pb210 10 
   minuit.FixParameter(104); // CuBox Sx th232 001
   minuit.FixParameter(105); // CuBox Sx th232 01
   minuit.FixParameter(106); // CuBox Sx th232 1
   minuit.FixParameter(107); // CuBox Sx th232 10
   minuit.FixParameter(108); // CuBox Sx u238 001
   minuit.FixParameter(109); // CuBox Sx u238 01
   minuit.FixParameter(110); // CuBox Sx u238 1
   minuit.FixParameter(111); // CuBox Sx u238 10
   */
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
  minuit.GetParameter(61,  fParameters[61],   fParError[61]);
  minuit.GetParameter(62,  fParameters[62],   fParError[62]);
  minuit.GetParameter(63,  fParameters[63],   fParError[63]);
/*
  minuit.GetParameter(64,  fParameters[64],   fParError[64]);
  minuit.GetParameter(65,  fParameters[65],   fParError[65]);
  minuit.GetParameter(66,  fParameters[66],   fParError[66]);
  minuit.GetParameter(67,  fParameters[67],   fParError[67]);
  minuit.GetParameter(68,  fParameters[68],   fParError[68]);
  minuit.GetParameter(69,  fParameters[69],   fParError[69]);
  minuit.GetParameter(70,  fParameters[70],   fParError[70]);
  minuit.GetParameter(71,  fParameters[71],   fParError[71]);
  minuit.GetParameter(72,  fParameters[72],   fParError[72]);
  minuit.GetParameter(73,  fParameters[73],   fParError[73]);
  minuit.GetParameter(74,  fParameters[74],   fParError[74]);
  minuit.GetParameter(75,  fParameters[75],   fParError[75]);
  minuit.GetParameter(76,  fParameters[76],   fParError[76]);
  minuit.GetParameter(77,  fParameters[77],   fParError[77]);
  minuit.GetParameter(78,  fParameters[78],   fParError[78]);
  minuit.GetParameter(79,  fParameters[79],   fParError[79]);
  minuit.GetParameter(80,  fParameters[80],   fParError[80]);
  minuit.GetParameter(81,  fParameters[81],   fParError[81]);
  minuit.GetParameter(82,  fParameters[82],   fParError[82]);
  minuit.GetParameter(83,  fParameters[83],   fParError[83]);
  minuit.GetParameter(84,  fParameters[84],   fParError[84]);
  minuit.GetParameter(85,  fParameters[85],   fParError[85]);
  minuit.GetParameter(86,  fParameters[86],   fParError[86]);
  minuit.GetParameter(87,  fParameters[87],   fParError[87]);
  minuit.GetParameter(88,  fParameters[88],   fParError[88]);
  minuit.GetParameter(89,  fParameters[89],   fParError[89]);
  minuit.GetParameter(90,  fParameters[90],   fParError[90]);
  minuit.GetParameter(91,  fParameters[91],   fParError[91]);
  minuit.GetParameter(92,  fParameters[92],   fParError[92]);
  minuit.GetParameter(93,  fParameters[93],   fParError[93]);
  minuit.GetParameter(94,  fParameters[94],   fParError[94]);
  minuit.GetParameter(95,  fParameters[95],   fParError[95]);
  minuit.GetParameter(96,  fParameters[96],   fParError[96]);
  minuit.GetParameter(97,  fParameters[97],   fParError[97]);
  minuit.GetParameter(98,  fParameters[98],   fParError[98]);
  minuit.GetParameter(99,  fParameters[99],   fParError[99]);
  minuit.GetParameter(100,  fParameters[100],   fParError[100]);
  minuit.GetParameter(101,  fParameters[101],   fParError[101]);
  minuit.GetParameter(102,  fParameters[102],   fParError[102]);
  minuit.GetParameter(103,  fParameters[103],   fParError[103]);
  minuit.GetParameter(104,  fParameters[104],   fParError[104]);
  minuit.GetParameter(105,  fParameters[105],   fParError[105]);
  minuit.GetParameter(106,  fParameters[106],   fParError[106]);
  minuit.GetParameter(107,  fParameters[107],   fParError[107]);
  minuit.GetParameter(108,  fParameters[108],   fParError[108]);
  minuit.GetParameter(109,  fParameters[109],   fParError[109]);
  minuit.GetParameter(110,  fParameters[110],   fParError[110]);
  minuit.GetParameter(111,  fParameters[111],   fParError[111]);
  */
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
  fModelTotthM1->Add( hSIth232M1,       fParameters[62] );

  fModelTotuM1->Add( hTeO2u238M1,       fParameters[11] );
  fModelTotuM1->Add( hCuFrameu238M1,    fParameters[21] );
  fModelTotuM1->Add( hCuBoxu238M1,      fParameters[29] );
  fModelTotuM1->Add( h50mKu238M1,       fParameters[37] );
  fModelTotuM1->Add( h600mKu238M1,      fParameters[41] );
  fModelTotuM1->Add( hPbRomu238M1,      fParameters[48] );
  fModelTotuM1->Add( hMBu238M1,         fParameters[52] );
  fModelTotuM1->Add( hIVCu238M1,        fParameters[56] );
  fModelTotuM1->Add( hOVCu238M1,        fParameters[60] );
  fModelTotuM1->Add( hSIu238M1,         fParameters[63] );


  fModelTotkM1->Add( hTeO2k40M1,        fParameters[3]  );
  fModelTotkM1->Add( hCuFramek40M1,     fParameters[17] );
  fModelTotkM1->Add( hCuBoxk40M1,       fParameters[25] );
  fModelTotkM1->Add( h50mKk40M1,        fParameters[33] );
  fModelTotkM1->Add( h600mKk40M1,       fParameters[39] );
  fModelTotkM1->Add( hPbRomk40M1,       fParameters[45] );
  fModelTotkM1->Add( hMBk40M1,          fParameters[50] );
  fModelTotkM1->Add( hIVCk40M1,         fParameters[54] );
  fModelTotkM1->Add( hOVCk40M1,         fParameters[58] ); 
  fModelTotkM1->Add( hSIk40M1,          fParameters[61] );


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
  fModelTotpbM1->Add( hCuFramepb210M1,  fParameters[19] );
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

  fModelTotSpbM1->Add( hTeO2Spb210M1_01,      fParameters[64] );
  fModelTotSpbM1->Add( hTeO2Sxpb210M1_001,    fParameters[69] );
  fModelTotSpbM1->Add( hTeO2Sxpb210M1_01,     fParameters[70] );
  fModelTotSpbM1->Add( hTeO2Sxpb210M1_1,      fParameters[71] );
  fModelTotSpbM1->Add( hTeO2Sxpb210M1_10,     fParameters[72] );
  fModelTotSpbM1->Add( hCuFrameSxpb210M1_001, fParameters[86] );
  fModelTotSpbM1->Add( hCuFrameSxpb210M1_01,  fParameters[87] );
  fModelTotSpbM1->Add( hCuFrameSxpb210M1_1,   fParameters[88] );
  fModelTotSpbM1->Add( hCuFrameSxpb210M1_10,  fParameters[89] );
  fModelTotSpbM1->Add( hCuBoxSxpb210M1_001,   fParameters[100] );
  fModelTotSpbM1->Add( hCuBoxSxpb210M1_01,    fParameters[101] );
  fModelTotSpbM1->Add( hCuBoxSxpb210M1_1,     fParameters[102] );
  fModelTotSpbM1->Add( hCuBoxSxpb210M1_10,    fParameters[103] );

  fModelTotSpoM1->Add( hTeO2Spo210M1_001,     fParameters[65] );
  fModelTotSpoM1->Add( hTeO2Spo210M1_01,      fParameters[66] );
  fModelTotSpoM1->Add( hTeO2Sxpo210M1_001,    fParameters[73] );
  fModelTotSpoM1->Add( hTeO2Sxpo210M1_01,     fParameters[74] );
  fModelTotSpoM1->Add( hTeO2Sxpo210M1_1,      fParameters[75] );

  fModelTotSthM1->Add( hTeO2Sth232M1_01,      fParameters[67] );
  fModelTotSthM1->Add( hTeO2Sxth232M1_001,    fParameters[76] );
  fModelTotSthM1->Add( hTeO2Sxth232M1_01,     fParameters[77] );
  fModelTotSthM1->Add( hTeO2Sxth232M1_1,      fParameters[78] );
  fModelTotSthM1->Add( hTeO2Sxth232M1_10,     fParameters[79] );
  fModelTotSthM1->Add( hCuFrameSth232M1_1,    fParameters[84] );
  fModelTotSthM1->Add( hCuFrameSxth232M1_001, fParameters[90] );
  fModelTotSthM1->Add( hCuFrameSxth232M1_01,  fParameters[91] );
  fModelTotSthM1->Add( hCuFrameSxth232M1_1,   fParameters[92] );
  fModelTotSthM1->Add( hCuFrameSxth232M1_10,  fParameters[93] );
  fModelTotSthM1->Add( hCuBoxSth232M1_1,      fParameters[98] );
  fModelTotSthM1->Add( hCuBoxSxth232M1_001,   fParameters[104] );
  fModelTotSthM1->Add( hCuBoxSxth232M1_01,    fParameters[105] );
  fModelTotSthM1->Add( hCuBoxSxth232M1_1,     fParameters[106] );
  fModelTotSthM1->Add( hCuBoxSxth232M1_10,    fParameters[107] );

  fModelTotSuM1->Add( hTeO2Su238M1_01,        fParameters[68] );
  fModelTotSuM1->Add( hTeO2Sxu238M1_001,      fParameters[80] );
  fModelTotSuM1->Add( hTeO2Sxu238M1_01,       fParameters[81] );
  fModelTotSuM1->Add( hTeO2Sxu238M1_1,        fParameters[82] );
  fModelTotSuM1->Add( hTeO2Sxu238M1_10,       fParameters[83] );
  fModelTotSuM1->Add( hCuFrameSu238M1_1,      fParameters[85] );
  fModelTotSuM1->Add( hCuFrameSxu238M1_001,   fParameters[94] );
  fModelTotSuM1->Add( hCuFrameSxu238M1_01,    fParameters[95] );
  fModelTotSuM1->Add( hCuFrameSxu238M1_1,     fParameters[96] );
  fModelTotSuM1->Add( hCuFrameSxu238M1_10,    fParameters[97] );
  fModelTotSuM1->Add( hCuBoxSu238M1_1,        fParameters[99] );
  fModelTotSuM1->Add( hCuBoxSxu238M1_001,     fParameters[108] );
  fModelTotSuM1->Add( hCuBoxSxu238M1_01,      fParameters[109] );
  fModelTotSuM1->Add( hCuBoxSxu238M1_1,       fParameters[110] );
  fModelTotSuM1->Add( hCuBoxSxu238M1_10,      fParameters[111] );


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
  fModelTotthM2->Add( hSIth232M2,       fParameters[62] );

  fModelTotuM2->Add( hTeO2u238M2,       fParameters[11] );
  fModelTotuM2->Add( hCuFrameu238M2,    fParameters[21] );
  fModelTotuM2->Add( hCuBoxu238M2,      fParameters[29] );
  fModelTotuM2->Add( h50mKu238M2,       fParameters[37] );
  fModelTotuM2->Add( h600mKu238M2,      fParameters[41] );
  fModelTotuM2->Add( hPbRomu238M2,      fParameters[48] );
  fModelTotuM2->Add( hMBu238M2,         fParameters[52] );
  fModelTotuM2->Add( hIVCu238M2,        fParameters[56] );
  fModelTotuM2->Add( hOVCu238M2,        fParameters[60] );
  fModelTotuM2->Add( hSIu238M2,         fParameters[63] );

  fModelTotkM2->Add( hTeO2k40M2,        fParameters[3]  );
  fModelTotkM2->Add( hCuFramek40M2,     fParameters[17] );
  fModelTotkM2->Add( hCuBoxk40M2,       fParameters[25] );
  fModelTotkM2->Add( h50mKk40M2,        fParameters[33] );
  fModelTotkM2->Add( h600mKk40M2,       fParameters[39] );
  fModelTotkM2->Add( hPbRomk40M2,       fParameters[45] );
  fModelTotkM2->Add( hMBk40M2,          fParameters[50] );
  fModelTotkM2->Add( hIVCk40M2,         fParameters[54] );
  fModelTotkM2->Add( hOVCk40M2,         fParameters[58] ); 
  fModelTotkM2->Add( hSIk40M2,          fParameters[61] );

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
  fModelTotpbM2->Add( hCuFramepb210M2,  fParameters[19] );
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

  fModelTotSpbM2->Add( hTeO2Spb210M2_01,      fParameters[64] );
  fModelTotSpbM2->Add( hTeO2Sxpb210M2_001,    fParameters[69] );
  fModelTotSpbM2->Add( hTeO2Sxpb210M2_01,     fParameters[70] );
  fModelTotSpbM2->Add( hTeO2Sxpb210M2_1,      fParameters[71] );
  fModelTotSpbM2->Add( hTeO2Sxpb210M2_10,     fParameters[72] );
  fModelTotSpbM2->Add( hCuFrameSxpb210M2_001, fParameters[86] );
  fModelTotSpbM2->Add( hCuFrameSxpb210M2_01,  fParameters[87] );
  fModelTotSpbM2->Add( hCuFrameSxpb210M2_1,   fParameters[88] );
  fModelTotSpbM2->Add( hCuFrameSxpb210M2_10,  fParameters[89] );
  fModelTotSpbM2->Add( hCuBoxSxpb210M2_001,   fParameters[100] );
  fModelTotSpbM2->Add( hCuBoxSxpb210M2_01,    fParameters[101] );
  fModelTotSpbM2->Add( hCuBoxSxpb210M2_1,     fParameters[102] );
  fModelTotSpbM2->Add( hCuBoxSxpb210M2_10,    fParameters[103] );

  fModelTotSpoM2->Add( hTeO2Spo210M2_001,     fParameters[65] );
  fModelTotSpoM2->Add( hTeO2Spo210M2_01,      fParameters[66] );
  fModelTotSpoM2->Add( hTeO2Sxpo210M2_001,    fParameters[73] );
  fModelTotSpoM2->Add( hTeO2Sxpo210M2_01,     fParameters[74] );
  fModelTotSpoM2->Add( hTeO2Sxpo210M2_1,      fParameters[75] );

  fModelTotSthM2->Add( hTeO2Sth232M2_01,      fParameters[67] );
  fModelTotSthM2->Add( hTeO2Sxth232M2_001,    fParameters[76] );
  fModelTotSthM2->Add( hTeO2Sxth232M2_01,     fParameters[77] );
  fModelTotSthM2->Add( hTeO2Sxth232M2_1,      fParameters[78] );
  fModelTotSthM2->Add( hTeO2Sxth232M2_10,     fParameters[79] );
  fModelTotSthM2->Add( hCuFrameSth232M2_1,    fParameters[84] );
  fModelTotSthM2->Add( hCuFrameSxth232M2_001, fParameters[90] );
  fModelTotSthM2->Add( hCuFrameSxth232M2_01,  fParameters[91] );
  fModelTotSthM2->Add( hCuFrameSxth232M2_1,   fParameters[92] );
  fModelTotSthM2->Add( hCuFrameSxth232M2_10,  fParameters[93] );
  fModelTotSthM2->Add( hCuBoxSth232M2_1,      fParameters[98] );
  fModelTotSthM2->Add( hCuBoxSxth232M2_001,   fParameters[104] );
  fModelTotSthM2->Add( hCuBoxSxth232M2_01,    fParameters[105] );
  fModelTotSthM2->Add( hCuBoxSxth232M2_1,     fParameters[106] );
  fModelTotSthM2->Add( hCuBoxSxth232M2_10,    fParameters[107] );

  fModelTotSuM2->Add( hTeO2Su238M2_01,        fParameters[68] );
  fModelTotSuM2->Add( hTeO2Sxu238M2_001,      fParameters[80] );
  fModelTotSuM2->Add( hTeO2Sxu238M2_01,       fParameters[81] );
  fModelTotSuM2->Add( hTeO2Sxu238M2_1,        fParameters[82] );
  fModelTotSuM2->Add( hTeO2Sxu238M2_10,       fParameters[83] );
  fModelTotSuM2->Add( hCuFrameSu238M2_1,      fParameters[85] );
  fModelTotSuM2->Add( hCuFrameSxu238M2_001,   fParameters[94] );
  fModelTotSuM2->Add( hCuFrameSxu238M2_01,    fParameters[95] );
  fModelTotSuM2->Add( hCuFrameSxu238M2_1,     fParameters[96] );
  fModelTotSuM2->Add( hCuFrameSxu238M2_10,    fParameters[97] );
  fModelTotSuM2->Add( hCuBoxSu238M2_1,        fParameters[99] );
  fModelTotSuM2->Add( hCuBoxSxu238M2_001,     fParameters[108] );
  fModelTotSuM2->Add( hCuBoxSxu238M2_01,      fParameters[109] );
  fModelTotSuM2->Add( hCuBoxSxu238M2_1,       fParameters[110] );
  fModelTotSuM2->Add( hCuBoxSxu238M2_10,      fParameters[111] );

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
  // fDataHistoM1->GetXaxis()->SetRange(1, 2650/dBinSize+5);
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
  // legfit1->AddEntry(fModelTotNDBDM1, "NDBD", "l");
  // legfit1->AddEntry(fModelTot2NDBDM1, "2NDBD", "l");
  // legfit1->AddEntry(fModelTotbiM1, "bi-207", "l");
  // legfit1->AddEntry(fModelTotmnM1, "mn-54", "l");
  // legfit1->AddEntry(fModelTotpbM1 , "pb-210", "l");    
  legfit1->Draw();





  TCanvas *c2 = new TCanvas("c2", "c2", 1200, 800);
  c2->SetLogy();

  ///// Draw Data M2
  fDataHistoM2->SetLineColor(1);
  fDataHistoM2->SetLineWidth(2);
  fDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
  fDataHistoM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)/yr", dBinSize));
  fDataHistoM2->SetMaximum(9000);
  // fDataHistoM2->GetXaxis()->SetRange(1/dBinSize-5, 2650/dBinSize+5);
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
  // legfit2->AddEntry(fModelTotNDBDM2, "NDBD", "l");
  // legfit2->AddEntry(fModelTot2NDBDM2, "2NDBD", "l");
  // legfit2->AddEntry(fModelTotbiM2, "bi-207", "l");
  // legfit2->AddEntry(fModelTotmnM2, "mn-54", "l");
  // legfit2->AddEntry(fModelTotpbM2 , "bi-210", "l");
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

  TCanvas *cBkgAdap1 = new TCanvas("cBkgAdap1", "cBkgAdap1", 1200, 800);
  cBkgAdap1->SetLogy();
  fAdapDataHistoM1->SetLineColor(1);
  fAdapDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)/yr", dBinSize));
  fAdapDataHistoM1->Draw();

  TCanvas *cBkgAdap2 = new TCanvas("cBkgAdap2", "cBkgAdap2", 1200, 800);
  cBkgAdap2->SetLogy();
  fAdapDataHistoM2->SetLineColor(1);
  fAdapDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM2->GetYaxis()->SetTitle(Form("Counts/(%d keV)/yr", dBinSize));
  fAdapDataHistoM2->Draw();

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

    datam1_i = fAdapDataHistoM1->GetBinContent(i)*fAdapDataHistoM1->GetBinWidth(i)/2.0; // For real data

    modelm1_i = fModelTotAdapM1->GetBinContent(i)*fAdapDataHistoM1->GetBinWidth(i)/2.0;

    if(modelm1_i != 0 && datam1_i != 0)
    {
      chiSquare += 2 * (modelm1_i - datam1_i + datam1_i * TMath::Log(datam1_i/modelm1_i));
    }
  }

  for(int i = dFitMinBinM2; i <= dFitMaxBinM2; i++)
  {
    datam2_i = fAdapDataHistoM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i)/2.0; // For real data

    modelm2_i = fModelTotAdapM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i)/2.0;

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
  fFile2 = new TFile("/Users/brian/macros/Simulations/Production/MCProduction_Surface_nonnormalized.root");

  // fFile = new TFile("/Users/brian/macros/Simulations/Production/MCProduction_Bulk.root"); 
  // fFile2 = new TFile("/Users/brian/macros/Simulations/Production/MCProduction_Surface.root");

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
  qtree->Add("/Users/brian/macros/Simulations/Toy/combi1/combi1.root"); 
  qtree->Project("fDataHistoTot", "Ener2");
  qtree->Project("fDataHistoM1",  "Ener2", "Multiplicity == 1");
  qtree->Project("fDataHistoM2",  "Ener2", "Multiplicity == 2");
*/
  qtree->Add("/Users/brian/macros/CUOREZ/Bkg/Q0_DR2_BackgroundSignalData.root"); 
  qtree->Project("fDataHistoTot", "Energy");
  qtree->Project("fDataHistoM1",  "Energy", "Multiplicity == 1");
  qtree->Project("fDataHistoM2",  "Energy", "Multiplicity == 2");

  cout << "Loaded Data" << endl;
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
  cout<< "Par26 = "  << fParameters[26] << " +/- " << fParError[26] << endl;
  cout<< "Par27 = "  << fParameters[27] << " +/- " << fParError[27] << endl;
  cout<< "Par28 = "  << fParameters[28] << " +/- " << fParError[28] << endl;
  cout<< "Par29 = "  << fParameters[29] << " +/- " << fParError[29] << endl;  
  cout<< "Par30 = "  << fParameters[30] << " +/- " << fParError[30] << endl;
  cout<< "Par31 = "  << fParameters[31] << " +/- " << fParError[31] << endl;
  cout<< "Par32 = "  << fParameters[32] << " +/- " << fParError[32] << endl;
  cout<< "Par33 = "  << fParameters[33] << " +/- " << fParError[33] << endl;
  cout<< "Par34 = "  << fParameters[34] << " +/- " << fParError[34] << endl;
  cout<< "Par35 = "  << fParameters[35] << " +/- " << fParError[35] << endl;
  cout<< "Par36 = "  << fParameters[36] << " +/- " << fParError[36] << endl;
  cout<< "Par37 = "  << fParameters[37] << " +/- " << fParError[37] << endl;
  cout<< "Par38 = "  << fParameters[38] << " +/- " << fParError[38] << endl;
  cout<< "Par39 = "  << fParameters[39] << " +/- " << fParError[39] << endl;  
  cout<< "Par40 = "  << fParameters[40] << " +/- " << fParError[40] << endl;
  cout<< "Par41 = "  << fParameters[41] << " +/- " << fParError[41] << endl;
  cout<< "Par42 = "  << fParameters[42] << " +/- " << fParError[42] << endl;
  cout<< "Par43 = "  << fParameters[43] << " +/- " << fParError[43] << endl;
  cout<< "Par44 = "  << fParameters[44] << " +/- " << fParError[44] << endl;
  cout<< "Par45 = "  << fParameters[45] << " +/- " << fParError[45] << endl;
  cout<< "Par46 = "  << fParameters[46] << " +/- " << fParError[46] << endl;
  cout<< "Par47 = "  << fParameters[47] << " +/- " << fParError[47] << endl;
  cout<< "Par48 = "  << fParameters[48] << " +/- " << fParError[48] << endl;
  cout<< "Par49 = "  << fParameters[49] << " +/- " << fParError[49] << endl;
  cout<< "Par50 = "  << fParameters[50] << " +/- " << fParError[50] << endl;
  cout<< "Par51 = "  << fParameters[51] << " +/- " << fParError[51] << endl;
  cout<< "Par52 = "  << fParameters[52] << " +/- " << fParError[52] << endl;
  cout<< "Par53 = "  << fParameters[53] << " +/- " << fParError[53] << endl;
  cout<< "Par54 = "  << fParameters[54] << " +/- " << fParError[54] << endl;
  cout<< "Par55 = "  << fParameters[55] << " +/- " << fParError[55] << endl;
  cout<< "Par56 = "  << fParameters[56] << " +/- " << fParError[56] << endl;
  cout<< "Par57 = "  << fParameters[57] << " +/- " << fParError[57] << endl;
  cout<< "Par58 = "  << fParameters[58] << " +/- " << fParError[58] << endl;
  cout<< "Par59 = "  << fParameters[59] << " +/- " << fParError[59] << endl;
  cout<< "Par60 = "  << fParameters[60] << " +/- " << fParError[60] << endl;
  cout<< "Par61 = "  << fParameters[61] << " +/- " << fParError[61] << endl;
  cout<< "Par62 = "  << fParameters[62] << " +/- " << fParError[62] << endl;
  cout<< "Par63 = "  << fParameters[63] << " +/- " << fParError[63] << endl;
  cout<< "Par64 = "  << fParameters[64] << " +/- " << fParError[64] << endl;
  cout<< "Par65 = "  << fParameters[65] << " +/- " << fParError[65] << endl;
  cout<< "Par66 = "  << fParameters[66] << " +/- " << fParError[66] << endl;
  cout<< "Par67 = "  << fParameters[67] << " +/- " << fParError[67] << endl;
  cout<< "Par68 = "  << fParameters[68] << " +/- " << fParError[68] << endl;
  cout<< "Par69 = "  << fParameters[69] << " +/- " << fParError[69] << endl;
  cout<< "Par70 = "  << fParameters[70] << " +/- " << fParError[70] << endl;
  cout<< "Par71 = "  << fParameters[71] << " +/- " << fParError[71] << endl;
  cout<< "Par72 = "  << fParameters[72] << " +/- " << fParError[72] << endl;
  cout<< "Par73 = "  << fParameters[73] << " +/- " << fParError[73] << endl;
  cout<< "Par74 = "  << fParameters[74] << " +/- " << fParError[74] << endl;
  cout<< "Par75 = "  << fParameters[75] << " +/- " << fParError[75] << endl;
  cout<< "Par76 = "  << fParameters[76] << " +/- " << fParError[76] << endl;
  cout<< "Par77 = "  << fParameters[77] << " +/- " << fParError[77] << endl;
  cout<< "Par78 = "  << fParameters[78] << " +/- " << fParError[78] << endl;
  cout<< "Par79 = "  << fParameters[79] << " +/- " << fParError[79] << endl;
  cout<< "Par80 = "  << fParameters[80] << " +/- " << fParError[80] << endl;
  cout<< "Par81 = "  << fParameters[81] << " +/- " << fParError[81] << endl;
  cout<< "Par82 = "  << fParameters[82] << " +/- " << fParError[82] << endl;
  cout<< "Par83 = "  << fParameters[83] << " +/- " << fParError[83] << endl;
  cout<< "Par84 = "  << fParameters[84] << " +/- " << fParError[84] << endl;
  cout<< "Par85 = "  << fParameters[85] << " +/- " << fParError[85] << endl;
  cout<< "Par86 = "  << fParameters[86] << " +/- " << fParError[86] << endl;
  cout<< "Par87 = "  << fParameters[87] << " +/- " << fParError[87] << endl;
  cout<< "Par88 = "  << fParameters[88] << " +/- " << fParError[88] << endl;
  cout<< "Par89 = "  << fParameters[89] << " +/- " << fParError[89] << endl;
  cout<< "Par90 = "  << fParameters[90] << " +/- " << fParError[90] << endl;
  cout<< "Par91 = "  << fParameters[91] << " +/- " << fParError[91] << endl;
  cout<< "Par92 = "  << fParameters[92] << " +/- " << fParError[92] << endl;
  cout<< "Par93 = "  << fParameters[93] << " +/- " << fParError[93] << endl;
  cout<< "Par94 = "  << fParameters[94] << " +/- " << fParError[94] << endl;
  cout<< "Par95 = "  << fParameters[95] << " +/- " << fParError[95] << endl;
  cout<< "Par96 = "  << fParameters[96] << " +/- " << fParError[96] << endl;
  cout<< "Par97 = "  << fParameters[97] << " +/- " << fParError[97] << endl;
  cout<< "Par98 = "  << fParameters[98] << " +/- " << fParError[98] << endl;
  cout<< "Par99 = "  << fParameters[99] << " +/- " << fParError[99] << endl;
  cout<< "Par100 = "  << fParameters[100] << " +/- " << fParError[100] << endl;
  cout<< "Par101 = "  << fParameters[101] << " +/- " << fParError[101] << endl;
  cout<< "Par102 = "  << fParameters[102] << " +/- " << fParError[102] << endl;
  cout<< "Par103 = "  << fParameters[103] << " +/- " << fParError[103] << endl;
  cout<< "Par104 = "  << fParameters[104] << " +/- " << fParError[104] << endl;
  cout<< "Par105 = "  << fParameters[105] << " +/- " << fParError[105] << endl;
  cout<< "Par106 = "  << fParameters[106] << " +/- " << fParError[106] << endl;
  cout<< "Par107 = "  << fParameters[107] << " +/- " << fParError[107] << endl;
  cout<< "Par108 = "  << fParameters[108] << " +/- " << fParError[108] << endl;
}


// Set Parameters in Model
void TBackgroundModel::SetParameters(int index, double value)
{
	// Change the index max depending on model
	if(index > 112) cout << "Index too large" << endl;
	else fParameters[index] = value;
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
  if(dNumCalls%100==0)
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

  fModelTotM1->Add( hSIk40M1,        fParameters[61]);
  fModelTotM1->Add( hSIth232M1,      fParameters[62]);
  fModelTotM1->Add( hSIu238M1,       fParameters[63]);

/*
  fModelTotM1->Add( hTeO2Spb210M1_01,      fParameters[64]);
  fModelTotM1->Add( hTeO2Spo210M1_001,     fParameters[65]);
  fModelTotM1->Add( hTeO2Spo210M1_01,      fParameters[66]);
  fModelTotM1->Add( hTeO2Sth232M1_01,      fParameters[67]);
  fModelTotM1->Add( hTeO2Su238M1_01,       fParameters[68]);
  fModelTotM1->Add( hTeO2Sxpb210M1_001,    fParameters[69]);
  fModelTotM1->Add( hTeO2Sxpb210M1_01,     fParameters[70]);
  fModelTotM1->Add( hTeO2Sxpb210M1_1,      fParameters[71]);
  fModelTotM1->Add( hTeO2Sxpb210M1_10,     fParameters[72]);
  fModelTotM1->Add( hTeO2Sxpo210M1_001,    fParameters[73]);
  fModelTotM1->Add( hTeO2Sxpo210M1_01,     fParameters[74]);
  fModelTotM1->Add( hTeO2Sxpo210M1_1,      fParameters[75]);
  fModelTotM1->Add( hTeO2Sxth232M1_001,    fParameters[76]);
  fModelTotM1->Add( hTeO2Sxth232M1_01,     fParameters[77]);
  fModelTotM1->Add( hTeO2Sxth232M1_1,      fParameters[78]);
  fModelTotM1->Add( hTeO2Sxth232M1_10,     fParameters[79]);
  fModelTotM1->Add( hTeO2Sxu238M1_001,     fParameters[80]);
  fModelTotM1->Add( hTeO2Sxu238M1_01,      fParameters[81]);
  fModelTotM1->Add( hTeO2Sxu238M1_1,       fParameters[82]);
  fModelTotM1->Add( hTeO2Sxu238M1_10,      fParameters[83]);

  fModelTotM1->Add( hCuFrameSth232M1_1,    fParameters[84]);
  fModelTotM1->Add( hCuFrameSu238M1_1,     fParameters[85]);
  fModelTotM1->Add( hCuFrameSxpb210M1_001, fParameters[86]);
  fModelTotM1->Add( hCuFrameSxpb210M1_01,  fParameters[87]);
  fModelTotM1->Add( hCuFrameSxpb210M1_1,   fParameters[88]);
  fModelTotM1->Add( hCuFrameSxpb210M1_10,  fParameters[89]);
  fModelTotM1->Add( hCuFrameSxth232M1_001, fParameters[90]);
  fModelTotM1->Add( hCuFrameSxth232M1_01,  fParameters[91]);
  fModelTotM1->Add( hCuFrameSxth232M1_1,   fParameters[92]);
  fModelTotM1->Add( hCuFrameSxth232M1_10,  fParameters[93]);
  fModelTotM1->Add( hCuFrameSxu238M1_001,  fParameters[94]);
  fModelTotM1->Add( hCuFrameSxu238M1_01,   fParameters[95]);
  fModelTotM1->Add( hCuFrameSxu238M1_1,    fParameters[96]);
  fModelTotM1->Add( hCuFrameSxu238M1_10,   fParameters[97]);

  fModelTotM1->Add( hCuBoxSth232M1_1,      fParameters[98]);
  fModelTotM1->Add( hCuBoxSu238M1_1,       fParameters[99]);
  fModelTotM1->Add( hCuBoxSxpb210M1_001,   fParameters[100]);
  fModelTotM1->Add( hCuBoxSxpb210M1_01,    fParameters[101]);
  fModelTotM1->Add( hCuBoxSxpb210M1_1,     fParameters[102]);
  fModelTotM1->Add( hCuBoxSxpb210M1_10,    fParameters[103]);
  fModelTotM1->Add( hCuBoxSxth232M1_001,   fParameters[104]);
  fModelTotM1->Add( hCuBoxSxth232M1_01,    fParameters[105]);
  fModelTotM1->Add( hCuBoxSxth232M1_1,     fParameters[106]);
  fModelTotM1->Add( hCuBoxSxth232M1_10,    fParameters[107]);
  fModelTotM1->Add( hCuBoxSxu238M1_001,    fParameters[108]);
  fModelTotM1->Add( hCuBoxSxu238M1_01,     fParameters[109]);
  fModelTotM1->Add( hCuBoxSxu238M1_1,      fParameters[110]);
  fModelTotM1->Add( hCuBoxSxu238M1_10,     fParameters[111]);    
*/
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

  fModelTotM2->Add( hSIk40M2,        fParameters[61]);
  fModelTotM2->Add( hSIth232M2,      fParameters[62]);
  fModelTotM2->Add( hSIu238M2,       fParameters[63]);  
/*
  fModelTotM2->Add( hTeO2Spb210M2_01,      fParameters[64]);
  fModelTotM2->Add( hTeO2Spo210M2_001,     fParameters[65]);
  fModelTotM2->Add( hTeO2Spo210M2_01,      fParameters[66]);
  fModelTotM2->Add( hTeO2Sth232M2_01,      fParameters[67]);
  fModelTotM2->Add( hTeO2Su238M2_01,       fParameters[68]);
  fModelTotM2->Add( hTeO2Sxpb210M2_001,    fParameters[69]);
  fModelTotM2->Add( hTeO2Sxpb210M2_01,     fParameters[70]);
  fModelTotM2->Add( hTeO2Sxpb210M2_1,      fParameters[71]);
  fModelTotM2->Add( hTeO2Sxpb210M2_10,     fParameters[72]);
  fModelTotM2->Add( hTeO2Sxpo210M2_001,    fParameters[73]);
  fModelTotM2->Add( hTeO2Sxpo210M2_01,     fParameters[74]);
  fModelTotM2->Add( hTeO2Sxpo210M2_1,      fParameters[75]);
  fModelTotM2->Add( hTeO2Sxth232M2_001,    fParameters[76]);
  fModelTotM2->Add( hTeO2Sxth232M2_01,     fParameters[77]);
  fModelTotM2->Add( hTeO2Sxth232M2_1,      fParameters[78]);
  fModelTotM2->Add( hTeO2Sxth232M2_10,     fParameters[79]);
  fModelTotM2->Add( hTeO2Sxu238M2_001,     fParameters[80]);
  fModelTotM2->Add( hTeO2Sxu238M2_01,      fParameters[81]);
  fModelTotM2->Add( hTeO2Sxu238M2_1,       fParameters[82]);
  fModelTotM2->Add( hTeO2Sxu238M2_10,      fParameters[83]);

  fModelTotM2->Add( hCuFrameSth232M2_1,    fParameters[84]);
  fModelTotM2->Add( hCuFrameSu238M2_1,     fParameters[85]);
  fModelTotM2->Add( hCuFrameSxpb210M2_001, fParameters[86]);
  fModelTotM2->Add( hCuFrameSxpb210M2_01,  fParameters[87]);
  fModelTotM2->Add( hCuFrameSxpb210M2_1,   fParameters[88]);
  fModelTotM2->Add( hCuFrameSxpb210M2_10,  fParameters[89]);
  fModelTotM2->Add( hCuFrameSxth232M2_001, fParameters[90]);
  fModelTotM2->Add( hCuFrameSxth232M2_01,  fParameters[91]);
  fModelTotM2->Add( hCuFrameSxth232M2_1,   fParameters[92]);
  fModelTotM2->Add( hCuFrameSxth232M2_10,  fParameters[93]);
  fModelTotM2->Add( hCuFrameSxu238M2_001,  fParameters[94]);
  fModelTotM2->Add( hCuFrameSxu238M2_01,   fParameters[95]);
  fModelTotM2->Add( hCuFrameSxu238M2_1,    fParameters[96]);
  fModelTotM2->Add( hCuFrameSxu238M2_10,   fParameters[97]);

  fModelTotM2->Add( hCuBoxSth232M2_1,      fParameters[98]);
  fModelTotM2->Add( hCuBoxSu238M2_1,       fParameters[99]);
  fModelTotM2->Add( hCuBoxSxpb210M2_001,   fParameters[100]);
  fModelTotM2->Add( hCuBoxSxpb210M2_01,    fParameters[101]);
  fModelTotM2->Add( hCuBoxSxpb210M2_1,     fParameters[102]);
  fModelTotM2->Add( hCuBoxSxpb210M2_10,    fParameters[103]);
  fModelTotM2->Add( hCuBoxSxth232M2_001,   fParameters[104]);
  fModelTotM2->Add( hCuBoxSxth232M2_01,    fParameters[105]);
  fModelTotM2->Add( hCuBoxSxth232M2_1,     fParameters[106]);
  fModelTotM2->Add( hCuBoxSxth232M2_10,    fParameters[107]);
  fModelTotM2->Add( hCuBoxSxu238M2_001,    fParameters[108]);
  fModelTotM2->Add( hCuBoxSxu238M2_01,     fParameters[109]);
  fModelTotM2->Add( hCuBoxSxu238M2_1,      fParameters[110]);
  fModelTotM2->Add( hCuBoxSxu238M2_10,     fParameters[111]);     
*/
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
  if(dNumCalls%100==0)
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


   TMinuit minuit(112);

   // Reduce Minuit Output
   minuit.SetPrintLevel(1);
//   minuit.Command("SET MINImize 1000 0.001");
   minuit.Command("SET STRategy 2");
  //minuit.Command("SET IMProve 1000 ");
   minuit.SetMaxIterations(10000);
   minuit.SetObjectFit(this); //see the external FCN  above
   
   //define the parameters and set the ranges and initial guesses see ROOTs TMinuit documentation

   // Close Th and Close Ra now split into its various sections, far Th and Ra still the same
   // This step after previous step full fit converges, just to see if any differences show up
   ////////////////////////////////////////////////
   // Using more parameters
   ////////////////////////////////////////////////
   minuit.DefineParameter(0, "TeO2 0nu",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(1, "TeO2 2nu",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(2, "TeO2 co60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(3, "TeO2 k40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(4, "TeO2 pb210",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(5, "TeO2 po210",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(6, "TeO2 te125",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(7, "TeO2 th232",  0.1, 0.1, 0., 1000000);
   minuit.DefineParameter(8, "TeO2 th228",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(9, "TeO2 ra226",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(10, "TeO2 rn222",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(11, "TeO2 u238",  0.1, 0.1, 0., 1000000);
   minuit.DefineParameter(12, "TeO2 th230",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(13, "TeO2 u234",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(14, "CuFrame co58",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(15, "CuFrame co60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(16, "CuFrame cs137",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(17, "CuFrame k40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(18, "CuFrame mn54",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(19, "CuFrame pb210",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(20, "CuFrame th232",  0.1, 0.1, 0., 1000000);
   minuit.DefineParameter(21, "CuFrame u238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(22, "CuBoxco58",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(23, "CuBoxco60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(24, "CuBoxcs137",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(25, "CuBoxk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(26, "CuBoxmn54",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(27, "CuBoxpb210",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(28, "CuBoxth232",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(29, "CuBoxu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(30, "50mKco58",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(31, "50mKco60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(32, "50mKcs137",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(33, "50mKk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(34, "50mKmn54",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(35, "50mKpb210",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(36, "50mKth232",  0.1, 0.1, 0., 1000000);
   minuit.DefineParameter(37, "50mKu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(38, "600mKco60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(39, "600mKk40",  0., 0.1, 0., 1000000); 
   minuit.DefineParameter(40, "600mKth232",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(41, "600mKu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(42, "PbRombi207",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(43, "PbRomco60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(44, "PbRomcs137",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(45, "PbRomk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(46, "PbRompb210",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(47, "PbRomth232",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(48, "PbRomu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(49, "MBco60",  0., 0.1, 0., 1000000); 
   minuit.DefineParameter(50, "MBk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(51, "MBth232",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(52, "MBu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(53, "IVCco60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(54, "IVCk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(55, "IVCth232",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(56, "IVCu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(57, "OVCco60",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(58, "OVCk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(59, "OVCth232",  0., 0.1, 0., 1000000);    
   minuit.DefineParameter(60, "OVCu238",  0., 0.1, 0., 1000000);

   minuit.DefineParameter(61, "SIk40",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(62, "SIth232",  0., 0.1, 0., 1000000);    
   minuit.DefineParameter(63, "SIu238",  0., 0.1, 0., 1000000);

/*
   minuit.DefineParameter(64, "TeO2Spb210_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(65, "TeO2Spo210_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(66, "TeO2Spo210_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(67, "TeO2Sth232_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(68, "TeO2Su238_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(69, "TeO2Sxpb210_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(70, "TeO2Sxpb210_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(71, "TeO2Sxpb210_1",  0.03, 0.1, 0., 1000000);
   minuit.DefineParameter(72, "TeO2Sxpb210_10",  0., 0.1, 0., 1000000);    
   minuit.DefineParameter(73, "TeO2Sxpo210_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(74, "TeO2Sxpo210_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(75, "TeO2Sxpo210_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(76, "TeO2Sxth232_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(77, "TeO2Sxth232_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(78, "TeO2Sxth232_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(79, "TeO2Sxth232_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(80, "TeO2Sxu238_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(81, "TeO2Sxu238_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(82, "TeO2Sxu238_1",  0., 0.1, 0., 1000000);   
   minuit.DefineParameter(83, "TeO2Sxu238_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(84, "CuFrameSth232_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(85, "CuFrameSu238_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(86, "CuFrameSxpb210_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(87, "CuFrameSxpb210_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(88, "CuFrameSxpb210_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(89, "CuFrameSxpb210_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(90, "CuFrameSxth232_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(91, "CuFrameSxth232_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(92, "CuFrameSxth232_1",  0., 0.1, 0., 1000000);   
   minuit.DefineParameter(93, "CuFrameSxth232_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(94, "CuFrameSxu238_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(95, "CuFrameSxu238_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(96, "CuFrameSxu238_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(97, "CuFrameSxu238_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(98, "CuBoxSth232_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(99, "CuBoxSu238_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(100, "CuBoxSxpb210_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(101, "CuBoxSxpb210_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(102, "CuBoxSxpb210_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(103, "CuBoxSxpb210_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(104, "CuBoxSxth232_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(105, "CuBoxSxth232_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(106, "CuBoxSxth232_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(107, "CuBoxSxth232_10",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(108, "CuBoxSxu238_001",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(109, "CuBoxSxu238_01",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(110, "CuBoxSxu238_1",  0., 0.1, 0., 1000000);
   minuit.DefineParameter(111, "CuBoxSxu238_10",  0., 0.1, 0., 1000000);
*/

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
   minuit.FixParameter(61); // SI k40
   minuit.FixParameter(62); // SI th232
   minuit.FixParameter(63); // SI u238   
   /*
   minuit.FixParameter(64); // TeO2 S pb210 01
   minuit.FixParameter(65); // TeO2 S po210 001
   minuit.FixParameter(66); // TeO2 S po210 01
   minuit.FixParameter(67); // TeO2 S th232 01
   minuit.FixParameter(68); // TeO2 S u238 01
   minuit.FixParameter(69); // TeO2 Sx pb210 001
   minuit.FixParameter(70); // TeO2 Sx pb210 01
   // minuit.FixParameter(71); // TeO2 Sx pb210 1
   minuit.FixParameter(72); // TeO2 S pb210 10
   minuit.FixParameter(73); // TeO2 Sx po210 001
   minuit.FixParameter(74); // TeO2 Sx po210 01
   minuit.FixParameter(75); // TeO2 Sx po210 1
   minuit.FixParameter(76); // TeO2 Sx th232 001
   minuit.FixParameter(77); // TeO2 Sx th232 01
   minuit.FixParameter(78); // TeO2 Sx th232 1
   minuit.FixParameter(79); // TeO2 Sx th232 10
   minuit.FixParameter(80); // TeO2 Sx u238 001
   minuit.FixParameter(81); // TeO2 Sx u238 01
   minuit.FixParameter(82); // TeO2 Sx u238 1
   minuit.FixParameter(83); // TeO2 Sx u238 10
   minuit.FixParameter(84); // Frame S th232 1
   minuit.FixParameter(85); // Frame S u238 1
   minuit.FixParameter(86); // Frame Sx pb210 001
   minuit.FixParameter(87); // Frame Sx pb210 01
   minuit.FixParameter(88); // Frame Sx pb210 1
   minuit.FixParameter(89); // Frame Sx pb210 10
   minuit.FixParameter(90); // Frame Sx th232 001
   minuit.FixParameter(91); // Frame Sx th232 01
   minuit.FixParameter(92); // Frame Sx th232 1
   minuit.FixParameter(93); // Frame Sx th232 10
   minuit.FixParameter(94); // Frame Sx u238 001
   minuit.FixParameter(95); // Frame Sx u238 01
   minuit.FixParameter(96); // Frame Sx u238 1
   minuit.FixParameter(97); // Frame Sx u238 10
   minuit.FixParameter(98); // CuBox S th232 1
   minuit.FixParameter(99); // CuBox S u238 1
   minuit.FixParameter(100); // CuBox Sx pb210 001
   minuit.FixParameter(101); // CuBox Sx pb210 01
   minuit.FixParameter(102); // CuBox Sx pb210 1
   minuit.FixParameter(103); // CuBox Sx pb210 10 
   minuit.FixParameter(104); // CuBox Sx th232 001
   minuit.FixParameter(105); // CuBox Sx th232 01
   minuit.FixParameter(106); // CuBox Sx th232 1
   minuit.FixParameter(107); // CuBox Sx th232 10
   minuit.FixParameter(108); // CuBox Sx u238 001
   minuit.FixParameter(109); // CuBox Sx u238 01
   minuit.FixParameter(110); // CuBox Sx u238 1
   minuit.FixParameter(111); // CuBox Sx u238 10
   */
   // Number of Parameters (for Chi-squared/NDF calculation)
   int dNumParameters = 4;
   //Tell minuit what external function to use 
   minuit.SetFCN(myExternal_FCNAdap);
   
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
  minuit.GetParameter(61,  fParameters[61],   fParError[61]);
  minuit.GetParameter(62,  fParameters[62],   fParError[62]);
  minuit.GetParameter(63,  fParameters[63],   fParError[63]);  
/*
  minuit.GetParameter(64,  fParameters[64],   fParError[64]);
  minuit.GetParameter(65,  fParameters[65],   fParError[65]);
  minuit.GetParameter(66,  fParameters[66],   fParError[66]);
  minuit.GetParameter(67,  fParameters[67],   fParError[67]);
  minuit.GetParameter(68,  fParameters[68],   fParError[68]);
  minuit.GetParameter(69,  fParameters[69],   fParError[69]);
  minuit.GetParameter(70,  fParameters[70],   fParError[70]);
  minuit.GetParameter(71,  fParameters[71],   fParError[71]);
  minuit.GetParameter(72,  fParameters[72],   fParError[72]);
  minuit.GetParameter(73,  fParameters[73],   fParError[73]);
  minuit.GetParameter(74,  fParameters[74],   fParError[74]);
  minuit.GetParameter(75,  fParameters[75],   fParError[75]);
  minuit.GetParameter(76,  fParameters[76],   fParError[76]);
  minuit.GetParameter(77,  fParameters[77],   fParError[77]);
  minuit.GetParameter(78,  fParameters[78],   fParError[78]);
  minuit.GetParameter(79,  fParameters[79],   fParError[79]);
  minuit.GetParameter(80,  fParameters[80],   fParError[80]);
  minuit.GetParameter(81,  fParameters[81],   fParError[81]);
  minuit.GetParameter(82,  fParameters[82],   fParError[82]);
  minuit.GetParameter(83,  fParameters[83],   fParError[83]);
  minuit.GetParameter(84,  fParameters[84],   fParError[84]);
  minuit.GetParameter(85,  fParameters[85],   fParError[85]);
  minuit.GetParameter(86,  fParameters[86],   fParError[86]);
  minuit.GetParameter(87,  fParameters[87],   fParError[87]);
  minuit.GetParameter(88,  fParameters[88],   fParError[88]);
  minuit.GetParameter(89,  fParameters[89],   fParError[89]);
  minuit.GetParameter(90,  fParameters[90],   fParError[90]);
  minuit.GetParameter(91,  fParameters[91],   fParError[91]);
  minuit.GetParameter(92,  fParameters[92],   fParError[92]);
  minuit.GetParameter(93,  fParameters[93],   fParError[93]);
  minuit.GetParameter(94,  fParameters[94],   fParError[94]);
  minuit.GetParameter(95,  fParameters[95],   fParError[95]);
  minuit.GetParameter(96,  fParameters[96],   fParError[96]);
  minuit.GetParameter(97,  fParameters[97],   fParError[97]);
  minuit.GetParameter(98,  fParameters[98],   fParError[98]);
  minuit.GetParameter(99,  fParameters[99],   fParError[99]);
  minuit.GetParameter(100,  fParameters[100],   fParError[100]);
  minuit.GetParameter(101,  fParameters[101],   fParError[101]);
  minuit.GetParameter(102,  fParameters[102],   fParError[102]);
  minuit.GetParameter(103,  fParameters[103],   fParError[103]);
  minuit.GetParameter(104,  fParameters[104],   fParError[104]);
  minuit.GetParameter(105,  fParameters[105],   fParError[105]);
  minuit.GetParameter(106,  fParameters[106],   fParError[106]);
  minuit.GetParameter(107,  fParameters[107],   fParError[107]);
  minuit.GetParameter(108,  fParameters[108],   fParError[108]);
  minuit.GetParameter(109,  fParameters[109],   fParError[109]);
  minuit.GetParameter(110,  fParameters[110],   fParError[110]);
  minuit.GetParameter(111,  fParameters[111],   fParError[111]);
  */
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
  ////////// Only for testing
  // Correction for M2 spectra, it's the M1 spectra but scaled down by N_M1*1-Exp(R*T)
  // fTotCorrection->Add(fCorrectionM2, 180197*(1-TMath::Exp(-0.05*0.1)));


  TCanvas *cadap1 = new TCanvas("cadap1", "cadap1", 1200, 800);
  cadap1->SetLogy();

  ///// Draw Data M1
  fAdapDataHistoM1->SetLineColor(1);
  fAdapDataHistoM1->SetLineWidth(2);
  fAdapDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM1->GetYaxis()->SetTitle("Counts/yr");
  fAdapDataHistoM1->SetMaximum(90000);
  // fAdapDataHistoM1->GetXaxis()->SetRange(1, 2650/dBinSize+5);
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
  // legfit1->AddEntry(fModelTotAdapNDBDM1, "NDBD", "l");
  // legfit1->AddEntry(fModelTotAdap2NDBDM1, "2NDBD", "l");
  // legfit1->AddEntry(fModelTotAdapbiM1, "bi-207", "l");
  // legfit1->AddEntry(fModelTotAdapmnM1, "mn-54", "l");
  legfit1->Draw();





  TCanvas *cadap2 = new TCanvas("cadap2", "cadap2", 1200, 800);
  cadap2->SetLogy();

  ///// Draw Data M2
  fAdapDataHistoM2->SetLineColor(1);
  fAdapDataHistoM2->SetLineWidth(2);
  fAdapDataHistoM2->GetXaxis()->SetTitle("Energy (keV)");
  fAdapDataHistoM2->GetYaxis()->SetTitle("Counts/yr");
  fAdapDataHistoM2->SetMaximum(9000);
  // fAdapDataHistoM2->GetXaxis()->SetRange(1/dBinSize-5, 2650/dBinSize+5);
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
  // legfit2->AddEntry(fModelTotAdapNDBDM2, "NDBD", "l");
  // legfit2->AddEntry(fModelTotAdap2NDBDM2, "2NDBD", "l");
  // legfit2->AddEntry(fModelTotAdapbiM2, "bi-207", "l");
  // legfit2->AddEntry(fModelTotAdapmnM2, "mn-54", "l");

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

void myExternal_FCNAdap(int &n, double *grad, double &fval, double x[], int code)
{
  // Required External Wrapper for function to be minimized by Minuit 
 
  // This gets called for each value of the parameters minuit tries
  // here the x array contains the parameters you are trying to fit
  
  // here myClass should inherit from TObject
  TBackgroundModel* Obj = (TBackgroundModel*)gMinuit->GetObjectFit(); 

  // implement a method in your class for setting the parameters and thus update the parameters of your fitter class 
  Obj->SetParameters(0, x[0]);   
  Obj->SetParameters(1, x[1]);  
  Obj->SetParameters(2, x[2]);
  Obj->SetParameters(3, x[3]);
  Obj->SetParameters(4, x[4]);
  Obj->SetParameters(5, x[5]);   
  Obj->SetParameters(6, x[6]);  
  Obj->SetParameters(7, x[7]);
  Obj->SetParameters(8, x[8]);
  Obj->SetParameters(9, x[9]);
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
  Obj->SetParameters(61, x[61]);
  Obj->SetParameters(62, x[62]);
  Obj->SetParameters(63, x[63]);
  // Surface
  /*
  Obj->SetParameters(64, x[64]);
  Obj->SetParameters(65, x[65]);
  Obj->SetParameters(66, x[66]);
  Obj->SetParameters(67, x[67]);
  Obj->SetParameters(68, x[68]);
  Obj->SetParameters(69, x[69]);  
  Obj->SetParameters(70, x[70]);
  Obj->SetParameters(71, x[71]);
  Obj->SetParameters(72, x[72]);
  Obj->SetParameters(73, x[73]);
  Obj->SetParameters(74, x[74]);
  Obj->SetParameters(75, x[75]);
  Obj->SetParameters(76, x[76]);
  Obj->SetParameters(77, x[77]);
  Obj->SetParameters(78, x[78]);
  Obj->SetParameters(79, x[79]);
  Obj->SetParameters(80, x[80]);
  Obj->SetParameters(81, x[81]);
  Obj->SetParameters(82, x[82]);
  Obj->SetParameters(83, x[83]);
  Obj->SetParameters(84, x[84]);
  Obj->SetParameters(85, x[85]);
  Obj->SetParameters(86, x[86]);
  Obj->SetParameters(87, x[87]);
  Obj->SetParameters(88, x[88]);
  Obj->SetParameters(89, x[89]);
  Obj->SetParameters(90, x[90]);
  Obj->SetParameters(91, x[91]);
  Obj->SetParameters(92, x[92]);
  Obj->SetParameters(93, x[93]);
  Obj->SetParameters(94, x[94]);
  Obj->SetParameters(95, x[95]);
  Obj->SetParameters(96, x[96]);
  Obj->SetParameters(97, x[97]);
  Obj->SetParameters(98, x[98]);
  Obj->SetParameters(99, x[99]);
  Obj->SetParameters(100, x[100]);
  Obj->SetParameters(101, x[101]);
  Obj->SetParameters(102, x[102]);
  Obj->SetParameters(103, x[103]);
  Obj->SetParameters(104, x[104]);
  Obj->SetParameters(105, x[105]);
  Obj->SetParameters(106, x[106]);
  Obj->SetParameters(107, x[107]);
  Obj->SetParameters(108, x[108]);
  Obj->SetParameters(109, x[109]);
  Obj->SetParameters(110, x[110]);
  Obj->SetParameters(111, x[111]);
*/
  // Implement a method in your class that calculates the quantity you want to minimize, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimise fval
    Obj->UpdateModelAdaptive();
    fval = Obj->GetChiSquareAdaptive();
}



void TBackgroundModel::DrawMC()
{
// Draws all MC spectra, must Initialize first!

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

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

  hAdapTeO2th232M1->SetLineColor(1);
  hAdapCuFrameth232M1->SetLineColor(2);
  hAdapCuBoxth232M1->SetLineColor(3);
  hAdap50mKth232M1->SetLineColor(4);
  hAdap600mKth232M1->SetLineColor(5);
  hAdapIVCth232M1->SetLineColor(6);
  hAdapOVCth232M1->SetLineColor(7);
  hAdapPbRomth232M1->SetLineColor(8);
  hAdapMBth232M1->SetLineColor(9);
  hAdapSIth232M1->SetLineColor(11);


  hAdapTeO2th232M1->DrawNormalized();
  hAdapCuFrameth232M1->DrawNormalized("SAME");
  hAdapCuBoxth232M1->DrawNormalized("SAME");
  hAdap50mKth232M1->DrawNormalized("SAME");
  hAdap600mKth232M1->DrawNormalized("SAME");
  hAdapIVCth232M1->DrawNormalized("SAME");
  hAdapOVCth232M1->DrawNormalized("SAME");
  hAdapPbRomth232M1->DrawNormalized("SAME");
  hAdapMBth232M1->DrawNormalized("SAME");
  hAdapSIth232M1->DrawNormalized("SAME");

  legth1->AddEntry(hAdapTeO2th232M1, "TeO2", "l");  
  legth1->AddEntry(hAdapCuFrameth232M1, "CuFrame", "l");
  legth1->AddEntry(hAdapCuBoxth232M1, "CuBox", "l");
  legth1->AddEntry(hAdap50mKth232M1, "50mK", "l");
  legth1->AddEntry(hAdap600mKth232M1, "600mK", "l");
  legth1->AddEntry(hAdapIVCth232M1, "IVC", "l");
  legth1->AddEntry(hAdapOVCth232M1, "OVC", "l");
  legth1->AddEntry(hAdapPbRomth232M1, "PbRom", "l");
  legth1->AddEntry(hAdapMBth232M1, "MB", "l");
  legth1->AddEntry(hAdapSIth232M1, "SI", "l");
  legth1->Draw();


  TCanvas *cTh2322 = new TCanvas("cTh2322", "cTh2322", 1200, 800);
  cTh2322->SetLogy();

  hAdapTeO2th232M2->SetLineColor(1);
  hAdapCuFrameth232M2->SetLineColor(2);
  hAdapCuBoxth232M2->SetLineColor(3);
  hAdap50mKth232M2->SetLineColor(4);
  hAdap600mKth232M2->SetLineColor(5);
  hAdapIVCth232M2->SetLineColor(6);
  hAdapOVCth232M2->SetLineColor(7);
  hAdapPbRomth232M2->SetLineColor(8);
  hAdapMBth232M2->SetLineColor(9);
  hAdapSIth232M2->SetLineColor(11);


  hAdapTeO2th232M2->DrawNormalized();
  hAdapCuFrameth232M2->DrawNormalized("SAME");
  hAdapCuBoxth232M2->DrawNormalized("SAME");
  hAdap50mKth232M2->DrawNormalized("SAME");
  hAdap600mKth232M2->DrawNormalized("SAME");
  hAdapIVCth232M2->DrawNormalized("SAME");
  hAdapOVCth232M2->DrawNormalized("SAME");
  hAdapPbRomth232M2->DrawNormalized("SAME");
  hAdapMBth232M2->DrawNormalized("SAME");
  hAdapSIth232M2->DrawNormalized("SAME");

  legth2->AddEntry(hAdapTeO2th232M2, "TeO2", "l");  
  legth2->AddEntry(hAdapCuFrameth232M2, "CuFrame", "l");
  legth2->AddEntry(hAdapCuBoxth232M2, "CuBox", "l");
  legth2->AddEntry(hAdap50mKth232M2, "50mK", "l");
  legth2->AddEntry(hAdap600mKth232M2, "600mK", "l");
  legth2->AddEntry(hAdapIVCth232M2, "IVC", "l");
  legth2->AddEntry(hAdapOVCth232M2, "OVC", "l");
  legth2->AddEntry(hAdapPbRomth232M2, "PbRom", "l");
  legth2->AddEntry(hAdapMBth232M2, "MB", "l");
  legth2->AddEntry(hAdapSIth232M2, "SI", "l");
  legth2->Draw();

  TCanvas *cU2381 = new TCanvas("cU2381", "cU2381", 1200, 800);
  cU2381->SetLogy();

  hAdapTeO2u238M1->SetLineColor(1);
  hAdapCuFrameu238M1->SetLineColor(2);
  hAdapCuBoxu238M1->SetLineColor(3);
  hAdap50mKu238M1->SetLineColor(4);
  hAdap600mKu238M1->SetLineColor(5);
  hAdapIVCu238M1->SetLineColor(6);
  hAdapOVCu238M1->SetLineColor(7);
  hAdapPbRomu238M1->SetLineColor(8);
  hAdapMBu238M1->SetLineColor(9);
  hAdapSIu238M1->SetLineColor(11);


  hAdapTeO2u238M1->DrawNormalized();
  hAdapCuFrameu238M1->DrawNormalized("SAME");
  hAdapCuBoxu238M1->DrawNormalized("SAME");
  hAdap50mKu238M1->DrawNormalized("SAME");
  hAdap600mKu238M1->DrawNormalized("SAME");
  hAdapIVCu238M1->DrawNormalized("SAME");
  hAdapOVCu238M1->DrawNormalized("SAME");
  hAdapPbRomu238M1->DrawNormalized("SAME");
  hAdapMBu238M1->DrawNormalized("SAME");
  hAdapSIu238M1->DrawNormalized("SAME");

  legu1->AddEntry(hAdapTeO2u238M2, "TeO2", "l");  
  legu1->AddEntry(hAdapCuFrameu238M2, "CuFrame", "l");
  legu1->AddEntry(hAdapCuBoxu238M2, "CuBox", "l");
  legu1->AddEntry(hAdap50mKu238M2, "50mK", "l");
  legu1->AddEntry(hAdap600mKu238M2, "600mK", "l");
  legu1->AddEntry(hAdapIVCu238M2, "IVC", "l");
  legu1->AddEntry(hAdapOVCu238M2, "OVC", "l");
  legu1->AddEntry(hAdapPbRomu238M2, "PbRom", "l");
  legu1->AddEntry(hAdapMBu238M2, "MB", "l");
  legu1->AddEntry(hAdapSIu238M2, "SI", "l");
  legu1->Draw();


  TCanvas *cU2382 = new TCanvas("cU2382", "cU2382", 1200, 800);
  cU2382->SetLogy();

  hAdapTeO2u238M2->SetLineColor(1);
  hAdapCuFrameu238M2->SetLineColor(2);
  hAdapCuBoxu238M2->SetLineColor(3);
  hAdap50mKu238M2->SetLineColor(4);
  hAdap600mKu238M2->SetLineColor(5);
  hAdapIVCu238M2->SetLineColor(6);
  hAdapOVCu238M2->SetLineColor(7);
  hAdapPbRomu238M2->SetLineColor(8);
  hAdapMBu238M2->SetLineColor(9);
  hAdapSIu238M2->SetLineColor(11);


  hAdapTeO2u238M2->DrawNormalized();
  hAdapCuFrameu238M2->DrawNormalized("SAME");
  hAdapCuBoxu238M2->DrawNormalized("SAME");
  hAdap50mKu238M2->DrawNormalized("SAME");
  hAdap600mKu238M2->DrawNormalized("SAME");
  hAdapIVCu238M2->DrawNormalized("SAME");
  hAdapOVCu238M2->DrawNormalized("SAME");
  hAdapPbRomu238M2->DrawNormalized("SAME");
  hAdapMBu238M2->DrawNormalized("SAME");
  hAdapSIu238M2->DrawNormalized("SAME");

  legu2->AddEntry(hAdapTeO2u238M2, "TeO2", "l");  
  legu2->AddEntry(hAdapCuFrameu238M2, "CuFrame", "l");
  legu2->AddEntry(hAdapCuBoxu238M2, "CuBox", "l");
  legu2->AddEntry(hAdap50mKu238M2, "50mK", "l");
  legu2->AddEntry(hAdap600mKu238M2, "600mK", "l");
  legu2->AddEntry(hAdapIVCu238M2, "IVC", "l");
  legu2->AddEntry(hAdapOVCu238M2, "OVC", "l");
  legu2->AddEntry(hAdapPbRomu238M2, "PbRom", "l");
  legu2->AddEntry(hAdapMBu238M2, "MB", "l");
  legu2->AddEntry(hAdapSIu238M2, "SI", "l");
  legu2->Draw();




  TCanvas *cK401 = new TCanvas("cK401", "cK401", 1200, 800);
  cK401->SetLogy();

  hAdapTeO2k40M1->SetLineColor(1);
  hAdapCuFramek40M1->SetLineColor(2);
  hAdapCuBoxk40M1->SetLineColor(3);
  hAdap50mKk40M1->SetLineColor(4);
  hAdap600mKk40M1->SetLineColor(5);
  hAdapIVCk40M1->SetLineColor(6);
  hAdapOVCk40M1->SetLineColor(7);
  hAdapPbRomk40M1->SetLineColor(8);
  hAdapMBk40M1->SetLineColor(9);
  hAdapSIk40M1->SetLineColor(11);


  hAdapTeO2k40M1->DrawNormalized();
  hAdapCuFramek40M1->DrawNormalized("SAME");
  hAdapCuBoxk40M1->DrawNormalized("SAME");
  hAdap50mKk40M1->DrawNormalized("SAME");
  hAdap600mKk40M1->DrawNormalized("SAME");
  hAdapIVCk40M1->DrawNormalized("SAME");
  hAdapOVCk40M1->DrawNormalized("SAME");
  hAdapPbRomk40M1->DrawNormalized("SAME");
  hAdapMBk40M1->DrawNormalized("SAME");
  hAdapSIk40M1->DrawNormalized("SAME");

  legk1->AddEntry(hAdapTeO2k40M1, "TeO2", "l");  
  legk1->AddEntry(hAdapCuFramek40M1, "CuFrame", "l");
  legk1->AddEntry(hAdapCuBoxk40M1, "CuBox", "l");
  legk1->AddEntry(hAdap50mKk40M1, "50mK", "l");
  legk1->AddEntry(hAdap600mKk40M1, "600mK", "l");
  legk1->AddEntry(hAdapIVCk40M1, "IVC", "l");
  legk1->AddEntry(hAdapOVCk40M1, "OVC", "l");
  legk1->AddEntry(hAdapPbRomk40M1, "PbRom", "l");
  legk1->AddEntry(hAdapMBk40M1, "MB", "l");
  legk1->AddEntry(hAdapSIk40M1, "SI", "l");
  legk1->Draw();



  TCanvas *cK402 = new TCanvas("cK402", "cK402", 1200, 800);
  cK402->SetLogy();

  hAdapTeO2k40M2->SetLineColor(1);
  hAdapCuFramek40M2->SetLineColor(2);
  hAdapCuBoxk40M2->SetLineColor(3);
  hAdap50mKk40M2->SetLineColor(4);
  hAdap600mKk40M2->SetLineColor(5);
  hAdapIVCk40M2->SetLineColor(6);
  hAdapOVCk40M2->SetLineColor(7);
  hAdapPbRomk40M2->SetLineColor(8);
  hAdapMBk40M2->SetLineColor(9);
  hAdapSIk40M2->SetLineColor(11);


  hAdapTeO2k40M2->DrawNormalized();
  hAdapCuFramek40M2->DrawNormalized("SAME");
  hAdapCuBoxk40M2->DrawNormalized("SAME");
  hAdap50mKk40M2->DrawNormalized("SAME");
  hAdap600mKk40M2->DrawNormalized("SAME");
  hAdapIVCk40M2->DrawNormalized("SAME");
  hAdapOVCk40M2->DrawNormalized("SAME");
  hAdapPbRomk40M2->DrawNormalized("SAME");
  hAdapMBk40M2->DrawNormalized("SAME");
  hAdapSIk40M2->DrawNormalized("SAME");

  legk2->AddEntry(hAdapTeO2k40M2, "TeO2", "l");  
  legk2->AddEntry(hAdapCuFramek40M2, "CuFrame", "l");
  legk2->AddEntry(hAdapCuBoxk40M2, "CuBox", "l");
  legk2->AddEntry(hAdap50mKk40M2, "50mK", "l");
  legk2->AddEntry(hAdap600mKk40M2, "600mK", "l");
  legk2->AddEntry(hAdapIVCk40M2, "IVC", "l");
  legk2->AddEntry(hAdapOVCk40M2, "OVC", "l");
  legk2->AddEntry(hAdapPbRomk40M2, "PbRom", "l");
  legk2->AddEntry(hAdapMBk40M2, "MB", "l");
  legk2->AddEntry(hAdapSIk40M2, "SI", "l");
  legk2->Draw();


  TCanvas *cCo601 = new TCanvas("cCo601", "cCo601", 1200, 800);
  cCo601->SetLogy();

  hAdapTeO2co60M1->SetLineColor(1);
  hAdapCuFrameco60M1->SetLineColor(2);
  hAdapCuBoxco60M1->SetLineColor(3);
  hAdap50mKco60M1->SetLineColor(4);
  hAdap600mKco60M1->SetLineColor(5);
  hAdapIVCco60M1->SetLineColor(6);
  hAdapOVCco60M1->SetLineColor(7);
  hAdapPbRomco60M1->SetLineColor(8);
  hAdapMBco60M1->SetLineColor(9);


  hAdapTeO2co60M1->DrawNormalized();
  hAdapCuFrameco60M1->DrawNormalized("SAME");
  hAdapCuBoxco60M1->DrawNormalized("SAME");
  hAdap50mKco60M1->DrawNormalized("SAME");
  hAdap600mKco60M1->DrawNormalized("SAME");
  hAdapIVCco60M1->DrawNormalized("SAME");
  hAdapOVCco60M1->DrawNormalized("SAME");
  hAdapPbRomco60M1->DrawNormalized("SAME");
  hAdapMBco60M1->DrawNormalized("SAME");

  legco1->AddEntry(hAdapTeO2co60M1, "TeO2", "l");  
  legco1->AddEntry(hAdapCuFrameco60M1, "CuFrame", "l");
  legco1->AddEntry(hAdapCuBoxco60M1, "CuBox", "l");
  legco1->AddEntry(hAdap50mKco60M1, "50mK", "l");
  legco1->AddEntry(hAdap600mKco60M1, "600mK", "l");
  legco1->AddEntry(hAdapIVCco60M1, "IVC", "l");
  legco1->AddEntry(hAdapOVCco60M1, "OVC", "l");
  legco1->AddEntry(hAdapPbRomco60M1, "PbRom", "l");
  legco1->AddEntry(hAdapMBco60M1, "MB", "l");
  legco1->Draw();



  TCanvas *cCo602 = new TCanvas("cCo602", "cCo602", 1200, 800);
  cCo602->SetLogy();

  hAdapTeO2co60M2->SetLineColor(1);
  hAdapCuFrameco60M2->SetLineColor(2);
  hAdapCuBoxco60M2->SetLineColor(3);
  hAdap50mKco60M2->SetLineColor(4);
  hAdap600mKco60M2->SetLineColor(5);
  hAdapIVCco60M2->SetLineColor(6);
  hAdapOVCco60M2->SetLineColor(7);
  hAdapPbRomco60M2->SetLineColor(8);
  hAdapMBco60M2->SetLineColor(9);


  hAdapTeO2co60M2->DrawNormalized();
  hAdapCuFrameco60M2->DrawNormalized("SAME");
  hAdapCuBoxco60M2->DrawNormalized("SAME");
  hAdap50mKco60M2->DrawNormalized("SAME");
  hAdap600mKco60M2->DrawNormalized("SAME");
  hAdapIVCco60M2->DrawNormalized("SAME");
  hAdapOVCco60M2->DrawNormalized("SAME");
  hAdapPbRomco60M2->DrawNormalized("SAME");
  hAdapMBco60M2->DrawNormalized("SAME");

  legco2->AddEntry(hAdapTeO2co60M2, "TeO2", "l");  
  legco2->AddEntry(hAdapCuFrameco60M2, "CuFrame", "l");
  legco2->AddEntry(hAdapCuBoxco60M2, "CuBox", "l");
  legco2->AddEntry(hAdap50mKco60M2, "50mK", "l");
  legco2->AddEntry(hAdap600mKco60M2, "600mK", "l");
  legco2->AddEntry(hAdapIVCco60M2, "IVC", "l");
  legco2->AddEntry(hAdapOVCco60M2, "OVC", "l");
  legco2->AddEntry(hAdapPbRomco60M2, "PbRom", "l");
  legco2->AddEntry(hAdapMBco60M2, "MB", "l");
  legco2->Draw();


}
