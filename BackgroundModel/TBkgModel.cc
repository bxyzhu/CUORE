#include "TMinuit.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TBkgModel.hh"
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

ClassImp(TBkgModel)
  
//first set up a global function that calls your classes method that calculates the quantity to minimize
void myExternal_FCNAdap(int &n, double *grad, double &fval, double x[], int code)
{
  // Required External Wrapper for function to be minimized by Minuit 
 
  // This gets called for each value of the parameters minuit tries
  // here the x array contains the parameters you are trying to fit
  
  // here myClass should inherit from TObject
  TBkgModel* Obj = (TBkgModel*)gMinuit->GetObjectFit(); 

  // implement a method in your class for setting the parameters and thus update the parameters of your fitter class 
  for(int i = 0; i < 44; i++ )
  {
    Obj->SetParameters(i, x[i]);
  }
  // Implement a method in your class that calculates the quantity you want to minimize, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimize fval
    Obj->UpdateModelAdaptive();
    fval = Obj->GetChiSquareAdaptive();
}

TBkgModel::TBkgModel()
{

}

TBkgModel::TBkgModel(double fFitMin, double fFitMax, int fBinBase, int fDataSet, bool fSave)
{

  TBkgModelSource *fSource = new TBkgModelSource(fFitMin, fFitMax, fBinBase, fDataSet);

  bSave = fSave;
  bToyData = false;

  tTime = new TDatime();
  dNParam = 44; // number of fitting parameters
  dNumCalls = 0;
  dSecToYears = 1./(60*60*24*365);

  minuit = new TMinuit(dNParam);

  TF1 *fEff = new TF1("fEff", "[0]+[1]/(1+[2]*exp(-[3]*x)) + [4]/(1+[5]*exp(-[6]*x))", dMinEnergy, dMaxEnergy);
  fEff->SetParameters(-4.71e-2, 1.12e-1, 2.29, -8.81e-5, 9.68e-1, 2.09, 1.58e-2);

  hEfficiency = new TH1D("hEfficiency", "", dNBins, dMinEnergy, dMaxEnergy);

  for(int i = 1; i <= hEfficiency->GetNbinsX(); i++)
  {
    hEfficiency->SetBinContent(i, fEff->Eval(hEfficiency->GetBinCenter(i)));
  }

  fDataHistoM1->Divide( hEfficiency );
  fDataHistoM2->Divide( hEfficiency );
  fDataHistoM2Sum->Divide( hEfficiency );


  // Loads all of the PDFs from file
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

  dBestChiSq = 0; // Chi-Squared from best fit (for ProfileNLL calculation)
  // Do the fit now if no other tests are needed 
  nLoop = 0;
  DoTheFitAdaptive(0,0);
  // DoTheFitAdaptive(0.0674202742, 0.0263278758);  
  // ProfileNLL(0.0685222152, 3968.95); 
  // ProfileNLL2D(0.0674202742, 0.0000003189, 3754);

  // Number of Toy fits
  if(bToyData)ToyFit(1);

}

// Probably needs updating  
TBkgModel::~TBkgModel()
{

}

double TBkgModel::GetChiSquareAdaptive()
{
  double chiSquare = 0.;
  double datam1_i, errm1_i;
  double datam2_i, errm2_i;
  double datam2sum_i, errm2sum_i;
  double modelm1_i, modelm2_i, modelm2sum_i;


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


/*
  for(int i = dFitMinBinM2Sum; i < dFitMaxBinM2Sum; i++)
  {
    if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 3150 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 3400)continue;
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


// Prints parameters, needs update 11-06-2014
void TBkgModel::PrintParameters()
{
  for(int i = 0; i < TBkgModel::dNParam; i++)
  {
    cout << i << " " << minuit->fCpnam[i] << ": " << fParameters[i] << " +/- " << fParError[i] << endl;
  }
}

void TBkgModel::PrintParActivity()
{ 
  // This only gets the number of counts corrected for detector efficiency
  for(int i = 0; i < TBkgModel::dNParam; i++)
  {
  cout << i << " " << minuit->fCpnam[i] << " activity: " << fParActivity[i] << endl;
  }

}

// Resets all parameters to 0
void TBkgModel::ResetParameters()
{
  dNumCalls = 0;
  dChiSquare = 0;
  
  fModelTotM1->Reset();
  fModelTotthM1->Reset();
  fModelTotuM1->Reset();
  fModelTotkM1->Reset();
  fModelTotcoM1->Reset();
  fModelTotmnM1->Reset();
  fModelTotNDBDM1->Reset();
  fModelTot2NDBDM1->Reset();
  fModelTotbiM1->Reset();
  fModelTotbi2M1->Reset();
  fModelTotptM1->Reset();
  fModelTotpbM1->Reset();
  fModelTotcsM1->Reset();
  fModelTotco2M1->Reset();
  fModelTotteo2M1->Reset();

  fModelTotSthM1->Reset();
  fModelTotSuM1->Reset();
  fModelTotSpoM1->Reset();
  fModelTotSpbM1->Reset();
  fModelTotExtM1->Reset();

  fModelTotAdapM1->Reset();
  fModelTotAdapthM1->Reset();
  fModelTotAdapuM1->Reset();
  fModelTotAdapkM1->Reset();
  fModelTotAdapcoM1->Reset();
  fModelTotAdapmnM1->Reset();
  fModelTotAdapNDBDM1->Reset();
  fModelTotAdap2NDBDM1->Reset();
  fModelTotAdapbiM1->Reset();
  fModelTotAdapbi2M1->Reset();
  fModelTotAdapptM1->Reset();
  fModelTotAdappbM1->Reset();
  fModelTotAdapcsM1->Reset();
  fModelTotAdapco2M1->Reset();
  fModelTotAdapteo2M1->Reset();

  fModelTotAdapSthM1->Reset();
  fModelTotAdapSuM1->Reset();
  fModelTotAdapSpoM1->Reset();
  fModelTotAdapSpbM1->Reset();
  fModelTotAdapExtM1->Reset();

  fModelTotAdapFudgeM1->Reset();

  // Total PDFs M2
  fModelTotM2->Reset();
  fModelTotthM2->Reset();
  fModelTotuM2->Reset();
  fModelTotkM2->Reset();
  fModelTotcoM2->Reset();
  fModelTotmnM2->Reset();
  fModelTotNDBDM2->Reset();
  fModelTot2NDBDM2->Reset();
  fModelTotbiM2->Reset();
  fModelTotbi2M2->Reset();
  fModelTotptM2->Reset();
  fModelTotpbM2->Reset();
  fModelTotcsM2->Reset();
  fModelTotco2M2->Reset();
  fModelTotteo2M2->Reset();

  fModelTotSthM2->Reset();
  fModelTotSuM2->Reset();
  fModelTotSpoM2->Reset();
  fModelTotSpbM2->Reset();
  fModelTotExtM2->Reset();

  fModelTotAdapM2->Reset();
  fModelTotAdapthM2->Reset();
  fModelTotAdapuM2->Reset();
  fModelTotAdapkM2->Reset();
  fModelTotAdapcoM2->Reset();
  fModelTotAdapmnM2->Reset();
  fModelTotAdapNDBDM2->Reset();
  fModelTotAdap2NDBDM2->Reset();
  fModelTotAdapbiM2->Reset();
  fModelTotAdapbi2M2->Reset();
  fModelTotAdapptM2->Reset();
  fModelTotAdappbM2->Reset();
  fModelTotAdapcsM2->Reset();
  fModelTotAdapco2M2->Reset();
  fModelTotAdapteo2M2->Reset();

  fModelTotAdapSthM2->Reset();
  fModelTotAdapSuM2->Reset();
  fModelTotAdapSpoM2->Reset();
  fModelTotAdapSpbM2->Reset();
  fModelTotAdapExtM2->Reset();

  fModelTotAdapFudgeM2->Reset();


  // Total PDFs M2Sum
  fModelTotM2Sum->Reset();
  fModelTotthM2Sum->Reset();
  fModelTotuM2Sum->Reset();
  fModelTotkM2Sum->Reset();
  fModelTotcoM2Sum->Reset();
  fModelTotmnM2Sum->Reset();
  fModelTotNDBDM2Sum->Reset();
  fModelTot2NDBDM2Sum->Reset();
  fModelTotbiM2Sum->Reset();
  fModelTotbi2M2Sum->Reset();
  fModelTotptM2Sum->Reset();
  fModelTotpbM2Sum->Reset();
  fModelTotcsM2Sum->Reset();
  fModelTotco2M2Sum->Reset();
  fModelTotteo2M2Sum->Reset();

  fModelTotSthM2Sum->Reset();
  fModelTotSuM2Sum->Reset();
  fModelTotSpoM2Sum->Reset();
  fModelTotSpbM2Sum->Reset();
  fModelTotExtM2Sum->Reset();

  fModelTotAdapM2Sum->Reset();
  fModelTotAdapthM2Sum->Reset();
  fModelTotAdapuM2Sum->Reset();
  fModelTotAdapkM2Sum->Reset();
  fModelTotAdapcoM2Sum->Reset();
  fModelTotAdapmnM2Sum->Reset();
  fModelTotAdapNDBDM2Sum->Reset();
  fModelTotAdap2NDBDM2Sum->Reset();
  fModelTotAdapbiM2Sum->Reset();
  fModelTotAdapbi2M2Sum->Reset();
  fModelTotAdapptM2Sum->Reset();
  fModelTotAdappbM2Sum->Reset();
  fModelTotAdapcsM2Sum->Reset();
  fModelTotAdapco2M2Sum->Reset();
  fModelTotAdapteo2M2Sum->Reset();

  fModelTotAdapSthM2Sum->Reset();
  fModelTotAdapSuM2Sum->Reset();
  fModelTotAdapSpoM2Sum->Reset();
  fModelTotAdapSpbM2Sum->Reset();
  fModelTotAdapExtM2Sum->Reset();

  fModelTotAdapFudgeM2Sum->Reset();

  for(int i = 0; i < TBkgModel::dNParam; i++)
  {
    fParameters[i] = 0;
    fParError[i] = 0;
  }
}

// Set Parameters in Model
void TBkgModel::SetParameters(int index, double value)
{
	// Change the index max depending on model
	if(index > dNParam) cout << "Index too large" << endl;
	else fParameters[index] = value;
}

void TBkgModel::SetParEfficiency()
{
   // fParEfficiencyM1[0] = 0.8897514; // TeO2 0nu
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

   // fParMass[0] = 39000/1000.; // TeO2 0nu
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
void TBkgModel::UpdateModelAdaptive()
{
  if(fModelTotAdapM1 == NULL || fModelTotAdapM2 == NULL || fModelTotAdapM2Sum == NULL) 
  {
    cout << "Model Histogram Not Created" << endl;
    return;
  }

  // Reset all bins in model histograms in every loop
  fModelTotAdapM1->Reset();
  fModelTotAdapM2->Reset();
  fModelTotAdapM2Sum->Reset();

  dNumCalls++;
  if(dNumCalls%1000==0)
  {
    cout << "Call #: "<< dNumCalls << endl;
  }

/////// Create model
/////// M1
  fModelTotAdapM1->Add( hAdapTeO22nuM1,              dDataIntegralM1*fParameters[0]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFrameco60M1,             dDataIntegralM1*fParameters[1]);
  fModelTotAdapM1->Add( hAdapTeO2th232onlyM1,        dDataIntegralM1*fParameters[2]);
  fModelTotAdapM1->Add( hAdapTeO2th230onlyM1,        dDataIntegralM1*fParameters[3]);
  fModelTotAdapM1->Add( hAdapTeO2Sxth232onlyM1_001,      dDataIntegralM1*fParameters[4]);
  fModelTotAdapM1->Add( hAdapTeO2Sxra228pb208M1_001, dDataIntegralM1*fParameters[5]);
  fModelTotAdapM1->Add( hAdapTeO2Sxu238th230M1_001,  dDataIntegralM1*fParameters[6]);
  fModelTotAdapM1->Add( hAdapTeO2Sxth230onlyM1_001,  dDataIntegralM1*fParameters[7]);
  fModelTotAdapM1->Add( hAdapTeO2Sxra226pb210M1_001, dDataIntegralM1*fParameters[8]);
  fModelTotAdapM1->Add( hAdapTeO2Sxpb210M1_1,        dDataIntegralM1*fParameters[9]);
  fModelTotAdapM1->Add( hAdapTeO2Sxpb210M1_001,      dDataIntegralM1*fParameters[10]);

  fModelTotAdapM1->Add( hAdapCuBox_CuFrameth232M1_10,   dDataIntegralM1*fParameters[11]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFrameu238M1_10,    dDataIntegralM1*fParameters[12]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFramepb210M1_01,   dDataIntegralM1*fParameters[13]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFramepb210M1_001,  dDataIntegralM1*fParameters[14]);

  fModelTotAdapM1->Add( hAdapPbRomk40M1,       dDataIntegralM1*fParameters[15]);
  fModelTotAdapM1->Add( hAdapOVCth232M1,     dDataIntegralM1*fParameters[16]);
  fModelTotAdapM1->Add( hAdapOVCu238M1,      dDataIntegralM1*fParameters[17]);
  fModelTotAdapM1->Add( hAdapOVCco60M1,      dDataIntegralM1*fParameters[18]);
  fModelTotAdapM1->Add( hAdapOVCk40M1,       dDataIntegralM1*fParameters[19]);
  fModelTotAdapM1->Add( hAdapExtPbbi210M1,   dDataIntegralM1*fParameters[20]);

  fModelTotAdapM1->Add( hAdapCuBox_CuFrameth232M1,     dDataIntegralM1*fParameters[21]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFrameu238M1,      dDataIntegralM1*fParameters[22]);
  fModelTotAdapM1->Add( hAdapPbRomcs137M1,             dDataIntegralM1*fParameters[23]);

  fModelTotAdapM1->Add( hAdapTeO2Sxth232M1_1,    dDataIntegralM1*fParameters[24]);
  fModelTotAdapM1->Add( hAdapTeO2Sxth232M1_10,    dDataIntegralM1*fParameters[25]);
  fModelTotAdapM1->Add( hAdapTeO2Sxu238M1_1,     dDataIntegralM1*fParameters[26]);
  fModelTotAdapM1->Add( hAdapTeO2Sxu238M1_10,     dDataIntegralM1*fParameters[27]);
  fModelTotAdapM1->Add( hAdapTeO2Sxpb210M1_10,     dDataIntegralM1*fParameters[28]);

  fModelTotAdapM1->Add( hAdapCuBox_CuFrameth232M1_1,    dDataIntegralM1*fParameters[29]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFrameth232M1_01,    dDataIntegralM1*fParameters[30]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFrameth232M1_001,    dDataIntegralM1*fParameters[31]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFrameu238M1_1,     dDataIntegralM1*fParameters[32]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFrameu238M1_01,     dDataIntegralM1*fParameters[33]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFrameu238M1_001,     dDataIntegralM1*fParameters[34]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFramepb210M1_10,     dDataIntegralM1*fParameters[35]);
  fModelTotAdapM1->Add( hAdapCuBox_CuFramepb210M1_1,     dDataIntegralM1*fParameters[36]);

  // fModelTotAdapM1->Add( hAdapCuBox_CuFrameth232M1_5,     dDataIntegralM1*fParameters[37]);
  // fModelTotAdapM1->Add( hAdapCuBox_CuFrameu238M1_5,      dDataIntegralM1*fParameters[38]);
  // fModelTotAdapM1->Add( hAdapCuBox_CuFramepb210M1_5,      dDataIntegralM1*fParameters[39]);

  fModelTotAdapM1->Add( hAdapInternalth232M1,     dDataIntegralM1*fParameters[37]);
  fModelTotAdapM1->Add( hAdapInternalu238M1,      dDataIntegralM1*fParameters[38]);
  fModelTotAdapM1->Add( hAdapInternalco60M1,      dDataIntegralM1*fParameters[39]);
  fModelTotAdapM1->Add( hAdapInternalk40M1,       dDataIntegralM1*fParameters[40]);

  fModelTotAdapM1->Add( hAdapPbRomth232M1,     dDataIntegralM1*fParameters[41]);
  fModelTotAdapM1->Add( hAdapPbRomu238M1,      dDataIntegralM1*fParameters[42]);
  fModelTotAdapM1->Add( hAdapPbRomco60M1,      dDataIntegralM1*fParameters[43]);

/////// M2
  fModelTotAdapM2->Add( hAdapTeO22nuM2,              dDataIntegralM1*fParameters[0]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFrameco60M2,             dDataIntegralM1*fParameters[1]);
  fModelTotAdapM2->Add( hAdapTeO2th232onlyM2,        dDataIntegralM1*fParameters[2]);
  fModelTotAdapM2->Add( hAdapTeO2th230onlyM2,        dDataIntegralM1*fParameters[3]);
  fModelTotAdapM2->Add( hAdapTeO2Sxth232onlyM2_001,      dDataIntegralM1*fParameters[4]);
  fModelTotAdapM2->Add( hAdapTeO2Sxra228pb208M2_001, dDataIntegralM1*fParameters[5]);
  fModelTotAdapM2->Add( hAdapTeO2Sxu238th230M2_001,  dDataIntegralM1*fParameters[6]);
  fModelTotAdapM2->Add( hAdapTeO2Sxth230onlyM2_001,  dDataIntegralM1*fParameters[7]);
  fModelTotAdapM2->Add( hAdapTeO2Sxra226pb210M2_001, dDataIntegralM1*fParameters[8]);
  fModelTotAdapM2->Add( hAdapTeO2Sxpb210M2_1,        dDataIntegralM1*fParameters[9]);
  fModelTotAdapM2->Add( hAdapTeO2Sxpb210M2_001,      dDataIntegralM1*fParameters[10]);

  fModelTotAdapM2->Add( hAdapCuBox_CuFrameth232M2_10,   dDataIntegralM1*fParameters[11]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFrameu238M2_10,    dDataIntegralM1*fParameters[12]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFramepb210M2_01,   dDataIntegralM1*fParameters[13]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFramepb210M2_001,  dDataIntegralM1*fParameters[14]);


  fModelTotAdapM2->Add( hAdapPbRomk40M2,       dDataIntegralM1*fParameters[15]);
  fModelTotAdapM2->Add( hAdapOVCth232M2,     dDataIntegralM1*fParameters[16]);
  fModelTotAdapM2->Add( hAdapOVCu238M2,      dDataIntegralM1*fParameters[17]);
  fModelTotAdapM2->Add( hAdapOVCco60M2,      dDataIntegralM1*fParameters[18]);
  fModelTotAdapM2->Add( hAdapOVCk40M2,       dDataIntegralM1*fParameters[19]);
  fModelTotAdapM2->Add( hAdapExtPbbi210M2,   dDataIntegralM1*fParameters[20]);

  fModelTotAdapM2->Add( hAdapCuBox_CuFrameth232M2,     dDataIntegralM1*fParameters[21]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFrameu238M2,      dDataIntegralM1*fParameters[22]);
  fModelTotAdapM2->Add( hAdapPbRomcs137M2,             dDataIntegralM1*fParameters[23]);

  fModelTotAdapM2->Add( hAdapTeO2Sxth232M2_1,    dDataIntegralM1*fParameters[24]);
  fModelTotAdapM2->Add( hAdapTeO2Sxth232M2_10,    dDataIntegralM1*fParameters[25]);
  fModelTotAdapM2->Add( hAdapTeO2Sxu238M2_1,     dDataIntegralM1*fParameters[26]);
  fModelTotAdapM2->Add( hAdapTeO2Sxu238M2_10,     dDataIntegralM1*fParameters[27]);
  fModelTotAdapM2->Add( hAdapTeO2Sxpb210M2_10,     dDataIntegralM1*fParameters[28]);

  fModelTotAdapM2->Add( hAdapCuBox_CuFrameth232M2_1,    dDataIntegralM1*fParameters[29]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFrameth232M2_01,    dDataIntegralM1*fParameters[30]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFrameth232M2_001,    dDataIntegralM1*fParameters[31]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFrameu238M2_1,     dDataIntegralM1*fParameters[32]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFrameu238M2_01,     dDataIntegralM1*fParameters[33]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFrameu238M2_001,     dDataIntegralM1*fParameters[34]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFramepb210M2_10,     dDataIntegralM1*fParameters[35]);
  fModelTotAdapM2->Add( hAdapCuBox_CuFramepb210M2_1,     dDataIntegralM1*fParameters[36]);

  // fModelTotAdapM2->Add( hAdapCuBox_CuFrameth232M2_5,     dDataIntegralM1*fParameters[37]);
  // fModelTotAdapM2->Add( hAdapCuBox_CuFrameu238M2_5,      dDataIntegralM1*fParameters[38]);
  // fModelTotAdapM2->Add( hAdapCuBox_CuFramepb210M2_5,      dDataIntegralM1*fParameters[39]);

  fModelTotAdapM2->Add( hAdapInternalth232M2,     dDataIntegralM1*fParameters[37]);
  fModelTotAdapM2->Add( hAdapInternalu238M2,      dDataIntegralM1*fParameters[38]);
  fModelTotAdapM2->Add( hAdapInternalco60M2,      dDataIntegralM1*fParameters[39]);
  fModelTotAdapM2->Add( hAdapInternalk40M2,       dDataIntegralM1*fParameters[40]);

  fModelTotAdapM2->Add( hAdapPbRomth232M2,     dDataIntegralM1*fParameters[41]);
  fModelTotAdapM2->Add( hAdapPbRomu238M2,      dDataIntegralM1*fParameters[42]);
  fModelTotAdapM2->Add( hAdapPbRomco60M2,      dDataIntegralM1*fParameters[43]);

/////// M2Sum
  // fModelTotAdapM2Sum->Add( hAdapTeO20nuM2Sum,              dDataIntegralM1*fParameters[0]);
  fModelTotAdapM2Sum->Add( hAdapTeO22nuM2Sum,              dDataIntegralM1*fParameters[0]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFrameco60M2Sum,             dDataIntegralM1*fParameters[1]);
  fModelTotAdapM2Sum->Add( hAdapTeO2th232onlyM2Sum,        dDataIntegralM1*fParameters[2]);
  fModelTotAdapM2Sum->Add( hAdapTeO2th230onlyM2Sum,        dDataIntegralM1*fParameters[3]);
  fModelTotAdapM2Sum->Add( hAdapTeO2Sxth232onlyM2Sum_001,      dDataIntegralM1*fParameters[4]);
  fModelTotAdapM2Sum->Add( hAdapTeO2Sxra228pb208M2Sum_001, dDataIntegralM1*fParameters[5]);
  fModelTotAdapM2Sum->Add( hAdapTeO2Sxu238th230M2Sum_001,  dDataIntegralM1*fParameters[6]);
  fModelTotAdapM2Sum->Add( hAdapTeO2Sxth230onlyM2Sum_001,  dDataIntegralM1*fParameters[7]);
  fModelTotAdapM2Sum->Add( hAdapTeO2Sxra226pb210M2Sum_001, dDataIntegralM1*fParameters[8]);
  fModelTotAdapM2Sum->Add( hAdapTeO2Sxpb210M2Sum_1,        dDataIntegralM1*fParameters[9]);
  fModelTotAdapM2Sum->Add( hAdapTeO2Sxpb210M2Sum_001,      dDataIntegralM1*fParameters[10]);

  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFrameth232M2Sum_10,   dDataIntegralM1*fParameters[11]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFrameu238M2Sum_10,    dDataIntegralM1*fParameters[12]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFramepb210M2Sum_01,   dDataIntegralM1*fParameters[13]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFramepb210M2Sum_001,  dDataIntegralM1*fParameters[14]);


  fModelTotAdapM2Sum->Add( hAdapPbRomk40M2Sum,       dDataIntegralM1*fParameters[15]);
  fModelTotAdapM2Sum->Add( hAdapOVCth232M2Sum,     dDataIntegralM1*fParameters[16]);
  fModelTotAdapM2Sum->Add( hAdapOVCu238M2Sum,      dDataIntegralM1*fParameters[17]);
  fModelTotAdapM2Sum->Add( hAdapOVCco60M2Sum,      dDataIntegralM1*fParameters[18]);
  fModelTotAdapM2Sum->Add( hAdapOVCk40M2Sum,       dDataIntegralM1*fParameters[19]);
  fModelTotAdapM2Sum->Add( hAdapExtPbbi210M2Sum,   dDataIntegralM1*fParameters[20]);

  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFrameth232M2Sum,     dDataIntegralM1*fParameters[21]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFrameu238M2Sum,      dDataIntegralM1*fParameters[22]);
  fModelTotAdapM2Sum->Add( hAdapPbRomcs137M2Sum,             dDataIntegralM1*fParameters[23]);

  fModelTotAdapM2Sum->Add( hAdapTeO2Sxth232M2Sum_1,    dDataIntegralM1*fParameters[24]);
  fModelTotAdapM2Sum->Add( hAdapTeO2Sxth232M2Sum_10,    dDataIntegralM1*fParameters[25]);
  fModelTotAdapM2Sum->Add( hAdapTeO2Sxu238M2Sum_1,     dDataIntegralM1*fParameters[26]);
  fModelTotAdapM2Sum->Add( hAdapTeO2Sxu238M2Sum_10,     dDataIntegralM1*fParameters[27]);
  fModelTotAdapM2Sum->Add( hAdapTeO2Sxpb210M2Sum_10,     dDataIntegralM1*fParameters[28]);

  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFrameth232M2Sum_1,    dDataIntegralM1*fParameters[29]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFrameth232M2Sum_01,    dDataIntegralM1*fParameters[30]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFrameth232M2Sum_001,    dDataIntegralM1*fParameters[31]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFrameu238M2Sum_1,     dDataIntegralM1*fParameters[32]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFrameu238M2Sum_01,     dDataIntegralM1*fParameters[33]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFrameu238M2Sum_001,     dDataIntegralM1*fParameters[34]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFramepb210M2Sum_10,     dDataIntegralM1*fParameters[35]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFramepb210M2Sum_1,     dDataIntegralM1*fParameters[36]);

  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFrameth232M2Sum_100,     dDataIntegralM1*fParameters[37]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFrameu238M2Sum_100,     dDataIntegralM1*fParameters[38]);
  fModelTotAdapM2Sum->Add( hAdapCuBox_CuFramepb210M2Sum_100,     dDataIntegralM1*fParameters[39]);

}
  
// This method sets up minuit and does the fit
bool TBkgModel::DoTheFitAdaptive(double f2nuValue, double fVariableValue)
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

  //// 50 counts
  minuit->DefineParameter(0, "TeO2 2nu",  f2nuValue, 1E-7, 0, 1.0);
  minuit->DefineParameter(1, "CuBox + CuFrame co60",  0, 1E-7, 0, 1.0);

  minuit->DefineParameter(2, "TeO2 th232 only", 0, 1E-7, 0, 1.0);
  minuit->DefineParameter(3, "TeO2 th230 only", 3.07203e-04, 1E-7, 0, 1.0);  
  minuit->DefineParameter(4, "TeO2 Sx th232 only 0.01", 1.43894e-04, 1E-7, 0, 1.0);
  minuit->DefineParameter(5, "TeO2 Sx ra228 to pb208 0.01", 3.07208e-03, 1E-7, 0, 1.0);
  minuit->DefineParameter(6, "TeO2 Sx u238 to th230 0.01", 1.71025e-03, 1E-7, 0, 1.0);
  minuit->DefineParameter(7, "TeO2 Sx th230 only 0.01", 7.23534e-04, 1E-7, 0, 1.0);
  minuit->DefineParameter(8, "TeO2 Sx ra226 to pb210 0.01", 2.98297e-03, 1E-7, 0, 1.0);
  minuit->DefineParameter(9, "TeO2 Sx pb210 1", 5.30974e-03, 1E-7, 0, 1.0);
  minuit->DefineParameter(10, "TeO2 Sx pb210 0.01", 4.11621e-02, 1E-7, 0, 1.0);
  minuit->DefineParameter(11, "CuBox + CuFrame Sx th232 10", 3.53539e-03, 1E-7, 0, 1.0);
  minuit->DefineParameter(12, "CuBox + CuFrame Sx u238 10 ", 5.80166e-03, 1E-7, 0, 1.0);
  minuit->DefineParameter(13, "CuBox + CuFrame Sx pb210 0.1", 5.91139e-03, 1E-7, 0, 1.0);
  minuit->DefineParameter(14, "CuBox + CuFrame Sx pb210 0.01", 1.79798e-02, 1E-7, 0, 1.0);

  // minuit->DefineParameter(2, "TeO2 th232 only", 0, 1E-7, 0, 1.0);
  // minuit->DefineParameter(3, "TeO2 th230 only", 0, 1E-7, 0, 1.0);  
  // minuit->DefineParameter(4, "TeO2 Sx th232 only 0.01", 0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(5, "TeO2 Sx ra228 to pb208 0.01", 0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(6, "TeO2 Sx u238 to th230 0.01", 0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(7, "TeO2 Sx th230 only 0.01", 0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(8, "TeO2 Sx ra226 to pb210 0.01", 0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(9, "TeO2 Sx pb210 1", 0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(10, "TeO2 Sx pb210 0.01", 0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(11, "CuBox + CuFrame Sx th232 10", 0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(12, "CuBox + CuFrame Sx u238 10 ", 0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(13, "CuBox + CuFrame Sx pb210 0.1", 0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(14, "CuBox + CuFrame Sx pb210 0.01", 0., 1E-7, 0, 1.0);

  minuit->DefineParameter(15, "PbRom k40",  0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(16, "OVC th232",  0., 1E-7, 0, 1.0);  
  // minuit->DefineParameter(17, "OVC u238",  0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(18, "OVC co60",  0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(19, "OVC k40",  0., 1E-7, 0, 1.0);
  minuit->DefineParameter(16, "OVC th232",  9.25179e-02, 1E-7, 0, 1.0);
  minuit->DefineParameter(17, "OVC u238",  1.33486e-01, 1E-7, 0, 1.0);
  minuit->DefineParameter(18, "OVC co60",  1.94100e-02, 1E-7, 0, 1.0);
  minuit->DefineParameter(19, "OVC k40",  9.90257e-02, 1E-7, 0, 1.0);  
  minuit->DefineParameter(20, "External Lead bi210", fVariableValue, 1E-7, 0, 1.0);

  minuit->DefineParameter(21, "CuBox + CuFrame th232",  0, 1E-7, 0, 1.0);
  // minuit->DefineParameter(21, "CuBox + CuFrame th232",  2.66019e-02, 1E-7, 0, 1.0);
  minuit->DefineParameter(22, "CuBox + CuFrame u238",  0, 1E-7, 0, 1.0);
  // minuit->DefineParameter(22, "CuBox + CuFrame u238",  4.98365e-03, 1E-7, 0, 1.0);
  minuit->DefineParameter(23, "PbRom cs137",  0, 1E-7, 0, 1.0);

  minuit->DefineParameter(24, "TeO2 Sx th232 1", 1.12618e-03, 1E-7, 0, 1.0);
  // minuit->DefineParameter(24, "TeO2 Sx th232 1", 0, 1E-7, 0, 1.0);  
  minuit->DefineParameter(25, "TeO2 Sx th232 10", 1.15477e-03, 1E-7, 0, 1.0);
  // minuit->DefineParameter(25, "TeO2 Sx th232 10", 0, 1E-7, 0, 1.0);
  minuit->DefineParameter(26, "TeO2 Sx u238 1", 0., 1E-7, 0, 1.0);
  minuit->DefineParameter(27, "TeO2 Sx u238 10", 0., 1E-7, 0, 1.0);
  minuit->DefineParameter(28, "TeO2 Sx pb210 10", 0., 1E-7, 0, 1.0);

  minuit->DefineParameter(29, "CuBox + CuFrame Sx th232 1", 0., 1E-7, 0, 1.0);
  minuit->DefineParameter(30, "CuBox + CuFrame Sx th232 0.1", 0., 1E-7, 0, 1.0);
  minuit->DefineParameter(31, "CuBox + CuFrame Sx th232 0.01", 5.75719e-04, 1E-7, 0, 1.0);
  // minuit->DefineParameter(31, "CuBox + CuFrame Sx th232 0.01", 0, 1E-7, 0, 1.0);
  minuit->DefineParameter(32, "CuBox + CuFrame Sx u238 1", 0., 1E-7, 0, 1.0);
  minuit->DefineParameter(33, "CuBox + CuFrame Sx u238 0.1", 5.81734e-04, 1E-7, 0, 1.0);
  // minuit->DefineParameter(33, "CuBox + CuFrame Sx u238 0.1", 0, 1E-7, 0, 1.0);
  minuit->DefineParameter(34, "CuBox + CuFrame Sx u238 0.01", 0., 1E-7, 0, 1.0);
  minuit->DefineParameter(35, "CuBox + CuFrame Sx pb210 10", 6.09448e-03, 1E-7, 0, 1.0);
  // minuit->DefineParameter(35, "CuBox + CuFrame Sx pb210 10", 0, 1E-7, 0, 1.0);
  // minuit->DefineParameter(36, "CuBox + CuFrame Sx pb210 1", 5.21599e-05, 1E-7, 0, 1.0);
  minuit->DefineParameter(36, "CuBox + CuFrame Sx pb210 1", 0, 1E-7, 0, 1.0);

  // minuit->DefineParameter(37, "CuBox + CuFrame Sx th232 5", 0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(38, "CuBox + CuFrame Sx u238 5", 0., 1E-7, 0, 1.0);
  // minuit->DefineParameter(39, "CuBox + CuFrame Sx pb210 5", 0, 1E-7, 0, 1.0);


  minuit->DefineParameter(37, "Internal th232",  0., 1E-7, 0, 1.0);
  minuit->DefineParameter(38, "Internal u238",  0., 1E-7, 0, 1.0);
  minuit->DefineParameter(39, "Internal co60",  0., 1E-7, 0, 1.0);    
  minuit->DefineParameter(40, "Internal k40",  0., 1E-7, 0, 1.0);

  minuit->DefineParameter(41, "PbRom th232",  0., 1E-7, 0, 1.0);
  minuit->DefineParameter(42, "PbRom u238",  0., 1E-7, 0, 1.0);
  minuit->DefineParameter(43, "PbRom co60",  0., 1E-7, 0, 1.0);



//////////////////////////////////////

  // Loop to fix parameters
  for(int i = 0; i < dNParam; i++)
  {
    if(bFixedArray[i]) 
    {
      minuit->FixParameter(i);
    // For debugging
    // cout << "Parameter " << fParameterName[i] << " is fixed" << endl; 
    }
  }

   // Uncommend to fix parameters here
   // minuit->FixParameter(0); // TeO2 2nu
   // minuit->FixParameter(1); // CuBox+Frame co60
   // minuit->FixParameter(2); // TeO2 th232 only
   // minuit->FixParameter(3); // TeO2 th230 only
   // minuit->FixParameter(4); // TeO2 Sx th232 only 0.01
   // minuit->FixParameter(5); // TeO2 Sx ra228 to pb208 0.01 
   // minuit->FixParameter(6); // TeO2 Sx u238 to th230 0.01  
   // minuit->FixParameter(7); // TeO2 Sx th230 only 0.01
   // minuit->FixParameter(8); // TeO2 Sx ra226 to pb210 0.01
   // minuit->FixParameter(9); // TeO2 Sx pb210 1 ==> necessary for bin below Po210 peak in M2
   // minuit->FixParameter(10); // TeO2 Sx pb210 0.01 ==> completely necessary for M2 spectrum
   // minuit->FixParameter(11); // CuBox+CuFrame Sx th232 10
   // minuit->FixParameter(12); // CuBox+CuFrame Sx u238 10
   // minuit->FixParameter(13); // CuBox+CuFrame Sx pb210 0.1 ==> useful for below Po210 peak in M1 but doesn't seem absolutely necessary
   // minuit->FixParameter(14); // CuBox+CuFrame Sx pb210 0.01 => necessary for below Po210 peak in M1
   // minuit->FixParameter(15); // PbRom k40 
   // minuit->FixParameter(16); // OVC th232
   // minuit->FixParameter(17); // OVC u238
   // minuit->FixParameter(18); // OVC co60
   // minuit->FixParameter(19); // OVC k40
   // minuit->FixParameter(20); // External Lead bi210 
   // minuit->FixParameter(21); // CuBox+Frame th232 
   // minuit->FixParameter(22); // CuBox+Frame u238
   // minuit->FixParameter(23); // PbRom cs137

   // minuit->FixParameter(24); // TeO2 Sx th232 1
   // minuit->FixParameter(25); // TeO2 Sx th232 10 
   // minuit->FixParameter(26); // TeO2 Sx u238 1
   // minuit->FixParameter(27); // TeO2 Sx u238 10
   // minuit->FixParameter(28); // TeO2 Sx pb210 10

   // minuit->FixParameter(29); // CuBox + CuFrame Sx th232 1
   // minuit->FixParameter(30); // CuBox + CuFrame Sx th232 0.1 
   // minuit->FixParameter(31); // CuBox + CuFrame Sx th232 0.01
   // minuit->FixParameter(32); // CuBox + CuFrame Sx u238 1
   // minuit->FixParameter(33); // CuBox + CuFrame Sx u238 0.1
   // minuit->FixParameter(34); // CuBox + CuFrame Sx u238 0.01
   // minuit->FixParameter(35); // CuBox + CuFrame Sx pb210 10
   // minuit->FixParameter(36); // CuBox + CuFrame Sx pb210 1

   minuit->FixParameter(37); // 
   minuit->FixParameter(38); // 
   minuit->FixParameter(39); //
   minuit->FixParameter(40); //

   minuit->FixParameter(41); // 
   minuit->FixParameter(42); //
   minuit->FixParameter(43); //

   // Number of free Parameters (for Chi-squared/NDF calculation only)
   dNumFreeParameters = minuit->GetNumPars() - minuit->GetNumFixedPars();
   // dNDF = (dFitMaxBinM1-dFitMinBinM1-dNumFreeParameters); // M1 only
   // dNDF = (dFitMaxBinM2-dFitMinBinM2-dNumFreeParameters); // M2 only
   dNDF = (dFitMaxBinM1+dFitMaxBinM2-dFitMinBinM1-dFitMinBinM2-dNumFreeParameters); // M1 and M2
   // dNDF = (dFitMaxBinM1+dFitMaxBinM2+dFitMaxBinM2Sum - dFitMinBinM1-dFitMinBinM2-dFitMinBinM2Sum - dNumFreeParameters);
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


// M1
  fModelTotAdap2NDBDM1->Add( hAdapTeO22nuM1,              dDataIntegralM1*fParameters[0]);
  // ///// Draw Data M1
  fAdapDataHistoM1->SetLineColor(kBlack);
  fAdapDataHistoM1->SetLineWidth(2);
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

  fModelTotAdapSpoM1->SetLineStyle(2);
  fModelTotAdapSpoM1->SetLineColor(44);

  fModelTotAdapSpbM1->SetLineStyle(2);
  fModelTotAdapSpbM1->SetLineColor(8);

  // fModelTotAdapteo2M1->Draw("SAME");
  // fModelTotAdapPbM1->Draw("SAME");
  fModelTotAdapFudgeM1->SetLineStyle(2);
  fModelTotAdapFudgeM1->SetLineColor(kRed+1);


  TLegend *legfit1 = new TLegend(0.7,0.7,0.95,0.92);
  legfit1->SetFillColor(0);
  legfit1->SetTextSize(0.02);
  legfit1->AddEntry(fModelTotAdapM1, "Total PDF", "l");
  legfit1->AddEntry(fModelTotAdapthM1, "Total th-232", "l");
  legfit1->AddEntry(fModelTotAdapuM1, "Total u-238", "l");
  legfit1->AddEntry(fModelTotAdapkM1, "Total k-40", "l");
  legfit1->AddEntry(fModelTotAdapcoM1, "Total co-60", "l");
  // legfit1->AddEntry(fModelTotAdapNDBDM1, "NDBD", "l");
  legfit1->AddEntry(fModelTotAdap2NDBDM1, "2NDBD", "l");
  legfit1->AddEntry(fModelTotAdapSpbM1, "Surface pb-210", "l");
  legfit1->AddEntry(fModelTotAdapbiM1, "Total External bi-210", "l");  
  legfit1->AddEntry(fModelTotAdapFudgeM1, "Total cs-137", "l");

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
  LineM1->DrawLine(dFitMin, 1, dFitMax-1, 1);

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
  // fModelTotAdapthM1->Draw("SAME");
  // fModelTotAdapuM1->Draw("SAME");
  // fModelTotAdapkM1->Draw("SAME");
  // fModelTotAdapcoM1->Draw("SAME");
  // fModelTotAdapNDBDM1->Draw("SAME");
  // fModelTotAdap2NDBDM1->Draw("SAME");
  // fModelTotAdapbiM1->Draw("SAME");
  // fModelTotAdapmnM1->Draw("SAME");
  // fModelTotAdapSpbM1->Draw("SAME");
  // fModelTotAdapSpoM1->Draw("SAME");
  // fModelTotAdappbM1->Draw("SAME");
  // fModelTotAdapFudgeM1->Draw("SAME");
  // legfit1->Draw();

  // fDataHistoM1->SetLineColor(17);
  // fDataHistoM1->Draw("SAME");

  ///// Draw Data M2
  fAdapDataHistoM2->SetLineColor(kBlack);
  fAdapDataHistoM2->SetLineWidth(2);
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

  fModelTotAdapSpoM2->SetLineStyle(2);
  fModelTotAdapSpoM2->SetLineColor(44);

  fModelTotAdapSpbM2->SetLineStyle(2);
  fModelTotAdapSpbM2->SetLineColor(8);

  fModelTotAdapFudgeM2->SetLineStyle(2);
  fModelTotAdapFudgeM2->SetLineColor(kRed+1);

  TLegend *legfit2 = new TLegend(0.7,0.7,0.95,0.92);
  legfit2->SetFillColor(0);
  legfit2->SetTextSize(0.02);
  legfit2->AddEntry(fModelTotAdapM2, "Total PDF", "l");
  legfit2->AddEntry(fModelTotAdapthM2, "Total th-232", "l");
  legfit2->AddEntry(fModelTotAdapuM2, "Total u-238", "l");
  legfit2->AddEntry(fModelTotAdapkM2, "Total k-40", "l");
  legfit2->AddEntry(fModelTotAdapcoM2, "Total co-60", "l");
  // legfit2->AddEntry(fModelTotAdapNDBDM2, "NDBD", "l");
  legfit2->AddEntry(fModelTotAdap2NDBDM2, "2NDBD", "l");
  legfit2->AddEntry(fModelTotAdapSpbM2, "Surface pb-210", "l");
  legfit2->AddEntry(fModelTotAdapbiM2, "Total External bi-210", "l");  
  legfit2->AddEntry(fModelTotAdapFudgeM2, "Total cs-137", "l");

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
  LineM1->DrawLine(dFitMin, 1, dFitMax-1, 1);



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
  // fModelTotAdapthM2->Draw("SAME");
  // fModelTotAdapuM2->Draw("SAME");
  // fModelTotAdapkM2->Draw("SAME");
  // fModelTotAdapcoM2->Draw("SAME");
  // fModelTotAdapNDBDM2->Draw("SAME");
  // fModelTotAdap2NDBDM2->Draw("SAME");
  // fModelTotAdapbiM2->Draw("SAME");
  // fModelTotAdapmnM2->Draw("SAME");
  // fModelTotAdapSpbM2->Draw("SAME");
  // fModelTotAdapSpoM2->Draw("SAME");
  // fModelTotAdappbM2->Draw("SAME");
  // fModelTotAdapFudgeM2->Draw("SAME");
  // legfit2->Draw();

  // fDataHistoM2->SetLineColor(17);
  // fDataHistoM2->Draw("SAME");


  ///// Draw Data M2Sum
  fAdapDataHistoM2Sum->SetLineColor(kBlack);
  fAdapDataHistoM2Sum->SetLineWidth(2);
  // fAdapDataHistoM2Sum->SetMaximum(9000);
  // fAdapDataHistoM2Sum->GetXaxis()->SetRange(1, fAdapDataHistoM2Sum->FindBin(3000));
  
  fModelTotAdapM2Sum->SetLineColor(2);
  fModelTotAdapM2Sum->SetLineWidth(1);
  fModelTotAdapthM2Sum->SetLineColor(3);
  fModelTotAdapthM2Sum->SetLineStyle(2);
  fModelTotAdapuM2Sum->SetLineColor(4);
  fModelTotAdapuM2Sum->SetLineStyle(2);
  fModelTotAdapkM2Sum->SetLineColor(6);
  fModelTotAdapkM2Sum->SetLineStyle(2);
  fModelTotAdapcoM2Sum->SetLineColor(7);
  fModelTotAdapcoM2Sum->SetLineStyle(2);
  fModelTotAdapNDBDM2Sum->SetLineColor(42);
  fModelTotAdapNDBDM2Sum->SetLineStyle(2);
  fModelTotAdap2NDBDM2Sum->SetLineColor(46);
  fModelTotAdap2NDBDM2Sum->SetLineStyle(2);
  fModelTotAdapbiM2Sum->SetLineColor(5);
  fModelTotAdapbiM2Sum->SetLineStyle(2);
  fModelTotAdapmnM2Sum->SetLineColor(40);
  fModelTotAdapmnM2Sum->SetLineStyle(2);

  fModelTotAdapSpoM2Sum->SetLineStyle(2);
  fModelTotAdapSpoM2Sum->SetLineColor(44);

  fModelTotAdapSpbM2Sum->SetLineStyle(2);
  fModelTotAdapSpbM2Sum->SetLineColor(8);

  fModelTotAdapFudgeM2Sum->SetLineStyle(2);
  fModelTotAdapFudgeM2Sum->SetLineColor(kRed+1);

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
  // fModelTotAdapthM2Sum->Draw("SAME");
  // fModelTotAdapuM2Sum->Draw("SAME");
  // fModelTotAdapkM2Sum->Draw("SAME");
  // fModelTotAdapcoM2Sum->Draw("SAME");
  // fModelTotAdapNDBDM2Sum->Draw("SAME");
  // fModelTotAdap2NDBDM2Sum->Draw("SAME");
  // fModelTotAdapbiM2Sum->Draw("SAME");
  // fModelTotAdapmnM2Sum->Draw("SAME");
  // fModelTotAdapSpbM2Sum->Draw("SAME");
  // fModelTotAdapSpoM2Sum->Draw("SAME");
  // fModelTotAdappbM2Sum->Draw("SAME");
  // fModelTotAdapFudgeM2Sum->Draw("SAME");

  // fDataHistoM2Sum->SetLineColor(17);
  // fDataHistoM2Sum->Draw("SAME");


  // Residuals
  // TCanvas *cResidual1 = new TCanvas("cResidual1", "cResidual1", 1200, 800);

  // TCanvas *cResidual2 = new TCanvas("cResidual2", "cResidual2", 1200, 800);

  // TCanvas *cResidual2Sum = new TCanvas("cResidual2Sum", "cResidual2Sum", 1200, 800);

  // TCanvas *cres1 = new TCanvas("cres1", "cres1", 1600, 600);

/*
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
  for (int i = dFitMinBinM2Sum; i < dFitMaxBinM2Sum; i++)
  {
    dResidualRMSM2Sum += hResidualDistM2Sum->GetBinContent(i)*hResidualDistM2Sum->GetBinContent(i);
  }


  dResidualRMSTot = TMath::Sqrt( (dResidualRMSM1 + dResidualRMSM2)/ (dNDF + dNumFreeParameters) );


  dResidualRMSM1 = TMath::Sqrt(dResidualRMSM1/(dFitMaxBinM1-dFitMinBinM1));
  dResidualRMSM2 = TMath::Sqrt(dResidualRMSM2/(dFitMaxBinM2-dFitMinBinM2));
  dResidualRMSM2Sum = TMath::Sqrt(dResidualRMSM2Sum/(dFitMaxBinM2Sum-dFitMinBinM2Sum));
*/


  //
  double dROIRange = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2570))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2570)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2470)); 
  double d2nbbRange = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2000))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2000)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(500));

  TH1D *fModelAlphaBulk = new TH1D("fModelAlphaBulk", "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelAlphaBulk->Add( hAdapTeO2th232onlyM1,        dDataIntegralM1*fParameters[2]);
  fModelAlphaBulk->Add( hAdapTeO2th230onlyM1,        dDataIntegralM1*fParameters[3]);


  TH1D *fModelAlphaTeSurface1 = new TH1D("fModelAlphaTeSurface1", "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelAlphaTeSurface1->Add( hAdapTeO2Sxpb210M1_1,        dDataIntegralM1*fParameters[9]);
  fModelAlphaTeSurface1->Add( hAdapTeO2Sxth232M1_1,    dDataIntegralM1*fParameters[24]);
  fModelAlphaTeSurface1->Add( hAdapTeO2Sxth232M1_10,    dDataIntegralM1*fParameters[25]);
  fModelAlphaTeSurface1->Add( hAdapTeO2Sxu238M1_1,     dDataIntegralM1*fParameters[26]);
  fModelAlphaTeSurface1->Add( hAdapTeO2Sxu238M1_10,     dDataIntegralM1*fParameters[27]);
  fModelAlphaTeSurface1->Add( hAdapTeO2Sxpb210M1_10,     dDataIntegralM1*fParameters[28]);

  TH1D *fModelAlphaTeSurface2 = new TH1D("fModelAlphaTeSurface2", "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelAlphaTeSurface2->Add( hAdapTeO2Sxth232onlyM1_001,      dDataIntegralM1*fParameters[4]);
  fModelAlphaTeSurface2->Add( hAdapTeO2Sxra228pb208M1_001, dDataIntegralM1*fParameters[5]);
  fModelAlphaTeSurface2->Add( hAdapTeO2Sxu238th230M1_001,  dDataIntegralM1*fParameters[6]);
  fModelAlphaTeSurface2->Add( hAdapTeO2Sxth230onlyM1_001,  dDataIntegralM1*fParameters[7]);
  fModelAlphaTeSurface2->Add( hAdapTeO2Sxra226pb210M1_001, dDataIntegralM1*fParameters[8]);
  fModelAlphaTeSurface2->Add( hAdapTeO2Sxpb210M1_001,      dDataIntegralM1*fParameters[10]);

  TH1D *fModelAlphaCuSurface1 = new TH1D("fModelAlphaCuSurface1", "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelAlphaCuSurface1->Add( hAdapCuBox_CuFrameth232M1_10,   dDataIntegralM1*fParameters[11]);
  fModelAlphaCuSurface1->Add( hAdapCuBox_CuFrameu238M1_10,    dDataIntegralM1*fParameters[12]);
  fModelAlphaCuSurface1->Add( hAdapCuBox_CuFrameth232M1_1,    dDataIntegralM1*fParameters[29]);
  fModelAlphaCuSurface1->Add( hAdapCuBox_CuFrameu238M1_1,     dDataIntegralM1*fParameters[32]);
  fModelAlphaCuSurface1->Add( hAdapCuBox_CuFramepb210M1_10,     dDataIntegralM1*fParameters[35]);
  fModelAlphaCuSurface1->Add( hAdapCuBox_CuFramepb210M1_1,     dDataIntegralM1*fParameters[36]);

  TH1D *fModelAlphaCuSurface2 = new TH1D("fModelAlphaCuSurface2", "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelAlphaCuSurface2->Add( hAdapCuBox_CuFramepb210M1_01,   dDataIntegralM1*fParameters[13]);
  fModelAlphaCuSurface2->Add( hAdapCuBox_CuFramepb210M1_001,  dDataIntegralM1*fParameters[14]);
  fModelAlphaCuSurface2->Add( hAdapCuBox_CuFrameth232M1_01,    dDataIntegralM1*fParameters[30]);
  fModelAlphaCuSurface2->Add( hAdapCuBox_CuFrameth232M1_001,    dDataIntegralM1*fParameters[31]);
  fModelAlphaCuSurface2->Add( hAdapCuBox_CuFrameu238M1_01,     dDataIntegralM1*fParameters[33]);
  fModelAlphaCuSurface2->Add( hAdapCuBox_CuFrameu238M1_001,     dDataIntegralM1*fParameters[34]);


  TH1D *fModelGammaClose = new TH1D("fModelGammaClose", "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelGammaClose->Add( hAdapCuBox_CuFrameco60M1,             dDataIntegralM1*fParameters[1]);
  fModelGammaClose->Add( hAdapCuBox_CuFrameth232M1,     dDataIntegralM1*fParameters[21]);
  fModelGammaClose->Add( hAdapCuBox_CuFrameu238M1,      dDataIntegralM1*fParameters[22]);

  TH1D *fModelGammaMid = new TH1D("fModelGammaMid", "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelGammaMid->Add( hAdapPbRomk40M1,       dDataIntegralM1*fParameters[15]);
  fModelGammaMid->Add( hAdapPbRomcs137M1,             dDataIntegralM1*fParameters[23]);
  fModelGammaMid->Add( hAdapInternalth232M1,     dDataIntegralM1*fParameters[37]);
  fModelGammaMid->Add( hAdapInternalu238M1,      dDataIntegralM1*fParameters[38]);
  fModelGammaMid->Add( hAdapInternalco60M1,      dDataIntegralM1*fParameters[39]);
  fModelGammaMid->Add( hAdapInternalk40M1,       dDataIntegralM1*fParameters[40]);
  fModelGammaMid->Add( hAdapPbRomth232M1,     dDataIntegralM1*fParameters[41]);
  fModelGammaMid->Add( hAdapPbRomu238M1,      dDataIntegralM1*fParameters[42]);
  fModelGammaMid->Add( hAdapPbRomco60M1,      dDataIntegralM1*fParameters[43]);

  TH1D *fModelGammaFar = new TH1D("fModelGammaFar", "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  fModelGammaFar->Add( hAdapOVCth232M1,     dDataIntegralM1*fParameters[16]);
  fModelGammaFar->Add( hAdapOVCu238M1,      dDataIntegralM1*fParameters[17]);
  fModelGammaFar->Add( hAdapOVCco60M1,      dDataIntegralM1*fParameters[18]);
  fModelGammaFar->Add( hAdapOVCk40M1,       dDataIntegralM1*fParameters[19]);
  fModelGammaFar->Add( hAdapExtPbbi210M1,   dDataIntegralM1*fParameters[20]);


  double dAlpha2nbbBulk = fModelAlphaBulk->Integral( fModelAlphaBulk->FindBin(500), fModelAlphaBulk->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double dAlpha2nbbBulkErr = TMath::Sqrt(fModelAlphaBulk->Integral( fModelAlphaBulk->FindBin(500), fModelAlphaBulk->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);
  double dAlpha2nbbTeSurface1 = fModelAlphaTeSurface1->Integral( fModelAlphaTeSurface1->FindBin(500), fModelAlphaTeSurface1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double dAlpha2nbbTeSurface1Err = TMath::Sqrt(fModelAlphaTeSurface1->Integral( fModelAlphaTeSurface1->FindBin(500), fModelAlphaTeSurface1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);
  double dAlpha2nbbTeSurface2 = fModelAlphaTeSurface2->Integral( fModelAlphaTeSurface2->FindBin(500), fModelAlphaTeSurface2->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double dAlpha2nbbTeSurface2Err = TMath::Sqrt(fModelAlphaTeSurface2->Integral( fModelAlphaTeSurface2->FindBin(500), fModelAlphaTeSurface2->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);  
  double dAlpha2nbbCuSurface1 = fModelAlphaCuSurface1->Integral( fModelAlphaCuSurface1->FindBin(500), fModelAlphaCuSurface1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double dAlpha2nbbCuSurface1Err = TMath::Sqrt(fModelAlphaCuSurface1->Integral( fModelAlphaCuSurface1->FindBin(500), fModelAlphaCuSurface1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);
  double dAlpha2nbbCuSurface2 = fModelAlphaCuSurface2->Integral( fModelAlphaCuSurface2->FindBin(500), fModelAlphaCuSurface2->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double dAlpha2nbbCuSurface2Err = TMath::Sqrt(fModelAlphaTeSurface2->Integral( fModelAlphaTeSurface2->FindBin(500), fModelAlphaTeSurface2->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);

  double dGammaClose2nbb = fModelGammaClose->Integral( fModelGammaClose->FindBin(500), fModelGammaClose->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double dGammaClose2nbbErr = TMath::Sqrt(fModelGammaClose->Integral( fModelGammaClose->FindBin(500), fModelGammaClose->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);
  double dGammaMid2nbb = fModelGammaMid->Integral( fModelGammaMid->FindBin(500), fModelGammaMid->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double dGammaMid2nbbErr = TMath::Sqrt(fModelGammaMid->Integral( fModelGammaMid->FindBin(500), fModelGammaMid->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);
  double dGammaFar2nbb = fModelGammaFar->Integral( fModelGammaFar->FindBin(500), fModelGammaFar->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double dGammaFar2nbbErr = TMath::Sqrt(fModelGammaFar->Integral( fModelGammaFar->FindBin(500), fModelGammaFar->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);


  double dAlpha0nbbBulk = fModelAlphaBulk->Integral( fModelAlphaBulk->FindBin(2470), fModelAlphaBulk->FindBin(2570), "width" )/(dROIRange*dLivetimeYr);
  double dAlpha0nbbBulkErr = TMath::Sqrt(fModelAlphaBulk->Integral( fModelAlphaBulk->FindBin(2470), fModelAlphaBulk->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr);
  double dAlpha0nbbTeSurface1 = fModelAlphaTeSurface1->Integral( fModelAlphaTeSurface1->FindBin(2470), fModelAlphaTeSurface1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr);
  double dAlpha0nbbTeSurface1Err = TMath::Sqrt(fModelAlphaTeSurface1->Integral( fModelAlphaTeSurface1->FindBin(2470), fModelAlphaTeSurface1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr);
  double dAlpha0nbbTeSurface2 = fModelAlphaTeSurface2->Integral( fModelAlphaTeSurface2->FindBin(2470), fModelAlphaTeSurface2->FindBin(2570), "width" )/(dROIRange*dLivetimeYr);
  double dAlpha0nbbTeSurface2Err = TMath::Sqrt(fModelAlphaTeSurface2->Integral( fModelAlphaTeSurface2->FindBin(2470), fModelAlphaTeSurface2->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr);  
  double dAlpha0nbbCuSurface1 = fModelAlphaCuSurface1->Integral( fModelAlphaCuSurface1->FindBin(2470), fModelAlphaCuSurface1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr);
  double dAlpha0nbbCuSurface1Err = TMath::Sqrt(fModelAlphaCuSurface1->Integral( fModelAlphaCuSurface1->FindBin(2470), fModelAlphaCuSurface1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr);
  double dAlpha0nbbCuSurface2 = fModelAlphaCuSurface2->Integral( fModelAlphaCuSurface2->FindBin(2470), fModelAlphaCuSurface2->FindBin(2570), "width" )/(dROIRange*dLivetimeYr);
  double dAlpha0nbbCuSurface2Err = TMath::Sqrt(fModelAlphaTeSurface2->Integral( fModelAlphaTeSurface2->FindBin(2470), fModelAlphaTeSurface2->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr);

  double dGammaClose0nbb = fModelGammaClose->Integral( fModelGammaClose->FindBin(2470), fModelGammaClose->FindBin(2570), "width" )/(dROIRange*dLivetimeYr);
  double dGammaClose0nbbErr = TMath::Sqrt(fModelGammaClose->Integral( fModelGammaClose->FindBin(2470), fModelGammaClose->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr);
  double dGammaMid0nbb = fModelGammaMid->Integral( fModelGammaMid->FindBin(2470), fModelGammaMid->FindBin(2570), "width" )/(dROIRange*dLivetimeYr);
  double dGammaMid0nbbErr = TMath::Sqrt(fModelGammaMid->Integral( fModelGammaMid->FindBin(2470), fModelGammaMid->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr);
  double dGammaFar0nbb = fModelGammaFar->Integral( fModelGammaFar->FindBin(2470), fModelGammaFar->FindBin(2570), "width" )/(dROIRange*dLivetimeYr);
  double dGammaFar0nbbErr = TMath::Sqrt(fModelGammaFar->Integral( fModelGammaFar->FindBin(2470), fModelGammaFar->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr);


  double d2nbbData = fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double d2nbbDataErr = TMath::Sqrt(fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);
  double d2nbbModel = fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double d2nbbModelErr = TMath::Sqrt(fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);
  double d2nbbPDF = fModelTotAdap2NDBDM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr);
  double d2nbbPDFErr = TMath::Sqrt(fModelTotAdap2NDBDM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr);
  double dROIData = fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr);
  double dROIDataErr = TMath::Sqrt(fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr);
  double dROIModel = fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr);
  double dROIModelErr = TMath::Sqrt(fModelTotAdapM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr);
  
  // Output integrals of stuff for limits
  cout << "ROI range: " << fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2470)) << " " << fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2570))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2570)) << " keV" << endl; // 2470 to 2572
  cout << "Integral Data in ROI: " << fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt( fAdapDataHistoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << "Integral Total PDF in ROI: " << fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width") << " +/- " << sqrt( fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << "Integral Total Th-232 PDF in ROI: " << fModelTotAdapthM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt( fModelTotAdapthM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << "Integral Total U-238 PDF in ROI: " << fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << "Integral Total Co PDF in ROI: " << fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << "Integral Total Pb-210 PDF in ROI: " << fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << "Integral Total Po-210 PDF in ROI: " << fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;  
  cout << "Integral Total 0NDBD PDF in ROI: " << fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  cout << endl;
  cout << endl;
  cout << "Integral Data in ROI (counts/keV/y): " << fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " +/- " << sqrt( fAdapDataHistoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << endl;
  cout << "Integral Total PDF in ROI (counts/keV/y): " << fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " +/- " << sqrt( fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << endl;
  cout << "Integral Total Th-232 PDF in ROI (counts/keV/y): " << fModelTotAdapthM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " +/- " << sqrt( fModelTotAdapthM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << endl;
  cout << "Integral Total U-238 PDF in ROI (counts/keV/y): " << fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " +/- " << sqrt(fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << endl;
  cout << "Integral Total Co PDF in ROI (counts/keV/y): " << fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " +/- " << sqrt(fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << endl;
  cout << "Integral Total Pb-210 PDF in ROI (counts/keV/y): " << fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " +/- " << sqrt(fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << endl;
  // cout << "Integral Total Po-210 PDF in ROI (counts/keV): " << fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " +/- " << sqrt(fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << endl;  
  cout << "Integral Total 0NDBD PDF in ROI (counts/keV/y): " << fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " +/- " << sqrt(fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << endl;
  cout << "Data in 2nbb region (c/keV/y): " << d2nbbData << " $\\pm$ " << d2nbbDataErr << endl;  
  // cout << "Number of 2nbb: " << fParameters[0]*dDataIntegralM1 << " +/- " << fParError[0]*dDataIntegralM1 << "\t 2nbb half life: " << (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1) << " +/- " << (fParError[0]/fParameters[0]) * (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1) << endl;
  // cout << "Counts in 2nbb (M1 + M2): " << fModelTotAdap2NDBDM1->Integral(1, fAdapDataHistoM1->FindBin(3000), "width") + fModelTotAdap2NDBDM2->Integral(1, fAdapDataHistoM2->FindBin(3000) , "width")/2 << "\t Half-Life " << (0.69314718056)*(4.726e25 * dLivetimeYr)/(fModelTotAdap2NDBDM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width") + fModelTotAdap2NDBDM2->Integral(1, fAdapDataHistoM2->FindBin(2700) , "width")/2) << endl;
  cout << "Model in 2nbb region (c/keV/y): " << d2nbbModel << " $\\pm$ " << d2nbbModelErr << endl;
  cout << "Model of 2nbb PDF (c/keV/y): " << d2nbbPDF << " $\\pm$ " << d2nbbPDFErr << endl;
  cout << "Model of PDFs in 2nbb region without 2nbb (c/keV/y): " << d2nbbModel - d2nbbPDF << " $\\pm$ " << d2nbbModelErr << endl;
  cout << "Data in 0nbb region (c/keV/y): " << dROIData << " $\\pm$ " << dROIDataErr << endl;
  cout << "Model in 0nbb region (c/keV/y): " << dROIModel << " $\\pm$ " << dROIModelErr << endl;
  cout << "Model/Data in 2nbb region: " << d2nbbModel/d2nbbData << endl;
  cout << "Model/Data in 0nbb region: " << dROIModel/dROIData << endl;
  cout << "Model of PDFs in 2nbb region without 2nbb/Data : " << (d2nbbModel - d2nbbPDF)/d2nbbData << endl;
  cout << endl;
  cout << endl;
  cout << "Crystal Bulk in 2nbb region (c/keV/y): " << dAlpha2nbbBulk << " $\\pm$ " << dAlpha2nbbBulkErr << endl;
  cout << "Crystal Surface (depth ~alpha range) in 2nbb PDF (c/keV/y): " << dAlpha2nbbTeSurface1 << " $\\pm$ " << dAlpha2nbbTeSurface1Err << endl;
  cout << "Crystal Surface (depth << alpha range) 2nbb (c/keV/y): " << dAlpha2nbbTeSurface2 << " $\\pm$ " << dAlpha2nbbTeSurface2Err << endl;
  cout << "Copper Surface (depth ~alpha range) in 2nbb PDF (c/keV/y): " << dAlpha2nbbCuSurface1 << " $\\pm$ " << dAlpha2nbbCuSurface1Err << endl;
  cout << "Copper Surface (depth << alpha range) 2nbb (c/keV/y): " << dAlpha2nbbCuSurface2 << " $\\pm$ " << dAlpha2nbbCuSurface2Err << endl;
  cout << endl;
  cout << "Crystal Bulk in 0nbb region (c/keV/y): " << dAlpha0nbbBulk << " $\\pm$ " << dAlpha0nbbBulkErr << endl;
  cout << "Crystal Surface (depth ~alpha range) in 0nbb PDF (c/keV/y): " << dAlpha0nbbTeSurface1 << " $\\pm$ " << dAlpha0nbbTeSurface1Err << endl;
  cout << "Crystal Surface (depth << alpha range) 0nbb (c/keV/y): " << dAlpha0nbbTeSurface2 << " $\\pm$ " << dAlpha0nbbTeSurface2Err << endl;
  cout << "Copper Surface (depth ~alpha range) in 0nbb PDF (c/keV/y): " << dAlpha0nbbCuSurface1 << " $\\pm$ " << dAlpha0nbbCuSurface1Err << endl;
  cout << "Copper Surface (depth << alpha range) 0nbb (c/keV/y): " << dAlpha0nbbCuSurface2 << " $\\pm$ " << dAlpha0nbbCuSurface2Err << endl;  
  cout << endl;
  cout << "Gamma Close in 2nbb region (c/keV/y): " << dGammaClose2nbb << " $\\pm$ " << dGammaClose2nbbErr << endl;
  cout << "Gamma Mid in 2nbb region (c/keV/y): " << dGammaMid2nbb << " $\\pm$ " << dGammaMid2nbbErr << endl;
  cout << "Gamma Far in 2nbb region (c/keV/y): " << dGammaFar2nbb << " $\\pm$ " << dGammaFar2nbbErr << endl;
  cout << endl;
  cout << "Gamma Close in 0nbb region (c/keV/y): " << dGammaClose0nbb << " $\\pm$ " << dGammaClose0nbbErr << endl;
  cout << "Gamma Mid in 0nbb region (c/keV/y): " << dGammaMid0nbb << " $\\pm$ " << dGammaMid0nbbErr << endl;
  cout << "Gamma Far in 0nbb region (c/keV/y): " << dGammaFar0nbb << " $\\pm$ " << dGammaFar0nbbErr << endl;
  cout << endl;  
  cout << endl;
  cout << "Residual RMS (Tot): " << dResidualRMSTot << endl;
  cout << "Residual RMS (M1): " << dResidualRMSM1 << "\t" << "Residual RMS (M2): " << dResidualRMSM2 << "\t Residual RMS (M2Sum): "  << dResidualRMSM2Sum << endl;

  dChiSquare = GetChiSquareAdaptive();

  cout << "Total number of calls = " << dNumCalls << "\t" << "ChiSq/NDF = " << dChiSquare/dNDF << endl; // for M1 and M2
  cout << "ChiSq = " << dChiSquare << "\t" << "NDF = " << dNDF << endl;
  cout << "Probability = " << TMath::Prob(dChiSquare, dNDF ) << endl;

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


  // This gets the number of counts corrected for detector efficiency
  for(int i = 0; i < TBkgModel::dNParam; i++)
  {
    fParActivity[i] = fParameters[i]*dDataIntegralM1/fParEfficiencyM1[i];
    fParActivityErr[i] = fParError[i]*dDataIntegralM1/fParEfficiencyM1[i];
  }


  // Correlation Matrix section
  TMatrixT<double> mCorrMatrix;
  TMatrixT<double> mCorrMatrixInverse;
  mCorrMatrix.ResizeTo(TBkgModel::dNParam, TBkgModel::dNParam);
  mCorrMatrixInverse.ResizeTo(TBkgModel::dNParam, TBkgModel::dNParam);
  minuit->mnemat(mCorrMatrix.GetMatrixArray(), TBkgModel::dNParam);

  for(int i = mCorrMatrix.GetRowLwb(); i <= mCorrMatrix.GetRowUpb(); i++)
    for(int j = mCorrMatrix.GetColLwb(); j <= mCorrMatrix.GetColUpb(); j++)
      mCorrMatrix(i,j) = mCorrMatrix(i,j)/(fParError[i]*fParError[j]);

  for(int i = mCorrMatrix.GetRowLwb(); i <= mCorrMatrix.GetRowUpb(); i++)
    for(int j = mCorrMatrix.GetColLwb(); j <= mCorrMatrix.GetColUpb(); j++)
      mCorrMatrixInverse(i,j) = mCorrMatrix(TBkgModel::dNParam-i-1, j); 

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
  // TH2D *h2Dummy = new TH2D("h2Dummy","", TBkgModel::dNParam, 0, TBkgModel::dNParam, TBkgModel::dNParam, 0, TBkgModel::dNParam);
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
  for(int i=0; i < TBkgModel::dNParam; i++)
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
    if( fAdapDataHistoM1->GetBinCenter(i) >= 3150 && fAdapDataHistoM1->GetBinCenter(i) <= 3400)continue;

    // Skipping unknown peaks
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 800 && fAdapDataHistoM1->GetBinCenter(i) <= 808)continue;
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 1060 && fAdapDataHistoM1->GetBinCenter(i) <= 1068)continue; 

    // if( fAdapDataHistoM1->GetBinCenter(i) >= 506 && fAdapDataHistoM1->GetBinCenter(i) <= 515)continue;
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 579 && fAdapDataHistoM1->GetBinCenter(i) <= 589)continue;
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 605 && fAdapDataHistoM1->GetBinCenter(i) <= 615)continue;
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 906 && fAdapDataHistoM1->GetBinCenter(i) <= 917)continue;
    // if( fAdapDataHistoM1->GetBinCenter(i) >= 1450 && fAdapDataHistoM1->GetBinCenter(i) <= 1475)continue;  
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
    if( fAdapDataHistoM2->GetBinCenter(i) >= 3150 && fAdapDataHistoM2->GetBinCenter(i) <= 3400)continue;

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
    if( fAdapDataHistoM2Sum->GetBinCenter(i) >= 3150 && fAdapDataHistoM2Sum->GetBinCenter(i) <= 3400)continue;

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
  cProgress->Divide(1, 3);
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

  cProgress->cd(3);
  hChiSquaredProgressM2Sum->SetTitle("Increase in #chi^{2} (M2Sum)");  
  hChiSquaredProgressM2Sum->GetXaxis()->SetTitle("Energy (keV)");
  hChiSquaredProgressM2Sum->GetYaxis()->SetTitle("#chi^{2}");  
  hChiSquaredProgressM2Sum->Draw();

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
  cadap1->SaveAs(Form("%s/FitResults/Test/FitM1_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  cadap2->SaveAs(Form("%s/FitResults/Test/FitM2_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  cResidual1->SaveAs(Form("%s/FitResults/Test/FitM1Residual_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  cResidual2->SaveAs(Form("%s/FitResults/Test/FitM2Residual_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  cres1->SaveAs(Form("%s/FitResults/Test/FitResidualDist_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  cMatrix->SaveAs(Form("%s/FitResults/Test/FitCovMatrix_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
 
  cProgress->SaveAs(Form("%s/FitResults/Test/ChiSquareProgress_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));

  // Save histograms to file
  fSaveResult = new TFile(Form("%s/FitResults/Test/FitResult_%d_%d.root", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime()), "RECREATE");
  fSaveResult->Add(fAdapDataHistoM1);
  fSaveResult->Add(fAdapDataHistoM2);
  fSaveResult->Add(fModelTotAdapM1);
  fSaveResult->Add(fModelTotAdapM2);

  // Scale histograms for saving
  // hAdapTeO20nuM1->Scale( dDataIntegralM1*fParameters[0]);
  hAdapTeO22nuM1->Scale( dDataIntegralM1*fParameters[0]);
  hAdapCuBox_CuFrameco60M1->Scale( dDataIntegralM1*fParameters[1]);
  hAdapTeO2th232onlyM1->Scale( dDataIntegralM1*fParameters[2]);
  hAdapTeO2th230onlyM1->Scale( dDataIntegralM1*fParameters[3]);
  hAdapTeO2Sxth232M1_001->Scale( dDataIntegralM1*fParameters[4]);
  hAdapTeO2Sxra228pb208M1_001->Scale( dDataIntegralM1*fParameters[5]);
  hAdapTeO2Sxu238th230M1_001->Scale( dDataIntegralM1*fParameters[6]);
  hAdapTeO2Sxth230onlyM1_001->Scale( dDataIntegralM1*fParameters[7]);
  hAdapTeO2Sxra226pb210M1_001->Scale( dDataIntegralM1*fParameters[8]);
  hAdapTeO2Sxpb210M1_1->Scale( dDataIntegralM1*fParameters[9]);
  hAdapTeO2Sxpb210M1_001->Scale( dDataIntegralM1*fParameters[10]);
  hAdapCuBox_CuFrameth232M1_10->Scale( dDataIntegralM1*fParameters[11]);
  hAdapCuBox_CuFrameu238M1_10->Scale( dDataIntegralM1*fParameters[12]);
  hAdapCuBox_CuFramepb210M1_01->Scale( dDataIntegralM1*fParameters[13]);
  hAdapCuBox_CuFramepb210M1_001->Scale( dDataIntegralM1*fParameters[14]);
  hAdapPbRomk40M1->Scale( dDataIntegralM1*fParameters[15]);
  hAdapOVCth232M1->Scale( dDataIntegralM1*fParameters[16]);
  hAdapOVCu238M1->Scale( dDataIntegralM1*fParameters[17]);
  hAdapOVCco60M1->Scale( dDataIntegralM1*fParameters[18]);
  hAdapOVCk40M1->Scale( dDataIntegralM1*fParameters[19]);
  hAdapExtPbbi210M1->Scale( dDataIntegralM1*fParameters[20]);
  hAdapCuBox_CuFrameth232M1->Scale( dDataIntegralM1*fParameters[21]);
  hAdapCuBox_CuFrameu238M1->Scale( dDataIntegralM1*fParameters[22]);
  hAdapPbRomcs137M1->Scale( dDataIntegralM1*fParameters[23]);

  // hAdapTeO20nuM2->Scale( dDataIntegralM2*fParameters[0]);
  hAdapTeO22nuM2->Scale( dDataIntegralM2*fParameters[0]);
  hAdapCuBox_CuFrameco60M2->Scale( dDataIntegralM2*fParameters[1]);
  hAdapTeO2th232onlyM2->Scale( dDataIntegralM2*fParameters[2]);
  hAdapTeO2th230onlyM2->Scale( dDataIntegralM2*fParameters[3]);
  hAdapTeO2Sxth232M2_001->Scale( dDataIntegralM2*fParameters[4]);
  hAdapTeO2Sxra228pb208M2_001->Scale( dDataIntegralM2*fParameters[5]);
  hAdapTeO2Sxu238th230M2_001->Scale( dDataIntegralM2*fParameters[6]);
  hAdapTeO2Sxth230onlyM2_001->Scale( dDataIntegralM2*fParameters[7]);
  hAdapTeO2Sxra226pb210M2_001->Scale( dDataIntegralM2*fParameters[8]);
  hAdapTeO2Sxpb210M2_1->Scale( dDataIntegralM2*fParameters[9]);
  hAdapTeO2Sxpb210M2_001->Scale( dDataIntegralM2*fParameters[10]);
  hAdapCuBox_CuFrameth232M2_10->Scale( dDataIntegralM2*fParameters[11]);
  hAdapCuBox_CuFrameu238M2_10->Scale( dDataIntegralM2*fParameters[12]);
  hAdapCuBox_CuFramepb210M2_01->Scale( dDataIntegralM2*fParameters[13]);
  hAdapCuBox_CuFramepb210M2_001->Scale( dDataIntegralM2*fParameters[14]);
  hAdapPbRomk40M2->Scale( dDataIntegralM2*fParameters[15]);
  hAdapOVCth232M2->Scale( dDataIntegralM2*fParameters[16]);
  hAdapOVCu238M2->Scale( dDataIntegralM2*fParameters[17]);
  hAdapOVCco60M2->Scale( dDataIntegralM2*fParameters[18]);
  hAdapOVCk40M2->Scale( dDataIntegralM2*fParameters[19]);
  hAdapExtPbbi210M2->Scale( dDataIntegralM2*fParameters[20]);
  hAdapCuBox_CuFrameth232M2->Scale( dDataIntegralM2*fParameters[21]);
  hAdapCuBox_CuFrameu238M2->Scale( dDataIntegralM2*fParameters[22]);
  hAdapPbRomcs137M2->Scale( dDataIntegralM2*fParameters[23]);

  fSaveResult->Add(hAdapTeO20nuM1);
  fSaveResult->Add(hAdapTeO22nuM1);
  // fSaveResult->Add(hAdapTeO2co60M1);
  // fSaveResult->Add(hAdapTeO2k40M1);
  fSaveResult->Add(hAdapCuBox_CuFrameco60M1);
  fSaveResult->Add(hAdapTeO2th232onlyM1);
  fSaveResult->Add(hAdapTeO2th230onlyM1);
  fSaveResult->Add(hAdapTeO2Sxth232M1_001);
  fSaveResult->Add(hAdapTeO2Sxra228pb208M1_001);
  fSaveResult->Add(hAdapTeO2Sxu238th230M1_001);
  fSaveResult->Add(hAdapTeO2Sxth230onlyM1_001);
  fSaveResult->Add(hAdapTeO2Sxra226pb210M1_001);
  fSaveResult->Add(hAdapTeO2Sxpb210M1_1);
  fSaveResult->Add(hAdapTeO2Sxpb210M1_001);
  fSaveResult->Add(hAdapCuBox_CuFrameth232M1_10);
  fSaveResult->Add(hAdapCuBox_CuFrameu238M1_10);
  fSaveResult->Add(hAdapCuBox_CuFramepb210M1_01);
  fSaveResult->Add(hAdapCuBox_CuFramepb210M1_001);
  fSaveResult->Add(hAdapPbRomk40M1);
  fSaveResult->Add(hAdapOVCth232M1);
  fSaveResult->Add(hAdapOVCu238M1);
  fSaveResult->Add(hAdapOVCco60M1);
  fSaveResult->Add(hAdapOVCk40M1);
  fSaveResult->Add(hAdapExtPbbi210M1);
  fSaveResult->Add(hAdapCuBox_CuFrameth232M1);
  fSaveResult->Add(hAdapCuBox_CuFrameu238M1);

  fSaveResult->Add(hAdapTeO20nuM2);
  fSaveResult->Add(hAdapTeO22nuM2);
  // fSaveResult->Add(hAdapTeO2co60M2);
  // fSaveResult->Add(hAdapTeO2k40M2);
  fSaveResult->Add(hAdapCuBox_CuFrameco60M2);  
  fSaveResult->Add(hAdapTeO2th232onlyM2);
  fSaveResult->Add(hAdapTeO2th230onlyM2);
  fSaveResult->Add(hAdapTeO2Sxth232M2_001);
  fSaveResult->Add(hAdapTeO2Sxra228pb208M2_001);
  fSaveResult->Add(hAdapTeO2Sxu238th230M2_001);
  fSaveResult->Add(hAdapTeO2Sxth230onlyM2_001);
  fSaveResult->Add(hAdapTeO2Sxra226pb210M2_001);
  fSaveResult->Add(hAdapTeO2Sxpb210M2_1);
  fSaveResult->Add(hAdapTeO2Sxpb210M2_001);
  fSaveResult->Add(hAdapCuBox_CuFrameth232M2_10);
  fSaveResult->Add(hAdapCuBox_CuFrameu238M2_10);
  fSaveResult->Add(hAdapCuBox_CuFramepb210M2_01);
  fSaveResult->Add(hAdapCuBox_CuFramepb210M2_001);
  fSaveResult->Add(hAdapPbRomk40M2);
  fSaveResult->Add(hAdapOVCth232M2);
  fSaveResult->Add(hAdapOVCu238M2);
  fSaveResult->Add(hAdapOVCco60M2);
  fSaveResult->Add(hAdapOVCk40M2);
  fSaveResult->Add(hAdapExtPbbi210M2);
  fSaveResult->Add(hAdapCuBox_CuFrameth232M2);
  fSaveResult->Add(hAdapCuBox_CuFrameu238M2);

  fSaveResult->Add(&mCorrMatrix);

  fSaveResult->Write(); 

  // cadap1->SaveAs(Form("%s/FitResults/CovMatrix/FitM1_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cadap2->SaveAs(Form("%s/FitResults/CovMatrix/FitM2_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cResidual1->SaveAs(Form("%s/FitResults/CovMatrix/FitM1Residual_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cResidual2->SaveAs(Form("%s/FitResults/CovMatrix/FitM2Residual_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cres1->SaveAs(Form("%s/FitResults/CovMatrix/FitResidualDist_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cMatrix->SaveAs(Form("%s/FitResults/CovMatrix/FitCovMatrix_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));

  LatexResultTable(0);
  } // end bSave

  // Kernal Convolution
  // TH1D *hKernalConvM1 = new TH1D("hKernalConvM1", "", dAdaptiveBinsM1, dAdaptiveArrayM1);
  // Kernal(fModelTotAdapM1, hKernalConvM1);

  // TCanvas *cKernal1 = new TCanvas("cKernal1", "cKernal1", 1200, 1200);
  // hKernalConvM1->Draw();



  return true;
}

// Txt file with useful stuff
void TBkgModel::LatexResultTable(double fValue)
{

  OutFile.open(Form("%s/FitResults/Test/FitOutputTable_%d_%d_%d.txt", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop ));
  OutFile << "Initial Value of 2nbb: " << fValue << endl;
  OutFile << "Fit Range: " << dFitMin << " to " << dFitMax << endl;
  OutFile << "Base binning: " << dBinBase << endl;
  OutFile << "Events in background spectrum (M1): " << dDataIntegralM1 << endl;
  OutFile << "Events in background spectrum (M2): " << dDataIntegralM2 << endl;
  OutFile << "Events in background spectrum (M2Sum): " << dDataIntegralM2Sum << endl;
  OutFile << "Livetime of background: " << dLivetimeYr << endl;
  OutFile << "Total number of calls = " << dNumCalls << "\t" << "ChiSq/NDF = " << dChiSquare/dNDF << endl; // for M1 and M2
  OutFile << "ChiSq = " << dChiSquare << "\t" << "NDF = " << dNDF << endl;
  OutFile << "Probability = " << TMath::Prob(dChiSquare, dNDF ) << endl;

  OutFile << endl;
  OutFile << endl;

  double dROIRange = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2570))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2570)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2470)); 
  double d2nbbRange = fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2000))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2000)) - fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(500)); 
  // Output integrals of stuff for limits
  OutFile << "ROI range: " << fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2470)) << " " << fAdapDataHistoM1->GetBinLowEdge(fAdapDataHistoM1->FindBin(2570))+fAdapDataHistoM1->GetBinWidth(fAdapDataHistoM1->FindBin(2570)) << " keV" << endl; // 2470 to 2572
  OutFile << "Integral Data in ROI: " << fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt( fAdapDataHistoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total PDF in ROI: " << fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width") << " +/- " << sqrt( fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total Th-232 PDF in ROI: " << fModelTotAdapthM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt( fModelTotAdapthM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total U-238 PDF in ROI: " << fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total Co PDF in ROI: " << fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total Pb-210 PDF in ROI: " << fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total Po-210 PDF in ROI: " << fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;  
  // OutFile << "Integral Total 2NDBD PDF in ROI: " << fModelTotAdap2NDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdap2NDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << "Integral Total 0NDBD PDF in ROI: " << fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ) << " +/- " << sqrt(fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )) << endl;
  OutFile << endl;
  OutFile << "Integral Data in ROI (counts/keV): " << fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt( fAdapDataHistoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Integral Total PDF in ROI (counts/keV): " << fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt( fModelTotAdapM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Integral Total Th-232 PDF in ROI (counts/keV): " << fModelTotAdapthM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt( fModelTotAdapthM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Integral Total U-238 PDF in ROI (counts/keV): " << fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapuM1->Integral( fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Integral Total Co PDF in ROI (counts/keV): " << fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapcoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Integral Total Pb-210 PDF in ROI (counts/keV): " << fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapSpbM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  // OutFile << "Integral Total Po-210 PDF in ROI (counts/keV): " << fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapSpoM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;  
  // OutFile << "Integral Total 2NDBD PDF in ROI (counts/keV): " << fModelTotAdap2NDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdap2NDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  OutFile << "Integral Total 0NDBD PDF in ROI (counts/keV): " << fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" )/dROIRange << " +/- " << sqrt(fModelTotAdapNDBDM1->Integral(fAdapDataHistoM1->FindBin(2470),fAdapDataHistoM1->FindBin(2570), "width" ))/dROIRange << endl;
  // OutFile << "Number of 2nbb: " << fParameters[0]*dDataIntegralM1 << " +/- " << fParError[0]*dDataIntegralM1 << "\t 2nbb half life: " << (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1) << " +/- " << fParError[0]/fParameters[0] * (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1) << endl;
  // OutFile << "Counts in 2nbb (M1 + M2): " << fModelTotAdap2NDBDM1->Integral(1, fAdapDataHistoM1->FindBin(3000), "width") + fModelTotAdap2NDBDM2->Integral(1, fAdapDataHistoM2->FindBin(3000) , "width")/2 << "\t Half-Life " << (0.69314718056)*(4.726e25 * dLivetimeYr)/(fModelTotAdap2NDBDM1->Integral(1, fAdapDataHistoM1->FindBin(2700), "width") + fModelTotAdap2NDBDM2->Integral(1, fAdapDataHistoM2->FindBin(2700) , "width")/2) << endl;
  
  OutFile << "Residual RMS (Tot): " << dResidualRMSTot << endl;
  OutFile << "Residual RMS (M1): " << dResidualRMSM1 << "\t" << "Residual RMS (M2): " << dResidualRMSM2 << endl;
  OutFile << endl;
  OutFile << endl;

  // Outputs table of best fit values
  OutFile << "// Name - BestFit - BestFit Err - Integral - Integral Err" << endl;
  for(int i = 0; i < TBkgModel::dNParam; i++)
  {
    OutFile << minuit->fCpnam[i] << Form(" & %.4e$\\pm$%.4e \t\t %.4e$\\pm$%.4e \\\\ ", fParameters[i], fParError[i], dDataIntegralM1*fParameters[i], dDataIntegralM1*fParError[i] ) << endl;
  }

  OutFile << endl;
  OutFile << endl;

  // Outputs table of activities of each cryostat element
  OutFile <<  "// Name - Activity(events) - Activity (Bq/kg) - Activity (Bq/cm2)" << endl;
  for(int i = 0; i < TBkgModel::dNParam; i++)
  {
    OutFile << minuit->fCpnam[i] << Form(" & %.4e \\pm %.4e", fParActivity[i], fParActivityErr[i]) << "\t \t" <<  Form("& %.4e \\pm %.4e", fParActivity[i]/fParMass[i]/dLivetime, fParActivityErr[i]/fParMass[i]/dLivetime) << "\t \t" << Form("& %.4e \\pm %.4e \\\\", fParActivity[i]/fParSurfaceArea[i]/dLivetime, fParActivityErr[i]/fParSurfaceArea[i]/dLivetime) << endl;
  }



  OutFile << endl;
  OutFile << endl;
  OutFile << endl;
  
  OutFile << "Contributions to ROI (M1) (counts/keV/y) +/- Err (counts/keV/y)" << endl;
  OutFile << minuit->fCpnam[0] << "& " << hAdapTeO22nuM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO22nuM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[1] << "& " << hAdapCuBox_CuFrameco60M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFrameco60M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[2] << "& " << hAdapTeO2th232onlyM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2th232onlyM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[3] << "& " << hAdapTeO2th230onlyM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2th230onlyM1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[4] << "& " << hAdapTeO2Sxth232M1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxth232M1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[5] << "& " << hAdapTeO2Sxra228pb208M1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxra228pb208M1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[6] << "& " << hAdapTeO2Sxu238th230M1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxu238th230M1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[7] << "& " << hAdapTeO2Sxth230onlyM1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxth230onlyM1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[8] << "& " << hAdapTeO2Sxra226pb210M1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxra226pb210M1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[9] << "& " << hAdapTeO2Sxpb210M1_1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxpb210M1_1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[10] << "& " << hAdapTeO2Sxpb210M1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxpb210M1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[11] << "& " << hAdapCuBox_CuFrameth232M1_10->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFrameth232M1_10->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[12] << "& " << hAdapCuBox_CuFrameu238M1_10->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFrameu238M1_10->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[13] << "& " << hAdapCuBox_CuFramepb210M1_01->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFramepb210M1_01->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[14] << "& " << hAdapCuBox_CuFramepb210M1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFramepb210M1_001->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[15] << "& " << hAdapPbRomk40M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapPbRomk40M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[16] << "& " << hAdapOVCth232M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapOVCth232M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[17] << "& " << hAdapOVCu238M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapOVCu238M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[18] << "& " << hAdapOVCco60M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapOVCco60M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[19] << "& " << hAdapOVCk40M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapOVCk40M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[20] << "& " << hAdapExtPbbi210M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapExtPbbi210M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[21] << "& " << hAdapCuBox_CuFrameth232M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFrameth232M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[22] << "& " << hAdapCuBox_CuFrameu238M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFrameu238M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[23] << "& " << hAdapPbRomcs137M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" )/(dROIRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapPbRomcs137M1->Integral( fAdapDataHistoM1->FindBin(2470), fAdapDataHistoM1->FindBin(2570), "width" ))/(dROIRange*dLivetimeYr) << " \\\\" << endl;

  OutFile << endl;
  OutFile << endl;
  OutFile << endl;

  OutFile << "Contributions to 2nbb (M1) (counts/keV/y) +/- Err (counts/keV/y)" << endl;
  OutFile << "Number in 2nbb region" <<fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(fAdapDataHistoM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << endl;
  OutFile << minuit->fCpnam[0] << " & " << hAdapTeO22nuM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO22nuM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[1] << " & " << hAdapCuBox_CuFrameco60M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFrameco60M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[2] << " & " << hAdapTeO2th232onlyM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2th232onlyM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[3] << " & " << hAdapTeO2th230onlyM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2th230onlyM1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[4] << " & " << hAdapTeO2Sxth232M1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxth232M1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[5] << " & " << hAdapTeO2Sxra228pb208M1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxra228pb208M1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[6] << " & " << hAdapTeO2Sxu238th230M1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxu238th230M1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[7] << " & " << hAdapTeO2Sxth230onlyM1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxth230onlyM1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[8] << " & " << hAdapTeO2Sxra226pb210M1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxra226pb210M1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[9] << " & " << hAdapTeO2Sxpb210M1_1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxpb210M1_1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[10] << " & " << hAdapTeO2Sxpb210M1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapTeO2Sxpb210M1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[11] << " & " << hAdapCuBox_CuFrameth232M1_10->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFrameth232M1_10->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[12] << " & " << hAdapCuBox_CuFrameu238M1_10->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFrameu238M1_10->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[13] << " & " << hAdapCuBox_CuFramepb210M1_01->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFramepb210M1_01->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[14] << " & " << hAdapCuBox_CuFramepb210M1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFramepb210M1_001->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[15] << " & " << hAdapPbRomk40M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapPbRomk40M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[16] << " & " << hAdapOVCth232M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapOVCth232M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[17] << " & " << hAdapOVCu238M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapOVCu238M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[18] << " & " << hAdapOVCco60M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapOVCco60M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[19] << " & " << hAdapOVCk40M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapOVCk40M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[20] << " & " << hAdapExtPbbi210M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapExtPbbi210M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[21] << " & " << hAdapCuBox_CuFrameth232M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFrameth232M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[22] << " & " << hAdapCuBox_CuFrameu238M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapCuBox_CuFrameu238M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;
  OutFile << minuit->fCpnam[23] << " & " << hAdapPbRomcs137M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" )/(d2nbbRange*dLivetimeYr) << " $\\pm$ " << TMath::Sqrt(hAdapPbRomcs137M1->Integral( fAdapDataHistoM1->FindBin(500), fAdapDataHistoM1->FindBin(2000), "width" ))/(d2nbbRange*dLivetimeYr) << " \\\\" << endl;

  OutFile << endl;
  OutFile << endl;
  // for(int i = dFitMinBinM1; i < dFitMaxBinM1; i++)
  // {
  //   OutFile << fAdapDataHistoM1->GetBinContent(i) << "\t" << fAdapDataHistoM1->GetBinWidth(i) << "\t" << fModelTotAdapM1->GetBinContent(i) << "\t" << fModelTotAdapM1->GetBinWidth(i) << endl;
  // }

  OutFile << endl;
  OutFile << endl;
  OutFile << endl;
  OutFile << endl;  

  OutFile.close();
}

// Adds a random percentage of events of 2nbb into spectrum
void TBkgModel::SanityCheck()
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
void TBkgModel::ProfileNLL(double fBestFitInit, double fBestFitChiSq)
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
    DoTheFitAdaptive(*iter, 0);
    // LatexResultTable(*iter);
    cout << "delta ChiSq = " << dChiSquare - dBestChiSq << endl; // Needs to be entered, otherwise just 0
    OutPNLL << Form("dX.push_back(%f); dT.push_back(%f);", dChiSquare-dBestChiSq, (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[0]*dDataIntegralM1) ) << endl;

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

void TBkgModel::ProfileNLL2D(double fBestFitInit, double fBestFitInit2, double fBestFitChiSq)
{
  dBestChiSq = fBestFitChiSq; // Chi-Squared from best fit (for ProfileNLL calculation)
  // Do the fit now if no other tests are needed 
  nLoop = 0;
  for(int i = -5; i < 5; i++)
  {
    fInitValues.push_back(fBestFitInit + fBestFitInit/100*i);
  }
  for(int j = -5; j < 5; j++)
  {
    fInitValues2.push_back(fBestFitInit2 + fBestFitInit2/100*j);
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
    OutPNLL << Form("dX.push_back(%f); dY.push_back(%.10f); dT.push_back(%f);", dChiSquare-dBestChiSq, *iter2, (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[1]*dDataIntegralM1) ) << endl;

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
  OutPNLL << "g1->GetYaxis()->SetTitle(\"External Bi-207\");" << endl;
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


void TBkgModel::ToyFit(int fNumFits)
{
    OutToy.open(Form("%s/FitResults/ToyMC/Toy_%d.C", dSaveDir.c_str(), tTime->GetDate() ));

    OutToy << "{" << endl;
    OutToy << endl;

    OutToy << "hPullDist = new TH1D(\"hPullDist\", \"Pull Distribution\", 20, -5, 5);" << endl;
    OutToy << "hToy2nbbRate = new TH1D(\"hToy2nbbRate\", \"Toy Monte Carlo half-life fit values\", 100, 6.7e+20, 6.9e+20);" << endl;
    OutToy << "hToy2nbbError = new TH1D(\"hToy2nbbError\", \"Toy Monte Carlo half-life error values\", 50, 1.e+19, 1.2e+19);" << endl;
    OutToy << "hChiSquared = new TH1D(\"hChiSquared\", \"Distribution of Chi-Squared values\", 60, 0, 20);" << endl;


    cout << "Number of Loops " << fNumFits << endl;
    // Number of toy fits
    for(int i = 1; i <= 10; i++)
    {
      cout << "Toy: " << i << endl;
      fAdapDataHistoM1->Reset();
      fAdapDataHistoM2->Reset();
      
      fToyData = new TFile(Form("%s/Toy/ToyTest_p%d.root", dMCDir.c_str(), i));
      fAdapDataHistoM1 = (TH1D*)fToyData->Get("fAdapDataHistoM1");
      fAdapDataHistoM2 = (TH1D*)fToyData->Get("fAdapDataHistoM2");
    
    // fAdapDataHistoM1->Scale(1.41872984571782959e+05/fAdapDataHistoM1->Integral("width"));
    // fAdapDataHistoM2->Scale(2.66538587693367845e+04/fAdapDataHistoM2->Integral("width"));
      fAdapDataHistoM1->Scale(274875./fAdapDataHistoM1->Integral());
      fAdapDataHistoM2->Scale(97635.6/fAdapDataHistoM2->Integral());

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
      OutToy << Form("hToy2nbbRate->Fill(%.5e);", (0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[1]*dDataIntegralM1) ) << endl;
      OutToy << Form("hToy2nbbError->Fill(%.5e);", fParError[1]/fParameters[1]*(0.69314718056)*(4.726e25*dLivetimeYr)/(fParameters[1]*dDataIntegralM1) ) << endl;
      OutToy << Form("hPullDist->Fill(%5e);", ((0.69314718056)*(4.726e25 * dLivetimeYr)/(fParameters[1]*dDataIntegralM1) - 6.80668e+20)/(fParError[1]/fParameters[1]*(0.69314718056)*(4.726e25*dLivetimeYr)/(fParameters[1]*dDataIntegralM1))  ) << endl;
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

TH1D *TBkgModel::Kernal(TH1D *hMC, TH1D *hSMC)
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
