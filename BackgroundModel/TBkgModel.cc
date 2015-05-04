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
#include "TF1.h"
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
  for(int i = 0; i < Obj->ShowNParameters(); i++)
  {
    Obj->SetParameters(i, x[i]);
  }

  // Implement a method in your class that calculates the quantity you want to minimize, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimize fval
  Obj->UpdateModelAdaptive();
  fval = Obj->GetChiSquareAdaptive();
}

TBkgModel::TBkgModel()
{}

// Probably needs to be done...
TBkgModel::~TBkgModel()
{}

void TBkgModel::CorrectForEfficiency()
{
  // Apply efficiency
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
}

// Shows total number of parameters
int TBkgModel::ShowNParameters()
{
  return dNParam;
}

// Initialize parameters
void TBkgModel::GenerateParameters()
{
  
  // Initialization (Name, Index, Initial Value, Min Limit, Max Limit, pointer to histogram)

  // M1
  BkgParM1[0] = new TBkgModelParameter("TeO2 2nbb", 0, 0, 1E-7, 0.0, 1.0, hAdapTeO22nuM1, hAdapTeO22nuM2, hAdapTeO22nuM2Sum);
  BkgParM1[1] = new TBkgModelParameter("CuBox + CuFrame co60", 1, 0, 1E-7, 0.0, 1.0, hAdapCuBox_CuFrameco60M1, hAdapCuBox_CuFrameco60M2, hAdapCuBox_CuFrameco60M2Sum);
  BkgParM1[2] = new TBkgModelParameter("TeO2 th232 only", 2, 0, 1E-7, 0, 1.0, hAdapTeO2th232onlyM1, hAdapTeO2th232onlyM2, hAdapTeO2th232onlyM2Sum);
  BkgParM1[3] = new TBkgModelParameter("TeO2 th230 only", 3, 3.07203e-04, 1E-7, 0, 1.0, hAdapTeO2th230onlyM1, hAdapTeO2th230onlyM2, hAdapTeO2th230onlyM2Sum);  
  BkgParM1[4] = new TBkgModelParameter("TeO2 Sx th232 only 0.01", 4, 1.43894e-04, 1E-7, 0, 1.0, hAdapTeO2Sxth232onlyM1_001, hAdapTeO2Sxth232onlyM2_001, hAdapTeO2Sxth232onlyM2Sum_001);
  BkgParM1[5] = new TBkgModelParameter("TeO2 Sx ra228 to pb208 0.01", 5, 3.07208e-03, 1E-7, 0, 1.0, hAdapTeO2Sxra228pb208M1_001, hAdapTeO2Sxra228pb208M2_001, hAdapTeO2Sxra228pb208M2Sum_001);
  BkgParM1[6] = new TBkgModelParameter("TeO2 Sx u238 to th230 0.01", 6, 1.71025e-03, 1E-7, 0, 1.0, hAdapTeO2Sxu238th230M1_001, hAdapTeO2Sxu238th230M2_001, hAdapTeO2Sxu238th230M2Sum_001);
  BkgParM1[7] = new TBkgModelParameter("TeO2 Sx th230 only 0.01", 7, 7.23534e-04, 1E-7, 0, 1.0, hAdapTeO2Sxth230onlyM1_001, hAdapTeO2Sxth230onlyM2_001, hAdapTeO2Sxth230onlyM2Sum_001);
  BkgParM1[8] = new TBkgModelParameter("TeO2 Sx ra226 to pb210 0.01", 8, 2.98297e-03, 1E-7, 0, 1.0, hAdapTeO2Sxra226pb210M1_001, hAdapTeO2Sxra226pb210M2_001, hAdapTeO2Sxra226pb210M2Sum_001);
  BkgParM1[9] = new TBkgModelParameter("TeO2 Sx pb210 1", 9, 5.30974e-03, 1E-7, 0, 1.0, hAdapTeO2Sxpb210M1_1, hAdapTeO2Sxpb210M2_1, hAdapTeO2Sxpb210M2Sum_1);
  BkgParM1[10] = new TBkgModelParameter("TeO2 Sx pb210 0.01", 10, 4.11621e-02, 1E-7, 0, 1.0, hAdapTeO2Sxpb210M1_001, hAdapTeO2Sxpb210M2_001, hAdapTeO2Sxpb210M2Sum_001);
  BkgParM1[11] = new TBkgModelParameter("CuBox + CuFrame Sx th232 10", 11, 3.53539e-03, 1E-7, 0, 1.0, hAdapCuBox_CuFrameth232M1_10, hAdapCuBox_CuFrameth232M2_10, hAdapCuBox_CuFrameth232M2Sum_10);
  BkgParM1[12] = new TBkgModelParameter("CuBox + CuFrame Sx u238 10 ", 12, 5.80166e-03, 1E-7, 0, 1.0, hAdapCuBox_CuFrameu238M1_10, hAdapCuBox_CuFrameu238M2_10, hAdapCuBox_CuFrameu238M2Sum_10);
  BkgParM1[13] = new TBkgModelParameter("CuBox + CuFrame Sx pb210 0.1", 13, 5.91139e-03, 1E-7, 0, 1.0, hAdapCuBox_CuFramepb210M1_01, hAdapCuBox_CuFramepb210M2_01, hAdapCuBox_CuFramepb210M2Sum_01);
  BkgParM1[14] = new TBkgModelParameter("CuBox + CuFrame Sx pb210 0.01", 14, 1.79798e-02, 1E-7, 0, 1.0, hAdapCuBox_CuFramepb210M1_001, hAdapCuBox_CuFramepb210M2_001, hAdapCuBox_CuFramepb210M2Sum_001);
  BkgParM1[15] = new TBkgModelParameter("PbRom k40",  15, 0., 1E-7, 0, 1.0, hAdapPbRomk40M1, hAdapPbRomk40M2, hAdapPbRomk40M2Sum);
  BkgParM1[16] = new TBkgModelParameter("OVC th232",  16, 9.25179e-02, 1E-7, 0, 1.0, hAdapOVCth232M1, hAdapOVCth232M2, hAdapOVCth232M2Sum);
  BkgParM1[17] = new TBkgModelParameter("OVC u238",  17, 1.33486e-01, 1E-7, 0, 1.0, hAdapOVCu238M1, hAdapOVCu238M2, hAdapOVCu238M2Sum);
  BkgParM1[18] = new TBkgModelParameter("OVC co60",  18, 1.94100e-02, 1E-7, 0, 1.0, hAdapOVCco60M1, hAdapOVCco60M2, hAdapOVCco60M2Sum);
  BkgParM1[19] = new TBkgModelParameter("OVC k40",  19, 9.90257e-02, 1E-7, 0, 1.0, hAdapOVCk40M1, hAdapOVCk40M2, hAdapOVCk40M2Sum);  
  BkgParM1[20] = new TBkgModelParameter("External Lead bi210", 20, 0, 1E-7, 0, 1.0, hAdapExtPbbi210M1, hAdapExtPbbi210M2, hAdapExtPbbi210M2Sum);
  BkgParM1[21] = new TBkgModelParameter("CuBox + CuFrame th232", 21, 2.66019e-02, 1E-7, 0, 1.0, hAdapCuBox_CuFrameth232M1, hAdapCuBox_CuFrameth232M2, hAdapCuBox_CuFrameth232M2Sum);
  BkgParM1[22] = new TBkgModelParameter("CuBox + CuFrame u238",  22, 4.98365e-03, 1E-7, 0, 1.0, hAdapCuBox_CuFrameu238M1, hAdapCuBox_CuFrameu238M2, hAdapCuBox_CuFrameu238M2Sum);
  BkgParM1[23] = new TBkgModelParameter("PbRom cs137",  23, 0, 1E-7, 0, 1.0, hAdapPbRomcs137M1, hAdapPbRomcs137M2, hAdapPbRomcs137M2Sum);
  BkgParM1[24] = new TBkgModelParameter("TeO2 Sx th232 1", 24, 1.12618e-03, 1E-7, 0, 1.0, hAdapTeO2Sxth232M1_1, hAdapTeO2Sxth232M2_1, hAdapTeO2Sxth232M2Sum_1);
  BkgParM1[25] = new TBkgModelParameter("TeO2 Sx th232 10", 25, 1.15477e-03, 1E-7, 0, 1.0, hAdapTeO2Sxth232M1_10, hAdapTeO2Sxth232M2_10, hAdapTeO2Sxth232M2Sum_10);
  BkgParM1[26] = new TBkgModelParameter("TeO2 Sx u238 1", 26, 0., 1E-7, 0, 1.0, hAdapTeO2Sxu238M1_1, hAdapTeO2Sxu238M2_1, hAdapTeO2Sxu238M2Sum_1);
  BkgParM1[27] = new TBkgModelParameter("TeO2 Sx u238 10", 27, 0., 1E-7, 0, 1.0, hAdapTeO2Sxu238M1_10, hAdapTeO2Sxu238M2_10, hAdapTeO2Sxu238M2Sum_10);
  BkgParM1[28] = new TBkgModelParameter("TeO2 Sx pb210 10", 28, 0., 1E-7, 0, 1.0, hAdapTeO2Sxpb210M1_10, hAdapTeO2Sxpb210M2_10, hAdapTeO2Sxpb210M2Sum_10);
  BkgParM1[29] = new TBkgModelParameter("CuBox + CuFrame Sx th232 1", 29, 0., 1E-7, 0, 1.0, hAdapCuBox_CuFrameth232M1_1, hAdapCuBox_CuFrameth232M2_1, hAdapCuBox_CuFrameth232M2Sum_1);
  BkgParM1[30] = new TBkgModelParameter("CuBox + CuFrame Sx th232 0.1", 30, 0., 1E-7, 0, 1.0, hAdapCuBox_CuFrameth232M1_01, hAdapCuBox_CuFrameth232M2_01, hAdapCuBox_CuFrameth232M2Sum_01);
  BkgParM1[31] = new TBkgModelParameter("CuBox + CuFrame Sx th232 0.01", 31, 5.75719e-04, 1E-7, 0, 1.0, hAdapCuBox_CuFrameth232M1_001, hAdapCuBox_CuFrameth232M2_001, hAdapCuBox_CuFrameth232M2Sum_001);
  BkgParM1[32] = new TBkgModelParameter("CuBox + CuFrame Sx u238 1", 32, 0., 1E-7, 0, 1.0, hAdapCuBox_CuFrameu238M1_1, hAdapCuBox_CuFrameu238M2_1, hAdapCuBox_CuFrameu238M2Sum_1);
  BkgParM1[33] = new TBkgModelParameter("CuBox + CuFrame Sx u238 0.1", 33, 5.81734e-04, 1E-7, 0, 1.0, hAdapCuBox_CuFrameu238M1_01, hAdapCuBox_CuFrameu238M2_01, hAdapCuBox_CuFrameu238M2Sum_01);
  BkgParM1[34] = new TBkgModelParameter("CuBox + CuFrame Sx u238 0.01", 34, 0., 1E-7, 0, 1.0, hAdapCuBox_CuFrameu238M1_001, hAdapCuBox_CuFrameu238M2_001, hAdapCuBox_CuFrameu238M2Sum_001);
  BkgParM1[35] = new TBkgModelParameter("CuBox + CuFrame Sx pb210 10", 35, 6.09448e-03, 1E-7, 0, 1.0, hAdapCuBox_CuFramepb210M1_10, hAdapCuBox_CuFramepb210M2_10, hAdapCuBox_CuFramepb210M2Sum_10);
  BkgParM1[36] = new TBkgModelParameter("CuBox + CuFrame Sx pb210 1", 36, 5.21599e-05, 1E-7, 0, 1.0, hAdapCuBox_CuFramepb210M1_1, hAdapCuBox_CuFramepb210M2_1, hAdapCuBox_CuFramepb210M2Sum_1);
  BkgParM1[37] = new TBkgModelParameter("Internal th232",  37, 0., 1E-7, 0, 1.0, hAdapInternalth232M1, hAdapInternalth232M2, hAdapInternalth232M2Sum);
  BkgParM1[38] = new TBkgModelParameter("Internal u238",  38, 0., 1E-7, 0, 1.0, hAdapInternalu238M1, hAdapInternalu238M2, hAdapInternalu238M2Sum);
  BkgParM1[39] = new TBkgModelParameter("Internal co60",  39, 0., 1E-7, 0, 1.0, hAdapInternalco60M1, hAdapInternalco60M2, hAdapInternalco60M2Sum);
  BkgParM1[40] = new TBkgModelParameter("Internal k40",  40, 0., 1E-7, 0, 1.0, hAdapInternalk40M1, hAdapInternalk40M2, hAdapInternalk40M2Sum);
  BkgParM1[41] = new TBkgModelParameter("PbRom th232",  41, 0., 1E-7, 0, 1.0, hAdapPbRomth232M1, hAdapPbRomth232M2, hAdapPbRomth232M2Sum);
  BkgParM1[42] = new TBkgModelParameter("PbRom u238",  42, 0., 1E-7, 0, 1.0, hAdapPbRomu238M1, hAdapPbRomu238M2, hAdapPbRomu238M2Sum);
  BkgParM1[43] = new TBkgModelParameter("PbRom co60",  43, 0., 1E-7, 0, 1.0, hAdapPbRomco60M1, hAdapPbRomco60M2, hAdapPbRomco60M2Sum);

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
  for(int i = 0; i < dNParam; i++)
  {
    cout << i << " " << minuit->fCpnam[i] << ": " << fParameters[i] << " +/- " << fParError[i] << endl;
  }
}

void TBkgModel::PrintParActivity()
{ 
  // This only gets the number of counts corrected for detector efficiency
  for(int i = 0; i < dNParam; i++)
  {
  cout << i << " " << minuit->fCpnam[i] << " activity: " << fParActivity[i] << endl;
  }

}

// Resets all parameters and total histograms (ones with addition) to 0
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

  for(int i = 0; i < dNParam; i++)
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

  // Create model
  for(int i = 0; i < dNParam; i++)
  {
    fModelTotAdapM1->Add( BkgParM1[i]->GetHistM1(),              dDataIntegralM1*fParameters[i]);
    fModelTotAdapM2->Add( BkgParM1[i]->GetHistM2(),              dDataIntegralM1*fParameters[i]);
    fModelTotAdapM2Sum->Add( BkgParM1[i]->GetHistM2Sum(),        dDataIntegralM1*fParameters[i]);
  }
  
}  
// This method sets up minuit and does the fit
bool TBkgModel::DoTheFitAdaptive()
{ 
  // Reset initial parameter/error values
  TBkgModel::ResetParameters();

  gStyle->SetOptStat(0);
  gStyle->SetOptFit();

  // Reduce Minuit Output
  minuit->SetPrintLevel(0); // Print level -1 (Quiet), 0 (Normal), 1 (Verbose)
  minuit->Command("SET STRategy 2"); // Sets strategy of fit
  minuit->SetMaxIterations(10000);
  minuit->SetObjectFit(this); //see the external FCN  above

  // Create and fix parameters
  for(int i = 0; i < dNParam; i++)
  {
    minuit->DefineParameter(BkgParM1[i]->GetParIndex(), BkgParM1[i]->GetParName(), BkgParM1[i]->GetParInitial(), BkgParM1[i]->GetParInitErr(), BkgParM1[i]->GetParMin(), BkgParM1[i]->GetParMax());
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

   // Do the minimization
   int status = minuit->Command("MINImize 500000 0.1");

  // Get final parameters from fit
  for(int i = 0; i < dNParam; i++)
  {
    minuit->GetParameter(i, fParameters[i], fParError[i]);
  }

  // Update model with final parameters
  TBkgModel::UpdateModelAdaptive();
  
  dChiSquare = TBkgModel::GetChiSquareAdaptive();

  cout << "Total number of calls = " << dNumCalls << "\t" << "ChiSq/NDF = " << dChiSquare/dNDF << endl; // for M1 and M2
  cout << "ChiSq = " << dChiSquare << "\t" << "NDF = " << dNDF << endl;
  cout << "Probability = " << TMath::Prob(dChiSquare, dNDF) << endl;



  // ///// Draw Data M1
  fAdapDataHistoM1->SetLineColor(kBlack);
  fAdapDataHistoM1->SetLineWidth(2);
  // fAdapDataHistoM1->SetMaximum(90000);
  // fAdapDataHistoM1->GetXaxis()->SetRange(1, fAdapDataHistoM1->FindBin(3000));
  fModelTotAdapM1->SetLineColor(2);
  fModelTotAdapM1->SetLineWidth(1);

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

  ///// Draw Data M2
  fAdapDataHistoM2->SetLineColor(kBlack);
  fAdapDataHistoM2->SetLineWidth(2);
  // fAdapDataHistoM2->SetMaximum(9000);
  // fAdapDataHistoM2->GetXaxis()->SetRange(1, fAdapDataHistoM2->FindBin(3000));
  
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


  // Correlation Matrix section
  // TMatrixT<double> mCorrMatrix;
  // mCorrMatrix.ResizeTo(dNParam, dNParam);
  mCorrMatrix = new TMatrixT<double>(dNParam, dNParam);
  minuit->mnemat(mCorrMatrix->GetMatrixArray(), dNParam);

  for(int i = mCorrMatrix->GetRowLwb(); i <= mCorrMatrix->GetRowUpb(); i++)
    for(int j = mCorrMatrix->GetColLwb(); j <= mCorrMatrix->GetColUpb(); j++)
      mCorrMatrix(i,j) = mCorrMatrix(i,j)/(fParError[i]*fParError[j]);

  // For inverse
  // for(int i = mCorrMatrix->GetRowLwb(); i <= mCorrMatrix->GetRowUpb(); i++)
  //   for(int j = mCorrMatrix->GetColLwb(); j <= mCorrMatrix->GetColUpb(); j++)
  //     mCorrMatrixInverse(i,j) = mCorrMatrix(dNParam-i-1, j); 

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
  mCorrMatrix->Draw("colz");

  // Text on right panel of matrix
  pM2->cd();
  pPave = new TPaveText(0,0,1,1, "NB"); // Text for matrix
  pPave->SetTextSize(0.04);
  pPave->SetFillColor(0);
  pPave->SetBorderSize(0);
  for(int i=0; i < dNParam; i++)
  {
    pPave->AddText(Form("%d: %s", i, minuit->fCpnam[i].Data() ) );
  }
  pPave->Draw();

/*
    How to calculate 2nbb:
     - TeO2 molar mass: 159.6 g/mol
     - half life is 9.81 * 10^20 years
     - how many in Q0 data so far? 1/rate = half life/ln(2) -> rate = ln(2)/half life = 7.066*10^-22 decays/year (Laura's thesis)
     - Moles = 750g * 49 crystals * 0.3408 abundance/159.6 g/mol = 78.474 mol total
     - N_TeO2 = 78.474 * N_A = 4.726*10^25 nuclei of Te130
     - N_2nbb = N_TeO2 * rate * livetime = 1.551*10^4 events
     - half life = rate * ln(2) = ln(2) * N_TeO2 * livetime / N_2nbb
*/

  // This gets the number of counts corrected for detector efficiency
  for(int i = 0; i < dNParam; i++)
  {
    fParActivity[i] = fParameters[i]*dDataIntegralM1/fParEfficiencyM1[i];
    fParActivityErr[i] = fParError[i]*dDataIntegralM1/fParEfficiencyM1[i];
  }


  if(bSave)
  {
  tTime = new TDatime();
  // // Saving plots
  cadap1->SaveAs(Form("%s/FitResults/FitM1_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  cadap2->SaveAs(Form("%s/FitResults/FitM2_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));

  // cResidual1->SaveAs(Form("%s/FitResults/FitM1Residual_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cResidual2->SaveAs(Form("%s/FitResults/FitM2Residual_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cres1->SaveAs(Form("%s/FitResults/FitResidualDist_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cMatrix->SaveAs(Form("%s/FitResults/FitCovMatrix_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cProgress->SaveAs(Form("%s/FitResults/ChiSquareProgress_%d_%d_%d.pdf", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));

  // Save histograms to file
  fSaveResult = new TFile(Form("%s/FitResults/FitResult_%d_%d.root", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime()), "RECREATE");
  fSaveResult->Add(fAdapDataHistoM1);
  fSaveResult->Add(fAdapDataHistoM2);
  fSaveResult->Add(fModelTotAdapM1);
  fSaveResult->Add(fModelTotAdapM2);

  // Scale histograms by parameter for saving (so that actual fit value is represented)
  for(int i = 0; i < dNParam; i++)
  {
    BkgParM1[i]->GetHistM1()->Scale( dDataIntegralM1*fParameters[i]);
    BkgParM1[i]->GetHistM2()->Scale( dDataIntegralM2*fParameters[i]);
    BkgParM1[i]->GetHistM2Sum()->Scale( dDataIntegralM2Sum*fParameters[i]);    

    fSaveResult->Add( BkgParM1[i]->GetHistM1() );
    fSaveResult->Add( BkgParM1[i]->GetHistM2() );    
    fSaveResult->Add( BkgParM1[i]->GetHistM2Sum() );    
  }

  // fSaveResult->Add(&mCorrMatrix);

  fSaveResult->Write(); 

  // cadap1->SaveAs(Form("%s/FitResults/CovMatrix/FitM1_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cadap2->SaveAs(Form("%s/FitResults/CovMatrix/FitM2_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cResidual1->SaveAs(Form("%s/FitResults/CovMatrix/FitM1Residual_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cResidual2->SaveAs(Form("%s/FitResults/CovMatrix/FitM2Residual_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cres1->SaveAs(Form("%s/FitResults/CovMatrix/FitResidualDist_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  // cMatrix->SaveAs(Form("%s/FitResults/CovMatrix/FitCovMatrix_%d_%d_%d.C", dSaveDir.c_str(), tTime->GetDate(), tTime->GetTime(), nLoop));
  } // end bSave



  return true;
}
