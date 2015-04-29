#include "TTest.hh"
// Gotta load TBackgroundModel shared object first
// gSystem->Load("TBackgroundModel_cc.so")


ClassImp(TTest)

TTest::TTest()
{
	TBackgroundModel *f1 = new TBackgroundModel(500, 8000, 50, 0, false);
}

TTest::TTest(double fFitMin, double fFitMax, int fBinBase, int fDataSet, bool fSave)
{
	
}


TTest::~TTest()
{}


double TTest::GetChiSquareM1(TH1D *fModelM1, TH1D* fDataM1)
{
  double chiSquare = 0.;
  double datam1_i, errm1_i;
  double modelm1_i;

  int dBinMinM1 = fDataM1->GetBinLowEdge(fDataM1->FindBin(fFitMin));
  int dBinMaxM1 = fDataM1->GetBinLowEdge(fDataM1->FindBin(fFitMin));

  for(int i = dBinMinM1; i < dBinMaxM1; i++)
  {
    if( fDataM1->GetBinCenter(i) >= 3150 && fDataM1->GetBinCenter(i) <= 3400)continue;
    datam1_i = fDataM1->GetBinContent(i)*fDataM1->GetBinWidth(i);
    modelm1_i = fModelM1->GetBinContent(i)*fDataM1->GetBinWidth(i);
    
    if(modelm1_i != 0 && datam1_i != 0)
    {
      // Log-likelihood
      chiSquare += 2 * (modelm1_i - datam1_i + datam1_i * TMath::Log(datam1_i/modelm1_i));
    }
  }


  return chiSquare;
}

double TTest::GetChiSquareM2()
{
  double chiSquare = 0.;
  double datam2_i, errm2_i;
  double modelm2_i;

  for(int i = dFitMinBinM2; i < dFitMaxBinM2; i++)
  {
    if( fAdapDataHistoM2->GetBinCenter(i) >= 3150 && fAdapDataHistoM2->GetBinCenter(i) <= 3400)continue;
    datam2_i = fAdapDataHistoM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i);
    modelm2_i = fModelTotAdapM2->GetBinContent(i)*fAdapDataHistoM2->GetBinWidth(i);
    
    if(modelm2_i != 0 && datam2_i != 0)
    {
      // Log-likelihood
      chiSquare += 2 * (modelm2_i - datam2_i + datam2_i * TMath::Log(datam2_i/modelm2_i));
    }
  }

  return chiSquare;

}
