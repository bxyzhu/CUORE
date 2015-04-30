// Fitter

void TBkgFitter(double fFitMin, double fFitMax, int fBinBase, int fDataSet, bool fSave)
{
  gSystem->Load("TBkgModelSource_cc.so");

  // Initialize fitter
  // TBkgModel *fBkgModel = new TBkgModel(double fFitMin, double fFitMax, int fBinBase, int fDataSet, bool fSave);	
  // fBkgModel->GenerateParameters()
  TBkgModelSource *fBkgModel = new TBkgModelSource(fFitMin, fFitMax, fBinBase, fDataSet);  

  hAdapTeO22nuM1->Draw();

}
