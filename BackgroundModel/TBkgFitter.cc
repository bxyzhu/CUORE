// Fitter

void TBkgFitter(double fFitMin, double fFitMax, int fBinBase, int fDataSet, bool fSave)
{

  // Initialize fitter
  // TBkgModel *fBkgModel = new TBkgModel(double fFitMin, double fFitMax, int fBinBase, int fDataSet, bool fSave);	
  // fBkgModel->GenerateParameters()


  TBkgModelSource *fBkgModel = new TBkgModelSource(double fFitMin, double fFitMax, int fBinBase, int fDataSet);  


}
