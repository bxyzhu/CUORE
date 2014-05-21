#include "TMinuit.h"
#include "TBackgroundModel.hh"
#include "TRandom.h"
#include <cmath>
#include <iostream>

using namespace std;

//ClassImp(TBackgroundModel)
  
//first set up a global function that calls your classes method that calculates the quantity to minimise
void myExternal_FCN(int &n, double *grad, double &fval, double x[], int code){
	// Required External Wrapper for function to be minimized by Minuit 
 
 	// This gets called for each value of the parameters minuit tries
	// here the x array contains the parameters you are trying to fit
  
	// here myClass should inherit from TObject
	TBackgroundModel* Obj = (TBackgroundModel*)gMinuit->GetObjectFit(); 

	// implement a method in your class for setting the parameters and thus update the parameters of your fitter class 
	Obj->SetParameters(0,x[0]);   
	Obj->SetParameters(1,x[1]);  
	Obj->SetParameters(2,x[2]);
	Obj->UpdateModel();

//  Obj->PrintParameters();
	//implement a method in your class that calculates the quantity you want to minimise, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimise fval
	fval = Obj->GetChiSquare();
}


TBackgroundModel::TBackgroundModel()
{
	// Data
	qtree = new TChain("qtree");

    base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "abs(BaselineSlope)<0.1";
    base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";

	fDataHistoTot = new TH1D("fDataHistoTot", "Data Histogram", 3000, 0., 3000.);
	fDataHistoM1 = new TH1D("fDataHistoM1", "Data Histogram", 3000, 0., 3000.);
	fDataHistoM2 = new TH1D("fDataHistoM2", "Data Histogram", 3000, 0., 3000.);

	// Model histograms
	fModel50mKTh = new TH1D("fModel50mKTh", "50mK", 3000, 0., 3000.);
	fModelMixingTh = new TH1D("fModelMixingTh", "Mixing Chamber", 3000, 0., 3000.);
	fModel600mKTh = new TH1D("fModel600mKTh", "600mK", 3000, 0., 3000.);
	fModelIVCTh = new TH1D("fModelIVCKTh", "IVC", 3000, 0., 3000.);
	fModelOVCTh = new TH1D("fModelOVCKTh", "OVC", 3000, 0., 3000.);

	fRandomGenerator = new TRandom3(0);
	// fParameters[0] = 0.;
	// fParameters[1] = 0.;

	// Currently filling data in the constructor, good idea?
	LoadData();	

}
  

TBackgroundModel::~TBackgroundModel()
{
	delete	fDataHistoTot;
	delete	fDataHistoM1;
	delete	fDataHistoM2;

	delete	fModelHisto;

}


// Prints parameters, make sure to update
void TBackgroundModel::PrintParameters()
{	
	cout<< "Par0 = "<< fParameters[0]<<endl;
	cout<< "Par1 = "<< fParameters[1]<<endl;

}


// Loads the data
void TBackgroundModel::LoadData()
{
	if(fDataHisto == NULL) {
		cout << "Data Histograms Not Created" << endl;
		return;
	}
	else{
		cout << "Data Histograms Created" << endl;
	}

/*
	// Test
	for(int i =0; i<10000; i++)
	{
		fDataHisto->Fill(gRandom->Gaus(2.0,1.0));
	}
*/	

    qtree->Add("/Users/brian/macros/CUOREZ/Bkg/ReducedCorrectedBkg.root");	
    qtree->Project("fDataHisto", "Energy", base_cut);


	cout << "Loaded Data" << endl;
}


// Creates/updates the background model
void TBackgroundModel::UpdateModel()
{
	if(fModelPDF == NULL) {
		cout << "Model Histogram Not Created" << endl;
		return;
	}

	// Reset all bins in model histogram(s)
	fModel50mKTh->Reset();
	fModelMixingTh->Reset();
	fModel600mKTh->Reset();
	fModelIVCTh->Reset();
	fModelOVCTh->Reset();





	// Smearing -> Energy resolution depends on energy of event
	fSmearedMC = TRandom1->Gaus(Energy,EnergyRes);


	fModel50mKTh = ;
	fModelMixingTh = ;
	fModel600mKTh = ;
	fModelIVCTh = ;
	fModelOVCTh = ;


/*
	// Test
	TF1	gaus("mygaus","gaus(0)",-10,10);
	gaus.SetParameters(1,fParameters[0],fParameters[1]);

	for(int i =1;i<fModelHisto->GetNbinsX(); i++)
	{
		fModelHisto->SetBinContent(i,gaus.Eval(fModelHisto->GetBinCenter(i)));
	}
	
	fModelHisto->Scale(fDataHisto->Integral()/fModelHisto->Integral());
	
*/	


	// cout << "Filled model" << endl;	

}


// Custom ChiSquare function
double TBackgroundModel::GetChiSquare()
{
	//cout<<"Calling GetChiSquare()"<<endl;
	double chiSquare = 0.;
	double data_i, err_i, model_i;	

	for(int i = 1; i<fModelHisto->GetNbinsX()-1; i++)
	{
		data_i = fDataHistoTot->GetBinContent(i);

		model_i = fModelPDF->GetBinContent(i);

//		err_i	=1000.;	
		//if(data_i>0)
		//{
			err_i = sqrt(data_i + model_i); // in real MC, include error for model and other errors
		//}

		//model_i = fModelHisto->GetBinContent(i);
		if(err_i>0)
		{
			chiSquare+=pow(data_i - model_i,2)/pow(err_i,2);
		}
		else
		{

			chiSquare+=pow(data_i - model_i,2);	
		}

	}

	return chiSquare;
}


// Set Parameters in Model
void TBackgroundModel::SetParameters(int index, double value)
{
	// Change the index max depending on model
	if(index > 1) cout << "Index too large" << endl;
	else fParameters[index] = value;

}


bool TBackgroundModel::DoTheFit(){
   // This method actually sets up minuit and does the fit

 
   TMinuit minuit(2); //initialise minuit, n is the number of parametes

   // Reduce Minuit Output
   minuit.SetPrintLevel(1);
//   minuit.Command("SET MINImize 1000 0.001");
   minuit.Command("SET STRategy 2");
  //minuit.Command("SET IMProve 1000 ");

//   minuit.SetMaxIterations(1000);	
   minuit.SetObjectFit(this); //see the external FCN  above
   
   //define the parameters and set the ranges and initial guesses see ROOTs TMinuit documentation
   minuit.DefineParameter(0, "Mean", -2, 100, -2.0, 2.0);
   minuit.DefineParameter(1, "Sigma", 0.4, 1, -2, 2);
   
   //Tell minuit what external function to use 
   minuit.SetFCN(myExternal_FCN);
   
   int status=minuit.Migrad(); // this actually does the minimisation
   

   //If you want to access the fitted parameters, supposing your class has a container member called fParValues and fParErrors
/* 
   for(Int_t i=0;i<n;i++){
     minuit.GetParameter(i,fParValues[i],fParErrors[i]);
   }
*/   

	fDataHistoTot->Draw();
	double 	dummy;
	minuit.GetParameter(0,fParameters[0],dummy);
	minuit.GetParameter(1,fParameters[1],dummy);

	UpdateModel();
	
	cout<<"At the end ChiSq = "<<GetChiSquare()<<endl;
	
	fModelPDF->SetLineColor(kRed);
	fModelPDF->Draw("SAME");

	return true;
   
 }


// 
int main(int argc, char** argv)
{
	TBackgroundModel *Fitter = new TBackgroundModel();

	Fitter->FillData();
	Fitter->DoTheFit();

	return 0;	
}