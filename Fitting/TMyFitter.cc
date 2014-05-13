#include "TMinuit.h"
#include "TMyFitter.hh"
#include "TRandom.h"
#include <cmath>
#include <iostream>

using namespace std;

//ClassImp(TMyFitter)
  
//first set up a global function that calls your classes method that calculates the quantity to minimise
void myExternal_FCN(int &n, double *grad, double &fval, double x[], int code){
  // Required External Wrapper for function to be minimized by Minuit 
 
  // This gets called for each value of the parameters minuit tries
  // here the x array contains the parameters you are trying to fit
  
  // here myClass should inherit from TObject
  TMyFitter* Obj = (TMyFitter*)gMinuit->GetObjectFit(); 

  // implement a method in your class for setting the parameters and thus update the parameters of your fitter class 
  Obj->SetParameters(0,x[0]);   
  Obj->SetParameters(1,x[1]);  

  Obj->UpdateModel();

//  Obj->PrintParameters();
  //implement a method in your class that calculates the quantity you want to minimise, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimise fval
  fval = Obj->GetChiSquare();   
}


TMyFitter::TMyFitter()
{
	fParameters[0] = 0.;
	fParameters[1] = 0.;
	fDataHisto = new TH1D("fDataHisto", "Data Histogram", 100, -10., 10.);
	fModelHisto = new TH1D("fModelHisto", "Model Histogram", 100, -10., 10.);
	FillData();	

}
  

TMyFitter::~TMyFitter()
{
	delete fDataHisto;
	delete fModelHisto;

}


void TMyFitter::PrintParameters()
{	
cout<< "Par0 = "<< fParameters[0]<<endl;
cout<< "Par1 = "<< fParameters[1]<<endl;

}

void TMyFitter::FillData()
{
	if(fDataHisto==NULL) {
		cout << "Data Histogram Not Created" << endl;
		return;
	}
	else{
		cout << "Data Histogram Created" << endl;
	}

	for(int i =0;i<10000; i++)
	{
		fDataHisto->Fill(gRandom->Gaus(2.0,1.0));
	}
	

	cout << "Filled fake data" << endl;
}

void TMyFitter::UpdateModel()
{
	if(fModelHisto == NULL) {
		cout << "Model Histogram Not Created" << endl;
		return;
	}

	fModelHisto->Reset();//resets all bins
// fModelHisto->Clear();
TF1	gaus("mygaus","gaus(0)",-10,10);
gaus.SetParameters(1,fParameters[0],fParameters[1]);

	for(int i =1;i<fModelHisto->GetNbinsX(); i++)
	{
		fModelHisto->SetBinContent(i,gaus.Eval(fModelHisto->GetBinCenter(i)));
	}
fModelHisto->Scale(fDataHisto->Integral()/fModelHisto->Integral());
	
	// cout << "Filled model" << endl;	

}

double TMyFitter::GetChiSquare()
{
//cout<<"Calling GetChiSquare()"<<endl;
	double chiSquare=0.;
	double	data_i,err_i,model_i;	

	for(int i = 1;i<fModelHisto->GetNbinsX()-1; i++)
	{
		data_i = fDataHisto->GetBinContent(i);
		model_i = fModelHisto->GetBinContent(i);

//		err_i	=1000.;	
		//if(data_i>0)
		//{
			err_i = sqrt(data_i); // in real MC, include error for model and other errors
		//}

		//model_i = fModelHisto->GetBinContent(i);
		if(err_i>0)
		{
			chiSquare+=pow(data_i - model_i,2)/pow(err_i,2);
		}
		else{

						chiSquare+=pow(data_i - model_i,2);	
		}

	}

	return chiSquare;
}

void TMyFitter::SetParameters(int index, double value)
{
	if(index > 1) cout << "Index too large" << endl;
	else fParameters[index] = value;

}


bool TMyFitter::DoTheFit(){
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
   minuit.DefineParameter(0, "Mean",-2,100,-2.0,2.0);
   minuit.DefineParameter(1, "Sigma",0.4,1,-2,2);
   
   //Tell minuit what external function to use 
   minuit.SetFCN(myExternal_FCN);
   
   Int_t status=minuit.Migrad(); // this actually does the minimisation
   

   //If you want to access the fitted parameters, supposing your class has a container member called fParValues and fParErrors
/* 
   for(Int_t i=0;i<n;i++){
     minuit.GetParameter(i,fParValues[i],fParErrors[i]);
   }
*/   

fDataHisto->Draw();
double 	dummy;
minuit.GetParameter(0,fParameters[0],dummy);
minuit.GetParameter(1,fParameters[1],dummy);
UpdateModel();
cout<<"At the end ChiSq = "<<GetChiSquare()<<endl;
fModelHisto->SetLineColor(kRed);
fModelHisto->Draw("SAME");
   return true;
   
 }

int main(int argc, char** argv)
{
	TMyFitter *Fitter = new TMyFitter();

	Fitter->FillData();
	Fitter->DoTheFit();

	return 0;	
}