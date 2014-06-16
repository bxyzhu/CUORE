#include "TMinuit.h"
#include "TBackgroundModel.hh"
#include "TRandom.h"
#include <cmath>
#include <iostream>
#include <string>

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
	Obj->SetParameters(3,x[3]);
	Obj->SetParameters(4,x[4]);

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

	fDataHistoTot = new TH1D("fDataHistoTot", "Data Histogram", 3500, 0., 3500.);
	fDataHistoM1 = new TH1D("fDataHistoM1", "Data Histogram", 3500, 0., 3500.);
	fDataHistoM2 = new TH1D("fDataHistoM2", "Data Histogram", 3500, 0., 3500.);

	// Model histograms
	fModelFrameTh = new TH1D("fModelFrameTh", "Frame", 3500, 0., 3500.);
	fModel50mKTh = new TH1D("fModel50mKTh", "50mK", 3500, 0., 3500.);
	fModel600mKTh = new TH1D("fModel600mKTh", "600mK", 3500, 0., 3500.);
	fModelIVCTh = new TH1D("fModelIVCTh", "IVC", 3500, 0., 3500.);
	fModelOVCTh = new TH1D("fModelOVCTh", "OVC", 3500, 0., 3500.);

	// fRandomGenerator = new TRandom3(0);
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

	delete	fModel50mKTh;
	delete	fModel600mKTh;
	delete 	fModelIVCTh;
	delete	fModelOVCTh;
	delete 	fModelFrameTh;

}


// Prints parameters, make sure to update
void TBackgroundModel::PrintParameters()
{	
	cout<< "Par0 = "<< fParameters[0]<<endl;
	cout<< "Par1 = "<< fParameters[1]<<endl;
	cout<< "Par2 = "<< fParameters[2]<<endl;
	cout<< "Par3 = "<< fParameters[3]<<endl;
	cout<< "Par4 = "<< fParameters[4]<<endl;
	// cout<< "Par5 = "<< fParameters[5]<<endl;
	// cout<< "Par6 = "<< fParameters[6]<<endl;

}


// Loads MC files
TChain *TBackgroundModel::LoadMC(std::string dLocation, std::string dSource, int dMult)
{
    TChain *outTree = new TChain("outTree");
    outTree->Add(Form("/Users/brian/macros/Simulations/Bkg/%s-%s-B-M%d-T50-r0.0425.root", dLocation.c_str(), dSource.c_str(), dMult));

    return outTree;
}

// Loads the data
void TBackgroundModel::LoadData()
{
	if(fDataHistoTot == NULL) {
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
    qtree->Project("fDataHistoTot", "Energy", base_cut);
    qtree->Project("fDataHistoM1", "Energy", base_cut && "Multiplicity_OFTime==1");
    qtree->Project("fDataHistoM2", "Energy", base_cut && "Multiplicity_OFTime==2");

	cout << "Loaded Data" << endl;
}

void TBackgroundModel::ReadMC()
{
    // Loads MC data
    outTree50mKTh = LoadMC("50mK", "Th232", 1);
    outTree600mKTh = LoadMC("600mK", "Th232", 1);
    outTreeIVCTh = LoadMC("IVC", "Th232", 1);
    outTreeOVCTh = LoadMC("OVC", "Th232", 1);
    outTreeFrameTh = LoadMC("Frame", "Th232", 1);

	outTreeFrameTh->Project("fModelFrameTh", "Ener1");
    outTree50mKTh->Project("fModel50mKTh", "Ener1");
    outTree600mKTh->Project("fModel600mKTh", "Ener1");
    outTreeIVCTh->Project("fModelIVCTh", "Ener1");
    outTreeOVCTh->Project("fModelOVCTh", "Ener1");

	cout << "Loaded MC" << endl;

}

// Creates/updates the background model
void TBackgroundModel::UpdateModel()
{
	if(fModelTot == NULL) {
		cout << "Model Histogram Not Created" << endl;
		return;
	}

	// Reset all bins in model histogram(s)
	// Necessary for MC?
	fModelTot->Reset();

	// Create model
	fModelTot->Add(fModelFrameTh,	fParameters[0]);
	fModelTot->Add(fModel50mKTh,	fParameters[1]);
	fModelTot->Add(fModel600mKTh,	fParameters[2]);
	fModelTot->Add(fModelIVCTh,		fParameters[3]);
	fModelTot->Add(fModelOVCTh,		fParameters[4]);





	// Smearing -> Energy resolution depends on energy of event
	// fSmearedMC = TRandom3->Gaus(Energy,EnergyRes);

	// Create total model here or elsewhere?
    // This way of adding doesn't work, should be bin-by-bin?
	// fModelTot = fParameters[0]*fModelFrameTh + fParameters[1]*fModel50mKTh + fParameters[2]*fModel600mKTh + fParameters[3]*fModelIVCTh 
			// + fParameters[4]*fModelOVCTh;

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

	for(int i = 1; i<fModelTot->GetNbinsX()-1; i++)
	{
		data_i = fDataHistoM1->GetBinContent(i);

		model_i = fModelTot->GetBinContent(i);

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
	if(index > 7) cout << "Index too large" << endl;
	else fParameters[index] = value;

}


bool TBackgroundModel::DoTheFit(){
   // This method actually sets up minuit and does the fit

 
   TMinuit minuit(2); //initialise minuit, n is the number of parameters

   // Reduce Minuit Output
   minuit.SetPrintLevel(1);
//   minuit.Command("SET MINImize 1000 0.001");
   minuit.Command("SET STRategy 2");
  //minuit.Command("SET IMProve 1000 ");

//   minuit.SetMaxIterations(1000);	
   minuit.SetObjectFit(this); //see the external FCN  above
   
   //define the parameters and set the ranges and initial guesses see ROOTs TMinuit documentation
   // How to set ranges for various parameters? 0 to ?
   minuit.DefineParameter(0, "Frame Th", 0.4, 1, 0., 10.0);
   minuit.DefineParameter(1, "50 mK Th", 0.4, 1, 0., 10.0);
   minuit.DefineParameter(2, "600 mK Th", 0.4, 1, 0., 10.);
   minuit.DefineParameter(3, "IVC Th", 0.4, 1, 0., 10.);
   minuit.DefineParameter(4, "OVC Th", 0.4, 1, 0., 10.);
   // minuit.DefineParameter(5, "", 0.4, 1, -2, 2);
   // minuit.DefineParameter(6, "", 0.4, 1, -2, 2);

   
   //Tell minuit what external function to use 
   minuit.SetFCN(myExternal_FCN);
   
   int status=minuit.Migrad(); // this actually does the minimisation
   

   //If you want to access the fitted parameters, supposing your class has a container member called fParValues and fParErrors
/* 
   for(Int_t i=0;i<n;i++){
     minuit.GetParameter(i,fParValues[i],fParErrors[i]);
   }
*/   

	fDataHistoM1->Draw();
	double 	dummy;
	minuit.GetParameter(0,fParameters[0],dummy);
	minuit.GetParameter(1,fParameters[1],dummy);
	minuit.GetParameter(2,fParameters[2],dummy);
	minuit.GetParameter(3,fParameters[3],dummy);
	minuit.GetParameter(4,fParameters[4],dummy);

	UpdateModel();
	
	cout<<"At the end ChiSq = "<<GetChiSquare()<<endl;
	
	fModelTot->SetLineColor(kRed);
	fModelTot->Draw("SAME");

	return true;
   
 }

/*
// 
int main(int argc, char** argv)
{
	TBackgroundModel *Fitter = new TBackgroundModel();

	// Fitter->FillData();
	Fitter->DoTheFit();

	return 0;	
}
*/