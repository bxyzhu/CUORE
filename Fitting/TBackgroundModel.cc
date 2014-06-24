#include "TMinuit.h"
#include "TBackgroundModel.hh"
#include "TRandom3.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

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
	Obj->SetParameters(10,	x[10]);
	// Obj->SetParameters(11,	x[11]);


	Obj->UpdateModel();

	//implement a method in your class that calculates the quantity you want to minimise, here I call it GetChiSquare. set its output equal to fval. minuit tries to minimise fval
	fval = Obj->GetChiSquare();
}


TBackgroundModel::TBackgroundModel()
{
	Initialize();
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

	delete	fModel50mKRa;
	delete	fModel600mKRa;
	delete 	fModelIVCRa;
	delete	fModelOVCRa;
	delete 	fModelFrameRa;

	delete	fModelFrameK;
	delete	fModelFrameCo;

	delete 	fModelTotTh;
	delete	fModelTotRa;
	delete	fModelTotK;
	delete	fModelTotCo;
}

void TBackgroundModel::Initialize()
{
	// Bin Size in KeV
	dBinSize = 10;
	dMaxEnergy = 2700.;
	dMinEnergy = 0.;
	dNBins = (dMaxEnergy - dMinEnergy)/ dBinSize;
	// Data
	qtree = new TChain("qtree");

    base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "abs(BaselineSlope)<0.1";
    base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";
    base_cut = base_cut && Form("Energy > %f && Energy < %f", dMinEnergy, dMaxEnergy);

    ener_cut = ener_cut && Form("Ener1 > %f && Ener1 < %f", dMinEnergy, dMaxEnergy);

	fDataHistoTot	 = new TH1D("fDataHistoTot", 	"", dNBins, dMinEnergy, dMaxEnergy);
	fDataHistoM1	 = new TH1D("fDataHistoM1", 	"", dNBins, dMinEnergy, dMaxEnergy);
	fDataHistoM2	 = new TH1D("fDataHistoM2", 	"", dNBins, dMinEnergy, dMaxEnergy);

	// Model histograms

	fModelFrameTh	 = new TH1D("fModelFrameTh", 	"Frame",		dNBins, dMinEnergy, dMaxEnergy);
	fModel50mKTh	 = new TH1D("fModel50mKTh", 	"50mK",			dNBins, dMinEnergy, dMaxEnergy);
	fModel600mKTh	 = new TH1D("fModel600mKTh", 	"600mK",		dNBins, dMinEnergy, dMaxEnergy);
	fModelIVCTh		 = new TH1D("fModelIVCTh", 		"IVC", 			dNBins, dMinEnergy, dMaxEnergy);
	fModelOVCTh		 = new TH1D("fModelOVCTh",		"OVC",	 		dNBins, dMinEnergy, dMaxEnergy);

	fModelFrameRa	 = new TH1D("fModelFrameRa", 	"Frame", 		dNBins, dMinEnergy, dMaxEnergy);
	fModel50mKRa	 = new TH1D("fModel50mKRa", 	"50mK", 		dNBins, dMinEnergy, dMaxEnergy);
	fModel600mKRa	 = new TH1D("fModel600mKRa", 	"600mK", 		dNBins, dMinEnergy, dMaxEnergy);
	fModelIVCRa		 = new TH1D("fModelIVCRa", 		"IVC", 			dNBins, dMinEnergy, dMaxEnergy);
	fModelOVCRa		 = new TH1D("fModelOVCRa", 		"OVC", 			dNBins, dMinEnergy, dMaxEnergy);

	fModelFrameK	 = new TH1D("fModelFrameK", 	"Frame", 		dNBins, dMinEnergy, dMaxEnergy);
	fModelFrameCo	 = new TH1D("fModelFrameCo", 	"Frame", 		dNBins, dMinEnergy, dMaxEnergy);

	// Total model histograms
	fModelTot 		 = new TH1D("fModelTot", 		"Frame", 		dNBins, dMinEnergy, dMaxEnergy);	
	fModelTotTh		 = new TH1D("fModelTotTh", 		"Total Th232", 	dNBins, dMinEnergy, dMaxEnergy);
	fModelTotRa		 = new TH1D("fModelTotRa", 		"Total Ra226", 	dNBins, dMinEnergy, dMaxEnergy);
	fModelTotK		 = new TH1D("fModelTotK", 		"Total K40",	dNBins, dMinEnergy, dMaxEnergy);
	fModelTotCo		 = new TH1D("fModelTotCo", 		"Total Co60", 	dNBins, dMinEnergy, dMaxEnergy);

	// Initial Parameters
	// fRandomGenerator = new TRandom3(0);
	fParameters[0] 	= 0.;
	fParameters[1] 	= 0.;
	fParameters[2] 	= 0.;
	fParameters[3] 	= 0.;
	fParameters[4] 	= 0.;
	fParameters[5] 	= 0.;
	fParameters[6] 	= 0.;
	fParameters[7] 	= 0.;
	fParameters[8]	= 0.;
	fParameters[9] 	= 0.;
	fParameters[10] = 0.;
	// fParameters[11] = 0.;

	// Loading all data in Initialize, correct or no?
	LoadData();	
	ReadMC();

}


// Prints parameters, make sure to update
void TBackgroundModel::PrintParameters()
{	
	cout<< "Par0 = "	<< fParameters[0]	<< endl;
	cout<< "Par1 = "	<< fParameters[1]	<< endl;
	cout<< "Par2 = "	<< fParameters[2]	<< endl;
	cout<< "Par3 = "	<< fParameters[3]	<< endl;
	cout<< "Par4 = "	<< fParameters[4]	<< endl;

	cout<< "Par5 = "	<< fParameters[5]	<< endl;
	cout<< "Par6 = "	<< fParameters[6]	<< endl;
	cout<< "Par7 = "	<< fParameters[7]	<< endl;
	cout<< "Par8 = "	<< fParameters[8]	<< endl;
	cout<< "Par9 = "	<< fParameters[9]	<< endl;
	cout<< "Par10 = "	<< fParameters[10]	<< endl;
	// cout<< "Par11 = "	<< fParameters[11]	<< endl;

	cout << "Sum = " << endl; 

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
	if(fDataHistoTot == NULL) 
	{
		cout << "Data Histograms Not Created" << endl;
		return;
	}
	else
	{
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
    qtree->Project("fDataHistoM1", 	"Energy", base_cut && "Multiplicity_OFTime==1");
    qtree->Project("fDataHistoM2", 	"Energy", base_cut && "Multiplicity_OFTime==2");

	cout << "Loaded Data" << endl;

	NormalizePDF(fDataHistoM1);

	cout << "Normalized Data" << endl;
}

void TBackgroundModel::ReadMC()
{
    // Loads MC data
    outTree50mKTh 	= LoadMC("50mK",	"Th232", 1);
    outTree600mKTh 	= LoadMC("600mK", 	"Th232", 1);
    outTreeIVCTh 	= LoadMC("IVC", 	"Th232", 1);
    outTreeOVCTh 	= LoadMC("OVC", 	"Th232", 1);
    outTreeFrameTh 	= LoadMC("Frame", 	"Th232", 1);

    outTree50mKRa 	= LoadMC("50mK", 	"Ra226", 1);
    outTree600mKRa	= LoadMC("600mK", 	"Ra226", 1);
    outTreeIVCRa 	= LoadMC("IVC", 	"Ra226", 1);
    outTreeOVCRa 	= LoadMC("OVC", 	"Ra226", 1);
    outTreeFrameRa 	= LoadMC("Frame", 	"Ra226", 1);

    outTreeFrameK 	= LoadMC("Frame", 	"K40",	 1);
    outTreeFrameCo 	= LoadMC("Frame", 	"Co60",	 1);

	outTreeFrameTh->Project("fModelFrameTh", 	"Ener1", ener_cut);
    outTree50mKTh->Project("fModel50mKTh", 		"Ener1", ener_cut);
    outTree600mKTh->Project("fModel600mKTh", 	"Ener1", ener_cut);
    outTreeIVCTh->Project("fModelIVCTh", 		"Ener1", ener_cut);
    outTreeOVCTh->Project("fModelOVCTh", 		"Ener1", ener_cut);

	outTreeFrameRa->Project("fModelFrameRa", 	"Ener1", ener_cut);
    outTree50mKRa->Project("fModel50mKRa", 		"Ener1", ener_cut);
    outTree600mKRa->Project("fModel600mKRa", 	"Ener1", ener_cut);
    outTreeIVCRa->Project("fModelIVCRa", 		"Ener1", ener_cut);
    outTreeOVCRa->Project("fModelOVCRa", 		"Ener1", ener_cut);

	outTreeFrameK->Project("fModelFrameK", 		"Ener1", ener_cut);
	outTreeFrameCo->Project("fModelFrameCo", 	"Ener1", ener_cut);


	cout << "Loaded MC" << endl;

	// Normalize all MC histograms
	NormalizePDF(fModelFrameTh);
	NormalizePDF(fModel50mKTh);
	NormalizePDF(fModel600mKTh);
	NormalizePDF(fModelIVCTh);
	NormalizePDF(fModelOVCTh);

	NormalizePDF(fModelFrameRa);
	NormalizePDF(fModel50mKRa);
	NormalizePDF(fModel600mKRa);
	NormalizePDF(fModelIVCRa);
	NormalizePDF(fModelOVCRa);

	NormalizePDF(fModelFrameK);
	NormalizePDF(fModelFrameCo);

	cout << "Normalized MC PDFs" << endl;
}

// Creates/updates the background model
void TBackgroundModel::UpdateModel()
{
	if(fModelTot == NULL) 
	{
		cout << "Model Histogram Not Created" << endl;
		return;
	}

	// Reset all bins in model histogram(s)
	// Necessary for MC?
	fModelTot->Reset();

	// Create model
	fModelTot->Add(fModelFrameTh,	1-fParameters[0] + fParameters[1] + fParameters[2] + fParameters[3] + fParameters[4] + fParameters[5]
									+ fParameters[6] + fParameters[7] + fParameters[8] + fParameters[9] + fParameters[10]);

	fModelTot->Add(fModel50mKTh,	fParameters[0]);
	fModelTot->Add(fModel600mKTh,	fParameters[1]);
	fModelTot->Add(fModelIVCTh,		fParameters[2]);
	fModelTot->Add(fModelOVCTh,		fParameters[3]);

	fModelTot->Add(fModelFrameRa,	fParameters[4]);
	fModelTot->Add(fModel50mKRa,	fParameters[5]);
	fModelTot->Add(fModel600mKRa,	fParameters[6]);
	fModelTot->Add(fModelIVCRa,		fParameters[7]);
	fModelTot->Add(fModelOVCRa,		fParameters[8]);

	fModelTot->Add(fModelFrameK,	fParameters[9]);
	fModelTot->Add(fModelFrameCo, 	fParameters[10]);


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

	// Start from bin 60 (5 keV bins -> 300 keV)
	for(int i = 2550/dBinSize; i < fModelTot->GetNbinsX(); i++)
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


void TBackgroundModel::NormalizePDF(TH1D *h1)
{
	double dIntegral;

	// bin 0 = underflow, bin dNBins = last bin with upper-edge xup Excluded
	dIntegral = h1->Integral(1, dNBins);

	h1->Scale(1/dIntegral);
}


TGraph *TBackgroundModel::CalculateResiduals(TH1D *h1, TH1D *h2)
{

	// Clone histograms for rebinning
	TH1D 	*hCloneBkg 		= (TH1D*)h1->Clone("hCloneBkg");
	TH1D 	*hCloneMC		= (TH1D*)h2->Clone("hCloneMC");
	TGraph	*g1				= new TGraph();
	// hResidualDist = new TH1D("hResidualDist","Residual Distribution", 39, -8, 11);


	// Variables used in Residual calculations
	double dResidualX, dResidualY, dResidualXErr = 0, dResidualYErr;

	// Residual plot and distribution
	for (int j = 300/dBinSize; j < hCloneMC->GetNbinsX(); j++)
	{
		dResidualX 		= hCloneBkg->GetBinCenter(j);
		dResidualY 		= (hCloneBkg->GetBinContent(j) - hCloneMC->GetBinContent(j))/
							TMath::Sqrt(hCloneBkg->GetBinContent(j) + hCloneMC->GetBinContent(j)); // Sqrt(MC + data) = sigma for poisson distribution

		g1->SetPoint(j, dResidualX, dResidualY);
		// hResidualDist->Fill(dResidualY);
	}

	return g1;
}


// Set Parameters in Model
void TBackgroundModel::SetParameters(int index, double value)
{
	// Change the index max depending on model
	if(index > 11) cout << "Index too large" << endl;
	else fParameters[index] = value;

}



bool TBackgroundModel::DoTheFit()
{
	gStyle->SetOptStat(0);
   // This method actually sets up minuit and does the fit

 
   TMinuit minuit(2); //initialise minuit, n is the number of parameters

   // Reduce Minuit Output
   minuit.SetPrintLevel(1);
//   minuit.Command("SET MINImize 1000 0.001");
   minuit.Command("SET STRategy 2");
  //minuit.Command("SET IMProve 1000 ");

   minuit.SetMaxIterations(10000);
   minuit.SetObjectFit(this); //see the external FCN  above
   
   //define the parameters and set the ranges and initial guesses see ROOTs TMinuit documentation
   // How to set ranges for various parameters? 0 to ?
   // minuit.DefineParameter(0, "Frame Th", 	0., 1., 0., 1.0);
   minuit.DefineParameter(0, "50 mK Th", 	0.5, 1., 0., 1.0);
   minuit.DefineParameter(1, "600 mK Th",	0., 1., 0., 1.0);
   minuit.DefineParameter(2, "IVC Th", 		0., 1., 0., 1.0);
   minuit.DefineParameter(3, "OVC Th", 		0., 1., 0., 1.0);

   minuit.DefineParameter(4, "Frame Ra", 	0., 1., 0., 1.0);
   minuit.DefineParameter(5, "50 mK Ra", 	0., 1., 0., 1.0);
   minuit.DefineParameter(6, "600 mK Ra", 	0., 1., 0., 1.0);
   minuit.DefineParameter(7, "IVC Ra", 		0., 1., 0., 1.0);
   minuit.DefineParameter(8, "OVC Ra", 		0., 1., 0., 1.0);

   minuit.DefineParameter(9, "Frame K", 	0.5, 1., 0., 1.0);
   minuit.DefineParameter(10, "Frame Co", 	0.5, 1., 0., 1.0);
   
   // Fix parameters for testing
   // minuit.FixParameter(0);
   minuit.FixParameter(1);
   minuit.FixParameter(2);
   minuit.FixParameter(3);

   minuit.FixParameter(5);
   minuit.FixParameter(6);
   minuit.FixParameter(7);
   minuit.FixParameter(8);


   //Tell minuit what external function to use 
   minuit.SetFCN(myExternal_FCN);
   
   int status = minuit.Migrad(); // this actually does the minimisation
   

   //If you want to access the fitted parameters, supposing your class has a container member called fParValues and fParErrors
/* 
   for(Int_t i=0;i<n;i++){
     minuit.GetParameter(i,fParValues[i],fParErrors[i]);
   }
*/  
	fModelTotTh->Add(fModelFrameTh,		1 - fParameters[0] + fParameters[1] + fParameters[2] + fParameters[3] + fParameters[4] + fParameters[5]
										+ fParameters[6] + fParameters[7] + fParameters[8] + fParameters[9] + fParameters[10]);
	fModelTotTh->Add(fModel50mKTh,		fParameters[0]);
	fModelTotTh->Add(fModel600mKTh,		fParameters[1]);
	fModelTotTh->Add(fModelIVCTh,		fParameters[2]);
	fModelTotTh->Add(fModelOVCTh,		fParameters[3]);

	fModelTotRa->Add(fModelFrameRa,		fParameters[4]);
	fModelTotRa->Add(fModel50mKRa,		fParameters[5]);
	fModelTotRa->Add(fModel600mKRa,		fParameters[6]);
	fModelTotRa->Add(fModelIVCRa,		fParameters[7]);
	fModelTotRa->Add(fModelOVCRa,		fParameters[8]);

	fModelTotK->Add(fModelFrameK,		fParameters[9]);
	fModelTotCo->Add(fModelFrameCo,		fParameters[10]);


    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);

   	fDataHistoM1->SetLineColor(1);
   	fDataHistoM1->SetLineWidth(2);
	fDataHistoM1->Draw();

	double 	dummy;

	minuit.GetParameter(0,	fParameters[0],		dummy);
	minuit.GetParameter(1,	fParameters[1],		dummy);
	minuit.GetParameter(2,	fParameters[2],		dummy);
	minuit.GetParameter(3,	fParameters[3],		dummy);
	minuit.GetParameter(4,	fParameters[4],		dummy);
	minuit.GetParameter(5,	fParameters[5],		dummy);
	minuit.GetParameter(6,	fParameters[6],		dummy);
	minuit.GetParameter(7,	fParameters[7],		dummy);
	minuit.GetParameter(8,	fParameters[8],		dummy);
	minuit.GetParameter(9,	fParameters[9],		dummy);
	minuit.GetParameter(10,	fParameters[10],	dummy);	
	// minuit.GetParameter(11,	fParameters[11],	dummy);

	UpdateModel();
	
	cout << "At the end ChiSq = " << GetChiSquare()<<endl;
	
	fModelTot->SetLineColor(2);
	fModelTot->Draw("SAME");

	fModelTotTh->SetLineColor(3);
	fModelTotTh->SetLineStyle(2);
	fModelTotRa->SetLineColor(4);
	fModelTotRa->SetLineStyle(2);
	fModelTotK->SetLineColor(6);
	fModelTotK->SetLineStyle(2);
	fModelTotCo->SetLineColor(7);
	fModelTotCo->SetLineStyle(2);

	fModelTotTh->Draw("SAME");
	fModelTotRa->Draw("SAME");
	fModelTotK->Draw("SAME");
	fModelTotCo->Draw("SAME");

	TCanvas *cResidual = new TCanvas("cResidual", "cResidual", 1200, 800);
	gResidual = CalculateResiduals(fModelTot, fDataHistoM1);

	TAxis *axis = gResidual->GetXaxis();
	axis->SetLimits(0, 3000);
	gResidual->SetName("Residuals");
	gResidual->GetXaxis()->SetTitle("Energy (keV)");
	gResidual->GetXaxis()->SetTitleSize(0.04);
	gResidual->GetXaxis()->SetLabelSize(0.05);
	gResidual->GetYaxis()->SetLabelSize(0.05);
	gResidual->GetYaxis()->SetTitleSize(0.04);
	gResidual->SetMarkerStyle(4);
	gResidual->SetMarkerSize(2);	
	gResidual->GetYaxis()->SetTitle("Residuals (#sigma)");


	gResidual->Draw("AP");

	return true;
   
 }
