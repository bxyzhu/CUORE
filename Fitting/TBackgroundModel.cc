#include "TMinuit.h"
#include "TBackgroundModel.hh"
#include "TRandom3.h"
#include "TPaveText.h"
#include "TAxis.h"
#include <cmath>
#include <iostream>
#include <string>

// Bookmarks: f2 = DoTheFit, GetChiSquare, Initialize


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

/*
	Obj->SetParameters(9,	x[9]);
	Obj->SetParameters(10,	x[10]);
	Obj->SetParameters(11,	x[11]);
	Obj->SetParameters(12,	x[12]);   
	Obj->SetParameters(13,	x[13]);  
	Obj->SetParameters(14,	x[14]);
	Obj->SetParameters(15,	x[15]);
	Obj->SetParameters(16,	x[16]);
	Obj->SetParameters(17,	x[17]);   
	Obj->SetParameters(18,	x[18]);  
	Obj->SetParameters(19,	x[19]);
	Obj->SetParameters(20,	x[20]);
	Obj->SetParameters(21,	x[21]);
	Obj->SetParameters(22,	x[22]);
	Obj->SetParameters(23,	x[23]);
*/
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
	delete	fToyData;

	delete 	fModelFrameTh;
	delete	fModelTShieldTh;
	delete	fModel50mKTh;
	delete	fModel600mKTh;
	delete 	fModelIVCTh;
	delete	fModelOVCTh;

	delete	fModelFrameRa
	delete 	fModelTShieldRa;
	delete	fModel50mKRa;
	delete	fModel600mKRa;
	delete 	fModelIVCRa;
	delete	fModelOVCRa;

	delete	fModelFrameK;
	delete 	fModelTShieldK;
	delete	fModel50mKK;
	delete	fModel600mKK;
	delete 	fModelIVCK;
	delete	fModelOVCK;

	delete	fModelFrameCo;
	delete 	fModelTShieldCo;
	delete	fModel50mKCo;
	delete	fModel600mKCo;
	delete 	fModelIVCCo;
	delete	fModelOVCCo;

	delete 	fModelTotTh;
	delete	fModelTotRa;
	delete	fModelTotK;
	delete	fModelTotCo;
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
		dResidualY 		= (hCloneBkg->GetBinContent(j) - hCloneMC->GetBinContent(j)) /
							TMath::Sqrt(hCloneBkg->GetBinContent(j)); // Sqrt(MC + data) = sigma for poisson distribution

		g1->SetPoint(j, dResidualX, dResidualY);
		// hResidualDist->Fill(dResidualY);
	}

	return g1;
}



bool TBackgroundModel::DoTheFit()
{
	gStyle->SetOptStat(0);
   // This method actually sets up minuit and does the fit

 
   TMinuit minuit(2); //initialize minuit, n is the number of parameters

   // Reduce Minuit Output
   minuit.SetPrintLevel(1);
//   minuit.Command("SET MINImize 1000 0.001");
   minuit.Command("SET STRategy 2");
  //minuit.Command("SET IMProve 1000 ");
   minuit.SetMaxIterations(50000);
   minuit.SetObjectFit(this); //see the external FCN  above
   
   //define the parameters and set the ranges and initial guesses see ROOTs TMinuit documentation
   // Range is from 0 to integral of the data
   // Around 60000 events in background spectrum


   ////////////////////////////////////////////////
   // Using less parameters
   ////////////////////////////////////////////////
   minuit.DefineParameter(0, "Close Th", 	2000., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(1, "Far Th",	 	2000., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(2, "Close Ra", 	1000., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(3, "Far Ra",		0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(4, "Close K", 	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(5, "Far K", 		0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(6, "Close Co", 	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(7, "Far Co",	 	0., 100.0, 0., dDataIntegral);  
   minuit.DefineParameter(8, "#sigma",	 	4., 1.0, 0., 10);  

   // Fix parameters for testing
   // minuit.FixParameter(0);
   // minuit.FixParameter(1);
   // minuit.FixParameter(2);
   minuit.FixParameter(3);
   minuit.FixParameter(4);
   minuit.FixParameter(5);
   minuit.FixParameter(6);
   minuit.FixParameter(7);



	//////////////////////////////////////////////
	// All parameters
   ///////////////////////////////////////////////
/*
   minuit.DefineParameter(0, "Frame Th", 	1000., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(1, "TShield Th", 	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(2, "50 mK Th", 	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(3, "600 mK Th",	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(4, "IVC Th", 		0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(5, "OVC Th", 		0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(6, "Frame Ra", 	500., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(7, "TShield Ra", 	0., 100.0, 0., dDataIntegral);   
   minuit.DefineParameter(8, "50 mK Ra", 	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(9, "600 mK Ra", 	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(10, "IVC Ra", 	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(11, "OVC Ra", 	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(12, "Frame K", 	500., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(13, "TShield K", 	0., 100.0, 0., dDataIntegral);   
   minuit.DefineParameter(14, "50 mK K", 	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(15, "600 mK K", 	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(16, "IVC K", 		0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(17, "OVC K", 		0., 100.0, 0., dDataIntegral);   
   minuit.DefineParameter(18, "Frame Co", 	500., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(19, "TShield Co", 0., 100.0, 0., dDataIntegral);   
   minuit.DefineParameter(20, "50 mK Co", 	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(21, "600 mK Co", 	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(22, "IVC Co", 	0., 100.0, 0., dDataIntegral);
   minuit.DefineParameter(23, "OVC Co", 	0., 100.0, 0., dDataIntegral);   
*/

   // Fix parameters for testing
   // minuit.FixParameter(0);
   // minuit.FixParameter(1);
   // minuit.FixParameter(2);
   // minuit.FixParameter(3);
   // minuit.FixParameter(4);
   // minuit.FixParameter(5);
   // minuit.FixParameter(6);
   // minuit.FixParameter(7);
   // minuit.FixParameter(8);
   // minuit.FixParameter(9);
   // minuit.FixParameter(10);
   // minuit.FixParameter(11);
   // minuit.FixParameter(12);
   // minuit.FixParameter(13);
   // minuit.FixParameter(14);
   // minuit.FixParameter(15);
   // minuit.FixParameter(16);
   // minuit.FixParameter(17);
   // minuit.FixParameter(18);
   // minuit.FixParameter(19);
   // minuit.FixParameter(20);
   // minuit.FixParameter(21);
   // minuit.FixParameter(22);
   // minuit.FixParameter(23);


   //Tell minuit what external function to use 
   minuit.SetFCN(myExternal_FCN);
   
   int status = minuit.Migrad(); // this actually does the minimisation
   

   //If you want to access the fitted parameters, supposing your class has a container member called fParValues and fParErrors
/* 
   for(Int_t i=0;i<n;i++){
     minuit.GetParameter(i,fParValues[i],fParErrors[i]);
   }
*/


   ///////////////////////////////////////////
   //// Few Parameters
   ///////////////////////////////////////////

	fModelTotTh->Add(fModelFrameTh,		fParameters[0]);
	fModelTotTh->Add(fModelTShieldTh,	fParameters[0]);
	fModelTotTh->Add(fModel50mKTh,		fParameters[0]);
	fModelTotTh->Add(fModel600mKTh,		fParameters[0]);
	fModelTotTh->Add(fModelIVCTh,		fParameters[1]);
	fModelTotTh->Add(fModelOVCTh,		fParameters[1]);

	fModelTotRa->Add(fModelFrameRa,		fParameters[2]);
	fModelTotRa->Add(fModelTShieldRa,	fParameters[2]);
	fModelTotRa->Add(fModel50mKRa,		fParameters[2]);
	fModelTotRa->Add(fModel600mKRa,		fParameters[2]);
	fModelTotRa->Add(fModelIVCRa,		fParameters[3]);
	fModelTotRa->Add(fModelOVCRa,		fParameters[3]);

	fModelTotK->Add(fModelFrameK,		fParameters[4]);
	fModelTotK->Add(fModelTShieldK,		fParameters[4]);
	fModelTotK->Add(fModel50mKK,		fParameters[4]);
	fModelTotK->Add(fModel600mKK,		fParameters[4]);
	fModelTotK->Add(fModelIVCK,			fParameters[5]);
	fModelTotK->Add(fModelOVCK,			fParameters[5]);

	fModelTotCo->Add(fModelFrameCo,		fParameters[6]);
	fModelTotCo->Add(fModelTShieldCo,	fParameters[6]);
	fModelTotCo->Add(fModel50mKCo,		fParameters[6]);
	fModelTotCo->Add(fModel600mKCo,		fParameters[6]);
	fModelTotCo->Add(fModelIVCCo,		fParameters[7]);
	fModelTotCo->Add(fModelOVCCo,		fParameters[7]);


	///////////////////////////////////////////
	//// All Parameters
	///////////////////////////////////////////
/*
	fModelTotTh->Add(fModelFrameTh,		fParameters[0]);
	fModelTotTh->Add(fModelTShieldTh,	fParameters[1]);
	fModelTotTh->Add(fModel50mKTh,		fParameters[2]);
	fModelTotTh->Add(fModel600mKTh,		fParameters[3]);
	fModelTotTh->Add(fModelIVCTh,		fParameters[4]);
	fModelTotTh->Add(fModelOVCTh,		fParameters[5]);

	fModelTotRa->Add(fModelFrameRa,		fParameters[6]);
	fModelTotRa->Add(fModelTShieldRa,	fParameters[7]);
	fModelTotRa->Add(fModel50mKRa,		fParameters[8]);
	fModelTotRa->Add(fModel600mKRa,		fParameters[9]);
	fModelTotRa->Add(fModelIVCRa,		fParameters[10]);
	fModelTotRa->Add(fModelOVCRa,		fParameters[11]);

	fModelTotK->Add(fModelFrameK,		fParameters[12]);
	fModelTotK->Add(fModelTShieldK,		fParameters[13]);
	fModelTotK->Add(fModel50mKK,		fParameters[14]);
	fModelTotK->Add(fModel600mKK,		fParameters[15]);
	fModelTotK->Add(fModelIVCK,			fParameters[16]);
	fModelTotK->Add(fModelOVCK,			fParameters[17]);

	fModelTotCo->Add(fModelFrameCo,		fParameters[18]);
	fModelTotCo->Add(fModelTShieldCo,	fParameters[19]);
	fModelTotCo->Add(fModel50mKCo,		fParameters[20]);
	fModelTotCo->Add(fModel600mKCo,		fParameters[21]);
	fModelTotCo->Add(fModelIVCCo,		fParameters[22]);
	fModelTotCo->Add(fModelOVCCo,		fParameters[23]);
*/
	// fModelTotCo->Add(fModelFrameCo,		1 - (fParameters[0] + fParameters[1] + fParameters[2] + fParameters[3] + fParameters[4] + fParameters[5]
										// + fParameters[6] + fParameters[7] + fParameters[8] + fParameters[9] + fParameters[10]));


    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->SetLogy();


    if(bToyFit)
    {
		////// Draw Toy Data
   		fToyData->SetLineColor(1);
   		fToyData->SetLineWidth(2);
   		fToyData->GetXaxis()->SetTitle("Energy (keV)");
   		fToyData->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));
		fToyData->Draw();
	}
	else
	{
   		///// Draw Data
  	 	fDataHistoM1->SetLineColor(1);
  	 	fDataHistoM1->SetLineWidth(2);
  	 	fDataHistoM1->GetXaxis()->SetTitle("Energy (keV)");
   		fDataHistoM1->GetYaxis()->SetTitle(Form("Counts/(%d keV)", dBinSize));
		fDataHistoM1->Draw();
	}

	// Dummy variable for error of parameter to throw away
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
/*	
	minuit.GetParameter(9,	fParameters[9],		dummy);
	minuit.GetParameter(10,	fParameters[10],	dummy);	
	minuit.GetParameter(11,	fParameters[11],	dummy);
	minuit.GetParameter(12,	fParameters[12],	dummy);
	minuit.GetParameter(13,	fParameters[13],	dummy);
	minuit.GetParameter(14,	fParameters[14],	dummy);
	minuit.GetParameter(15,	fParameters[15],	dummy);
	minuit.GetParameter(16,	fParameters[16],	dummy);
	minuit.GetParameter(17,	fParameters[17],	dummy);
	minuit.GetParameter(18,	fParameters[18],	dummy);
	minuit.GetParameter(19,	fParameters[19],	dummy);
	minuit.GetParameter(20,	fParameters[20],	dummy);
	minuit.GetParameter(21,	fParameters[21],	dummy);
	minuit.GetParameter(22,	fParameters[22],	dummy);	
	minuit.GetParameter(23,	fParameters[23],	dummy);
*/

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


 	TLegend *legfit = new TLegend(0.82,0.82,0.95,0.95);
 	legfit->AddEntry(fModelTotTh, "Total Th-232", "l");
  	legfit->AddEntry(fModelTotRa, "Total Ra-226", "l");
 	legfit->AddEntry(fModelTotK, "Total K-40", "l");
 	legfit->AddEntry(fModelTotCo, "Total Co-60", "l");

 	legfit->Draw();


 	TCanvas *ctable = new TCanvas("ctable", "ctable", 800, 1200);
 	// TPaveText *pt = new TPaveText(.0,.0,1.,1.);
 	// pt->AddText("Fit Parameters");
 	// pt->AddBox(.0, .0, 1.0, .6);
 	// pt->AddText("Blah");
 	// pt->AddText("Blah");
 	// pt->AddText("Blah");
 	// pt->AddText("Blah");

 	// pt->AddLine(.0,.65,1.,.65);

 	if(bToyFit)
 	{
 		pt->AddText(Form("Toy Data Integral (Tl-208 peak): %.2f", fToyData->Integral(2600/dBinSize, 2700/dBinSize))); 	
 		pt->AddText(Form("Toy Data Integral (Fit Range %.0f to %.0f): %.2f", dFitMin, dFitMax, fToyData->Integral(dFitMin/dBinSize, dFitMax/dBinSize)));
 	}
 	else
 	{
 		pt->AddText(Form("Data Integral (Tl-208 peak): %.2f", fDataHistoM1->Integral(2600/dBinSize, 2700/dBinSize)));
	 	pt->AddText(Form("Data Integral (Fit Range %.0f to %.0f): %.2f", dFitMin, dFitMax, fDataHistoM1->Integral(dFitMin/dBinSize, dFitMax/dBinSize)));
 	}

 	pt->AddText(Form("Model Integral (Tl-208 peak): %.2f", fModelTot->Integral(2600/dBinSize, 2700/dBinSize)));
 	pt->AddText(Form("Model Integral (Fit Range %.0f to %.0f): %.2f", dFitMin, dFitMax, fModelTot->Integral(dFitMin/dBinSize, dFitMax/dBinSize)));
 	pt->Draw();


/*
	// Residuals
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
*/

	// gResidual->Draw("AP");

	return true;
   
 }


// Draws all MC spectra
 void TBackgroundModel::DrawMC()
 {
 	gStyle->SetOptStat(0);
 	gStyle->SetOptTitle(0);

 	TLegend *legth = new TLegend(0.75,0.80,0.925,0.925);
 	TLegend *legra = new TLegend(0.75,0.80,0.925,0.925);
 	TLegend *legk = new TLegend(0.75,0.80,0.925,0.925);
 	TLegend *legco = new TLegend(0.75,0.80,0.925,0.925);


    TCanvas *cTh232 = new TCanvas("cTh232", "cTh232", 1200, 800);
    cTh232->SetLogy();

    fModelFrameTh->SetLineColor(1);
    fModelTShieldTh->SetLineColor(2);
    fModel50mKTh->SetLineColor(3);
    fModel600mKTh->SetLineColor(4);
    fModelIVCTh->SetLineColor(6);
    fModelOVCTh->SetLineColor(7);


    fModelFrameTh->Draw();
    fModelTShieldTh->Draw("SAME");
    fModel50mKTh->Draw("SAME");
    fModel600mKTh->Draw("SAME");
    fModelIVCTh->Draw("SAME");
    fModelOVCTh->Draw("SAME");
    fModelFrameTh->GetXaxis()->SetRange(2200/dBinSize, 2700/dBinSize);


    legth->AddEntry(fModelFrameTh, "Frame" ,"l");
    legth->AddEntry(fModelTShieldTh, "TShield", "l");
    legth->AddEntry(fModel50mKTh, "50mK" ,"l");
    legth->AddEntry(fModel600mKTh, "600mK" ,"l");
    legth->AddEntry(fModelIVCTh, "IVC" ,"l");
    legth->AddEntry(fModelOVCTh, "OVC" ,"l");
    legth->Draw();

    TCanvas *cRa226 = new TCanvas("cRa226", "cRa226", 1200, 800);
    cRa226->SetLogy();

    fModelFrameRa->SetLineColor(1);
    fModelTShieldRa->SetLineColor(2);
    fModel50mKRa->SetLineColor(3);
    fModel600mKRa->SetLineColor(4);
    fModelIVCRa->SetLineColor(6);
    fModelOVCRa->SetLineColor(7);

    fModelFrameRa->Draw();
    fModelTShieldRa->Draw("SAME");
    fModel50mKRa->Draw("SAME");
    fModel600mKRa->Draw("SAME");
    fModelIVCRa->Draw("SAME");
    fModelOVCRa->Draw("SAME");
    fModelFrameRa->GetXaxis()->SetRange(2200/dBinSize, 2700/dBinSize);


    legra->AddEntry(fModelFrameRa, "Frame" ,"l");
    legra->AddEntry(fModelTShieldRa, "TShield", "l");
    legra->AddEntry(fModel50mKRa, "50mK" ,"l");
    legra->AddEntry(fModel600mKRa, "600mK" ,"l");
    legra->AddEntry(fModelIVCRa, "IVC" ,"l");
    legra->AddEntry(fModelOVCRa, "OVC" ,"l");
    legra->Draw();

    TCanvas *cK40 = new TCanvas("cK40", "cK40", 1200, 800);
    cK40->SetLogy();

    fModelFrameK->SetLineColor(1);
    fModelTShieldK->SetLineColor(2);
    fModel50mKK->SetLineColor(3);
    fModel600mKK->SetLineColor(4);
    fModelIVCK->SetLineColor(6);
    fModelOVCK->SetLineColor(7);

    fModelFrameK->Draw();
    fModelTShieldK->Draw("SAME");
    fModel50mKK->Draw("SAME");
    fModel600mKK->Draw("SAME");
    fModelIVCK->Draw("SAME");
    fModelOVCK->Draw("SAME");
    fModelFrameK->GetXaxis()->SetRange(1000/dBinSize, 1600/dBinSize);


    legk->AddEntry(fModelFrameK, "Frame" ,"l");
    legk->AddEntry(fModelTShieldK, "TShield", "l");
    legk->AddEntry(fModel50mKK, "50mK" ,"l");
    legk->AddEntry(fModel600mKK, "600mK" ,"l");
    legk->AddEntry(fModelIVCK, "IVC" ,"l");
    legk->AddEntry(fModelOVCK, "OVC" ,"l");    
    legk->Draw();

    TCanvas *cCo60 = new TCanvas("cCo60", "cCo60", 1200, 800);
    cCo60->SetLogy();

    fModelFrameCo->SetLineColor(1);
    fModelTShieldCo->SetLineColor(2);
    fModel50mKCo->SetLineColor(3);
    fModel600mKCo->SetLineColor(4);
    fModelIVCCo->SetLineColor(6);
    fModelOVCCo->SetLineColor(7);

    fModelFrameCo->Draw();
    fModelTShieldCo->Draw("SAME");
    fModel50mKCo->Draw("SAME");
    fModel600mKCo->Draw("SAME");
    fModelIVCCo->Draw("SAME");
    fModelOVCCo->Draw("SAME");
    fModelFrameCo->GetXaxis()->SetRange(2200/dBinSize, 2700/dBinSize);


    legco->AddEntry(fModelFrameCo, "Frame" ,"l");
    legco->AddEntry(fModelTShieldCo, "TShield", "l");
    legco->AddEntry(fModel50mKCo, "50mK" ,"l");
    legco->AddEntry(fModel600mKCo, "600mK" ,"l");
    legco->AddEntry(fModelIVCCo, "IVC" ,"l");
    legco->AddEntry(fModelOVCCo, "OVC" ,"l");
    legco->Draw();

 }


// Generates toy data using MC histogram
void TBackgroundModel::GenerateToyData()
{

	bToyFit = true;
	// Create some RNG for weights?
	// Currently just put in by hand
	std::string dToyDir = "/Users/brian/macros/Simulations/Bkg/";

	// Load g2tas smeared data:
    outTreeToyTh 		= LoadMC(dToyDir.c_str(),	"Frame", 	"Th232", 1);
    outTreeToyRa 		= LoadMC(dToyDir.c_str(),	"Frame", 	"Ra226", 1);
    outTreeToyCo 		= LoadMC(dToyDir.c_str(),	"Frame", 	"Co60", 1);
    outTreeToyK 		= LoadMC(dToyDir.c_str(),	"Frame", 	"K40", 1);

	outTreeToyTh->Project("fToyDataTh", 	"Ener1", ener_cut);
	outTreeToyRa->Project("fToyDataRa", 	"Ener1", ener_cut);
	outTreeToyCo->Project("fToyDataCo", 	"Ener1", ener_cut);
	outTreeToyK->Project("fToyDataK",	 	"Ener1", ener_cut);

	NormalizePDF(fToyDataTh, 	dFitMin, dFitMax);
	NormalizePDF(fToyDataRa, 	dFitMin, dFitMax);
	NormalizePDF(fToyDataCo, 	dFitMin, dFitMax);
	NormalizePDF(fToyDataK, 	dFitMin, dFitMax);



	fToyData->Add(fToyDataTh,	2000);
	fToyData->Add(fToyDataRa,	550);
	fToyData->Add(fToyDataCo,	600);
	fToyData->Add(fToyDataK,	700);

}


// Calculates ChiSquare... model parameters not set here!
double TBackgroundModel::GetChiSquare()
{
	//cout<<"Calling GetChiSquare()"<<endl;
	double chiSquare = 0.;
	double data_i, err_i, model_i;	

	for(int i = dFitMin/dBinSize; i < dFitMax/dBinSize; i++)
	{
		if(bToyFit)
		{
			data_i = fToyData->GetBinContent(i); // For Toy data
		}
		else 
		{
			data_i = fDataHistoM1->GetBinContent(i); // For real data
		}

		model_i = fModelTot->GetBinContent(i);


		// Log-likelihood Chi-Squared

		if(model_i != 0 && data_i != 0)
		{
			chiSquare += 2 * (model_i - data_i + data_i * TMath::Log(data_i/model_i));
		}


		// Neyman chi-squared
/*		
		err_i = sqrt(data_i);	// Assuming no data bins are 0!
		if(err_i>0)
		{
			chiSquare += pow(data_i - model_i,2)/pow(err_i,2);
		}
		else
		{

			chiSquare+=pow(data_i - model_i,2);	
		}
*/


	}

	return chiSquare;
}


void TBackgroundModel::Initialize()
{	

	dDataDir = 	"/Users/brian/macros/Simulations/Bkg/Unsmeared/";
	dDataIntegral = 0;
	bToyFit = false;

	// Bin size (keV)
	dBinSize = 10;
	// Histogram range
	dMaxEnergy = 3500.;
	dMinEnergy = 0.;

	// Fitting range
	dFitMin = 0.;
	dFitMax = 2700.;


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
	fToyData		 = new TH1D("fToyData",			"", dNBins, dMinEnergy, dMaxEnergy);


	// Model histograms

	fModelFrameTh	 = new TH1D("fModelFrameTh", 	"Frame",		dNBins, dMinEnergy, dMaxEnergy);
	fModelTShieldTh	 = new TH1D("fModelTShieldTh", 	"TShield",		dNBins, dMinEnergy, dMaxEnergy);
	fModel50mKTh	 = new TH1D("fModel50mKTh", 	"50mK",			dNBins, dMinEnergy, dMaxEnergy);
	fModel600mKTh	 = new TH1D("fModel600mKTh", 	"600mK",		dNBins, dMinEnergy, dMaxEnergy);
	fModelIVCTh		 = new TH1D("fModelIVCTh", 		"IVC", 			dNBins, dMinEnergy, dMaxEnergy);
	fModelOVCTh		 = new TH1D("fModelOVCTh",		"OVC",	 		dNBins, dMinEnergy, dMaxEnergy);

	fModelFrameRa	 = new TH1D("fModelFrameRa", 	"Frame", 		dNBins, dMinEnergy, dMaxEnergy);
	fModelTShieldRa	 = new TH1D("fModelTShieldRa", 	"TShield",		dNBins, dMinEnergy, dMaxEnergy);	
	fModel50mKRa	 = new TH1D("fModel50mKRa", 	"50mK", 		dNBins, dMinEnergy, dMaxEnergy);
	fModel600mKRa	 = new TH1D("fModel600mKRa", 	"600mK", 		dNBins, dMinEnergy, dMaxEnergy);
	fModelIVCRa		 = new TH1D("fModelIVCRa", 		"IVC", 			dNBins, dMinEnergy, dMaxEnergy);
	fModelOVCRa		 = new TH1D("fModelOVCRa", 		"OVC", 			dNBins, dMinEnergy, dMaxEnergy);

	fModelFrameK	 = new TH1D("fModelFrameK", 	"Frame", 		dNBins, dMinEnergy, dMaxEnergy);
	fModelTShieldK	 = new TH1D("fModelTShieldK", 	"TShield",		dNBins, dMinEnergy, dMaxEnergy);
	fModel50mKK	 	 = new TH1D("fModel50mKK",	 	"50mK",			dNBins, dMinEnergy, dMaxEnergy);
	fModel600mKK	 = new TH1D("fModel600mKK", 	"600mK",		dNBins, dMinEnergy, dMaxEnergy);
	fModelIVCK		 = new TH1D("fModelIVCK", 		"IVC", 			dNBins, dMinEnergy, dMaxEnergy);
	fModelOVCK		 = new TH1D("fModelOVCK",		"OVC",	 		dNBins, dMinEnergy, dMaxEnergy);

	fModelFrameCo	 = new TH1D("fModelFrameCo", 	"Frame", 		dNBins, dMinEnergy, dMaxEnergy);
	fModelTShieldCo	 = new TH1D("fModelTShieldCo", 	"TShield",		dNBins, dMinEnergy, dMaxEnergy);
	fModel50mKCo	 = new TH1D("fModel50mKCo", 	"50mK",			dNBins, dMinEnergy, dMaxEnergy);
	fModel600mKCo	 = new TH1D("fModel600mKCo", 	"600mK",		dNBins, dMinEnergy, dMaxEnergy);
	fModelIVCCo		 = new TH1D("fModelIVCCo", 		"IVC", 			dNBins, dMinEnergy, dMaxEnergy);
	fModelOVCCo		 = new TH1D("fModelOVCCo",		"OVC",	 		dNBins, dMinEnergy, dMaxEnergy);

	// Total model histograms
	fModelTot 		 = new TH1D("fModelTot", 		"Frame", 		dNBins, dMinEnergy, dMaxEnergy);	
	fModelTotTh		 = new TH1D("fModelTotTh", 		"Total Th232", 	dNBins, dMinEnergy, dMaxEnergy);
	fModelTotRa		 = new TH1D("fModelTotRa", 		"Total Ra226", 	dNBins, dMinEnergy, dMaxEnergy);
	fModelTotK		 = new TH1D("fModelTotK", 		"Total K40",	dNBins, dMinEnergy, dMaxEnergy);
	fModelTotCo		 = new TH1D("fModelTotCo", 		"Total Co60", 	dNBins, dMinEnergy, dMaxEnergy);


	// Smearing
	gaus = new TF1("gaus","gaus(0)", dMinEnergy, dMaxEnergy);

	// Clear Initial Parameters
	fParameters[0] 	= 0.;
	fParameters[1] 	= 0.;
	fParameters[2] 	= 0.;
	fParameters[3] 	= 0.;
	fParameters[4] 	= 0.;
	fParameters[5] 	= 0.;
	fParameters[6] 	= 0.;
	fParameters[7] 	= 0.;
	fParameters[8]	= 0.;
/*	
	fParameters[8]	= 0.;
	fParameters[9] 	= 0.;
	fParameters[10] = 0.;
	fParameters[11] = 0.;
	fParameters[12] = 0.;
	fParameters[13] = 0.;
	fParameters[14] = 0.;
	fParameters[15] = 0.;
	fParameters[16] = 0.;
	fParameters[17] = 0.;
	fParameters[18] = 0.;
	fParameters[19] = 0.;
	fParameters[20]	= 0.;
	fParameters[21] = 0.;
	fParameters[22] = 0.;
	fParameters[23] = 0.;
*/

	// Loading all data in Initialize, correct or no?
	LoadData();	
	ReadMC();

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

	// Normalizing data (don't!)
	// bin 0 = underflow, bin dNBins = last bin with upper-edge xup Excluded
	dDataIntegral = fDataHistoM1->Integral(1, dNBins);
	cout << "Events in background spectrum: " << dDataIntegral << endl;
	// cout << "Normalized Data" << endl;
}


// Loads MC files into Trees
TChain *TBackgroundModel::LoadMC(std::string dDir, std::string dLocation, std::string dSource, int dMult)
{
    TChain *outTree = new TChain("outTree");
    outTree->Add(Form("%s%s-%s-B-M%d-T50-r0.0425.root", dDir.c_str(), dLocation.c_str(), dSource.c_str(), dMult));

    return outTree;
}


// Normalize histogram
void TBackgroundModel::NormalizePDF(TH1D *h1, int minE, int maxE)
{
	double dIntegral;

	// bin 0 = underflow, bin dNBins = last bin with upper-edge xup Excluded
	dIntegral = h1->Integral(minE/dBinSize, maxE/dBinSize);
	// cout << "Integral for " << h1->GetTitle() << " :" << dIntegral << endl;

	// Make sure integral isn't 0 --> Need to double check if this is the right thing to do!
	if(dIntegral != 0)
	{
		h1->Scale(1/dIntegral);
	}
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
/*	
	cout<< "Par8 = "	<< fParameters[8]	<< endl;
	cout<< "Par9 = "	<< fParameters[9]	<< endl;
	cout<< "Par10 = "	<< fParameters[10]	<< endl;
	cout<< "Par11 = "	<< fParameters[11] 	<< endl;
	cout<< "Par12 = "	<< fParameters[12]	<< endl;
	cout<< "Par13 = "	<< fParameters[13]	<< endl;
	cout<< "Par14 = "	<< fParameters[14]	<< endl;
	cout<< "Par15 = "	<< fParameters[15]	<< endl;
	cout<< "Par16 = "	<< fParameters[16]	<< endl;
	cout<< "Par17 = "	<< fParameters[17]	<< endl;
	cout<< "Par18 = "	<< fParameters[18]	<< endl;
	cout<< "Par19 = "	<< fParameters[19]	<< endl;
	cout<< "Par20 = "	<< fParameters[20]	<< endl;
	cout<< "Par21 = "	<< fParameters[21]	<< endl;
	cout<< "Par22 = "	<< fParameters[22]	<< endl;
	cout<< "Par23 = "	<< fParameters[23] 	<< endl;
*/
	double dSum = fParameters[0] + fParameters[1] + fParameters[2] + fParameters[3]
					+ fParameters[4] + fParameters[5] + fParameters[6] + fParameters[7];
					// + fParameters[8] + fParameters[9] + fParameters[10] + fParameters[11]
					// + fParameters[12] + fParameters[13] + fParameters[14] + fParameters[15]
					// + fParameters[16] + fParameters[17] + fParameters[18] + fParameters[19]
					// + fParameters[20] + fParameters[21] + fParameters[22] + fParameters[23];

	// cout << "Par11 (1 - Sum) = " << 1 - dSum << endl;
	cout << "Sum = " << dSum << endl; 


}

// Fills MC histograms
void TBackgroundModel::ReadMC()
{

    // Loads MC data
    outTreeFrameTh 		= LoadMC(dDataDir.c_str(),	"Frame", 	"Th232", 1);
    outTreeTShieldTh 	= LoadMC(dDataDir.c_str(),	"TShield", 	"Th232", 1);
    outTree50mKTh 		= LoadMC(dDataDir.c_str(),	"50mK",		"Th232", 1);
    outTree600mKTh 		= LoadMC(dDataDir.c_str(),	"600mK", 	"Th232", 1);
    outTreeIVCTh 		= LoadMC(dDataDir.c_str(),	"IVC", 		"Th232", 1);
    outTreeOVCTh 		= LoadMC(dDataDir.c_str(),	"OVC", 		"Th232", 1);

    outTreeFrameRa	 	= LoadMC(dDataDir.c_str(),	"Frame", 	"Ra226", 1);
    outTreeTShieldRa 	= LoadMC(dDataDir.c_str(),	"TShield", 	"Ra226", 1);    
    outTree50mKRa	 	= LoadMC(dDataDir.c_str(),	"50mK", 	"Ra226", 1);
    outTree600mKRa		= LoadMC(dDataDir.c_str(),	"600mK", 	"Ra226", 1);
    outTreeIVCRa	 	= LoadMC(dDataDir.c_str(),	"IVC", 		"Ra226", 1);
    outTreeOVCRa 		= LoadMC(dDataDir.c_str(),	"OVC", 		"Ra226", 1);

    outTreeFrameK 		= LoadMC(dDataDir.c_str(),	"Frame", 	"K40", 1);
    outTreeTShieldK 	= LoadMC(dDataDir.c_str(),	"TShield", 	"K40", 1);    
    outTree50mKK	 	= LoadMC(dDataDir.c_str(),	"50mK", 	"K40", 1);
    outTree600mKK		= LoadMC(dDataDir.c_str(),	"600mK", 	"K40", 1);
    outTreeIVCK		 	= LoadMC(dDataDir.c_str(),	"IVC", 		"K40", 1);
    outTreeOVCK 		= LoadMC(dDataDir.c_str(),	"OVC", 		"K40", 1);


    outTreeFrameCo 		= LoadMC(dDataDir.c_str(),	"Frame", 	"Co60",	1);
    outTreeTShieldCo 	= LoadMC(dDataDir.c_str(),	"TShield", 	"Co60", 1);    
    outTree50mKCo	 	= LoadMC(dDataDir.c_str(),	"50mK", 	"Co60", 1);
    outTree600mKCo		= LoadMC(dDataDir.c_str(),	"600mK", 	"Co60", 1);
    outTreeIVCCo	 	= LoadMC(dDataDir.c_str(),	"IVC", 		"Co60", 1);
    outTreeOVCCo 		= LoadMC(dDataDir.c_str(),	"OVC", 		"Co60", 1);


	outTreeFrameTh->Project("fModelFrameTh", 		"Ener1", ener_cut);
	outTreeTShieldTh->Project("fModelTShieldTh",	"Ener1", ener_cut);
    outTree50mKTh->Project("fModel50mKTh", 			"Ener1", ener_cut);
    outTree600mKTh->Project("fModel600mKTh", 		"Ener1", ener_cut);
    outTreeIVCTh->Project("fModelIVCTh", 			"Ener1", ener_cut);
    outTreeOVCTh->Project("fModelOVCTh", 			"Ener1", ener_cut);

	outTreeFrameRa->Project("fModelFrameRa", 		"Ener1", ener_cut);
	outTreeTShieldRa->Project("fModelTShieldRa",	"Ener1", ener_cut);	
    outTree50mKRa->Project("fModel50mKRa", 			"Ener1", ener_cut);
    outTree600mKRa->Project("fModel600mKRa", 		"Ener1", ener_cut);
    outTreeIVCRa->Project("fModelIVCRa", 			"Ener1", ener_cut);
    outTreeOVCRa->Project("fModelOVCRa", 			"Ener1", ener_cut);

	outTreeFrameK->Project("fModelFrameK", 			"Ener1", ener_cut);
	outTreeTShieldK->Project("fModelTShieldK",		"Ener1", ener_cut);	
    outTree50mKK->Project("fModel50mKK", 			"Ener1", ener_cut);
    outTree600mKK->Project("fModel600mKK", 			"Ener1", ener_cut);
    outTreeIVCK->Project("fModelIVCK", 				"Ener1", ener_cut);
    outTreeOVCK->Project("fModelOVCK", 				"Ener1", ener_cut);	

	outTreeFrameCo->Project("fModelFrameCo", 		"Ener1", ener_cut);
	outTreeTShieldCo->Project("fModelTShieldCo",	"Ener1", ener_cut);	
    outTree50mKCo->Project("fModel50mKCo", 			"Ener1", ener_cut);
    outTree600mKCo->Project("fModel600mKCo", 		"Ener1", ener_cut);
    outTreeIVCCo->Project("fModelIVCCo", 			"Ener1", ener_cut);
    outTreeOVCCo->Project("fModelOVCCo", 			"Ener1", ener_cut);

	cout << "Loaded MC" << endl;

	// Normalize all MC histograms
	NormalizePDF(fModelFrameTh, 	dFitMin, dFitMax);
	NormalizePDF(fModelTShieldTh, 	dFitMin, dFitMax);
	NormalizePDF(fModel50mKTh, 		dFitMin, dFitMax);
	NormalizePDF(fModel600mKTh, 	dFitMin, dFitMax);
	NormalizePDF(fModelIVCTh, 		dFitMin, dFitMax);
	NormalizePDF(fModelOVCTh, 		dFitMin, dFitMax);

	NormalizePDF(fModelFrameRa, 	dFitMin, dFitMax);
	NormalizePDF(fModelTShieldRa, 	dFitMin, dFitMax);	
	NormalizePDF(fModel50mKRa, 		dFitMin, dFitMax);
	NormalizePDF(fModel600mKRa, 	dFitMin, dFitMax);
	NormalizePDF(fModelIVCRa, 		dFitMin, dFitMax);
	NormalizePDF(fModelOVCRa, 		dFitMin, dFitMax);

	// Normalizing K-40 for full range since no peaks above 1500 keV
	NormalizePDF(fModelFrameK, 		dFitMin, dFitMax);
	NormalizePDF(fModelTShieldK, 	dFitMin, dFitMax);	
	NormalizePDF(fModel50mKK, 		dFitMin, dFitMax);
	NormalizePDF(fModel600mKK, 		dFitMin, dFitMax);
	NormalizePDF(fModelIVCK, 		dFitMin, dFitMax);
	NormalizePDF(fModelOVCK, 		dFitMin, dFitMax);

	NormalizePDF(fModelFrameCo, 	dFitMin, dFitMax);
	NormalizePDF(fModelTShieldCo, 	dFitMin, dFitMax);	
	NormalizePDF(fModel50mKCo, 		dFitMin, dFitMax);
	NormalizePDF(fModel600mKCo, 	dFitMin, dFitMax);
	NormalizePDF(fModelIVCCo, 		dFitMin, dFitMax);
	NormalizePDF(fModelOVCCo, 		dFitMin, dFitMax);


	cout << "Normalized MC PDFs" << endl;
}

// Set Parameters in Model
void TBackgroundModel::SetParameters(int index, double value)
{
	// Change the index max depending on model
	if(index > 24) cout << "Index too large" << endl;
	else fParameters[index] = value;

}


// For custom smearing
TH1D *TBackgroundModel::SmearMC(TH1D *hMC, double resolution)
{

	TH1D *hSmeared = new TH1D("hSmeared", "", dNBins, dMinEnergy, dMaxEnergy);

	double dArea;
	double dSmearedValue;

	// If i only goes through fit range, saves time?
	for(int i = 0; i<dNBins; i++)
	{
		for(int j = 0; j<dNBins; j++)
		{
			// Normalization of gaussian = (bin size * Area of bin j in MC) / Sigma of bin j (fit function evaluated at bin center)
			dArea = dBinSize*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*resolution);

			// Set parameters of gaussian ... resolution floating in fit
			gaus->SetParameters(dArea, hMC->GetBinCenter(j), resolution);

			// Smeared contribution from gaussian centered at bin j for bin i 
			dSmearedValue = gaus->Eval(hSmeared->GetBinCenter(i));

			// Fill bin i with contribution from gaussian centered at bin j
			hSmeared->Fill(hSmeared->GetBinCenter(i), dSmearedValue);
		}
	}
	
	return hSmeared;
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
	fModelTot->Reset();

	// Create model


	////////////////////////////////////////
	// Few parameters
	////////////////////////////////////////
	fModelTot->Add(SmearMC(fModelFrameTh,	fParameters[8]), fParameters[0]);
	fModelTot->Add(SmearMC(fModelTShieldTh,	fParameters[8]), fParameters[0]);	
	fModelTot->Add(SmearMC(fModel50mKTh,	fParameters[8]), fParameters[0]);
	fModelTot->Add(SmearMC(fModel600mKTh,	fParameters[8]), fParameters[0]);
	fModelTot->Add(SmearMC(fModelIVCTh, 	fParameters[8]), fParameters[1]);
	fModelTot->Add(SmearMC(fModelOVCTh, 	fParameters[8]), fParameters[1]);

	fModelTot->Add(SmearMC(fModelFrameRa,	fParameters[8]), fParameters[2]);
	fModelTot->Add(SmearMC(fModelTShieldRa,	fParameters[8]), fParameters[2]);	
	fModelTot->Add(SmearMC(fModel50mKRa, 	fParameters[8]), fParameters[2]);
	fModelTot->Add(SmearMC(fModel600mKRa, 	fParameters[8]), fParameters[2]);
	fModelTot->Add(SmearMC(fModelIVCRa, 	fParameters[8]), fParameters[3]);
	fModelTot->Add(SmearMC(fModelOVCRa, 	fParameters[8]), fParameters[3]);

	fModelTot->Add(SmearMC(fModelFrameK, 	fParameters[8]), fParameters[4]);
	fModelTot->Add(SmearMC(fModelTShieldK, 	fParameters[8]), fParameters[4]);
	fModelTot->Add(SmearMC(fModel50mKK, 	fParameters[8]), fParameters[4]);
	fModelTot->Add(SmearMC(fModel600mKK, 	fParameters[8]), fParameters[4]);
	fModelTot->Add(SmearMC(fModelIVCK, 		fParameters[8]), fParameters[5]);
	fModelTot->Add(SmearMC(fModelOVCK, 		fParameters[8]), fParameters[5]);

	fModelTot->Add(SmearMC(fModelFrameCo, 	fParameters[8]), fParameters[6]);
	fModelTot->Add(SmearMC(fModelTShieldCo, fParameters[8]), fParameters[6]);
	fModelTot->Add(SmearMC(fModel50mKCo, 	fParameters[8]), fParameters[6]);
	fModelTot->Add(SmearMC(fModel600mKCo, 	fParameters[8]), fParameters[6]);
	fModelTot->Add(SmearMC(fModelIVCCo, 	fParameters[8]), fParameters[7]);
	fModelTot->Add(SmearMC(fModelOVCCo, 	fParameters[8]), fParameters[7]);	


	////////////////////////////////////////
	// All Parameters ... probably won't use this
	////////////////////////////////////////
/*
	fModelTot->Add(fModelFrameTh,	fParameters[0]);
	fModelTot->Add(fModelTShieldTh,	fParameters[1]);	
	fModelTot->Add(fModel50mKTh,	fParameters[2]);
	fModelTot->Add(fModel600mKTh,	fParameters[3]);
	fModelTot->Add(fModelIVCTh,		fParameters[4]);
	fModelTot->Add(fModelOVCTh,		fParameters[5]);

	fModelTot->Add(fModelFrameRa,	fParameters[6]);
	fModelTot->Add(fModelTShieldRa,	fParameters[7]);	
	fModelTot->Add(fModel50mKRa,	fParameters[8]);
	fModelTot->Add(fModel600mKRa,	fParameters[9]);
	fModelTot->Add(fModelIVCRa,		fParameters[10]);
	fModelTot->Add(fModelOVCRa,		fParameters[11]);

	fModelTot->Add(fModelFrameK,	fParameters[12]);
	fModelTot->Add(fModelTShieldK,	fParameters[13]);
	fModelTot->Add(fModel50mKK,		fParameters[14]);
	fModelTot->Add(fModel600mKK,	fParameters[15]);
	fModelTot->Add(fModelIVCK,		fParameters[16]);
	fModelTot->Add(fModelOVCK,		fParameters[17]);

	fModelTot->Add(fModelFrameCo,	fParameters[18]);
	fModelTot->Add(fModelTShieldCo,	fParameters[19]);
	fModelTot->Add(fModel50mKCo,	fParameters[20]);
	fModelTot->Add(fModel600mKCo,	fParameters[21]);
	fModelTot->Add(fModelIVCCo,		fParameters[22]);
	fModelTot->Add(fModelOVCCo,		fParameters[23]);	
*/


	// Don't use this
	// fModelTot->Add(fModelFrameCo, 	1 - (fParameters[0] + fParameters[1] + fParameters[2] + fParameters[3] + fParameters[4] + fParameters[5]
									// + fParameters[6] + fParameters[7] + fParameters[8] + fParameters[9] + fParameters[10]));

/*
	// Test gaussian
	TF1	gaus("mygaus","gaus(0)",-10,10);
	gaus.SetParameters(1,fParameters[0],fParameters[1]);

	for(int i =1;i<fModelHisto->GetNbinsX(); i++)
	{
		fModelHisto->SetBinContent(i,gaus.Eval(fModelHisto->GetBinCenter(i)));
	}
	
	fModelHisto->Scale(fDataHisto->Integral()/fModelHisto->Integral());
	
*/	


}