/*
	This macro compares the Calibration data with the MC spectra
	This macro is for g2tas ROOT output files (with outTree)
*/

#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TChain.h"
#include "TCut.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TAxis.h"

#include <iostream>
#include <vector>


TH1D *SmearMC(TH1D *hMC, TH1D *hSMC, double resolution1, double resolution2)
{
	// Reset previously Modeled histogram
	hSMC->Reset();

	double dNorm;
    double dNorm2;
	double dSmearedValue;

	TF1 *gaus1 = new TF1("gaus","gaus(0)", 0, 3500);
	TF1 *gaus2 = new TF1("gaus","gaus(0)", 0, 3500);


	for(int i = 0; i<1750; i++)
	{
		for(int j = 0; j<1750; j++)
		{
			// Normalization of gaussian = (bsin size * Area of bin j in MC) / Sigma of bin j (fit function evaluated at bin center)
			dNorm = (853.3/1215.8)*2*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*resolution1);
    	    dNorm2 = (362.5/1215.8)*2*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*resolution2);

			// Set parameters of gaussian ... 2nd gaussian *slightly* shifted... not sure if this works
			gaus->SetParameters(dNorm, hMC->GetBinCenter(j), resolution1);
     	    gaus2->SetParameters(dNorm2, hMC->GetBinCenter(j)-1, resolution2);

			// Smeared contribution from gaussian centered at bin j for bin i 
			dSmearedValue = gaus->Eval(hSMC->GetBinCenter(i)) + gaus2->Eval(hSMC->GetBinCenter(i));

			// Fill bin i with contribution from gaussian centered at bin j
			hSMC->Fill(hSMC->GetBinCenter(i), dSmearedValue);
		}
	}
	
	return hSMC;
}



void FitCalib()
{

	gStyle->SetOptFit();

	TChain *qtree = new TChain("qredtree", "qredtree");
	qtree->Add("/Users/brian/macros/CUOREZ/Bkg/Q0_DR2_CalibrationSignalData.root");


	TH1D *hCal = new TH1D("hCal", "", 15, 2600, 2630);

	qtree->Project("hCal","Energy", "Multiplicity == 1");


	// TF1 *f1 = new TF1("f1", "gaus(0)+gaus(3)");
	TF1 *f1 = new TF1("f1", "gaus(0)");
	f1->SetParameter(0, 100);
	f1->SetParameter(1, 2615);
	f1->SetParameter(2, 2);
	// f1->SetParameter(3, 10);
	// f1->SetParameter(4, 2612);
	// f1->SetParameter(5, 4);



	// hCal->SetTitle("Double Gaussian Fit");
	hCal->GetXaxis()->SetTitle("Energy (keV)");
	hCal->GetYaxis()->SetTitle("Counts/keV");
	hCal->Fit("f1");
	hCal->Draw();

}


// Choose multiplicity, saving plots, and energy range for plots
void CalibrationBenchmark(int dMult = 1, bool bSavePlot = false, double dEMin = 0, double dEMax = 3500)
{
//	gStyle->Reset();
	gStyle->SetOptStat("emr");
	gStyle->SetOptFit();
	gStyle->SetStatW(0.15);
	gStyle->SetStatH(0.2);
//	gStyle->SetOptTitle(0); 	// Gets rid of titles


	TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
	int nHisto = 2;  
	double width1 = 0.025;
	double width2 = 0.975;
	double canBotMargin = 0.025;
	double canTopMargin = 0.025;
	double padHeight = (1.-canTopMargin-canBotMargin)/nHisto;
  
	TPad *p1 = new TPad("p1","p1",width1,canBotMargin,width2,canBotMargin+padHeight,0,0);
	p1->SetTopMargin(0);
	p1->SetLeftMargin(0.075);
	p1->SetRightMargin(0.075);
	p1->SetFillColor(0);
	p1->SetBorderMode(0);
	p1->SetBorderSize(0);
	p1->Draw();
  
	// p2 is on top  
	TPad* p2 = new TPad("p2","p2",width1,canBotMargin+padHeight,width2,canBotMargin+2*padHeight,0,0);
	p2->SetBottomMargin(0);
	p2->SetTopMargin(0.1);
	p2->SetLeftMargin(0.075);
	p2->SetRightMargin(0.075);
	p2->SetFillColor(0);
	p2->SetBorderMode(0);
	p2->SetBorderSize(0);
	p2->Draw();


	// Set Histogram min/max and bin size
	int dBinSize = 5;
	int dBin = (dEMax - dEMin)/dBinSize;

	// Set up range for residuals plot
	int dRebin = 5; // Rebin value (for residuals)
	int dEMinRes, dEMaxRes;
	if (dEMin < 300)
	{
		dEMinRes = 300;
	}
	else 
	{
		dEMinRes = dEMin;
	}
	if (dEMax > 2700)
	{
		dEMaxRes = 2700;
	}
	else
	{
		dEMaxRes = dEMax;
	}

	// Cuts and stuff
	double dMC, dMC2;
	double dCal, dCal2;
	double dMCFull;
	double dCalFull;
	double dCutMin 		= 	2590; // first normalization
	double dCutMax 		= 	2630;
	double dCutMin2 	= 	1000; // second normalization
	double dCutMax2 	= 	2000;
	double Normalization;
	double Normalization2;
	double dRatio, dRatio2;

	int dThres = 50;



	TChain *qtree;
	TChain *qtree_bkg;
	TChain *outTree;
	TH1D *hMC;
	TH1D *hCal;
	TH1D *hBkg;
	TLegend *leg;
	leg = new TLegend(0.72,0.55,0.925,0.9);

///////////////////////////////////////////////////////////////////////
// Load Calibration Data
///////////////////////////////////////////////////////////////////////

	// Select 1st pad
	p2->cd();
	p2->SetLogy();

	qtree = new TChain("qredtree", "qredtree");
	qtree->Add("/Users/brian/macros/CUOREZ/Bkg/Q0_DR2_CalibrationSignalData.root");
	// qtree->Add("/Users/brian/macros/CUOREZ/Bkg/ReducedCalib-ds2061.root"); // for 50 keV
	// qtree->Add(Form("/Users/brian/data/CUOREZ/BlindedReduced_200730_C_p001_T%d.root", dThres));
	// qtree->Add(Form("/Users/brian/data/CUOREZ/BlindedReduced_200757_C_p001_T%d.root", dThres));


	// qtree_bkg = new TChain("qtree","qtree");
	// qtree_bkg->Add("/Users/brian/data/CUOREZ/BlindedReduced_200713_B_p001.root");

	// TCut base_cut("IsSignal");
	TCut base_cut;
	// base_cut = base_cut && "Filter_RejectBadIntervals";
	// base_cut = base_cut && "NumberOfPulses==1";
	// base_cut = base_cut && "TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0";
	// base_cut = base_cut && "TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0";
	// base_cut = base_cut && "abs(BaselineSlope)<0.1";
	// base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";
	// base_cut = base_cut && "Run == 200730 || Run == 200757";

	// TCut mult_cut(Form("Multiplicity_OFTime==%d", dMult));
	TCut mult_cut(Form("Multiplicity==%d", dMult));

	// Calibration histogram
	hCal = new TH1D("hCal", Form("Energy Spectrum (M%d)", dMult), dBin, dEMin, dEMax);
	qtree->Project("hCal","Energy",base_cut && mult_cut);

	// Background histogram
	// hBkg = new TH1D("hBkg", Form("Energy Spectrum (M%d) - Background", dMult), dBin, dEMin, dEMax);
	// qtree_bkg->Project("hBkg","Energy", base_cut && mult_cut);

	// Subtract off background histogram to get an estimation of "pure" calibration events
	// Currently not scaled by livetime, runs chosen were similar in livetime
	// Doesn't really make a difference
	// hCal->Add(hBkg,-1);


	dCal = qtree->GetEntries(base_cut && mult_cut && Form("Energy > %f && Energy < %f", dCutMin, dCutMax));
	dCal2 = qtree->GetEntries(base_cut && mult_cut && Form("Energy > %f && Energy < %f", dCutMin2, dCutMax2));
	dCalFull = qtree->GetEntries(base_cut && mult_cut && Form("Energy > %f && Energy < %f", dEMin, dEMax));


	hCal->SetLineColor(1);
	hCal->GetXaxis()->SetTitle("Energy (keV)");
	hCal->SetTitleSize(1);
	hCal->GetYaxis()->SetTitle("Counts/(5 keV)");
	hCal->GetXaxis()->SetLabelSize(0.04);
	hCal->GetYaxis()->SetLabelSize(0.04);
	hCal->GetXaxis()->SetTitleSize(0.04);
	hCal->GetYaxis()->SetTitleSize(0.04);
	hCal->SetStats(0);
	hCal->Draw();


	leg->AddEntry(hCal,Form("Calibration (%d keV threshold)",dThres),"l");

///////////////////////////////////////////////////////////////////////
// Load Montecarlo Data
///////////////////////////////////////////////////////////////////////

	TCut MC_cut = Form("Ener1 > %f && Ener1 < %f", dEMin, dEMax);

	// Multiplicity
	for (int k = 1; k <= 1; k++)
	{

					
		outTree = new TChain("outTree");

		// outTreeM1->Add(Form("/Users/brian/macros/Simulations/Calib/Jon/Calib-M%d-T50.root", dMult));
		outTree->Add("/Users/brian/macros/Simulations/Calib_new/Calib-p*-r2.root");

				
		// Old
		/*
		hMC = new TH1D(Form("hC%d",*angle_itr), Form("MC Energy Spectrum (M%d)", dMult), dBin, dEMin, dEMax);
		outTreeM1->Project(Form("hC%d",*angle_itr), "Ener1", MC_cut);
		// hMC = SmearSpectra();
		// hMC->Rebin(dBinSize); // To match the bin size
		*/

		hMC = new TH1D("hMC", Form("MC Energy Spectrum (M%d)", dMult), dBin, dEMin, dEMax);

		outTree->Project("hMC", "Ener1", Form("Multiplicity == %d", dMult));
		dMCFull = outTree->GetEntries(Form("Ener1 > %f && Ener1 < %f", dEMin, dEMax));

/*
		double resolution1 = 2.1512;
		double resolution2 = 4.9617;


		hSmeared = new TH1D("hSmeared", Form("MC Energy Spectrum (M%d)", dMult), dBin, dEMin, dEMax);


		double dNorm;
	    double dNorm2;
		double dSmearedValue;
		TF1 *gaus1 = new TF1("gaus","gaus(0)", 0, 3500);
		TF1 *gaus2 = new TF1("gaus","gaus(0)", 0, 3500);
		for(int i = 0; i<dBin; i++)
		{
			for(int j = 0; j<dBin; j++)
			{
			// Normalization of gaussian = (bsin size * Area of bin j in MC) / Sigma of bin j (fit function evaluated at bin center)
			dNorm = (57.584/80.296)*dBinSize*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*resolution1);
    		dNorm2 = (22.712/80.296)*dBinSize*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*resolution2);

			// Set parameters of gaussian ... 2nd gaussian *slightly* shifted... not sure if this works
			gaus->SetParameters(dNorm, hMC->GetBinCenter(j), resolution1);
     		gaus2->SetParameters(dNorm2, hMC->GetBinCenter(j)-2, resolution2);

			// Smeared contribution from gaussian centered at bin j for bin i 
			dSmearedValue = gaus->Eval(hSmeared->GetBinCenter(i)) + gaus2->Eval(hSmeared->GetBinCenter(i));

			// Fill bin i with contribution from gaussian centered at bin j
			hSmeared->Fill(hSmeared->GetBinCenter(i), dSmearedValue);
			}
		}
*/
	
		// Normalization calculation
		// dMC = outTree->GetEntries(Form("Ener1 > %f && Ener1 < %f", dCutMin, dCutMax)); // 2600 to 2630 keV
		// dMC = hMC->Integral(2600/dBinSize, 2630/dBinSize);


		dMC = hMC->Integral(2600/dBinSize, 2630/dBinSize);
		Normalization = dCal/dMC; // 2600 to 2630 keV

		// Normalization = 184.394;

		// Normalization = 148.555; // r2.5
		Normalization = 129.293; // r2.0

		cout << "Normalization (2600 to 2630 keV): " << Normalization << endl;

		cout << "Full Integral Cal: " << dCalFull << endl;
		cout << "Full Integral MC: " << dMCFull << endl;
		dRatio = dCalFull/dMCFull/Normalization;
		cout << "Ratio: " << dRatio << endl;

		hCal->Scale(1/Normalization);
		// hMC->Scale(Normalization);
		hMC->GetXaxis()->SetTitle("Energy (keV)");
		hMC->GetYaxis()->SetTitle("Counts/(5 keV)");
		hMC->GetXaxis()->SetLabelSize(0.04);
		hMC->GetYaxis()->SetLabelSize(0.04);
		hMC->SetStats(0);
		hMC->SetLineColor(kRed);	
		hMC->Draw("SAME");

		leg->AddEntry(hMC, "MC", "l");

	}

//	leg->SetTextSize(0.02);
	leg->Draw();






///////////////////////////////////////////////////////////////////////
// Residuals plot
///////////////////////////////////////////////////////////////////////

	// Select 2nd pad for residuals
	p1->cd();

	// Clone histograms for rebinning
	TH1D *hCloneCal 	= (TH1D*)hCal->Clone("hCloneCal");
	TH1D *hCloneMC		= (TH1D*)hMC->Clone("hCloneMC");

	hCloneCal->Rebin(dRebin);
	hCloneMC->Rebin(dRebin);

	// Variables used in Residual calculations
	double dResidualX, dResidualY, dResidualXErr = 0, dResidualYErr;

	// Residual plot and distribution
	TGraph *gResidual = new TGraph();
	TH1D *h1 = new TH1D("h1",Form("Residual Distribution (M%d) (%d keV threshold)", dMult, dThres), 39, -8, 11);

	for (int j = dEMinRes/(dBinSize*dRebin); j<dEMaxRes/(dBinSize*dRebin); j++)
	{
		dResidualX 		= hCloneCal->GetBinCenter(j);
		dResidualY 		= (hCloneCal->GetBinContent(j) - hCloneMC->GetBinContent(j))/
							TMath::Sqrt(hCloneCal->GetBinContent(j) + hCloneMC->GetBinContent(j)); // Sqrt(MC + data) = sigma for poisson distribution

		// dResidualY 		= (1/dBin)*TMath::Power((hCloneCal->GetBinContent(j) - hCloneMC->GetBinContent(j)),2)/
							// TMath::Sqrt(hCloneCal->GetBinContent(j)*hCloneCal->GetBinContent(j) + hCloneMC->GetBinContent(j)*hCloneMC->GetBinContent(j));

		gResidual->SetPoint(j, dResidualX, dResidualY);
		h1->Fill(dResidualY);
	}

	TAxis *axis = gResidual->GetXaxis();
	axis->SetLimits(dEMin, dEMax);
	gResidual->SetName("Residuals");
	gResidual->GetXaxis()->SetTitle("Energy (keV)");
	gResidual->GetXaxis()->SetTitleSize(0.04);
	gResidual->GetXaxis()->SetLabelSize(0.05);
	gResidual->GetYaxis()->SetLabelSize(0.05);
	gResidual->GetYaxis()->SetTitleSize(0.04);
	gResidual->GetYaxis()->SetTitle("Residuals (#sigma)");
	gResidual->SetMarkerStyle(2);
	gResidual->SetMarkerSize(2);
	gResidual->Draw("AP");

	TLine *line = new TLine();
	line->DrawLine(300,3,3000,3);
	line->DrawLine(300,-3,3000,-3);


/*
	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
//	TStyle *style2 = (TStyle)gStyle->Clone("new style");
//	style2->SetOptStat(1111);
//	h1->SetStats(1);
	h1->Fit("gaus");
	h1->Draw();



///////////////////////////////////////////////////////////////////////
// Fitting
///////////////////////////////////////////////////////////////////////



	TCanvas *c3 = new TCanvas("c3","c3", 1600, 1200);
	c3->Divide(1,2);

	// Clone histograms for fitting (I don't want them drawn on the main plot)
	TH1D *hFitCloneCal	 	= (TH1D*)hCal->Clone("hFitCloneCal");
	TH1D *hFitCloneMC		= (TH1D*)hMC->Clone("hFitCloneMC");

	// Set ranges for fits
	double dFitMinTl 	=	2595;
	double dFitMaxTl 	=	2630;
	double dFitMinSE 	=	2070;
	double dFitMaxSE 	=	2130;
	double dFitMinDE 	=	1560;
	double dFitMaxDE 	=	1610;
	double dFitMinAc1 	=	940;
	double dFitMaxAc1 	=	995;
	double dFitMinAc2 	=	875;
	double dFitMaxAc2 	=	940;
	double dFitMaxE 	=	530;
	double dFitMinE 	=	495;


	TF1 *fFitTl = new TF1("fFitTl","gaus(0) + pol3(3)", dFitMinTl, dFitMaxTl);
	TF1 *fFitSE = new TF1("fFitSE","gaus(0) + pol3(3)", dFitMinSE, dFitMaxSE);
	TF1 *fFitDE = new TF1("fFitDE","gaus(0) + pol3(3)", dFitMinDE, dFitMaxDE);
	TF1 *fFitAc1 = new TF1("fFitAc1","gaus(0) + pol3(3)", dFitMinAc1, dFitMaxAc1);
	TF1 *fFitAc2 = new TF1("fFitAc2","gaus(0) + pol3(3)", dFitMinAc2, dFitMaxAc2);
	TF1 *fFitE = new TF1("fFitE", "gaus(0) + pol3(3)", dFitMinE, dFitMaxE);


	// Initialize Parameters for the fit functions
	fFitTl->SetParameters(5000,2615,3,0,0,0,0);
	fFitTl->SetLineColor(1);
	fFitSE->SetParameters(5000,2104,5,0,0,0,0);
	fFitSE->SetLineColor(2);
	fFitDE->SetParameters(5000,1593,3,0,0,0,0);
	fFitDE->SetLineColor(3);
	fFitAc1->SetParameters(5000,968,3,0,0,0,0);
	fFitAc1->SetLineColor(4);
	fFitAc2->SetParameters(500,911,3,0,0,0,0);
	fFitAc2->SetLineColor(5);
	fFitE->SetParameters(500,511,3,0,0,0,0);
	fFitE->SetLineColor(6);

	// Fit the calibration run
	c3->cd(1);
	hFitCloneCal->Fit("fFitTl","R");
	hFitCloneCal->Fit("fFitSE","R+");
	hFitCloneCal->Fit("fFitDE","R+");
	hFitCloneCal->Fit("fFitAc1","R+");
	hFitCloneCal->Fit("fFitAc2","R+");
	hFitCloneCal->Fit("fFitE","R+");

	// Get fit parameters for Calib
	double dCalTlParam[7], dCalSEParam[7], dCalDEParam[7], dCalAc1Param[7], dCalAc2Param[7], dCalEParam[7];
	double dCalTlErr[7], dCalSEErr[7], dCalDEErr[7], dCalAc1Err[7], dCalAc2Err[7], dCalEErr[7];
	for (int j=0; j<7; j++)
	{
		dCalTlParam[j] 	=	fFitTl->GetParameter(j);
		dCalTlErr[j]	=	fFitTl->GetParError(j);
		dCalSEParam[j]	=	fFitSE->GetParameter(j);
		dCalSEErr[j]	=	fFitSE->GetParError(j);
		dCalDEParam[j]	=	fFitDE->GetParameter(j);
		dCalDEErr[j]	=	fFitDE->GetParError(j);
		dCalAc1Param[j]	=	fFitAc1->GetParameter(j);
		dCalAc1Err[j]	=	fFitAc1->GetParError(j);
		dCalAc2Param[j]	=	fFitAc2->GetParameter(j);
		dCalAc2Err[j]	=	fFitAc2->GetParError(j);
		dCalEParam[j]	=	fFitE->GetParameter(j);
		dCalEErr[j]		=	fFitE->GetParError(j);
	}


	// Fit the MC run
	c3->cd(2);
	hFitCloneMC->SetLineColor(1); // Make the cloned histogram black
	hFitCloneMC->Fit("fFitTl","R");
	hFitCloneMC->Fit("fFitSE","R+");
	hFitCloneMC->Fit("fFitDE","R+");
	hFitCloneMC->Fit("fFitAc1","R+");
	hFitCloneMC->Fit("fFitAc2","R+");
	hFitCloneMC->Fit("fFitE","R+");

	// Get fit parameters for MC
	double dMCTlParam[7], dMCSEParam[7], dMCDEParam[7], dMCAc1Param[7], dMCAc2Param[7], dMCEParam[7];
	double dMCTlErr[7], dMCSEErr[7], dMCDEErr[7], dMCAc1Err[7], dMCAc2Err[7], dMCEErr[7];
	for (int j=0; j<7; j++)
	{
		dMCTlParam[j] 	=	fFitTl->GetParameter(j);
		dMCTlErr[j]		=	fFitTl->GetParError(j);
		dMCSEParam[j]	=	fFitSE->GetParameter(j);
		dMCSEErr[j]		=	fFitSE->GetParError(j);
		dMCDEParam[j]	=	fFitDE->GetParameter(j);
		dMCDEErr[j]		=	fFitDE->GetParError(j);
		dMCAc1Param[j]	=	fFitAc1->GetParameter(j);
		dMCAc1Err[j]	=	fFitAc1->GetParError(j);
		dMCAc2Param[j]	=	fFitAc2->GetParameter(j);
		dMCAc2Err[j]	=	fFitAc2->GetParError(j);
		dMCEParam[j]	=	fFitE->GetParameter(j);
		dMCEErr[j]		=	fFitE->GetParError(j);
	}

	// Note for the Gaussian Norm: [0] = Area/([3]sqrt(2*pi))
	// Area = [0]*[3]*sqrt(2*pi)/5 <- 5 is for bin size


	// Saving all the output
	if(bSavePlot)
	{
		c1->SaveAs(Form("SpectraComparison-M%d-T%d.pdf", dMult, dThres));
		c1->SaveAs(Form("SpectraComparison-M%d-T%d.png", dMult, dThres));

		c2->SaveAs(Form("SpectraResiduals-M%d-T%d.pdf",	 dMult, dThres));
		c2->SaveAs(Form("SpectraResiduals-M%d-T%d.png",	 dMult, dThres));

		c3->SaveAs(Form("SpectraFitting-M%d-T%d.pdf",	 dMult, dThres));
		c3->SaveAs(Form("SpectraFitting-M%d-T%d.png",	 dMult, dThres));

		// Save fit parameters to a file
		ofstream FitOutput(Form("FitParams-M%d-T%d.txt", dMult, dThres));

		FitOutput << "Parameters from Tl, SE, DE, Ac1(969), Ac2(911), and Electron (511)" << endl;
		FitOutput << "1st block [0] = Cal Norm - Err - MC Norm - Err" << endl;
		FitOutput << "2nd block [1] = Cal Mean - Err - MC Mean - Err" << endl;
		FitOutput << "3rd block [2] = Cal Sigma - Err - MC Sigma - Err" << endl;
		FitOutput << "4th block [3] = Cal p0 - Err - MC p0 - Err" << endl;
		FitOutput << "5th block [4] = Cal p1 - Err - MC p1 - Err" << endl;
		FitOutput << "6th block [6] = Cal p2 - Err - MC p2 - Err" << endl;
		FitOutput << "7th block [7] = Cal p3 - Err - MC p3 - Err" << endl;
		FitOutput << "8th block [8] = Cal p3 - Err - MC p3 - Err \n" << endl;

		for (int j = 0; j < 7; j++)
		{
			FitOutput << dCalTlParam[j] << "\t" << dCalTlErr[j] << "\t" << dMCTlParam[j] << "\t" << dMCTlErr[j] << endl;
			FitOutput << dCalSEParam[j] << "\t" << dCalSEErr[j] << "\t" << dMCSEParam[j] << "\t" << dMCSEErr[j] << endl;
			FitOutput << dCalDEParam[j] << "\t" << dCalDEErr[j] << "\t" << dMCDEParam[j] << "\t" << dMCDEErr[j] << endl;
			FitOutput << dCalAc1Param[j] << "\t" << dCalAc1Err[j] << "\t" << dMCAc1Param[j] << "\t" << dMCAc1Err[j] << endl;
			FitOutput << dCalAc2Param[j] << "\t" << dCalAc2Err[j] << "\t" << dMCAc2Param[j] << "\t" << dMCAc2Err[j] << endl;
			FitOutput << dCalEParam[j] << "\t" << dCalEErr[j] << "\t" << dMCEParam[j] << "\t" << dMCEErr[j] << endl;
			FitOutput << "\n \n" << endl;

		}
		FitOutput.close();

		ofstream FitResolutions(Form("FitResolutions-M%d-T%d.txt", dMult, dThres));

		FitResolutions << "Tl208 (2615 keV): " << dCalTlParam[2] << "\t" << dCalTlErr[2] << endl;
		FitResolutions << "Single Escape (2104 keV): "<< dCalSEParam[2] << "\t" << dCalSEErr[2] << endl;
		FitResolutions << "Double Escape (1593 keV): "<< dCalDEParam[2] << "\t" << dCalDEErr[2] << endl;
		FitResolutions << "Ac228 (968 keV): "<< dCalAc1Param[2] << "\t" << dCalAc1Err[2] << endl;
		FitResolutions << "Ac228 (911 keV): "<< dCalAc2Param[2] << "\t" << dCalAc2Err[2] << endl;
		FitResolutions << "Electron (511 keV): "<< dCalEParam[2] << "\t" << dCalEErr[2] << endl;
		FitResolutions << "\n \n" << endl;

		FitResolutions.close();

	}
*/


}
