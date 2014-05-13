/*
	This macro compares the Calibration data with the multiplicity of manually calculated spectra
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

/*
TCut SetBaseCuts()
{
	TCut base_cut("IsSignal");
	base_cut = base_cut && "Filter_RejectBadIntervals";
	base_cut = base_cut && "NumberOfPulses==1";
	base_cut = base_cut && "TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0";
	base_cut = base_cut && "TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0";
	base_cut = base_cut && "abs(BaselineSlope)<0.1";
	base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";

	return base_cuts;
}

TChain *LoadCalib()
{
	TChain *qtree_calib = new TChain("qtree","qtree");
	qtree_calib->Add("");
	return qtree_calib;
}


TChain *LoadBkg()
{
	TChain *qtree_bkg = new TChain("qtree","qtree");
	qtree_bkg->Add("");
	return qtree_bkg;
}
*/

// Choose multiplicity, saving plots, and energy range for plots
void DrawSpecManual(int dMult = 1, bool bSavePlot = false, double dEMin = 0, double dEMax = 3500)
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
	double canBotMargin =0.025;
	double canTopMargin =0.025;
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


	// Histogram min/max and bin size
	int dBinSize = 5;
	int dBin = (dEMax - dEMin)/dBinSize;

	// Set up range for residuals plot
	int dRebin = 2; // Rebin value (for residuals)
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
	double dCutMin 		= 	2600; // first normalization
	double dCutMax 		= 	2630;
	double dCutMin2 	= 	1000; // second normalization
	double dCutMax2 	= 	2000;
	double Normalization;
	double Normalization2;
	double dRatio, dRatio2;

	int angleList[] = {0};
	vector<int> angle;
	for (unsigned int i = 0; i<sizeof(angleList)/sizeof(angleList[0]); ++i)
	{
		angle.push_back(angleList[i]);
	}

	int threshList[] = {300};
	vector<int> thresh;
	for (unsigned int i = 0; i<sizeof(threshList)/sizeof(threshList[0]); ++i)
	{
		thresh.push_back(threshList[i]);
	}
	int dThres = 300;

//	float coinList[] = {0.10000, 0.01000, 0.00100, 0.00010, 0.00001};
	float coinList[] = {0.10000};
	vector<double> coin;
	for (unsigned int i = 0; i<sizeof(coinList)/sizeof(coinList[0]); ++i)
	{
		coin.push_back(coinList[i]);
	}

	float rateList[] = {3.2537};
	vector<double> rate;
	for (unsigned int i = 0; i<sizeof(rateList)/sizeof(rateList[0]); ++i)
	{
		rate.push_back(rateList[i]);
	}

	float pileupList[] = {0.0};
	vector<double> pileUp;
	for (unsigned int i = 0; i<sizeof(pileupList)/sizeof(pileupList[0]); ++i)
	{
		pileUp.push_back(pileupList[i]);
	}

	TChain *qtree;
	TChain *qtree_bkg;
	TChain *outTreeM1;
	TH1D *hMC;
	TH1D *hCal;
	TH1D *hBkg;
	TLegend *leg;
	leg = new TLegend(0.74,0.55,0.925,0.9);

///////////////////////////////////////////////////////////////////////
// Load Calibration Data
///////////////////////////////////////////////////////////////////////

	// Select 1st pad
	p2->cd();
	p2->SetLogy();

	qtree = new TChain("qtree", "qtree");
	qtree->Add("/Users/brian/data/CUOREZ/BlindedReduced_200730_C_p001.root");
	qtree->Add("/Users/brian/data/CUOREZ/BlindedReduced_200757_C_p001.root");


	qtree_bkg = new TChain("qtree","qtree");
	qtree_bkg->Add("/Users/brian/data/CUOREZ/BlindedReduced_200713_B_p001.root");

	TCut base_cut("IsSignal");
	base_cut = base_cut && "Filter_RejectBadIntervals";
	base_cut = base_cut && "NumberOfPulses==1";
	base_cut = base_cut && "TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0";
	base_cut = base_cut && "TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0";
	base_cut = base_cut && "abs(BaselineSlope)<0.1";
	base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";

	TCut mult_cut(Form("Multiplicity_OFTime==%d", dMult));

	TCut mult_cut_mc(Form("Multiplicity==%d", dMult));
	// Calibration histogram
	hCal = new TH1D("hCal", Form("Energy Spectrum (M%d)", dMult), dBin, dEMin, dEMax);
	qtree->Project("hCal","Energy",base_cut && mult_cut);

	// Background histogram
	hBkg = new TH1D("hBkg", Form("Energy Spectrum (M%d) - Background", dMult), dBin, dEMin, dEMax);
	qtree_bkg->Project("hBkg","Energy", base_cut && mult_cut);

	// Subtract off background histogram to get an estimation of "pure" calibration events
	// Currently not scaled by livetime, runs chosen were similar in livetime 
	// hCal->Add(hBkg,-1);


	dCal = qtree->GetEntries(base_cut && mult_cut && Form("Energy > %f && Energy < %f", dCutMin, dCutMax));
	dCal2 = qtree->GetEntries(base_cut && mult_cut && Form("Energy > %f && Energy < %f", dCutMin2, dCutMax2));
	dCalFull = qtree->GetEntries(base_cut && mult_cut && Form("Energy > %f && Energy < %f", dEMin, dEMax));


	hCal->SetLineColor(1);
	hCal->GetXaxis()->SetTitle("Energy (keV)");
	hCal->SetTitleSize(1);
	hCal->GetYaxis()->SetTitle("Counts/(5 keV)");
	hCal->SetStats(0);
	hCal->Draw();


	leg->AddEntry(hCal,"Calibration (300 keV threshold)","l");

///////////////////////////////////////////////////////////////////////
// Load Montecarlo Data
///////////////////////////////////////////////////////////////////////

	TCut MC_cut = Form("Heat > %f && Heat < %f", dEMin, dEMax);

	// Int to set number for histogram color (when iterating over multiple values)
	int j = 2;

	for (vector<int>::const_iterator angle_itr = angle.begin(); angle_itr != angle.end(); angle_itr++)
	{
		for (vector<int>::const_iterator thresh_itr = thresh.begin(); thresh_itr != thresh.end(); thresh_itr++)
		{
			for (vector<double>::const_iterator coin_itr = coin.begin(); coin_itr != coin.end(); coin_itr++)
			{
				for (vector<double>::const_iterator rate_itr = rate.begin(); rate_itr != rate.end(); rate_itr++)
				{

				outTreeM1 = new TChain("qtree");
				// outTreeM1->Add(Form("/Users/brian/macros/Simulations/Calib2/Left-209cm-C%d-M%d-T%d-r%.4f-c%.5f.root", 
										// *angle_itr, dMult, *thresh_itr, *rate_itr, *coin_itr));
				// outTreeM1->Add(Form("/Users/brian/macros/Simulations/Calib2/Right-207cm-C%d-M%d-T%d-r%.4f-c%.5f.root", 
										// *angle_itr, dMult, *thresh_itr, *rate_itr, *coin_itr));

				outTreeM1->Add(Form("/Users/brian/macros/Simulations/Calib2/qtree-C%d-T%d-r%.4f-c%.5f-test.root",
										*angle_itr, *thresh_itr, *rate_itr, *coin_itr));


					for(vector<double>::const_iterator pileup_itr = pileUp.begin(); pileup_itr != pileUp.end(); pileup_itr++)
					{
					hMC = new TH1D(Form("hC%d",*angle_itr), Form("MC Energy Spectrum (M%d)", dMult), dBin, dEMin, dEMax);
					outTreeM1->Project(Form("hC%d",*angle_itr), "Heat", MC_cut && mult_cut_mc && Form("TimeSinceEvent >= %f", *pileup_itr));
					// outTreeM1->Project(Form("hC%d",*angle_itr), "Heat", MC_cut && mult_cut_mc); // Without pile-up

					dMCFull = outTreeM1->GetEntries(Form("Heat > %f && Heat < %f", dEMin, dEMax) && mult_cut_mc && Form("TimeSinceEvent >= %f", *pileup_itr));
					// dMCFull = outTreeM1->GetEntries(Form("Heat > %f && Heat < %f", dEMin, dEMax) && mult_cut_mc); // without pile-up
				
					// Normalization calculation
					dMC = outTreeM1->GetEntries(Form("Heat > %f && Heat < %f", dCutMin, dCutMax) && mult_cut_mc && Form("TimeSinceEvent >= %f", *pileup_itr)); // 2600 to 2630 keV
					// dMC = outTreeM1->GetEntries(Form("Heat > %f && Heat < %f", dCutMin, dCutMax) && mult_cut_mc); // 2600 to 2630 keV without pile-up
					dMC2 = outTreeM1->GetEntries(Form("Heat > %f && Heat < %f", dCutMin2, dCutMax2) && mult_cut_mc && Form("TimeSinceEvent >= %f", *pileup_itr)); // 1000 to 2000 keV
					// dMC2 = outTreeM1->GetEntries(Form("Heat > %f && Heat < %f", dCutMin2, dCutMax2) && mult_cut_mc); // 1000 to 2000 keV without pile-up

					// Normalization = dCal/dMC; // 2600 to 2630 keV
					// Normalization2 = dCal2/dMC2; // 1000 to 2000 keV

					// Normalization of higher statistic MC
					Normalization = 0.299925; // Manual multiplicity Paralyzed coincidence window
					Normalization2 = 0.318198;
					// Normalization = 0.24185; // Manual multiplicity Paralyzed coincidence window of 0.5 sec
					// Normalization2 = 0.256583;
					// Normalization = 0.309526; // Manual multiplicity + pile-up of 1 sec
					// Normalization2 = 0.329368; 
					// Normalization = 0.346524; // Manual multiplicity + pile-up of 4 secs
					// Normalization2 = 0.3694; 

					cout << "Normalization (2600 to 2630 keV): " << Normalization << endl;
					cout << "Normalization (1000 to 2000 keV): " << Normalization2 << endl;				

					cout << "Full Integral Cal: " << dCalFull << endl;
					cout << "Full Integral MC: " << dMCFull << endl;
					dRatio = dCalFull/dMCFull/Normalization;
					dRatio2 = dCalFull/dMCFull/Normalization2;

					hMC->Scale(Normalization);
					hMC->GetXaxis()->SetTitle("Energy (keV)");
					hMC->GetYaxis()->SetTitle("Counts/(5 keV)");
					hMC->SetStats(0);
					hMC->SetLineColor(j);	
					hMC->Draw("SAME");

					leg->AddEntry(hMC, Form("hMC Calib/MC Ratio = %f", dRatio), "l");

					++j;
					}	// Pile-up
				}	// Rate
			}	// Coincidence
		}	// Threshold
	}	// Angle
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
	TGraph *gResidual = new TGraphErrors();
	TH1D *h1 = new TH1D("h1",Form("Residual Distribution (M%d)", dMult), 39, -8, 11);

	for (int j = dEMinRes/(dBinSize*dRebin); j<dEMaxRes/(dBinSize*dRebin); j++)
	{
		dResidualX 		= hCloneCal->GetBinCenter(j);
		dResidualY 		= (hCloneCal->GetBinContent(j) - hCloneMC->GetBinContent(j))/
							TMath::Sqrt(hCloneCal->GetBinContent(j) + hCloneMC->GetBinContent(j));

		gResidual->SetPoint(j, dResidualX, dResidualY);
		h1->Fill(dResidualY);
	}

	TAxis *axis = gResidual->GetXaxis();
	axis->SetLimits(dEMin, dEMax);
	gResidual->SetName("Residuals");
	gResidual->GetXaxis()->SetTitle("Energy (keV)");
	gResidual->GetXaxis()->SetTitleSize(0.045);
	gResidual->GetYaxis()->SetTitle("Residuals (#sigma)");
	gResidual->SetMarkerStyle(2);
	gResidual->SetMarkerSize(2);
	gResidual->Draw("AP");

	TLine *line = new TLine();
	line->DrawLine(300,3,3000,3);
	line->DrawLine(300,-3,3000,-3);

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
		c1->SaveAs(Form("SpectraComparison-M%d-P0.pdf", dMult));
		c2->SaveAs(Form("SpectraResiduals-M%d-P0.pdf",	 dMult));
		c3->SaveAs(Form("SpectraFitting-M%d-P0.pdf",	 dMult));

		// Save fit parameters to a file
		ofstream FitOutput(Form("FitParams-M%d-P0.txt", dMult));

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
	}

}