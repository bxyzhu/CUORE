/*
	This macro compares the Calibration data with the MC spectra
	This macro is for g2tas output files (with outTree)
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

// Choose multiplicity, saving plots, and energy range for plots
void DrawSpec_Old(int dMult = 1, bool bSavePlot = false, double dEMin = 0, double dEMax = 3500)
{
//	gStyle->Reset();
	gStyle->SetOptStat("emr");
	gStyle->SetOptFit();
	gStyle->SetStatW(0.15);
	gStyle->SetStatH(0.2);
//	gStyle->SetOptTitle(0); 	// Gets rid of titles


	TCanvas *c1 = new TCanvas("c1","c1",1200,700);

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
	float coinList[] = {0.10};
	vector<double> coin;
	for (unsigned int i = 0; i<sizeof(coinList)/sizeof(coinList[0]); ++i)
	{
		coin.push_back(coinList[i]);
	}

	float rateList[] = {2.3259};
	vector<double> rate;
	for (unsigned int i = 0; i<sizeof(rateList)/sizeof(rateList[0]); ++i)
	{
		rate.push_back(rateList[i]);
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

	c1->SetLogy();

	qtree = new TChain("qtree", "qtree");
	// qtree->Add("/Users/brian/data/CUOREZ/noThresh/BlindedReduced_200730_C_p001.root");
	qtree->Add("/Users/brian/data/CUOREZ/BlindedReduced_200730_C_p001.root");

//	qtree->Add("/Users/brian/data/CUOREZ/BlindedReduced_200757_C_p001.root");


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
				outTreeM1 = new TChain("outTree");
				outTreeM1->Add(Form("/Users/brian/macros/Simulations/Calib2/rootfile_rate2/Left-209cm-C%d-M%d-T%d-r%.4f-c%.2f.root", 
										*angle_itr, dMult, *thresh_itr, *rate_itr, *coin_itr));
				outTreeM1->Add(Form("/Users/brian/macros/Simulations/Calib2/rootfile_rate2/Right-207cm-C%d-M%d-T%d-r%.4f-c%.2f.root", 
										*angle_itr, dMult, *thresh_itr, *rate_itr, *coin_itr));


				hMC = new TH1D(Form("hC%d",*angle_itr), Form("MC Energy Spectrum (M%d)", dMult), dBin, dEMin, dEMax);
				outTreeM1->Project(Form("hC%d",*angle_itr), "Heat", MC_cut);
				dMCFull = outTreeM1->GetEntries(Form("Heat > %f && Heat < %f", dEMin, dEMax));
				
				// Normalization calculation
				dMC = outTreeM1->GetEntries(Form("Heat > %f && Heat < %f", dCutMin, dCutMax)); // 2600 to 2630 keV
				// dMC2 = outTreeM1->GetEntries(Form("Heat > %f && Heat < %f", dCutMin2, dCutMax2)); // 1000 to 2000 keV
				// Normalization = dCal/dMC; // 2600 to 2630 keV
				// Normalization2 = dCal2/dMC2; // 1000 to 2000 keV



				// Old normalization calculations
				// Normalization = 1.44326; // without threshold
//				Normalization = 1.33944; // With PSA
				Normalization = 1.44676; // w/out PSA w/ threshold
				cout << "Normalization (2600 to 2630 keV): " << Normalization << endl;
				// cout << "Normalization (1000 to 2000 keV): " << Normalization2 << endl;				

				cout << "Full Integral Cal: " << dCalFull << endl;
				cout << "Full Integral MC: " << dMCFull << endl;
				dRatio = dCalFull/dMCFull/Normalization;
				// dRatio2 = dCalFull/dMCFull/Normalization2;

				hMC->Scale(Normalization);
				hMC->GetXaxis()->SetTitle("Energy (keV)");
				hMC->GetYaxis()->SetTitle("Counts/(5 keV)");
				hMC->SetStats(0);
				hMC->SetLineColor(j);	
				hMC->Draw("SAME");

				// leg->AddEntry(hMC, Form("hMC Calib/MC Ratio = %f", dRatio), "l");
				leg->AddEntry(hMC, "hMC", "l");

				++j;
				}
			}
		}
	}
//	leg->SetTextSize(0.02);
	leg->Draw();



	// Saving all the output
	if(bSavePlot)
	{
		// c1->SaveAs(Form("SpectraComparison-M%d.pdf", dMult));
		c1->SaveAs(Form("SpectraComparison-M%d_oldT300.png", dMult));
	}

}
