// This macro compares ratios of M2 and M1 (as well as outputs the M2 and M1 quantities) of both calibration and MC

void CalibRatio()
{
	gStyle->SetOptStat(0);
	// Get rid of title
//	gStyle->SetOptTitle(0);

	TH1D *hCal;
	TH1D *hCal2;
	TLegend *leg;
	leg = new TLegend(0.72,0.65,0.9,0.9);
	TCanvas *c1 = new TCanvas("c1","c1",1600,800);
	c1->SetLogy();

	// Histogram min/max and threshold
//	double dEMin = 2000;
//	double dEMax = 3000;
	double dEMin = 0;
	double dEMax = 3500;
//	double dEMax = 2000;
	int bin = (dEMax - dEMin)/5;
	
	// Cuts and stuff
	double dCalM1;
	double dCalM2;
	double dCalFullM1;
	double dCalFullM2;
	double dCutMin = 2600;
	double dCutMax = 2630;
	double dRatio;
	double dRatioTl;

///////////// Calibration Data

	TChain *qtree = new TChain("qtree", "qtree");
	qtree->Add("/Users/brian/data/CUOREZ/BlindedReduced_200730_C_p001_T300.root");

	TCut base_cut("IsSignal");
	base_cut = base_cut && "Filter_RejectBadIntervals";
//	base_cut = base_cut && "NumberOfPulses==1";		
//	base_cut = base_cut && "TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0";
//	base_cut = base_cut && "TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0";
//	base_cut = base_cut && "abs(BaselineSlope)<0.1";
//	base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";

	TCut mult_cut("Multiplicity_OFTime==1");

	hCal = new TH1D("hCal", "",bin, dEMin, dEMax);
	qtree->Project("hCal","Energy",base_cut && mult_cut);

	hCal2 = new TH1D("hCal2", "", bin, dEMin, dEMax);
	qtree->Project("hCal2","Energy",base_cut && "Multiplicity_OFTime==2");

	dCalM1 = qtree->GetEntries(base_cut && mult_cut && Form("Energy > %f && Energy < %f", dCutMin, dCutMax));
	dCalFullM1 = qtree->GetEntries(base_cut && mult_cut && Form("Energy > %f && Energy < %f", dEMin, dEMax));

	dCalM2 = qtree->GetEntries(base_cut && "Multiplicity_OFTime==2" && Form("Energy > %f && Energy < %f", dCutMin, dCutMax));
	dCalFullM2 = qtree->GetEntries(base_cut && "Multiplicity_OFTime==2" && Form("Energy > %f && Energy < %f", dEMin, dEMax));

	dRatio = dCalFullM2/dCalFullM1;
	dRatioTl = dCalM2/dCalM1;
	cout << "M1: " << dCalFullM1 << " " << "M2: " << dCalFullM2 << endl;
	cout << "M2/M1 Ratio: " << dRatio << endl;
	cout << "M2/M1 Ratio Tl-208: " << dRatioTl << endl;

	hCal->SetLineColor(1);
	hCal->GetXaxis()->SetTitle("Energy (keV)");
	hCal->GetYaxis()->SetTitle("Counts/(5 keV)");
//	hCal->Sumw2();
	hCal->Draw();

	hCal2->SetLineColor(4);
	hCal2->SetLineStyle(2);
	hCal2->Draw("SAME");
	leg->AddEntry(hCal,Form("M1 - M2/M1 Ratio: %f",dRatio),"l");
	leg->AddEntry(hCal2,Form("M2 - M2/M1 Ratio(Tl): %f", dRatioTl), "l");


	leg->Draw();



}



void MCRatio(bool bDraw = false)
{
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	ofstream Outfile("MCRatio-t4.txt");

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

//	float coinList[] = {0.10000, 0.01000, 0.00100, 0.00010, 0.00001};
	float coinList[] = {0.10000};
	vector<double> coin;
	for (unsigned int i = 0; i<sizeof(coinList)/sizeof(coinList[0]); ++i)
	{
		coin.push_back(coinList[i]);
	}

	float rateList[] = {2.7500, 3.0000, 3.2537, 3.5000, 3.7500};
	vector<double> rate;
	for (unsigned int i = 0; i<sizeof(rateList)/sizeof(rateList[0]); ++i)
	{
		rate.push_back(rateList[i]);
	}

	TChain *outTreeM1;
	TChain *outTreeM2;
	TH1D *hMC;
	TH1D *hMC2;

	if (bDraw)
	{
		leg = new TLegend(0.72,0.65,0.9,0.9);
		TCanvas *c1 = new TCanvas("c1","c1",1600,800);
		c1->SetLogy();
	}
///// Histogram min/max and threshold /////
//	double dEMin = 2000;
//	double dEMax = 3000;
	double dEMin = 0;
	double dEMax = 3500;
//	double dEMax = 2000;
	int bin = (dEMax - dEMin)/5;
	int dThres = 300;


////// Cuts and stuff /////
	double dMCM1;
	double dMCM2;
	double dMCM1Full;
	double dMCM2Full;
	double dCutMin = 2600;
	double dCutMax = 2630;
	double dRatio;
	double dRatioTl;

	TCut MC_cut = Form("Heat > %f && Heat < %f", dEMin, dEMax);


	int j = 0;
	for (vector<int>::const_iterator angle_itr = angle.begin(); angle_itr != angle.end(); angle_itr++)
	{
		for (vector<int>::const_iterator thresh_itr = thresh.begin(); thresh_itr != thresh.end(); thresh_itr++)
		{
			for (vector<double>::const_iterator coin_itr = coin.begin(); coin_itr != coin.end(); coin_itr++)
			{
				for (vector<double>::const_iterator rate_itr = rate.begin(); rate_itr != rate.end(); rate_itr++)
				{
				outTreeM1 = new TChain("outTree");
				outTreeM1->Add(Form("/Users/brian/macros/Simulations/Calib2/Left-209cm-C%d-M1-T%d-r%.4f-c%.5f.root", *angle_itr, *thresh_itr, *rate_itr, *coin_itr));
				outTreeM1->Add(Form("/Users/brian/macros/Simulations/Calib2/Right-207cm-C%d-M1-T%d-r%.4f-c%.5f.root", *angle_itr, *thresh_itr, *rate_itr, *coin_itr));

				outTreeM2 = new TChain("outTree");
				outTreeM2->Add(Form("/Users/brian/macros/Simulations/Calib2/Left-209cm-C%d-M2-T%d-r%.4f-c%.5f.root", *angle_itr, *thresh_itr, *rate_itr, *coin_itr));
				outTreeM2->Add(Form("/Users/brian/macros/Simulations/Calib2/Right-207cm-C%d-M2-T%d-r%.4f-c%.5f.root", *angle_itr, *thresh_itr, *rate_itr, *coin_itr));	

				hMC = new TH1D(Form("hC%d-T%d-r%.4f-c%.5f-M1", *angle_itr, *thresh_itr, *rate_itr, *coin_itr), "", bin, dEMin, dEMax);
				outTreeM1->Project(Form("hC%d-T%d-r%.4f-c%.5f-M1", *angle_itr, *thresh_itr, *rate_itr, *coin_itr), "Heat", MC_cut);
				hMC2 = new TH1D(Form("hC%d-T%d-r%.4f-c%.5f-M2", *angle_itr, *thresh_itr, *rate_itr, *coin_itr), "", bin, dEMin, dEMax);
				outTreeM2->Project(Form("hC%d-T%d-r%.4f-c%.5f-M2", *angle_itr, *thresh_itr, *rate_itr, *coin_itr), "Heat", MC_cut);


				dMCM1 = outTreeM1->GetEntries(Form("Heat > %f && Heat < %f", dCutMin, dCutMax));
				dMCM1Full = outTreeM1->GetEntries(Form("Heat > %f && Heat < %f", dEMin, dEMax));

				dMCM2 = outTreeM2->GetEntries(Form("Heat > %f && Heat < %f", dCutMin, dCutMax));
				dMCM2Full = outTreeM2->GetEntries(Form("Heat > %f && Heat < %f", dEMin, dEMax));

				dRatio = dMCM2Full/dMCM1Full;
				dRatioTl = dMCM2/dMCM1;

				Outfile << "Threshold: " << *thresh_itr << " " << "Coincidence: " << *coin_itr << " " << "Rate: " << *rate_itr << endl;
				Outfile << "M1: " << dMCM1Full << " " << "M2: " << dMCM2Full << endl;
				Outfile << "Ratio M2/M1: " << dRatio << endl;
				Outfile << "RatioTl M2/M1: " << dRatioTl << endl;
				Outfile << " " << endl;

				if(bDraw)
				{
					hMC->SetLineColor(j+1);	
					hMC->Draw("SAME");
					hMC2->SetLineColor(j+2);	
					hMC2->Draw("SAME");
					leg->AddEntry(hMC, Form("hC%d-M1 - M2/M1 Ratio: %f", *angle_itr, dRatio), "l");
					leg->AddEntry(hMC2, Form("hC%d-M2 - M2/M1 Ratio(Tl): %f", *angle_itr, dRatioTl), "l");		
				}
				j++;
				}
			}
		}
	}
	Outfile.close();

	if(bDraw)
	{
		leg->Draw();
//		c1->SaveAs(Form("RatioMC.pdf", dThres));

	}

}