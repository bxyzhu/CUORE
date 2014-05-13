/*
	This macro compares spectra of Calibration runs with the thresholds
*/
{
	gStyle->SetOptStat(0);

	TH1D *hCal;
	TH1D *hCal2;
	TLegend *leg;
	leg = new TLegend(0.75,0.55,0.9,0.9);
	TCanvas *c1 = new TCanvas("c1","c1",1600,800);
	c1->SetLogy();

	// Histogram min/max and threshold
	double dEMin = 0;
	double dEMax = 3500;
	int bin = (dEMax - dEMin)/5;
	int dThres = 300;

	TChain *qtree = new TChain("qtree", "qtree");
	qtree->Add("/Users/brian/data/CUOREZ/BlindedReduced_200730_C_p001_T300.root");
	
	TChain *qtree2 = new TChain("qtree","qtree");
	qtree2->Add("/Users/brian/data/CUOREZ/Blinded_200730_C_p*.root");

	TCut base_cut("IsSignal");
	base_cut = base_cut && "Filter_RejectBadIntervals";
	base_cut = base_cut && "NumberOfPulses==1";
	base_cut = base_cut && "TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0";
	base_cut = base_cut && "TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0";
	base_cut = base_cut && "abs(BaselineSlope)<0.1";
	base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";
	base_cut = base_cut && "Multiplicity_OFTime==2";
//	base_cut = base_cut && Form("Energy > %d && (TotalEnergy_OFTime - Energy) > %d", dThres, dThres);
	base_cut = base_cut && Form("Energy>%f && Energy<%f", dEMin, dEMax);

	hCal = new TH1D("hCal", "Energy Spectrum (M2)",bin, dEMin, dEMax);
	hCal2 = new TH1D("hCal2", "", bin, dEMin, dEMax);


	qtree->Project("hCal","Energy",base_cut);
	qtree2->Project("hCal2", "Energy", base_cut);

	hCal->SetLineColor(1);
	hCal->GetXaxis()->SetTitle("Energy (keV)");
	hCal->GetYaxis()->SetTitle("Counts/(5 keV)");
	hCal->Draw();

	hCal2->SetLineColor(2);
	hCal2->Draw("SAME");

	leg->AddEntry(hCal,"Calib 300 keV threshold","l");
	leg->AddEntry(hCal2,"Calib no threshold","l");
	leg->Draw();

//	c1->SaveAs("SpectraComparison_M2.pdf");
}