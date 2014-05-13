// This macro compares the M2Tot spectra and other asc2tas txt files (prefer to use the ROOT files)

void DrawSpec2(const int dMult=2, bool fSavePlot = false)
{

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


	TH1D *hMC;
	TH1D *hMC2;
	TLegend *leg;
	leg = new TLegend(0.75,0.55,0.9,0.9);


	// Histogram min/max and threshold
	int dBinSize = 5;
	double dEMin = 0;
	double dEMax = 3500;
	int dBin = (dEMax - dEMin)/dBinSize;

	// Bins 
	int dRebin = 2; // Rebin value (for residuals)
	int dEMinRes = 300;
	int dEMaxRes = 2700;

	// Cuts and stuff
	double dMC, dMC2;
	double dMCFull, dMCFull2;
	double dCutMin 		= 	2600; // first normalization
	double dCutMax 		= 	2630;
	double Normalization;

	int dThres = 100;

	qtree = new TChain("qtree", "qtree");
	qtree->Add(Form("/Users/brian/data/CUOREZ/BlindedReduced_200730_C_p001_T%d.root", dThres));
	qtree->Add(Form("/Users/brian/data/CUOREZ/BlindedReduced_200757_C_p001_T%d.root", dThres));


	// qtree_bkg = new TChain("qtree","qtree");
	// qtree_bkg->Add("/Users/brian/data/CUOREZ/BlindedReduced_200713_B_p001.root");

	// TCut base_cut("IsSignal");
	TCut base_cut;
	// base_cut = base_cut && "Filter_RejectBadIntervals";
	base_cut = base_cut && "NumberOfPulses==1";
	base_cut = base_cut && "TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0";
	base_cut = base_cut && "TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0";
	base_cut = base_cut && "abs(BaselineSlope)<0.1";
	base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";
	base_cut = base_cut && "OrderInMultiple_OFTime == 1";

	TCut mult_cut(Form("Multiplicity_OFTime==%d", dMult));

	// Calibration histogram
	hCal = new TH1D("hCal", Form("Total Energy Spectrum (M%dTot)", dMult), dBin, dEMin, dEMax);
	qtree->Project("hCal","TotalEnergy_OFTime",base_cut && mult_cut);
	// qtree->Project("hCal","Energy", base_cut && mult_cut);

	// Background histogram
	// hBkg = new TH1D("hBkg", Form("Energy Spectrum (M%d) - Background", dMult), dBin, dEMin, dEMax);
	// qtree_bkg->Project("hBkg","Energy", base_cut && mult_cut);

	// Subtract off background histogram to get an estimation of "pure" calibration events
	// Currently not scaled by livetime, runs chosen were similar in livetime
	// Doesn't really make a difference
	// hCal->Add(hBkg,-1);


///////////// Montecarlo Data

	hMC = new TH1D("hMC", "hMC", dBin, dEMin, dEMax);


	ifstream input;
	input.open(Form("/Users/brian/macros/Simulations/Calib/Calib-C0-T%d-M2A-sum.txt", dThres));
	int k = 0;
	double number;
	while(input)
	{
		input >> number; 
		hMC->Fill(k,number);
		k++;
	}


	// Select 1st pad
	p2->cd();
	p2->SetLogy();

	// dCal = qtree->GetEntries(base_cut && mult_cut && Form("Energy > %f && Energy < %f", dCutMin, dCutMax));
	// dMC = hMC->Integral();
	// double Normalization = 0.314588;
	double Normalization = 0.313728;
	hMC->Scale(Normalization);


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

	hMC->SetLineColor(2);	
	hMC->Draw("SAME");		

	TLegend *leg;
	leg = new TLegend(0.72,0.70,0.925,0.9);

	leg->AddEntry(hCal,Form("Calibration (%d keV threshold)",dThres),"l");
	leg->AddEntry(hMC, "MC", "l");
	leg->Draw();


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
	TH1D *h1 = new TH1D("h1",Form("Residual Distribution (M%dTot) (%d keV threshold)", dMult, dThres), 39, -8, 11);

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

	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
//	TStyle *style2 = (TStyle)gStyle->Clone("new style");
//	style2->SetOptStat(1111);
//	h1->SetStats(1);
	h1->Fit("gaus");
	h1->Draw();

	if(fSavePlot)
	{
		c1->SaveAs(Form("SpectraComparison-M2Tot-T%d.png", dThres));
		c1->SaveAs(Form("SpectraComparison-M2Tot-T%d.pdf", dThres));
		c2->SaveAs(Form("SpectraResiduals-M2Tot-T%d.png", dThres));
		c2->SaveAs(Form("SpectraResiduals-M2Tot-T%d.pdf", dThres));
	}

}
