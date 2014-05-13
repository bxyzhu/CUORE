void DrawBkg()
{
    gStyle->SetOptStat(0);
	gStyle->SetOptFit();	

    int dMult = 1;
    float dEMin = 0, dEMax = 8000;

    // Load Data files

    // const int datasets[] = {2061, 2064, 2067, 2070, 2073, 2076};
    const int datasets[] = {2061};
    const int nDataset = sizeof(datasets)/sizeof(int);

    TChain *qtree = new TChain("qtree");
    qtree->Add("/Users/brian/macros/CUOREZ/Bkg/ReducedBkg-ds*.root");

    TCut base_cut;
    base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "abs(BaselineSlope)<0.1";
    base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";

    TH1D *hfull;
    TH1D *hzoomed;
    TH1D *hzoomedm1;
    TH1D *hzoomedm2;

    TCanvas *cfull = new TCanvas("cfull", "cfull", 1200, 800);
    // gPad->SetLogy();
    hzoomed = new TH1D("hzoomed", "Zoomed Spectra", 1250, 3000, 8000);
    hzoomedm1 = new TH1D("hzoomedm1", "Zoomed Spectra", 1250, 3000, 8000);
    hzoomedm2 = new TH1D("hzoomedm2", "Zoomed Spectra", 1250, 3000, 8000);

    hzoomed->SetFillColor(kBlue);
    hzoomed->GetXaxis()->SetTitle("Energy (keV)");
    hzoomed->GetYaxis()->SetTitle("Counts/(4 keV)");

    hzoomedm1->SetFillColor(2);
    hzoomedm2->SetFillColor(5);

    qtree->Project("hzoomed", "Energy", base_cut);
    qtree->Project("hzoomedm1", "Energy", base_cut && "Multiplicity_OFTime == 1");
    qtree->Project("hzoomedm2", "Energy", base_cut && "Multiplicity_OFTime == 2");

    hzoomed->Draw();
    // hzoomedm1->Draw("SAME");
    // hzoomedm2->Draw("SAME");

    TLegend *leg;
    leg = new TLegend(0.60,0.75,0.925,0.9);

    leg->AddEntry(hzoomed,"Total","f");
    // leg->AddEntry(hzoomedm1,"Multiplicity 1","f");
    // leg->AddEntry(hzoomedm2,"Multiplicity 2","f");

    // leg->Draw();


// Drawing out spectra for each channel individually    
/*


    TCanvas *c1 = new TCanvas("c1","c1",1600,600);
    c1->Divide(2);

    for(int i = 2; i<53; i++)
    // for(int i = 2)
    {
        if(i == 1 || i == 10 || i == 49) continue;

    // int i = 2;

        hfull = new TH1D(Form("hfull-ch%d", i), Form("Full Spectra Channel %d", i), 1000, 0, 10000);
        hzoomed = new TH1D(Form("hzoomed-ch%d", i), Form("Zoomed Spectra Channel %d", i), 200, 5200, 5600);

        hfull->SetFillColor(kBlue);
        hfull->GetXaxis()->SetTitle("Energy (keV)");
        hfull->GetYaxis()->SetTitle("Counts/(10 keV)");        
        hzoomed->SetFillColor(kBlue);
        hzoomed->GetXaxis()->SetTitle("Energy (keV)");
        hzoomed->GetYaxis()->SetTitle("Counts/(2 keV)");


        qtree->Project(Form("hfull-ch%d", i), "Energy", base_cut && Form("Channel == %d", i));
        qtree->Project(Form("hzoomed-ch%d", i), "Energy", base_cut && Form("Channel == %d", i));

        c1->cd(1);
        hfull->Draw();
        c1->cd(2);
        hzoomed->Draw();

        c1->SaveAs(Form("BkgSpectra-Ch%d.png",i));
    }
*/
}


// Adds new branch with alpha energies to ROOT TTree
void AddAlphaEnergies()
{
    TChain *qtree = new TChain("qtree");
    // qtree->Add(Form("/Users/brian/macros/CUOREZ/Bkg/ReducedBkg-ds%d.root", dataset));
    qtree->Add("/Users/brian/macros/CUOREZ/Bkg/ReducedBkg-ds*.root");
    // qtree->Add("/Users/brian/macros/CUOREZ/Bkg/ReducedCorrectedBkg.root");
    // QTree variables
    int Run, Channel, Multiplicity_OFTime, OrderInMultiple;
    double Energy, OF_TVR, OF_TVL, BaselineSlope, TotalEnergy_OFTime, CorrectedEnergy, CorrectedEnergyPt;
    double TimeSinceSignalEvent_SameChannel, TimeUntilSignalEvent_SameChannel;

    // QTree branches
    qtree->SetBranchAddress("Run",                                  &Run);
    qtree->SetBranchAddress("Channel",                              &Channel);
    qtree->SetBranchAddress("Energy",                               &Energy);
    qtree->SetBranchAddress("OF_TVL",                               &OF_TVL);
    qtree->SetBranchAddress("OF_TVR",                               &OF_TVR);
    qtree->SetBranchAddress("BaselineSlope",                        &BaselineSlope);
    qtree->SetBranchAddress("OrderInMultiple",                      &OrderInMultiple);
    qtree->SetBranchAddress("Multiplicity_OFTime",                  &Multiplicity_OFTime);
    qtree->SetBranchAddress("TotalEnergy_OFTime",                   &TotalEnergy_OFTime);
    // qtree->SetBranchAddress("CorrectedEnergy",                      &CorrectedEnergy);
    qtree->SetBranchAddress("TimeUntilSignalEvent_SameChannel",     &TimeUntilSignalEvent_SameChannel);
    qtree->SetBranchAddress("TimeSinceSignalEvent_SameChannel",     &TimeSinceSignalEvent_SameChannel);



    // New Tree variables
    // TFile *f1 = new TFile(Form("ReducedCorrectedBkg-ds%d.root",dataset), "recreate");
    TFile *f1 = new TFile("ReducedCorrectedBkg.root", "recreate");
    TTree *outTree = new TTree("qtree", "Reduced Data Tree");

    outTree->Branch("Run",                  &Run,                   "Run/I");
    outTree->Branch("Channel",              &Channel,               "Channel/I");
    outTree->Branch("Energy",               &Energy,                "Energy/D");
    outTree->Branch("OF_TVL",               &OF_TVL,                "OF_TVL/D");
    outTree->Branch("OF_TVR",               &OF_TVR,                "OF_TVR/D");
    outTree->Branch("BaselineSlope",        &BaselineSlope,         "BaselineSlope/D");
    outTree->Branch("OrderInMultiple",      &OrderInMultiple,       "OrderInMultiple/I");
    outTree->Branch("Multiplicity_OFTime",  &Multiplicity_OFTime,   "Multiplicity_OFTime/I");
    outTree->Branch("TotalEnergy_OFTime",   &TotalEnergy_OFTime,    "TotalEnergy_OFTime/D");
    outTree->Branch("CorrectedEnergy",      &CorrectedEnergy,       "CorrectedEnergy/D");
    // outTree->Branch("CorrectedEnergyPt",    &CorrectedEnergyPt,     "CorrectedEnergyPt/D");
    outTree->Branch("TimeSinceSignalEvent_SameChannel",     &TimeSinceSignalEvent_SameChannel,  "TimeSinceSignalEvent_SameChannel/D");
    outTree->Branch("TimeUntilSignalEvent_SameChannel",     &TimeUntilSignalEvent_SameChannel,  "TimeUntilSignalEvent_SameChannel/D");


    int nEvent = qtree->GetEntries();
    // std::cout <<nEvent << " Events" << std::endl;

    for(int i = 0; i < nEvent; i++)
    {
        qtree->GetEntry(i);

        if(Energy > 3000)
        {
            CorrectedEnergy = Energy*1.00576 - 63.5934;
            // CorrectedEnergyPt = Energy*1.00845 - 77.6892;
        }
        else if(Energy <= 3000)
        {
            CorrectedEnergy = Energy;
        }
        outTree->Fill();
    }

    outTree->Write();

}

void EnergyResolutions(int dMult = 1, bool bSavePlots = false)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit();    

    // int dMult = 1;
    float dEMin = 0, dEMax = 8000;    
    TChain *qtree = new TChain("qtree");
    qtree->Add("/Users/brian/macros/CUOREZ/Bkg/ReducedCorrectedBkg.root");

    TCut base_cut;
    base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "abs(BaselineSlope)<0.1";
    base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";
    // base_cut = base_cut && "Run == 2"; // Test out individual channels

    TH1D *hzoomed;
    TH1D *hzoomed2;
    TH1D *hgamma;
    TH1D *halpha;

    // gPad->SetLogy();
    hzoomed = new TH1D("hzoomed", "", 8000, 0000, 8000);  
    // hzoomed = new TH1D("hzoomed", "", 4000, 0000, 8000);    
    hzoomed->SetLineColor(kBlack);
    hzoomed2 = new TH1D("hzoomed2", "", 10000, 0000, 8000);
    // hzoomed2 = new TH1D("hzoomed2", "", 4000, 0000, 8000);
    hzoomed2->SetLineColor(kBlack);

    hgamma = new TH1D("hgamma","", 3000, 0, 3000);

    // hzoomedm1 = new TH1D("hzoomedm1", "Zoomed Spectra", 1250, 3000, 8000);
    // hzoomedm2 = new TH1D("hzoomedm2", "Zoomed Spectra", 1250, 3000, 8000);


    qtree->Project("hzoomed", "Energy", base_cut && Form("Multiplicity_OFTime==%d", dMult));
    qtree->Project("hzoomed2", "Energy", base_cut && Form("Multiplicity_OFTime==%d", dMult));
    qtree->Project("hgamma", "Energy", base_cut && Form("Multiplicity_OFTime==%d", dMult));

    // Po210 - 5407 - 5304
    double dFitMaxPo    =   5520;
    double dFitMinPo    =   5280;
 
    // Th232 - 4083
    double dFitMaxTh    =   4065;
    double dFitMinTh    =   4088;
    
    // Pt190 - 3249
    double dFitMaxPt    =   3320;
    double dFitMinPt    =   3260;

    // Tl208 - 2615
    double dFitMinTl    =   2600;
    double dFitMaxTl    =   2630;
    // Single Escape - 2104
    double dFitMinSE    =   2090;
    double dFitMaxSE    =   2115;

    // Double Escape - 1593
    double dFitMinDE    =   1560;
    double dFitMaxDE    =   1610;

    // K40 - 1460
    double dFitMaxK     =   1480;
    double dFitMinK     =   1440;

    // Ac228 - 968
    double dFitMinAc1   =   960;
    double dFitMaxAc1   =   976;
    // Ac228 - 911
    double dFitMinAc2   =   875;
    double dFitMaxAc2   =   940;

    double dFitMaxE     =   530;
    double dFitMinE     =   495;


    TF1 *fFitPo = new TF1("fFitPo", "gaus(0)+gaus(3)", dFitMinPo, dFitMaxPo);
    TF1 *fFitTh = new TF1("fFitTh", "gaus(0)", dFitMinTh, dFitMaxTh);
    TF1 *fFitPt = new TF1("fFitPt", "gaus(0)", dFitMinPt, dFitMaxPt);
    TF1 *fFitK = new TF1("fFitK", "gaus(0)+pol3(3)",dFitMinK, dFitMaxK);

    TF1 *fFitTl = new TF1("fFitTl","gaus(0)", dFitMinTl, dFitMaxTl);
    TF1 *fFitSE = new TF1("fFitSE","gaus(0)", dFitMinSE, dFitMaxSE);
    TF1 *fFitDE = new TF1("fFitDE","gaus(0) + pol3(3)", dFitMinDE, dFitMaxDE);
    TF1 *fFitAc1 = new TF1("fFitAc1","[0]*exp(-pow(x-[1],2)/(2*[2]*[2])) + 3*[0]*exp(-pow(x-[3],2)/(2*[2]*[2])) + pol0(4)", dFitMinAc1, dFitMaxAc1); // Fit 968 with 2 gaussians
    TF1 *fFitAc2 = new TF1("fFitAc2","gaus(0) + pol1(3)", dFitMinAc2, dFitMaxAc2);
    TF1 *fFitE = new TF1("fFitE", "gaus(0) + pol3(3)", dFitMinE, dFitMaxE);


    // Initialize Parameters for the fit functions
    fFitPo->SetParameters(5000,5304,5,0,5407,5,0);
    fFitPo->SetLineColor(2);
    fFitTh->SetParameters(3,4078,5);
    fFitTh->SetLineColor(2);
    fFitPt->SetParameters(17,3249,7);
    fFitPt->SetLineColor(2);
    fFitK->SetParameters(5000,1460,5,0,0,0,0);
    fFitK->SetLineColor(2);

    fFitTl->SetParameters(100,2615,2.8);
    fFitTl->SetLineColor(2);
    fFitSE->SetParameters(5000,2104,5,0,0,0,0);
    fFitSE->SetLineColor(2);
    fFitDE->SetParameters(100,1593,3,0,0,0,0);
    fFitDE->SetLineColor(2);
    fFitAc1->SetParameters(20,965,1,968,0);
    fFitAc1->SetLineColor(2);
    fFitAc2->SetParameters(500,911,2,0,0);
    fFitAc2->SetLineColor(2);
    fFitE->SetParameters(500,511,3,0,0,0,0);
    fFitE->SetLineColor(2);


    vector<double> dEnergy, dEnergyErr, dResolution, dResolutionErr;

    // If parameters in fit change, must change the input here to make sure the graph is OK
    TCanvas *ce = new TCanvas("ce", "ce", 1100, 750);
    TH1D *he     = (TH1D*)hzoomed->Clone("he");
    he->SetTitle("Electron (511 keV)");
    he->SetAxisRange(dFitMinE, dFitMaxE);
    he->Fit("fFitE","R");
    dEnergy.push_back(fFitE->GetParameter(1));
    dEnergyErr.push_back(fFitE->GetParError(1));
    dResolution.push_back(TMath::Abs(fFitE->GetParameter(2)));
    dResolutionErr.push_back(TMath::Abs(fFitE->GetParError(2)));


    TCanvas *cac2 = new TCanvas("cac2", "cac2", 1100, 750);    
    TH1D *hac2     = (TH1D*)hzoomed->Clone("hac2");
    hac2->SetTitle("Ac228 (911 keV)");
    hac2->SetAxisRange(dFitMinAc2, dFitMaxAc2);
    hac2->Fit("fFitAc2","R");
    dEnergy.push_back(fFitAc2->GetParameter(1));
    dEnergyErr.push_back(fFitAc2->GetParError(1));
    dResolution.push_back(TMath::Abs(fFitAc2->GetParameter(2)));
    dResolutionErr.push_back(TMath::Abs(fFitAc2->GetParError(2)));

/*
    TCanvas *cac1 = new TCanvas("cac1", "cac1", 1100, 750);    
    // TH1D *hac1     = (TH1D*)hzoomed->Clone("hac1");
    hzoomed2->SetTitle("Ac228 (968 keV)");
    hzoomed2->SetAxisRange(dFitMinAc1-5, dFitMaxAc1+20);
    hzoomed2->Draw();
    hzoomed2->Fit("fFitAc1","R");
    dEnergy.push_back(fFitAc1->GetParameter(3));
    dEnergyErr.push_back(fFitAc1->GetParError(3));
    dResolution.push_back(TMath::Abs(fFitAc1->GetParameter(2)));
    dResolutionErr.push_back(TMath::Abs(fFitAc1->GetParError(2)));
*/


    TCanvas *ck = new TCanvas("ck", "ck", 1100, 750);
    TH1D *hk     = (TH1D*)hzoomed->Clone("hk");
    hk->SetTitle("K40 (1460 keV)");
    hk->SetAxisRange(dFitMinK, dFitMaxK);
    hk->Fit("fFitK","R");
    dEnergy.push_back(fFitK->GetParameter(1));
    dEnergyErr.push_back(fFitK->GetParError(1));
    dResolution.push_back(TMath::Abs(fFitK->GetParameter(2)));
    dResolutionErr.push_back(TMath::Abs(fFitK->GetParError(2)));


    TCanvas *ctl = new TCanvas("ctl", "ctl", 1100, 750);
    TH1D *htl     = (TH1D*)hzoomed->Clone("htl");
    htl->SetTitle("Tl208 (2615 keV)");
    htl->SetAxisRange(dFitMinTl, dFitMaxTl);
    htl->Fit("fFitTl","R");
    dEnergy.push_back(fFitTl->GetParameter(1));
    dEnergyErr.push_back(fFitTl->GetParError(1));
    dResolution.push_back(TMath::Abs(fFitTl->GetParameter(2)));
    dResolutionErr.push_back(TMath::Abs(fFitTl->GetParError(2)));


/*
    TCanvas *cpt = new TCanvas("cpt", "cpt", 1100, 750);
    TH1D *hpt     = (TH1D*)hzoomed->Clone("hpt");
    hpt->SetTitle("Pt190 (3249 keV)");
    // hpt->Rebin(2);
    hpt->SetAxisRange(dFitMinPt, dFitMaxPt);
    hpt->Fit("fFitPt","R");
    dEnergy.push_back(fFitPt->GetParameter(1));
    dEnergyErr.push_back(fFitPt->GetParError(1));
    dResolution.push_back(TMath::Abs(fFitPt->GetParameter(2)));
    dResolutionErr.push_back(TMath::Abs(fFitPt->GetParError(2)));


    TCanvas *cpo = new TCanvas("cpo", "cpo", 1100, 750);
    TH1D *hpo     = (TH1D*)hzoomed->Clone("hpo");
    hpo->SetTitle("Po210 (5304 and 5407 keV)");
    hpo->Rebin(2);
    hpo->SetAxisRange(dFitMinPo, dFitMaxPo+30);
    hpo->Fit("fFitPo","R");
    dEnergy.push_back(fFitPo->GetParameter(1));
    dEnergyErr.push_back(fFitPo->GetParError(1));
    dResolution.push_back(TMath::Abs(fFitPo->GetParameter(2)));
    dResolutionErr.push_back(TMath::Abs(fFitPo->GetParError(2)));

    dEnergy.push_back(fFitPo->GetParameter(4));
    dEnergyErr.push_back(fFitPo->GetParError(4));
    dResolution.push_back(TMath::Abs(fFitPo->GetParameter(5)));
    dResolutionErr.push_back(TMath::Abs(fFitPo->GetParError(5)));    
*/

    TCanvas *resolution = new TCanvas("resolution", "resolution", 1100, 750);
    int n = dEnergy.size();
    double *x = &dEnergy[0];
    double *errx = &dEnergyErr[0];
    // Resolution (FWHM as percentage of Energy)
    double *y = &dResolution[0];
    double *erry = &dResolutionErr[0];


    TF1 *fFitRes = new TF1("fFitRes","pol1(0)",300,5500);
    TGraphErrors *g1 = new TGraphErrors(n,x,y,errx, erry);
    g1->SetTitle("Resolution vs Energy");
    g1->GetXaxis()->SetTitle("Energy (keV)");
    g1->GetYaxis()->SetTitle("Resolution (keV)");
    // g1->SetLineColor(1);
    g1->Draw("A*");  
    g1->Fit("fFitRes","R");
    resolution->Update();



    TCanvas *cgamma = new TCanvas("cgamma", "cgamma", 1100, 750);
    hgamma->SetTitle("Energy Spectra (#gamma region)");
    hgamma->Draw();


    if(bSavePlots)
    {
        ctl->SaveAs("cTl208.png");
        cac1->SaveAs("cAc228-1.png");
        cac2->SaveAs("cAc228-2.png");
        ce->SaveAs("cE.png");
        cpo->SaveAs("cPo210.png");
        cpt->SaveAs("cPt190.png");
        ck->SaveAs("cK40.png");
    }
}


void AlphaRecalibration()
{
    // gStyle->SetOptFit();
    // Theoretical energy
    double x[7] = {3249, 4083, 4770, 4871, 5304, 5407, 5788};
    // double x[6] = {4083, 4770, 4871, 5304, 5407, 5788};

    // Energy in data
    double y[7] = {3291, 4116, 4826, 4898, 5349, 5441, 5802};
    // double y[6] = {4116, 4826, 4898, 5349, 5441, 5802};

    TGraph *g1 = new TGraph(7, y, x);
    g1->SetTitle("Alpha Recalibration");
    g1->GetXaxis()->SetTitle("Energy (Data) (keV)");
    g1->GetYaxis()->SetTitle("Theoretical Energy (keV)");
    g1->SetLineColor(1);
    g1->Fit("pol1");
    g1->Draw("A*");
}