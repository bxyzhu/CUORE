#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooRealVar.h"
#include "TChain.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"

using namespace RooFit;

TH1D *LoadMCHist(char* dSource, char* dLocation, int dMult)
{
    // Fixed range and bin number right now...
    // TH1D *MCHist = new TH1D(Form("h%s-%s-M%d",dSource, dLocation, dMult), Form("%s-%s-M%d", dSource, dLocation, dMult), 800, 0, 8000);
    TH1D *MCHist = new TH1D("MCHist", Form("%s-%s-M%d", dSource, dLocation, dMult), 800, 0, 8000);

    // Load MC files
    TChain *outTree = new TChain("outTree");
    outTree->Add(Form("/Users/brian/macros/Simulations/Bkg/%s-%s-M%d-T50-r0.0425-c0.10000.root", dSource, dLocation, dMult));

    float dEMin = 0, dEMax = 8000;
    TCut MC_cut = Form("Ener1 > %f && Ener1 < %f", dEMin, dEMax);

    outTree->Project("MCHist","Ener1",MC_cut);

    return MCHist;
}

void RooFitBkg()
{
    gStyle->SetOptStat();
	gStyle->SetOptFit();	

    int dMult = 1;
    float dEMin = 0, dEMax = 8000;

    // Load Data files

    // const int datasets[] = {2061, 2064, 2067, 2070, 2073, 2076};
    const int datasets[] = {2061};
    const int nDataset = sizeof(datasets)/sizeof(int);

    TChain *qtree = new TChain("qtree");
    qtree->Add("/Users/brian/macros/CUOREZ/Bkg/ReducedBkg-ds*.root");

/*
    for(int i = 0; i < nDataset; i++)
    {
        qtree->Add(Form("/Users/brian/macros/CUOREZ/Bkg/ReducedBkg-ds%d.root", datasets[i]));
    }
*/


    TCut base_cut;
    base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "abs(BaselineSlope)<0.1";
    base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";
    // base_cut = base_cut && Form("Multiplicity_OFTime==%d", dMult);
    // base_cut = base_cut && Form("Energy > %f && Energy < %f", dEMin, dEMax);


    TH1D *dataHistM1 = new TH1D("dataHistM1", "", 800, 0, 8000);
    qtree->Project("dataHistM1","Energy",base_cut && "Multiplicity_OFTime==1");


    RooRealVar Energy("Energy", "Energy", 0, 8000);
/*
    // This way doesn't seem to work well...
    RooRealVar Multiplicity_OFTime("Multiplicity_OFTime", "Multiplicity_OFTime", 0, 5);
    RooRealVar BaselineSlope("BaselineSlope", "BaselineSlope", -10.0, 10.0);
    RooRealVar OF_TVR("OF_TVR", "OF_TVR", -10.0, 10.0);
    RooRealVar OF_TVL("OF_TVL", "OF_TVL", -10.0, 10.0);
    RooRealVar TimeUntilSignalEvent_SameChannel("TimeUntilSignalEvent_SameChannel","TimeUntilSignalEvent_SameChannel",-1000.0,100000.0);
    RooRealVar TimeSinceSignalEvent_SameChannel("TimeSinceSignalEvent_SameChannel","TimeSinceSignalEvent_SameChannel",-1000.0,100000.0);

    RooDataSet *mult1Data = new RooDataSet("mult1Data", "Multiplicity 1 Data", 
            // RooArgSet(Energy, Multiplicity_OFTime, BaselineSlope, OF_TVR, OF_TVL),
            RooArgSet(Energy, Multiplicity_OFTime, BaselineSlope, OF_TVR, OF_TVL, TimeUntilSignalEvent_SameChannel, TimeSinceSignalEvent_SameChannel), 
            Import(*qtree), 
            Cut("Multiplicity_OFTime==1 && abs(BaselineSlope<0.1) && OF_TVR < 1.75 && OF_TVL < 2.05 && (TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0) && (TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)")
        );
*/


    RooDataHist *dataBkgM1 = new RooDataHist("dataBkgM1", "Data from Bkg M1", Energy, dataHistM1);

    // Load MC histograms
    TH1D *hM1_1 = LoadMCHist("Crystal-U238","B",1);
    TH1D *hM1_2 = LoadMCHist("Crystal-U238","S100",1);
    TH1D *hM1_3 = LoadMCHist("Crystal-U238","S10",1);
    TH1D *hM1_4 = LoadMCHist("Crystal-U238","S1",1);
    TH1D *hM1_5 = LoadMCHist("Crystal-U238","S01",1);

/*
    TH1D *hM2_1 = LoadMCHist("Crystal-U238","B",2);
    TH1D *hM2_2 = LoadMCHist("Crystal-U238","S100",2);
    TH1D *hM2_3 = LoadMCHist("Crystal-U238","S10",2);
    TH1D *hM2_4 = LoadMCHist("Crystal-U238","S1",2);
    TH1D *hM2_5 = LoadMCHist("Crystal-U238","S01",2);
*/
    // 
    RooRealVar n1("n1","norm1", 0.1, 0, 1);
    RooRealVar n2("n2","norm2", 0.1, 0, 1);
    RooRealVar n3("n3","norm3", 0.1, 0, 1);
    RooRealVar n4("n4","norm4", 0.1, 0, 1);
    RooRealVar n5("n5","norm5", 0.1, 0, 1);

    // Make histograms into RooFit data
    RooDataHist *dM1_1 = new RooDataHist("dM1_1", "Data from Montecarlo", Energy, hM1_1);
    RooDataHist *dM1_2 = new RooDataHist("dM1_2", "Data from Montecarlo", Energy, hM1_2);
    RooDataHist *dM1_3 = new RooDataHist("dM1_3", "Data from Montecarlo", Energy, hM1_3);
    RooDataHist *dM1_4 = new RooDataHist("dM1_4", "Data from Montecarlo", Energy, hM1_4);
    RooDataHist *dM1_5 = new RooDataHist("dM1_5", "Data from Montecarlo", Energy, hM1_5);

    // Make data into PDFs
    RooHistPdf *pdfMC_1 = new RooHistPdf("pdfMC", "PDF from Montecarlo", Energy, *dM1_1);
    RooHistPdf *pdfMC_2 = new RooHistPdf("pdfMC", "PDF from Montecarlo", Energy, *dM1_2);
    RooHistPdf *pdfMC_3 = new RooHistPdf("pdfMC", "PDF from Montecarlo", Energy, *dM1_3);
    RooHistPdf *pdfMC_4 = new RooHistPdf("pdfMC", "PDF from Montecarlo", Energy, *dM1_4);
    RooHistPdf *pdfMC_5 = new RooHistPdf("pdfMC", "PDF from Montecarlo", Energy, *dM1_5);


    // RooAddPdf *model = new RooAddPdf("model", "Background Model", RooArgList(*pdfMC_1, *pdfMC_2, *pdfMC_3, *pdfMC_4, *pdfMC_5), 
                                                            // RooArgList(n1, n2, n3, n4, n5));

    RooAddPdf *model = new RooAddPdf("model", "Background Model", RooArgList(*pdfMC_1, *pdfMC_2), RooArgList(n1, n2));

/*
    RooDataSet *mult2Data = new RooDataSet("mult2Data", "Multiplicity 2 Data", 
        RooArgSet(Energy, Multiplicity_OFTime, BaselineSlope, OF_TVR, OF_TVL, TimeUntilSignalEvent_SameChannel, TimeSinceSignalEvent_SameChannel), 
            Import(*qtree) 
            // Cut("Multiplicity_OFTime==2 && abs(BaselineSlope<0.1) && OF_TVR < 1.75 && OF_TVL < 2.05 && (TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0) && (TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)")
        );

    // Combined data
    RooDataSet *combinedData = new RooDataSet("combinedData", "Combined Data", Energy, Index(*c), 
                                    Import("Multiplicity1", *mult1Data), Import("Multiplicity2", *mult2Data));
*/


    // data statistics

	TCanvas *c1 = new TCanvas("c1","c1",1000,600);
	gPad->SetLogy();
	RooPlot *frame = Energy.frame(Title("Background Fit"));
	dataBkgM1->plotOn(frame);
	// model->fitTo(*dataBkgM1, Range(4500, 6000));
	// model->plotOn(frame, LineColor(kRed));
    // model->plotOn(frame, Components(RooArgSet(*pdfMC_1, *pdfMC_2, *pdfMC_3, *pdfMC_4, *pdfMC_5)), LineStyle(kDashed));
    // model->plotOn(frame, Components(RooArgSet(*pdfMC_1,*pdfMC_2)), LineStyle(kDashed));

	// pdfMC->Print("t");

    // Add stats to plot
    // pdfMC->paramOn(frame,Layout(0.65,99,0.95));
    // dataCalib->statOn(frame,Layout(0.65,0.99,0.95));


	frame->Draw();

/*
    // Testing stuff
    TCanvas *c2 = new TCanvas("c2","c2",1000,1200);
    c2->Divide(1,2);

    // Residuals
    RooHist *hresid = frame->residHist();

    //
    RooHist *hpull = frame->pullHist();

    RooPlot *frame2 = Energy.frame(Title("Residual Distribution"));
    frame2->addPlotable(hresid,"P");

    RooPlot *frame3 = Energy.frame(Title("Pull Distribution"));
    frame3->addPlotable(hpull,"P");

    c2->cd(1);
    frame2->Draw();
    c2->cd(2);
    frame3->Draw();

    cout << "chi^2 = " << frame->chiSquare() << endl;
*/


}

void TestHists()
{

    TH1D *hM1_1 = LoadMCHist("Crystal-U238","B",1);
    TH1D *hM1_2 = LoadMCHist("Crystal-U238","S100",1);

    TCanvas *c1 = new TCanvas("c1");
    c1->Divide(2,1);
    c1->cd(1);
    hM1_1->Draw();
    c1->cd(2);
    hM1_2->Draw();
}
