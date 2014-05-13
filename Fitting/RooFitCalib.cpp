#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooPlot.h"
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




void RooFitCalib()
{
	gStyle->SetOptFit();	

    int dMult = 1;
    float dEMin = 0, dEMax = 2700;

    TH1D *calibHist = new TH1D("calibHist", "", 350, 0, 3500);
    TH1D *MCHist = new TH1D("MCHist", "", 350, 0, 3500);


    // Load MC files
    TChain *outTree = new TChain("outTree");
    outTree->Add(Form("/Users/brian/macros/Simulations/Calib/Calib-C0-M%d-T300-r3.2537-c0.10000.root", dMult));


    // Load Calibration files
    TChain *qtree = new TChain("qtree");
    // qtree->Add("/Users/brian/data/CUOREZ/BlindedReduced_200730_C_p001.root");
    // qtree->Add("/Users/brian/data/CUOREZ/BlindedReduced_200757_C_p001.root");
    qtree->Add("/Users/brian/macros/CUOREZ/Bkg/ReducedCalib-ds2061.root");

    TCut base_cut;
    // base_cut = base_cut && "Filter_RejectBadIntervals";
    // base_cut = base_cut && "NumberOfPulses==1";
    base_cut = base_cut && "(TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "(TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)";
    base_cut = base_cut && "abs(BaselineSlope)<0.1";
    base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";
    base_cut = base_cut && Form("Multiplicity_OFTime==%d", dMult);
    base_cut = base_cut && Form("Energy > %f && Energy < %f", dEMin, dEMax);

    TCut MC_cut = Form("Ener1 > %f && Ener1 < %f", dEMin, dEMax);


    qtree->Project("calibHist","Energy",base_cut);
    outTree->Project("MCHist","Ener1",MC_cut);

    // RooFit stuff
    RooRealVar Energy("Energy", "Energy", 0, 2700);
    RooRealVar Multiplicity_OFTime("Multiplicity_OFTime", "Multiplicity_OFTime", 0, 5);
    RooRealVar BaselineSlope("BaselineSlope", "BaselineSlope", -5, 5);
    RooRealVar OF_TVR("OF_TVR", "OF_TVR", -5, 5);
    RooRealVar OF_TVL("OF_TVL", "OF_TVL", -5, 5);
    RooRealVar TimeUntilSignalEvent_SameChannel("TimeUntilSignalEvent_SameChannel","TimeUntilSignalEvent_SameChannel",-5,5);
    RooRealVar TimeSinceSignalEvent_SameChannel("TimeSinceSignalEvent_SameChannel","TimeSinceSignalEvent_SameChannel",-5,5);

    RooDataSet *dataCalib = new RooDataSet("dataCalib", "Data from Calibration", 
        RooArgSet(Energy, Multiplicity_OFTime, BaselineSlope, OF_TVR, OF_TVL, TimeUntilSignalEvent_SameChannel, TimeSinceSignalEvent_SameChannel), 
            Import(*qtree), 
            // Cut("Multiplicity_OFTime==1","abs(BaselineSlope<0.1)") // Doesn't work
            Cut("Multiplicity_OFTime==1 && abs(BaselineSlope<0.1) && OF_TVR < 1.75 && OF_TVL < 2.05 && (TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0) && (TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0)")
        );


    dataCalib->Print();
    // RooDataHist *dataCalib = new RooDataHist("dataCalib", "Data from Calibration", x, calibHist);

    RooDataHist *dataMC = new RooDataHist("dataMC", "Data from Montecarlo", Energy, MCHist);
    RooHistPdf *pdfMC = new RooHistPdf("pdfMC", "PDF from Montecarlo", Energy, *dataMC);

    // RooRealVar para("para", "parameter", 0.0, 1.0);

    // data statistics

	TCanvas *c1 = new TCanvas("c1","c1",1000,600);
	gPad->SetLogy();
	RooPlot *frame = Energy.frame(Title("Calibration Fit"));
	dataCalib->plotOn(frame,DataError(RooAbsData::SumW2));
	pdfMC->fitTo(*dataCalib, Range(300, 2700));
	pdfMC->plotOn(frame);

	pdfMC->Print("t");

    // Add stats to plot
    pdfMC->paramOn(frame,Layout(0.65,99,0.95));
    // dataCalib->statOn(frame,Layout(0.65,0.99,0.95));


	frame->Draw();

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


}
