// This macro should load at least 2 MC histogram files and try to plot them
#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "TChain.h"
#include "TCut.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include <cmath>
#include <iostream>
#include <map>
#include <vector>


using namespace RooFit;

//void LoadMCData(int numFiles, const char* fileName, ...)
void LoadDataAndPlot()
{
    // Problem... only loading 1 MC file at a time?
//    va_list files;
//    va_start (files, numFiles);

    ifstream input; 

/////////////// Calibration files
    


    
/////////////////////////////////////////////////////////////////////////////////



    // Create variables
    RooRealVar x("x","x",0, 10000);

    // Inport U238 data from histograms
    RooDataHist *datau1 = new RooDataHist("datau1", "Data U238 Crystal Bulk",x, hu1);
    RooDataHist *datau2 = new RooDataHist("datau2", "Data U238 Crystal Surf 10",x, hu2);
    RooDataHist *datau3 = new RooDataHist("datau3", "Data U238 Crystal Surf 1",x, hu3);
    RooDataHist *datau4 = new RooDataHist("datau4", "Data U238 Crystal Surf 0.1",x, hu4);
    RooDataHist *datau5 = new RooDataHist("datau5", "Data U238 Frame Bulk",x, hu5);
    RooDataHist *datau6 = new RooDataHist("datau6", "Data U238 Frame Surf 10",x, hu6);
    RooDataHist *datau7 = new RooDataHist("datau7", "Data U238 Frame Surf 1",x, hu7);
    RooDataHist *datau8 = new RooDataHist("datau8", "Data U238 Frame Surf 0.1",x, hu8);

    // Convert histograms to PDF models
    RooHistPdf *pdfMCu1 = new RooHistPdf("pdfMCu1", "PDF U238 Crystal Bulk", x, *datau1);
    RooHistPdf *pdfMCu2 = new RooHistPdf("pdfMCu2", "PDF U238 Crystal Surf 10", x, *datau2);
    RooHistPdf *pdfMCu3 = new RooHistPdf("pdfMCu3", "PDF U238 Crystal Surf 1", x, *datau3);
    RooHistPdf *pdfMCu4 = new RooHistPdf("pdfMCu4", "PDF U238 Crystal Surf 0.1", x, *datau4);
    RooHistPdf *pdfMCu5 = new RooHistPdf("pdfMCu5", "PDF U238 Frame Bulk", x, *datau5);
    RooHistPdf *pdfMCu6 = new RooHistPdf("pdfMCu6", "PDF U238 Frame Surf 10", x, *datau6);
    RooHistPdf *pdfMCu7 = new RooHistPdf("pdfMCu7", "PDF U238 Frame Surf 1", x, *datau7);
    RooHistPdf *pdfMCu8 = new RooHistPdf("pdfMCu8", "PDF U238 Frame Surf 0.1", x, *datau8);

    // Create fractional variables (with initial vars for now)
    RooRealVar fbkgu1("fbkgu1", "background fraction 1", 0., 1.);
    RooRealVar fbkgu2("fbkgu2", "background fraction 2", 0., 1.);
    RooRealVar fbkgu3("fbkgu3", "background fraction 3", 0., 1.);
    RooRealVar fbkgu4("fbkgu4", "background fraction 4", 0., 1.);
    RooRealVar fbkgu5("fbkgu5", "background fraction 5", 0., 1.);
    RooRealVar fbkgu6("fbkgu6", "background fraction 6", 0., 1.);
    RooRealVar fbkgu7("fbkgu7", "background fraction 7", 0., 1.);

    // Add together PDFs
    RooAddPdf modelUTot("modelUTot", "Total U238 Model", 
	RooArgList(*pdfMCu1, *pdfMCu2, *pdfMCu3, *pdfMCu4, *pdfMCu5, *pdfMCu6, *pdfMCu7, *pdfMCu8), 
	RooArgList(fbkgu1, fbkgu2, fbkgu3, fbkgu4, fbkgu5, fbkgu6, fbkgu7), kTRUE);




    // Inport Th232 data from histograms
    RooDataHist *datat1 = new RooDataHist("datat1", "Data Th232 Crystal Bulk",x, ht1);
    RooDataHist *datat2 = new RooDataHist("datat2", "Data Th232 Crystal Surf 10",x, ht2);
    RooDataHist *datat3 = new RooDataHist("datat3", "Data Th232 Crystal Surf 1",x, ht3);
    RooDataHist *datat5 = new RooDataHist("datat5", "Data Th232 Frame Bulk",x, ht5);
    RooDataHist *datat6 = new RooDataHist("datat6", "Data Th232 Frame Surf 10",x, ht6);
    RooDataHist *datat7 = new RooDataHist("datat7", "Data Th232 Frame Surf 1",x, ht7);

    // Convert histograms to PDF models
    RooHistPdf *pdfMCt1 = new RooHistPdf("pdfMCt1", "PDF Th232 Crystal Bulk", x, *datat1);
    RooHistPdf *pdfMCt2 = new RooHistPdf("pdfMCt2", "PDF Th232 Crystal Surf 10", x, *datat2);
    RooHistPdf *pdfMCt3 = new RooHistPdf("pdfMCt3", "PDF Th232 Crystal Surf 1", x, *datat3);
    RooHistPdf *pdfMCt5 = new RooHistPdf("pdfMCt5", "PDF Th232 Frame Bulk", x, *datat5);
    RooHistPdf *pdfMCt6 = new RooHistPdf("pdfMCt6", "PDF Th232 Frame Surf 10", x, *datat6);
    RooHistPdf *pdfMCt7 = new RooHistPdf("pdfMCt7", "PDF Th232 Frame Surf 1", x, *datat7);

    // Create fractional variables (with initial vars for now)
    RooRealVar fbkgt1("fbkgt1", "background fraction 1", 0., 1.);
    RooRealVar fbkgt2("fbkgt2", "background fraction 2", 0., 1.);
    RooRealVar fbkgt3("fbkgt3", "background fraction 3", 0., 1.);
    RooRealVar fbkgt4("fbkgt4", "background fraction 4", 0., 1.);
    RooRealVar fbkgt5("fbkgt5", "background fraction 5", 0., 1.);

    // Add together PDFs
    RooAddPdf modelThTot("modelThTot", "Total Th232 Model", 
	RooArgList(*pdfMCt1, *pdfMCt2, *pdfMCt3, *pdfMCt5, *pdfMCt6, *pdfMCt7), 
	RooArgList(fbkgt1, fbkgt2, fbkgt3, fbkgt4, fbkgt5), kTRUE);

    RooRealVar fbkgTot("fbkgTot", "Total background fraction", 0., 1.);

    RooAddPdf modelTot("modelTot", "Total Model (Th + U)",
	RooArgList(modelUTot, modelThTot),
	RooArgList(fbkgTot), kTRUE);



//////////////////////////////////////////////////////////



    const char* inDir = "/home/brian/data/Parylene/";
    int maxE = 4000;
    int minE = 1000;
    // Load files
    TChain *qtree = new TChain("qtree");
    qtree->Add(Form("%s/Reduced_*_B_*.root",inDir));
    
/*    vector<int> ch;
    int chanList[] = {2, 4, 6, 7, 8, 10, 12, 13, 15, 16};
    for (unsigned int i = 0; i<sizeof(chanList)/sizeof(chanList[0]); i++)
    {
	ch.push_back(chanList[i]);
	cout << "channels[" << i << "] = " << ch[i] << endl;
    }
*/

    int bin = (maxE - minE)/5;
    TH1F *dataHist = new TH1F("dataHist","dataHist", bin, minE, maxE);

    TCut cuts("IsSignal");
    cuts = cuts && "Filter_ReTrigger";
    cuts = cuts && "Filter_RejectBadIntervals";
    cuts = cuts && "NumberOfPulses == 1";
//    cuts = cuts && Form("Channel==%d",chan);
    cuts = cuts && "OF_SecondAmplitude==0";
    cuts = cuts && "abs(BaselineSlope<0.1)";
//    cuts = cuts && Form("OF_TVR<%f",TVR[chan]);
//    cuts = cuts && Form("OF_TVL<%f",TVL[chan]);
    cuts = cuts && "Multiplicity==1";
    cuts = cuts && Form("Energy>%d", minE);
    cuts = cuts && Form("Energy<%d", maxE);
    cuts = cuts && "PulserFlagByRegularTiming==0";
    qtree->Project("dataHist","Energy",cuts);

    // Create variables (not needed if I load model?)
//    RooRealVar x("x","x",0, 10000);

    // Convert data histogram into RooFit variable
    RooDataHist *dataP = new RooDataHist("dataP", "Data from Parylene Run",x,dataHist);

    TCanvas * c1 = new TCanvas("c1","PDFs converted from Geant4 data", 10,10,700,700);
    gPad->SetLogy();    
    RooPlot* frame = x.frame();
    dataP->plotOn(frame);
    modelTot.fitTo(*dataP, Range(1000, 4000));
    modelTot.plotOn(frame);
    frame->Draw();

}
