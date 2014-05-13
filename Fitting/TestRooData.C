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

/////////////// Input U238 files
    input.open("Crystal-U238-Bulk-M1.txt");
    double number;
    int i1 = 0;
    TH1D *hu1 = new TH1D("hu1", "U238 Crystal Bulk", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	hu1->Fill(i1,number);
	i1++;	
    }
    input.close(); 

    input.open("Crystal-U238-Surf10-M1.txt");
    int i2 = 0;
    TH1D *hu2 = new TH1D("hu2", "U238 Crystal Surf 10", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	hu2->Fill(i2,number);
	i2++;	
    }
    input.close(); 

    input.open("Crystal-U238-Surf1-M1.txt");
    int i3 = 0;
    TH1D *hu3 = new TH1D("hu3", "U238 Crystal Surf 1", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	hu3->Fill(i3,number);
	i3++;	
    }
    input.close();   

    input.open("Crystal-U238-Surf01-M1.txt");
    int i4 = 0;
    TH1D *hu4 = new TH1D("hu4", "U238 Crystal Surf 0.1", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	hu4->Fill(i4,number);
	i4++;	
    }
    input.close();   

    input.open("Frame-U238-Bulk-M1.txt");
    int i5 = 0;
    TH1D *hu5 = new TH1D("hu5", "U238 Frame Bulk", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	hu5->Fill(i5,number);
	i5++;	
    }
    input.close();  

    input.open("Frame-U238-Surf10-M1.txt");
    int i6 = 0;
    TH1D *hu6 = new TH1D("hu6", "U238 Frame Surf 10", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	hu6->Fill(i6,number);
	i6++;	
    }
    input.close();  

    input.open("Frame-U238-Surf1-M1.txt");
    int i7 = 0;
    TH1D *hu7 = new TH1D("hu7", "U238 Frame Surf 1", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	hu7->Fill(i7,number);
	i7++;	
    }
    input.close();  

    input.open("Frame-U238-Surf01-M1.txt");
    int i8 = 0;
    TH1D *hu8 = new TH1D("hu8", "U238 Frame Surf 0.1", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	hu8->Fill(i8,number);
	i8++;	
    }
    input.close(); 


///////////// Input Th232 files

    input.open("Crystal-Th232-Bulk-M1.txt");
    int ii1 = 0;
    TH1D *ht1 = new TH1D("ht1", "Th232 Crystal Bulk", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	ht1->Fill(ii1,number);
	ii1++;	
    }
    input.close(); 

    input.open("Crystal-Th232-Surf10-M1.txt");
    int ii2 = 0;
    TH1D *ht2 = new TH1D("ht2", "Th232 Crystal Surf 10", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	ht2->Fill(ii2,number);
	ii2++;	
    }
    input.close(); 

    input.open("Crystal-Th232-Surf1-M1.txt");
    int ii3 = 0;
    TH1D *ht3 = new TH1D("ht3", "Th232 Crystal Surf 1", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	ht3->Fill(ii3,number);
	ii3++;	
    }
    input.close();   

    input.open("Crystal-Th232-Surf01-M1.txt");
    int ii4 = 0;
    TH1D *ht4 = new TH1D("ht4", "Th232 Crystal Surf 0.1", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	ht4->Fill(ii4,number);
	ii4++;	
    }
    input.close();   

    input.open("Frame-Th232-Bulk-M1.txt");
    int ii5 = 0;
    TH1D *ht5 = new TH1D("ht5", "Th232 Frame Bulk", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	ht5->Fill(ii5,number);
	ii5++;	
    }
    input.close();  

    input.open("Frame-Th232-Surf10-M1.txt");
    int ii6 = 0;
    TH1D *ht6 = new TH1D("ht6", "Th232 Frame Surf 10", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	ht6->Fill(ii6,number);
	ii6++;	
    }
    input.close();  

    input.open("Frame-Th232-Surf1-M1.txt");
    int ii7 = 0;
    TH1D *ht7 = new TH1D("ht7", "Th232 Frame Surf 1", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	ht7->Fill(ii7,number);
	ii7++;
    }
    input.close();  

    input.open("Frame-Th232-Surf01-M1.txt");
    int ii8 = 0;
    TH1D *ht8 = new TH1D("ht8", "Th232 Frame Surf 0.1", 10000, 0, 10000);
    while(input)
    {
	input >> number;
	ht8->Fill(ii8,number);
	ii8++;
    }
    input.close(); 

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


/*
void LoadThData()
{
    // Problem... only loading 1 MC file at a time?
//    va_list files;
//    va_start (files, numFiles);

    vector<char> inputList;

    ifstream input;
    input.open("Crystal-Th232-Bulk-M1.txt");
    double number1;
    int i1 = 0;
    TH1D *hu1 = new TH1D("hu1", "Th232 Crystal Bulk", 10000, 0, 10000);
    while(input)
    {
	input >> number1;
	hu1->Fill(i1,number1);
	i1++;	
    }
    input.close(); 

    input.open("Crystal-Th232-Surf10-M1.txt");
    double number2;
    int i2 = 0;
    TH1D *hu2 = new TH1D("hu2", "Th232 Crystal Surf 10", 10000, 0, 10000);
    while(input)
    {
	input >> number2;
	hu2->Fill(i2,number2);
	i2++;	
    }
    input.close(); 


    input.open("Crystal-Th232-Surf1-M1.txt");
    double number3;
    int i3 = 0;
    TH1D *hu3 = new TH1D("hu3", "Th232 Crystal Surf 1", 10000, 0, 10000);
    while(input)
    {
	input >> number3;
	hu3->Fill(i3,number3);
	i3++;	
    }
    input.close();   

    input.open("Crystal-Th232-Surf01-M1.txt");
    double number4;
    int i4 = 0;
    TH1D *hu4 = new TH1D("hu4", "Th232 Crystal Surf 0.1", 10000, 0, 10000);
    while(input)
    {
	input >> number4;
	hu4->Fill(i4,number4);
	i4++;	
    }
    input.close();   


    input.open("Frame-Th232-Bulk-M1.txt");
    double number5;
    int i5 = 0;
    TH1D *h5 = new TH1D("h5", "Th232 Frame Bulk", 10000, 0, 10000);
    while(input)
    {
	input >> number5;
	h5->Fill(i5,number5);
	i5++;	
    }
    input.close();  


    input.open("Frame-Th232-Surf10-M1.txt");
    double number6;
    int i6 = 0;
    TH1D *h6 = new TH1D("h6", "Th232 Frame Surf 10", 10000, 0, 10000);
    while(input)
    {
	input >> number6;
	h6->Fill(i6,number6);
	i6++;	
    }
    input.close();  


    input.open("Frame-Th232-Surf1-M1.txt");
    double number7;
    int i7 = 0;
    TH1D *h7 = new TH1D("h7", "Th232 Frame Surf 1", 10000, 0, 10000);
    while(input)
    {
	input >> number7;
	h7->Fill(i7,number7);
	i7++;
    }
    input.close();  


    input.open("Frame-Th232-Surf01-M1.txt");
    double number8;
    int i8 = 0;
    TH1D *h8 = new TH1D("h8", "Th232 Frame Surf 0.1", 10000, 0, 10000);
    while(input)
    {
	input >> number8;
	h8->Fill(i8,number8);
	i8++;
    }
    input.close(); 

    // Create variables
    RooRealVar x("x","x",0, 10000);

    // Inport data from histograms
    RooDataHist *data1 = new RooDataHist("data1", "Data Th232 Crystal Bulk",x, h1);
    RooDataHist *data2 = new RooDataHist("data2", "Data Th232 Crystal Surf 10",x, h2);
    RooDataHist *data3 = new RooDataHist("data3", "Data Th232 Crystal Surf 1",x, h3);
//    RooDataHist *data4 = new RooDataHist("data4", "Data Th232 Crystal Surf 0.1",x, h4);
    RooDataHist *data5 = new RooDataHist("data5", "Data Th232 Frame Bulk",x, h5);
    RooDataHist *data6 = new RooDataHist("data6", "Data Th232 Frame Surf 10",x, h6);
    RooDataHist *data7 = new RooDataHist("data7", "Data Th232 Frame Surf 1",x, h7);
//    RooDataHist *data8 = new RooDataHist("data8", "Data Th232 Frame Surf 0.1",x, h8);

    // Convert histograms to PDF models
    RooHistPdf *pdfMC1 = new RooHistPdf("pdfMC1", "PDF Th232 Crystal Bulk", x, *data1);
    RooHistPdf *pdfMC2 = new RooHistPdf("pdfMC2", "PDF Th232 Crystal Surf 10", x, *data2);
    RooHistPdf *pdfMC3 = new RooHistPdf("pdfMC3", "PDF Th232 Crystal Surf 1", x, *data3);
//    RooHistPdf *pdfMC4 = new RooHistPdf("pdfMC4", "PDF Th232 Crystal Surf 0.1", x, *data4);
    RooHistPdf *pdfMC5 = new RooHistPdf("pdfMC5", "PDF Th232 Frame Bulk", x, *data5);
    RooHistPdf *pdfMC6 = new RooHistPdf("pdfMC6", "PDF Th232 Frame Surf 10", x, *data6);
    RooHistPdf *pdfMC7 = new RooHistPdf("pdfMC7", "PDF Th232 Frame Surf 1", x, *data7);
//    RooHistPdf *pdfMC8 = new RooHistPdf("pdfMC8", "PDF Th232 Frame Surf 0.1", x, *data8);

    // Create fractional variables (with initial vars for now)
    RooRealVar fbkg1("fbkg1", "background fraction 1", 0., 1.);
    RooRealVar fbkg2("fbkg2", "background fraction 2", 0., 1.);
    RooRealVar fbkg3("fbkg3", "background fraction 3", 0., 1.);
    RooRealVar fbkg4("fbkg4", "background fraction 4", 0., 1.);
    RooRealVar fbkg5("fbkg5", "background fraction 5", 0., 1.);
//    RooRealVar fbkg6("fbkg6", "background fraction 6", 0., 1.);
//    RooRealVar fbkg7("fbkg7", "background fraction 7", 0., 1.);

    RooAbsReal *intMC1 = pdfMC1->createIntegral(x);
    RooAbsReal *intMC2 = pdfMC2->createIntegral(x);
    RooAbsReal *intMC3 = pdfMC3->createIntegral(x);
//    RooAbsReal *intMC4 = pdfMC4->createIntegral(x);
    RooAbsReal *intMC5 = pdfMC5->createIntegral(x);
    RooAbsReal *intMC6 = pdfMC6->createIntegral(x);
    RooAbsReal *intMC7 = pdfMC7->createIntegral(x);
//    RooAbsReal *intMC8 = pdfMC8->createIntegral(x);

    // Add together PDFs (not recursively for now)
    RooAddPdf modelTot("modelTot", "Total U238 Model", 
	RooArgList(*pdfMC1, *pdfMC2, *pdfMC3, *pdfMC5, *pdfMC6, *pdfMC7), 
	RooArgList(fbkg1, fbkg2, fbkg3, fbkg4, fbkg5), kTRUE);

    RooAbsReal *intTotMC = modelTot.createIntegral(x);



    cout << "Th232 Crystal Bulk: " << intMC1->getVal() << endl;
    cout << "Th232 Crystal Surf 10: " << intMC2->getVal() << endl;
    cout << "Th232 Crystal Surf 1: " << intMC3->getVal() << endl;
//    cout << "Th232 Crystal Surf 0.1: " << intMC4->getVal() << endl;
    cout << "Th232 Frame Bulk: " << intMC5->getVal() << endl;
    cout << "Th232 Frame Surf 10: " << intMC6->getVal() << endl;
    cout << "Th232 Frame Surf 1: " << intMC7->getVal() << endl;
//    cout << "Th232 Frame Surf 0.1: " << intMC8->getVal() << endl;
    cout << "Total Integral: " << intTotMC->getVal() << endl;

    // Draw Model
    TCanvas * c1 = new TCanvas("c1","PDFs converted from Geant4 data", 10,10,700,700);
    c1->Divide(1,2);
    c1->cd(1);
    gPad->SetLogy();
    // Plotting the model
    RooPlot* frame1 = x.frame();
    modelTot.plotOn(frame1, LineColor(kBlack));
    frame1->Draw();

    c1->cd(2);
    gPad->SetLogy();
    RooPlot* frame2 = x.frame();
    modelTot.plotOn(frame2, Components("pdfMC1"), LineColor(2), LineStyle(2));
    modelTot.plotOn(frame2, Components("pdfMC2"), LineColor(3), LineStyle(3));
    modelTot.plotOn(frame2, Components("pdfMC3"), LineColor(4), LineStyle(4));
//    modelTot.plotOn(frame2, Components("pdfMC4"), LineColor(5), LineStyle(5));
    modelTot.plotOn(frame2, Components("pdfMC5"), LineColor(6), LineStyle(6));
    modelTot.plotOn(frame2, Components("pdfMC6"), LineColor(7), LineStyle(7));
    modelTot.plotOn(frame2, Components("pdfMC7"), LineColor(8), LineStyle(8));
//    modelTot.plotOn(frame2, Components("pdfMC8"), LineColor(9), LineStyle(9));
    frame2->Draw();

}
*/