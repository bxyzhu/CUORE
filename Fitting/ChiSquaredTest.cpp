#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TCut.h"
#include "THStack.h"
#include "TMath.h"
#include "TH2F.h"
#include "TLegend.h"

#include "TRandom1.h"

#include <TMinuitMinimizer.h>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <cmath>
#include <string>

using namespace std;


double CalibChiSquare(const double* para)
{
    // string MCDir = "/Users/brian/macros/Simulations/Calib/";
    // string CalibDir = "/Users/brian/data/CUOREZ/";

    TH1D *calibHist = new TH1D("calibHist", "", 2700, 0, 2700);
    TH1D *MCHist = new TH1D("MCHist", "", 2700, 0, 2700);

    int dMult = 1;
    int dEMin=0, dEMax=2700;
    // const int dTotBinNbr = calibHist->GetSize();

    float CalibData[2700];
    float CalibPDF[2700];
    float TotalPDF[2700];

    // Load MC files
    TChain *outTree = new TChain("outTree");
    outTree->Add(Form("/Users/brian/macros/Simulations/Calib/Calib-C0-M%d-T300-r3.2537-c0.10000.root", dMult));



    // Load Calibration files
    TChain *qtree = new TChain("qtree");
    qtree->Add("/Users/brian/data/CUOREZ/BlindedReduced_200730_C_p001.root");
    qtree->Add("/Users/brian/data/CUOREZ/BlindedReduced_200757_C_p001.root");


    TCut base_cut("IsSignal");
    base_cut = base_cut && "Filter_RejectBadIntervals";
    base_cut = base_cut && "NumberOfPulses==1";
    base_cut = base_cut && "TimeUntilSignalEvent_SameChannel > 4.0 || TimeUntilSignalEvent_SameChannel < 0";
    base_cut = base_cut && "TimeSinceSignalEvent_SameChannel > 3.1 || TimeSinceSignalEvent_SameChannel < 0";
    base_cut = base_cut && "abs(BaselineSlope)<0.1";
    base_cut = base_cut && "OF_TVR < 1.75 && OF_TVL < 2.05";
    base_cut = base_cut && Form("Multiplicity_OFTime==%d", dMult);
    base_cut = base_cut && Form("Energy > %f && Energy < %f", dEMin, dEMax);

    TCut MC_cut = Form("Ener1 > %f && Ener1 < %f", dEMin, dEMax);



    qtree->Project("calibHist","Energy",base_cut);
    outTree->Project("MCHist","Ener1",MC_cut);


    double ChiSquare = 0;


    for(int i=2000; i<=2600; i++)
    {
        // Get value from histograms
        CalibData[i] = calibHist->GetBinContent(i);
        CalibPDF[i] = MCHist->GetBinContent(i);

        TotalPDF[i] = para[0]*CalibPDF[i];

        ChiSquare += TMath::Power(
            CalibData[i] - TotalPDF[i], 2)
            /(CalibData[i] + TotalPDF[i]);
    }

    delete calibHist;
    delete MCHist;
    
    return ChiSquare;
}

int ChiSquaredTest()
{
    
    ofstream Output("Minuit_Fitting_Result.txt");



   ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
	// ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Simplex");
//    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Minuit2");

	min->SetMaxFunctionCalls(50);
	min->SetMaxIterations(50);

	// Adjust tolerance to make fit easier/harder
	min->SetTolerance(50.);
    min->SetPrintLevel(1);

	const int ParaNbr = 1;

	ROOT::Math::Functor f(&CalibChiSquare, ParaNbr);

    // Step Sizes
	double step[ParaNbr] = {0.001., //  
				};

    // Starting points
	double variable[ParaNbr] = {0.5, // 
				};


    // Sets function to be minimized and variables to be minimized
	min->SetFunction(f);
	min->SetVariable(0, "Calib Norm", variable[0], step[0]);

	min->Minimize();


	cout << "Minimized" << endl;

	const double *Result = min->X();

	Output << min->MinValue() << "    " << Result[0] << " " << Result[1] << " " << Result[2] << " " << Result[3] << "  " << endl;


    Output.close();

    return 0;

}
