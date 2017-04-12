// Run the low energy fitter... do whatever
#include "WenqinFitter.hh"

#include "TStyle.h"
#include "TChain.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"

using namespace RooFit;

int main(int argc, char** argv)
{
	if(argc <= 2)
	{
		std::cout << "Usage: " << argv[0] << " [Fit Min] [Fit Max]" << std::endl;
	}
	int fitMin = atoi(argv[1]);
	int fitMax = atoi(argv[2]);

	//// For drawing pretty shit
	gStyle->SetOptStat(0);
	gStyle->SetPalette(kRainBow);

	//// Set cuts here
    std::string theCut = "trapENFCal > 0.8 && gain==0 && mHClean==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin < 0 && D != 0";
    std::string noiseTime = "&& !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s < 7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3)";
    std::string noisyRuns = "&& run!=9648 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12745 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004";
    std::string PSACuts = "&& waveS5/TMath::Power(trapENFCal, 1/4) < 1200 && tOffset < 10 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 320";
    std::string SlowCut = "&& (((kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1) || (channel==580 && ((kvorrT/trapENFCal) >= 0.9 && (kvorrT/trapENFCal) <= 2.1)) || (channel==664 && ((kvorrT/trapENFCal) >= 1.1 && (kvorrT/trapENFCal) <= 2.1)))";
    // std::string SlowCut = "&& ( ((kvorrT/trapENFCal) >= 1.2) && ((kvorrT/trapENFCal) <= 2.1) )"; // DS0 ToE

    theCut += noiseTime;
    theCut += noisyRuns;
    theCut += PSACuts;
    theCut += SlowCut;
    theCut += "&& isEnr && threshkeV > 0.5 && (threshkeV+threshSig) < 2.";

	WenqinFitter *fitter = new WenqinFitter(fitMin, fitMax);
	fitter->SetSavePrefix(Form("DS1_Enr_Test_%d_%d", fitMin, fitMax));
	
	// Load data from TChain with a cut string
	TChain *skimTree = new TChain("skimTree");
	skimTree->Add("~/project/wavelet-skim/waveletSkimDS1_*.root");
	fitter->LoadChainData(skimTree, theCut);

	// Construct PDF and do fit
	fitter->ConstructPDF();
	fitter->DoFit();

	// Output stuff
	fitter->DrawBasicShit();

	//// Both of these take a long time!
	//// Will literally run forever if the parameter is close to 0 or limited
	//// Make Contour before ProfileNLL -- contour gets fked up if profile NLL is done first, don't know why
	// fitter->DrawContour();
	// fitter->DrawContour("Tritium", "Bkg");
	// fitter->DrawContour("Tritium", "Co57");
	// fitter->DrawContour("Tritium", "Mn54");
	// fitter->ProfileNLL();
	// fitter->ProfileNLL("Bkg");
	// fitter->ProfileNLL("Co57");
	// fitter->ProfileNLL("Mn54");

	// Generate Toy MC study
	// fitter->GenerateMCStudy("Tritium", 5000);

	//// Print fit result again at the end
	// RooWorkspace *fitWorkspace = fitter->GetWorkspace();
	// fitWorkspace->Print();

	RooFitResult *fitResult = fitter->GetFitResult();
	std::cout << "Fit Range: " << fitMin <<  " " << fitMax << std::endl;
	fitResult->Print();

	return 0;
}