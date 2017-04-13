// Run the low energy fitter... do whatever
#include "WenqinFitter.hh"

#include "TStyle.h"
#include "TChain.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"

using namespace RooFit;

/*
  Results (8 Mar 2016)
  DS-0:
  Total live time (days): 47.6217
  Reduction from veto (sec): 3683.46
  Reduced live time (days): 47.5791
  Active mass:  14.6000 kg   enriched 10.6900   natural 3.9100
  Exposure:  694.6549 kg-days total  natural 186.0343   enriched 508.6206
  Tritium (counts):  natural: 8.9170e+02 +/-  3.84e+01  enriched: 1.5540e+02 +/- 1.92e+01

  	(Lower by about 10%)
  	Rate: natural: 4.79320 +/- 0.206414    enriched: 0.305532 +/- 0.0377492

  DS-1:
  Total live time (days): 60.1582
  Reduction from veto (sec): 3435.08
  Reduced live time (days): 60.1184
  Active mass:  12.4300 kg   enriched 11.3100   natural 1.1200
  Exposure:  747.2720 kg-days total  natural 67.3326   enriched 679.9394
  Tritium (counts):  natural: 4.7636e+02 +/- 2.94e+01  enriched: 2.0131e+02 +/- 2.01e+01

	Rate: natural: 7.07473 +/- 0.436638		enriched: 0.296071 +/- 0.0295615
*/

int main(int argc, char** argv)
{
	if(argc <= 2)
	{
		std::cout << "Usage: " << argv[0] << " [Fit Min] [Fit Max]" << std::endl;
		return 0;
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
    theCut += "&& isEnr";

	WenqinFitter *fitter = new WenqinFitter(fitMin, fitMax);
	fitter->SetSavePrefix(Form("DS1_Enr_WithPb_%d_%d", fitMin, fitMax));
	
	// Load data from TChain with a cut string
	TChain *skimTree = new TChain("skimTree");
	skimTree->Add("~/project/wavelet-skim/waveletSkimDS1_*.root");
	fitter->LoadChainData(skimTree, theCut);

	// Construct PDF and do fit
	fitter->ConstructPDF();
	fitter->DoFit();

	// Output stuff
	// fitter->DrawBasicShit();
	// fitter->DrawBasicShit(1.0);

	//// Both of these take a long time!
	//// Will literally run forever if the parameter is close to 0 or limited
	//// Make Contour before ProfileNLL -- contour gets fked up if profile NLL is done first, don't know why
	// fitter->DrawContour();
	// fitter->DrawContour("Tritium", "Bkg");
	// fitter->DrawContour("Tritium", "Co57");
	// fitter->DrawContour("Tritium", "Mn54");
	// fitter->ProfileNLL();
	fitter->ProfileNLL("Ge68");
	// fitter->ProfileNLL("Co57");
	// fitter->ProfileNLL("Mn54");

	// Generate Toy MC study
	fitter->GenerateMCStudy("Ge68", 5000);

	//// Print fit result again at the end
	// RooWorkspace *fitWorkspace = fitter->GetWorkspace();
	// fitWorkspace->Print();

	// RooFitResult *fitResult = fitter->GetFitResult();
	// std::cout << "Fit Range: " << fitMin <<  " " << fitMax << std::endl;
	// fitResult->Print();
	// std::cout << "Chi Square: " << fitter->fChiSquare << std::endl;

	return 0;
}