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

  DS-1:
  Total live time (days): 60.1582
  Reduction from veto (sec): 3435.08
  Reduced live time (days): 60.1184
  Active mass:  12.4300 kg   enriched 11.3100   natural 1.1200
  Exposure:  747.2720 kg-days total  natural 67.3326   enriched 679.9394

  DS-3:
  Total live time (days): 29.9128
  Reduction from veto (sec): 1117.99
  Reduced live time (days): 29.8999
  Active mass:  15.4120 kg   enriched 12.6310   natural 2.7810
  Exposure:  460.8174 kg-days total  natural 83.1516   enriched 377.6657
  
  DS-4:
  Total live time (days): 23.6908
  Reduction from veto (sec): 9986.62
  Reduced live time (days): 23.5752
  Active mass:  9.4212 kg   enriched 5.4712   natural 3.9500
  Exposure:  222.1064 kg-days total  natural 93.1219   enriched 128.9845


	Note: DS-0 might be lower by about 10% due to 60 hours of data being incorrect
	DS0 Enriched:				Counts 							Rate (Counts/kg/year)
                   Bkg    1.2645e+03 +/-  3.75e+01			2.48614 +/- 0.0737288
                  Co57    4.4537e-06 +/-  6.59e+00			8.756e-09 +/- 0.01295
                  Fe55    2.5077e-06 +/-  6.22e+00			4.930e-09 +/- 0.01223
                  Ge68    1.1885e+01 +/-  5.47e+00			0.023367 +/- 0.010755
                  Mn54    7.4401e+00 +/-  6.21e+00			0.014628 +/- 0.012209
                 Pb210    1.2458e+01 +/-  4.89e+00			0.024494 +/- 0.009614
               Tritium    1.5544e+02 +/-  1.92e+01			0.305532 +/- 0.037749
                  Zn65    2.2986e+00 +/-  4.72e+00			0.004519 +/- 0.009280

	DS1 Enriched:				Counts 							Rate (Counts/kg/year)
                   Bkg    6.2617e+02 +/-  2.66e+01			0.920920 +/- 0.0391211
                  Co57    3.5904e+00 +/-  5.28e+00			0.005280 +/- 0.0077654
                  Fe55    1.6264e+01 +/-  6.70e+00			0.023919 +/- 0.0098538
                  Ge68    7.5651e+00 +/-  4.80e+00			0.011126 +/- 0.0070594
                  Mn54    4.6256e-05 +/-  3.14e+01			6.803e-08 +/- 0.046181
                 Pb210    9.2066e+00 +/-  3.80e+00			0.013540 +/- 0.0055887
               Tritium    2.0133e+02 +/-  2.01e+01			0.296071 +/- 0.0295615
                  Zn65    5.8570e+00 +/-  4.98e+00			0.008614 +/- 0.0073242
	
	DS0 Natural:				Counts 							Rate (Counts/kg/year)
                   Bkg    6.6178e+02 +/-  2.75e+01			3.55730 +/- 0.147822
                  Co57    8.8618e-08 +/-  7.37e+00			4.76e-10 +/- 0.0396164
                  Fe55    1.1705e+02 +/-  1.52e+01			0.62919 +/- 0.081705
                  Ge68    2.1266e+02 +/-  1.68e+01			1.14312 +/- 0.090306
                  Mn54    1.3809e+01 +/-  1.11e+01			0.07423 +/- 0.059666
                 Pb210    4.8461e-01 +/-  2.79e+00			0.00260 +/- 0.014997
               Tritium    8.9171e+02 +/-  3.84e+01			4.79320 +/- 0.206414
                  Zn65    4.3513e+01 +/-  1.08e+01			0.23389 +/- 0.058054

	DS1 Natural:				Counts 							Rate (Counts/kg/year)
                   Bkg    1.6964e+02 +/-  1.41e+01			2.51943 +/- 0.209408
                  Co57    1.0509e+00 +/-  1.31e+01			0.01561 +/- 0.194557
                  Fe55    6.1692e+01 +/-  1.12e+01			0.91623 +/- 0.166338
                  Ge68    1.0079e+02 +/-  1.17e+01			1.49690 +/- 0.173764
                  Mn54    6.5398e+00 +/-  8.04e+00			0.09713 +/- 0.119407
                 Pb210    4.8488e+00 +/-  2.61e+00			0.07201 +/- 0.038763
               Tritium    4.7636e+02 +/-  2.94e+01			7.07473 +/- 0.436638
                  Zn65    1.4089e+01 +/-  7.27e+00			0.20925 +/- 0.107971
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
    
    // Noisy cuts for DS1
    std::string noiseTime = "&& !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s < 7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3)";
    std::string noisyRuns = "&& run!=9648 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12745 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004";
    std::string PSACuts = "&& waveS5/TMath::Power(trapENFCal, 1/4) < 1200 && tOffset < 10 && (waveS3-waveS2)/trapENFCal > 100 && (waveS3-waveS2)/trapENFCal < 320";
    // ToE for DS1
    // std::string SlowCut = "&& (((kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1) || (channel==580 && ((kvorrT/trapENFCal) >= 0.9 && (kvorrT/trapENFCal) <= 2.1)) || (channel==664 && ((kvorrT/trapENFCal) >= 1.1 && (kvorrT/trapENFCal) <= 2.1)))";
    // ToE for DS0
    // std::string SlowCut = "&& ( ((kvorrT/trapENFCal) >= 1.2) && ((kvorrT/trapENFCal) <= 2.1) )"; // DS0 ToE

    //// Add all the cuts to "theCut"
    // theCut += noiseTime;
    // theCut += noisyRuns;
    theCut += PSACuts;
    // theCut += SlowCut;
    theCut += "&& isEnr"; // Set Enriched or Natural

	WenqinFitter *fitter = new WenqinFitter(fitMin, fitMax);
	// This is just a string for output files
	fitter->SetSavePrefix(Form("DS0_Enr_WithPb_%d_%d", fitMin, fitMax));
	
	// Load data from TChain with a cut string
	TChain *skimTree = new TChain("skimTree");
	skimTree->Add("~/project/wavelet-skim/waveletSkimDS0_*.root");
	fitter->LoadChainData(skimTree, theCut);

	// Construct PDF and do fit
	fitter->ConstructPDF();
	fitter->DoFit();

	// Output stuff -- argument sets the bin size for the output files
	// This draws the spectrum as well as the covariance matrix
	fitter->DrawBasicShit();
	// fitter->DrawBasicShit(1.0);

	//// Both of these take a long time!
	// The contour plot is kinda fked up, 
	// fitter->DrawContour();
	// fitter->DrawContour("Tritium", "Bkg");
	// fitter->DrawContour("Tritium", "Co57");
	// fitter->DrawContour("Tritium", "Mn54");
	// fitter->ProfileNLL();
	// fitter->ProfileNLL("Ge68");
	// fitter->ProfileNLL("Co57");
	// fitter->ProfileNLL("Mn54");
	// fitter->ProfileNLL("Zn65");
	// fitter->ProfileNLL("Pb210");
	// fitter->ProfileNLL("Bkg");
	// fitter->ProfileNLL("Fe55");

	// Generate Toy MC study -- this shit takes a long time
	// fitter->GenerateMCStudy("Ge68", 5000);

	//// Print fit result again at the end
	// RooWorkspace *fitWorkspace = fitter->GetWorkspace();
	// fitWorkspace->Print();
	RooFitResult *fitResult = fitter->GetFitResult();
	// std::cout << "Fit Range: " << fitMin <<  " " << fitMax << std::endl;
	fitResult->Print();
	// std::cout << "Chi Square: " << fitter->fChiSquare << std::endl;

	return 0;
}