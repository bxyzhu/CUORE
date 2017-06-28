// Run the low energy fitter... do whatever
#include "WenqinFitter.hh"
#include "TStyle.h"
#include "TChain.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooAbsArg.h"

using namespace RooFit;
/*
  Old Results For Comparison (8 Mar 2016)
	Note: DS-0 might be lower by about 10% due to 60 hours of data being incorrect
	DS0 Enriched:				Counts 							             Rate (Counts/kg/year)
                   Bkg    1.2645e+03 +/-  3.75e+01			2.48614 +/- 0.0737288
                  Co57    4.4537e-06 +/-  6.59e+00			8.756e-09 +/- 0.01295
                  Fe55    2.5077e-06 +/-  6.22e+00			4.930e-09 +/- 0.01223
                  Ge68    1.1885e+01 +/-  5.47e+00			0.023367 +/- 0.010755
                  Mn54    7.4401e+00 +/-  6.21e+00			0.014628 +/- 0.012209
                 Pb210    1.2458e+01 +/-  4.89e+00			0.024494 +/- 0.009614
               Tritium    1.5544e+02 +/-  1.92e+01			0.305532 +/- 0.037749
                  Zn65    2.2986e+00 +/-  4.72e+00			0.004519 +/- 0.009280

	DS1 Enriched:				Counts 							             Rate (Counts/kg/year)
                   Bkg    6.2617e+02 +/-  2.66e+01			0.920920 +/- 0.0391211
                  Co57    3.5904e+00 +/-  5.28e+00			0.005280 +/- 0.0077654
                  Fe55    1.6264e+01 +/-  6.70e+00			0.023919 +/- 0.0098538
                  Ge68    7.5651e+00 +/-  4.80e+00			0.011126 +/- 0.0070594
                  Mn54    4.6256e-05 +/-  3.14e+01			6.803e-08 +/- 0.046181
                 Pb210    9.2066e+00 +/-  3.80e+00			0.013540 +/- 0.0055887
               Tritium    2.0133e+02 +/-  2.01e+01			0.296071 +/- 0.0295615
                  Zn65    5.8570e+00 +/-  4.98e+00			0.008614 +/- 0.0073242
	
	DS0 Natural:				Counts 							             Rate (Counts/kg/year)
                   Bkg    6.6178e+02 +/-  2.75e+01			3.55730 +/- 0.147822
                  Co57    8.8618e-08 +/-  7.37e+00			4.76e-10 +/- 0.0396164
                  Fe55    1.1705e+02 +/-  1.52e+01			0.62919 +/- 0.081705
                  Ge68    2.1266e+02 +/-  1.68e+01			1.14312 +/- 0.090306
                  Mn54    1.3809e+01 +/-  1.11e+01			0.07423 +/- 0.059666
                 Pb210    4.8461e-01 +/-  2.79e+00			0.00260 +/- 0.014997
               Tritium    8.9171e+02 +/-  3.84e+01			4.79320 +/- 0.206414
                  Zn65    4.3513e+01 +/-  1.08e+01			0.23389 +/- 0.058054

	DS1 Natural:				Counts 							             Rate (Counts/kg/year)
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
	if(argc <= 4)
	{
		std::cout << "Usage: " << argv[0] << " [DS] [Fit Min] [Fit Max] [Nat/Enr]" << std::endl;
		return 0;
	}
	int fDS = atoi(argv[1]);
  float fitMin = atof(argv[2]);
	float fitMax = atof(argv[3]);
	std::string ftype = argv[4];
  //// For drawing pretty shit
	gStyle->SetOptStat(0);
	gStyle->SetPalette(kRainBow);

	//// Set cuts here
  // std::string theCut = "gain==0 && mH==1 && isGood && !muVeto && !wfDCBits && !isLNFill1 && !isLNFill2 && trapETailMin<0. && P!=0 && D!=0";
  // std::string DS0runCut = "&& run!=3341 && run!=3523 && run!=3524 && run!=3529 && run!=3530 && run!=4057 && run!=4125 && run!=4473 && run!=4554 && run!=5248 && run!=5249 && run!=5889 && run!=5890 && run!=5894 && run!=5895 && run!=6775";
  // Noisy cuts for DS1
  // std::string noiseTime = "&& !(time_s > 2192e3 && time_s < 2195e3) && !(time_s > 7370e3 && time_s < 7371e3) && !(time_s > 7840e3 && time_s < 7860e3) && !(time_s > 8384e3 && time_s < 8387e3) && !(time_s > 8984e3 && time_s < 8985e3) && !(time_s > 9002e3 && time_s < 9005e3)";
  // std::string noisyRuns = "&& run!=9648 && run!=10663 && run!=10745 && run!=11175 && run!=12445 && run!=12723 && run!=12735 && run!=12745 && run!=12746 && run!=12765 && run!=12766 && run!=12767 && run!=13004";
  // std::string PSACuts = "&& bcMax < 599";
  // ToE for DS1
  // std::string SlowCut = "&& (((kvorrT/trapENFCal) >= 1.2 && (kvorrT/trapENFCal) <= 2.1) || (channel==580 && ((kvorrT/trapENFCal) >= 0.9 && (kvorrT/trapENFCal) <= 2.1)) || (channel==664 && ((kvorrT/trapENFCal) >= 1.1 && (kvorrT/trapENFCal) <= 2.1)))";
  // std::string ThreshCut = "&& threshKeV > 0 && threshKeV < 1.5";
  // ToE for DS0
  // std::string SlowCut = "&& ( ((kvorrT/trapENFCal) >= 1.2) && ((kvorrT/trapENFCal) <= 2.1) )"; // DS0 ToE

  std::string theCut = "";
  if(fDS == 0) theCut += "!(run==3523 && (channel==690)) && channel!=656";
  
  // else if(fDS == 1) theCut += "!(run==10663 && (channel==578)) && !(run==11175 && (channel==578)) && !(run==12746 && (channel==692)) && !(run==12765 && (channel==692)) && !(run==12766 && (channel==692)) && !(run==13004 && (channel==598|| channel==600))";
  else if (fDS==1) theCut += "!(run==10663&&(channel==578||channel==610||channel==692))&&!(run==10745&&(channel==578||channel==664||channel==692))&&!(run==11175&&(channel==578))&&!(run==12735&&(channel==580||channel==598||channel==600||channel==664||channel==672||channel==692))&&!(run==12745&&(channel==580||channel==582||channel==664||channel==692))&&!(run==12746&&(channel==692))&&!(run==12765&&(channel==608||channel==692))&&!(run==12766&&(channel==598||channel==600||channel==648||channel==690||channel==692))&&!(run==12767&&(channel==600||channel==692))&&!(run==13004&&(channel==580||channel==598||channel==600||channel==640||channel==648||channel==664||channel==672)) && bcMax < 400 && fitSlo < 15";

  else if(fDS == 3) theCut += "!(run==17169 && (channel==614)) && !(run==17178 && (channel==614)) && !(run==17505 && (channel==598)) && !(run==17875 && (channel==614)) && !(run==17878 && (channel==614)) && !(run==17935 && (channel==614)) && !(run==17936 && (channel==614)) && !(run==17940 && (channel==614)) && !(run==17971 && (channel==614)) && !(run==17978 && (channel==614)) && channel!=592";
  
  else if(fDS == 4) theCut += "!(run==60001001 && (channel==1106)) && !(run==60001002 && (channel==1106)) && !(run==60001003 && (channel==1106)) && !(run==60001011 && (channel==1106)) && !(run==60001078 && (channel==1144)) && !(run==60001121 && (channel==1144)) && !(run==60001171 && (channel==1106)) && !(run==60001177 && (channel==1144)) && !(run==60001178 && (channel==1106|| channel==1144)) && !(run==60001184 && (channel==1144)) && !(run==60001463 && (channel==1106)) && !(run==60001469 && (channel==1106)) && !(run==60001471 && (channel==1106)) && !(run==60001477 && (channel==1106)) && !(run==60001492 && (channel==1106)) && !(run==60001498 && (channel==1106|| channel==1298)) && !(run==60001692 && (channel==1144)) && !(run==60001698 && (channel==1106)) && !(run==60001797 && (channel==1144)) && !(run==60001827 && (channel==1144)) && !(run==60001838 && (channel==1106)) && !(run==60001847 && (channel==1106)) && !(run==60001848 && (channel==1106)) && !(run==60001881 && (channel==1106)) && !(run==60001882 && (channel==1106)) && channel!=1332";

  else if(fDS == 5) theCut += "!(run==22326&&(channel==614))&&!(run==22328 && (channel==1124))&&!(run==22406&&(channel==614))&&!(run==22428&&(channel==1124))&&!(run==22444&&(channel==614))&&!(run==22461 && (channel==1124)) && !(run==22462 && (channel==1124)) && !(run==22480 && (channel==614)) && !(run==22493 && (channel==1124)) && !(run==22494 && (channel==1124)) && !(run==22666 && (channel==614)) && !(run==22667 && (channel==614)) && !(run==22939 && (channel==614)) && !(run==23533 && (channel==614)) && !(run==23705 && (channel==1298|| channel==1302)) && !(run==23706 && (channel==1302)) && !(run==23707 && (channel==1302)) && !(run==23708 && (channel==1298|| channel==1302)) && !(run==23709 && (channel==1298|| channel==1302)) && !(run==23710 && (channel==1298|| channel==1302)) && !(run==23711 && (channel==1298|| channel==1302)) && !(run==23712 && (channel==1298|| channel==1302)) && !(run==23713 && (channel==1298|| channel==1302)) && !(run==23714 && (channel==1298|| channel==1302)) && !(run==23715 && (channel==1298|| channel==1302)) && !(run==23718 && (channel==1302)) && !(run==23719 && (channel==1298|| channel==1302)) && !(run==23721 && (channel==1298|| channel==1302)) && !(run==23725 && (channel==1330|| channel==1332)) && !(run==23726 && (channel==1330|| channel==1332)) && !(run==23728 && (channel==1330|| channel==1332)) && !(run==23729 && (channel==1330|| channel==1332)) && !(run==23730 && (channel==1330|| channel==1332)) && !(run==23731 && (channel==1330|| channel==1332)) && !(run==23732 && (channel==1330|| channel==1332)) && !(run==23733 && (channel==1330|| channel==1332)) && !(run==23734 && (channel==1330|| channel==1332)) && !(run==23735 && (channel==1330|| channel==1332)) && !(run==23736 && (channel==1330|| channel==1332)) && !(run==23737 && (channel==1330|| channel==1332)) && !(run==23738 && (channel==1330|| channel==1332)) && !(run==23739 && (channel==1330|| channel==1332)) && !(run==23740 && (channel==1330|| channel==1332)) && !(run==23741 && (channel==1330|| channel==1332)) && !(run==23742 && (channel==1330|| channel==1332)) && !(run==23743 && (channel==1330|| channel==1332)) && !(run==23744 && (channel==1330|| channel==1332)) && !(run==23745 && (channel==1330|| channel==1332)) && !(run==23746 && (channel==1330|| channel==1332)) && !(run==23747 && (channel==1330|| channel==1332)) && !(run==23748 && (channel==1330|| channel==1332)) && !(run==23749 && (channel==1330|| channel==1332)) && !(run==23750 && (channel==1330|| channel==1332)) && !(run==23751 && (channel==1330|| channel==1332)) && !(run==23752 && (channel==1330|| channel==1332)) && !(run==23753 && (channel==1330|| channel==1332)) && !(run==23754 && (channel==1330|| channel==1332)) && !(run==23755 && (channel==1330|| channel==1332)) && !(run==23756 && (channel==1330|| channel==1332)) && !(run==23757 && (channel==1332)) && !(run==23758 && (channel==1330|| channel==1332)) && !(run==23759 && (channel==1330|| channel==1332))";

  if(ftype == "Nat") theCut += "&&isNat"; // Set Enriched or Natural
  else if(ftype == "Enr") theCut += "&&isEnr";

  theCut += Form("&& trapENFCal>=%.2f && trapENFCal<=%.2f", fitMin, fitMax); // Energy cut for fit range

	WenqinFitter *fitter = new WenqinFitter(fDS, fitMin, fitMax);
	// This is just a string for output files
	fitter->SetSavePrefix(Form("LowE_DS%d_%s_%.1f_%.1f", fDS, ftype.c_str(), fitMin, fitMax));

	// Load data from TChain with a cut string
	TChain *skimTree = new TChain("skimTree");
  skimTree->Add(Form("~/project/latskim/latSkimDS%d_*.root", fDS) );
  fitter->LoadChainData(skimTree, theCut);

	// Construct PDF and do fit
	fitter->ConstructPDF();
	fitter->DoFit();

	// This draws the spectrum as well as the covariance matrix and residuals if you want
	fitter->DrawBasicShit(0.2, false, false);

  // List of parameters we want to do more studies on
  // std::vector<std::string> argList = {"Tritium", "Ge68"};
	//// Profile likelihood calculation
	// std::map<std::string, std::vector<double>> LimitMap = fitter->ProfileNLL(argList);
	
  // Generate Toy MC study -- input a list of parameters to look at from toy MC
  // This shit takes a long time
  // fitter->GenerateMCStudy(argList, 1000);

  // Print limits
  // for(auto &kv : LimitMap) std::cout << kv.first << " limit: " << kv.second[0] << " " << kv.second[1] << std::endl;

	//// Print fit result again at the end
	// RooWorkspace *fitWorkspace = fitter->GetWorkspace();
	// fitWorkspace->Print();
	// RooFitResult *fitResult = fitter->GetFitResult();
	// std::cout << "Fit Range: " << fitMin <<  " " << fitMax << std::endl;
	// fitResult->Print();
	// RooArgList floatPars = fitResult->floatParsFinal();
	// floatPars.Print();

	return 0;
}