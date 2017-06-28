// Run the low energy fitter... do whatever
#include "GPXFitter.hh"
#include "TStyle.h"
#include "TChain.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooAbsArg.h"
#include "RooMsgService.h"
#include "TROOT.h"

using namespace RooFit;

int main(int argc, char** argv)
{
  gROOT->ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);");

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

	GPXFitter *fitter = new GPXFitter(fDS, fitMin, fitMax);
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

	return 0;
}