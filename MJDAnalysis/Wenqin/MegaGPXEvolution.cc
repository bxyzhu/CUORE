// Run the evolved version of Wenqin

#include "GPXFitter.hh"
#include "TStyle.h"
#include "TChain.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooAbsArg.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "TROOT.h"

using namespace RooFit;

// std::vector<double> RunWenqin(std::string argN);

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

  	if(ftype == "Nat" || ftype == "nat") theCut += "isNat"; // Set Enriched or Natural
  	else if(ftype == "Enr" || ftype == "enr") theCut += "isEnr";

  	if(fDS == 0) theCut += "&&!(run==6811&&(channel==600||channel==696)) && channel!=656";
  	else if(fDS == 3) theCut += "&&channel!=592 && channel!=692";
  	else if(fDS == 4) theCut += "&&!(run==60001692 && (channel==1144))&&channel!=1332";
  	else if(fDS == 5) theCut += "&&channel!=1124";
  	else if(fDS == 6) theCut += "&&!(run==6811&&(channel==600||channel==696)) && !(run==60001692 && (channel==1144))";

  	theCut += Form("&& trapENFCal>=%.2f && trapENFCal<=%.2f", fitMin, fitMax); // Energy cut for fit range

	GPXFitter *fitter = new GPXFitter(fDS, fitMin, fitMax);
	// This is just a string for output files
	fitter->SetSavePrefix(Form("LowE_DS%d_%s_%.1f_%.1f", fDS, ftype.c_str(), fitMin, fitMax));

	// Load data from TChain with a cut string
	TChain *skimTree = new TChain("skimTree");
  	if(fDS == 6) skimTree->Add("~/project/latskim/latSkimDS*.root" );
  	else skimTree->Add(Form("~/project/latskim/latSkimDS%d_*.root", fDS) );
  	fitter->LoadChainData(skimTree, theCut);

	// Construct PDF and do fit
	fitter->ConstructPDF();
	fitter->DoFit();

	// This draws the spectrum as well as the covariance matrix and residuals if you want
	// fitter->DrawBasicShit(0.2, false, false);

  // List of parameters we want to do more studies on
  // std::vector<std::string> argList = {"Tritium", "Ge68"};
	//// Profile likelihood calculation
	// std::map<std::string, std::vector<double>> LimitMap = fitter->ProfileNLL(argList);
	
  // Generate Toy MC study -- input a list of parameters to look at from toy MC
  // This shit takes a long time
  // fitter->GenerateMCStudy(argList, 1000);

  // Print limits
  // for(auto &kv : LimitMap) std::cout << kv.first << " limit: " << kv.second[0] << " " << kv.second[1] << std::endl;

	RooFitResult *fitResult = fitter->GetFitResult();
	fitResult->Print();

	// RooWorkspace *fitWorkspace = fitter->GetWorkspace();
	// RooRealVar *Ge68 = dynamic_cast<RooRealVar*>(fitResult->floatParsFinal().find("Ge68"));
	// double geVal = dynamic_cast<RooRealVar*>(fitResult->floatParsFinal().find("Ge68"))->getValV();
	std::vector<double> GeVals = fitter->GetVar("Ge68");
	std::cout << "Ge68: " << GeVals[0] << " " <<  GeVals[1] << " " << GeVals[2] << std::endl;

	std::vector<double> FeVals = fitter->GetVar("Fe55");
	std::cout << "Fe55: " << FeVals[0] << " " <<  FeVals[1] << " " << FeVals[2] << std::endl;

	return 0;
}