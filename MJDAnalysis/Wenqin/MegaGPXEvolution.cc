// Run the evolved version of Wenqin

#include "GPXFitter.hh"
#include "TStyle.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooAbsArg.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "TROOT.h"
// #include <Python/Python.h>

using namespace RooFit;
using namespace std;

map<string, vector<double>> RunWenqin(vector<string> argS, int fDS, float fitMin, float fitMax, string theCut);
map<string, vector<double>> RunWenqin(vector<string> argS, int fDS, float fitMin, float fitMax, string ftype, string theCut);

int main(int argc, char** argv)
{
  	gROOT->ProcessLine("gErrorIgnoreLevel = 3001;");
  	gROOT->ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);");
	if(argc <= 4)
	{
		cout << "Usage: " << argv[0] << " [DS] [Fit Min] [Fit Max] [Nat/Enr/Cut]" << endl;
		return 0;
	}

	int fDS = atoi(argv[1]);
  	float fitMin = atof(argv[2]);
	float fitMax = atof(argv[3]);
	string ftype = argv[4];
	gStyle->SetOptStat(0);
	gStyle->SetPalette(kRainBow);

  // Prelim exposure
  std::map<int, std::vector<double>> liveTime;
	liveTime[0] = {460.054, 171.021};
	liveTime[1] = {661.811, 63.2937};
  liveTime[2] = {106.286, 10.6791};
	liveTime[3] = {368.52, 81.7408};
	liveTime[4] = {102.858, 73.8446};
	liveTime[5] = {492.158+182.193, 138.461+197.769};
  liveTime[6] = {460.054+661.811+106.286+368.52+102.858+492.158+182.193, 171.021+63.2937+10.6791+81.7408+73.8446+138.461+197.769};
  int bNat = 0;

	//// Set cuts here
    string theCut = "";
  	theCut += Form("trapENFCal>=%.2f&&trapENFCal<=%.2f", fitMin, fitMax); // Energy cut for fit range

  	if(ftype == "Nat" || ftype == "nat"){
  		theCut += "&&isNat"; // Set Enriched or Natural
  		bNat = 1;
  	}
  	else if(ftype == "Enr" || ftype == "enr") {
  		theCut += "&&isEnr";
  		bNat = 0;
  	}
  	else theCut = ftype;

    std::string tritDir = "/Users/brianzhu/macros/code/MJDAnalysis/Axion";
    TFile *tritFile = new TFile(Form("%s/TritSpec.root", tritDir.c_str()));
    TH1D *tritSpec = dynamic_cast<TH1D*>(tritFile->Get("tritHist"));
    double tritScale = 1./tritSpec->Integral(tritSpec->FindBin(fitMin), tritSpec->FindBin(fitMax));
	vector<string> argList = {"Tritium", "Ge68", "Fe55", "Bkg"};

	map<string, vector<double>> ValMap;
	if(ftype == "Nat" || ftype == "nat" || ftype == "Enr" || ftype == "enr") ValMap = RunWenqin(argList, fDS, fitMin, fitMax, ftype, theCut);
	// else  ValMap = RunWenqin(argList, fDS, fitMin, fitMax, theCut);
	cout << "LiveTime: " << liveTime[fDS][bNat] << endl;

  // Calculate Rates
  for(auto &kv : ValMap)
	{
		if(kv.first == "Tritium") {
			cout << kv.first << ":[" << kv.second[0]*tritScale << "," << kv.second[1]*tritScale << "," << kv.second[2]*tritScale << "]"<< endl;
			cout << "Rates (c/keV/kg/day):" << kv.second[0]*tritScale/liveTime[fDS][bNat]/(fitMax-fitMin) << " +" << kv.second[1]*tritScale/liveTime[fDS][bNat]/(fitMax-fitMin) << " " << kv.second[2]*tritScale/liveTime[fDS][bNat]/(fitMax-fitMin) << endl;
			cout << "Rates (2-4 keV) (c/keV/kg/day):" << kv.second[0]*tritScale*0.224487/liveTime[fDS][bNat]/(fitMax-fitMin) << " +" << kv.second[1]*tritScale*0.224487/liveTime[fDS][bNat]/(fitMax-fitMin) << " " << kv.second[2]*tritScale*0.224487/liveTime[fDS][bNat]/(fitMax-fitMin) << endl;
		}
		else {
			cout << kv.first << ":[" << kv.second[0] << "," << kv.second[1] << "," << kv.second[2] << "]"<< endl;
			cout << "Rates (c/keV/kg/day):" << kv.second[0]/liveTime[fDS][bNat]/(fitMax-fitMin) << " +" << kv.second[1]/liveTime[fDS][bNat]/(fitMax-fitMin) << " " << kv.second[2]/liveTime[fDS][bNat]/(fitMax-fitMin) << endl;
		}
	}

	return 0;
}

map<string, vector<double>> RunWenqin(vector<string> argS, int fDS, float fitMin, float fitMax, string theCut, int idx)
{
	string ftype = Form("Ch%d",idx);
	map<string, vector<double>> WenqinMap = RunWenqin(argS, fDS, fitMin, fitMax, ftype, theCut);

	return WenqinMap;
}

map<string, vector<double>> RunWenqin(vector<string> argS, int fDS, float fitMin, float fitMax, string ftype, string theCut)
{
	GPXFitter *fitter = new GPXFitter(fDS, fitMin, fitMax);
	// This is just a string for output files
	fitter->SetSavePrefix(Form("DS%d_%s_%.1f_%.1f", fDS, ftype.c_str(), fitMin, fitMax));

	// Load data from TChain with a cut string
	TChain *skimTree = new TChain("skimTree");
    if(fDS == 6) {
      skimTree->Add("/Users/brianzhu/project/cuts/corrfs_rn/corrfs_rn-DS1-*.root" );
      skimTree->Add("/Users/brianzhu/project/cuts/corrfs_rn/corrfs_rn-DS2-*.root" );
      skimTree->Add("/Users/brianzhu/project/cuts/corrfs_rn/corrfs_rn-DS3-*.root" );
      skimTree->Add("/Users/brianzhu/project/cuts/corrfs_rn/corrfs_rn-DS4-*.root" );
      skimTree->Add("/Users/brianzhu/project/cuts/corrfs_rn/corrfs_rn-DS5-*.root" );
    }
  	else skimTree->Add(Form("/Users/brianzhu/project/cuts/corrfs_rn/corrfs_rn-DS%d-*.root", fDS) );
  	fitter->LoadChainData(skimTree, theCut);

	// Construct PDF and do fit
	fitter->ConstructPDF();
  fitter->DoFit();
	fitter->DrawBasicShit(0.2, false, false);
	fitter->GetFitResult()->Print("v");
	map<string, vector<double>> WenqinMap;
	for(auto &argN: argS){
		vector<double> WenqinParameter = fitter->GetVar( Form("%s", argN.c_str()) );
		WenqinMap[argN.c_str()] = WenqinParameter;
	}
	return WenqinMap;
}
