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
#include <Python/Python.h>

using namespace RooFit;
using namespace std;

map<string, vector<double>> RunWenqin(vector<string> argS, int fDS, float fitMin, float fitMax, string theCut, int idx);
map<string, vector<double>> RunWenqin(vector<string> argS, int fDS, float fitMin, float fitMax, string ftype, string theCut);

int main(int argc, char** argv)
{
  	gROOT->ProcessLine("gErrorIgnoreLevel = 3001;");
  	gROOT->ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);");
	if(argc <= 4)
	{
		cout << "Usage: " << argv[0] << " [DS] [Fit Min] [Fit Max] [Nat/Enr/Cut] [idx]" << endl;
		return 0;
	}
	
	// Py_SetProgramName(argv[0]);  /* optional but recommended */
 // 	Py_Initialize();
 //  	PyRun_SimpleString("from time import time,ctime\n"
 //                     "print 'Today is',ctime(time())\n");
 //  	Py_Finalize();

	int fDS = atoi(argv[1]);
  	float fitMin = atof(argv[2]);
	float fitMax = atof(argv[3]);
	string ftype = argv[4];
	int fidx = atoi(argv[5]);
  	//// For drawing pretty shit
	gStyle->SetOptStat(0);
	gStyle->SetPalette(kRainBow);

	std::map<int, std::vector<double>> liveTime; 
	// liveTime[0] = {449.241, 185.208};
	// liveTime[1] = {661.597, 66.0657};
	// liveTime[3] = {331.405, 66.4533};
	// liveTime[4] = {105.232, 93.5911};
	liveTime[0] = {332.269, 64.8862};
	liveTime[1] = {646.728, 65.6867};
	liveTime[3] = {327.951, 44.8998};
	liveTime[4] = {89.6705, 67.3392};
	liveTime[5] = {662.177+319.476, 154.725+189.524};
	int bNat = 0;

	//// Set cuts here
  	string theCut = "";
  	theCut += Form("trapENFCal>=%.2f&&trapENFCal<=%.2f", fitMin, fitMax); // Energy cut for fit range
  	if(fDS == 0) theCut += "&&!(run==6811&&(channel==600||channel==696)) && channel!=656";
  	else if(fDS == 3) theCut += "&&channel!=592 && channel!=692";
  	else if(fDS == 4) theCut += "&&!(run==60001692 && (channel==1144))&&channel!=1332";
  	else if(fDS == 5) theCut += "&&channel!=1124";
  	else if(fDS == 6) theCut += "&&!(run==6811&&(channel==600||channel==696)) && !(run==60001692 && (channel==1144))";

  	if(ftype == "Nat" || ftype == "nat"){ 
  		theCut += "&&isNat"; // Set Enriched or Natural
  		bNat = 1;
  	}
  	else if(ftype == "Enr" || ftype == "enr") {
  		theCut += "&&isEnr && fitSlo > 0.4";
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
	else  ValMap = RunWenqin(argList, fDS, fitMin, fitMax, theCut, fidx);
	cout << "LiveTime: " << liveTime[fDS][bNat] << endl;
	for(auto &kv : ValMap)
	{
		if(kv.first == "Tritium") {
			cout << kv.first << fidx << ":[" << kv.second[0]*tritScale << "," << kv.second[1]*tritScale << "," << kv.second[2]*tritScale << "]"<< endl;
			cout << "Rates (c/kg/day):" << kv.second[0]*tritScale/liveTime[fDS][bNat] << " +" << kv.second[1]*tritScale/liveTime[fDS][bNat] << " " << kv.second[2]*tritScale/liveTime[fDS][bNat] << endl;
			cout << "Rates (2-4 keV) (c/kg/day):" << kv.second[0]*tritScale*0.224487/liveTime[fDS][bNat] << " +" << kv.second[1]*tritScale*0.224487/liveTime[fDS][bNat] << " " << kv.second[2]*tritScale*0.224487/liveTime[fDS][bNat] << endl;
			cout << "Rates (c/keV/kg/day):" << kv.second[0]*tritScale/liveTime[fDS][bNat]/(fitMax-fitMin) << " +" << kv.second[1]*tritScale/liveTime[fDS][bNat]/(fitMax-fitMin) << " " << kv.second[2]*tritScale/liveTime[fDS][bNat]/(fitMax-fitMin) << endl;
			cout << "Rates (2-4 keV) (c/keV/kg/day):" << kv.second[0]*tritScale*0.224487/liveTime[fDS][bNat]/(fitMax-fitMin) << " +" << kv.second[1]*tritScale*0.224487/liveTime[fDS][bNat]/(fitMax-fitMin) << " " << kv.second[2]*tritScale*0.224487/liveTime[fDS][bNat]/(fitMax-fitMin) << endl;		
		}
		else {
			cout << kv.first << fidx << ":[" << kv.second[0] << "," << kv.second[1] << "," << kv.second[2] << "]"<< endl;
			cout << "Rates (c/kg/day):" << kv.second[0]/liveTime[fDS][bNat] << " +" << kv.second[1]/liveTime[fDS][bNat] << " " << kv.second[2]/liveTime[fDS][bNat] << endl;
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
	fitter->SetSavePrefix(Form("FloatMean_DS%d_%s_%.1f_%.1f", fDS, ftype.c_str(), fitMin, fitMax));

	// Load data from TChain with a cut string
	TChain *skimTree = new TChain("skimTree");
    if(fDS == 6) skimTree->Add("~/project/lat/latSkimDS*.root" );
  	else skimTree->Add(Form("~/project/lat/latSkimDS%d_*.root", fDS) );
  	// if(fDS == 6){ 
  		// skimTree->Add("~/project/latskimv4/latSkimDS0*.root" );
  		// skimTree->Add("~/project/latskimv4/latSkimDS1*.root" );
  		// skimTree->Add("~/project/latskimv4/latSkimDS3*.root" );
  		// skimTree->Add("~/project/latskimv4/latSkimDS4*.root" );
  	// }
  	// else skimTree->Add(Form("~/project/latskimv4/latSkimDS%d_*.root", fDS) );


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
