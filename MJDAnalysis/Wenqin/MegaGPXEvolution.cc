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
#include <Python/Python.h>

using namespace RooFit;
using namespace std;

map<string, vector<double>> RunWenqin(vector<string> argS, int fDS, float fitMin, float fitMax, string theCut, int idx);
map<string, vector<double>> RunWenqin(vector<string> argS, int fDS, float fitMin, float fitMax, string ftype, string theCut);

int main(int argc, char** argv)
{
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

	//// Set cuts here
  	string theCut = "";
  	theCut += Form("trapENFCal>=%.2f&&trapENFCal<=%.2f", fitMin, fitMax); // Energy cut for fit range
  	if(fDS == 0) theCut += "&&!(run==6811&&(channel==600||channel==696)) && channel!=656";
  	else if(fDS == 3) theCut += "&&channel!=592 && channel!=692";
  	else if(fDS == 4) theCut += "&&!(run==60001692 && (channel==1144))&&channel!=1332";
  	else if(fDS == 5) theCut += "&&channel!=1124";
  	else if(fDS == 6) theCut += "&&!(run==6811&&(channel==600||channel==696)) && !(run==60001692 && (channel==1144))";

  	if(ftype == "Nat" || ftype == "nat") theCut += "&&isNat"; // Set Enriched or Natural
  	else if(ftype == "Enr" || ftype == "enr") theCut += "&&isEnr";
  	else theCut = ftype;

	vector<string> argList = {"Tritium", "Ge68", "Fe55", "Bkg"};

	map<string, vector<double>> ValMap;
	if(ftype == "Nat" || ftype == "nat" || ftype == "Enr" || ftype == "enr") ValMap = RunWenqin(argList, fDS, fitMin, fitMax, ftype, theCut);
	else  ValMap = RunWenqin(argList, fDS, fitMin, fitMax, theCut, fidx);
	for(auto &kv : ValMap)
	{
		cout << kv.first << fidx << ":[" << kv.second[0] << "," << kv.second[1] << "," << kv.second[2] << "]"<< endl;
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
    if(fDS == 6) skimTree->Add("~/project/lat/latSkimDS*.root" );
  	else skimTree->Add(Form("~/project/lat/latSkimDS%d_*.root", fDS) );
  	// if(fDS == 6) skimTree->Add("~/project/latskimv3/latSkimDS*.root" );
  	// else skimTree->Add(Form("~/project/latskimv3/latSkimDS%d_*.root", fDS) );


  	fitter->LoadChainData(skimTree, theCut);

	// Construct PDF and do fit
	fitter->ConstructPDF();
	fitter->DoFit();
	fitter->DrawBasicShit(0.2, false, false);

	map<string, vector<double>> WenqinMap;
	for(auto &argN: argS){
		vector<double> WenqinParameter = fitter->GetVar( Form("%s", argN.c_str()) );
		WenqinMap[argN.c_str()] = WenqinParameter;
	}
	return WenqinMap;
}
