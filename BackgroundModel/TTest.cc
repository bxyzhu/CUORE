#include "TTest.hh"
#include "TCanvas.h"
// gSystem->Load("TBkgModelSource_cc.so")
ClassImp(TTest)

TTest::TTest()
{
}

TTest::~TTest()
{}

void TTest::TestStuff()
{

	TBkgModelParameter *BkgParM1[10];

	BkgParM1[0] = new TBkgModelParameter("TeO2 2nbb", 0, 0, 0, 0.0, 1.0, hAdapTeO22nuM1, hAdapTeO22nuM2, hAdapTeO22nuM2Sum);
	// hAdapTeO22nuM1->Draw(); // works

	TCanvas *c1 = new TCanvas("c1","c1");
	BkgParM1[0]->GetHistM1()->Draw(); 
	TCanvas *c2 = new TCanvas("c2","c2");
	BkgParM1[0]->GetHistM2()->Draw(); 
	TCanvas *c2Sum = new TCanvas("c2Sum","c2Sum");
	BkgParM1[0]->GetHistM2Sum()->Draw(); 


	cout << "Test name: " << BkgParM1[0]->GetParName() << endl;
}
