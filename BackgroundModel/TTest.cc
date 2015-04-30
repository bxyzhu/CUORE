#include "TTest.hh"
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
	BkgParM1[0] = new TBkgModelParameter("TeO2 2nbb", 0, 0, 0.0, 1.0, hAdapTeO22nuM1);
	hAdapTeO22nuM1->Draw();

	// BkgParM1[0]->GetHist()->Draw();

	cout << "Test name: " << BkgParM1[0]->GetParName() << endl;
}