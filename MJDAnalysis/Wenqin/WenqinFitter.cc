#include "WenqinFitter.hh"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooMinimizer.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooMCStudy.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"
#include "TStyle.h"
#include "TEntryList.h"
#include "TBranch.h"
#include <iostream>

using namespace RooFit;

WenqinFitter::WenqinFitter() : 
    fFitMin(2.), 
    fFitMax(100.), 
    fEnergy(nullptr), 
    fRealData(nullptr),
    fMCData(nullptr),
    fMCStudy(nullptr),
    fModelPDF(nullptr),
    fNLL(nullptr),
    fProfileNLL(nullptr),
    fMinimizer(nullptr),
    fFitResult(nullptr),
    fFitWorkspace(nullptr),
    fSavePrefix("FitResult")
{ }

WenqinFitter::~WenqinFitter()
{
	delete fEnergy;
	fEnergy = nullptr;
	
    delete fRealData;
	fRealData = nullptr;

	delete fModelPDF;
	fModelPDF = nullptr;

    delete fMCData;
    fMCData = nullptr;

    delete fMCStudy;
    fMCStudy = nullptr;

    delete fMinimizer;
    fMinimizer = nullptr;

    delete fNLL;
    fNLL = nullptr;

    delete fProfileNLL;
    fProfileNLL = nullptr;

	delete fFitResult;
	fFitResult = nullptr;

    delete fFitWorkspace;
    fFitWorkspace = nullptr;
}

// Constructs model PDF, use only after LoadData or else!
void WenqinFitter::ConstructPDF()
{
	if(fRealData == nullptr)
	{
		std::cout << "Error: Data not loaded! Will Segfault!" << std::endl;
		return;
	}

    // Energy shift -- 
    double fDeltaE = -0.05;

	std::string tritDir = "/Users/brianzhu/macros/code/MJDAnalysis/Axion";
    std::shared_ptr<TFile> tritFile( std::make_shared<TFile>(Form("%s/TritSpec.root", tritDir.c_str())) );
    std::shared_ptr<TH1D> tritSpec(dynamic_cast<TH1D*>(tritFile->Get("tritHist")));

    RooDataHist tritRooHist("trit", "Tritium Histogram", *fEnergy, Import(*tritSpec));
    // Because of Steve's histogram
    // The range of the histogram is maxed out at 50 keV, so need to reset range after loading histogram
    fEnergy->setRange(fFitMin, fFitMax);
    RooHistPdf tritPdf("tritPdf", "TritiumPdf", *fEnergy, tritRooHist, 1);

    RooRealVar polySlope("polySlope", "Background Slope", 0.00002, -0.2, 0.2);
    RooArgList polyList(polySlope);
    RooPolynomial BkgPoly("Background", "Linear Background function", *fEnergy, polyList);
    // RooPolynomial BkgPoly("Background", "Linear Background function", *fEnergy, RooArgList() );

    RooRealVar Mn54_mean("Mn54_mean", "Mn54_mean", 5.99 + fDeltaE);
    RooRealVar Mn54_sigma("Mn54_sigma", "Mn54_sigma", GetSigma(5.99 + fDeltaE));
    RooGaussian Mn54_gauss("Mn54_gauss", "Mn54 Gaussian", *fEnergy, Mn54_mean, Mn54_sigma);

    RooRealVar Fe55_mean("Fe55_mean", "Fe55_mean", 6.54 + fDeltaE);
    RooRealVar Fe55_sigma("Fe55_sigma", "Fe55_sigma", GetSigma(6.54 + fDeltaE));
    RooGaussian Fe55_gauss("Fe55_gauss", "Fe55 Gaussian", *fEnergy, Fe55_mean, Fe55_sigma);

    RooRealVar Co57_mean("Co57_mean", "Co57_mean", 7.11 + fDeltaE);
    RooRealVar Co57_sigma("Co57_sigma", "Co57_sigma", GetSigma(7.11 + fDeltaE));
    RooGaussian Co57_gauss("Co57_gauss", "Co57 Gaussian", *fEnergy, Co57_mean, Co57_sigma);

    RooRealVar Zn65_mean("Zn65_mean", "Zn65_mean", 8.98 + fDeltaE);
    RooRealVar Zn65_sigma("Zn65_sigma", "Zn65_sigma", GetSigma(8.98 + fDeltaE));
    RooGaussian Zn65_gauss("Zn65_gauss", "Zn65 Gaussian", *fEnergy, Zn65_mean, Zn65_sigma);

    RooRealVar Ge68_mean("Ge68_mean", "Ge68_mean", 10.37 + fDeltaE, 10.2, 10.5);
    RooRealVar Ge68_sigma("Ge68_sigma", "Ge68_sigma", GetSigma(10.37 + fDeltaE));
    RooGaussian Ge68_gauss("Ge68_gauss", "Ge68 Gaussian", *fEnergy, Ge68_mean, Ge68_sigma);

    // Normalization parameters
    // Make names pretty for plots
    RooRealVar num_trit("Tritium", "Tritium", 20.0, 0.0, 10000.);
    RooRealVar num_bkg("Bkg", "Background", 20.0, 0.0, 10000.);
    RooRealVar num_Mn54("Mn54", "Mn54", 1.0, 0.0, 5000.);
    RooRealVar num_Fe55("Fe55", "Fe55", 1.0, 0.0, 5000.);
    RooRealVar num_Co57("Co57", "Co57", 1.0, 0.0, 5000.);
    RooRealVar num_Zn65("Zn65", "Zn65", 1.0, 0.0, 5000.);
    RooRealVar num_Ge68("Ge68", "Ge68", 1.0, 0.0, 5000.);

    // Non-extended model
    // RooArgList shapes(tritPdf, BkgPoly, Mn54_gauss, Fe55_gauss, Co57_gauss, Zn65_gauss, Ge68_gauss);
    // RooArgList yields(num_trit, num_bkg, num_Mn54, num_Fe55, num_Co57, num_Zn65, num_Ge68);
    // RooAddPdf model("model", "total pdf", shapes, yields);

    // Extended PDF model
    RooExtendPdf tritPdfe("tritPdfe", "Extended trit", tritPdf, num_trit);
    RooExtendPdf BkgPolye("BkgPolye", "Extended BkgPoly", BkgPoly, num_bkg);
    RooExtendPdf Mn54_gausse("Mn54_gausse", "Extended Mn54_gauss", Mn54_gauss, num_Mn54);
    RooExtendPdf Fe55_gausse("Fe55_gausse", "Extended Fe55_gauss", Fe55_gauss, num_Fe55);
    RooExtendPdf Co57_gausse("Co57_gausse", "Extended Co57_gauss", Co57_gauss, num_Co57);
    RooExtendPdf Zn65_gausse("Zn65_gausse", "Extended Zn65_gauss", Zn65_gauss, num_Zn65);
    RooExtendPdf Ge68_gausse("Ge68_gausse", "Extended Ge68_gauss", Ge68_gauss, num_Ge68);

    RooArgList shapes(tritPdfe, BkgPolye, Mn54_gausse, Fe55_gausse, Co57_gausse, Zn65_gausse, Ge68_gausse);
    RooAddPdf model("model", "total pdf", shapes);

    // RooWorkspace is necessary for the model and parameters to be persistent
    // this is necessary because we created a bunch of objects that aren't persistent here
    fFitWorkspace = new RooWorkspace("fFitWorkspace", "Fit Workspace");
    
    // Add model to workspace -- also adds all of the normalization constants
    // If this step isn't done, a lot of the later functions won't work!
    fFitWorkspace->import(RooArgSet(model));
    fModelPDF = fFitWorkspace->pdf("model");
}

void WenqinFitter::DoFit(std::string Minimizer)
{
    // Create NLL (This is not a profile! When you draw it to one axis, it's just a projection!)
    fNLL = fModelPDF->createNLL(*fRealData, Extended(), NumCPU(4));
    
    // Create minimizer, fit model to data and save result
    fMinimizer = new RooMinimizer(*fNLL);
    fMinimizer->setMinimizerType(Form("%s", Minimizer.c_str()));
    fMinimizer->setStrategy(2);
    fMinimizer->migrad();
    fMinimizer->hesse();
    // Use all these, fk if I know they're any good
    fMinimizer->improve();
    fMinimizer->minos();

    fFitResult = fMinimizer->save();
    fFitWorkspace->import(*fFitResult);
}

void WenqinFitter::DrawBasicShit(double binSize)
{
    std::shared_ptr<TCanvas> cSpec( std::make_shared<TCanvas>("cSpec", "cSpec", 1100, 800) );
	RooPlot* frameFit = fEnergy->frame(Range(fFitMin, fFitMax), Bins((fFitMax - fFitMin)*1.0/binSize + 0.5));
    fRealData->plotOn(frameFit);
    fModelPDF->plotOn(frameFit, LineColor(kRed));
	frameFit->SetTitle("");
    frameFit->Draw();

    std::shared_ptr<TCanvas> cMatrix( std::make_shared<TCanvas>("cMatrix", "cMatrix", 1100, 800) );
    TH2D *fCorrMatrix = dynamic_cast<TH2D*>(fFitResult->correlationHist("Correlation Matrix"));
    fCorrMatrix->Draw("colz");
    
    // Save plots into workspace and pdf
    cSpec->SaveAs(Form("./Results/%s_Spec.pdf", fSavePrefix.c_str()) );
    // cMatrix->SaveAs(Form("./Results/%s_CorrMatrix.pdf", fSavePrefix.c_str()) );
    fFitWorkspace->import(*fCorrMatrix);
    fFitWorkspace->import(*frameFit);
}

void WenqinFitter::DrawContour(std::string argN1, std::string argN2)
{
    std::shared_ptr<TCanvas> cContour( std::make_shared<TCanvas>("cContour", "cContour", 1100, 800) );
    RooPlot *frameContour = fMinimizer->contour( *fFitWorkspace->var(Form("%s", argN1.c_str())), *fFitWorkspace->var(Form("%s", argN2.c_str())), 1, 2, 3);
    frameContour->SetTitle(Form("Contour of %s vs %s", argN2.c_str(), argN1.c_str()) );
    
    // Get range for plot -- make stuff pretty
    double meanx = fFitWorkspace->var(Form("%s", argN1.c_str()))->getValV();
    double uperrx = fFitWorkspace->var(Form("%s", argN1.c_str()))->getErrorHi();
    double lowerrx = fFitWorkspace->var(Form("%s", argN1.c_str()))->getErrorLo(); // This value is negative!
    double meany = fFitWorkspace->var(Form("%s", argN2.c_str()))->getValV();
    double uperry = fFitWorkspace->var(Form("%s", argN2.c_str()))->getErrorHi();
    double lowerry = fFitWorkspace->var(Form("%s", argN2.c_str()))->getErrorLo(); // This value is negative!

    frameContour->GetXaxis()->SetRangeUser(meanx + 4*lowerrx, meanx + 4*uperrx);
    frameContour->GetYaxis()->SetRangeUser(meany + 4*lowerry, meany + 4*uperry);
    frameContour->Draw();
    
    // Save plots into workspace and pdf
    cContour->SaveAs(Form("./Results/%s_Contour_%svs%s.pdf", fSavePrefix.c_str(), argN2.c_str(), argN1.c_str()));
    fFitWorkspace->import(*frameContour);
}


// Use after constructing the model and minimization!
// This is a simple MC study only generating parameter information and pulls, the toy MC data can be saved
void WenqinFitter::GenerateMCStudy(std::string argN, int nMC)
{
    fMCStudy = new RooMCStudy(*fModelPDF, *fEnergy, Extended(), FitOptions(Save()) );
    // Generate 2000 
    fMCStudy->generateAndFit(nMC);

    // Make test plots
    std::shared_ptr<TCanvas> cMCStudy( std::make_shared<TCanvas>("cMCStudy", "cMCStudy", 1100, 800) );
    RooPlot *frame1 = fMCStudy->plotParam( *fFitWorkspace->var(Form("%s", argN.c_str())), Bins(40) );
    RooPlot *frame2 = fMCStudy->plotError( *fFitWorkspace->var(Form("%s", argN.c_str())), Bins(40) );
    RooPlot *frame3 = fMCStudy->plotPull( *fFitWorkspace->var(Form("%s", argN.c_str())), Bins(40), FitGauss() );
    RooPlot* frame4 = fMCStudy->plotNLL(Bins(40)) ;

   cMCStudy->Divide(2,2) ;
   cMCStudy->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
   cMCStudy->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
   cMCStudy->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;
   cMCStudy->cd(4) ; gPad->SetLeftMargin(0.15) ; frame4->GetYaxis()->SetTitleOffset(1.4) ; frame4->Draw() ;
   cMCStudy->SaveAs(Form("./Results/%s_%s_MCStudy.pdf", fSavePrefix.c_str(), argN.c_str()) );

    // Save MC Study to workspace
    // fFitWorkspace->import(*fMCStudy);
}

// Use after constructing the model and minimization!
void WenqinFitter::GenerateToyMC(std::string fileName)
{
    std::shared_ptr<TFile> fOut( std::make_shared<TFile>( Form("./Data/%s_%s.root", fSavePrefix.c_str(), fileName.c_str()), "RECREATE" ) );
    fMCData = fModelPDF->generate( RooArgSet(*fEnergy), Name("Toy_dataset"), Extended());

    // Save data to workspace and write to a file
    RooWorkspace *workspace = new RooWorkspace("workspace", "workspace");
    workspace->import(*fMCData);
    fOut->cd();
    workspace->Write();
    fOut->Close();
}


// Gets resolution, function and parameters from BDM PRL paper
// https://arxiv.org/abs/1612.00886
double WenqinFitter::GetSigma(double energy)
{
  	double sig0 = 0.16, F=0.11, eps=0.00296;
  	double sig = std::sqrt(std::pow(sig0,2) + eps * F * energy);
	
	return sig;
}

// Load data from file
void WenqinFitter::LoadData(std::string fileName, std::string treeName, std::string parName)
{
    std::shared_ptr<TFile> input_file( std::make_shared<TFile>(fileName.c_str()));
    std::shared_ptr<TTree> skimTree( dynamic_cast<TTree*>(input_file->Get(Form("%s", treeName.c_str()) )) );

    // Can and perhaps should split the data up by channel in a more complicated fit
    fEnergy = new RooRealVar(Form("%s", parName.c_str()), Form("%s", parName.c_str()), fFitMin, fFitMax, "keV");
    // Stupid constructor of RooDataSet needs pointer
    fRealData = new RooDataSet("data", "data", skimTree.get(), RooArgSet(*fEnergy));
}

// Assumes standard skim format -- converts stuff from vector<double> to scalar
// Also assumes trapENFCal is the parameter of choice
void WenqinFitter::LoadChainData(TChain *skimTree, std::string theCut)
{
    // First get TEntryList with TCut
    skimTree->Draw(">> elist", Form("%s", theCut.c_str()), "entrylist goff");
    std::shared_ptr<TEntryList> elist( dynamic_cast<TEntryList*>(gDirectory->Get("elist")) );
    skimTree->SetEntryList(elist.get());

    // I found it easier to work like this rather than with a TTreeReader... Ian probably hates me
    std::vector<double> *ftrapENFCal = nullptr;
    std::vector<int> *fchannel = nullptr;
    skimTree->SetBranchAddress("trapENFCal", &ftrapENFCal);
    skimTree->SetBranchAddress("channel", &fchannel);

    // Create and fill a dummy TTree to load into the RooDataSet
    double trapENFCal = 0;
    int channel = 0;
    int treeNum = 0;
    TTree *dummyTree = new TTree("dummyTree", "Tree for RooDataSet");
    dummyTree->Branch("trapENFCal", &trapENFCal, "trapENFCal/D");
    dummyTree->Branch("channel", &channel, "channel/I");
    for(int i = 0; i < elist->GetN(); i++)
    {
        int treeEntry = elist->GetEntryAndTree(i, treeNum);
        skimTree->GetEntry( treeEntry + skimTree->GetTreeOffset()[treeNum] );
        
        for(int j = 0; j < ftrapENFCal->size(); j++)
        {
            trapENFCal = ftrapENFCal->at(j);
            channel = fchannel->at(j);
        }
        dummyTree->Fill();
    }

    std::cout << "Dummy Tree filled entries: " << dummyTree->GetEntries() << std::endl;

    // Can and perhaps should split the data up by channel in a more complicated fit
    fEnergy = new RooRealVar("trapENFCal", "trapENFCal", fFitMin, fFitMax, "keV");
    fRealData = new RooDataSet("data", "data", dummyTree, RooArgSet(*fEnergy));
}


void WenqinFitter::ProfileNLL(std::string argN)
{
    // Create Profile NLL (This should take a while since the actual profiling goes on here)
    // Variable name here must match name in "ConstructPDF()"
    fProfileNLL = fNLL->createProfile(RooArgSet(*fFitWorkspace->var(Form("%s", argN.c_str()) )));

    // Draw both NLL and Profile NLL and save as PDF
    std::shared_ptr<TCanvas> cNLL( std::make_shared<TCanvas>("cNLL", "cNLL", 900, 600) );
    RooPlot* frameProfileNLL = fFitWorkspace->var(Form("%s", argN.c_str()))->frame(Title(Form("%s Profile Likelihood", argN.c_str()) ));
    fProfileNLL->plotOn(frameProfileNLL, LineColor(kRed));
    double mean = fFitWorkspace->var(Form("%s", argN.c_str()))->getValV();
    double uperr = fFitWorkspace->var(Form("%s", argN.c_str()))->getErrorHi();
    double lowerr = fFitWorkspace->var(Form("%s", argN.c_str()))->getErrorLo(); // This value is negative!
    frameProfileNLL->GetXaxis()->SetRangeUser(mean + 10*lowerr, mean + 10*uperr);
    frameProfileNLL->GetYaxis()->SetRangeUser(0, 50); // 0 to 50 should be enough for profile likelihood
    frameProfileNLL->Draw();

    // Save plot into workspace and pdf
    cNLL->SaveAs(Form("./Results/%s_%sNLL.pdf", fSavePrefix.c_str(), argN.c_str()) );
    fFitWorkspace->import(*fProfileNLL);
    fFitWorkspace->import(*frameProfileNLL);
}

void WenqinFitter::SaveShit(std::string outfileName)
{
    std::shared_ptr<TFile> fOut( std::make_shared<TFile>( Form("./Results/%s_%s", fSavePrefix.c_str(), outfileName.c_str()), "RECREATE" ) );
    fOut->cd();
    fFitWorkspace->Write();
    fOut->Close();
}

void WenqinFitter::SetFitRange(double fitMin, double fitMax)
{
    fFitMin = fitMin;
    fFitMax = fitMax;
}
