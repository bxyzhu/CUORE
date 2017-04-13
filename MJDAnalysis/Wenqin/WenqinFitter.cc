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
#include "RooHist.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"
#include "TLine.h"
#include "TTree.h"
#include "TStyle.h"
#include "TEntryList.h"
#include "TBranch.h"
#include "TPaveText.h"
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
    fChiSquare(0),
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

    // RooRealVar polySlope("polySlope", "Background Slope", 0.00002, -0.2, 0.2);
    // RooArgList polyList(polySlope);
    // RooPolynomial BkgPoly("Background", "Linear Background function", *fEnergy, polyList);
    RooPolynomial BkgPoly("Background", "Linear Background function", *fEnergy, RooArgList() );

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
    
    RooRealVar Pb210_mean("Pb210_mean", "Pb210_mean", 46.54 + fDeltaE);
    RooRealVar Pb210_sigma("Pb210_sigma", "Pb210_sigma", GetSigma(46.54 + fDeltaE));
    RooGaussian Pb210_gauss("Pb210_gauss", "Pb210 Gaussian", *fEnergy, Pb210_mean, Pb210_sigma);


    // Normalization parameters
    // Make names pretty for plots
    RooRealVar num_trit("Tritium", "Tritium", 20.0, 0.0, 10000.);
    RooRealVar num_bkg("Bkg", "Background", 20.0, 0.0, 10000.);
    RooRealVar num_Mn54("Mn54", "Mn54", 1.0, 0.0, 5000.);
    RooRealVar num_Fe55("Fe55", "Fe55", 1.0, 0.0, 5000.);
    RooRealVar num_Co57("Co57", "Co57", 1.0, 0.0, 5000.);
    RooRealVar num_Zn65("Zn65", "Zn65", 1.0, 0.0, 5000.);
    RooRealVar num_Ge68("Ge68", "Ge68", 1.0, 0.0, 5000.);
    RooRealVar num_Pb210("Pb210", "Pb210", 1.0, 0.0, 5000.);
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
    RooExtendPdf Pb210_gausse("Pb210_gausse", "Extended Pb210_gauss", Pb210_gauss, num_Pb210);

    RooArgList shapes(tritPdfe, BkgPolye, Mn54_gausse, Fe55_gausse, Co57_gausse, Zn65_gausse, Ge68_gausse, Pb210_gausse);
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
    
    // Add chi-square to the plot
    // The chi-square is fairly meaningless as it's an unbinned fit...
    fChiSquare = frameFit->chiSquare(9);
    TPaveText *leg = new TPaveText(0.65, 0.78, 0.88, .88, "NDC");
    leg->SetTextFont(133);
    leg->SetFillColor(0);
    leg->SetBorderSize(1);
    leg->SetTextSize(22);
    leg->AddText(Form("#chi^{2}/NDF = %.3f" ,fChiSquare ) );
    frameFit->addObject(leg);
    frameFit->Draw();

    std::shared_ptr<TCanvas> cMatrix( std::make_shared<TCanvas>("cMatrix", "cMatrix", 1100, 800) );
    TH2D *fCorrMatrix = dynamic_cast<TH2D*>(fFitResult->correlationHist("Correlation Matrix"));
    fCorrMatrix->Draw("colz");
    
    // Save plots into workspace and pdf
    cSpec->SaveAs(Form("./Results/%s_Spec.pdf", fSavePrefix.c_str()) );
    cMatrix->SaveAs(Form("./Results/%s_CorrMatrix.pdf", fSavePrefix.c_str()) );
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
    // Right now I'm saving the fit output
    fMCStudy = new RooMCStudy(*fModelPDF, *fEnergy, Extended(), Silence(), FitOptions(Save()));
    fMCStudy->generateAndFit(nMC);

    // Get parameter values from first fit... these methods suck
    double parVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find(Form("%s", argN.c_str())))->getValV();
    double parErr = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find(Form("%s", argN.c_str())))->getError();

    // Make test plots
    std::shared_ptr<TCanvas> cMCStudy( std::make_shared<TCanvas>("cMCStudy", "cMCStudy", 1100, 800) );
    RooPlot *frame1 = fMCStudy->plotParam( *fFitWorkspace->var(Form("%s", argN.c_str())), Bins(50) );
    RooPlot *frame2 = fMCStudy->plotError( *fFitWorkspace->var(Form("%s", argN.c_str())), FrameRange(parErr-0.5*parErr, parErr+0.5*parErr), Bins(50));
    RooPlot *frame3 = fMCStudy->plotPull( *fFitWorkspace->var(Form("%s", argN.c_str())), FrameRange(-5, 5), Bins(50));
    RooPlot *frame4 = fMCStudy->plotNLL(Bins(50));

    // Add PaveTexts with values and such
    TPaveText *legParam = new TPaveText(0.60, 0.78, 0.89, 0.89, "NDC");
    legParam->SetTextFont(133);
    legParam->SetFillColor(0);
    legParam->SetBorderSize(1);
    legParam->SetTextSize(14);
    legParam->AddText(Form("Best Fit: %.3f #pm %.3f", parVal, parErr));
    frame1->addObject(legParam);

    // Workaround because fitting built into plotPull is terrible...
    RooHist *hist = frame3->getHist();
    hist->Fit("gaus", "ME");
    TF1 *gaus = hist->GetFunction("gaus");
    TPaveText *legpull = new TPaveText(0.60, 0.75, 0.89, 0.89, "NDC");
    legpull->SetTextFont(133);
    legpull->SetFillColor(0);
    legpull->SetBorderSize(1);
    legpull->SetTextSize(14);
    legpull->AddText(Form("Pull Mean: %.3f #pm %.3f", gaus->GetParameter(1), gaus->GetParError(1)) );
    legpull->AddText(Form("Pull Sigma: %.3f #pm %.3f", gaus->GetParameter(2), gaus->GetParError(2)) );
    frame3->addObject(legpull);

    TPaveText *legNLL = new TPaveText(0.60, 0.78, 0.89, 0.89, "NDC");
    legNLL->SetTextFont(133);
    legNLL->SetFillColor(0);
    legNLL->SetBorderSize(1);
    legNLL->SetTextSize(14);
    legNLL->AddText(Form("Best Fit NLL: %.3f", fFitResult->minNll()));
    frame4->addObject(legNLL);

    TLine l1;
    l1.SetLineColor(kBlue);
    l1.SetLineWidth(2);
    l1.SetLineStyle(3);

    cMCStudy->Divide(2,2);
    cMCStudy->cd(1); gPad->SetLeftMargin(0.15); frame1->GetYaxis()->SetTitleOffset(1.4); frame1->Draw();
    // Draw a line at best fit position
    l1.DrawLine(parVal, frame1->GetMinimum(), parVal, frame1->GetMaximum());
    cMCStudy->cd(2); gPad->SetLeftMargin(0.15); frame2->GetYaxis()->SetTitleOffset(1.4); frame2->Draw();
    l1.DrawLine(parErr, frame2->GetMinimum(), parErr, frame2->GetMaximum());
    cMCStudy->cd(3); gPad->SetLeftMargin(0.15); frame3->GetYaxis()->SetTitleOffset(1.4); frame3->Draw();
    cMCStudy->cd(4); gPad->SetLeftMargin(0.15); frame4->GetYaxis()->SetTitleOffset(1.4); frame4->Draw();
    // Draw a line at minimum NLL position
    l1.DrawLine(fFitResult->minNll(), frame4->GetMinimum(), fFitResult->minNll(), frame4->GetMaximum());
    
    // Save MC Study to plot and workspace
    // fFitWorkspace->import(*fMCStudy);
    cMCStudy->SaveAs(Form("./Results/%s_%s_MCStudy.pdf", fSavePrefix.c_str(), argN.c_str()) );
}

// Use after constructing the model and minimization!
// How useful is this when there's RooMCStudy?
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
// In the future it should just be a convolution with all the other PDFs
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
    fRealData = new RooDataSet("data", "data", &*skimTree, RooArgSet(*fEnergy));
}

// Assumes standard skim format -- converts stuff from vector<double> to scalar
// Also assumes trapENFCal is the parameter of choice
void WenqinFitter::LoadChainData(TChain *skimTree, std::string theCut)
{
    // First get TEntryList with TCut
    skimTree->Draw(">> elist", Form("%s", theCut.c_str()), "entrylist goff");
    std::shared_ptr<TEntryList> elist( dynamic_cast<TEntryList*>(gDirectory->Get("elist")) );
    // This works
    skimTree->SetEntryList(&*elist);

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

// Implemented now in RooStats rather than RooFit
void WenqinFitter::ProfileNLL(std::string argN)
{
    // Draw both NLL and Profile NLL and save as PDF
    std::shared_ptr<TCanvas> cNLL( std::make_shared<TCanvas>("cNLL", "cNLL", 900, 600) );

    // Best fit value -- just using this to set range
    double parVal = dynamic_cast<RooRealVar*>(fFitResult->floatParsFinal().find(Form("%s", argN.c_str())))->getValV();

    RooStats::ProfileLikelihoodCalculator plc(*fRealData, *fModelPDF, RooArgSet(*fFitWorkspace->var(Form("%s", argN.c_str()))) );
    // Set 1 sigma interval
    plc.SetConfidenceLevel(0.683);

    RooStats::LikelihoodInterval *interval = plc.GetInterval();
    double lowerLimit = interval->LowerLimit(*fFitWorkspace->var(Form("%s", argN.c_str())));
    double upperLimit = interval->UpperLimit(*fFitWorkspace->var(Form("%s", argN.c_str())));

    RooStats::LikelihoodIntervalPlot plot(interval);
    plot.SetRange(parVal - 1.5*(parVal - lowerLimit), parVal + 1.5*(upperLimit-parVal) );
    plot.Draw();
    cNLL->SaveAs(Form("./Results/%s_%sNLL.pdf", fSavePrefix.c_str(), argN.c_str()) );
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
