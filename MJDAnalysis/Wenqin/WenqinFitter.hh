#ifndef _WENQIN_FITTER_HH_
#define _WENQIN_FITTER_HH_

// Class that fits low energy crap

#include <string>

class RooDataSet;
class RooFitResult;
class RooPlot;
class RooRealVar;
class RooAbsReal;
class RooWorkspace;
class RooAbsPdf;
class RooMinimizer;
class RooMCStudy;
class TChain;

class WenqinFitter {
	
	public: 
		WenqinFitter();
		
                WenqinFitter(double fitMin, double fitMax) {fFitMin = fitMin; fFitMax = fitMax;}

		virtual ~WenqinFitter();

        // Constructs model PDF -- MUST be called after LoadData()!
        virtual void ConstructPDF();
        
        // Do the fit
        // Set minimizer type here also... there's like Minuit, Minuit2, probably more
        // Honestly Minuit2 and Minuit are the same shit, one's just potentially slightly faster
        // virtual void DoFit(std::string Minimizer = "Minuit2");
        virtual void DoFit(std::string Minimizer = "Minuit");

        // Draws and saves a plot of the fit as well as correlation matrix -- default binning is 0.2 keV
        // Binning is simply for visualization! 
        void DrawBasicShit(double binSize = 0.2);
        
        // Draws and saves contour plot -- arguments must have same name as in ConstructPDF()!
        // Parameters that become limited will take forever (as in never finish)
        void DrawContour(std::string argN1 = "Tritium", std::string argN2 = "Ge68");

        // This function calculates the energy resolution as a function of energy
        // According to the BDM PRL paper -- https://arxiv.org/abs/1612.00886
        double GetSigma(double energy);

        // This function uses RooMCStudy to generate toy MC and then fit the results
        void GenerateMCStudy(std::string argN = "Tritium", int nMC = 1000);

        // This function generates Toy MC data according to the best fit model and saves to a file
        void GenerateToyMC(std::string fileName);

        // Get Fit Result
        RooFitResult *GetFitResult() {return fFitResult;}

        // Get Workspace
        RooWorkspace *GetWorkspace() {return fFitWorkspace;}

        // Load data from skim TChain with a TCut
        // This assumes standard skimTree format
        void LoadChainData(TChain *skimTree, std::string theCut);

        // Creates, draws, and saves Profile Likelihood -- argument must have same name as in ConstructPDF()!
        // This is the ProfileNLL built into RooFit
        void ProfileNLL(std::string argN = "Tritium");

        // Creates, draws, and saves Profile Likelihood
        // This is a custom one to get a finer scan, probably takes much longer
        void ProfileNLLCustom(std::string argN = "Tritium");

        // Saves fit results into file
        void SaveShit(std::string outfileName = "Test.root");

        // Sets range for fit
        void SetFitRange(double fitMin, double fitMax);

        // Sets prefix of output files
        void SetSavePrefix(std::string savePrefix) {fSavePrefix = savePrefix;}

        // String to append to beginning of saved files
        std::string fSavePrefix;

        double fChiSquare;

	private:
                // Fit range -- in keV
		double fFitMin;
		double fFitMax;

        // Energy
        RooRealVar *fEnergy;

        // Real dataset (as in real data)
        RooDataSet *fRealData;

        // Toy MC data -- This is for studying systematics...
        RooDataSet *fMCData;
        RooMCStudy *fMCStudy;

        // Total PDF -- should change to RooSimultaneous for simultaneous fits
        RooAbsPdf *fModelPDF;

        // Minimizer
        RooMinimizer *fMinimizer;

        // NLL and ProfileNLL
        RooAbsReal *fNLL;
        RooAbsReal *fProfileNLL;

        // Saved fit result
        RooFitResult *fFitResult;

        // Fit workspace 
        RooWorkspace *fFitWorkspace;

};

#endif
