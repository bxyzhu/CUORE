// Macro takes outputted g2root files and saves the histograms (M1 and M2 for now) into a ROOT file for easier usage
// File names are subject to change (unfortunately) - Nov 5, 2014

TChain *LoadMC(std::string dDir, std::string dLocation, std::string dSource)
{
    TChain *outTree = new TChain("outTree");
    outTree->Add(Form("%s%s-%s.root", dDir.c_str(), dLocation.c_str(), dSource.c_str()));

    cout << "Loading: " << Form("%s-%s", dLocation.c_str(), dSource.c_str()) << endl;

    return outTree;

}


TH1D *SmearMC(TH1D *hMC, TH1D *hSMC, double dSigma)
{
	// Reset previously Modeled histogram
	hSMC->Reset();
	cout << "Smearing: " << hMC->GetName() << endl;

	double dNorm;
  	// double dNorm2;
	double dSmearedValue;

	double dBinSize = 1;

	TF1 *gaus = new TF1("gaus","gaus(0)", 0, 10000);
  	// TF1 *gaus2 = new TF1("gaus2","gaus(0)", 0, 10000);


	for(int i = 499; i < hMC->GetNbinsX(); i++)
	{
		// cout << "bin " << i << endl;
		for(int j = 499; j < hMC->GetNbinsX(); j++)
		{
			// Normalization of gaussian = (bsin size * Area of bin j in MC) / Sigma of bin j (fit function evaluated at bin center)
			// dNorm = (853.3/1215.8)*dBinSize*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*dRes1);
      		// dNorm2 = (362.5/1215.8)*dBinSize*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*dRes2);

			// dNorm = (304.5/338.74)*dBinSize*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*dRes1);
   //    		dNorm2 = (34.24/338.74)*dBinSize*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*dRes2);

			dNorm = dBinSize*hMC->GetBinContent(j)/(sqrt(2*TMath::Pi())*dSigma);


			// Set parameters of gaussian ... 2nd gaussian *slightly* shifted... not sure if this works
			gaus->SetParameters(dNorm, hMC->GetBinCenter(j), dSigma);
      		// gaus2->SetParameters(dNorm2, hMC->GetBinCenter(j)-1, dRes2);

			// Smeared contribution from gaussian centered at bin j for bin i 
			// dSmearedValue = gaus->Eval(hSMC->GetBinCenter(i)) + gaus2->Eval(hSMC->GetBinCenter(i));
			dSmearedValue = gaus->Eval(hSMC->GetBinCenter(i));

			// Fill bin i with contribution from gaussian centered at bin j
			// cout << "Filling Bin: " << i << " with " << dSmearedValue << endl;
			hSMC->Fill(hSMC->GetBinCenter(i), dSmearedValue);
		}
	}
	
	return hSMC;
}


void NormalizePDF(TH1D *h1, int minE, int maxE)
{

  double dIntegral = 0;
  int dBinSize = 1;

  // bin 0 = underflow, bin dNBins = last bin with upper-edge xup Excluded
  dIntegral = h1->Integral(minE/dBinSize, maxE/dBinSize);
  // cout << "Integral for " << h1->GetTitle() << " :" << dIntegral << endl;

  cout << "Normalizing " << h1->GetName() << endl;

  // Make sure integral isn't 0
  // If it is 0, clear model... 
  if(dIntegral == 0)
  {
    cout << Form("Integral of %s is 0, resetting histogram", h1->GetName()) << endl;
    h1->Reset();
  }

  if(dIntegral != 0)
  {
    h1->Scale(1.0/dIntegral);
  }
}

void NormalizePDFs(TH1D *h1, TH1D *h2, TH1D *h3, int minE, int maxE)
{
  double dIntegral1 = 0;
  double dIntegral2 = 0;
  double dIntegral3 = 0;

  int dBinSize = 1;

  // bin 0 = underflow, bin dNBins = last bin with upper-edge xup Excluded
  dIntegral1 = h1->Integral(minE/dBinSize, maxE/dBinSize);
  dIntegral2 = h2->Integral(minE/dBinSize, maxE/dBinSize);
  dIntegral3 = h3->Integral(minE/dBinSize, maxE/dBinSize);
  // cout << "Integral for " << h1->GetTitle() << " :" << dIntegral << endl;

  cout << "Normalizing " << h1->GetName() << " : " << h2->GetName() << " : " << h3->GetName() << endl;

  // Make sure integral isn't 0
  // If it is 0, clear model... 
  if(dIntegral1 == 0)
  {
    cout << Form("Integral of %s is 0, resetting histogram", h1->GetName()) << endl;
    h1->Reset();
    h2->Reset();
    h3->Reset();
  }

  if(dIntegral1 != 0)
  {
    h1->Scale(1.0/dIntegral1);
    h2->Scale(1.0/dIntegral1);
    h3->Scale(1.0/dIntegral1);
  }
}

void NormalizePDFPair(TH1D *h1, TH1D *h2, int minE, int maxE)
{
  double dIntegral1 = 0;
  double dIntegral2 = 0;
  int dBinSize = 1;

  // bin 0 = underflow, bin dNBins = last bin with upper-edge xup Excluded
  dIntegral1 = h1->Integral(minE/dBinSize, maxE/dBinSize);
  dIntegral2 = h2->Integral(minE/dBinSize, maxE/dBinSize);
  // cout << "Integral for " << h1->GetTitle() << " :" << dIntegral << endl;

  cout << "Normalizing " << h1->GetName() << " : " << h2->GetName() << endl;


  // Make sure integral isn't 0
  // If it is 0, clear model... 
  if(dIntegral1 == 0)
  {
    cout << Form("Integral of %s is 0, resetting histogram", h1->GetName()) << endl;
    h1->Reset();
    h2->Reset();
  }

  if(dIntegral1 != 0)
  {
    h1->Scale(1.0/dIntegral1);
    h2->Scale(1.0/dIntegral1);
  }
}


void SaveHistogramsReducedBulk()
{
	std::string sDataDir = "/cuore/user/zhubrian/MC/Bkg/Unsmeared/";

	double dMinEnergy = 0;
	double dMaxEnergy = 10000;
	int dBinSize = 1;
	int dNBins = (dMaxEnergy - dMinEnergy)/dBinSize;

	//2614.53
	// double dB = 0.53
	double dSigma[52] = {3.6309, 2.23914, 1.85379, 2.12793, 4.41549, 1.74575, 2.34558, 2.37479, 1.94436, 8.20342, 2.17113, 2.10644, 2.2757, 2.7982, 4.02611, 2.87497, 1.78139, 1.92231, 1.87181, 2.95697, 1.9673, 4.89725, 2.06469, 2.33388, 2.18577, 2.36421, 7.53857, 2.03118, 2.86947, 2.57709, 2.43256, 2.01672, 3.1208, 2.67247, 5.11125, 2.20974, 3.177, 2.43366, 2.55086, 2.23907, 1.94238, 2.14217, 2.74023, 2.05892, 1.97427, 2.4128, 2.48429, 3.50036, 3, 3.12728, 2.50535, 3.07077};

	// Unsmeared

	TH1D *uTeO20nuM1 = new TH1D("uTeO20nuM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO22nuM1 = new TH1D("uTeO22nuM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2co60M1 = new TH1D("uTeO2co60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2k40M1 = new TH1D("uTeO2k40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2po210M1 = new TH1D("uTeO2po210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2ra228pb208M1 = new TH1D("uTeO2ra228pb208M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2th230onlyM1 = new TH1D("uTeO2th230onlyM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2th232onlyM1 = new TH1D("uTeO2th232onlyM1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *uCuFrameco60M1 = new TH1D("uCuFrameco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uCuFramek40M1 = new TH1D("uCuFramek40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uCuFrameth232M1 = new TH1D("uCuFrameth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uCuFrameu238M1 = new TH1D("uCuFrameu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *u600mKco60M1 = new TH1D("u600mKco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *u600mKk40M1 = new TH1D("u600mKk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *u600mKth232M1 = new TH1D("u600mKth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *u600mKu238M1 = new TH1D("u600mKu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *uPbRomco60M1 = new TH1D("uPbRomco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uPbRomk40M1 = new TH1D("uPbRomk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uPbRomth232M1 = new TH1D("uPbRomth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uPbRomu238M1 = new TH1D("uPbRomu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *uOVCco60M1 = new TH1D("uOVCco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uOVCk40M1 = new TH1D("uOVCk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uOVCth232M1 = new TH1D("uOVCth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uOVCu238M1 = new TH1D("uOVCu238M1", "", dNBins, dMinEnergy, dMaxEnergy);	

	TH1D *uExtPbbi210M1 = new TH1D("uExtPbbi210M1", "", dNBins, dMinEnergy, dMaxEnergy);


	TH1D *uTeO20nuM2 = new TH1D("uTeO20nuM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO22nuM2 = new TH1D("uTeO22nuM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2co60M2 = new TH1D("uTeO2co60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2k40M2 = new TH1D("uTeO2k40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2po210M2 = new TH1D("uTeO2po210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2ra228pb208M2 = new TH1D("uTeO2ra228pb208M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2th230onlyM2 = new TH1D("uTeO2th230onlyM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2th232onlyM2 = new TH1D("uTeO2th232onlyM2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *uCuFrameco60M2 = new TH1D("uCuFrameco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uCuFramek40M2 = new TH1D("uCuFramek40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uCuFrameth232M2 = new TH1D("uCuFrameth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uCuFrameu238M2 = new TH1D("uCuFrameu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *u600mKco60M2 = new TH1D("u600mKco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *u600mKk40M2 = new TH1D("u600mKk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *u600mKth232M2 = new TH1D("u600mKth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *u600mKu238M2 = new TH1D("u600mKu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *uPbRomco60M2 = new TH1D("uPbRomco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uPbRomk40M2 = new TH1D("uPbRomk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uPbRomth232M2 = new TH1D("uPbRomth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uPbRomu238M2 = new TH1D("uPbRomu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *uOVCco60M2 = new TH1D("uOVCco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uOVCk40M2 = new TH1D("uOVCk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uOVCth232M2 = new TH1D("uOVCth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uOVCu238M2 = new TH1D("uOVCu238M2", "", dNBins, dMinEnergy, dMaxEnergy);	

	TH1D *uExtPbbi210M2 = new TH1D("uExtPbbi210M2", "", dNBins, dMinEnergy, dMaxEnergy);



	// Fill these

	TH1D *hTeO20nuM1 = new TH1D("hTeO20nuM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO22nuM1 = new TH1D("hTeO22nuM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2co60M1 = new TH1D("hTeO2co60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2k40M1 = new TH1D("hTeO2k40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2po210M1 = new TH1D("hTeO2po210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra228pb208M1 = new TH1D("hTeO2ra228pb208M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th230onlyM1 = new TH1D("hTeO2th230onlyM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232onlyM1 = new TH1D("hTeO2th232onlyM1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuFrameco60M1 = new TH1D("hCuFrameco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramek40M1 = new TH1D("hCuFramek40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameth232M1 = new TH1D("hCuFrameth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameu238M1 = new TH1D("hCuFrameu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *h600mKco60M1 = new TH1D("h600mKco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKk40M1 = new TH1D("h600mKk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKth232M1 = new TH1D("h600mKth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKu238M1 = new TH1D("h600mKu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hPbRomco60M1 = new TH1D("hPbRomco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomk40M1 = new TH1D("hPbRomk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomth232M1 = new TH1D("hPbRomth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomu238M1 = new TH1D("hPbRomu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hOVCco60M1 = new TH1D("hOVCco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCk40M1 = new TH1D("hOVCk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCth232M1 = new TH1D("hOVCth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCu238M1 = new TH1D("hOVCu238M1", "", dNBins, dMinEnergy, dMaxEnergy);	

	TH1D *hExtPbbi210M1 = new TH1D("hExtPbbi210M1", "", dNBins, dMinEnergy, dMaxEnergy);


	TH1D *hTeO20nuM2 = new TH1D("hTeO20nuM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO22nuM2 = new TH1D("hTeO22nuM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2co60M2 = new TH1D("hTeO2co60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2k40M2 = new TH1D("hTeO2k40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2po210M2 = new TH1D("hTeO2po210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra228pb208M2 = new TH1D("hTeO2ra228pb208M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th230onlyM2 = new TH1D("hTeO2th230onlyM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232onlyM2 = new TH1D("hTeO2th232onlyM2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuFrameco60M2 = new TH1D("hCuFrameco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramek40M2 = new TH1D("hCuFramek40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameth232M2 = new TH1D("hCuFrameth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameu238M2 = new TH1D("hCuFrameu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *h600mKco60M2 = new TH1D("h600mKco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKk40M2 = new TH1D("h600mKk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKth232M2 = new TH1D("h600mKth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKu238M2 = new TH1D("h600mKu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hPbRomco60M2 = new TH1D("hPbRomco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomk40M2 = new TH1D("hPbRomk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomth232M2 = new TH1D("hPbRomth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomu238M2 = new TH1D("hPbRomu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hOVCco60M2 = new TH1D("hOVCco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCk40M2 = new TH1D("hOVCk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCth232M2 = new TH1D("hOVCth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCu238M2 = new TH1D("hOVCu238M2", "", dNBins, dMinEnergy, dMaxEnergy);	

	TH1D *hExtPbbi210M2 = new TH1D("hExtPbbi210M2", "", dNBins, dMinEnergy, dMaxEnergy);


	TChain *tTeO20nu = LoadMC(sDataDir.c_str(), "TeO2", "0nu_1um");
	TChain *tTeO22nu = LoadMC(sDataDir.c_str(), "TeO2", "2nu");
	// TChain *tTeO2co60 = LoadMC(sDataDir.c_str(), "TeO2", "co60");
	// TChain *tTeO2k40 = LoadMC(sDataDir.c_str(), "TeO2", "k40");
	// TChain *tTeO2po210 = LoadMC(sDataDir.c_str(), "TeO2", "po210");
	// TChain *tTeO2ra228pb208 = LoadMC(sDataDir.c_str(), "TeO2", "ra228-pb208");
	// TChain *tTeO2th230only = LoadMC(sDataDir.c_str(), "TeO2", "th230_only");
	// TChain *tTeO2th232only = LoadMC(sDataDir.c_str(), "TeO2", "th232_only");

	// TChain *tCuFrameco60 = LoadMC(sDataDir.c_str(), "CuFrame", "co60");
	// TChain *tCuFramek40 = LoadMC(sDataDir.c_str(), "CuFrame", "k40");
	// TChain *tCuFrameth232 = LoadMC(sDataDir.c_str(), "CuFrame", "th232");
	// TChain *tCuFrameu238 = LoadMC(sDataDir.c_str(), "CuFrame", "u238");

	// TChain *t600mKco60 = LoadMC(sDataDir.c_str(), "600mK", "co60");
	// TChain *t600mKk40 = LoadMC(sDataDir.c_str(), "600mK", "k40");
	// TChain *t600mKth232 = LoadMC(sDataDir.c_str(), "600mK", "th232");
	// TChain *t600mKu238 = LoadMC(sDataDir.c_str(), "600mK", "u238");

	// TChain *tOVCco60 = LoadMC(sDataDir.c_str(), "OVC", "co60");
	// TChain *tOVCk40 = LoadMC(sDataDir.c_str(), "OVC", "k40");
	// TChain *tOVCth232 = LoadMC(sDataDir.c_str(), "OVC", "th232");
	// TChain *tOVCu238 = LoadMC(sDataDir.c_str(), "OVC", "u238");

	// TChain *tPbRomco60 = LoadMC(sDataDir.c_str(), "PbRom", "co60");
	// TChain *tPbRomk40 = LoadMC(sDataDir.c_str(), "PbRom", "k40");
	// TChain *tPbRomth232 = LoadMC(sDataDir.c_str(), "PbRom", "th232");
	// TChain *tPbRomu238 = LoadMC(sDataDir.c_str(), "PbRom", "u238");	

	// TChain *tExtPbbi210 = LoadMC(sDataDir.c_str(), "ExtPb4cm_210BiBStot", "1.23T");


	// tTeO20nu->Project("uTeO20nuM1", "Ener2", "Multiplicity==1");
	// tTeO22nu->Project("uTeO22nuM1", "Ener2", "Multiplicity==1");
	// tTeO2co60->Project("uTeO2co60M1", "Ener2", "Multiplicity==1");
	// tTeO2k40->Project("uTeO2k40M1", "Ener2", "Multiplicity==1");
	// tTeO2po210->Project("uTeO2po210M1", "Ener2", "Multiplicity==1");
	// tTeO2ra228pb208->Project("uTeO2ra228pb208M1", "Ener2", "Multiplicity==1");
	// tTeO2th230only->Project("uTeO2th230onlyM1", "Ener2", "Multiplicity==1");
	// tTeO2th232only->Project("uTeO2th232onlyM1", "Ener2", "Multiplicity==1");

	// tCuFrameco60->Project("uCuFrameco60M1", "Ener2", "Multiplicity==1");
	// tCuFramek40->Project("uCuFramek40M1", "Ener2", "Multiplicity==1");
	// tCuFrameth232->Project("uCuFrameth232M1", "Ener2", "Multiplicity==1");
	// tCuFrameu238->Project("uCuFrameu238M1", "Ener2", "Multiplicity==1");

	// t600mKco60->Project("u600mKco60M1", "Ener2", "Multiplicity==1");
	// t600mKk40->Project("u600mKk40M1", "Ener2", "Multiplicity==1");
	// t600mKth232->Project("u600mKth232M1", "Ener2", "Multiplicity==1");
	// t600mKu238->Project("u600mKu238M1", "Ener2", "Multiplicity==1");

	// tOVCco60->Project("uOVCco60M1", "Ener2", "Multiplicity==1");
	// tOVCk40->Project("uOVCk40M1", "Ener2", "Multiplicity==1");
	// tOVCth232->Project("uOVCth232M1", "Ener2", "Multiplicity==1");
	// tOVCu238->Project("uOVCu238M1", "Ener2", "Multiplicity==1");

	// tPbRomco60->Project("uPbRomco60M1", "Ener2", "Multiplicity==1");
	// tPbRomk40->Project("uPbRomk40M1", "Ener2", "Multiplicity==1");
	// tPbRomth232->Project("uPbRomth232M1", "Ener2", "Multiplicity==1");
	// tPbRomu238->Project("uPbRomu238M1", "Ener2", "Multiplicity==1");

	// tExtPbbi210->Project("uExtPbbi210M1", "Ener2", "Multiplicity==1");


	// tTeO20nu->Project("uTeO20nuM2", "Ener2", "Multiplicity==2");
	// tTeO22nu->Project("uTeO22nuM2", "Ener2", "Multiplicity==2");
	// tTeO2co60->Project("uTeO2co60M2", "Ener2", "Multiplicity==2");
	// tTeO2k40->Project("uTeO2k40M2", "Ener2", "Multiplicity==2");
	// tTeO2po210->Project("uTeO2po210M2", "Ener2", "Multiplicity==2");
	// tTeO2ra228pb208->Project("uTeO2ra228pb208M2", "Ener2", "Multiplicity==2");
	// tTeO2th230only->Project("uTeO2th230onlyM2", "Ener2", "Multiplicity==2");
	// tTeO2th232only->Project("uTeO2th232onlyM2", "Ener2", "Multiplicity==2");

	// tCuFrameco60->Project("uCuFrameco60M2", "Ener2", "Multiplicity==2");
	// tCuFramek40->Project("uCuFramek40M2", "Ener2", "Multiplicity==2");
	// tCuFrameth232->Project("uCuFrameth232M2", "Ener2", "Multiplicity==2");
	// tCuFrameu238->Project("uCuFrameu238M2", "Ener2", "Multiplicity==2");

	// t600mKco60->Project("u600mKco60M2", "Ener2", "Multiplicity==2");
	// t600mKk40->Project("u600mKk40M2", "Ener2", "Multiplicity==2");
	// t600mKth232->Project("u600mKth232M2", "Ener2", "Multiplicity==2");
	// t600mKu238->Project("u600mKu238M2", "Ener2", "Multiplicity==2");

	// tOVCco60->Project("uOVCco60M2", "Ener2", "Multiplicity==2");
	// tOVCk40->Project("uOVCk40M2", "Ener2", "Multiplicity==2");
	// tOVCth232->Project("uOVCth232M2", "Ener2", "Multiplicity==2");
	// tOVCu238->Project("uOVCu238M2", "Ener2", "Multiplicity==2");

	// tPbRomco60->Project("uPbRomco60M2", "Ener2", "Multiplicity==2");
	// tPbRomk40->Project("uPbRomk40M2", "Ener2", "Multiplicity==2");
	// tPbRomth232->Project("uPbRomth232M2", "Ener2", "Multiplicity==2");
	// tPbRomu238->Project("uPbRomu238M2", "Ener2", "Multiplicity==2");

	// tExtPbbi210->Project("uExtPbbi210M2", "Ener2", "Multiplicity==2");

	CombineAndSmear(tTeO20nu, uTeO20nuM1, hTeO20nuM1, 1, dSigma);
	CombineAndSmear(tTeO22nu, uTeO22nuM1, hTeO22nuM1, 1, dSigma);
	// CombineAndSmear(tTeO2co60, uTeO2co60M1, hTeO2co60M1, 1, dSigma);
	// CombineAndSmear(tTeO2k40, uTeO2k40M1, hTeO2k40M1, 1, dSigma);
	// CombineAndSmear(tTeO2po210, uTeO2po210M1, hTeO2po210M1, 1, dSigma);
	// CombineAndSmear(tTeO2ra228pb208, uTeO2ra228pb208M1, hTeO2ra228pb208M1, 1, dSigma);
	// CombineAndSmear(tTeO2th230only, uTeO2th230onlyM1, hTeO2th230onlyM1, 1, dSigma);
	// CombineAndSmear(tTeO2th232only, uTeO2th232onlyM1, hTeO2th232onlyM1, 1, dSigma);
	// CombineAndSmear(tCuFrameco60, uCuFrameco60M1, hCuFrameco60M1, 1, dSigma);
	// CombineAndSmear(tCuFrameth232, uCuFrameth232M1, hCuFrameth232M1, 1, dSigma);
	// CombineAndSmear(tCuFrameu238, uCuFrameu238M1, hCuFrameu238M1, 1, dSigma);
	// CombineAndSmear(tCuFramek40, uCuFramek40M1, hCuFramek40M1, 1, dSigma);
	// CombineAndSmear(t600mKco60, u600mKco60M1, h600mKco60M1, 1, dSigma);
	// CombineAndSmear(t600mKth232, u600mKth232M1, h600mKth232M1, 1, dSigma);
	// CombineAndSmear(t600mKu238, u600mKu238M1, h600mKu238M1, 1, dSigma);
	// CombineAndSmear(t600mKk40, u600mKk40M1, h600mKk40M1, 1, dSigma);
	// CombineAndSmear(tOVCco60, uOVCco60M1, hOVCco60M1, 1, dSigma);
	// CombineAndSmear(tOVCth232, uOVCth232M1, hOVCth232M1, 1, dSigma);
	// CombineAndSmear(tOVCu238, uOVCu238M1, hOVCu238M1, 1, dSigma);
	// CombineAndSmear(tOVCk40, uOVCk40M1, hOVCk40M1, 1, dSigma);
	// CombineAndSmear(tPbRomco60, uPbRomco60M1, hPbRomco60M1, 1, dSigma);
	// CombineAndSmear(tPbRomth232, uPbRomth232M1, hPbRomth232M1, 1, dSigma);
	// CombineAndSmear(tPbRomu238, uPbRomu238M1, hPbRomu238M1, 1, dSigma);
	// CombineAndSmear(tPbRomk40, uPbRomk40M1, hPbRomk40M1, 1, dSigma);
	// CombineAndSmear(tExtPbbi210, uExtPbbi210M1, hExtPbbi210M1, 1, dSigma);


	CombineAndSmear(tTeO20nu, uTeO20nuM2, hTeO20nuM2, 2, dSigma);
	CombineAndSmear(tTeO22nu, uTeO22nuM2, hTeO22nuM2, 2, dSigma);
	// CombineAndSmear(tTeO2co60, uTeO2co60M2, hTeO2co60M2, 2, dSigma);
	// CombineAndSmear(tTeO2k40, uTeO2k40M2, hTeO2k40M2, 2, dSigma);
	// CombineAndSmear(tTeO2po210, uTeO2po210M2, hTeO2po210M2, 2, dSigma);
	// CombineAndSmear(tTeO2ra228pb208, uTeO2ra228pb208M2, hTeO2ra228pb208M2, 2, dSigma);
	// CombineAndSmear(tTeO2th230only, uTeO2th230onlyM2, hTeO2th230onlyM2, 2, dSigma);
	// CombineAndSmear(tTeO2th232only, uTeO2th232onlyM2, hTeO2th232onlyM2, 2, dSigma);
	// CombineAndSmear(tCuFrameco60, uCuFrameco60M2, hCuFrameco60M2, 2, dSigma);
	// CombineAndSmear(tCuFrameth232, uCuFrameth232M2, hCuFrameth232M2, 2, dSigma);
	// CombineAndSmear(tCuFrameu238, uCuFrameu238M2, hCuFrameu238M2, 2, dSigma);
	// CombineAndSmear(tCuFramek40, uCuFramek40M2, hCuFramek40M2, 2, dSigma);
	// CombineAndSmear(t600mKco60, u600mKco60M2, h600mKco60M2, 2, dSigma);
	// CombineAndSmear(t600mKth232, u600mKth232M2, h600mKth232M2, 2, dSigma);
	// CombineAndSmear(t600mKu238, u600mKu238M2, h600mKu238M2, 2, dSigma);
	// CombineAndSmear(t600mKk40, u600mKk40M2, h600mKk40M2, 2, dSigma);
	// CombineAndSmear(tOVCco60, uOVCco60M2, hOVCco60M2, 2, dSigma);
	// CombineAndSmear(tOVCth232, uOVCth232M2, hOVCth232M2, 2, dSigma);
	// CombineAndSmear(tOVCu238, uOVCu238M2, hOVCu238M2, 2, dSigma);
	// CombineAndSmear(tOVCk40, uOVCk40M2, hOVCk40M2, 2, dSigma);
	// CombineAndSmear(tPbRomco60, uPbRomco60M2, hPbRomco60M2, 2, dSigma);
	// CombineAndSmear(tPbRomth232, uPbRomth232M2, hPbRomth232M2, 2, dSigma);
	// CombineAndSmear(tPbRomu238, uPbRomu238M2, hPbRomu238M2, 2, dSigma);
	// CombineAndSmear(tPbRomk40, uPbRomk40M2, hPbRomk40M2, 2, dSigma);
	// CombineAndSmear(tExtPbbi210, uExtPbbi210M2, hExtPbbi210M2, 2, dSigma);


	NormalizePDFPair( hTeO20nuM1, hTeO20nuM2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO22nuM1, hTeO22nuM2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hTeO2co60M1, hTeO2co60M2,  dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hTeO2k40M1, hTeO2k40M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hTeO2po210M1, hTeO2po210M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hTeO2ra228pb208M1, hTeO2ra228pb208M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hTeO2th230onlyM1, hTeO2th230onlyM2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hTeO2th232onlyM1, hTeO2th232onlyM2, dMinEnergy, dMaxEnergy);

	// NormalizePDFPair( hCuFrameco60M1, hCuFrameco60M2,  dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hCuFramek40M1, hCuFramek40M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hCuFrameth232M1, hCuFrameth232M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hCuFrameu238M1, hCuFrameu238M2, dMinEnergy, dMaxEnergy);

	// NormalizePDFPair( h600mKco60M1, h600mKco60M2,  dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( h600mKk40M1, h600mKk40M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( h600mKth232M1, h600mKth232M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( h600mKu238M1, h600mKu238M2, dMinEnergy, dMaxEnergy);

	// NormalizePDFPair( hOVCco60M1, hOVCco60M2,  dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hOVCk40M1, hOVCk40M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hOVCth232M1, hOVCth232M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hOVCu238M1, hOVCu238M2, dMinEnergy, dMaxEnergy);

	// NormalizePDFPair( hPbRomco60M1, hPbRomco60M2,  dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hPbRomk40M1, hPbRomk40M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hPbRomth232M1, hPbRomth232M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hPbRomu238M1, hPbRomu238M2, dMinEnergy, dMaxEnergy);

	TFile *file1 = new TFile("MCProduction_ReducedBulk.root", "RECREATE");
	hTeO20nuM1->Write();
	hTeO22nuM1->Write();
	hTeO2co60M1->Write();
	hTeO2k40M1->Write();
	hTeO2po210M1->Write();
	hTeO2ra228pb208M1->Write();
	hTeO2th230onlyM1->Write();
	hTeO2th232onlyM1->Write();

	hCuFrameco60M1->Write();
	hCuFramek40M1->Write();
	hCuFrameth232M1->Write();
	hCuFrameu238M1->Write();

	h600mKco60M1->Write();
	h600mKk40M1->Write();
	h600mKth232M1->Write();
	h600mKu238M1->Write();

	hPbRomco60M1->Write();
	hPbRomk40M1->Write();
	hPbRomth232M1->Write();
	hPbRomu238M1->Write();

	hOVCco60M1->Write();
	hOVCk40M1->Write();
	hOVCth232M1->Write();
	hOVCu238M1->Write();

	hExtPbbi210M1->Write();


	hTeO20nuM2->Write();
	hTeO22nuM2->Write();
	hTeO2co60M2->Write();
	hTeO2k40M2->Write();
	hTeO2po210M2->Write();
	hTeO2ra228pb208M2->Write();
	hTeO2th230onlyM2->Write();
	hTeO2th232onlyM2->Write();

	hCuFrameco60M2->Write();
	hCuFramek40M2->Write();
	hCuFrameth232M2->Write();
	hCuFrameu238M2->Write();

	h600mKco60M2->Write();
	h600mKk40M2->Write();
	h600mKth232M2->Write();
	h600mKu238M2->Write();

	hPbRomco60M2->Write();
	hPbRomk40M2->Write();
	hPbRomth232M2->Write();
	hPbRomu238M2->Write();

	hOVCco60M2->Write();
	hOVCk40M2->Write();
	hOVCth232M2->Write();
	hOVCu238M2->Write();

	hExtPbbi210M2->Write();

	file1->Write();

}

void SaveHistogramsReducedSurface()
{
	std::string sDataDir = "/cuore/user/zhubrian/MC/Bkg/Unsmeared/";

	double dMinEnergy = 0;
	double dMaxEnergy = 10000;
	int dBinSize = 1;
	int dNBins = (dMaxEnergy - dMinEnergy)/dBinSize;

	// Unsmeared
	TH1D *uTeO2Sxpb210M1_001 = new TH1D("uTeO2Sxpb210M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxpb210M1_1 = new TH1D("uTeO2Sxpb210M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxth230onlyM1_001 = new TH1D("uTeO2Sxth230onlyM1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxth232onlyM1_001 = new TH1D("uTeO2Sxth232onlyM1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxra228pb208M1_001 = new TH1D("uTeO2Sxra228pb208M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxu238th230M1_001 = new TH1D("uTeO2Sxu238th230M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxra226pb210M1_001 = new TH1D("uTeO2Sxra226pb210M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxth232M1_001 = new TH1D("uTeO2Sxth232M1_001", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *uCuBoxSxpb210M1_01 = new TH1D("uCuBoxSxpb210M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uCuBoxSxpb210M1_10 = new TH1D("uCuBoxSxpb210M1_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uCuBoxSxth232M1_10 = new TH1D("uCuBoxSxth232M1_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uCuBoxSxu238M1_10 = new TH1D("uCuBoxSxu238M1_10", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *uTeO2Sxpb210M2_001 = new TH1D("uTeO2Sxpb210M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxpb210M2_1 = new TH1D("uTeO2Sxpb210M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxth230onlyM2_001 = new TH1D("uTeO2Sxth230onlyM2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxth232onlyM2_001 = new TH1D("uTeO2Sxth232onlyM2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxra228pb208M2_001 = new TH1D("uTeO2Sxra228pb208M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxu238th230M2_001 = new TH1D("uTeO2Sxu238th230M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxra226pb210M2_001 = new TH1D("uTeO2Sxra226pb210M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uTeO2Sxth232M2_001 = new TH1D("uTeO2Sxth232M2_001", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *uCuBoxSxpb210M2_01 = new TH1D("uCuBoxSxpb210M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uCuBoxSxpb210M2_10 = new TH1D("uCuBoxSxpb210M2_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uCuBoxSxth232M2_10 = new TH1D("uCuBoxSxth232M2_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *uCuBoxSxu238M2_10 = new TH1D("uCuBoxSxu238M2_10", "", dNBins, dMinEnergy, dMaxEnergy);


	// Fill these
	TH1D *hTeO2Sxpb210M1_001 = new TH1D("hTeO2Sxpb210M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M1_1 = new TH1D("hTeO2Sxpb210M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth230onlyM1_001 = new TH1D("hTeO2Sxth230onlyM1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth232onlyM1_001 = new TH1D("hTeO2Sxth232onlyM1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxra228pb208M1_001 = new TH1D("hTeO2Sxra228pb208M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238th230M1_001 = new TH1D("hTeO2Sxu238th230M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxra226pb210M1_001 = new TH1D("hTeO2Sxra226pb210M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth232M1_001 = new TH1D("hTeO2Sxth232M1_001", "", dNBins, dMinEnergy, dMaxEnergy);


	TH1D *hCuBoxSxpb210M1_01 = new TH1D("hCuBoxSxpb210M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxpb210M1_10 = new TH1D("hCuBoxSxpb210M1_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M1_10 = new TH1D("hCuBoxSxth232M1_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M1_10 = new TH1D("hCuBoxSxu238M1_10", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hTeO2Sxpb210M2_001 = new TH1D("hTeO2Sxpb210M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M2_1 = new TH1D("hTeO2Sxpb210M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth230onlyM2_001 = new TH1D("hTeO2Sxth230onlyM2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth232onlyM2_001 = new TH1D("hTeO2Sxth232onlyM2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxra228pb208M2_001 = new TH1D("hTeO2Sxra228pb208M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238th230M2_001 = new TH1D("hTeO2Sxu238th230M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxra226pb210M2_001 = new TH1D("hTeO2Sxra226pb210M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth232M2_001 = new TH1D("hTeO2Sxth232M2_001", "", dNBins, dMinEnergy, dMaxEnergy);


	TH1D *hCuBoxSxpb210M2_01 = new TH1D("hCuBoxSxpb210M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxpb210M2_10 = new TH1D("hCuBoxSxpb210M2_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M2_10 = new TH1D("hCuBoxSxth232M2_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M2_10 = new TH1D("hCuBoxSxu238M2_10", "", dNBins, dMinEnergy, dMaxEnergy);


	TChain *tTeO2Sxpb210_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "pb210-.01");
	TChain *tTeO2Sxpb210_1 = LoadMC(sDataDir.c_str(), "TeO2Sx", "pb210-1");
	TChain *tTeO2Sxth232only_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "th232_only-.01");
	TChain *tTeO2Sxra228pb208_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "ra228-pb208-.01");
	TChain *tTeO2Sxu238th230_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "u238-th230-.01");
	TChain *tTeO2Sxth230only_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "th230_only-.01");
	TChain *tTeO2Sxra226pb210_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "ra226-pb210-.01");
	TChain *tTeO2Sxth232_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "th232-.01");

	TChain *tCuBoxSxpb210_01 = LoadMC(sDataDir.c_str(), "CuBoxSx", "pb210-.1");
	TChain *tCuBoxSxpb210_10 = LoadMC(sDataDir.c_str(), "CuBoxSx", "pb210-10");
	TChain *tCuBoxSxth232_10 = LoadMC(sDataDir.c_str(), "CuBoxSx", "th232-10");
	TChain *tCuBoxSxu238_10 = LoadMC(sDataDir.c_str(), "CuBoxSx", "u238-10");

	tTeO2Sxpb210_001->Project("uTeO2Sxpb210M1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxpb210_1->Project("uTeO2Sxpb210M1_1", "Ener2", "Multiplicity==1");
	tTeO2Sxth232only_001->Project("uTeO2Sxth232onlyM1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxra228pb208_001->Project("uTeO2Sxra228pb208M1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxu238th230_001->Project("uTeO2Sxu238th230M1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxth230only_001->Project("uTeO2Sxth230onlyM1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxra226pb210_001->Project("uTeO2Sxra226pb210M1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxth232_001->Project("uTeO2Sxth232M1_001", "Ener2", "Multiplicity==1");

	tCuBoxSxpb210_01->Project("uCuBoxSxpb210M1_01", "Ener2", "Multiplicity==1");
	tCuBoxSxpb210_10->Project("uCuBoxSxpb210M1_10", "Ener2", "Multiplicity==1");
	tCuBoxSxth232_10->Project("uCuBoxSxth232M1_10", "Ener2", "Multiplicity==1");
	tCuBoxSxu238_10->Project("uCuBoxSxu238M1_10", "Ener2", "Multiplicity==1");


	tTeO2Sxpb210_001->Project("uTeO2Sxpb210M2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxpb210_1->Project("uTeO2Sxpb210M2_1", "Ener2", "Multiplicity==2");
	tTeO2Sxth232only_001->Project("uTeO2Sxth232onlyM2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxra228pb208_001->Project("uTeO2Sxra228pb208M2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxu238th230_001->Project("uTeO2Sxu238th230M2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxth230only_001->Project("uTeO2Sxth230onlyM2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxra226pb210_001->Project("uTeO2Sxra226pb210M2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxth232_001->Project("uTeO2Sxth232M2_001", "Ener2", "Multiplicity==2");

	tCuBoxSxpb210_01->Project("uCuBoxSxpb210M2_01", "Ener2", "Multiplicity==2");
	tCuBoxSxpb210_10->Project("uCuBoxSxpb210M2_10", "Ener2", "Multiplicity==2");
	tCuBoxSxth232_10->Project("uCuBoxSxth232M2_10", "Ener2", "Multiplicity==2");
	tCuBoxSxu238_10->Project("uCuBoxSxu238M2_10", "Ener2", "Multiplicity==2");




	SmearMC(uTeO2Sxpb210M1_001, hTeO2Sxpb210M1_001);
	SmearMC(uTeO2Sxpb210M1_1, hTeO2Sxpb210M1_1);
	SmearMC(uTeO2Sxth232onlyM1_001, hTeO2Sxth232onlyM1_001);
	SmearMC(uTeO2Sxra228pb208M1_001, hTeO2Sxra228pb208M1_001);
	SmearMC(uTeO2Sxu238th230M1_001, hTeO2Sxu238th230M1_001);
	SmearMC(uTeO2Sxth230onlyM1_001, hTeO2Sxth230onlyM1_001);
	SmearMC(uTeO2Sxra226pb210M1_001, hTeO2Sxra226pb210M1_001);
	SmearMC(uTeO2Sxth232M1_001, hTeO2Sxth232M1_001);

	SmearMC(uCuBoxSxpb210M1_01, hCuBoxSxpb210M1_01);
	SmearMC(uCuBoxSxpb210M1_10, hCuBoxSxpb210M1_10);
	SmearMC(uCuBoxSxth232M1_10, hCuBoxSxth232M1_10);
	SmearMC(uCuBoxSxu238M1_10, hCuBoxSxu238M1_10);

	SmearMC(uTeO2Sxpb210M2_001, hTeO2Sxpb210M2_001);
	SmearMC(uTeO2Sxpb210M2_1, hTeO2Sxpb210M2_1);
	SmearMC(uTeO2Sxth232onlyM2_001, hTeO2Sxth232onlyM2_001);
	SmearMC(uTeO2Sxra228pb208M2_001, hTeO2Sxra228pb208M2_001);
	SmearMC(uTeO2Sxu238th230M2_001, hTeO2Sxu238th230M2_001);
	SmearMC(uTeO2Sxth230onlyM2_001, hTeO2Sxth230onlyM2_001);
	SmearMC(uTeO2Sxra226pb210M2_001, hTeO2Sxra226pb210M2_001);
	SmearMC(uTeO2Sxth232M2_001, hTeO2Sxth232M2_001);

	SmearMC(uCuBoxSxpb210M2_01, hCuBoxSxpb210M2_01);
	SmearMC(uCuBoxSxpb210M2_10, hCuBoxSxpb210M2_10);
	SmearMC(uCuBoxSxth232M2_10, hCuBoxSxth232M2_10);
	SmearMC(uCuBoxSxu238M2_10, hCuBoxSxu238M2_10);


	NormalizePDFPair(hTeO2Sxpb210M1_001, hTeO2Sxpb210M2_001, dMinEnergy, dMaxEnergy);
	NormalizePDFPair(hTeO2Sxpb210M1_1, hTeO2Sxpb210M2_1, dMinEnergy, dMaxEnergy);
	NormalizePDFPair(hTeO2Sxth232onlyM1_001, hTeO2Sxth232onlyM2_001, dMinEnergy, dMaxEnergy);
	NormalizePDFPair(hTeO2Sxra228pb208M1_001, hTeO2Sxra228pb208M2_001, dMinEnergy, dMaxEnergy);
	NormalizePDFPair(hTeO2Sxu238th230M1_001, hTeO2Sxu238th230M2_001, dMinEnergy, dMaxEnergy);
	NormalizePDFPair(hTeO2Sxth230onlyM1_001, hTeO2Sxth230onlyM2_001, dMinEnergy, dMaxEnergy);
	NormalizePDFPair(hTeO2Sxra226pb210M1_001, hTeO2Sxra226pb210M2_001, dMinEnergy, dMaxEnergy);
	NormalizePDFPair(hTeO2Sxth232M1_001, hTeO2Sxth232M2_001, dMinEnergy, dMaxEnergy);

	NormalizePDFPair(hCuBoxSxpb210M1_01, hCuBoxSxpb210M2_01, dMinEnergy, dMaxEnergy);
	NormalizePDFPair(hCuBoxSxpb210M1_10, hCuBoxSxpb210M2_10, dMinEnergy, dMaxEnergy);
	NormalizePDFPair(hCuBoxSxth232M1_10, hCuBoxSxth232M2_10, dMinEnergy, dMaxEnergy);
	NormalizePDFPair(hCuBoxSxu238M1_10, hCuBoxSxu238M2_10, dMinEnergy, dMaxEnergy);


	TFile *file1 = new TFile("MCProduction_ReducedSurface.root", "RECREATE");

	hTeO2Sxpb210M1_001->Write();
	hTeO2Sxpb210M1_1->Write();
	hTeO2Sxth232onlyM1_001->Write();
	hTeO2Sxra228pb208M1_001->Write();
	hTeO2Sxu238th230M1_001->Write();
	hTeO2Sxth230onlyM1_001->Write();
	hTeO2Sxra226pb210M1_001->Write();
	hTeO2Sxth232M1_001->Write();

	hCuBoxSxpb210M1_01->Write();
	hCuBoxSxpb210M1_10->Write();
	hCuBoxSxth232M1_10->Write();
	hCuBoxSxu238M1_10->Write();

	hTeO2Sxpb210M2_001->Write();
	hTeO2Sxpb210M2_1->Write();
	hTeO2Sxth232onlyM2_001->Write();
	hTeO2Sxra228pb208M2_001->Write();
	hTeO2Sxu238th230M2_001->Write();
	hTeO2Sxth230onlyM2_001->Write();
	hTeO2Sxra226pb210M2_001->Write();
	hTeO2Sxth232M2_001->Write();

	hCuBoxSxpb210M2_01->Write();
	hCuBoxSxpb210M2_10->Write();
	hCuBoxSxth232M2_10->Write();
	hCuBoxSxu238M2_10->Write();

}


void SaveHistogramsFudge()
{
	std::string sDataDir = "/cuore/data/simulation/CUORE0/t14.08/production_g2root-r363/ntp/Bulk/";

	double dMinEnergy = 0;
	double dMaxEnergy = 10000;
	int dBinSize = 1;
	int dNBins = (dMaxEnergy - dMinEnergy)/dBinSize;

	TH1D *h50mKk40M1 = new TH1D("h50mKk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKk40M1 = new TH1D("h600mKk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCk40M1 = new TH1D("hIVCk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalk40M1 = new TH1D("hInternalk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCk40M1 = new TH1D("hOVCk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomk40M1 = new TH1D("hPbRomk40M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hFudge661M1 = new TH1D("hFudge661M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hFudge803M1 = new TH1D("hFudge803M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hFudge1063M1 = new TH1D("hFudge1063M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *h50mKk40M2 = new TH1D("h50mKk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKk40M2 = new TH1D("h600mKk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCk40M2 = new TH1D("hIVCk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalk40M2 = new TH1D("hInternalk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCk40M2 = new TH1D("hOVCk40M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hFudge661M2 = new TH1D("hFudge661M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hFudge803M2 = new TH1D("hFudge803M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hFudge1063M2 = new TH1D("hFudge1063M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TChain *t50mKk40 = LoadMC(sDataDir.c_str(), "50mK", "k40");
	TChain *t600mKk40 = LoadMC(sDataDir.c_str(), "600mK", "k40");
	TChain *tIVCk40 = LoadMC(sDataDir.c_str(), "IVC", "k40");
	TChain *tOVCk40 = LoadMC(sDataDir.c_str(), "OVC", "k40");
	TChain *tPbRomk40 = LoadMC(sDataDir.c_str(), "PbRom", "k40");

	t50mKk40->Project("h50mKk40M1", "Ener2", "Multiplicity==1");
	t600mKk40->Project("h600mKk40M1", "Ener2", "Multiplicity==1");
	tIVCk40->Project("hIVCk40M1", "Ener2", "Multiplicity==1");
	tOVCk40->Project("hOVCk40M1", "Ener2", "Multiplicity==1");
	tPbRomk40->Project("hPbRomk40M1", "Ener2", "Multiplicity==1");

	t50mKk40->Project("h50mKk40M2", "Ener2", "Multiplicity==2");
	t600mKk40->Project("h600mKk40M2", "Ener2", "Multiplicity==2");
	tIVCk40->Project("hIVCk40M2", "Ener2", "Multiplicity==2");
	tOVCk40->Project("hOVCk40M2", "Ener2", "Multiplicity==2");
	tPbRomk40->Project("hPbRomk40M2", "Ener2", "Multiplicity==2");

	hInternalk40M1->Add(h600mKk40M1);
	hInternalk40M1->Add(h50mKk40M1, 0.8625620249);
	hInternalk40M1->Add(hIVCk40M1, 2.3168024764);	

	hInternalk40M2->Add(h600mKk40M2);
	hInternalk40M2->Add(h50mKk40M2, 0.8625620249);
	hInternalk40M2->Add(hIVCk40M2, 2.3168024764);		

	// for (int i = 1; i <= dNBins; i++)
	// {
	// 	// hFudge661M1->SetBinContent(i, hInternalk40M1->GetBinContent(0.4527397*i) );

	// }
	for (int i = 1; i < 10000000; i++)
	{
		hFudge661M1->Fill(hInternalk40M1->GetRandom()*0.4527397);
		hFudge803M1->Fill(hInternalk40M1->GetRandom()*0.55);
		hFudge1063M1->Fill(hInternalk40M1->GetRandom()*0.72808);

		hFudge661M2->Fill(hInternalk40M2->GetRandom()*0.4527397);
		hFudge803M2->Fill(hInternalk40M2->GetRandom()*0.55);
		hFudge1063M2->Fill(hInternalk40M2->GetRandom()*0.72808);
	}

	double dIntegral1 = hInternalk40M1->Integral(1, dNBins);
	double dEntries1 = hInternalk40M1->GetEntries();
	double dIntegral2 = hInternalk40M2->Integral(1, dNBins);
	double dEntries2 = hInternalk40M2->GetEntries();

	hFudge661M1->Scale(1/hFudge661M1->Integral(1, dNBins));
	hFudge661M2->Scale(1/hFudge661M2->Integral(1, dNBins) * dIntegral2/dIntegral1);

	hFudge803M1->Scale(1/hFudge803M1->Integral(1, dNBins));
	hFudge803M2->Scale(1/hFudge803M2->Integral(1, dNBins) * dIntegral2/dIntegral1);

	hFudge1063M1->Scale(1/hFudge1063M1->Integral(1, dNBins));
	hFudge1063M2->Scale(1/hFudge1063M2->Integral(1, dNBins) * dIntegral2/dIntegral1);

	TFile *file1 = new TFile("MCProduction_FudgeFactor_1keV.root", "RECREATE");
	hFudge661M1->Write();
	hFudge661M2->Write();

	hFudge803M1->Write();
	hFudge803M2->Write();

	hFudge1063M1->Write();
	hFudge1063M2->Write();
	file1->Write();
}


void SaveHistogramsBulkInner()
{
	std::string sDataDir = "/cuore/user/zhubrian/MC/Bkg/Production/";

	double dMinEnergy = 0;
	double dMaxEnergy = 10000;
	int dBinSize = 1;
	int dNBins = (dMaxEnergy - dMinEnergy)/dBinSize;

	TH1D *hCuBoxco58M1 = new TH1D("hCuBoxco58M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxco60M1 = new TH1D("hCuBoxco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxcs137M1 = new TH1D("hCuBoxcs137M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxk40M1 = new TH1D("hCuBoxk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxmn54M1 = new TH1D("hCuBoxmn54M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxpb210M1 = new TH1D("hCuBoxpb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxth232M1 = new TH1D("hCuBoxth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxu238M1 = new TH1D("hCuBoxu238M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxeu152M1 = new TH1D("hCuBoxeu152M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxeu154M1 = new TH1D("hCuBoxeu154M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuFrameco58M1 = new TH1D("hCuFrameco58M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameco60M1 = new TH1D("hCuFrameco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramecs137M1 = new TH1D("hCuFramecs137M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramek40M1 = new TH1D("hCuFramek40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramemn54M1 = new TH1D("hCuFramemn54M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramepb210M1 = new TH1D("hCuFramepb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameth232M1 = new TH1D("hCuFrameth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameu238M1 = new TH1D("hCuFrameu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuBox_CuFrameco60M1 = new TH1D("hCuBox_CuFrameco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFramek40M1 = new TH1D("hCuBox_CuFramek40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M1 = new TH1D("hCuBox_CuFrameth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M1 = new TH1D("hCuBox_CuFrameu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	// TH1D *hPTFEth232M1 = new TH1D("hPTFEth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hPTFEu238M1 = new TH1D("hPTFEu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hTeO20nuM1 = new TH1D("hTeO20nuM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO22nuM1 = new TH1D("hTeO22nuM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2co60M1 = new TH1D("hTeO2co60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2k40M1 = new TH1D("hTeO2k40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2pb210M1 = new TH1D("hTeO2pb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2po210M1 = new TH1D("hTeO2po210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226onlyM1 = new TH1D("hTeO2ra226onlyM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226pb206M1 = new TH1D("hTeO2ra226pb206M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226pb210M1 = new TH1D("hTeO2ra226pb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra228pb208M1 = new TH1D("hTeO2ra228pb208M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2rn222onlyM1 = new TH1D("hTeO2rn222onlyM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2rn222pb206M1 = new TH1D("hTeO2rn222pb206M1", "", dNBins, dMinEnergy, dMaxEnergy);	
	TH1D *hTeO2sb125M1 = new TH1D("hTeO2sb125M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2te125M1 = new TH1D("hTeO2te125M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th230onlyM1 = new TH1D("hTeO2th230onlyM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th230pb206M1 = new TH1D("hTeO2th230pb206M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232onlyM1 = new TH1D("hTeO2th232onlyM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232M1 = new TH1D("hTeO2th232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232th228M1 = new TH1D("hTeO2th232th228M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238ra226M1 = new TH1D("hTeO2u238ra226M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238rn222M1 = new TH1D("hTeO2u238rn222M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238M1 = new TH1D("hTeO2u238M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238th230M1 = new TH1D("hTeO2u238th230M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238u234M1 = new TH1D("hTeO2u238u234M1", "", dNBins, dMinEnergy, dMaxEnergy);



////////////////////////

	TH1D *hCuBoxco58M2 = new TH1D("hCuBoxco58M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxco60M2 = new TH1D("hCuBoxco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxcs137M2 = new TH1D("hCuBoxcs137M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxk40M2 = new TH1D("hCuBoxk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxmn54M2 = new TH1D("hCuBoxmn54M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxpb210M2 = new TH1D("hCuBoxpb210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxth232M2 = new TH1D("hCuBoxth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxu238M2 = new TH1D("hCuBoxu238M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxeu152M2 = new TH1D("hCuBoxeu152M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxeu154M2 = new TH1D("hCuBoxeu154M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuFrameco58M2 = new TH1D("hCuFrameco58M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameco60M2 = new TH1D("hCuFrameco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramecs137M2 = new TH1D("hCuFramecs137M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramek40M2 = new TH1D("hCuFramek40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramemn54M2 = new TH1D("hCuFramemn54M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramepb210M2 = new TH1D("hCuFramepb210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameth232M2 = new TH1D("hCuFrameth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameu238M2 = new TH1D("hCuFrameu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuBox_CuFrameco60M2 = new TH1D("hCuBox_CuFrameco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFramek40M2 = new TH1D("hCuBox_CuFramek40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M2 = new TH1D("hCuBox_CuFrameth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M2 = new TH1D("hCuBox_CuFrameu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	// TH1D *hPTFEth232M2 = new TH1D("hPTFEth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hPTFEu238M2 = new TH1D("hPTFEu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hTeO20nuM2 = new TH1D("hTeO20nuM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO22nuM2 = new TH1D("hTeO22nuM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2co60M2 = new TH1D("hTeO2co60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2k40M2 = new TH1D("hTeO2k40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2pb210M2 = new TH1D("hTeO2pb210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2po210M2 = new TH1D("hTeO2po210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226onlyM2 = new TH1D("hTeO2ra226onlyM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226pb206M2 = new TH1D("hTeO2ra226pb206M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226pb210M2 = new TH1D("hTeO2ra226pb210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra228pb208M2 = new TH1D("hTeO2ra228pb208M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2rn222onlyM2 = new TH1D("hTeO2rn222onlyM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2rn222pb206M2 = new TH1D("hTeO2rn222pb206M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2sb125M2 = new TH1D("hTeO2sb125M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2te125M2 = new TH1D("hTeO2te125M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th230onlyM2 = new TH1D("hTeO2th230onlyM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th230pb206M2 = new TH1D("hTeO2th230pb206M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232onlyM2 = new TH1D("hTeO2th232onlyM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232M2 = new TH1D("hTeO2th232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232th228M2 = new TH1D("hTeO2th232th228M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238ra226M2 = new TH1D("hTeO2u238ra226M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238rn222M2 = new TH1D("hTeO2u238rn222M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238M2 = new TH1D("hTeO2u238M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238th230M2 = new TH1D("hTeO2u238th230M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238u234M2 = new TH1D("hTeO2u238u234M2", "", dNBins, dMinEnergy, dMaxEnergy);



///////////////////////////////

	TChain *tCuBoxco58 = LoadMC(sDataDir.c_str(), "CuBox", "co58");
	TChain *tCuBoxco60 = LoadMC(sDataDir.c_str(), "CuBox", "co60");
	TChain *tCuBoxcs137 = LoadMC(sDataDir.c_str(), "CuBox", "cs137");
	TChain *tCuBoxk40 = LoadMC(sDataDir.c_str(), "CuBox", "k40");
	TChain *tCuBoxmn54 = LoadMC(sDataDir.c_str(), "CuBox", "mn54");
	TChain *tCuBoxpb210 = LoadMC(sDataDir.c_str(), "CuBox", "pb210");
	TChain *tCuBoxth232 = LoadMC(sDataDir.c_str(), "CuBox", "th232");
	TChain *tCuBoxu238 = LoadMC(sDataDir.c_str(), "CuBox", "u238");
	TChain *tCuBoxeu152 = LoadMC(sDataDir.c_str(), "CuBox", "eu152");
	TChain *tCuBoxeu154 = LoadMC(sDataDir.c_str(), "CuBox", "eu154");

	TChain *tCuFrameco58 = LoadMC(sDataDir.c_str(), "CuFrame", "co58");
	TChain *tCuFrameco60 = LoadMC(sDataDir.c_str(), "CuFrame", "co60");
	TChain *tCuFramecs137 = LoadMC(sDataDir.c_str(), "CuFrame", "cs137");
	TChain *tCuFramek40 = LoadMC(sDataDir.c_str(), "CuFrame", "k40");
	TChain *tCuFramemn54 = LoadMC(sDataDir.c_str(), "CuFrame", "mn54");
	TChain *tCuFramepb210 = LoadMC(sDataDir.c_str(), "CuFrame", "pb210");
	TChain *tCuFrameth232 = LoadMC(sDataDir.c_str(), "CuFrame", "th232");
	TChain *tCuFrameu238 = LoadMC(sDataDir.c_str(), "CuFrame", "u238");

	// TChain *tCuBox_CuFrameco60 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "co60");
	// TChain *tCuBox_CuFramek40 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "k40");
	// TChain *tCuBox_CuFrameth232 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "th232");
	// TChain *tCuBox_CuFrameu238 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "u238");

	// TChain *tInternalco60 = LoadMC(sDataDir.c_str(), "InternalShields", "co60");
	// TChain *tInternalk40 = LoadMC(sDataDir.c_str(), "InternalShields", "k40");
	// TChain *tInternalth232 = LoadMC(sDataDir.c_str(), "InternalShields", "th232");
	// TChain *tInternalu238 = LoadMC(sDataDir.c_str(), "InternalShields", "u238");

	// TChain *tPTFEth232 = LoadMC(sDataDir.c_str(), "PTFE", "th232");
	// TChain *tPTFEu238 = LoadMC(sDataDir.c_str(), "PTFE", "u238");

	TChain *tTeO20nu = LoadMC(sDataDir.c_str(), "TeO2", "0nu_1um");
	TChain *tTeO22nu = LoadMC(sDataDir.c_str(), "TeO2", "2nu");
	TChain *tTeO2co60 = LoadMC(sDataDir.c_str(), "TeO2", "co60");
	TChain *tTeO2k40 = LoadMC(sDataDir.c_str(), "TeO2", "k40");
	TChain *tTeO2pb210 = LoadMC(sDataDir.c_str(), "TeO2", "pb210");
	TChain *tTeO2po210 = LoadMC(sDataDir.c_str(), "TeO2", "po210");
	TChain *tTeO2ra226only = LoadMC(sDataDir.c_str(), "TeO2", "ra226_only");
	TChain *tTeO2ra226pb206 = LoadMC(sDataDir.c_str(), "TeO2", "ra226-pb206");
	TChain *tTeO2ra226pb210 = LoadMC(sDataDir.c_str(), "TeO2", "ra226-pb210");
	TChain *tTeO2ra228pb208 = LoadMC(sDataDir.c_str(), "TeO2", "ra228-pb208");
	TChain *tTeO2rn222only = LoadMC(sDataDir.c_str(), "TeO2", "rn222_only");
	TChain *tTeO2rn222pb206 = LoadMC(sDataDir.c_str(), "TeO2", "rn222-pb206");
	TChain *tTeO2sb125 = LoadMC(sDataDir.c_str(), "TeO2", "sb125");
	TChain *tTeO2te125 = LoadMC(sDataDir.c_str(), "TeO2", "te125m");
	TChain *tTeO2th230only = LoadMC(sDataDir.c_str(), "TeO2", "th230_only");
	TChain *tTeO2th230pb206 = LoadMC(sDataDir.c_str(), "TeO2", "th230-pb206");
	TChain *tTeO2th232only = LoadMC(sDataDir.c_str(), "TeO2", "th232_only");
	TChain *tTeO2th232 = LoadMC(sDataDir.c_str(), "TeO2", "th232");
	TChain *tTeO2th232th228 = LoadMC(sDataDir.c_str(), "TeO2", "th232-th228");
	TChain *tTeO2u238ra226 = LoadMC(sDataDir.c_str(), "TeO2", "u238-ra226");
	TChain *tTeO2u238rn222 = LoadMC(sDataDir.c_str(), "TeO2", "u238-rn222");
	TChain *tTeO2u238 = LoadMC(sDataDir.c_str(), "TeO2", "u238");
	TChain *tTeO2u238th230 = LoadMC(sDataDir.c_str(), "TeO2", "u238-th230");
	TChain *tTeO2u238u234 = LoadMC(sDataDir.c_str(), "TeO2", "u238-u234");

////////////////////////////////////////


	tCuBoxco58->Project("hCuBoxco58M1", "Ener2", "Multiplicity==1");
	tCuBoxco60->Project("hCuBoxco60M1", "Ener2", "Multiplicity==1");
	tCuBoxcs137->Project("hCuBoxcs137M1", "Ener2", "Multiplicity==1");
	tCuBoxk40->Project("hCuBoxk40M1", "Ener2", "Multiplicity==1");
	tCuBoxmn54->Project("hCuBoxmn54M1", "Ener2", "Multiplicity==1");
	tCuBoxpb210->Project("hCuBoxpb210M1", "Ener2", "Multiplicity==1");
	tCuBoxth232->Project("hCuBoxth232M1", "Ener2", "Multiplicity==1");
	tCuBoxu238->Project("hCuBoxu238M1", "Ener2", "Multiplicity==1");	
	tCuBoxeu152->Project("hCuBoxeu152M1", "Ener2", "Multiplicity==1");
	tCuBoxeu154->Project("hCuBoxeu154M1", "Ener2", "Multiplicity==1");

	tCuFrameco58->Project("hCuFrameco58M1", "Ener2", "Multiplicity==1");
	tCuFrameco60->Project("hCuFrameco60M1", "Ener2", "Multiplicity==1");
	tCuFramecs137->Project("hCuFramecs137M1", "Ener2", "Multiplicity==1");
	tCuFramek40->Project("hCuFramek40M1", "Ener2", "Multiplicity==1");
	tCuFramemn54->Project("hCuFramemn54M1", "Ener2", "Multiplicity==1");
	tCuFramepb210->Project("hCuFramepb210M1", "Ener2", "Multiplicity==1");
	tCuFrameth232->Project("hCuFrameth232M1", "Ener2", "Multiplicity==1");
	tCuFrameu238->Project("hCuFrameu238M1", "Ener2", "Multiplicity==1");

	// tCuBox_CuFrameco60->Project("hCuBox_CuFrameco60M1", "Ener2", "Multiplicity==1");
	// tCuBox_CuFramek40->Project("hCuBox_CuFramek40M1", "Ener2", "Multiplicity==1");
	// tCuBox_CuFrameth232->Project("hCuBox_CuFrameth232M1", "Ener2", "Multiplicity==1");
	// tCuBox_CuFrameu238->Project("hCuBox_CuFrameu238M1", "Ener2", "Multiplicity==1");

	// tInternalco60->Project("hInternalco60M1", "Ener2", "Multiplicity==1");
	// tInternalk40->Project("hInternalk40M1", "Ener2", "Multiplicity==1");
	// tInternalth232->Project("hInternalth232M1", "Ener2", "Multiplicity==1");
	// tInternalu238->Project("hInternalu238M1", "Ener2", "Multiplicity==1");

	// tPTFEth232->Project("hPTFEth232M1", "Ener2", "Multiplicity==1");
	// tPTFEu238->Project("hPTFEu238M1", "Ener2", "Multiplicity==1");


	tTeO20nu->Project("hTeO20nuM1", "Ener2", "Multiplicity==1");
	tTeO22nu->Project("hTeO22nuM1", "Ener2", "Multiplicity==1");
	tTeO2co60->Project("hTeO2co60M1", "Ener2", "Multiplicity==1");
	tTeO2k40->Project("hTeO2k40M1", "Ener2", "Multiplicity==1");
	tTeO2pb210->Project("hTeO2pb210M1", "Ener2", "Multiplicity==1");
	tTeO2po210->Project("hTeO2po210M1", "Ener2", "Multiplicity==1");
	tTeO2ra226only->Project("hTeO2ra226onlyM1", "Ener2", "Multiplicity==1");
	tTeO2ra226pb206->Project("hTeO2ra226pb206M1", "Ener2", "Multiplicity==1");
	tTeO2ra226pb210->Project("hTeO2ra226pb210M1", "Ener2", "Multiplicity==1");
	tTeO2ra228pb208->Project("hTeO2ra228pb208M1", "Ener2", "Multiplicity==1");
	tTeO2rn222only->Project("hTeO2rn222onlyM1", "Ener2", "Multiplicity==1");
	tTeO2rn222pb206->Project("hTeO2rn222pb206M1", "Ener2", "Multiplicity==1");
	tTeO2sb125->Project("hTeO2sb125M1", "Ener2", "Multiplicity==1");
	tTeO2te125->Project("hTeO2te125M1", "Ener2", "Multiplicity==1");
	tTeO2th230only->Project("hTeO2th230onlyM1", "Ener2", "Multiplicity==1");
	tTeO2th230pb206->Project("hTeO2th230pb206M1", "Ener2", "Multiplicity==1");
	tTeO2th232only->Project("hTeO2th232onlyM1", "Ener2", "Multiplicity==1");
	tTeO2th232->Project("hTeO2th232M1", "Ener2", "Multiplicity==1");
	tTeO2th232th228->Project("hTeO2th232th228M1", "Ener2", "Multiplicity==1");
	tTeO2u238ra226->Project("hTeO2u238ra226M1", "Ener2", "Multiplicity==1");
	tTeO2u238rn222->Project("hTeO2u238rn222M1", "Ener2", "Multiplicity==1");
	tTeO2u238->Project("hTeO2u238M1", "Ener2", "Multiplicity==1");
	tTeO2u238th230->Project("hTeO2u238th230M1", "Ener2", "Multiplicity==1");
	tTeO2u238u234->Project("hTeO2u238u234M1", "Ener2", "Multiplicity==1");


////////////////////////////////////////

	tCuBoxco58->Project("hCuBoxco58M2", "Ener2", "Multiplicity==2");
	tCuBoxco60->Project("hCuBoxco60M2", "Ener2", "Multiplicity==2");
	tCuBoxcs137->Project("hCuBoxcs137M2", "Ener2", "Multiplicity==2");
	tCuBoxk40->Project("hCuBoxk40M2", "Ener2", "Multiplicity==2");
	tCuBoxmn54->Project("hCuBoxmn54M2", "Ener2", "Multiplicity==2");
	tCuBoxpb210->Project("hCuBoxpb210M2", "Ener2", "Multiplicity==2");
	tCuBoxth232->Project("hCuBoxth232M2", "Ener2", "Multiplicity==2");
	tCuBoxu238->Project("hCuBoxu238M2", "Ener2", "Multiplicity==2");	
	tCuBoxeu152->Project("hCuBoxeu152M2", "Ener2", "Multiplicity==2");
	tCuBoxeu154->Project("hCuBoxeu154M2", "Ener2", "Multiplicity==2");

	tCuFrameco58->Project("hCuFrameco58M2", "Ener2", "Multiplicity==2");
	tCuFrameco60->Project("hCuFrameco60M2", "Ener2", "Multiplicity==2");
	tCuFramecs137->Project("hCuFramecs137M2", "Ener2", "Multiplicity==2");
	tCuFramek40->Project("hCuFramek40M2", "Ener2", "Multiplicity==2");
	tCuFramemn54->Project("hCuFramemn54M2", "Ener2", "Multiplicity==2");
	tCuFramepb210->Project("hCuFramepb210M2", "Ener2", "Multiplicity==2");
	tCuFrameth232->Project("hCuFrameth232M2", "Ener2", "Multiplicity==2");
	tCuFrameu238->Project("hCuFrameu238M2", "Ener2", "Multiplicity==2");

	// tCuBox_CuFrameco60->Project("hCuBox_CuFrameco60M2", "Ener2", "Multiplicity==2");
	// tCuBox_CuFramek40->Project("hCuBox_CuFramek40M2", "Ener2", "Multiplicity==2");
	// tCuBox_CuFrameth232->Project("hCuBox_CuFrameth232M2", "Ener2", "Multiplicity==2");
	// tCuBox_CuFrameu238->Project("hCuBox_CuFrameu238M2", "Ener2", "Multiplicity==2");

	// tInternalco60->Project("hInternalco60M2", "Ener2", "Multiplicity==2");
	// tInternalk40->Project("hInternalk40M2", "Ener2", "Multiplicity==2");
	// tInternalth232->Project("hInternalth232M2", "Ener2", "Multiplicity==2");
	// tInternalu238->Project("hInternalu238M2", "Ener2", "Multiplicity==2");

	// tPTFEth232->Project("hPTFEth232M2", "Ener2", "Multiplicity==2");
	// tPTFEu238->Project("hPTFEu238M2", "Ener2", "Multiplicity==2");

	tTeO20nu->Project("hTeO20nuM2", "Ener2", "Multiplicity==2");
	tTeO22nu->Project("hTeO22nuM2", "Ener2", "Multiplicity==2");
	tTeO2co60->Project("hTeO2co60M2", "Ener2", "Multiplicity==2");
	tTeO2k40->Project("hTeO2k40M2", "Ener2", "Multiplicity==2");
	tTeO2pb210->Project("hTeO2pb210M2", "Ener2", "Multiplicity==2");
	tTeO2po210->Project("hTeO2po210M2", "Ener2", "Multiplicity==2");
	tTeO2ra226only->Project("hTeO2ra226onlyM2", "Ener2", "Multiplicity==2");
	tTeO2ra226pb206->Project("hTeO2ra226pb206M2", "Ener2", "Multiplicity==2");
	tTeO2ra226pb210->Project("hTeO2ra226pb210M2", "Ener2", "Multiplicity==2");
	tTeO2ra228pb208->Project("hTeO2ra228pb208M2", "Ener2", "Multiplicity==2");
	tTeO2rn222only->Project("hTeO2rn222onlyM2", "Ener2", "Multiplicity==2");
	tTeO2rn222pb206->Project("hTeO2rn222pb206M2", "Ener2", "Multiplicity==2");
	tTeO2sb125->Project("hTeO2sb125M2", "Ener2", "Multiplicity==2");
	tTeO2te125->Project("hTeO2te125M2", "Ener2", "Multiplicity==2");
	tTeO2th230only->Project("hTeO2th230onlyM2", "Ener2", "Multiplicity==2");
	tTeO2th230pb206->Project("hTeO2th230pb206M2", "Ener2", "Multiplicity==2");
	tTeO2th232only->Project("hTeO2th232onlyM2", "Ener2", "Multiplicity==2");
	tTeO2th232->Project("hTeO2th232M2", "Ener2", "Multiplicity==2");
	tTeO2th232th228->Project("hTeO2th232th228M2", "Ener2", "Multiplicity==2");
	tTeO2u238ra226->Project("hTeO2u238ra226M2", "Ener2", "Multiplicity==2");
	tTeO2u238rn222->Project("hTeO2u238rn222M2", "Ener2", "Multiplicity==2");
	tTeO2u238->Project("hTeO2u238M2", "Ener2", "Multiplicity==2");
	tTeO2u238th230->Project("hTeO2u238th230M2", "Ener2", "Multiplicity==2");
	tTeO2u238u234->Project("hTeO2u238u234M2", "Ener2", "Multiplicity==2");


	//////////////////////////////////
	// CuBox+CuFrame here
	// Total = CuBox + 0.2444254783*CuFrame
	hCuBox_CuFrameco60M1->Add(hCuBoxco60M1, hCuFrameco60M1, 1.0, 0.2444254783);
	hCuBox_CuFramek40M1->Add(hCuBoxk40M1, hCuFramek40M1, 1.0, 0.2444254783);
	hCuBox_CuFrameth232M1->Add(hCuBoxth232M1, hCuFrameth232M1, 1.0, 0.2444254783);
	hCuBox_CuFrameu238M1->Add(hCuBoxu238M1, hCuFrameu238M1, 1.0, 0.2444254783);

	hCuBox_CuFrameco60M2->Add(hCuBoxco60M2, hCuFrameco60M2, 1.0, 0.2444254783);
	hCuBox_CuFramek40M2->Add(hCuBoxk40M2, hCuFramek40M2, 1.0, 0.2444254783);
	hCuBox_CuFrameth232M2->Add(hCuBoxth232M2, hCuFrameth232M2, 1.0, 0.2444254783);
	hCuBox_CuFrameu238M2->Add(hCuBoxu238M2, hCuFrameu238M2, 1.0, 0.2444254783);

	//////////////////////////////////

	NormalizePDFPair( hCuBoxco58M1, hCuBoxco58M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxco60M1, hCuBoxco60M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxcs137M1, hCuBoxcs137M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxk40M1, hCuBoxk40M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxmn54M1, hCuBoxmn54M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxpb210M1, hCuBoxpb210M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxth232M1, hCuBoxth232M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxu238M1, hCuBoxu238M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxeu152M1, hCuBoxeu152M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxeu154M1, hCuBoxeu154M2, dMinEnergy, dMaxEnergy);


	NormalizePDFPair( hCuFrameco58M1, hCuFrameco58M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFrameco60M1, hCuFrameco60M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFramecs137M1, hCuFramecs137M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFramek40M1, hCuFramek40M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFramemn54M1, hCuFramemn54M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFramepb210M1, hCuFramepb210M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFrameth232M1, hCuFrameth232M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFrameu238M1, hCuFrameu238M2, dMinEnergy, dMaxEnergy);

	NormalizePDFPair( hCuBox_CuFrameco60M1, hCuBox_CuFrameco60M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBox_CuFramek40M1, hCuBox_CuFramek40M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBox_CuFrameth232M1, hCuBox_CuFrameth232M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBox_CuFrameu238M1, hCuBox_CuFrameu238M2,  dMinEnergy, dMaxEnergy);

	// NormalizePDFPair( hPTFEth232M1, hPTFEth232M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hPTFEu238M1, hPTFEu238M2, dMinEnergy, dMaxEnergy);	

	NormalizePDFPair( hTeO20nuM1, hTeO20nuM2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO22nuM1, hTeO22nuM2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2co60M1, hTeO2co60M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2k40M1, hTeO2k40M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2pb210M1, hTeO2pb210M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2po210M1, hTeO2po210M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2ra226onlyM1, hTeO2ra226onlyM2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2ra226pb206M1, hTeO2ra226pb206M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2ra226pb210M1, hTeO2ra226pb210M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2ra228pb208M1, hTeO2ra228pb208M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2rn222onlyM1, hTeO2rn222onlyM2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2rn222pb206M1, hTeO2rn222pb206M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2sb125M1, hTeO2sb125M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2te125M1, hTeO2te125M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2th230onlyM1, hTeO2th230onlyM2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2th230pb206M1, hTeO2th230pb206M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2th232onlyM1, hTeO2th232onlyM2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2th232M1, hTeO2th232M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2th232th228M1, hTeO2th232th228M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2u238ra226M1, hTeO2u238ra226M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2u238rn222M1, hTeO2u238rn222M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2u238M1,  hTeO2u238M2,   dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2u238th230M1, hTeO2u238th230M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2u238u234M1,  hTeO2u238u234M2,   dMinEnergy, dMaxEnergy);


	TFile *file1 = new TFile("MCProduction_BulkInner_1keV.root", "RECREATE");

	hCuBoxco58M1->Write();
	hCuBoxco60M1->Write();
	hCuBoxcs137M1->Write();
	hCuBoxeu152M1->Write();
	hCuBoxeu154M1->Write();
	hCuBoxk40M1->Write();
	hCuBoxmn54M1->Write();
	hCuBoxpb210M1->Write();
	hCuBoxth232M1->Write();
	hCuBoxu238M1->Write();	

	hCuFrameco58M1->Write();
	hCuFrameco60M1->Write();
	hCuFramecs137M1->Write();
	hCuFramek40M1->Write();
	hCuFramemn54M1->Write();
	hCuFramepb210M1->Write();
	hCuFrameth232M1->Write();
	hCuFrameu238M1->Write();

	hCuBox_CuFrameco60M1->Write();
	hCuBox_CuFramek40M1->Write();
	hCuBox_CuFrameth232M1->Write();
	hCuBox_CuFrameu238M1->Write();

	// hPTFEth232M1->Write();
	// hPTFEu238M1->Write();

	hTeO20nuM1->Write();
	hTeO22nuM1->Write();
	hTeO2co60M1->Write();
	hTeO2k40M1->Write();
	hTeO2pb210M1->Write();
	hTeO2po210M1->Write();
	hTeO2ra226onlyM1->Write();
	hTeO2ra226pb206M1->Write();
	hTeO2ra226pb210M1->Write();
	hTeO2ra228pb208M1->Write();
	hTeO2rn222onlyM1->Write();
	hTeO2rn222pb206M1->Write();
	hTeO2sb125M1->Write();
	hTeO2te125M1->Write();
	hTeO2th230onlyM1->Write();
	hTeO2th230pb206M1->Write();	
	hTeO2th232onlyM1->Write();	
	hTeO2th232M1->Write();
	hTeO2th232th228M1->Write();
	hTeO2u238ra226M1->Write();
	hTeO2u238rn222M1->Write();
	hTeO2u238M1->Write();
	hTeO2u238th230M1->Write();
	hTeO2u238u234M1->Write();

//////////////////////////////

	hCuBoxco58M2->Write();
	hCuBoxco60M2->Write();
	hCuBoxcs137M2->Write();
	hCuBoxeu152M2->Write();
	hCuBoxeu154M2->Write();
	hCuBoxk40M2->Write();
	hCuBoxmn54M2->Write();
	hCuBoxpb210M2->Write();
	hCuBoxth232M2->Write();
	hCuBoxu238M2->Write();	

	hCuFrameco58M2->Write();
	hCuFrameco60M2->Write();
	hCuFramecs137M2->Write();
	hCuFramek40M2->Write();
	hCuFramemn54M2->Write();
	hCuFramepb210M2->Write();
	hCuFrameth232M2->Write();
	hCuFrameu238M2->Write();

	hCuBox_CuFrameco60M2->Write();
	hCuBox_CuFramek40M2->Write();
	hCuBox_CuFrameth232M2->Write();
	hCuBox_CuFrameu238M2->Write();

	// hPTFEth232M2->Write();
	// hPTFEu238M2->Write();

	hTeO20nuM2->Write();
	hTeO22nuM2->Write();
	hTeO2co60M2->Write();
	hTeO2k40M2->Write();
	hTeO2pb210M2->Write();
	hTeO2po210M2->Write();
	hTeO2ra226onlyM2->Write();
	hTeO2ra226pb206M2->Write();
	hTeO2ra226pb210M2->Write();
	hTeO2ra228pb208M2->Write();
	hTeO2rn222onlyM2->Write();
	hTeO2rn222pb206M2->Write();
	hTeO2sb125M2->Write();
	hTeO2te125M2->Write();
	hTeO2th230onlyM2->Write();
	hTeO2th230pb206M2->Write();	
	hTeO2th232onlyM2->Write();
	hTeO2th232M2->Write();
	hTeO2th232th228M2->Write();
	hTeO2u238ra226M2->Write();
	hTeO2u238rn222M2->Write();
	hTeO2u238M2->Write();
	hTeO2u238th230M2->Write();
	hTeO2u238u234M2->Write();

	file1->Write();

}


void SaveHistogramsBulkOuter()
{

	std::string sDataDir = "/cuore/user/zhubrian/MC/Bkg/Production/";

	double dMinEnergy = 0;
	double dMaxEnergy = 10000;
	int dBinSize = 1;
	int dNBins = (dMaxEnergy - dMinEnergy)/dBinSize;


	TH1D *h50mKco58M1 = new TH1D("h50mKco58M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKco60M1 = new TH1D("h50mKco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKcs137M1 = new TH1D("h50mKcs137M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKk40M1 = new TH1D("h50mKk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKmn54M1 = new TH1D("h50mKmn54M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKpb210M1 = new TH1D("h50mKpb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKth232M1 = new TH1D("h50mKth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKu238M1 = new TH1D("h50mKu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *h600mKco60M1 = new TH1D("h600mKco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKk40M1 = new TH1D("h600mKk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKth232M1 = new TH1D("h600mKth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKu238M1 = new TH1D("h600mKu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hExtPbbi210M1 = new TH1D("hExtPbbi210M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hIVCco60M1 = new TH1D("hIVCco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCk40M1 = new TH1D("hIVCk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCth232M1 = new TH1D("hIVCth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCu238M1 = new TH1D("hIVCu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hMBco60M1 = new TH1D("hMBco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBk40M1 = new TH1D("hMBk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBth232M1 = new TH1D("hMBth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBu238M1 = new TH1D("hMBu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hInternalco60M1 = new TH1D("hInternalco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalk40M1 = new TH1D("hInternalk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalth232M1 = new TH1D("hInternalth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalu238M1 = new TH1D("hInternalu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hMCco60M1 = new TH1D("hMCco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMCk40M1 = new TH1D("hMCk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMCth232M1 = new TH1D("hMCth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMCu238M1 = new TH1D("hMCu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hOVCbi207M1 = new TH1D("hOVCbi207M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCco60M1 = new TH1D("hOVCco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCk40M1 = new TH1D("hOVCk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCth232M1 = new TH1D("hOVCth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCu238M1 = new TH1D("hOVCu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	// TH1D *hPbRombi207M1 = new TH1D("hPbRombi207M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomco60M1 = new TH1D("hPbRomco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomcs137M1 = new TH1D("hPbRomcs137M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomeu152M1 = new TH1D("hPbRomeu152M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomeu154M1 = new TH1D("hPbRomeu154M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomk40M1 = new TH1D("hPbRomk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRompb210M1 = new TH1D("hPbRompb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomth232M1 = new TH1D("hPbRomth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomu238M1 = new TH1D("hPbRomu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hSIk40M1 = new TH1D("hSIk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hSIth232M1 = new TH1D("hSIth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hSIu238M1 = new TH1D("hSIu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

////////////////////////

	TH1D *h50mKco58M2 = new TH1D("h50mKco58M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKco60M2 = new TH1D("h50mKco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKcs137M2 = new TH1D("h50mKcs137M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKk40M2 = new TH1D("h50mKk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKmn54M2 = new TH1D("h50mKmn54M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKpb210M2 = new TH1D("h50mKpb210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKth232M2 = new TH1D("h50mKth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKu238M2 = new TH1D("h50mKu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *h600mKco60M2 = new TH1D("h600mKco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKk40M2 = new TH1D("h600mKk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKth232M2 = new TH1D("h600mKth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKu238M2 = new TH1D("h600mKu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hExtPbbi210M2 = new TH1D("hExtPbbi210M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hIVCco60M2 = new TH1D("hIVCco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCk40M2 = new TH1D("hIVCk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCth232M2 = new TH1D("hIVCth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCu238M2 = new TH1D("hIVCu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hInternalco60M2 = new TH1D("hInternalco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalk40M2 = new TH1D("hInternalk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalth232M2 = new TH1D("hInternalth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalu238M2 = new TH1D("hInternalu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hMBco60M2 = new TH1D("hMBco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBk40M2 = new TH1D("hMBk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBth232M2 = new TH1D("hMBth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBu238M2 = new TH1D("hMBu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hMCco60M2 = new TH1D("hMCco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMCk40M2 = new TH1D("hMCk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMCth232M2 = new TH1D("hMCth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMCu238M2 = new TH1D("hMCu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hOVCbi207M2 = new TH1D("hOVCbi207M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCco60M2 = new TH1D("hOVCco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCk40M2 = new TH1D("hOVCk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCth232M2 = new TH1D("hOVCth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCu238M2 = new TH1D("hOVCu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	// TH1D *hPbRombi207M2 = new TH1D("hPbRombi207M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomco60M2 = new TH1D("hPbRomco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomcs137M2 = new TH1D("hPbRomcs137M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomeu152M2 = new TH1D("hPbRomeu152M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomeu154M2 = new TH1D("hPbRomeu154M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomk40M2 = new TH1D("hPbRomk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRompb210M2 = new TH1D("hPbRompb210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomth232M2 = new TH1D("hPbRomth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomu238M2 = new TH1D("hPbRomu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hSIk40M2 = new TH1D("hSIk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hSIth232M2 = new TH1D("hSIth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hSIu238M2 = new TH1D("hSIu238M2", "", dNBins, dMinEnergy, dMaxEnergy);



///////////////////////////////

	TChain *t50mKco58 = LoadMC(sDataDir.c_str(), "50mK", "co58");
	TChain *t50mKco60 = LoadMC(sDataDir.c_str(), "50mK", "co60");
	TChain *t50mKcs137 = LoadMC(sDataDir.c_str(), "50mK", "cs137");
	TChain *t50mKk40 = LoadMC(sDataDir.c_str(), "50mK", "k40");
	TChain *t50mKmn54 = LoadMC(sDataDir.c_str(), "50mK", "mn54");
	TChain *t50mKpb210 = LoadMC(sDataDir.c_str(), "50mK", "pb210");
	TChain *t50mKth232 = LoadMC(sDataDir.c_str(), "50mK", "th232");
	TChain *t50mKu238 = LoadMC(sDataDir.c_str(), "50mK", "u238");

	TChain *t600mKco60 = LoadMC(sDataDir.c_str(), "600mK", "co60");
	TChain *t600mKk40 = LoadMC(sDataDir.c_str(), "600mK", "k40");
	TChain *t600mKth232 = LoadMC(sDataDir.c_str(), "600mK", "th232");
	TChain *t600mKu238 = LoadMC(sDataDir.c_str(), "600mK", "u238");

	TChain *tExtPbbi210 = LoadMC(sDataDir.c_str(), "ExtPb4cm_210BiBStot", "1.23T");

	TChain *tIVCco60 = LoadMC(sDataDir.c_str(), "IVC", "co60");
	TChain *tIVCk40 = LoadMC(sDataDir.c_str(), "IVC", "k40");
	TChain *tIVCth232 = LoadMC(sDataDir.c_str(), "IVC", "th232");
	TChain *tIVCu238 = LoadMC(sDataDir.c_str(), "IVC", "u238");

	// TChain *tOVCbi207 = LoadMC(sDataDir.c_str(), "OVC", "bi207");
	TChain *tOVCco60 = LoadMC(sDataDir.c_str(), "OVC", "co60");
	TChain *tOVCk40 = LoadMC(sDataDir.c_str(), "OVC", "k40");
	TChain *tOVCth232 = LoadMC(sDataDir.c_str(), "OVC", "th232");
	TChain *tOVCu238 = LoadMC(sDataDir.c_str(), "OVC", "u238");

	// TChain *tSIk40 = LoadMC(sDataDir.c_str(), "SI", "k40");
	// TChain *tSIth232 = LoadMC(sDataDir.c_str(), "SI", "th232");
	// TChain *tSIu238 = LoadMC(sDataDir.c_str(), "SI", "u238");

	// TChain *tMBco60 = LoadMC(sDataDir.c_str(), "MB", "co60");
	// TChain *tMBk40 = LoadMC(sDataDir.c_str(), "MB", "k40");
	// TChain *tMBth232 = LoadMC(sDataDir.c_str(), "MB", "th232");
	// TChain *tMBu238 = LoadMC(sDataDir.c_str(), "MB", "u238");

	// TChain *tMCco60 = LoadMC(sDataDir.c_str(), "MC", "co60");
	// TChain *tMCk40 = LoadMC(sDataDir.c_str(), "MC", "k40");
	// TChain *tMCth232 = LoadMC(sDataDir.c_str(), "MC", "th232");
	// TChain *tMCu238 = LoadMC(sDataDir.c_str(), "MC", "u238");

	// TChain *tPbRombi207 = LoadMC(sDataDir.c_str(), "PbRom", "bi207");
	TChain *tPbRomco60 = LoadMC(sDataDir.c_str(), "PbRom", "co60");
	TChain *tPbRomcs137 = LoadMC(sDataDir.c_str(), "PbRom", "cs137");
	TChain *tPbRomeu152 = LoadMC(sDataDir.c_str(), "PbRom", "eu152");
	TChain *tPbRomeu154 = LoadMC(sDataDir.c_str(), "PbRom", "eu154");
	TChain *tPbRomk40 = LoadMC(sDataDir.c_str(), "PbRom", "k40");
	TChain *tPbRompb210 = LoadMC(sDataDir.c_str(), "PbRom", "pb210");
	TChain *tPbRomth232 = LoadMC(sDataDir.c_str(), "PbRom", "th232");
	TChain *tPbRomu238 = LoadMC(sDataDir.c_str(), "PbRom", "u238");

////////////////////////////////////////

	t50mKco58->Project("h50mKco58M1", "Ener2", "Multiplicity==1");
	t50mKco60->Project("h50mKco60M1", "Ener2", "Multiplicity==1");
	t50mKcs137->Project("h50mKcs137M1", "Ener2", "Multiplicity==1");
	t50mKk40->Project("h50mKk40M1", "Ener2", "Multiplicity==1");
	t50mKmn54->Project("h50mKmn54M1", "Ener2", "Multiplicity==1");
	t50mKpb210->Project("h50mKpb210M1", "Ener2", "Multiplicity==1");
	t50mKth232->Project("h50mKth232M1", "Ener2", "Multiplicity==1");
	t50mKu238->Project("h50mKu238M1", "Ener2", "Multiplicity==1");

	t600mKco60->Project("h600mKco60M1", "Ener2", "Multiplicity==1");
	t600mKk40->Project("h600mKk40M1", "Ener2", "Multiplicity==1");
	t600mKth232->Project("h600mKth232M1", "Ener2", "Multiplicity==1");
	t600mKu238->Project("h600mKu238M1", "Ener2", "Multiplicity==1");

	tExtPbbi210->Project("hExtPbbi210M1", "Ener2", "Multiplicity==1");

	tIVCco60->Project("hIVCco60M1", "Ener2", "Multiplicity==1");
	tIVCk40->Project("hIVCk40M1", "Ener2", "Multiplicity==1");
	tIVCth232->Project("hIVCth232M1", "Ener2", "Multiplicity==1");
	tIVCu238->Project("hIVCu238M1", "Ener2", "Multiplicity==1");

	// tOVCbi207->Project("hOVCbi207M1", "Ener2", "Multiplicity==1");
	tOVCco60->Project("hOVCco60M1", "Ener2", "Multiplicity==1");
	tOVCk40->Project("hOVCk40M1", "Ener2", "Multiplicity==1");
	tOVCth232->Project("hOVCth232M1", "Ener2", "Multiplicity==1");
	tOVCu238->Project("hOVCu238M1", "Ener2", "Multiplicity==1");

	// tMBco60->Project("hMBco60M1", "Ener2", "Multiplicity==1");
	// tMBk40->Project("hMBk40M1", "Ener2", "Multiplicity==1");
	// tMBth232->Project("hMBth232M1", "Ener2", "Multiplicity==1");
	// tMBu238->Project("hMBu238M1", "Ener2", "Multiplicity==1");

	// tMCco60->Project("hMCco60M1", "Ener2", "Multiplicity==1");
	// tMCk40->Project("hMCk40M1", "Ener2", "Multiplicity==1");
	// tMCth232->Project("hMCth232M1", "Ener2", "Multiplicity==1");
	// tMCu238->Project("hMCu238M1", "Ener2", "Multiplicity==1");

	// tSIk40->Project("hSIk40M1", "Ener2", "Multiplicity==1");
	// tSIth232->Project("hSIth232M1", "Ener2", "Multiplicity==1");
	// tSIu238->Project("hSIu238M1", "Ener2", "Multiplicity==1");

	// tPbRombi207->Project("hPbRombi207M1", "Ener2", "Multiplicity==1");
	tPbRomco60->Project("hPbRomco60M1", "Ener2", "Multiplicity==1");
	tPbRomcs137->Project("hPbRomcs137M1", "Ener2", "Multiplicity==1");
	tPbRomeu152->Project("hPbRomeu152M1", "Ener2", "Multiplicity==1");
	tPbRomeu154->Project("hPbRomeu154M1", "Ener2", "Multiplicity==1");
	tPbRomk40->Project("hPbRomk40M1", "Ener2", "Multiplicity==1");
	tPbRompb210->Project("hPbRompb210M1", "Ener2", "Multiplicity==1");
	tPbRomth232->Project("hPbRomth232M1", "Ener2", "Multiplicity==1");
	tPbRomu238->Project("hPbRomu238M1", "Ener2", "Multiplicity==1");


////////////////////////////////////////

	t50mKco58->Project("h50mKco58M2", "Ener2", "Multiplicity==2");
	t50mKco60->Project("h50mKco60M2", "Ener2", "Multiplicity==2");
	t50mKcs137->Project("h50mKcs137M2", "Ener2", "Multiplicity==2");
	t50mKk40->Project("h50mKk40M2", "Ener2", "Multiplicity==2");
	t50mKmn54->Project("h50mKmn54M2", "Ener2", "Multiplicity==2");
	t50mKpb210->Project("h50mKpb210M2", "Ener2", "Multiplicity==2");
	t50mKth232->Project("h50mKth232M2", "Ener2", "Multiplicity==2");
	t50mKu238->Project("h50mKu238M2", "Ener2", "Multiplicity==2");

	t600mKco60->Project("h600mKco60M2", "Ener2", "Multiplicity==2");
	t600mKk40->Project("h600mKk40M2", "Ener2", "Multiplicity==2");
	t600mKth232->Project("h600mKth232M2", "Ener2", "Multiplicity==2");
	t600mKu238->Project("h600mKu238M2", "Ener2", "Multiplicity==2");

	tExtPbbi210->Project("hExtPbbi210M2", "Ener2", "Multiplicity==2");

	tIVCco60->Project("hIVCco60M2", "Ener2", "Multiplicity==2");
	tIVCk40->Project("hIVCk40M2", "Ener2", "Multiplicity==2");
	tIVCth232->Project("hIVCth232M2", "Ener2", "Multiplicity==2");
	tIVCu238->Project("hIVCu238M2", "Ener2", "Multiplicity==2");

	// tOVCbi207->Project("hOVCbi207M2", "Ener2", "Multiplicity==2");
	tOVCco60->Project("hOVCco60M2", "Ener2", "Multiplicity==2");
	tOVCk40->Project("hOVCk40M2", "Ener2", "Multiplicity==2");
	tOVCth232->Project("hOVCth232M2", "Ener2", "Multiplicity==2");
	tOVCu238->Project("hOVCu238M2", "Ener2", "Multiplicity==2");

	// tMBco60->Project("hMBco60M2", "Ener2", "Multiplicity==2");
	// tMBk40->Project("hMBk40M2", "Ener2", "Multiplicity==2");
	// tMBth232->Project("hMBth232M2", "Ener2", "Multiplicity==2");
	// tMBu238->Project("hMBu238M2", "Ener2", "Multiplicity==2");

	// tMCco60->Project("hMCco60M2", "Ener2", "Multiplicity==2");
	// tMCk40->Project("hMCk40M2", "Ener2", "Multiplicity==2");
	// tMCth232->Project("hMCth232M2", "Ener2", "Multiplicity==2");
	// tMCu238->Project("hMCu238M2", "Ener2", "Multiplicity==2");

	// tSIk40->Project("hSIk40M2", "Ener2", "Multiplicity==2");
	// tSIth232->Project("hSIth232M2", "Ener2", "Multiplicity==2");
	// tSIu238->Project("hSIu238M2", "Ener2", "Multiplicity==2");

	// tPbRombi207->Project("hPbRombi207M2", "Ener2", "Multiplicity==2");
	tPbRomco60->Project("hPbRomco60M2", "Ener2", "Multiplicity==2");
	tPbRomcs137->Project("hPbRomcs137M2", "Ener2", "Multiplicity==2");
	tPbRomeu152->Project("hPbRomeu152M2", "Ener2", "Multiplicity==2");
	tPbRomeu154->Project("hPbRomeu154M2", "Ener2", "Multiplicity==2");
	tPbRomk40->Project("hPbRomk40M2", "Ener2", "Multiplicity==2");
	tPbRompb210->Project("hPbRompb210M2", "Ener2", "Multiplicity==2");
	tPbRomth232->Project("hPbRomth232M2", "Ener2", "Multiplicity==2");
	tPbRomu238->Project("hPbRomu238M2", "Ener2", "Multiplicity==2");

	//////////////////////////////////
	// Internal Shielding here:
	// Total = 600mK + 0.8625620249*50mK + 2.3168024764*IVC
	hInternalco60M1->Add(h600mKco60M1);
	hInternalco60M1->Add(h50mKco60M1, 0.8625620249);
	hInternalco60M1->Add(hIVCco60M1, 2.3168024764);

	hInternalk40M1->Add(h600mKk40M1);
	hInternalk40M1->Add(h50mKk40M1, 0.8625620249);
	hInternalk40M1->Add(hIVCk40M1, 2.3168024764);

	hInternalth232M1->Add(h600mKth232M1);
	hInternalth232M1->Add(h50mKth232M1, 0.8625620249);
	hInternalth232M1->Add(hIVCth232M1, 2.3168024764);

	hInternalu238M1->Add(h600mKu238M1);
	hInternalu238M1->Add(h50mKu238M1, 0.8625620249);
	hInternalu238M1->Add(hIVCu238M1, 2.3168024764);

	hInternalco60M2->Add(h600mKco60M2);
	hInternalco60M2->Add(h50mKco60M2, 0.8625620249);
	hInternalco60M2->Add(hIVCco60M2, 2.3168024764);

	hInternalk40M2->Add(h600mKk40M2);
	hInternalk40M2->Add(h50mKk40M2, 0.8625620249);
	hInternalk40M2->Add(hIVCk40M2, 2.3168024764);

	hInternalth232M2->Add(h600mKth232M2);
	hInternalth232M2->Add(h50mKth232M2, 0.8625620249);
	hInternalth232M2->Add(hIVCth232M2, 2.3168024764);

	hInternalu238M2->Add(h600mKu238M2);
	hInternalu238M2->Add(h50mKu238M2, 0.8625620249);
	hInternalu238M2->Add(hIVCu238M2, 2.3168024764);

	//////////////////////////////////


	NormalizePDFPair( h50mKco58M1, h50mKco58M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKco60M1, h50mKco60M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKcs137M1, h50mKcs137M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKk40M1, h50mKk40M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKmn54M1, h50mKmn54M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKpb210M1, h50mKpb210M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKth232M1, h50mKth232M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKu238M1, h50mKu238M2, dMinEnergy, dMaxEnergy);

	NormalizePDFPair( h600mKco60M1, h600mKco60M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h600mKk40M1, h600mKk40M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h600mKth232M1, h600mKth232M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h600mKu238M1, h600mKu238M2, dMinEnergy, dMaxEnergy);

	NormalizePDFPair( hIVCco60M1, hIVCco60M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hIVCk40M1, hIVCk40M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hIVCth232M1, hIVCth232M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hIVCu238M1, hIVCu238M2, dMinEnergy, dMaxEnergy);

	NormalizePDFPair( hInternalco60M1, hInternalco60M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hInternalk40M1, hInternalk40M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hInternalth232M1, hInternalth232M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hInternalu238M1, hInternalu238M2, dMinEnergy, dMaxEnergy);


	// NormalizePDFPair( hOVCbi207M1, hOVCbi207M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hOVCco60M1, hOVCco60M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hOVCk40M1, hOVCk40M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hOVCth232M1, hOVCth232M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hOVCu238M1, hOVCu238M2, dMinEnergy, dMaxEnergy);

	// NormalizePDFPair( hMBco60M1, hMBco60M2,  dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hMBk40M1, hMBk40M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hMBth232M1, hMBth232M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hMBu238M1, hMBu238M2, dMinEnergy, dMaxEnergy);

	// NormalizePDFPair( hMCco60M1, hMCco60M2,  dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hMCk40M1, hMCk40M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hMCth232M1, hMCth232M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hMCu238M1, hMCu238M2, dMinEnergy, dMaxEnergy);

	// NormalizePDFPair( hSIk40M1, hSIk40M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hSIth232M1, hSIth232M2, dMinEnergy, dMaxEnergy);
	// NormalizePDFPair( hSIu238M1, hSIu238M2, dMinEnergy, dMaxEnergy);


	// NormalizePDFPair( hPbRombi207M1, hPbRombi207M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomco60M1, hPbRomco60M2,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomcs137M1, hPbRomcs137M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomeu152M1, hPbRomeu152M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomeu154M1, hPbRomeu154M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomk40M1, hPbRomk40M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRompb210M1, hPbRompb210M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomth232M1, hPbRomth232M2, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomu238M1, hPbRomu238M2, dMinEnergy, dMaxEnergy);



	TFile *file1 = new TFile("MCProduction_BulkOuter_1keV.root", "RECREATE");

	h50mKco58M1->Write();
	h50mKco60M1->Write();
	h50mKcs137M1->Write();
	h50mKk40M1->Write();
	h50mKmn54M1->Write();
	h50mKpb210M1->Write();
	h50mKth232M1->Write();
	h50mKu238M1->Write();

	h600mKco60M1->Write();
	h600mKk40M1->Write();
	h600mKth232M1->Write();
	h600mKu238M1->Write();

	hInternalco60M1->Write();
	hInternalk40M1->Write();
	hInternalth232M1->Write();
	hInternalu238M1->Write();	

	hExtPbbi210M1->Write();


	hIVCco60M1->Write();
	hIVCk40M1->Write();
	hIVCth232M1->Write();
	hIVCu238M1->Write();

	// hOVCbi207M1->Write();
	hOVCco60M1->Write();
	hOVCk40M1->Write();
	hOVCth232M1->Write();
	hOVCu238M1->Write();

	// hMBco60M1->Write();
	// hMBk40M1->Write();
	// hMBth232M1->Write();
	// hMBu238M1->Write();

	// hMCco60M1->Write();
	// hMCk40M1->Write();
	// hMCth232M1->Write();
	// hMCu238M1->Write();

	// hSIk40M1->Write();
	// hSIth232M1->Write();
	// hSIu238M1->Write();

	// hPbRombi207M1->Write();
	hPbRomco60M1->Write();
	hPbRomcs137M1->Write();
	hPbRomeu152M1->Write();
	hPbRomeu154M1->Write();
	hPbRomk40M1->Write();
	hPbRompb210M1->Write();
	hPbRomth232M1->Write();
	hPbRomu238M1->Write();

//////////////////////////////

	h50mKco58M2->Write();
	h50mKco60M2->Write();
	h50mKcs137M2->Write();
	h50mKk40M2->Write();
	h50mKmn54M2->Write();
	h50mKpb210M2->Write();
	h50mKth232M2->Write();
	h50mKu238M2->Write();

	h600mKco60M2->Write();
	h600mKk40M2->Write();
	h600mKth232M2->Write();
	h600mKu238M2->Write();

	hExtPbbi210M2->Write();

	hIVCco60M2->Write();
	hIVCk40M2->Write();
	hIVCth232M2->Write();
	hIVCu238M2->Write();

	// hOVCbi207M2->Write();
	hOVCco60M2->Write();
	hOVCk40M2->Write();
	hOVCth232M2->Write();
	hOVCu238M2->Write();

	hInternalco60M2->Write();
	hInternalk40M2->Write();
	hInternalth232M2->Write();
	hInternalu238M2->Write();	

	// hMBco60M2->Write();
	// hMBk40M2->Write();
	// hMBth232M2->Write();
	// hMBu238M2->Write();

	// hMCco60M2->Write();
	// hMCk40M2->Write();
	// hMCth232M2->Write();
	// hMCu238M2->Write();

	// hSIk40M2->Write();
	// hSIth232M2->Write();
	// hSIu238M2->Write();

	// hPbRombi207M2->Write();
	hPbRomco60M2->Write();
	hPbRomcs137M2->Write();
	hPbRomeu152M2->Write();
	hPbRomeu154M2->Write();
	hPbRomk40M2->Write();
	hPbRompb210M2->Write();
	hPbRomth232M2->Write();
	hPbRomu238M2->Write();

	file1->Write();

}



void SaveHistogramsBulkM2SumInner()
{
	std::string sDataDir = "/cuore/user/zhubrian/MC/Bkg/Production/";

	double dMinEnergy = 0;
	double dMaxEnergy = 10000;
	int dBinSize = 1;
	int dNBins = (dMaxEnergy - dMinEnergy)/dBinSize;

	TH1D *hCuBoxco58M1 = new TH1D("hCuBoxco58M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxco60M1 = new TH1D("hCuBoxco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxcs137M1 = new TH1D("hCuBoxcs137M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxk40M1 = new TH1D("hCuBoxk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxmn54M1 = new TH1D("hCuBoxmn54M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxpb210M1 = new TH1D("hCuBoxpb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxth232M1 = new TH1D("hCuBoxth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxu238M1 = new TH1D("hCuBoxu238M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxeu152M1 = new TH1D("hCuBoxeu152M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxeu154M1 = new TH1D("hCuBoxeu154M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuFrameco58M1 = new TH1D("hCuFrameco58M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameco60M1 = new TH1D("hCuFrameco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramecs137M1 = new TH1D("hCuFramecs137M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramek40M1 = new TH1D("hCuFramek40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramemn54M1 = new TH1D("hCuFramemn54M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramepb210M1 = new TH1D("hCuFramepb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameth232M1 = new TH1D("hCuFrameth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameu238M1 = new TH1D("hCuFrameu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuBox_CuFrameco60M1 = new TH1D("hCuBox_CuFrameco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFramek40M1 = new TH1D("hCuBox_CuFramek40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M1 = new TH1D("hCuBox_CuFrameth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M1 = new TH1D("hCuBox_CuFrameu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hInternalco60M1 = new TH1D("hInternalco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalk40M1 = new TH1D("hInternalk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalth232M1 = new TH1D("hInternalth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalu238M1 = new TH1D("hInternalu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hPTFEth232M1 = new TH1D("hPTFEth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPTFEu238M1 = new TH1D("hPTFEu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hTeO20nuM1 = new TH1D("hTeO20nuM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO22nuM1 = new TH1D("hTeO22nuM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2co60M1 = new TH1D("hTeO2co60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2k40M1 = new TH1D("hTeO2k40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2pb210M1 = new TH1D("hTeO2pb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2po210M1 = new TH1D("hTeO2po210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226onlyM1 = new TH1D("hTeO2ra226onlyM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226pb206M1 = new TH1D("hTeO2ra226pb206M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226pb210M1 = new TH1D("hTeO2ra226pb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra228pb208M1 = new TH1D("hTeO2ra228pb208M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2rn222onlyM1 = new TH1D("hTeO2rn222onlyM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2rn222pb206M1 = new TH1D("hTeO2rn222pb206M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2sb125M1 = new TH1D("hTeO2sb125M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2te125M1 = new TH1D("hTeO2te125M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th230onlyM1 = new TH1D("hTeO2th230onlyM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th230pb206M1 = new TH1D("hTeO2th230pb206M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232onlyM1 = new TH1D("hTeO2th232onlyM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232M1 = new TH1D("hTeO2th232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232th228M1 = new TH1D("hTeO2th232th228M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238ra226M1 = new TH1D("hTeO2u238ra226M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238rn222M1 = new TH1D("hTeO2u238rn222M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238M1 = new TH1D("hTeO2u238M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238th230M1 = new TH1D("hTeO2u238th230M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238u234M1 = new TH1D("hTeO2u238u234M1", "", dNBins, dMinEnergy, dMaxEnergy);



////////////////////////

	TH1D *hCuBoxco58M2Sum = new TH1D("hCuBoxco58M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxco60M2Sum = new TH1D("hCuBoxco60M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxcs137M2Sum = new TH1D("hCuBoxcs137M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxk40M2Sum = new TH1D("hCuBoxk40M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxmn54M2Sum = new TH1D("hCuBoxmn54M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxpb210M2Sum = new TH1D("hCuBoxpb210M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxth232M2Sum = new TH1D("hCuBoxth232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxu238M2Sum = new TH1D("hCuBoxu238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxeu152M2Sum = new TH1D("hCuBoxeu152M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxeu154M2Sum = new TH1D("hCuBoxeu154M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuFrameco58M2Sum = new TH1D("hCuFrameco58M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameco60M2Sum = new TH1D("hCuFrameco60M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramecs137M2Sum = new TH1D("hCuFramecs137M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramek40M2Sum = new TH1D("hCuFramek40M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramemn54M2Sum = new TH1D("hCuFramemn54M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramepb210M2Sum = new TH1D("hCuFramepb210M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameth232M2Sum = new TH1D("hCuFrameth232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameu238M2Sum = new TH1D("hCuFrameu238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuBox_CuFrameco60M2Sum = new TH1D("hCuBox_CuFrameco60M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFramek40M2Sum = new TH1D("hCuBox_CuFramek40M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M2Sum = new TH1D("hCuBox_CuFrameth232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M2Sum = new TH1D("hCuBox_CuFrameu238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hInternalco60M2Sum = new TH1D("hInternalco60M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalk40M2Sum = new TH1D("hInternalk40M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalth232M2Sum = new TH1D("hInternalth232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hInternalu238M2Sum = new TH1D("hInternalu238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hPTFEth232M2Sum = new TH1D("hPTFEth232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPTFEu238M2Sum = new TH1D("hPTFEu238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hTeO20nuM2Sum = new TH1D("hTeO20nuM2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO22nuM2Sum = new TH1D("hTeO22nuM2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2co60M2Sum = new TH1D("hTeO2co60M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2k40M2Sum = new TH1D("hTeO2k40M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2pb210M2Sum = new TH1D("hTeO2pb210M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2po210M2Sum = new TH1D("hTeO2po210M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226onlyM2Sum = new TH1D("hTeO2ra226onlyM2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226pb206M2Sum = new TH1D("hTeO2ra226pb206M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226pb210M2Sum = new TH1D("hTeO2ra226pb210M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra228pb208M2Sum = new TH1D("hTeO2ra228pb208M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2rn222onlyM2Sum = new TH1D("hTeO2rn222onlyM2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2rn222pb206M2Sum = new TH1D("hTeO2rn222pb206M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2sb125M2Sum = new TH1D("hTeO2sb125M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2te125M2Sum = new TH1D("hTeO2te125M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th230onlyM2Sum = new TH1D("hTeO2th230onlyM2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th230pb206M2Sum = new TH1D("hTeO2th230pb206M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232onlyM2Sum = new TH1D("hTeO2th232onlyM2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232M2Sum = new TH1D("hTeO2th232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232th228M2Sum = new TH1D("hTeO2th232th228M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238ra226M2Sum = new TH1D("hTeO2u238ra226M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238rn222M2Sum = new TH1D("hTeO2u238rn222M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238M2Sum = new TH1D("hTeO2u238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238th230M2Sum = new TH1D("hTeO2u238th230M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238u234M2Sum = new TH1D("hTeO2u238u234M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);



///////////////////////////////

	TChain *tCuBoxco58 = LoadMC(sDataDir.c_str(), "CuBox", "co58");
	TChain *tCuBoxco60 = LoadMC(sDataDir.c_str(), "CuBox", "co60");
	TChain *tCuBoxcs137 = LoadMC(sDataDir.c_str(), "CuBox", "cs137");
	TChain *tCuBoxk40 = LoadMC(sDataDir.c_str(), "CuBox", "k40");
	TChain *tCuBoxmn54 = LoadMC(sDataDir.c_str(), "CuBox", "mn54");
	TChain *tCuBoxpb210 = LoadMC(sDataDir.c_str(), "CuBox", "pb210");
	TChain *tCuBoxth232 = LoadMC(sDataDir.c_str(), "CuBox", "th232");
	TChain *tCuBoxu238 = LoadMC(sDataDir.c_str(), "CuBox", "u238");
	TChain *tCuBoxeu152 = LoadMC(sDataDir.c_str(), "n7_CuBox", "eu152");
	TChain *tCuBoxeu154 = LoadMC(sDataDir.c_str(), "n7_CuBox", "eu154");

	TChain *tCuFrameco58 = LoadMC(sDataDir.c_str(), "CuFrame", "co58");
	TChain *tCuFrameco60 = LoadMC(sDataDir.c_str(), "CuFrame", "co60");
	TChain *tCuFramecs137 = LoadMC(sDataDir.c_str(), "CuFrame", "cs137");
	TChain *tCuFramek40 = LoadMC(sDataDir.c_str(), "CuFrame", "k40");
	TChain *tCuFramemn54 = LoadMC(sDataDir.c_str(), "CuFrame", "mn54");
	TChain *tCuFramepb210 = LoadMC(sDataDir.c_str(), "CuFrame", "pb210");
	TChain *tCuFrameth232 = LoadMC(sDataDir.c_str(), "CuFrame", "th232");
	TChain *tCuFrameu238 = LoadMC(sDataDir.c_str(), "CuFrame", "u238");

	TChain *tCuBox_CuFrameco60 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "co60");
	TChain *tCuBox_CuFramek40 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "k40");
	TChain *tCuBox_CuFrameth232 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "th232");
	TChain *tCuBox_CuFrameu238 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "u238");

	TChain *tInternalco60 = LoadMC(sDataDir.c_str(), "InternalShields", "co60");
	TChain *tInternalk40 = LoadMC(sDataDir.c_str(), "InternalShields", "k40");
	TChain *tInternalth232 = LoadMC(sDataDir.c_str(), "InternalShields", "th232");
	TChain *tInternalu238 = LoadMC(sDataDir.c_str(), "InternalShields", "u238");

	TChain *tPTFEth232 = LoadMC(sDataDir.c_str(), "PTFE", "th232");
	TChain *tPTFEu238 = LoadMC(sDataDir.c_str(), "PTFE", "u238");

	TChain *tTeO20nu = LoadMC(sDataDir.c_str(), "TeO2", "0nu_1um");
	TChain *tTeO22nu = LoadMC(sDataDir.c_str(), "TeO2", "2nu");
	TChain *tTeO2co60 = LoadMC(sDataDir.c_str(), "TeO2", "co60");
	TChain *tTeO2k40 = LoadMC(sDataDir.c_str(), "TeO2", "k40");
	TChain *tTeO2pb210 = LoadMC(sDataDir.c_str(), "TeO2", "pb210");
	TChain *tTeO2po210 = LoadMC(sDataDir.c_str(), "TeO2", "po210");
	TChain *tTeO2ra226only = LoadMC(sDataDir.c_str(), "TeO2", "ra226_only");
	TChain *tTeO2ra226pb206 = LoadMC(sDataDir.c_str(), "TeO2", "ra226-pb206");
	TChain *tTeO2ra226pb210 = LoadMC(sDataDir.c_str(), "TeO2", "ra226-pb210");
	TChain *tTeO2ra228pb208 = LoadMC(sDataDir.c_str(), "TeO2", "ra228-pb208");
	TChain *tTeO2rn222only = LoadMC(sDataDir.c_str(), "TeO2", "rn222_only");
	TChain *tTeO2rn222pb206 = LoadMC(sDataDir.c_str(), "TeO2", "rn222-pb206");
	TChain *tTeO2sb125 = LoadMC(sDataDir.c_str(), "TeO2", "sb125");
	TChain *tTeO2te125 = LoadMC(sDataDir.c_str(), "TeO2", "te125m");
	TChain *tTeO2th230only = LoadMC(sDataDir.c_str(), "TeO2", "th230_only");
	TChain *tTeO2th230pb206 = LoadMC(sDataDir.c_str(), "TeO2", "th230-pb206");
	TChain *tTeO2th232only = LoadMC(sDataDir.c_str(), "TeO2", "th232_only");
	TChain *tTeO2th232 = LoadMC(sDataDir.c_str(), "TeO2", "th232");
	TChain *tTeO2th232th228 = LoadMC(sDataDir.c_str(), "TeO2", "th232-th228");
	TChain *tTeO2u238ra226 = LoadMC(sDataDir.c_str(), "TeO2", "u238-ra226");
	TChain *tTeO2u238rn222 = LoadMC(sDataDir.c_str(), "TeO2", "u238-rn222");
	TChain *tTeO2u238 = LoadMC(sDataDir.c_str(), "TeO2", "u238");
	TChain *tTeO2u238th230 = LoadMC(sDataDir.c_str(), "TeO2", "u238-th230");
	TChain *tTeO2u238u234 = LoadMC(sDataDir.c_str(), "TeO2", "u238-u234");

////////////////////////////////////////


	tCuBoxco58->Project("hCuBoxco58M1", "Ener2", "Multiplicity==1");
	tCuBoxco60->Project("hCuBoxco60M1", "Ener2", "Multiplicity==1");
	tCuBoxcs137->Project("hCuBoxcs137M1", "Ener2", "Multiplicity==1");
	tCuBoxk40->Project("hCuBoxk40M1", "Ener2", "Multiplicity==1");
	tCuBoxmn54->Project("hCuBoxmn54M1", "Ener2", "Multiplicity==1");
	tCuBoxpb210->Project("hCuBoxpb210M1", "Ener2", "Multiplicity==1");
	tCuBoxth232->Project("hCuBoxth232M1", "Ener2", "Multiplicity==1");
	tCuBoxu238->Project("hCuBoxu238M1", "Ener2", "Multiplicity==1");	
	tCuBoxeu152->Project("hCuBoxeu152M1", "Ener2", "Multiplicity==1");
	tCuBoxeu154->Project("hCuBoxeu154M1", "Ener2", "Multiplicity==1");

	tCuFrameco58->Project("hCuFrameco58M1", "Ener2", "Multiplicity==1");
	tCuFrameco60->Project("hCuFrameco60M1", "Ener2", "Multiplicity==1");
	tCuFramecs137->Project("hCuFramecs137M1", "Ener2", "Multiplicity==1");
	tCuFramek40->Project("hCuFramek40M1", "Ener2", "Multiplicity==1");
	tCuFramemn54->Project("hCuFramemn54M1", "Ener2", "Multiplicity==1");
	tCuFramepb210->Project("hCuFramepb210M1", "Ener2", "Multiplicity==1");
	tCuFrameth232->Project("hCuFrameth232M1", "Ener2", "Multiplicity==1");
	tCuFrameu238->Project("hCuFrameu238M1", "Ener2", "Multiplicity==1");

	tCuBox_CuFrameco60->Project("hCuBox_CuFrameco60M1", "Ener2", "Multiplicity==1");
	tCuBox_CuFramek40->Project("hCuBox_CuFramek40M1", "Ener2", "Multiplicity==1");
	tCuBox_CuFrameth232->Project("hCuBox_CuFrameth232M1", "Ener2", "Multiplicity==1");
	tCuBox_CuFrameu238->Project("hCuBox_CuFrameu238M1", "Ener2", "Multiplicity==1");

	tInternalco60->Project("hInternalco60M1", "Ener2", "Multiplicity==1");
	tInternalk40->Project("hInternalk40M1", "Ener2", "Multiplicity==1");
	tInternalth232->Project("hInternalth232M1", "Ener2", "Multiplicity==1");
	tInternalu238->Project("hInternalu238M1", "Ener2", "Multiplicity==1");

	tPTFEth232->Project("hPTFEth232M1", "Ener2", "Multiplicity==1");
	tPTFEu238->Project("hPTFEu238M1", "Ener2", "Multiplicity==1");


	tTeO20nu->Project("hTeO20nuM1", "Ener2", "Multiplicity==1");
	tTeO22nu->Project("hTeO22nuM1", "Ener2", "Multiplicity==1");
	tTeO2co60->Project("hTeO2co60M1", "Ener2", "Multiplicity==1");
	tTeO2k40->Project("hTeO2k40M1", "Ener2", "Multiplicity==1");
	tTeO2pb210->Project("hTeO2pb210M1", "Ener2", "Multiplicity==1");
	tTeO2po210->Project("hTeO2po210M1", "Ener2", "Multiplicity==1");
	tTeO2ra226only->Project("hTeO2ra226onlyM1", "Ener2", "Multiplicity==1");
	tTeO2ra226pb206->Project("hTeO2ra226pb206M1", "Ener2", "Multiplicity==1");
	tTeO2ra226pb210->Project("hTeO2ra226pb210M1", "Ener2", "Multiplicity==1");
	tTeO2ra228pb208->Project("hTeO2ra228pb208M1", "Ener2", "Multiplicity==1");
	tTeO2rn222only->Project("hTeO2rn222onlyM1", "Ener2", "Multiplicity==1");
	tTeO2rn222pb206->Project("hTeO2rn222pb206M1", "Ener2", "Multiplicity==1");
	tTeO2sb125->Project("hTeO2sb125M1", "Ener2", "Multiplicity==1");
	tTeO2te125->Project("hTeO2te125M1", "Ener2", "Multiplicity==1");
	tTeO2th230only->Project("hTeO2th230onlyM1", "Ener2", "Multiplicity==1");
	tTeO2th230pb206->Project("hTeO2th230pb206M1", "Ener2", "Multiplicity==1");
	tTeO2th232only->Project("hTeO2th232onlyM1", "Ener2", "Multiplicity==1");
	tTeO2th232->Project("hTeO2th232M1", "Ener2", "Multiplicity==1");
	tTeO2th232th228->Project("hTeO2th232th228M1", "Ener2", "Multiplicity==1");
	tTeO2u238ra226->Project("hTeO2u238ra226M1", "Ener2", "Multiplicity==1");
	tTeO2u238rn222->Project("hTeO2u238rn222M1", "Ener2", "Multiplicity==1");
	tTeO2u238->Project("hTeO2u238M1", "Ener2", "Multiplicity==1");
	tTeO2u238th230->Project("hTeO2u238th230M1", "Ener2", "Multiplicity==1");
	tTeO2u238u234->Project("hTeO2u238u234M1", "Ener2", "Multiplicity==1");


////////////////////////////////////////

	tCuBoxco58->Project("hCuBoxco58M2Sum", "Esum2", "Multiplicity==2");
	tCuBoxco60->Project("hCuBoxco60M2Sum", "Esum2", "Multiplicity==2");
	tCuBoxcs137->Project("hCuBoxcs137M2Sum", "Esum2", "Multiplicity==2");
	tCuBoxk40->Project("hCuBoxk40M2Sum", "Esum2", "Multiplicity==2");
	tCuBoxmn54->Project("hCuBoxmn54M2Sum", "Esum2", "Multiplicity==2");
	tCuBoxpb210->Project("hCuBoxpb210M2Sum", "Esum2", "Multiplicity==2");
	tCuBoxth232->Project("hCuBoxth232M2Sum", "Esum2", "Multiplicity==2");
	tCuBoxu238->Project("hCuBoxu238M2Sum", "Esum2", "Multiplicity==2");	
	tCuBoxeu152->Project("hCuBoxeu152M2Sum", "Esum2", "Multiplicity==2");
	tCuBoxeu154->Project("hCuBoxeu154M2Sum", "Esum2", "Multiplicity==2");

	tCuFrameco58->Project("hCuFrameco58M2Sum", "Esum2", "Multiplicity==2");
	tCuFrameco60->Project("hCuFrameco60M2Sum", "Esum2", "Multiplicity==2");
	tCuFramecs137->Project("hCuFramecs137M2Sum", "Esum2", "Multiplicity==2");
	tCuFramek40->Project("hCuFramek40M2Sum", "Esum2", "Multiplicity==2");
	tCuFramemn54->Project("hCuFramemn54M2Sum", "Esum2", "Multiplicity==2");
	tCuFramepb210->Project("hCuFramepb210M2Sum", "Esum2", "Multiplicity==2");
	tCuFrameth232->Project("hCuFrameth232M2Sum", "Esum2", "Multiplicity==2");
	tCuFrameu238->Project("hCuFrameu238M2Sum", "Esum2", "Multiplicity==2");

	tCuBox_CuFrameco60->Project("hCuBox_CuFrameco60M2Sum", "Esum2", "Multiplicity==2");
	tCuBox_CuFramek40->Project("hCuBox_CuFramek40M2Sum", "Esum2", "Multiplicity==2");
	tCuBox_CuFrameth232->Project("hCuBox_CuFrameth232M2Sum", "Esum2", "Multiplicity==2");
	tCuBox_CuFrameu238->Project("hCuBox_CuFrameu238M2Sum", "Esum2", "Multiplicity==2");

	tInternalco60->Project("hInternalco60M2Sum", "Esum2", "Multiplicity==2");
	tInternalk40->Project("hInternalk40M2Sum", "Esum2", "Multiplicity==2");
	tInternalth232->Project("hInternalth232M2Sum", "Esum2", "Multiplicity==2");
	tInternalu238->Project("hInternalu238M2Sum", "Esum2", "Multiplicity==2");

	tPTFEth232->Project("hPTFEth232M2Sum", "Esum2", "Multiplicity==2");
	tPTFEu238->Project("hPTFEu238M2Sum", "Esum2", "Multiplicity==2");

	tTeO20nu->Project("hTeO20nuM2Sum", "Esum2", "Multiplicity==2");
	tTeO22nu->Project("hTeO22nuM2Sum", "Esum2", "Multiplicity==2");
	tTeO2co60->Project("hTeO2co60M2Sum", "Esum2", "Multiplicity==2");
	tTeO2k40->Project("hTeO2k40M2Sum", "Esum2", "Multiplicity==2");
	tTeO2pb210->Project("hTeO2pb210M2Sum", "Esum2", "Multiplicity==2");
	tTeO2po210->Project("hTeO2po210M2Sum", "Esum2", "Multiplicity==2");
	tTeO2ra226only->Project("hTeO2ra226onlyM2Sum", "Esum2", "Multiplicity==2");
	tTeO2ra226pb206->Project("hTeO2ra226pb206M2Sum", "Esum2", "Multiplicity==2");
	tTeO2ra226pb210->Project("hTeO2ra226pb210M2Sum", "Esum2", "Multiplicity==2");
	tTeO2ra228pb208->Project("hTeO2ra228pb208M2Sum", "Esum2", "Multiplicity==2");
	tTeO2rn222only->Project("hTeO2rn222onlyM2Sum", "Esum2", "Multiplicity==2");
	tTeO2rn222pb206->Project("hTeO2rn222pb206M2Sum", "Esum2", "Multiplicity==2");
	tTeO2sb125->Project("hTeO2sb125M2Sum", "Esum2", "Multiplicity==2");
	tTeO2te125->Project("hTeO2te125M2Sum", "Esum2", "Multiplicity==2");
	tTeO2th230only->Project("hTeO2th230onlyM2Sum", "Esum2", "Multiplicity==2");
	tTeO2th230pb206->Project("hTeO2th230pb206M2Sum", "Esum2", "Multiplicity==2");
	tTeO2th232only->Project("hTeO2th232onlyM2Sum", "Esum2", "Multiplicity==2");
	tTeO2th232->Project("hTeO2th232M2Sum", "Esum2", "Multiplicity==2");
	tTeO2th232th228->Project("hTeO2th232th228M2Sum", "Esum2", "Multiplicity==2");
	tTeO2u238ra226->Project("hTeO2u238ra226M2Sum", "Esum2", "Multiplicity==2");
	tTeO2u238rn222->Project("hTeO2u238rn222M2Sum", "Esum2", "Multiplicity==2");
	tTeO2u238->Project("hTeO2u238M2Sum", "Esum2", "Multiplicity==2");
	tTeO2u238th230->Project("hTeO2u238th230M2Sum", "Esum2", "Multiplicity==2");
	tTeO2u238u234->Project("hTeO2u238u234M2Sum", "Esum2", "Multiplicity==2");


	//////////////////////////////////


	NormalizePDFPair( hCuBoxco58M1, hCuBoxco58M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxco60M1, hCuBoxco60M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxcs137M1, hCuBoxcs137M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxk40M1, hCuBoxk40M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxmn54M1, hCuBoxmn54M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxpb210M1, hCuBoxpb210M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxth232M1, hCuBoxth232M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxu238M1, hCuBoxu238M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxeu152M1, hCuBoxeu152M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBoxeu154M1, hCuBoxeu154M2Sum, dMinEnergy, dMaxEnergy);


	NormalizePDFPair( hCuFrameco58M1, hCuFrameco58M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFrameco60M1, hCuFrameco60M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFramecs137M1, hCuFramecs137M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFramek40M1, hCuFramek40M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFramemn54M1, hCuFramemn54M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFramepb210M1, hCuFramepb210M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFrameth232M1, hCuFrameth232M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuFrameu238M1, hCuFrameu238M2Sum, dMinEnergy, dMaxEnergy);

	NormalizePDFPair( hCuBox_CuFrameco60M1, hCuBox_CuFrameco60M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBox_CuFramek40M1, hCuBox_CuFramek40M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBox_CuFrameth232M1, hCuBox_CuFrameth232M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hCuBox_CuFrameu238M1, hCuBox_CuFrameu238M2Sum,  dMinEnergy, dMaxEnergy);

	NormalizePDFPair( hInternalco60M1, hInternalco60M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hInternalk40M1, hInternalk40M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hInternalth232M1, hInternalth232M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hInternalu238M1, hInternalu238M2Sum, dMinEnergy, dMaxEnergy);

	NormalizePDFPair( hPTFEth232M1, hPTFEth232M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPTFEu238M1, hPTFEu238M2Sum, dMinEnergy, dMaxEnergy);	

	NormalizePDFPair( hTeO20nuM1, hTeO20nuM2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO22nuM1, hTeO22nuM2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2co60M1, hTeO2co60M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2k40M1, hTeO2k40M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2pb210M1, hTeO2pb210M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2po210M1, hTeO2po210M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2ra226onlyM1, hTeO2ra226onlyM2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2ra226pb206M1, hTeO2ra226pb206M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2ra226pb210M1, hTeO2ra226pb210M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2ra228pb208M1, hTeO2ra228pb208M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2rn222onlyM1, hTeO2rn222onlyM2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2rn222pb206M1, hTeO2rn222pb206M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2sb125M1, hTeO2sb125M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2te125M1, hTeO2te125M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2th230onlyM1, hTeO2th230onlyM2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2th230pb206M1, hTeO2th230pb206M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2th232onlyM1, hTeO2th232onlyM2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2th232M1, hTeO2th232M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2th232th228M1, hTeO2th232th228M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2u238ra226M1, hTeO2u238ra226M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2u238rn222M1, hTeO2u238rn222M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2u238M1,  hTeO2u238M2Sum,   dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2u238th230M1, hTeO2u238th230M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hTeO2u238u234M1,  hTeO2u238u234M2Sum,   dMinEnergy, dMaxEnergy);


	TFile *file1 = new TFile("MCProduction_BulkInnerM2Sum_1keV.root", "RECREATE");

//////////////////////////////

	hCuBoxco58M2Sum->Write();
	hCuBoxco60M2Sum->Write();
	hCuBoxcs137M2Sum->Write();
	hCuBoxeu152M2Sum->Write();
	hCuBoxeu154M2Sum->Write();
	hCuBoxk40M2Sum->Write();
	hCuBoxmn54M2Sum->Write();
	hCuBoxpb210M2Sum->Write();
	hCuBoxth232M2Sum->Write();
	hCuBoxu238M2Sum->Write();	

	hCuFrameco58M2Sum->Write();
	hCuFrameco60M2Sum->Write();
	hCuFramecs137M2Sum->Write();
	hCuFramek40M2Sum->Write();
	hCuFramemn54M2Sum->Write();
	hCuFramepb210M2Sum->Write();
	hCuFrameth232M2Sum->Write();
	hCuFrameu238M2Sum->Write();

	hCuBox_CuFrameco60M2Sum->Write();
	hCuBox_CuFramek40M2Sum->Write();
	hCuBox_CuFrameth232M2Sum->Write();
	hCuBox_CuFrameu238M2Sum->Write();

	hInternalco60M2Sum->Write();
	hInternalk40M2Sum->Write();
	hInternalth232M2Sum->Write();
	hInternalu238M2Sum->Write();	

	hPTFEth232M2Sum->Write();
	hPTFEu238M2Sum->Write();

	hTeO20nuM2Sum->Write();
	hTeO22nuM2Sum->Write();
	hTeO2co60M2Sum->Write();
	hTeO2k40M2Sum->Write();
	hTeO2pb210M2Sum->Write();
	hTeO2po210M2Sum->Write();
	hTeO2ra226onlyM2Sum->Write();
	hTeO2ra226pb206M2Sum->Write();
	hTeO2ra226pb210M2Sum->Write();
	hTeO2ra228pb208M2Sum->Write();
	hTeO2rn222onlyM2Sum->Write();
	hTeO2rn222pb206M2Sum->Write();
	hTeO2sb125M2Sum->Write();
	hTeO2te125M2Sum->Write();
	hTeO2th230onlyM2Sum->Write();
	hTeO2th230pb206M2Sum->Write();	
	hTeO2th232onlyM2Sum->Write();
	hTeO2th232M2Sum->Write();
	hTeO2th232th228M2Sum->Write();
	hTeO2u238ra226M2Sum->Write();
	hTeO2u238rn222M2Sum->Write();
	hTeO2u238M2Sum->Write();
	hTeO2u238th230M2Sum->Write();
	hTeO2u238u234M2Sum->Write();

	file1->Write();	

	
}


void SaveHistogramsBulkM2SumOuter()
{
	std::string sDataDir = "/cuore/user/zhubrian/MC/Bkg/Production/";

	double dMinEnergy = 0;
	double dMaxEnergy = 10000;
	int dBinSize = 1;
	int dNBins = (dMaxEnergy - dMinEnergy)/dBinSize;


	TH1D *h50mKco58M1 = new TH1D("h50mKco58M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKco60M1 = new TH1D("h50mKco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKcs137M1 = new TH1D("h50mKcs137M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKk40M1 = new TH1D("h50mKk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKmn54M1 = new TH1D("h50mKmn54M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKpb210M1 = new TH1D("h50mKpb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKth232M1 = new TH1D("h50mKth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKu238M1 = new TH1D("h50mKu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *h600mKco60M1 = new TH1D("h600mKco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKk40M1 = new TH1D("h600mKk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKth232M1 = new TH1D("h600mKth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKu238M1 = new TH1D("h600mKu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hExtPbbi210M1 = new TH1D("hExtPbbi210M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hIVCco60M1 = new TH1D("hIVCco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCk40M1 = new TH1D("hIVCk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCth232M1 = new TH1D("hIVCth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCu238M1 = new TH1D("hIVCu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hMBco60M1 = new TH1D("hMBco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBk40M1 = new TH1D("hMBk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBth232M1 = new TH1D("hMBth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBu238M1 = new TH1D("hMBu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hMCco60M1 = new TH1D("hMCco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMCk40M1 = new TH1D("hMCk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMCth232M1 = new TH1D("hMCth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMCu238M1 = new TH1D("hMCu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hOVCbi207M1 = new TH1D("hOVCbi207M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCco60M1 = new TH1D("hOVCco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCk40M1 = new TH1D("hOVCk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCth232M1 = new TH1D("hOVCth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCu238M1 = new TH1D("hOVCu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	// TH1D *hPbRombi207M1 = new TH1D("hPbRombi207M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomco60M1 = new TH1D("hPbRomco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomcs137M1 = new TH1D("hPbRomcs137M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomeu152M1 = new TH1D("hPbRomeu152M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomeu154M1 = new TH1D("hPbRomeu154M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomk40M1 = new TH1D("hPbRomk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRompb210M1 = new TH1D("hPbRompb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomth232M1 = new TH1D("hPbRomth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomu238M1 = new TH1D("hPbRomu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hSIk40M1 = new TH1D("hSIk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hSIth232M1 = new TH1D("hSIth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hSIu238M1 = new TH1D("hSIu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

////////////////////////

	TH1D *h50mKco58M2Sum = new TH1D("h50mKco58M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKco60M2Sum = new TH1D("h50mKco60M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKcs137M2Sum = new TH1D("h50mKcs137M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKk40M2Sum = new TH1D("h50mKk40M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKmn54M2Sum = new TH1D("h50mKmn54M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKpb210M2Sum = new TH1D("h50mKpb210M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKth232M2Sum = new TH1D("h50mKth232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h50mKu238M2Sum = new TH1D("h50mKu238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *h600mKco60M2Sum = new TH1D("h600mKco60M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKk40M2Sum = new TH1D("h600mKk40M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKth232M2Sum = new TH1D("h600mKth232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *h600mKu238M2Sum = new TH1D("h600mKu238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hExtPbbi210M2Sum = new TH1D("hExtPbbi210M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hIVCco60M2Sum = new TH1D("hIVCco60M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCk40M2Sum = new TH1D("hIVCk40M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCth232M2Sum = new TH1D("hIVCth232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCu238M2Sum = new TH1D("hIVCu238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hMBco60M2Sum = new TH1D("hMBco60M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBk40M2Sum = new TH1D("hMBk40M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBth232M2Sum = new TH1D("hMBth232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBu238M2Sum = new TH1D("hMBu238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hMCco60M2Sum = new TH1D("hMCco60M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMCk40M2Sum = new TH1D("hMCk40M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMCth232M2Sum = new TH1D("hMCth232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMCu238M2Sum = new TH1D("hMCu238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hOVCbi207M2Sum = new TH1D("hOVCbi207M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCco60M2Sum = new TH1D("hOVCco60M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCk40M2Sum = new TH1D("hOVCk40M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCth232M2Sum = new TH1D("hOVCth232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCu238M2Sum = new TH1D("hOVCu238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);

	// TH1D *hPbRombi207M2Sum = new TH1D("hPbRombi207M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomco60M2Sum = new TH1D("hPbRomco60M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomcs137M2Sum = new TH1D("hPbRomcs137M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomeu152M2Sum = new TH1D("hPbRomeu152M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomeu154M2Sum = new TH1D("hPbRomeu154M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomk40M2Sum = new TH1D("hPbRomk40M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRompb210M2Sum = new TH1D("hPbRompb210M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomth232M2Sum = new TH1D("hPbRomth232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomu238M2Sum = new TH1D("hPbRomu238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hSIk40M2Sum = new TH1D("hSIk40M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hSIth232M2Sum = new TH1D("hSIth232M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hSIu238M2Sum = new TH1D("hSIu238M2Sum", "", dNBins, dMinEnergy, dMaxEnergy);



///////////////////////////////

	TChain *t50mKco58 = LoadMC(sDataDir.c_str(), "50mK", "co58");
	TChain *t50mKco60 = LoadMC(sDataDir.c_str(), "50mK", "co60");
	TChain *t50mKcs137 = LoadMC(sDataDir.c_str(), "50mK", "cs137");
	TChain *t50mKk40 = LoadMC(sDataDir.c_str(), "50mK", "k40");
	TChain *t50mKmn54 = LoadMC(sDataDir.c_str(), "50mK", "mn54");
	TChain *t50mKpb210 = LoadMC(sDataDir.c_str(), "50mK", "pb210");
	TChain *t50mKth232 = LoadMC(sDataDir.c_str(), "50mK", "th232");
	TChain *t50mKu238 = LoadMC(sDataDir.c_str(), "50mK", "u238");

	TChain *t600mKco60 = LoadMC(sDataDir.c_str(), "n7_600mK", "co60");
	TChain *t600mKk40 = LoadMC(sDataDir.c_str(), "600mK", "k40");
	TChain *t600mKth232 = LoadMC(sDataDir.c_str(), "n7_600mK", "th232");
	TChain *t600mKu238 = LoadMC(sDataDir.c_str(), "600mK", "u238");

	TChain *tExtPbbi210 = LoadMC(sDataDir.c_str(), "ExtPb4cm_210BiBStot", "1.23T");

	TChain *tIVCco60 = LoadMC(sDataDir.c_str(), "IVC", "co60");
	TChain *tIVCk40 = LoadMC(sDataDir.c_str(), "IVC", "k40");
	TChain *tIVCth232 = LoadMC(sDataDir.c_str(), "IVC", "th232");
	TChain *tIVCu238 = LoadMC(sDataDir.c_str(), "IVC", "u238");

	TChain *tOVCbi207 = LoadMC(sDataDir.c_str(), "OVC", "bi207");
	TChain *tOVCco60 = LoadMC(sDataDir.c_str(), "OVC", "co60");
	TChain *tOVCk40 = LoadMC(sDataDir.c_str(), "OVC", "k40");
	TChain *tOVCth232 = LoadMC(sDataDir.c_str(), "OVC", "th232");
	TChain *tOVCu238 = LoadMC(sDataDir.c_str(), "OVC", "u238");

	TChain *tSIk40 = LoadMC(sDataDir.c_str(), "SI", "k40");
	TChain *tSIth232 = LoadMC(sDataDir.c_str(), "SI", "th232");
	TChain *tSIu238 = LoadMC(sDataDir.c_str(), "SI", "u238");

	TChain *tMBco60 = LoadMC(sDataDir.c_str(), "MB", "co60");
	TChain *tMBk40 = LoadMC(sDataDir.c_str(), "MB", "k40");
	TChain *tMBth232 = LoadMC(sDataDir.c_str(), "MB", "th232");
	TChain *tMBu238 = LoadMC(sDataDir.c_str(), "MB", "u238");

	TChain *tMCco60 = LoadMC(sDataDir.c_str(), "MC", "co60");
	TChain *tMCk40 = LoadMC(sDataDir.c_str(), "MC", "k40");
	TChain *tMCth232 = LoadMC(sDataDir.c_str(), "MC", "th232");
	TChain *tMCu238 = LoadMC(sDataDir.c_str(), "MC", "u238");

	// TChain *tPbRombi207 = LoadMC(sDataDir.c_str(), "PbRom", "bi207");
	TChain *tPbRomco60 = LoadMC(sDataDir.c_str(), "PbRom", "co60");
	TChain *tPbRomcs137 = LoadMC(sDataDir.c_str(), "PbRom", "cs137");
	TChain *tPbRomeu152 = LoadMC(sDataDir.c_str(), "PbRom", "eu152");
	TChain *tPbRomeu154 = LoadMC(sDataDir.c_str(), "PbRom", "eu154");
	TChain *tPbRomk40 = LoadMC(sDataDir.c_str(), "PbRom", "k40");
	TChain *tPbRompb210 = LoadMC(sDataDir.c_str(), "PbRom", "pb210");
	TChain *tPbRomth232 = LoadMC(sDataDir.c_str(), "PbRom", "th232");
	TChain *tPbRomu238 = LoadMC(sDataDir.c_str(), "PbRom", "u238");

////////////////////////////////////////

	t50mKco58->Project("h50mKco58M1", "Ener2", "Multiplicity==1");
	t50mKco60->Project("h50mKco60M1", "Ener2", "Multiplicity==1");
	t50mKcs137->Project("h50mKcs137M1", "Ener2", "Multiplicity==1");
	t50mKk40->Project("h50mKk40M1", "Ener2", "Multiplicity==1");
	t50mKmn54->Project("h50mKmn54M1", "Ener2", "Multiplicity==1");
	t50mKpb210->Project("h50mKpb210M1", "Ener2", "Multiplicity==1");
	t50mKth232->Project("h50mKth232M1", "Ener2", "Multiplicity==1");
	t50mKu238->Project("h50mKu238M1", "Ener2", "Multiplicity==1");

	t600mKco60->Project("h600mKco60M1", "Ener2", "Multiplicity==1");
	t600mKk40->Project("h600mKk40M1", "Ener2", "Multiplicity==1");
	t600mKth232->Project("h600mKth232M1", "Ener2", "Multiplicity==1");
	t600mKu238->Project("h600mKu238M1", "Ener2", "Multiplicity==1");

	tExtPbbi210->Project("hExtPbbi210M1", "Ener2", "Multiplicity==1");

	tIVCco60->Project("hIVCco60M1", "Ener2", "Multiplicity==1");
	tIVCk40->Project("hIVCk40M1", "Ener2", "Multiplicity==1");
	tIVCth232->Project("hIVCth232M1", "Ener2", "Multiplicity==1");
	tIVCu238->Project("hIVCu238M1", "Ener2", "Multiplicity==1");

	tOVCbi207->Project("hOVCbi207M1", "Ener2", "Multiplicity==1");
	tOVCco60->Project("hOVCco60M1", "Ener2", "Multiplicity==1");
	tOVCk40->Project("hOVCk40M1", "Ener2", "Multiplicity==1");
	tOVCth232->Project("hOVCth232M1", "Ener2", "Multiplicity==1");
	tOVCu238->Project("hOVCu238M1", "Ener2", "Multiplicity==1");

	tMBco60->Project("hMBco60M1", "Ener2", "Multiplicity==1");
	tMBk40->Project("hMBk40M1", "Ener2", "Multiplicity==1");
	tMBth232->Project("hMBth232M1", "Ener2", "Multiplicity==1");
	tMBu238->Project("hMBu238M1", "Ener2", "Multiplicity==1");

	tMCco60->Project("hMCco60M1", "Ener2", "Multiplicity==1");
	tMCk40->Project("hMCk40M1", "Ener2", "Multiplicity==1");
	tMCth232->Project("hMCth232M1", "Ener2", "Multiplicity==1");
	tMCu238->Project("hMCu238M1", "Ener2", "Multiplicity==1");

	tSIk40->Project("hSIk40M1", "Ener2", "Multiplicity==1");
	tSIth232->Project("hSIth232M1", "Ener2", "Multiplicity==1");
	tSIu238->Project("hSIu238M1", "Ener2", "Multiplicity==1");

	// tPbRombi207->Project("hPbRombi207M1", "Ener2", "Multiplicity==1");
	tPbRomco60->Project("hPbRomco60M1", "Ener2", "Multiplicity==1");
	tPbRomcs137->Project("hPbRomcs137M1", "Ener2", "Multiplicity==1");
	tPbRomeu152->Project("hPbRomeu152M1", "Ener2", "Multiplicity==1");
	tPbRomeu154->Project("hPbRomeu154M1", "Ener2", "Multiplicity==1");
	tPbRomk40->Project("hPbRomk40M1", "Ener2", "Multiplicity==1");
	tPbRompb210->Project("hPbRompb210M1", "Ener2", "Multiplicity==1");
	tPbRomth232->Project("hPbRomth232M1", "Ener2", "Multiplicity==1");
	tPbRomu238->Project("hPbRomu238M1", "Ener2", "Multiplicity==1");


////////////////////////////////////////

	t50mKco58->Project("h50mKco58M2Sum", "Esum2", "Multiplicity==2");
	t50mKco60->Project("h50mKco60M2Sum", "Esum2", "Multiplicity==2");
	t50mKcs137->Project("h50mKcs137M2Sum", "Esum2", "Multiplicity==2");
	t50mKk40->Project("h50mKk40M2Sum", "Esum2", "Multiplicity==2");
	t50mKmn54->Project("h50mKmn54M2Sum", "Esum2", "Multiplicity==2");
	t50mKpb210->Project("h50mKpb210M2Sum", "Esum2", "Multiplicity==2");
	t50mKth232->Project("h50mKth232M2Sum", "Esum2", "Multiplicity==2");
	t50mKu238->Project("h50mKu238M2Sum", "Esum2", "Multiplicity==2");

	t600mKco60->Project("h600mKco60M2Sum", "Esum2", "Multiplicity==2");
	t600mKk40->Project("h600mKk40M2Sum", "Esum2", "Multiplicity==2");
	t600mKth232->Project("h600mKth232M2Sum", "Esum2", "Multiplicity==2");
	t600mKu238->Project("h600mKu238M2Sum", "Esum2", "Multiplicity==2");

	tExtPbbi210->Project("hExtPbbi210M2Sum", "Esum2", "Multiplicity==2");

	tIVCco60->Project("hIVCco60M2Sum", "Esum2", "Multiplicity==2");
	tIVCk40->Project("hIVCk40M2Sum", "Esum2", "Multiplicity==2");
	tIVCth232->Project("hIVCth232M2Sum", "Esum2", "Multiplicity==2");
	tIVCu238->Project("hIVCu238M2Sum", "Esum2", "Multiplicity==2");

	tOVCbi207->Project("hOVCbi207M2Sum", "Esum2", "Multiplicity==2");
	tOVCco60->Project("hOVCco60M2Sum", "Esum2", "Multiplicity==2");
	tOVCk40->Project("hOVCk40M2Sum", "Esum2", "Multiplicity==2");
	tOVCth232->Project("hOVCth232M2Sum", "Esum2", "Multiplicity==2");
	tOVCu238->Project("hOVCu238M2Sum", "Esum2", "Multiplicity==2");

	tMBco60->Project("hMBco60M2Sum", "Esum2", "Multiplicity==2");
	tMBk40->Project("hMBk40M2Sum", "Esum2", "Multiplicity==2");
	tMBth232->Project("hMBth232M2Sum", "Esum2", "Multiplicity==2");
	tMBu238->Project("hMBu238M2Sum", "Esum2", "Multiplicity==2");

	tMCco60->Project("hMCco60M2Sum", "Esum2", "Multiplicity==2");
	tMCk40->Project("hMCk40M2Sum", "Esum2", "Multiplicity==2");
	tMCth232->Project("hMCth232M2Sum", "Esum2", "Multiplicity==2");
	tMCu238->Project("hMCu238M2Sum", "Esum2", "Multiplicity==2");

	tSIk40->Project("hSIk40M2Sum", "Esum2", "Multiplicity==2");
	tSIth232->Project("hSIth232M2Sum", "Esum2", "Multiplicity==2");
	tSIu238->Project("hSIu238M2Sum", "Esum2", "Multiplicity==2");

	// tPbRombi207->Project("hPbRombi207M2Sum", "Esum2", "Multiplicity==2");
	tPbRomco60->Project("hPbRomco60M2Sum", "Esum2", "Multiplicity==2");
	tPbRomcs137->Project("hPbRomcs137M2Sum", "Esum2", "Multiplicity==2");
	tPbRomeu152->Project("hPbRomeu152M2Sum", "Esum2", "Multiplicity==2");
	tPbRomeu154->Project("hPbRomeu154M2Sum", "Esum2", "Multiplicity==2");
	tPbRomk40->Project("hPbRomk40M2Sum", "Esum2", "Multiplicity==2");
	tPbRompb210->Project("hPbRompb210M2Sum", "Esum2", "Multiplicity==2");
	tPbRomth232->Project("hPbRomth232M2Sum", "Esum2", "Multiplicity==2");
	tPbRomu238->Project("hPbRomu238M2Sum", "Esum2", "Multiplicity==2");


	//////////////////////////////////


	NormalizePDFPair( h50mKco58M1, h50mKco58M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKco60M1, h50mKco60M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKcs137M1, h50mKcs137M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKk40M1, h50mKk40M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKmn54M1, h50mKmn54M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKpb210M1, h50mKpb210M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKth232M1, h50mKth232M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h50mKu238M1, h50mKu238M2Sum, dMinEnergy, dMaxEnergy);

	NormalizePDFPair( h600mKco60M1, h600mKco60M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h600mKk40M1, h600mKk40M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h600mKth232M1, h600mKth232M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( h600mKu238M1, h600mKu238M2Sum, dMinEnergy, dMaxEnergy);

	NormalizePDFPair( hIVCco60M1, hIVCco60M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hIVCk40M1, hIVCk40M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hIVCth232M1, hIVCth232M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hIVCu238M1, hIVCu238M2Sum, dMinEnergy, dMaxEnergy);

	NormalizePDFPair( hOVCbi207M1, hOVCbi207M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hOVCco60M1, hOVCco60M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hOVCk40M1, hOVCk40M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hOVCth232M1, hOVCth232M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hOVCu238M1, hOVCu238M2Sum, dMinEnergy, dMaxEnergy);

	NormalizePDFPair( hMBco60M1, hMBco60M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hMBk40M1, hMBk40M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hMBth232M1, hMBth232M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hMBu238M1, hMBu238M2Sum, dMinEnergy, dMaxEnergy);

	NormalizePDFPair( hMCco60M1, hMCco60M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hMCk40M1, hMCk40M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hMCth232M1, hMCth232M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hMCu238M1, hMCu238M2Sum, dMinEnergy, dMaxEnergy);

	NormalizePDFPair( hSIk40M1, hSIk40M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hSIth232M1, hSIth232M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hSIu238M1, hSIu238M2Sum, dMinEnergy, dMaxEnergy);


	// NormalizePDFPair( hPbRombi207M1, hPbRombi207M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomco60M1, hPbRomco60M2Sum,  dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomcs137M1, hPbRomcs137M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomeu152M1, hPbRomeu152M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomeu154M1, hPbRomeu154M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomk40M1, hPbRomk40M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRompb210M1, hPbRompb210M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomth232M1, hPbRomth232M2Sum, dMinEnergy, dMaxEnergy);
	NormalizePDFPair( hPbRomu238M1, hPbRomu238M2Sum, dMinEnergy, dMaxEnergy);



	TFile *file1 = new TFile("MCProduction_BulkOuterM2Sum_1keV.root", "RECREATE");

//////////////////////////////

	h50mKco58M2Sum->Write();
	h50mKco60M2Sum->Write();
	h50mKcs137M2Sum->Write();
	h50mKk40M2Sum->Write();
	h50mKmn54M2Sum->Write();
	h50mKpb210M2Sum->Write();
	h50mKth232M2Sum->Write();
	h50mKu238M2Sum->Write();

	h600mKco60M2Sum->Write();
	h600mKk40M2Sum->Write();
	h600mKth232M2Sum->Write();
	h600mKu238M2Sum->Write();

	hExtPbbi210M2Sum->Write();

	hIVCco60M2Sum->Write();
	hIVCk40M2Sum->Write();
	hIVCth232M2Sum->Write();
	hIVCu238M2Sum->Write();

	hOVCbi207M2Sum->Write();
	hOVCco60M2Sum->Write();
	hOVCk40M2Sum->Write();
	hOVCth232M2Sum->Write();
	hOVCu238M2Sum->Write();

	hMBco60M2Sum->Write();
	hMBk40M2Sum->Write();
	hMBth232M2Sum->Write();
	hMBu238M2Sum->Write();

	hMCco60M2Sum->Write();
	hMCk40M2Sum->Write();
	hMCth232M2Sum->Write();
	hMCu238M2Sum->Write();

	hSIk40M2Sum->Write();
	hSIth232M2Sum->Write();
	hSIu238M2Sum->Write();

	// hPbRombi207M2Sum->Write();
	hPbRomco60M2Sum->Write();
	hPbRomcs137M2Sum->Write();
	hPbRomeu152M2Sum->Write();
	hPbRomeu154M2Sum->Write();
	hPbRomk40M2Sum->Write();
	hPbRompb210M2Sum->Write();
	hPbRomth232M2Sum->Write();
	hPbRomu238M2Sum->Write();

	file1->Write();
	
}





void SaveHistogramsSurfaceCrystal()
{

	std::string sDataDir = "/cuore/data/simulation/CUORE0/t14.08/production_g2root-r363/ntp/Sup/";

	double dMinEnergy = 0;
	double dMaxEnergy = 10000;
	int dBinSize = 1;
	int dNBins = (dMaxEnergy - dMinEnergy)/dBinSize;

	// TH1D *hTeO2Spb210M1_01 = new TH1D("hTeO2Spb210M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Spo210M1_001 = new TH1D("hTeO2Spo210M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Spo210M1_01 = new TH1D("hTeO2Spo210M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Sth232M1_01 = new TH1D("hTeO2Sth232M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Su238M1_01 = new TH1D("hTeO2Su238M1_01", "", dNBins, dMinEnergy, dMaxEnergy);


	TH1D *hTeO2Sxpb210M1_001 = new TH1D("hTeO2Sxpb210M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M1_01 = new TH1D("hTeO2Sxpb210M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M1_1 = new TH1D("hTeO2Sxpb210M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M1_10 = new TH1D("hTeO2Sxpb210M1_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpo210M1_001 = new TH1D("hTeO2Sxpo210M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpo210M1_01 = new TH1D("hTeO2Sxpo210M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpo210M1_1 = new TH1D("hTeO2Sxpo210M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Sxpo210M1_10 = new TH1D("hTeO2Sxpo210M1_10", "", dNBins, dMinEnergy, dMaxEnergy);	
	TH1D *hTeO2Sxth232M1_001 = new TH1D("hTeO2Sxth232M1_001", "", dNBins, dMinEnergy, dMaxEnergy);	
	TH1D *hTeO2Sxth232M1_01 = new TH1D("hTeO2Sxth232M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth232M1_1 = new TH1D("hTeO2Sxth232M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth232M1_10 = new TH1D("hTeO2Sxth232M1_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238M1_001 = new TH1D("hTeO2Sxu238M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238M1_01 = new TH1D("hTeO2Sxu238M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238M1_1 = new TH1D("hTeO2Sxu238M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238M1_10 = new TH1D("hTeO2Sxu238M1_10", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hTeO2Sxth232onlyM1_001 = new TH1D("hTeO2Sxth232onlyM1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxra228pb208M1_001 = new TH1D("hTeO2Sxra228pb208M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238th230M1_001 = new TH1D("hTeO2Sxu238th230M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth230onlyM1_001 = new TH1D("hTeO2Sxth230onlyM1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxra226pb210M1_001 = new TH1D("hTeO2Sxra226pb210M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M1_0001 = new TH1D("hTeO2Sxpb210M1_0001", "", dNBins, dMinEnergy, dMaxEnergy);

//////////////////////////////

	// TH1D *hTeO2Spb210M2_01 = new TH1D("hTeO2Spb210M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Spo210M2_001 = new TH1D("hTeO2Spo210M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Spo210M2_01 = new TH1D("hTeO2Spo210M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Sth232M2_01 = new TH1D("hTeO2Sth232M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Su238M2_01 = new TH1D("hTeO2Su238M2_01", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hTeO2Sxpb210M2_001 = new TH1D("hTeO2Sxpb210M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M2_01 = new TH1D("hTeO2Sxpb210M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M2_1 = new TH1D("hTeO2Sxpb210M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M2_10 = new TH1D("hTeO2Sxpb210M2_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpo210M2_001 = new TH1D("hTeO2Sxpo210M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpo210M2_01 = new TH1D("hTeO2Sxpo210M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpo210M2_1 = new TH1D("hTeO2Sxpo210M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Sxpo210M2_10 = new TH1D("hTeO2Sxpo210M2_10", "", dNBins, dMinEnergy, dMaxEnergy);	
	TH1D *hTeO2Sxth232M2_001 = new TH1D("hTeO2Sxth232M2_001", "", dNBins, dMinEnergy, dMaxEnergy);	
	TH1D *hTeO2Sxth232M2_01 = new TH1D("hTeO2Sxth232M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth232M2_1 = new TH1D("hTeO2Sxth232M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth232M2_10 = new TH1D("hTeO2Sxth232M2_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238M2_001 = new TH1D("hTeO2Sxu238M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238M2_01 = new TH1D("hTeO2Sxu238M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238M2_1 = new TH1D("hTeO2Sxu238M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238M2_10 = new TH1D("hTeO2Sxu238M2_10", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hTeO2Sxth232onlyM2_001 = new TH1D("hTeO2Sxth232onlyM2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxra228pb208M2_001 = new TH1D("hTeO2Sxra228pb208M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238th230M2_001 = new TH1D("hTeO2Sxu238th230M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth230onlyM2_001 = new TH1D("hTeO2Sxth230onlyM2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxra226pb210M2_001 = new TH1D("hTeO2Sxra226pb210M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M2_0001 = new TH1D("hTeO2Sxpb210M2_0001", "", dNBins, dMinEnergy, dMaxEnergy);

//////////////////////////////

	// TH1D *hTeO2Spb210M2Sum_01 = new TH1D("hTeO2Spb210M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Spo210M2Sum_001 = new TH1D("hTeO2Spo210M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Spo210M2Sum_01 = new TH1D("hTeO2Spo210M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Sth232M2Sum_01 = new TH1D("hTeO2Sth232M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Su238M2Sum_01 = new TH1D("hTeO2Su238M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hTeO2Sxpb210M2Sum_001 = new TH1D("hTeO2Sxpb210M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M2Sum_01 = new TH1D("hTeO2Sxpb210M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M2Sum_1 = new TH1D("hTeO2Sxpb210M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M2Sum_10 = new TH1D("hTeO2Sxpb210M2Sum_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpo210M2Sum_001 = new TH1D("hTeO2Sxpo210M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpo210M2Sum_01 = new TH1D("hTeO2Sxpo210M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpo210M2Sum_1 = new TH1D("hTeO2Sxpo210M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);
	// TH1D *hTeO2Sxpo210M2Sum_10 = new TH1D("hTeO2Sxpo210M2Sum_10", "", dNBins, dMinEnergy, dMaxEnergy);	
	TH1D *hTeO2Sxth232M2Sum_001 = new TH1D("hTeO2Sxth232M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);	
	TH1D *hTeO2Sxth232M2Sum_01 = new TH1D("hTeO2Sxth232M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth232M2Sum_1 = new TH1D("hTeO2Sxth232M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth232M2Sum_10 = new TH1D("hTeO2Sxth232M2Sum_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238M2Sum_001 = new TH1D("hTeO2Sxu238M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238M2Sum_01 = new TH1D("hTeO2Sxu238M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238M2Sum_1 = new TH1D("hTeO2Sxu238M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238M2Sum_10 = new TH1D("hTeO2Sxu238M2Sum_10", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hTeO2Sxth232onlyM2Sum_001 = new TH1D("hTeO2Sxth232onlyM2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxra228pb208M2Sum_001 = new TH1D("hTeO2Sxra228pb208M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxu238th230M2Sum_001 = new TH1D("hTeO2Sxu238th230M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxth230onlyM2Sum_001 = new TH1D("hTeO2Sxth230onlyM2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxra226pb210M2Sum_001 = new TH1D("hTeO2Sxra226pb210M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sxpb210M2Sum_0001 = new TH1D("hTeO2Sxpb210M2Sum_0001", "", dNBins, dMinEnergy, dMaxEnergy);

	// TChain *tTeO2Spb210_01 = LoadMC(sDataDir.c_str(), "TeO2S", "pb210-.1");
	// TChain *tTeO2Spo210_001 = LoadMC(sDataDir.c_str(), "TeO2S", "po210-.01");
	// TChain *tTeO2Spo210_01 = LoadMC(sDataDir.c_str(), "TeO2S", "po210-.1");
	// TChain *tTeO2Sth232_01 = LoadMC(sDataDir.c_str(), "TeO2S", "th232-.1");
	// TChain *tTeO2Su238_01 = LoadMC(sDataDir.c_str(), "TeO2S", "u238-.1");


	TChain *tTeO2Sxpb210_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "pb210-.01");
	TChain *tTeO2Sxpb210_01 = LoadMC(sDataDir.c_str(), "TeO2Sx", "pb210-.1");
	TChain *tTeO2Sxpb210_1 = LoadMC(sDataDir.c_str(), "TeO2Sx", "pb210-1");
	TChain *tTeO2Sxpb210_10 = LoadMC(sDataDir.c_str(), "TeO2Sx", "pb210-10");
	TChain *tTeO2Sxpo210_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "po210-.01");
	TChain *tTeO2Sxpo210_01 = LoadMC(sDataDir.c_str(), "TeO2Sx", "po210-.1");
	TChain *tTeO2Sxpo210_1 = LoadMC(sDataDir.c_str(), "TeO2Sx", "po210-1");
	// TChain *tTeO2Sxpo210_10 = LoadMC(sDataDir.c_str(), "TeO2Sx", "po210-10");	
	TChain *tTeO2Sxth232_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "th232-.01");
	TChain *tTeO2Sxth232_01 = LoadMC(sDataDir.c_str(), "TeO2Sx", "th232-.1");
	TChain *tTeO2Sxth232_1 = LoadMC(sDataDir.c_str(), "TeO2Sx", "th232-1");
	TChain *tTeO2Sxth232_10 = LoadMC(sDataDir.c_str(), "TeO2Sx", "th232-10");
	TChain *tTeO2Sxu238_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "u238-.01");
	TChain *tTeO2Sxu238_01 = LoadMC(sDataDir.c_str(), "TeO2Sx", "u238-.1");
	TChain *tTeO2Sxu238_1 = LoadMC(sDataDir.c_str(), "TeO2Sx", "u238-1");
	TChain *tTeO2Sxu238_10 = LoadMC(sDataDir.c_str(), "TeO2Sx", "u238-10");	

	TChain *tTeO2Sxth232only_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "th232_only-.01");
	TChain *tTeO2Sxra228pb208_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "ra228-pb208-.01");
	TChain *tTeO2Sxu238th230_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "u238-th230-.01");
	TChain *tTeO2Sxth230only_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "th230_only-.01");
	TChain *tTeO2Sxra226pb210_001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "ra226-pb210-.01");
	TChain *tTeO2Sxpb210_0001 = LoadMC(sDataDir.c_str(), "TeO2Sx", "pb210-.001");


///////////////////

	// tTeO2Spb210_01->Project("hTeO2Spb210M1_01", "Ener2", "Multiplicity==1");
	// tTeO2Spo210_001->Project("hTeO2Spo210M1_001", "Ener2", "Multiplicity==1");
	// tTeO2Spo210_01->Project("hTeO2Spo210M1_01", "Ener2", "Multiplicity==1");
	// tTeO2Sth232_01->Project("hTeO2Sth232M1_01", "Ener2", "Multiplicity==1");
	// tTeO2Su238_01->Project("hTeO2Su238M1_01", "Ener2", "Multiplicity==1");

	tTeO2Sxpb210_001->Project("hTeO2Sxpb210M1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxpb210_01->Project("hTeO2Sxpb210M1_01", "Ener2", "Multiplicity==1");
	tTeO2Sxpb210_1->Project("hTeO2Sxpb210M1_1", "Ener2", "Multiplicity==1");
	tTeO2Sxpb210_10->Project("hTeO2Sxpb210M1_10", "Ener2", "Multiplicity==1");
	tTeO2Sxpo210_001->Project("hTeO2Sxpo210M1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxpo210_01->Project("hTeO2Sxpo210M1_01", "Ener2", "Multiplicity==1");
	tTeO2Sxpo210_1->Project("hTeO2Sxpo210M1_1", "Ener2", "Multiplicity==1");
	// tTeO2Sxpo210_10->Project("hTeO2Sxpo210M1_10", "Ener2", "Multiplicity==1");
	tTeO2Sxth232_001->Project("hTeO2Sxth232M1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxth232_01->Project("hTeO2Sxth232M1_01", "Ener2", "Multiplicity==1");
	tTeO2Sxth232_1->Project("hTeO2Sxth232M1_1", "Ener2", "Multiplicity==1");
	tTeO2Sxth232_10->Project("hTeO2Sxth232M1_10", "Ener2", "Multiplicity==1");
	tTeO2Sxu238_001->Project("hTeO2Sxu238M1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxu238_01->Project("hTeO2Sxu238M1_01", "Ener2", "Multiplicity==1");
	tTeO2Sxu238_1->Project("hTeO2Sxu238M1_1", "Ener2", "Multiplicity==1");
	tTeO2Sxu238_10->Project("hTeO2Sxu238M1_10", "Ener2", "Multiplicity==1");

	tTeO2Sxth232only_001->Project("hTeO2Sxth232onlyM1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxra228pb208_001->Project("hTeO2Sxra228pb208M1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxu238th230_001->Project("hTeO2Sxu238th230M1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxth230only_001->Project("hTeO2Sxth230onlyM1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxra226pb210_001->Project("hTeO2Sxra226pb210M1_001", "Ener2", "Multiplicity==1");
	tTeO2Sxpb210_0001->Project("hTeO2Sxpb210M1_0001", "Ener2", "Multiplicity==1");

//////////////////////////////////

	// tTeO2Spb210_01->Project("hTeO2Spb210M2_01", "Ener2", "Multiplicity==2");
	// tTeO2Spo210_001->Project("hTeO2Spo210M2_001", "Ener2", "Multiplicity==2");
	// tTeO2Spo210_01->Project("hTeO2Spo210M2_01", "Ener2", "Multiplicity==2");
	// tTeO2Sth232_01->Project("hTeO2Sth232M2_01", "Ener2", "Multiplicity==2");
	// tTeO2Su238_01->Project("hTeO2Su238M2_01", "Ener2", "Multiplicity==2");

	tTeO2Sxpb210_001->Project("hTeO2Sxpb210M2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxpb210_01->Project("hTeO2Sxpb210M2_01", "Ener2", "Multiplicity==2");
	tTeO2Sxpb210_1->Project("hTeO2Sxpb210M2_1", "Ener2", "Multiplicity==2");
	tTeO2Sxpb210_10->Project("hTeO2Sxpb210M2_10", "Ener2", "Multiplicity==2");
	tTeO2Sxpo210_001->Project("hTeO2Sxpo210M2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxpo210_01->Project("hTeO2Sxpo210M2_01", "Ener2", "Multiplicity==2");
	tTeO2Sxpo210_1->Project("hTeO2Sxpo210M2_1", "Ener2", "Multiplicity==2");
	// tTeO2Sxpo210_10->Project("hTeO2Sxpo210M2_10", "Ener2", "Multiplicity==2");
	tTeO2Sxth232_001->Project("hTeO2Sxth232M2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxth232_01->Project("hTeO2Sxth232M2_01", "Ener2", "Multiplicity==2");
	tTeO2Sxth232_1->Project("hTeO2Sxth232M2_1", "Ener2", "Multiplicity==2");
	tTeO2Sxth232_10->Project("hTeO2Sxth232M2_10", "Ener2", "Multiplicity==2");
	tTeO2Sxu238_001->Project("hTeO2Sxu238M2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxu238_01->Project("hTeO2Sxu238M2_01", "Ener2", "Multiplicity==2");
	tTeO2Sxu238_1->Project("hTeO2Sxu238M2_1", "Ener2", "Multiplicity==2");
	tTeO2Sxu238_10->Project("hTeO2Sxu238M2_10", "Ener2", "Multiplicity==2");

	tTeO2Sxth232only_001->Project("hTeO2Sxth232onlyM2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxra228pb208_001->Project("hTeO2Sxra228pb208M2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxu238th230_001->Project("hTeO2Sxu238th230M2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxth230only_001->Project("hTeO2Sxth230onlyM2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxra226pb210_001->Project("hTeO2Sxra226pb210M2_001", "Ener2", "Multiplicity==2");
	tTeO2Sxpb210_0001->Project("hTeO2Sxpb210M2_0001", "Ener2", "Multiplicity==2");

//////////////////////////////////

	// tTeO2Spb210_01->Project("hTeO2Spb210M2Sum_01", "Esum2", "Multiplicity==2");
	// tTeO2Spo210_001->Project("hTeO2Spo210M2Sum_001", "Esum2", "Multiplicity==2");
	// tTeO2Spo210_01->Project("hTeO2Spo210M2Sum_01", "Esum2", "Multiplicity==2");
	// tTeO2Sth232_01->Project("hTeO2Sth232M2Sum_01", "Esum2", "Multiplicity==2");
	// tTeO2Su238_01->Project("hTeO2Su238M2Sum_01", "Esum2", "Multiplicity==2");

	// tTeO2Sxpb210_001->Project("hTeO2Sxpb210M2Sum_001", "Esum2", "Multiplicity==2");
	// tTeO2Sxpb210_01->Project("hTeO2Sxpb210M2Sum_01", "Esum2", "Multiplicity==2");
	// tTeO2Sxpb210_1->Project("hTeO2Sxpb210M2Sum_1", "Esum2", "Multiplicity==2");
	// tTeO2Sxpb210_10->Project("hTeO2Sxpb210M2Sum_10", "Esum2", "Multiplicity==2");
	// tTeO2Sxpo210_001->Project("hTeO2Sxpo210M2Sum_001", "Esum2", "Multiplicity==2");
	// tTeO2Sxpo210_01->Project("hTeO2Sxpo210M2Sum_01", "Esum2", "Multiplicity==2");
	// tTeO2Sxpo210_1->Project("hTeO2Sxpo210M2Sum_1", "Esum2", "Multiplicity==2");
	// // tTeO2Sxpo210_10->Project("hTeO2Sxpo210M2Sum_10", "Esum2", "Multiplicity==2");
	// tTeO2Sxth232_001->Project("hTeO2Sxth232M2Sum_001", "Esum2", "Multiplicity==2");
	// tTeO2Sxth232_01->Project("hTeO2Sxth232M2Sum_01", "Esum2", "Multiplicity==2");
	// tTeO2Sxth232_1->Project("hTeO2Sxth232M2Sum_1", "Esum2", "Multiplicity==2");
	// tTeO2Sxth232_10->Project("hTeO2Sxth232M2Sum_10", "Esum2", "Multiplicity==2");
	// tTeO2Sxu238_001->Project("hTeO2Sxu238M2Sum_001", "Esum2", "Multiplicity==2");
	// tTeO2Sxu238_01->Project("hTeO2Sxu238M2Sum_01", "Esum2", "Multiplicity==2");
	// tTeO2Sxu238_1->Project("hTeO2Sxu238M2Sum_1", "Esum2", "Multiplicity==2");
	// tTeO2Sxu238_10->Project("hTeO2Sxu238M2Sum_10", "Esum2", "Multiplicity==2");

	// tTeO2Sxth232only_001->Project("hTeO2Sxth232onlyM2Sum_001", "Esum2", "Multiplicity==2");
	// tTeO2Sxra228pb208_001->Project("hTeO2Sxra228pb208M2Sum_001", "Esum2", "Multiplicity==2");
	// tTeO2Sxu238th230_001->Project("hTeO2Sxu238th230M2Sum_001", "Esum2", "Multiplicity==2");
	// tTeO2Sxth230only_001->Project("hTeO2Sxth230onlyM2Sum_001", "Esum2", "Multiplicity==2");
	// tTeO2Sxra226pb210_001->Project("hTeO2Sxra226pb210M2Sum_001", "Esum2", "Multiplicity==2");
	// tTeO2Sxpb210_0001->Project("hTeO2Sxpb210M2Sum_0001", "Esum2", "Multiplicity==2");


////////////////////////////

	// NormalizePDFs(hTeO2Spb210M1_01, hTeO2Spb210M2_01, hTeO2Spb210M2Sum_01, 0, 10000);
	// NormalizePDFs(hTeO2Spo210M1_001, hTeO2Spo210M2_001, hTeO2Spo210M2Sum_001, 0, 10000);
	// NormalizePDFs(hTeO2Spo210M1_01, hTeO2Spo210M2_01, hTeO2Spo210M2Sum_01, 0, 10000);
	// NormalizePDFs(hTeO2Sth232M1_01, hTeO2Sth232M2_01, hTeO2Sth232M2Sum_01, 0, 10000);
	// NormalizePDFs(hTeO2Su238M1_01, hTeO2Su238M2_01, hTeO2Su238M2Sum_01, 0, 10000);

	NormalizePDFs(hTeO2Sxpb210M1_001, hTeO2Sxpb210M2_001, hTeO2Sxpb210M2Sum_001, 0, 10000);
	NormalizePDFs(hTeO2Sxpb210M1_01, hTeO2Sxpb210M2_01, hTeO2Sxpb210M2Sum_01, 0, 10000);
	NormalizePDFs(hTeO2Sxpb210M1_1, hTeO2Sxpb210M2_1, hTeO2Sxpb210M2Sum_1, 0, 10000);
	NormalizePDFs(hTeO2Sxpb210M1_10, hTeO2Sxpb210M2_10, hTeO2Sxpb210M2Sum_10, 0, 10000);
	NormalizePDFs(hTeO2Sxpo210M1_001, hTeO2Sxpo210M2_001, hTeO2Sxpo210M2Sum_001, 0, 10000);
	NormalizePDFs(hTeO2Sxpo210M1_01, hTeO2Sxpo210M2_01, hTeO2Sxpo210M2Sum_01, 0, 10000);
	NormalizePDFs(hTeO2Sxpo210M1_1, hTeO2Sxpo210M2_1, hTeO2Sxpo210M2Sum_1, 0, 10000);
	// NormalizePDFs(hTeO2Sxpo210M1_10, hTeO2Sxpo210M2_10, , 0, 10000);
	NormalizePDFs(hTeO2Sxth232M1_001, hTeO2Sxth232M2_001, hTeO2Sxth232M2Sum_001, 0, 10000);
	NormalizePDFs(hTeO2Sxth232M1_01, hTeO2Sxth232M2_01, hTeO2Sxth232M2Sum_01, 0, 10000);
	NormalizePDFs(hTeO2Sxth232M1_1, hTeO2Sxth232M2_1, hTeO2Sxth232M2Sum_1, 0, 10000);
	NormalizePDFs(hTeO2Sxth232M1_10, hTeO2Sxth232M2_10, hTeO2Sxth232M2Sum_10, 0, 10000);
	NormalizePDFs(hTeO2Sxu238M1_001, hTeO2Sxu238M2_001, hTeO2Sxu238M2Sum_001, 0, 10000);
	NormalizePDFs(hTeO2Sxu238M1_01, hTeO2Sxu238M2_01, hTeO2Sxu238M2Sum_01, 0, 10000);
	NormalizePDFs(hTeO2Sxu238M1_1, hTeO2Sxu238M2_1, hTeO2Sxu238M2Sum_1, 0, 10000);
	NormalizePDFs(hTeO2Sxu238M1_10, hTeO2Sxu238M2_10, hTeO2Sxu238M2Sum_10, 0, 10000);	

	NormalizePDFs(hTeO2Sxth232onlyM1_001, hTeO2Sxth232onlyM2_001, hTeO2Sxth232onlyM2Sum_001, 0, 10000);
	NormalizePDFs(hTeO2Sxra228pb208M1_001, hTeO2Sxra228pb208M2_001, hTeO2Sxra228pb208M2Sum_001, 0, 10000);
	NormalizePDFs(hTeO2Sxu238th230M1_001, hTeO2Sxu238th230M2_001, hTeO2Sxu238th230M2Sum_001, 0, 10000);
	NormalizePDFs(hTeO2Sxra226pb210M1_001, hTeO2Sxra226pb210M2_001, hTeO2Sxra226pb210M2Sum_001, 0, 10000);
	NormalizePDFs(hTeO2Sxth230onlyM1_001, hTeO2Sxth230onlyM2_001, hTeO2Sxth230onlyM2Sum_001, 0, 10000);
	NormalizePDFs(hTeO2Sxpb210M1_0001, hTeO2Sxpb210M2_0001, hTeO2Sxpb210M2Sum_0001, 0, 10000);


	TFile *file2 = new TFile("MCProduction_SurfaceCrystal_1keV.root", "RECREATE");


	// hTeO2Spb210M1_01->Write();
	// hTeO2Spo210M1_001->Write();
	// hTeO2Spo210M1_01->Write();
	// hTeO2Sth232M1_01->Write();
	// hTeO2Su238M1_01->Write();


	hTeO2Sxpb210M1_001->Write();
	hTeO2Sxpb210M1_01->Write();
	hTeO2Sxpb210M1_1->Write();
	hTeO2Sxpb210M1_10->Write();
	hTeO2Sxpo210M1_001->Write();
	hTeO2Sxpo210M1_01->Write();
	hTeO2Sxpo210M1_1->Write();
	// hTeO2Sxpo210M1_10->Write();	
	hTeO2Sxth232M1_001->Write();
	hTeO2Sxth232M1_01->Write();
	hTeO2Sxth232M1_1->Write();
	hTeO2Sxth232M1_10->Write();	
	hTeO2Sxu238M1_001->Write();
	hTeO2Sxu238M1_01->Write();
	hTeO2Sxu238M1_1->Write();
	hTeO2Sxu238M1_10->Write();

	hTeO2Sxth232onlyM1_001->Write();
	hTeO2Sxra228pb208M1_001->Write();
	hTeO2Sxu238th230M1_001->Write();
	hTeO2Sxth230onlyM1_001->Write();
	hTeO2Sxra226pb210M1_001->Write();
	hTeO2Sxpb210M1_0001->Write();

///////////////////////////////

	// hTeO2Spb210M2_01->Write();
	// hTeO2Spo210M2_001->Write();
	// hTeO2Spo210M2_01->Write();
	// hTeO2Sth232M2_01->Write();
	// hTeO2Su238M2_01->Write();


	hTeO2Sxpb210M2_001->Write();
	hTeO2Sxpb210M2_01->Write();
	hTeO2Sxpb210M2_1->Write();
	hTeO2Sxpb210M2_10->Write();
	hTeO2Sxpo210M2_001->Write();
	hTeO2Sxpo210M2_01->Write();
	hTeO2Sxpo210M2_1->Write();
	// hTeO2Sxpo210M2_10->Write();	
	hTeO2Sxth232M2_001->Write();
	hTeO2Sxth232M2_01->Write();
	hTeO2Sxth232M2_1->Write();
	hTeO2Sxth232M2_10->Write();	
	hTeO2Sxu238M2_001->Write();
	hTeO2Sxu238M2_01->Write();
	hTeO2Sxu238M2_1->Write();
	hTeO2Sxu238M2_10->Write();


	hTeO2Sxth232onlyM2_001->Write();
	hTeO2Sxra228pb208M2_001->Write();
	hTeO2Sxu238th230M2_001->Write();
	hTeO2Sxth230onlyM2_001->Write();
	hTeO2Sxra226pb210M2_001->Write();
	hTeO2Sxpb210M2_0001->Write();

///////////////////////////////

	// hTeO2Spb210M2Sum_01->Write();
	// hTeO2Spo210M2Sum_001->Write();
	// hTeO2Spo210M2Sum_01->Write();
	// hTeO2Sth232M2Sum_01->Write();
	// hTeO2Su238M2Sum_01->Write();


	// hTeO2Sxpb210M2Sum_001->Write();
	// hTeO2Sxpb210M2Sum_01->Write();
	// hTeO2Sxpb210M2Sum_1->Write();
	// hTeO2Sxpb210M2Sum_10->Write();
	// hTeO2Sxpo210M2Sum_001->Write();
	// hTeO2Sxpo210M2Sum_01->Write();
	// hTeO2Sxpo210M2Sum_1->Write();
	// // hTeO2Sxpo210M2Sum_10->Write();	
	// hTeO2Sxth232M2Sum_001->Write();
	// hTeO2Sxth232M2Sum_01->Write();
	// hTeO2Sxth232M2Sum_1->Write();
	// hTeO2Sxth232M2Sum_10->Write();	
	// hTeO2Sxu238M2Sum_001->Write();
	// hTeO2Sxu238M2Sum_01->Write();
	// hTeO2Sxu238M2Sum_1->Write();
	// hTeO2Sxu238M2Sum_10->Write();

	// hTeO2Sxth232onlyM2Sum_001->Write();
	// hTeO2Sxra228pb208M2Sum_001->Write();
	// hTeO2Sxu238th230M2Sum_001->Write();
	// hTeO2Sxth230onlyM2Sum_001->Write();
	// hTeO2Sxra226pb210M2Sum_001->Write();
	// hTeO2Sxpb210M2Sum_0001->Write();

	file2->Write();

}


void SaveHistogramsSurfaceOther()
{


	std::string sDataDir = "/cuore/data/simulation/CUORE0/t14.08/production_g2root-r363/ntp/Sup/";

	double dMinEnergy = 0;
	double dMaxEnergy = 10000;
	int dBinSize = 1;
	int dNBins = (dMaxEnergy - dMinEnergy)/dBinSize;

	// TH1D *hCuBoxSth232M1_1 = new TH1D("hCuBoxSth232M1_1", "", dNBins, dMinEnergy, dMaxEnergy);	
	// TH1D *hCuBoxSu238M1_1 = new TH1D("hCuBoxSu238M1_1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuBoxSxpb210M1_001 = new TH1D("hCuBoxSxpb210M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxpb210M1_01 = new TH1D("hCuBoxSxpb210M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxpb210M1_1 = new TH1D("hCuBoxSxpb210M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxpb210M1_10 = new TH1D("hCuBoxSxpb210M1_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M1_001 = new TH1D("hCuBoxSxth232M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M1_01 = new TH1D("hCuBoxSxth232M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M1_1 = new TH1D("hCuBoxSxth232M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M1_10 = new TH1D("hCuBoxSxth232M1_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M1_001 = new TH1D("hCuBoxSxu238M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M1_01 = new TH1D("hCuBoxSxu238M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M1_1 = new TH1D("hCuBoxSxu238M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M1_10 = new TH1D("hCuBoxSxu238M1_10", "", dNBins, dMinEnergy, dMaxEnergy);

	// TH1D *hCuFrameSth232M1_1 = new TH1D("hCuFrameSth232M1_1", "", dNBins, dMinEnergy, dMaxEnergy);	
	// TH1D *hCuFrameSu238M1_1 = new TH1D("hCuFrameSu238M1_1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuFrameSxpb210M1_001 = new TH1D("hCuFrameSxpb210M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxpb210M1_01 = new TH1D("hCuFrameSxpb210M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxpb210M1_1 = new TH1D("hCuFrameSxpb210M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxpb210M1_10 = new TH1D("hCuFrameSxpb210M1_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxth232M1_001 = new TH1D("hCuFrameSxth232M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxth232M1_01 = new TH1D("hCuFrameSxth232M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxth232M1_1 = new TH1D("hCuFrameSxth232M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxth232M1_10 = new TH1D("hCuFrameSxth232M1_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxu238M1_001 = new TH1D("hCuFrameSxu238M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxu238M1_01 = new TH1D("hCuFrameSxu238M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxu238M1_1 = new TH1D("hCuFrameSxu238M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxu238M1_10 = new TH1D("hCuFrameSxu238M1_10", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuBox_CuFramepb210M1_001 = new TH1D("hCuBox_CuFramepb210M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFramepb210M1_01 = new TH1D("hCuBox_CuFramepb210M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFramepb210M1_1 = new TH1D("hCuBox_CuFramepb210M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFramepb210M1_10 = new TH1D("hCuBox_CuFramepb210M1_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M1_001 = new TH1D("hCuBox_CuFrameth232M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M1_01 = new TH1D("hCuBox_CuFrameth232M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M1_1 = new TH1D("hCuBox_CuFrameth232M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M1_10 = new TH1D("hCuBox_CuFrameth232M1_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M1_001 = new TH1D("hCuBox_CuFrameu238M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M1_01 = new TH1D("hCuBox_CuFrameu238M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M1_1 = new TH1D("hCuBox_CuFrameu238M1_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M1_10 = new TH1D("hCuBox_CuFrameu238M1_10", "", dNBins, dMinEnergy, dMaxEnergy);

//////////////////////////////
	// TH1D *hCuBoxSth232M2_1 = new TH1D("hCuBoxSth232M2_1", "", dNBins, dMinEnergy, dMaxEnergy);	
	// TH1D *hCuBoxSu238M2_1 = new TH1D("hCuBoxSu238M2_1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuBoxSxpb210M2_001 = new TH1D("hCuBoxSxpb210M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxpb210M2_01 = new TH1D("hCuBoxSxpb210M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxpb210M2_1 = new TH1D("hCuBoxSxpb210M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxpb210M2_10 = new TH1D("hCuBoxSxpb210M2_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M2_001 = new TH1D("hCuBoxSxth232M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M2_01 = new TH1D("hCuBoxSxth232M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M2_1 = new TH1D("hCuBoxSxth232M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M2_10 = new TH1D("hCuBoxSxth232M2_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M2_001 = new TH1D("hCuBoxSxu238M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M2_01 = new TH1D("hCuBoxSxu238M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M2_1 = new TH1D("hCuBoxSxu238M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M2_10 = new TH1D("hCuBoxSxu238M2_10", "", dNBins, dMinEnergy, dMaxEnergy);

	// TH1D *hCuFrameSth232M2_1 = new TH1D("hCuFrameSth232M2_1", "", dNBins, dMinEnergy, dMaxEnergy);	
	// TH1D *hCuFrameSu238M2_1 = new TH1D("hCuFrameSu238M2_1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuFrameSxpb210M2_001 = new TH1D("hCuFrameSxpb210M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxpb210M2_01 = new TH1D("hCuFrameSxpb210M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxpb210M2_1 = new TH1D("hCuFrameSxpb210M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxpb210M2_10 = new TH1D("hCuFrameSxpb210M2_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxth232M2_001 = new TH1D("hCuFrameSxth232M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxth232M2_01 = new TH1D("hCuFrameSxth232M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxth232M2_1 = new TH1D("hCuFrameSxth232M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxth232M2_10 = new TH1D("hCuFrameSxth232M2_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxu238M2_001 = new TH1D("hCuFrameSxu238M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxu238M2_01 = new TH1D("hCuFrameSxu238M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxu238M2_1 = new TH1D("hCuFrameSxu238M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxu238M2_10 = new TH1D("hCuFrameSxu238M2_10", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuBox_CuFramepb210M2_001 = new TH1D("hCuBox_CuFramepb210M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFramepb210M2_01 = new TH1D("hCuBox_CuFramepb210M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFramepb210M2_1 = new TH1D("hCuBox_CuFramepb210M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFramepb210M2_10 = new TH1D("hCuBox_CuFramepb210M2_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M2_001 = new TH1D("hCuBox_CuFrameth232M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M2_01 = new TH1D("hCuBox_CuFrameth232M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M2_1 = new TH1D("hCuBox_CuFrameth232M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M2_10 = new TH1D("hCuBox_CuFrameth232M2_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M2_001 = new TH1D("hCuBox_CuFrameu238M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M2_01 = new TH1D("hCuBox_CuFrameu238M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M2_1 = new TH1D("hCuBox_CuFrameu238M2_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M2_10 = new TH1D("hCuBox_CuFrameu238M2_10", "", dNBins, dMinEnergy, dMaxEnergy);
//////////////////////////////
	// TH1D *hCuBoxSth232M2Sum_1 = new TH1D("hCuBoxSth232M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);	
	// TH1D *hCuBoxSu238M2Sum_1 = new TH1D("hCuBoxSu238M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuBoxSxpb210M2Sum_001 = new TH1D("hCuBoxSxpb210M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxpb210M2Sum_01 = new TH1D("hCuBoxSxpb210M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxpb210M2Sum_1 = new TH1D("hCuBoxSxpb210M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxpb210M2Sum_10 = new TH1D("hCuBoxSxpb210M2Sum_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M2Sum_001 = new TH1D("hCuBoxSxth232M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M2Sum_01 = new TH1D("hCuBoxSxth232M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M2Sum_1 = new TH1D("hCuBoxSxth232M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxth232M2Sum_10 = new TH1D("hCuBoxSxth232M2Sum_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M2Sum_001 = new TH1D("hCuBoxSxu238M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M2Sum_01 = new TH1D("hCuBoxSxu238M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M2Sum_1 = new TH1D("hCuBoxSxu238M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxSxu238M2Sum_10 = new TH1D("hCuBoxSxu238M2Sum_10", "", dNBins, dMinEnergy, dMaxEnergy);

	// TH1D *hCuFrameSth232M2Sum_1 = new TH1D("hCuFrameSth232M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);	
	// TH1D *hCuFrameSu238M2Sum_1 = new TH1D("hCuFrameSu238M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuFrameSxpb210M2Sum_001 = new TH1D("hCuFrameSxpb210M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxpb210M2Sum_01 = new TH1D("hCuFrameSxpb210M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxpb210M2Sum_1 = new TH1D("hCuFrameSxpb210M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxpb210M2Sum_10 = new TH1D("hCuFrameSxpb210M2Sum_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxth232M2Sum_001 = new TH1D("hCuFrameSxth232M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxth232M2Sum_01 = new TH1D("hCuFrameSxth232M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxth232M2Sum_1 = new TH1D("hCuFrameSxth232M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxth232M2Sum_10 = new TH1D("hCuFrameSxth232M2Sum_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxu238M2Sum_001 = new TH1D("hCuFrameSxu238M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxu238M2Sum_01 = new TH1D("hCuFrameSxu238M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxu238M2Sum_1 = new TH1D("hCuFrameSxu238M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameSxu238M2Sum_10 = new TH1D("hCuFrameSxu238M2Sum_10", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuBox_CuFramepb210M2Sum_001 = new TH1D("hCuBox_CuFramepb210M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFramepb210M2Sum_01 = new TH1D("hCuBox_CuFramepb210M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFramepb210M2Sum_1 = new TH1D("hCuBox_CuFramepb210M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFramepb210M2Sum_10 = new TH1D("hCuBox_CuFramepb210M2Sum_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M2Sum_001 = new TH1D("hCuBox_CuFrameth232M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M2Sum_01 = new TH1D("hCuBox_CuFrameth232M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M2Sum_1 = new TH1D("hCuBox_CuFrameth232M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameth232M2Sum_10 = new TH1D("hCuBox_CuFrameth232M2Sum_10", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M2Sum_001 = new TH1D("hCuBox_CuFrameu238M2Sum_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M2Sum_01 = new TH1D("hCuBox_CuFrameu238M2Sum_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M2Sum_1 = new TH1D("hCuBox_CuFrameu238M2Sum_1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBox_CuFrameu238M2Sum_10 = new TH1D("hCuBox_CuFrameu238M2Sum_10", "", dNBins, dMinEnergy, dMaxEnergy);

	// TChain *tCuBoxSth232_1 = LoadMC(sDataDir.c_str(), "CuBoxS", "th232-1");
	// TChain *tCuBoxSu238_1 = LoadMC(sDataDir.c_str(), "CuBoxS", "u238-1");

	TChain *tCuBoxSxpb210_001 = LoadMC(sDataDir.c_str(), "CuBoxSx", "pb210-.01");
	TChain *tCuBoxSxpb210_01 = LoadMC(sDataDir.c_str(), "CuBoxSx", "pb210-.1");
	TChain *tCuBoxSxpb210_1 = LoadMC(sDataDir.c_str(), "CuBoxSx", "pb210-1");
	TChain *tCuBoxSxpb210_10 = LoadMC(sDataDir.c_str(), "CuBoxSx", "pb210-10");
	TChain *tCuBoxSxth232_001 = LoadMC(sDataDir.c_str(), "CuBoxSx", "th232-.01");
	TChain *tCuBoxSxth232_01 = LoadMC(sDataDir.c_str(), "CuBoxSx", "th232-.1");
	TChain *tCuBoxSxth232_1 = LoadMC(sDataDir.c_str(), "CuBoxSx", "th232-1");
	TChain *tCuBoxSxth232_10 = LoadMC(sDataDir.c_str(), "CuBoxSx", "th232-10");
	TChain *tCuBoxSxu238_001 = LoadMC(sDataDir.c_str(), "CuBoxSx", "u238-.01");
	TChain *tCuBoxSxu238_01 = LoadMC(sDataDir.c_str(), "CuBoxSx", "u238-.1");
	TChain *tCuBoxSxu238_1 = LoadMC(sDataDir.c_str(), "CuBoxSx", "u238-1");
	TChain *tCuBoxSxu238_10 = LoadMC(sDataDir.c_str(), "CuBoxSx", "u238-10");

	// TChain *tCuFrameSth232_1 = LoadMC(sDataDir.c_str(), "CuFrameS", "th232-1");
	// TChain *tCuFrameSu238_1 = LoadMC(sDataDir.c_str(), "CuFrameS", "u238-1");

	TChain *tCuFrameSxpb210_001 = LoadMC(sDataDir.c_str(), "CuFrameSx", "pb210-.01");
	TChain *tCuFrameSxpb210_01 = LoadMC(sDataDir.c_str(), "CuFrameSx", "pb210-.1");
	TChain *tCuFrameSxpb210_1 = LoadMC(sDataDir.c_str(), "CuFrameSx", "pb210-1");
	TChain *tCuFrameSxpb210_10 = LoadMC(sDataDir.c_str(), "CuFrameSx", "pb210-10");
	TChain *tCuFrameSxth232_001 = LoadMC(sDataDir.c_str(), "CuFrameSx", "th232-.01");
	TChain *tCuFrameSxth232_01 = LoadMC(sDataDir.c_str(), "CuFrameSx", "th232-.1");
	TChain *tCuFrameSxth232_1 = LoadMC(sDataDir.c_str(), "CuFrameSx", "th232-1");
	TChain *tCuFrameSxth232_10 = LoadMC(sDataDir.c_str(), "CuFrameSx", "th232-10");
	TChain *tCuFrameSxu238_001 = LoadMC(sDataDir.c_str(), "CuFrameSx", "u238-.01");
	TChain *tCuFrameSxu238_01 = LoadMC(sDataDir.c_str(), "CuFrameSx", "u238-.1");
	TChain *tCuFrameSxu238_1 = LoadMC(sDataDir.c_str(), "CuFrameSx", "u238-1");
	TChain *tCuFrameSxu238_10 = LoadMC(sDataDir.c_str(), "CuFrameSx", "u238-10");

	// TChain *tCuBox_CuFramepb210_001 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "pb210-.01");
	// TChain *tCuBox_CuFramepb210_01 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "pb210-.1");
	// TChain *tCuBox_CuFramepb210_1 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "pb210-1");
	// TChain *tCuBox_CuFramepb210_10 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "pb210-10");
	// TChain *tCuBox_CuFrameth232_001 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "th232-.01");
	// TChain *tCuBox_CuFrameth232_01 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "th232-.1");
	// TChain *tCuBox_CuFrameth232_1 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "th232-1");
	// TChain *tCuBox_CuFrameth232_10 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "th232-10");
	// TChain *tCuBox_CuFrameu238_001 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "u238-.01");
	// TChain *tCuBox_CuFrameu238_01 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "u238-.1");
	// TChain *tCuBox_CuFrameu238_1 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "u238-1");
	// TChain *tCuBox_CuFrameu238_10 = LoadMC(sDataDir.c_str(), "CuBox+CuFrame", "u238-10");

///////////////////
	// tCuBoxSth232_1->Project("hCuBoxSth232M1_1", "Ener2", "Multiplicity==1");
	// tCuBoxSu238_1->Project("hCuBoxSu238M1_1", "Ener2", "Multiplicity==1");

	tCuBoxSxpb210_001->Project("hCuBoxSxpb210M1_001", "Ener2", "Multiplicity==1");
	tCuBoxSxpb210_01->Project("hCuBoxSxpb210M1_01", "Ener2", "Multiplicity==1");
	tCuBoxSxpb210_1->Project("hCuBoxSxpb210M1_1", "Ener2", "Multiplicity==1");
	tCuBoxSxpb210_10->Project("hCuBoxSxpb210M1_10", "Ener2", "Multiplicity==1");
	tCuBoxSxth232_001->Project("hCuBoxSxth232M1_001", "Ener2", "Multiplicity==1");
	tCuBoxSxth232_01->Project("hCuBoxSxth232M1_01", "Ener2", "Multiplicity==1");
	tCuBoxSxth232_1->Project("hCuBoxSxth232M1_1", "Ener2", "Multiplicity==1");
	tCuBoxSxth232_10->Project("hCuBoxSxth232M1_10", "Ener2", "Multiplicity==1");
	tCuBoxSxu238_001->Project("hCuBoxSxu238M1_001", "Ener2", "Multiplicity==1");
	tCuBoxSxu238_01->Project("hCuBoxSxu238M1_01", "Ener2", "Multiplicity==1");
	tCuBoxSxu238_1->Project("hCuBoxSxu238M1_1", "Ener2", "Multiplicity==1");
	tCuBoxSxu238_10->Project("hCuBoxSxu238M1_10", "Ener2", "Multiplicity==1");

	// tCuFrameSth232_1->Project("hCuFrameSth232M1_1", "Ener2", "Multiplicity==1");
	// tCuFrameSu238_1->Project("hCuFrameSu238M1_1", "Ener2", "Multiplicity==1");

	tCuFrameSxpb210_001->Project("hCuFrameSxpb210M1_001", "Ener2", "Multiplicity==1");
	tCuFrameSxpb210_01->Project("hCuFrameSxpb210M1_01", "Ener2", "Multiplicity==1");
	tCuFrameSxpb210_1->Project("hCuFrameSxpb210M1_1", "Ener2", "Multiplicity==1");
	tCuFrameSxpb210_10->Project("hCuFrameSxpb210M1_10", "Ener2", "Multiplicity==1");
	tCuFrameSxth232_001->Project("hCuFrameSxth232M1_001", "Ener2", "Multiplicity==1");
	tCuFrameSxth232_01->Project("hCuFrameSxth232M1_01", "Ener2", "Multiplicity==1");
	tCuFrameSxth232_1->Project("hCuFrameSxth232M1_1", "Ener2", "Multiplicity==1");
	tCuFrameSxth232_10->Project("hCuFrameSxth232M1_10", "Ener2", "Multiplicity==1");
	tCuFrameSxu238_001->Project("hCuFrameSxu238M1_001", "Ener2", "Multiplicity==1");
	tCuFrameSxu238_01->Project("hCuFrameSxu238M1_01", "Ener2", "Multiplicity==1");
	tCuFrameSxu238_1->Project("hCuFrameSxu238M1_1", "Ener2", "Multiplicity==1");
	tCuFrameSxu238_10->Project("hCuFrameSxu238M1_10", "Ener2", "Multiplicity==1");

	// tCuBox_CuFramepb210_001->Project("hCuBox_CuFramepb210M1_001", "Ener2", "Multiplicity==1");
	// tCuBox_CuFramepb210_01->Project("hCuBox_CuFramepb210M1_01", "Ener2", "Multiplicity==1");
	// tCuBox_CuFramepb210_1->Project("hCuBox_CuFramepb210M1_1", "Ener2", "Multiplicity==1");
	// tCuBox_CuFramepb210_10->Project("hCuBox_CuFramepb210M1_10", "Ener2", "Multiplicity==1");
	// tCuBox_CuFrameth232_001->Project("hCuBox_CuFrameth232M1_001", "Ener2", "Multiplicity==1");
	// tCuBox_CuFrameth232_01->Project("hCuBox_CuFrameth232M1_01", "Ener2", "Multiplicity==1");
	// tCuBox_CuFrameth232_1->Project("hCuBox_CuFrameth232M1_1", "Ener2", "Multiplicity==1");
	// tCuBox_CuFrameth232_10->Project("hCuBox_CuFrameth232M1_10", "Ener2", "Multiplicity==1");
	// tCuBox_CuFrameu238_001->Project("hCuBox_CuFrameu238M1_001", "Ener2", "Multiplicity==1");
	// tCuBox_CuFrameu238_01->Project("hCuBox_CuFrameu238M1_01", "Ener2", "Multiplicity==1");
	// tCuBox_CuFrameu238_1->Project("hCuBox_CuFrameu238M1_1", "Ener2", "Multiplicity==1");
	// tCuBox_CuFrameu238_10->Project("hCuBox_CuFrameu238M1_10", "Ener2", "Multiplicity==1");


//////////////////////////////////
	// tCuBoxSth232_1->Project("hCuBoxSth232M2_1", "Ener2", "Multiplicity==2");
	// tCuBoxSu238_1->Project("hCuBoxSu238M2_1", "Ener2", "Multiplicity==2");

	tCuBoxSxpb210_001->Project("hCuBoxSxpb210M2_001", "Ener2", "Multiplicity==2");
	tCuBoxSxpb210_01->Project("hCuBoxSxpb210M2_01", "Ener2", "Multiplicity==2");
	tCuBoxSxpb210_1->Project("hCuBoxSxpb210M2_1", "Ener2", "Multiplicity==2");
	tCuBoxSxpb210_10->Project("hCuBoxSxpb210M2_10", "Ener2", "Multiplicity==2");
	tCuBoxSxth232_001->Project("hCuBoxSxth232M2_001", "Ener2", "Multiplicity==2");
	tCuBoxSxth232_01->Project("hCuBoxSxth232M2_01", "Ener2", "Multiplicity==2");
	tCuBoxSxth232_1->Project("hCuBoxSxth232M2_1", "Ener2", "Multiplicity==2");
	tCuBoxSxth232_10->Project("hCuBoxSxth232M2_10", "Ener2", "Multiplicity==2");
	tCuBoxSxu238_001->Project("hCuBoxSxu238M2_001", "Ener2", "Multiplicity==2");
	tCuBoxSxu238_01->Project("hCuBoxSxu238M2_01", "Ener2", "Multiplicity==2");
	tCuBoxSxu238_1->Project("hCuBoxSxu238M2_1", "Ener2", "Multiplicity==2");
	tCuBoxSxu238_10->Project("hCuBoxSxu238M2_10", "Ener2", "Multiplicity==2");


	// tCuFrameSth232_1->Project("hCuFrameSth232M2_1", "Ener2", "Multiplicity==2");
	// tCuFrameSu238_1->Project("hCuFrameSu238M2_1", "Ener2", "Multiplicity==2");

	tCuFrameSxpb210_001->Project("hCuFrameSxpb210M2_001", "Ener2", "Multiplicity==2");
	tCuFrameSxpb210_01->Project("hCuFrameSxpb210M2_01", "Ener2", "Multiplicity==2");
	tCuFrameSxpb210_1->Project("hCuFrameSxpb210M2_1", "Ener2", "Multiplicity==2");
	tCuFrameSxpb210_10->Project("hCuFrameSxpb210M2_10", "Ener2", "Multiplicity==2");
	tCuFrameSxth232_001->Project("hCuFrameSxth232M2_001", "Ener2", "Multiplicity==2");
	tCuFrameSxth232_01->Project("hCuFrameSxth232M2_01", "Ener2", "Multiplicity==2");
	tCuFrameSxth232_1->Project("hCuFrameSxth232M2_1", "Ener2", "Multiplicity==2");
	tCuFrameSxth232_10->Project("hCuFrameSxth232M2_10", "Ener2", "Multiplicity==2");
	tCuFrameSxu238_001->Project("hCuFrameSxu238M2_001", "Ener2", "Multiplicity==2");
	tCuFrameSxu238_01->Project("hCuFrameSxu238M2_01", "Ener2", "Multiplicity==2");
	tCuFrameSxu238_1->Project("hCuFrameSxu238M2_1", "Ener2", "Multiplicity==2");
	tCuFrameSxu238_10->Project("hCuFrameSxu238M2_10", "Ener2", "Multiplicity==2");

	// tCuBox_CuFramepb210_001->Project("hCuBox_CuFramepb210M2_001", "Ener2", "Multiplicity==2");
	// tCuBox_CuFramepb210_01->Project("hCuBox_CuFramepb210M2_01", "Ener2", "Multiplicity==2");
	// tCuBox_CuFramepb210_1->Project("hCuBox_CuFramepb210M2_1", "Ener2", "Multiplicity==2");
	// tCuBox_CuFramepb210_10->Project("hCuBox_CuFramepb210M2_10", "Ener2", "Multiplicity==2");
	// tCuBox_CuFrameth232_001->Project("hCuBox_CuFrameth232M2_001", "Ener2", "Multiplicity==2");
	// tCuBox_CuFrameth232_01->Project("hCuBox_CuFrameth232M2_01", "Ener2", "Multiplicity==2");
	// tCuBox_CuFrameth232_1->Project("hCuBox_CuFrameth232M2_1", "Ener2", "Multiplicity==2");
	// tCuBox_CuFrameth232_10->Project("hCuBox_CuFrameth232M2_10", "Ener2", "Multiplicity==2");
	// tCuBox_CuFrameu238_001->Project("hCuBox_CuFrameu238M2_001", "Ener2", "Multiplicity==2");
	// tCuBox_CuFrameu238_01->Project("hCuBox_CuFrameu238M2_01", "Ener2", "Multiplicity==2");
	// tCuBox_CuFrameu238_1->Project("hCuBox_CuFrameu238M2_1", "Ener2", "Multiplicity==2");
	// tCuBox_CuFrameu238_10->Project("hCuBox_CuFrameu238M2_10", "Ener2", "Multiplicity==2");


//////////////////////////////////
	// tCuBoxSth232_1->Project("hCuBoxSth232M2Sum_1", "Esum2", "Multiplicity==2");
	// tCuBoxSu238_1->Project("hCuBoxSu238M2Sum_1", "Esum2", "Multiplicity==2");

	// tCuBoxSxpb210_001->Project("hCuBoxSxpb210M2Sum_001", "Esum2", "Multiplicity==2");
	// tCuBoxSxpb210_01->Project("hCuBoxSxpb210M2Sum_01", "Esum2", "Multiplicity==2");
	// tCuBoxSxpb210_1->Project("hCuBoxSxpb210M2Sum_1", "Esum2", "Multiplicity==2");
	// tCuBoxSxpb210_10->Project("hCuBoxSxpb210M2Sum_10", "Esum2", "Multiplicity==2");
	// tCuBoxSxth232_001->Project("hCuBoxSxth232M2Sum_001", "Esum2", "Multiplicity==2");
	// tCuBoxSxth232_01->Project("hCuBoxSxth232M2Sum_01", "Esum2", "Multiplicity==2");
	// tCuBoxSxth232_1->Project("hCuBoxSxth232M2Sum_1", "Esum2", "Multiplicity==2");
	// tCuBoxSxth232_10->Project("hCuBoxSxth232M2Sum_10", "Esum2", "Multiplicity==2");
	// tCuBoxSxu238_001->Project("hCuBoxSxu238M2Sum_001", "Esum2", "Multiplicity==2");
	// tCuBoxSxu238_01->Project("hCuBoxSxu238M2Sum_01", "Esum2", "Multiplicity==2");
	// tCuBoxSxu238_1->Project("hCuBoxSxu238M2Sum_1", "Esum2", "Multiplicity==2");
	// tCuBoxSxu238_10->Project("hCuBoxSxu238M2Sum_10", "Esum2", "Multiplicity==2");


	// tCuFrameSth232_1->Project("hCuFrameSth232M2Sum_1", "Esum2", "Multiplicity==2");
	// tCuFrameSu238_1->Project("hCuFrameSu238M2Sum_1", "Esum2", "Multiplicity==2");

	// tCuFrameSxpb210_001->Project("hCuFrameSxpb210M2Sum_001", "Esum2", "Multiplicity==2");
	// tCuFrameSxpb210_01->Project("hCuFrameSxpb210M2Sum_01", "Esum2", "Multiplicity==2");
	// tCuFrameSxpb210_1->Project("hCuFrameSxpb210M2Sum_1", "Esum2", "Multiplicity==2");
	// tCuFrameSxpb210_10->Project("hCuFrameSxpb210M2Sum_10", "Esum2", "Multiplicity==2");
	// tCuFrameSxth232_001->Project("hCuFrameSxth232M2Sum_001", "Esum2", "Multiplicity==2");
	// tCuFrameSxth232_01->Project("hCuFrameSxth232M2Sum_01", "Esum2", "Multiplicity==2");
	// tCuFrameSxth232_1->Project("hCuFrameSxth232M2Sum_1", "Esum2", "Multiplicity==2");
	// tCuFrameSxth232_10->Project("hCuFrameSxth232M2Sum_10", "Esum2", "Multiplicity==2");
	// tCuFrameSxu238_001->Project("hCuFrameSxu238M2Sum_001", "Esum2", "Multiplicity==2");
	// tCuFrameSxu238_01->Project("hCuFrameSxu238M2Sum_01", "Esum2", "Multiplicity==2");
	// tCuFrameSxu238_1->Project("hCuFrameSxu238M2Sum_1", "Esum2", "Multiplicity==2");
	// tCuFrameSxu238_10->Project("hCuFrameSxu238M2Sum_10", "Esum2", "Multiplicity==2");

	// tCuBox_CuFramepb210_001->Project("hCuBox_CuFramepb210M2Sum_001", "Esum2", "Multiplicity==2");
	// tCuBox_CuFramepb210_01->Project("hCuBox_CuFramepb210M2Sum_01", "Esum2", "Multiplicity==2");
	// tCuBox_CuFramepb210_1->Project("hCuBox_CuFramepb210M2Sum_1", "Esum2", "Multiplicity==2");
	// tCuBox_CuFramepb210_10->Project("hCuBox_CuFramepb210M2Sum_10", "Esum2", "Multiplicity==2");
	// tCuBox_CuFrameth232_001->Project("hCuBox_CuFrameth232M2Sum_001", "Esum2", "Multiplicity==2");
	// tCuBox_CuFrameth232_01->Project("hCuBox_CuFrameth232M2Sum_01", "Esum2", "Multiplicity==2");
	// tCuBox_CuFrameth232_1->Project("hCuBox_CuFrameth232M2Sum_1", "Esum2", "Multiplicity==2");
	// tCuBox_CuFrameth232_10->Project("hCuBox_CuFrameth232M2Sum_10", "Esum2", "Multiplicity==2");
	// tCuBox_CuFrameu238_001->Project("hCuBox_CuFrameu238M2Sum_001", "Esum2", "Multiplicity==2");
	// tCuBox_CuFrameu238_01->Project("hCuBox_CuFrameu238M2Sum_01", "Esum2", "Multiplicity==2");
	// tCuBox_CuFrameu238_1->Project("hCuBox_CuFrameu238M2Sum_1", "Esum2", "Multiplicity==2");
	// tCuBox_CuFrameu238_10->Project("hCuBox_CuFrameu238M2Sum_10", "Esum2", "Multiplicity==2");

	//////////////////////////////////
	// CuBox+CuFrame here
	// Total = CuBox + 0.2444254783*CuFrame
	hCuBox_CuFramepb210M1_10->Add(hCuBoxSxpb210M1_10, hCuFrameSxpb210M1_10, 1.0, 0.2444254783);
	hCuBox_CuFramepb210M1_1->Add(hCuBoxSxpb210M1_1, hCuFrameSxpb210M1_1, 1.0, 0.2444254783);
	hCuBox_CuFramepb210M1_01->Add(hCuBoxSxpb210M1_01, hCuFrameSxpb210M1_01, 1.0, 0.2444254783);
	hCuBox_CuFramepb210M1_001->Add(hCuBoxSxpb210M1_001, hCuFrameSxpb210M1_001, 1.0, 0.2444254783);

	hCuBox_CuFrameth232M1_10->Add(hCuBoxSxth232M1_10, hCuFrameSxth232M1_10, 1.0, 0.2444254783);
	hCuBox_CuFrameu238M1_10->Add(hCuBoxSxu238M1_10, hCuFrameSxu238M1_10, 1.0, 0.2444254783);

	hCuBox_CuFramepb210M2_10->Add(hCuBoxSxpb210M2_10, hCuFrameSxpb210M2_10, 1.0, 0.2444254783);
	hCuBox_CuFramepb210M2_1->Add(hCuBoxSxpb210M2_1, hCuFrameSxpb210M2_1, 1.0, 0.2444254783);
	hCuBox_CuFramepb210M2_01->Add(hCuBoxSxpb210M2_01, hCuFrameSxpb210M2_01, 1.0, 0.2444254783);
	hCuBox_CuFramepb210M2_001->Add(hCuBoxSxpb210M2_001, hCuFrameSxpb210M2_001, 1.0, 0.2444254783);

	hCuBox_CuFrameth232M2_10->Add(hCuBoxSxth232M2_10, hCuFrameSxth232M2_10, 1.0, 0.2444254783);
	hCuBox_CuFrameu238M2_10->Add(hCuBoxSxu238M2_10, hCuFrameSxu238M2_10, 1.0, 0.2444254783);

////////////////////////////



	// NormalizePDFs(hCuBoxSth232M1_1, hCuBoxSth232M2_1, hCuBoxSth232M2Sum_1, 0, 10000);
	// NormalizePDFs(hCuBoxSu238M1_1, hCuBoxSu238M2_1, hCuBoxSu238M2Sum_1, 0, 10000);

	NormalizePDFs(hCuBoxSxpb210M1_001, hCuBoxSxpb210M2_001, hCuBoxSxpb210M2Sum_001, 0, 10000);
	NormalizePDFs(hCuBoxSxpb210M1_01, hCuBoxSxpb210M2_01, hCuBoxSxpb210M2Sum_01, 0, 10000);
	NormalizePDFs(hCuBoxSxpb210M1_1, hCuBoxSxpb210M2_1, hCuBoxSxpb210M2Sum_1, 0, 10000);
	NormalizePDFs(hCuBoxSxpb210M1_10, hCuBoxSxpb210M2_10, hCuBoxSxpb210M2Sum_10, 0, 10000);
	NormalizePDFs(hCuBoxSxth232M1_001, hCuBoxSxth232M2_001, hCuBoxSxth232M2Sum_001, 0, 10000);
	NormalizePDFs(hCuBoxSxth232M1_01, hCuBoxSxth232M2_01, hCuBoxSxth232M2Sum_01, 0, 10000);
	NormalizePDFs(hCuBoxSxth232M1_1, hCuBoxSxth232M2_1, hCuBoxSxth232M2Sum_1, 0, 10000);
	NormalizePDFs(hCuBoxSxth232M1_10, hCuBoxSxth232M2_10, hCuBoxSxth232M2Sum_10, 0, 10000);
	NormalizePDFs(hCuBoxSxu238M1_001, hCuBoxSxu238M2_001, hCuBoxSxu238M2Sum_001, 0, 10000);
	NormalizePDFs(hCuBoxSxu238M1_01, hCuBoxSxu238M2_01, hCuBoxSxu238M2Sum_01, 0, 10000);
	NormalizePDFs(hCuBoxSxu238M1_1, hCuBoxSxu238M2_1, hCuBoxSxu238M2Sum_1, 0, 10000);
	NormalizePDFs(hCuBoxSxu238M1_10, hCuBoxSxu238M2_10, hCuBoxSxu238M2Sum_10, 0, 10000);

	// NormalizePDFs(hCuFrameSth232M1_1, hCuFrameSth232M2_1, hCuFrameSth232M2Sum_1, 0, 10000);
	// NormalizePDFs(hCuFrameSu238M1_1, hCuFrameSu238M2_1, hCuFrameSu238M2Sum_1, 0, 10000);

	NormalizePDFs(hCuFrameSxpb210M1_001, hCuFrameSxpb210M2_001, hCuFrameSxpb210M2Sum_001, 0, 10000);
	NormalizePDFs(hCuFrameSxpb210M1_01, hCuFrameSxpb210M2_01, hCuFrameSxpb210M2Sum_01, 0, 10000);
	NormalizePDFs(hCuFrameSxpb210M1_1, hCuFrameSxpb210M2_1, hCuFrameSxpb210M2Sum_1, 0, 10000);
	NormalizePDFs(hCuFrameSxpb210M1_10, hCuFrameSxpb210M2_10, hCuFrameSxpb210M2Sum_10, 0, 10000);
	NormalizePDFs(hCuFrameSxth232M1_001, hCuFrameSxth232M2_001, hCuFrameSxth232M2Sum_001, 0, 10000);
	NormalizePDFs(hCuFrameSxth232M1_01, hCuFrameSxth232M2_01, hCuFrameSxth232M2Sum_01, 0, 10000);
	NormalizePDFs(hCuFrameSxth232M1_1, hCuFrameSxth232M2_1, hCuFrameSxth232M2Sum_1, 0, 10000);
	NormalizePDFs(hCuFrameSxth232M1_10, hCuFrameSxth232M2_10, hCuFrameSxth232M2Sum_10, 0, 10000);
	NormalizePDFs(hCuFrameSxu238M1_001, hCuFrameSxu238M2_001, hCuFrameSxu238M2Sum_001, 0, 10000);
	NormalizePDFs(hCuFrameSxu238M1_01, hCuFrameSxu238M2_01, hCuFrameSxu238M2Sum_01, 0, 10000);
	NormalizePDFs(hCuFrameSxu238M1_1, hCuFrameSxu238M2_1, hCuFrameSxu238M2Sum_1, 0, 10000);
	NormalizePDFs(hCuFrameSxu238M1_10, hCuFrameSxu238M2_10, hCuFrameSxu238M2Sum_10, 0, 10000);	

	NormalizePDFs(hCuBox_CuFramepb210M1_001 , hCuBox_CuFramepb210M2_001 , hCuBox_CuFramepb210M2Sum_001, 0, 10000); 
	NormalizePDFs(hCuBox_CuFramepb210M1_01 , hCuBox_CuFramepb210M2_01 , hCuBox_CuFramepb210M2Sum_01, 0, 10000); 
	NormalizePDFs(hCuBox_CuFramepb210M1_1 , hCuBox_CuFramepb210M2_1 , hCuBox_CuFramepb210M2Sum_1, 0, 10000); 
	NormalizePDFs(hCuBox_CuFramepb210M1_10 , hCuBox_CuFramepb210M2_10 , hCuBox_CuFramepb210M2Sum_10, 0, 10000); 
	// NormalizePDFs(hCuBox_CuFrameth232M1_001 , hCuBox_CuFrameth232M2_001 , hCuBox_CuFrameth232M2Sum_001, 0, 10000); 
	// NormalizePDFs(hCuBox_CuFrameth232M1_01 , hCuBox_CuFrameth232M2_01 , hCuBox_CuFrameth232M2Sum_01, 0, 10000); 
	// NormalizePDFs(hCuBox_CuFrameth232M1_1 , hCuBox_CuFrameth232M2_1 , hCuBox_CuFrameth232M2Sum_1, 0, 10000); 
	NormalizePDFs(hCuBox_CuFrameth232M1_10 , hCuBox_CuFrameth232M2_10 , hCuBox_CuFrameth232M2Sum_10, 0, 10000); 
	// NormalizePDFs(hCuBox_CuFrameu238M1_001 , hCuBox_CuFrameu238M2_001 , hCuBox_CuFrameu238M2Sum_001, 0, 10000); 
	// NormalizePDFs(hCuBox_CuFrameu238M1_01 , hCuBox_CuFrameu238M2_01 , hCuBox_CuFrameu238M2Sum_01, 0, 10000); 
	// NormalizePDFs(hCuBox_CuFrameu238M1_1 , hCuBox_CuFrameu238M2_1 , hCuBox_CuFrameu238M2Sum_1, 0, 10000); 
	NormalizePDFs(hCuBox_CuFrameu238M1_10 , hCuBox_CuFrameu238M2_10 , hCuBox_CuFrameu238M2Sum_10, 0, 10000); 

	TFile *file2 = new TFile("MCProduction_SurfaceOther_1keV.root", "RECREATE");

	// hCuBoxSth232M1_1->Write();
	// hCuBoxSu238M1_1->Write();

	hCuBoxSxpb210M1_001->Write();
	hCuBoxSxpb210M1_01->Write();
	hCuBoxSxpb210M1_1->Write();
	hCuBoxSxpb210M1_10->Write();
	hCuBoxSxth232M1_001->Write();
	hCuBoxSxth232M1_01->Write();
	hCuBoxSxth232M1_1->Write();
	hCuBoxSxth232M1_10->Write();	
	hCuBoxSxu238M1_001->Write();
	hCuBoxSxu238M1_01->Write();
	hCuBoxSxu238M1_1->Write();
	hCuBoxSxu238M1_10->Write();

	// hCuFrameSth232M1_1->Write();
	// hCuFrameSu238M1_1->Write();

	hCuFrameSxpb210M1_001->Write();
	hCuFrameSxpb210M1_01->Write();
	hCuFrameSxpb210M1_1->Write();
	hCuFrameSxpb210M1_10->Write();
	hCuFrameSxth232M1_001->Write();
	hCuFrameSxth232M1_01->Write();
	hCuFrameSxth232M1_1->Write();
	hCuFrameSxth232M1_10->Write();	
	hCuFrameSxu238M1_001->Write();
	hCuFrameSxu238M1_01->Write();
	hCuFrameSxu238M1_1->Write();
	hCuFrameSxu238M1_10->Write();


	hCuBox_CuFramepb210M1_001->Write();
	hCuBox_CuFramepb210M1_01->Write();
	hCuBox_CuFramepb210M1_1->Write();
	hCuBox_CuFramepb210M1_10->Write();
	// hCuBox_CuFrameth232M1_001->Write();
	// hCuBox_CuFrameth232M1_01->Write();
	// hCuBox_CuFrameth232M1_1->Write();
	hCuBox_CuFrameth232M1_10->Write();
	// hCuBox_CuFrameu238M1_001->Write();
	// hCuBox_CuFrameu238M1_01->Write();
	// hCuBox_CuFrameu238M1_1->Write();
	hCuBox_CuFrameu238M1_10->Write();

///////////////////////////////
	// hCuBoxSth232M2_1->Write();
	// hCuBoxSu238M2_1->Write();

	hCuBoxSxpb210M2_001->Write();
	hCuBoxSxpb210M2_01->Write();
	hCuBoxSxpb210M2_1->Write();
	hCuBoxSxpb210M2_10->Write();
	hCuBoxSxth232M2_001->Write();
	hCuBoxSxth232M2_01->Write();
	hCuBoxSxth232M2_1->Write();
	hCuBoxSxth232M2_10->Write();	
	hCuBoxSxu238M2_001->Write();
	hCuBoxSxu238M2_01->Write();
	hCuBoxSxu238M2_1->Write();
	hCuBoxSxu238M2_10->Write();

	// hCuFrameSth232M2_1->Write();
	// hCuFrameSu238M2_1->Write();

	hCuFrameSxpb210M2_001->Write();
	hCuFrameSxpb210M2_01->Write();
	hCuFrameSxpb210M2_1->Write();
	hCuFrameSxpb210M2_10->Write();
	hCuFrameSxth232M2_001->Write();
	hCuFrameSxth232M2_01->Write();
	hCuFrameSxth232M2_1->Write();
	hCuFrameSxth232M2_10->Write();	
	hCuFrameSxu238M2_001->Write();
	hCuFrameSxu238M2_01->Write();
	hCuFrameSxu238M2_1->Write();
	hCuFrameSxu238M2_10->Write();

	hCuBox_CuFramepb210M2_001->Write();
	hCuBox_CuFramepb210M2_01->Write();
	hCuBox_CuFramepb210M2_1->Write();
	hCuBox_CuFramepb210M2_10->Write();
	// hCuBox_CuFrameth232M2_001->Write();
	// hCuBox_CuFrameth232M2_01->Write();
	// hCuBox_CuFrameth232M2_1->Write();
	hCuBox_CuFrameth232M2_10->Write();
	// hCuBox_CuFrameu238M2_001->Write();
	// hCuBox_CuFrameu238M2_01->Write();
	// hCuBox_CuFrameu238M2_1->Write();
	hCuBox_CuFrameu238M2_10->Write();

///////////////////////////////
	// hCuBoxSth232M2Sum_1->Write();
	// hCuBoxSu238M2Sum_1->Write();

	// hCuBoxSxpb210M2Sum_001->Write();
	// hCuBoxSxpb210M2Sum_01->Write();
	// hCuBoxSxpb210M2Sum_1->Write();
	// hCuBoxSxpb210M2Sum_10->Write();
	// hCuBoxSxth232M2Sum_001->Write();
	// hCuBoxSxth232M2Sum_01->Write();
	// hCuBoxSxth232M2Sum_1->Write();
	// hCuBoxSxth232M2Sum_10->Write();	
	// hCuBoxSxu238M2Sum_001->Write();
	// hCuBoxSxu238M2Sum_01->Write();
	// hCuBoxSxu238M2Sum_1->Write();
	// hCuBoxSxu238M2Sum_10->Write();

	// hCuFrameSth232M2Sum_1->Write();
	// hCuFrameSu238M2Sum_1->Write();

	// hCuFrameSxpb210M2Sum_001->Write();
	// hCuFrameSxpb210M2Sum_01->Write();
	// hCuFrameSxpb210M2Sum_1->Write();
	// hCuFrameSxpb210M2Sum_10->Write();
	// hCuFrameSxth232M2Sum_001->Write();
	// hCuFrameSxth232M2Sum_01->Write();
	// hCuFrameSxth232M2Sum_1->Write();
	// hCuFrameSxth232M2Sum_10->Write();	
	// hCuFrameSxu238M2Sum_001->Write();
	// hCuFrameSxu238M2Sum_01->Write();
	// hCuFrameSxu238M2Sum_1->Write();
	// hCuFrameSxu238M2Sum_10->Write();

	// hCuBox_CuFramepb210M2Sum_001->Write();
	// hCuBox_CuFramepb210M2Sum_01->Write();
	// hCuBox_CuFramepb210M2Sum_1->Write();
	// hCuBox_CuFramepb210M2Sum_10->Write();
	// hCuBox_CuFrameth232M2Sum_001->Write();
	// hCuBox_CuFrameth232M2Sum_01->Write();
	// hCuBox_CuFrameth232M2Sum_1->Write();
	// hCuBox_CuFrameth232M2Sum_10->Write();
	// hCuBox_CuFrameu238M2Sum_001->Write();
	// hCuBox_CuFrameu238M2Sum_01->Write();
	// hCuBox_CuFrameu238M2Sum_1->Write();
	// hCuBox_CuFrameu238M2Sum_10->Write();

	file2->Write();











}