// Macro takes outputted g2root files and saves the histograms (M1 and M2 for now) into a ROOT file for easier usage
// File names are subject to change (unfortunately) - Nov 5, 2014

TChain *LoadMC(std::string dDir, std::string dLocation, std::string dSource)
{
    TChain *outTree = new TChain("outTree");
    outTree->Add(Form("%s%s-%s.root", dDir.c_str(), dLocation.c_str(), dSource.c_str()));

    return outTree;

}


void NormalizePDFPair(TH1D *h1, TH1D *h2, int minE, int maxE)
{
  double dIntegral = 0;
  int dBinSize = 2;

  // bin 0 = underflow, bin dNBins = last bin with upper-edge xup Excluded
  dIntegral = h1->Integral(minE/dBinSize, maxE/dBinSize);
  // cout << "Integral for " << h1->GetTitle() << " :" << dIntegral << endl;


  // Make sure integral isn't 0
  // If it is 0, clear model... 
  if(dIntegral == 0)
  {
    cout << Form("Integral of %s is 0, resetting histogram", h1->GetName()) << endl;
    h1->Reset();
    h2->Reset();
  }

  if(dIntegral != 0)
  {
    h1->Scale(1.0/dIntegral);
    h2->Scale(1.0/dIntegral);
  }
}

void SaveHistogramsBulk()
{
	std::string sDataDir = "/cuore/data/simulation/CUORE0/t14.08/production_g2root-r349/ntp/Bulk/";

	double dMinEnergy = 0;
	double dMaxEnergy = 10000;
	int dBinSize = 2;
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

	TH1D *hCuBoxco58M1 = new TH1D("hCuBoxco58M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxco60M1 = new TH1D("hCuBoxco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxcs137M1 = new TH1D("hCuBoxcs137M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxk40M1 = new TH1D("hCuBoxk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxmn54M1 = new TH1D("hCuBoxmn54M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxpb210M1 = new TH1D("hCuBoxpb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxth232M1 = new TH1D("hCuBoxth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxu238M1 = new TH1D("hCuBoxu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuFrameco58M1 = new TH1D("hCuFrameco58M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameco60M1 = new TH1D("hCuFrameco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramecs137M1 = new TH1D("hCuFramecs137M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramek40M1 = new TH1D("hCuFramek40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramemn54M1 = new TH1D("hCuFramemn54M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramepb210M1 = new TH1D("hCuFramepb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameth232M1 = new TH1D("hCuFrameth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameu238M1 = new TH1D("hCuFrameu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hIVCco60M1 = new TH1D("hIVCco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCk40M1 = new TH1D("hIVCk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCth232M1 = new TH1D("hIVCth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCu238M1 = new TH1D("hIVCu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hMBco60M1 = new TH1D("hMBco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBk40M1 = new TH1D("hMBk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBth232M1 = new TH1D("hMBth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBu238M1 = new TH1D("hMBu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hMCu238M1 = new TH1D("hMCu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hOVCco60M1 = new TH1D("hOVCco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCk40M1 = new TH1D("hOVCk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCth232M1 = new TH1D("hOVCth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCu238M1 = new TH1D("hOVCu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hPbRombi207M1 = new TH1D("hPbRombi207M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomco60M1 = new TH1D("hPbRomco60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomcs137M1 = new TH1D("hPbRomcs137M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomk40M1 = new TH1D("hPbRomk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRompb210M1 = new TH1D("hPbRompb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomth232M1 = new TH1D("hPbRomth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomu238M1 = new TH1D("hPbRomu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hSIk40M1 = new TH1D("hSIk40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hSIth232M1 = new TH1D("hSIth232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hSIu238M1 = new TH1D("hSIu238M1", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hTeO20nuM1 = new TH1D("hTeO20nuM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO22nuM1 = new TH1D("hTeO22nuM1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2co60M1 = new TH1D("hTeO2co60M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2k40M1 = new TH1D("hTeO2k40M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2pb210M1 = new TH1D("hTeO2pb210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2po210M1 = new TH1D("hTeO2po210M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2te125M1 = new TH1D("hTeO2te125M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232M1 = new TH1D("hTeO2th232M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th228M1 = new TH1D("hTeO2th228M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226M1 = new TH1D("hTeO2ra226M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2rn222M1 = new TH1D("hTeO2rn222M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238M1 = new TH1D("hTeO2u238M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th230M1 = new TH1D("hTeO2th230M1", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u234M1 = new TH1D("hTeO2u234M1", "", dNBins, dMinEnergy, dMaxEnergy);


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

	TH1D *hCuBoxco58M2 = new TH1D("hCuBoxco58M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxco60M2 = new TH1D("hCuBoxco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxcs137M2 = new TH1D("hCuBoxcs137M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxk40M2 = new TH1D("hCuBoxk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxmn54M2 = new TH1D("hCuBoxmn54M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxpb210M2 = new TH1D("hCuBoxpb210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxth232M2 = new TH1D("hCuBoxth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuBoxu238M2 = new TH1D("hCuBoxu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hCuFrameco58M2 = new TH1D("hCuFrameco58M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameco60M2 = new TH1D("hCuFrameco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramecs137M2 = new TH1D("hCuFramecs137M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramek40M2 = new TH1D("hCuFramek40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramemn54M2 = new TH1D("hCuFramemn54M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFramepb210M2 = new TH1D("hCuFramepb210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameth232M2 = new TH1D("hCuFrameth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hCuFrameu238M2 = new TH1D("hCuFrameu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hIVCco60M2 = new TH1D("hIVCco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCk40M2 = new TH1D("hIVCk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCth232M2 = new TH1D("hIVCth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hIVCu238M2 = new TH1D("hIVCu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hMBco60M2 = new TH1D("hMBco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBk40M2 = new TH1D("hMBk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBth232M2 = new TH1D("hMBth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hMBu238M2 = new TH1D("hMBu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hMCu238M2 = new TH1D("hMCu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hOVCco60M2 = new TH1D("hOVCco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCk40M2 = new TH1D("hOVCk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCth232M2 = new TH1D("hOVCth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hOVCu238M2 = new TH1D("hOVCu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hPbRombi207M2 = new TH1D("hPbRombi207M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomco60M2 = new TH1D("hPbRomco60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomcs137M2 = new TH1D("hPbRomcs137M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomk40M2 = new TH1D("hPbRomk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRompb210M2 = new TH1D("hPbRompb210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomth232M2 = new TH1D("hPbRomth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hPbRomu238M2 = new TH1D("hPbRomu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hSIk40M2 = new TH1D("hSIk40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hSIth232M2 = new TH1D("hSIth232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hSIu238M2 = new TH1D("hSIu238M2", "", dNBins, dMinEnergy, dMaxEnergy);

	TH1D *hTeO20nuM2 = new TH1D("hTeO20nuM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO22nuM2 = new TH1D("hTeO22nuM2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2co60M2 = new TH1D("hTeO2co60M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2k40M2 = new TH1D("hTeO2k40M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2pb210M2 = new TH1D("hTeO2pb210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2po210M2 = new TH1D("hTeO2po210M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2te125M2 = new TH1D("hTeO2te125M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th232M2 = new TH1D("hTeO2th232M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th228M2 = new TH1D("hTeO2th228M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2ra226M2 = new TH1D("hTeO2ra226M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2rn222M2 = new TH1D("hTeO2rn222M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u238M2 = new TH1D("hTeO2u238M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2th230M2 = new TH1D("hTeO2th230M2", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2u234M2 = new TH1D("hTeO2u234M2", "", dNBins, dMinEnergy, dMaxEnergy);

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

	TChain *tCuBoxco58 = LoadMC(sDataDir.c_str(), "CuBox", "co58");
	TChain *tCuBoxco60 = LoadMC(sDataDir.c_str(), "CuBox", "co60");
	TChain *tCuBoxcs137 = LoadMC(sDataDir.c_str(), "CuBox", "cs137");
	TChain *tCuBoxk40 = LoadMC(sDataDir.c_str(), "CuBox", "k40");
	TChain *tCuBoxmn54 = LoadMC(sDataDir.c_str(), "CuBox", "mn54");
	TChain *tCuBoxpb210 = LoadMC(sDataDir.c_str(), "CuBox", "pb210");
	TChain *tCuBoxth232 = LoadMC(sDataDir.c_str(), "CuBox", "th232");
	TChain *tCuBoxu238 = LoadMC(sDataDir.c_str(), "CuBox", "u238");

	TChain *tCuFrameco58 = LoadMC(sDataDir.c_str(), "CuFrame", "co58");
	TChain *tCuFrameco60 = LoadMC(sDataDir.c_str(), "CuFrame", "co60");
	TChain *tCuFramecs137 = LoadMC(sDataDir.c_str(), "CuFrame", "cs137");
	TChain *tCuFramek40 = LoadMC(sDataDir.c_str(), "CuFrame", "k40");
	TChain *tCuFramemn54 = LoadMC(sDataDir.c_str(), "CuFrame", "mn54");
	TChain *tCuFramepb210 = LoadMC(sDataDir.c_str(), "CuFrame", "pb210");
	TChain *tCuFrameth232 = LoadMC(sDataDir.c_str(), "CuFrame", "th232");
	TChain *tCuFrameu238 = LoadMC(sDataDir.c_str(), "CuFrame", "u238");

	TChain *tIVCco60 = LoadMC(sDataDir.c_str(), "IVC", "co60");
	TChain *tIVCk40 = LoadMC(sDataDir.c_str(), "IVC", "k40");
	TChain *tIVCth232 = LoadMC(sDataDir.c_str(), "IVC", "th232");
	TChain *tIVCu238 = LoadMC(sDataDir.c_str(), "IVC", "u238");

	TChain *tOVCco60 = LoadMC(sDataDir.c_str(), "OVC", "co60");
	TChain *tOVCk40 = LoadMC(sDataDir.c_str(), "OVC", "k40");
	TChain *tOVCth232 = LoadMC(sDataDir.c_str(), "OVC", "th232");
	TChain *tOVCu238 = LoadMC(sDataDir.c_str(), "OVC", "u238");

	TChain *tMBco60 = LoadMC(sDataDir.c_str(), "MB", "co60");
	TChain *tMBk40 = LoadMC(sDataDir.c_str(), "MB", "k40");
	TChain *tMBth232 = LoadMC(sDataDir.c_str(), "MB", "th232");
	TChain *tMBu238 = LoadMC(sDataDir.c_str(), "MB", "u238");

	TChain *tPbRombi207 = LoadMC(sDataDir.c_str(), "PbRom", "bi207");
	TChain *tPbRomco60 = LoadMC(sDataDir.c_str(), "PbRom", "co60");
	TChain *tPbRomcs137 = LoadMC(sDataDir.c_str(), "PbRom", "cs137");
	TChain *tPbRomk40 = LoadMC(sDataDir.c_str(), "PbRom", "k40");
	TChain *tPbRompb210 = LoadMC(sDataDir.c_str(), "PbRom", "pb210");
	TChain *tPbRomth232 = LoadMC(sDataDir.c_str(), "PbRom", "th232");
	TChain *tPbRomu238 = LoadMC(sDataDir.c_str(), "PbRom", "u238");

	TChain *tTeO20nu = LoadMC(sDataDir.c_str(), "TeO2", "0nu");
	TChain *tTeO22nu = LoadMC(sDataDir.c_str(), "TeO2", "2nu");
	TChain *tTeO2co60 = LoadMC(sDataDir.c_str(), "TeO2", "co60");
	TChain *tTeO2k40 = LoadMC(sDataDir.c_str(), "TeO2", "k40");
	TChain *tTeO2pb210 = LoadMC(sDataDir.c_str(), "TeO2", "pb210");
	TChain *tTeO2po210 = LoadMC(sDataDir.c_str(), "TeO2", "po210");
	TChain *tTeO2te125 = LoadMC(sDataDir.c_str(), "TeO2", "te125m");
	TChain *tTeO2th232 = LoadMC(sDataDir.c_str(), "TeO2", "th232");
	TChain *tTeO2th228 = LoadMC(sDataDir.c_str(), "TeO2", "th232-th228");
	TChain *tTeO2ra226 = LoadMC(sDataDir.c_str(), "TeO2", "u238-ra226");
	TChain *tTeO2rn222 = LoadMC(sDataDir.c_str(), "TeO2", "u238-rn222");
	TChain *tTeO2u238 = LoadMC(sDataDir.c_str(), "TeO2", "u238");
	TChain *tTeO2th230 = LoadMC(sDataDir.c_str(), "TeO2", "u238-th230");
	TChain *tTeO2u234 = LoadMC(sDataDir.c_str(), "TeO2", "u238-u234");

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

	tCuBoxco58->Project("hCuBoxco58M1", "Ener2", "Multiplicity==1");
	tCuBoxco60->Project("hCuBoxco60M1", "Ener2", "Multiplicity==1");
	tCuBoxcs137->Project("hCuBoxcs137M1", "Ener2", "Multiplicity==1");
	tCuBoxk40->Project("hCuBoxk40M1", "Ener2", "Multiplicity==1");
	tCuBoxmn54->Project("hCuBoxmn54M1", "Ener2", "Multiplicity==1");
	tCuBoxpb210->Project("hCuBoxpb210M1", "Ener2", "Multiplicity==1");
	tCuBoxth232->Project("hCuBoxth232M1", "Ener2", "Multiplicity==1");
	tCuBoxu238->Project("hCuBoxu238M1", "Ener2", "Multiplicity==1");	

	tCuFrameco58->Project("hCuFrameco58M1", "Ener2", "Multiplicity==1");
	tCuFrameco60->Project("hCuFrameco60M1", "Ener2", "Multiplicity==1");
	tCuFramecs137->Project("hCuFramecs137M1", "Ener2", "Multiplicity==1");
	tCuFramek40->Project("hCuFramek40M1", "Ener2", "Multiplicity==1");
	tCuFramemn54->Project("hCuFramemn54M1", "Ener2", "Multiplicity==1");
	tCuFramepb210->Project("hCuFramepb210M1", "Ener2", "Multiplicity==1");
	tCuFrameth232->Project("hCuFrameth232M1", "Ener2", "Multiplicity==1");
	tCuFrameu238->Project("hCuFrameu238M1", "Ener2", "Multiplicity==1");

	tIVCco60->Project("hIVCco60M1", "Ener2", "Multiplicity==1");
	tIVCk40->Project("hIVCk40M1", "Ener2", "Multiplicity==1");
	tIVCth232->Project("hIVCth232M1", "Ener2", "Multiplicity==1");
	tIVCu238->Project("hIVCu238M1", "Ener2", "Multiplicity==1");

	tOVCco60->Project("hOVCco60M1", "Ener2", "Multiplicity==1");
	tOVCk40->Project("hOVCk40M1", "Ener2", "Multiplicity==1");
	tOVCth232->Project("hOVCth232M1", "Ener2", "Multiplicity==1");
	tOVCu238->Project("hOVCu238M1", "Ener2", "Multiplicity==1");

	tMBco60->Project("hMBco60M1", "Ener2", "Multiplicity==1");
	tMBk40->Project("hMBk40M1", "Ener2", "Multiplicity==1");
	tMBth232->Project("hMBth232M1", "Ener2", "Multiplicity==1");
	tMBu238->Project("hMBu238M1", "Ener2", "Multiplicity==1");

	tPbRombi207->Project("hPbRombi207M1", "Ener2", "Multiplicity==1");
	tPbRomco60->Project("hPbRomco60M1", "Ener2", "Multiplicity==1");
	tPbRomcs137->Project("hPbRomcs137M1", "Ener2", "Multiplicity==1");
	tPbRomk40->Project("hPbRomk40M1", "Ener2", "Multiplicity==1");
	tPbRompb210->Project("hPbRompb210M1", "Ener2", "Multiplicity==1");
	tPbRomth232->Project("hPbRomth232M1", "Ener2", "Multiplicity==1");
	tPbRomu238->Project("hPbRomu238M1", "Ener2", "Multiplicity==1");


	tTeO20nu->Project("hTeO20nuM1", "Ener2", "Multiplicity==1");
	tTeO22nu->Project("hTeO22nuM1", "Ener2", "Multiplicity==1");
	tTeO2co60->Project("hTeO2co60M1", "Ener2", "Multiplicity==1");
	tTeO2k40->Project("hTeO2k40M1", "Ener2", "Multiplicity==1");
	tTeO2pb210->Project("hTeO2pb210M1", "Ener2", "Multiplicity==1");
	tTeO2po210->Project("hTeO2po210M1", "Ener2", "Multiplicity==1");
	tTeO2te125->Project("hTeO2te125M1", "Ener2", "Multiplicity==1");
	tTeO2th232->Project("hTeO2th232M1", "Ener2", "Multiplicity==1");
	tTeO2th228->Project("hTeO2th228M1", "Ener2", "Multiplicity==1");
	tTeO2ra226->Project("hTeO2ra226M1", "Ener2", "Multiplicity==1");
	tTeO2rn222->Project("hTeO2rn222M1", "Ener2", "Multiplicity==1");
	tTeO2u238->Project("hTeO2u238M1", "Ener2", "Multiplicity==1");
	tTeO2th230->Project("hTeO2th230M1", "Ener2", "Multiplicity==1");
	tTeO2u234->Project("hTeO2u234M1", "Ener2", "Multiplicity==1");

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

	tCuBoxco58->Project("hCuBoxco58M2", "Ener2", "Multiplicity==2");
	tCuBoxco60->Project("hCuBoxco60M2", "Ener2", "Multiplicity==2");
	tCuBoxcs137->Project("hCuBoxcs137M2", "Ener2", "Multiplicity==2");
	tCuBoxk40->Project("hCuBoxk40M2", "Ener2", "Multiplicity==2");
	tCuBoxmn54->Project("hCuBoxmn54M2", "Ener2", "Multiplicity==2");
	tCuBoxpb210->Project("hCuBoxpb210M2", "Ener2", "Multiplicity==2");
	tCuBoxth232->Project("hCuBoxth232M2", "Ener2", "Multiplicity==2");
	tCuBoxu238->Project("hCuBoxu238M2", "Ener2", "Multiplicity==2");	

	tCuFrameco58->Project("hCuFrameco58M2", "Ener2", "Multiplicity==2");
	tCuFrameco60->Project("hCuFrameco60M2", "Ener2", "Multiplicity==2");
	tCuFramecs137->Project("hCuFramecs137M2", "Ener2", "Multiplicity==2");
	tCuFramek40->Project("hCuFramek40M2", "Ener2", "Multiplicity==2");
	tCuFramemn54->Project("hCuFramemn54M2", "Ener2", "Multiplicity==2");
	tCuFramepb210->Project("hCuFramepb210M2", "Ener2", "Multiplicity==2");
	tCuFrameth232->Project("hCuFrameth232M2", "Ener2", "Multiplicity==2");
	tCuFrameu238->Project("hCuFrameu238M2", "Ener2", "Multiplicity==2");

	tIVCco60->Project("hIVCco60M2", "Ener2", "Multiplicity==2");
	tIVCk40->Project("hIVCk40M2", "Ener2", "Multiplicity==2");
	tIVCth232->Project("hIVCth232M2", "Ener2", "Multiplicity==2");
	tIVCu238->Project("hIVCu238M2", "Ener2", "Multiplicity==2");

	tOVCco60->Project("hOVCco60M2", "Ener2", "Multiplicity==2");
	tOVCk40->Project("hOVCk40M2", "Ener2", "Multiplicity==2");
	tOVCth232->Project("hOVCth232M2", "Ener2", "Multiplicity==2");
	tOVCu238->Project("hOVCu238M2", "Ener2", "Multiplicity==2");

	tMBco60->Project("hMBco60M2", "Ener2", "Multiplicity==2");
	tMBk40->Project("hMBk40M2", "Ener2", "Multiplicity==2");
	tMBth232->Project("hMBth232M2", "Ener2", "Multiplicity==2");
	tMBu238->Project("hMBu238M2", "Ener2", "Multiplicity==2");

	tPbRombi207->Project("hPbRombi207M2", "Ener2", "Multiplicity==2");
	tPbRomco60->Project("hPbRomco60M2", "Ener2", "Multiplicity==2");
	tPbRomcs137->Project("hPbRomcs137M2", "Ener2", "Multiplicity==2");
	tPbRomk40->Project("hPbRomk40M2", "Ener2", "Multiplicity==2");
	tPbRompb210->Project("hPbRompb210M2", "Ener2", "Multiplicity==2");
	tPbRomth232->Project("hPbRomth232M2", "Ener2", "Multiplicity==2");
	tPbRomu238->Project("hPbRomu238M2", "Ener2", "Multiplicity==2");


	tTeO20nu->Project("hTeO20nuM2", "Ener2", "Multiplicity==2");
	tTeO22nu->Project("hTeO22nuM2", "Ener2", "Multiplicity==2");
	tTeO2co60->Project("hTeO2co60M2", "Ener2", "Multiplicity==2");
	tTeO2k40->Project("hTeO2k40M2", "Ener2", "Multiplicity==2");
	tTeO2pb210->Project("hTeO2pb210M2", "Ener2", "Multiplicity==2");
	tTeO2po210->Project("hTeO2po210M2", "Ener2", "Multiplicity==2");
	tTeO2te125->Project("hTeO2te125M2", "Ener2", "Multiplicity==2");
	tTeO2th232->Project("hTeO2th232M2", "Ener2", "Multiplicity==2");
	tTeO2th228->Project("hTeO2th228M2", "Ener2", "Multiplicity==2");
	tTeO2ra226->Project("hTeO2ra226M2", "Ener2", "Multiplicity==2");
	tTeO2rn222->Project("hTeO2rn222M2", "Ener2", "Multiplicity==2");
	tTeO2u238->Project("hTeO2u238M2", "Ener2", "Multiplicity==2");
	tTeO2th230->Project("hTeO2th230M2", "Ener2", "Multiplicity==2");
	tTeO2u234->Project("hTeO2u234M2", "Ener2", "Multiplicity==2");



	//////////////////////////////////

	NormalizePDFPair(h50mKco58M1, h50mKco58M2, 0, 10000);
	NormalizePDFPair(h50mKco60M1, h50mKco60M2, 0, 10000);
	NormalizePDFPair(h50mKcs137M1, h50mKcs137M2, 0, 10000);
	NormalizePDFPair(h50mKk40M1, h50mKk40M2, 0, 10000);
	NormalizePDFPair(h50mKmn54M1, h50mKmn54M2, 0, 10000);
	NormalizePDFPair(h50mKpb210M1, h50mKpb210M2, 0, 10000);
	NormalizePDFPair(h50mKth232M1, h50mKth232M2, 0, 10000);
	NormalizePDFPair(h50mKu238M1, h50mKu238M2, 0, 10000);

	NormalizePDFPair(h600mKco60M1, h600mKco60M2, 0, 10000);
	NormalizePDFPair(h600mKk40M1, h600mKk40M2, 0, 10000);
	NormalizePDFPair(h600mKth232M1, h600mKth232M2, 0, 10000);
	NormalizePDFPair(h600mKu238M1, h600mKu238M2, 0, 10000);

	NormalizePDFPair(hCuBoxco58M1, hCuBoxco58M2, 0, 10000);
	NormalizePDFPair(hCuBoxco60M1, hCuBoxco60M2, 0, 10000);
	NormalizePDFPair(hCuBoxcs137M1, hCuBoxcs137M2, 0, 10000);
	NormalizePDFPair(hCuBoxk40M1, hCuBoxk40M2, 0, 10000);
	NormalizePDFPair(hCuBoxmn54M1, hCuBoxmn54M2, 0, 10000);
	NormalizePDFPair(hCuBoxpb210M1, hCuBoxpb210M2, 0, 10000);
	NormalizePDFPair(hCuBoxth232M1, hCuBoxth232M2, 0, 10000);
	NormalizePDFPair(hCuBoxu238M1, hCuBoxu238M2, 0, 10000);

	NormalizePDFPair(hCuFrameco58M1, hCuFrameco58M2, 0, 10000);
	NormalizePDFPair(hCuFrameco60M1, hCuFrameco60M2, 0, 10000);
	NormalizePDFPair(hCuFramecs137M1, hCuFramecs137M2, 0, 10000);
	NormalizePDFPair(hCuFramek40M1, hCuFramek40M2, 0, 10000);
	NormalizePDFPair(hCuFramemn54M1, hCuFramemn54M2, 0, 10000);
	NormalizePDFPair(hCuFramepb210M1, hCuFramepb210M2, 0, 10000);
	NormalizePDFPair(hCuFrameth232M1, hCuFrameth232M2, 0, 10000);
	NormalizePDFPair(hCuFrameu238M1, hCuFrameu238M2, 0, 10000);

	NormalizePDFPair(hIVCco60M1, hIVCco60M2, 0, 10000);
	NormalizePDFPair(hIVCk40M1, hIVCk40M2, 0, 10000);
	NormalizePDFPair(hIVCth232M1, hIVCth232M2, 0, 10000);
	NormalizePDFPair(hIVCu238M1, hIVCu238M2, 0, 10000);

	NormalizePDFPair(hOVCco60M1, hOVCco60M2, 0, 10000);
	NormalizePDFPair(hOVCk40M1, hOVCk40M2, 0, 10000);
	NormalizePDFPair(hOVCth232M1, hOVCth232M2, 0, 10000);
	NormalizePDFPair(hOVCu238M1, hOVCu238M2, 0, 10000);

	NormalizePDFPair(hMBco60M1, hMBco60M2, 0, 10000);
	NormalizePDFPair(hMBk40M1, hMBk40M2, 0, 10000);
	NormalizePDFPair(hMBth232M1, hMBth232M2, 0, 10000);
	NormalizePDFPair(hMBu238M1, hMBu238M2, 0, 10000);

	NormalizePDFPair(hPbRombi207M1, hPbRombi207M2, 0, 10000);
	NormalizePDFPair(hPbRomco60M1, hPbRomco60M2, 0, 10000);
	NormalizePDFPair(hPbRomcs137M1, hPbRomcs137M2, 0, 10000);
	NormalizePDFPair(hPbRomk40M1, hPbRomk40M2, 0, 10000);
	NormalizePDFPair(hPbRompb210M1, hPbRompb210M2, 0, 10000);
	NormalizePDFPair(hPbRomth232M1, hPbRomth232M2, 0, 10000);
	NormalizePDFPair(hPbRomu238M1, hPbRomu238M2, 0, 10000);

	NormalizePDFPair(hTeO20nuM1, hTeO20nuM2, 0, 10000);
	NormalizePDFPair(hTeO22nuM1, hTeO22nuM2, 0, 10000);
	NormalizePDFPair(hTeO2co60M1, hTeO2co60M2, 0, 10000);
	NormalizePDFPair(hTeO2k40M1, hTeO2k40M2, 0, 10000);
	NormalizePDFPair(hTeO2pb210M1, hTeO2pb210M2, 0, 10000);
	NormalizePDFPair(hTeO2po210M1, hTeO2po210M2, 0, 10000);
	NormalizePDFPair(hTeO2te125M1, hTeO2te125M2, 0, 10000);
	NormalizePDFPair(hTeO2th232M1, hTeO2th232M2, 0, 10000);
	NormalizePDFPair(hTeO2th228M1, hTeO2th228M2, 0, 10000);
	NormalizePDFPair(hTeO2ra226M1, hTeO2ra226M2, 0, 10000);
	NormalizePDFPair(hTeO2rn222M1, hTeO2rn222M2, 0, 10000);
	NormalizePDFPair(hTeO2u238M1, hTeO2u238M2, 0, 10000);
	NormalizePDFPair(hTeO2th230M1, hTeO2th230M2, 0, 10000);
	NormalizePDFPair(hTeO2u234M1, hTeO2u234M2, 0, 10000);


	TFile *file1 = new TFile("MCProduction_Bulk.root", "RECREATE");

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

	hCuBoxco58M1->Write();
	hCuBoxco60M1->Write();
	hCuBoxcs137M1->Write();
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

	hIVCco60M1->Write();
	hIVCk40M1->Write();
	hIVCth232M1->Write();
	hIVCu238M1->Write();

	hOVCco60M1->Write();
	hOVCk40M1->Write();
	hOVCth232M1->Write();
	hOVCu238M1->Write();

	hMBco60M1->Write();
	hMBk40M1->Write();
	hMBth232M1->Write();
	hMBu238M1->Write();

	hPbRombi207M1->Write();
	hPbRomco60M1->Write();
	hPbRomcs137M1->Write();
	hPbRomk40M1->Write();
	hPbRompb210M1->Write();
	hPbRomth232M1->Write();
	hPbRomu238M1->Write();


	hTeO20nuM1->Write();
	hTeO22nuM1->Write();
	hTeO2co60M1->Write();
	hTeO2k40M1->Write();
	hTeO2pb210M1->Write();
	hTeO2po210M1->Write();
	hTeO2te125M1->Write();
	hTeO2th232M1->Write();
	hTeO2th228M1->Write();
	hTeO2ra226M1->Write();
	hTeO2rn222M1->Write();
	hTeO2u238M1->Write();
	hTeO2th230M1->Write();
	hTeO2u234M1->Write();

///////////////////

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

	hCuBoxco58M2->Write();
	hCuBoxco60M2->Write();
	hCuBoxcs137M2->Write();
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

	hIVCco60M2->Write();
	hIVCk40M2->Write();
	hIVCth232M2->Write();
	hIVCu238M2->Write();

	hOVCco60M2->Write();
	hOVCk40M2->Write();
	hOVCth232M2->Write();
	hOVCu238M2->Write();

	hMBco60M2->Write();
	hMBk40M2->Write();
	hMBth232M2->Write();
	hMBu238M2->Write();

	hPbRombi207M2->Write();
	hPbRomco60M2->Write();
	hPbRomcs137M2->Write();
	hPbRomk40M2->Write();
	hPbRompb210M2->Write();
	hPbRomth232M2->Write();
	hPbRomu238M2->Write();


	hTeO20nuM2->Write();
	hTeO22nuM2->Write();
	hTeO2co60M2->Write();
	hTeO2k40M2->Write();
	hTeO2pb210M2->Write();
	hTeO2po210M2->Write();
	hTeO2te125M2->Write();
	hTeO2th232M2->Write();
	hTeO2th228M2->Write();
	hTeO2ra226M2->Write();
	hTeO2rn222M2->Write();
	hTeO2u238M2->Write();
	hTeO2th230M2->Write();
	hTeO2u234M2->Write();



	file1->Write();

}





void SaveHistogramsSurface()
{

	std::string sDataDir = "/cuore/data/simulation/CUORE0/t14.08/production_g2root-r349/ntp/Sup/";

	double dMinEnergy = 0;
	double dMaxEnergy = 10000;
	int dBinSize = 2;
	int dNBins = (dMaxEnergy - dMinEnergy)/dBinSize;

	TH1D *hCuBoxSth232M1_1 = new TH1D("hCuBoxSth232M1_1", "", dNBins, dMinEnergy, dMaxEnergy);	
	TH1D *hCuBoxSu238M1_1 = new TH1D("hCuBoxSu238M1_1", "", dNBins, dMinEnergy, dMaxEnergy);

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

	TH1D *hCuFrameSth232M1_1 = new TH1D("hCuFrameSth232M1_1", "", dNBins, dMinEnergy, dMaxEnergy);	
	TH1D *hCuFrameSu238M1_1 = new TH1D("hCuFrameSu238M1_1", "", dNBins, dMinEnergy, dMaxEnergy);

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


	TH1D *hTeO2Spb210M1_01 = new TH1D("hTeO2Spb210M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Spo210M1_001 = new TH1D("hTeO2Spo210M1_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Spo210M1_01 = new TH1D("hTeO2Spo210M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sth232M1_01 = new TH1D("hTeO2Sth232M1_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Su238M1_01 = new TH1D("hTeO2Su238M1_01", "", dNBins, dMinEnergy, dMaxEnergy);


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

//////////////////////////////
	TH1D *hCuBoxSth232M2_1 = new TH1D("hCuBoxSth232M2_1", "", dNBins, dMinEnergy, dMaxEnergy);	
	TH1D *hCuBoxSu238M2_1 = new TH1D("hCuBoxSu238M2_1", "", dNBins, dMinEnergy, dMaxEnergy);

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

	TH1D *hCuFrameSth232M2_1 = new TH1D("hCuFrameSth232M2_1", "", dNBins, dMinEnergy, dMaxEnergy);	
	TH1D *hCuFrameSu238M2_1 = new TH1D("hCuFrameSu238M2_1", "", dNBins, dMinEnergy, dMaxEnergy);

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

	TH1D *hTeO2Spb210M2_01 = new TH1D("hTeO2Spb210M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Spo210M2_001 = new TH1D("hTeO2Spo210M2_001", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Spo210M2_01 = new TH1D("hTeO2Spo210M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Sth232M2_01 = new TH1D("hTeO2Sth232M2_01", "", dNBins, dMinEnergy, dMaxEnergy);
	TH1D *hTeO2Su238M2_01 = new TH1D("hTeO2Su238M2_01", "", dNBins, dMinEnergy, dMaxEnergy);

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


	TChain *tCuBoxSth232_1 = LoadMC(sDataDir.c_str(), "CuBoxS", "th232-1");
	TChain *tCuBoxSu238_1 = LoadMC(sDataDir.c_str(), "CuBoxS", "u238-1");


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

	TChain *tCuFrameSth232_1 = LoadMC(sDataDir.c_str(), "CuFrameS", "th232-1");
	TChain *tCuFrameSu238_1 = LoadMC(sDataDir.c_str(), "CuFrameS", "u238-1");

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

	TChain *tTeO2Spb210_01 = LoadMC(sDataDir.c_str(), "TeO2S", "pb210-.1");
	TChain *tTeO2Spo210_001 = LoadMC(sDataDir.c_str(), "TeO2S", "po210-.01");
	TChain *tTeO2Spo210_01 = LoadMC(sDataDir.c_str(), "TeO2S", "po210-.1");
	TChain *tTeO2Sth232_01 = LoadMC(sDataDir.c_str(), "TeO2S", "th232-.1");
	TChain *tTeO2Su238_01 = LoadMC(sDataDir.c_str(), "TeO2S", "u238-.1");


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


///////////////////
	tCuBoxSth232_1->Project("hCuBoxSth232M1_1", "Ener2", "Multiplicity==1");
	tCuBoxSu238_1->Project("hCuBoxSu238M1_1", "Ener2", "Multiplicity==1");

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

	tCuFrameSth232_1->Project("hCuFrameSth232M1_1", "Ener2", "Multiplicity==1");
	tCuFrameSu238_1->Project("hCuFrameSu238M1_1", "Ener2", "Multiplicity==1");

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

	tTeO2Spb210_01->Project("hTeO2Spb210M1_01", "Ener2", "Multiplicity==1");
	tTeO2Spo210_001->Project("hTeO2Spo210M1_001", "Ener2", "Multiplicity==1");
	tTeO2Spo210_01->Project("hTeO2Spo210M1_01", "Ener2", "Multiplicity==1");
	tTeO2Sth232_01->Project("hTeO2Sth232M1_01", "Ener2", "Multiplicity==1");
	tTeO2Su238_01->Project("hTeO2Su238M1_01", "Ener2", "Multiplicity==1");

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

//////////////////////////////////
	tCuBoxSth232_1->Project("hCuBoxSth232M2_1", "Ener2", "Multiplicity==2");
	tCuBoxSu238_1->Project("hCuBoxSu238M2_1", "Ener2", "Multiplicity==2");

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


	tCuFrameSth232_1->Project("hCuFrameSth232M2_1", "Ener2", "Multiplicity==2");
	tCuFrameSu238_1->Project("hCuFrameSu238M2_1", "Ener2", "Multiplicity==2");

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

	tTeO2Spb210_01->Project("hTeO2Spb210M2_01", "Ener2", "Multiplicity==2");
	tTeO2Spo210_001->Project("hTeO2Spo210M2_001", "Ener2", "Multiplicity==2");
	tTeO2Spo210_01->Project("hTeO2Spo210M2_01", "Ener2", "Multiplicity==2");
	tTeO2Sth232_01->Project("hTeO2Sth232M2_01", "Ener2", "Multiplicity==2");
	tTeO2Su238_01->Project("hTeO2Su238M2_01", "Ener2", "Multiplicity==2");

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


////////////////////////////
	NormalizePDFPair(hCuBoxSth232M1_1, hCuBoxSth232M2_1, 0, 10000);
	NormalizePDFPair(hCuBoxSu238M1_1, hCuBoxSu238M2_1, 0, 10000);

	NormalizePDFPair(hCuBoxSxpb210M1_001, hCuBoxSxpb210M2_001, 0, 10000);
	NormalizePDFPair(hCuBoxSxpb210M1_01, hCuBoxSxpb210M2_01, 0, 10000);
	NormalizePDFPair(hCuBoxSxpb210M1_1, hCuBoxSxpb210M2_1, 0, 10000);
	NormalizePDFPair(hCuBoxSxpb210M1_10, hCuBoxSxpb210M2_10, 0, 10000);
	NormalizePDFPair(hCuBoxSxth232M1_001, hCuBoxSxth232M2_001, 0, 10000);
	NormalizePDFPair(hCuBoxSxth232M1_01, hCuBoxSxth232M2_01, 0, 10000);
	NormalizePDFPair(hCuBoxSxth232M1_1, hCuBoxSxth232M2_1, 0, 10000);
	NormalizePDFPair(hCuBoxSxth232M1_10, hCuBoxSxth232M2_10, 0, 10000);
	NormalizePDFPair(hCuBoxSxu238M1_001, hCuBoxSxu238M2_001, 0, 10000);
	NormalizePDFPair(hCuBoxSxu238M1_01, hCuBoxSxu238M2_01, 0, 10000);
	NormalizePDFPair(hCuBoxSxu238M1_1, hCuBoxSxu238M2_1, 0, 10000);
	NormalizePDFPair(hCuBoxSxu238M1_10, hCuBoxSxu238M2_10, 0, 10000);

	NormalizePDFPair(hCuFrameSth232M1_1, hCuFrameSth232M2_1, 0, 10000);
	NormalizePDFPair(hCuFrameSu238M1_1, hCuFrameSu238M2_1, 0, 10000);

	NormalizePDFPair(hCuFrameSxpb210M1_001, hCuFrameSxpb210M2_001, 0, 10000);
	NormalizePDFPair(hCuFrameSxpb210M1_01, hCuFrameSxpb210M2_01, 0, 10000);
	NormalizePDFPair(hCuFrameSxpb210M1_1, hCuFrameSxpb210M2_1, 0, 10000);
	NormalizePDFPair(hCuFrameSxpb210M1_10, hCuFrameSxpb210M2_10, 0, 10000);
	NormalizePDFPair(hCuFrameSxth232M1_001, hCuFrameSxth232M2_001, 0, 10000);
	NormalizePDFPair(hCuFrameSxth232M1_01, hCuFrameSxth232M2_01, 0, 10000);
	NormalizePDFPair(hCuFrameSxth232M1_1, hCuFrameSxth232M2_1, 0, 10000);
	NormalizePDFPair(hCuFrameSxth232M1_10, hCuFrameSxth232M2_10, 0, 10000);
	NormalizePDFPair(hCuFrameSxu238M1_001, hCuFrameSxu238M2_001, 0, 10000);
	NormalizePDFPair(hCuFrameSxu238M1_01, hCuFrameSxu238M2_01, 0, 10000);
	NormalizePDFPair(hCuFrameSxu238M1_1, hCuFrameSxu238M2_1, 0, 10000);
	NormalizePDFPair(hCuFrameSxu238M1_10, hCuFrameSxu238M2_10, 0, 10000);	

	NormalizePDFPair(hTeO2Spb210M1_01, hTeO2Spb210M2_01, 0, 10000);
	NormalizePDFPair(hTeO2Spo210M1_001, hTeO2Spo210M2_001, 0, 10000);
	NormalizePDFPair(hTeO2Spo210M1_01, hTeO2Spo210M2_01, 0, 10000);
	NormalizePDFPair(hTeO2Sth232M1_01, hTeO2Sth232M2_01, 0, 10000);
	NormalizePDFPair(hTeO2Su238M1_01, hTeO2Su238M2_01, 0, 10000);

	NormalizePDFPair(hTeO2Sxpb210M1_001, hTeO2Sxpb210M2_001, 0, 10000);
	NormalizePDFPair(hTeO2Sxpb210M1_01, hTeO2Sxpb210M2_01, 0, 10000);
	NormalizePDFPair(hTeO2Sxpb210M1_1, hTeO2Sxpb210M2_1, 0, 10000);
	NormalizePDFPair(hTeO2Sxpb210M1_10, hTeO2Sxpb210M2_10, 0, 10000);
	NormalizePDFPair(hTeO2Sxpo210M1_001, hTeO2Sxpo210M2_001, 0, 10000);
	NormalizePDFPair(hTeO2Sxpo210M1_01, hTeO2Sxpo210M2_01, 0, 10000);
	NormalizePDFPair(hTeO2Sxpo210M1_1, hTeO2Sxpo210M2_1, 0, 10000);
	// NormalizePDFPair(hTeO2Sxpo210M1_10, hTeO2Sxpo210M2_10, 0, 10000);
	NormalizePDFPair(hTeO2Sxth232M1_001, hTeO2Sxth232M2_001, 0, 10000);
	NormalizePDFPair(hTeO2Sxth232M1_01, hTeO2Sxth232M2_01, 0, 10000);
	NormalizePDFPair(hTeO2Sxth232M1_1, hTeO2Sxth232M2_1, 0, 10000);
	NormalizePDFPair(hTeO2Sxth232M1_10, hTeO2Sxth232M2_10, 0, 10000);
	NormalizePDFPair(hTeO2Sxu238M1_001, hTeO2Sxu238M2_001, 0, 10000);
	NormalizePDFPair(hTeO2Sxu238M1_01, hTeO2Sxu238M2_01, 0, 10000);
	NormalizePDFPair(hTeO2Sxu238M1_1, hTeO2Sxu238M2_1, 0, 10000);
	NormalizePDFPair(hTeO2Sxu238M1_10, hTeO2Sxu238M2_10, 0, 10000);	


	TFile *file2 = new TFile("MCProduction_Surface.root", "RECREATE");

	hCuBoxSth232M1_1->Write();
	hCuBoxSu238M1_1->Write();

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

	hCuFrameSth232M1_1->Write();
	hCuFrameSu238M1_1->Write();

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

	hTeO2Spb210M1_01->Write();
	hTeO2Spo210M1_001->Write();
	hTeO2Spo210M1_01->Write();
	hTeO2Sth232M1_01->Write();
	hTeO2Su238M1_01->Write();


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

///////////////////////////////
	hCuBoxSth232M2_1->Write();
	hCuBoxSu238M2_1->Write();

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

	hCuFrameSth232M2_1->Write();
	hCuFrameSu238M2_1->Write();

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

	hTeO2Spb210M2_01->Write();
	hTeO2Spo210M2_001->Write();
	hTeO2Spo210M2_01->Write();
	hTeO2Sth232M2_01->Write();
	hTeO2Su238M2_01->Write();


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


	file2->Write();

}


