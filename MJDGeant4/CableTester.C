// This macro tests the generation of the source locations in ROOT
// Instead of using G4UniformRand, use TRandom3 here
// Using TVector3 here instead of G4ThreeVector

void Test()
{
	gStyle->SetOptStat(0);

  // std::random_device rd;     // only used once to initialise (seed) engine
  // std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
  // std::uniform_int_distribution<int> uni(0,max); // guaranteed unbiased

  // auto random_integer = uni(rng);
  TRandom3 *fRand = new TRandom3(0);

  double fCableRadius = 0.5; // Just a placeholder so far

  double fPositionX, fPositionY;

  // TH1D *hTest = new TH1D("hTest", "hTest", 10, 0, 10);
  TH2D *hTest = new TH2D("hTest", "hTest", 500, -1, 1, 500, -1, 1);
  
  double randRadius;
  double randAngle;

  for(int i = 0; i < 10000000; i ++)
  {
    randRadius = fCableRadius*fCableRadius*fRand->Rndm();
  	randAngle = 2*TMath::Pi()*fRand->Rndm();
  	
  	fPositionX = TMath::Sqrt( randRadius ) * TMath::Cos( randAngle );
  	fPositionY = TMath::Sqrt( randRadius ) * TMath::Sin( randAngle );

  	hTest->Fill(fPositionX, fPositionY);

  }
  hTest->Draw("colz");
}

// All units are in cm here
void CableTester()
{
	gStyle->SetOptStat(0);

  double fPositionX, fPositionY, fPositionZ;
  TVector3 *fCableOffset[14]; // Centers of cables wrt center of cold plate, even is signal and odd is HV


  double fCableRadius = 0.05; // Just a placeholder so far
  double fCablePosition[4]; // Position o
  
  // These aren't good so far...
  double fCableCenter[4] = {-9.0/2,-9.0/2,-9.0/2,-9.0/2}; // Centers of signal cables
  double fHVCenter[4] = {-9.0/2,-9.0/2,-9.0/2,-9.0/2}; // Centers of HV cables, assuming the same as signal for now
  
  double fCableLength[4] = {9.0/2, 9.0/2, 9.0/2, 9.0/2}; // half length of cable... isn't this the same as the center?
  double fHVLength[4] = {9.0/2, 9.0/2, 9.0/2, 9.0/2};


  TRandom3 *fRand = new TRandom3(0);

  // Offsets are listed with signal before HV (signal is left of string)
  // Units are in inches?
  fCableOffset[0] = new TVector3(-2.54*1.839, 2.54*0.560, fCableCenter[0]); // P1
  fCableOffset[1] = new TVector3(-2.54*1.839, -2.54*0.560, fHVCenter[0]);  

  fCableOffset[2] = new TVector3(2.54*3.453, 2.54*3.846, fCableCenter[0]); // P2
  fCableOffset[3] = new TVector3(2.54*2.350, 2.54*3.620, fHVCenter[0]);

  fCableOffset[4] = new TVector3(-2.54*0.089, 2.54*5.576, fCableCenter[0]); // P3
  fCableOffset[5] = new TVector3(-2.54*1.186, 2.54*5.350, fHVCenter[0]);

  fCableOffset[6] = new TVector3(-2.54*3.356, 2.54*3.846, fCableCenter[0]); // P4
  fCableOffset[7] = new TVector3(-2.54*4.387, 2.54*3.408, fHVCenter[0]);

  fCableOffset[8] = new TVector3(-2.54*4.387, -2.54*3.408, fCableCenter[0]); // P5
  fCableOffset[9] = new TVector3(-2.54*3.356, -2.54*3.846, fHVCenter[0]);

  fCableOffset[10] = new TVector3(-2.54*1.186, -2.54*5.350, fCableCenter[0]); // P6
  fCableOffset[11] = new TVector3(-2.54*0.089, -2.54*5.576, fHVCenter[0]);

  fCableOffset[12] = new TVector3(2.54*2.350, -2.54*3.620, fCableCenter[0]); // P7
  fCableOffset[13] = new TVector3(2.54*3.453, -2.54*3.846, fHVCenter[0]);

  // Works
  // for(int i = 0; i < 14; i++)
  // {
  	// cout << fCableOffset[i](0) << "\t" << fCableOffset[i](1) << "\t" << fCableOffset(2) << endl;
  // }


  TH3D *hTest = new TH3D("hTest", "hTest", 200, -15, 10, 200, -20, 20, 200, -10, 10); 
  TH2D *hTestXY = new TH2D("hTestXY", "hTestXY", 1000, -15, 10, 1000, -20, 20);
  TH2D *hTestXZ = new TH2D("hTestXZ", "hTestXZ", 1000, -15, 10, 200, -10, 10);
  TH2D *hTestYZ = new TH2D("hTestYZ", "hTestYZ", 1000, -15, 15, 200, -10, 10);

  TH1D *hTestX = new TH1D("hTestX", "hTestX", 1000, -15, 10);
  TH1D *hTestY = new TH1D("hTestY", "hTestY", 1000, -20, 20);
  TH1D *hTestZ = new TH1D("hTestZ", "hTestZ", 200, -10, 10);


  // TGraph2D *gTest = new TGraph2D();

  int j;
  double randAngle;
  double randRadius;

  // Fill histogram
  for(int i = 0; i < 5000000; i++)
  {
	// Choose random integer for cable
  	j = fRand->Integer(14);

  	randRadius = fCableRadius*fCableRadius*fRand->Rndm();
  	randAngle = 2*TMath::Pi()*fRand->Rndm();

  	// Choose a random XY point along a disk
  	fPositionX = TMath::Sqrt( randRadius ) * TMath::Cos( randAngle ) + fCableOffset[j](0);
  	fPositionY = TMath::Sqrt( randRadius ) * TMath::Sin( randAngle ) + fCableOffset[j](1);
	
	fPositionZ =  (1. - 2.*fRand->Rndm())*fCableLength[0] + fCableOffset[j](2);
  	

  	hTest->Fill(fPositionX, fPositionY, fPositionZ);
  	hTestXY->Fill(fPositionX, fPositionY);
  	hTestXZ->Fill(fPositionX, fPositionZ);
  	hTestYZ->Fill(fPositionY, fPositionZ);

  	hTestX->Fill(fPositionX);
  	hTestY->Fill(fPositionY);
	hTestZ->Fill(fPositionZ);
  }

  TCanvas *cTest = new TCanvas("cTest", "cTest", 1200, 800);
  hTest->GetXaxis()->SetTitle("X");
  hTest->GetYaxis()->SetTitle("Y");
  hTest->GetZaxis()->SetTitle("Z");
  hTest->Draw("colz");


  TCanvas *cTestXY = new TCanvas("cTestXY", "cTestXY", 1200, 800);
  hTestXY->GetXaxis()->SetTitle("X");
  hTestXY->GetYaxis()->SetTitle("Y");
  hTestXY->Draw("colz");


  TCanvas *cTestXZ = new TCanvas("cTestXZ", "cTestXZ", 1200, 800);
  hTestXZ->GetXaxis()->SetTitle("X");
  hTestXZ->GetYaxis()->SetTitle("Z");
  hTestXZ->Draw("colz");

  TCanvas *cTestYZ = new TCanvas("cTestYZ", "cTestYZ", 1200, 800);
  hTestYZ->GetXaxis()->SetTitle("Y");
  hTestYZ->GetYaxis()->SetTitle("Z");
  hTestYZ->Draw("colz");


  TCanvas *cSingleAxis = new TCanvas("cSingleAxis", "cSingleAxis", 1200, 800);
  cSingleAxis->Divide(1,3);
  cSingleAxis->cd(1);
  hTestX->Draw();
  cSingleAxis->cd(2);
  hTestY->Draw();
  cSingleAxis->cd(3);
  hTestZ->Draw();
}