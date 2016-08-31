// This macro tests the generation of the source locations in ROOT
// Instead of using G4UniformRand, use TRandom3 here
// Using TVector3 here instead of G4ThreeVector

  // Calculating cable locations:
  // 99 to -187 => Middle = -44 (simulation) = -120 (real) => Add 76
  // 287.0 mm total length currently => 280
  // Signal
  // 91+189, -126 center;  91+117, -89 center; 91+50, -56 center; 91+15, -38 => For the 4 detectors
  // 91+183, 91+128, 91+77, 91+24, 91-27 => For the 5 detectors, P
  // HV
  // 91+139, -100 center; 91+74, -67.5; 91+14, -37.5, 91-58, -1.5 => For the 4 detectors
  // 91+143, 91+91, 91+39, 91-5, 91-56 => For the 5 detectors

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
  TVector3 *fCableOffset[4][14]; // Centers of cables wrt center of cold plate, even is signal and odd is HV
  // TVector3 *fCableOffset[1][14];
  // TVector3 *fCableOffset3[14];
  // TVector3 *fCableOffset4[14];

  double fCableRadius = 0.5; // Just a placeholder so far
  double fCablePosition[4]; // Position o
  
  // These aren't good so far...
  double fCableCenter[4] = {-2.54*10.0/2, -2.54*8.0/2, -2.54*5.0/2, -2.54*2.5/2}; // Centers of signal cables
  double fHVCenter[4] = {-2.54*10.0/2, -2.54*8.0/2, -2.54*5.0/2, -2.54*2.5/2}; // Centers of HV cables, assuming the same as signal for now
  
  double fCableLength[4] = {2.54*10.0/2, 2.54*8.0/2, 2.54*5.0/2, 2.54*2.5/2}; // Half length of cable
  double fHVLength[4] = {2.54*10.0/2, 2.54*8.0/2, 2.54*5.0/2, 2.54*2.5/2};

  // Treat diving board as a rectangular box probably...
  double fDiveBoardL = 2.54*6.15;
  double fDiveBoardW = ;
  double fDiveBoardH = 2*fCableRadius; // should be basically same as cable radius right?


  TRandom3 *fRand = new TRandom3(0);

  // Offsets are listed with signal before HV (signal is left of string)
  // Units are in inches?
  fCableOffset[0][0] = new TVector3(-2.54*1.839, 2.54*0.560, fCableCenter[0]); // P1
  fCableOffset[0][1] = new TVector3(-2.54*1.839, -2.54*0.560, fHVCenter[0]);  
  fCableOffset[0][2] = new TVector3(2.54*3.453, 2.54*3.846, fCableCenter[0]); // P2
  fCableOffset[0][3] = new TVector3(2.54*2.350, 2.54*3.620, fHVCenter[0]);
  fCableOffset[0][4] = new TVector3(-2.54*0.089, 2.54*5.576, fCableCenter[0]); // P3
  fCableOffset[0][5] = new TVector3(-2.54*1.186, 2.54*5.350, fHVCenter[0]);
  fCableOffset[0][6] = new TVector3(-2.54*3.356, 2.54*3.846, fCableCenter[0]); // P4
  fCableOffset[0][7] = new TVector3(-2.54*4.387, 2.54*3.408, fHVCenter[0]);
  fCableOffset[0][8] = new TVector3(-2.54*4.387, -2.54*3.408, fCableCenter[0]); // P5
  fCableOffset[0][9] = new TVector3(-2.54*3.356, -2.54*3.846, fHVCenter[0]);
  fCableOffset[0][10] = new TVector3(-2.54*1.186, -2.54*5.350, fCableCenter[0]); // P6
  fCableOffset[0][11] = new TVector3(-2.54*0.089, -2.54*5.576, fHVCenter[0]);
  fCableOffset[0][12] = new TVector3(2.54*2.350, -2.54*3.620, fCableCenter[0]); // P7
  fCableOffset[0][13] = new TVector3(2.54*3.453, -2.54*3.846, fHVCenter[0]);

  fCableOffset[1][0] = new TVector3(-2.54*1.839, 2.54*0.560, fCableCenter[1]); // P1
  fCableOffset[1][1] = new TVector3(-2.54*1.839, -2.54*0.560, fHVCenter[1]);  
  fCableOffset[1][2] = new TVector3(2.54*3.453, 2.54*3.846, fCableCenter[1]); // P2
  fCableOffset[1][3] = new TVector3(2.54*2.350, 2.54*3.620, fHVCenter[1]);
  fCableOffset[1][4] = new TVector3(-2.54*0.089, 2.54*5.576, fCableCenter[1]); // P3
  fCableOffset[1][5] = new TVector3(-2.54*1.186, 2.54*5.350, fHVCenter[1]);
  fCableOffset[1][6] = new TVector3(-2.54*3.356, 2.54*3.846, fCableCenter[1]); // P4
  fCableOffset[1][7] = new TVector3(-2.54*4.387, 2.54*3.408, fHVCenter[1]);
  fCableOffset[1][8] = new TVector3(-2.54*4.387, -2.54*3.408, fCableCenter[1]); // P5
  fCableOffset[1][9] = new TVector3(-2.54*3.356, -2.54*3.846, fHVCenter[1]);
  fCableOffset[1][10] = new TVector3(-2.54*1.186, -2.54*5.350, fCableCenter[1]); // P6
  fCableOffset[1][11] = new TVector3(-2.54*0.089, -2.54*5.576, fHVCenter[1]);
  fCableOffset[1][12] = new TVector3(2.54*2.350, -2.54*3.620, fCableCenter[1]); // P7
  fCableOffset[1][13] = new TVector3(2.54*3.453, -2.54*3.846, fHVCenter[1]);


  fCableOffset[2][0] = new TVector3(-2.54*1.839, 2.54*0.560, fCableCenter[2]); // P1
  fCableOffset[2][1] = new TVector3(-2.54*1.839, -2.54*0.560, fHVCenter[2]);  
  fCableOffset[2][2] = new TVector3(2.54*3.453, 2.54*3.846, fCableCenter[2]); // P2
  fCableOffset[2][3] = new TVector3(2.54*2.350, 2.54*3.620, fHVCenter[2]);
  fCableOffset[2][4] = new TVector3(-2.54*0.089, 2.54*5.576, fCableCenter[2]); // P3
  fCableOffset[2][5] = new TVector3(-2.54*1.186, 2.54*5.350, fHVCenter[2]);
  fCableOffset[2][6] = new TVector3(-2.54*3.356, 2.54*3.846, fCableCenter[2]); // P4
  fCableOffset[2][7] = new TVector3(-2.54*4.387, 2.54*3.408, fHVCenter[2]);
  fCableOffset[2][8] = new TVector3(-2.54*4.387, -2.54*3.408, fCableCenter[2]); // P5
  fCableOffset[2][9] = new TVector3(-2.54*3.356, -2.54*3.846, fHVCenter[2]);
  fCableOffset[2][10] = new TVector3(-2.54*1.186, -2.54*5.350, fCableCenter[2]); // P6
  fCableOffset[2][11] = new TVector3(-2.54*0.089, -2.54*5.576, fHVCenter[2]);
  fCableOffset[2][12] = new TVector3(2.54*2.350, -2.54*3.620, fCableCenter[2]); // P7
  fCableOffset[2][13] = new TVector3(2.54*3.453, -2.54*3.846, fHVCenter[2]);

  fCableOffset[3][0] = new TVector3(-2.54*1.839, 2.54*0.560, fCableCenter[3]); // P1
  fCableOffset[3][1] = new TVector3(-2.54*1.839, -2.54*0.560, fHVCenter[3]);  
  fCableOffset[3][2] = new TVector3(2.54*3.453, 2.54*3.846, fCableCenter[3]); // P2
  fCableOffset[3][3] = new TVector3(2.54*2.350, 2.54*3.620, fHVCenter[3]);
  fCableOffset[3][4] = new TVector3(-2.54*0.089, 2.54*5.576, fCableCenter[3]); // P3
  fCableOffset[3][5] = new TVector3(-2.54*1.186, 2.54*5.350, fHVCenter[3]);
  fCableOffset[3][6] = new TVector3(-2.54*3.356, 2.54*3.846, fCableCenter[3]); // P4
  fCableOffset[3][7] = new TVector3(-2.54*4.387, 2.54*3.408, fHVCenter[3]);
  fCableOffset[3][8] = new TVector3(-2.54*4.387, -2.54*3.408, fCableCenter[3]); // P5
  fCableOffset[3][9] = new TVector3(-2.54*3.356, -2.54*3.846, fHVCenter[3]);
  fCableOffset[3][10] = new TVector3(-2.54*1.186, -2.54*5.350, fCableCenter[3]); // P6
  fCableOffset[3][11] = new TVector3(-2.54*0.089, -2.54*5.576, fHVCenter[3]);
  fCableOffset[3][12] = new TVector3(2.54*2.350, -2.54*3.620, fCableCenter[3]); // P7
  fCableOffset[3][13] = new TVector3(2.54*3.453, -2.54*3.846, fHVCenter[3]);

  // Works
  // for(int i = 0; i < 14; i++)
  // {
  	// cout << fCableOffset[i](0) << "\t" << fCableOffset[i](1) << "\t" << fCableOffset(2) << endl;
  // }


  TH3D *hTest = new TH3D("hTest", "hTest", 200, -15, 10, 200, -20, 20, 200, -30, 10); 
  TH2D *hTestXY = new TH2D("hTestXY", "hTestXY", 1000, -15, 10, 1000, -20, 20);
  TH2D *hTestXZ = new TH2D("hTestXZ", "hTestXZ", 1000, -15, 10, 200, -30, 10);
  TH2D *hTestYZ = new TH2D("hTestYZ", "hTestYZ", 1000, -15, 15, 200, -30, 10);

  TH1D *hTestX = new TH1D("hTestX", "hTestX", 1000, -15, 10);
  TH1D *hTestY = new TH1D("hTestY", "hTestY", 1000, -20, 20);
  TH1D *hTestZ = new TH1D("hTestZ", "hTestZ", 200, -30, 10);

  // TGraph2D *gTest = new TGraph2D();

  int j,k;
  double randAngle;
  double randRadius;

  // Fill histogram
  for(int i = 0; i < 10000000; i++)
  {
	// Choose random integer for cable
  	k = fRand->Integer(14); // random string position
  	j = fRand->Integer(4); // random cable position

  	randRadius = fCableRadius*fCableRadius*fRand->Rndm();
  	randAngle = 2*TMath::Pi()*fRand->Rndm();

  	// Choose a random XY point along a disk
  	fPositionX = TMath::Sqrt( randRadius ) * TMath::Cos( randAngle ) + fCableOffset[j][k](0);
  	fPositionY = TMath::Sqrt( randRadius ) * TMath::Sin( randAngle ) + fCableOffset[j][k](1);
	
	fPositionZ =  (1. - 2.*fRand->Rndm())*fCableLength[j] + fCableOffset[j][k](2);
  	

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
  hTest->Draw();


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