// This macro tests the generation of the source locations in ROOT
// Instead of using G4UniformRand, use TRandom3 here

void Test()
{
	double fCableOffset[10][3];
	fCableOffset[0]

}


// All units are in cm here
void CableTester()
{

  double fPositionX, fPositionY, fPositionZ;
  double fCableOffset[14][3];
  TVector3 *fCableOffset[14];
  double fCableRadius = 0.005; // Just a placeholder so far
  double fCablePosition[4]; // Position of cables
  double fCableCenter[4] = {0,0,0,0}; // Centers of signal cables
  double fHVCenter[4]; // Centers of HV cables

  TRandom3 *fRand = new TRandom3(0);

  // Offsets are listed with signal before HV (signal is left of string)
  // Units are in inches?
  fCableOffset[0] = new TVector3() {2.54*1.839, 2.54*0.560, fCableCenter[0]}; // P1
  fCableOffset[1] = {2.54*1.839, -2.54*0.560, fCableCenter[1]};  

  fCableOffset[2] = {2.54*3.453, 2.54*3.846, fCableCenter[0]}; // P2
  fCableOffset[3] = {2.54*2.350, 2.54*3.620, fCableCenter[1]};

  fCableOffset[4] = {2.54*0.089, 2.54*5.576, fCableCenter[0]}; // P2
  fCableOffset[5] = {2.54*1.186, 2.54*5.350, fCableCenter[1]};

  fCableOffset[6] = {2.54*3.356, 2.54*3.846, fCableCenter[0]}; // P4
  fCableOffset[7] = {2.54*4.387, 2.54*3.408, fCableCenter[1]};

  fCableOffset[8] = {2.54*4.387, -2.54*3.408, fCableCenter[0]}; // P5
  fCableOffset[9] = {2.54*3.356, -2.54*3.846, fCableCenter[1]};

  fCableOffset[10] = {2.54*1.186, -2.54*5.350, fCableCenter[0]}; // P6
  fCableOffset[11] = {2.54*0.089, -2.54*5.576, fCableCenter[1]};

  fCableOffset[12] = {2.54*2.350, -2.54*3.620, fCableCenter[0]}; // P7
  fCableOffset[13] = {2.54*3.453, -2.54*3.846, fCableCenter[1]};

  TH3D *hTest = new TH3D("hTest", "hTest", 20, -10, 10, 20, -10, 10, 20, -10, 10); 
  TGraph2D *gTest = new TGraph2D();


  // Fill histogram
  for(int i = 0; i < 100000; i++)
  {
	// This chooses a random cable and generates a random point along the cable (cables are infinitely thin here)
  	fPositionZ =  (1. - 2.*fRand->Rndm())*fCableLength[i];
  
  	// Choose a random XY point along a disk
  	fPositionX = TMath::Sqrt( fCableRadius*fRand->Rndm() ) * TMath::Cos( 2*TMath::Pi()*fRand->Rndm() );
  	fPositionY = TMath::Sqrt( fCableRadius*fRand->Rndm() ) * TMath::Sin( 2*TMath::Pi()*fRand->Rndm() );
  	
  	hTest->Fill(fPositionX, fPositionY, fPositionZ);
  }

  TCanvas *cTest = new TCanvas("cTest", "cTest", 1200, 800);



  hTest->Draw();

}