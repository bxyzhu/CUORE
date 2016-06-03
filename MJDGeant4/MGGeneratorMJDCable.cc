//---------------------------------------------------------------------------//
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                                                                           //
//                         MaGe Simulation                                   //
//                                                                           //
//      This code implementation is the intellectual property of the         //
//      MAJORANA and Gerda Collaborations. It is based on Geant4, an         //
//      intellectual property of the RD44 GEANT4 collaboration.              //
//                                                                           //
//                        *********************                              //
//                                                                           //
//    Neither the authors of this software system, nor their employing       //
//    institutes, nor the agencies providing financial support for this      //
//    work  make  any representation or  warranty, express or implied,       //
//    regarding this software system or assume any liability for its use.    //
//    By copying, distributing or modifying the Program (or any work based   //
//    on on the Program) you indicate your acceptance of this statement,     //
//    and all its terms.                                                     //
//                                                                           //
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//---------------------------------------------------------------------------//
//                                                          
// $Id: MGGeneratorMJDCable.cc,v 1.4 2014-07-14 09:24:09 mgmarino Exp $ 
//      
// CLASS IMPLEMENTATION:  MGGeneratorMJDCable.cc
//
//---------------------------------------------------------------------------//
/**
 * SPECIAL NOTES:
 */
// 
//---------------------------------------------------------------------------//
/**
 * AUTHOR: K. Vorren
 * CONTACT: 
 * FIRST SUBMISSION:
 * 
 * REVISION:
 *
 * 02-24-2016, Modified generator to read calibration source parameters
 *             from messenger commands, T. Caldwell 
 * 07-14-2014, Created, K. Vorren
 *
 */
//---------------------------------------------------------------------------//
//

#include <sstream>
#include "Randomize.hh"

#include "G4ios.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4RandomDirection.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Run.hh"
#include "G4ThreeVector.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandTree.hh"

#include "G4AffineTransform.hh"
#include "G4Geantino.hh"
#include "generators/MGGeneratorMJDCable.hh"
#include "generators/MGGeneratorMJDCableMessenger.hh"
// #include <CLHEP/Random/RandFlat.h>
#include "io/MGLogger.hh"

#include <math.h>
#define pi 3.14159265358979323846

//---------------------------------------------------------------------------//

#include "generators/MGGeneratorMJDCable.hh" 

//---------------------------------------------------------------------------//


using namespace CLHEP;

MGGeneratorMJDCable::MGGeneratorMJDCable()
{
  fGeneratorName = "MJDCable";
  fG4Messenger = new MGGeneratorMJDCableMessenger(this);
  fParticleGun = new G4ParticleGun(1);

  // Cryostat global position, from MJGeometryDemonstrator
  G4ThreeVector fCryo1Pos = G4ThreeVector(-8.1417 * 25.4 * mm, 0.0, 4.4265 * 25.4 * mm);
  G4double fCryo1Rot = pi / 2;
  G4RotationMatrix* rotationC = new G4RotationMatrix();
  rotationC->rotateZ(fCryo1Rot);

  G4ThreeVector fCryo2Pos = G4ThreeVector(-fCryo1Pos.x(), fCryo1Pos.y(), fCryo1Pos.z());
  G4double fCryo2Rot = 0.0;
  G4RotationMatrix* rotationD = new G4RotationMatrix();
  rotationD->rotateZ(fCryo2Rot);
  G4double eps = 0.01*mm;

  G4AffineTransform *assemAffine1 = new G4AffineTransform(rotationC,fCryo1Pos);
  G4AffineTransform *assemAffine2 = new G4AffineTransform(rotationD,fCryo2Pos);

  // Cold plate position w.r.t Cryostat
  G4ThreeVector *CPlocalPos = new G4ThreeVector(0, 0, -1.05*25.4*mm-eps);
  G4RotationMatrix *CPlocalRot = new G4RotationMatrix();
  CPlocalRot->rotateZ(pi);
  
  G4AffineTransform *CPaffine1 = new G4AffineTransform(CPlocalRot,*CPlocalPos);
  *CPaffine1 *= *assemAffine1;  
  fColdPlateOffset[0] = CPaffine1->NetTranslation();
  G4RotationMatrix *CPglobalRot1= new G4RotationMatrix(CPaffine1->NetRotation());
  fColdPlateOffset[0] *= *CPglobalRot1;

  G4AffineTransform *CPaffine2 = new G4AffineTransform(CPlocalRot,*CPlocalPos);
  *CPaffine2 *= *assemAffine2;  
  fColdPlateOffset[1] = CPaffine2->NetTranslation();
  G4RotationMatrix *CPglobalRot2= new G4RotationMatrix(CPaffine2->NetRotation());
  fColdPlateOffset[1] *= *CPglobalRot2; 

  fCableRadius = 0.5*mm;
  // fCableCenter[4] = {0.,0.,0.,0.};
  // fCableLength[4] = {2.54*10.0/2*cm, 2.54*8.0/2*cm, 2.54*5.0/2*cm, 2.54*2.5/2*cm};

  // Offsets are listed with signal before HV (signal is left of string)
  // Units are converted to cm
  fCableOffset[0] = G4ThreeVector(-2.54*1.839*cm, 2.54*0.560*cm, 0.); // P1
  fCableOffset[1] = G4ThreeVector(-2.54*1.839*cm, -2.54*0.560*cm, 0.);  
  fCableOffset[2] = G4ThreeVector(2.54*3.453*cm, 2.54*3.846*cm, 0.); // P2
  fCableOffset[3] = G4ThreeVector(2.54*2.350*cm, 2.54*3.620*cm, 0.);
  fCableOffset[4] = G4ThreeVector(-2.54*0.089*cm, 2.54*5.576*cm, 0.); // P3
  fCableOffset[5] = G4ThreeVector(-2.54*1.186*cm, 2.54*5.350*cm, 0.);
  fCableOffset[6] = G4ThreeVector(-2.54*3.356*cm, 2.54*3.846*cm, 0.); // P4
  fCableOffset[7] = G4ThreeVector(-2.54*4.387*cm, 2.54*3.408*cm, 0.);
  fCableOffset[8] = G4ThreeVector(-2.54*4.387*cm, -2.54*3.408*cm, 0.); // P5
  fCableOffset[9] = G4ThreeVector(-2.54*3.356*cm, -2.54*3.846*cm, 0.);
  fCableOffset[10] = G4ThreeVector(-2.54*1.186*cm, -2.54*5.350*cm, 0.); // P6
  fCableOffset[11] = G4ThreeVector(-2.54*0.089*cm, -2.54*5.576*cm, 0.);
  fCableOffset[12] = G4ThreeVector(2.54*2.350*cm, -2.54*3.620*cm, 0.); // P7
  fCableOffset[13] = G4ThreeVector(2.54*3.453*cm, -2.54*3.846*cm, 0.);

}

//---------------------------------------------------------------------------//

MGGeneratorMJDCable::MGGeneratorMJDCable(const MGGeneratorMJDCable & other) : MGVGenerator(other)
{;}

//---------------------------------------------------------------------------//

MGGeneratorMJDCable::~MGGeneratorMJDCable()
{
  delete fG4Messenger;
  delete fParticleGun;
}

//---------------------------------------------------------------------------//

void MGGeneratorMJDCable::BeginOfRunAction(G4Run const*)
{;}  

//---------------------------------------------------------------------------//

void MGGeneratorMJDCable::EndOfRunAction(G4Run const*)
{;}

//---------------------------------------------------------------------------//

void MGGeneratorMJDCable::GeneratePrimaryVertex(G4Event *event)
{

  // Generate random variables
  fRandString = G4RandFlat::shootInt(14);
  fRandPos = G4RandFlat::shootInt(4);
  fRandRadiusSq = fCableRadius*fCableRadius*G4UniformRand();
  fRandAngle = 2*pi*G4UniformRand();
  
  // Choose a random XY point along a disk
  fPositionX = sqrt( fRandRadiusSq ) * cos( fRandAngle );
  fPositionY = sqrt( fRandRadiusSq ) * sin( fRandAngle );
  // Choose a random cable and generate a random point along the cable
  fPositionZ = (1. - 2.*G4UniformRand())*fCableLength[fRandPos];

  // Set source position
  fPosition = fColdPlateOffset[0] + fCableOffset[fRandString] + G4ThreeVector(fPositionX, fPositionY, fPositionZ + fCableCenter[fRandPos]);

  // G4IonTable *theIonTable = (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
  // G4ParticleDefinition *aIon = theIonTable->GetIon(fZ, fA);

  // fParticleGun->SetParticleDefinition(aIon);
  // fParticleGun->SetParticleEnergy(0.0);
  fParticleGun->SetParticleMomentumDirection(G4RandomDirection()); // Geantinos for testing
  fParticleGun->SetParticleDefinition(G4Geantino::Geantino());
  fParticleGun->SetParticleEnergy(1.0*GeV); 
  fParticleGun->SetParticlePosition(fPosition);
  fParticleGun->GeneratePrimaryVertex(event);
}

//---------------------------------------------------------------------------//
