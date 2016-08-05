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
 * AUTHOR: B. Zhu
 * CONTACT: 
 * FIRST SUBMISSION:
 * 
 * REVISION:
 *
 * 07-2016, Created, B. Zhu
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

  // No implementation for "E" yet
  std::string sourcePos = "W";
  G4UIcommandTree* cmdTree = G4UImanager::GetUIpointer()->GetTree()->GetTree("/MG/");
  cmdTree = cmdTree->GetTree(G4String("/MG/demonstrator/"));
  for(int i=0; i<cmdTree->GetCommandEntry(); i++)
  {
    std::string param = cmdTree->GetCommand(i+1)->GetCommandName();
    if(param == "cryo1Pos" || param == "cryo2Pos")
    {
      G4UIcommand* cmd = cmdTree->GetCommand(i+1);
      G4double x, y, z;
      std::stringstream(cmd->GetParameter(0)->GetDefaultValue()) >> x;
      std::stringstream(cmd->GetParameter(1)->GetDefaultValue()) >> y;
      std::stringstream(cmd->GetParameter(2)->GetDefaultValue()) >> z;
      if(param == "cryo1Pos" && sourcePos == "W")
  fColdPlateOffset[0] = G4ThreeVector(x, y, z);
      else if(param == "cryo2Pos" && sourcePos == "E")
  fColdPlateOffset[0] = G4ThreeVector(x, y, z);
    }
    else if(param == "cryo1Rot" && sourcePos == "W")
    {
      std::string val =
  cmdTree->GetCommand(i+1)->GetParameter(0)->GetDefaultValue();
      // std::stringstream(val) >> fZrotation;
    }
    else if(param == "cryo2Rot" && sourcePos == "E")
    {
      std::string val =
  cmdTree->GetCommand(i+1)->GetParameter(0)->GetDefaultValue();
      // std::stringstream(val) >> fZrotation;
    }
  }
  cmdTree = G4UImanager::GetUIpointer()->GetTree()->GetTree("/MG/");
  cmdTree = cmdTree->GetTree(G4String("/MG/mjdemocryoassembly"+sourcePos+"/"));
  for(int i=0; i<cmdTree->GetCommandEntry(); i++)
  {
    std::string param = cmdTree->GetCommand(i+1)->GetCommandName();
    if(param == "calAssemblyPos"){
      G4UIcommand* cmd = cmdTree->GetCommand(i+1);
      G4double x, y, z;
      std::stringstream(cmd->GetParameter(0)->GetDefaultValue()) >> x;
      std::stringstream(cmd->GetParameter(1)->GetDefaultValue()) >> y;
      std::stringstream(cmd->GetParameter(2)->GetDefaultValue()) >> z;
      fColdPlateOffset[0] += G4ThreeVector(x, y, z/3); // use 0 here
    }
  }
  cmdTree = G4UImanager::GetUIpointer()->GetTree()->GetTree("/MG/");
  cmdTree = cmdTree->GetTree(G4String("/MG/mjdemocalassembly"+sourcePos+"/"));
  for(int i=0; i<cmdTree->GetCommandEntry(); i++)
  {
    std::string param = cmdTree->GetCommand(i+1)->GetCommandName();
    if(param == "sourceOffset")
    {
      G4UIcommand* cmd = cmdTree->GetCommand(i+1);
      G4double x, y, z;
      std::stringstream(cmd->GetParameter(0)->GetDefaultValue()) >> x;
      std::stringstream(cmd->GetParameter(1)->GetDefaultValue()) >> y;
      std::stringstream(cmd->GetParameter(2)->GetDefaultValue()) >> z;
      fColdPlateOffset[0] += G4ThreeVector(x, y, z);
    }
  }

  fCableRadius = 0.05*mm;

  // Offsets are listed with signal before HV (signal is left of string)
  // Units were originally in inches and then converted to cm
  // The drawing and simulation geometries aren't one-to-one so I made some slight adjustments
  fCableOffset[0] = G4ThreeVector(-2.54*1.839*cm, 2.54*0.560*cm, 0.); // P1
  fCableOffset[1] = G4ThreeVector(-2.54*1.839*cm, -2.54*0.560*cm, 0.);  

  fCableOffset[2] = G4ThreeVector(2.54*3.425*cm, 2.54*3.75*cm, 0.); // P2
  fCableOffset[3] = G4ThreeVector(2.54*2.350*cm, 2.54*3.6*cm, 0.);

  fCableOffset[4] = G4ThreeVector(-2.54*0.090*cm, 2.54*5.59*cm, 0.); // P3
  fCableOffset[5] = G4ThreeVector(-2.54*1.15*cm, 2.54*5.25*cm, 0.);

  fCableOffset[6] = G4ThreeVector(-2.54*3.356*cm, 2.54*3.81*cm, 0.); // P4
  fCableOffset[7] = G4ThreeVector(-2.54*4.387*cm, 2.54*3.408*cm, 0.);
  
  fCableOffset[8] = G4ThreeVector(-2.54*4.387*cm, -2.54*3.408*cm, 0.); // P5
  fCableOffset[9] = G4ThreeVector(-2.54*3.356*cm, -2.54*3.81*cm, 0.);

  fCableOffset[10] = G4ThreeVector(-2.54*1.15 *cm, -2.54*5.25*cm, 0.); // P6
  fCableOffset[11] = G4ThreeVector(-2.54*0.090*cm, -2.54*5.59*cm, 0.);
  
  fCableOffset[12] = G4ThreeVector(2.54*2.350*cm, -2.54*3.6*cm, 0.); // P7
  fCableOffset[13] = G4ThreeVector(2.54*3.425*cm, -2.54*3.75*cm, 0.);

  // Haven't fixed for E yet
  if(sourcePos == "W")
  {
    for(int i = 0; i < 14; i++) fCableOffset[i].rotateZ(-pi/2);
  }

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

  G4IonTable *theIonTable = (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
  G4ParticleDefinition *aIon = theIonTable->GetIon(fZ, fA);

  fParticleGun->SetParticleDefinition(aIon);
  fParticleGun->SetParticleEnergy(0.0);
  fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
  fParticleGun->SetParticlePosition(fPosition);
  fParticleGun->GeneratePrimaryVertex(event);
}

//---------------------------------------------------------------------------//

