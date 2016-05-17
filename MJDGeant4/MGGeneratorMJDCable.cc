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

#include "generators/MGGeneratorMJDCable.hh"
#include "generators/MGGeneratorMJDCableMessenger.hh"
#include "io/MGLogger.hh"

//---------------------------------------------------------------------------//

#include "generators/MGGeneratorMJDCable.hh" 

//---------------------------------------------------------------------------//

using namespace CLHEP;

MGGeneratorMJDCable::MGGeneratorMJDCable()
{
  fGeneratorName = "MJDCable";
  fG4Messenger = new MGGeneratorMJDCableMessenger(this);
  fParticleGun = new G4ParticleGun(1);
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

void MGGeneratorMJDCable::SetCableDimensions()
{;}

//---------------------------------------------------------------------------//

void MGGeneratorMJDCable::GeneratePrimaryVertex(G4Event *event)
{
  G4IonTable *theIonTable = (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
  G4ParticleDefinition *aIon = theIonTable->GetIon(fZ, fA);

  fParticleGun->SetParticleDefinition(aIon);
  fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
  fParticleGun->SetParticleEnergy(0.0);

  fPosition =  fStringCenter[i] + fStringOffset[j] + (1. - 2.*G4UniformRand())*fCableLength;  

  fParticleGun->SetParticlePosition(fPosition);
  fParticleGun->GeneratePrimaryVertex(event);
}
  
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
