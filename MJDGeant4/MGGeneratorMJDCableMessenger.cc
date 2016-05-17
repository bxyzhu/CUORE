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
// $Id: MGGeneratorMJDCableMessenger.cc,v 1.2 2007-02-21 09:31:33 mgmarino Exp $ 
//      
// CLASS IMPLEMENTATION:  MGGeneratorMJDCableMessenger.cc
//
//---------------------------------------------------------------------------//
/**
 * SPECIAL NOTES:
 * 
 *
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
 * 02-24-2016, Removed most commands from messenger since they are now
 *             read directly from the geometry, T. Caldwell
 * 07-14-2014, Created, K. Vorren
 */
//---------------------------------------------------------------------------//
//

#include "globals.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

#include "generators/MGGeneratorMJDCable.hh"

//---------------------------------------------------------------------------//

#include "generators/MGGeneratorMJDCableMessenger.hh"

//---------------------------------------------------------------------------//

MGGeneratorMJDCableMessenger::MGGeneratorMJDCableMessenger(MGGeneratorMJDCable *generator) : fMJDCableGenerator(generator)
{
  // /MG/generator/MJDCable
  fMJDCableDirectory = new G4UIdirectory("/MG/generator/MJDCable/");
  fMJDCableDirectory->SetGuidance("MJD Cable generator");

  // MG/generator/MJDCable/dump
  // fDumpCmd = new G4UIcmdWithoutParameter("/MG/generator/MJDCable/dump", this);
  // fDumpCmd->SetGuidance("Dump parameters of generator to screen");

  // /MG/generator/MJDCable/setZ
  fZCmd = new G4UIcmdWithAnInteger("/MG/generator/MJDCable/setZ", this);
  fZCmd->SetGuidance("The Z value of the ion");
  generator->SetIonZ(90);
  fZCmd->SetDefaultValue(generator->GetIonZ());
  
  // /MG/generator/MJDCable/setA
  fACmd = new G4UIcmdWithAnInteger("/MG/generator/MJDCable/setA", this);
  fACmd->SetGuidance("The A value of the ion");
  generator->SetIonA(228);
  fACmd->SetDefaultValue(generator->GetIonZ());
  
}

//---------------------------------------------------------------------------//

MGGeneratorMJDCableMessenger::MGGeneratorMJDCableMessenger(const MGGeneratorMJDCableMessenger & other) : G4UImessenger(other)
{;}

//---------------------------------------------------------------------------//

MGGeneratorMJDCableMessenger::~MGGeneratorMJDCableMessenger()
{
  delete fACmd;
  delete fZCmd;
  delete fDumpCmd;
  delete fMJDCableDirectory;
}

//---------------------------------------------------------------------------//

G4String MGGeneratorMJDCableMessenger::GetCurrentValue(G4UIcommand *cmd)
{
  G4String cv = " ";
  if(cmd == fDumpCmd)
    fMJDCableGenerator->Dump();

  return cv;
}

//---------------------------------------------------------------------------//

void MGGeneratorMJDCableMessenger::SetNewValue(G4UIcommand *cmd, G4String newValues)
{
  if(cmd == fZCmd)
    fMJDCableGenerator->SetIonZ(fZCmd->GetNewIntValue(newValues));
  else if(cmd == fACmd)
    fMJDCableGenerator->SetIonA(fACmd->GetNewIntValue(newValues));
  else if(cmd == fSourcePosCmd)
    fMJDCableGenerator->SetSourcePos(newValues);
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
