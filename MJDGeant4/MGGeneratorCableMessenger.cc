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
// $Id: MGGeneratorCableMessenger.cc,v 1.2 2007-02-21 09:31:33 mgmarino Exp $ 
//      
// CLASS IMPLEMENTATION:  MGGeneratorCableMessenger.cc
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
 * AUTHOR: B. Zhu
 * CONTACT: 
 * FIRST SUBMISSION:
 * 
 * REVISION:
 *
 * 06-2016, Created, B. Zhu
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

#include "generators/MGGeneratorCableMessenger.hh"

//---------------------------------------------------------------------------//

MGGeneratorCableMessenger::MGGeneratorCableMessenger(MGGeneratorMJDCable *generator) : fMJDCableGenerator(generator)
{
  // /MG/generator/Cable
  fCableDirectory = new G4UIdirectory("/MG/generator/Cable/");
  fCableDirectory->SetGuidance("MJD Cable generator");

  // /MG/generator/Cable/setZ
  fZCmd = new G4UIcmdWithAnInteger("/MG/generator/Cable/setZ", this);
  fZCmd->SetGuidance("The Z value of the ion");
  generator->SetIonZ(90);
  fZCmd->SetDefaultValue(generator->GetIonZ());
  
  // /MG/generator/Cable/setA
  fACmd = new G4UIcmdWithAnInteger("/MG/generator/Cable/setA", this);
  fACmd->SetGuidance("The A value of the ion");
  generator->SetIonA(228);
  fACmd->SetDefaultValue(generator->GetIonZ());
  
  // /MG/generator/Cable/setSourceType
  fSourceTypeCmd =
    new G4UIcmdWithAString("/MG/generator/Cable/setSourceType", this);
  fSourceTypeCmd->SetGuidance("Set source type: H(igh voltage), S(ignal), (Cold) P(late), or C(ross arm)");
  fSourceTypeCmd->SetCandidates("H S C P");
  fSourceTypeCmd->SetDefaultValue("S");
  fCableGenerator->SetSourceType("S");

  // /MG/generator/Cable/setSourcePos
  fSourcePosCmd =
    new G4UIcmdWithAString("/MG/generator/Cable/setSourcePos", this);
  fSourcePosCmd->SetGuidance("West or East cryostat position.");
  fSourcePosCmd->SetCandidates("W E");
  fSourcePosCmd->SetDefaultValue("W");
  fMJDCableGenerator->SetSourcePos("W");

}

//---------------------------------------------------------------------------//

MGGeneratorCableMessenger::MGGeneratorCableMessenger(const MGGeneratorCableMessenger & other) : G4UImessenger(other)
{;}

//---------------------------------------------------------------------------//

MGGeneratorCableMessenger::~MGGeneratorCableMessenger()
{
  delete fACmd;
  delete fZCmd;
  delete fMJDCableDirectory;
}

//---------------------------------------------------------------------------//

void MGGeneratorCableMessenger::SetNewValue(G4UIcommand *cmd, G4String newValues)
{
  if(cmd == fZCmd)
    fMJDCableGenerator->SetIonZ(fZCmd->GetNewIntValue(newValues));
  else if(cmd == fACmd)
    fMJDCableGenerator->SetIonA(fACmd->GetNewIntValue(newValues));
  else if(cmd == fSourcePosCmd)
    fMJDCableGenerator->SetSourcePos(newValues);
  else if(cmd == fSourceTypeCmd)
    fMJDCableGenerator->SetSourceType(newValues);
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
