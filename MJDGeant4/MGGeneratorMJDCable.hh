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
/**                                                            
 * $Id: MGGeneratorMJDCable.hh, $
 *      
 * CLASS DECLARATION:  MGGeneratorMJDCable.hh
 *
 *---------------------------------------------------------------------------//
 *
 * DESCRIPTION: 
 *
 */ 
// Begin description of class here
/**
 *
 * Generate Beam of Gammas as from the TUNL FEL
 * Currently the beam is simulated with a Gaussian energy and
 * width distribution in x and y.
 * The beam is created with the major axis along the x-axis and the minor
 * axis along the y axis. It is then rotated by fRho degrees on the x-y plane
 * fRho is relative to the x-axis. Particles are then generated parallel to 
 * the z-axis.
 *
 */
// End class description
//
/**  
 * SPECIAL NOTES:
 *
 */
// 
// --------------------------------------------------------------------------//
/** 
 * AUTHOR: B. Zhu
 * CONTACT: 
 * FIRST SUBMISSION:
 * 
 * REVISION:
 * 
 * 06-2016, Created, B. Zhu
 *
 */
// --------------------------------------------------------------------------//

#ifndef _MGGENERATORMJDCABLE_HH
#define _MGGENERATORMJDCABLE_HH

//---------------------------------------------------------------------------//

#include "globals.hh"
#include "G4ThreeVector.hh"

#include "generators/MGVGenerator.hh"

//---------------------------------------------------------------------------//

class G4Event;
class G4Messenger;
class G4ParticleGun;
class G4Run;

class MGGeneratorMJDCable : public MGVGenerator
{
public:

  //default constructor
  MGGeneratorMJDCable();

  //copy constructor
  MGGeneratorMJDCable(const MGGeneratorMJDCable &);

  //destructor
  ~MGGeneratorMJDCable();

  //public interface
  void BeginOfEventAction(G4Event *event);
  void BeginOfRunAction(G4Run const *run);
  void EndOfRunAction(G4Run const *run);
  void GeneratePrimaryVertex(G4Event *event);

  // G4ThreeVector GetCurrentPosition() { return fCurrentPosition; }
  
  //This method is not used but it is necessary 
  //because it is purely virtual in MGVGenerator
  void SetParticlePosition(G4ThreeVector) {;}

  // Sets dimensions of cables
  // void SetCableOffset();

  void SetIonZ(G4int z) {fZ = z;}
  void SetIonA(G4int a) {fA = a;}
  G4double GetIonZ() const {return fZ;}
  G4double GetIonA() const {return fA;}

  // Single particle
  // void SetParticle(G4int i) {fI = i};
  
  // Ion
  // void SetParticle(G4int z, G4int a) {fZ = z, fA = a;};
  
  //protected members
protected:

  //private  members
private:
  G4ParticleGun *fParticleGun;

  // Particle types
  G4int fI;
  G4int fZ;
  G4int fA;

  // Type of source? Should be radioactive probably, maybe just leave the A and Z settings?
  // To generate a line, find center point and half-length
  G4ThreeVector fPosition; // position of particle generated
  G4double fPositionX;
  G4double fPositionY;
  G4double fPositionZ; // Z position of particle generated

  G4double fRandRadiusSq;
  G4double fRandAngle;
  G4int fRandString; // Random integer to choose string position
  G4int fRandPos; // Random integer to choose cable position

  G4double fCableRadius; // Radius of a bundle of cables
  
  G4ThreeVector fStringCenter[14]; // center of strings, 7 strings per module
  G4ThreeVector fCableOffset[14]; // XY location of cables wrt center of cold plate, even is signal and odd is HV
  
  // G4double fStringOffset[14]; // offset of cable position from string center
  
  G4double fCableLength[4] = {2.54*10.0/2*cm, 2.54*8.0/2*cm, 2.54*5.0/2*cm, 2.54*2.5/2*cm}; // Half length of signal cable, one for each detector
  G4double fCableCenter[4] = {-2.54*10.0/2*cm, -2.54*8.0/2*cm, -2.54*5.0/2*cm, -2.54*2.5/2*cm}; // Centers of signal cables, one for each detector
  G4double fHVLength[4]; // Half length of HV cable, one for each detector
  G4double fHVCenter[4]; // Centers of HV cables, one for each detector
  
  G4ThreeVector fColdPlateOffset[2]; // offset of cold plate to origin in world. 0 for Module 1, 1 for Module 2

};
#endif
