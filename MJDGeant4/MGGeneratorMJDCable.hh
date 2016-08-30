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
  void SetSourcePos(std::string sourcePos);
  void SetSourceType(std::string sourceType) {fSourceType = sourceType;}
  void SetIonZ(G4int z) {fZ = z;}
  void SetIonA(G4int a) {fA = a;}

  G4double GetIonZ() const {return fZ;}
  G4double GetIonA() const {return fA;}
  
  //protected members
protected:

  //private  members
private:
  G4ParticleGun *fParticleGun;

  // Particle types
  G4int fZ;
  G4int fA;

  std::string fSourceType;

  // Type of source? Should be radioactive probably, maybe just leave the A and Z settings?
  // To generate a line, find center point and half-length
  G4ThreeVector fPosition; // position of particle generated
  G4double fPositionX;
  G4double fPositionY;
  G4double fPositionZ; // Z position of particle generated
  G4double fZrotation;   // rotation about the z axis for the source


  G4double fRandRadiusSq;
  G4double fRandAngle;
  G4int fRandString; // Random integer to choose string position
  G4int fRandPos; // Random integer to choose cable position

  G4double fCableRadius; // Radius of a bundle of cables
  G4ThreeVector fCableOffset[7]; // XY location of cables wrt center of cold plate, even is signal and odd is HV
  G4ThreeVector fHVOffset[7]; // XY location of cables wrt center of cold plate, even is signal and odd is HV

  G4double fColdPlateRadius;
  G4double fColdPlateZ;    
  G4double fCrossArmLength;
  G4double fCrossArmWidth;
  G4double fCrossArmIRadius;
  G4double fCrossArmT;

  // 99 to -187 => Middle = -44 (simulation) = -120 (real) => Add 76
  // 287.0 mm total length currently => 280
  // Signal
  // 91+189, -126 center;  91+117, -89 center; 91+50, -56 center; 91+15, -38 => For the 4 detectors
  // 91+183, 91+128, 91+77, 91+24, 91-27 => For the 5 detectors, P
  // HV
  // 91+139, -100 center; 91+74, -67.5; 91+14, -37.5, 91-58, -1.5 => For the 4 detectors
  // 91+143, 91+91, 91+39, 91-5, 91-56 => For the 5 detectors
  G4double fCableLength[4] = {28.0/2*cm, 20.8/2*cm, 14.1/2*cm, 10.6/2*cm}; // Half length of signal cable, one for each detector 
  G4double fCableCenter[4] = {-12.6*cm, -8.9*cm, -5.6*cm, -3.8*cm}; // Centers of signal cables, one for each detector
  G4double fHVLength[4] = {23.0/2*cm, 16.5/2*cm, 10.5/2*cm, 3.3/2*cm}; // Half length of HV cable, one for each detector
  G4double fHVCenter[4] = {-10.0*cm, -6.75*cm, -3.75*cm, -0.15*cm}; // Centers of HV cables, one for each detector

  G4ThreeVector fColdPlateOffset[2]; // offset of cold plate to origin in world. 0 for Module 1, 1 for Module 2

};
#endif
