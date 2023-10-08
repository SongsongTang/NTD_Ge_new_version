//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file DBDecay/include/TrackingAction.hh
/// \brief Definition of the TrackingAction class
//
//
// $Id: TrackingAction.hh 78307 2013-12-11 10:55:57Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"
#include "HistoManager.hh"
#include "SteppingAction.hh"
#include <iomanip>
//class SteppingVerbose;
class RunAction;
class EventAction;
class TrackingMessenger;
class ParticleInfo;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackingAction : public G4UserTrackingAction {

  public:  
  TrackingAction(RunAction*,EventAction*,HistoManager*);
   ~TrackingAction();
   
    virtual void  PreUserTrackingAction(const G4Track*);
    virtual void  PostUserTrackingAction(const G4Track*);
    
    void SetFullChain(G4bool flag) { fFullChain = flag;};

    void SetVolumeFlag1(G4bool flag){ fInScoringVolume1 = flag;};
    G4bool GetVolumeFlag1(){return fInScoringVolume1;};
    void SetVolumeFlag2(G4bool flag){ fInScoringVolume2 = flag;};
    G4bool GetVolumeFlag2(){return fInScoringVolume2;};
    // G4double zmax;
    G4double MaxPosition[3];
    // G4double zmin;
    G4double MinPosition[3];
    G4double MaxEdep;
    G4double MaxEdepPos;
    G4double MaxEdepPosZ;
    std::vector<G4double> TrackStartPos;
    std::vector<G4double> EventStartPos;
    G4double E_primary;         //primary energy of the injected particle

  
  //ParticleInfo fParticleInfo_Track;
 /*
    void AddEdep_CuCollimation(G4double edep){
      fTrackEdepInCuCollimation +=edep;
      return;
    }
    void AddEdep_Hole1(G4double edep){
      fTrackEdepInHole1 +=edep;
      return;
    }
*/
    // void AddEdep_PET(G4double edep){
    //   fTrackEdepInPET +=edep;
    //   return;
    // }
    // void AddEdep_Source(G4double edep){
    //   fTrackEdepInSource +=edep;
    //   return;
    // }
    // void AddEdep_Gas(G4double edep){
    //   fTrackEdepInGas +=edep;
    //   return;
    // }
    void AddTracklen_ScoringVolume(G4double len){
      fTracklenInSV += len;
      return;
    }
    void AddTracklen_MMVolume(G4double len){
      fTracklenInMM += len;
      return;
    }
    void AddStepInfo(G4ThreeVector stepVertexPos, G4double dE_dx, G4double length){
      fTrackInfo_Stepping.fStepVertexPosX.push_back(stepVertexPos.x());
      fTrackInfo_Stepping.fStepVertexPosY.push_back(stepVertexPos.y());
      fTrackInfo_Stepping.fStepVertexPosZ.push_back(stepVertexPos.z());
      fTrackInfo_Stepping.fStepdE_dx.push_back(dE_dx);
      fTrackInfo_Stepping.fStepTrackLen.push_back(length);
    }
    void AddTrackEdep_SV(G4double edep){
      fTrackEnergyInSV += edep;
      return;
    }
    void AddTrackEdep_MM(G4double edep){
      fTrackEnergyInMM += edep;
      return;
    }
    G4double GetTrackLenInSV(){return fTracklenInSV;}
    G4double GetTrackLenInMM(){return fTracklenInMM;}
  


  private:
  RunAction* fRun;
  EventAction*        fEvent;
  HistoManager*       fHistoManager_Track;
  TrackingMessenger*  fTrackMessenger;
    
  G4double fCharge, fBaryonNo, fMass;
  G4double  fParticleEnCode;        
  //G4double fTrackEdepInCuCollimation;
  // G4double fTrackEdepInSource;
  //G4double fTrackEdepInHole1;
  // G4double fTrackEdepInPET;
  // G4double fTrackEdepInGas;
  G4double fTracklenInSV;
  G4double fTracklenInMM;
  G4double fTrackEnergyInSV;
  G4double fTrackEnergyInMM;
  G4bool   fFullChain;
  G4bool fInScoringVolume1;
  G4bool fInScoringVolume2;
  G4int nCounts;
  ParticleInfo fParticleInfo_Tracking;
  TrackInfo fTrackInfo_Stepping;
  G4String ParentTrackParticleName;
  G4int ParentID;
  G4String volume_name;
  G4String creator_process;
  G4double KinEnergy_start;
  //SteppingAction fSteppingAction_Tracking;
 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
