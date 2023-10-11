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
#include "Constant.h"

using namespace TPCsystem;

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

    //-----added on 2023.10.09, some variables for digitization-----
    bool fDigitization = false;   //whether digitization process is used
    double v_drift = 3.5;         //drift velocity in cm/us
    const int nch = 128;          //total channel numbers for x and y each
    const double chnwidth = 1.3;        //diagonal of one pad, unit is mm
    const double E_ion = 26.0;          //ionization energy of Ar, in eV
    const double timestep = 40.0;       //the timestep of one ADC sampling point
    std::vector<std::vector<double>> charge_X(nch);       //charge collected in each X channel, in MeV
    std::vector<std::vector<double>> time_X(nch);         //charge collection time in each X channel, in ns
    std::vector<std::vector<double>> charge_Y(nch);       //charge collected in each Y channel, in MeV
    std::vector<std::vector<double>> time_Y(nch);         //charge collection time in each Y channel, in ns
    //response function, bin width is 40ns
    const double response_func = {4.118, 7.618, -2.782, -1.982, 8.318, 7.818, 7.018,
        25.418, 23.418, -5.782, -3.882, 17.918, 12.718, 7.218, 11.718, 17.918, 22.418,
        18.918, 7.718, 17.118, 24.918, 18.918, 17.318, 32.518, 31.818, 10.918, 5.218,
        14.818, 10.818, 7.918, 15.118, 29.818, 128.618, 425.318, 941.018, 1568.718,
        2208.018, 2814.318, 3355.718, 3809.818, 4170.118, 4435.618, 4612.418, 4706.418,
        4718.318, 4672.318, 4582.418, 4445.818, 4267.718, 4070.018, 3840.418, 3587.518,
        3355.418, 3127.018, 2882.218, 2634.518, 2404.418, 2178.018, 1971.918, 1792.618,
        1604.918, 1423.118, 1271.618, 1132.818, 983.518, 850.718, 747.718, 649.018, 560.618,
        495.818, 431.318, 367.018, 315.018, 267.218, 224.018, 187.918, 156.018, 120.718,
        94.618, 86.818, 74.718, 57.418, 49.018, 39.418, 22.518, 16.818, 11.118, -5.482, -11.482,
        -8.082, -9.082, -17.082, -4.382, 1.518, -8.182, -18.782, -29.882, -24.982, -12.482,
        -11.982, -13.582, -10.382, -15.882, -14.882, 2.718, 9.418, -4.382, -11.482, -12.982,
        -23.582, -26.782, -23.382, -23.982, -21.882, -13.082, -18.582, -31.082, -33.782, -19.682,
        -11.782, -18.782, -32.482, -33.282, -25.882, -27.382, -27.482, -29.282, -29.082, -30.782,
        -34.382, -32.482, -28.382, -27.782, -21.982, -18.282, -27.382, -30.782, -21.082, -19.782,
        -28.582, -40.182, -37.282, -22.882, -30.182, -43.882, -42.182, -22.982, -28.882, -50.082, -42.882};
    //gain of the TPC system
    const double gain = 1;

    double waveform_X[Tch][Nsp]={0};
    double waveform_Y[Tch][Nsp]={0};

    int waveform_X_int[Tch][Nsp]={0};
    int waveform_Y_int[Tch][Nsp]={0};

    //--------------------------------------------------------------

  
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

  //----added on 2023.10.09, for digitization process----
  void ClearChannelBuffer(){
    charge_X.clear();
    time_X.clear();
    charge_Y.clear();
    time_Y.clear();
    memset(waveform_X,0,sizeof(waveform_X));
    memset(waveform_Y,0,sizeof(waveform_Y));
  }
  //calculate and save the waveforms of each channel
  void FillChannelWaveforms();
  //from the waveforms, get the data to save into the root file (should have the same format as the raw root file in testing)

  //------------------------------------------
 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
