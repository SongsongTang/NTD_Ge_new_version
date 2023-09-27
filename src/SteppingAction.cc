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
/// \file analysis/shared/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
// $Id: SteppingAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "TrackingAction.hh"

#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//
SteppingAction* SteppingAction::fgInstance =0;
SteppingAction* SteppingAction::Instance()
{
  // G4cout<<"<<------------SteppingAction::Instance()-------------------->>"<<G4endl;
  return fgInstance;
}
//
SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* event, TrackingAction* tracking) 
  : G4UserSteppingAction(),
    fDetector(det),
    fEventAction(event), 
    fTrackingAction(tracking)                                         
{
  // G4cout<<"<<------------SteppingAction::SteppingAction(DetectorConstruction* det,TrackingAction* tracking) -------------------->>"<<G4endl;
  fgInstance = this;
  fEdep = 0;
  fStepLen = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ 
//  G4cout<<"<<------------SteppingAction::~SteppingAction()-------------------->>"<<G4endl;
  fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // G4cout<<"<<------------SteppingAction::UserSteppingAction(const G4Step* aStep)-------------------->>"<<G4endl;
  // get volume of the current step
  G4String preVolumeName 
    = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  // if(preVolumeName != "collimation" && preVolumeName != "source" && preVolumeName != "Hole1" && preVolumeName != "PET" && preVolumeName != "Gas")
  //   return;
  
  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double steplen = aStep->GetStepLength();



  /*
  if(edep>0 && preVolumeName == "collimation"){
    fTrackingAction->AddEdep_CuCollimation(edep);
  }
  
	
	if(edep>0 && preVolumeName == "Hole1"){
    fTrackingAction->AddEdep_Hole1(edep);
  }
*/

	// if(edep>0 && preVolumeName == "PET"){
  //   fTrackingAction->AddEdep_PET(edep);
  // }
	// if(edep>0 && preVolumeName == "source"){
  //   fTrackingAction->AddEdep_Source(edep);
  // } 

	// if(edep>0 && preVolumeName == "Gas"){
  //   fTrackingAction->AddEdep_Gas(edep);
  // }

  //add this steps' energy lost if it is in the Effective Gas volume
  if(preVolumeName == "GasEff"){
    fEventAction->AddEdep_ScoringVolume(edep);
  }

  //NEW PART: if this step does NOT begin in the scoring volume, drops it
  if(fTrackingAction->GetVolumeFlag1()){ 
    //if the previous tracks are ALL NOT in the periphery volume and the anticoincident MM region
    if(preVolumeName == "Gas" || preVolumeName == "GasEff2")  fTrackingAction->SetVolumeFlag1(false);   //this means the track get into the frame
    if(preVolumeName == "GasEff")  {
      fTrackingAction->SetVolumeFlag2(true);        //this means this track get into the gas volume

      // //Add secondary particles' energy
      // G4TrackVector* trackvector = astep->GetSecondary();
      // for(int i=0;i<trackvector.size();i++){
      //   fTrackingAction->AddEdep_ScoringVolume(trackvector[i]->GetKineticEnergy());
      // }

      

      //record the step points of this track
      G4ThreeVector stepVertexPos = aStep->GetPreStepPoint()->GetPosition();
      //calculate the total energy transfer in this step, including the energy transferred to secondaries
      G4double TotalEneTransfer = aStep->GetPreStepPoint()->GetKineticEnergy() - aStep->GetPostStepPoint()->GetKineticEnergy();
      
      fTrackingAction->AddTracklen_ScoringVolume(steplen);
      fTrackingAction->AddTrackEdep_SV(TotalEneTransfer);

      //this part sum the edep until the length is long enough, then fill
      fEdep += TotalEneTransfer;
      fStepLen += steplen;
      if(fStepLen > 0.5*mm || aStep->GetPostStepPoint()->GetKineticEnergy()==0){
        G4double dE_dx = fEdep/fStepLen;
        fEdep = 0;        //reset
        fStepLen = 0;     //reset
        //Get the total track length until this step
        // G4double TotalTrackLen = aStep->GetTrack()->GetTrackLength();

        if(stepVertexPos.x()>fTrackingAction->MaxPosition[0]) fTrackingAction->MaxPosition[0] = stepVertexPos.x();
        if(stepVertexPos.x()<fTrackingAction->MinPosition[0]) fTrackingAction->MinPosition[0] = stepVertexPos.x();
        if(stepVertexPos.y()>fTrackingAction->MaxPosition[1]) fTrackingAction->MaxPosition[1] = stepVertexPos.y();
        if(stepVertexPos.y()<fTrackingAction->MinPosition[1]) fTrackingAction->MinPosition[1] = stepVertexPos.y();
        if(stepVertexPos.z()>fTrackingAction->MaxPosition[2]) fTrackingAction->MaxPosition[2] = stepVertexPos.z();
        if(stepVertexPos.z()<fTrackingAction->MinPosition[2]) fTrackingAction->MinPosition[2] = stepVertexPos.z();

        if(dE_dx>(fTrackingAction->MaxEdep)){ 
          fTrackingAction->MaxEdep = dE_dx;
          fTrackingAction->MaxEdepPos = fTrackingAction->GetTrackLenInSV();
          fTrackingAction->MaxEdepPosZ = stepVertexPos.z();
        }

        //record step informations: the step points of this track, the tracklen & dE/dx, for plotting
        fTrackingAction->AddStepInfo(stepVertexPos, dE_dx, fTrackingAction->GetTrackLenInSV());
      }
    }
  }


}
 
void SteppingAction::Reset()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
