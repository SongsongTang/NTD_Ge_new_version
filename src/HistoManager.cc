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
//
// $Id: HistoManager.cc 72249 2013-07-12 08:57:49Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include <TTree.h>
#include <TFile.h>
#include <TPolyLine3D.h>
#include <TGraph.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("NTD_Ge"),fRootFile(0),fNtuple(0)
{
  G4cout<<"<<------------HistoManager::HistoManager()-------------------->>"<<G4endl;
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  G4cout<<"<<------------HistoManager::~HistoManager()-------------------->>"<<G4endl;
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  G4cout<<"<<------------HistoManager::Book()-------------------->>"<<G4endl;
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);     //enable inactivation of histograms
  
  // Define histograms start values
  // const G4int kMaxHisto = 16;
  const G4String id[] = {"0",
                        "1",
                        "2",
                        "3",
                        "4",
                        "5",
                        "6",
                        "7",
                        "8",
                        "energy_spectrum",
                        "track_length",
                        "primary_spectrum",
                        "secondary_particle_spectrum",
                        "trackstartpos_x_selected",
                        "trackstartpos_y_selected",
                        "primary_spectrum_selected",
                        "primary_spectrum_selected_PCB",
                        "primary_spectrum_selected_Gas",
                        "primary_spectrum_selected_film",
                        "primary_spectrum_selected_Cubrd",
                        "primary_spectrum_selected_other",
                        "primary_spectrum_selected_compt",
                        "primary_spectrum_selected_phont",
                        "primary_spectrum_selected_eIoni",
                        "primary_spectrum_selected_other_process"};
  const G4String title[] = 
          { "dummy",                                    //0
            "energy spectrum (%): e+ e-",               //1
            "energy spectrum (%): nu_e anti_nu_e",      //2
            "energy spectrum (%): gamma",               //3                  
            "energy spectrum (%): alpha",               //4
            "energy spectrum (%): ions",                //5
            "total kinetic energy (Q)",                 //6                        
            "momentum balance",                         //7
            "total time of life of decay chain",         //8
            "energy spectrum of interested particle in Gas volume",     //9
            "track length (cm) of interested particle in Gas volume",      //10
            "primary particle spectrum",                   //11
            "secondary_particle_spectrum",                //12
            "track start position in x dimension",        //13
            "track start position in y dimension",         //14
            "primary particle energy after selection",       //15
            "primary particle energy after selection from PCB",       //16
            "primary particle energy after selection from Gas",       //17
            "primary particle energy after selection from film",      //18
            "primary particle energy after selection from Cu board",   //19
            "primary particle energy after selection from other places",     //20
            "primary particle energy after selection from compton process",    //21
            "primary particle energy after selection from photon electric process",   //22
            "primary particle energy after selection from eIoni process",           //23
            "primary particle energy after selection from other processes"        //24
          };  

    const G4int kMaxHisto = sizeof(id)/sizeof(id[0]);

    if(kMaxHisto!=sizeof(title)/sizeof(title[0])) {
      std::cout << "Error: the length of the id and the title is not compatible!!" << std::endl;
      exit(0);
    }

  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }
  
}

void HistoManager::book()
{
  G4cout<<"<<------------HistoManager::book()-------------------->>"<<G4endl;
  G4cout<<"------>create Particle info rootfile"<<G4endl;
  
  fFileName = G4AnalysisManager::Instance()->GetFileName();

  // to satisfy the condition that the file is in a sub directory
  size_t pos = fFileName.find("/");
  G4String datafile_name, trackfile_name;
  if(pos!=std::string::npos){
    G4String firstpart  = fFileName.substr(0,pos);
    G4String secondpart = fFileName.substr(pos+1);
    datafile_name = firstpart+"/tree_"+secondpart+".root";
    trackfile_name = firstpart+"/tracks_"+secondpart+".root";
  }
  else{
    datafile_name = "tree_"+fFileName+".root";
    trackfile_name = "tracks_"+fFileName+".root";
  }
  
  fRootFile = new TFile(datafile_name,"RECREATE");
  fNtuple = new TTree("T","data of decay");
  fNtuple->Branch("nTrack",                 &fParticleInfo.nTrack,                  "nTrack/I");
  // fNtuple->Branch("ParticleZ",              &fParticleInfo.fParticleZ);
  // fNtuple->Branch("ParticleA",              &fParticleInfo.fParticleA);
  // fNtuple->Branch("ParticleMass",           &fParticleInfo.fParticleMass);
  // fNtuple->Branch("ParticleCode",           &fParticleInfo.fParticleCode);
  // fNtuple->Branch("ParticleKineticEnergy",  &fParticleInfo.fTrackKineticEnergy);
  // fNtuple->Branch("TrackVertexPosX",        &fParticleInfo.fTrackVertexPosX);
  // fNtuple->Branch("TrackVertexPosY",        &fParticleInfo.fTrackVertexPosY);
  // fNtuple->Branch("TrackVertexPosZ",        &fParticleInfo.fTrackVertexPosZ);
  // fNtuple->Branch("TrackVertexDirX",        &fParticleInfo.fTrackVertexDirX);
  // fNtuple->Branch("TrackVertexDirY",        &fParticleInfo.fTrackVertexDirY);
  // fNtuple->Branch("TrackVertexDirZ",        &fParticleInfo.fTrackVertexDirZ);
  // fNtuple->Branch("TrackID",                &fParticleInfo.fTrackID);
  // fNtuple->Branch("TrackParentID",          &fParticleInfo.fParentID);
  // // fNtuple->Branch("TrackEdepInSource",    &fParticleInfo.fTrackEdepInSource);
  // //fNtuple->Branch("TrackEdepInCuCollimation",     &fParticleInfo.fTrackEdepInCuCollimation);
  // //fNtuple->Branch("TrackEdepInHole1",    &fParticleInfo.fTrackEdepInHole1);
  // // fNtuple->Branch("TrackEdepInPET",    &fParticleInfo.fTrackEdepInPET);
  // // fNtuple->Branch("TrackEdepInGas",    &fParticleInfo.fTrackEdepInGas);
  // fNtuple->Branch("TrackTime",              &fParticleInfo.fTrackTime);

  //new branches
  fNtuple->Branch("StartPosX",              &fParticleInfo.fStartPosX);
  fNtuple->Branch("StartPosY",              &fParticleInfo.fStartPosY);
  fNtuple->Branch("TrackLength",              &fParticleInfo.fTrackLength);
  fNtuple->Branch("XTrackLength",              &fParticleInfo.fXTrackLength);
  fNtuple->Branch("YTrackLength",              &fParticleInfo.fYTrackLength);
  fNtuple->Branch("DriftDistance",              &fParticleInfo.fDriftDistance);
  fNtuple->Branch("CotTheta",              &fParticleInfo.fCotTheta);
  fNtuple->Branch("TrackEnergy",              &fParticleInfo.fTrackEnergy);
  fNtuple->Branch("MaxEdepPosition",              &fParticleInfo.fMaxEdepPosition);
  fNtuple->Branch("MaxEdepPositionZ",             &fParticleInfo.fMaxEdepPositionZ);

  fNtuple->Branch("TrackStartPos",      &fParticleInfo.fTrackStartPos);
  fNtuple->Branch("EventStartPos",      &fParticleInfo.fEventStartPos);
  fNtuple->Branch("TrackVertexPosX_per_step",      &fParticleInfo.fTrackVertexPosX_per_step);
  fNtuple->Branch("TrackVertexPosY_per_step",      &fParticleInfo.fTrackVertexPosY_per_step);
  fNtuple->Branch("TrackVertexPosZ_per_step",      &fParticleInfo.fTrackVertexPosZ_per_step);
  fNtuple->Branch("TrackdE_dx_per_step",           &fParticleInfo.fTrackdE_dx_per_step);
  fNtuple->Branch("TrackLength_per_step",          &fParticleInfo.fTrackLength_per_step);

  //=========================================================
  G4cout<<"------>create track graph rootfile"<<G4endl;

  // G4String trackfile_name = "tracks_"+fFileName+".root";
  G4cout<<"The output track file name is : "<<trackfile_name<<"\n"<< G4endl;
  fGraphRootFile = new TFile(trackfile_name,"RECREATE");

}

void HistoManager::save()
{
  G4cout<<"<<------------HistoManager::save()-------------------->>"<<G4endl;
  fRootFile->cd();
  fNtuple->Write();
 // fRootFile->Write();
  fRootFile->Close();
  fGraphRootFile->Close();
  G4cout<<"------>close rootfiles"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillTrackGraph(TrackInfo* fTrackInfo, G4int nEvents, G4int nTracks){
    
  // if(nEvents>100)  return;
  G4double x, y, z;
  G4double dE_dx, tracklen;
  G4int npoints = fTrackInfo->fStepVertexPosX.size();
  TPolyLine3D* fTrackGraph = new TPolyLine3D(npoints);
  TGraph* Track_XZ = new TGraph();
  TGraph* Track_YZ = new TGraph();
  TGraph* dE_dx_graph = new TGraph();
  char grtitle[100];

  Track_XZ->SetMarkerStyle(7);
  memset(grtitle, 0, sizeof(grtitle));
  sprintf(grtitle, "Event%d_Track%d_XZ;X(mm);Z(mm)", nEvents, nTracks);
  Track_XZ->SetTitle(grtitle);
  Track_XZ->SetName(grtitle);

  Track_YZ->SetMarkerStyle(7);
  memset(grtitle, 0, sizeof(grtitle));
  sprintf(grtitle, "Event%d_Track%d_YZ;Y(mm);Z(mm)", nEvents, nTracks);
  Track_YZ->SetTitle(grtitle);
  Track_YZ->SetName(grtitle);

  dE_dx_graph->SetMarkerStyle(7);
  memset(grtitle, 0, sizeof(grtitle));
  sprintf(grtitle, "Event%d_Track%d_dE/dx;track length(mm);dE/dx(keV/mm)", nEvents, nTracks);
  dE_dx_graph->SetTitle(grtitle);
  dE_dx_graph->SetName(grtitle);

  // G4cout << "-----X   " << "Y   " << "Z   "<< G4endl;

  for(G4int i = 0; i<npoints; i++){
    x = fTrackInfo->fStepVertexPosX[i];
    y = fTrackInfo->fStepVertexPosY[i];
    z = fTrackInfo->fStepVertexPosZ[i];
    dE_dx = fTrackInfo->fStepdE_dx[i];
    tracklen = fTrackInfo->fStepTrackLen[i];

    // G4cout << "     " << x <<"   " << y << "   " << z << G4endl;

    fTrackGraph->SetPoint(i,x,y,z);
    Track_XZ->SetPoint(i,x,z);
    Track_YZ->SetPoint(i,y,z);
    dE_dx_graph->SetPoint(i,tracklen,dE_dx*1000);     //MeV/mm---->keV/mm
  }

  fGraphRootFile->cd();
  // fTrackGraph->Write();
  Track_XZ->Write();
  Track_YZ->Write();
  dE_dx_graph->Write();

}
