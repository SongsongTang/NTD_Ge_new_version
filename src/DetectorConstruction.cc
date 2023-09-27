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
/// \file radioactivedecay/rdecay01/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
// $Id: DetectorConstruction.cc 78307 2013-12-11 10:55:57Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "SteppingAction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4UserLimits.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4AutoDelete.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0;

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{G4cout<<"<<------------DetectorConstruction::DetectorConstruction()-------------------->>"<<G4endl;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{G4cout<<"<<------------DetectorConstruction::~DetectorConstruction()-------------------->>"<<G4endl;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	G4cout<<"<<------------DetectorConstruction::Construct()-------------------->>"<<G4endl;
	//
	// World volume
	//   

	// Material ---> Vacuum  
	G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
	G4Material* Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");  
	G4Material* Ge = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ge"); 
	G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
	//Material --->CO2
	G4Material* CO2 = G4NistManager::Instance()->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
	//Material --->Ar
	G4Material* Ar = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ar");
	//Material --->Ne
	G4Material* Ne = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ne");
	G4Material* Xe = G4NistManager::Instance()->FindOrBuildMaterial("G4_Xe");

	//Materials defined by self
	G4String symbol;
	G4double a;                                              // atomic mass
	G4double z;                                              // atomic number
	G4double density;
	G4int ncomponents, natoms;
	G4Element* C  = new G4Element("Carbon",     symbol= "C",  z= 6.,  a= 12.00*g/mole);
	G4Element* H  = new G4Element("Hydrogen",   symbol= "H",  z= 1.,  a= 1.00*g/mole);
	G4Element* O  = new G4Element("Oxygen",     symbol= "O" , z= 8. , a= 16.00*g/mole);
	G4Element* N  = new G4Element("Nitrogen",     symbol= "N" , z= 7. , a= 14.00*g/mole);
	G4Element* Si = new G4Element("Si",			symbol= "Si", z= 14., a= 28.00*g/mole);
	//Material--->PET C10H8O4
	G4Material* PET = new G4Material("PET", density=1.38*g/cm3, ncomponents=3); //PET 200um
	PET->AddElement(C, natoms=10);
	PET->AddElement(H, natoms=8);
	PET->AddElement(O, natoms=4);

	//Define FR4, from https://agenda.infn.it/event/14179/contributions/23405/attachments/16712/18966/Geant_FR4_Raffaella_v1.pdf
	//epoxy
	G4Material* Epoxy = new G4Material("Epoxy", density = 1.2*g/cm3, ncomponents = 2);
	Epoxy->AddElement(H, natoms = 2);
	Epoxy->AddElement(C, natoms = 2);
	//SiO2
	G4Material* SiO2 = new G4Material("SiO2", density = 2.200*g/cm3, ncomponents = 2);
	SiO2->AddElement(Si, natoms = 1);
	SiO2->AddElement(O, natoms = 2);
	//Material--->FR4
	G4Material* FR4 = new G4Material("FR4", density = 1.86*g/cm3, ncomponents = 2);
	FR4->AddMaterial(Epoxy, 47.2*perCent);
	FR4->AddMaterial(SiO2, 52.8*perCent);

	G4Material* Polyacrylate = new G4Material("Polyacrylate", density=1.39*g/cm3, ncomponents=3); //PET 200um
	Polyacrylate->AddElement(C, natoms=3);
	Polyacrylate->AddElement(H, natoms=4);
	Polyacrylate->AddElement(O, natoms=2);
	//Material ---> Atlasgas
	G4Material* Atlasgas = new G4Material("Atlasgas", density=1.79e-3*g/cm3, ncomponents=2);   //atlas gas 10mm
	Atlasgas->AddMaterial(CO2, 7.0*perCent);
	Atlasgas->AddMaterial(Ar,  93.0*perCent);

	//Material ---> C4H10
	G4Material* iC4H10 = new G4Material("iC4H10",  density=2.487e-3*g/cm3, ncomponents=2);  
	//density should be checked!!!
	iC4H10->AddElement(C, natoms=4);
	iC4H10->AddElement(H, natoms=10);

	//Material ---> Ne+C4H10
	G4Material* NeiC4H10 = new G4Material("NeiC4H10", density=0.98e-3*g/cm3, ncomponents=2);   
	NeiC4H10->AddMaterial(iC4H10,  5.0*perCent);
	NeiC4H10->AddMaterial(Ne,     95.0*perCent);

	//Material ---> Ar+iC4H10				//density is calculated by volume fraction 96.5/3.5
	G4Material* AriC4H10 = new G4Material("AriC4H10", density = (0.001782 * 0.965 + 0.00251 * 0.035)*g/cm3, ncomponents = 2);
	AriC4H10->AddMaterial(iC4H10, 	3.5*perCent);
	AriC4H10->AddMaterial(Ar,	  	96.5*perCent);

	//Material ---> Kapton
	G4Material* Kapton = new G4Material("Kapton", density = 1.43*g/cm3, ncomponents = 4);
	Kapton->AddElement(C, natoms = 22);
	Kapton->AddElement(H, natoms = 10);
	Kapton->AddElement(N, natoms = 2);
	Kapton->AddElement(O, natoms = 5);

	//---------------construct detector---------------
	// Full sphere shape
	G4double solidWorld_rmax = 100*cm;
	G4Orb* solidWorld = new G4Orb(
			"World",                   // its name
			solidWorld_rmax);                 // its size 

	G4LogicalVolume* logicWorld = new G4LogicalVolume(
			solidWorld,             // its solid
			Air,                 // its material
			"World");               // its name
	G4VPhysicalVolume* physicalWorld = new G4PVPlacement(
			0,                        // no rotation
			G4ThreeVector(),          // at (0,0,0)
			logicWorld,               // its logical volume
			"World",                  // its name
			0,                        // its mother  volume
			false,                    // no boolean operation
			0);                       // copy number


	//common parameters of the aluminum frame
	G4double FrameSizeX = 20.6*cm;
	G4double FrameSizeY = 20.6*cm;
	G4double InnerFrameSizeX = 18.2*cm;
	G4double InnerFrameSizeY = 18.2*cm;
	G4double DetectorSizeX = 15*cm;
	G4double DetectorSizeY = 15*cm;
	G4double FrameStepHeight = 0.05*cm;
	G4double WindowSizeX = 16.6*cm;
	G4double WindowSizeY = 16.6*cm;

	// parameters of PET
	G4double petsizeX = InnerFrameSizeX;
	G4double petsizeY = InnerFrameSizeY;
	G4double petsizeZ = (1e-4)*cm;

	// Double-sided adhesive
	// parameters of Polyacrylate
	G4double polysizeX = InnerFrameSizeX;
	G4double polysizeY = InnerFrameSizeY;
	G4double polysizeZ = (1e-4)*cm;

	// parameters of TPC
	G4double TPCSizeZ = 7*cm;	
	G4double TPCGasThickness = TPCSizeZ-FrameStepHeight-petsizeZ;
	G4double EffectiveSizeX = 12.*cm;
	G4double EffectiveSizeY = 12.*cm;

	// parameters of the micromegas
	G4double MMSizeZ = 0.9*cm;
	G4double MMGasThickness = 0.75*cm;
	G4double MMEffectiveSizeX = 15.*cm;
	G4double MMEffectiveSizeY = 15.*cm;

	// parameters of the PCB board
	G4double PCBsizeX = 24*cm;
	G4double PCBsizeY = 24*cm;
	G4double PCBdeltaY = 0*cm;
	G4double PCBthickness = 0.273*cm;

	// parameters of Cu board
	G4double CusizeX = PCBsizeX;
	G4double CusizeY = PCBsizeY;
	G4double CudeltaY = PCBdeltaY;
	// G4double Cuthickness = 10e-4*cm;
	G4double Cuthickness = 40e-4*cm;
	
	// parameter of the Cu board below TPC
	G4double gap = 2*cm;
	G4double Cuthickness2 = 1*cm;

	// parameter of the field cage
	G4double FieldCageSize = 17*cm;
	G4double FieldCageThickness = 0.22*cm;
	G4double FieldCageHeight = 6.2*cm;

	G4double source_thickness = 0.2*cm;
	G4double source_radius = 1*cm;

	G4double ShellSize = 30*cm;
	G4double ShellThickness = 1*cm;
	G4double ShellHeight = 7*cm;

	G4double Cuthickness3 = 1*cm;
	
/*
*/
//=================Added part: a Aluminum plate to simulate the beta source ==========================

	// G4ThreeVector positionsource = G4ThreeVector(0., 0., -TPCSizeZ-0.5*source_thickness);

	// G4Tubs* solidsource = new G4Tubs(
	// 		"beta_source",                   // its name
	// 		0,                 // r_min 
	// 		source_radius,                 // r_max
	// 		source_thickness,				//height
	// 		0.,
	// 		360.);

	// G4LogicalVolume* logicsource = new G4LogicalVolume(
	// 		solidsource,                                    // its solid
	// 		Al,                                    // its material
	// 		"beta_source");                                      // its name

	// new G4PVPlacement(
	// 		0,                                           // no rotation
	// 		positionsource,                                 // at (0,0,0)
	// 		logicsource,                                    // its logical volume
	// 		"beta_source",                                       // its name
	// 		logicWorld,                                  // its mother  volume
	// 		false,                                       // no boolean operation
	// 		0);                                          // copy number

//=============================================

//=======Copper board below TPC (if exists)=====================

	// G4ThreeVector positionCubrd2 = G4ThreeVector(0., CudeltaY, -TPCSizeZ-gap-0.5*Cuthickness2);	

	// G4Box* solidCubrd2 = new G4Box("Cu_board2",                                    // its name
	// 		0.5*CusizeX, 0.5*CusizeY, 0.5*Cuthickness2);                      // its size

	// G4LogicalVolume* logicCubrd2 = new G4LogicalVolume(
	// 		solidCubrd2,                                    // its solid
	// 		Cu,                                    // its material
	// 		"Cu_board2");                                      // its name

	// new G4PVPlacement(
	// 		0,                                           // no rotation
	// 		positionCubrd2,                                 // at (0,0,0)
	// 		logicCubrd2,                                    // its logical volume
	// 		"Cu_board2",                                       // its name
	// 		logicWorld,                                  // its mother  volume
	// 		false,                                       // no boolean operation
	// 		0);                                          // copy number

//=============================================

//=======Copper board above anti-coincident detector (if exists)=====================

	// G4ThreeVector positionCubrd3 = G4ThreeVector(0., CudeltaY, PCBthickness+Cuthickness+PCBthickness+MMGasThickness+gap);	

	// G4Box* solidCubrd3 = new G4Box("Cu_board3",                                    // its name
	// 		0.5*CusizeX*1.2, 0.5*CusizeY*1.2, 0.5*Cuthickness3);                      // its size

	// G4LogicalVolume* logicCubrd3 = new G4LogicalVolume(
	// 		solidCubrd3,                                    // its solid
	// 		Cu,                                    // its material
	// 		"Cu_board3");                                      // its name

	// new G4PVPlacement(
	// 		0,                                           // no rotation
	// 		positionCubrd3,                                 // at (0,0,0)
	// 		logicCubrd3,                                    // its logical volume
	// 		"Cu_board3",                                       // its name
	// 		logicWorld,                                  // its mother  volume
	// 		false,                                       // no boolean operation
	// 		0);                                          // copy number

//=============================================

//=============Added part: a copper shell around the TPC ==================

	// G4ThreeVector positionShell = G4ThreeVector(0., 0., -0.5*ShellHeight+(PCBthickness+Cuthickness+PCBthickness+MMGasThickness+gap-0.5*Cuthickness3));		//(0,0,-3.35)cm center, top at z=0

	// G4Box* solidShellOut = new G4Box("ShellOut",                                    // its name
	// 		0.5*CusizeX*1.2, 0.5*CusizeY*1.2, 0.5*ShellHeight);                      // its size
	
	// G4Box* solidShellIn = new G4Box("ShellIn",                                    // its name
	// 		0.5*CusizeX*1.2-ShellThickness, 0.5*CusizeY*1.2-ShellThickness, 0.5*ShellHeight);                      // its size

	// G4SubtractionSolid* solidShell = new G4SubtractionSolid("CuShell",solidShellOut,solidShellIn);
	
	// G4LogicalVolume* logicShell = new G4LogicalVolume(
	// 		solidShell,                                    // its solid
	// 		Cu,                                    // its material
	// 		"CuShell");                                      // its name

	// new G4PVPlacement(
	// 		0,                                           // no rotation
	// 		positionShell,                                 // at (0,0,0)
	// 		logicShell,                                    // its logical volume
	// 		"CuShell",                                       // its name
	// 		logicWorld,                                  // its mother  volume
	// 		false,                                       // no boolean operation
	// 		0);                                          // copy number

//==============================================================

//=============PART1: The TPC detector=========

	// Gas chamber
	G4ThreeVector positionTPC = G4ThreeVector(0., 0., -0.5*TPCSizeZ);		//(0,0,-3.35)cm center, top at z=0

	G4Box* solidframe = new G4Box("AlFrame",                                    // its name
			0.5*FrameSizeX, 0.5*FrameSizeY, 0.5*TPCSizeZ);                      // its size

	G4LogicalVolume* logicframe = new G4LogicalVolume(
			solidframe,                                    // its solid
			Al,                                    // its material
			"AlFrame");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionTPC,                                 // at (0,0,0)
			logicframe,                                    // its logical volume
			"AlFrame",                                       // its name
			logicWorld,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number


	G4ThreeVector positionStep = G4ThreeVector(0., 0., -TPCSizeZ+0.5*FrameStepHeight);		
	
	G4Box* solidstep = new G4Box("Step",
			0.5*WindowSizeX, 0.5*WindowSizeY, 0.5*FrameStepHeight);

	G4LogicalVolume* logicstep = new G4LogicalVolume(
			solidstep,                                    // its solid
			Air,                                    // its material
			"Step");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionStep-positionTPC,                                 // at (0,0,0)
			logicstep,                                    // its logical volume
			"Step",                                       // its name
			logicframe,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	G4ThreeVector positionPET = G4ThreeVector(0., 0., -TPCSizeZ+FrameStepHeight+0.5*petsizeZ);				//-0.20835cm

	G4Box* solidPET = new G4Box("PET",                                    // its name
			0.5*petsizeX, 0.5*petsizeY, 0.5*petsizeZ);                      // its size


	G4LogicalVolume* logicPET = new G4LogicalVolume(
			solidPET,                                     // its solid
			PET,                                          // its material
			"PET");                                       // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionPET-positionTPC,                                 // at (0,0,0)
			logicPET,                                    // its logical volume
			"PET",                                       // its name
			logicframe,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number


	// full gas volume
	G4ThreeVector positionGas0 = G4ThreeVector(0., 0., -0.5*TPCGasThickness);

	G4Box* solidGas0 = new G4Box("Gas0",                                    // its name
			0.5*InnerFrameSizeX, 0.5*InnerFrameSizeY, 0.5*TPCGasThickness);                      // its size

	G4LogicalVolume* logicGas0 = new G4LogicalVolume(
			solidGas0,                                    // its solid
			AriC4H10,                                    // its material
			"Gas0");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionGas0-positionTPC,                                 // at (0,0,0)
			logicGas0,                                    // its logical volume
			"Gas0",                                       // its name
			logicframe,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	// field cage

	G4Box* solidFieldCageOut = new G4Box("FieldCageOut",                                    // its name
			0.5*FieldCageSize+FieldCageThickness, 0.5*FieldCageSize+FieldCageThickness, 0.5*FieldCageHeight);                      // its size
	
	G4Box* solidFieldCageIn = new G4Box("FieldCageIn",                                    // its name
			0.5*FieldCageSize, 0.5*FieldCageSize, 0.5*FieldCageHeight);                      // its size

	G4SubtractionSolid* solidFieldCage = new G4SubtractionSolid("FieldCage",solidFieldCageOut,solidFieldCageIn);
	
	G4LogicalVolume* logicFieldCage = new G4LogicalVolume(
			solidFieldCage,                                    // its solid
			FR4,                                    // its material
			"FieldCage");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			G4ThreeVector(),                                 // at (0,0,0)
			logicFieldCage,                                    // its logical volume
			"FieldCage",                                       // its name
			logicGas0,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number


	// detector sensitive volume
	G4Box* solidGas = new G4Box("Gas",                                    // its name
			0.5*DetectorSizeX, 0.5*DetectorSizeY, 0.5*TPCGasThickness);                      // its size

	G4LogicalVolume* logicGas = new G4LogicalVolume(
			solidGas,                                    // its solid
			AriC4H10,                                    // its material
			"Gas");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			G4ThreeVector(),                                 // at (0,0,0)
			logicGas,                                    // its logical volume
			"Gas",                                       // its name
			logicGas0,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number


	//selected TPC sensitive volume
	G4Box* solidGasEff = new G4Box("GasEff",                                    // its name
			0.5*EffectiveSizeX, 0.5*EffectiveSizeY, 0.5*TPCGasThickness);                      // its size

	G4LogicalVolume* logicGasEff = new G4LogicalVolume(
			solidGasEff,                                    // its solid
			AriC4H10,                                    // its material
			"GasEff");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			G4ThreeVector(),                                 // at (0,0,0)
			logicGasEff,                                    // its logical volume
			"GasEff",                                       // its name
			logicGas,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	
	

	//	PCB board

	G4ThreeVector positionPCB = G4ThreeVector(0., PCBdeltaY, 0.5*PCBthickness);				//at (0,0,-6.8365)		
	G4Box* solidPCB = new G4Box("PCB",                                    // its name
			0.5*PCBsizeX, 0.5*PCBsizeY, 0.5*PCBthickness);                      // its size

	G4LogicalVolume* logicPCB = new G4LogicalVolume(
			solidPCB,                                    // its solid
			Kapton,                                    // its material
			"PCB");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionPCB,                                 // at (0,0,0)
			logicPCB,                                    // its logical volume
			"PCB",                                       // its name
			logicWorld,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	//	Copper board

	G4ThreeVector positionCubrd = G4ThreeVector(0., CudeltaY, PCBthickness+0.5*Cuthickness);	

	G4Box* solidCubrd = new G4Box("Cu_board",                                    // its name
			0.5*CusizeX, 0.5*CusizeY, 0.5*Cuthickness);                      // its size

	G4LogicalVolume* logicCubrd = new G4LogicalVolume(
			solidCubrd,                                    // its solid
			Cu,                                    // its material
			"Cu_board");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionCubrd,                                 // at (0,0,0)
			logicCubrd,                                    // its logical volume
			"Cu_board",                                       // its name
			logicWorld,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number
			
	//	PCB board

	G4ThreeVector positionPCB2 = G4ThreeVector(0., PCBdeltaY, PCBthickness+Cuthickness+0.5*PCBthickness);				//at (0,0,-6.8365)		
	G4Box* solidPCB2 = new G4Box("PCB2",                                    // its name
			0.5*PCBsizeX, 0.5*PCBsizeY, 0.5*PCBthickness);                      // its size

	G4LogicalVolume* logicPCB2 = new G4LogicalVolume(
			solidPCB2,                                    // its solid
			Kapton,                                    // its material
			"PCB2");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionPCB2,                                 // at (0,0,0)
			logicPCB2,                                    // its logical volume
			"PCB2",                                       // its name
			logicWorld,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	
	// anticoincident MM

	G4ThreeVector positionMM = G4ThreeVector(0., 0., PCBthickness+Cuthickness+PCBthickness+0.5*MMSizeZ);

	G4Box* solidMM = new G4Box("MM",                                    // its name
			0.5*FrameSizeX, 0.5*FrameSizeY, 0.5*MMSizeZ);                      // its size

	G4LogicalVolume* logicMM = new G4LogicalVolume(
			solidMM,                                    // its solid
			Al,                                    // its material
			"MM");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionMM,                                 // at (0,0,0)
			logicMM,                                    // its logical volume
			"MM",                                       // its name
			logicWorld,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	// full gas volume of MM

	G4ThreeVector positionMMGas0 = G4ThreeVector(0., 0., PCBthickness+Cuthickness+PCBthickness+0.5*MMGasThickness);		
	
	G4Box* solidMMGas0 = new G4Box("MMGas0",
			0.5*InnerFrameSizeX, 0.5*InnerFrameSizeY, 0.5*MMGasThickness);

	G4LogicalVolume* logicMMGas0 = new G4LogicalVolume(
			solidMMGas0,                                    // its solid
			AriC4H10,                                    // its material
			"MMGas0");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionMMGas0-positionMM,                                 // at (0,0,0)
			logicMMGas0,                                    // its logical volume
			"MMGas0",                                       // its name
			logicMM,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	// sensitive volume of the MM

	G4ThreeVector positionMMGas = G4ThreeVector(0., 0., PCBthickness+Cuthickness+PCBthickness+0.5*MMGasThickness);		
	
	G4Box* solidMMGas = new G4Box("GasEff2",
			0.5*MMEffectiveSizeX, 0.5*MMEffectiveSizeX, 0.5*MMGasThickness);

	G4LogicalVolume* logicMMGas = new G4LogicalVolume(
			solidMMGas,                                    // its solid
			AriC4H10,                                    // its material
			"GasEff2");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			G4ThreeVector(),                                 // at (0,0,0)
			logicMMGas,                                    // its logical volume
			"GasEff2",                                       // its name
			logicMMGas0,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number





// ==========================================================



// ===============================================================


	//-----Set the step limits in the Gas volume-------------
	G4double maxStep = 1.0*mm;
	fStepLimits = new G4UserLimits(maxStep);
	logicGasEff->SetUserLimits(fStepLimits);
	//-------------------------------



	G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0.9,0.0,0.0));
	visAttributes->SetVisibility(false);
	logicWorld->SetVisAttributes(visAttributes);

	// visAttributes = new G4VisAttributes(G4Colour(0.0,0.0,1.0)); // red
	// logicCollimation->SetVisAttributes(visAttributes);

	visAttributes = new G4VisAttributes(G4Colour(1.0,0.0,1.0)); 
	logicGas->SetVisAttributes(visAttributes);
	visAttributes = new G4VisAttributes(G4Colour(0.0,1.0,0.0)); 
	logicGasEff->SetVisAttributes(visAttributes);
	visAttributes = new G4VisAttributes(G4Colour(1.0,0.0,1.0)); 
	logicMMGas->SetVisAttributes(visAttributes);


	//
	//always return the physical World
	//  
	return physicalWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField(){
	// Create global magnetic field messenger.
  	// Uniform magnetic field is then created automatically if
  	// the field value is not zero.
  	G4ThreeVector fieldValue = G4ThreeVector();
  	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  	fMagFieldMessenger->SetVerboseLevel(1);
	
  	// Register the field messenger for deleting
  	G4AutoDelete::Register(fMagFieldMessenger);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......