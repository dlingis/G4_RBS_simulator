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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 71404 2013-06-14 16:56:38Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
                           
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* event, DetectorConstruction* DET)
: G4UserSteppingAction(), fEventAction(event), detector(DET)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

   //apsauga nuo neigiamu verciu
   G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep <= 0.) return;


	G4VPhysicalVolume* preMAT = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
	G4VPhysicalVolume* posMAT = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume();


 	G4double charge = aStep->GetTrack()->GetDynamicParticle()->GetCharge(); //GetDefinition()->GetPDGCharge()
  	G4String parName = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();

  	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  	G4double niel = aStep->GetNonIonizingEnergyDeposit();


	G4double detsize = detector->GetLength(0)/2;
	G4ThreeVector end = aStep->GetPostStepPoint()->GetPosition();
	G4ThreeVector start = aStep->GetPreStepPoint()->GetPosition();
	G4ThreeVector point = start + G4UniformRand()*(end - start);

	G4double endPos = end.z()+detsize;
	G4double startPos = start.z()+detsize;

    	G4double x  = startPos + G4UniformRand()*(endPos-startPos);
    	analysisManager->FillH1(39, x, niel);

	// primary particles 
  	G4int IDp =  aStep->GetTrack()->GetParentID();
  	Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
             
	// main volume
	if(preMAT == detector->GetIntAbsorber(0) && posMAT == detector->GetIntAbsorber(0) )
	{
		if(IDp == 0)
		{
			run->absLEN(aStep->GetStepLength());
			run->absION(aStep->GetTotalEnergyDeposit());
			run->absNON(aStep->GetNonIonizingEnergyDeposit());
			fEventAction->absSTP();
		}
	}
	
	// layer1
	if(preMAT == detector->GetIntAbsorber(1) && posMAT == detector->GetIntAbsorber(1) )
	{
		if(IDp == 0)
		{
			run->abs1LEN(aStep->GetStepLength());
			run->abs1ION(aStep->GetTotalEnergyDeposit());
			run->abs1NON(aStep->GetNonIonizingEnergyDeposit());
			fEventAction->abs1STP();
		}
	}
	
	// layer2
	if(preMAT == detector->GetIntAbsorber(2) && posMAT == detector->GetIntAbsorber(2) )
	{
		if(IDp == 0)
		{
			run->abs2LEN(aStep->GetStepLength());
			run->abs2ION(aStep->GetTotalEnergyDeposit());
			run->abs2NON(aStep->GetNonIonizingEnergyDeposit());
			fEventAction->abs2STP();
		}
	}		
	// layer3
	if(preMAT == detector->GetIntAbsorber(3) && posMAT == detector->GetIntAbsorber(3) )
	{
		if(IDp == 0)
		{
			run->abs3LEN(aStep->GetStepLength());
			run->abs3ION(aStep->GetTotalEnergyDeposit());
			run->abs3NON(aStep->GetNonIonizingEnergyDeposit());
			fEventAction->abs3STP();
		}
	}		
	// layer4
	if(preMAT == detector->GetIntAbsorber(4) && posMAT == detector->GetIntAbsorber(4) )
	{
		if(IDp == 0)
		{
			run->abs4LEN(aStep->GetStepLength());
			run->abs4ION(aStep->GetTotalEnergyDeposit());
			run->abs4NON(aStep->GetNonIonizingEnergyDeposit());
			fEventAction->abs4STP();
		}
	}

  	// primary particles
   	if (IDp==0 )
   	{//*****************************************
   	
		run->AddTotStep();
	   	run->AddTrakLenPrim(aStep->GetStepLength());
	   	fEventAction->CountStepsPrim();
	   	run->AddEnergy(aStep->GetTotalEnergyDeposit());
	   	run->AddNonIonEnergy(aStep->GetNonIonizingEnergyDeposit());
	} //**********************************
   	if (IDp > 0 && charge > 0)
   	{
   		run->AddTrakLenSec(aStep->GetStepLength()); 
   		//fEventAction->AddTrakLenSec(aStep->GetStepLength()); 
	}

  	const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  	const G4VProcess* process   = endPoint->GetProcessDefinedStep();
  	run->CountProcesses(process);
  
   	if (IDp == 0)
   	{// primary energy deposition
  		G4double edepStep = aStep->GetTotalEnergyDeposit();
  		G4double nielStep = aStep->GetNonIonizingEnergyDeposit();

		if (nielStep > 0.)
		{ 
			G4ThreeVector prePoint  = aStep->GetPreStepPoint()->GetPosition();
 			G4ThreeVector postPoint = aStep->GetPostStepPoint()->GetPosition();
 			G4double r = point.z();
			// NIEL
 			analysisManager->FillH1(17, r+detsize, nielStep/keV);
			analysisManager->FillH1(18, nielStep);
			// IONiZING energy loss
 			analysisManager->FillH1(2, r+detsize, edepStep);
 			analysisManager->FillH1(19, edepStep);
		}
	}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


