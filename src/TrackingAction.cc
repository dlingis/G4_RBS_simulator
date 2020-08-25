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
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
// $Id: TrackingAction.cc 69099 2013-04-18 12:25:19Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4StepStatus.hh"
#include "G4ParticleTypes.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(EventAction* event, DetectorConstruction* detector,PrimaryGeneratorAction* prim)
:G4UserTrackingAction(), fEventAction(event), fDetector(detector),fPrimary(prim)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{  
  
	

	// counts secondary particles in both materials
if (track->GetTouchableHandle()->GetVolume() == fDetector->GetAbsorber() || track->GetTouchableHandle()->GetVolume() == fDetector->GetAbsorber2() ) {
  //count secondary particles
  if (track->GetTrackID() == 1) return;  
  G4String name   = track->GetDefinition()->GetParticleName();
  G4double energy = track->GetKineticEnergy();
  Run* run = static_cast<Run*>(
      G4RunManager::GetRunManager()->GetNonConstCurrentRun());    
  run->ParticleCount(name,energy);

  G4double charge = track->GetDynamicParticle()->GetDefinition()->GetPDGCharge();
	if (charge > 0) {
		run->addsec();
	G4ThreeVector pos = track->GetVertexPosition(); // nustato pradine pozicija treko
/*	G4cout << " pozicija vertex " << G4BestUnit(pos,"Length") << " ir energija: " << G4BestUnit(energy,"Energy") << G4endl;
			const G4StepPoint* startPoint = track->GetStep()->GetPreStepPoint();
			const G4StepPoint* endPoint = track->GetStep()->GetPostStepPoint();
			G4ThreeVector end = endPoint->GetPosition();
			G4ThreeVector start = startPoint->GetPosition();
	G4cout << " pozicija step " << start << G4endl;
	G4cout << " pozicija position " << track->GetPosition() << G4endl;} */

   FILE * output=fopen("secondaries.txt","a");
   fprintf(output, "%s\t %f\t %f\t %f\t %f\t %f\n", name.c_str(), charge, energy/keV, pos.x()/um, pos.y()/um, pos.z()/um);
   fclose(output); 	
	}

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{

 Run* run = static_cast<Run*>(
              G4RunManager::GetRunManager()->GetNonConstCurrentRun());

	//remove("backscattered.txt");
	// pradzioje trekinimas pasaulyjejei pirmines daleles neiseina is absorberio, treko nera. Jeigu cia uzregistruoja, reiskia vyksta atgaline sklaida. 
if (track->GetTouchableHandle()->GetVolume() == fDetector->GetWorld()) {
	if (track->GetParentID() == 0) {
		G4ThreeVector momentas = track->GetMomentumDirection();
		G4double momz = momentas.z();
			if (momz <= 0) {
			// printinti backscattered particles i faila
				G4String name   = track->GetDefinition()->GetParticleName();
 				G4double energy = track->GetKineticEnergy();


				//remove("backscattered.txt");
   				FILE * output=fopen("backscattered.txt","a");
   				fprintf(output, "%s\t %f\t %f\n", name.c_str(), energy/keV, std::acos(momz)/degree);
   				fclose(output); 
				//G4cout << " Vyksta atgaline sklaida. Momentas " << momz << G4endl;
				run->addrbs();
					}
		}
}

  if (track->GetTouchableHandle()->GetVolume() == fDetector->GetAbsorber() || track->GetTouchableHandle()->GetVolume() == fDetector->GetAbsorber2() ) 
	{
  		G4String parName = track->GetDefinition()->GetParticleName();

	//pliusas, nes -z kryptimi pluostelis sklinda
  G4double Trleng = track->GetTrackLength()+fPrimary->GetParticleGun()->GetParticlePosition().z()+0.5*fDetector->GetLength(0);

  if (track->GetParentID() == 0) {
    fEventAction->AddTrueTrakLen(Trleng);
    G4ThreeVector vertex = track->GetVertexPosition();
    G4ThreeVector position = track->GetPosition()+0.5*fDetector->GetDimensions(0);// - vertex;    //pakeista is   track->GetPosition() - vertex;
    fEventAction->AddProjTrakLen(position.z());
	} 

 if (track->GetTrackID() == 1) {
   G4double z = track->GetPosition().z() + 0.5*fDetector->GetLength(0);
   run->AddProjectedRange(z);
 	G4AnalysisManager* analysis = G4AnalysisManager::Instance();
   analysis->FillH1(16, z);
 }


  // buvo vx, keičiu į z ašį pagal mano geometriją
  //scattering angle of  primary particle
  if (track->GetParentID() == 0 ){

    G4double theta = 0.0;
    G4double    vz = (track->GetMomentumDirection()).z(); 
    //G4cout<<"Angle_Lab: "<<std::acos(vx)<<" Position_x: "<<px<<G4endl;
    if(vz <= -1.0)    { theta = -CLHEP::pi; }
    else if(vz < 1.0) { theta = std::acos(vz); }
    fEventAction->AddTheta(theta);   
  }
  }
 // keep only outgoing particle
 G4StepStatus status = track->GetStep()->GetPostStepPoint()->GetStepStatus();
 if (status != fWorldBoundary) return; 

	if (track->GetParentID() == 0) 
		{ run->countEmerging();}
 
 const G4ParticleDefinition* particle = track->GetParticleDefinition();
 G4String name   = particle->GetParticleName();
 G4double energy = track->GetKineticEnergy();

 fEventAction->AddEflow(energy);   
 run->ParticleFlux(name,energy);               
 
 // histograms: enery flow
 //
 G4AnalysisManager* analysis = G4AnalysisManager::Instance();
 
 G4int ih = 0; 
 G4String type   = particle->GetParticleType();      
 G4double charge = particle->GetPDGCharge();
 if (charge > 3.)  ih = 10; 
 else if (particle == G4Gamma::Gamma())       ih = 4;
 else if (particle == G4Electron::Electron()) ih = 5;
 else if (particle == G4Positron::Positron()) ih = 5;  
 else if (particle == G4Neutron::Neutron())   ih = 6;
 else if (particle == G4Proton::Proton())     ih = 7;
 else if (particle == G4Deuteron::Deuteron()) ih = 8;
 else if (particle == G4Alpha::Alpha())       ih = 9;       
 else if (type == "nucleus")                  ih = 10;
 else if (type == "baryon")                   ih = 11;         
 else if (type == "meson")                    ih = 12;
 else if (type == "lepton")                   ih = 13;        
 if (ih > 0) analysis->FillH1(ih,energy);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

