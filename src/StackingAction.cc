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
// StackingAction.cc
// 
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include <iostream>
#include <iomanip>

#include "RunAction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "StackingMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Material.hh"
#include "G4RunManager.hh"

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(EventAction* EA, DetectorConstruction* DE )
:eventaction(EA), detector(DE)
{
  killSecondary  = 0;
  stackMessenger = new StackingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{
  delete stackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
	// Lindhard - Robinson partition function for any ion
G4double StackingAction::PartitionFunction(G4double T, G4double M1, G4double M2, G4double Z1, G4double Z2)
{

  G4double E_nu;

  G4double e_ch_squared = 1.4399764e-6;//*(MeV*nm);
  //screening distance
  G4double a;
  G4double a_0 = 0.0529*nm;   // Bohr radius
  a = 0.8853*a_0/(std::pow(std::pow(Z1,(2./3.))+std::pow(Z2,(2./3.)),(1./2.)))/nm;  // Thomas-Fermi
  //G4double x = 0.23;
  //a = 0.8853*a_0/(std::pow(Z1,x)+std::pow(Z2,x));  // ZBL
  G4double E_l = (Z1*Z2*e_ch_squared/a)*(M1+M2)/M2;
  G4double first_term = std::pow(M1+M2,(3./2.))*std::pow(Z1,(2./3.))*std::pow(Z2,(1./2.));
  G4double second_term = std::pow(M1,(3./2.))*std::pow(std::pow(Z1,(2./3.))+std::pow(Z2,(2./3.)),(3./4.));
  G4double m_0 = 0.510998950; 		// electron rest mass
  G4double M2_x = 931.49410242*M2;	// target atom mass
  G4double k_l = (32/(3*pi))*std::pow((m_0/M2_x),(1./2.))*first_term/second_term;

  G4double T_1 = T/eV;	// kinetic energy of projectile in eV
  G4double T_El = T_1/(E_l*1000000);

  G4double g_l = T_El + 0.40244*std::pow(T_El,(3./4.))+3.4008*std::pow(T_El,(1./6.));


  E_nu = 1./(1.+k_l*g_l);

  return E_nu;
}

	// Akerman partition function
	//https://ieeexplore.ieee.org/document/4033183
G4double StackingAction::PartitionFunction2(G4double T, G4double M1, G4double M2, G4double Z1, G4double Z2)
{
  G4double E_nu;
  G4double e_ch_squared = 1.4399764e-6;//*(MeV*nm);
  //screening distance
  G4double a;
  G4double a_0 = 0.0529*nm;   // Bohr radius
  // Thomas-Fermi
  a = (0.8853*a_0/(std::pow(std::pow(Z1,(2./3.))+std::pow(Z2,(2./3.)),(0.5))))/nm;
  // ZBL
  //G4double x = 0.23;
  //a = 0.8853*a_0/(std::pow(Z1,x)+std::pow(Z2,x));

  G4double E_l = (Z1*Z2*e_ch_squared/a)*(M1+M2)/M2;
  G4double first_term = std::pow(M1+M2,(3./2.))*std::pow(Z1,(1./3.))*std::pow(Z2,(0.5));
  G4double second_term = std::pow(M1,(3./2.))*std::pow(std::pow(Z1,(2./3.))+std::pow(Z2,(2./3.)),(3./4.));
  G4double m_0 = 0.510998950; 		// electron rest mass
  G4double M2_x = 931.49410242*M2;	// target atom mass
  G4double k_l = (32/(3*pi))*std::pow((m_0/M2_x),(0.5))*first_term/second_term;
  

  G4double T_1 = T/eV;
  G4double T_El = T_1/(E_l*1000000);

  G4double g_n = 0.74422*T_El + 1.6812*std::pow(T_El,(0.75))+0.90565*std::pow(T_El,(1./6.));


  E_nu = 1./(1.+k_l*g_n);

  return E_nu;
}


	// old Lindhard - Robinson partition function for protons
G4double  StackingAction::DamageEnergy(G4double T,G4double A, G4double Z)
{
  //.................. T in  eV!!!!!!!!!!!!!
  G4double Z2= Z;
  G4double M2= A;
  G4double k_d;
  G4double epsilon_d;
  G4double g_epsilon_d;
  G4double E_nu;

	// original, for protons
  k_d=0.1334*std::pow(Z2,(2./3.))*std::pow(M2,(-1./2.));
  epsilon_d=0.01014*std::pow(Z2,(-7./3.))*(T/eV);
  g_epsilon_d= epsilon_d+0.40244*std::pow(epsilon_d,(3./4.))+3.4008*std::pow(epsilon_d,(1./6.));

  E_nu=1./(1.+ k_d*g_epsilon_d);

  return E_nu;//partition fraction!!!
}

//...................................................................

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
    G4int IDp= aTrack->GetParentID();
  //G4int TrID = aTrack->GetTrackID();

    G4double energy = aTrack->GetKineticEnergy();
    G4double charge = aTrack->GetDefinition()->GetPDGCharge();
    G4double u = 931.49410242;

    G4double A1 = aTrack->GetParticleDefinition()->GetPDGMass()/u;
    G4double Z1 = aTrack->GetParticleDefinition()->GetPDGCharge();

    G4double Spoint  = (aTrack->GetPosition()).mag();

    G4Material*  material  = detector->GetMaterialM(0);
    G4Material*  material1 = detector->GetMaterialM(1);
    G4Material*  material2 = detector->GetMaterialM(2);
    G4Material*  material3 = detector->GetMaterialM(3);
    G4Material*  material4 = detector->GetMaterialM(4);   
    
    //G4cout << " check " << G4endl;
    
    G4double mA2   = material->GetDensity()/(material->GetTotNbOfAtomsPerVolume()*amu);
    G4double mZ2   = material->GetTotNbOfElectPerVolume()/material->GetTotNbOfAtomsPerVolume();   
    
    G4double m1A2  = material1->GetDensity()/(material1->GetTotNbOfAtomsPerVolume()*amu);
    G4double m1Z2  = material1->GetTotNbOfElectPerVolume()/material1->GetTotNbOfAtomsPerVolume();       
    
    G4double m2A2  = material2->GetDensity()/(material2->GetTotNbOfAtomsPerVolume()*amu);
    G4double m2Z2  = material2->GetTotNbOfElectPerVolume()/material2->GetTotNbOfAtomsPerVolume();           
         
    G4double m3A2  = material3->GetDensity()/(material3->GetTotNbOfAtomsPerVolume()*amu);
    G4double m3Z2  = material3->GetTotNbOfElectPerVolume()/material3->GetTotNbOfAtomsPerVolume();   

    G4double m4A2  = material4->GetDensity()/(material4->GetTotNbOfAtomsPerVolume()*amu);
    G4double m4Z2  = material4->GetTotNbOfElectPerVolume()/material4->GetTotNbOfAtomsPerVolume();   

  //keep primary particle
  if (IDp == 0) {return fUrgent;}
  Run* run = static_cast<Run*>(
  G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  
  G4double partition = 0.;
  G4double dam_energy = 0.;
  	G4VPhysicalVolume* touchable_vol = aTrack->GetTouchableHandle()->GetVolume();
	if (IDp > 0 && charge > 0) 
	{

	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

		// MAIN LAYER
		if (aTrack->GetTouchableHandle()->GetVolume() == detector->GetIntAbsorber(0))
		{
			partition = PartitionFunction2(energy, A1,mA2,Z1,mZ2);
			dam_energy = partition * energy;
			
			run->abssumTL(dam_energy);	// damage energy
			run->absnbrec(1);		// recoil
			run->abssumT(energy);		// kinetic energy
			analysisManager->FillH1(14, Spoint, dam_energy);	
			run->NumberRec(1);	
		}
		// FIRST LAYER
		else if(touchable_vol == detector->GetIntAbsorber(1))
		{
			partition = PartitionFunction2(energy, A1,m1A2,Z1,m1Z2);
			dam_energy = partition * energy;
			
			run->abs1sumTL(dam_energy);	// damage energy
			run->abs1nbrec(1);		// recoil
			run->abs1sumT(energy);		// kinetic energy	
			analysisManager->FillH1(14, Spoint, dam_energy);	
			run->NumberRec(1);				
		}
		// SECOND LAYER
		else if(touchable_vol == detector->GetIntAbsorber(2))
		{
			partition = PartitionFunction2(energy, A1,m2A2,Z1,m2Z2);
			dam_energy = partition * energy;
			
			run->abs2sumTL(dam_energy);	// damage energy
			run->abs2nbrec(1);		// recoil
			run->abs2sumT(energy);		// kinetic energy		
			analysisManager->FillH1(14, Spoint, dam_energy);
			run->NumberRec(1);	
		}
		// THIRD LAYER
		else if(touchable_vol == detector->GetIntAbsorber(3))
		{
			partition = PartitionFunction2(energy, A1,m3A2,Z1,m3Z2);
			dam_energy = partition * energy;
			
			run->abs3sumTL(dam_energy);	// damage energy
			run->abs3nbrec(1);		// recoil
			run->abs3sumT(energy);		// kinetic energy		
			analysisManager->FillH1(14, Spoint, dam_energy);
			run->NumberRec(1);	
		}	
		// FOURTH LAYER
		else if(touchable_vol == detector->GetIntAbsorber(4))
		{
			partition = PartitionFunction2(energy, A1,m4A2,Z1,m4Z2);
			dam_energy = partition * energy;
			
			run->abs4sumTL(dam_energy);	// damage energy
			run->abs4nbrec(1);		// recoil
			run->abs4sumT(energy);		// kinetic energy		
			analysisManager->FillH1(14, Spoint, dam_energy);
			run->NumberRec(1);	
		}						
	}

  //stack or delete secondaries
  G4ClassificationOfNewTrack status = fUrgent;
  if (killSecondary) 
    {if (killSecondary == 1) {
      status = fKill;
        }
    }

  return status;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
