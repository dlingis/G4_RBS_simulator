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
/// \file HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 67909 2013-03-12 18:51:09Z vnivanch $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("Hadr06")
{

G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysisManager->GetType() << G4endl;
    
    //** Create directories **//
    analysisManager->SetVerboseLevel(1);
    //analysisManager->SetNtupleMerging(true);

  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  //analysisManager->SetFirstHistoId(1);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);     //enable inactivation of histograms


  // Define histograms start values
  const G4int kMaxHisto = 51; // 47
  const G4String id[] = {"0",
  			"Total_energy_deposit",
  			"EDEP_(MeV/mm)",
  			"3",
  			"4",
  			"5",
  			"6",
  			"7",
  			"8",
  			"9",
                       "10",
                       "11",
                       "12",
                       "13",
                       "14",
                       "15",
                       "Projected_range",
                       "17",
                       "18",
			"19",
			"RBS_total_en",
			"RBS_total_dist",
			"RBS_subs_en",
			"RBS_subs_dist",
			"RBS_subs_el1_en",
			"RBS_subs_el2_en",
			"RBS_subs_el1_dist",
			"RBS_subs_el2_dist",
			"RBS_la1_en",
			"RBS_la1_dist",
			"RBS_la1_el1_en",
			"RBS_la1_el2_en",
			"RBS_la1_el1_dist",
			"RBS_la1_el2_dist",
			"RBS_la2_en",
			"RBS_la2_dist",
			"RBS_la2_el1_en",
			"RBS_la2_el2_en",
			"RBS_la2_el1_dist",
			"RBS_la2_el2_dist",
			"RBS_la3_en",
			"RBS_la3_dist",
			"RBS_la3_el1_en",
			"RBS_la3_el2_en",
			"RBS_la3_el1_dist",
			"RBS_la3_el2_dist",
			"RBS_la4_en",
			"RBS_la4_dist",
			"RBS_la4_el1_en",
			"RBS_la4_el2_en",
			"RBS_la4_el1_dist",
			"RBS_la4_el2_dist",
			
};
  const G4String title[] = 
                { "dummy",                                        		//0
                  "total energy deposit",                         		//1
                  "Edep (MeV/mm) along absorber",                 		//2
                  "total kinetic energy flow",                    		//3
                  "gamma flux (dN/dE) at exit",                   		//4
                  "e+- flux (dN/dE) at exit",                     		//5
                  "neutrons flux (dN/dE) at exit",                		//6
                  "protons flux (dN/dE) at exit",                 		//7
                  "deuterons flux (dN/dE) at exit",               		//8
                  "alphas flux (dN/dE) at exit",                  		//9
                  "all others ions flux (dN/dE) at exit",         		//10
                  "all others baryons flux (dN/dE) at exit",      		//11
                  "all others mesons flux (dN/dE) at exit",       		//12
                  "all others leptons flux (dN/dE) at exit",       		//13 
 		  "Frenkelio poru koncentracija nuo saltinio(#/cm3)",         //14
		  " Sklaidos kampas ",						//15
		  " Projected range distribution ",				//16	
		  " NIEL pagal atstuma, keV	",				//17 
		  " NIEL pagal energija",					//18
		  " EDEP pagal energija ",					//19
		  " Total RBS Spectrum ",					//20
		  " Total RBS Spectrum distance", 				//21
		  " RBS spectrum substrate en",				//22
		  " RBS spectrum substrate distance",				//23
		  " RBS spectrum substrate el 1",				//24
		  " RBS spectrum substrate el 2",				//25
		  " RBS spectrum substrate el 1 distance",			//26
		  " RBS spectrum substrate el 2 distance",			//27
		  " RBS Layer 1 en",						//28
		  " RBS Layer 1 distance",					//29
		  " RBS Layer 1 el 1 ",					//30
		  " RBS Layer 1 el 2 ",					//31
		  " RBS Layer 1 el 1 distance",				//32
		  " RBS Layer 1 el 2 distance",				//33
		  " RBS Layer 2 en",						//34
		  " RBS Layer 2 distance",					//35
		  " RBS Layer 2 el 1 ",					//36
		  " RBS Layer 2 el 2 ",					//37
		  " RBS Layer 2 el 1 distance",				//38
		  " RBS Layer 2 el 2 distance",				//39
		  " RBS Layer 3 en",						//40
		  " RBS Layer 3 distance",					//41
		  " RBS Layer 3 el 1 ",					//42
		  " RBS Layer 3 el 2 ",					//43
		  " RBS Layer 3 el 1 distance",				//44
		  " RBS Layer 3 el 2 distance",				//45		  
		  " RBS Layer 4 en",						//46
		  " RBS Layer 4 distance",					//47
		  " RBS Layer 4 el 1 ",					//48
		  " RBS Layer 4 el 2 ",					//49
		  " RBS Layer 4 el 1 distance",				//50
		  " RBS Layer 4 el 2 distance"				//51		  
		};  

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
