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
/// \file Hadr06.cc
/// \brief Main program of the hadronic/Hadr06 example
//
//
// $Id: TestEm1.cc,v 1.16 2010-04-06 11:11:24 maire Exp $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
//#else
#include "G4RunManager.hh"
//#endif

#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "SteppingVerbose.hh"

#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BIC.hh"

#ifdef G4VIS_USE
 #include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif


// channe
#include "G4Types.hh"
#include "G4ScoringManager.hh"
#include "G4EmStandardPhysics_option4_channeling.hh"
#include "G4GenericBiasingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"

#include "G4EmStandardPhysics_option4.hh"




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
int main(int argc,char** argv) {
 
  //choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  //G4Random::setTheEngine(new CLHEP::RanshiEngine);
 
  // Construct the default run manager
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  G4int nThreads = G4Threading::G4GetNumberOfCores();
  if (argc==3) nThreads = G4UIcommand::ConvertToInt(argv[2]);
  runManager->SetNumberOfThreads(nThreads-2);
#else
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
  G4RunManager* runManager = new G4RunManager;
#endif

  //****************
  // channeling bullshit
    G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
    scManager->SetVerboseLevel(0);
    
    // Set mandatory initialization classes
    //G4VModularPhysicsList* physlist= new FTFP_BERT();
    G4VModularPhysicsList* physlist= new QGSP_BIC();
    //G4VModularPhysicsList* physlist= new G4EmStandardPhysics_option4();
    //G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
    //physlist->RegisterPhysics(new G4ChannelingPhysics());
    physlist->ReplacePhysics(new G4EmStandardPhysics_option4_channeling());
    //physlist->ReplacePhysics(new G4EmStandardPhysics_option4());
    //biasingPhysics->PhysicsBiasAllCharged();
    //physlist->RegisterPhysics(biasingPhysics);
    runManager->SetUserInitialization(physlist);

  // set mandatory initialization classes
  DetectorConstruction* det= new DetectorConstruction;
  runManager->SetUserInitialization(det);
  //runManager->SetUserInitialization(chan);
    
  runManager->SetUserInitialization(new ActionInitialization(det));    
     
  // get the pointer to the User Interface manager 
    G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc!=1)   // batch mode  
    {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
    }
    
  else           //define visualization and UI terminal for interactive mode
    { 
#ifdef G4VIS_USE
      G4VisManager* visManager = new G4VisExecutive;
      visManager->Initialize();
#endif    
     
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);      
      ui->SessionStart();
      delete ui;
#endif
          
#ifdef G4VIS_USE
     delete visManager;
#endif     
    }

  // job termination 
  //
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
