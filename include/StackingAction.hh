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
// StackingAction.hh
// 
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef StackingAction_h
#define StackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"

class Run;
class RunAction;
class EventAction;
class StackingMessenger;
class DetectorConstruction; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class StackingAction : public G4UserStackingAction
{
  public:
    StackingAction(EventAction*,DetectorConstruction* );
   ~StackingAction();
   
    void SetKillStatus(G4int value) {killSecondary = value;};
    
    
    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);
    
    G4double DamageEnergy(G4double T, G4double A, G4double Z);
  
    G4double PartitionFunction(G4double T, G4double M1, G4double M2, G4double Z1, G4double Z2);
    G4double PartitionFunction2(G4double T, G4double M1, G4double M2, G4double Z1, G4double Z2);


    
  private:
    //Run*          run;
    EventAction*        eventaction;    
    DetectorConstruction* detector;
    G4int               killSecondary;
    StackingMessenger*  stackMessenger;
    G4int rec;
    G4int SecRec;
    G4int absnbrec;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

