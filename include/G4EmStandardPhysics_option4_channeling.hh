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
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmStandardPhysics_option4_channeling
//
// Author:      V.Ivanchenko 28.09.2012
//
// Modified:
//
//----------------------------------------------------------------------------
//
// This class provides construction of EM physics using the best models
// of standard and low-energy packages and set of 
// the most adavced options allowing precise simulation at low
// and intermediate energies
//

#ifndef G4EmStandardPhysics_option4_channeling_h
#define G4EmStandardPhysics_option4_channeling_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include "G4PhysicsFreeVector.hh"
#include "G4ExtDEDXTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




class G4EmStandardPhysics_option4_channeling : public G4VPhysicsConstructor
{
public:

  explicit G4EmStandardPhysics_option4_channeling(G4int ver=1, const G4String& name="");

  virtual ~G4EmStandardPhysics_option4_channeling();

  virtual void ConstructParticle();
  virtual void ConstructProcess();
  
  //G4ExtDEDXTable* srim_table = new G4ExtDEDXTable();

private:
  G4int  verbose;
  /*
    G4PhysicsFreeVector* fVectorLiHf;
    G4PhysicsFreeVector* fVectorLiO;
    G4PhysicsFreeVector* fVectorLiAr;
    G4PhysicsFreeVector* fVectorLiSi;     
    G4PhysicsFreeVector* fVectorLiAu;
    G4PhysicsFreeVector* fVectorLiHfO2;   
    G4PhysicsFreeVector* fVectorLiSiO2;                
    

    G4PhysicsFreeVector* FillStopVector(G4String filename, G4int atomic_no);    
*/
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
