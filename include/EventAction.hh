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
/// \file EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 76293 2013-11-08 13:11:23Z gcosmo $
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"


//chan
#include <iostream>
#include <fstream>
#include "G4PhysicsFreeVector.hh"
#include "G4ElementVector.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"

#include "G4UAtomicDeexcitation.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class RunAction;
class DetectorConstruction; 
class PrimaryGeneratorAction;

class EventAction : public G4UserEventAction
{
  public:
    EventAction(DetectorConstruction*, PrimaryGeneratorAction*);
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void AddEdep (G4double Edep);
    void AddEflow(G4double Eflow);    

    void AddEnergy      (G4double edep)   {EnergyDeposit  += edep;};
    void AddNonIonEnergy(G4double enondep){NonIonEnergyDeposit  += enondep;};
    void AddTheta(G4double tet){theta+=tet;};	
    void AddTrakLenPrim(G4double length) {TrakLenPrim += length;};
    void CountStepsPrim()		 {nbStepsPrim++ ;};
    void absSTP()				{absStepPrim++ ;};
    void abs2STP()				{abs2StepPrim++ ;};
    void abs3STP()				{abs3StepPrim++ ;};
    void AddTrakLenSec(G4double length)  {TrakLenSec  += length;};

    void AddTrueTrakLen(G4double trueLength) {TrueTrakLen += trueLength;};
    void AddProjTrakLen(G4double projLength) {ProjTrakLen += projLength;};
    
    
    

	// functions for RBS evaluations
    G4double RecoilEnergy(G4double E, G4double angle, G4double M1, G4double M2);	//evaluates recoiled energy
    G4double RandomEnLoss(G4double E, G4double dedx, G4double position, G4double angle);	//calculates the energy lost when particle reaches outside volume after the scattering event
    G4double CalcDiffRuthXsec(G4double E, G4double M1, G4double M2, G4double angle,G4double Z1, G4double Z2);	// calculates the differential Rutherford xsec in laboratory system
    G4double CalcRBSYield(G4double xsec, G4double dist, G4double solidAngle, G4double atomDens);
    //G4double CalcRBSYield(G4double xsec, G4double dist, G4double solidAngle, G4double atomDens, G4double angle);
	//function to calculate Bohr energy straggling
    G4double CalcBohrStrag(G4double Z1, G4double Z2, G4double atomDens, G4double dist);
	//function to calculate modified rutherford scattering xsec due to screening corrections, see (E. Huttel, et.al., Screening corrections to the rutherford cross section, 1985, NIMB 12)
    //G4double CalcRuthXsecMod(G4double xsec, G4double Z1, G4double Z2, G4double energy, G4double angle);
	// see https://www.sciencedirect.com/science/article/abs/pii/0168583X85900503?via%3Dihub
    G4double CalcRuthXsecMod(G4double xsec, G4double Z1, G4double Z2, G4double energy);
	//function to calculate scattering angle in the CM reference frame
    G4double CalcAngleCMFrame(G4double angle, G4double M1, G4double M2);
	//function to calculate energy in the CM reference frame
    G4double CalcEnergyCMFrame(G4double energy, G4double M1, G4double M2);	       
	// RBS xsec in the CM reference fram
    G4double CalcDiffRuthXsecCM(G4double E, G4double angleCM ,G4double Z1, G4double Z2); 
	// RBS xsec in the Lab frame from the CM reference fram
    G4double CalcDiffRuthXsecLAB(G4double M1, G4double M2, G4double angle, G4double xsection);
	// integrates energy loss on the "particle way out" 
    G4double CalcTotEnLoss(G4double E, G4double distance, G4int steps, G4ParticleDefinition* fParticle, G4Material* mat);
    
	// calculates stopping number
    G4double CalcStoppingNumb(G4double energy, G4double mass, G4double ioniz); 
	// calculates chi for determination of which straggling regime to choose
    G4double CalcChi(G4double energy, G4double mass, G4double Z2);
    G4PhysicsFreeVector* FillRTRVector(G4String filename);
    
    	// function for Total RBS yield, combining other functions into single one
    G4double CalculateTotalRBSYield(G4double energy, G4double M1, G4double M2, G4double Z1, G4double Z2, G4double angle, G4double dist,G4double solidAngle, G4double xsecmod, G4double atomDensity);

    G4double GetRTRValue(G4double energy, G4String name);
    //G4double CalculateTotalBohrStraggling(G4double energy, G4double A1_PDG, G4double Z1, G4double Z2,G4double atomdensity, G4double distance, G4double ionisation);
    
    G4double CalculateTotalBohrStraggling(G4double energy, G4ParticleDefinition* particle, G4Material* mat, G4double distance);
    	//https://www.sciencedirect.com/science/article/abs/pii/0168583X9195454L
    G4double CalcHeavyIonStragglingFactor(G4double energy, G4double Z1, G4double Z2);
    
    G4double CalculateDeadLayerEffect(G4double energy, const G4Material* dead_mat, G4double thickness,G4ParticleDefinition* particle);
    
    void SetFreeSurf(G4int a) 	{free_surf = a;}
    G4int GetFreeSurf()		{return free_surf;}
    
    
    

    G4double CalcErrorFunction(G4double a);
    G4double GenerateGaussian(G4double x, G4double y, G4double sigma_sq);


  private:
    DetectorConstruction* detector;
    PrimaryGeneratorAction* fPrimary;


    G4double fTotalEnergyDeposit;
    G4double fTotalEnergyFlow;   

    G4double EnergyDeposit;
    G4double NonIonEnergyDeposit;
    G4double theta;	

    G4double TrakLenPrim;
    G4int nbStepsPrim;
    G4double TrakLenSec;

    G4int absStepPrim;
    G4int abs2StepPrim;
    G4int abs3StepPrim;

    
    G4double TrueTrakLen;
    G4double ProjTrakLen;

    G4double En, Angle, Ma1, Ma2;
    G4double ang2;

    int hit;
    
    G4int free_surf;




    //chan bullshit
    G4int sdht_ID;
    G4int sdct_ID;
	// for second material
    G4int sdxt_ID;
    
    
    G4int sd0_ID,sd1_ID,sd2_ID,sd3_ID,sd4_ID;


	// RTR databases for elements
    G4PhysicsFreeVector* fVectorSi;
    G4PhysicsFreeVector* fVectorO;
    G4PhysicsFreeVector* fVectorB;
    G4PhysicsFreeVector* fVectorC;
    G4PhysicsFreeVector* fVectorF;
    G4PhysicsFreeVector* fVectorN;
    G4PhysicsFreeVector* fVectorNa;
    G4PhysicsFreeVector* fVectorMg;
    G4PhysicsFreeVector* fVectorP;
    G4PhysicsFreeVector* fVectorS;
    G4PhysicsFreeVector* fVectorK;
    G4PhysicsFreeVector* fVectorCa;
    G4PhysicsFreeVector* fVectorTi;
    G4PhysicsFreeVector* fVectorFe;
    G4PhysicsFreeVector* fVectorCr;
    G4PhysicsFreeVector* fVectorNi;

	const G4ElementVector* ElementVector;
	const G4Element* Element;
	const G4Isotope* Isotope;


	G4UAtomicDeexcitation* deExcitation = new G4UAtomicDeexcitation();
	//G4double *TOTALXsec;
public:

	G4Material* sample_material[5];
	/*G4int NoOfElements[5];
	G4double Znumb[5][4];
	G4double Mnumb[5][4];
	G4double Adens[5][4];
	*/
	
	//G4Material* sample_material	= new G4Material* [5];
	G4int* NoOfElements 		= new G4int [5];
	G4double** Znumb		= new G4double* [5];	//[4]
	G4double** Mnumb		= new G4double* [5];	//[4]
	G4double** Adens		= new G4double* [5];	// [4]
	
	G4double* layer_thickness	= new G4double [4];
	G4ThreeVector* layer_pos	= new G4ThreeVector[4];
	G4double* layer_pos_z		= new G4double [4];
	G4double* neighbour_dist	= new G4double [3];
	

	//G4double ionization[5];
	
	//G4double KinEnergy[5];

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
