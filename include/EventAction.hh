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
#include "G4Physics2DVector.hh"
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
    void abs1STP()				{abs1StepPrim++ ;};
    void abs2STP()				{abs2StepPrim++ ;};
    void abs3STP()				{abs3StepPrim++ ;};
    void abs4STP()				{abs4StepPrim++ ;};        
    void AddTrakLenSec(G4double length)  {TrakLenSec  += length;};

    void AddTrueTrakLen(G4double trueLength) {TrueTrakLen += trueLength;};
    void AddProjTrakLen(G4double projLength) {ProjTrakLen += projLength;};
    
	// functions for RBS evaluations
    G4double RecoilEnergy(G4double E, G4double angle, G4double M1, G4double M2);	//evaluates recoiled energy
    G4double KinematicFactor(G4double angle, G4double M1, G4double M2);
    G4double RandomEnLoss(G4double E, G4double dedx, G4double position, G4double angle);	//calculates the energy lost when particle reaches outside volume after the scattering event
    G4double CalcDiffRuthXsec(G4double E, G4double M1, G4double M2, G4double angle,G4double Z1, G4double Z2);	// calculates the differential Rutherford xsec in laboratory system
    //G4double CalcRBSYield(G4double xsec, G4double dist, G4double solidAngle, G4double atomDens);
    G4double CalcRBSYield(G4double xsec, G4double dist, G4double solidAngle, G4double atomDens, G4double inc_angle);
	//function to calculate Bohr energy straggling
    G4double CalcBohrStrag(G4double Z1, G4double Z2, G4double atomDens, G4double dist);
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
    	// function for Total RBS yield, combining other functions into single one
    G4double CalculateTotalRBSYield(G4double energy, G4double M1, G4double M2, G4double Z1, G4double Z2, G4double angle, G4double dist,G4double solidAngle, G4double xsecmod, G4double atomDensity,G4double inc_angle);
    	// Total energy straggling, both electronic and nuclear included
    G4double CalculateTotalBohrStraggling(G4double energy, G4ParticleDefinition* particle, G4Material* mat, G4double distance);
    // based on yang theory https://www.sciencedirect.com/science/article/abs/pii/0168583X9195454L
    //G4double CalculateTotalStraggling(G4double energy, G4ParticleDefinition* particle, G4Material* mat, G4double distance, G4double q);
    
    
    	//https://www.sciencedirect.com/science/article/abs/pii/0168583X9195454L
    //G4double CalcEnergyStragglingFactor(G4double energy, G4double Z1, G4double Z2,G4double mass);
    	// dead layer energy loss
    G4double CalculateDeadLayerEffect(G4double energy, const G4Material* dead_mat, G4double thickness,G4ParticleDefinition* particle);
    	//https://www.sciencedirect.com/science/article/pii/S0168583X97006642
    G4double CalcDetectorFWHM(G4double energy, G4double Z1);
    
    G4double CalcAndersenScreening(G4double energy_cm, G4double angle_cm, G4double Z1, G4double Z2);
    	// from SIMNRA user's guide
    G4double CalcNuclEnStraggling(G4double Z1, G4double Z2, G4double M1, G4double M2, G4double atdens,G4double distance);
    
    G4double CalcScreening_TF(G4double Z1, G4double Z2);
    G4double CalcScreening_ZBL(G4double Z1, G4double Z2);
    
    void SetFreeSurf(G4int a) 	{free_surf = a;}
    G4int GetFreeSurf()		{return free_surf;}
    
    // Multiple scattering functions
    G4double CalcScallingFactorS(G4double a, G4double b, G4double v);
    G4double CalcScallingFactorNu(G4double tau);
    G4double CalcScallingFactorDeltaPhi(G4double tau);
    G4double CalcDeltaF_K(G4double M1, G4double energy, G4double angle, G4double M2, G4double kinematic_factor);
    G4double CalcScallingFactorMiu(G4double energy, G4double Z1, G4double Z2,G4double screening_radius);   
    //G4double CalcVectorAB_out(G4String vec, G4double delta_Fk, G4double C_out, G4double beta, G4double K_out);

    // rewrite of previous function to only account on the way out
    G4double CalculateMS_spread_out(G4double E_mid, G4double E_recoil,G4ParticleDefinition* fParticle, G4Material* mat, G4double distance_out, G4double exit_angle, G4double scat_angle,G4int mat_no,G4double Mm); 
    	
    	
    G4double CalculateMSandFillArrays(G4ParticleDefinition* fParticle, G4double exit_angle, G4double scattering_angle);   
    G4double CalcTauFromPhi(G4double delta_phi);    
    	
    G4double Fill_MS_Vector(G4double primary_energy);
    
    
    void SetShapeFactorM2(G4double m)		{shape_factor_m2 = m;}
    G4double GetShapeFactorM2()			{return shape_factor_m2;}
    
    void SetWidthFactorF2(G4double w)		{width_factor_w2 = w;}
    G4double GetWidthFactorF2()			{return width_factor_w2;}
	
	//G4double FindFWHMandM(G4double m_2, G4double f_2);
	void FindFWHMandM(G4double m_2, G4double f_2);
	
	void SetConvM(G4double a)			{conv_shape_factor = a;}
	G4double GetConvM()				{return conv_shape_factor;}
	
	void SetConvFWHM(G4double a)			{conv_fwhm = a;}
	G4double GetConvFWHM()				{return conv_fwhm;}
	
	void SetLastDphi(G4double a)			{last_dphi = a;}
	G4double GetLastDphi()				{return last_dphi;}
	
	void SetLastKout(G4double a)			{last_K_out = a;}
	G4double GetLastKout()				{return last_K_out;}	
	
	void SetLastTau(G4double a)			{last_tau = a;}
	G4double GetLastTau()				{return last_tau;}
	
	void SetLastMiuUp(G4double a)			{last_miu_up = a;}
	G4double GetLastMiuUp()			{return last_miu_up;}

	void SetLastMiuDown(G4double a)		{last_miu_down = a;}
	G4double GetLastMiuDown()			{return last_miu_down;}	
	
	void SetLastCout(G4double a)			{last_C_out = a;}
	G4double GetLastCout()				{return last_C_out;}
	
	void SetDepthLimit(G4double a)			{depth_lim = a;}
	G4double GetDepthLimit()			{return depth_lim;}
	
	void ResetMSParameters();
	void UpdateMSParameters();
	
	void ResetMSFactors();
	
	G4Physics2DVector* Get2DRTRVector(G4String element, G4int Z);
	G4double Get2DRTRValue(G4double energy, G4String elname, G4double angle);
	
	
	// distribution functions
	G4double CalcErrorFunction(G4double a);
	G4double GenerateGaussian(G4double x, G4double y, G4double sigma_sq);
    G4double CalcPearsonVII(G4double value/*, G4double average*/, G4double deviation, G4double shape);
    G4double CalcPearsonVII_new(G4double value/*, G4double average*/, G4double deviation, G4double shape);
    G4double CalcNewGauss(G4double value, G4double average, G4double deviation);	
    	
    G4double* CalcEnergyLeft(G4double depth, G4double exit_angle, G4double energy, G4double elementZ, G4double steps);
    	
    G4double FillTauPhiVector();
    void Construct2DVectors();
    	
    G4Physics2DVector* ReadMSFiles(G4String name);
    	
    	
  private:
  
    G4double max_permisible_energy;
  
    G4double last_dphi;
    G4double last_K_out;
    G4double last_tau;
    G4double last_miu_up;
    G4double last_miu_down;
    G4double last_C_out;
    G4double depth_lim;
    
  
    G4int nn_calc;
    G4double conv_shape_factor;
    G4double conv_fwhm;
    G4double shape_factor_m2;
    G4double width_factor_w2;
  
    G4double sum_part_phi_T_M;
    G4double sum_part_phi;
    G4double tot_sum_part_tau;
    G4double tau_part_hat;
    G4double last_delta_phi;
  
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
    G4int abs1StepPrim;
    G4int abs2StepPrim;
    G4int abs3StepPrim;
    G4int abs4StepPrim;
    
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

    //2D vectors of rtr values
    G4Physics2DVector* fVectorSi_total;
    G4Physics2DVector* fVectorO_total;   
    G4Physics2DVector* fVectorNa_total;       
    G4Physics2DVector* fVectorN_total;
    G4Physics2DVector* fVectorC_total;
    G4Physics2DVector* fVectorF_total;    
    G4Physics2DVector* fVectorB_total;      
    G4Physics2DVector* fVectorNi_total;  
    G4Physics2DVector* fVectorCu_total;     
    
    // 2D physics vectors for MS values
    G4Physics2DVector* MS_shape_surf;     
    G4Physics2DVector* MS_width_surf;           

    G4Physics2DVector* MS_shape_l1;     
    G4Physics2DVector* MS_width_l1;       

    G4Physics2DVector* MS_shape_l2;     
    G4Physics2DVector* MS_width_l2;       
    
    G4Physics2DVector* MS_shape_l3;     
    G4Physics2DVector* MS_width_l3;       

    G4Physics2DVector* MS_shape_l4;     
    G4Physics2DVector* MS_width_l4;            

    G4Physics2DVector* MS_shape_l12;     
    G4Physics2DVector* MS_width_l12;     
   
    G4Physics2DVector* MS_shape_l23;     
    G4Physics2DVector* MS_width_l23;                                                                     
    
    G4Physics2DVector* MS_shape_l34;     
    G4Physics2DVector* MS_width_l34;         
    
    G4Physics2DVector* MS_shape_l4e;     
    G4Physics2DVector* MS_width_l4e;         
    
    
    //2D vectors for FWHM and M retrieval of convoluted distros
    G4Physics2DVector* fVectorFWHM;//	= new G4Physics2DVector(8,7);
    G4Physics2DVector* fVectorM;//	= new G4Physics2DVector(8,7);
    
    //G4PhysicsFreeVector* fVectorSiXsec;
    
    G4PhysicsFreeVector* tauphi;
    G4PhysicsFreeVector* phitau;
    G4PhysicsFreeVector* taunu;

    
	const G4ElementVector* ElementVector;
	const G4Element* Element;
	const G4Isotope* Isotope;


	G4UAtomicDeexcitation* deExcitation = new G4UAtomicDeexcitation();
	//G4double *TOTALXsec;

public:

	G4Material* sample_material[5];
	
	G4int* NoOfElements 		= new G4int [5];
	G4double** Znumb		= new G4double* [5];	//[4]
	G4double** Mnumb		= new G4double* [5];	//[4]
	G4double** Adens		= new G4double* [5];	// [4]
	
	G4double* layer_thickness	= new G4double [4];
	G4ThreeVector* layer_pos	= new G4ThreeVector[4];
	G4double* layer_pos_z		= new G4double [4];
	G4double* neighbour_dist	= new G4double [3];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
