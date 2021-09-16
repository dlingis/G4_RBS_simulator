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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 76293 2013-11-08 13:11:23Z gcosmo $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "RunAction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"

// is chanell

#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "CrystalDetectorHit.hh"

#include "Analysis.hh"

// new addition

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Track.hh"
#include "DetectorConstruction.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4EmCalculator.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"

#include <cmath>

#include "G4NistManager.hh"


//test fluctuations
#include "G4IonYangFluctuationModel.hh"
#include "G4IonChuFluctuationModel.hh"

#include "G4hIonEffChargeSquare.hh"


#include "G4Physics2DVector.hh"
#include "G4Pow.hh"


//using std::vector;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* DE, PrimaryGeneratorAction* PRA )
:G4UserEventAction(),detector(DE),fPrimary(PRA),
 fTotalEnergyDeposit(0.), fTotalEnergyFlow(0.), EnergyDeposit(0.), NonIonEnergyDeposit(0.),
  theta(0), TrakLenPrim(0.), nbStepsPrim(0.), TrakLenSec(0.), absStepPrim(0.), abs2StepPrim(0.), TrueTrakLen(0.), ProjTrakLen(0.), sdht_ID(-1),sdct_ID(-1),sdxt_ID(-1), // paskutiniai is chan
  sd0_ID(-1),sd1_ID(-1),sd2_ID(-1),sd3_ID(-1),sd4_ID(-1)
{


	//G4UAtomicDeexcitation* deExcitation = new G4UAtomicDeexcitation();
	//deExcitation->InitialiseAtomicDeexcitation();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
		for (int i = 0; i<5; i++)
		{
			delete[] Znumb[i];
			delete[] Mnumb[i];
			delete[] Adens[i];
		}
	//delete[] sample_material;
	delete[] NoOfElements;
	delete[] Znumb;
	delete[] Mnumb;
	delete[] Adens;
	delete[] layer_thickness;
	delete[] layer_pos;
	delete[] layer_pos_z;
	delete[] neighbour_dist;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
	fTotalEnergyDeposit = 0.;
  	fTotalEnergyFlow = 0.;
 	// initialisation per event
 	EnergyDeposit  = 0.;
 	NonIonEnergyDeposit=0.;
 	theta=0.;
 	TrakLenPrim=0.;
 	TrakLenSec=0.;
 	nbStepsPrim=0.;
 	TrueTrakLen=0.;
	ProjTrakLen=0.;
 	hit = 0.;
 	absStepPrim = 0;
 	abs1StepPrim = 0;
 	abs2StepPrim = 0; 
 	abs3StepPrim = 0;
 	abs4StepPrim = 0;

 	SetFreeSurf(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
	// start of functions relevant to RBS spectrum
G4double EventAction::RecoilEnergy(G4double E, G4double angle, G4double M1, G4double M2)
{
	//kinematic factor
	// M1 - incident particle, M2 - target atom, E - incident energy, angle - scattering angle
	G4double k;
	G4double square = std::sqrt(std::pow(M2/M1,2.)-std::pow(sin(angle),2.));
	G4double M1co = cos(angle);
	G4double denominator = 1+(M2/M1);
	k = std::pow((M1co + square)/denominator, 2.);
	return  E*k;
}

// kinematic factor for multiple scattering evaluations
G4double EventAction::KinematicFactor(G4double angle, G4double M1, G4double M2)
{
	G4double k = 0.;
	G4double square = std::sqrt(std::pow(M2,2.)-std::pow(M1*sin(angle),2.));
	G4double M1co = M1*cos(angle);
	G4double denominator = M1+M2;
	k = std::pow((M1co + square)/denominator, 2.);
	return k;
}
G4double EventAction::CalcDiffRuthXsec(G4double E, G4double M1, G4double M2, G4double angle, G4double Z1, G4double Z2)
{
	G4double cose = cos(angle);
	G4double sine = sin(angle);
	G4double epsilon = 55.26349406/1000; // in units of e^2/(MeV*fm)
	G4double FirstTermDenom = 8*pi*epsilon*E;
	G4double FirstTerm = std::pow((Z1*Z2)/(FirstTermDenom),2.);
	G4double SecondTerm = 1/(std::pow(sine,4.));
	G4double ThirdTerm = std::pow(M2*cose+std::pow(std::pow(M2,2)-std::pow(M1*sine,2),0.5),2);
	G4double ThirdTermDenom = M2*std::pow(std::pow(M2,2)-std::pow(M1*sine,2),0.5);
	G4double fullTerm = FirstTerm*SecondTerm*ThirdTerm/ThirdTermDenom;
	return fullTerm*10;//fm^2 conversion to milibarns
}

G4double EventAction::CalcRBSYield(G4double xsec, G4double dist, G4double solidAngle, G4double atomDens, G4double inc_angle)
{
	G4double crossSect = xsec * 1e-27; 	// mbarn to cm2 = 1e−27
	G4double z = dist/cm;
	G4double atom_thickness = atomDens*z;
	G4double yield = (solidAngle*crossSect*atom_thickness/cos(inc_angle));
	return yield;
}

	//energy straggling in MeV^2
G4double EventAction::CalcBohrStrag(G4double Z1, G4double Z2, G4double atomDens, G4double dist)
{
	G4double elmcharge_squared = 1.4399764;	// in units of MeV*fm
	G4double e2 = std::pow(elmcharge_squared,2.)*1e-26;	// in units of MeV^2*cm^2
	G4double fourpi = 4*3.14159265358979;
	G4double fBohr = fourpi * e2;
	G4double Bohr_t = Z1*Z1*Z2*fBohr*atomDens*(dist/cm);
	return Bohr_t;
}

G4double EventAction::CalcRuthXsecMod(G4double xsec, G4double Z1, G4double Z2, G4double energy)
{
	G4double a = 0.049*Z1*std::pow(Z2,4/3);
	G4double b = a/(energy/keV);
	G4double c = 1-b;
	G4double result = xsec*c;
	return result;
}
	//from http://atlas.physics.arizona.edu/~shupe/Indep_Studies_2015/Notes_Goethe_Univ/L4_Scattering_Basic.pdf
G4double EventAction::CalcAngleCMFrame(G4double angle, G4double M1, G4double M2)
{
	// old version

	G4double CMang=0., x= 0., y=0.;
	x = (M1/M2)*sin(angle);
	y = asin(x);
	CMang = y+angle;

	return CMang;
}
	//energy in the CM reference frame
G4double EventAction::CalcEnergyCMFrame(G4double energy, G4double M1, G4double M2)
{
	G4double en=0.;
	en = energy*M2/(M1+M2);
	return en;
}

	//xsec in the CM reference system
G4double EventAction::CalcDiffRuthXsecCM(G4double E, G4double angleCM, G4double Z1, G4double Z2)
{
	G4double sine = sin(angleCM/2);
	G4double epsilon = 55.26349406/1000; // in units of e^2/(MeV*fm)
	G4double FirstTermDenom = 16*pi*epsilon*E;
	G4double FirstTerm = std::pow((Z1*Z2)/(FirstTermDenom),2.);
	G4double SecondTerm = 1/(std::pow(sine,4.));
	G4double fullTerm = FirstTerm*SecondTerm;
	return fullTerm*10;
}
	//xsec in the Lab reference system from the CM reference system
G4double EventAction::CalcDiffRuthXsecLAB(G4double M1, G4double M2, G4double angle, G4double xsection)
{
	G4double ratio = M1/M2;
	G4double multiplier = std::pow(1+(ratio*ratio)+(2*ratio*cos(angle)),3./2.)/(1+(ratio*cos(angle)));
	//G4cout << " multiplier " << multiplier << G4endl;
	G4double labXSEC = multiplier * xsection;
	return labXSEC;
}

	// integrates energy loss in distance
G4double EventAction::CalcTotEnLoss(G4double E, G4double distance, G4int steps, G4ParticleDefinition* fParticle, G4Material* mat)
{
	G4EmCalculator emCalculator;
	// step distance
	G4double stp = (distance/steps)/cm;
	// stopping power more or less accurate compared to pstar
	for (int i=1; i<=steps; i++)
	{
		G4double stop = emCalculator.ComputeTotalDEDX(E,fParticle,mat)/(MeV/cm);
		//G4double stop = emCalculator.GetDEDX(E,fParticle,mat)/(MeV/cm);
		E -= stop*stp;
		if (E <= 0.01/MeV)
		{
			break;
		}
	}
	return E;
}

	// end of functions relevant to RBS spectrum
void EventAction::AddEdep(G4double Edep)
{
  fTotalEnergyDeposit += Edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddEflow(G4double Eflow)
{
  fTotalEnergyFlow += Eflow;
}

void EventAction::EndOfEventAction(const G4Event* evt)
{

  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  run->AddEdep (fTotalEnergyDeposit);
  run->AddEflow(fTotalEnergyFlow);
  run->AddTrakLenPrim(TrakLenPrim);
  run->AddNumberOfSteps(nbStepsPrim);
  run->AddTheta(theta);
  run->AddTrakLenSec(TrakLenSec);
  run->AddTrueRange(TrueTrakLen);
  run->AddProjRange(ProjTrakLen);
  
  run->absSTP(absStepPrim);
  run->abs1STP(abs1StepPrim);
  run->abs2STP(abs2StepPrim);
  run->abs3STP(abs3StepPrim);
  run->abs4STP(abs4StepPrim);    
  
  G4AnalysisManager::Instance()->FillH1(1,fTotalEnergyDeposit);
  G4AnalysisManager::Instance()->FillH1(3,fTotalEnergyFlow);  
    
  G4int rbs_eval = detector->GetRBSCalc();
  
  // ******************************

  	for (int i = 0; i<5; i++) {
		Znumb[i] = new G4double[4];
		Mnumb[i] = new G4double[4];
		Adens[i] = new G4double[4];
	}

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	//detector energy resolution
	G4double det_FWHM = detector->GetDetectorResolution();
	G4int gauss_counter = detector->GetGaussCounter();

	G4ParticleDefinition*  fParticle = fPrimary->GetParticleGun()->GetParticleDefinition();
    G4double A1 = fPrimary->GetParticleGun()->GetParticleDefinition()->GetAtomicMass();
	//G4double A1_PDG = fPrimary->GetParticleGun()->GetParticleDefinition()->GetPDGMass();
    G4double Z1 = fPrimary->GetParticleGun()->GetParticleDefinition()->GetAtomicNumber();
    G4double primary_energy = fPrimary->GetParticleGun()->GetParticleEnergy();
    run->PrimaryEnergy(primary_energy);
	//**************************************************************************
	// Incidence angle
	G4double angle_of_incidence = 0.;
	angle_of_incidence = atan(fPrimary->GetParticleGun()->GetParticleMomentumDirection().x()/fPrimary->GetParticleGun()->GetParticleMomentumDirection().z());
	run->Inc_angle(angle_of_incidence);
	

	//scattering angle
	Angle 				= detector->GetRBSAngle();
	ang2 				= 3.14159265358979; //pi radians

	//G4double exit_angle 		= pi - Angle - angle_of_incidence;
	G4double exit_angle 		= pi - Angle;
	
  	G4EmCalculator emCalculator;


    free_surf = 0.;
    // detector dead layer


    G4double RBS_norm_dist 	= 0.1*nm;

    G4double step_number 		= detector->GetEnLossStep();
    G4double dead_thickness 	= detector->GetDeadLayerThickness();
    G4String dead_material_name = detector->GetDeadLayer();
	G4double solidAngle 		= detector->GetSolidAngle();
   	G4int fwhm_calc 			= detector->GetCalcFWHM();
   	G4int ms_calc				= detector->GetMSCalc(); 	
   	G4double ROI_region			= detector->GetRBSROImin();
   	G4int use_const_angle 		= detector->GetConstAngle();
    		
    // MS relevant parameters
    last_dphi = 0.;
    last_K_out = 0.;
    last_tau = 0.;
    last_miu_up = 0.;
    last_miu_down = 0.;

    const G4Material* dead_material = G4NistManager::Instance()->FindOrBuildMaterial(dead_material_name);
    G4Material* d_mat = G4NistManager::Instance()->FindOrBuildMaterial(dead_material_name);


	for (int i=0;i<5;i++) {
		sample_material[i] 	= detector->GetMaterialM(i);
		NoOfElements[i]    	= sample_material[i]->GetNumberOfElements();
	}
	G4int max_M = 0;

	for (int i=0;i<5;i++) { 
		for (int j = 0; j<NoOfElements[i]; j++) {
			Znumb[i][j] 				= sample_material[i]->GetElement(j)->GetZ();
			Mnumb[i][j] 				= sample_material[i]->GetElement(j)->GetA()/(g/mole);
			if(Mnumb[i][j] > max_M)
				max_M = Mnumb[i][j];
			const G4double *atomDensVector	= sample_material[i]->GetVecNbOfAtomsPerVolume();
			Adens[i][j] 				= atomDensVector[j]/(1/cm3);
		}
	}

	max_permisible_energy = 1.05*RecoilEnergy(primary_energy, Angle, A1, max_M);

	// RTR vector fillup
	if (detector->GetSigmaCalc() == 1 && rbs_eval == 1) {
  		G4String mat_name;
  		for(int i=0;i<5;i++) {
  			for(int j=0;j<NoOfElements[i]; j++)	{
				mat_name = sample_material[i]->GetElement(j)->GetName();
			  	if 	(mat_name == "Si" && fVectorSi_total == 0)	{ 
  				  	fVectorSi_total = Get2DRTRVector("Si", A1);
  	  				fVectorSi_total->SetBicubicInterpolation(true);
  				} else if (mat_name == "O" && fVectorO_total == 0) { 
  				  	fVectorO_total = Get2DRTRVector("O",A1);
	  				fVectorO_total->SetBicubicInterpolation(true);
  				} else if (mat_name == "B"  && fVectorB_total == 0)	{ 
  				  	fVectorB_total = Get2DRTRVector("B",A1);
	  				fVectorB_total->SetBicubicInterpolation(true);  					
  				} else if (mat_name == "C"  && fVectorC_total == 0) { 
  				  	fVectorC_total = Get2DRTRVector("C",A1);
	  				fVectorC_total->SetBicubicInterpolation(true);  	  						
  				} else if (mat_name == "F" && fVectorF_total == 0)	{
  				  	fVectorF_total = Get2DRTRVector("F",A1);
	  				fVectorF_total->SetBicubicInterpolation(true);  	  					
  				} else if (mat_name == "N"  && fVectorN_total == 0)	{
  				  	fVectorN_total = Get2DRTRVector("N",A1);
	  				fVectorN_total->SetBicubicInterpolation(true);  	  					
  				} else if (mat_name == "Na" && fVectorNa_total == 0) {
  				  	fVectorNa_total = Get2DRTRVector("Na",A1);
	  				fVectorNa_total->SetBicubicInterpolation(true);    					
  				} else if (mat_name == "Ni" && fVectorNi_total == 0) { 
  				  	fVectorNi_total = Get2DRTRVector("Ni",A1);
	  				fVectorNi_total->SetBicubicInterpolation(true);  	  		
  				} else if (mat_name == "Cu" && fVectorCu_total == 0) { 
  				  	fVectorCu_total = Get2DRTRVector("Cu",A1);
	  				fVectorCu_total->SetBicubicInterpolation(true);  	  		
  				}  				
  			} // end of cycle throu elements
  		} // end of cycle throu materials
  	}	// end of sigma calc


  	// sensitive detectors
  	G4SDManager* SDman = G4SDManager::GetSDMpointer();
  	G4String sdName;
	if(sd0_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector",0))
		{
            		sd0_ID = SDman->GetCollectionID(sdName="crystaldetector/collection");
		}
   	 }
   	 if(sd1_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector2",0))
		{
            		sd1_ID = SDman->GetCollectionID(sdName="crystaldetector2/collection");
		}
   	 }
   	 if(sd2_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector3",0))
		{
            		sd2_ID = SDman->GetCollectionID(sdName="crystaldetector3/collection");
		}
   	 }
   	 if(sd3_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector4",0))
		{
            		sd3_ID = SDman->GetCollectionID(sdName="crystaldetector4/collection");
		}
   	 }
   	 if(sd4_ID == -1) {
		if(SDman->FindSensitiveDetector(sdName="crystaldetector5",0))
		{
            		sd4_ID = SDman->GetCollectionID(sdName="crystaldetector5/collection");
		}
   	 }


    	CrystalDetectorHitsCollection* sd0 = 0;
    	G4HCofThisEvent *hx0 = evt->GetHCofThisEvent();

    	if(hx0){
        	if(sd0_ID != -1){
            		G4VHitsCollection* aHCSD0 = hx0->GetHC(sd0_ID);
            		sd0 = (CrystalDetectorHitsCollection*)(aHCSD0);
        		}
    		}
    	CrystalDetectorHitsCollection* sd1 = 0;
    	G4HCofThisEvent *hx1 = evt->GetHCofThisEvent();

    	if(hx1){
        	if(sd1_ID != -1){
            		G4VHitsCollection* aHCSD1 = hx1->GetHC(sd1_ID);
            		sd1 = (CrystalDetectorHitsCollection*)(aHCSD1);
        		}
    		}
    	CrystalDetectorHitsCollection* sd2 = 0;
    	G4HCofThisEvent *hx2 = evt->GetHCofThisEvent();

    	if(hx2){
        	if(sd2_ID != -1){
            		G4VHitsCollection* aHCSD2 = hx2->GetHC(sd2_ID);
            		sd2 = (CrystalDetectorHitsCollection*)(aHCSD2);
        		}
    		}
    	CrystalDetectorHitsCollection* sd3 = 0;
    	G4HCofThisEvent *hx3 = evt->GetHCofThisEvent();

    	if(hx3){
        	if(sd3_ID != -1){
            		G4VHitsCollection* aHCSD3 = hx3->GetHC(sd3_ID);
            		sd3 = (CrystalDetectorHitsCollection*)(aHCSD3);
        		}
    		}
    	CrystalDetectorHitsCollection* sd4 = 0;
    	G4HCofThisEvent *hx4 = evt->GetHCofThisEvent();

    	if(hx4){
        	if(sd4_ID != -1){
            		G4VHitsCollection* aHCSD4 = hx4->GetHC(sd4_ID);
            		sd4 = (CrystalDetectorHitsCollection*)(aHCSD4);
        		}
    		}

   for (int i=0; i<4; i++) {
  	layer_thickness[i] 	= detector->GetLength(i+1);
  	layer_pos[i]		= detector->GetPosition(i);
  	layer_pos_z[i]		= layer_pos[i].z();
   }
   for (int i=0; i< 3; i++) {
  	neighbour_dist[i]	=(layer_pos_z[i+1]-layer_thickness[i+1]/2)-(layer_pos_z[i]+layer_thickness[i]/2);
   	if (neighbour_dist[i]/um < 0.0001) 
   		neighbour_dist[i] =0.;
   }

	G4double surf_distance = detector->GetLength(0)/2+layer_pos_z[0]-(detector->GetLength(1)/2);
	if (surf_distance/um < 0.0001) { surf_distance = 0.;}
	
	G4double steps;
	G4ThreeVector position;
	G4double sample_energy =0.;
	G4double xsecRTR;
	G4double RBS_yield =0.;
	G4double energy_left = 0.;
	G4double depth =0.;

	G4double trueWorldPosition = 0.;
	G4double bohr_straggling = 0.;
	free_surf = 0;

	G4double tot_step =0.;
	G4double sigma_det = 10.*keV;	
	G4double tot_sigma;	
	
	G4double MS_M_conv = 0.;
	G4double MS_FWHM_conv = 0.;
	
	G4double FWHM_G(1.), f2_ratio(1.), tot_shape(1.), tot_fwhm(1.);
	
	G4int GK = gauss_counter;
	G4double santykis;
	G4double Gauss_sum	= 0;
	G4double Pearson_sum   = 0;
	G4double Pearson_value;
	G4double Gauss_value;
	G4double mod_P_value;
	G4double mod_G_value;
	G4double Gauss_en;	
	
	if(rbs_eval == 1 && ms_calc == 1) {
		if(fVectorFWHM == 0 || fVectorM == 0) {
			fVectorFWHM = new G4Physics2DVector;
			fVectorM = new G4Physics2DVector;
			std::ifstream vFileInF("MS/MS_FWHM.DAT");
			if(!vFileInF.is_open()) { 
				G4cout << " No MS_FWHM.DAT file " << G4endl;
				exit(1);
			}	
			fVectorFWHM->Retrieve(vFileInF);
			fVectorFWHM->SetBicubicInterpolation(true);
			
			std::ifstream vFileInM("MS/MS_M.DAT");
			if(!vFileInM.is_open())  { 
				G4cout << " No MS_M.DAT file " << G4endl;
				exit(1);
			}	
			G4cout << " Reading values from MS/MS_M.dat and MS/MS_FWHM.DAT for interpolation " << G4endl;
			fVectorM->Retrieve(vFileInM);
			fVectorM->SetBicubicInterpolation(true);
		}
		if(tauphi == 0 || phitau == 0 || taunu == 0)
			FillTauPhiVector();
		// test to see whether vectors are set up
		if(MS_shape_surf == 0 || MS_width_l1 == 0)
			Fill_MS_Vector(primary_energy);
	}

	// mother volume sensitive detector
    if(sd0)
    {
        int n_hit_sd = sd0->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++) {
           	CrystalDetectorHit* aHit = (*sd0)[i1];
           	steps = aHit->GetStep();
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();

	    	tot_step +=steps;

	    	trueWorldPosition = position.z()+detector->GetLength(0)/2;
			G4ThreeVector momDir = aHit->GetWorldMomentumDirection();
		
			G4double rand_angle_x = atan(momDir.x()/momDir.z());	
			//G4double rand_angle_y = atan(momDir.y()/momDir.z());

			//G4cout << " angles = " << rand_angle_x/degree << " " << rand_angle_y/degree << G4endl;
			G4double curr_angle = rand_angle_x;
			G4double diff_angle;
		
			if(use_const_angle == 1)
				diff_angle = Angle;
			else
				diff_angle = Angle-curr_angle;

    		for(int i=0;i<NoOfElements[0];i++) {
    			G4double RecEn = RecoilEnergy(sample_energy,diff_angle,A1,Mnumb[0][i]);
    			G4String el_name = sample_material[0]->GetElement(i)->GetName();
    				
    			if(detector->GetSigmaCalc() == 1)
					xsecRTR = Get2DRTRValue(sample_energy, el_name,diff_angle);
    			else 
    				xsecRTR = 1.;
	
    			if (RecEn/MeV > 0.1) {
    				RBS_yield = CalculateTotalRBSYield(sample_energy, A1, Mnumb[0][i], Z1, Znumb[0][i], diff_angle, RBS_norm_dist, solidAngle, xsecRTR, Adens[0][i],angle_of_incidence);
	    			G4double *main_l_parameters = CalcEnergyLeft(trueWorldPosition, exit_angle, RecEn,Znumb[0][i],step_number);	
    				energy_left 	= main_l_parameters[0];
    				bohr_straggling = main_l_parameters[1];
    				MS_M_conv 	= main_l_parameters[2];
					MS_FWHM_conv 	= main_l_parameters[3];
    				
					for(int h = 0; h<4; h++) {
						main_l_parameters[h] = 0;
					}	
					// filling of RBS histograms
					if (energy_left > ROI_region)
					{
					//G4cout << " DEPTH = " << trueWorldPosition/nm << " [nm] " << G4endl;
					//G4cout << " C_out " << last_C_out << G4endl;
					//G4cout << " Depth " << trueWorldPosition/nm << " Element " <<  sample_material[0]->GetElement(i)->GetName() << " dphi " << last_dphi*57.2957795 << " espread " << GetWidthFactorF2() << G4endl;
					//G4cout << " Dphi = " << GetLastDphi()*57.2957795 << " [degree] " << G4endl; 
					//G4cout << " Dphi " << GetLastDphi()*57.2957795 << " espread " << GetWidthFactorF2() << G4endl;
					//G4cout << " Energy spread = " << GetWidthFactorF2() << G4endl;
					
						// include detector dead layer energy loss
						energy_left = CalculateDeadLayerEffect(energy_left, dead_material, dead_thickness, fParticle);
						G4double dead_layer_straggling = CalculateTotalBohrStraggling(energy_left, fParticle, d_mat, dead_thickness);
						bohr_straggling += dead_layer_straggling;
	
						G4double rbs_depth = depth;// - detector->GetLength(0)/2;
	
						// rbs dist
						analysisManager->FillH1(21, rbs_depth, RBS_yield); // total
						analysisManager->FillH1(23, rbs_depth, RBS_yield); // substrate
						// el1
						if (i == 0)
							analysisManager->FillH1(26, rbs_depth, RBS_yield); //el1
						// el2
						if (i == 1)
							analysisManager->FillH1(27, rbs_depth, RBS_yield); // el2
	
						if (Z1 == 1)
							sigma_det = (det_FWHM/MeV)/2.355;
						if (Z1 == 3) {
							if (fwhm_calc == 1)
								sigma_det = (CalcDetectorFWHM(energy_left,Z1)/1000*MeV)/2.355;
							else if (fwhm_calc == 0)
								sigma_det = (det_FWHM/MeV)/2.355;
						}
						if (Z1 == 2)
							sigma_det = (det_FWHM/MeV)/2.355;
	
						if  (sigma_det/keV < 10.*keV)
							sigma_det = (10.*keV)/2.355;
	
						G4double sigma_det_sq = std::pow(sigma_det,2);
	
						tot_sigma = sigma_det_sq + bohr_straggling;
						if(ms_calc == 1) {
							//Gauss convolution with Pearson VII
							FWHM_G = sqrt(tot_sigma)*1000; // Gauss fwhm due to straggling and detector
							f2_ratio  = MS_FWHM_conv/FWHM_G;	// ratio of MS fwhm to Gauss fwhm
							FindFWHMandM(MS_M_conv, f2_ratio);// evaluation of convoluted fwhm and shape factor					
							tot_shape = GetConvM();		// shape factor of convoluted
							tot_fwhm = GetConvFWHM()*FWHM_G;	// fwhm of convoluted pearson vii distro
						}
	
						Gauss_sum	= 0;
						Pearson_sum   	= 0;
						
						for (int k=-GK; k<=GK; k++)	{
							santykis	= k;
							Gauss_en 	= energy_left*(1-(santykis/1000));
							if(ms_calc == 0) {
								Gauss_value 	= (GenerateGaussian(Gauss_en,energy_left,tot_sigma))/GenerateGaussian(0,0,tot_sigma);
								Gauss_sum 	+= Gauss_value;
							} else if(ms_calc == 1)	{
								Pearson_value	= CalcPearsonVII(k, tot_fwhm,tot_shape)/CalcPearsonVII(0,tot_fwhm,tot_shape);
								Pearson_sum 	+= Pearson_value;
							}
						}
						for (int k=-GK; k<=GK; k++)	{
							santykis	= k;					
							Gauss_en 	= energy_left*(1-(santykis/1000));
							if(ms_calc == 0) {
								Gauss_value 	= (GenerateGaussian(Gauss_en,energy_left,tot_sigma))/GenerateGaussian(0,0,tot_sigma);
								mod_G_value 	= Gauss_value/Gauss_sum;
							} else if(ms_calc == 1)	{
								Pearson_value 	= CalcPearsonVII(k,tot_fwhm,tot_shape)/CalcPearsonVII(0,tot_fwhm,tot_shape);
								mod_P_value 	= Pearson_value/Pearson_sum;
							}
						
							if (Gauss_en/MeV > ROI_region)
							{
								run->add_entry_reach(1);
								if(ms_calc == 1) {
									if(Gauss_en/MeV < max_permisible_energy) {									
										analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));
										analysisManager->FillH1(22, Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	//substrate								
										if (i == 0)
											analysisManager->FillH1(24,Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	//el1						
										else if (i == 1)
											analysisManager->FillH1(25,Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));
									}									
								} else if (ms_calc == 0) {
									analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	// total
									analysisManager->FillH1(22, Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	//substrate								
									if (i == 0)
										analysisManager->FillH1(24,Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	//el1						
									else if (i == 1)
										analysisManager->FillH1(25,Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));
								}	
							}
						}
					}	
				}		
	    	}
		}
		run->add_total_step(tot_step);
  }
	// end of mother sensitive detector

//======================================================
// LAYER 1
//======================================================

	// layer1 sensitive detector
    if(sd1){
    	tot_step =0.;
        int n_hit_sd = sd1->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++){
            CrystalDetectorHit* aHit = (*sd1)[i1];
            steps = aHit->GetStep();
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();

			tot_step +=steps;
			G4ThreeVector momDir = aHit->GetWorldMomentumDirection();
			G4double rand_angle_x = atan(momDir.x()/momDir.z());	
			G4double curr_angle = rand_angle_x;
			G4double diff_angle;
		
			if(use_const_angle == 1)
				diff_angle = Angle;
			else
				diff_angle = Angle-curr_angle;

    		G4double layer_to_surf = std::abs(-(detector->GetLength(0)/2)-layer_pos_z[0]+(detector->GetLength(1)/2));
    		if (layer_to_surf/um < 0.0001)
	    		layer_to_surf = 0.;
	   		// position in the sample
	    	trueWorldPosition = position.z()+detector->GetLength(0)/2;
	    	// position in the layer
	    	//G4double position_in_layer = trueWorldPosition - layer_to_surf;		
	    	G4double newWorldPosition = (position.z()+detector->GetLength(0)/2);	    	


	    	for(int i=0;i<NoOfElements[1];i++) {
				G4double RecEn = RecoilEnergy(sample_energy, diff_angle, A1, Mnumb[1][i]);
	   			G4String el_name = sample_material[1]->GetElement(i)->GetName();

    			if(detector->GetSigmaCalc() == 1)
					xsecRTR = Get2DRTRValue(sample_energy, el_name,diff_angle);
    			else 
    				xsecRTR = 1.;
    			if (RecEn/MeV > 0.1) {
	    			RBS_yield = CalculateTotalRBSYield(sample_energy, A1, Mnumb[1][i], Z1, Znumb[1][i], diff_angle, RBS_norm_dist, solidAngle, xsecRTR, Adens[1][i],angle_of_incidence);	
    				G4double *l1_parameters = CalcEnergyLeft(newWorldPosition, exit_angle, RecEn,Znumb[1][i],step_number);		
    				energy_left 	= l1_parameters[0];
    				bohr_straggling = l1_parameters[1];
    				MS_M_conv 	= l1_parameters[2];
					MS_FWHM_conv 	= l1_parameters[3];
				
   					for(int h = 0; h<4; h++) {
						l1_parameters[h] = 0;
					}
					// filling of RBS histograms
					if (energy_left > ROI_region) {
						// include detector dead layer energy loss
						energy_left = CalculateDeadLayerEffect(energy_left, dead_material,dead_thickness,fParticle);

						G4double dead_layer_straggling = CalculateTotalBohrStraggling(energy_left, fParticle, d_mat, dead_thickness);
						bohr_straggling += dead_layer_straggling;
					
						G4double rbs_depth = depth;
						// rbs dist
						analysisManager->FillH1(21, rbs_depth, RBS_yield); // total
						analysisManager->FillH1(29, rbs_depth, RBS_yield); // layer1
						// el1
						if (i == 0)
							analysisManager->FillH1(32, rbs_depth, RBS_yield); //el1
						// el2
						if (i == 1)
							analysisManager->FillH1(33, rbs_depth, RBS_yield); // el2
	
						if (Z1 == 1)
							sigma_det = (det_FWHM/MeV)/2.355;
						if (Z1 == 3) {
							if (fwhm_calc == 1)
								sigma_det = (CalcDetectorFWHM(energy_left,Z1)/1000*MeV)/2.355;
							else if (fwhm_calc == 0)
								sigma_det = (det_FWHM/MeV)/2.355;
						}
						if (Z1 == 2)
						{sigma_det = (det_FWHM/MeV)/2.355;}
					/*
					else if (Z1 == 2 && fwhm_calc == 1 && energy_left/MeV > 0.2)
						{
							G4double energy_for_cal 	= 3.0417*MeV;
							G4double test_calc = CalcDetectorFWHM(energy_for_cal,Z1)/1000*MeV;
							G4double exp_vs_theor = det_FWHM/test_calc;

							sigma_det = (CalcDetectorFWHM(energy_left,Z1)/1000*MeV)*exp_vs_theor;
							//G4cout << " det eff " << sigma_det/keV << " prie [MeV] " << energy_left/MeV << G4endl;

						}
					*/
					//G4cout << " sigma det " << sigma_det << G4endl;
					//G4cout << " en left " << energy_left/MeV << " det eff " << sigma_det/keV << G4endl;

						if (sigma_det/keV < 10.*keV)
							sigma_det = (10.*keV)/2.355;


						G4double sigma_det_sq = std::pow(sigma_det,2);
						tot_sigma = sigma_det_sq + bohr_straggling;

						if(ms_calc == 1) {
							//Gauss convolution with Pearson VII
							FWHM_G = sqrt(tot_sigma)*1000; // Gauss fwhm due to straggling and detector
							f2_ratio  = MS_FWHM_conv/FWHM_G;	// ratio of MS fwhm to Gauss fwhm
							FindFWHMandM(MS_M_conv, f2_ratio);// evaluation of convoluted fwhm and shape factor					
							tot_shape = GetConvM();		// shape factor of convoluted
							tot_fwhm = GetConvFWHM()*FWHM_G;	// fwhm of convoluted pearson vii distro
						}						
						Gauss_sum	= 0;
						Pearson_sum   = 0;

						for (int k=-GK; k<=GK; k++) {
							santykis	= k;
							Gauss_en 	= energy_left*(1-(santykis/1000));
							if(ms_calc == 0) {
								Gauss_value 	= (GenerateGaussian(Gauss_en,energy_left,tot_sigma))/GenerateGaussian(0,0,tot_sigma);
								Gauss_sum 	+= Gauss_value;
							} else if(ms_calc == 1) {
						 		Pearson_value	= CalcPearsonVII(k,tot_fwhm,tot_shape)/CalcPearsonVII(0,tot_fwhm,tot_shape);
						 		Pearson_sum 	+= Pearson_value;
							}							
						}

						for (int k=-GK; k<=GK; k++)	{
							santykis	= k;					
							Gauss_en 	= energy_left*(1-(santykis/1000));
							if(ms_calc == 0) {
								Gauss_value 	= (GenerateGaussian(Gauss_en,energy_left,tot_sigma))/GenerateGaussian(0,0,tot_sigma);
								mod_G_value 	= Gauss_value/Gauss_sum;
							} else if(ms_calc == 1) {
								Pearson_value = CalcPearsonVII(k,tot_fwhm,tot_shape)/CalcPearsonVII(0,tot_fwhm,tot_shape);
								mod_P_value 	= Pearson_value/Pearson_sum;
							}
							if (Gauss_en/MeV > ROI_region) {
								run->add_entry_reach(1);
								if(ms_calc == 1) {
									if(Gauss_en/MeV < max_permisible_energy) {
										analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	
										analysisManager->FillH1(28, Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	//layer1								
										if (i == 0)
											analysisManager->FillH1(30,Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	//el1			
										if (i == 1)
											analysisManager->FillH1(31,Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	
									}									
								} else if (ms_calc == 0) {
									analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	// total
									analysisManager->FillH1(28, Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	//layer1								
									if (i == 0)
										analysisManager->FillH1(30,Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	//el1					
									if (i == 1)
										analysisManager->FillH1(31,Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));				
								}					
							}
						}
					}
				}
    		}//end of elements
    	}	
    	run->add_total_step(tot_step);
    }// end of layer1 sensitive detector

//======================================================
// LAYER 2
//======================================================

	// layer2 sensitive detector
    if(sd2){
    	tot_step =0.;
        int n_hit_sd = sd2->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++){
            	CrystalDetectorHit* aHit = (*sd2)[i1];
            	steps = aHit->GetStep();
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();

		tot_step +=steps;

	    	G4double layer_to_surf = std::abs(-(detector->GetLength(0)/2)-layer_pos_z[1]+(detector->GetLength(2)/2));

	    	trueWorldPosition = (position.z()+detector->GetLength(0)/2)-layer_to_surf;
	    	G4double newWorldPosition = (position.z()+detector->GetLength(0)/2);
	    	
		G4ThreeVector momDir = aHit->GetWorldMomentumDirection();
		G4double rand_angle_x = atan(momDir.x()/momDir.z());	
		G4double curr_angle = rand_angle_x;
		G4double diff_angle;
		
		if(use_const_angle == 1)
		{
			diff_angle = Angle;
		}
		else
		{
			diff_angle = Angle-curr_angle;
		}


	    	for(int i=0;i<NoOfElements[2];i++)
		{
			G4double RecEn = RecoilEnergy(sample_energy,diff_angle,A1,Mnumb[2][i]);
	    		G4String el_name = sample_material[2]->GetElement(i)->GetName();
	    		
    			if(detector->GetSigmaCalc() == 1)
    			{
				xsecRTR = Get2DRTRValue(sample_energy, el_name,diff_angle);
			}
    			else 
    			{ 
    				xsecRTR = 1.;
    			}
	    		if (RecEn/MeV > 0.1)
    			{
    				RBS_yield = CalculateTotalRBSYield(sample_energy, A1, Mnumb[2][i], Z1, Znumb[2][i], diff_angle, RBS_norm_dist, solidAngle, xsecRTR, Adens[2][i],angle_of_incidence);
    				
    					G4double *l2_parameters = CalcEnergyLeft(newWorldPosition, exit_angle, RecEn,Znumb[2][i],step_number);
    					
    					energy_left 	= l2_parameters[0];
    					bohr_straggling = l2_parameters[1];
    					MS_M_conv 	= l2_parameters[2];
					MS_FWHM_conv 	= l2_parameters[3];
    				
					for(int h = 0; h<4; h++)
					{
						l2_parameters[h] = 0;
					}			

				// filling of RBS histograms
				if (energy_left > ROI_region)
				{
					// include detector dead layer energy loss
					energy_left = CalculateDeadLayerEffect(energy_left, dead_material,dead_thickness,fParticle);

					G4double dead_layer_straggling = CalculateTotalBohrStraggling(energy_left, fParticle, d_mat, dead_thickness);
					bohr_straggling += dead_layer_straggling;

					G4double rbs_depth = depth;
										
					// rbs dist
					analysisManager->FillH1(21, rbs_depth, RBS_yield); // total
					analysisManager->FillH1(35, rbs_depth, RBS_yield); // layer1
					// el1
					if (i == 0)
					{
						analysisManager->FillH1(38, rbs_depth, RBS_yield); //el1
					}
					// el2
					if (i == 1)
					{
						analysisManager->FillH1(39, rbs_depth, RBS_yield); // el2
					}

					//G4double sigma_det;

					if (Z1 == 1)
					{
						sigma_det = (det_FWHM/MeV)/2.355;
					}
					if (Z1 == 3)
					{
						if (fwhm_calc == 1)
						{
							sigma_det = (CalcDetectorFWHM(energy_left,Z1)/1000*MeV)/2.355;
						}
						else if (fwhm_calc == 0)
						{
							sigma_det = (det_FWHM/MeV)/2.355;
						}
					}
					if (Z1 == 2)
					{
						sigma_det = (det_FWHM/MeV)/2.355;
					}

					if  (sigma_det/keV < 10.*keV)
					{
						sigma_det = (10.*keV)/2.355;
					}
					G4double sigma_det_sq = std::pow(sigma_det,2);
					tot_sigma = sigma_det_sq + bohr_straggling;

					if(ms_calc == 1)
					{
						//Gauss convolution with Pearson VII
						FWHM_G = sqrt(tot_sigma)*1000; // Gauss fwhm due to straggling and detector
						f2_ratio  = MS_FWHM_conv/FWHM_G;	// ratio of MS fwhm to Gauss fwhm
						FindFWHMandM(MS_M_conv, f2_ratio);// evaluation of convoluted fwhm and shape factor					
						tot_shape = GetConvM();		// shape factor of convoluted
						tot_fwhm = GetConvFWHM()*FWHM_G;	// fwhm of convoluted pearson vii distro
						//G4cout << " shape = " << tot_shape << " fwhm = " << tot_fwhm << G4endl;
					}						
									
					Gauss_sum	= 0;
					Pearson_sum  	= 0;
					
					for (int k=-GK; k<=GK; k++)
					{
						santykis	= k;
						Gauss_en 	= energy_left*(1-(santykis/1000));
						//if (Gauss_en/MeV > ROI_region)
						//{
							if(ms_calc == 0)
							{
								Gauss_value 	= (GenerateGaussian(Gauss_en,energy_left,tot_sigma))/GenerateGaussian(0,0,tot_sigma);
								Gauss_sum 	+= Gauss_value;
							}
							else if(ms_calc == 1)
							{
								//Pearson_value 	= CalcPearsonVII_new(k, tot_fwhm,tot_shape)/CalcPearsonVII_new(0, tot_fwhm,tot_shape);
							 	Pearson_value	= CalcPearsonVII(k,tot_fwhm,tot_shape)/CalcPearsonVII(0,tot_fwhm,tot_shape);
							 	Pearson_sum 	+= Pearson_value;
							}							
						//}
					}

					for (int k=-GK; k<GK; k++)
					{
						santykis	= k;					
						Gauss_en 	= energy_left*(1-(santykis/1000));
						if(ms_calc == 0)
						{
							Gauss_value 	= (GenerateGaussian(Gauss_en,energy_left,tot_sigma))/GenerateGaussian(0,0,tot_sigma);
							mod_G_value 	= Gauss_value/Gauss_sum;
						}
						else if(ms_calc == 1)
						{
							//Pearson_value 	= CalcPearsonVII_new(k, tot_fwhm,tot_shape)/CalcPearsonVII_new(0, tot_fwhm,tot_shape);
							Pearson_value = CalcPearsonVII(k,tot_fwhm,tot_shape)/CalcPearsonVII(0,tot_fwhm,tot_shape);
							mod_P_value 	= Pearson_value/Pearson_sum;
						}
						if (Gauss_en/MeV > ROI_region)
						{
							run->add_entry_reach(1);
							if(ms_calc == 1)
							{
								if(Gauss_en/MeV < max_permisible_energy)
								{
									analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));
									analysisManager->FillH1(34, Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	//layer2								
									if (i == 0)
									{
										analysisManager->FillH1(36,Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	//el1
									}							
									if (i == 1)
									{
										analysisManager->FillH1(37,Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));
									}
								}									
							}
							else if (ms_calc == 0)
							{
								analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	// total
								analysisManager->FillH1(34, Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	//substrate								
								if (i == 0)
								{
									analysisManager->FillH1(36,Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	//el1
								}							
								if (i == 1)
								{
									analysisManager->FillH1(37,Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));
								}						
							}						
						}
					}
				}
			}//Rec en
    		} //end of elements
    	}run->add_total_step(tot_step); //all hits		
    }// end of layer2 sensitive detector

//======================================================
// LAYER 3
//======================================================

	// layer3 sensitive detector
    if(sd3){
    	tot_step =0.;
        int n_hit_sd = sd3->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++){
            	CrystalDetectorHit* aHit = (*sd3)[i1];
            	steps = aHit->GetStep();
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();

		tot_step +=steps;
		
	    	G4double layer_to_surf = std::abs(-(detector->GetLength(0)/2)-layer_pos_z[2]+(detector->GetLength(3)/2));

	    	trueWorldPosition = (position.z()+detector->GetLength(0)/2)-layer_to_surf;
	    	G4double newWorldPosition = (position.z()+detector->GetLength(0)/2);
	    	//G4cout << " true world position " << trueWorldPosition/um << G4endl;

		G4ThreeVector momDir = aHit->GetWorldMomentumDirection();
		G4double rand_angle_x = atan(momDir.x()/momDir.z());	
		G4double curr_angle = rand_angle_x;
		G4double diff_angle;
		
		if(use_const_angle == 1)
		{
			diff_angle = Angle;
		}
		else
		{
			diff_angle = Angle-curr_angle;
		}

	    	for(int i=0;i<NoOfElements[3];i++)
		{
	    		G4double RecEn = RecoilEnergy(sample_energy,diff_angle,A1,Mnumb[3][i]);
	    		//G4double RecEn = RecoilEnergy(sample_energy, Angle, A1, Mnumb[3][i]);
	    		G4String el_name = sample_material[3]->GetElement(i)->GetName();

    			if(detector->GetSigmaCalc() == 1)
    			{
				xsecRTR = Get2DRTRValue(sample_energy, el_name,diff_angle);
			}
    			else 
    			{ 
    				xsecRTR = 1.;
    			}
    			
			//G4double RBS_norm_dist = steps;
	    		if (RecEn/MeV > 0.1)
    			{	   			
    				RBS_yield = CalculateTotalRBSYield(sample_energy, A1, Mnumb[3][i], Z1, Znumb[3][i], diff_angle, RBS_norm_dist, solidAngle, xsecRTR, Adens[3][i],angle_of_incidence);    		

    				G4double *l3_parameters = CalcEnergyLeft(newWorldPosition, exit_angle, RecEn,Znumb[3][i],step_number);
    					
    				energy_left 	= l3_parameters[0];
    				bohr_straggling = l3_parameters[1];
    				MS_M_conv 	= l3_parameters[2];
				MS_FWHM_conv 	= l3_parameters[3];
    				
				for(int h = 0; h<4; h++)
				{
					l3_parameters[h] = 0;
				}
				delete[] l3_parameters;
	    			// filling of RBS histograms
				if (energy_left > ROI_region)
				{
					// include detector dead layer energy loss
					energy_left = CalculateDeadLayerEffect(energy_left, dead_material,dead_thickness,fParticle);
					G4double dead_layer_straggling = CalculateTotalBohrStraggling(energy_left, fParticle, d_mat, dead_thickness);
					bohr_straggling += dead_layer_straggling;
												
					G4double rbs_depth = depth;
					// rbs dist
					analysisManager->FillH1(21, rbs_depth, RBS_yield); // total
					analysisManager->FillH1(41, rbs_depth, RBS_yield); // layer3
					// el1
					if (i == 0)
					{
						analysisManager->FillH1(44, rbs_depth, RBS_yield); //el1
					}
					// el2
					if (i == 1)
					{
							analysisManager->FillH1(45, rbs_depth, RBS_yield); // el2
						}

					//G4double sigma_det;

					if (Z1 == 1)
					{
						sigma_det = (det_FWHM/MeV)/2.355;
					}
					if (Z1 == 3)
					{
						if (fwhm_calc == 1)
						{
							sigma_det = (CalcDetectorFWHM(energy_left,Z1)/1000*MeV)/2.355;
						}
						else if (fwhm_calc == 0)
						{
							sigma_det = (det_FWHM/MeV)/2.355;
						}
					}
					if (Z1 == 2)
					{sigma_det = (det_FWHM/MeV)/2.355;}

					if  (sigma_det/keV < 10.*keV)
					{
						sigma_det = (10.*keV)/2.355;
					}

					G4double sigma_det_sq = std::pow(sigma_det,2);
					
					tot_sigma = sigma_det_sq + bohr_straggling;

					if(ms_calc == 1)
					{
						//Gauss convolution with Pearson VII
						FWHM_G = sqrt(tot_sigma)*1000; // Gauss fwhm due to straggling and detector
						f2_ratio  = MS_FWHM_conv/FWHM_G;	// ratio of MS fwhm to Gauss fwhm
						FindFWHMandM(MS_M_conv, f2_ratio);// evaluation of convoluted fwhm and shape factor					
						tot_shape = GetConvM();		// shape factor of convoluted
						tot_fwhm = GetConvFWHM()*FWHM_G;	// fwhm of convoluted pearson vii distro
						//G4cout << " shape = " << tot_shape << " fwhm = " << tot_fwhm << G4endl;
					}								

					Gauss_sum	= 0;
					Pearson_sum   	= 0;
					
					for (int k=-GK; k<=GK; k++)
					{
						santykis	= k;
						Gauss_en 	= energy_left*(1-(santykis/1000));
						//if (Gauss_en/MeV > ROI_region)
						//{
							if(ms_calc == 0)
							{
								Gauss_value 	= (GenerateGaussian(Gauss_en,energy_left,tot_sigma))/GenerateGaussian(0,0,tot_sigma);
								Gauss_sum 	+= Gauss_value;
							}
							else if(ms_calc == 1)
							{
								//Pearson_value 	= CalcPearsonVII_new(k, tot_fwhm,tot_shape)/CalcPearsonVII_new(0, tot_fwhm,tot_shape);
							 	Pearson_value	= CalcPearsonVII(k,tot_fwhm,tot_shape)/CalcPearsonVII(0,tot_fwhm,tot_shape);
							 	Pearson_sum 	+= Pearson_value;
							}							
						//}
					}

					for (int k=-GK; k<GK; k++)
					{
						santykis	= k;					
						Gauss_en 	= energy_left*(1-(santykis/1000));
						if(ms_calc == 0)
						{
							Gauss_value 	= (GenerateGaussian(Gauss_en,energy_left,tot_sigma))/GenerateGaussian(0,0,tot_sigma);
							mod_G_value 	= Gauss_value/Gauss_sum;
						}
						else if(ms_calc == 1)
						{
							//Pearson_value 	= CalcPearsonVII_new(k, tot_fwhm,tot_shape)/CalcPearsonVII_new(0, tot_fwhm,tot_shape);
							Pearson_value = CalcPearsonVII(k,tot_fwhm,tot_shape)/CalcPearsonVII(0,tot_fwhm,tot_shape);
							mod_P_value 	= Pearson_value/Pearson_sum;
						}
						if (Gauss_en/MeV > ROI_region)
						{
							run->add_entry_reach(1);
							if(ms_calc == 1)
							{
								if(Gauss_en/MeV < max_permisible_energy)
								{
									analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	// total
									analysisManager->FillH1(40, Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	//layer3						
									if (i == 0)
									{
										analysisManager->FillH1(42,Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	//el1
									}							
									if (i == 1)
									{
										analysisManager->FillH1(43,Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	//el2
									}
								}									
							}
							else if (ms_calc == 0)
							{
								analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	// total
								analysisManager->FillH1(40, Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	//layer3								
								if (i == 0)
								{
									analysisManager->FillH1(42,Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	//el1
								}							
								if (i == 1)
								{
									analysisManager->FillH1(43,Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	//el2
								}						
							}						
						}
					}
				}
			}
	    	}// end of elements
		}
		run->add_total_step(tot_step);
    }


//======================================================
// LAYER 4
//======================================================

	// layer4 sensitive detector
    if(sd4){
    	tot_step =0.;
        int n_hit_sd = sd4->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++){
            	CrystalDetectorHit* aHit = (*sd4)[i1];
            	steps = aHit->GetStep();
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();

		tot_step +=steps;
	    	G4double layer_to_surf = std::abs(-(detector->GetLength(0)/2)-layer_pos_z[3]+(detector->GetLength(4)/2));

	    	trueWorldPosition = (position.z()+detector->GetLength(0)/2)-layer_to_surf;
	    	//G4cout << " true world position " << trueWorldPosition/um << G4endl;
	    	G4double newWorldPosition = (position.z()+detector->GetLength(0)/2);
	    	
		G4ThreeVector momDir = aHit->GetWorldMomentumDirection();
		G4double rand_angle_x = atan(momDir.x()/momDir.z());	
		G4double curr_angle = rand_angle_x;
		G4double diff_angle;
		
		if(use_const_angle == 1)
		{
			diff_angle = Angle;
		}
		else
		{
			diff_angle = Angle-curr_angle;
		}

	    	for(int i=0;i<NoOfElements[4];i++)
		{
			G4double RecEn = RecoilEnergy(sample_energy,diff_angle,A1,Mnumb[4][i]);
			

	    		//G4double RecEn = RecoilEnergy(sample_energy, Angle, A1, Mnumb[4][i]);
	    		G4String el_name = sample_material[4]->GetElement(i)->GetName();

    			if(detector->GetSigmaCalc() == 1)
    			{
				xsecRTR = Get2DRTRValue(sample_energy, el_name,diff_angle);
			}
    			else 
    			{ 
    				xsecRTR = 1.;
    			}

	    		if (RecEn/MeV > 0.1)
    			{
    				RBS_yield = CalculateTotalRBSYield(sample_energy, A1, Mnumb[4][i], Z1, Znumb[4][i], diff_angle, RBS_norm_dist, solidAngle, xsecRTR, Adens[4][i],angle_of_incidence);

    					G4double *l4_parameters = CalcEnergyLeft(newWorldPosition, exit_angle, RecEn,Znumb[4][i],step_number);
    					
    					energy_left 	= l4_parameters[0];
    					bohr_straggling = l4_parameters[1];
    					MS_M_conv 	= l4_parameters[2];
					MS_FWHM_conv 	= l4_parameters[3];
    				
					for(int h = 0; h<4; h++)
					{
						l4_parameters[h] = 0;
					}
	    				// filling of RBS histograms
				if (energy_left > ROI_region)
				{
					// include detector dead layer energy loss
					energy_left = CalculateDeadLayerEffect(energy_left, dead_material,dead_thickness,fParticle);
					G4double dead_layer_straggling = CalculateTotalBohrStraggling(energy_left, fParticle, d_mat, dead_thickness);
					bohr_straggling += dead_layer_straggling;
					G4double rbs_depth = depth;
					// rbs dist
					analysisManager->FillH1(21, rbs_depth, RBS_yield); // total
					analysisManager->FillH1(47, rbs_depth, RBS_yield); // layer4
					// el1
					if (i == 0)
					{
						analysisManager->FillH1(50, rbs_depth, RBS_yield); //el1
					}
					// el2
					if (i == 1)
					{
						analysisManager->FillH1(51, rbs_depth, RBS_yield); // el2
					}


					//G4double sigma_det;

					if (Z1 == 1)
					{
						sigma_det = (det_FWHM/MeV)/2.355;
					}
					if (Z1 == 3)
					{
						if (fwhm_calc == 1)
						{
							sigma_det = (CalcDetectorFWHM(energy_left,Z1)/1000*MeV)/2.355;
						}
						else if (fwhm_calc == 0)
						{
							sigma_det = (det_FWHM/MeV)/2.355;
						}
					}
					if (Z1 == 2)
					{sigma_det = (det_FWHM/MeV)/2.355;}

					if  (sigma_det/keV < 10.*keV)
					{
						sigma_det = (10.*keV)/2.355;
					}

					G4double sigma_det_sq = std::pow(sigma_det,2);
					tot_sigma = sigma_det_sq + bohr_straggling;

					if(ms_calc == 1)
					{
						//Gauss convolution with Pearson VII
						FWHM_G = sqrt(tot_sigma)*1000; // Gauss fwhm due to straggling and detector
						f2_ratio  = MS_FWHM_conv/FWHM_G;	// ratio of MS fwhm to Gauss fwhm
						FindFWHMandM(MS_M_conv, f2_ratio);// convoluted fwhm and shape					
						tot_shape = GetConvM();		// shape factor of convoluted
						tot_fwhm = GetConvFWHM()*FWHM_G;	// fwhm of convoluted pearson vii distro
					}								

					Gauss_sum	= 0;
					Pearson_sum   	= 0;

					for (int k=-GK; k<=GK; k++)
					{
						santykis	= k;
						Gauss_en 	= energy_left*(1-(santykis/1000));
						//if (Gauss_en/MeV > ROI_region)
						//{
							if(ms_calc == 0)
							{
								Gauss_value 	= (GenerateGaussian(Gauss_en,energy_left,tot_sigma))/GenerateGaussian(0,0,tot_sigma);
								Gauss_sum 	+= Gauss_value;
							}
							else if(ms_calc == 1)
							{
							 	Pearson_value	= CalcPearsonVII(k,tot_fwhm,tot_shape)/CalcPearsonVII(0,tot_fwhm,tot_shape);
							 	Pearson_sum 	+= Pearson_value;
							}							
						//}
					}

					for (int k=-GK; k<GK; k++)
					{
						santykis	= k;					
						Gauss_en 	= energy_left*(1-(santykis/1000));
						if(ms_calc == 0)
						{
							Gauss_value 	= (GenerateGaussian(Gauss_en,energy_left,tot_sigma))/GenerateGaussian(0,0,tot_sigma);
							mod_G_value 	= Gauss_value/Gauss_sum;
						}
						else if(ms_calc == 1)
						{
							Pearson_value = CalcPearsonVII(k,tot_fwhm,tot_shape)/CalcPearsonVII(0,tot_fwhm,tot_shape);
							mod_P_value 	= Pearson_value/Pearson_sum;
						}
						if (Gauss_en/MeV > ROI_region)
						{
						
							run->add_entry_reach(1);
							if(ms_calc == 1)
							{
								if(Gauss_en/MeV < max_permisible_energy)
								{
									analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	// total
									analysisManager->FillH1(46, Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	//layer4						
								
									if (i == 0)
									{
										analysisManager->FillH1(48,Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	//el1
									}							
									if (i == 1)
									{
										analysisManager->FillH1(49,Gauss_en, (RBS_yield*mod_P_value*(steps/RBS_norm_dist)));	//el2
									}
								}									
							}
							else if (ms_calc == 0)
							{
							
								analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	// total
								analysisManager->FillH1(46, Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	//layer4								
								if (i == 0)
								{
									analysisManager->FillH1(48,Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	//el1
								}							
								if (i == 1)
								{
									analysisManager->FillH1(49,Gauss_en, (RBS_yield*mod_G_value*(steps/RBS_norm_dist)));	//el2
								}						
							}						
						}
					}
				}
			}
			}// end of elements
		}		run->add_total_step(tot_step);
	}
	//**************************************************************************
	//**************************************************************************
	//**************************************************************************
	//**************************************************************************
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EventAction::GenerateGaussian(G4double x, G4double y, G4double sigma_sq)
{
	G4double inv_sqrt_2pi = 0.398942280401433;
	G4double a = (x-y);

	return (inv_sqrt_2pi/std::sqrt(sigma_sq)) * std::exp(-0.5*a*a/sigma_sq);

}

void EventAction::FindFWHMandM(G4double m_2, G4double f_2)
//G4double EventAction::FindFWHMandM(G4double m_2, G4double f_2)
{
	G4double FWHM = 0;
	G4double M = 0;

	std::ofstream exportname("testas.txt");
	
	MS_width_l1->Store(exportname);
	// separate cases
	// no 1
	if(f_2 == 0.01) {
		FWHM = 1;
		M = m_2;
		SetConvM(M);
		SetConvFWHM(FWHM);
		//return 1;
	} else if(f_2 == 100) {// no 2
		FWHM = 100;
		M = 100;
		SetConvM(M);
		SetConvFWHM(FWHM);
		//return 1;			
	} else {
		size_t fIndx(0);
		size_t fIndy(0);
		M = fVectorM->Value(m_2,f_2,fIndx,fIndy);
		FWHM = fVectorFWHM->Value(m_2,f_2,fIndx,fIndy);
		SetConvM(M);
		SetConvFWHM(FWHM);
	}
	//return 1;
}

	// function for Total RBS yield, combining other functions into single one
G4double EventAction::CalculateTotalRBSYield(G4double energy, G4double M1, G4double M2, G4double Z1, G4double Z2, G4double angle, G4double dist,G4double solidAngle, G4double xsecmod, G4double atomDensity,G4double inc_angle)
{
	G4double CMangleX 	= CalcAngleCMFrame(angle,M1,M2);
	G4double CMenergyX 	= CalcEnergyCMFrame(energy,M1,M2);
	G4double CMxsecX 	= CalcDiffRuthXsecCM(CMenergyX,CMangleX,Z1,Z2);
	G4double andersen_factor = CalcAndersenScreening(CMenergyX,CMangleX,Z1,Z2);
	G4double CMxsecModX	= CMxsecX*andersen_factor;
	G4double CMtoLABxsecX	= CalcDiffRuthXsecLAB(M1,M2,CMangleX,CMxsecModX);
	G4double RBSyieldX 	= CalcRBSYield((CMtoLABxsecX*xsecmod), dist, solidAngle, atomDensity, inc_angle);

	return RBSyieldX;
}

// yang+chu
G4double EventAction::CalculateTotalBohrStraggling(G4double energy, G4ParticleDefinition* particle, G4Material* mat, G4double distance)
{

	G4double straggling = 0;
	G4double nucl_strag = 0.;
	G4double elec_strag = 0.;
	
	G4int elements    	= mat->GetNumberOfElements();
    G4double *Z2 		= new G4double[elements];
    G4double *M2 		= new G4double[elements];
    G4double *aDensity 	= new G4double[elements];	

	if (energy > 0) {

	    //functions for en loss straggling evaluation
    	G4hIonEffChargeSquare* eff_charge 	= new G4hIonEffChargeSquare(""); // gets effective charge square of ion based on energy
    	G4IonYangFluctuationModel* model_yang	= new G4IonYangFluctuationModel("");
    	G4IonChuFluctuationModel* model_chu 	= new G4IonChuFluctuationModel("");

    	G4double q_squared	= 0.;
    	G4double yang_sq	= 0.;
		G4double chu_sq	= 0.;

		// Particle parameters
		q_squared 		= eff_charge->TheValue(particle,mat,energy);		// effective charge squared
		G4double Z1 		= particle->GetAtomicNumber();
		G4double M1 		= particle->GetAtomicMass();
		G4double q_sqrt 	= sqrt(q_squared);
		G4double q_sqrt_z 	= q_sqrt/Z1;
		
		yang_sq 		= model_yang->TheValue(particle,mat,energy);
		chu_sq			= model_chu->TheValue(particle,mat,energy);

		G4double strag_factor 	= q_sqrt_z*chu_sq+yang_sq;

    	const G4double *atomDensVector = mat->GetVecNbOfAtomsPerVolume();

    	for (int i=0; i<elements; i++) {
    		Z2[i] = mat->GetElement(i)->GetZ();
    		M2[i] = mat->GetElement(i)->GetA()/(g/mole);
    		aDensity[i] = atomDensVector[i]/(1/cm3);
    	}

		for (int j=0; j<elements;j++) {
			elec_strag += CalcBohrStrag(Z1,Z2[j],aDensity[j],distance)*strag_factor;	// electronic straggling straggling, in MeV2
			nucl_strag += CalcNuclEnStraggling(Z1,Z2[j],M1,M2[j],aDensity[j],distance);	// nuclear straggling, in MeV2

			if (elec_strag < 0) 
				elec_strag = 0.;
			if (nucl_strag < 0) 
				nucl_strag = 0.;
		}
	}
	straggling = elec_strag+nucl_strag;

    delete []Z2;
    delete []M2;
    delete []aDensity;

	return straggling;
}
// energy loss in detector dead layer
G4double EventAction::CalculateDeadLayerEffect(G4double energy, const G4Material* material, G4double thickness,G4ParticleDefinition* particle)
{
	G4double steps = 100.;
	G4double stp = (thickness/steps)/cm;
	G4EmCalculator emCalculator;
	
	for (int i=1; i<=steps; i++) {
		G4double stop = emCalculator.ComputeTotalDEDX(energy,particle,material)/(MeV/cm);
		energy -= stop*stp;
		if (energy <= 0.01/MeV)
			break;
	}	
	return energy;
}

// detector fwhm calculation
G4double EventAction::CalcDetectorFWHM(G4double energy, G4double Z1)
{
	// energy must be in keV
	G4double E = energy/keV;
	G4double fwhm;
	G4double C1 = 0.0999288;
	G4double C2 = 1.1871;
	G4double C3 = 1.94699;
	G4double C4 = 0.18;
	G4double C5 = 2.70004;
	G4double C6 = 9.29965;

	G4double first_term = C1*std::pow(Z1,C2);
	G4double second_term = std::pow(log(E),C3);
	G4double third_term = C4*std::pow(Z1,C5);
	G4double fourth_term = std::pow(log(E),C6);

	fwhm = first_term*second_term - (third_term/fourth_term);

	if (fwhm < 10.)
		fwhm = 10.;
	return fwhm;
}

// andersen screening factor
G4double EventAction::CalcAndersenScreening(G4double energy_cm, G4double angle_cm, G4double Z1, G4double Z2)
{
	G4double screening_factor;
	G4double V1 = (0.04873*Z1*Z2*std::sqrt(std::pow(Z1,2./3.)+std::pow(Z2,2./3.)))*keV;	// in keV
	G4double first_term = std::pow(1+(0.5*(V1/energy_cm)),2.);
	G4double second_term = V1/(2*energy_cm*sin(angle_cm/2));
	G4double second_term_sq = std::pow(second_term,2);
	G4double denominator = std::pow(1+(V1/energy_cm)+second_term_sq,2);
	screening_factor = first_term/denominator;

	return screening_factor;
}

// nuclear energy loss straggling
G4double EventAction::CalcNuclEnStraggling(G4double Z1, G4double Z2, G4double M1, G4double M2, G4double atdens,G4double distance)
{
	G4double nucl_strag = 0.;
	nucl_strag = 0.26*std::pow(Z1*Z2,2)*std::pow(M1/(M1+M2),2)*(distance/cm)*(atdens)/(1e+24); // in MeV2

	return nucl_strag;
}


G4double EventAction::CalcScreening_TF(G4double Z1, G4double Z2)
{
	G4double a_bohr = 52900; // in fm
	G4double zz1 = std::pow(Z1,2./3.);
	G4double zz2 = std::pow(Z2,2./3.);
	G4double screen_radius = a_bohr*0.8853*std::pow(zz1+zz2,-0.5);

	return screen_radius;
}
G4double EventAction::CalcScreening_ZBL(G4double Z1, G4double Z2)
{
	G4double a_bohr = 52900; // in fm
	G4double zz1 = std::pow(Z1,0.23);
	G4double zz2 = std::pow(Z2,0.23);
	G4double screen_radius = a_bohr*0.8853/(zz1+zz2);

	return screen_radius;
}

G4double EventAction::CalculateMS_spread_out(G4double E_mid, G4double E_recoil, G4ParticleDefinition* fParticle, G4Material* mat, G4double distance_out, G4double exit_angle, G4double scat_angle,G4int mat_no,G4double Mm)
{
	G4EmCalculator emCalculator;
	G4double Z1  =  fParticle->GetAtomicNumber();
	G4double A1  =  fParticle->GetAtomicMass();
	//2021-09-20 add
	//G4double primary_energy = fPrimary->GetParticleGun()->GetParticleEnergy();
	//best for alpha ~100 nm
	//best for protons ~500 nm
	G4double step_size;
	if(Z1 == 1)
		step_size = 500*nm;
	else if(Z1 == 2)
		step_size = 100*nm;
	else 
		step_size = 100*nm;	
	G4int step 	= 10;
	G4double m_tau_out = 0;
	G4double layers_part_out = distance_out/(step_size);
	// conversion to iterations
	G4int no_of_iterations_out 	= ((int)layers_part_out)+1;
	G4double distance_out_k = distance_out/(no_of_iterations_out);
	
	G4double energy_k_rec_start, energy_k_rec_finish;
	// constants
	G4double C_out;
	G4double K_out;
	//G4double K_in = 0.;
	G4double dFk, kin_fact;
	// angles
	G4double beta = (pi/2)-std::abs(exit_angle);	// exit angle

	// material parameters
	G4int elements    		= mat->GetNumberOfElements();
    const G4double *atomDensVector  = mat->GetVecNbOfAtomsPerVolume();
    G4double *Z2 			= new G4double[elements];
    G4double *M2	  		= new G4double[elements];
    G4double *aDensity 		= new G4double[elements];

	G4double miu;
    // material components for multielemental materials
    G4Material* test_mat[3];

    G4double a_bohr = 52900; // in fm
	G4double zz1 = std::pow(Z1,2./3.);	
	G4double *screen		= new G4double[elements];
	G4double *tau			= new G4double[elements];
	
    for (int i=0; i<elements; i++) {
    	Z2[i] 		= mat->GetElement(i)->GetZ();
    	M2[i] 		= mat->GetElement(i)->GetA()/(g/mole);
    	aDensity[i] 	= atomDensVector[i]/(1/cm3);
    	test_mat[i] 	= detector->GetMaterialComponents(mat_no,i);
		G4double zz2 	= std::pow(Z2[i],2./3.);
		screen[i] 	= a_bohr*0.8853*std::pow(zz1+zz2,-0.5);    
		tau[i] 	= pi * screen[i]*screen[i]*1e-26*aDensity[i]*(distance_out_k/cm);			
    }
    
    G4double miu_hat_out = 0;
	G4double tau_hat_out = 0;
	G4double tau_asterisk_out = 0.;
	// initialising primary energies
	energy_k_rec_start = E_recoil;	// recoil energy

	G4double sum_tau_out;
	G4double a_out		= 0;
	G4double b_out 	= 0;
	
	G4double sum_dphi_out;
	G4double test_miu_up;
	G4double test_miu_down;
	
	//Update parameters based on last value	
	UpdateMSParameters();	
	//Retrieve updated parameters	
	C_out		= last_C_out;
	sum_dphi_out 	= last_dphi;
	K_out 		= last_K_out;
	test_miu_up 	= last_miu_up;
	test_miu_down 	= last_miu_down;
	sum_tau_out 	= last_tau;

	for(int k = 0; k< no_of_iterations_out; k++) {
		for (int i=0;i<elements;i++) {
			energy_k_rec_finish = CalcTotEnLoss(energy_k_rec_start, distance_out_k, step, fParticle, test_mat[i]);
			G4double avg_energy_out = (energy_k_rec_start+energy_k_rec_finish)/2;
			miu 		= CalcScallingFactorMiu(avg_energy_out,Z1,Z2[i],screen[i]);		
		
			test_miu_up 	+= tau[i]/std::pow(miu,0.65);
			test_miu_down 	+= tau[i]/std::pow(miu,1.65);
	
			if (sum_dphi_out == 0)
				sum_dphi_out = tauphi->Value(tau[i])/miu;
			else {
				tau_hat_out = phitau->Value(miu*sum_dphi_out);
				sum_dphi_out = tauphi->Value(tau_hat_out+tau[i])/miu;
			}
			energy_k_rec_start	= energy_k_rec_finish;
			sum_tau_out += tau[i];
		}
	}
	miu_hat_out = test_miu_up/test_miu_down;
	
	G4int iters_for_hat = (int)(distance_out/step_size)+1;
	G4double dist_int_out = distance_out/iters_for_hat;
	G4double en_st_out	= E_recoil;

	for(int j = 0; j<iters_for_hat; j++) {
		G4double en_fin_out = CalcTotEnLoss(en_st_out,dist_int_out,10,fParticle,mat);
		G4double end_stp_pwr_out = emCalculator.ComputeTotalDEDX(en_fin_out,fParticle,mat)/(keV/nm);
		G4double start_stp_pwr_out = emCalculator.ComputeTotalDEDX(en_st_out,fParticle,mat)/(keV/nm);
		C_out *= end_stp_pwr_out/start_stp_pwr_out;	
		K_out += (dist_int_out/nm)*end_stp_pwr_out;
		en_st_out = en_fin_out;	
	}
	// Updates parameters to latest
	SetLastCout(C_out);
	SetLastKout(K_out);
	SetLastDphi(sum_dphi_out);
	SetLastTau(sum_tau_out);
	SetLastMiuUp(test_miu_up);
	SetLastMiuDown(test_miu_down);
	//
	tau_asterisk_out = phitau->Value(miu_hat_out*sum_dphi_out);
	G4double ln_aste_out = std::log(tau_asterisk_out);
	m_tau_out = 0.05754*std::pow(ln_aste_out,2.)+0.62037*ln_aste_out+2.0887;

	SetShapeFactorM2(m_tau_out);
	
	G4double energy_fwhm_out = 0.;
	
	kin_fact = KinematicFactor(scat_angle, A1, Mm);
	dFk = CalcDeltaF_K(A1,E_mid/keV,scat_angle,Mm,kin_fact);	
	a_out = -dFk*C_out;
	b_out = - (cos(beta)/sin(beta))*K_out;
	G4double nu_out = taunu->Value(sum_tau_out);
	G4double s_out = CalcScallingFactorS(a_out,b_out, nu_out);	
	
	energy_fwhm_out = std::abs((s_out*sum_dphi_out));	
	SetWidthFactorF2(energy_fwhm_out);	
	
    delete []Z2;
    delete []M2;
    delete []aDensity;	
    delete []screen;
    delete []tau;
	
	return 1;
} 

G4double EventAction::CalcDeltaF_K(G4double M1, G4double energy, G4double angle, G4double M2, G4double kinematic_factor)
{
	G4double delta_fk = 0.;
	G4double delta_fk1 = -2*M1*kinematic_factor*energy*sin(angle);
	G4double delta_fk2 = sqrt(M2*M2-(M1*M1*sin(angle)*sin(angle)));
	delta_fk = delta_fk1/delta_fk2;
	
	return delta_fk;
}
G4double EventAction::CalcScallingFactorS(G4double a, G4double b, G4double v)
{
	G4double sf		= 0;
	G4double k 		= b/a;
	G4double k_v	 	= k*(v+1);
	G4double one_plus_k 	= std::abs(1+k);
	G4double one_plus_k_sq = std::pow(one_plus_k,v+1);

	G4double plus_one 	= (one_plus_k_sq+1)/k_v;
	G4double minus_one 	= (one_plus_k_sq-1)/k_v;

	G4double v_plus 	= v+1;
	G4double inv_v		= 1/v;

	if (b==0)
		sf = a;
	else if (a==0)
		sf = b/(std::pow(v_plus,inv_v));
	else if (a*b != 0 && k < -1)
		sf = a * std::pow(std::abs(plus_one),inv_v);
	else if (a*b != 0 && k >= -1)
		sf = a * std::pow(std::abs(minus_one),inv_v);
	return sf;
}

G4double EventAction::CalcScallingFactorMiu(G4double energy, G4double Z1, G4double Z2,G4double screening_radius)
{
	G4double miu = 0.;
	G4double elmcharge_squared = 1.4399764;
	G4double miu_up = (energy/MeV)*screening_radius;
	G4double miu_down = 2*Z1*Z2*elmcharge_squared;
	miu = miu_up/miu_down;

	return miu;
}

G4double EventAction::CalcPearsonVII(G4double value,/* G4double average,*/ G4double deviation, G4double shape)
{
	// calculation of norm const
	G4double inv_sq_pi = 0.564189583547;
	G4double upper_gamma = std::exp(gamma((shape+1)/2));
	G4double lower_gamma = std::exp(gamma((shape)/2));

	G4double first_term = upper_gamma/lower_gamma*inv_sq_pi;
	G4double second_term = 2/deviation;
	G4double third_term = std::sqrt(std::pow(2,(2/(shape+1)))-1);
	G4double norm_const = first_term*second_term*third_term;
	// Pearson type VII function
	G4double first_term_p =  4*std::pow((value/*-average*/),2)*(std::pow(2,(2/(shape+1)))-1);
	G4double second_term_p = 1+(first_term_p/(deviation*deviation));
	G4double third_term_p = std::pow(second_term_p,((-(shape+1))/2));
	G4double full_term_p = norm_const*third_term_p;

	return full_term_p;
}

G4double EventAction::CalcPearsonVII_new(G4double value, G4double deviation, G4double shape)
{
	// calculation of norm const
	G4double alpha = deviation*std::sqrt(2*shape-3);
	G4double inv_sq_pi = 0.564189583547;
	G4double upper_gamma = std::exp(gamma(shape));
	G4double lower_gamma = std::exp(gamma(shape-0.5));

	G4double first_term = ((upper_gamma/lower_gamma)*inv_sq_pi)/alpha;
	// Pearson type VII function
	G4double first_term_p = 1+std::pow((value/alpha),2);
	G4double third_term_p = std::pow(first_term_p,(-shape));
	G4double full_term_p = first_term*third_term_p;

	return full_term_p;
}

void EventAction::ResetMSParameters()
{
	last_dphi = last_K_out = last_tau = last_miu_up = last_miu_down = 0.;
	last_C_out = 1.;
	SetWidthFactorF2(0.);
	SetShapeFactorM2(0.);
}

void EventAction::UpdateMSParameters()
{
	last_dphi 	= GetLastDphi();
	last_K_out 	= GetLastKout();
	last_tau 	= GetLastTau();
	last_miu_up 	= GetLastMiuUp();
	last_miu_down 	= GetLastMiuDown();
	last_C_out	= GetLastCout();
}


G4Physics2DVector* EventAction::Get2DRTRVector(G4String element, G4int A)
{
	G4Physics2DVector* vector = new G4Physics2DVector;
	G4String start_file = "RTR_values/";
	G4String end_file = "_total.txt";
	
	G4String folder_def;
	G4String part_def;
	if (A == 1) {
		folder_def = "Hydrogen/";
		part_def = "H_";
	} else if (A == 4) {
		folder_def = "Helium/";
		part_def = "He_";
	}
	G4String tot_filename = start_file+folder_def+part_def+element+end_file;
	G4cout << " Reading RTR values from file: " << tot_filename << G4endl;
	std::ifstream vFile(tot_filename);
	if(!vFile.is_open()) {
		G4cout << " No RTR file " << tot_filename << G4endl;
		exit(1);	
	}
	vector->Retrieve(vFile);

	return vector;
}

G4double EventAction::Get2DRTRValue(G4double energy, G4String elname, G4double angle)
{
	size_t fIndx(0);
	size_t fIndy(0);

	G4double value = 0.;
	if(elname == "Si")
		value = fVectorSi_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else if(elname == "O")
		value = fVectorO_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else if(elname == "Na")
		value = fVectorNa_total->Value(angle/degree,energy/keV,fIndx,fIndy);		
	else if(elname == "N")
		value = fVectorN_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else if(elname == "C")
		value = fVectorC_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else if(elname == "F")
		value = fVectorF_total->Value(angle/degree,energy/keV,fIndx,fIndy);	
	else if(elname == "B")
		value = fVectorB_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else if(elname == "Ni")
		value = fVectorNi_total->Value(angle/degree,energy/keV,fIndx,fIndy);
	else if(elname == "Cu")
		value = fVectorCu_total->Value(angle/degree,energy/keV,fIndx,fIndy);					
	else 
		value = 1.;
	return value;
}


// calculates energy loss, energy fwhm, multiple scattering
G4double* EventAction::CalcEnergyLeft(G4double depth, G4double exit_angle, G4double energy, G4double elementZ, G4double steps)
{
	G4double* final_parameters 	= new G4double[4];

	G4double ROI_region		= detector->GetRBSROImin();
	G4int MS_enabled 		= detector->GetMSCalc();
	G4double shape_factor_m		= 0.;
	G4double fwhm_factor_f		= 0.;
	G4double final_energy 		= 0.;
	G4double mid_energy 		= 0.;
	G4double final_energy_fwhm 	= 0.;
	G4double mid_energy_fwhm 	= 0.;
	G4double path_to_exit;

	// layer thicknesses
	G4double t1 = detector->GetLength(1);
	G4double t2 = detector->GetLength(2);	
	G4double t3 = detector->GetLength(3);
	G4double t4 = detector->GetLength(4);
	// Main layer shift
	G4double a = detector->GetLength(0)/2;
	// Mid positions of layers
	G4double p1 = a+detector->GetPosition(0).z();
	G4double p2 = a+detector->GetPosition(1).z();
	G4double p3 = a+detector->GetPosition(2).z();
	G4double p4 = a+detector->GetPosition(3).z();

	// positions of start and end of layers
	G4double p1s = p1-t1/2;	// pos of start
	G4double p1e = p1+t1/2;	// pos of end
	G4double p2s = p2-t2/2;
	G4double p2e = p2+t2/2;
	G4double p3s = p3-t3/2;
	G4double p3e = p3+t3/2;
	G4double p4s = p4-t4/2;
	G4double p4e = p4+t4/2;			
	
	// distance between layers
	G4double d12 = p2s-p1e;
	G4double d23 = p3s-p2e;
	G4double d34 = p4s-p3e;
	
	G4int location = 0;

	G4Material* mat;
	G4ParticleDefinition* particle = fPrimary->GetParticleGun()->GetParticleDefinition();	

	if( depth<p1s) 		 			 location = 1;	// surface
	else if (depth>p1s && depth<p1e) location = 2;  // first layer
	else if (depth>p1e && depth<p2s) location = 3;	// l1 - l2
	else if (depth>p2s && depth<p2e) location = 4;	// layer 2
	else if (depth>p2e && depth<p3s) location = 5;	// l2-l3
	else if (depth>p3s && depth<p3e) location = 6;	// layer 3
	else if (depth>p3e && depth<p4s) location = 7;	// l3-l4
	else if (depth>p4s && depth<p4e) location = 8;	// layer 4
	else if (depth>p4e) 		 	 location = 9;	// > layer4	

	G4double surf_layer = (detector->GetLength(0)/2)+(detector->GetPosition(0).z())-(detector->GetLength(1)/2);
	if (surf_layer/um < 0.0001) 
		surf_layer = 0.;	
	
	//for 2d MS vectors
	size_t fIndx(0);
	size_t fIndy(0);
	
	// calculation based on particle location
	switch(location) {
		case 1: // particle in substrate surface layer
			mat = detector->GetMaterialM(0);
			path_to_exit = depth/cos(exit_angle);
	    	final_energy 	= CalcTotEnLoss(energy, path_to_exit, steps, particle, mat);
	    	if(final_energy < ROI_region) {
	    		final_energy = 0.;
	    		final_energy_fwhm =0.;
	    		break;
	    	} else 	{
	    		final_energy_fwhm = CalculateTotalBohrStraggling(energy, particle, mat, path_to_exit);	
	    		if(MS_enabled == 1) {
	 	   			shape_factor_m = MS_shape_surf->Value(elementZ,depth/um,fIndx,fIndy);
	    			fwhm_factor_f = MS_width_surf->Value(elementZ,depth/um,fIndx,fIndy);	
	    		}
	    	}
	  		break;  			
		case 2: // particle in layer 1
			mat = detector->GetMaterialM(1);
			path_to_exit = (depth-p1s)/cos(exit_angle);		
			mid_energy 	= CalcTotEnLoss(energy,path_to_exit,steps,particle,mat);
	    	mid_energy_fwhm = CalculateTotalBohrStraggling(energy, particle, mat, path_to_exit);
	    	if(MS_enabled == 1)	{
	 	   		shape_factor_m = MS_shape_l1->Value(elementZ,depth/um,fIndx,fIndy);
	    		fwhm_factor_f = MS_width_l1->Value(elementZ,depth/um,fIndx,fIndy);		    		
	   		}		    			
			if(surf_layer != 0)	{// surface layer exists
				mat = detector->GetMaterialM(0);
				path_to_exit = surf_layer/cos(exit_angle);
	    		final_energy 	= CalcTotEnLoss(mid_energy, path_to_exit, steps, particle, mat);
	    		final_energy_fwhm = mid_energy_fwhm + CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);	
	  			break;  	    			
			} else { // no surface layer
				final_energy = mid_energy;
				final_energy_fwhm = mid_energy_fwhm;
				break;
			}
		case 3: // particle between layers 2 and 1
			// l1-l2
			mat = detector->GetMaterialM(0); 
			path_to_exit = (depth-p1e)/cos(exit_angle);
	    	if(MS_enabled == 1) {
	 	   		shape_factor_m = MS_shape_l12->Value(elementZ,depth/um,fIndx,fIndy);
	    		fwhm_factor_f = MS_width_l12->Value(elementZ,depth/um,fIndx,fIndy);		
	    	}	
			mid_energy 	= CalcTotEnLoss(energy,path_to_exit,steps,particle,mat);
	    	mid_energy_fwhm = CalculateTotalBohrStraggling(energy, particle, mat, path_to_exit);
	    	if(mid_energy > ROI_region)	{
	    			// layer 1
				mat = detector->GetMaterialM(1); 
				path_to_exit = t1/cos(exit_angle);			
	    		mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);					
				mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);				
				if(surf_layer != 0)	{// surface layer exists
					mat = detector->GetMaterialM(0);
					path_to_exit = surf_layer/cos(exit_angle);				
	    			final_energy 	= CalcTotEnLoss(mid_energy, path_to_exit, steps, particle, mat);
	    			final_energy_fwhm = mid_energy_fwhm + CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);		
	    			break;
				} else {// no surface layer
					final_energy = mid_energy;
					final_energy_fwhm = mid_energy_fwhm;
					break;	
				}				
	    	} else {
	    		final_energy = 0.;
	    		final_energy_fwhm = 0.;
	    		break;
	    	}		
		case 4: // particle in layer2
			// layer2
			mat = detector->GetMaterialM(2); 
			path_to_exit = (depth-p2s)/cos(exit_angle);
	    	if(MS_enabled == 1) {
	 	   		shape_factor_m = MS_shape_l2->Value(elementZ,depth/um,fIndx,fIndy);
	    		fwhm_factor_f = MS_width_l2->Value(elementZ,depth/um,fIndx,fIndy);		
	    	}					
			mid_energy 	= CalcTotEnLoss(energy,path_to_exit,steps,particle,mat);
	    	mid_energy_fwhm = CalculateTotalBohrStraggling(energy, particle, mat, path_to_exit);		    		
	    	if(mid_energy > ROI_region) {
	    			// l2-l1
				mat = detector->GetMaterialM(0); 
				path_to_exit = (d12)/cos(exit_angle);
				mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
				mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);
				if(mid_energy > ROI_region) {
					// layer1
					mat = detector->GetMaterialM(1);
					path_to_exit = t1/cos(exit_angle);						
					mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
					mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);					
					if(mid_energy > ROI_region) {
						if(surf_layer!= 0) {// surface layer exists
							mat = detector->GetMaterialM(0);
							path_to_exit = surf_layer/cos(exit_angle);							
	    					final_energy 	= CalcTotEnLoss(mid_energy, path_to_exit, steps, particle, mat);
	    					final_energy_fwhm = mid_energy_fwhm + CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);	
	    					break;
	    				} else {	// no surface layer, but energy > than roi
							final_energy = mid_energy;
							final_energy_fwhm = mid_energy_fwhm;
							break;						
						}    					
					} else {	// energy less than ROI
						final_energy = 0;
						final_energy_fwhm = 0;
						break;						
					}
				} else { // energy less than ROI
					final_energy = 0;
					final_energy_fwhm = 0;
					break;	
				}			
	    	} else { 	// energy less than ROI
	    		final_energy = 0.;
	    		final_energy_fwhm = 0.;
	    		break;
	    	}
		case 5: // particle between layers 3 and 2
			mat = detector->GetMaterialM(0); 
			path_to_exit = (depth-p2e)/cos(exit_angle);
	    	if(MS_enabled == 1) {
	 	   		shape_factor_m = MS_shape_l23->Value(elementZ,depth/um,fIndx,fIndy);
	    		fwhm_factor_f = MS_width_l23->Value(elementZ,depth/um,fIndx,fIndy);		
	    	}			
			mid_energy 	= CalcTotEnLoss(energy,path_to_exit,steps,particle,mat);
	    	mid_energy_fwhm = CalculateTotalBohrStraggling(energy, particle, mat, path_to_exit);	
	    	if(mid_energy > ROI_region) {
	    			// l2
				mat = detector->GetMaterialM(2); 
				path_to_exit = t2/cos(exit_angle);				
				mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
				mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);
				if(mid_energy > ROI_region) {
					// l2-l1
					mat = detector->GetMaterialM(0);
					path_to_exit = (d12)/cos(exit_angle);					
					mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
					mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);					
					if(mid_energy > ROI_region) {
						// l1
						mat = detector->GetMaterialM(1); 
						path_to_exit = t1/cos(exit_angle);						
						mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
						mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);
						if(mid_energy > ROI_region) {														
							// surface
							if(surf_layer!= 0) {// surface layer exists
								mat = detector->GetMaterialM(0);
								path_to_exit = surf_layer/cos(exit_angle);								
	    						final_energy 	= CalcTotEnLoss(mid_energy, path_to_exit, steps, particle, mat);
	    						final_energy_fwhm = mid_energy_fwhm + CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);	
	    						break;
	    					} else	{// no surface layer, but energy > than roi
								final_energy = mid_energy;
								final_energy_fwhm = mid_energy_fwhm;
								break;						
							}
						} else {	// energy less than ROI
							final_energy = 0;
							final_energy_fwhm = 0;
							break;						
						}
					} else {	// energy less than ROI
						final_energy = 0;
						final_energy_fwhm = 0;
						break;						
					}
				} else { // energy less than ROI
					final_energy = 0;
					final_energy_fwhm = 0;
					break;	
				}			
	    	} else { 	// energy less than ROI
	    		final_energy = 0.;
	    		final_energy_fwhm = 0.;
	    		break;
	    	}			    		
		case 6: // particle in layer3
			mat = detector->GetMaterialM(3); 
			path_to_exit = (depth-p3s)/cos(exit_angle);
			if(MS_enabled == 1) {
	 	   		shape_factor_m = MS_shape_l3->Value(elementZ,depth/um,fIndx,fIndy);
	    		fwhm_factor_f = MS_width_l3->Value(elementZ,depth/um,fIndx,fIndy);		
	    	}					
			mid_energy 	= CalcTotEnLoss(energy,path_to_exit,steps,particle,mat);
	    	mid_energy_fwhm = CalculateTotalBohrStraggling(energy, particle, mat, path_to_exit);	
	    	if(mid_energy > ROI_region) {
	    		// l3-l2
				mat = detector->GetMaterialM(0); 
				path_to_exit = (d23)/cos(exit_angle);				
				mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
				mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);
				if(mid_energy > ROI_region) {
					// l2
					mat = detector->GetMaterialM(2);
					path_to_exit = t2/cos(exit_angle);						
					mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
					mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);					
					if(mid_energy > ROI_region) {
						// l2-l1
						mat = detector->GetMaterialM(0); 
						path_to_exit = (d12)/cos(exit_angle);						
						mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
						mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);
						if(mid_energy > ROI_region) {							
							// layer1
							mat = detector->GetMaterialM(1);
							path_to_exit = t1/cos(exit_angle);								
							mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
							mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);							
							if(mid_energy > ROI_region) {						
								// surface
								if(surf_layer!= 0) {// surface layer exists
									mat = detector->GetMaterialM(0);
									path_to_exit = surf_layer/cos(exit_angle);
										
	    							final_energy 	= CalcTotEnLoss(mid_energy, path_to_exit, steps, particle, mat);
	    							final_energy_fwhm = mid_energy_fwhm + CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);	
	    							break;
	    						} else	{// no surface layer, but energy > than roi
									final_energy = mid_energy;
									final_energy_fwhm = mid_energy_fwhm;
									break;						
								}
							} else {	// energy less than ROI
								final_energy = 0;
								final_energy_fwhm = 0;
								break;						
							}    										
						} else {	// energy less than ROI
							final_energy = 0;
							final_energy_fwhm = 0;
							break;						
						}    					
					} else {	// energy less than ROI
						final_energy = 0;
						final_energy_fwhm = 0;
						break;						
					}
				} else { // energy less than ROI
					final_energy = 0;
					final_energy_fwhm = 0;
					break;	
				}
	    	} else { 	// energy less than ROI
    			final_energy = 0.;
    			final_energy_fwhm = 0.;
    			break;
    		}		    		
		case 7: // particle between l3-l4
			mat = detector->GetMaterialM(0); 
			path_to_exit = (depth-p3e)/cos(exit_angle);
	    	if(MS_enabled == 1) {
	 	   		shape_factor_m = MS_shape_l34->Value(elementZ,depth/um,fIndx,fIndy);
	    		fwhm_factor_f = MS_width_l34->Value(elementZ,depth/um,fIndx,fIndy);		
	    	}								
			mid_energy 	= CalcTotEnLoss(energy,path_to_exit,steps,particle,mat);
	    	mid_energy_fwhm = CalculateTotalBohrStraggling(energy, particle, mat, path_to_exit);	 		
	    	if(mid_energy > ROI_region) {
	    		// l3
				mat = detector->GetMaterialM(3); 
				path_to_exit = (t3)/cos(exit_angle);				
				mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
				mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);
				if(mid_energy > ROI_region) {
					// l3-l2
					mat = detector->GetMaterialM(0);
					path_to_exit = (d23)/cos(exit_angle);					
					mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
					mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);					
					if(mid_energy > ROI_region) {
						//l2
						mat = detector->GetMaterialM(2); 
						path_to_exit = (t2)/cos(exit_angle);							
						mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
						mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);
						if(mid_energy > ROI_region) {				
							// l2-l1
							mat = detector->GetMaterialM(0);
							path_to_exit = d12/cos(exit_angle);							
							mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
							mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);							
							if(mid_energy > ROI_region) {
								//l1
								mat = detector->GetMaterialM(1);
								path_to_exit = t1/cos(exit_angle);									
								mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
								mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);									
								if(mid_energy > ROI_region) {				
									// surface
									if(surf_layer!= 0) {// surface layer exists
										mat = detector->GetMaterialM(0);
										path_to_exit = surf_layer/cos(exit_angle);										
	    								final_energy 	= CalcTotEnLoss(mid_energy, path_to_exit, steps, particle, mat);
	    								final_energy_fwhm = mid_energy_fwhm + CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);	
	    								break;
	    							} else {	// no surface layer, but energy > than roi
										final_energy = mid_energy;
										final_energy_fwhm = mid_energy_fwhm;
										break;						
									}
								} else {	// energy less than ROI
									final_energy = 0;
									final_energy_fwhm = 0;
									break;						
								} 										
							} else {	// energy less than ROI
								final_energy = 0;
								final_energy_fwhm = 0;
								break;						
							}    									
						} else {	// energy less than ROI
							final_energy = 0;
							final_energy_fwhm = 0;
							break;						
						}    					
					} else {	// energy less than ROI
						final_energy = 0;
						final_energy_fwhm = 0;
						break;						
					}
				} else {// energy less than ROI
					final_energy = 0;
					final_energy_fwhm = 0;
					break;	
				}				
	    	} else { 	// energy less than ROI
	    		final_energy = 0.;
	    		final_energy_fwhm = 0.;
	    		break;
	    	}
		case 8: // particle in layer4
			mat = detector->GetMaterialM(4); 
			path_to_exit = (depth-p4s)/cos(exit_angle);
	    	if(MS_enabled == 1) {
	 	   		shape_factor_m = MS_shape_l4->Value(elementZ,depth/um,fIndx,fIndy);
	    		fwhm_factor_f = MS_width_l4->Value(elementZ,depth/um,fIndx,fIndy);		
	    	}			
			mid_energy 	= CalcTotEnLoss(energy,path_to_exit,steps,particle,mat);
	    	mid_energy_fwhm = CalculateTotalBohrStraggling(energy, particle, mat, path_to_exit);	    		
	    	if(mid_energy > ROI_region) {
	    		// l4-l3
				mat = detector->GetMaterialM(0); 
				path_to_exit = (d34)/cos(exit_angle);				
				mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
				mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);
				if(mid_energy > ROI_region) {
					// l3
					mat = detector->GetMaterialM(3);
					path_to_exit = (t3)/cos(exit_angle);					
					mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
					mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);					
					if(mid_energy > ROI_region) {
						//l3-l2
						mat = detector->GetMaterialM(0); 
						path_to_exit = (d23)/cos(exit_angle);						
						mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
						mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);
						if(mid_energy > ROI_region) {				
							// l2
							mat = detector->GetMaterialM(2);
							path_to_exit = t2/cos(exit_angle);								
							mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
							mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);							
							if(mid_energy > ROI_region) {
								// l2-l1
								mat = detector->GetMaterialM(0);
								path_to_exit = d12/cos(exit_angle);								
								mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
								mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);									
								if(mid_energy > ROI_region) {
									// l1
									mat = detector->GetMaterialM(1);
									path_to_exit = t1/cos(exit_angle);									
									mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
									mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);									
									if(mid_energy > ROI_region) {															
										// surface
										if(surf_layer!= 0) {// surface layer exists
											mat = detector->GetMaterialM(0);
											path_to_exit = surf_layer/cos(exit_angle);												
	    									final_energy 	= CalcTotEnLoss(mid_energy, path_to_exit, steps, particle, mat);
	    									final_energy_fwhm = mid_energy_fwhm + CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);	
	    									break;
	    								} else {	// no surface layer, but energy > than roi
											final_energy = mid_energy;
											final_energy_fwhm = mid_energy_fwhm;
											break;						
										}
									} else {	// energy less than ROI
										final_energy = 0;
										final_energy_fwhm = 0;
										break;						
									} 													
								} else {	// energy less than ROI
									final_energy = 0;
									final_energy_fwhm = 0;
									break;						
								} 										
							} else {	// energy less than ROI
								final_energy = 0;
								final_energy_fwhm = 0;
								break;						
							}    									
						} else {	// energy less than ROI
							final_energy = 0;
							final_energy_fwhm = 0;
							break;						
						}    					
					} else {	// energy less than ROI
						final_energy = 0;
						final_energy_fwhm = 0;
						break;						
					}
				} else { // energy less than ROI
					final_energy = 0;
					final_energy_fwhm = 0;
					break;	
				}				
	    	} else { 	// energy less than ROI
    			final_energy = 0.;
    			final_energy_fwhm = 0.;
    			break;
    		}	    		
		case 9: // particle in outside layer 4
			mat = detector->GetMaterialM(0); 
			path_to_exit = (depth-p4e)/cos(exit_angle);
	    	if(MS_enabled == 1) {
	 	   		shape_factor_m = MS_shape_l4e->Value(elementZ,depth/um,fIndx,fIndy);
	    		fwhm_factor_f = MS_width_l4e->Value(elementZ,depth/um,fIndx,fIndy);		
	    	}		
			mid_energy 	= CalcTotEnLoss(energy,path_to_exit,steps,particle,mat);
	    	mid_energy_fwhm = CalculateTotalBohrStraggling(energy, particle, mat, path_to_exit);		
	    	if(mid_energy > ROI_region) {
	    		// l4
				mat = detector->GetMaterialM(4); 
				path_to_exit = (t4)/cos(exit_angle);				
				mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
				mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);
				if(mid_energy > ROI_region) {
					// l4-l3
					mat = detector->GetMaterialM(0);
					path_to_exit = (d34)/cos(exit_angle);					
					mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
					mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);					
					if(mid_energy > ROI_region) {
						//l3
						mat = detector->GetMaterialM(3); 
						path_to_exit = (t3)/cos(exit_angle);						
						mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
						mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);
						if(mid_energy > ROI_region) {			
							// l3-l2
							mat = detector->GetMaterialM(0);
							path_to_exit = d23/cos(exit_angle);							
							mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
							mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);							
							if(mid_energy > ROI_region) {		
								// l2
								mat = detector->GetMaterialM(2);
								path_to_exit = t2/cos(exit_angle);								
								mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
								mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);									
								if(mid_energy > ROI_region) {
									// l2-l1
									mat = detector->GetMaterialM(0);
									path_to_exit = d12/cos(exit_angle);										
									mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
									mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);									
									if(mid_energy > ROI_region) {
										//l1
										mat = detector->GetMaterialM(1);
										path_to_exit = t1/cos(exit_angle);										
										mid_energy_fwhm += CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);
										mid_energy 	= CalcTotEnLoss(mid_energy,path_to_exit,steps,particle,mat);												
										if(mid_energy > ROI_region) {														
											// surface
											if(surf_layer!= 0) {// surface layer exists
												mat = detector->GetMaterialM(0);
												path_to_exit = surf_layer/cos(exit_angle);													
	    										final_energy 	= CalcTotEnLoss(mid_energy, path_to_exit, steps, particle, mat);
	    										final_energy_fwhm = mid_energy_fwhm + CalculateTotalBohrStraggling(mid_energy, particle, mat, path_to_exit);	
	    										break;
	    									} else {	// no surface layer, but energy > than roi
												final_energy = mid_energy;
												final_energy_fwhm = mid_energy_fwhm;
												break;						
											}
										} else {	// energy less than ROI
											final_energy = 0;
											final_energy_fwhm = 0;
											break;						
										} 												
									} else {	// energy less than ROI
										final_energy = 0;
										final_energy_fwhm = 0;
										break;						
									} 													
								} else {	// energy less than ROI
									final_energy = 0;
									final_energy_fwhm = 0;
									break;						
								} 										
							} else {	// energy less than ROI
								final_energy = 0;
								final_energy_fwhm = 0;
								break;						
							}    										
						} else {	// energy less than ROI
							final_energy = 0;
							final_energy_fwhm = 0;
							break;						
						}    					
					} else {	// energy less than ROI
						final_energy = 0;
						final_energy_fwhm = 0;
						break;						
					}
				} else { // energy less than ROI
					final_energy = 0;
					final_energy_fwhm = 0;
					break;	
				}				
	    	} else { 	// energy less than ROI
	    		final_energy = 0.;
	    		final_energy_fwhm = 0.;
	    		break;
	    	}
		default:
	    	break;
	}


	if(final_energy < 0) 
		final_energy = 0;
	if(final_energy_fwhm < 0) 
		final_energy_fwhm = 0.;

	if(final_energy == 0) {
		shape_factor_m = 0.;
		fwhm_factor_f = 0.;
	}
	final_parameters[0] = final_energy;
	final_parameters[1] = final_energy_fwhm;
	final_parameters[2] = shape_factor_m;
	final_parameters[3] = fwhm_factor_f;

	// returns final energy and fwhm as array
	// if ms calc is enabled, shape factor m and fwhm factor f is also returned
	// if ms disabled, m and f returns 0
	return final_parameters;
}

G4double EventAction::FillTauPhiVector()
{
	std::ifstream vFileIn;
	G4String filename = "MS/TAU_PHI_NU.DAT";
	G4cout << " Reading TAU_PHI values from " << filename << G4endl;
    vFileIn.open(filename);
    G4String line;
    if(!vFileIn) {
    	G4cout <<" No TAU-PHI file " << G4endl;
        exit(1); 
    }
    G4int linenum = 0;

	vFileIn >> linenum;
	
	G4double val1, val2, val3;
	
	tauphi = new G4PhysicsFreeVector(linenum);	
	phitau = new G4PhysicsFreeVector(linenum);
	taunu = new G4PhysicsFreeVector(linenum);
	
	// filling of vectors
	for(int i = 0; i< linenum; i++) {
		vFileIn >> val1 >> val2 >> val3;
		tauphi->PutValue(i,val1,val2);
		phitau->PutValue(i,val2,val1);
		taunu->PutValue(i,val1,val3);
	}
	vFileIn.close();
	return 1;
}

G4double EventAction::Fill_MS_Vector(G4double primary_energy)
{
	G4double incidence_angle = atan(fPrimary->GetParticleGun()->GetParticleMomentumDirection().x()/fPrimary->GetParticleGun()->GetParticleMomentumDirection().z());
    G4double incidence_modifier = cos(incidence_angle);
    	
	G4ParticleDefinition* particle = fPrimary->GetParticleGun()->GetParticleDefinition();	
    G4double A1 = particle->GetAtomicMass();
    	//G4double Z1 = particle->GetAtomicNumber();	
	G4double scat_angle = detector->GetRBSAngle();

	G4int use_ms_corrections	= detector->GetUseMSCorrections();

	G4double path_to_exit;
	G4double step_interval = 2.*nm;
	G4double standard_size = detector->GetLength(0);
	G4int number_of_intervals = round(standard_size/step_interval);
	
	G4double end_en = 0;
	G4double start_en = 0;
	G4double recoil_st, recoil_end;
	G4double exit_angle = pi-detector->GetRBSAngle()-incidence_angle;

	G4double pseudo_depth;

	// layer thicknesses
	G4double t1 = detector->GetLength(1);
	G4double t2 = detector->GetLength(2);	
	G4double t3 = detector->GetLength(3);
	G4double t4 = detector->GetLength(4);

	// Main layer shift
	G4double a = detector->GetLength(0)/2;
	// Mid positions of layers
	G4double p1 = a+detector->GetPosition(0).z();
	G4double p2 = a+detector->GetPosition(1).z();
	G4double p3 = a+detector->GetPosition(2).z();
	G4double p4 = a+detector->GetPosition(3).z();

	// positions of start and end of layers
	G4double p1s = p1-t1/2;	// pos of start
	G4double p1e = p1+t1/2;	// pos of end
	G4double p2s = p2-t2/2;
	G4double p2e = p2+t2/2;
	G4double p3s = p3-t3/2;
	G4double p3e = p3+t3/2;
	G4double p4s = p4-t4/2;
	G4double p4e = p4+t4/2;			
	
	// distance between layers
	G4double d12 = p2s-p1e;
	G4double d23 = p3s-p2e;
	G4double d34 = p4s-p3e;
	
	G4int location = 0;
	//G4double Z2;
	G4double M2;

	G4Material* mat;
	G4int no_of_el;
	ResetMSParameters();

  	G4double* shape_m 		= new G4double [5];
	G4double* width_f 		= new G4double [5];
  	
  	std::ofstream myfile_shape_surf;
	std::ofstream myfile_width_surf;
  	std::ofstream myfile_shape_l1;
	std::ofstream myfile_width_l1;	
  	std::ofstream myfile_shape_l1_l2;
	std::ofstream myfile_width_l1_l2;
  	std::ofstream myfile_shape_l2;
	std::ofstream myfile_width_l2;	
  	std::ofstream myfile_shape_l2_l3;
	std::ofstream myfile_width_l2_l3;
  	std::ofstream myfile_shape_l3;
	std::ofstream myfile_width_l3;	
  	std::ofstream myfile_shape_l3_l4;
	std::ofstream myfile_width_l3_l4;	
  	std::ofstream myfile_shape_l4;
	std::ofstream myfile_width_l4;	
  	std::ofstream myfile_shape_l4_end;
	std::ofstream myfile_width_l4_end;				
	
  	myfile_shape_surf.open ("TEMP/surf_temp_shape.txt");	
  	myfile_width_surf.open ("TEMP/surf_temp_width.txt");   	
  	myfile_shape_l1.open ("TEMP/l1_temp_shape.txt");
  	myfile_width_l1.open ("TEMP/l1_temp_width.txt");    
  	myfile_shape_l1_l2.open ("TEMP/l1-l2_temp_shape.txt");
  	myfile_width_l1_l2.open ("TEMP/l1-l2_temp_width.txt");      	
  	myfile_shape_l2.open ("TEMP/l2_temp_shape.txt");
  	myfile_width_l2.open ("TEMP/l2_temp_width.txt");      	  
  	myfile_shape_l2_l3.open ("TEMP/l2-l3_temp_shape.txt");
  	myfile_width_l2_l3.open ("TEMP/l2-l3_temp_width.txt");  
  	myfile_shape_l3.open ("TEMP/l3_temp_shape.txt");
  	myfile_width_l3.open ("TEMP/l3_temp_width.txt");      	  	
  	myfile_shape_l3_l4.open ("TEMP/l3-l4_temp_shape.txt");
  	myfile_width_l3_l4.open ("TEMP/l3-l4_temp_width.txt");        
  	myfile_shape_l4.open ("TEMP/l4_temp_shape.txt");
  	myfile_width_l4.open ("TEMP/l4_temp_width.txt");      	  	
  	myfile_shape_l4_end.open ("TEMP/l4-end_temp_shape.txt");
  	myfile_width_l4_end.open ("TEMP/l4-end_temp_width.txt");        		      				
		
	start_en = primary_energy;
	

	G4double surf_layer = (detector->GetLength(0)/2)+(detector->GetPosition(0).z())-(detector->GetLength(1)/2);
	if (surf_layer/um < 0.0001) { surf_layer = 0.;}	
	
	// calculate number of lines
	G4int no_of_lines_surf 	= round(surf_layer/step_interval);
	G4int no_of_lines_l1   	= round(t1/step_interval);
	G4int no_of_lines_l1_l2   	= round(d12/step_interval);
	G4int no_of_lines_l2   	= round(t2/step_interval);
	G4int no_of_lines_l2_l3   	= round(d23/step_interval);
	G4int no_of_lines_l3   	= round(t3/step_interval);
	G4int no_of_lines_l3_l4   	= round(d34/step_interval);
	G4int no_of_lines_l4   	= round(t4/step_interval);
	G4int no_of_lines_l4_end   	= round((2*a-p4e)/step_interval);	
	// print no of lines into file
	// number of elements 
	G4int number_el = detector->GetMaterialM(0)->GetNumberOfElements();
	G4double temp;
	G4double temp2;
	G4double Z;
	// **************************************************
	// setting up of files for 2d phys vector read
	// surface layer
	// needs to sort elements Z by increasing order, because otherwise the physics vectors doesn't work as intended
	G4int number_el_l1 = detector->GetMaterialM(1)->GetNumberOfElements();
	G4int number_el_l2 = detector->GetMaterialM(2)->GetNumberOfElements();
	G4int number_el_l3 = detector->GetMaterialM(3)->GetNumberOfElements();
	G4int number_el_l4 = detector->GetMaterialM(4)->GetNumberOfElements();			
	//************************
	G4double* mat0_Z = new G4double[number_el];
	G4double* mat1_Z = new G4double[number_el_l1];
	G4double* mat2_Z = new G4double[number_el_l2];
	G4double* mat3_Z = new G4double[number_el_l3];
	G4double* mat4_Z = new G4double[number_el_l4];
	
	G4double* mat0_M = new G4double[number_el];
	G4double* mat1_M = new G4double[number_el_l1];
	G4double* mat2_M = new G4double[number_el_l2];
	G4double* mat3_M = new G4double[number_el_l3];
	G4double* mat4_M = new G4double[number_el_l4];	
	// substrate sorting
	for(int i = 0; i<number_el; i++) {
		mat0_Z[i] = detector->GetMaterialM(0)->GetElement(i)->GetZ();
		mat0_M[i] = detector->GetMaterialM(0)->GetElement(i)->GetA()/(g/mole);
	}
	for(int i = 0; i<number_el; i++) {
		for(int j=i+1; j<number_el;j++) {	
			if(mat0_Z[i]>mat0_Z[j]) {
				temp = mat0_Z[i];
				mat0_Z[i] = mat0_Z[j];
				mat0_Z[j] = temp;
			}	
			if(mat0_M[i]>mat0_M[j]) {
				temp2 = mat0_M[i];
				mat0_M[i] = mat0_M[j];
				mat0_M[j] = temp2;
			}			
		}
	}
	// layer 1 sorting
	for(int i = 0; i<number_el_l1; i++) {
		mat1_Z[i] = detector->GetMaterialM(1)->GetElement(i)->GetZ();
		mat1_M[i] = detector->GetMaterialM(1)->GetElement(i)->GetA()/(g/mole);
	}
	for(int i = 0; i<number_el_l1; i++) {
		for(int j=i+1; j<number_el_l1;j++) {	
			if(mat1_Z[i]>mat1_Z[j]) {
				temp = mat1_Z[i];
				mat1_Z[i] = mat1_Z[j];
				mat1_Z[j] = temp;
			}	
			if(mat1_M[i]>mat1_M[j])	{
				temp2 = mat1_M[i];
				mat1_M[i] = mat1_M[j];
				mat1_M[j] = temp2;
			}			
		}
	}
	// layer 2 sorting
	for(int i = 0; i<number_el_l2; i++) {
		mat2_Z[i] = detector->GetMaterialM(2)->GetElement(i)->GetZ();
		mat2_M[i] = detector->GetMaterialM(2)->GetElement(i)->GetA()/(g/mole);
	}
	for(int i = 0; i<number_el_l2; i++)	{
		for(int j=i+1; j<number_el_l2;j++)	{	
			if(mat2_Z[i]>mat2_Z[j])	{
				temp = mat2_Z[i];
				mat2_Z[i] = mat2_Z[j];
				mat2_Z[j] = temp;
			}	
			if(mat2_M[i]>mat2_M[j])	{
				temp2 = mat2_M[i];
				mat2_M[i] = mat2_M[j];
				mat2_M[j] = temp2;
			}			
		}
	}
	// layer 3 sorting
	for(int i = 0; i<number_el_l3; i++)	{
		mat3_Z[i] = detector->GetMaterialM(3)->GetElement(i)->GetZ();
		mat3_M[i] = detector->GetMaterialM(3)->GetElement(i)->GetA()/(g/mole);
	}
	for(int i = 0; i<number_el_l3; i++)	{
		for(int j=i+1; j<number_el_l3;j++) {	
			if(mat3_Z[i]>mat3_Z[j])	{
				temp = mat3_Z[i];
				mat3_Z[i] = mat3_Z[j];
				mat3_Z[j] = temp;
			}	
			if(mat3_M[i]>mat3_M[j])	{
				temp2 = mat3_M[i];
				mat3_M[i] = mat3_M[j];
				mat3_M[j] = temp2;
			}			
		}
	}
	// layer 4 sorting
	for(int i = 0; i<number_el_l4; i++)	{
		mat4_Z[i] = detector->GetMaterialM(4)->GetElement(i)->GetZ();
		mat4_M[i] = detector->GetMaterialM(4)->GetElement(i)->GetA()/(g/mole);
	}
	for(int i = 0; i<number_el_l4; i++)	{
		for(int j=i+1; j<number_el_l4;j++) {	
			if(mat4_Z[i]>mat4_Z[j])	{
				temp = mat4_Z[i];
				mat4_Z[i] = mat4_Z[j];
				mat4_Z[j] = temp;
			}	
			if(mat4_M[i]>mat4_M[j])	{
				temp2 = mat4_M[i];
				mat4_M[i] = mat4_M[j];
				mat4_M[j] = temp2;
			}			
		}
	}
	// end of sorting
	//******************
	// fill of output files
	myfile_shape_surf << 4 << " " << number_el << " " << no_of_lines_surf << "\n";
	myfile_width_surf << 4 << " " << number_el << " " << no_of_lines_surf << "\n";

	for(int i = 0; i<number_el;i++)	{
		Z = mat0_Z[i];
		myfile_shape_surf << Z << " ";
		myfile_width_surf << Z << " ";		
	}
	myfile_shape_surf << "\n";
	myfile_width_surf << "\n";	

	for(int i = 0; i<no_of_lines_surf; i++) {
		myfile_shape_surf << (step_interval*(i+1))/um << " ";
		myfile_width_surf << (step_interval*(i+1))/um << " ";		
	}	
	myfile_shape_surf << "\n";
	myfile_width_surf << "\n";	
	// layer 1-layer2
	myfile_shape_l1_l2 << 4 << " " << number_el << " " << no_of_lines_l1_l2 << "\n";
	myfile_width_l1_l2 << 4 << " " << number_el << " " << no_of_lines_l1_l2 << "\n";
	for(int i = 0; i<number_el;i++)	{
		Z = mat0_Z[i];
		myfile_shape_l1_l2 << Z << " ";
		myfile_width_l1_l2 << Z << " ";		
	}
	myfile_shape_l1_l2 << "\n";
	myfile_width_l1_l2 << "\n";	

	for(int i = 0; i<no_of_lines_l1_l2; i++) {
		myfile_shape_l1_l2 << (p1e/um)+(step_interval*(i+1))/um << " ";
		myfile_width_l1_l2 << (p1e/um)+(step_interval*(i+1))/um << " ";		
	}	
	myfile_shape_l1_l2 << "\n";
	myfile_width_l1_l2 << "\n";	
	// layer2-layer3
	myfile_shape_l2_l3 << 4 << " " << number_el << " " << no_of_lines_l2_l3 << "\n";
	myfile_width_l2_l3 << 4 << " " << number_el << " " << no_of_lines_l2_l3 << "\n";
	for(int i = 0; i<number_el;i++)	{
		Z = mat0_Z[i];
		myfile_shape_l2_l3 << Z << " ";
		myfile_width_l2_l3 << Z << " ";		
	}
	myfile_shape_l2_l3 << "\n";
	myfile_width_l2_l3 << "\n";	

	for(int i = 0; i<no_of_lines_l2_l3; i++) {
		myfile_shape_l2_l3 << (p2e/um)+(step_interval*(i+1))/um << " ";
		myfile_width_l2_l3 << (p2e/um)+(step_interval*(i+1))/um << " ";		
	}	
	myfile_shape_l2_l3 << "\n";
	myfile_width_l2_l3 << "\n";	
	// layer3-layer4
	myfile_shape_l3_l4 << 4 << " " << number_el << " " << no_of_lines_l3_l4 << "\n";
	myfile_width_l3_l4 << 4 << " " << number_el << " " << no_of_lines_l3_l4 << "\n";
	for(int i = 0; i<number_el;i++)	{
		Z = mat0_Z[i];
		myfile_shape_l3_l4 << Z << " ";
		myfile_width_l3_l4 << Z << " ";		
	}
	myfile_shape_l3_l4 << "\n";
	myfile_width_l3_l4 << "\n";	

	for(int i = 0; i<no_of_lines_l3_l4; i++) {
		myfile_shape_l3_l4 << (p3e/um)+(step_interval*(i+1))/um << " ";
		myfile_width_l3_l4 << (p3e/um)+(step_interval*(i+1))/um << " ";		
	}	
	myfile_shape_l3_l4 << "\n";
	myfile_width_l3_l4 << "\n";		
	// layer4-end
	myfile_shape_l4_end << 4 << " " << number_el << " " << no_of_lines_l4_end << "\n";
	myfile_width_l4_end << 4 << " " << number_el << " " << no_of_lines_l4_end << "\n";
	for(int i = 0; i<number_el;i++)	{
		Z = mat0_Z[i];
		myfile_shape_l4_end << Z << " ";
		myfile_width_l4_end << Z << " ";	
	}
	myfile_shape_l4_end << "\n";
	myfile_width_l4_end << "\n";	

	for(int i = 0; i<no_of_lines_l4_end; i++) {
		myfile_shape_l4_end << (p4e/um)+(step_interval*(i+1))/um << " ";
		myfile_width_l4_end << (p4e/um)+(step_interval*(i+1))/um << " ";		
	}	
	myfile_shape_l4_end << "\n";
	myfile_width_l4_end << "\n";		
	//layer1
	number_el = detector->GetMaterialM(1)->GetNumberOfElements();
	myfile_shape_l1 << 4 << " " << number_el << " " << no_of_lines_l1 << "\n";
	myfile_width_l1 << 4 << " " << number_el << " " << no_of_lines_l1 << "\n";
	for(int i = 0; i<number_el;i++)	{
		Z = mat1_Z[i];
		myfile_shape_l1 << Z << " ";
		myfile_width_l1 << Z << " ";		
	}
	myfile_shape_l1 << "\n";
	myfile_width_l1 << "\n";	

	for(int i = 0; i<no_of_lines_l1; i++) {
		myfile_shape_l1 << (p1s/um)+(step_interval*(i+1))/um << " ";
		myfile_width_l1 << (p1s/um)+(step_interval*(i+1))/um << " ";		
	}	
	myfile_shape_l1 << "\n";
	myfile_width_l1 << "\n";
	//layer2
	number_el = detector->GetMaterialM(2)->GetNumberOfElements();
	myfile_shape_l2 << 4 << " " << number_el << " " << no_of_lines_l2 << "\n";
	myfile_width_l2 << 4 << " " << number_el << " " << no_of_lines_l2 << "\n";
	
	for(int i = 0; i<number_el;i++)	{
		Z = mat2_Z[i];
		myfile_shape_l2 << Z << " ";
		myfile_width_l2 << Z << " ";			
	}
	myfile_shape_l2 << "\n";
	myfile_width_l2 << "\n";	

	for(int i = 0; i<no_of_lines_l2; i++) {
		myfile_shape_l2 << (p1s/um)+(step_interval*(i+1))/um << " ";
		myfile_width_l2 << (p1s/um)+(step_interval*(i+1))/um << " ";		
	}	
	myfile_shape_l2 << "\n";
	myfile_width_l2 << "\n";
	//layer3
	number_el = detector->GetMaterialM(3)->GetNumberOfElements();
	myfile_shape_l3 << 4 << " " << number_el << " " << no_of_lines_l3 << "\n";
	myfile_width_l3 << 4 << " " << number_el << " " << no_of_lines_l3 << "\n";
	for(int i = 0; i<number_el;i++)	{
		Z = mat3_Z[i];
		myfile_shape_l3 << Z << " ";
		myfile_width_l3 << Z << " ";	
	}
	myfile_shape_l3 << "\n";
	myfile_width_l3 << "\n";	

	for(int i = 0; i<no_of_lines_l3; i++) {
		myfile_shape_l3 << (p1s/um)+(step_interval*(i+1))/um << " ";
		myfile_width_l3 << (p1s/um)+(step_interval*(i+1))/um << " ";		
	}	
	myfile_shape_l3 << "\n";
	myfile_width_l3 << "\n";	
	//layer4
	number_el = detector->GetMaterialM(4)->GetNumberOfElements();
	myfile_shape_l4 << 4 << " " << number_el << " " << no_of_lines_l4 << "\n";
	myfile_width_l4 << 4 << " " << number_el << " " << no_of_lines_l4 << "\n";
	for(int i = 0; i<number_el;i++)	{
		Z = mat4_Z[i];
		myfile_shape_l4 << Z << " ";
		myfile_width_l4 << Z << " ";			
	}
	myfile_shape_l4 << "\n";
	myfile_width_l4 << "\n";	

	for(int i = 0; i<no_of_lines_l4; i++) {
		myfile_shape_l4 << (p1s/um)+(step_interval*(i+1))/um << " ";
		myfile_width_l4 << (p1s/um)+(step_interval*(i+1))/um << " ";		
	}	
	myfile_shape_l4 << "\n";
	myfile_width_l4 << "\n";				
	// **************************************************	

 	//G4cout << " mod " << incidence_modifier << G4endl;
	
	G4double calc_percent = 100./number_of_intervals;
	G4cout << " Calculating Multiple Scattering Energy spread and shape " << G4endl;	
	//G4cout << " incidence " << incidence_angle/deg << " exit " << exit_angle/deg << " scat " << scat_angle/deg << G4endl;
	//G4cout << " Progress: " << G4endl;
	for(int i = 1; i<=number_of_intervals; i++) {
		// Progress
		if(i % 5000 == 0)
			G4cout << " " << calc_percent*i << " % " << G4endl;
		
		pseudo_depth = step_interval*i;

		if( pseudo_depth<p1s) 		 		location = 1;	// surface
		else if (pseudo_depth>p1s && pseudo_depth<p1e)	location = 2;  // first layer
		else if (pseudo_depth>p1e && pseudo_depth<p2s)	location = 3;	// l1 - l2
		else if (pseudo_depth>p2s && pseudo_depth<p2e)	location = 4;	// layer 2
		else if (pseudo_depth>p2e && pseudo_depth<p3s)	location = 5;	// l2-l3
		else if (pseudo_depth>p3s && pseudo_depth<p3e)	location = 6;	// layer 3
		else if (pseudo_depth>p3e && pseudo_depth<p4s)	location = 7;	// l3-l4
		else if (pseudo_depth>p4s && pseudo_depth<p4e)	location = 8;	// layer 4
		else if (pseudo_depth>p4e) 				location = 9;	// > layer4	
		
		// calculation based on particle location
		switch(location) {
			case 1: // particle in substrate surface layer
				mat = detector->GetMaterialM(0);
				no_of_el = mat->GetNumberOfElements();
				path_to_exit = pseudo_depth/cos(exit_angle);
				end_en = CalcTotEnLoss(start_en, step_interval/incidence_modifier, 10, particle, mat);				
				for(int j = 0; j<no_of_el;j++) {
					ResetMSParameters();
	    			M2 = mat0_M[j];
	    			recoil_st = RecoilEnergy(end_en,scat_angle,A1,M2);
	    			recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
	    			if(recoil_end > 100.*keV) {
						CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
	    				shape_m[j]		= GetShapeFactorM2();
						
						if(use_ms_corrections == 1)
							width_f[j] 	= GetWidthFactorF2()*last_C_out;
						else
							width_f[j]	= GetWidthFactorF2();	
						/*
						if(j == 0 && (int)(pseudo_depth/nm)%50 == 0)  {
							G4cout << " depth = " << pseudo_depth/nm << " deg = " << GetLastDphi()*57.2957795 << " deg, energy =  " <<  width_f[j] << G4endl;

						}*/
	    			} else { 
	    				shape_m[j]		= 0;
						width_f[j]		= 0;		    				
	    			}
	    			myfile_shape_surf << shape_m[j] << "\t";
					myfile_width_surf << width_f[j] << "\t";	
		
	    		}
	    		myfile_shape_surf << "\n";
				myfile_width_surf << "\n";		
	  			break;  		
	  			// end of case 1	
	  			
			case 2: // particle in layer 1
				// layer1
				mat = detector->GetMaterialM(1);
				no_of_el = mat->GetNumberOfElements();	
				end_en = CalcTotEnLoss(start_en, step_interval/incidence_modifier, 10, particle, mat);	
				for(int j = 0; j<no_of_el;j++) {
					mat = detector->GetMaterialM(1);
					path_to_exit = (pseudo_depth-p1s)/cos(exit_angle);					
					ResetMSParameters();
	    			M2 = mat1_M[j];
	    			recoil_st = RecoilEnergy(end_en,scat_angle,A1,M2);
	    			recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
	    			
	    			if(recoil_end > 100.*keV) {
	   					CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,1,M2);	
						shape_m[j]		= GetShapeFactorM2();
						width_f[j]		= GetWidthFactorF2();
						// surface layer substrate, if it exists	
						if(surf_layer != 0) {
							mat = detector->GetMaterialM(0);
							path_to_exit = surf_layer/cos(exit_angle);	
							recoil_st = recoil_end;
							recoil_end = CalcTotEnLoss(recoil_end, path_to_exit,100, particle, mat);

							if(recoil_end > 100*keV) {																
		    					CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
		    					if(use_ms_corrections == 1)
									width_f[j] 	= GetWidthFactorF2()*last_C_out;
								else
									width_f[j]	= GetWidthFactorF2();												
							} else { 
								shape_m[j]		= 0;
								width_f[j]		= 0;	
							}
						}
	    			} else {
	    				shape_m[j]		= 0;
						width_f[j]		= 0;		    				
	    				break;
	    			}
	    			myfile_shape_l1 << shape_m[j] << "\t";
					myfile_width_l1 << width_f[j] << "\t";
	    			}
				myfile_shape_l1 << "\n";
				myfile_width_l1 << "\n";
	  			break;  
	  			// end of case 2
	  			
			case 3: // particle between layers l2-l1
				// substrate
				mat = detector->GetMaterialM(0);
				no_of_el = mat->GetNumberOfElements();	
				end_en = CalcTotEnLoss(start_en, step_interval/incidence_modifier, 10, particle, mat);	
				for(int j = 0; j<no_of_el;j++) {
					mat = detector->GetMaterialM(0);
					path_to_exit = (pseudo_depth-p1e)/cos(exit_angle);					
					ResetMSParameters();
	    			M2 = mat0_M[j];
	    			recoil_st = RecoilEnergy(end_en,scat_angle,A1,M2);
	    			recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
	    			if(recoil_end > 100.*keV) {
	    				CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
						shape_m[j]		= GetShapeFactorM2();
						width_f[j]		= GetWidthFactorF2();	
						
						// layer 1
						mat = detector->GetMaterialM(1);
						path_to_exit = t1/cos(exit_angle);
						recoil_st = recoil_end;
						recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
						
						if(recoil_end > 100.*keV) {
	    					CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,1,M2);	
							shape_m[j]		= GetShapeFactorM2();
							width_f[j]		= GetWidthFactorF2();		
											
							// surface layer substrate, if it exists	
							if(surf_layer != 0) {
								mat = detector->GetMaterialM(0);
								path_to_exit = surf_layer/cos(exit_angle);	
								recoil_st = recoil_end;
								recoil_end = CalcTotEnLoss(recoil_end, path_to_exit,100, particle, mat);

								if(recoil_end > 100*keV) {																	
		    						CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
		    							
									if(use_ms_corrections == 1)
										width_f[j] 	= GetWidthFactorF2()*last_C_out;
									else
										width_f[j]	= GetWidthFactorF2();													
								} else { 
									shape_m[j]		= 0;
									width_f[j]		= 0;	
								}
							}
						} else { 
							shape_m[j]		= 0;
							width_f[j]		= 0;	
						}
					} else { 
						shape_m[j]		= 0;
						width_f[j]		= 0;	
					}
					myfile_shape_l1_l2 << shape_m[j] << "\t";
					myfile_width_l1_l2 << width_f[j] << "\t";	
				}					
				myfile_shape_l1_l2 << "\n";
				myfile_width_l1_l2 << "\n";
	  			break;  	  			
				// end of case 3
	  				    
			case 4: // particle in layer2
				// l2
				mat = detector->GetMaterialM(2);
				no_of_el = mat->GetNumberOfElements();	
				end_en = CalcTotEnLoss(start_en, step_interval/incidence_modifier, 10, particle, mat);	
				for(int j = 0; j<no_of_el;j++) {
					mat = detector->GetMaterialM(2);
					path_to_exit = (pseudo_depth-p2s)/cos(exit_angle);					
					ResetMSParameters();
	    			M2 = mat2_M[j];
	    			recoil_st = RecoilEnergy(end_en,scat_angle,A1,M2);
	    			recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
	    				
	    			if(recoil_end > 100.*keV) {
	    				CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,2,M2);	
						shape_m[j]		= GetShapeFactorM2();
						width_f[j]		= GetWidthFactorF2();	
						
						// substrate l2-l1
						mat = detector->GetMaterialM(0);
						path_to_exit = d12/cos(exit_angle);
						recoil_st = recoil_end;
						recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
						if(recoil_end > 100.*keV) {
	    					CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
							shape_m[j]		= GetShapeFactorM2();
							width_f[j]		= GetWidthFactorF2();	
						
							// l1
							mat = detector->GetMaterialM(1);
							path_to_exit = t1/cos(exit_angle);
							recoil_st = recoil_end;
							recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);							
							if(recoil_end > 100.*keV) {
	    						CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,1,M2);	
								shape_m[j]		= GetShapeFactorM2();
								width_f[j]		= GetWidthFactorF2();					
								// surface layer substrate, if it exists	
								if(surf_layer != 0) {
									mat = detector->GetMaterialM(0);
									path_to_exit = surf_layer/cos(exit_angle);	
									recoil_st = recoil_end;
									recoil_end = CalcTotEnLoss(recoil_end, path_to_exit,100, particle, mat);
									if(recoil_end > 100*keV) {												
		    							CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);		
										if(use_ms_corrections == 1)
											width_f[j] 	= GetWidthFactorF2()*last_C_out;
										else
											width_f[j]	= GetWidthFactorF2();														
									} else { 
										shape_m[j]		= 0;
										width_f[j]		= 0;	
									}
								}
							} else { 
								shape_m[j]		= 0;
								width_f[j]		= 0;	
							}
						} else { 
							shape_m[j]		= 0;
							width_f[j]		= 0;	
						}
					} else { 
						shape_m[j]		= 0;
						width_f[j]		= 0;	
					}
					myfile_shape_l2 << shape_m[j] << "\t";
					myfile_width_l2 << width_f[j] << "\t";
				}					
				myfile_shape_l2 << "\n";
				myfile_width_l2 << "\n";
	  			break;  	  			
				// end of case 4
	  
			case 5: // particle in between l3-l2
				// substrate
				mat = detector->GetMaterialM(0);
				no_of_el = mat->GetNumberOfElements();	
				end_en = CalcTotEnLoss(start_en, step_interval/incidence_modifier, 10, particle, mat);	
				for(int j = 0; j<no_of_el;j++) {
					mat = detector->GetMaterialM(0);
					path_to_exit = (pseudo_depth-p2e)/cos(exit_angle);					
					ResetMSParameters();
	    			M2 = mat0_M[j];
	    			recoil_st = RecoilEnergy(end_en,scat_angle,A1,M2);
	    			recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
	    				
	    			if(recoil_end > 100.*keV) { 
	    				CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
						shape_m[j]		= GetShapeFactorM2();
						width_f[j]		= GetWidthFactorF2();	
						
						// l2
						mat = detector->GetMaterialM(2);
						path_to_exit = t2/cos(exit_angle); 
						recoil_st = recoil_end;
						recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
						if(recoil_end > 100.*keV) {
    						CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,2,M2);	
							shape_m[j]		= GetShapeFactorM2();
							width_f[j]		= GetWidthFactorF2();	
						
							// l2-l1
							mat = detector->GetMaterialM(0);
							path_to_exit = d12/cos(exit_angle);
							recoil_st = recoil_end;
							recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);							
							if(recoil_end > 100.*keV) {
    							CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
								shape_m[j]		= GetShapeFactorM2();
								width_f[j]		= GetWidthFactorF2();		
											
								//l1
								mat = detector->GetMaterialM(1);
								path_to_exit = t1/cos(exit_angle);
								recoil_st = recoil_end;
								recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);											
								if(recoil_end > 100.*keV) {
    								CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,1,M2);	
									shape_m[j]		= GetShapeFactorM2();
									width_f[j]		= GetWidthFactorF2();									
									
									// surface layer substrate, if it exists	
									if(surf_layer != 0) {
										mat = detector->GetMaterialM(0);
										path_to_exit = surf_layer/cos(exit_angle);	
										recoil_st = recoil_end;
										recoil_end = CalcTotEnLoss(recoil_end, path_to_exit,100, particle, mat);
										if(recoil_end > 100*keV) {
											CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
		    									
											if(use_ms_corrections == 1)
												width_f[j] 	= GetWidthFactorF2()*last_C_out;
											else
												width_f[j]	= GetWidthFactorF2();													
										} else {
											shape_m[j]		= 0;
											width_f[j]		= 0;	
										}
									}
								} else {
									shape_m[j]		= 0;
									width_f[j]		= 0;	
								}	
							} else { 
								shape_m[j]		= 0;
								width_f[j]		= 0;	
							}						
						} else { 
							shape_m[j]		= 0;
							width_f[j]		= 0;	
						}
					} else { 
						shape_m[j]		= 0;
						width_f[j]		= 0;	
					}
					myfile_shape_l2_l3 << shape_m[j] << "\t";
					myfile_width_l2_l3 << width_f[j] << "\t";
				}					
				myfile_shape_l2_l3 << "\n";
				myfile_width_l2_l3 << "\n";
	  			break;  	  			
				// end of case 5	  				  
				
			case 6: // particle in layer3
				// l3
				mat = detector->GetMaterialM(3);
				no_of_el = mat->GetNumberOfElements();	
				end_en = CalcTotEnLoss(start_en, step_interval/incidence_modifier, 10, particle, mat);	
				for(int j = 0; j<no_of_el;j++) {
					mat = detector->GetMaterialM(3);
					path_to_exit = (pseudo_depth-p3s)/cos(exit_angle);					
					ResetMSParameters();
	    			M2 = mat3_M[j];
	    			recoil_st = RecoilEnergy(end_en,scat_angle,A1,M2);
	    			recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
	    				
	    			if(recoil_end > 100.*keV) {
	   					CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,3,M2);	
						shape_m[j]		= GetShapeFactorM2();
						width_f[j]		= GetWidthFactorF2();	
						
						// l3-l2
						mat = detector->GetMaterialM(0);
						path_to_exit = d23/cos(exit_angle);
						recoil_st = recoil_end;
						recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);					
						if(recoil_end > 100.*keV) {
    						CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
							shape_m[j]		= GetShapeFactorM2();
							width_f[j]		= GetWidthFactorF2();	
						
							// l2
							mat = detector->GetMaterialM(2);
							path_to_exit = t2/cos(exit_angle);
							recoil_st = recoil_end;
							recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);							
							if(recoil_end > 100.*keV) {
	 							CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,2,M2);	
								shape_m[j]		= GetShapeFactorM2();
								width_f[j]		= GetWidthFactorF2();		
											
								//l2-l1
								mat = detector->GetMaterialM(0);
								path_to_exit = d12/cos(exit_angle);
								recoil_st = recoil_end;
								recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);											
								if(recoil_end > 100.*keV) {
    								CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
									shape_m[j]		= GetShapeFactorM2();
									width_f[j]		= GetWidthFactorF2();									
									
									// l1
									mat = detector->GetMaterialM(1);
									path_to_exit = t1/cos(exit_angle);
									recoil_st = recoil_end;
									recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);											
									
									if(recoil_end > 100.*keV) {
	    								CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,1,M2);	
										shape_m[j]		= GetShapeFactorM2();
										width_f[j]		= GetWidthFactorF2();											
									
										// surface layer substrate, if it exists	
										if(surf_layer != 0) {
											mat = detector->GetMaterialM(0);
											path_to_exit = surf_layer/cos(exit_angle);	
											recoil_st = recoil_end;
											recoil_end = CalcTotEnLoss(recoil_end, path_to_exit,100, particle, mat);
											if(recoil_end > 100*keV) {																
		    									CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);		
												if(use_ms_corrections == 1)
													width_f[j] 	= GetWidthFactorF2()*last_C_out;
												else
													width_f[j]	= GetWidthFactorF2();													
											} else { 
												shape_m[j]		= 0;
												width_f[j]		= 0;	
											}
										}
									} else { 
										shape_m[j]		= 0;
										width_f[j]		= 0;	
									}											
								} else { 
									shape_m[j]		= 0;
									width_f[j]		= 0;	
								}	
							} else { 
								shape_m[j]		= 0;
								width_f[j]		= 0;	
							}						
						} else { 
							shape_m[j]		= 0;
							width_f[j]		= 0;	
						}
					} else { 
						shape_m[j]		= 0;
						width_f[j]		= 0;	
					}
					myfile_shape_l3 << shape_m[j] << "\t";
					myfile_width_l3 << width_f[j] << "\t";
				}					
				myfile_shape_l3 << "\n";
				myfile_width_l3 << "\n";
	  			break;  	  			
				// end of case 6	

			case 7: // particle between layers 4 and 3
				// substrate
				mat = detector->GetMaterialM(0);
				no_of_el = mat->GetNumberOfElements();	
				end_en = CalcTotEnLoss(start_en, step_interval/incidence_modifier, 10, particle, mat);	
				for(int j = 0; j<no_of_el;j++) {
					mat = detector->GetMaterialM(0);
					path_to_exit = (pseudo_depth-p3e)/cos(exit_angle);					
					ResetMSParameters();
	    			M2 = mat0_M[j];
	    			recoil_st = RecoilEnergy(end_en,scat_angle,A1,M2);
	    			recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
	    				
	    			if(recoil_end > 100.*keV) {
    					CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
						shape_m[j]		= GetShapeFactorM2();
						width_f[j]		= GetWidthFactorF2();	
						
						// l3
						mat = detector->GetMaterialM(3);
						path_to_exit = t3/cos(exit_angle);
						recoil_st = recoil_end;
						recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
						if(recoil_end > 100.*keV) {
    						CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,3,M2);	
							shape_m[j]		= GetShapeFactorM2();
							width_f[j]		= GetWidthFactorF2();	
						
							// l3-l2
							mat = detector->GetMaterialM(0);
							path_to_exit = d23/cos(exit_angle);
							recoil_st = recoil_end;
							recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);							
							if(recoil_end > 100.*keV) {
    							CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
								shape_m[j]		= GetShapeFactorM2();
								width_f[j]		= GetWidthFactorF2();		
											
								//l2
								mat = detector->GetMaterialM(2);
								path_to_exit = t2/cos(exit_angle);
								recoil_st = recoil_end;
								recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);											
								if(recoil_end > 100.*keV) {
    								CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,2,M2);	
									shape_m[j]		= GetShapeFactorM2();
									width_f[j]		= GetWidthFactorF2();									
									
									// l2-l1
									mat = detector->GetMaterialM(0);
									path_to_exit = d12/cos(exit_angle);
									recoil_st = recoil_end;
									recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);											
									if(recoil_end > 100.*keV) {
    									CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
										shape_m[j]		= GetShapeFactorM2();
										width_f[j]		= GetWidthFactorF2();											
									
										//l1
										mat = detector->GetMaterialM(1);
										path_to_exit = t1/cos(exit_angle);
										recoil_st = recoil_end;
										recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);										
										if(recoil_end > 100.*keV) {		
	    									CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,1,M2);	
											shape_m[j]		= GetShapeFactorM2();
											width_f[j]		= GetWidthFactorF2();										

											// surface layer substrate, if it exists	
											if(surf_layer != 0) {
												mat = detector->GetMaterialM(0);
												path_to_exit = surf_layer/cos(exit_angle);	
												recoil_st = recoil_end;
												recoil_end = CalcTotEnLoss(recoil_end, path_to_exit,100, particle, mat);
												if(recoil_end > 100*keV) {												
		    										CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
													if(use_ms_corrections == 1)
														width_f[j] 	= GetWidthFactorF2()*last_C_out;
													else
														width_f[j]	= GetWidthFactorF2();												
												} else { 
													shape_m[j]		= 0;
													width_f[j]		= 0;	
												}
											}
										} else { 
											shape_m[j]		= 0;
											width_f[j]		= 0;	
										}														
									} else { 
										shape_m[j]		= 0;
										width_f[j]		= 0;	
									}											
								} else {
									shape_m[j]		= 0;
									width_f[j]		= 0;	
								}	
							} else {
								shape_m[j]		= 0;
								width_f[j]		= 0;	
							}						
						} else {
							shape_m[j]		= 0;
							width_f[j]		= 0;	
						}
					} else {
						shape_m[j]		= 0;
						width_f[j]		= 0;	
					}
					myfile_shape_l3_l4 << shape_m[j] << "\t";
					myfile_width_l3_l4 << width_f[j] << "\t";
				}					
				myfile_shape_l3_l4 << "\n";
				myfile_width_l3_l4 << "\n";
	  			break;  	  			
				// end of case 7	

			case 8: // particle in layer4
				// l4
				mat = detector->GetMaterialM(4);
				no_of_el = mat->GetNumberOfElements();	
				end_en = CalcTotEnLoss(start_en, step_interval/incidence_modifier, 10, particle, mat);	
				for(int j = 0; j<no_of_el;j++) {
					mat = detector->GetMaterialM(4);
					path_to_exit = (pseudo_depth-p4s)/cos(exit_angle);					
					ResetMSParameters();
	    			M2 = mat4_M[j];
	    			recoil_st = RecoilEnergy(end_en,scat_angle,A1,M2);
	    			recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
	    			if(recoil_end > 100.*keV) {
	   					CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,4,M2);	
						shape_m[j]		= GetShapeFactorM2();
						width_f[j]		= GetWidthFactorF2();	
						
						// l4-l3
						mat = detector->GetMaterialM(0);
						path_to_exit = d34/cos(exit_angle);
						recoil_st = recoil_end;
						recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
						if(recoil_end > 100.*keV) {
    						CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
							shape_m[j]		= GetShapeFactorM2();
							width_f[j]		= GetWidthFactorF2();	
						
							// l3
							mat = detector->GetMaterialM(3);
							path_to_exit = t3/cos(exit_angle);
							recoil_st = recoil_end;
							recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);							
							if(recoil_end > 100.*keV) {
    							CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,3,M2);	
								shape_m[j]		= GetShapeFactorM2();
								width_f[j]		= GetWidthFactorF2();		
											
								//l3-l2
								mat = detector->GetMaterialM(0);
								path_to_exit = d23/cos(exit_angle);
								recoil_st = recoil_end;
								recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);											
								if(recoil_end > 100.*keV) {
    								CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
									shape_m[j]		= GetShapeFactorM2();
									width_f[j]		= GetWidthFactorF2();									
									
									// l2
									mat = detector->GetMaterialM(2);
									path_to_exit = t2/cos(exit_angle);
									recoil_st = recoil_end;
									recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);											
									if(recoil_end > 100.*keV) {
    									CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,2,M2);	
										shape_m[j]		= GetShapeFactorM2();
										width_f[j]		= GetWidthFactorF2();											
									
										//l2-l1
										mat = detector->GetMaterialM(0);
										path_to_exit = d12/cos(exit_angle);
										recoil_st = recoil_end;
										recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);										
										if(recoil_end > 100.*keV) {
    										CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
											shape_m[j]		= GetShapeFactorM2();
											width_f[j]		= GetWidthFactorF2();										

											//l1
											mat = detector->GetMaterialM(1);
											path_to_exit = t1/cos(exit_angle);
											recoil_st = recoil_end;
											recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
											if(recoil_end > 100.*keV) {
    											CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,1,M2);	
												shape_m[j]		= GetShapeFactorM2();
												width_f[j]		= GetWidthFactorF2();		
												// surface layer substrate, if it exists	
												if(surf_layer != 0)	{
													mat = detector->GetMaterialM(0);
													path_to_exit = surf_layer/cos(exit_angle);	
													recoil_st = recoil_end;
													recoil_end = CalcTotEnLoss(recoil_end, path_to_exit,100, particle, mat);
													if(recoil_end > 100*keV) {
		    											CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
														if(use_ms_corrections == 1)
															width_f[j] 	= GetWidthFactorF2()*last_C_out;
														else
															width_f[j]	= GetWidthFactorF2();													
													} else {
														shape_m[j]		= 0;
														width_f[j]		= 0;	
													}
												}
											} else {
												shape_m[j]		= 0;
												width_f[j]		= 0;	
											}													
										} else {
											shape_m[j]		= 0;
											width_f[j]		= 0;	
										}														
									} else {
										shape_m[j]		= 0;
										width_f[j]		= 0;	
									}											
								} else {
									shape_m[j]		= 0;
									width_f[j]		= 0;	
								}	
							} else {
								shape_m[j]		= 0;
								width_f[j]		= 0;	
							}						
						} else {
							shape_m[j]		= 0;
							width_f[j]		= 0;	
						}
					} else {
						shape_m[j]		= 0;
						width_f[j]		= 0;	
					}
					myfile_shape_l4 << shape_m[j] << "\t";
					myfile_width_l4 << width_f[j] << "\t";
				}					
				myfile_shape_l4 << "\n";
				myfile_width_l4 << "\n";
	  			break;  	  			
				// end of case 8	

			case 9: // particle after layer4
				// substrate
				mat = detector->GetMaterialM(0);
				no_of_el = mat->GetNumberOfElements();	
				end_en = CalcTotEnLoss(start_en, step_interval/incidence_modifier, 10, particle, mat);	
				for(int j = 0; j<no_of_el;j++) {
					mat = detector->GetMaterialM(0);
					path_to_exit = (pseudo_depth-p4s)/cos(exit_angle);					
					ResetMSParameters();
	    			M2 = mat0_M[j];
	    			recoil_st = RecoilEnergy(end_en,scat_angle,A1,M2);
	    			recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);				
	    			if(recoil_end > 100.*keV) {				
	    				CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
						shape_m[j]		= GetShapeFactorM2();
						width_f[j]		= GetWidthFactorF2();	
						// layer4				
						mat = detector->GetMaterialM(4);
						path_to_exit = (t4)/cos(exit_angle);					
	    				recoil_st = RecoilEnergy(end_en,scat_angle,A1,M2);
	    				recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
    					if(recoil_end > 100.*keV) {
    						CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,4,M2);	
							shape_m[j]		= GetShapeFactorM2();
							width_f[j]		= GetWidthFactorF2();	
							// l4-l3
							mat = detector->GetMaterialM(0);
							path_to_exit = d34/cos(exit_angle);
							recoil_st = recoil_end;
							recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
							if(recoil_end > 100.*keV) {
    							CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
								shape_m[j]		= GetShapeFactorM2();
								width_f[j]		= GetWidthFactorF2();	
								// l3
								mat = detector->GetMaterialM(3);
								path_to_exit = t3/cos(exit_angle);
								recoil_st = recoil_end;
								recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);							
								if(recoil_end > 100.*keV) {
    								CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,3,M2);	
									shape_m[j]		= GetShapeFactorM2();
									width_f[j]		= GetWidthFactorF2();		
									//l3-l2
									mat = detector->GetMaterialM(0);
									path_to_exit = d23/cos(exit_angle);
									recoil_st = recoil_end;
									recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);											
									if(recoil_end > 100.*keV) {
    									CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
										shape_m[j]		= GetShapeFactorM2();
										width_f[j]		= GetWidthFactorF2();									
										// l2
										mat = detector->GetMaterialM(2);
										path_to_exit = t2/cos(exit_angle);
										recoil_st = recoil_end;
										recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);											
										if(recoil_end > 100.*keV) {
    										CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,2,M2);	
											shape_m[j]		= GetShapeFactorM2();
											width_f[j]		= GetWidthFactorF2();											
											//l2-l1
											mat = detector->GetMaterialM(0);
											path_to_exit = d12/cos(exit_angle);
											recoil_st = recoil_end;
											recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);										
											if(recoil_end > 100.*keV) {
    											CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
												shape_m[j]		= GetShapeFactorM2();
												width_f[j]		= GetWidthFactorF2();									
												//l1
												mat = detector->GetMaterialM(1);
												path_to_exit = t1/cos(exit_angle);
												recoil_st = recoil_end;
												recoil_end = CalcTotEnLoss(recoil_st,path_to_exit, 100, particle, mat);
												if(recoil_end > 100.*keV) {
    												CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,1,M2);	
													shape_m[j]		= GetShapeFactorM2();
													width_f[j]		= GetWidthFactorF2();		
													// surface layer substrate, if it exists	
													if(surf_layer != 0) {
														mat = detector->GetMaterialM(0);
														path_to_exit = surf_layer/cos(exit_angle);	
														recoil_st = recoil_end;
														recoil_end = CalcTotEnLoss(recoil_end, path_to_exit,100, particle, mat);
														if(recoil_end > 100*keV) {
		    												CalculateMS_spread_out(end_en, recoil_st, particle, mat, path_to_exit,exit_angle, scat_angle,0,M2);	
															if(use_ms_corrections == 1)
																width_f[j] 	= GetWidthFactorF2()*last_C_out;
															else
																width_f[j]	= GetWidthFactorF2();													
														} else { 
															shape_m[j]		= 0;
															width_f[j]		= 0;	
														}
													}
												} else { 
													shape_m[j]		= 0;
													width_f[j]		= 0;	
												}													
											} else { 
												shape_m[j]		= 0;
												width_f[j]		= 0;	
											}														
										} else { 
											shape_m[j]		= 0;
											width_f[j]		= 0;	
										}											
									} else { 
										shape_m[j]		= 0;
										width_f[j]		= 0;	
									}	
								} else { 
									shape_m[j]		= 0;
									width_f[j]		= 0;	
								}						
							} else {
								shape_m[j]		= 0;
								width_f[j]		= 0;	
							}
						} else { 
							shape_m[j]		= 0;
							width_f[j]		= 0;	
						}
					} else {
						shape_m[j]		= 0;
						width_f[j]		= 0;	
					}	
					myfile_shape_l4_end << shape_m[j] << "\t";
					myfile_width_l4_end << width_f[j] << "\t";						
				}					
				myfile_shape_l4_end << "\n";
				myfile_width_l4_end << "\n";
	  			break;  	  			
				// end of case 9	
		}
		// sets the new energy to the end energy obtained after loop iteration
		start_en = end_en;
	}
	myfile_shape_surf.close();
	myfile_width_surf.close();
	myfile_shape_l1.close();
	myfile_width_l1.close();
	myfile_shape_l1_l2.close();
	myfile_width_l1_l2.close();						
	myfile_shape_l2.close();
	myfile_width_l2.close();										
	myfile_shape_l2_l3.close();
	myfile_width_l2_l3.close();						
	myfile_shape_l3.close();
	myfile_width_l3.close();
	myfile_shape_l3_l4.close();
	myfile_width_l3_l4.close();						
	myfile_shape_l4.close();
	myfile_width_l4.close();				
	myfile_shape_l4_end.close();
	myfile_width_l4_end.close();	
	
	G4cout << " Calculation and print to file is finished " << G4endl;	
	
	Construct2DVectors();
	
	G4cout << " Vector set up is finished " << G4endl;
	
	// deleting of arrays
	delete [] mat0_Z;
	delete [] mat1_Z;
	delete [] mat2_Z;
	delete [] mat3_Z;
	delete [] mat4_Z;		

	delete [] mat0_M;
	delete [] mat1_M;
	delete [] mat2_M;
	delete [] mat3_M;
	delete [] mat4_M;		

	return 1;
}

void EventAction::Construct2DVectors()
{
	G4cout << " Setting up MS 2D vectors " << G4endl;	
	MS_shape_surf = new G4Physics2DVector;
	MS_width_surf = new G4Physics2DVector;
	MS_shape_l1 = new G4Physics2DVector;
   	MS_width_l1 = new G4Physics2DVector;       
	MS_shape_l2 = new G4Physics2DVector;     
	MS_width_l2 = new G4Physics2DVector;       
	MS_shape_l3 = new G4Physics2DVector;     
	MS_width_l3 = new G4Physics2DVector;       
	MS_shape_l4 = new G4Physics2DVector;     
	MS_width_l4 = new G4Physics2DVector;            
	MS_shape_l12 = new G4Physics2DVector;     
	MS_width_l12 = new G4Physics2DVector;     
	MS_shape_l23 = new G4Physics2DVector;     
	MS_width_l23 = new G4Physics2DVector;                                                                     
	MS_shape_l34 = new G4Physics2DVector;     
	MS_width_l34 = new G4Physics2DVector;         
	MS_shape_l4e = new G4Physics2DVector;     
	MS_width_l4e = new G4Physics2DVector;   

	MS_shape_surf = ReadMSFiles("surf_temp_shape.txt");
	MS_shape_surf->SetBicubicInterpolation(true);
	
	MS_width_surf = ReadMSFiles("surf_temp_width.txt");
	MS_width_surf->SetBicubicInterpolation(true);		

	MS_shape_l1 = ReadMSFiles("l1_temp_shape.txt");
	MS_shape_l1->SetBicubicInterpolation(true);
	
	MS_width_l1 = ReadMSFiles("l1_temp_width.txt");
	MS_width_l1->SetBicubicInterpolation(true);
	
	MS_shape_l2 = ReadMSFiles("l2_temp_shape.txt");
	MS_shape_l2->SetBicubicInterpolation(true);
	
	MS_width_l2 = ReadMSFiles("l2_temp_width.txt");
	MS_width_l2->SetBicubicInterpolation(true);		

	MS_shape_l3 = ReadMSFiles("l3_temp_shape.txt");
	MS_shape_l3->SetBicubicInterpolation(true);
	
	MS_width_l3 = ReadMSFiles("l3_temp_width.txt");
	MS_width_l3->SetBicubicInterpolation(true);
	
	MS_shape_l4 = ReadMSFiles("l4_temp_shape.txt");
	MS_shape_l4->SetBicubicInterpolation(true);
	
	MS_width_l4 = ReadMSFiles("l4_temp_width.txt");
	MS_width_l4->SetBicubicInterpolation(true);		

	MS_shape_l12 = ReadMSFiles("l1-l2_temp_shape.txt");
	MS_shape_l12->SetBicubicInterpolation(true);
	
	MS_width_l12 = ReadMSFiles("l1-l2_temp_width.txt");
	MS_width_l12->SetBicubicInterpolation(true);
	
	MS_shape_l23 = ReadMSFiles("l2-l3_temp_shape.txt");
	MS_shape_l23->SetBicubicInterpolation(true);
	
	MS_width_l23 = ReadMSFiles("l2-l3_temp_width.txt");
	MS_width_l23->SetBicubicInterpolation(true);	
	
	MS_shape_l34 = ReadMSFiles("l3-l4_temp_shape.txt");
	MS_shape_l34->SetBicubicInterpolation(true);
	
	MS_width_l34 = ReadMSFiles("l3-l4_temp_width.txt");
	MS_width_l34->SetBicubicInterpolation(true);		
	
	MS_shape_l4e = ReadMSFiles("l4-end_temp_shape.txt");
	MS_shape_l4e->SetBicubicInterpolation(true);
	
	MS_width_l4e = ReadMSFiles("l4-end_temp_width.txt");
	MS_width_l4e->SetBicubicInterpolation(true);		
}

G4Physics2DVector* EventAction::ReadMSFiles(G4String name)
{
	G4Physics2DVector* vector = new G4Physics2DVector;
	G4String folder_def = "TEMP/";
	G4String tot_filename = folder_def+name;
	std::ifstream vFile(tot_filename);
	if(!vFile.is_open()) {
		G4cout << " No MS file present " << tot_filename << G4endl;
		exit(1);	
	}
	vector->Retrieve(vFile);
	return vector;
}
