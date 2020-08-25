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
#include "SensitiveDetectorHit.hh"

#include "Analysis.hh"

// new addition

#include "G4Material.hh"
#include "G4Track.hh"
#include "DetectorConstruction.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4EmCalculator.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"

#include <cmath>

#include "G4Material.hh"
#include "G4NistManager.hh"


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
 abs2StepPrim = 0;
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

G4double EventAction::RandomEnLoss(G4double E, G4double dedx, G4double position, G4double angle)
{
	// depth position in channeling system in cm
	//G4double z = position;
	//G4double energy = E;
	// path to outside according to scattering angle 
	G4double pathOutside = position/cos(angle);
	G4double energyLost = dedx * pathOutside/cm;
	G4double energyLeft = E - energyLost; 
	return energyLeft;
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


//G4double EventAction::CalcRBSYield(G4double xsec, G4double dist, G4double solidAngle, G4double atomDens, G4double angle)
G4double EventAction::CalcRBSYield(G4double xsec, G4double dist, G4double solidAngle, G4double atomDens)
{
	G4double crossSect = xsec * 1e-27; 	// mbarn to cm2 = 1e−27
	G4double z = dist/cm;
	G4double yield = (atomDens * crossSect * solidAngle * z);///cos(angle);
	return yield;
}

	//energy straggling in MeV^2
G4double EventAction::CalcBohrStrag(G4double Z1, G4double Z2, G4double atomDens, G4double dist)
{
	G4double elmcharge_squared = 1.4399764;	// in units of MeV*fm
	G4double e2 = std::pow(elmcharge_squared,2.)*1e-26;	// in units of MeV^2*cm^2
	G4double fourpi = 4*3.14159265358979;
	G4double fBohr = fourpi * e2;	
	//G4cout << " atom dens " << atomDens << " distance " << dist/cm << G4endl;
	G4double Bohr_t = Z1*Z1*Z2*fBohr*atomDens*(dist/cm);
	//G4cout << " bohr_t " << Bohr_t << G4endl;
	return Bohr_t;
}

	//https://doi.org/10.1016/0168-583X(85)90050-3
	// 50-150 degrees
/*
G4double EventAction::CalcRuthXsecMod(G4double xsec, G4double Z1, G4double Z2, G4double energy, G4double angle)
{
	G4double b = 3.3e-5*MeV;
	G4double c = 0.6;
	G4double result = xsec*(1- (b/energy)*Z1*Z2*std::pow(std::pow(Z1,2/3)+std::pow(Z2,2/3),0.5)*((1-c)+c/(sin(angle/2))));
	return result;
}
*/
G4double EventAction::CalcRuthXsecMod(G4double xsec, G4double Z1, G4double Z2, G4double energy)
{
	G4double a = 0.049*Z1*std::pow(Z2,4/3);
	G4double b = a/(energy/keV);
	G4double c = 1-b;
	//G4cout << " C " << c << G4endl;
	G4double result = xsec*c;
	return result;
}
	//from http://atlas.physics.arizona.edu/~shupe/Indep_Studies_2015/Notes_Goethe_Univ/L4_Scattering_Basic.pdf
G4double EventAction::CalcAngleCMFrame(G4double angle, G4double M1, G4double M2)
	{
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
	G4double multiplier = std::pow(1+(ratio*ratio)+(2*ratio*cos(angle)),1.5)/(1+(ratio*cos(angle)));
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


G4double EventAction::CalcStoppingNumb(G4double energy, G4double mass, G4double ioniz)
{
	G4double light_speed = 3e+10; // cm/s
	G4double light_speed_sq = light_speed*light_speed;
	G4double el_mass = 0.511; // in units of MeV/c^2 -> [MeV*s/m]
	G4double el_mass_2 = el_mass/(light_speed*light_speed); //[MeV*s^2/m^2]
	G4double part_mass = mass/(light_speed*light_speed); //[MeV*s^2/m^2]
	G4double velocity_sq = (2*(energy))/part_mass;	// in units of (cm/s)^2]
	G4double ion_const = ioniz/MeV;
	G4double beta_sq = velocity_sq/light_speed_sq;
	G4double stop_numb = std::log((2*el_mass_2*light_speed_sq*beta_sq)/(ion_const*(1-beta_sq)))-beta_sq;
	//G4cout << " stop nb " << stop_numb << G4endl;
	//not correct for ions, must be fixed later
	return std::abs(stop_numb);
}


G4double EventAction::CalcChi(G4double energy, G4double mass, G4double Z2)
{
	G4double elmcharge_squared = 1.4399764;	// in units of MeV*fm
	G4double e_sq = elmcharge_squared * 1e-7;	// in units of eV*cm
	G4double planck_bar = 6.58E-16; // in units of eV*s
	G4double Bohr_velocity = e_sq/planck_bar;
	G4double B_vel_sq = Bohr_velocity*Bohr_velocity;
	G4double light_speed = 3e+10; // cm/s
	G4double part_mass = mass/(light_speed*light_speed); //[MeV*s^2/m^2]
	G4double velocity_sq = (2*(energy))/part_mass;	// in units of (cm/s)^2]

	G4double chi_t = velocity_sq/(B_vel_sq*Z2);
	return chi_t;
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
               
  G4AnalysisManager::Instance()->FillH1(1,fTotalEnergyDeposit);
  G4AnalysisManager::Instance()->FillH1(3,fTotalEnergyFlow); 

  //run->AddEnergy(EnergyDeposit);
  //run->AddNonIonEnergy(NonIonEnergyDeposit);
  run->AddTrakLenPrim(TrakLenPrim);
  run->AddNumberOfSteps(nbStepsPrim); 
  run->AddTheta(theta);
  run->AddTrakLenSec(TrakLenSec);
  run->AddTrueRange(TrueTrakLen);
  run->AddProjRange(ProjTrakLen);
  run->absSTP(absStepPrim);
  run->abs2STP(abs2StepPrim);
  run->abs3STP(abs3StepPrim);
  // ******************************
  
  
  	for (int i = 0; i<5; i++)
		{
			Znumb[i] = new G4double[4];
			Mnumb[i] = new G4double[4];
			Adens[i] = new G4double[4];
		}
  
  
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	//functions to obtain parameters Z, M
	//particle parameters

	//detector energy resolution
	G4double det_FWHM = detector->GetDetectorResolution();
	G4int gauss_counter = detector->GetGaussCounter();
	//G4cout << " det fwhm " << det_FWHM/keV << G4endl;
	// check whether RBS evaluation is enabled 

	G4ParticleDefinition*  fParticle = fPrimary->GetParticleGun()->GetParticleDefinition();
    	G4double A1 = fPrimary->GetParticleGun()->GetParticleDefinition()->GetAtomicMass();
	//G4double A1_PDG = fPrimary->GetParticleGun()->GetParticleDefinition()->GetPDGMass();
    	G4double Z1 = fPrimary->GetParticleGun()->GetParticleDefinition()->GetAtomicNumber();
	//target atom parameters
	/*
    	G4Material*  material		= detector->GetMaterialM(0);
	G4Material*  material2 	= detector->GetMaterialM(1);
	*/
	
	//**************************************************************************
	//**************************************************************************
	//**************************************************************************
	//**************************************************************************
	
		// transliuoja gerai
	G4double angle_of_incidence = atan(fPrimary->GetParticleGun()->GetParticleMomentumDirection().x()/fPrimary->GetParticleGun()->GetParticleMomentumDirection().z());
	
	
	//G4cout << " " << angle_of_incidence/degree << G4endl;
	

	
	
	
	
	CalcHeavyIonStragglingFactor(1, Z1, 14);
	
	
	// test of new RBS code
	
		//scattering angle
	Angle = detector->GetRBSAngle();
	ang2 = 3.14159265358979; //pi radians
	G4double distance = detector->GetLength(0);
	//G4double distance2 = detector->GetLength(0);
	G4double solidAngle = detector->GetSolidAngle();
	//G4double solidAngle = 1.; // sterradian
	//G4double pathToOutside = 0.;
	//G4double pathToOutside2 = 0.;
  	G4EmCalculator emCalculator;
    	G4double step_number = detector->GetEnLossStep();
    	free_surf = 0.;
    	// detector dead layer
    	G4double dead_thickness 	= detector->GetDeadLayerThickness();
    	G4String dead_material_name 	= detector->GetDeadLayer();
    	
    	
    	//G4cout << " A1 " << Angle << G4endl;
    	Angle = Angle-angle_of_incidence;
	
	//G4cout << " A2 " << Angle << G4endl;
    	
    	
    	//G4cout << " dead mat name " << dead_material_name << G4endl;
    	
    	
    	const G4Material* dead_material = G4NistManager::Instance()->FindOrBuildMaterial(dead_material_name);
	
	
	
	for (int i=0;i<5;i++)
		{
		sample_material[i] 	= detector->GetMaterialM(i);
		NoOfElements[i]    	= sample_material[i]->GetNumberOfElements();
		}
	
	for (int i=0;i<5;i++)
		{ for (int j = 0; j<NoOfElements[i]; j++)
			{
				Znumb[i][j] 				= sample_material[i]->GetElement(j)->GetZ();
				Mnumb[i][j] 				= sample_material[i]->GetElement(j)->GetA()/(g/mole);
				const G4double *atomDensVector	= sample_material[i]->GetVecNbOfAtomsPerVolume();
				Adens[i][j] 				= atomDensVector[j]/(1/cm3);
				//G4cout << " el name " << sample_material[i]->GetElement(j)->GetName() << " density " << Adens[i][j] << G4endl;
			}
		}
	/*
	for (int i=0; i<NoOfElements[1];i++)
		{	
		G4cout << " atom density " << Adens[1][i] << " atom " << sample_material[1]->GetElement(i)->GetName() << G4endl;
		}
		*/
			// RTR vector fillup
	if (detector->GetSigmaCalc() == 1)
		{

  			int det_angle = (int)(detector->GetRBSAngle()/degree);
  			G4String det_angle_s = std::to_string(det_angle);
  			G4String start = "RTR_values/";
  			G4String mid = "_";
  			G4String end = ".txt";
  			
  			G4String part_designation;
  			
  			if (A1 == 1)
  				{ part_designation = "Hydrogen/H_";}
  			if (A1 == 4)
  				{part_designation = "Helium/He_";}
  			G4String mat_name;
  			G4String failas;
  				for(int i=0;i<5;i++)
  					{
  					for(int j=0;j<NoOfElements[i]; j++)
  						{
						mat_name = sample_material[i]->GetElement(j)->GetName();  	
						failas = start+part_designation+mat_name+mid+det_angle_s+end;
						  	if 	(mat_name == "Si" && fVectorSi == 0)
  								{ fVectorSi = FillRTRVector(failas); }
  							else if (mat_name == "O" && fVectorO == 0)
  								{ fVectorO = FillRTRVector(failas); }
   							else if (mat_name == "B" && fVectorB == 0)
  								{ fVectorB = FillRTRVector(failas); }
   							else if (mat_name == "C" && fVectorC == 0)
  								{ fVectorC = FillRTRVector(failas); }
   							else if (mat_name == "F" && fVectorF == 0)
  								{fVectorF = FillRTRVector(failas);  }
   							else if (mat_name == "N" && fVectorN == 0)
  								{fVectorN = FillRTRVector(failas); }
   							else if (mat_name == "Na" && fVectorNa == 0)
  								{fVectorNa = FillRTRVector(failas); }
   							else if (mat_name == "Mg" && fVectorMg == 0)
  								{fVectorMg = FillRTRVector(failas); }
   							else if (mat_name == "P" && fVectorP == 0)
  								{ fVectorP = FillRTRVector(failas); }
   							else if (mat_name == "S" && fVectorS == 0)
  								{ fVectorS = FillRTRVector(failas); }
   							else if (mat_name == "K" && fVectorK == 0)
  								{ fVectorK = FillRTRVector(failas); }
   							else if (mat_name == "Ca" && fVectorCa == 0)
  								{ fVectorCa = FillRTRVector(failas); }
   							else if (mat_name == "Ti" && fVectorTi == 0)
  								{ fVectorTi = FillRTRVector(failas); }
   							else if (mat_name == "Fe" && fVectorFe == 0)
  								{ fVectorFe = FillRTRVector(failas); }
   							else if (mat_name == "Cr" && fVectorCr == 0)
  								{ fVectorCr = FillRTRVector(failas); }
   							else if (mat_name == "Ni" && fVectorNi == 0)
  								{ fVectorNi = FillRTRVector(failas); }
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

  /*G4double layer_thickness[4]*/

  /*G4double *///layer_thickness[] = {detector->GetLength(1), detector->GetLength(2), detector->GetLength(3), detector->GetLength(4)};  
  for (int i=0; i<4; i++)
  	{
  		layer_thickness[i] 	= detector->GetLength(i+1);
  		layer_pos[i]		= detector->GetPosition(i);
  		layer_pos_z[i]		= layer_pos[i].z();
	}
   for (int i=0; i< 3; i++)
   	{
   	  		neighbour_dist[i]	=(layer_pos_z[i+1]-layer_thickness[i+1]/2)-(layer_pos_z[i]+layer_thickness[i]/2);
   			if (neighbour_dist[i]/um < 0.0001) {neighbour_dist[i] =0.;}
   	}
  		
  /*
  G4ThreeVector layer_pos[4] = {detector->GetPosition(0),detector->GetPosition(1),detector->GetPosition(2),detector->GetPosition(3)};
  G4double layer_pos_z[4] = {layer_pos[0].z(), layer_pos[1].z(), layer_pos[2].z(), layer_pos[3].z()};
  
  G4double dist_all[4] = {0}; 
  for (int i=0; i<4; i++)
  	{
  		dist_all[i] = -(detector->GetLength(0)/2)-layer_pos_z[i]+layer_thickness[i]/2; 
  		//G4cout << " dist " << dist_all[i]/um << G4endl;
  	}
  	
   G4double neighbour_dist[3];
   for (int i=0; i<3; i++)
   	{
   		//neighbour_dist[i] = std::abs(layer_pos_z[i]-layer_pos_z[i+1])+(layer_thickness[i]/2)+(layer_thickness[i+1]/2);
   		neighbour_dist[i] = (layer_pos_z[i+1]-layer_thickness[i+1]/2)-(layer_pos_z[i]+layer_thickness[i]/2);
   		if (neighbour_dist[i]/um < 0.0001) {neighbour_dist[i] =0.;}
   		//G4cout << " atstumas " << i << " " << neighbour_dist[i]/um << G4endl;	
   		
   	}
   	*/
   	
   	
	G4double surf_distance = detector->GetLength(0)/2+layer_pos_z[0]-(detector->GetLength(1)/2);
	if (surf_distance/um < 0.0001) { surf_distance = 0.;}
	//G4cout << " surf dist " << surf_distance/um << G4endl;

	G4double steps;
	G4ThreeVector position;
	G4double sample_energy =0.;
	G4double xsecRTR;
	G4double RBS_yield =0.;
	G4double path =0.;
	G4double energy_left = 0.;
	G4double depth =0.;
	
	G4double trueWorldPosition = 0.;
	G4double bohr_straggling = 0.;
	free_surf = 0;
	G4ThreeVector worldMomentum;
	G4double transv_momentum;
	G4double momentum_z;
	
	G4double tot_step =0.;
	
	
	
	
	

	
	


	// mother volume sensitive detector
    if(sd0){
        int n_hit_sd = sd0->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++){
            	CrystalDetectorHit* aHit = (*sd0)[i1];
            	steps = aHit->GetStep();	
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();
	
	    	
	    	trueWorldPosition = position.z()+detector->GetLength(0)/2;
	    	
	    	for(int i=0;i<NoOfElements[0];i++)
			{
	    			G4double RecEn = RecoilEnergy(sample_energy, Angle, A1, Mnumb[0][i]);
	    			//G4cout << " rec E subs surf " << RecEn/MeV << G4endl;
	    			G4String el_name = sample_material[0]->GetElement(i)->GetName();
	    			if (detector->GetSigmaCalc() == 1)
					{
						xsecRTR = GetRTRValue(sample_energy, el_name);
	    				}
	    			else { xsecRTR = 1.;}
	    			
	    			
	    			
	    			// xsec obtain
	    			
	    			/*
	    			G4double CMangleX	= CalcAngleCMFrame(Angle,A1,Mnumb[0][i]);
	    			G4double CMenergyX 	= CalcEnergyCMFrame(sample_energy,A1,Mnumb[0][i]);
				G4double CMxsecX 	= CalcDiffRuthXsecCM(CMenergyX,CMangleX,Z1,Znumb[0][i]);
				G4double CMxsecModX	= CalcRuthXsecMod(CMxsecX,Z1,Znumb[0][i],CMenergyX);
				G4double CMtoLABxsecX	= CalcDiffRuthXsecLAB(A1,Mnumb[0][i],CMangleX,CMxsecModX);	// cia buvo ModX
	    			
	    			
	    			
	    			
	    		
	    				if (i == 1){
	    				//G4cout << " CM angle " << CMangleX/degree << G4endl;
					//G4cout << " CM ANGLE " << CMangle[i]/degree << G4endl;
					if (sample_energy/MeV > 1.8) {
					//G4cout << " energy [MeV] " << sample_energy/MeV << " xsec " << (CMtoLABxsecX*xsecRTR) << G4endl;
					
					
					
					 //FILE * output=fopen("XSEC_O.txt","a");
					 //FILE * output=fopen("XSEC_Si_LAB.txt","a");
					  FILE * output=fopen("XSEC_O_LAB.txt","a");
  
  					 fprintf(output, "%f\t %f\n", sample_energy/MeV, (CMtoLABxsecX*xsecRTR));
   					//flush(output);
   					fclose(output); 
					
					
					
					
					
					}
					}
	*/

	
				// end of xsec obtain
	    			
	    			
	    			
	    			
	    			
	    			RBS_yield = CalculateTotalRBSYield(sample_energy, A1, Mnumb[0][i], Z1, Znumb[0][i], Angle, distance, solidAngle, xsecRTR, Adens[0][i]);
				//G4cout << " RBS Yield subs " << RBS_yield << G4endl;
//******************************************************
// surface region	    			
//******************************************************
	    			// energy loss if particle is in mother material, and first layer is below the current particle position
	    			if (position.z() < layer_pos_z[0])
	    				{
	    				//G4cout << " free surf " << G4endl;
	    				// set variable free_surf to use in energy loss evaluations from other layers
	    				SetFreeSurf(1); 
	    				//G4cout << " pos surf " << trueWorldPosition/um << G4endl;
	    				path = trueWorldPosition/(cos(ang2-Angle));
	    				depth = trueWorldPosition;
	    				energy_left = CalcTotEnLoss(RecEn, path, step_number, fParticle, sample_material[0]);
	    				bohr_straggling = CalculateTotalBohrStraggling(RecEn, fParticle, sample_material[0], path);
	    				if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    				{
	    					energy_left = 0.;
	    					run->MaxRBSDepth(trueWorldPosition);
	    					run->AddCount();
	    				}
	    			}
//******************************************************
// region between l1 and l2	    			
//******************************************************
	    			
	    			// energy loss if particle is in mother material, between layers 1 and 2
	    			if (position.z() > layer_pos_z[0] && position.z() < layer_pos_z[1])
	    			{	
	    				//G4cout << " surface not free " << G4endl;
	    				// calculate energy loss to layer1
	    				G4double dist_to_layer1 = std::abs((layer_pos_z[0] + layer_thickness[0]/2) - position.z());
	    				//G4cout << " dist to layer " << dist_to_layer1/um << G4endl;
	    				depth = dist_to_layer1;
					path = dist_to_layer1/(cos(ang2-Angle));
					energy_left = CalcTotEnLoss(RecEn, path, step_number, fParticle, sample_material[0]);
					bohr_straggling = CalculateTotalBohrStraggling(RecEn, fParticle, sample_material[0], path);
					
					// calculate energy loss in layer 1
	    				if (energy_left/MeV > 0.01)
	    				{
	    					// layer 1 thickness
	    					path = layer_thickness[0]/(cos(ang2-Angle));
	    					depth += layer_thickness[0];
	    					energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[1]);
	    					bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[1], path);
	    					
	    					if (energy_left/MeV > 0.01){
	    						if (GetFreeSurf() == 1)
	    						{
	    							// calculate en loss in surface, after layer1
	    							depth += surf_distance;
	    							path = surf_distance/(cos(ang2-Angle));
	    							energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    							bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    					
	    						}
	    					}
	    					else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    					{
	    						energy_left = 0.;
	    						run->MaxRBSDepth(trueWorldPosition);
	    						run->AddCount();
	    					}
	    					
	    				}
	    				else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    				{
	    					energy_left = 0.;
	    					run->MaxRBSDepth(trueWorldPosition);
	    					run->AddCount();
	    				}	
	    			}
//******************************************************
// region between l2 and l3	    			
//******************************************************
	    			// energy loss if particle is in mother material, between layers 2 and 3
				if (position.z() > layer_pos_z[1] && position.z() < layer_pos_z[2]) 
				{
	    				// calculate energy loss to layer2
	    				G4double dist_to_layer2 = std::abs((layer_pos_z[1] + layer_thickness[1]/2) - position.z());
	    				//G4cout << " dist to l2 " << dist_to_layer2/um << " pos " << position.z()/um << G4endl;
	    				depth = dist_to_layer2;
	    				path = dist_to_layer2/(cos(ang2-Angle));
					energy_left = CalcTotEnLoss(RecEn, path, step_number, fParticle, sample_material[0]);
					bohr_straggling = CalculateTotalBohrStraggling(RecEn, fParticle, sample_material[0], path);
				
					// calculate energy loss in layer 2
	    				if (energy_left/MeV > 0.01)
	    				{
	    					// layer 2 thickness
	    					path = layer_thickness[1]/(cos(ang2-Angle));
	    					depth += layer_thickness[1];
	    					energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[2]);
	    					bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[2], path);
	    					
	    					if (energy_left/MeV > 0.01)
	    					{
	    						// calculate en loss from l2 to l1 
	    						depth += neighbour_dist[0];
	    						//G4cout << " neig " << neighbour_dist[0]/um << G4endl;
	    						path = neighbour_dist[0]/(cos(ang2-Angle));
	    						energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    						bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    						
	    						// calculate en loss in l1
	    						if(energy_left/MeV > 0.01)
	    						{
	    							
	    							path = layer_thickness[0]/(cos(ang2-Angle));
	    							depth += layer_thickness[0];
	    							energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[1]);
	    							bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[1], path);
	    							
	    							// calculate en loss in surface, after layer1
	    							if (energy_left/MeV > 0.01) {
	    								if (GetFreeSurf() == 1 )
	    								{
	    									depth += surf_distance;
	    									path = surf_distance/(cos(ang2-Angle));
	    									energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    									bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    								}
	    							}
	    							else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    							{
	    								energy_left = 0.;
	    								run->MaxRBSDepth(trueWorldPosition);
	    								run->AddCount();
	    							}
	    						
	    						}
	    						else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    						{
	    							energy_left = 0.;
	    							run->MaxRBSDepth(trueWorldPosition);
	    							run->AddCount();
	    						}
	    					
	    					}
	    					else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    					{
	    						energy_left = 0.;
	    						run->MaxRBSDepth(trueWorldPosition);
	    						run->AddCount();
	    					}
	    					
	    				}
	    				else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    				{
	    					energy_left = 0.;
	    					run->MaxRBSDepth(trueWorldPosition);
	    					run->AddCount();
	    				}	
				
				}
				// end of region
				
//******************************************************
// region between l3 and l4	    			
//******************************************************
				// energy loss if particle is in mother material, between layers 3 and 4
				if (position.z() > layer_pos_z[2] && position.z() < layer_pos_z[3]) 
				{
					// calculate energy loss to layer3
	    				G4double dist_to_layer3 = std::abs((layer_pos_z[2] + layer_thickness[2]/2) - position.z());
	    				//G4cout << " dist to l3 " << dist_to_layer3/um << " pos " << position.z()/um << G4endl;
	    				depth = dist_to_layer3;
	    				path = dist_to_layer3/(cos(ang2-Angle));
					energy_left = CalcTotEnLoss(RecEn, path, step_number, fParticle, sample_material[0]);
					bohr_straggling = CalculateTotalBohrStraggling(RecEn, fParticle, sample_material[0], path);
					// calculate energy loss in layer 3
	    				if (energy_left/MeV > 0.01)
	    				{
	    					// layer 3 thickness
	    					path = layer_thickness[2]/(cos(ang2-Angle));;
	    					depth += layer_thickness[2];
	    					energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[3]);
	    					bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[3], path);
	    					
	    					// calculate energy loss between layers 2 and 3 
	    					if (energy_left/MeV > 0.01)
	    					{
	    						depth += neighbour_dist[1];
	    						path = neighbour_dist[1]/(cos(ang2-Angle));;
	    						energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    						bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    						
	    						// calculate energy loss in layer 2
	    						if (energy_left/MeV > 0.01)
	    						{
	    							// layer 2 thickness
	    							path = layer_thickness[1]/(cos(ang2-Angle));;
	    							depth += layer_thickness[1];
	    							energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[2]);
	    							bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[2], path);
	    							
	    							// calculate energy loss in between layers 2 and 1 
	    							if (energy_left/MeV > 0.01)
	    							{
	    								depth += neighbour_dist[0];
	    								path = neighbour_dist[0]/(cos(ang2-Angle));;
	    								energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    								bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    								
	    								// calculate energy loss in layer 1
	    								if(energy_left/MeV > 0.01)
	    								{
	    									path = layer_thickness[0]/(cos(ang2-Angle));;
	    									depth += layer_thickness[0];
	    									energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[1]);
	    									bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[1], path);
	    							
	    									// calculate en loss in surface, after layer1
	    									if (GetFreeSurf() == 1 && energy_left/MeV > 0.01)
	    									{
	    										depth += surf_distance;
	    										path = surf_distance/(cos(ang2-Angle));;
	    										energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    										bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    										
	    										// kolkas blogas, reikia tikrinti.
	    										//G4cout << " palyginimas, world " << trueWorldPosition/um << " suma " << depth/um << G4endl;
	    					
	    									}
	    									else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    									{
	    										run->MaxRBSDepth(trueWorldPosition);
	    										run->AddCount();
	    									}
	    								}
	    								else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    								{
	    									run->MaxRBSDepth(trueWorldPosition);
	    									run->AddCount();
	    								}
	    							}
	    							else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    							{
	    								run->MaxRBSDepth(trueWorldPosition);
	    								run->AddCount();
	    							}
	    						}
	    						else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    						{
	    							run->MaxRBSDepth(trueWorldPosition);
	    							run->AddCount();
	    						}
	    					}	
	    					else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    					{
	    						run->MaxRBSDepth(trueWorldPosition);
	    						run->AddCount();
	    					}
	    				}
					else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    				{
	    					run->MaxRBSDepth(trueWorldPosition);
	    					run->AddCount();
	    				}
				
				} // end of region
				
				
//******************************************************
// region after layer 4	    			
//******************************************************

				// energy loss if particle is in mother material, after layer4
				if (position.z() > layer_pos_z[3]) 
				{
					// calculate energy loss to layer4
	    				G4double dist_to_layer4 = std::abs((layer_pos_z[3] + layer_thickness[3]/2) - position.z());
	    				depth = dist_to_layer4;
	    				path = dist_to_layer4/(cos(ang2-Angle));
					energy_left = CalcTotEnLoss(RecEn, path, step_number, fParticle, sample_material[0]);
					bohr_straggling = CalculateTotalBohrStraggling(RecEn, fParticle, sample_material[0], path);	
					
					if(energy_left/MeV > 0.01)	
					// calculate energy loss in layer 4
					{
	    					path = layer_thickness[3]/(cos(ang2-Angle));;
	    					depth += layer_thickness[3];
	    					energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[4]);
	    					bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[4], path);
						// calculate energy loss in layer 3

						if(energy_left/MeV > 0.01)
						{	
							// calculate energy loss to layer3
	    						G4double dist_to_layer3 = std::abs((layer_pos_z[2] + layer_thickness[2]/2) - position.z());
	    						//G4cout << " dist to l3 " << dist_to_layer3/um << " pos " << position.z()/um << G4endl;
	    						depth = dist_to_layer3;
	    						path = dist_to_layer3/(cos(ang2-Angle));
							energy_left = CalcTotEnLoss(RecEn, path, step_number, fParticle, sample_material[0]);
							bohr_straggling = CalculateTotalBohrStraggling(RecEn, fParticle, sample_material[0], path);
							// calculate energy loss in layer 3
	    						if (energy_left/MeV > 0.01)
	    						{
	    							// layer 3 thickness
	    							path = layer_thickness[2]/(cos(ang2-Angle));;
	    							depth += layer_thickness[2];
	    							energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[3]);
	    							bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[3], path);
	    							
	    							// calculate energy loss between layers 2 and 3 
	    							if (energy_left/MeV > 0.01)
	    							{
	    								depth += neighbour_dist[1];
	    								path = neighbour_dist[1]/(cos(ang2-Angle));;
	    								energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    								bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    								
	    								// calculate energy loss in layer 2
	    								if (energy_left/MeV > 0.01)
	    								{
	    									// layer 2 thickness
	    									path = layer_thickness[1]/(cos(ang2-Angle));;
	    									depth += layer_thickness[1];
	    									energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[2]);
	    									bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[2], path);
	    									
	    									// calculate energy loss in between layers 2 and 1 
	    									if (energy_left/MeV > 0.01)
	    									{
	    										depth += neighbour_dist[0];
	    										path = neighbour_dist[0]/(cos(ang2-Angle));;
	    										energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    										bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    										
	    										// calculate energy loss in layer 1
	    										if(energy_left/MeV > 0.01)
	    										{
	    											path = layer_thickness[0]/(cos(ang2-Angle));;
	    											depth += layer_thickness[0];
	    											energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[1]);
	    											bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[1], path);
	    									
	    											// calculate en loss in surface, after layer1
	    											if (GetFreeSurf() == 1 && energy_left/MeV > 0.01)
	    											{
	    												depth += surf_distance;
	    												path = surf_distance/(cos(ang2-Angle));;
	    												energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    												bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    												
	    												// kolkas blogas, reikia tikrinti.
	    												//G4cout << " palyginimas, world " << trueWorldPosition/um << " suma " << depth/um << G4endl;
	    							
	    											}
	    											else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    											{
	    												run->MaxRBSDepth(trueWorldPosition);
	    												run->AddCount();
	    											}
	    										}
	    										else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	   	 									{
	    											run->MaxRBSDepth(trueWorldPosition);
	    											run->AddCount();
	    										}
	    									}
	    									else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    									{
	    										run->MaxRBSDepth(trueWorldPosition);
	   	 									run->AddCount();
	    									}
	    								}
	    								else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	   	 							{
	    									run->MaxRBSDepth(trueWorldPosition);
	    									run->AddCount();
	    								}
	    							}	
	    							else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    							{
	    								run->MaxRBSDepth(trueWorldPosition);
	    								run->AddCount();
	    							}
	    						}
							else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    						{
	    							run->MaxRBSDepth(trueWorldPosition);
	    							run->AddCount();
	    						}
						}
						else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    					{
	    						run->MaxRBSDepth(trueWorldPosition);
	    						run->AddCount();
	    					}
					}
					
				}	//end of region
				
				// filling of RBS histograms
				if (energy_left/MeV > 0.01)
					{
					
					tot_step +=steps;

					// include detector dead layer energy loss
					energy_left = CalculateDeadLayerEffect(energy_left, dead_material, dead_thickness, fParticle);

					
					G4double rbs_depth = depth;// - detector->GetLength(0)/2;
					//G4cout << " rbs depth " << rbs_depth/um << G4endl;
					
					
					// rbs dist 
					analysisManager->FillH1(21, rbs_depth, RBS_yield*steps); // total
					analysisManager->FillH1(23, rbs_depth, RBS_yield*steps); // substrate
					// el1
					if (i == 0) 
						{
						analysisManager->FillH1(26, rbs_depth, RBS_yield*steps); //el1
						}
					// el2
					if (i == 1)
						{
						analysisManager->FillH1(27, rbs_depth, RBS_yield*steps); // el2
						}
					
					G4double sigma_det = (det_FWHM/MeV)/2.355;
					G4double sigma_det_sq = std::pow(sigma_det,2);
					
					G4double tot_sigma = sigma_det_sq + bohr_straggling;
					
					
					//G4cout << " Boras " << bohr_straggling << G4endl;
				
					G4int GK = gauss_counter;
					G4double Gauss_sum	= 0;
					for (int k=-GK; k<=GK; k++)	
						{
							G4double santykis	= k;
							G4double Gauss_value 	= (GenerateGaussian(energy_left*(1-(santykis/1000)),energy_left,tot_sigma))/GenerateGaussian(energy_left,energy_left,tot_sigma);
							Gauss_sum += Gauss_value;
						}
	
					for (int k=-GK; k<=GK; k++)	
						{
							G4double santykis	= k;
							G4double Gauss_value 	= (GenerateGaussian(energy_left*(1-(santykis/1000)),energy_left,tot_sigma))/GenerateGaussian(energy_left,energy_left,tot_sigma);
							//G4double mod_G_value 	= Gauss_value/Gauss_sum;			// original Gauss_value/Gauss_sum;
							G4double mod_G_value 	= (Gauss_value/Gauss_sum)/(2*GK+1);
							G4double Gauss_en 	= energy_left*(1-(santykis/1000)); 
							
							//G4cout << " G energy " << Gauss_en/MeV << G4endl;
							//G4cout << " mod G value " << mod_G_value << G4endl;
							//G4cout << " Rbs yield " << RBS_yield << G4endl;
							
							
							//G4cout << " G; k= " << k << " verte " << mod_G_value << G4endl;
							
							if (Gauss_en/MeV > 0.01)
								{
								//G4cout << " RBS YIELD SUB " << RBS_yield*mod_G_value*steps << G4endl;
								analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_G_value*steps));	// total
								analysisManager->FillH1(22, Gauss_en, (RBS_yield*mod_G_value*steps));	//substrate
								if (i == 0) 
									{ 
										analysisManager->FillH1(24,Gauss_en, (RBS_yield*mod_G_value*steps));	//el1
									}
								if (i == 1) 
									{ 
										analysisManager->FillH1(25,Gauss_en, (RBS_yield*mod_G_value*steps)); //el2
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
        int n_hit_sd = sd1->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++){
            	CrystalDetectorHit* aHit = (*sd1)[i1];
            	steps = aHit->GetStep();	
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();
	    	
	    	//G4cout << " hit obtained " << G4endl;
	
	    	
	    	
	    	G4double layer_to_surf = std::abs(-(detector->GetLength(0)/2)-layer_pos_z[0]+(detector->GetLength(1)/2));
	    	if (layer_to_surf/um < 0.0001)
	    		{
	    		layer_to_surf = 0.;
	    		}
	    	// position in the sample
	    	trueWorldPosition = position.z()+detector->GetLength(0)/2;
	    	//G4cout << " true world position " << trueWorldPosition/um << G4endl;
	    	// position in the layer
	    	G4double position_in_layer = trueWorldPosition - layer_to_surf;
	    	//G4cout << " pos " << position_in_layer/um << G4endl;
 
	    	
	    	//G4cout << " pos in layer " << position_in_layer/um << G4endl;
	    	
	    	for(int i=0;i<NoOfElements[1];i++)
			{
	    			G4double RecEn = RecoilEnergy(sample_energy, Angle, A1, Mnumb[1][i]);
	    			G4String el_name = sample_material[1]->GetElement(i)->GetName();
	    			if (detector->GetSigmaCalc() == 1)
					{
						xsecRTR = GetRTRValue(sample_energy, el_name);
	    				}
	    			else { xsecRTR = 1.;}
	    			RBS_yield = CalculateTotalRBSYield(sample_energy, A1, Mnumb[1][i], Z1, Znumb[1][i], Angle, distance, solidAngle, xsecRTR, Adens[1][i]); 
	    			//G4cout << " RB Y " << RBS_yield << G4endl;
	    			//G4cout << " layer 1 RBS Y " << RBS_yield << G4endl;
//******************************************************
// surface region	    			
//******************************************************
	    			// energy loss if layer1 is on the surface
	    			if (surf_distance == 0)
	    			{ 
	    				// GOOD
	    				//G4cout << " no free surface " << G4endl;
	    				path = position_in_layer/(cos(ang2-Angle));
	    				depth = position_in_layer;
	    				energy_left = CalcTotEnLoss(RecEn, path, step_number, fParticle, sample_material[1]);
	    				bohr_straggling = CalculateTotalBohrStraggling(RecEn, fParticle, sample_material[1],path);
	    				
	    				if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    				{
	    					energy_left = 0.;
	    					run->MaxRBSDepth(trueWorldPosition);
	    					run->AddCount();
	    				}
	    			}
	    			// energy loss if layer1 is below the surface
	    			if (surf_distance != 0)
	    			{
	    				//G4cout << " free surface " << G4endl;
					//G4cout << " true pos " << trueWorldPosition/nm << G4endl;
					

	    				path = position_in_layer/(cos(ang2-Angle));
	    				depth = position_in_layer;
	    				energy_left = CalcTotEnLoss(RecEn, path, step_number, fParticle, sample_material[1]);
	    				bohr_straggling = CalculateTotalBohrStraggling(RecEn, fParticle, sample_material[1], path);
	    				
	    				//if (path < 0 || path/um > 0.32) { G4cout << " xujovas atstumas " << G4endl;}
	    				if (energy_left/MeV > 0.01)
	    				{
	    					//if (energy_left < 1*MeV) { G4cout << " xuine " << " " << i << G4endl;}
	    					path = surf_distance/(cos(ang2-Angle));
	    					depth += surf_distance;
	    					energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    					bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    				}
	    				
	    				else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    				{
	    					energy_left = 0.;
	    					run->MaxRBSDepth(trueWorldPosition);
	    					run->AddCount();
	    				}
				}
				
				// filling of RBS histograms
				if (energy_left/MeV > 0.01)
					{
					
					// include detector dead layer energy loss
					energy_left = CalculateDeadLayerEffect(energy_left, dead_material,dead_thickness,fParticle);
					
					
					//G4cout << " boras " << bohr_straggling << G4endl;
					G4double rbs_depth = depth;
					// rbs dist 
					analysisManager->FillH1(21, rbs_depth, RBS_yield*steps); // total
					analysisManager->FillH1(29, rbs_depth, RBS_yield*steps); // layer1
					// el1
					if (i == 0) 
						{
						analysisManager->FillH1(32, rbs_depth, RBS_yield*steps); //el1
						}
					// el2
					if (i == 1)
						{
						analysisManager->FillH1(33, rbs_depth, RBS_yield*steps); // el2
						}
					
					G4double sigma_det = (det_FWHM/MeV)/2.355;
					G4double sigma_det_sq = std::pow(sigma_det,2);
					G4double tot_sigma = sigma_det_sq + bohr_straggling;
					
					
					//G4cout << " en left " << energy_left/MeV << "  " << i <<  G4endl;
					
					
					G4int GK = gauss_counter;
					G4double Gauss_sum	= 0;
					for (int k=-GK; k<=GK; k++)	
						{
							G4double santykis	= k;
							G4double Gauss_value 	= (GenerateGaussian(energy_left*(1-(santykis/1000)),energy_left,tot_sigma))/GenerateGaussian(energy_left,energy_left,tot_sigma);
							Gauss_sum += Gauss_value;
						}
	
					for (int k=-GK; k<=GK; k++)	
						{
							G4double santykis	= k;
							G4double Gauss_value 	= (GenerateGaussian(energy_left*(1-(santykis/1000)),energy_left,tot_sigma))/GenerateGaussian(energy_left,energy_left,tot_sigma);
							//G4double mod_G_value 	= Gauss_value/Gauss_sum;
							G4double mod_G_value 	= (Gauss_value/Gauss_sum)/(2*GK+1);
							G4double Gauss_en 	= energy_left*(1-(santykis/1000)); 
							if (Gauss_en/MeV > 0.01)
								{
								
								//G4cout << " RBS YIELD L 1 " << RBS_yield*mod_G_value*steps << G4endl;
								analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_G_value*steps));	// total
								analysisManager->FillH1(28, Gauss_en, (RBS_yield*mod_G_value*steps));	//layer1
								if (i == 0) 
									{ 
										analysisManager->FillH1(30,Gauss_en, (RBS_yield*mod_G_value*steps));	//el1
									}
								if (i == 1) 
									{ 
										analysisManager->FillH1(31,Gauss_en, (RBS_yield*mod_G_value*steps)); //el2
									}
								}
						}					
				
					}

    			}//end of elements
    		}		run->add_total_step(tot_step);
    	}// end of layer1 sensitive detector
    
//======================================================
// LAYER 2 	    			
//======================================================

	// layer2 sensitive detector
    if(sd2){
        int n_hit_sd = sd2->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++){
            	CrystalDetectorHit* aHit = (*sd2)[i1];
            	steps = aHit->GetStep();	
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();
	    	
	    	
	    	
	    	
	    	
	    	    
	    	G4double layer_to_surf = std::abs(-(detector->GetLength(0)/2)-layer_pos_z[1]+(detector->GetLength(2)/2));	    	    
		//G4cout << " D " << layer_to_surf/um << G4endl;

	    	trueWorldPosition = (position.z()+detector->GetLength(0)/2)-layer_to_surf;
	    	//G4cout << " true world position " << trueWorldPosition/um << G4endl;

	    	for(int i=0;i<NoOfElements[2];i++)
			{
	    			G4double RecEn = RecoilEnergy(sample_energy, Angle, A1, Mnumb[2][i]);
	    			G4String el_name = sample_material[2]->GetElement(i)->GetName();
	    			if (detector->GetSigmaCalc() == 1)
					{
						xsecRTR = GetRTRValue(sample_energy, el_name);
	    				}
	    			else { xsecRTR = 1.;}
	    			RBS_yield = CalculateTotalRBSYield(sample_energy, A1, Mnumb[2][i], Z1, Znumb[2][i], Angle, distance, solidAngle, xsecRTR, Adens[2][i]);	    	
	    	    
	    				path = trueWorldPosition/(cos(ang2-Angle));
	    				depth = trueWorldPosition;
	    				energy_left = CalcTotEnLoss(RecEn, path, step_number, fParticle, sample_material[2]);
	    				bohr_straggling = CalculateTotalBohrStraggling(RecEn, fParticle, sample_material[2], path);	    	 
	    				
	    				if (energy_left/MeV > 0.01)
	    					{
	    						// energy loss between layers 2 and 1
	    						//path = neighbour_dist[0];
	    						//G4cout << " path " << path/um << G4endl;
	    						path = neighbour_dist[0]/(cos(ang2-Angle));
	    						depth += neighbour_dist[0];
	    						energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    						bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);	 	
	    						if (energy_left/MeV > 0.01)
	    						{
	    							// energy loss in layer1
	    							path = layer_thickness[0]/(cos(ang2-Angle));
	    							//G4cout << " path " << path/um << G4endl;
	    							depth += layer_thickness[0];
	    							energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[1]);
	    							bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[1], path);	
	    							
	    							if (energy_left/MeV > 0.01 && GetFreeSurf() == 1)
	    							{
	    								path = surf_distance/(cos(ang2-Angle));
	    								depth += surf_distance;
	    								energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    								bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    							}
	    							else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    							{
	    								run->MaxRBSDepth(trueWorldPosition);
	    								run->AddCount();
	    							}   
	    							
	    						}
	    						else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    						{
	    							energy_left = 0.;
	    							run->MaxRBSDepth(trueWorldPosition);
	    							run->AddCount();
	    						}   
	    					}
	    				else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    				{
	    					energy_left = 0.;
	    					run->MaxRBSDepth(trueWorldPosition);
	    					run->AddCount();
	    				}   
	    				   
	    				   
	    	    		else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    			{
	    				energy_left = 0.;
	    				run->MaxRBSDepth(trueWorldPosition);
	    				run->AddCount();
	    			}
	    	    
				// filling of RBS histograms
				if (energy_left/MeV > 0.01)
					{
					
					
					// include detector dead layer energy loss
					energy_left = CalculateDeadLayerEffect(energy_left, dead_material,dead_thickness,fParticle);
					
					G4double rbs_depth = depth;
					// rbs dist 
					analysisManager->FillH1(21, rbs_depth, RBS_yield*steps); // total
					analysisManager->FillH1(35, rbs_depth, RBS_yield*steps); // layer1
					// el1
					if (i == 0) 
						{
						analysisManager->FillH1(38, rbs_depth, RBS_yield*steps); //el1
						}
					// el2
					if (i == 1)
						{
						analysisManager->FillH1(39, rbs_depth, RBS_yield*steps); // el2
						}
					
					G4double sigma_det = (det_FWHM/MeV)/2.355;
					G4double sigma_det_sq = std::pow(sigma_det,2);
					G4double tot_sigma = sigma_det_sq + bohr_straggling;
					
					
					G4int GK = gauss_counter;
					G4double Gauss_sum	= 0;
					for (int k=-GK; k<GK; k++)	
						{
							G4double santykis	= k;
							G4double Gauss_value 	= (GenerateGaussian(energy_left*(1-(santykis/1000)),energy_left,tot_sigma))/GenerateGaussian(energy_left,energy_left,tot_sigma);
							Gauss_sum += Gauss_value;
						}
	
					for (int k=-GK; k<GK; k++)	
						{
							G4double santykis	= k;
							G4double Gauss_value 	= (GenerateGaussian(energy_left*(1-(santykis/1000)),energy_left,tot_sigma))/GenerateGaussian(energy_left,energy_left,tot_sigma);
							//G4double mod_G_value 	= Gauss_value/Gauss_sum;
							G4double mod_G_value 	= (Gauss_value/Gauss_sum)/(2*GK+1);
							G4double Gauss_en 	= energy_left*(1-(santykis/1000)); 
							if (Gauss_en/MeV > 0.01)
								{
								analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_G_value*steps));	// total
								analysisManager->FillH1(34, Gauss_en, (RBS_yield*mod_G_value*steps));	//layer1
								if (i == 0) 
									{ 
										analysisManager->FillH1(36,Gauss_en, (RBS_yield*mod_G_value*steps));	//el1
									}
								if (i == 1) 
									{ 
										analysisManager->FillH1(37,Gauss_en, (RBS_yield*mod_G_value*steps)); //el2
									}
								}
						}					
					}	    	    
    			} //end of elements
    		}		run->add_total_step(tot_step);
    	}// end of layer2 sensitive detector

//======================================================
// LAYER 3 	    			
//======================================================

	// layer3 sensitive detector
    if(sd3){
        int n_hit_sd = sd3->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++){
            	CrystalDetectorHit* aHit = (*sd3)[i1];
            	steps = aHit->GetStep();	
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();
	    	
	    	
	    	
	    	
	    	
	    	
	    	    
	    	G4double layer_to_surf = std::abs(-(detector->GetLength(0)/2)-layer_pos_z[2]+(detector->GetLength(3)/2));	    	    

	    	trueWorldPosition = (position.z()+detector->GetLength(0)/2)-layer_to_surf;
	    	//G4cout << " true world position " << trueWorldPosition/um << G4endl;

	    	for(int i=0;i<NoOfElements[3];i++)
			{
	    			G4double RecEn = RecoilEnergy(sample_energy, Angle, A1, Mnumb[3][i]);
	    			G4String el_name = sample_material[3]->GetElement(i)->GetName();
	    			if (detector->GetSigmaCalc() == 1)
					{
						xsecRTR = GetRTRValue(sample_energy, el_name);
	    				}
	    			else { xsecRTR = 1.;}
	    			RBS_yield = CalculateTotalRBSYield(sample_energy, A1, Mnumb[3][i], Z1, Znumb[3][i], Angle, distance, solidAngle, xsecRTR, Adens[3][i]);	    	
	    	    
	    	    			// energy loss in layer3
	    				path = trueWorldPosition/(cos(ang2-Angle));
	    				depth = trueWorldPosition;
	    				energy_left = CalcTotEnLoss(RecEn, path, step_number, fParticle, sample_material[3]);
	    				bohr_straggling = CalculateTotalBohrStraggling(RecEn, fParticle, sample_material[3], path);	    	 
	    				
	    				// calculate energy loss between layers 3 and 2 
	    				if (energy_left/MeV > 0.01)
	    					{
	    						path = neighbour_dist[1]/(cos(ang2-Angle));
	    						depth += neighbour_dist[1];
	    						energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    						bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);	 	
	    						if (energy_left/MeV > 0.01)
	    						{
	    							//calculate energy loss in layer 2
	    							path = layer_thickness[1]/(cos(ang2-Angle));
	    							//G4cout << " path " << path/um << G4endl;
	    							depth += layer_thickness[1];
	    							energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[2]);
	    							bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[2], path);	
	    							
	    							if (energy_left/MeV > 0.01)
	    							{
	    								// calculate energy loss between layers 2 and 1
	    								path = neighbour_dist[0]/(cos(ang2-Angle));
	    								depth += neighbour_dist[0];
	    								energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    								bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);	 	
	    								if (energy_left/MeV > 0.01)
	    								{
	    									//calculate energy loss in layer1 
	    										    						
	    									path = layer_thickness[0]/(cos(ang2-Angle));
	    									//G4cout << " path " << path/um << G4endl;
	    									depth += layer_thickness[0];
	    									energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[1]);
	    									bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[1], path);	
	    							
	    									if (energy_left/MeV > 0.01)
	    									{	
	    										// check whether energy loss in surface calculation is needed
	    										if(GetFreeSurf() == 1)
	    										{
	    											path = surf_distance/(cos(ang2-Angle));
	    											depth += surf_distance;
	    											energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    											bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    										}
	    									
	    									}
	    									else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    									{
	    										energy_left = 0.;
	    										run->MaxRBSDepth(trueWorldPosition);
	    										run->AddCount();
	    									}  
	    									
	    								}
	    								else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    								{
	    									energy_left = 0.;
	    									run->MaxRBSDepth(trueWorldPosition);
	    									run->AddCount();
	    								}  
	    							}
	    							else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    							{
	    								energy_left = 0.;
	    								run->MaxRBSDepth(trueWorldPosition);
	    								run->AddCount();
	    							}     	    	
	    						}
	    						else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    						{
	    							energy_left = 0.;
	    							run->MaxRBSDepth(trueWorldPosition);
	    							run->AddCount();
	    						}     	    
	    	   				}
	    	   				
	    	   			else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    					{
	    						energy_left = 0.;
	    						run->MaxRBSDepth(trueWorldPosition);
	    						run->AddCount();
	    					}
	    					
	    					
	    									// filling of RBS histograms
					if (energy_left/MeV > 0.01)
					{
					
					
						// include detector dead layer energy loss
						energy_left = CalculateDeadLayerEffect(energy_left, dead_material,dead_thickness,fParticle);
						
						G4double rbs_depth = depth;
						// rbs dist 
						analysisManager->FillH1(21, rbs_depth, RBS_yield*steps); // total
						analysisManager->FillH1(41, rbs_depth, RBS_yield*steps); // layer3
						// el1
						if (i == 0) 
							{
							analysisManager->FillH1(44, rbs_depth, RBS_yield*steps); //el1
							}
						// el2
						if (i == 1)
							{
							analysisManager->FillH1(45, rbs_depth, RBS_yield*steps); // el2
							}
					
						G4double sigma_det = (det_FWHM/MeV)/2.355;
						G4double sigma_det_sq = std::pow(sigma_det,2);
						G4double tot_sigma = sigma_det_sq + bohr_straggling;
					
					
						G4int GK = gauss_counter;
						G4double Gauss_sum	= 0;
						for (int k=-GK; k<GK; k++)	
							{
								G4double santykis	= k;
								G4double Gauss_value 	= (GenerateGaussian(energy_left*(1-(santykis/1000)),energy_left,tot_sigma))/GenerateGaussian(energy_left,energy_left,tot_sigma);
								Gauss_sum += Gauss_value;
							}
			
						for (int k=-GK; k<GK; k++)	
							{
								G4double santykis	= k;
								G4double Gauss_value 	= (GenerateGaussian(energy_left*(1-(santykis/1000)),energy_left,tot_sigma))/GenerateGaussian(energy_left,energy_left,tot_sigma);
								//G4double mod_G_value 	= Gauss_value/Gauss_sum;
								G4double mod_G_value 	= (Gauss_value/Gauss_sum)/(2*GK+1);
								G4double Gauss_en 	= energy_left*(1-(santykis/1000)); 
								if (Gauss_en/MeV > 0.01)
									{
									analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_G_value*steps));	// total
									analysisManager->FillH1(40, Gauss_en, (RBS_yield*mod_G_value*steps));	//layer1
									if (i == 0) 
										{ 
											analysisManager->FillH1(42,Gauss_en, (RBS_yield*mod_G_value*steps));	//el1
										}
									if (i == 1) 
										{ 
											analysisManager->FillH1(43,Gauss_en, (RBS_yield*mod_G_value*steps)); //el2
										}
									}
							}					
						}
	    								   
	    	   	}// end of elements
		}		run->add_total_step(tot_step);
	}


//======================================================
// LAYER 4    			
//======================================================

	// layer4 sensitive detector
    if(sd4){
        int n_hit_sd = sd4->entries();
        run->add_entry_sd(n_hit_sd);
        for(int i1=0;i1<n_hit_sd;i1++){
            	CrystalDetectorHit* aHit = (*sd4)[i1];
            	steps = aHit->GetStep();	
	    	position = aHit->GetWorldPos();
	    	sample_energy = aHit->GetKinECR();
	    	
	    	
	    	
	    	
	    	
	    	
	    	
	    	    
	    	G4double layer_to_surf = std::abs(-(detector->GetLength(0)/2)-layer_pos_z[3]+(detector->GetLength(4)/2));	    	    

	    	trueWorldPosition = (position.z()+detector->GetLength(0)/2)-layer_to_surf;
	    	//G4cout << " true world position " << trueWorldPosition/um << G4endl;

	    	for(int i=0;i<NoOfElements[4];i++)
			{
	    			G4double RecEn = RecoilEnergy(sample_energy, Angle, A1, Mnumb[4][i]);
	    			G4String el_name = sample_material[4]->GetElement(i)->GetName();
	    			if (detector->GetSigmaCalc() == 1)
					{
						xsecRTR = GetRTRValue(sample_energy, el_name);
	    				}
	    			else { xsecRTR = 1.;}
	    			RBS_yield = CalculateTotalRBSYield(sample_energy, A1, Mnumb[4][i], Z1, Znumb[4][i], Angle, distance, solidAngle, xsecRTR, Adens[4][i]);	    	
	    	    		
	    	    			// energy loss in layer4
	    				path = trueWorldPosition/(cos(ang2-Angle));
	    				depth = trueWorldPosition;
	    				energy_left = CalcTotEnLoss(RecEn, path, step_number, fParticle, sample_material[4]);
	    				bohr_straggling = CalculateTotalBohrStraggling(RecEn, fParticle, sample_material[4], path);	    	 
	    			
	    				// calculate energy loss between layers 4 and 3 
	    				if (energy_left/MeV > 0.01)
	    					{
	    					
	    						path = neighbour_dist[2]/(cos(ang2-Angle));
	    						depth += neighbour_dist[2];
	    						energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    						bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);	 	
	    						
	    						if (energy_left/MeV > 0.01)
	    							{
	    								//calculate energy loss in layer 3 
	    							path = layer_thickness[2]/(cos(ang2-Angle));
	    							//G4cout << " path " << path/um << G4endl;
	    							depth += layer_thickness[2];
	    							energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[3]);
	    							bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[3], path);	
	    						
	    						
	    							if (energy_left/MeV > 0.01)
	    								{	    						
	    									// calculate energy loss between layers 3 and 2
	    									path = neighbour_dist[1]/(cos(ang2-Angle));
	    									depth += neighbour_dist[1];
	    									energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    									bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);	 	
	    						
	    									if (energy_left/MeV > 0.01)
	    										{
	    											// calculate energy loss in layer 2
	    											path = layer_thickness[1]/(cos(ang2-Angle));
	    											//G4cout << " path " << path/um << G4endl;
	    											depth += layer_thickness[1];
	    											energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[2]);
	    											bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[2], path);	
	    											
	    											if (energy_left/MeV > 0.01)
	    											{
	    												// calculate energy between layers 2 and 1
	    												path = neighbour_dist[0]/(cos(ang2-Angle));
	    												depth += neighbour_dist[0];
	    												energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    												bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    												
	    												if (energy_left/MeV > 0.01)
	    													{
	    														//calculate energy loss in layer1
	    															    									
	    														path = layer_thickness[0]/(cos(ang2-Angle));
	    														//G4cout << " path " << path/um << G4endl;
	    														depth += layer_thickness[0];
	    														energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[1]);
	    														bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[1], path);	
	    													    																if (energy_left/MeV > 0.01)
	    															{
 																		// check whether energy loss in surface calculation is needed
	    																if(GetFreeSurf() == 1)
	    																	{
	    																	path = surf_distance/(cos(ang2-Angle));
	    																	depth += surf_distance;
	    																	energy_left = CalcTotEnLoss(energy_left, path, step_number, fParticle, sample_material[0]);
	    																	bohr_straggling += CalculateTotalBohrStraggling(energy_left, fParticle, sample_material[0], path);
	    																	}
	    															
	    															}
	    															else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    															{
	    																energy_left = 0.;
	    																run->MaxRBSDepth(trueWorldPosition);
	    																run->AddCount();
	    															}  				
	    													}
	    													else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    													{
	    														energy_left = 0.;
	    														run->MaxRBSDepth(trueWorldPosition);
	    														run->AddCount();
	    													}  
	    											}
	    											else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    											{
	    												energy_left = 0.;
	    												run->MaxRBSDepth(trueWorldPosition);
	    												run->AddCount();
	    											}  

	    										}
	    										else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    										{
	    											energy_left = 0.;
	    											run->MaxRBSDepth(trueWorldPosition);
	    											run->AddCount();
	    										}  
	    								}
	    								else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    									{
	    										energy_left = 0.;
	    										run->MaxRBSDepth(trueWorldPosition);
	    										run->AddCount();
	    									}  
	    							}
	    							else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    								{
	    									energy_left = 0.;
	    									run->MaxRBSDepth(trueWorldPosition);
	    									run->AddCount();
	    								}  
	    	    				}
	    	    				else if (energy_left/MeV > 0. && energy_left/MeV < 0.01)
	    						{
	    							energy_left = 0.;
	    							run->MaxRBSDepth(trueWorldPosition);
	    							run->AddCount();
	    						}  
	    	    		
	    					
	    				// filling of RBS histograms
					if (energy_left/MeV > 0.01)
					{
					
					
						// include detector dead layer energy loss
						energy_left = CalculateDeadLayerEffect(energy_left, dead_material,dead_thickness,fParticle);
					
						G4double rbs_depth = depth;
						// rbs dist 
						analysisManager->FillH1(21, rbs_depth, RBS_yield*steps); // total
						analysisManager->FillH1(47, rbs_depth, RBS_yield*steps); // layer4
						// el1
						if (i == 0) 
							{
							analysisManager->FillH1(50, rbs_depth, RBS_yield*steps); //el1
							}
						// el2
						if (i == 1)
							{
							analysisManager->FillH1(51, rbs_depth, RBS_yield*steps); // el2
							}
					
						G4double sigma_det = (det_FWHM/MeV)/2.355;
						G4double sigma_det_sq = std::pow(sigma_det,2);
						G4double tot_sigma = sigma_det_sq + bohr_straggling;
					
					
						G4int GK = gauss_counter;
						G4double Gauss_sum	= 0;
						for (int k=-GK; k<GK; k++)	
							{
								G4double santykis	= k;
								G4double Gauss_value 	= (GenerateGaussian(energy_left*(1-(santykis/1000)),energy_left,tot_sigma))/GenerateGaussian(energy_left,energy_left,tot_sigma);
								Gauss_sum += Gauss_value;
							}
			
						for (int k=-GK; k<GK; k++)	
							{
								G4double santykis	= k;
								G4double Gauss_value 	= (GenerateGaussian(energy_left*(1-(santykis/1000)),energy_left,tot_sigma))/GenerateGaussian(energy_left,energy_left,tot_sigma);
								//G4double mod_G_value 	= Gauss_value/Gauss_sum;
								G4double mod_G_value 	= (Gauss_value/Gauss_sum)/(2*GK+1);
								G4double Gauss_en 	= energy_left*(1-(santykis/1000)); 
								if (Gauss_en/MeV > 0.01)
									{
									analysisManager->FillH1(20, Gauss_en, (RBS_yield*mod_G_value*steps));	// total
									analysisManager->FillH1(46, Gauss_en, (RBS_yield*mod_G_value*steps));	//layer1
									if (i == 0) 
										{ 
											analysisManager->FillH1(48,Gauss_en, (RBS_yield*mod_G_value*steps));	//el1
										}
									if (i == 1) 
										{ 
											analysisManager->FillH1(49,Gauss_en, (RBS_yield*mod_G_value*steps)); //el2
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

    
    G4ThreeVector ssd[4];
    ssd[0]= G4ThreeVector(0.,0.,0.);
    ssd[1]= G4ThreeVector(0.,0.,0.);
    ssd[2]= G4ThreeVector(0.,0.,0.);
    ssd[3]= G4ThreeVector(0.,0.,0.);
    G4double energy = 0.;

    G4double kinen = 0.;

	//G4cout << " step " << step_number << G4endl;
   

    if(sdht_ID == -1) {
        //G4String sdName;
        if(SDman->FindSensitiveDetector(sdName="telescope",0)){
            sdht_ID = SDman->GetCollectionID(sdName="telescope/collection");
        }
    }
    
    SensitiveDetectorHitsCollection* sdht = 0;
    G4HCofThisEvent *hce = evt->GetHCofThisEvent();
    
    if(hce){
        if(sdht_ID != -1){
            G4VHitsCollection* aHCSD = hce->GetHC(sdht_ID);
            sdht = (SensitiveDetectorHitsCollection*)(aHCSD);
        }
    }
    int bTotalHits = 0;
    if(sdht){
        int n_hit_sd = sdht->entries();

	for(int i2=0;i2<4;i2++){
            for(int i1=0;i1<n_hit_sd;i1++)
            {
                SensitiveDetectorHit* aHit = (*sdht)[i1];
                if(aHit->GetLayerID()==i2) {
                    ssd[i2] = aHit->GetWorldPos();
                    bTotalHits++;	// checks, how many detectors got a hit
                }
                if(aHit->GetLayerID()==2) {	// info from the 3rd detector
                    energy  = aHit->GetKinE();
                }
		// particle hit in the last detector [no 4]
		if(aHit->GetLayerID()==3) {
		run->addh(); // adds hits
		kinen = aHit->GetKinE(); // takes kinetinic energy
		run->addkinen(kinen); // add kinetic energy 
		}
            }
        }
    }



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EventAction::CalcErrorFunction(G4double a)
	{
	G4double probability = std::erf(a);
	return probability;
	}


G4double EventAction::GenerateGaussian(G4double x, G4double y, G4double sigma_sq)
{
	G4double inv_sqrt_2pi = 0.398942280401433;
	G4double a = (x-y);
	
	return (inv_sqrt_2pi/std::sqrt(sigma_sq)) * std::exp(-0.5*a*a/sigma_sq);

}

G4PhysicsFreeVector* EventAction::FillRTRVector(G4String filename)
{
	//delete fPhysicsVector;
	std::ifstream vFileIn;
    	vFileIn.open(filename);
    	//G4cout << " RTR vector " << G4endl; 
    	G4String start = "RTR";
    	G4String line;
    	if(!vFileIn) 
    	{
    	G4cout <<" NO SIGMA CALC VALUES AVAILABLE " << G4endl;
        G4cout <<"Couldn't open the file"<< G4endl; //checks for input
        exit(1); 
        }
    	G4int linenum = 0;

	// 2 steps, first one checks for number of lines
	// second files up the physics vector
	
while (getline(vFileIn, line)) {
        linenum++;
        }
    	vFileIn.close();
    	
    	//G4cout << " number of lines " << linenum << G4endl;
	G4double value1;
	G4double value2;
	//G4cout << " physics vector size " << linenum-1 << G4endl;
	
	G4PhysicsFreeVector* fPhysicsVector = new G4PhysicsFreeVector(linenum-1);
	fPhysicsVector->PutValue(0,0*CLHEP::keV,1);

    	vFileIn.open(filename);
while (getline(vFileIn, line)) 
{
	// searches for line RTR
    if (line.find(start) != G4String::npos) 
    	{
    	// takes values of lines from the first one after "RTR" to the last on in the file
        for (int i=linenum-3; i; --i) 
        	{       
        		vFileIn >> value1 >> value2;
            		fPhysicsVector->PutValue(linenum-2-i,value1*CLHEP::keV,value2);
       	}     
        }                            
}
	fPhysicsVector->PutValue(linenum-2,1.01*value1*CLHEP::keV,1);

    	vFileIn.close();
	return fPhysicsVector;
	
}

	// function for Total RBS yield, combining other functions into single one
G4double EventAction::CalculateTotalRBSYield(G4double energy, G4double M1, G4double M2, G4double Z1, G4double Z2, G4double angle, G4double dist,G4double solidAngle, G4double xsecmod, G4double atomDensity)
{
	G4double CMangleX 	= CalcAngleCMFrame(angle,M1,M2);
	G4double CMenergyX 	= CalcEnergyCMFrame(energy,M1,M2);
	G4double CMxsecX 	= CalcDiffRuthXsecCM(CMenergyX,CMangleX,Z1,Z2);
	G4double CMxsecModX	= CalcRuthXsecMod(CMxsecX,Z1,Z2,CMenergyX);
	G4double CMtoLABxsecX	= CalcDiffRuthXsecLAB(M1,M2,CMangleX,CMxsecModX);
	G4double RBSyieldX 	= CalcRBSYield((CMtoLABxsecX*xsecmod), dist, solidAngle, atomDensity);
	
	return RBSyieldX;
}

G4double EventAction::GetRTRValue(G4double energy, G4String name)
{
	G4double RTR;
		if 	(name == "Si")
			{ RTR = fVectorSi->Value(energy);	}
		else if (name == "O")
			{ RTR = fVectorO->Value(energy);	}
		else if (name == "B")
			{ RTR = fVectorB->Value(energy);	} 		
		else if (name == "C")
			{ RTR = fVectorC->Value(energy);	}
		else if (name == "F")
			{ RTR = fVectorF->Value(energy);	}
		else if (name == "N")
			{ RTR = fVectorN->Value(energy);	}
		else if (name == "Na")
			{ RTR = fVectorNa->Value(energy);	}
		else if (name == "Mg")
			{ RTR = fVectorMg->Value(energy);	}
		else if (name == "P")
			{ RTR = fVectorP->Value(energy);	}
		else if (name == "S")
			{ RTR = fVectorS->Value(energy);	}
		else if (name == "K")
			{ RTR = fVectorK->Value(energy);	}
		else if (name == "Ca")
			{ RTR = fVectorCa->Value(energy);	}
		else if (name == "Ti")
			{ RTR = fVectorTi->Value(energy);	}
		else if (name == "Fe")
			{ RTR = fVectorFe->Value(energy);	}
		else if (name == "Ni")
			{ RTR = fVectorNi->Value(energy);	}
		else	{ RTR = 1.;				}

	return RTR;
}

G4double EventAction::CalculateTotalBohrStraggling(G4double energy, G4ParticleDefinition* particle, G4Material* mat, G4double distance)
{
	
	G4int elements    	= mat->GetNumberOfElements();
	G4double ionisation 	= mat->GetIonisation()->GetMeanExcitationEnergy();
	
	G4double A1_PDG = particle->GetPDGMass();
    	G4double Z1 = particle->GetAtomicNumber();
    	
    	const G4double *atomDensVector = mat->GetVecNbOfAtomsPerVolume();
    	G4double *Z2 = new G4double[elements];
    	G4double *aDensity = new G4double[elements];
    	
    	for (int i=0; i<elements; i++)
    		{
    			Z2[i] = mat->GetElement(i)->GetZ();
    			aDensity[i] = atomDensVector[i]/(1/cm3);
    		}
	G4double straggling = 0;
	
	for (int j=0; j<elements;j++)
	{
		G4double chic = CalcChi(energy, A1_PDG, Z2[j]);
		if (chic > 3)
			{ 
			straggling += CalcBohrStrag(Z1,Z2[j],aDensity[j],distance);
			}
		else 
			{
			straggling += CalcBohrStrag(Z1,Z2[j],aDensity[j],distance) * CalcStoppingNumb(energy, A1_PDG,ionisation)/2;
			}
	}
	return straggling;
}

G4double EventAction::CalculateDeadLayerEffect(G4double energy, const G4Material* material, G4double thickness,G4ParticleDefinition* particle)
{
		G4EmCalculator emCalculator;

		//G4double stop = emCalculator.GetDEDX(energy,particle,material)/(MeV/cm);
		G4double stop = emCalculator.ComputeTotalDEDX(energy,particle,material)/(MeV/cm);
		G4double en_loss = stop*(thickness/cm);

		energy -= en_loss;

	return energy;


}

G4double EventAction::CalcHeavyIonStragglingFactor(G4double energy, G4double Z1, G4double Z2)
{
	G4double scalling = 0.;
	
	G4double C1=0.01273;
	G4double C2=0.03458;
	G4double C3=0.3931;
	G4double C4=3.812;
	
	G4double kin_en = energy;
	
	G4double epsilon = kin_en/(std::pow(Z1,1.5)*sqrt(Z2));
	G4double gamma = C3*(1-exp(-C4*epsilon));
	
	//G4cout << " epsilon " << epsilon << G4endl;
	//G4cout << " gamma " << gamma << G4endl;
	
	G4double a = (std::pow(Z1,1.3333)/std::pow(Z2,0.3333))*C1*gamma;
	G4double b = std::pow((epsilon-C2),2)+std::pow(gamma,2);
	
	scalling = a/b;

	//G4cout << " scalling " << scalling << G4endl;

	return scalling;
}

