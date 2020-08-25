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

#include "DetectorConstruction.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4CrystalExtension.hh"
#include "G4ExtendedMaterial.hh"
#include "G4LogicalCrystalVolume.hh"

#include "G4ChannelingMaterialData.hh"
#include "G4ChannelingOptrMultiParticleChangeCrossSection.hh"

#include "CrystalDetector.hh"
#include "SensitiveDetector.hh"

#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():
stepLimit(0),
fMaterialName{"G4_Si","G4_Si","G4_Si","G4_Si","G4_Si"},
fSizes{G4ThreeVector(1.,1.,1.),G4ThreeVector(1.,1.,1.),G4ThreeVector(1.,1.,1.),G4ThreeVector(1.,1.,1.),G4ThreeVector(1.,1.,1.)},
fAngles(G4ThreeVector(0.,0.,0.)),
fWorldMaterial("G4_Galactic"),
fDetectorMaterialName(""),
fDetectorSizes(G4ThreeVector(50. * CLHEP::mm,50. * CLHEP::mm,1 * CLHEP::mm)),
fDetectorDistance{-20. * CLHEP::cm,-19. * CLHEP::cm,+19. * CLHEP::cm,+20. * CLHEP::cm,+40. * CLHEP::cm},maxStep(1.*nm),
position{G4ThreeVector(0.,0.,0.),G4ThreeVector(0.,0.,0.),G4ThreeVector(0.,0.,0.),G4ThreeVector(0.,0.,0.)},
rbs_angle(160.*degree),sigma_calc(0),rbs_calc(0),rbs_step(0),material_mixing(0.), sec_material_ratio(0.),
detector_resolution(10.*keV),gauss_counter(5),
/*custom material definition*/enable_custom_material(0),element1(""),element2(""),custom_density(0.),
part1(0),part2(0),dead_thickness(0.),solidAngle(1.),
material_for_mix(0.),mix_material_name("G4_Si") //dead_material_name(""),
{
    fMessenger = new DetectorConstructionMessenger(this); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::DefineMaterials(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){
    
	

    //** World **//
    G4Material* worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial(fWorldMaterial);
    
    G4double worldSizeXY = 2. * CLHEP::meter;
    G4double worldSizeZ = 30. * CLHEP::meter; // buvo 30 metru
    
    G4Box* worldSolid = new G4Box("world.solid",
                                  worldSizeXY/2.,
                                  worldSizeXY/2.,
                                  worldSizeZ/2.);
    
    G4LogicalVolume* worldLogic = new G4LogicalVolume(worldSolid,
                                                      worldMaterial,
                                                      "world.logic");
    
    /*G4PVPlacement* */worldPhysical = new G4PVPlacement(0,
                                                     G4ThreeVector(),
                                                     worldLogic,
                                                     "world.physic",
                                                     0,
                                                     false,
                                                     0);
    
    
    	if (position[0]>position[1] || position[1] > position[2] || position[2] > position[3])
		{
			G4cout << " ************************ " << G4endl;
			G4cout << " ERROR in layer placement " << G4endl;
			G4cout << " ************************ " << G4endl;
			exit(1);
		}
    
    
    //** Detectors instantiation **//
    G4Box* ssdSolid = new G4Box("ssd.solid",
                                fDetectorSizes.x()/2.,
                                fDetectorSizes.y()/2.,
                                fDetectorSizes.z()/2.);
    
    G4Material* detectorMaterial;
    if(fDetectorMaterialName == ""){
        detectorMaterial = worldMaterial;
    }
    else{
        detectorMaterial = G4NistManager::Instance()->FindOrBuildMaterial(fDetectorMaterialName);
    }
    
    G4LogicalVolume* ssdLogic = new G4LogicalVolume(ssdSolid,
                                                    detectorMaterial,
                                                    "ssd.logic");
    
    for(size_t i1=0;i1<4;i1++){
        new G4PVPlacement(0,
                          G4ThreeVector(0.,0.,fDetectorDistance[i1]),
                          ssdLogic,
                          "ssd.physic",
                          worldLogic,
                          false,
                          i1);
        
    }

//    for(int i = 0;i<5;i++)
//    {
//    	G4String solidName = "crystal.solid" + std::to_string(i+1);

   	G4Box* crystalSolid = new G4Box("crystal.Solid",
                                    fSizes[0].x()/2.,
                                    fSizes[0].y()/2.,
                                    fSizes[0].z()/2.);
                                         
    // another crystal by me
    G4Box* crystalSolid2 = new G4Box("crystal.solid2",
                                    fSizes[1].x()/2.,
                                    fSizes[1].y()/2.,
                                    fSizes[1].z()/2.);

    G4Box* crystalSolid3 = new G4Box("crystal.solid3",
                                    fSizes[2].x()/2.,
                                    fSizes[2].y()/2.,
                                    fSizes[2].z()/2.);

    // another crystal by me
    G4Box* crystalSolid4 = new G4Box("crystal.solid4",
                                    fSizes[3].x()/2.,
                                    fSizes[3].y()/2.,
                                    fSizes[3].z()/2.);
                                    
    // another crystal by me
    G4Box* crystalSolid5 = new G4Box("crystal.solid5",
                                    fSizes[4].x()/2.,
                                    fSizes[4].y()/2.,
                                    fSizes[4].z()/2.);




	// niobium pentoxide
  	//G4double z;                            
  	G4String name, symbol;
  	//G4double density1 = 92.9*g/mole;
  	//G4Element* elNb = new G4Element(name="Niobium", symbol="Nb", z=41., density1);
  	//G4double density2 = 16.00*g/mole;
  	//G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8., density2);
    
	G4Element* elHfn = G4NistManager::Instance()->FindOrBuildElement(72);
    	G4Element* elOx = G4NistManager::Instance()->FindOrBuildElement(8);
    	G4Element* elNb = G4NistManager::Instance()->FindOrBuildElement(41);
    	G4Element* elIn = G4NistManager::Instance()->FindOrBuildElement(49);
    	G4Element* elP = G4NistManager::Instance()->FindOrBuildElement(15);
    
    
  	G4int ncomponents, natoms;

  	G4double density3 = 4.60*g/cm3;
  	G4Material* Nb2O5 = new G4Material(name="Nb2O5", density3, ncomponents=2);
  	Nb2O5->AddElement(elNb, natoms=2);
  	Nb2O5->AddElement(elOx, natoms=5);
  	
  	
  	//G4cout << " nb2o5 density " << Nb2O5->GetDensity()/(g/cm3) << G4endl;


	G4Element* elC = new G4Element(name="Carbon", symbol="C", 6., 12.*g/mole);
	G4Element* elN = new G4Element(name="Nitrogen", symbol="N", 7., 14.*g/mole);

	G4Material* CN = new G4Material(name="CN", 2.25*g/cm3, ncomponents=2);
	CN->AddElement(elC, natoms=9);
	CN->AddElement(elN, natoms=1);


    	G4Material* InP = new G4Material("InP",4.81 *CLHEP::g/CLHEP::cm3,2);
    	InP->AddElement(elIn,1);
    	InP->AddElement(elP,1);
    	
    	G4Material* HfO2 = new G4Material("HfO2",9.68 *CLHEP::g/CLHEP::cm3,2);
    	HfO2->AddElement(elHfn,1);
    	HfO2->AddElement(elOx,2);

	

	G4String plus = " + ";
	G4String gap  = " ";
	
	

    //** Crystal Definition Start **//
    G4cout << "Material Name: " << fMaterialName[0] << G4endl;
    if(fMaterialName[0] == ""){
        mat[0] = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    }
    else if (fMaterialName[0] == "Nb2O5"){
	mat[0] = Nb2O5;
	}

    else if (fMaterialName[0] == "CN"){
	mat[0] = CN;
	}

    else if (fMaterialName[0] == "InP"){
	mat[0] = InP;
	}
    else if (fMaterialName[0] == "HfO2"){
	mat[0] = HfO2;
	}
    else{
        mat[0] = G4NistManager::Instance()->FindOrBuildMaterial(fMaterialName[0]);
    }
    
    //** Crystal Definition Start **//
    G4cout << "Material Name2: " << fMaterialName[1] << G4endl;
    if(fMaterialName[1] == ""){
        mat[1] = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    }
    else if (fMaterialName[1] == "Nb2O5"){
	mat[1] = Nb2O5;
	}

    else if (fMaterialName[1] == "CN"){
	mat[1] = CN;
	}

    else if (fMaterialName[1] == "InP"){
	mat[1] = InP;
	}
    else if (fMaterialName[1] == "HfO2"){
	mat[1] = HfO2;
	}
    else{
        mat[1] = G4NistManager::Instance()->FindOrBuildMaterial(fMaterialName[1]);
    }
    
    //** Crystal Definition Start **//
    G4cout << "Material Name3: " << fMaterialName[2] << G4endl;
    if(fMaterialName[2] == ""){
        mat[2] = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    }
    else if (fMaterialName[2] == "Nb2O5"){
	mat[2] = Nb2O5;
	}

    else if (fMaterialName[2] == "CN"){
	mat[2] = CN;
	}

    else if (fMaterialName[2] == "InP"){
	mat[2] = InP;
	}
    else if (fMaterialName[2] == "HfO2"){
	mat[2] = HfO2;
	}
    else{
        mat[2] = G4NistManager::Instance()->FindOrBuildMaterial(fMaterialName[2]);
    }
    
    //** Crystal Definition Start **//
    G4cout << "Material Name4: " << fMaterialName[3] << G4endl;
    if(fMaterialName[3] == ""){
        mat[3] = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    }
    else if (fMaterialName[3] == "Nb2O5"){
	mat[3] = Nb2O5;
	}

    else if (fMaterialName[3] == "CN"){
	mat[3] = CN;
	}

    else if (fMaterialName[3] == "InP"){
	mat[3] = InP;
	}
    else if (fMaterialName[3] == "HfO2"){
	mat[3] = HfO2;
	}
    else{
        mat[3] = G4NistManager::Instance()->FindOrBuildMaterial(fMaterialName[3]);
    }
    
    //** Crystal Definition Start **//
    G4cout << "Material Name5: " << fMaterialName[4] << G4endl;
    if(fMaterialName[4] == ""){
        mat[4] = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    }
    else if (fMaterialName[4] == "Nb2O5"){
	mat[4] = Nb2O5;
	}

    else if (fMaterialName[4] == "CN"){
	mat[4] = CN;
	}

    else if (fMaterialName[4] == "InP"){
	mat[4] = InP;
	}
    else if (fMaterialName[4] == "HfO2"){
	mat[4] = HfO2;
	}
    else{
        mat[4] = G4NistManager::Instance()->FindOrBuildMaterial(fMaterialName[4]);
    }
    

    
    
    if(dead_material_name == ""){
        dead_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    }
    else {    dead_material = G4NistManager::Instance()->FindOrBuildMaterial(dead_material_name); }
    
    
    G4cout << "Detector dead layer: " << dead_thickness/nm << " nm of " << dead_material_name << " Z " << dead_material->GetZ() << G4endl;   
	// end
	
	
	
    G4double a = 1.*nm;
    G4Box* phantom_layer = new G4Box("phantom.layer",
    					a,a,dead_thickness);
    	
    	G4LogicalVolume* phantom_log = new G4LogicalVolume (phantom_layer, dead_material, "phantom_logic");				
    					
    	G4double b = 20.*cm;
    G4VPhysicalVolume* phantom   = new G4PVPlacement(0,
                      G4ThreeVector(b,b,b),
                      phantom_log,
                      "phantom_physic",
                      worldLogic,
                      false,
                      0,
		      true);
    
	
	
	
	
	if (material_mixing != 0)
		{	
		G4int no_of_material_mix = GetMaterialForMix();	
		G4cout << " mix no " << no_of_material_mix << G4endl;
		mixing_material = G4NistManager::Instance()->FindOrBuildMaterial(mix_material_name);
		G4double den1 = mixing_material->GetDensity()/(g/cm3);
		G4double den2 = mat[no_of_material_mix]->GetDensity()/(g/cm3);
		G4String ratio = std::to_string(sec_material_ratio);
		G4String ratio2 = std::to_string(1-sec_material_ratio);
		//G4String plus = " + ";
		//G4String gap  = " ";
		G4String mix_name = ratio + gap + fMaterialName[no_of_material_mix] + plus + ratio2 + gap + mix_material_name;
		G4Material* material_mix = new G4Material(name = mix_name,(den1*(1-sec_material_ratio)+den2*sec_material_ratio)*g/cm3,ncomponents = 2);
		
		material_mix->AddMaterial(mixing_material, 1-sec_material_ratio);
		material_mix->AddMaterial(mat[no_of_material_mix], sec_material_ratio);
		
		//if (
		
		mat[no_of_material_mix] = material_mix;
		//mat2=material_mix;
		}

	if (enable_custom_material == 1){
	//G4cout << " enabled " << G4endl;
	//part2 = 1-part1;
	G4String ratio1 = std::to_string(part1);
	G4String ratio2 = std::to_string(part2);
	G4String custom_material_name = ratio1 + element1 + plus + ratio2 + element2;
	// addition of custom material
	G4Element* el1 = G4NistManager::Instance()->FindOrBuildElement(element1);
	G4cout << " el 1: name " << el1->GetName() << " Z " << el1->GetZ() << " symb " << el1->GetSymbol() << G4endl;
	G4Element* el2 = G4NistManager::Instance()->FindOrBuildElement(element2);
	G4cout << " el 2: name " << el2->GetName() << " Z " << el2->GetZ() << " symb " << el2->GetSymbol() << G4endl;
	G4Material* custom_material = new G4Material(custom_material_name,custom_density,2);
	custom_material->AddElement(el1,part1);
	custom_material->AddElement(el2,part2);
	
	
	//mat2=custom_material;
	}


   //step limiter
  //maxStep = 1.*nm;
// **************************************************************
// bandymas primesti amorfini suda




	// test rotation
    G4RotationMatrix* rot = new G4RotationMatrix;
    if(fAngles.x()!=0.){
        rot->rotateX(fAngles.x());
    }
    if(fAngles.y()!=0.){
        rot->rotateY(fAngles.y());
    }
    if(fAngles.z()!=0.){
        rot->rotateZ(fAngles.z());
    }
	//


	G4LogicalVolume* crystalLogic_amo = new G4LogicalVolume (crystalSolid, mat[0], "crystal.logic_amo");

	
	// check is main crystal is amorphous

	//G4cout << "\n test fCrystalAmorphous == 1 \n" << G4endl;
   	physAbsor = new G4PVPlacement(0, //
                      G4ThreeVector(),
                      crystalLogic_amo,
                      "crystal.physic_amo",
                      worldLogic,
                      false,
                      0,
		      true);

  stepLimit = new G4UserLimits(maxStep);
  crystalLogic_amo->SetUserLimits(stepLimit);
    // end of test


/*
    for(int i = 1;i<5;i++)
    {
    	G4String logicName = "crystal.logic" + std::to_string(i+1);

    	crystalLogic_amo_mid[i] = new G4LogicalVolume(crystalSolid[i],mat[i],logicName);
	crystalLogic_amo_mid[i]->SetUserLimits(stepLimit);
    	
    }*/
    

	G4LogicalVolume* crystalLogic2_amo = new G4LogicalVolume (crystalSolid2, mat[1], "crystal.logic2_amo");
	G4LogicalVolume* crystalLogic3_amo = new G4LogicalVolume (crystalSolid3, mat[2], "crystal.logic3_amo");
	G4LogicalVolume* crystalLogic4_amo = new G4LogicalVolume (crystalSolid4, mat[3], "crystal.logic4_amo");	
	G4LogicalVolume* crystalLogic5_amo = new G4LogicalVolume (crystalSolid5, mat[4], "crystal.logic5_amo");




  		stepLimit = new G4UserLimits(maxStep);
		crystalLogic2_amo->SetUserLimits(stepLimit);
		crystalLogic3_amo->SetUserLimits(stepLimit);
		crystalLogic4_amo->SetUserLimits(stepLimit);
		crystalLogic5_amo->SetUserLimits(stepLimit);

// if second material is outside first one, the check on whether the first one is amorphous or not is not made
//if (fCrystalInside == 0){
/*	
		   	physAbsor2 = new G4PVPlacement(0,
                	position[0],
                	crystalLogic2_amo,
                	"crystal.physic2_amo",
                	worldLogic,	
                	false,
                	0,
			true);
*/
//}
// if second material is inside the first one, the check for first material must be made 
//else if(fCrystalInside == 1)
//{

		   	physAbsor2 = new G4PVPlacement(0,
                      	position[0],
                      	crystalLogic2_amo,
                      	"crystal.physic2_amo",
                     	crystalLogic_amo,	
                      	false,
                      	0,
		      	true);
//}

		   	physAbsor3 = new G4PVPlacement(0,
                      	position[1],
                      	crystalLogic3_amo,
                      	"crystal.physic3_amo",
                     	crystalLogic_amo,	
                      	false,
                      	0,
		      	true);
		      	
		      	
		   	physAbsor4 = new G4PVPlacement(0,
                      	position[2],
                      	crystalLogic4_amo,
                      	"crystal.physic4_amo",
                     	crystalLogic_amo,	
                      	false,
                      	0,
		      	true);

		   	physAbsor5 = new G4PVPlacement(0,
                      	position[3],
                      	crystalLogic5_amo,
                      	"crystal.physic5_amo",
                     	crystalLogic_amo,	
                      	false,
                      	0,
		      	true);



    return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField(){

        G4LogicalVolume* crystalLogic_amo = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal.logic_amo");
        G4VSensitiveDetector* crystaldetector = new CrystalDetector("/crystaldetector");
        G4SDManager::GetSDMpointer()->AddNewDetector(crystaldetector);
        crystalLogic_amo->SetSensitiveDetector(crystaldetector);
        
        G4LogicalVolume* crystalLogic2_amo = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal.logic2_amo");
        G4VSensitiveDetector* crystaldetector2 = new CrystalDetector("/crystaldetector2");
        G4SDManager::GetSDMpointer()->AddNewDetector(crystaldetector2);
        crystalLogic2_amo->SetSensitiveDetector(crystaldetector2);


        G4LogicalVolume* crystalLogic3_amo = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal.logic3_amo");
        G4VSensitiveDetector* crystaldetector3 = new CrystalDetector("/crystaldetector3");
        G4SDManager::GetSDMpointer()->AddNewDetector(crystaldetector3);
        crystalLogic3_amo->SetSensitiveDetector(crystaldetector3);

        G4LogicalVolume* crystalLogic4_amo = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal.logic4_amo");
        G4VSensitiveDetector* crystaldetector4 = new CrystalDetector("/crystaldetector4");
        G4SDManager::GetSDMpointer()->AddNewDetector(crystaldetector4);
        crystalLogic4_amo->SetSensitiveDetector(crystaldetector4);

        G4LogicalVolume* crystalLogic5_amo = G4LogicalVolumeStore::GetInstance()->GetVolume("crystal.logic5_amo");
        G4VSensitiveDetector* crystaldetector5 = new CrystalDetector("/crystaldetector5");
        G4SDManager::GetSDMpointer()->AddNewDetector(crystaldetector5);
        crystalLogic5_amo->SetSensitiveDetector(crystaldetector5);
        

        
    G4LogicalVolume* ssdLogic = G4LogicalVolumeStore::GetInstance()->GetVolume("ssd.logic");
    G4VSensitiveDetector* telescope = new SensitiveDetector("/telescope");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    for(unsigned int i1=0;i1<4;i1++){
        ssdLogic->SetSensitiveDetector(telescope);
    }
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

