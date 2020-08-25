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
// --------------------------------------------------------------
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1
#endif

#include "G4VUserDetectorConstruction.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "DetectorConstructionMessenger.hh"

#include "G4LogicalCrystalVolume.hh"

#include "globals.hh"


#include "G4UserLimits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


class DetectorConstruction : public G4VUserDetectorConstruction
{



public:
    
    DetectorConstruction();
    ~DetectorConstruction();
    
    void DefineMaterials();
    G4VPhysicalVolume* Construct();

private:
    void ConstructSDandField();

private:
    DetectorConstructionMessenger* fMessenger;
    G4VPhysicalVolume* physAbsor;
    G4VPhysicalVolume* physAbsor2;
    G4VPhysicalVolume* physAbsor3;
    G4VPhysicalVolume* physAbsor4;
    G4VPhysicalVolume* physAbsor5;
    G4VPhysicalVolume* worldPhysical; 

    G4UserLimits* stepLimit;             // pointer to user step limits
    
    G4String fMaterialName[5];

    G4ThreeVector fSizes[5];
 

    G4ThreeVector fAngles;
    G4String fWorldMaterial;

    //G4LogicalVolume*				crystalLogic_amo;
    //G4LogicalVolume*				crystalLogic_amo_mid[4];

public:
    
    G4String GetMaterial(G4int i) 			{return fMaterialName[i];}
    void SetMaterial(G4int i,G4String aString) 	{fMaterialName[i] = aString;}

    
    G4ThreeVector GetSize(G4int i) 			{return fSizes[i];}
    void SetSize(G4int i, G4ThreeVector a3vec) 	{fSizes[i] = a3vec;}


    G4ThreeVector GetAngles() 			{return fAngles;}
    void SetAngles(G4ThreeVector a3vec) 		{fAngles = a3vec;}


// mano
     G4VPhysicalVolume* GetAbsorber()   	{return physAbsor;};
     G4VPhysicalVolume* GetAbsorber2()   	{return physAbsor2;};
     G4VPhysicalVolume* GetAbsorber3()   	{return physAbsor3;};
     G4VPhysicalVolume* GetAbsorber4()   	{return physAbsor4;};
     G4VPhysicalVolume* GetAbsorber5()   	{return physAbsor5;};
     
     const G4VPhysicalVolume* GetWorld()   	{return worldPhysical;};
     
     G4Material* GetMaterialM(G4int i)	{return mat[i];};

     G4double GetLength(G4int i)		{return fSizes[i].z();};

     G4ThreeVector GetDimensions(G4int i )		{return fSizes[i];};
/*     G4ThreeVector GetDimensions2()		{return fSizes[1];};
     G4ThreeVector GetDimensions3()		{return fSizes[2];};
     G4ThreeVector GetDimensions4()		{return fSizes[3];};
     G4ThreeVector GetDimensions5()		{return fSizes[4];};
*/

     void SetPosition(G4int i, G4ThreeVector a3vec) 	{position[i] = a3vec;};
     G4ThreeVector GetPosition(G4int i) 		{return position[i];};	

    

public:
    void SetWorldMaterial(G4String aString) {fWorldMaterial = aString;}
    G4String GetWorldMaterial() 		{return fWorldMaterial;}

    void SetDetectorMaterial(G4String aString) 	{fDetectorMaterialName = aString;}
    G4String GetDetectorMaterial() 		{return fDetectorMaterialName;}
    
    void SetDetectorSizes(G4ThreeVector a3vec) 	{fDetectorSizes = a3vec;}
    G4ThreeVector GetDetectorSizes() 		{return fDetectorSizes;}

    void SetDetectorDistance(G4int aInt,G4double aDouble) {fDetectorDistance[aInt] = aDouble;}
    G4double GetDetectorDistance(G4int aInt) 	{return fDetectorDistance[aInt];}

    //pridetas amorfizavimas, testui
    //void SetInside(bool a) 			{fCrystalInside = a;};
    //G4int GetCrystalInside() 			{return fCrystalInside;}

    G4double GetMaxStep()			{return maxStep;}
    void SetMaxStep(G4double stp)		{maxStep = stp;}

    G4double GetRBSAngle()			{return rbs_angle;}
    void SetRBSAngle(G4double ang)		{rbs_angle = ang;}

    G4int GetSigmaCalc()			{return sigma_calc;}
    void SetSigmaCal(bool a)			{sigma_calc = a;}

    G4int GetRBSCalc()				{return rbs_calc;}
    void SetRBSCalc(bool a)			{rbs_calc = a;}

    void SetEnLossStep(G4double a)		{rbs_step = a;}
    G4double GetEnLossStep()			{return rbs_step;}
    
    void SetMixing(bool a)			{material_mixing = a;}
    G4double GetMaterialMixing()		{return material_mixing;}
    
    void SetMixingRatio(G4double a)		{sec_material_ratio = a;}
    G4double GetMixingRatio()			{return sec_material_ratio;}
    
    void SetDetectorResolution(G4double a)	{detector_resolution = a;}
    G4double GetDetectorResolution()		{return detector_resolution;}
    
    void SetGaussCounter(G4int a)		{gauss_counter = a;}
    G4int GetGaussCounter()			{return gauss_counter;}
    
    
    // custom material
    
    void SetCustomMaterial(bool a)		{enable_custom_material = a;}
    G4int GetCustomMaterial()			{return enable_custom_material;}
    
    void SetCustomElement1(G4String a)	{element1 = a;}
    G4String GetCustomElement1()		{return element1;}
    
    void SetCustomElement2(G4String a)	{element2 = a;}
    G4String GetCustomElement2()		{return element2;}
    
    void SetCustomMaterialDensity(G4double a)	{custom_density = a;}
    G4String GetCustomMaterialDensity()	{return custom_density;}
    
    void SetCustomElement1Part(G4double a)	{part1 = a;}
    G4String GetCustomElement1Part()		{return part1;}
    
    void SetCustomElement2Part(G4double a)	{part2 = a;}
    G4String GetCustomElement2Part()		{return part2;}
    
    void SetDeadLayer(G4String a)		{dead_material_name = a;};
    G4String GetDeadLayer()			{return dead_material_name;};
    void SetDeadLayerThickness(G4double a)	{dead_thickness = a;};
    G4double GetDeadLayerThickness()		{return dead_thickness;};
    G4Material* GetDeadLayerMaterial()	{return dead_material;};
    
    void SetSolidAngle(G4double a)		{solidAngle = a;};
    G4double GetSolidAngle()			{return solidAngle;};
    
    // set the number of material which should be mixed
    void SetMaterialForMix(G4int a)		{material_for_mix = a;};
    G4int GetMaterialForMix()			{return material_for_mix;};
    
    // set the material 2 for mix
    G4String GetMixingMaterial() 		{return mix_material_name;}
    void SetMixingMaterial(G4String a) 	{mix_material_name = a;}
    //
	
private:
    G4String fDetectorMaterialName;
    G4ThreeVector fDetectorSizes;
    G4double fDetectorDistance[5];

    //G4bool fCrystalInside;

    G4double maxStep;
    G4ThreeVector position[4];

    G4double rbs_angle;
    G4int sigma_calc;
    G4int rbs_calc;
    G4double rbs_step;

    G4double material_mixing;    
    G4double sec_material_ratio;
    G4double detector_resolution;
    G4int gauss_counter;
    
    G4int				enable_custom_material;
    G4String				element1;
    G4String				element2;
    G4double				custom_density;
    G4double 				part1;
    G4double				part2;
    
    
    G4String				dead_material_name;

    G4double				dead_thickness;
    
    G4double				solidAngle;
    
    G4int 				material_for_mix;
    G4String				mix_material_name;
    
    //mano
    //G4Material* material;
    G4Material* mat[5];
    G4Material*			dead_material;
    G4Material*			mixing_material;
    //G4Box* 	crystalSolid[5];


};


