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

public:
    
    G4String GetMaterial(G4int i) 			{return fMaterialName[i];}
    void SetMaterial(G4int i,G4String aString) 	{fMaterialName[i] = aString;}

    
    G4ThreeVector GetSize(G4int i) 			{return fSizes[i];}
    void SetSize(G4int i, G4ThreeVector a3vec) 	{fSizes[i] = a3vec;}


    G4ThreeVector GetAngles() 				{return fAngles;}
    void SetAngles(G4ThreeVector a3vec) 		{fAngles = a3vec;}

     G4VPhysicalVolume* GetIntAbsorber(G4int i)	{
     						if(i == 0) return physAbsor;
     						else if(i == 1) return physAbsor2;
     						else if(i == 2) return physAbsor3;
     						else if(i == 3) return physAbsor4;
     						else if(i == 4) return physAbsor5;
     						else return physAbsor;
     						};
     
// mano
     G4VPhysicalVolume* GetAbsorber()   	{return physAbsor;};
     G4VPhysicalVolume* GetAbsorber2()   	{return physAbsor2;};
     G4VPhysicalVolume* GetAbsorber3()   	{return physAbsor3;};
     G4VPhysicalVolume* GetAbsorber4()   	{return physAbsor4;};
     G4VPhysicalVolume* GetAbsorber5()   	{return physAbsor5;};
     
     const G4VPhysicalVolume* GetWorld()   	{return worldPhysical;};
     
     G4Material* GetMaterialM(G4int i)		{return mat[i];};
     
     
     G4Material* GetMaterialComponents(G4int i, G4int j)	{return mat_components[i][j];};

     G4double GetLength(G4int i)		{return fSizes[i].z();};

     G4ThreeVector GetDimensions(G4int i )	{return fSizes[i];};

     void SetPosition(G4int i, G4ThreeVector a3vec) 	{position[i] = a3vec;};
     G4ThreeVector GetPosition(G4int i) 		{return position[i];};	

    

public:
    void SetWorldMaterial(G4String aString) {fWorldMaterial = aString;}
    G4String GetWorldMaterial() 		{return fWorldMaterial;}

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
    
    void SetCustomElement3(G4String a)	{element3 = a;}
    G4String GetCustomElement3()		{return element3;}
    
    void SetCustomMaterialDensity(G4double a)	{custom_density = a;}
    G4String GetCustomMaterialDensity()	{return custom_density;}
    
    void SetCustomElement1Part(G4double a)	{part1 = a;}
    G4String GetCustomElement1Part()		{return part1;}
    
    void SetCustomElement2Part(G4double a)	{part2 = a;}
    G4String GetCustomElement2Part()		{return part2;}
    
    void SetCustomElement3Part(G4double a)	{part3 = a;}
    G4String GetCustomElement3Part()		{return part3;}    

    
    void SetDeadLayer(G4String a)		{dead_material_name = a;};
    G4String GetDeadLayer()			{return dead_material_name;};
    void SetDeadLayerThickness(G4double a)	{dead_thickness = a;};
    G4double GetDeadLayerThickness()		{return dead_thickness;};
    G4Material* GetDeadLayerMaterial()		{return dead_material;};
    
    void SetSolidAngle(G4double a)		{solidAngle = a;};
    G4double GetSolidAngle()			{return solidAngle;};
    
    // set the number of material which should be mixed
    void SetMaterialForMix(G4int a)		{material_for_mix = a;};
    G4int GetMaterialForMix()			{return material_for_mix;};
    
    // set the material 2 for mix
    G4String GetMixingMaterial() 		{return mix_material_name;}
    void SetMixingMaterial(G4String a) 	{mix_material_name = a;}
    //
    
    // calculate detector fwhm
    G4int GetCalcFWHM()			{return enable_fwhm_calc;}
    void SetCalcFWHM(bool a)			{enable_fwhm_calc = a;}

    // Calculate MS induced energy spread    
    G4int GetMSCalc()				{return enable_MS_calc;}
    void SetMSCalc(bool a)			{enable_MS_calc = a;}        
    
    G4double GetRBSROImin()			{return rbs_roi_min;}
    void SetRBSROImin(G4double a)		{rbs_roi_min = a;}
    
    G4int GetConstAngle()			{return use_const_angle;}
    void UseConstAngle(bool a)			{use_const_angle = a;}
    
    G4int GetUseMSCorrections()		{return ms_corrections;}
    void SetUseMSCorrections(bool a)		{ms_corrections = a;}
    
	
private:

    G4double maxStep;
    G4ThreeVector position[4];

    G4double 				rbs_angle;
    G4int 				sigma_calc;
    G4int 				rbs_calc;
    G4double 				rbs_step;

    G4double 				material_mixing;    
    G4double 				sec_material_ratio;
    G4double 				detector_resolution;
    G4int 				gauss_counter;
    
    G4int				enable_custom_material;
    G4String				element1;
    G4String				element2;
    G4String				element3;
    G4double				custom_density;
    G4double 				part1;
    G4double				part2;
    G4double				part3;
    
    
    G4String				dead_material_name;

    G4double				dead_thickness;
    
    G4double				solidAngle;
    
    G4int 				material_for_mix;
    G4String				mix_material_name;
    
    G4int 				enable_fwhm_calc;
    G4int				enable_MS_calc;
    G4int				enable_NN_calc;
    G4double 				rbs_roi_min;
    G4int				use_const_angle;
    G4int				ms_corrections;
    
    G4Material* mat[5];
    G4Material*			dead_material;
    G4Material*			mixing_material;
    
    G4Material* 			mat_components[5][5];


};


