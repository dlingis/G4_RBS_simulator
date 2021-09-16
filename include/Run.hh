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
/// \file Run.hh
/// \brief Definition of the Run class
//
// $Id: Run.hh 71375 2013-06-14 07:39:33Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include "G4THitsMap.hh"
#include <map>



//test
class DetectorConstruction;
class G4ParticleDefinition;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(DetectorConstruction*);
   ~Run();

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);         
    void CountProcesses(const G4VProcess* process);
    void ParticleCount(G4String, G4double); 
    void AddEdep (G4double edep);
    void AddEflow (G4double eflow);                   
    void ParticleFlux(G4String, G4double);
    void AddEnergy (G4double);
    void AddNonIonEnergy (G4double);
    void SumTL (G4double);
    void SumTL2 (G4double energy)	{sum_TL_2 += energy; sum_TL_22 += energy*energy;}
    void SumTL3 (G4double energy)	{sum_TL_3 += energy; sum_TL_32 += energy*energy;}
    void SumT (G4double);
    void SumTt (G4double);
    void NumberRec (G4int);
    void NumberSecRec (G4int);
    void AddNumberOfSteps (G4int);
    void AddTheta (G4double);
    void AddTrakLenPrim (G4double);
    void AddTrakLenSec (G4double);
    //void AddHit() {tothit++;};
    void addh() {hit++;};	// sumuoja detektoriaus hitus
    void countEmerging() {partEmerging++;}; // sumuoja iseinusias is kristalo daleles
    void addkinen(G4double energy) { detKinEn += energy; detKinEn2 += energy*energy;}; // sumuoja kinetines energijas is 4 detektoriaus
    void addrbs() {rbs++;};
    void addsec() {second++;};
    void AddTotStep() {totstep++;};
    void AddNielStep() {nielstep++;};

    void AddProjectedRange(G4double z);
    
    void PrimaryEnergy(G4double en) 	{primary_energy += en;};

    virtual void Merge(const G4Run*);
    void EndOfRun();     


    void AddTrueRange (G4double l) { fTrueRange += l; fTrueRange2 += l*l;};
    void AddProjRange (G4double x) { fProjRange += x; fProjRange2 += x*x;};


   void AddNonIonisEnergy (G4double niel)   {totniel += niel; totniel2 += niel*niel;};

    //  absorber
   void absLEN (G4double length) {abslen += length; abslen2 = length*length;};
   void absSTP (G4int step) 	 {absstep += double(step); absstep2 += double(step)*double(step);};
   void absION (G4double edep)   {absfEdep += edep; absfEdep2 += edep*edep;};
   void absNON (G4double niel)   {absniel += niel; absniel2 += niel*niel;};
   void abssumTL(G4double energy){abssum_tl += energy; abssum_tl2 += energy*energy;};
   void abssumTL2(G4double energy){abssum_tl_2 += energy; abssum_tl_22 += energy*energy;};
   void abssumTL3(G4double energy){abssum_tl_3 += energy; abssum_tl_32 += energy*energy;};
   void abssumT(G4double energy) {abssum_t += energy; abssum_t2 += energy*energy;};
   void absnbrec(G4int rec)	 {absrec += rec;}
   //layer1 absorber
   void abs1LEN (G4double length){abs1len += length; abs1len2 = length*length;};
   void abs1STP (G4int step) 	 {abs1step += double(step); abs1step2 += double(step)*double(step);};
   void abs1ION (G4double edep)  {abs1fEdep += edep; abs1fEdep2 += edep*edep;};
   void abs1NON (G4double niel)  {abs1niel += niel; abs1niel2 += niel*niel;};
   void abs1sumTL(G4double energy){abs1sum_tl += energy; abs1sum_tl2 += energy*energy;};
   void abs1sumT(G4double energy){abs1sum_t += energy; abs1sum_t2 += energy*energy;};
   void abs1nbrec(G4int rec)	 {abs1rec += rec;}
   //layer2 absorber
   void abs2LEN (G4double length){abs2len += length; abs2len2 = length*length;};
   void abs2STP (G4int step) 	 {abs2step += double(step); abs2step2 += double(step)*double(step);};
   void abs2ION (G4double edep)  {abs2fEdep += edep; abs2fEdep2 += edep*edep;};
   void abs2NON (G4double niel)  {abs2niel += niel; abs2niel2 += niel*niel;};
   void abs2sumTL(G4double energy){abs2sum_tl += energy; abs2sum_tl2 += energy*energy;};
   void abs2sumT(G4double energy){abs2sum_t += energy; abs2sum_t2 += energy*energy;};
   void abs2nbrec(G4int rec)	 {abs2rec += rec;}   
   //layer3 absorber
   void abs3LEN (G4double length){abs3len += length; abs3len2 = length*length;};
   void abs3STP (G4int step) 	 {abs3step += double(step); abs3step2 += double(step)*double(step);};
   void abs3ION (G4double edep)  {abs3fEdep += edep; abs3fEdep2 += edep*edep;};
   void abs3NON (G4double niel)  {abs3niel += niel; abs3niel2 += niel*niel;};
   void abs3sumTL(G4double energy){abs3sum_tl += energy; abs3sum_tl2 += energy*energy;};
   void abs3sumT(G4double energy){abs3sum_t += energy; abs3sum_t2 += energy*energy;};
   void abs3nbrec(G4int rec)	 {abs3rec += rec;}      
   //layer4 absorber
   void abs4LEN (G4double length){abs4len += length; abs4len2 = length*length;};
   void abs4STP (G4int step) 	 {abs4step += double(step); abs4step2 += double(step)*double(step);};
   void abs4ION (G4double edep)  {abs4fEdep += edep; abs4fEdep2 += edep*edep;};
   void abs4NON (G4double niel)  {abs4niel += niel; abs4niel2 += niel*niel;};
   void abs4sumTL(G4double energy){abs4sum_tl += energy; abs4sum_tl2 += energy*energy;};
   void abs4sumT(G4double energy){abs4sum_t += energy; abs4sum_t2 += energy*energy;};
   void abs4nbrec(G4int rec)	 {abs4rec += rec;}         
   
   // check for entries
   void add_entry_sd(G4int en) 		{entry_sd += en;};
   void add_total_step(G4double en)		{total_step += en;};
   void MaxRBSDepth(G4double dist) 	{RBSDepth += dist; RBSDepth2 += dist*dist;};
   void AddCount()			{counts++;};		
   
   void Inc_angle(G4double a) 	{angle_of_incidence += a;};
   
   void add_entry_reach(G4double en)		{entry_reach += en;};


  
  private:
    struct ParticleData {
     ParticleData()
       : fCount(0), fEmean(0.), fEmin(0.), fEmax(0.) {}
     ParticleData(G4int count, G4double ekin, G4double emin, G4double emax)
       : fCount(count), fEmean(ekin), fEmin(emin), fEmax(emax) {}
     G4int     fCount;
     G4double  fEmean;
     G4double  fEmin;
     G4double  fEmax;
    };

  private:
    DetectorConstruction* 	fDetector;
    G4ParticleDefinition* 	fParticle;
    G4double              	fEkin;
    G4double			primary_energy;
    
    G4double fEnergyDeposit, fEnergyDeposit2;
    G4double fEnergyFlow,    fEnergyFlow2;            
    std::map<G4String,G4int>        fProcCounter;
    std::map<G4String,ParticleData> fParticleDataMap1;                    
    std::map<G4String,ParticleData> fParticleDataMap2;

    G4double        fTrueRange, fTrueRange2;             
    G4double        fProjRange, fProjRange2;

	// first material
    G4double abslen,abslen2;
    G4double absfEdep,absfEdep2;
    G4double absniel,absniel2;
    G4double absstep,absstep2;
    G4double abssum_tl,abssum_tl2;
    G4double abssum_tl_2,abssum_tl_22;
    G4double abssum_tl_3,abssum_tl_32;
    G4double abssum_t,abssum_t2;
    G4double    absrec;
	// other material
    G4double abs1len,abs1len2;
    G4double abs1fEdep,abs1fEdep2;
    G4double abs1niel,abs1niel2;
    G4double abs1step,abs1step2;
    G4double abs1sum_tl,abs1sum_tl2;
    G4double abs1sum_t,abs1sum_t2;
    G4double    abs1rec;    
	// other material
    G4double abs2len,abs2len2;
    G4double abs2fEdep,abs2fEdep2;
    G4double abs2niel,abs2niel2;
    G4double abs2step,abs2step2;
    G4double abs2sum_tl,abs2sum_tl2;
    G4double abs2sum_t,abs2sum_t2;
    G4double    abs2rec;
	// other material
    G4double abs3len,abs3len2;
    G4double abs3fEdep,abs3fEdep2;
    G4double abs3niel,abs3niel2;
    G4double abs3step,abs3step2;
    G4double abs3sum_tl,abs3sum_tl2;
    G4double abs3sum_t,abs3sum_t2;
    G4double    abs3rec;
	// other material
    G4double abs4len,abs4len2;
    G4double abs4fEdep,abs4fEdep2;
    G4double abs4niel,abs4niel2;
    G4double abs4step,abs4step2;
    G4double abs4sum_tl,abs4sum_tl2;
    G4double abs4sum_t,abs4sum_t2;
    G4double    abs4rec;    



    G4double vInitialTime;
    G4double totniel, totniel2;
    G4double RBSDepth, RBSDepth2;
    G4int counts;

  G4double EnergyDeposit,  EnergyDeposit2;
  G4double NonIonEnergyDeposit,  NonIonEnergyDeposit2;
  G4double sum_TL, sum_TL2;
  G4double sum_T, sum_T2, sum_Tt, sum_Tt2;
  G4double sum_TL_2, sum_TL_22;
  G4double sum_TL_3, sum_TL_32;

  G4double Th; 

  G4int N_rec;
  G4int N_Sec_Rec;
  G4double Nsteps, Nsteps2;
  G4double theta, theta2;

  G4double TrakLenPrim, TrakLenPrim2;
  G4double TrakLenSec, TrakLenSec2;
	

	G4int entry_sd;
	G4double entry_reach;
	G4double total_step;
	G4double rbs_angle;

	G4int hit;
	G4int partEmerging;
	G4double detKinEn, detKinEn2;
	G4int rbs;
	G4int second;
	G4double nielstep;
	G4double totstep;

    G4double projectedR, projectedR2;
    
    G4double angle_of_incidence;


 void Print(const std::vector<G4String>& title,
               const std::map< G4int, std::vector<G4double> >&out,
               G4String&) const;
    
    std::map<G4int, G4THitsMap<G4double>* > fMap;
    G4String fOutputFileSpec;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

