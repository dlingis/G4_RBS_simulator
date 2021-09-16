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
/// \file Run.cc
/// \brief Implementation of the Run class
//
// $Id: Run.cc 71376 2013-06-14 07:44:50Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4BetheBlochModel.hh"
#include "G4BraggModel.hh"
#include "G4Proton.hh"

#include "G4hIonisation.hh"

#include <ctime>


#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
: G4Run(),
  fDetector(det), fParticle(0), fEkin(0.),  fTrueRange(0.), fTrueRange2(0.),
  fProjRange(0.), fProjRange2(0.), fMap()
{
  //TotalCount = 0.;
  fEnergyDeposit = fEnergyDeposit2 = 0.;
  fEnergyFlow    = fEnergyFlow2    = 0.;  
  EnergyDeposit  = EnergyDeposit2  = 0.;
  NonIonEnergyDeposit = NonIonEnergyDeposit2 = 0.;
  sum_TL = sum_TL2 = 0.;
  sum_T = sum_T2 = 0.;
  sum_Tt = sum_Tt2 = 0.;
  Nsteps = Nsteps2 = 0.;
  theta = theta2 = 0.;
  TrakLenPrim = TrakLenPrim2 = 0.;
  TrakLenSec  = TrakLenSec2  = 0.;
  N_rec = 0.;
  N_Sec_Rec = 0.;
  Th=28.*eV;  
  fTrueRange=0.;
  fTrueRange2=0.;
  fProjRange=0.;
  fProjRange2=0.;
   hit = 0.;
	partEmerging = 0.;
	detKinEn = detKinEn2 = 0.;
	rbs = 0.;
	second = 0.;
	nielstep=0.;
	totstep=0.;
  //	first material

   abslen=abslen2=0.;
   absfEdep=absfEdep2=0.;
   absniel=absniel2=0.;
   absstep=absstep2=0.;
   abssum_tl=abssum_tl2=0.;
   abssum_t=abssum_t2=0.;
   absrec=0.;

   abs1len=abs1len2=0.;
   abs1fEdep=abs1fEdep2=0.;
   abs1niel=abs1niel2=0.;
   abs1step=abs1step2=0.;
   abs1sum_tl=abs1sum_tl2=0.;
   abs1sum_t=abs1sum_t2=0.;
   abs1rec=0.;

	// other material

   abs2len=abs2len2=0.;
   abs2fEdep=abs2fEdep2=0.;
   abs2niel=abs2niel2=0.;
   abs2step=abs2step2=0.;
   abs2sum_tl=abs2sum_tl2=0.;
   abs2sum_t=abs2sum_t2=0.;
   abs2rec=0.;

   abs3len=abs3len2=0.;
   abs3fEdep=abs3fEdep2=0.;
   abs3niel=abs3niel2=0.;
   abs3step=abs3step2=0.;
   abs3sum_tl=abs3sum_tl2=0.;
   abs3sum_t=abs3sum_t2=0.;
   abs3rec=0.;

   abs4len=abs4len2=0.;
   abs4fEdep=abs4fEdep2=0.;
   abs4niel=abs4niel2=0.;
   abs4step=abs4step2=0.;
   abs4sum_tl=abs4sum_tl2=0.;
   abs4sum_t=abs4sum_t2=0.;
   abs4rec=0.;


   totniel=totniel2=0.;


   projectedR=projectedR2=0.;

	entry_sd 	= 0.;
	total_step 	= 0.;
	entry_reach	= 0.;
   RBSDepth=RBSDepth2=0.;
   counts=0;
        

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ 
    // Important to clean up the map
    std::map<G4int, G4THitsMap<G4double>* >::iterator iter = fMap.begin();
    
    while (iter != fMap.end()) {
        delete iter->second;
        iter++;}
}
//========================================================
void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(const G4VProcess* process) 
{
  G4String procName = process->GetProcessName();
  std::map<G4String,G4int>::iterator it = fProcCounter.find(procName);
  if ( it == fProcCounter.end()) {
    fProcCounter[procName] = 1;
  }
  else {
    fProcCounter[procName]++; 
  }
}

void Run::AddProjectedRange (G4double z) 
{
  projectedR  += z;
  projectedR2 += z*z;
}


void Run::AddEnergy (G4double edep)
{ 
 EnergyDeposit += edep; 
 EnergyDeposit2 += edep*edep; 
}

void Run::AddNonIonEnergy (G4double enondep)
{ 
  NonIonEnergyDeposit += enondep; 
  NonIonEnergyDeposit2 += enondep*enondep; 
}

   // sum of secondary Kinetic energy*L(T):
void Run::SumTL(G4double energyL)
{
 sum_TL +=energyL;  
 sum_TL2 +=energyL*energyL;
}

  // sum of secondary Kinetic energy:
void Run::SumT(G4double energy)
{
  sum_T +=energy;  
  sum_T2 +=energy*energy; 
}

  // sum of tertiary particles energies
void Run::SumTt(G4double energy)
{
  sum_Tt +=energy;  
  sum_Tt2 +=energy*energy; 
}


  //number of recoils
void Run::NumberRec(G4int i)
{ 
  N_rec += i; 
}

  //number of tertiary particles
void Run::NumberSecRec(G4int i)
{
  N_Sec_Rec +=i;
}

  //number of steps
void Run::AddNumberOfSteps(G4int i )
{ 
  Nsteps+=double(i); 
  Nsteps2+=double(i)*double(i);
}
  //scattering angle	
void Run::AddTheta(G4double tet)
{ 
  theta+=tet; 
  theta2+=tet*tet;
}	
  // track length	
void Run::AddTrakLenPrim (G4double length)
{
  //TotalCount++;
  //
  TrakLenPrim += length; 
  TrakLenPrim2 += length*length;
}
void Run::AddTrakLenSec (G4double length)
{
  TrakLenSec += length; 
  TrakLenSec2 += length*length;
}
  
                  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ParticleCount(G4String name, G4double Ekin)
{
  std::map<G4String, ParticleData>::iterator it = fParticleDataMap1.find(name);
  if ( it == fParticleDataMap1.end()) {
    fParticleDataMap1[name] = ParticleData(1, Ekin, Ekin, Ekin);
  }
  else {
    ParticleData& data = it->second;
    data.fCount++;
    data.fEmean += Ekin;
    //update min max
    G4double emin = data.fEmin;
    if (Ekin < emin) data.fEmin = Ekin;
    G4double emax = data.fEmax;
    if (Ekin > emax) data.fEmax = Ekin; 
  }   
}
                 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddEdep(G4double edep)
{ 
  fEnergyDeposit += edep;
  fEnergyDeposit2 += edep*edep;
}
                 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddEflow(G4double eflow)
{ 
  fEnergyFlow += eflow;
  fEnergyFlow2 += eflow*eflow;
}                  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ParticleFlux(G4String name, G4double Ekin)
{
  std::map<G4String, ParticleData>::iterator it = fParticleDataMap2.find(name);
  if ( it == fParticleDataMap2.end()) {
    fParticleDataMap2[name] = ParticleData(1, Ekin, Ekin, Ekin);
  }
  else {
    ParticleData& data = it->second;
    data.fCount++;
    data.fEmean += Ekin;
    //update min max
    G4double emin = data.fEmin;
    if (Ekin < emin) data.fEmin = Ekin;
    G4double emax = data.fEmax;
    if (Ekin > emax) data.fEmax = Ekin; 
  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);
  
  //primary particle info
  //
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;
  
  // accumulate sums
  //
  fEnergyDeposit   += localRun->fEnergyDeposit;  
  fEnergyDeposit2  += localRun->fEnergyDeposit2;
  fEnergyFlow      += localRun->fEnergyFlow;
  fEnergyFlow2     += localRun->fEnergyFlow2;

  //additional parameters
  //TotalCount += localRun->TotalCount;
  //
  EnergyDeposit    += localRun->EnergyDeposit;  
  EnergyDeposit2   += localRun->EnergyDeposit2;  
  NonIonEnergyDeposit += localRun->NonIonEnergyDeposit;
  NonIonEnergyDeposit2 += localRun->NonIonEnergyDeposit2;
  sum_TL	   += localRun->sum_TL; 
  sum_TL2	   += localRun->sum_TL2;
  sum_T	  	   += localRun->sum_T; 
  sum_T2	   += localRun->sum_T2; 
  sum_Tt	   += localRun->sum_Tt; 
  sum_Tt2	   += localRun->sum_Tt2; 
  N_rec		   += localRun->N_rec; 
  N_Sec_Rec 	   += localRun->N_Sec_Rec;
  Nsteps	   += localRun->Nsteps;
  Nsteps2	   += localRun->Nsteps2;
  theta		   += localRun->theta;
  theta2	   += localRun->theta2;
  TrakLenPrim      += localRun->TrakLenPrim;
  TrakLenPrim2 	   += localRun->TrakLenPrim2;
  TrakLenSec       += localRun->TrakLenSec;
  TrakLenSec2 	   += localRun->TrakLenSec2;

  fTrueRange  += localRun->fTrueRange;
  fTrueRange2 += localRun->fTrueRange2;
  fProjRange  += localRun->fProjRange;
  fProjRange2 += localRun->fProjRange2;

  // ABSORBER
  abslen		+= localRun->abslen;
  abslen2 		+= localRun->abslen2;
  absstep		+= localRun->absstep;
  absstep2		+= localRun->absstep2;
  absfEdep      	+= localRun->absfEdep;
  absniel		+= localRun->absniel;
  absfEdep2     	+= localRun->absfEdep2;
  absniel2		+= localRun->absniel2;
  abssum_tl		+= localRun->abssum_tl;
  abssum_tl2		+= localRun->abssum_tl2;
  abssum_t		+= localRun->abssum_t;
  abssum_t2		+= localRun->abssum_t2;
  absrec 		+= localRun->absrec;
  
  abs1len		+= localRun->abs1len;
  abs1len2 		+= localRun->abs1len2;
  abs1step		+= localRun->abs1step;
  abs1step2		+= localRun->abs1step2;
  abs1fEdep      	+= localRun->abs1fEdep;
  abs1niel		+= localRun->abs1niel;
  abs1fEdep2     	+= localRun->abs1fEdep2;
  abs1niel2		+= localRun->abs1niel2;
  abs1sum_tl		+= localRun->abs1sum_tl;
  abs1sum_tl2		+= localRun->abs1sum_tl2;
  abs1sum_t		+= localRun->abs1sum_t;
  abs1sum_t2		+= localRun->abs1sum_t2;
  abs1rec 		+= localRun->abs1rec;  
  
  // other material
  abs2len		+= localRun->abs2len;
  abs2len2 		+= localRun->abs2len2;
  abs2step		+= localRun->abs2step;
  abs2step2		+= localRun->abs2step2;
  abs2fEdep      	+= localRun->abs2fEdep;
  abs2niel		+= localRun->abs2niel;
  abs2fEdep2     	+= localRun->abs2fEdep2;
  abs2niel2		+= localRun->abs2niel2;
  abs2sum_tl		+= localRun->abs2sum_tl;
  abs2sum_tl2		+= localRun->abs2sum_tl2;
  abs2sum_t		+= localRun->abs2sum_t;
  abs2sum_t2		+= localRun->abs2sum_t2;
  abs2rec 		+= localRun->abs2rec;
  // in between
  abs3len		+= localRun->abs3len;
  abs3len2 		+= localRun->abs3len2;
  abs3step		+= localRun->abs3step;
  abs3step2		+= localRun->abs3step2;
  abs3fEdep      	+= localRun->abs3fEdep;
  abs3niel		+= localRun->abs3niel;
  abs3fEdep2     	+= localRun->abs3fEdep2;
  abs3niel2		+= localRun->abs3niel2;
  abs3sum_tl		+= localRun->abs3sum_tl;
  abs3sum_tl2		+= localRun->abs3sum_tl2;
  abs3sum_t		+= localRun->abs3sum_t;
  abs3sum_t2		+= localRun->abs3sum_t2;
  abs3rec 		+= localRun->abs3rec;
  
  abs4len		+= localRun->abs4len;
  abs4len2 		+= localRun->abs4len2;
  abs4step		+= localRun->abs4step;
  abs4step2		+= localRun->abs4step2;
  abs4fEdep      	+= localRun->abs4fEdep;
  abs4niel		+= localRun->abs4niel;
  abs4fEdep2     	+= localRun->abs4fEdep2;
  abs4niel2		+= localRun->abs4niel2;
  abs4sum_tl		+= localRun->abs4sum_tl;
  abs4sum_tl2		+= localRun->abs4sum_tl2;
  abs4sum_t		+= localRun->abs4sum_t;
  abs4sum_t2		+= localRun->abs4sum_t2;
  abs4rec 		+= localRun->abs4rec;  
  
  primary_energy	+= localRun->primary_energy;
  
  angle_of_incidence 	+= localRun->angle_of_incidence;


  totniel       += localRun->totniel;
  totniel2      += localRun->totniel2;

  hit		+= localRun->hit;
  partEmerging  += localRun->partEmerging;

  detKinEn 	+= localRun->detKinEn;
  detKinEn2	+= localRun->detKinEn2;

  rbs 		+= localRun->rbs;
  second 	+= localRun->second;
  nielstep	+= localRun->nielstep;
  totstep	+= localRun->totstep;


  projectedR	+=localRun->projectedR;
  projectedR2	+=localRun->projectedR2;


  entry_sd 	+=localRun->entry_sd;
  entry_reach	+=localRun->entry_reach;
  total_step	+=localRun->total_step;
  RBSDepth	+=localRun->RBSDepth;
  counts 	+=localRun->counts;

  //tothit 	+= localRun->hit;
  // END

  //map: processes count
  std::map<G4String,G4int>::const_iterator itp;
  for ( itp = localRun->fProcCounter.begin();
        itp != localRun->fProcCounter.end(); ++itp ) {

    G4String procName = itp->first;
    G4int localCount = itp->second;
    if ( fProcCounter.find(procName) == fProcCounter.end()) {
      fProcCounter[procName] = localCount;
    }
    else {
      fProcCounter[procName] += localCount;
    }  
  }
  
  //map: created particles count    
  std::map<G4String,ParticleData>::const_iterator itc;
  for (itc = localRun->fParticleDataMap1.begin(); 
       itc != localRun->fParticleDataMap1.end(); ++itc) {
    
    G4String name = itc->first;
    const ParticleData& localData = itc->second;   
    if ( fParticleDataMap1.find(name) == fParticleDataMap1.end()) {
      fParticleDataMap1[name]
       = ParticleData(localData.fCount, 
                      localData.fEmean, 
                      localData.fEmin, 
                      localData.fEmax);
    }
    else {
      ParticleData& data = fParticleDataMap1[name];   
      data.fCount += localData.fCount;
      data.fEmean += localData.fEmean;
      G4double emin = localData.fEmin;
      if (emin < data.fEmin) data.fEmin = emin;
      G4double emax = localData.fEmax;
      if (emax > data.fEmax) data.fEmax = emax; 
    }   
  }
  
  //map: particles flux count       
  std::map<G4String,ParticleData>::const_iterator itn;
  for (itn = localRun->fParticleDataMap2.begin(); 
       itn != localRun->fParticleDataMap2.end(); ++itn) {
    
    G4String name = itn->first;
    const ParticleData& localData = itn->second;   
    if ( fParticleDataMap2.find(name) == fParticleDataMap2.end()) {
      fParticleDataMap2[name]
       = ParticleData(localData.fCount, 
                      localData.fEmean, 
                      localData.fEmin, 
                      localData.fEmax);
    }
    else {
      ParticleData& data = fParticleDataMap2[name];   
      data.fCount += localData.fCount;
      data.fEmean += localData.fEmean;
      G4double emin = localData.fEmin;
      if (emin < data.fEmin) data.fEmin = emin;
      G4double emax = localData.fEmax;
      if (emax > data.fEmax) data.fEmax = emax; 
    }   
  }
    const std::map< G4int, G4THitsMap<G4double>* >& localMap = localRun->fMap;
    std::map< G4int, G4THitsMap<G4double>* >::const_iterator iter = localMap.begin();
    for ( ; iter != localMap.end() ; ++iter)
        (*(fMap[iter->first])) += (*(iter->second));


  G4Run::Merge(run); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun() 
{



  G4int prec = 5, wid = prec + 2;  
  G4int dfprec = G4cout.precision(prec);
  //fDetector->GetMaterial()->GetName()

  G4Material* material = fDetector->GetMaterialM(0);
  G4double density = material->GetDensity();

  G4Material* material2 = fDetector->GetMaterialM(1);
  G4double density2 = material2->GetDensity();
  G4Material* material3 = fDetector->GetMaterialM(2);
  G4double density3 = material3->GetDensity();
  G4Material* material4 = fDetector->GetMaterialM(3);  
  G4double density4 = material4->GetDensity();
  G4Material* material5 = fDetector->GetMaterialM(4);  
  G4double density5 = material5->GetDensity();
  
  G4cout << " ********************************************************** " << G4endl;
  G4cout << " ******************** FINISH OF RUN *********************** " << G4endl;
     
  G4String Particle = fParticle->GetParticleName();    
  G4cout << "\n The run is " << numberOfEvent << " "<< Particle << " of "
         << G4BestUnit(fEkin,"Energy") << " through "  << G4BestUnit(fDetector->GetLength(0),"Length") << " length parallelpiped of " << material->GetName() << " (density: " << G4BestUnit(density,"Volumic Mass") << ")" << "\n and second material " << G4BestUnit(fDetector->GetLength(1),"Length") << " of " << material2->GetName() <<  " (density: " << G4BestUnit(density2,"Volumic Mass") << ")" << G4endl;

  if (numberOfEvent == 0) { G4cout.precision(dfprec);   return;}
             
  //frequency of processes
  //
  G4cout << "\n Process calls frequency :" << G4endl;
  G4int index = 0;
  std::map<G4String,G4int>::iterator it;    
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first;
     G4int    count    = it->second;
     G4String space = " "; if (++index%3 == 0) space = "\n";
     G4cout << " " << std::setw(20) << procName << "="<< std::setw(7) << count
            << space;
  }
  G4cout << G4endl;
  
  //particles count
  //
  G4cout << "\n List of generated particles:" << G4endl;
     
 std::map<G4String,ParticleData>::iterator itc;               
 for (itc = fParticleDataMap1.begin(); itc != fParticleDataMap1.end(); itc++) { 
    G4String name = itc->first;
    ParticleData data = itc->second;
    G4int count = data.fCount;
    G4double eMean = data.fEmean/count;
    G4double eMin = data.fEmin;
    G4double eMax = data.fEmax;    
         
    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") 
           << ")" << G4endl;           
 }
   
  // compute mean Energy deposited and rms
  //
  G4int TotNbofEvents = numberOfEvent;
  fEnergyDeposit /= TotNbofEvents; fEnergyDeposit2 /= TotNbofEvents;
  G4double rmsEdep = fEnergyDeposit2 - fEnergyDeposit*fEnergyDeposit;
  if (rmsEdep>0.) rmsEdep = std::sqrt(rmsEdep);
  else            rmsEdep = 0.;
  
  G4cout << "\n Mean energy deposit per event = "
         << G4BestUnit(fEnergyDeposit,"Energy") << ";  rms = "
         << G4BestUnit(rmsEdep,      "Energy") 
         << G4endl;
  
  // compute mean Energy flow and rms
  //
  fEnergyFlow /= TotNbofEvents; fEnergyFlow2 /= TotNbofEvents;
  G4double rmsEflow = fEnergyFlow2 - fEnergyFlow*fEnergyFlow;
  if (rmsEflow>0.) rmsEflow = std::sqrt(rmsEflow);
  else             rmsEflow = 0.;
  
  G4cout << " Mean energy flow per event    = "
         << G4BestUnit(fEnergyFlow,"Energy") << ";  rms = "
         << G4BestUnit(rmsEflow,   "Energy") 
         << G4endl;
                                
 //particles flux
 //
 G4cout << "\n List of particles emerging from the absorber :" << G4endl;
     
 std::map<G4String,ParticleData>::iterator itn;               
 for (itn = fParticleDataMap2.begin(); itn != fParticleDataMap2.end(); itn++) { 
    G4String name = itn->first;
    ParticleData data = itn->second;
    G4int count = data.fCount;
    G4double eMean = data.fEmean/count;
    G4double eMin = data.fEmin;
    G4double eMax = data.fEmax;
    G4double Eflow = data.fEmean/TotNbofEvents;        
         
    G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
           << "  Emean = " << std::setw(wid) << G4BestUnit(eMean, "Energy")
           << "\t( "  << G4BestUnit(eMin, "Energy")
           << " --> " << G4BestUnit(eMax, "Energy") 
           << ") \tEflow/event = " << G4BestUnit(Eflow, "Energy") << G4endl;
 }

  // compute mean and rms
  if (TotNbofEvents == 0) return;
  
  //total energy loss
  EnergyDeposit /= TotNbofEvents; EnergyDeposit2 /= TotNbofEvents;
  G4double FrmsEdep = EnergyDeposit2 - EnergyDeposit*EnergyDeposit;
  if(FrmsEdep >0.)FrmsEdep=std::sqrt(FrmsEdep/TotNbofEvents);
  else FrmsEdep =0;

  //nuclear energy loss
  NonIonEnergyDeposit /= TotNbofEvents; NonIonEnergyDeposit2 /= TotNbofEvents;
  G4double rmsEnondep = NonIonEnergyDeposit2 - NonIonEnergyDeposit*NonIonEnergyDeposit;
  if(rmsEnondep>0.) rmsEnondep= std::sqrt(rmsEnondep/TotNbofEvents);
  else rmsEnondep=0;

  //mean sum of T( kinetic energy of secondary)x L(T) partition energy 
  sum_TL/=TotNbofEvents;     sum_TL2/=TotNbofEvents;
  G4double rmssum_TL =sum_TL2- sum_TL*sum_TL;
  if(rmssum_TL>0.) rmssum_TL=std::sqrt(rmssum_TL/TotNbofEvents);
  else rmssum_TL =0;

  //mean kinetic energy of secondary particles (IDp==1) 
  G4double rmssum_T = 0.0;
  if(N_rec > 0) {
    sum_T/=N_rec;     sum_T2/=N_rec;
    rmssum_T =sum_T2- sum_T*sum_T;
    if(rmssum_T>0.) rmssum_T=std::sqrt(rmssum_T/N_rec);  }

  //mean kinetic energy of tertiary particles (IDp>1) 
  G4double rmssum_Tt = 0.0;
  if(N_Sec_Rec > 0) {
    sum_Tt/=N_Sec_Rec;     sum_Tt2/=N_Sec_Rec;
    rmssum_Tt =sum_Tt2- sum_Tt*sum_Tt;
    if(rmssum_Tt>0.) rmssum_Tt=std::sqrt(rmssum_Tt/N_Sec_Rec);  }

  //mean number of steps:
  Nsteps/=TotNbofEvents;  Nsteps2/=TotNbofEvents;
  G4double rmsSteps= Nsteps2 -Nsteps*Nsteps;
  if(rmsSteps>0) rmsSteps= std::sqrt(rmsSteps/TotNbofEvents);
  else rmsSteps=0;

  //scattering angle
  theta/=TotNbofEvents ; theta2/=TotNbofEvents;	
  G4double rmsTheta =theta2-theta*theta;
  if(rmsTheta>0.) rmsTheta =std::sqrt(rmsTheta/TotNbofEvents);
  else rmsTheta =0;
  
  //track length
  TrakLenPrim /= TotNbofEvents; TrakLenPrim2 /= TotNbofEvents;
  G4double rmsTLPrim = TrakLenPrim2 - TrakLenPrim*TrakLenPrim;
  if (rmsTLPrim>0.) rmsTLPrim = std::sqrt(rmsTLPrim/TotNbofEvents);
  else rmsTLPrim = 0.;
  // secondaries track length
  TrakLenSec /= N_rec; TrakLenSec2 /= N_rec;
  G4double rmsTLSec = TrakLenSec2 - TrakLenSec*TrakLenSec;
  if (rmsTLSec>0.) rmsTLSec = std::sqrt(rmsTLSec/N_rec);
  else rmsTLSec = 0.;

  // AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  //************************************************
  // --------------------------------------------
  // MAIN LAYER
  //  edep
  absfEdep /= TotNbofEvents; absfEdep2 /= TotNbofEvents;
  G4double absrmsEdep = absfEdep2 - absfEdep*absfEdep;
  if(absrmsEdep >0.)absrmsEdep=std::sqrt(absrmsEdep/TotNbofEvents);
  else absrmsEdep =0;
  // track length 
  abslen /= TotNbofEvents; abslen2 /= TotNbofEvents;
  G4double absrmsTLPrim = abslen2 - abslen*abslen;
  if (absrmsTLPrim>0.) absrmsTLPrim = std::sqrt(absrmsTLPrim/TotNbofEvents);
  else absrmsTLPrim = 0.;
  //nuclear energy loss
  absniel /= TotNbofEvents; absniel2 /= TotNbofEvents;
  G4double absrmsEnondep = absniel2 - absniel*absniel;
  if(absrmsEnondep>0.) absrmsEnondep= std::sqrt(absrmsEnondep/TotNbofEvents);
  else absrmsEnondep=0;
  //mean sum of T( kinetic energy of secondary)x L(T) partition energy 
  abssum_tl/=TotNbofEvents;     abssum_tl2/=TotNbofEvents;
  G4double absrmssum_tl =abssum_tl2- abssum_tl*abssum_tl;
  if(absrmssum_tl>0.) absrmssum_tl=std::sqrt(absrmssum_tl/TotNbofEvents);
  else absrmssum_tl =0;
  //mean kinetic energy of secondary particles (IDp==1) 
  G4double absrmssum_t = 0.0;
  if(absrec > 0) {
    abssum_t/=absrec;     abssum_t2/=absrec;
    absrmssum_t =abssum_t2- abssum_t*abssum_t;
    if(absrmssum_t>0.) absrmssum_t=std::sqrt(absrmssum_t/absrec);  }
  //mean number of steps:
  absstep/=TotNbofEvents;  absstep2/=TotNbofEvents;
  G4double absrmsSteps= absstep2 -absstep*absstep;
  if(absrmsSteps>0) absrmsSteps= std::sqrt(absrmsSteps/TotNbofEvents);
  else absrmsSteps=0;
	// mean dE/dx in absorber
  G4double dEdx_abs = absfEdep/abslen;
	// mean niel dE/dx in absorber
  G4double dEdx_niel_abs = absniel/abslen;

  // **************************************
  // FIRST LAYER
  //  edep
  abs1fEdep /= TotNbofEvents; abs1fEdep2 /= TotNbofEvents;
  G4double abs1rmsEdep = abs1fEdep2 - abs1fEdep*abs1fEdep;
  if(abs1rmsEdep >0.)abs1rmsEdep=std::sqrt(abs1rmsEdep/TotNbofEvents);
  else abs1rmsEdep =0;
  // track length 
  abs1len /= TotNbofEvents; abs1len2 /= TotNbofEvents;
  G4double abs1rmsTLPrim = abs1len2 - abs1len*abs1len;
  if (abs1rmsTLPrim>0.) abs1rmsTLPrim = std::sqrt(abs1rmsTLPrim/TotNbofEvents);
  else abs1rmsTLPrim = 0.;
  //nuclear energy loss
  abs1niel /= TotNbofEvents; abs1niel2 /= TotNbofEvents;
  G4double abs1rmsEnondep = abs1niel2 - abs1niel*abs1niel;
  if(abs1rmsEnondep>0.) abs1rmsEnondep= std::sqrt(abs1rmsEnondep/TotNbofEvents);
  else abs1rmsEnondep=0;
  //mean sum of T( kinetic energy of secondary)x L(T) partition energy 
  abs1sum_tl/=TotNbofEvents;     abs1sum_tl2/=TotNbofEvents;
  G4double abs1rmssum_tl =abs1sum_tl2- abs1sum_tl*abs1sum_tl;
  if(abs1rmssum_tl>0.) abs1rmssum_tl=std::sqrt(abs1rmssum_tl/TotNbofEvents);
  else abs1rmssum_tl =0;
  //mean kinetic energy of secondary particles (IDp==1) 
  G4double abs1rmssum_t = 0.0;
  if(abs1rec > 0) {
    abs1sum_t/=abs1rec;     abs1sum_t2/=abs1rec;
    abs1rmssum_t =abs1sum_t2- abs1sum_t*abs1sum_t;
    if(abs1rmssum_t>0.) abs1rmssum_t=std::sqrt(abs1rmssum_t/abs1rec);  }
  //mean number of steps:
  abs1step/=TotNbofEvents;  abs1step2/=TotNbofEvents;
  G4double abs1rmsSteps= abs1step2 -abs1step*abs1step;
  if(abs1rmsSteps>0) abs1rmsSteps= std::sqrt(abs1rmsSteps/TotNbofEvents);
  else abs1rmsSteps=0;
	// mean dE/dx in absorber 1
  G4double dEdx_abs1 = abs1fEdep/abs1len;
	// mean niel dE/dx in absorber 1 
  G4double dEdx_niel_abs1 = abs1niel/abs1len;


//*********************************
  // SECOND LAYER
  //  edep
  abs2fEdep /= TotNbofEvents; abs2fEdep2 /= TotNbofEvents;
  G4double abs2rmsEdep = abs2fEdep2 - abs2fEdep*abs2fEdep;
  if(abs2rmsEdep >0.)abs2rmsEdep=std::sqrt(abs2rmsEdep/TotNbofEvents);
  else abs2rmsEdep =0;
  // track length 
  abs2len /= TotNbofEvents; abs2len2 /= TotNbofEvents;
  G4double abs2rmsTLPrim = abs2len2 - abs2len*abs2len;
  if (abs2rmsTLPrim>0.) abs2rmsTLPrim = std::sqrt(abs2rmsTLPrim/TotNbofEvents);
  else abs2rmsTLPrim = 0.;
  //nuclear energy loss
  abs2niel /= TotNbofEvents; abs2niel2 /= TotNbofEvents;
  G4double abs2rmsEnondep = abs2niel2 - abs2niel*abs2niel;
  if(abs2rmsEnondep>0.) abs2rmsEnondep= std::sqrt(abs2rmsEnondep/TotNbofEvents);
  else abs2rmsEnondep=0;
  //mean sum of T( kinetic energy of secondary)x L(T) partition energy 
  abs2sum_tl/=TotNbofEvents;     abs2sum_tl2/=TotNbofEvents;
  G4double abs2rmssum_tl =abs2sum_tl2- abs2sum_tl*abs2sum_tl;
  if(abs2rmssum_tl>0.) abs2rmssum_tl=std::sqrt(abs2rmssum_tl/TotNbofEvents);
  else abs2rmssum_tl =0;
  //mean kinetic energy of secondary particles (IDp==1) 
  G4double abs2rmssum_t = 0.0;
  if(abs2rec > 0) {
    abs2sum_t/=abs2rec;     abs2sum_t2/=abs2rec;
    abs2rmssum_t =abs2sum_t2- abs2sum_t*abs2sum_t;
    if(abs2rmssum_t>0.) abs2rmssum_t=std::sqrt(abs2rmssum_t/abs2rec);  }
  //mean number of steps:
  abs2step/=TotNbofEvents;  abs2step2/=TotNbofEvents;
  G4double abs2rmsSteps= abs2step2 -abs2step*abs2step;
  if(abs2rmsSteps>0) abs2rmsSteps= std::sqrt(abs2rmsSteps/TotNbofEvents);
  else abs2rmsSteps=0;
	// mean dE/dx in absorber 2
  G4double dEdx_abs2 = abs2fEdep/abs2len;
	// mean niel dE/dx in absorber 2 
  G4double dEdx_niel_abs2 = abs2niel/abs2len;


//*********************************
  // THIRD LAYER
  //  edep
  abs3fEdep /= TotNbofEvents; abs3fEdep2 /= TotNbofEvents;
  G4double abs3rmsEdep = abs3fEdep2 - abs3fEdep*abs3fEdep;
  if(abs3rmsEdep >0.)abs3rmsEdep=std::sqrt(abs3rmsEdep/TotNbofEvents);
  else abs3rmsEdep =0;
  // track length 
  abs3len /= TotNbofEvents; abs3len2 /= TotNbofEvents;
  G4double abs3rmsTLPrim = abs3len2 - abs3len*abs3len;
  if (abs3rmsTLPrim>0.) abs3rmsTLPrim = std::sqrt(abs3rmsTLPrim/TotNbofEvents);
  else abs3rmsTLPrim = 0.;
  //nuclear energy loss
  abs3niel /= TotNbofEvents; abs3niel2 /= TotNbofEvents;
  G4double abs3rmsEnondep = abs3niel2 - abs3niel*abs3niel;
  if(abs3rmsEnondep>0.) abs3rmsEnondep= std::sqrt(abs3rmsEnondep/TotNbofEvents);
  else abs3rmsEnondep=0;
  //mean sum of T( kinetic energy of secondary)x L(T) partition energy 
  abs3sum_tl/=TotNbofEvents;     abs3sum_tl2/=TotNbofEvents;
  G4double abs3rmssum_tl =abs3sum_tl2- abs3sum_tl*abs3sum_tl;
  if(abs3rmssum_tl>0.) abs3rmssum_tl=std::sqrt(abs3rmssum_tl/TotNbofEvents);
  else abs3rmssum_tl =0;
  //mean kinetic energy of secondary particles (IDp==1) 
  G4double abs3rmssum_t = 0.0;
  if(abs3rec > 0) {
    abs3sum_t/=abs3rec;     abs3sum_t2/=abs3rec;
    abs3rmssum_t =abs3sum_t2- abs3sum_t*abs3sum_t;
    if(abs3rmssum_t>0.) abs3rmssum_t=std::sqrt(abs3rmssum_t/abs3rec);  }
  //mean number of steps:
  abs3step/=TotNbofEvents;  abs3step2/=TotNbofEvents;
  G4double abs3rmsSteps= abs3step2 -abs3step*abs3step;
  if(abs3rmsSteps>0) abs3rmsSteps= std::sqrt(abs3rmsSteps/TotNbofEvents);
  else abs3rmsSteps=0;
	// mean dE/dx in absorber 3
  G4double dEdx_abs3 = abs3fEdep/abs3len;
	// mean niel dE/dx in absorber 3 
  G4double dEdx_niel_abs3 = abs3niel/abs3len;


//*********************************
  // FOURTH LAYER
  //  edep
  abs4fEdep /= TotNbofEvents; abs4fEdep2 /= TotNbofEvents;
  G4double abs4rmsEdep = abs4fEdep2 - abs4fEdep*abs4fEdep;
  if(abs4rmsEdep >0.)abs4rmsEdep=std::sqrt(abs4rmsEdep/TotNbofEvents);
  else abs4rmsEdep =0;
  // track length 
  abs4len /= TotNbofEvents; abs4len2 /= TotNbofEvents;
  G4double abs4rmsTLPrim = abs4len2 - abs4len*abs4len;
  if (abs4rmsTLPrim>0.) abs4rmsTLPrim = std::sqrt(abs4rmsTLPrim/TotNbofEvents);
  else abs4rmsTLPrim = 0.;
  //nuclear energy loss
  abs4niel /= TotNbofEvents; abs4niel2 /= TotNbofEvents;
  G4double abs4rmsEnondep = abs4niel2 - abs4niel*abs4niel;
  if(abs4rmsEnondep>0.) abs4rmsEnondep= std::sqrt(abs4rmsEnondep/TotNbofEvents);
  else abs4rmsEnondep=0;
  //mean sum of T( kinetic energy of secondary)x L(T) partition energy 
  abs4sum_tl/=TotNbofEvents;     abs4sum_tl2/=TotNbofEvents;
  G4double abs4rmssum_tl =abs4sum_tl2- abs4sum_tl*abs4sum_tl;
  if(abs4rmssum_tl>0.) abs4rmssum_tl=std::sqrt(abs4rmssum_tl/TotNbofEvents);
  else abs4rmssum_tl =0;
  //mean kinetic energy of secondary particles (IDp==1) 
  G4double abs4rmssum_t = 0.0;
  if(abs4rec > 0) {
    abs4sum_t/=abs4rec;     abs4sum_t2/=abs4rec;
    abs4rmssum_t =abs4sum_t2- abs4sum_t*abs4sum_t;
    if(abs4rmssum_t>0.) abs4rmssum_t=std::sqrt(abs4rmssum_t/abs4rec);  }
  //mean number of steps:
  abs4step/=TotNbofEvents;  abs4step2/=TotNbofEvents;
  G4double abs4rmsSteps= abs4step2 -abs4step*abs4step;
  if(abs4rmsSteps>0) abs4rmsSteps= std::sqrt(abs4rmsSteps/TotNbofEvents);
  else abs4rmsSteps=0;
	// mean dE/dx in absorber 4
  G4double dEdx_abs4 = abs4fEdep/abs4len;
	// mean niel dE/dx in absorber 4 
  G4double dEdx_niel_abs4 = abs4niel/abs4len;

  //*************************************************

  //nuclear energy loss
  totniel /= TotNbofEvents; totniel2 /= TotNbofEvents;
  G4double totrmsEnondep = totniel2 - totniel*totniel;
  if(totrmsEnondep>0.) totrmsEnondep= std::sqrt(totrmsEnondep/TotNbofEvents);
  else totrmsEnondep=0;


  //RBSDepth /= TotNbofEvents; RBSDepth2 /= TotNbofEvents;
  G4double rmsRBSDepth = RBSDepth2 - RBSDepth*RBSDepth;
  if (rmsRBSDepth>0.) rmsRBSDepth = std::sqrt(rmsRBSDepth);
  else rmsRBSDepth=0.;





  //..............................................................

  //G4double thickness  = fDetector->GetRadius();
  //Stopping Power and NIEL from simulation.
  // effective length
  G4double length=TrakLenPrim;
  // total energy loss  
  G4double meandEdx  = EnergyDeposit/length;
  // nuclear energy loss
  G4double meandEdx_nucl  = NonIonEnergyDeposit/length;
  // NIEL 
  G4double meandEdx_sumTL=sum_TL/length;

  //G4double RealMFP = TrakLenPrim/TotalCount;

	//[MeVcm2/g]
  G4double stopPower = meandEdx/density;  
  G4double stopPower_nucl = meandEdx_nucl/density;
  G4double stopPower_sumTL=meandEdx_sumTL/density;
 
 	//mean free path & corss section 
  G4double freepath= TrakLenPrim/Nsteps;
  G4double er1=rmsTLPrim/Nsteps;	
  G4double er2=freepath*rmsSteps/Nsteps;
  G4double rmsFreepath=std::sqrt(er1*er1+er2*er2);


  G4double NA  = material->GetTotNbOfAtomsPerVolume();
  G4double CrossSection =1./(NA*freepath); 
  G4double rmsCrossSection=rmsFreepath*CrossSection/freepath;


  // true and projected range from TestEm1
  fTrueRange /= TotNbofEvents; fTrueRange2 /= TotNbofEvents;
  G4double trueRms = fTrueRange2 - fTrueRange*fTrueRange;        
  if (trueRms>0.) trueRms = std::sqrt(trueRms); else trueRms = 0.;
        
  fProjRange /= TotNbofEvents; fProjRange2 /= TotNbofEvents;
  G4double projRms = fProjRange2 - fProjRange*fProjRange;        
  if (projRms>0.) projRms = std::sqrt(projRms); else projRms = 0.;

	// projected range from testem11
  //compute projected range of primary track
  //
  G4double rmsPR;
  projectedR /= TotNbofEvents; projectedR2 /= TotNbofEvents;
  rmsPR = projectedR2 - projectedR*projectedR;        
  if (rmsPR>0.) rmsPR = std::sqrt(rmsPR); else rmsPR = 0.;




  G4cout << "\n ============= Simulation statistics ==============\n";
  G4cout << "\n Primary Total track length in absorber:\n "
         << G4BestUnit(TrakLenPrim,"Length") << " +- "
         << G4BestUnit(rmsTLPrim,       "Length") << G4endl;

  G4cout << "\n Secondaries total track length in absorber:\n "
         << G4BestUnit(TrakLenSec,"Length") << " +- "
         << G4BestUnit(rmsTLSec,       "Length") << G4endl;

  G4cout << "\n ==================================================\n ";
  G4cout << "\n Primary particle statistics\n ";
  G4cout << "\n Mean Number of Steps:\n "<<Nsteps<< " +/- "<< rmsSteps<<G4endl;
								
  G4cout << "\n Mean Free Path :\n "<<G4BestUnit(freepath,"Length")<<  
				" +/- "<< G4BestUnit(rmsFreepath,"Length")<<G4endl;

  G4cout << "\n Mean Cross Section :\n "<<G4BestUnit(CrossSection,"Surface")<<
			" +/- "<<G4BestUnit(rmsCrossSection,"Surface")<<G4endl;

  G4cout << "\n Mean scattering angle :\n "<<G4BestUnit(theta,"Angle")<< " +/- "
			<< G4BestUnit(rmsTheta,"Angle")<<G4endl;
  
  G4cout << "\n Total energy deposit in absorber:\n "
         << G4BestUnit(EnergyDeposit,"Energy") << " +/- "
         << G4BestUnit(FrmsEdep,      "Energy") 
         << G4endl;
  G4cout << "-----> dE/dx total= " << meandEdx/(MeV/cm) << " MeV/cm"
         << "\t(" << stopPower/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;


  G4cout << "\n Nuclear energy deposit in absorber:\n "
         << G4BestUnit(NonIonEnergyDeposit,"Energy") << " +/- "
         << G4BestUnit(rmsEnondep,      "Energy")
         << G4endl;
  G4cout << "-----> dE/dx  nucl = " << meandEdx_nucl/(MeV/cm) << " MeV/cm"
         << "\t(" << stopPower_nucl/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;

  G4cout <<"\n NIEL in absorber (Th>"<<Th/eV <<" eV):\n "
         << G4BestUnit(sum_TL,"Energy") << " +/- "
         << G4BestUnit(rmssum_TL,      "Energy")
         << G4endl;
  G4cout << "-----> NIEL = " << meandEdx_sumTL/(MeV/cm) << " MeV/cm"
         << "\t(" << stopPower_sumTL/(MeV*cm2/g) << " MeV*cm2/g)"
         << G4endl;

  G4cout << "\n ===========================================================\n";
  G4cout << "\n true Range = " << G4BestUnit(fTrueRange,"Length")
         << "   rms = "        << G4BestUnit(trueRms,  "Length");

  G4cout << "\n proj Range = " << G4BestUnit(fProjRange,"Length")
         << "   rms = "        << G4BestUnit(projRms,  "Length");
  G4cout << "\n ===========================================================\n";
  G4cout << " ====================== GEOMETRY ===========================";
  G4cout << "\n ===========================================================\n";
  G4cout << " \n"; 
  
    G4ThreeVector pozicija1 = fDetector->GetPosition(0);
    G4ThreeVector pozicija2 = fDetector->GetPosition(1);
    G4ThreeVector pozicija3 = fDetector->GetPosition(2);
    G4ThreeVector pozicija4 = fDetector->GetPosition(3);
    
    G4double ilgis1 = pozicija1.z();
    G4double ilgis2 = pozicija2.z();
    G4double ilgis3 = pozicija3.z();
    G4double ilgis4 = pozicija4.z();

  G4double atstums1 = (fDetector->GetLength(0)/2)+ilgis1-(fDetector->GetLength(1)/2);
  if(atstums1/nm < 1e-5)
  {
  	atstums1 = 0.;
  }
  G4double atstums2 = (fDetector->GetLength(0)/2)+ilgis2-(fDetector->GetLength(2)/2);  
  G4double atstums3 = (fDetector->GetLength(0)/2)+ilgis3-(fDetector->GetLength(3)/2);  
  G4double atstums4 = (fDetector->GetLength(0)/2)+ilgis4-(fDetector->GetLength(4)/2);  
  
  G4cout << " >>>>>>>>>>>>>>>>> MAIN LAYER <<<<<<<<<<<<<<<<<<<<<<<<\n";
  G4cout << " Thickness = " << G4BestUnit(fDetector->GetLength(0),"Length") << " or surface layer thickness " << (material->GetTotNbOfAtomsPerVolume()/(1/cm3)*atstums1/cm)/(1e+15) << "[TFU] " << G4endl;
  G4cout << " Total layer dimensions : " << G4BestUnit(fDetector->GetSize(0), "Length") << G4endl;	
  G4cout << " Material: " << fDetector->GetMaterial(0) << G4endl;
  G4cout << " Density = " << density/(g/cm3) << " g/cm3 "  << ", or " << material->GetTotNbOfAtomsPerVolume()/(1/cm3) << " [at3] "  << G4endl;  
  G4cout << " Composition: " << G4endl;
  const G4double *atomDensVector1	= material->GetVecNbOfAtomsPerVolume();
  G4double NoOfElements_material    	= material->GetNumberOfElements();
  for (int i=0; i<NoOfElements_material; i++){
  G4cout << " " << material->GetElement(i)->GetName() << " = " << atomDensVector1[i]/(1/cm3) << " [cm-3] " << G4endl;}
  G4cout << " -------------------- Particle parameters -------------------- " << G4endl;
  G4cout << " Ionizing energy loss : \t\t" << G4BestUnit(absfEdep, "Energy") << " +/- " << G4BestUnit(absrmsEdep, "Energy")  << G4endl;	
  G4cout << " Non ionizing energy loss [NIEL] : \t" << G4BestUnit(absniel, "Energy") << " +/- " << G4BestUnit(absrmsEnondep, "Energy")  << G4endl;	
  G4cout << " Mean number of steps : \t\t" << absstep << " +/- " << absrmsSteps << G4endl;
  G4cout << " Primary track length : \t\t" << G4BestUnit(abslen, "Length") << " +/- " << G4BestUnit(absrmsTLPrim,"Length") << G4endl; 
  G4cout << " Mean dE/dx : \t\t\t\t" << dEdx_abs/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean NIEL dE/dx: \t\t\t" << dEdx_niel_abs/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean free path: \t\t\t" <<G4BestUnit(abslen/absstep,"Length")<<  G4endl;
  G4cout << " Steps per length [nm-1]: \t\t" << absstep/(fDetector->GetLength(0)/nm) << G4endl;
  G4cout << " # of secondaries: \t\t\t" << absrec << " \t per primary : " << absrec/TotNbofEvents << G4endl;
  G4cout << " Mean E of secondaries: \t\t" << G4BestUnit(abssum_t, "Energy") << " +/- " << G4BestUnit(absrmssum_t,"Energy") << G4endl;  
  G4cout << " Mean damage en of secondaries: \t" << G4BestUnit(abssum_tl, "Energy") << " +/- " << G4BestUnit(absrmssum_tl,"Energy") << G4endl;    
  G4cout << " \n"; 
  
  G4double NoOfElements_material2    	= material2->GetNumberOfElements();
  G4cout << " >>>>>>>>>>>>>>>>> FIRST LAYER <<<<<<<<<<<<<<<<<<<<<<<<\n";
  G4cout << " Thickness = " << G4BestUnit(fDetector->GetLength(1),"Length")<< " or " << (material2->GetTotNbOfAtomsPerVolume()/(1/cm3)*fDetector->GetLength(1)/cm)/(1e+15) << " [TFU] "<< G4endl;
  G4cout << " Layer dimensions : " << G4BestUnit(fDetector->GetSize(1), "Length") << G4endl;	
  G4cout << " Material: " << fDetector->GetMaterial(1) << " or name: " << material2->GetName() << G4endl;
  G4cout << " Density: " << density2/(g/cm3) << " g/cm3 " << ", or " << material2->GetTotNbOfAtomsPerVolume()/(1/cm3)  << " [at3] "  << G4endl;
  const G4double *atomDensVector	= material2->GetVecNbOfAtomsPerVolume();
  G4cout << " Composition: " << G4endl;
  for (int i=0; i<NoOfElements_material2; i++){
  G4cout << " " << material2->GetElement(i)->GetName() << " = " << atomDensVector[i]/(1/cm3) << " [cm-3] " << G4endl;}
  G4cout << " Position of the layer: " << G4BestUnit(fDetector->GetPosition(0), "Length") << G4endl;	
  G4cout << " Distance from the surface: " << G4BestUnit(atstums1, "Length") << G4endl;	
  G4cout << " -------------------- Particle parameters -------------------- " << G4endl;
  G4cout << " Ionizing energy loss : \t\t" << G4BestUnit(abs1fEdep, "Energy") << " +/- " << G4BestUnit(abs1rmsEdep, "Energy")  << G4endl;	
  G4cout << " Non ionizing energy loss [NIEL] : \t" << G4BestUnit(abs1niel, "Energy") << " +/- " << G4BestUnit(abs1rmsEnondep, "Energy")  << G4endl;	
  G4cout << " Mean number of steps : \t\t" << abs1step << " +/- " << abs1rmsSteps << G4endl;
  G4cout << " Primary track length : \t\t" << G4BestUnit(abs1len, "Length") << " +/- " << G4BestUnit(abs2rmsTLPrim,"Length") << G4endl; 
  G4cout << " Mean dE/dx : \t\t\t\t" << dEdx_abs1/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean NIEL dE/dx: \t\t\t" << dEdx_niel_abs1/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean free path: \t\t\t" <<G4BestUnit(abs1len/abs1step,"Length")<<  G4endl;
  G4cout << " Steps per length [nm-1]: \t\t" << abs1step/(fDetector->GetLength(1)/nm) << G4endl;
  G4cout << " # of secondaries: \t\t\t" << abs1rec << " \t per primary : " << abs1rec/TotNbofEvents << G4endl;
  G4cout << " Mean E of secondaries: \t\t" << G4BestUnit(abs1sum_t, "Energy") << " +/- " << G4BestUnit(abs1rmssum_t,"Energy") << G4endl;  
  G4cout << " Mean damage en of secondaries: \t" << G4BestUnit(abs1sum_tl, "Energy") << " +/- " << G4BestUnit(abs1rmssum_tl,"Energy") << G4endl;    
  G4cout << " \n";  
  
  G4double NoOfElements_material3    	= material3->GetNumberOfElements();  
  G4cout << " >>>>>>>>>>>>>>>>> SECOND LAYER <<<<<<<<<<<<<<<<<<<<<<<<\n";

  G4cout << " Thickness = " << G4BestUnit(fDetector->GetLength(2),"Length")<< " or " << (material3->GetTotNbOfAtomsPerVolume()/(1/cm3)*fDetector->GetLength(2)/cm)/(1e+15) << " [TFU] "<< G4endl;
  G4cout << " Layer dimensions : " << G4BestUnit(fDetector->GetSize(2), "Length") << G4endl;	
  G4cout << " Material: " << fDetector->GetMaterial(2) << G4endl;
  G4cout << " Density: " << density3/(g/cm3) << " g/cm3 " << ", or " << material3->GetTotNbOfAtomsPerVolume()/(1/cm3)  << " [at3] "  << G4endl;
  const G4double *atomDensVector3	= material3->GetVecNbOfAtomsPerVolume();
  G4cout << " Composition: " << G4endl;
  for (int i=0; i<NoOfElements_material3; i++)
  	{
  G4cout << " " << material3->GetElement(i)->GetName() << " = " << atomDensVector3[i]/(1/cm3) << " [cm-3] " << G4endl;
  	}
  G4cout << " Position of the layer: " << G4BestUnit(fDetector->GetPosition(1), "Length") << G4endl;	
  G4cout << " Distance from the surface: " << G4BestUnit(atstums2, "Length") << G4endl;	
  G4cout << " -------------------- Particle parameters -------------------- " << G4endl;
  G4cout << " Ionizing energy loss : \t\t" << G4BestUnit(abs2fEdep, "Energy") << " +/- " << G4BestUnit(abs2rmsEdep, "Energy")  << G4endl;	
  G4cout << " Non ionizing energy loss [NIEL] : \t" << G4BestUnit(abs2niel, "Energy") << " +/- " << G4BestUnit(abs2rmsEnondep, "Energy")  << G4endl;	
  G4cout << " Mean number of steps : \t\t" << abs2step << " +/- " << abs2rmsSteps << G4endl;
  G4cout << " Primary track length : \t\t" << G4BestUnit(abs2len, "Length") << " +/- " << G4BestUnit(abs2rmsTLPrim,"Length") << G4endl; 
  G4cout << " Mean dE/dx : \t\t\t\t" << dEdx_abs2/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean NIEL dE/dx: \t\t\t" << dEdx_niel_abs2/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean free path: \t\t\t" <<G4BestUnit(abs2len/abs2step,"Length")<<  G4endl;
  G4cout << " Steps per length [nm-1]: \t\t" << abs2step/(fDetector->GetLength(2)/nm) << G4endl;  
  G4cout << " # of secondaries: \t\t\t" << abs2rec << " \t per primary : " << abs2rec/TotNbofEvents << G4endl;
  G4cout << " Mean E of secondaries: \t\t" << G4BestUnit(abs2sum_t, "Energy") << " +/- " << G4BestUnit(abs2rmssum_t,"Energy") << G4endl;  
  G4cout << " Mean damage en of secondaries: \t" << G4BestUnit(abs2sum_tl, "Energy") << " +/- " << G4BestUnit(abs2rmssum_tl,"Energy") << G4endl;    
  G4cout << " \n";    

  G4double NoOfElements_material4    	= material4->GetNumberOfElements();  
  G4cout << " >>>>>>>>>>>>>>>>> THIRD LAYER <<<<<<<<<<<<<<<<<<<<<<<<\n";

  G4cout << " Thickness = " << G4BestUnit(fDetector->GetLength(3),"Length")<< " or " << (material4->GetTotNbOfAtomsPerVolume()/(1/cm3)*fDetector->GetLength(3)/cm)/(1e+15) << " [TFU] "<< G4endl;
  G4cout << " Layer dimensions : " << G4BestUnit(fDetector->GetSize(3), "Length") << G4endl;	
  G4cout << " Material: " << fDetector->GetMaterial(3) << G4endl;
  G4cout << " Density: " << density4/(g/cm3) << " g/cm3 " << ", or " << material4->GetTotNbOfAtomsPerVolume()/(1/cm3)  << " [at3] "  << G4endl;
  const G4double *atomDensVector4	= material4->GetVecNbOfAtomsPerVolume();
  G4cout << " Composition: " << G4endl;
  for (int i=0; i<NoOfElements_material4; i++)
  	{
  G4cout << " " << material4->GetElement(i)->GetName() << " = " << atomDensVector4[i]/(1/cm3) << " [cm-3] " << G4endl;
  	}
  G4cout << " Position of the layer: " << G4BestUnit(fDetector->GetPosition(2), "Length") << G4endl;	
  G4cout << " Distance from the surface: " << G4BestUnit(atstums3, "Length") << G4endl;	
  G4cout << " -------------------- Particle parameters -------------------- " << G4endl;
  G4cout << " Ionizing energy loss : \t\t" << G4BestUnit(abs3fEdep, "Energy") << " +/- " << G4BestUnit(abs3rmsEdep, "Energy")  << G4endl;	
  G4cout << " Non ionizing energy loss [NIEL] : \t" << G4BestUnit(abs3niel, "Energy") << " +/- " << G4BestUnit(abs3rmsEnondep, "Energy")  << G4endl;	
  G4cout << " Mean number of steps : \t\t" << abs3step << " +/- " << abs2rmsSteps << G4endl;
  G4cout << " Primary track length : \t\t" << G4BestUnit(abs3len, "Length") << " +/- " << G4BestUnit(abs3rmsTLPrim,"Length") << G4endl; 
  G4cout << " Mean dE/dx : \t\t\t\t" << dEdx_abs3/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean NIEL dE/dx: \t\t\t" << dEdx_niel_abs3/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean free path: \t\t\t" <<G4BestUnit(abs3len/abs3step,"Length")<<  G4endl;
  G4cout << " Steps per length [nm-1]: \t\t" << abs3step/(fDetector->GetLength(3)/nm) << G4endl;    
  G4cout << " # of secondaries: \t\t\t" << abs3rec << " \t per primary : " << abs3rec/TotNbofEvents << G4endl;
  G4cout << " Mean E of secondaries: \t\t" << G4BestUnit(abs3sum_t, "Energy") << " +/- " << G4BestUnit(abs3rmssum_t,"Energy") << G4endl;  
  G4cout << " Mean damage en of secondaries: \t" << G4BestUnit(abs3sum_tl, "Energy") << " +/- " << G4BestUnit(abs3rmssum_tl,"Energy") << G4endl;    
  G4cout << " \n";  
  
  G4double NoOfElements_material5    	= material5->GetNumberOfElements();  
  G4cout << " >>>>>>>>>>>>>>>>> FOURTH LAYER <<<<<<<<<<<<<<<<<<<<<<<<\n";

  G4cout << " Thickness = " << G4BestUnit(fDetector->GetLength(4),"Length")<< " or " << (material5->GetTotNbOfAtomsPerVolume()/(1/cm3)*fDetector->GetLength(4)/cm)/(1e+15) << " [TFU] "<< G4endl;
  G4cout << " Layer dimensions : " << G4BestUnit(fDetector->GetSize(4), "Length") << G4endl;	
  G4cout << " Material: " << fDetector->GetMaterial(4) << G4endl;
  G4cout << " Density: " << density5/(g/cm3) << " g/cm3 " << ", or " << material5->GetTotNbOfAtomsPerVolume()/(1/cm3)  << " [at3] "  << G4endl;
  const G4double *atomDensVector5	= material4->GetVecNbOfAtomsPerVolume();
  G4cout << " Composition: " << G4endl;
  for (int i=0; i<NoOfElements_material5; i++)
  	{
  G4cout << " " << material5->GetElement(i)->GetName() << " = " << atomDensVector5[i]/(1/cm3) << " [cm-3] " << G4endl;
  	}
  G4cout << " Position of the layer: " << G4BestUnit(fDetector->GetPosition(3), "Length") << G4endl;	
  G4cout << " Distance from the surface: " << G4BestUnit(atstums4, "Length") << G4endl;	
  G4cout << " -------------------- Particle parameters -------------------- " << G4endl;
  G4cout << " Ionizing energy loss : \t\t" << G4BestUnit(abs4fEdep, "Energy") << " +/- " << G4BestUnit(abs4rmsEdep, "Energy")  << G4endl;	
  G4cout << " Non ionizing energy loss [NIEL] : \t" << G4BestUnit(abs4niel, "Energy") << " +/- " << G4BestUnit(abs4rmsEnondep, "Energy")  << G4endl;	
  G4cout << " Mean number of steps : \t\t" << abs4step << " +/- " << abs4rmsSteps << G4endl;
  G4cout << " Primary track length : \t\t" << G4BestUnit(abs4len, "Length") << " +/- " << G4BestUnit(abs4rmsTLPrim,"Length") << G4endl; 
  G4cout << " Mean dE/dx : \t\t\t\t" << dEdx_abs4/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean NIEL dE/dx: \t\t\t" << dEdx_niel_abs4/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean free path: \t\t\t" <<G4BestUnit(abs4len/abs4step,"Length")<<  G4endl;
  G4cout << " Steps per length [nm-1]: \t\t" << abs4step/(fDetector->GetLength(4)/nm) << G4endl;    
  G4cout << " # of secondaries: \t\t\t" << abs4rec << " \t per primary : " << abs4rec/TotNbofEvents << G4endl;
  G4cout << " Mean E of secondaries: \t\t" << G4BestUnit(abs4sum_t, "Energy") << " +/- " << G4BestUnit(abs4rmssum_t,"Energy") << G4endl;  
  G4cout << " Mean damage en of secondaries: \t" << G4BestUnit(abs4sum_tl, "Energy") << " +/- " << G4BestUnit(abs4rmssum_tl,"Energy") << G4endl;    
  G4cout << " \n";       

   G4cout << " **************************************************************\n"; 
  G4cout << " Max Step size : \t\t\t" <<G4BestUnit(fDetector->GetMaxStep(),"Length")<<  G4endl; 
	rbs_angle = fDetector->GetRBSAngle()/degree;
  G4cout << " RBS detector angle: \t\t\t" << rbs_angle << " degrees " << G4endl;
  G4cout << " Average Max RBS depth: \t\t" << G4BestUnit(RBSDepth/(counts), "Length") << " +/- " << G4BestUnit(rmsRBSDepth/counts, "Length") << G4endl;
  G4cout << " Use of sigmacalc: \t\t\t" ;
  if (fDetector->GetSigmaCalc() == 1) { G4cout << "ENABLED " << G4endl; };
  if (fDetector->GetSigmaCalc() == 0) { G4cout << "DISABLED " << G4endl; };
  G4cout << " RBS evaluation was: \t\t\t" ;
  if (fDetector->GetRBSCalc() == 1) { G4cout << "ENABLED " << G4endl; };
  if (fDetector->GetRBSCalc() == 0) { G4cout << "DISABLED " << G4endl; };
  G4cout << " MS evaluation was: \t\t\t" ;
  if (fDetector->GetMSCalc() == 1) { G4cout << "ENABLED " << G4endl; };
  if (fDetector->GetMSCalc() == 0) { G4cout << "DISABLED " << G4endl; };  
  G4cout << " MS stopping corrections were:   \t" ;
  if (fDetector->GetUseMSCorrections() == 1) { G4cout << "ENABLED " << G4endl; };
  if (fDetector->GetUseMSCorrections() == 0) { G4cout << "DISABLED " << G4endl; };  
  G4cout << " Constant scattering angle: \t\t" ;
  if (fDetector->GetConstAngle() == 1) { G4cout << "ENABLED " << G4endl; };
  if (fDetector->GetConstAngle() == 0) { G4cout << "DISABLED " << G4endl; };       
  G4cout << " RBS minimum ROI : \t\t\t" << G4BestUnit(fDetector->GetRBSROImin(), "Energy") << G4endl;
  G4cout << " \n";


  G4double chan = hit/4;
  G4double totchan = chan/TotNbofEvents;
  G4double totem = chan/partEmerging;
  // kinetic energy of particles reaching the 4th detector
  G4double rmssum_kinen = 0.0;
  if(chan > 0) {
    detKinEn/=hit;     detKinEn2/=hit;
    rmssum_kinen =detKinEn2- detKinEn*detKinEn;
    if(rmssum_kinen>0.) rmssum_kinen=std::sqrt(rmssum_kinen/chan);  }

  G4cout << " **************************************************** " << G4endl;
  G4cout 
    << " Projected range = " << G4BestUnit(projectedR,"Length")
    << " +- " << G4BestUnit( rmsPR,"Length")<< G4endl;
  G4cout << " Rotation of the sample: " << fDetector->GetAngles()/degree << " degrees " << G4endl;
    
  G4cout << " **************************************************** " << G4endl;
  G4cout << " # of particles that reaches last detector " << hit/4 << G4endl;
  G4cout << " part of total particles " << std::fixed << std::setprecision(2) << totchan*100 << " % " <<G4endl;
  G4cout << " part of emitted particles " << std::setprecision(2) << totem*100 << " % "<< G4endl;

  G4cout << " Mean KinEn of particles reaching 4th detector = "<<G4BestUnit(detKinEn,"Energy")
                <<" +/- "  <<G4BestUnit(rmssum_kinen,"Energy")<<G4endl; 
  G4cout << " ****************************************************** " << G4endl;
  G4double no_of_steps_per_particle = entry_sd/numberOfEvent;
  G4double no_of_reach_per_particle = entry_reach/numberOfEvent;
  G4cout << " number of steps : " << no_of_steps_per_particle << G4endl; 
  G4cout << " number of entries reached histo : " << no_of_reach_per_particle << G4endl;
  
  G4cout << " average number of entries per step : " << no_of_reach_per_particle/no_of_steps_per_particle << G4endl;
  
  G4double total_step_length =  total_step/numberOfEvent;
  G4cout << " total step : " << G4BestUnit(total_step_length, "Length") << G4endl;
  G4cout << " average step length: " << G4BestUnit(total_step_length/no_of_steps_per_particle, "Length") << G4endl;
  G4cout << " ****************************************************** " << G4endl;
  G4cout << " ******************** GEOMETRY ************************ " << G4endl;
  G4cout << " ****************************************************** " << G4endl;  




  G4cout << " Angle of incidence: " << (angle_of_incidence/numberOfEvent)/degree << " [degree] " << G4endl;
  G4cout << " Primary particles " << Particle << " with energy: " << G4BestUnit(primary_energy/TotNbofEvents, "Energy") << G4endl;

  G4double plotis[4] = {fDetector->GetLength(1), fDetector->GetLength(2), fDetector->GetLength(3), fDetector->GetLength(4)};
  
  
  G4ThreeVector pozicija_visu[4] = {fDetector->GetPosition(0),fDetector->GetPosition(1),fDetector->GetPosition(2),fDetector->GetPosition(3)};
  G4double ilgis_visu[4] = {pozicija_visu[0].z(), pozicija_visu[1].z(), pozicija_visu[2].z(), pozicija_visu[3].z()};
   G4double atstumas_tarp_gretimu[3];
   for (int i=0; i<3; i++)
   	{
   		atstumas_tarp_gretimu[i] = (ilgis_visu[i+1]-plotis[i+1]/2)-(ilgis_visu[i]+plotis[i]/2);
   	}

  
  G4cout << " Materials: \n " << G4endl;
  G4cout << " X = " << material->GetName() << G4endl;
  G4cout << " A = " << material2->GetName() << G4endl;
  G4cout << " B = " << material3->GetName() << G4endl;
  G4cout << " C = " << material4->GetName() << G4endl;
  G4cout << " D = " << material5->GetName() << G4endl;
  G4cout << " \n" << G4endl;

  G4double final_dist = fDetector->GetLength(0)/2 - ilgis_visu[3] - plotis[3]/2;
	//G4cout << " test dist " << (fDetector->GetLength(0)/2 - ilgis_visu[3] - plotis[3]/2)/nm << G4endl;

  G4double a = 1e+15;
  //G4double TFU0 = (material->GetTotNbOfAtomsPerVolume()/(1/cm3)*fDetector->GetLength(0)/cm)/a;
  G4double TFU1 = (material2->GetTotNbOfAtomsPerVolume()/(1/cm3)*fDetector->GetLength(1)/cm)/a;
  G4double TFU2 = (material3->GetTotNbOfAtomsPerVolume()/(1/cm3)*fDetector->GetLength(2)/cm)/a;
  G4double TFU3 = (material4->GetTotNbOfAtomsPerVolume()/(1/cm3)*fDetector->GetLength(3)/cm)/a;
  G4double TFU4 = (material5->GetTotNbOfAtomsPerVolume()/(1/cm3)*fDetector->GetLength(4)/cm)/a;  
  G4double TFU5 = (material->GetTotNbOfAtomsPerVolume()/(1/cm3)*final_dist/cm)/a; 


  G4double TFU01 = (material->GetTotNbOfAtomsPerVolume()/(1/cm3)*atstums1/cm)/a; 
  G4double TFU12 = (material->GetTotNbOfAtomsPerVolume()/(1/cm3)*atstumas_tarp_gretimu[0]/cm)/a;
  G4double TFU23 = (material->GetTotNbOfAtomsPerVolume()/(1/cm3)*atstumas_tarp_gretimu[1]/cm)/a;
  G4double TFU34 = (material->GetTotNbOfAtomsPerVolume()/(1/cm3)*atstumas_tarp_gretimu[2]/cm)/a;

  G4cout << " TFU UNITS " << G4endl;
  G4cout << " " << TFU01 << " || " << TFU1 << " || " << TFU12 << " || " << TFU2 << " || " << TFU23 << " || " << TFU3 << " || " << TFU34 << " || " << TFU4 << " || " << TFU5 << G4endl;
  
    G4cout << " " << G4BestUnit(atstums1, "Length") <<  "  " << G4BestUnit(plotis[0], "Length") << "  "<< G4BestUnit(atstumas_tarp_gretimu[0], "Length")  <<  "  " << G4BestUnit(plotis[1], "Length") << "  " << G4BestUnit(atstumas_tarp_gretimu[1], "Length") <<  "  " << G4BestUnit(plotis[2], "Length") << "  "<<G4BestUnit(atstumas_tarp_gretimu[2], "Length") <<  "  " << G4BestUnit(plotis[3], "Length") << "  " << G4BestUnit(final_dist, "Length") <<G4endl;
  
  
  G4cout << " XXXXXX->|<-AAAAAA->|<-XXXXXX->|<-BBBBBB->|<-XXXXXX->|<-CCCCCC->|<-XXXXXX->|<-DDDDDD->|<-XXXXXX" << G4endl;

  G4cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX " << G4endl;
  G4cout << " |<-------------- Substrate " << G4BestUnit(fDetector->GetLength(0), "Length") << " --------------------------------------------------->| " << G4endl;
  G4cout << "                                                        " << G4endl;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
  if (analysisManager->IsActive() ) 
  {      
	G4double fac; 
  	for (G4int ih=3; ih<15; ih++) 
  	{
	
		G4double binWidth = analysisManager->GetH1Width(ih);
    //G4double unit     = analysisManager->GetH1Unit(ih);  

    //fac = unit/binWidth;
    		fac = (1./(numberOfEvent*binWidth))*(mm/MeV);
			analysisManager->ScaleH1(2,fac);
    //fac = unit/binWidth;
    		fac = (1./totstep); 
			analysisManager->ScaleH1(15,fac);
    		fac = (1./(numberOfEvent*binWidth))*(mm/keV);
			analysisManager->ScaleH1(17,fac);
    		fac = (1./numberOfEvent);
			analysisManager->ScaleH1(18,fac);
    		fac = (1./(totstep*binWidth));
			analysisManager->ScaleH1(19,fac);
    	}
    	for (G4int ih=20; ih<53; ih++) 
    	{

    		G4double binW = analysisManager->GetH1Width(ih);
    
    		G4double ave_step = total_step_length/no_of_steps_per_particle;
    		G4double RBS_norm_dist = 0.1*nm;
    		G4double norm = ave_step/RBS_norm_dist;
    		G4double exponent = exp(1);
    		fac = (1/(exponent*norm*no_of_steps_per_particle*numberOfEvent*binW));	// latest, 2021-02-23
    		analysisManager->ScaleH1(ih,fac);
	} 
		G4double BW1 = analysisManager->GetH1Width(16);
        	G4double fcc = (1./(BW1*numberOfEvent));
		analysisManager->ScaleH1(16,fcc);


}
  //remove all contents in fProcCounter, fCount 
  fProcCounter.clear();
  fParticleDataMap2.clear();
                          
  //restore default format         
  G4cout.precision(dfprec);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
