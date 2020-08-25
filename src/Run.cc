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

//#include "G4MultiFunctionalDetector.hh"
//#include "G4SDManager.hh"
//#include "G4VPrimitiveScorer.hh"

//keisti neutronu srautui
//#include "RunMessenger.hh"

//#include <fstream>
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

	// other material

   abs2len=abs2len2=0.;
   abs2fEdep=abs2fEdep2=0.;
   abs2niel=abs2niel2=0.;
   abs2step=abs2step2=0.;
   abs2sum_tl=abs2sum_tl2=0.;
   abs2sum_t=abs2sum_t2=0.;
   abs2rec=0.;

   abs3niel=abs3niel2=0;



  abs3fEdep=abs3fEdep2=0.;
  abs3step =abs3step2=0.;
  abs3len  = abs3len2=0.;



   totniel=totniel2=0.;


   projectedR=projectedR2=0.;

	entry_sd 	= 0.;
	total_step 	= 0.;
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
  // absorber 
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
  // end of absorber


  // absorber2 
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
  //mean number of steps:
  abs3step/=TotNbofEvents;  abs3step2/=TotNbofEvents;
  G4double abs3rmsSteps= abs3step2 -abs3step*abs3step;
  if(abs3rmsSteps>0) abs3rmsSteps= std::sqrt(abs3rmsSteps/TotNbofEvents);
  else abs3rmsSteps=0;
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



//totniel
  // PLEASE END MY MYSERY
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
  G4cout << " ================= GEOMETRIJOS parametrai==================";
  G4cout << "\n ===========================================================\n";
  G4cout << " \n"; 
  G4cout << " >>>>>>>>>>>>>>>>> Pirmasis sluoksnis <<<<<<<<<<<<<<<<<<<<<<<<\n";
  G4cout << " \n"; 
  G4cout << " Sluoksnio storis = " << G4BestUnit(fDetector->GetLength(0),"Length")<<G4endl;
  G4cout << " Bendri absorberio matmenys : " << G4BestUnit(fDetector->GetSize(0), "Length") << G4endl;	
  G4cout << " Medziaga: " << fDetector->GetMaterial(0) << G4endl;
  G4cout << " Tankis: " << density/(g/cm3) << " g/cm3 " << G4endl;
  G4cout << " -------------------- Dalelių parametrai -------------------- " << G4endl;
  G4cout << " Ionizing energy loss : \t\t" << G4BestUnit(abs2fEdep, "Energy") << " +/- " << G4BestUnit(abs2rmsEdep, "Energy")  << G4endl;	
  G4cout << " Non ionizing energy loss [NIEL] : \t" << G4BestUnit(abs2niel, "Energy") << " +/- " << G4BestUnit(abs2rmsEnondep, "Energy")  << G4endl;	
  G4cout << " Mean number of steps : \t\t" << abs2step << " +/- " << abs2rmsSteps << G4endl;
  G4cout << " Primary track length : \t\t" << G4BestUnit(abs2len, "Length") << " +/- " << G4BestUnit(abs2rmsTLPrim,"Length") << G4endl; 
  G4cout << " Mean dE/dx : \t\t\t\t" << dEdx_abs2/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean NIEL dE/dx: \t\t\t" << dEdx_niel_abs2/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean free path: \t\t\t" <<G4BestUnit(abs2len/abs2step,"Length")<<  G4endl;
  G4cout << " Steps per length [nm-1]: \t\t" << abs2step/(fDetector->GetLength(0)/nm) << G4endl;
  G4cout << " \n"; 
  G4cout << " >>>>>>>>>>>>>>>>> Antrasis sluoksnis <<<<<<<<<<<<<<<<<<<<<<<<\n";
  G4cout << " \n"; 
  G4cout << " Antrojo sluoksnio storis = " << G4BestUnit(fDetector->GetLength(1),"Length")<<G4endl;
  G4cout << " Bendri antrojo absorberio matmenys : " << G4BestUnit(fDetector->GetSize(1), "Length") << G4endl;	
  G4cout << " Medziaga: " << fDetector->GetMaterial(1) << G4endl;
  G4cout << " Tankis: " << density2/(g/cm3) << " g/cm3 " << G4endl;
  G4cout << " Antrojo sluoksnio pozicija: " << G4BestUnit(fDetector->GetPosition(0), "Length") << G4endl;	
  G4ThreeVector pozicija = fDetector->GetPosition(0);
  
  G4double ilgis = pozicija.z();
  //G4cout << " pozicija " << ilgis/um << G4endl; //+(ilgis)
  //G4cout << " det ilgis " << (fDetector->GetLength()/2)/um << G4endl;
  G4double atstums = (fDetector->GetLength(0)/2)+ilgis-(fDetector->GetLength(1)/2);

  
  
  G4cout << " Atstumas nuo pavirsiaus: " << G4BestUnit(atstums, "Length") << G4endl;	
  G4cout << " -------------------- Dalelių parametrai -------------------- " << G4endl;
  G4cout << " Ionizing energy loss : \t\t" << G4BestUnit(absfEdep, "Energy") << " +/- " << G4BestUnit(absrmsEdep, "Energy")  << G4endl;	
  G4cout << " Non ionizing energy loss [NIEL] : \t" << G4BestUnit(absniel, "Energy") << " +/- " << G4BestUnit(absrmsEnondep, "Energy")  << G4endl;	
  G4cout << " Mean number of steps : \t\t" << absstep << " +/- " << absrmsSteps << G4endl;
  G4cout << " Primary track length : \t\t" << G4BestUnit(abslen, "Length") << " +/- " << G4BestUnit(absrmsTLPrim,"Length") << G4endl; 
  G4cout << " Mean dE/dx : \t\t\t\t" << dEdx_abs/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean NIEL dE/dx: \t\t\t" << dEdx_niel_abs/(MeV/cm) << " MeV/cm " << G4endl;
  G4cout << " Mean free path: \t\t\t" <<G4BestUnit(abslen/absstep,"Length")<<  G4endl;
  G4cout << " Steps per length [nm-1]: \t\t" << absstep/(fDetector->GetLength(1)/nm) << G4endl;
  G4cout << " **************************************************************\n";
  G4cout << " \t\t SUM " << G4endl;
  G4cout << " Ionizing energy loss : \t\t" << G4BestUnit(abs2fEdep+absfEdep, "Energy") << G4endl;	
  G4cout << " Non ionizing energy loss [NIEL] : \t" << G4BestUnit(abs2niel+absniel, "Energy") << G4endl;	
  G4cout << " Mean number of steps : \t\t" << abs2step+absstep << G4endl;
  G4cout << " Primary track length : \t\t" << G4BestUnit(abs2len+abslen, "Length") << G4endl; 
  G4cout << " Mean free path: \t\t\t" <<G4BestUnit((abs2len+abslen)/(abs2step+absstep),"Length")<<  G4endl;
  G4cout << " **************************************************************\n";
  G4cout << " \t\t In Between materials " << G4endl;
  G4cout << " Ionizing energy loss : \t\t" << G4BestUnit(abs3fEdep, "Energy") << " +/- " << G4BestUnit(abs3rmsEdep, "Energy")  << G4endl;	
  G4cout << " Non ionizing energy loss [NIEL] : \t" << G4BestUnit(abs3niel, "Energy") << " +/- " << G4BestUnit(abs3rmsEnondep, "Energy")  << G4endl;	
  G4cout << " Mean number of steps : \t\t" << abs3step << " +/- " << abs3rmsSteps << G4endl;
  G4cout << " Primary track length : \t\t" << G4BestUnit(abs3len, "Length") << " +/- " << G4BestUnit(abs3rmsTLPrim,"Length") << G4endl;
  G4cout << " Mean free path: \t\t\t" <<G4BestUnit(abs3len/abs3step,"Length")<<  G4endl; 
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
  G4cout << " \n";
  G4cout << " **************************************************************\n";
  G4cout << " ***************** Detektoriai ********************************\n";
  G4cout << " Detektoriu matmenys: "<< G4BestUnit(fDetector->GetDetectorSizes(), "Length") << G4endl;
  G4cout << " Detektoriu medziaga: "<< fDetector->GetDetectorMaterial() << G4endl;
  G4cout << " Detektoriu pozicijos: "<< G4BestUnit(fDetector->GetDetectorDistance(0), "Length") << "  " << G4BestUnit(fDetector->GetDetectorDistance(1), "Length") << "  " << G4BestUnit(fDetector->GetDetectorDistance(2), "Length") << "  " << G4BestUnit(fDetector->GetDetectorDistance(3), "Length") << G4endl;



  G4double chan = hit/4;
  G4double totchan = chan/TotNbofEvents;
  G4double totem = chan/partEmerging;
  // kinetic energy of particles reaching the 4th detector
  G4double rmssum_kinen = 0.0;
  if(chan > 0) {
    detKinEn/=hit;     detKinEn2/=hit;
    rmssum_kinen =detKinEn2- detKinEn*detKinEn;
    if(rmssum_kinen>0.) rmssum_kinen=std::sqrt(rmssum_kinen/chan);  }
  //


  G4cout << " **************************************************** " << G4endl;
  G4cout 
    << " Projected range = " << G4BestUnit(projectedR,"Length")
    << " +- " << G4BestUnit( rmsPR,"Length")<< G4endl;
  G4cout << " Kristalo pasukimo kampas x y ir z kryptimis: " << fDetector->GetAngles()/degree << " degrees " << G4endl;
    
  G4cout << " **************************************************** " << G4endl;
  G4cout << " # of particles that reaches last detector " << hit/4 << G4endl;
  G4cout << " part of total particles " << std::fixed << std::setprecision(2) << totchan*100 << " % " <<G4endl;
  G4cout << " part of emitted particles " << std::setprecision(2) << totem*100 << " % "<< G4endl;

  G4cout << " Mean KinEn of particles reaching 4th detector = "<<G4BestUnit(detKinEn,"Energy")
                <<" +/- "  <<G4BestUnit(rmssum_kinen,"Energy")<<G4endl; 
  G4cout << " ****************************************************** " << G4endl;
  G4cout << " # of backscattered primaries " << rbs << " see backscattered.txt for more info. " << G4endl;
  G4cout << " # of secondaries generated " << second << " see secondaries.txt for more info. " << G4endl;
  G4cout << " ****************************************************** " << G4endl;
  G4cout << " ****************************************************** " << G4endl;
  G4cout << " DON'T FORGET TO FLUSH secondaries AND backscattered FILES " << G4endl;
  G4cout << " ****************************************************** " << G4endl;

  G4cout << " ****************************************************** " << G4endl;

  G4cout << " number of entries : " << entry_sd << G4endl; 
  //G4double sum_entries = entry_sdx + entry_sd;
  //G4cout << " sum of entries " << sum_entries << " per particle " << sum_entries/numberOfEvent << G4endl;
  //G4cout << " entry density, 1 : " << entry_sdx/(abslen/um) << " and 2 : " << entry_sd/(abs2len/um) << G4endl;
  //........................................................
  G4cout << " ****************************************************** " << G4endl;
  G4cout << " ************ Eksperimento geometrija ***************** " << G4endl;
  G4cout << " ****************************************************** " << G4endl;  
  
  G4double plotis[4] = {fDetector->GetLength(1), fDetector->GetLength(2), fDetector->GetLength(3), fDetector->GetLength(4)};
  
  
  G4ThreeVector pozicija_visu[4] = {fDetector->GetPosition(0),fDetector->GetPosition(1),fDetector->GetPosition(2),fDetector->GetPosition(3)};
  G4double ilgis_visu[4] = {pozicija_visu[0].z(), pozicija_visu[1].z(), pozicija_visu[2].z(), pozicija_visu[3].z()};
  
  G4double atstumas_visu[4] = {0}; 
  for (int i=0; i<4; i++)
  	{
  		atstumas_visu[i] = -(fDetector->GetLength(0)/2)-ilgis_visu[i]+plotis[i]/2; 
  		//G4cout << " atstumas " << i << " " << abs(atstumas_visu[i])/um << G4endl;
  	}
  	
   G4double atstumas_tarp_gretimu[3];
   for (int i=0; i<3; i++)
   	{
   		atstumas_tarp_gretimu[i] = (ilgis_visu[i+1]-plotis[i+1]/2)-(ilgis_visu[i]+plotis[i]/2);

   		//G4cout << " atstumas " << i << " " << atstumas_tarp_gretimu[i]/um << G4endl;	
   		
   	}

  
  G4cout << " Materials: \n " << G4endl;
  G4cout << " X = " << fDetector->GetMaterialM(0)->GetName() << G4endl;
  G4cout << " A = " << fDetector->GetMaterialM(1)->GetName() << G4endl;
  G4cout << " B = " << fDetector->GetMaterialM(2)->GetName() << G4endl;
  G4cout << " C = " << fDetector->GetMaterialM(3)->GetName() << G4endl;
  G4cout << " D = " << fDetector->GetMaterialM(4)->GetName() << G4endl;
  G4cout << " \n" << G4endl;
  
    G4cout << " " << G4BestUnit(atstums, "Length") <<  " " << G4BestUnit(plotis[0], "Length") << " "<< G4BestUnit(atstumas_tarp_gretimu[0], "Length")  <<  "  " << G4BestUnit(plotis[1], "Length") << "  " << G4BestUnit(atstumas_tarp_gretimu[1], "Length") <<  "  " << G4BestUnit(plotis[2], "Length") << "  "<<G4BestUnit(atstumas_tarp_gretimu[2], "Length") <<  "  " << G4BestUnit(plotis[3], "Length") << "  " << G4endl;
  
  
  
  G4cout << " XXXXX->|<-AAAAA->|<-XXXXX->|<-BBBBB->|<-XXXXX->|<-CCCCC->|<-XXXXX->|<-DDDD->|<-XXXX" << G4endl;

  G4cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX " << G4endl;
  G4cout << " |<-------------- Substrate " << G4BestUnit(fDetector->GetLength(0), "Length") << " --------------------------------------------->| " << G4endl;
  G4cout << "                                                        " << G4endl;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
  if (analysisManager->IsActive() ) {      



	G4double total_step_length = total_step/numberOfEvent;


	//analysisManager->OpenFile();
  //normalize histograms
  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();       
	G4double fac; 
  for (G4int ih=3; ih<15; ih++) {

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
    for (G4int ih=20; ih<40; ih++) {
    //fac = (1./numberOfEvent);
    fac = (1./(numberOfEvent*total_step_length));
    analysisManager->ScaleH1(ih,fac);
	} 
	G4double BW1 = analysisManager->GetH1Width(16);
        G4double fcc = (1./(BW1*numberOfEvent));
	analysisManager->ScaleH1(16,fcc);


    //G4double BW = analysisManager->GetH1Width(39);  
    G4double sme = (mm/MeV)/(numberOfEvent);// * BW);
    analysisManager->ScaleH1(39, sme);

}
  //remove all contents in fProcCounter, fCount 
  fProcCounter.clear();
  fParticleDataMap2.clear();
                          
  //restore default format         
  G4cout.precision(dfprec);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
