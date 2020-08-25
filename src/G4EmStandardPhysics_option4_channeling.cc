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
//---------------------------------------------------------------------------
//
// ClassName:   G4EmStandardPhysics_option4_channeling
//
// Author:      V.Ivanchenko 28.09.2012 from Option3 physics constructor
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4EmStandardPhysics_option4_channeling.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4PEEffectFluoModel.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LowEPComptonModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4MscStepLimitType.hh"
#include "G4UrbanMscModel.hh"
#include "G4DummyModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4Generator2BS.hh"
#include "G4Generator2BN.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4eplusAnnihilation.hh"
#include "G4UAtomicDeexcitation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"
#include "G4ePairProduction.hh"

#include "G4MuBremsstrahlungModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4hBremsstrahlungModel.hh"
#include "G4hPairProductionModel.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmModelActivator.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"



// jomajo

#include "G4StepLimiter.hh"
#include "G4ProcessManager.hh"

#include "G4IonCoulombScatteringModel.hh"


// Wrapper

#include "G4HadronElasticProcess.hh"
#include "G4HadronElastic.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4ChipsProtonElasticXS.hh"
#include "G4ChipsElasticModel.hh"
#include "G4IonPhysics.hh"

#include "G4HadronInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4CrossSectionInelastic.hh"

#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4FTFBuilder.hh"
#include "G4HadronicInteraction.hh"
#include "G4BuilderType.hh"
#include "G4ProtonInelasticProcess.hh"

//
/*
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmStandardPhysics_option4_channeling);
	// new addition
G4ThreadLocal G4BinaryLightIonReaction* G4EmStandardPhysics_option4_channeling::theIonBC = 0;
G4ThreadLocal G4HadronicInteraction*  G4EmStandardPhysics_option4_channeling::theFTFP = 0;
G4ThreadLocal G4VCrossSectionDataSet*  G4EmStandardPhysics_option4_channeling::theNuclNuclData = 0;
G4ThreadLocal G4VComponentCrossSection* G4EmStandardPhysics_option4_channeling::theGGNuclNuclXS = 0;
G4ThreadLocal G4FTFBuilder* G4EmStandardPhysics_option4_channeling::theBuilder;
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmStandardPhysics_option4_channeling::G4EmStandardPhysics_option4_channeling(G4int ver,
                                                                               const G4String&)
  : G4VPhysicsConstructor("G4EmStandard_opt4_channeling"), verbose(ver)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetVerbose(verbose);
  param->SetMinEnergy(100*eV);// orgas
  param->SetMaxEnergy(10*MeV);
  param->SetLowestElectronEnergy(10*eV);
  //param->SetNumberOfBins(40000);//pridetas
  param->SetNumberOfBinsPerDecade(200);
  param->ActivateAngularGeneratorForIonisation(true);
  param->SetMscThetaLimit(0.0);
  param->SetFluo(true);
  param->SetAuger(true);
  param->SetPixe(true);
  param->SetBuildCSDARange(true); // nebuvo, mano pridetas 
  param->SetLossFluctuations(true);
  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmStandardPhysics_option4_channeling::~G4EmStandardPhysics_option4_channeling()
{    /*
    delete theBuilder; theBuilder = 0;
    delete theGGNuclNuclXS; theGGNuclNuclXS = 0;
    delete theNuclNuclData; theNuclNuclData = 0;
    delete theIonBC; theIonBC = 0;
    delete theFTFP; theFTFP = 0;*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmStandardPhysics_option4_channeling::ConstructParticle()
{
  // gamma
  G4Gamma::Gamma();

  // leptons
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();

  // mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();

  // barions
  G4Proton::Proton();
  G4AntiProton::AntiProton();

  // ions
  G4Deuteron::Deuteron();
  G4Triton::Triton();
  G4He3::He3();
  G4Alpha::Alpha();
  G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmStandardPhysics_option4_channeling::ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // muon & hadron bremsstrahlung and pair production
  G4MuBremsstrahlung* mub = new G4MuBremsstrahlung();
  G4MuPairProduction* mup = new G4MuPairProduction();
  G4hBremsstrahlung* pib = new G4hBremsstrahlung();
  G4hPairProduction* pip = new G4hPairProduction();
  G4hBremsstrahlung* kb = new G4hBremsstrahlung();
  G4hPairProduction* kp = new G4hPairProduction();
  //G4hBremsstrahlung* pb = new G4hBremsstrahlung();
  //G4hPairProduction* pp = new G4hPairProduction();
  G4ePairProduction* ee = new G4ePairProduction();

  // muon & hadron multiple scattering
  G4CoulombScattering* muss = new G4CoulombScattering();
  G4CoulombScattering* piss = new G4CoulombScattering();
  G4CoulombScattering* kss = new G4CoulombScattering();
  //G4CoulombScattering* pss = new G4CoulombScattering();

  // energy limits for e+- scattering models
  //G4double highEnergyLimit = 100*MeV;
  // energy limits for e+- ionisation models
  G4double penEnergyLimit = 1*MeV;

  // nuclear stopping
  G4NuclearStopping* pnuc = new G4NuclearStopping();
  


  // Add standard EM Processes
  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() ){
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4String particleName = particle->GetParticleName();
      	    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (particleName == "gamma") {

      // Photoelectric
      G4PhotoElectricEffect* pe = new G4PhotoElectricEffect();
      G4VEmModel* theLivermorePEModel = new G4LivermorePhotoElectricModel();
      pe->SetEmModel(theLivermorePEModel,1);
      ph->RegisterProcess(pe, particle);

      // Compton scattering
      G4ComptonScattering* cs = new G4ComptonScattering;
      cs->SetEmModel(new G4KleinNishinaModel(),1);
      G4VEmModel* theLowEPComptonModel = new G4LowEPComptonModel();
      theLowEPComptonModel->SetHighEnergyLimit(20*MeV);
      cs->AddEmModel(0, theLowEPComptonModel);
      ph->RegisterProcess(cs, particle);

      // Gamma conversion
      G4GammaConversion* gc = new G4GammaConversion();
      G4VEmModel* thePenelopeGCModel = new G4PenelopeGammaConversionModel();
      thePenelopeGCModel->SetHighEnergyLimit(1*GeV);
      gc->SetEmModel(thePenelopeGCModel,1);
      ph->RegisterProcess(gc, particle);

      // Rayleigh scattering
      ph->RegisterProcess(new G4RayleighScattering(), particle);
 
    } else if (particleName == "e-") {

      // multiple scattering
      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel();
      G4CoulombScattering* ss = new G4CoulombScattering();
      ss->SetEmModel(ssm, 1); 

      // ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);
      G4PenelopeIonisationModel* pen = new G4PenelopeIonisationModel();
      pen->SetHighEnergyLimit(penEnergyLimit);
      eIoni->AddEmModel(0, pen, new G4UniversalFluctuation());

      // bremsstrahlung
      G4eBremsstrahlung* brem = new G4eBremsstrahlung();
      G4SeltzerBergerModel* br1 = new G4SeltzerBergerModel();
      G4eBremsstrahlungRelModel* br2 = new G4eBremsstrahlungRelModel();
      br1->SetAngularDistribution(new G4Generator2BS());
      br2->SetAngularDistribution(new G4Generator2BS());
      brem->SetEmModel(br1,1);
      brem->SetEmModel(br2,2);
      br2->SetLowEnergyLimit(GeV);

      // register processes
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(brem, particle);
      ph->RegisterProcess(ee, particle);
      ph->RegisterProcess(ss, particle);

    } else if (particleName == "e+") {

      // multiple scattering
      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel();
      G4CoulombScattering* ss = new G4CoulombScattering();
      ss->SetEmModel(ssm, 1); 

      // ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);
      G4PenelopeIonisationModel* pen = new G4PenelopeIonisationModel();
      pen->SetHighEnergyLimit(penEnergyLimit);
      eIoni->AddEmModel(0, pen, new G4UniversalFluctuation());

      // bremsstrahlung
      G4eBremsstrahlung* brem = new G4eBremsstrahlung();
      G4SeltzerBergerModel* br1 = new G4SeltzerBergerModel();
      G4eBremsstrahlungRelModel* br2 = new G4eBremsstrahlungRelModel();
      br1->SetAngularDistribution(new G4Generator2BS());
      br2->SetAngularDistribution(new G4Generator2BS());
      brem->SetEmModel(br1,1);
      brem->SetEmModel(br2,2);
      br2->SetLowEnergyLimit(GeV);

      // register processes
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(brem, particle);
      ph->RegisterProcess(ee, particle);
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);
      ph->RegisterProcess(ss, particle);

    } else if (particleName == "mu+" ||
               particleName == "mu-"    ) {

      G4MuIonisation* muIoni = new G4MuIonisation();
      muIoni->SetStepFunction(0.2, 50*um);          

      ph->RegisterProcess(muIoni, particle);
      ph->RegisterProcess(mub, particle);
      ph->RegisterProcess(mup, particle);
      ph->RegisterProcess(muss, particle);

    } else if (particleName == "alpha" || particleName == "He3") {
// COPY FROM OLD CHANNELING
            // ionisation
            G4ionIonisation* ionIoni = new G4ionIonisation();
	    ionIoni->SetEmModel(new G4IonParametrisedLossModel());
            ionIoni->SetStepFunction(0.1, 0.02*mm); // org 10 um
            ph->RegisterProcess(ionIoni, particle);

            // coulomb scattering
            G4CoulombScattering* ecs = new G4CoulombScattering();
            G4IonCoulombScatteringModel* mod = new G4IonCoulombScatteringModel();
	    ecs->AddEmModel(0,mod);
            ecs->SetBuildTableFlag(false);
            ph->RegisterProcess(ecs, particle);

		ph->RegisterProcess(pnuc, particle);
		
      pmanager->AddProcess(new G4StepLimiter,       -1,-1,3); 
    } else if (particleName == "GenericIon") {

// is older channeling
            // ionisation
            G4ionIonisation* ionIoni = new G4ionIonisation();
      	    //ionIoni->SetEmModel(new G4IonParametrisedLossModel());
            ionIoni->SetStepFunction(0.1, 0.02*mm);
            ph->RegisterProcess(ionIoni, particle);

            // coulomb scattering
            G4CoulombScattering* ecs = new G4CoulombScattering();
            G4IonCoulombScatteringModel* mod = new G4IonCoulombScatteringModel();
            mod->SetLowEnergyLimit(100.*keV);
            mod->SetRecoilThreshold(28.*eV);
            ecs->AddEmModel(0,mod);
            ecs->SetBuildTableFlag(false);
            ph->RegisterProcess(ecs, particle);

	    //ph->RegisterProcess(pnuc,particle);

      // step limiteris
      //G4ProcessManager* pmanager = particle->GetProcessManager();
      pmanager->AddProcess(new G4StepLimiter,       -1,-1,3); 
      // sudu demono pabaiga


		// new addition, gali buti blogas 2020-01-28
/*
    	G4double emax = 100.*TeV;

    	G4ExcitationHandler* handler = new G4ExcitationHandler();
    	G4PreCompoundModel* thePreCompound = new G4PreCompoundModel(handler);
    
    	// Binary Cascade
    	G4BinaryLightIonReaction* theIonBC = new G4BinaryLightIonReaction(thePreCompound);
    	theIonBC->SetMinEnergy(0.0);
    	theIonBC->SetMaxEnergy(4*GeV);
    
    	// FTFP
    	G4FTFBuilder* theBuilder = new G4FTFBuilder("FTFP",thePreCompound);
    	G4HadronicInteraction* theFTFP = theBuilder->GetModel();
    	theFTFP->SetMinEnergy(2*GeV);
    	theFTFP->SetMaxEnergy(emax);
    
    	G4VComponentCrossSection* theGGNuclNuclXS = new G4ComponentGGNuclNuclXsc();
    	G4VCrossSectionDataSet* theNuclNuclData = new G4CrossSectionInelastic(theGGNuclNuclXS);

    	G4HadronInelasticProcess* hadi = new G4HadronInelasticProcess("ionInelastic", particle);    
    	hadi->AddDataSet(theNuclNuclData);
    
    	hadi->RegisterMe(theIonBC);
    	hadi->RegisterMe(theFTFP);
    
    	pmanager->AddDiscreteProcess(hadi);*/
    	      pmanager->AddProcess(new G4StepLimiter,       -1,-1,3); 


    } else if (particleName == "pi+" ||
               particleName == "pi-" ) {

      //G4hMultipleScattering* pimsc = new G4hMultipleScattering();
      G4hIonisation* hIoni = new G4hIonisation();
      hIoni->SetStepFunction(0.2, 50*um);

      ph->RegisterProcess(hIoni, particle);
      ph->RegisterProcess(pib, particle);
      ph->RegisterProcess(pip, particle);
      ph->RegisterProcess(piss, particle);

    } else if (particleName == "kaon+" ||
               particleName == "kaon-" ) {

      //G4hMultipleScattering* kmsc = new G4hMultipleScattering();
      G4hIonisation* hIoni = new G4hIonisation();
      hIoni->SetStepFunction(0.2, 50*um);

      ph->RegisterProcess(hIoni, particle);
      ph->RegisterProcess(kb, particle);
      ph->RegisterProcess(kp, particle);
      ph->RegisterProcess(kss, particle);

    } else if (particleName == "proton" ||
               particleName == "anti_proton") {

// is older channeling
            //G4cout << particleName << G4endl;
            // ionisation
            G4hIonisation* hIoni = new G4hIonisation();
            hIoni->SetStepFunction(0.1, 2*um);
            ph->RegisterProcess(hIoni, particle);

            // bremsstrahlung
            G4hBremsstrahlung* hBrem = new G4hBremsstrahlung();
            ph->RegisterProcess(hBrem, particle);

            // pair production
            G4hPairProduction* hPair = new G4hPairProduction();
            ph->RegisterProcess(hPair, particle);

            // coulomb scattering
            G4CoulombScattering* ecs = new G4CoulombScattering();
            ecs->SetBuildTableFlag(false);
            ph->RegisterProcess(ecs, particle);

	    //ph->RegisterProcess(pnuc, particle);

		// single scattering, vietoje multiple scattering, jei amorfinis

		// addition from old channeling source, elastic shit

    	    G4ChipsElasticModel* chipsp = new G4ChipsElasticModel();
	    G4HadronElasticProcess* hel = new G4HadronElasticProcess();

            hel->AddDataSet(G4CrossSectionDataSetRegistry::
                    Instance()->GetCrossSectionDataSet(
                    G4ChipsProtonElasticXS::Default_Name()));

            hel->RegisterMe(chipsp);

            pmanager->AddDiscreteProcess(hel);
            if(verbose > 1) {
                G4cout << "### HadronElasticPhysics: " << hel->GetProcessName()
                << " added for " << particle->GetParticleName() << G4endl;
            }


      // step limiteris
      //G4ProcessManager* pmanager = particle->GetProcessManager();
      pmanager->AddProcess(new G4StepLimiter,       -1,-1,3); 



	    // inelastic part
/*

    	std::vector<G4VProtonBuilder *>::iterator i;
    	for(i=theModelCollections.begin(); i!=theModelCollections.end(); i++)
    	{
    		    (*i)->Build(theProtonInelastic);
    	}


		//G4ProtonInelasticProcess* theProtonInelastic=new G4ProtonInelasticProcess;

    		G4ProcessManager * theProcMan = G4Proton::Proton()->GetProcessManager();
    		XWrapperDiscreteProcess* theProtonInelastic_wrapper =
    	  	  new XWrapperDiscreteProcess("prot_Inelastic");
    		theProtonInelastic_wrapper->RegisterProcess(theProtonInelastic,1);
    		theProcMan->AddDiscreteProcess(theProtonInelastic_wrapper);
*/

		// end of addition


    } else if (particleName == "B+" ||
               particleName == "B-" ||
               particleName == "D+" ||
               particleName == "D-" ||
               particleName == "Ds+" ||
               particleName == "Ds-" ||
               particleName == "anti_He3" ||
               particleName == "anti_alpha" ||
               particleName == "anti_deuteron" ||
               particleName == "anti_lambda_c+" ||
               particleName == "anti_omega-" ||
               particleName == "anti_sigma_c+" ||
               particleName == "anti_sigma_c++" ||
               particleName == "anti_sigma+" ||
               particleName == "anti_sigma-" ||
               particleName == "anti_triton" ||
               particleName == "anti_xi_c+" ||
               particleName == "anti_xi-" ||
               particleName == "deuteron" ||
               particleName == "lambda_c+" ||
               particleName == "omega-" ||
               particleName == "sigma_c+" ||
               particleName == "sigma_c++" ||
               particleName == "sigma+" ||
               particleName == "sigma-" ||
               particleName == "tau+" ||
               particleName == "tau-" ||
               particleName == "triton" ||
               particleName == "xi_c+" ||
               particleName == "xi-" ) {

      ph->RegisterProcess(new G4CoulombScattering(), particle);
      ph->RegisterProcess(new G4hIonisation(), particle);
      //ph->RegisterProcess(pnuc, particle);
    }
  }
    
  // Nuclear stopping
  pnuc->SetMaxKinEnergy(MeV);

  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  de->SetFluo(true);
  de->SetAuger(true);   
  de->SetPIXE(true);  
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);

  G4EmModelActivator mact(GetPhysicsName());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
