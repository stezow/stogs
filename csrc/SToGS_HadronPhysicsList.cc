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

#include "SToGS_HadronPhysicsList.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4Version.hh"

#include <iomanip>

#include "G4LossTableManager.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4LossTableManager.hh"
#include "G4PhysicsListHelper.hh"

// PARTICLE DEFINITIONS
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Neutron.hh"
// to deal with version dependant interfaces
#if G4VERSION_NUMBER < 1000
#else
#endif

// PROCESSES FOR GAMMA
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"
// PROCESSES FOR e-
#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh"
#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4KleinNishinaModel.hh"
#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"
#include "G4SystemOfUnits.hh"
// to deal with version dependant interfaces
#if G4VERSION_NUMBER < 1000
    #if G4VERSION_NUMBER >= 960
        #include "G4UrbanMscModel96.hh"
        #define URBANMSCMODEL G4UrbanMscModel96
    #else
        #include "G4UrbanMscModel95.hh"
        #define URBANMSCMODEL G4UrbanMscModel95
    #endif
#else
#include "G4UrbanMscModel.hh"
#define URBANMSCMODEL G4UrbanMscModel
#endif

// hadronic processes
#include "G4HadronElasticProcess.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronElastic.hh"
// to deal with version dependant interfaces
#if G4VERSION_NUMBER < 1000
#else
#include "FTFP_BERT_HP.hh"
#endif

//HPNeutron
#include "G4NeutronHPThermalScatteringData.hh"
#include "G4NeutronHPThermalScattering.hh"
//
#include "G4HadronFissionProcess.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPFission.hh"
//
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
//
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
//
#include "G4NeutronInelasticProcess.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
#if G4VERSION_NUMBER < 1000
#else
#endif

//c-s
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"

// RadioactiveDecay
#include "G4RadioactiveDecay.hh"
#include "G4GenericIon.hh"

SToGS::HadronPhysicsList::HadronPhysicsList(const G4String& name):  G4VPhysicsConstructor(name)
{
}

SToGS::HadronPhysicsList::~HadronPhysicsList()
{
  
}

void SToGS:: HadronPhysicsList::ConstructParticle()
{

  G4LeptonConstructor lepton;
  lepton.ConstructParticle();
 
  G4BosonConstructor boson;
  boson.ConstructParticle();

  G4MesonConstructor meson;
  meson.ConstructParticle();

  G4BaryonConstructor baryon;
  baryon.ConstructParticle();

  G4ShortLivedConstructor shortLived;
  shortLived.ConstructParticle();

  G4IonConstructor ion;
  ion.ConstructParticle();


  G4cout << " ------ INFO ------ PARTICLES DEFINITION FINISHED" << G4endl;
}


void SToGS::HadronPhysicsList::ConstructProcess()
{
    /* OS TODO - figure out
  theParticleIterator->reset();
  //
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager(); 
    G4String particleName = particle->GetParticleName();
	
    //define EM physics for gamma
    if (particleName == "gamma") {
			
      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      thePhotoElectricEffect->SetModel(new G4LivermorePhotoElectricModel());
			
      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      theComptonScattering->SetModel(new G4LivermoreComptonModel());
			
      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      theGammaConversion->SetModel(new G4LivermoreGammaConversionModel());
			
      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      theRayleigh->SetModel(new G4LivermoreRayleighModel());
			
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);
      pmanager->AddDiscreteProcess(theComptonScattering);
      pmanager->AddDiscreteProcess(theGammaConversion);
      pmanager->AddDiscreteProcess(theRayleigh);

      //define EM physics for electron			
    } else if (particleName == "e-") {
			
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc -> AddEmModel(0, new URBANMSCMODEL() );
		
      // Ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetEmModel(new G4LivermoreIonisationModel());
      eIoni->SetFluctModel(new G4UniversalFluctuation() );
			
      // Bremsstrahlung
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      eBrem->SetEmModel(new G4LivermoreBremsstrahlungModel());
			
      pmanager->AddProcess(msc   ,-1, 1,1);
      pmanager->AddProcess(eIoni ,-1, 2,2);
      pmanager->AddProcess(eBrem ,-1,-1,3);
     
    }
     */
 /*   // define EM physics for positron
    else if (particleName == "e+") {
			
      G4eMultipleScattering* msc = new G4eMultipleScattering();
			
      // Ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetEmModel(new G4LivermoreIonisationModel());
      eIoni->SetFluctModel(new G4UniversalFluctuation());
			
      // Bremsstrahlung
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      eBrem->SetEmModel(new G4LivermoreBremsstrahlungModel());
			
      //Annihilation
      G4eplusAnnihilation* eAnni = new G4eplusAnnihilation();
			
      pmanager->AddProcess(msc   ,-1, 1,1);
      pmanager->AddProcess(eIoni ,-1, 2,2);
      pmanager->AddProcess(eBrem ,-1,-1,3);
      pmanager->AddProcess(eAnni , 0,-1,4);

    }*/
//netronHP physics
  
/*  // define physics for neutron
  else if (particleName == "neutron") {
    // process: elastic
    //
    G4HadronElasticProcess* process1 = new G4HadronElasticProcess();
    pManager->AddDiscreteProcess(process1);   
    //
    // cross section data set
    G4NeutronHPElasticData* dataSet1a = new G4NeutronHPElasticData();
    G4NeutronHPThermalScatteringData* dataSet1b 
      = new G4NeutronHPThermalScatteringData();
    process1->AddDataSet(dataSet1a);                               
    if (fThermal) process1->AddDataSet(dataSet1b);
    //
    // models
    G4NeutronHPElastic*           model1a = new G4NeutronHPElastic();
    G4NeutronHPThermalScattering* model1b = new G4NeutronHPThermalScattering();
    if (fThermal)  model1a->SetMinEnergy(4*eV);
    process1->RegisterMe(model1a);    
    if (fThermal) process1->RegisterMe(model1b);
   
    // process: inelastic
    //
    G4NeutronInelasticProcess* process2 = new G4NeutronInelasticProcess();
    pManager->AddDiscreteProcess(process2);   
    //
    // cross section data set
    G4NeutronHPInelasticData* dataSet2 = new G4NeutronHPInelasticData();
    process2->AddDataSet(dataSet2);                               
   //
   // models
   G4NeutronHPInelastic* model2 = new G4NeutronHPInelastic();
   process2->RegisterMe(model2);    

   // process: nCapture   
   //
   G4HadronCaptureProcess* process3 = new G4HadronCaptureProcess();
   pManager->AddDiscreteProcess(process3);    
   //
   // cross section data set
   G4NeutronHPCaptureData* dataSet3 = new G4NeutronHPCaptureData();
   process3->AddDataSet(dataSet3);                               
   //
   // models
   G4NeutronHPCapture* model3 = new G4NeutronHPCapture();
   process3->RegisterMe(model3);
   
   // process: nFission   
   //
   G4HadronFissionProcess* process4 = new G4HadronFissionProcess();
   pManager->AddDiscreteProcess(process4);    
   //
  
   G4ProcessManager* pManager = G4Neutron::Neutron()->GetProcessManager();
   
   // process: elastic
   //
   G4HadronElasticProcess* process1 = new G4HadronElasticProcess();
   pManager->AddDiscreteProcess(process1);   
   //
   // cross section data set
   G4NeutronHPElasticData* dataSet1a = new G4NeutronHPElasticData();
   G4NeutronHPThermalScatteringData* dataSet1b 
                               = new G4NeutronHPThermalScatteringData();
   process1->AddDataSet(dataSet1a);                               
   if (fThermal) process1->AddDataSet(dataSet1b);
   //
   // models
   G4NeutronHPElastic*           model1a = new G4NeutronHPElastic();
   G4NeutronHPThermalScattering* model1b = new G4NeutronHPThermalScattering();
  if (fThermal)  model1a->SetMinEnergy(4*eV);
   process1->RegisterMe(model1a);    
   if (fThermal) process1->RegisterMe(model1b);
   
   // process: inelastic
   //
   G4NeutronInelasticProcess* process2 = new G4NeutronInelasticProcess();
   pManager->AddDiscreteProcess(process2);   
   //
   // cross section data set
   G4NeutronHPInelasticData* dataSet2 = new G4NeutronHPInelasticData();
   process2->AddDataSet(dataSet2);                               
   //
   // models
   G4NeutronHPInelastic* model2 = new G4NeutronHPInelastic();
   process2->RegisterMe(model2);    

   // process: nCapture   
   //
   G4HadronCaptureProcess* process3 = new G4HadronCaptureProcess();
   pManager->AddDiscreteProcess(process3);    
   //
   // cross section data set
   G4NeutronHPCaptureData* dataSet3 = new G4NeutronHPCaptureData();
   process3->AddDataSet(dataSet3);                               
   //
   // models
   G4NeutronHPCapture* model3 = new G4NeutronHPCapture();
   process3->RegisterMe(model3);
   
   // process: nFission   
   //
   G4HadronFissionProcess* process4 = new G4HadronFissionProcess();
   pManager->AddDiscreteProcess(process4);    
   //
   // cross section data set
   G4NeutronHPFissionData* dataSet4 = new G4NeutronHPFissionData();
   process4->AddDataSet(dataSet4);                               
   //
   // models
   G4NeutronHPFission* model4 = new G4NeutronHPFission();
   process4->RegisterMe(model4);   

 }*/
 // }
    // ***********************************************
    // *             HADRONIC PHYSICS                *  
    // *********************************************** 

  // define hadronic elastic processes
    /*
  G4HadronElasticPhysicsHP* hadronicElastic = new G4HadronElasticPhysicsHP();
  hadronicElastic->ConstructProcess();

  // define hadronic inelastic processes
  HadronPhysicsFTFP_BERT_HP* hadronicPhysics = new HadronPhysicsFTFP_BERT_HP();
  hadronicPhysics->ConstructProcess();
  */
   
  
}

/*
void SToGS:: HadronPhysicsList::ConstructEM()
{
//   G4VProcess *aprocess;
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    
    if (particleName == "gamma") {     
      pmanager->AddDiscreteProcess(new G4LowEnergyPhotoElectric);
      pmanager->AddDiscreteProcess(new G4LowEnergyCompton);
      pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion);
      pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh);
      
    } else if (particleName == "e-") {
      pmanager->AddProcess(new G4eMultipleScattering,        -1, 1,1);
      pmanager->AddProcess(new G4LowEnergyIonisation,       -1, 2,2);
      pmanager->AddProcess(new G4LowEnergyBremsstrahlung,   -1,-1,3);      

    } else if (particleName == "e+") { 
      pmanager->AddProcess(new G4eMultipleScattering,        -1, 1,1);
      pmanager->AddProcess(new G4eIonisation,               -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,           -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,          0,-1,4);

    } 
  }
}
*/

void SToGS::HadronPhysicsList::SetCuts()
{
	/*
  if ( verboseLevel > 0 ) {
    G4cout << " HadronPhysicsList::SetCuts:";
    G4cout << " CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  //
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");

  if (verboseLevel>0) DumpCutValuesTable();
	 */
}



 
