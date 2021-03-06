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

// to deal with version dependant interfaces
#if G4VERSION_NUMBER < 1000

//*******************************************************
//*                     GEANT4.9 VERSION               *
//*******************************************************

#include "G4LossTableManager.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
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
#include "G4HadronElasticPhysics.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronElastic.hh"
#include "G4HadronElasticProcess.hh"

//c-s
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"

// RadioactiveDecay
#include "G4RadioactiveDecay.hh"
#include "G4GenericIon.hh"

#include "G4UrbanMscModel96.hh"

// hadronic

#include "G4HadronElasticPhysicsHP.hh"
#include "HadronPhysicsFTFP_BERT_HP.hh"

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



/*



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
*/




SToGS::HadronPhysicsList::HadronPhysicsList(const G4String& name):  G4VPhysicsConstructor(name)
{
}

SToGS::HadronPhysicsList::~HadronPhysicsList()
{
    G4LossTableManager::Instance();
  //  defaultCutValue = 1.0*mm;

    G4cout << G4endl << " ------ INFO ------ You are working with an hadronic PhysicsList (transportation, electromagnetic and hadronic processes)" << G4endl;
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


G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();


theParticleIterator->reset();
//
while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();





    // define EM physics for gamma
    if (particleName == "gamma") {

      // activate photoelectric effect
      ph->RegisterProcess(new G4PhotoElectricEffect, particle);

      // activate Compton scattering
      G4ComptonScattering* cs   = new G4ComptonScattering;
      cs->SetEmModel(new G4KleinNishinaModel());
      ph->RegisterProcess(cs, particle);

      // activate gamma conversion
      ph->RegisterProcess(new G4GammaConversion, particle);

    }

    // define EM physics for electron
    else if (particleName == "e-") {

      // activate multiple scattering
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc -> AddEmModel(0, new G4UrbanMscModel96());
      ph->RegisterProcess(msc, particle);

      // activate ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.1, 100*um);
      ph->RegisterProcess(eIoni, particle);

      // activate bremsstrahlung
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);

    }

    // define EM physics for positron
    else if (particleName == "e+") {

      // activate multiple scattering
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc -> AddEmModel(0, new G4UrbanMscModel96());
      ph->RegisterProcess(msc, particle);

      // activate ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.1, 100*um);
      ph->RegisterProcess(eIoni, particle);

      // activate bremsstrahlung
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);

      // activate annihilation
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);

    }

    // define EM physics for muons
    else if (particleName == "mu+" || particleName == "mu-") {

      // activate multiple scattering
      ph->RegisterProcess(new G4MuMultipleScattering(), particle);

      // activate ionisation
      G4MuIonisation* muIoni = new G4MuIonisation();
      muIoni->SetStepFunction(0.1, 50*um);
      ph->RegisterProcess(muIoni, particle);

      // activate bremsstrahlung
      ph->RegisterProcess(new G4MuBremsstrahlung(), particle);

      // activate pair production
      ph->RegisterProcess(new G4MuPairProduction(), particle);

    }

    // define EM physics for protons and pions
    else if( particleName == "proton" ||
         particleName == "pi-" ||
         particleName == "pi+") {

      // activate multiple scattering
      ph->RegisterProcess(new G4hMultipleScattering(), particle);

      // activate ionisation
      G4hIonisation* hIoni = new G4hIonisation();
      hIoni->SetStepFunction(0.1, 20*um);
      ph->RegisterProcess(hIoni, particle);

      // activate bremsstrahlung
      ph->RegisterProcess(new G4hBremsstrahlung(), particle);

      // activate pair production
      ph->RegisterProcess(new G4hPairProduction(), particle);

    }

    // define EM physics for alpha and helium-3
    else if( particleName == "alpha" ||
         particleName == "He3") {

      // activate multiple scattering
      ph->RegisterProcess(new G4hMultipleScattering(), particle);

      // activate ionisation
      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetStepFunction(0.1, 1*um);
      ph->RegisterProcess(ionIoni, particle);

      // activate nuclear stopping
      ph->RegisterProcess(new G4NuclearStopping(), particle);

    }

    // define EM physics for all the other ions
    else if( particleName == "GenericIon" ) {

      // activate multiple scattering
      ph->RegisterProcess(new G4hMultipleScattering(), particle);

      // activate ionisation
      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ionIoni->SetStepFunction(0.1, 1*um);
      ph->RegisterProcess(ionIoni, particle);

      // activate nuclear stopping
      ph->RegisterProcess(new G4NuclearStopping(), particle);

    }

    // activate EM physics for all others charged particles except geantino
    else if ((!particle->IsShortLived()) &&
         (particle->GetPDGCharge() != 0.0) &&
         (particle->GetParticleName() != "chargedgeantino")) {

      // activate multiple scattering
      ph->RegisterProcess(new G4hMultipleScattering(), particle);

      // activate ionisation
      ph->RegisterProcess(new G4hIonisation(), particle);
    }
  }




  /* ***********************************************
     *                HADRONIC PHYSICS             *
     *********************************************** */


  // define hadronic elastic processes
  G4HadronElasticPhysicsHP* hadronicElastic = new G4HadronElasticPhysicsHP();
  hadronicElastic->ConstructProcess();

  // define hadronic inelastic processes
  HadronPhysicsFTFP_BERT_HP* hadronicPhysics = new HadronPhysicsFTFP_BERT_HP();
  hadronicPhysics->ConstructProcess();
}

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




#else

//*******************************************************
//*                     GEANT4.10 VERSION               *
//*******************************************************

//#if G4VERSION_NUMBER >= 1000
#include "G4LossTableManager.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
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
#include "G4HadronElasticPhysics.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronElastic.hh"
#include "G4HadronElasticProcess.hh"

//c-s
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"

// RadioactiveDecay
#include "G4RadioactiveDecay.hh"
#include "G4GenericIon.hh"

#include "G4UrbanMscModel.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "QGSP_BERT_HP.hh"

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



SToGS::HadronPhysicsList::HadronPhysicsList(const G4String& name):  G4VPhysicsConstructor(name)
{
}

SToGS::HadronPhysicsList::~HadronPhysicsList()
{
    G4LossTableManager::Instance();
    //  defaultCutValue = 1.0*CLHEP::mm;

    G4cout << G4endl << " ------ INFO ------ You are working with an hadronic PhysicsList (transportation, electromagnetic and hadronic processes)" << G4endl;
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


    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

    theParticleTable->GetIterator()->reset();

    //
    while( (*theParticleTable->GetIterator())() ){
        G4ParticleDefinition* particle = theParticleTable->GetIterator()->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

        /*theParticleIterator->reset();
    //
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

*/



        // define EM physics for gamma
        if (particleName == "gamma") {

            // activate photoelectric effect
            ph->RegisterProcess(new G4PhotoElectricEffect, particle);

            // activate Compton scattering
            G4ComptonScattering* cs   = new G4ComptonScattering;
            cs->SetEmModel(new G4KleinNishinaModel());
            ph->RegisterProcess(cs, particle);

            // activate gamma conversion
            ph->RegisterProcess(new G4GammaConversion, particle);

        }


        // define EM physics for electron
            else if (particleName == "e-") {

              // activate multiple scattering
              G4eMultipleScattering* msc = new G4eMultipleScattering();
              msc -> AddEmModel(0, new G4UrbanMscModel());
              ph->RegisterProcess(msc, particle);

              // activate ionisation
              G4eIonisation* eIoni = new G4eIonisation();
              eIoni->SetStepFunction(0.1, 100*CLHEP::um);
              ph->RegisterProcess(eIoni, particle);

              // activate bremsstrahlung
              ph->RegisterProcess(new G4eBremsstrahlung(), particle);

            }

        // define EM physics for positron
            else if (particleName == "e+") {

              // activate multiple scattering
              G4eMultipleScattering* msc = new G4eMultipleScattering();
              msc -> AddEmModel(0, new G4UrbanMscModel());
              ph->RegisterProcess(msc, particle);

              // activate ionisation
              G4eIonisation* eIoni = new G4eIonisation();
              eIoni->SetStepFunction(0.1, 100*CLHEP::um);
              ph->RegisterProcess(eIoni, particle);

              // activate bremsstrahlung
              ph->RegisterProcess(new G4eBremsstrahlung(), particle);

              // activate annihilation
              ph->RegisterProcess(new G4eplusAnnihilation(), particle);

            }


        // define EM physics for muons
            else if (particleName == "mu+" || particleName == "mu-") {

              // activate multiple scattering
              ph->RegisterProcess(new G4MuMultipleScattering(), particle);

              // activate ionisation
              G4MuIonisation* muIoni = new G4MuIonisation();
              muIoni->SetStepFunction(0.1, 50*CLHEP::um);
              ph->RegisterProcess(muIoni, particle);

              // activate bremsstrahlung
              ph->RegisterProcess(new G4MuBremsstrahlung(), particle);

              // activate pair production
              ph->RegisterProcess(new G4MuPairProduction(), particle);

            }


        // define EM physics for protons and pions
            else if( particleName == "proton" ||
                 particleName == "pi-" ||
                 particleName == "pi+") {

              // activate multiple scattering
              ph->RegisterProcess(new G4hMultipleScattering(), particle);

              // activate ionisation
              G4hIonisation* hIoni = new G4hIonisation();
              hIoni->SetStepFunction(0.1, 20*CLHEP::um);
              ph->RegisterProcess(hIoni, particle);

              // activate bremsstrahlung
              ph->RegisterProcess(new G4hBremsstrahlung(), particle);

              // activate pair production
              ph->RegisterProcess(new G4hPairProduction(), particle);

            }

            // define EM physics for alpha and helium-3
            else if( particleName == "alpha" ||
                 particleName == "He3") {

              // activate multiple scattering
              ph->RegisterProcess(new G4hMultipleScattering(), particle);

              // activate ionisation
              G4ionIonisation* ionIoni = new G4ionIonisation();
              ionIoni->SetStepFunction(0.1, 1*CLHEP::um);
              ph->RegisterProcess(ionIoni, particle);

              // activate nuclear stopping
              ph->RegisterProcess(new G4NuclearStopping(), particle);

            }

        // define EM physics for all the other ions
            else if( particleName == "GenericIon" ) {

              // activate multiple scattering
              ph->RegisterProcess(new G4hMultipleScattering(), particle);

              // activate ionisation
              G4ionIonisation* ionIoni = new G4ionIonisation();
              ionIoni->SetEmModel(new G4IonParametrisedLossModel());
              ionIoni->SetStepFunction(0.1, 1*CLHEP::um);
              ph->RegisterProcess(ionIoni, particle);

              // activate nuclear stopping
              ph->RegisterProcess(new G4NuclearStopping(), particle);

            }

            // activate EM physics for all others charged particles except geantino
            else if ((!particle->IsShortLived()) &&
                 (particle->GetPDGCharge() != 0.0) &&
                 (particle->GetParticleName() != "chargedgeantino")) {

              // activate multiple scattering
              ph->RegisterProcess(new G4hMultipleScattering(), particle);

              // activate ionisation
              ph->RegisterProcess(new G4hIonisation(), particle);
            }



    }

    /* ***********************************************
     *                HADRONIC PHYSICS             *
     *********************************************** */


    // define hadronic elastic processes
    G4HadronElasticPhysicsHP* hadronicElastic = new G4HadronElasticPhysicsHP();
    hadronicElastic->ConstructProcess();

    // define hadronic inelastic processes
    QGSP_BERT_HP* hadronicPhysics = new QGSP_BERT_HP();
    hadronicPhysics->ConstructProcess();
}

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

#endif
    
    
