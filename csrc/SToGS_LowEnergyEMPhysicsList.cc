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

#include "ParisLowEnergyEMPhysicsList.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"

// gamma
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"

#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"

// e-
#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh"

#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"


ParisLowEnergyEMPhysicsList::ParisLowEnergyEMPhysicsList(const G4String& name):  G4VPhysicsConstructor(name)
{ 
}

ParisLowEnergyEMPhysicsList::~ParisLowEnergyEMPhysicsList()
{
}

void ParisLowEnergyEMPhysicsList::ConstructParticle()
{
	// gamma
	G4Gamma::GammaDefinition();
	
	// electron
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();
	
	// pseudo-particles
	G4Geantino::GeantinoDefinition();
}


void ParisLowEnergyEMPhysicsList::ConstructProcess()
{
	theParticleIterator->reset();
	//
	while( (*theParticleIterator)() ){
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager(); 
		G4String particleName = particle->GetParticleName();
		
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
			
		} else if (particleName == "e-") {
			
			G4eMultipleScattering* msc = new G4eMultipleScattering();
			
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
			
		} else if (particleName == "e+") {
			
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

		}
	}
}

/*
void ParisLowEnergyEMPhysicsList::ConstructEM()
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

void ParisLowEnergyEMPhysicsList::SetCuts()
{
	/*
  if ( verboseLevel > 0 ) {
    G4cout << " ParisLowEnergyEMPhysicsList::SetCuts:";
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



 
