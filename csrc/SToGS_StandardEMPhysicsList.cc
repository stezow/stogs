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

#include "ParisStandardEMPhysicsList.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4ios.hh"

// gamma
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
//#include "G4PolarizedComptonScattering.hh"

// charged
#include "G4eMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"


ParisStandardEMPhysicsList::ParisStandardEMPhysicsList(const G4String& name):  G4VPhysicsConstructor(name)
{ 
}

ParisStandardEMPhysicsList::~ParisStandardEMPhysicsList()
{
}

void ParisStandardEMPhysicsList::ConstructParticle()
{
	// gamma
	G4Gamma::GammaDefinition();
	
	// electron
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();
	
	// pseudo-particles
	G4Geantino::GeantinoDefinition();
}

void ParisStandardEMPhysicsList::ConstructProcess()
{	
	// add processes to the manager for each kind of particles
	G4ProcessManager * pManager = 0;
	
	// Gamma Physics : all processes
	G4PhotoElectricEffect* thePhotoEffect = new G4PhotoElectricEffect();
	G4ComptonScattering* theComptonEffect = new G4ComptonScattering();
	G4GammaConversion* thePairProduction = new G4GammaConversion();
	// 
	pManager = G4Gamma::Gamma()->GetProcessManager();
	//
	pManager->AddDiscreteProcess(thePhotoEffect);
	pManager->AddDiscreteProcess(theComptonEffect);
	pManager->AddDiscreteProcess(thePairProduction);
	
	// Electron physics : all processes
	G4eMultipleScattering* theElectronMultipleScattering = new G4eMultipleScattering();
	G4eIonisation* theElectronIonisation = new G4eIonisation();
	G4eBremsstrahlung* theElectronBremsStrahlung = new G4eBremsstrahlung();
	//
	pManager = G4Electron::Electron()->GetProcessManager();
	//
	pManager->AddProcess(theElectronMultipleScattering, -1, 1, 1);
	pManager->AddProcess(theElectronIonisation,         -1, 2, 2);
	pManager->AddProcess(theElectronBremsStrahlung,     -1, 3, 3); 
	
	// Positron physics : all processes
	G4eMultipleScattering* thePositronMultipleScattering = new G4eMultipleScattering();
	G4eIonisation* thePositronIonisation = new G4eIonisation(); 
	G4eBremsstrahlung* thePositronBremsStrahlung = new G4eBremsstrahlung();  
	G4eplusAnnihilation* theAnnihilation = new G4eplusAnnihilation();
	//
	pManager = G4Positron::Positron()->GetProcessManager();
	//
	pManager->AddProcess(thePositronMultipleScattering, -1, 1, 1);
	pManager->AddProcess(thePositronIonisation,         -1, 2, 2);
	pManager->AddProcess(thePositronBremsStrahlung,     -1, 3, 3);  
	pManager->AddProcess(theAnnihilation,                0,-1, 4);  
}


void ParisStandardEMPhysicsList::SetCuts()
{
	/*
  if ( verboseLevel > 0 ) {
    G4cout << " ParisStandardEMPhysicsList::SetCuts:";
    G4cout << " CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  //
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");

  if (verboseLevel>0) 
		DumpCutValuesTable();
	 
	 */
}



 
