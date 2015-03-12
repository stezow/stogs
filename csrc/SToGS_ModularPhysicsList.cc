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
// $Id: SToGS::PhysicsList.cc,v 1.39 2010-06-04 15:42:23 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-04-patch-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SToGS_ModularPhysicsList.hh"
#include "SToGS_G4_GeneralPhysics.hh"

#include "G4Version.hh"

//Electromagnetic
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmConfigurator.hh"
#include "G4UnitsTable.hh"

// transport and decay
#include "G4ProcessManager.hh"
#include "G4Decay.hh"

// scintillation
#include "G4OpticalPhysics.hh"

// hadrons
#include "SToGS_HadronPhysicsList.hh"
#if G4VERSION_NUMBER < 1000
#include "HadronPhysicsQGSP_BIC_HP.hh"
#define QGSP_BIC_HP_MODEL HadronPhysicsQGSP_BIC_HP
#else
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#define QGSP_BIC_HP_MODEL G4HadronPhysicsQGSP_BIC_HP
#endif

/*
#include "SToGS_StandardEMPhysicsList.hh"
#include "SToGS_LowEnergyEMPhysicsList.hh"
#include "SToGS_PenelopeEMPhysicsList.hh"
#include "SToGS_HadronPhysicsList.hh"
#include "SToGS::Hadron1.hh"
 */


SToGS::ModularPhysicsList::ModularPhysicsList(const G4String &option) : G4VModularPhysicsList()
{
    G4cout << " ------ INF ------ from SToGS::ModularPhysicsList::ModularPhysicsList with " << option << G4endl;

	// parse option
	std::vector <G4String> all_opt; G4String lopt = option, tmp = ""; lopt += ';';
	//
	for (size_t i = 0; i < lopt.size(); i++) { // search for one sequence
		if ( lopt[i] == ';' && tmp.size() > 0 ) {
			all_opt.push_back(tmp);
			tmp = "";
		}
		else {
			tmp += G4String(lopt[i]); 
		}
	}
	
	// default cut value  (1.0mm) 
	defaultCutValue = 1.0*CLHEP::mm;
	
	// General Physics
	for (size_t i = 0; i < all_opt.size(); i++) {
		// Register general physics list
		if (all_opt[i].contains(G4String("general")) ) {
			//
			RegisterPhysics( new SToGS::GeneralPhysics(all_opt[i]) );
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}
		// G4 EM physics list
		if ( all_opt[i] == "emstandard_opt0" ) {
			RegisterPhysics( new G4EmStandardPhysics(1) );
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}
		if ( all_opt[i] == "emstandard_opt1" ) {
			RegisterPhysics( new G4EmStandardPhysics_option1() );
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}	
		if ( all_opt[i] == "emstandard_opt2" ) {
			RegisterPhysics( new G4EmStandardPhysics_option2() );
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}
		if ( all_opt[i] == "emstandard_opt3" ) {
			RegisterPhysics( new G4EmStandardPhysics_option3() );
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}			
		if ( all_opt[i] == "emlivermore" ) {
			RegisterPhysics( new G4EmLivermorePhysics() );
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}
		if ( all_opt[i] == "empenelope" ) {
			RegisterPhysics( new G4EmPenelopePhysics() );
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}

//SToGS:: hadron physics for list
        if ( all_opt[i] == "SToGS_hadron" ) {
			RegisterPhysics( new SToGS::HadronPhysicsList() );
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}
        if ( all_opt[i] == "QGSP_BIC_HP" ) {
            RegisterPhysics( new QGSP_BIC_HP_MODEL() );
		}
		// SToGS:: EM physics list
        /*
		if ( all_opt[i] == "SToGS::StandardEM" ) {
			RegisterPhysics( new SToGS::StandardEMPhysicsList() );
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}
		if ( all_opt[i] == "SToGS::LowEnergyEM" ) {
			RegisterPhysics( new SToGS::LowEnergyEMPhysicsList() );
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}
		if ( all_opt[i] == "SToGS::PenelopeEM" ) {
			RegisterPhysics( new SToGS::LowEnergyEMPhysicsList() );
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}
         */
		if ( all_opt[i] == "Optical" ) {
			
			G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
			RegisterPhysics( opticalPhysics );
			
			// in case Wavelength Shifting (WLS) fibers are used. It means the time between abs. and re-emission is constant equal to WLSTIMECONSTANT
			// other possibility could be exponential
			// opticalPhysics->SetWLSTimeProfile("delta");
			
			// to allow rising time ==> 	FASTSCINTILLATIONRISETIME and SLOWSCINTILLATIONRISETIME to be defined

			opticalPhysics->SetFiniteRiseTime(true);
			
			//opticalPhysics->SetScintillationYieldFactor(1.0);
			// 1.0 so that optical goes to two components in case of fast and slow components .... not clear why .... ??
			opticalPhysics->SetScintillationExcitationRatio(1.0);
			
			opticalPhysics->SetMaxNumPhotonsPerStep(100);
			opticalPhysics->SetMaxBetaChangePerStep(10.0);
			
			opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
			opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}

        /*
		if ( all_opt[i] == "SToGS::Hadron0" ) {
			RegisterPhysics( new SToGS::Hadron0() );
			
    
			RegisterPhysics( new HadronPhysicsQGSP_BIC());
			RegisterPhysics( new G4EmExtraPhysics());
			RegisterPhysics( new G4HadronElasticPhysics());
			RegisterPhysics( new G4QStoppingPhysics());
			RegisterPhysics( new G4IonBinaryCascadePhysics());
			RegisterPhysics( new G4NeutronTrackingCut());		
        
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}		
		if ( all_opt[i] == "SToGS::Hadron1" ) {
			RegisterPhysics( new SToGS::Hadron1() );
			
        
			 RegisterPhysics( new HadronPhysicsQGSP_BIC());
			 RegisterPhysics( new G4EmExtraPhysics());
			 RegisterPhysics( new G4HadronElasticPhysics());
			 RegisterPhysics( new G4QStoppingPhysics());
			 RegisterPhysics( new G4IonBinaryCascadePhysics());
			 RegisterPhysics( new G4NeutronTrackingCut());		
         
			//
			G4cout << " ==> Physics list " << all_opt[i] << " registered " << G4endl;
		}		
         */
	}
    G4cout <<  " ------ END ------ from SToGS::ModularPhysicsList::ModularPhysicsList " << G4endl ;
}


SToGS::ModularPhysicsList::~ModularPhysicsList()
{
}


/*
void SToGS::PhysicsList::AddIonGasModels()
{
	G4EmConfigurator* em_config = G4LossTableManager::Instance()->EmConfigurator();
	theParticleIterator->reset();
	while ((*theParticleIterator)())
	{
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4String partname = particle->GetParticleName();
		if(partname == "alpha" || partname == "He3" || partname == "GenericIon") {
			G4BraggIonGasModel* mod1 = new G4BraggIonGasModel();
			G4BetheBlochIonGasModel* mod2 = new G4BetheBlochIonGasModel();
			G4double eth = 2.*MeV*particle->GetPDGMass()/proton_mass_c2;
			em_config->SetExtraEmModel(partname,"ionIoni",mod1,"",0.0,eth,
												new G4IonFluctuations());
			em_config->SetExtraEmModel(partname,"ionIoni",mod2,"",eth,100*TeV,
												new G4UniversalFluctuation());
			
		}
	}
}
*/

