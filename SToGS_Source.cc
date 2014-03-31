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
// --------------------------------------------------------------
//      GEANT 4 -
//
//      For information related to this code contact:
//
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "SToGS_UserActionManager.hh"

// TMP to test MIGRATION
#include "SToGS_TwoShellsDetectorConstruction.hh"
#include "SToGS_ModularPhysicsList.hh"
#include "SToGS_PrintOut.hh"

//! This is the main program of the simulation for 'Sources'
/*!
 
    It is used for detector characterization/test but can be used also for 'simple' experiments
    In this case, 'simple' means the source of events does not have any temporal structure.
    It thus excludes the case of beams with event number / time stamp and treatments of radiactivity created by it
 
 \code
 # This is the default file used to configure the SToGS_Source program.
 # It is read at the beginning of the program to select
 #
 # The UserActionInitialization (which deals without outputs and the generator [thread local objects])
 #
 actions: printout run;event;track;step
 #actions: ascii setup/SToGS_ascii_actions.conf
 #
 # Detector geometry
 #
 setup:  
 #
 # The Physics list
 #
 physics: general0;emstandard_opt0;
 #physics: general0;emstandard_opt0;Optical;
 #physics: general0;emstandard_opt0;ParisHadron0
 #
 # generator
 #
 generator: GPS G4Macros/GPS_Cs137.mac
 #

 \endcode
 
 Without arguments, the program starts a G4 interactive session. It is also possible to run it
 in batch mode by adding on the command line the name of the G4 macro you would like to execute. Ex:
 
 \c SToGS_Source -c conf -b mymacro.mac
 
 By default, SToGS_Source loads the setup/SToGS_Source.default file to set up the session. 
 You can change this with the -c conf option where conf is the name of the configuration file you would like
 */
int main(int argc,char** argv)
{
	// check out the command line options. default is to read global.paris unless another file is givn on command line
	G4String gconf = "setup/SToGS_Source.conf", macro; G4bool is_interactive = true;
	//
	for( G4int i = 1; i < argc ; i++) {
		G4String arg = argv[i];
		// G4cout << arg << G4endl;
		if ( arg == "-c" && i < argc - 1 )
			gconf = argv[i+1]; // a different global conf file
		if ( arg == "-b" && i < argc - 1 ) {
			is_interactive = false;
			macro = argv[i+1]; // additional standard G4 macro, ex : definition of gps
		}
        // choose the Random engine
        // CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
	}
    
    SToGS::UserActionManager *user_action_manager = new SToGS::UserActionManager(gconf);
    
    // Construct the default run manager which is necessary
#ifdef G4MULTITHREADED
    G4MTRunManager* theRunManager = new G4MTRunManager();
    theRunManager->SetNumberOfThreads(user_action_manager->GetNbThread());
#else
    G4RunManager* theRunManager = new G4RunManager();
#endif
    G4VUserDetectorConstruction *setup = user_action_manager->GetDetectorConstruction();
    if ( setup == 0x0 ) {
        G4cout << " [ERR] Cannot init setup " << G4endl;
        exit(-1);
    }
    else theRunManager->SetUserInitialization ( setup );
    G4VUserPhysicsList *physics_list = user_action_manager->GetPhysicsList();
    if ( physics_list == 0x0 ) {
        G4cout << " [ERR] Cannot init Physics List " << G4endl;
        exit(-1);
    }
    else theRunManager->SetUserInitialization ( physics_list );
#ifdef G4MULTITHREADED
    theRunManager->SetUserInitialization( user_action_manager );
#else
    theRunManager->SetUserAction( user_action_manager->GetGun() );
    theRunManager->SetUserAction( user_action_manager->GetRunAction() );
    theRunManager->SetUserAction( user_action_manager->GetEventAction() );
    theRunManager->SetUserAction( user_action_manager->GetTrackingAction() );
    theRunManager->SetUserAction( user_action_manager->GetSteppingAction() );
#endif
	
	// Initialize G4 kernel
	theRunManager->Initialize();
	
	// Visualization manager
	G4VisManager* visManager = 0; G4UIsession *session = 0;
	
	if (is_interactive)   // Define UI terminal for interactive mode
	{
		visManager = new G4VisExecutive();
        visManager->SetVerboseLevel(G4VisManager::quiet);
        visManager->Initialize();
		
#ifdef G4UI_USE_TCSH
		session = new G4UIterminal(new G4UItcsh);
#else
		session = new G4UIterminal();
#endif
		session->SessionStart();
		delete session;
	}
	else  // Batch mode
	{
		G4cout << "UI interface is started" << G4endl;
		
		G4String command = "/control/execute ";
		G4String fileName = macro;
		G4UImanager::GetUIpointer()->ApplyCommand(command+fileName);
	}
	
	// job termination
	if ( visManager )
        delete visManager;
    delete theRunManager;
	
	return 0;
}





