
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

// G4 includes
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "SToGS_DetectorFactory.hh"
#include "SToGS_ModularPhysicsList.hh"
#include "SToGS_PrintOut.hh"
#include "SToGS_LoadFromDetectorFactory.hh"

// TMP to test MIGRATION
//#include "SToGS_ShellDetectorConstruction.hh"


/*! DetectorBuilder helps in building detectors/setup using the Detector Factories
 
 Construction is done through the Detector Factory using the
 */
int main(int argc, char** argv)
{
    // file to be given on the command line. If not, set a demo
    G4int what_detector = -1; G4String filedfb = "", macro_visu = "G4Macros/vis.mac";
    for( G4int i = 1; i < argc ; i++) {
		G4String arg = argv[i];
		if ( arg == "-dfb" && i < argc - 1 ) {
			filedfb = argv[i+1];
        }
        if ( arg == "-myd" ) {
            what_detector = 0;
        }
        if ( arg == "-vis" ) {
			macro_visu = argv[i+1];
        }
	}
    if ( filedfb == "" ) {
        filedfb = "default.dfb";
    }
    
    //Pure SToGS related
    // Make sure an output manager is set and the main factory is built
    SToGS::DetectorFactory::theMainFactory()->MakeStore();
    // simple printout manager enough at the level of detector construction
    // it shows also how it can be directy used without having to go through ActionManager
    SToGS::UserActionInitialization *user_action = new SToGS::PrintOut("run;event;track;step","GPS","G4Macros/GPSPointLike.mac");
    
	// Construct the default run manager which is necessary
    G4RunManager* theRunManager = new G4RunManager;
    //
    switch ( what_detector ) {
        case 0:
            break;
        default:
            theRunManager->SetUserInitialization ( new SToGS::BuildFromDetectorFactory(filedfb) );
            break;
    }
    theRunManager->SetUserInitialization ( new SToGS::ModularPhysicsList() );
#if G4VERSION_NUMBER < 1000
    theRunManager->SetUserAction( user_action->GetGun() );
    theRunManager->SetUserAction( user_action->GetRunAction() );
    theRunManager->SetUserAction( user_action->GetEventAction() );
    theRunManager->SetUserAction( user_action->GetTrackingAction() );
    theRunManager->SetUserAction( user_action->GetSteppingAction() );
#else
    theRunManager->SetUserInitialization( user_action );
#endif
    
    // Initialize G4 kernel
    theRunManager->Initialize();
    
    // Get the pointer to the User Interface manager
    G4UImanager * UImanager = G4UImanager::GetUIpointer();
    
    // And the Visualization manager
    G4VisManager* visManager = new G4VisExecutive();
    //
    if ( visManager ) {
        visManager->SetVerboseLevel(G4VisManager::quiet);
        visManager->Initialize();
    }
    
    G4UIsession *session = 0x0;
#ifdef G4UI_USE
    G4UIExecutive * ui = new G4UIExecutive(argc,argv);
    UImanager->ApplyCommand("/vis/ASCIITree/verbose 11");
    UImanager->ApplyCommand("/vis/drawTree");
    UImanager->ApplyCommand("/geometry/test/grid_test");
#ifdef G4VIS_USE
    G4String cmd = "/control/execute ";
    cmd += macro_visu;
    UImanager->ApplyCommand(cmd.data());
#endif
    ui->SessionStart();
    delete ui;
#endif
    // job termination
    delete session;
    
    if ( visManager ) delete visManager;
    delete theRunManager;
	
    return 0;
}









