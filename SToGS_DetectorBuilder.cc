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
#include "G4Version.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#include "SToGS_ActionInitialization.hh"
#else
#include "G4RunManager.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "SToGS_DetectorFactory.hh"

/*! DetectorBuilder helps in building detectors/setup using the Detector Factories
*/
int main(int argc, char** argv)
{
    // file to be given on the command line. If not, set a demo
    G4String filedfb = "";
    for( G4int i = 1; i < argc ; i++) {
		G4String arg = argv[i];
		if ( arg == "-b" && i < argc - 1 )
			filedfb = argv[i+1];
	}
    if ( filedfb == "" ) {
        filedfb = "default.dfb";
    }
    G4cout << G4VERSION_NUMBER << G4endl;
    
    //! Make sure an output manager is set and the main factory is built
    // ParisOutputManager::SetTheOutputManager( new ParisPrintOut("") );
    SToGS::DetectorFactory::theMainFactory()->MakeStore();
	
	// Construct the default run manager which is necessary
#ifdef G4MULTITHREADED
    G4MTRunManager* theRunManager = new G4MTRunManager;
    theRunManager->SetNumberOfThreads(2);
#else
    G4RunManager* theRunManager = new G4RunManager;
#endif
    
#ifdef G4MULTITHREADED
    theRunManager->SetUserInitialization( new SToGS::ActionInitialization() );
#else
#endif
    
	// mandatory initialization G4 classes provided by the user
	// detector construction
	// physics  list
	// initialisation of the generator

	// Init from setup file
    //G4RunManager::GetRunManager()->SetUserInitialization( new ParisLoadDetectorConstruction(filedfb) );
    //G4RunManager::GetRunManager()->SetUserInitialization( new ParisPhysicsList("general0;ParisStandardEM") );
    // G4RunManager::GetRunManager()->SetUserAction( new GPSPrimaryGeneratorAction() );

	//
	if ( G4RunManager::GetRunManager()->GetUserDetectorConstruction() == 0x0 || 
		  G4RunManager::GetRunManager()->GetUserPhysicsList() == 0x0 || 
		  G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction() == 0x0 ) {
		exit (-1);
	}

	// Initialize G4 kernel
	theRunManager->Initialize();  
	
	// Visualization manager
	G4VisManager* visManager = new G4VisExecutive();
    //
    if ( visManager ) {
        visManager->SetVerboseLevel(G4VisManager::quiet);
        visManager->Initialize();
    }
		
    G4UIsession *session = 0;
    //
#ifdef G4UI_USE_TCSH
    session = new G4UIterminal(new G4UItcsh);
#else
    session = new G4UIterminal();
#endif
    // a set of commands to check the geometry
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/ASCIITree/verbose 11");
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/drawTree");
    G4UImanager::GetUIpointer()->ApplyCommand("/geometry/test/grid_test");
    
    G4cout << "\n A detector/setup has been built using " << filedfb << "\n";
    G4cout << " To view it through openGL:\n    /control/execute macros/visGL.mac \n";
    G4cout << " Try also /DetectorFactory commands \n" << G4endl;

    session->SessionStart();
    delete session;

	// job termination
	if ( visManager )
        delete visManager;
    delete theRunManager;
	
	return 0; 
}





