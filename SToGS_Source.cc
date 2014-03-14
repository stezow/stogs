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

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "ParisBasicRunAction.hh"
#include "ParisBasicEventAction.hh"
#include "ParisBasicSteppingAction.hh"
#include "ParisBasicTrackingAction.hh"

#include <fstream>

using namespace std;

/** \file SToGS_Source Main program for running test / characterization / sources like experiment */

//! To select the different modules from an ascii file
/*!
 */
void InitGlobal(const char *filename = "setup/global.stogs");

//! This is the main program of the simulation.
/*!
 \code
 #
 # This file is used to configure the Paris program. It is read at the beginning
 # of the program to select
 #
 # Analysis (output manager)
 #	 0: ParisPrintOut
 #	 1: ParisBaseAscii
 #	 2: ParisAscii
 #
 0 setup/ascii.ana
 # Detector geometry
 #	 0: ParisShellDetectorConstruction
 #	 1: ParisSegmentedDetectorConstruction
 #
 0 setup/shells.geo
 #
 # physics list
 #	 0: ParisStandardEMPhysicsList
 #	 1: ParisLowEnergyEMPhysicsList
 #	 2: ParisPenelopeEMPhysicsList
 #
 0
 # generator
 #	 0: ParisBasicPrimaryGeneratorAction
 #	 1: ParisUniformPrimaryGeneratorAction
 #
 O setup/basic.gene
 #
 \endcode
 
 Without arguments, the program starts a G4 interactive session. It is also possible to run it
 in batch mode by adding on the command line the name of the G4 macro you would like to execute. Ex:
 
 \c Paris -c conf -b mymacro.mac
 
 By default, Paris loads the setup/global.paris file to set up the session. You can change this
 with the -c conf option where conf is the name of the configuration file you would like
 */
int main(int argc,char** argv)
{
	// check out the command line options. default is to read global.paris unless another file is givn on command line
	G4String gconf = "setup/global.stogs", macro; G4bool is_interactive = true;
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
	
	// Construct the default run manager which is necessary
	G4RunManager *theRunManager = new G4RunManager();
    
	// Init from setup file
	InitGlobal(gconf.data());
	//
	if ( G4RunManager::GetRunManager()->GetUserDetectorConstruction() == 0x0 ||
        G4RunManager::GetRunManager()->GetUserPhysicsList() == 0x0 ||
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction() == 0x0 ) {
		exit (-1);
	}
	// now set the others user actions
	theRunManager->SetUserAction( new ParisBasicRunAction() );
	theRunManager->SetUserAction( new ParisBasicEventAction() );
	theRunManager->SetUserAction( new ParisBasicTrackingAction() );
	theRunManager->SetUserAction( new ParisBasicSteppingAction() );
	
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

// function to be modified to add new geom/physics/gene
//
// geometries
#include "ParisShellDetectorConstruction.hh"
#include "ParisSegmentedDetectorConstruction.hh"
#include "ParisSegmentedShellDetectorConstruction.hh"
#include "ParisLoadDetectorConstruction.hh"
#include "ParisSDCRotated.hh"
#include "ParisSDCRotated_CsI.hh"
#include "ParisSemiSphe2LayersConstruction.hh"
#include "Paris3by3CubicConstruction.hh"

#ifdef HAS_MYDETECTOR
#include "MyDetectorConstruction.hh"
#endif

//
void SetDetectorConstruction(G4int detectorid, const G4String &configdet)
{
	G4VUserDetectorConstruction *detector = 0x0;
	switch( detectorid ){
		case -1:
			detector = new ParisLoadDetectorConstruction(configdet);
			G4cout << " You are working with detector # -1 (from factory using file " << configdet << ")"  << G4endl;
			break;
		case 0:
			detector = new ParisShellDetectorConstruction(configdet);
			G4cout << " You are working with detector # 0 (two perfect shells from " << configdet << ")"  << G4endl;
			break;
		case 1:
			detector = new ParisSegmentedDetectorConstruction(configdet);
			G4cout << " You are working with detector #1 (configuration from " << configdet << ")"<< G4endl;
			break;
		case 2:
			detector = new ParisSDCRotated(configdet);
			G4cout << " You are working with detector # 2 (configuration from " << configdet << ")"  << G4endl;
			break;
		case 3:
			detector = new ParisSDCRotated_CsI(configdet);
			G4cout << " You are working with detector #3 (configuration from " << configdet << ")"<< G4endl;
			break;
		case 4:
			detector = new ParisSegmentedShellDetectorConstruction(configdet);
			G4cout << " You are working with detector #4 (configuration from " << configdet << ")"<< G4endl;
			break;
		case 5:
			detector = new ParisSemiSphe2LayersConstruction(configdet);
			G4cout << " You are working with detector #5 (configuration from " << configdet << ")"<< G4endl;
			break;
		case 6:
			detector = new Paris3by3CubicConstruction(); G4cout << "You are working with detector #6 "<< G4endl;
			break;
#ifdef HAS_MYDETECTOR
		case 100:
			detector = new MyDetectorConstruction(configdet); G4cout << " You are working with detector #100 : MyDetectorConstruction "<< G4endl;
			break;
#endif
		default:
			G4cout << "  ==> ** ERROR IN ** SetDetectorConstruction # " << detectorid << " " << configdet << " is not defined " << G4endl;
			break;
	}
	// register detector
	if ( detector )
		G4RunManager::GetRunManager()->SetUserInitialization(detector);
	else
		exit(-1);
}
// physics list
#include "ParisPhysicsList.hh"
//
void SetPhysicsList(G4int physicsid, const G4String &opt)
{
	G4VUserPhysicsList *physicsList = 0x0;
	switch( physicsid ){
		case -1:
			physicsList = new ParisPhysicsList(opt) ;
			G4cout << " You are working with physics list # -1 (customisable physics list) " << G4endl;
			break;
		case 0:
			physicsList = new ParisPhysicsList("general0;ParisStandardEM") ;
			G4cout << " You are working with physics list # 0 (Standard EM) " << G4endl;
			break;
		case 1:
			physicsList = new ParisPhysicsList("general0;ParisLowEnergyEM") ;
			G4cout << " You are working with physics list # 1 (Low energy EM) " << G4endl;
			break;
		case 2:
			physicsList = new ParisPhysicsList("general0;ParisPenelopeEM");
			G4cout << " You are working with physics list # 2 (Penelope EM) " << G4endl;
			break;
            /*		case 3:
             physicsList = new ParisIonPhysList() ;
             G4cout << " You are working with physics list # 3 (Ion + EM) " << G4endl;
             break;
             case 4:
             physicsList = new ParisRadIonPhysList() ;
             G4cout << " You are working with physics list # 4 (Rad Ion + EM) " << G4endl;
             break;	*/
		default:
			G4cout << " ** WARNING ** physics list # " << physicsid << " is not defined " << G4endl;
			break;
	}
	// register physics list
	if ( physicsList )
		G4RunManager::GetRunManager()->SetUserInitialization(physicsList);
	else
		exit(-1);
}
// generators
#include "ParisBasicPrimaryGeneratorAction.hh"
#include "ParisUniformPrimaryGeneratorAction.hh"
#include "ParisRandomCascadeGeneratorAction.hh"
#include "ParisGDRPrimaryGeneratorAction.hh"
//#include "ParisCascadeGenerator.hh"
//#include "ParisCascadeGenerator_ttree.hh"
//#include "ParisIonPrimaryGeneratorAction.hh"
#include "GPSPrimaryGeneratorAction.hh"

//
void SetGenerator(G4int generatorid, G4String &configgen)
{
	G4VUserPrimaryGeneratorAction *generator = 0x0;
	switch( generatorid ){
		case -1:
			generator = new GPSPrimaryGeneratorAction(configgen);
			G4cout << " You are working with the basic generator # -1 ( GPS generator " << ")" << G4endl;
			break;
		case 0:
			generator = new ParisBasicPrimaryGeneratorAction(configgen);
			G4cout << " You are working with the basic generator # 0 (a cascade reads from " << configgen << ")" << G4endl;
			break;
		case 1:
			generator = new ParisUniformPrimaryGeneratorAction(configgen);
			G4cout << " You are working with the uniform generator # 1 (reads from " << configgen << ")" << G4endl;
			break;
		case 2:
			generator = new ParisGDRPrimaryGeneratorAction(configgen);
			G4cout << " You are working with the uniform generator # 2 (reads from " << configgen << ")" << G4endl;
			break;
            /*		case 3:
             generator = new ParisRandomCascade(configgen);
             G4cout << " You are working with the uniform generator # 3 (reads from " << configgen << ")" << G4endl;
             break;
             case 4:
             generator = new ParisCascadeGenerator(configgen);
             G4cout << " You are working with the cascade generator # 4 (reads from " << configgen << ")" << G4endl;
             break;
             case 5:
             generator = new ParisCascadeGenerator_ttree(configgen);
             G4cout << " You are working with the cascade generator # 5 (reads from " << configgen << ")" << G4endl;
             break; */
            /*		case 6:
             generator = new ParisIonPrimaryGeneratorAction(configgen);
             G4cout << " You are working with the Ion generator # 6 (reads from " << configgen << ")" << G4endl;
             break; */
		default:
			G4cout << " ** WARNING ** generator # " << generatorid << " is not defined " << G4endl;
			break;
	}
	// register physics list
	if ( generator )
		G4RunManager::GetRunManager()->SetUserAction(generator);
	else
		exit(-1);
}
// analysis manager
#include "ParisPrintOut.hh"
#include "ParisAscii.hh"
#include "ParisROOT.hh"
#include "ParisOpticalROOT.hh"

ParisOutputManager *SetOutputManager(G4int outputid, G4String &configana)
{
	// Construction the Paris output manager
	switch( outputid ){
		case 0:
			ParisOutputManager::SetTheOutputManager( new ParisPrintOut(configana) );
			G4cout << " Analysis consists in printing out information on standard output with init file "
            << configana << G4endl;
			break;
		case 1:
			ParisOutputManager::SetTheOutputManager( new ParisBaseAscii() );
			G4cout << " Analysis consists in printing informations in an ascii file " << G4endl;
			break;
		case 2:
			ParisOutputManager::SetTheOutputManager( new ParisAscii(configana) );
			G4cout << " Analysis consists in printing informations in an ascii file. Initialization from "
            << configana << G4endl;
			break;
		case 3:
			ParisOutputManager::SetTheOutputManager( new ParisROOT(configana) );
			G4cout << " Analysis consists in printing informations in a ROOT file. Initialization from "
            << configana << G4endl;
			break;
		case 10:
			ParisOutputManager::SetTheOutputManager( new ParisOpticalROOT() );
			G4cout << " Analysis consists in basic ROOT + information on optical photons "
			<< configana << G4endl;
			break;
		default:
			ParisOutputManager::SetTheOutputManager( new ParisPrintOut() );
			G4cout << " ** WARNING ** analysis # " << outputid << " is not defined " << G4endl;
			G4cout << " ==> DEFAULT INSTEAD: Analysis consists in printing out information on standard output " << G4endl;
			break;
	}
    
	return ParisOutputManager::GetTheOutputManager();
}

void InitGlobal(const char *filename)
{
	G4cout << G4endl << " ------ INFO ------ from InitGlobal " << G4endl << G4endl;
	
	// default parameters
	G4int ana,det,phy,gen; G4String configana, configdet, configphy, configgen;
    
	ana = det = phy = gen = 0; // default parameters
	
	// open the ascii file
	ifstream file(filename);
	if ( file.is_open() == false ) {
		G4cout << " ** WARNING ** cannot open file " << filename << " (Default parameters are used) "<< G4endl; return; }
	
	const G4int MAXWIDTH = 1000; G4int tmpi, which; char keys[30], tmps[MAXWIDTH]; char aline[MAXWIDTH]; file.getline(aline,MAXWIDTH);
	which = 0;
	while ( file.good() ) {
        
		if ( aline[0] == '#' ){
            file.getline(aline,MAXWIDTH);
            continue;
        } // this line is a comment
		
        ::memset(tmps, 0, MAXWIDTH); G4int format = sscanf(aline,"%s %d %s",keys,&tmpi,tmps);
		if ( format > 0 ) {
			G4String key = keys;
			if ( key == "analysis:" ) {
				ana = tmpi; configana = tmps; which++;
			}
			if ( key == "geom:" ) {
				det = tmpi; configdet = tmps; which++;
			}
			if ( key == "physics:" ) {
				phy = tmpi; configphy = tmps; which++;
			}
			if ( key == "gene:" ) {
				gen = tmpi; configgen = tmps; which++;
				
			}
		}
		file.getline(aline,MAXWIDTH);
	}
	G4cout << " (analysisID,detectorID,physicsID,generatorID) read from " << filename
    << " = (" << ana << "," << det << "," << phy << "," << gen << ")" << G4endl;
	
	SetOutputManager(ana, configana);
	SetDetectorConstruction(det, configdet);
	SetPhysicsList(phy, configphy);
	SetGenerator(gen, configgen);
    
	G4cout << G4endl << " ------ END ------ from InitGlobal " << G4endl << G4endl; 
	
	file.close();
}






