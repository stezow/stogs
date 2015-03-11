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
//----------------------------------------------------------------------------------


#include "SToGS_UserActionManager.hh"
#include "SToGSConfig.hh"

#include "G4ios.hh"
#include <fstream>

SToGS::UserActionManager::UserActionManager(G4String filename) :
    fImplementation(0x0),
    fWhichGenerator("GPS","G4Macros/GPSPointLike.mac"),
    fWhichGeometry("factory","DetectorFactory/Generics/Target"),
    fWhichPhysics("stogs_m","general0;emstandard_opt0;"),
    fWhichActionManager("PrintOut","run;event;track;step"),
    fNbThreads(2)
{
    G4cout << G4endl << " ------ INF ------ from UserActionManager::UserActionManager with " << filename << G4endl;
    
	// open the ascii file
    std::ifstream file(filename);
	if ( file.is_open() == false ) {
		G4cout << " ** SToGS WARNING ** Cannot open file, Default parameters to be used "<< G4endl;
    }
	else {
        std::string aline; getline(file,aline);
        while ( file.good() && !file.eof() ) {
            
            if ( aline[0] == '#' ){
                getline(file,aline);
                continue;
            } // this line is a comment
            
            std::string key, which, option; std::istringstream decode(aline); decode >> key;
            
            if ( key == "actions:" ) {
                decode >> which >> option;
                fWhichActionManager = std::pair<G4String, G4String> (which,option);
            }
            if ( key == "setup:" ) {
                decode >> which >> option;
                fWhichGeometry = std::pair<G4String, G4String> (which,option);
            }
            if ( key == "physics:" ) {
                decode >> which >> option;
                fWhichPhysics = std::pair<G4String, G4String> (which,option);
            }
            if ( key == "generator:" ) {
                decode >> which >> option;
                fWhichGenerator = std::pair<G4String, G4String> (which,option);
            }
            if ( key == "nbthread:" ) {
                decode >> fNbThreads;
            }
            getline(file,aline);
        }
        file.close();
    }
    

    fImplementation = ProvideUserActionInitialization();
    if ( fImplementation == 0x0 ) {
        G4cout << " *** SToGS ERROR *** Action Manager is not known *** SToGS ERROR *** " << G4endl;
    }
    else {
        fImplementation->SetWhichGenerator(fWhichGenerator.first,fWhichGenerator.second);
    }
    
    G4cout << " ActionManager has been initiated with " << G4endl;
	G4cout << "   + generator        " << fWhichGenerator.first << " " << fWhichGenerator.second << G4endl;
	G4cout << "   + physics list     " <<   fWhichPhysics.first << " " << fWhichPhysics.second   << G4endl;
	G4cout << "   + geometry used    " <<  fWhichGeometry.first << " " << fWhichGeometry.second  << G4endl;
	G4cout << "   + action manager   " << fWhichActionManager.first << " " << fWhichActionManager.second << G4endl;
    
	G4cout << " ------ END ------ from UserActionManager::UserActionManager " << G4endl;
}

#include "SToGS_PrintOut.hh"
#include "SToGS_Ascii.hh"
#ifdef HAS_STOGS_ROOT_EVENTS
#include "SToGS_BaseROOTEventsActions.hh"
#endif

SToGS::UserActionInitialization *SToGS::UserActionManager::ProvideUserActionInitialization()
{
    if ( fWhichActionManager.first == "printout" ) {
        fImplementation = new SToGS::PrintOut(fWhichActionManager.second);
    }
    if ( fWhichActionManager.first == "ascii" ) {
        fImplementation = new SToGS::Ascii(fWhichActionManager.second);
    }
#ifdef HAS_STOGS_ROOT_EVENTS
    if ( fWhichActionManager.first == "stogstree" ) {
        fImplementation = new SToGS::BaseROOTEvents(fWhichActionManager.second);
    }
#endif
    // based on My plugins defined in SToGSConfig, it build the user actio
#ifdef HAS_MYACT
    if ( fWhichActionManager.first.contains(MYACT_) ) {
       fImplementation = new MYACT_CLASSTYPE(fWhichActionManager.second);
    }
#endif
    return fImplementation;
}

#include "SToGS_LoadFromDetectorFactory.hh"

G4VUserDetectorConstruction *SToGS::UserActionManager::GetDetectorConstruction() const
{
    G4VUserDetectorConstruction *detectorconstruction = 0x0;
    if ( fWhichGeometry.first == "factory" ) {
        detectorconstruction = new SToGS::LoadFromDetectorFactory(fWhichGeometry.second);
    }
    
    return detectorconstruction;
}

#include "SToGS_ModularPhysicsList.hh"

G4VUserPhysicsList *SToGS::UserActionManager::GetPhysicsList() const
{
    G4VUserPhysicsList *physics_list = 0x0;
    
    if ( fWhichPhysics.first == "stogs_m" ) {
        physics_list = new SToGS::ModularPhysicsList(fWhichPhysics.second);
    }
    
    return physics_list;
}

