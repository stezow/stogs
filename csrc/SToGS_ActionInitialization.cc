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
#include "G4ios.hh"
#include "G4SDManager.hh"

#include <fstream>

/*
#include "SToGS_UserActionManager.hh"
#include "SToGS_UserActionManager.hh"
#include "SToGS_UserActionManager.hh"
#include "SToGS_UserActionManager.hh"
 */


SToGS::UserActionManager::UserActionManager(G4String filename) :
#if G4VERSION_NUMBER < 1000
#else
    G4VUserUserActionManager(),
#endif
    fWhichGeometry(),
    fWhichPhysics(),
    fWhichGenerator("GPS","G4Macros/GPS_Cs137.mac"),
    fWhichActionManager()
{
    G4cout << G4endl << " ------ INFO ------ from UserActionManager::UserActionManager with " << filename << G4endl << G4endl;
		
	// open the ascii file
    std::ifstream file(filename);
	if ( file.is_open() == false ) {
		G4cout << " ** WARNING ** cannot open file  (Default parameters are used) "<< G4endl;
		G4cout << " ** WARNING ** defaults are to be used " << G4endl;
    }
	else {
        std::string aline; getline(file,aline);
        while ( file.good() && !file.eof() ) {
            
            if ( aline[0] == '#' ){
                getline(file,aline);
                continue;
            } // this line is a comment
            
            std::string key, which, option; std::istringstream decode(aline); decode >> key >> which >> option;

            if ( decode.good() > 0 ) {
                if ( key == "analysis:" ) {
                    fWhichActionManager = std::pair<G4String, G4String> (which,option);
                }
                if ( key == "geom:" ) {
                    fWhichGeometry = std::pair<G4String, G4String> (which,option);
                }
                if ( key == "physics:" ) {
                    fWhichPhysics = std::pair<G4String, G4String> (which,option);
                }
                if ( key == "gene:" ) {
                    fWhichGenerator = std::pair<G4String, G4String> (which,option);
                }
            }
            getline(file,aline);
        }
        file.close();
    }

	G4cout << " ActionManager has been initiated with " << G4endl;
	G4cout << "   + generator      " << fWhichGenerator.first << " " << fWhichGenerator.second << G4endl;
	G4cout << "   + physics list   " <<   fWhichPhysics.first << " " << fWhichPhysics.second   << G4endl;
	G4cout << "   + geometry used  " <<  fWhichGeometry.first << " " << fWhichGeometry.second  << G4endl;
	G4cout << "   + Action manager " << fWhichActionManager.first << " " << fWhichActionManager.second << G4endl;
    
	G4cout << G4endl << " ------ END ------ from UserActionManager::UserActionManager " << G4endl << G4endl;
}

SToGS::UserActionManager::~UserActionManager()
{
    ;
}

// Primary Generators
#include "SToGS_G4_GPSPrimaryGeneratorAction.hh"

G4VUserPrimaryGeneratorAction *SToGS::UserActionManager::GetGun(G4String which, G4String opt)
{
    G4VUserPrimaryGeneratorAction *gun = 0x0;
    if ( which == "GPS" )
        gun = new SToGS::GPSPrimaryGeneratorAction(opt);
    
    return gun;
}

#include "SToGS_G4_TrackerSD.hh"
#include "SToGS_G4_CaloSD.hh"

G4VSensitiveDetector *SToGS::UserActionManager::GetTrackerSD( G4String name )
{
    G4cout << "[+[SToGS::UserActionManager::GetTrackerSD()]] Creating a tracker SD " << G4endl;
    G4VSensitiveDetector *aSD = 0x0;
    
    G4SDManager *SDman =
        G4SDManager::GetSDMpointer();
    
#ifdef G4MULTITHREADED
    // TO BE CHECKED IF SHOULD BE LIKE THAT !
    aSD = new SToGS::TrackerSD(name);
    if ( aSD )
        SDman->AddNewDetector(aSD);
#else
    aSD = SDman->FindSensitiveDetector(name);
    if ( aSD == 0x0 ) {
        aSD = new SToGS::TrackerSD(name);
        if ( aSD )
            SDman->AddNewDetector(aSD);
    }
#endif
    G4cout << "[_[SToGS::UserActionManager::GetTrackerSD()]] Creating a tracker SD " << G4endl;
    return aSD;
}

G4VSensitiveDetector *SToGS::UserActionManager::GetCaloSD( G4String name )
{
    G4cout << "[+[SToGS::UserActionManager::GetCaloSD()]] Creating a tracker SD " << G4endl;
    G4VSensitiveDetector *aSD = 0x0;
    
    G4SDManager *SDman =
        G4SDManager::GetSDMpointer();
    
#ifdef G4MULTITHREADED
    // TO BE CHECKED IF SHOULD BE LIKE THAT !
    aSD = new SToGS::CaloSD(name);
    if ( aSD )
        SDman->AddNewDetector(aSD);
#else
    aSD = SDman->FindSensitiveDetector(name);
    if ( aSD == 0x0 ) {
        aSD = new SToGS::CaloSD(name);
        if ( aSD )
            SDman->AddNewDetector(aSD);
    }
#endif
    G4cout << "[_[SToGS::UserActionManager::GetCaloSD()]] Creating a calo SD " << G4endl;
    return aSD;
}

void SToGS::UserActionManager::BuildForMaster() const
{
    //   SetUserAction(new MyRunAction);
    fImplementation->BuildForMaster();
}

void SToGS::UserActionManager::Build() const
{
    fImplementation->Build();

   // SetUserAction( new SToGS::GPSPrimaryGeneratorAction () );
    /*
    SetUserAction(new MyRunAction);
    SetUserAction(new MySteppingAction);
     */
}
