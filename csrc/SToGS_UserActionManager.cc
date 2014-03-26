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

#include <fstream>

/*
#include "SToGS_UserActionInitialization.hh"
#include "SToGS_UserActionInitialization.hh"
#include "SToGS_UserActionInitialization.hh"
#include "SToGS_UserActionInitialization.hh"
 */

// includes here all possible actions
#include "SToGS_PrintOut.hh"

SToGS::UserActionManager::UserActionManager(G4String filename) :
    UserActionInitialization(),
    fImplementation(0x0),
    fWhichGeometry(),
    fWhichPhysics(),
    fWhichActionManager("PrintOut","")
{
    G4cout << G4endl << " ------ INFO ------ from UserActionManager::UserActionManager with " << filename << G4endl << G4endl;
    
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
            
            std::string key, which, option; std::istringstream decode(aline); decode >> key >> which >> option;
            
      //      if ( decode.good() ) {
                if ( key == "actions:" ) {
                    fWhichActionManager = std::pair<G4String, G4String> (which,option);
                }
                if ( key == "setup:" ) {
                    fWhichGeometry = std::pair<G4String, G4String> (which,option);
                }
                if ( key == "physics:" ) {
                    fWhichPhysics = std::pair<G4String, G4String> (which,option);
                }
                if ( key == "generator:" ) {
                    fWhichGenerator = std::pair<G4String, G4String> (which,option);
       //         }
            }
            getline(file,aline);
        }
        file.close();
    }
    

    fImplementation = GetUserActionInitialization();
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
    
	G4cout << G4endl << " ------ END ------ from UserActionManager::UserActionManager " << G4endl << G4endl;
}

SToGS::UserActionInitialization *SToGS::UserActionManager::GetUserActionInitialization()
{
    if ( fWhichActionManager.first == "PrintOut" ) {
        fImplementation = new PrintOut(fWhichActionManager.second);
    }
    return fImplementation;
}
G4UserRunAction *SToGS::UserActionManager::GetRunAction() const
{
    return fImplementation->GetRunAction();
}
G4UserEventAction *SToGS::UserActionManager::GetEventAction() const
{
    return fImplementation->GetEventAction();
}
G4UserTrackingAction *SToGS::UserActionManager::GetTrackingAction() const
{
    return fImplementation->GetTrackingAction();
}
G4UserSteppingAction *SToGS::UserActionManager::GetSteppingAction() const
{
    return fImplementation->GetSteppingAction();
}
void SToGS::UserActionManager::BuildForMaster() const
{
    fImplementation->BuildForMaster();
}

void SToGS::UserActionManager::Build() const
{
    fImplementation->Build();
}

