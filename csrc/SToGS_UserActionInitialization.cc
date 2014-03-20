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
#include "SToGS_UserActionInitialization.hh"
#include "SToGS_UserActionInitialization.hh"
#include "SToGS_UserActionInitialization.hh"
#include "SToGS_UserActionInitialization.hh"
 */


SToGS::UserActionInitialization::UserActionInitialization()
#if G4VERSION_NUMBER < 1000
    :
#else
    : G4VUserActionInitialization(),
#endif
    fWhichGenerator("-","-")
{
}
SToGS::UserActionInitialization::~UserActionInitialization()
{
    ;
}

std::pair < G4String, G4String > SToGS::UserActionInitialization::SetWhichGenerator(G4String which_gene, G4String option)
{
    std::pair < G4String, G4String > current = fWhichGenerator, new_one(which_gene,option);
    
    fWhichGenerator = new_one;
    return current;
}


// Primary Generators
#include "SToGS_G4_GPSPrimaryGeneratorAction.hh"

G4VUserPrimaryGeneratorAction *SToGS::UserActionInitialization::GetGun(G4String which, G4String opt) const
{
    G4VUserPrimaryGeneratorAction *gun = 0x0;
    if ( which == "GPS" )
        gun = new SToGS::GPSPrimaryGeneratorAction(opt);
    
    return gun;
}
G4VUserPrimaryGeneratorAction *SToGS::UserActionInitialization::GetGun() const
{
    return GetGun(fWhichGenerator.first,fWhichGenerator.second);
}

#include "SToGS_G4_TrackerSD.hh"
#include "SToGS_G4_CaloSD.hh"

G4VSensitiveDetector *SToGS::UserActionInitialization::GetTrackerSD( G4String name )
{
    G4cout << "[+[SToGS::UserActionInitialization::GetTrackerSD()]] Creating a tracker SD " << G4endl;
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
    G4cout << "[_[SToGS::UserActionInitialization::GetTrackerSD()]] Creating a tracker SD " << G4endl;
    return aSD;
}

G4VSensitiveDetector *SToGS::UserActionInitialization::GetCaloSD( G4String name )
{
    G4cout << "[+[SToGS::UserActionInitialization::GetCaloSD()]] Creating a tracker SD " << G4endl;
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
    G4cout << "[_[SToGS::UserActionInitialization::GetCaloSD()]] Creating a calo SD " << G4endl;
    return aSD;
}


