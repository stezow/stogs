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


#include "SToGS_UserActionInitialization.hh"
#include "SToGS_G4_TrackInformation.hh"

#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4TrackingManager.hh"

#ifdef G4MULTITHREADED
#include "G4AutoLock.hh"
namespace { G4Mutex buildMutex = G4MUTEX_INITIALIZER; }
#endif

#include <fstream>

void SToGS::TrackingAction::PreUserTrackingAction(const G4Track* atrack)
{
    // if real primary, no parent, add user information
    if ( atrack->GetParentID() == 0 ) {
        if ( atrack->GetUserInformation() ) {
            SToGS::PrimaryTrackInformation *pinfo = new SToGS::PrimaryTrackInformation(atrack->GetTrackID());
            fpTrackingManager->SetUserTrackInformation(pinfo);
        }
    }
    if (theRealAction) {
        theRealAction->PreUserTrackingAction(atrack);
    }
}
void SToGS::TrackingAction::PostUserTrackingAction(const G4Track *atrack)
{
    G4TrackVector *secondaries = fpTrackingManager->GimmeSecondaries();
    if ( secondaries ) {
        size_t nbSecondaries = secondaries->size();
        if ( nbSecondaries > 0 ) {
            SToGS::PrimaryTrackInformation *pinfo = (SToGS::PrimaryTrackInformation *)atrack->GetUserInformation();
            G4int pid = atrack->GetTrackID();
            if ( pinfo ) {
                pid = pinfo->GetPrimaryID();
            }
            for (size_t i = 0; i < nbSecondaries; i++) {
                SToGS::PrimaryTrackInformation *sinfo = new SToGS::PrimaryTrackInformation(pid);
                (*secondaries)[i]->SetUserInformation(sinfo);
            }
        }
    }
    if (theRealAction) {
        theRealAction->PostUserTrackingAction(atrack);
    }
}

SToGS::UserActionInitialization::UserActionInitialization(G4String user_action_opt, G4String which_gene, G4String which_gene_opt)
#if G4VERSION_NUMBER < 1000
    :
#else
    : G4VUserActionInitialization(),
#endif
    fUserActionOption(user_action_opt),
    fWhichGenerator(which_gene,which_gene_opt)
{
}
SToGS::UserActionInitialization::~UserActionInitialization()
{
    for (size_t i = 0; i < fAllUserAction.size(); i++) {
        if ( fAllUserAction[i] ) {
            delete fAllUserAction[i];
            fAllUserAction[i] = 0x0;
        }
    };
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
#include "SToGS_G4_CopClusterSD.hh"

G4VSensitiveDetector *SToGS::UserActionInitialization::GetTrackerSD( G4String name )
{
    G4cout << "[+[SToGS::UserActionInitialization::GetTrackerSD()]] Creating a Tracker SD " << G4endl;
    G4VSensitiveDetector *aSD = 0x0;
    
    G4SDManager *SDman =
        G4SDManager::GetSDMpointer();
    
#ifdef G4MULTITHREADED
    // TO BE CHECKED IF SHOULD BE LIKE THAT !
    aSD = new SToGS::TrackerSD(name);
    if ( aSD )
        SDman->AddNewDetector(aSD);
#else
    aSD = SDman->FindSensitiveDetector(name,false);
    if ( aSD == 0x0 ) {
        aSD = new SToGS::TrackerSD(name);
        if ( aSD )
            SDman->AddNewDetector(aSD);
    }
#endif
    G4cout << "[_[SToGS::UserActionInitialization::GetTrackerSD()]] Creating a Tracker SD " << G4endl;
    return aSD;
}

G4VSensitiveDetector *SToGS::UserActionInitialization::GetCopClusterSD( G4String name )
{
    G4cout << "[+[SToGS::UserActionInitialization::GetCopClusterSD()]] Creating a CopCluster SD " << G4endl;
    G4VSensitiveDetector *aSD = 0x0;
    
    G4SDManager *SDman =
        G4SDManager::GetSDMpointer();
    
#ifdef G4MULTITHREADED
    // TO BE CHECKED IF SHOULD BE LIKE THAT !
    aSD = new SToGS::CopClusterSD(name);
    if ( aSD )
        SDman->AddNewDetector(aSD);
#else
    aSD = SDman->FindSensitiveDetector(name,false);
    if ( aSD == 0x0 ) {
        aSD = new SToGS::CopClusterSD(name);
        if ( aSD )
            SDman->AddNewDetector(aSD);
    }
#endif
    G4cout << "[_[SToGS::UserActionInitialization::GetCopClusterSD()]] Creating a CopCluster SD " << G4endl;
    return aSD;
}

void SToGS::UserActionInitialization::BuildAndRegister ( AllActions *actions ) const
{
#ifdef G4MULTITHREADED
    G4AutoLock lock(&buildMutex);
#endif
    G4cout << " ------ INF ------ from SToGS::UserActionInitialization::BuildAndRegister() " << G4endl;
#if G4VERSION_NUMBER < 1000
    G4cout << " *** ERROR *** SToGS::UserActionInitialization::BuildAndRegister() should never be called by this version of Geant4 " << G4endl;
#else
    SetUserAction( GetGun(fWhichGenerator.first,fWhichGenerator.second) );
    //
    fAllUserAction.push_back(actions);
    //
    SetUserAction( new SToGS::RunAction(actions) ) ;
    SetUserAction( new SToGS::EventAction(actions) ) ;
    SetUserAction( new SToGS::TrackingAction(actions) ) ;
    SetUserAction( new SToGS::SteppingAction(actions) ) ;
#endif
    G4cout << " ------ END ------  from SToGS::UserActionInitialization::BuildAndRegister() " << G4endl;
}

void SToGS::UserActionInitialization::Register ( AllActions *actions ) const
{
#ifdef G4MULTITHREADED
    G4AutoLock lock(&buildMutex);
#endif    
    G4cout << " ------ INF ------ from SToGS::UserActionInitialization::Register() " << G4endl;
#if G4VERSION_NUMBER < 1000
    G4cout << " *** ERROR *** SToGS::AllInOneUserActionInitialization::Build() should never be called by this version of Geant4 " << G4endl;
#else
    fAllUserAction.push_back(actions);
#endif
    G4cout << " ------ END ------  from SToGS::UserActionInitialization::Register() " << G4endl;
}


