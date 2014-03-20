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

#include "SToGS_PrintOut.hh"
#include "SToGS_G4_TrackerSD.hh"
#include "SToGS_G4_CaloSD.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

SToGS::PrintOutRun::PrintOutRun() :
    G4Run(),
    colltrackerID(-1),
    collcaloID(-1)
{
   	G4SDManager* SDman = G4SDManager::GetSDMpointer();
    if ( SDman ) {
        // keep the collection ID associated to this collection
        colltrackerID = SDman->GetCollectionID("/SToGS/Track");
        
        // keep the collection ID associated to this collection
        collcaloID = SDman->GetCollectionID("/SToGS/Calo");
    }
}
SToGS::PrintOutRun::~PrintOutRun()
{
    ;
}
void SToGS::PrintOutRun::RecordEvent(const G4Event* evt)
{
    SToGS::SingleHitsCollection *THC = 0x0; SToGS::CaloHitsCollection *CHC = 0x0;
	
	if( colltrackerID < 0 || collcaloID < 0 )
		return;
	
	G4HCofThisEvent *HCE = evt->GetHCofThisEvent();
    
	if(HCE)
	{
        THC = (SToGS::SingleHitsCollection *)(HCE->GetHC(colltrackerID));
		CHC = (SToGS::CaloHitsCollection*)(HCE->GetHC(collcaloID));
  	}
    
	if(THC)
	{
		int n_hit = THC->entries();
		G4cout  << "     " << n_hit
                << " hits are stored in tracker collection" << G4endl;
        
		for (int i = 0 ;i < n_hit; i++) {
            //			(*THC)[i]->SetPrimaryID(PrimaryID[(*THC)[i]->GetTrackID()]);
			
			G4cout  << " hit # " << i << G4endl;
            (*THC)[i]->Print();
		}
  	}
	
	if(CHC)
	{
		int n_hit = CHC->entries();
		G4cout  << "     " << n_hit
                << " hits are stored in calo collection" << G4endl;
        
		for(int i = 0; i < n_hit; i++) {
			G4cout  << " hit # " << i << G4endl;
            (*CHC)[i]->Print();
		}
	}
}

SToGS::PrintOutRunAction::PrintOutRunAction()
{
    ;
}
SToGS::PrintOutRunAction::~PrintOutRunAction()
{
    ;
}
void SToGS::PrintOutRunAction::BeginOfRunAction(const G4Run *aRun)
{
    G4cout  << ">>>> Begin of Run: " <<  aRun->GetRunID() << " " << aRun->GetNumberOfEventToBeProcessed() << G4endl ;
}
void SToGS::PrintOutRunAction::EndOfRunAction(const G4Run* aRun)
{
    G4cout  << ">>>> End of Run: " <<  aRun->GetRunID() << " " << aRun->GetNumberOfEvent() << G4endl ;
}

G4int printModuloEvt = 1;

void SToGS::PrintOutEventAction::BeginOfEventAction(const G4Event *evt)
{
    G4int evtNb = evt->GetEventID();
    if (evtNb%printModuloEvt == 0)
    {
        G4cout << ">>>> Begin of Event: " << evtNb + 1 << G4endl;
    }
}
void SToGS::PrintOutEventAction::EndOfEventAction(const G4Event *evt)
{
    G4int evtNb = evt->GetEventID();
    if (evtNb%printModuloEvt == 0)
    {
        G4cout << ">>>> End of Event: " << evtNb + 1 << G4endl;
    }
}

void SToGS::PrintOutTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
    G4cout << ">>>> Begin of Track: (PreUserTrackingAction) " << G4endl;
    G4cout << "  PARTICLE: " << aTrack->GetDefinition()->GetParticleName() << G4endl;
    G4cout << "  TRACK ID: " << aTrack->GetTrackID() << G4endl;
    G4cout << "  PARENT ID: " << aTrack->GetParentID() << G4endl;
    G4cout << "  TOTAL ENERGY: " << G4BestUnit(aTrack->GetTotalEnergy(), "Energy") << G4endl;
    G4cout << "  KINETIC ENERGY: " << G4BestUnit(aTrack->GetKineticEnergy(), "Energy") << G4endl;
    G4cout << "  VELOCITY: " << aTrack->GetVelocity() << G4endl;
}

void SToGS::PrintOutTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
    G4cout << ">>>> End of Track: (PostUserTrackingAction) " << G4endl;
    
    // verify if next volume exist -> if it does not exist it means that the particle went out of the world volume
    if(aTrack->GetNextVolume())
        G4cout << " VOLUME WHERE THE PARTICLE WAS KILLED: " << aTrack->GetTouchable()->GetVolume()->GetName() << G4endl;
    else
        G4cout << " THE PARTICLE WAS KILLED BECAUSE IT WENT OUT OF THE WORLD VOLUME" << G4endl;
}

void SToGS::PrintOutSteppingAction::UserSteppingAction(const G4Step* aStep)
{
    G4cout << ">>>> A Step: " << G4endl;

    // protection against out of world case
    if(!aStep->GetTrack()->GetTouchable()->GetVolume())
        return;
    
    // check if it is the same particle and the same event
    // if it is not -> clear vector of volumes of the previous particle and assign new values for the variables
    const G4int id_track = aStep->GetTrack()->GetTrackID();
    const G4int id_event = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    if( trackID != id_track || id_event != eventID ){
        trackID = id_track;
        eventID = id_event;
        volumesName.clear();
    }
    const G4String currentVolumeName = aStep->GetTrack()->GetTouchable()->GetVolume()->GetName();
    // check if it is the first time inside this volume
    if( std::find(volumesName.begin(), volumesName.end(), currentVolumeName) != volumesName.end())
        return;
    
    volumesName.push_back(currentVolumeName);
    
    const G4double stepLength = aStep->GetStepLength();
    const G4String volumeName = aStep->GetTrack()->GetTouchable()->GetVolume()->GetName();
    
    G4cout << " FIRST STEP IN VOLUME \"" << volumeName << "\" WITH A STEP LENGTH OF " << G4BestUnit(stepLength, "Length") << G4endl ;
}

SToGS::PrintOut::~PrintOut()
{
    ;
}
G4UserRunAction *SToGS::PrintOut::GetRunAction() const
{
    if ( fOption.contains("run") )
        return new SToGS::PrintOutRunAction();
    return 0x0;
}
G4UserEventAction *SToGS::PrintOut::GetEventAction() const
{
    if ( fOption.contains("event") )
        return new SToGS::PrintOutEventAction();
    return 0x0;
}
G4UserTrackingAction *SToGS::PrintOut::GetTrackingAction() const
{
    if ( fOption.contains("track") )
        return new PrintOutTrackingAction();
    return 0x0;
}
G4UserSteppingAction *SToGS::PrintOut::GetSteppingAction() const
{
    if ( fOption.contains("step") )
        return new SToGS::PrintOutSteppingAction();
    return 0x0;
}

void SToGS::PrintOut::BuildForMaster () const
{
#if G4VERSION_NUMBER < 1000
    G4cout << " *** ERROR *** SToGS::PrintOut::BuildForMaster() should never be called by this version of Geant4 " << G4endl;
#else
    SetUserAction( new SToGS::PrintOutRunAction() );
#endif
}
void SToGS::PrintOut::Build () const
{
#if G4VERSION_NUMBER < 1000
    G4cout << " *** ERROR *** SToGS::PrintOut::BuildForMaster() should never be called by this version of Geant4 " << G4endl;
#else
    SetUserAction( GetGun(fWhichGenerator.first,fWhichGenerator.second) );
    if ( fOption.contains("run") )
        SetUserAction( new SToGS::PrintOutRunAction() ) ;
    if ( fOption.contains("event") )
        SetUserAction( new SToGS::PrintOutEventAction() ) ;
    if ( fOption.contains("track") )
        SetUserAction( new SToGS::PrintOutTrackingAction() ) ;
    if ( fOption.contains("step") )
        SetUserAction( new SToGS::PrintOutSteppingAction() ) ;
#endif
}







