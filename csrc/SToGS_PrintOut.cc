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
#include "SToGS_G4_CopClusterSD.hh"
#include "SToGS_G4_TrackInformation.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#ifdef G4MULTITHREADED
#include "G4AutoLock.hh"
namespace { G4Mutex buildMutex = G4MUTEX_INITIALIZER; }
#endif

SToGS::PrintOutRun::PrintOutRun() :
    G4Run(),
    colltrackerID(-1),
    collcaloID(-1)
{
   	G4SDManager* SDman = G4SDManager::GetSDMpointer();
    if ( SDman ) {
        // keep the collection ID associated to this collection
        colltrackerID = SDman->GetCollectionID("TrackerHits");
        
        // keep the collection ID associated to this collection
        collcaloID = SDman->GetCollectionID("CopClusterHits");
    }
}
SToGS::PrintOutRun::~PrintOutRun()
{
    ;
}
void SToGS::PrintOutRun::RecordEvent(const G4Event* evt)
{
    SToGS::TrackerHitsCollection *THC = 0x0; SToGS::CopClusterHitsCollection *CHC = 0x0;
	
	if( colltrackerID < 0 && collcaloID < 0  ) {
		return;
    }
    G4cout << "  Beginning Record of Event: " << evt->GetEventID() + 1 << G4endl;

	
	G4HCofThisEvent *HCE = evt->GetHCofThisEvent();
    
	if(HCE)
	{
        THC = (SToGS::TrackerHitsCollection *)(HCE->GetHC(colltrackerID));
		CHC = (SToGS::CopClusterHitsCollection*)(HCE->GetHC(collcaloID));
  	}
    
	if(THC && colltrackerID > -1)
	{
		int n_hit = THC->entries();
		G4cout  << "  " << n_hit
                << " hits are stored in tracker collection" << G4endl;
        
		for (int i = 0 ;i < n_hit; i++) {
			G4cout  << "  hit # " << i << " " << G4BestUnit((*THC)[i]->GetEdep(), "Energy") << G4endl;
//            (*THC)[i]->Print();
		}
  	}
	
	if(CHC && collcaloID > -1)
	{
		int n_hit = CHC->entries();
		G4cout  << "  " << n_hit
                << " hits are stored in calo collection" << G4endl;
        
		for(int i = 0; i < n_hit; i++) {
			G4cout  << "  hit # " << i << G4BestUnit((*CHC)[i]->GetEdep(), "Energy") << G4endl;
//            (*CHC)[i]->Print();
		}
	}
    
    G4cout << "  End Record of Event: " << evt->GetEventID() + 1 << G4endl;
}

G4Run* SToGS::PrintOutAction::GenerateRun()
{
    G4Run* therun = 0x0; G4cout << " In PrintOutAction, Generate a new Run " << G4endl;
    
    // creates a new run. As the file is open bu AsciiRun no need to open something for master ... maybe one day to keep some globals ?
#if G4MULTITHREADED
    if ( G4Threading::IsWorkerThread() ) {
        //    if ( IsMaster() ) {
        SToGS::PrintOutRun *loc_therun = new SToGS::PrintOutRun();
        therun = loc_therun;
    }
    else therun = G4UserRunAction::GenerateRun();
#else
    therun = new SToGS::PrintOutRun();
#endif
    return therun;
}

void SToGS::PrintOutAction::BeginOfRunAction(const G4Run *aRun)
{
    G4cout  << "Begin of Run: " <<  aRun->GetRunID() << " " << aRun->GetNumberOfEventToBeProcessed() << G4endl ;
}
void SToGS::PrintOutAction::EndOfRunAction(const G4Run* aRun)
{
    G4cout  << "End of Run: " <<  aRun->GetRunID() << " " << aRun->GetNumberOfEvent() << G4endl ;
}

void SToGS::PrintOutAction::BeginOfEventAction(const G4Event *evt)
{
    G4int evtNb = evt->GetEventID();
    G4cout << "  Begin of Event: " << evtNb + 1 << G4endl;
}
void SToGS::PrintOutAction::EndOfEventAction(const G4Event *evt)
{
    G4int evtNb = evt->GetEventID();
    G4cout << "  End of Event: " << evtNb + 1 << G4endl;
}

void SToGS::PrintOutAction::PreUserTrackingAction(const G4Track* aTrack)
{
    //
    G4cout << "    Begin of Track: (PreUserTrackingAction) " << G4endl;
    G4cout << "     PARTICLE: " << aTrack->GetDefinition()->GetParticleName() << G4endl;
    G4cout << "     TRACK ID: " << aTrack->GetTrackID() << G4endl;
    G4cout << "     PARENT ID: " << aTrack->GetParentID() << G4endl;
    G4cout << "     TOTAL ENERGY: " << G4BestUnit(aTrack->GetTotalEnergy(), "Energy") << G4endl;
    G4cout << "     KINETIC ENERGY: " << G4BestUnit(aTrack->GetKineticEnergy(), "Energy") << G4endl;
    G4cout << "     VELOCITY: " << aTrack->GetVelocity() << G4endl;
}

void SToGS::PrintOutAction::PostUserTrackingAction(const G4Track* aTrack)
{
    // verify if next volume exist -> if it does not exist it means that the particle went out of the world volume
    if(aTrack->GetNextVolume())
        G4cout << "     VOLUME WHERE THE PARTICLE WAS KILLED: " << aTrack->GetTouchable()->GetVolume()->GetName() << G4endl;
    else
        G4cout << "     THE PARTICLE WAS KILLED BECAUSE IT WENT OUT OF THE WORLD VOLUME" << G4endl;

    G4cout << "    End of Track: (PostUserTrackingAction) " << G4endl;
}

void SToGS::PrintOutAction::UserSteppingAction(const G4Step* aStep)
{
    G4Track* theTrack = aStep->GetTrack();

    G4int primaryID = theTrack->GetParentID(); SToGS::PrimaryTrackInformation *pinfo = (SToGS::PrimaryTrackInformation *)theTrack->GetUserInformation();
    if ( pinfo )
       primaryID = pinfo->GetPrimaryID();
    
    G4cout << "        A Step: " << G4endl;
	G4cout << "         trackID: " << theTrack->GetTrackID()
        << ", parentID: " << theTrack->GetParentID()
        << ", primaryID: " << primaryID
        << ", particleName: " << theTrack->GetDefinition()->GetParticleName()
        << ", PDG: " << theTrack->GetDefinition()->GetPDGEncoding ()
        << ", processName: " << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
        << G4endl
        << "         detID: " << theTrack->GetVolume()->GetCopyNo()
        << ", detName: " << theTrack->GetVolume()->GetName()
        << ", motherID: " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(1)
        << ", motherDetName: " /* << motherDetName */ << G4endl;
    G4cout << "         edep: " <<  aStep->GetTotalEnergyDeposit() / CLHEP::keV << " keV"
        << ", pos: " << aStep->GetPostStepPoint()->GetPosition()
        << ", ToF: " << aStep->GetPostStepPoint()->GetGlobalTime()  / CLHEP::ns << " ns" << G4endl;
}

G4UserRunAction *SToGS::PrintOut::GetRunAction() const
{
    if ( fUserActionOption.contains("run") )
        return SToGS::AllInOneUserActionInitialization<PrintOutAction>::GetRunAction();
    return 0x0;
}
G4UserEventAction *SToGS::PrintOut::GetEventAction() const
{
    if ( fUserActionOption.contains("event") )
        return SToGS::AllInOneUserActionInitialization<PrintOutAction>::GetEventAction();
    return 0x0;
}
G4UserTrackingAction *SToGS::PrintOut::GetTrackingAction() const
{
    if ( fUserActionOption.contains("track") )
        return SToGS::AllInOneUserActionInitialization<PrintOutAction>::GetTrackingAction();
    return 0x0;
}
G4UserSteppingAction *SToGS::PrintOut::GetSteppingAction() const
{
    if ( fUserActionOption.contains("step") )
        return SToGS::AllInOneUserActionInitialization<PrintOutAction>::GetSteppingAction();
    return 0x0;
}

void SToGS::PrintOut::Build () const
{
#ifdef G4MULTITHREADED
    G4AutoLock lock(&buildMutex);
#endif
    
    G4cout << " ------ INF ------ from SToGS::PrintOut::Build() " << G4endl;
#if G4VERSION_NUMBER < 1000
    G4cout << " *** ERROR *** SToGS::PrintOut::Build() should never be called by this version of Geant4 " << G4endl;
#else
    SetUserAction( GetGun(fWhichGenerator.first,fWhichGenerator.second) );
    
    PrintOutAction *all_actions = new PrintOutAction(fUserActionOption); fAllUserAction.push_back(all_actions);
    SetUserAction( new SToGS::RunAction(all_actions) ) ;
    if ( fUserActionOption.contains("event") )
        SetUserAction( new SToGS::EventAction(all_actions) ) ;
    if ( fUserActionOption.contains("track") )
        SetUserAction( new SToGS::TrackingAction(all_actions) ) ;
    if ( fUserActionOption.contains("step") )
        SetUserAction( new SToGS::SteppingAction(all_actions) ) ;
#endif
    G4cout << " ------ END ------  from SToGS::PrintOut::Build() " << G4endl;
}






