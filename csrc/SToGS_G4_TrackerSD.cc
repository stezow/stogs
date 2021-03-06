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

#include "SToGS_G4_TrackerSD.hh"
#include "SToGS_G4_TrackerHit.hh"
#include "SToGS_G4_TrackInformation.hh"

#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include "G4OpticalPhoton.hh"

SToGS::TrackerSD::TrackerSD(G4String name): G4VSensitiveDetector(name)
{
	G4String HCname = "TrackerHits";
    collectionName.insert(HCname);
	
}
SToGS::TrackerSD::~TrackerSD()
{
    ;
}
void SToGS::TrackerSD::Initialize(G4HCofThisEvent* HCE)
{
//	G4cout << " In SToGS::TrackerSD::Initialize " << G4endl;
	static int HCID = -1;
	trackerCollection = new SToGS::TrackerHitsCollection(SensitiveDetectorName,collectionName[0]); 
	
	if ( HCID < 0 ) 
		HCID = GetCollectionID(0);
	HCE->AddHitsCollection(HCID,trackerCollection);
//	G4cout << " Out SToGS::TrackerSD::Initialize " << G4endl;
}

G4bool SToGS::TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory * /*touch*/)
{
	//	G4cout << " In SToGS::TrackerSD::ProcessHits" << G4endl;

	G4String tmp; G4int temp_code;
	
	// avoid keeping optical
	G4Track* theTrack = aStep->GetTrack();
	if ( theTrack && theTrack->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() )
		return false ;
	
	// nothing to be stored if no energy 
	G4double edep = aStep->GetTotalEnergyDeposit();  
	if ( edep == 0. ) {
		return false;
	}
	
	// a new hit is created
	SToGS::TrackerHit *newHit = new SToGS::TrackerHit();
	
	// set hit properties  
	newHit->SetTrackID(theTrack->GetTrackID());
	newHit->SetParentID(theTrack->GetParentID());
    
    SToGS::PrimaryTrackInformation *pinfo = (SToGS::PrimaryTrackInformation *)theTrack->GetUserInformation();
    if ( pinfo )
        newHit->SetPrimaryID(pinfo->GetPrimaryID());
    else
        newHit->SetPrimaryID(theTrack->GetParentID()); // in case primary is not know, at least parent
    
	newHit->SetEdep( edep );
	newHit->SetPos(  aStep->GetPostStepPoint()->GetPosition() );
	newHit->SetToF ( aStep->GetPostStepPoint()->GetGlobalTime() );	
	
	newHit->SetDetName  (aStep->GetTrack()->GetVolume()->GetName());
    
    const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int depth = touchable->GetHistoryDepth();
    
    G4int detector_number = touchable->GetCopyNumber(0); // get curent detector number and add offset of the on-top structure
    if ( depth > 1 ) {
        detector_number += touchable->GetCopyNumber(depth-1);
        /*
        G4cout << " 0 " << touchable->GetCopyNumber(0) << G4endl;;
        G4cout << " 1 " << touchable->GetCopyNumber(1) << G4endl;;
        G4cout << " 2 " << touchable->GetCopyNumber(2) << G4endl;;
        G4cout << " 3 " << touchable->GetCopyNumber(3) << G4endl;;
        */
    }
    
	newHit->SetDetID(detector_number);
	
//	newHit->SetMotherDetName(touch->GetVolume()->GetName());
   	newHit->SetMotherID(aStep->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(1)); 
	
	tmp = aStep->GetTrack()->GetDefinition()->GetParticleName();
	newHit->SetParticleName( tmp );

	temp_code = aStep->GetTrack()->GetDefinition()->GetPDGEncoding ();
	newHit->SetPDGcode(temp_code);
	
	tmp = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();	
	newHit->SetProcessName( tmp );
		
	// add this hit to the collection
	trackerCollection->insert( newHit );

//	G4cout << " Out SToGS::TrackerSD::ProcessHits" << G4endl;
	return true;
}

void SToGS::TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
}

void SToGS::TrackerSD::clear()
{
} 

void SToGS::TrackerSD::DrawAll()
{
} 

void SToGS::TrackerSD::PrintAll()
{
} 
