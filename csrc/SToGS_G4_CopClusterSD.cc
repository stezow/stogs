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

#include "SToGS_G4_CopClusterSD.hh"
#include "SToGS_G4_CopClusterHit.hh"

#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"

SToGS::CopClusterSD::CopClusterSD(G4String name): G4VSensitiveDetector(name)
{
	G4String HCname = "CopClusterHits";
    collectionName.insert(HCname);
}
SToGS::CopClusterSD::~CopClusterSD()
{
    ;
}
void SToGS::CopClusterSD::Initialize(G4HCofThisEvent* HCE)
{
//	G4cout << " In SToGS::CopClusterSD::Initialize " << G4endl;
	static int HCID = -1;
	caloCollection = new SToGS::CopClusterHitsCollection(SensitiveDetectorName,collectionName[0]); 
	
	if ( HCID < 0 ) 
		HCID = GetCollectionID(0);
		
	HCE->AddHitsCollection(HCID,caloCollection);
//	G4cout << " Out SToGS::CopClusterSD::Initialize " << G4endl;
}
G4bool SToGS::CopClusterSD::ProcessHits(G4Step* aStep, G4TouchableHistory *touch)
{
	G4String tmp;
	
//	G4cout << " In SToGS::CopClusterSD::ProcessHits" << G4endl;
	// nothing to be stored if no energy 
	
	G4double edep = aStep->GetTotalEnergyDeposit();  if ( edep == 0. ) return false;
	
	G4bool anew = true; G4int id = aStep->GetTrack()->GetVolume()->GetCopyNo(); SToGS::CopClusterHit *newHit; 
	
	for( unsigned int i = 0; i < caloCollection->GetSize(); i++ ){ // look 
		newHit = (SToGS::CopClusterHit *)caloCollection->GetHit(i); if ( id == newHit->GetDetID() ) { anew = false; break; }
	}

	if ( anew ) { // a new hit is created
		newHit = new SToGS::CopClusterHit();
	
		newHit->SetEdep( edep );
		newHit->SetPos(  aStep->GetPostStepPoint()->GetPosition() );
		newHit->SetToF ( aStep->GetPostStepPoint()->GetGlobalTime() );	
	
		newHit->SetDetName  (aStep->GetTrack()->GetVolume()->GetName());
		newHit->SetDetID    (id);	
	
		// newHit->SetMotherDetName(touch->GetVolume()->GetName());
   		newHit->SetMotherID(aStep->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(1)); 
	
		newHit->SetNbHits(0);
		// add this hit to the collection
		caloCollection->insert( newHit );
	}
	else { // add the informations concerning this hit to the existing one
		newHit->AddOneHit();
		newHit->AddEdep( edep );
		newHit->AddPos ( edep, aStep->GetPostStepPoint()->GetPosition() );
		newHit->AddToF ( edep, aStep->GetPostStepPoint()->GetGlobalTime() );
	}
	

//	G4cout << " Out SToGS::CopClusterSD::ProcessHits" << G4endl;
	return true;
}
void SToGS::CopClusterSD::EndOfEvent(G4HCofThisEvent*)
{
	for( unsigned int i = 0; i < caloCollection->GetSize(); i++ ){ // to calculate properly the mean values
		SToGS::CopClusterHit *hit = (SToGS::CopClusterHit *)caloCollection->GetHit(i); 
		hit->EndOfEvent();
	}	
}
void SToGS::CopClusterSD::clear()
{
}
void SToGS::CopClusterSD::DrawAll()
{
} 
void SToGS::CopClusterSD::PrintAll()
{
} 
