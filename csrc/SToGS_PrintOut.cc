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

#include "G4SDManager.hh"
#include "ParisTrackerSD.hh"
#include "ParisCaloSD.hh"
#include "ParisPrintOut.hh"
#include "G4Event.hh" 
  
ParisPrintOut::ParisPrintOut(G4String conf) : 
	ConcreteParisOutputManager(), 
	colltrackerID(-1), 
	collcaloID(-1)
	

 {
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	
	
	
	// keep the collection ID associated to this collection	
	colltrackerID = SDman->GetCollectionID(theTracker->GetCollectionName(0)); 
	
	
	// keep the collection ID associated to this collection	
	collcaloID = SDman->GetCollectionID(theCalo->GetCollectionName(0)); 
	
	SetPrintModulo(1);
	
	
	
}

ParisPrintOut::~ParisPrintOut() 
{  

}
 
 
 
 
 
void ParisPrintOut::EndOfEventAction(const G4Event *evt)
{
	ParisSingleHitsCollection *THC = NULL; ParisCaloHitsCollection *CHC = NULL;
	
	if( colltrackerID < 0 || collcaloID < 0 ) 
		return;
	
	G4HCofThisEvent * HCE = evt->GetHCofThisEvent();

	if(HCE)
	{
    		THC = (ParisSingleHitsCollection *)(HCE->GetHC(colltrackerID));
		CHC = (ParisCaloHitsCollection*)(HCE->GetHC(collcaloID));
  	}

	if(THC)
	{
		int n_hit = THC->entries();
		G4cout  << "     " << n_hit
			<< " hits are stored in tracker collection" << G4endl;
			
		for (int i = 0 ;i < n_hit; i++) { 
			(*THC)[i]->SetPrimaryID(PrimaryID[(*THC)[i]->GetTrackID()]);
			
			G4cout  << " hit # " << i << G4endl; (*THC)[i]->Print(); 
		}	
  	}
	
	if(CHC)
	{
		int n_hit = CHC->entries();
		G4cout  << "     " << n_hit
			<< " hits are stored in calo collection" << G4endl;

		for(int i = 0; i < n_hit; i++) { 
			G4cout  << " hit # " << i << G4endl; (*CHC)[i]->Print(); 
		}
	} 

	// call print out the event number from time to time
	ParisOutputManager::EndOfEventAction(evt);
	
	// reset information concerning the primaries
	ResetPrimaryID();
	
}



