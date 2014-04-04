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

#include "ParisUserAction.hh"
#include "SToGS_G4_TrackerSD.hh"
#include "SToGS_G4_CopClusterSD.hh"

#include "TTree.h"

#include "G4Event.hh"
#if G4MULTITHREADED
#include "G4Threading.hh"
#endif
#ifdef G4MULTITHREADED
#include "G4AutoLock.hh"
namespace { G4Mutex buildMutex = G4MUTEX_INITIALIZER; }
#endif

ParisEventRun::ParisEventRun(TTree *tree, PEvent *primaryevent, PEvent *event) :
    SToGS::BaseROOTTreeRun(),
    fTree(tree),
	fPrimaryEvent(primaryevent),
	fEvent(event)
{
}

void ParisEventRun::RecordEvent(const G4Event* evt)
{
    SToGS::TrackerHitsCollection *THC = 0x0; SToGS::CopClusterHitsCollection *CHC = 0x0; G4int nb_hits_tot = 0, nb_hits_trac = 0, nb_hits_calo = 0;
    
    // check there is something to write in files
	if( colltrackerID < 0 && collcaloID < 0 )
		return;
    
    fPrimaryEvent->Clear("");
    fEvent->Clear("");
    
    G4HCofThisEvent *HCE = evt->GetHCofThisEvent();
	if(HCE)
	{
        THC = (SToGS::TrackerHitsCollection *)(HCE->GetHC(colltrackerID));
        if ( colltrackerID > -1 ) {
            nb_hits_trac += THC->entries();
        }
		CHC = (SToGS::CopClusterHitsCollection*)(HCE->GetHC(collcaloID));
        if ( collcaloID > -1 ) {
            nb_hits_calo += CHC->entries();
        }
        nb_hits_tot = nb_hits_trac + nb_hits_calo;
  	}

	// write primary part
	G4int istart = 0, iend = evt->GetNumberOfPrimaryVertex();
    G4int K = 0; G4double H = 0.0;
	for( G4int i = istart; i < iend ; i++ ){
		
		G4int jstart = 0, jend = evt->GetPrimaryVertex(i)->GetNumberOfParticle(); // to get the next vertex
		for( G4int j = jstart; j < jend; j++ ) {
			
			PHit *hit = fPrimaryEvent->AddHit();
			
			K++;
			G4PrimaryParticle *prim = evt->GetPrimaryVertex(i)->GetPrimary(j); // get the next particle for vertex #i

			hit->fE = (std::sqrt( prim->GetMomentum().mag2() + prim->GetMass() * prim->GetMass()) - prim->GetMass())/keV;
			hit->fX = prim->GetPx();
			hit->fY = prim->GetPy();
			hit->fZ = prim->GetPz();
			
			hit->fFlag = prim->GetTrackID()-1;
			
			H += (std::sqrt( prim->GetMomentum().mag2() + prim->GetMass() * prim->GetMass()) - prim->GetMass())/MeV;
		}
	} 
	fPrimaryEvent->SetHK(H,K);
    
	if ( nb_hits_trac ) {
		int n_hit = nb_hits_trac;
		
		K = n_hit; H = 0.0;
		for (int i = 0 ; i < n_hit; i++) { 
			
			PHit *hit = fEvent->AddHit();
			
			if ( hit == 0x0 ) {
				G4cout << " Cannot add more that " << i << " hits in the collection " << G4endl;
				break;
			}
	
			hit->fID = PHit::kGamma ;

			hit->fE = (*THC)[i]->GetEdep()/keV;
			hit->fX = (*THC)[i]->GetPos().x()/cm;
			hit->fY = (*THC)[i]->GetPos().y()/cm;
			hit->fZ = (*THC)[i]->GetPos().z()/cm;
			hit->fT = (*THC)[i]->GetToF();
			
			hit->fFlag = (*THC)[i]->GetPrimaryID()-1 ;
			hit->fUID  = (*THC)[i]->GetDetID();		

			H += (*THC)[i]->GetEdep()/MeV;
		}
		fEvent->SetHK(H,K);
  	}  
	if ( nb_hits_calo )
	{
		int n_hit = nb_hits_calo;
		
		K = n_hit; H = 0.0;
		for (int i = 0 ; i < n_hit; i++) { 
			
			PHit *hit = fEvent->AddHit();
			if ( hit == 0x0 ) {
				G4cout << " Cannot add more that " << i << " hits in the collection " << G4endl;
				break;
			}			
			hit->fID = PHit::kUnknown ;		

			hit->fE = (*CHC)[i]->GetEdep()/keV;
			hit->fX = (*CHC)[i]->GetPos().x()/cm;
			hit->fY = (*CHC)[i]->GetPos().y()/cm;
			hit->fZ = (*CHC)[i]->GetPos().z()/cm;
			hit->fT = (*CHC)[i]->GetToF();
			
			hit->fFlag = (*CHC)[i]->GetNbHits() ;		
			hit->fUID  = (*CHC)[i]->GetDetID();		
			
			H += (*CHC)[i]->GetEdep()/MeV;
		}
		fEvent->SetHK(H,K);
	}
    /*
    if ( fRecordOption == 1 ) {
		if ( fEvent->GetK() > 0 ) fTree->Fill();
	}
	else fTree->Fill(); */
}

ParisUserAction::ParisUserAction(G4String conffile):
    SToGS::BaseROOTTreeAction(conffile),
    fPrimaryEvent(),
    fEvent()
{
    
}

G4Run* ParisUserAction::GenerateRun()
{
    G4Run* therun = 0x0; G4cout << " In ParisUserAction, Generate a new Run " << G4endl;
    
    // creates a new run. As the file is open bu AsciiRun no need to open something for master ... maybe one day to keep some globals ?
#if G4MULTITHREADED
    if ( G4Threading::IsWorkerThread() ) {
        //    if ( IsMaster() ) {
        ParisEventRun *loc_therun = new ParisEventRun(fTree,&fPrimaryEvent,&fEvent);
        therun = loc_therun;
    }
    else therun = G4UserRunAction::GenerateRun();
#else
    therun = new ParisEventRun(fTree,&fPrimaryEvent,&fEvent);
#endif
    return therun;
}


