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

#include "SToGS_BaseROOTEventsActions.hh"

#include "SToGS_G4_TrackerSD.hh"
#include "SToGS_G4_CopClusterSD.hh"
#include "SToGS_G4_TrackInformation.hh"

#include "TTree.h"

#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"
#include "G4Event.hh"
#if G4MULTITHREADED
#include "G4Threading.hh"
#endif
#ifdef G4MULTITHREADED
#include "G4AutoLock.hh"
namespace { G4Mutex buildMutex = G4MUTEX_INITIALIZER; }
#endif

SToGS::BaseROOTEventsRun::BaseROOTEventsRun(TTree *tree) :
    SToGS::BaseROOTTreeRun(),
    fTree(tree),
    fPrimaryEvent(),
    fEvent()
{
    fTree->Branch("Pr.", &fPrimaryEvent);
    fTree->Branch("Ev.", &fEvent);
}

void SToGS::BaseROOTEventsRun::RecordEvent(const G4Event* evt)
{
    SToGS::TrackerHitsCollection *THC = 0x0; SToGS::CopClusterHitsCollection *CHC = 0x0;
    G4int nb_hits_tot = 0, nb_hits_trac = 0, nb_hits_calo = 0;
    
    // check there is something to write in files
	if( colltrackerID < 0 && collcaloID < 0 )
		return;
    
    fPrimaryEvent.Clear("");
    fEvent.Clear("");
    
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
			
			SBRPHit *hit = fPrimaryEvent.AddHit();
			
			K++;
			G4PrimaryParticle *prim = evt->GetPrimaryVertex(i)->GetPrimary(j); // get the next particle for vertex #i
            
            hit->fPDG =  prim->GetPDGcode();
            //
            //hit->fE = (std::sqrt( prim->GetMomentum().mag2() + prim->GetMass() * prim->GetMass()) - prim->GetMass())/CLHEP::keV;
            hit->fE = (std::sqrt( prim->GetMomentum().mag2() + prim->GetMass() * prim->GetMass()) - prim->GetMass())/CLHEP::MeV;
            //
 			hit->fX = evt->GetPrimaryVertex(i)->GetX0()/CLHEP::cm;
			hit->fY = evt->GetPrimaryVertex(i)->GetY0()/CLHEP::cm;
			hit->fZ = evt->GetPrimaryVertex(i)->GetZ0()/CLHEP::cm;
            //
            hit->fT = evt->GetPrimaryVertex(i)->GetT0();
            //
			hit->fPX = prim->GetPx();
			hit->fPY = prim->GetPy();
			hit->fPZ = prim->GetPz();
			
			hit->fFlag = prim->GetTrackID()-1;
			
            H += (std::sqrt( prim->GetMomentum().mag2() + prim->GetMass() * prim->GetMass()) - prim->GetMass())/CLHEP::MeV;
            //H += (std::sqrt( prim->GetMomentum().mag2() + prim->GetMass() * prim->GetMass()) - prim->GetMass())/CLHEP::keV;
		}
	}
	fPrimaryEvent.SetEMult(H,K);
    
	if ( nb_hits_trac ) {
		int n_hit = nb_hits_trac;
		
		K = n_hit; H = 0.0;
		for (int i = 0 ; i < n_hit; i++) {
			
			SBRHit *hit = fEvent.AddHit();
			
			if ( hit == 0x0 ) {
				G4cout << " Cannot add more that " << i << " hits in the collection " << G4endl;
				break;
			}
            
			hit->fPDG =  (*THC)[i]->GetPDGcode();

            hit->fE = (*THC)[i]->GetEdep()/CLHEP::MeV;
            //hit->fE = (*THC)[i]->GetEdep()/CLHEP::keV;
			hit->fX = (*THC)[i]->GetPos().x()/CLHEP::cm;
			hit->fY = (*THC)[i]->GetPos().y()/CLHEP::cm;
			hit->fZ = (*THC)[i]->GetPos().z()/CLHEP::cm;
			hit->fT = (*THC)[i]->GetToF();
			
			hit->fFlag = (*THC)[i]->GetPrimaryID()-1 ;
			hit->fUID  = (*THC)[i]->GetDetID();
            
            H += (*THC)[i]->GetEdep()/CLHEP::MeV;
            //H += (*THC)[i]->GetEdep()/CLHEP::keV;
		}
		fEvent.SetEMult(H,K);
  	}
	if ( nb_hits_calo )
	{
		int n_hit = nb_hits_calo;
		
		K = n_hit; H = 0.0;
		for (int i = 0 ; i < n_hit; i++) {
			
			SBRHit *hit = fEvent.AddHit();
			if ( hit == 0x0 ) {
				G4cout << " Cannot add more that " << i << " hits in the collection " << G4endl;
				break;
			}
			hit->fPDG = -1 ;

            hit->fE = (*CHC)[i]->GetEdep()/CLHEP::MeV;
            //hit->fE = (*CHC)[i]->GetEdep()/CLHEP::keV;
			hit->fX = (*CHC)[i]->GetPos().x()/CLHEP::cm;
			hit->fY = (*CHC)[i]->GetPos().y()/CLHEP::cm;
			hit->fZ = (*CHC)[i]->GetPos().z()/CLHEP::cm;
			hit->fT = (*CHC)[i]->GetToF();
			
			hit->fFlag = (*CHC)[i]->GetNbHits() ;
			hit->fUID  = (*CHC)[i]->GetDetID();
			
            H += (*CHC)[i]->GetEdep()/CLHEP::MeV;
            //H += (*CHC)[i]->GetEdep()/CLHEP::keV;
		}
		fEvent.SetEMult(H,K);
	}
    if ( fRecordOption == 1 ) {
		if ( fEvent.GetMultTot() > 0 ) fTree->Fill();
	}
	else fTree->Fill();
}

SToGS::BaseROOTEventsUserAction::BaseROOTEventsUserAction(G4String conffile):
    SToGS::BaseROOTTreeAction(conffile),
    fisOptical(0),
    fOpticalEventBeg(0x0),
    fOpticalEventEnd(0x0)
{
    // open the ascii file
    std::ifstream file(conffile.data());
	if ( file.is_open() == false ) {
		G4cout << " ** SToGS::BaseROOTEventsUserAction WARNING ** Cannot open file " << conffile << ", default parameters to be used "<< G4endl;
    }
	else {
        std::string aline; getline(file,aline);
        while ( file.good() && !file.eof() ) {
            
            if ( aline[0] == '#' ){
                getline(file,aline);
                continue;
            } // this line is a comment
            
            std::string key, what, option; std::istringstream decode(aline); decode.clear();
            decode >> key;
            
            if ( key == "scintillation:" ) { // to switch on following scintillation photons
                decode >> fisOptical ;
            }
            getline(file,aline);
        }
        file.close();
    }
}

G4Run* SToGS::BaseROOTEventsUserAction::GenerateRun()
{
#ifdef G4MULTITHREADED
    G4AutoLock lock(&buildMutex);
#endif
    
    G4Run* therun = 0x0; G4cout << " In SToGS::BaseROOTEventsUserAction, Generate a new Run " << G4endl;
    
    if ( fTree == 0x0 ) {
        fTree = new TTree(fTreeName.data(),fTreeTitle.data());
        fTree->SetDirectory(0x0);
        //        fTree->SetMaxTreeSize(fMaxEvents);
        
        if (fisOptical > 0) {
            fOpticalEventBeg = new SBROpticalEvent(); fOpticalEventEnd = new SBROpticalEvent();
            //
            fTree->Branch("OpticalEvBegin.", fOpticalEventBeg);
            fTree->Branch("OpticalEvEnd.", fOpticalEventEnd);
        }
    }
    
    // creates a new run. As the file is open bu AsciiRun no need to open something for master ... maybe one day to keep some globals ?
#if G4MULTITHREADED
    if ( G4Threading::IsWorkerThread() ) {
        //    if ( IsMaster() ) {
        BaseROOTEventsRun *loc_therun = new BaseROOTEventsRun(fTree);
        therun = loc_therun;
    }
    else therun = G4UserRunAction::GenerateRun();
#else
    therun = new BaseROOTEventsRun(fTree);
#endif
    return therun;
}

void SToGS::BaseROOTEventsUserAction::EndOfRunAction(const G4Run *therun)
{
    SToGS::BaseROOTTreeAction::EndOfRunAction(therun);
    if (fisOptical > 0) {
        delete fOpticalEventBeg; delete fOpticalEventEnd;
    }
}

void SToGS::BaseROOTEventsUserAction::BeginOfEventAction(const G4Event *evt)
{
	SToGS::BaseROOTTreeAction::BeginOfEventAction(evt);
	
	// clear events
    if ( fisOptical > 0 ) {
        fOpticalEventBeg->Clear("");
        fOpticalEventEnd->Clear("");
    }
}

void SToGS::BaseROOTEventsUserAction::PreUserTrackingAction(const G4Track *atrack)
{
    if ( fisOptical < 1 ) {
        return;
    }
	//Count what process generated the optical photons
	if(atrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()){
		// particle is optical photon
		if(atrack->GetParentID()>0){
			// particle is secondary
			if(atrack->GetCreatorProcess()->GetProcessName()=="Scintillation") {
                //				G4cout << "Scintillation Pre " << atrack->GetTrackID() << " "  << atrack->GetTrackLength() << G4endl;
                
                SToGS::PrimaryTrackInformation *pinfo = (SToGS::PrimaryTrackInformation *)atrack->GetUserInformation();
                G4int pid = atrack->GetTrackID();
                if ( pinfo ) {
                    pid = pinfo->GetPrimaryID();
                }
                
				SBROpticalHit *hit = fOpticalEventBeg->AddHit();
				
				hit->fX = atrack->GetPosition().x();
				hit->fY = atrack->GetPosition().y();
				hit->fZ = atrack->GetPosition().z();
				
				hit->fTA = atrack->GetGlobalTime();
				hit->fTL = atrack->GetLocalTime();
				
				hit->fLength = atrack->GetTrackLength();
                hit->fNbSteps = atrack->GetCurrentStepNumber();
                
  				hit->fPrimaryID   = pid-1;
				hit->fSecondaryID = atrack->GetVolume()->GetCopyNo(); // hit->fUID
				
			}
		}
	}
}
void SToGS::BaseROOTEventsUserAction::PostUserTrackingAction(const G4Track *atrack)
{
    if ( fisOptical < 1 ) {
        return;
    }
	//Count what process generated the optical photons
	if(atrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()){
		// particle is optical photon
		if(atrack->GetParentID()>0){
			// particle is secondary
			if(atrack->GetCreatorProcess()->GetProcessName()=="Scintillation") {
                //	G4cout << "Scintillation Post " << atrack->GetTrackID() << " "  << atrack->GetTrackLength()  << G4endl;
				
                SToGS::PrimaryTrackInformation *pinfo = (SToGS::PrimaryTrackInformation *)atrack->GetUserInformation();
                G4int pid = atrack->GetTrackID();
                if ( pinfo ) {
                    pid = pinfo->GetPrimaryID();
                }
                
                SBROpticalHit *hit = fOpticalEventEnd->AddHit();
                
				hit->fX = atrack->GetPosition().x();
				hit->fY = atrack->GetPosition().y();
				hit->fZ = atrack->GetPosition().z();
				
				hit->fTA = atrack->GetGlobalTime();
				hit->fTL = atrack->GetLocalTime();
				
				hit->fLength = atrack->GetTrackLength();
                hit->fNbSteps = atrack->GetCurrentStepNumber();
                
  				hit->fPrimaryID   = pid-1; // PrimaryID[(*THC)[i]->GetTrackID()]-1 ;
				hit->fSecondaryID = atrack->GetVolume()->GetCopyNo(); // hit->fUID
            }
		}
	}
}



