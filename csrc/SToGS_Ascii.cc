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

#include "SToGS_Ascii.hh"
#include "SToGS_G4_TrackerSD.hh"
#include "SToGS_G4_CopClusterSD.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#if G4MULTITHREADED
#include "G4Threading.hh"
#endif

SToGS::AsciiRun::AsciiRun(std::ofstream &out) :
    G4Run(),
    colltrackerID(-1),
    collcaloID(-1),
    fOutputFile(out)
//    fEventMarkOff('$','\0')
{
   	G4SDManager* SDman = G4SDManager::GetSDMpointer();
    if ( SDman ) {
        // keep the collection ID associated to this collection
        colltrackerID = SDman->GetCollectionID("TrackerHits");
        
        // keep the collection ID associated to this collection
        collcaloID = SDman->GetCollectionID("CopClusterHits");
    }
}

SToGS::AsciiRun::~AsciiRun()
{
}

void SToGS::AsciiRun::RecordEvent(const G4Event* evt)
{
    SToGS::TrackerHitsCollection *THC = 0x0; SToGS::CopClusterHitsCollection *CHC = 0x0; G4int nb_hits = 0;

    // check there is something to write in files
	if( colltrackerID < 0 && collcaloID < 0 )
		return;
    
    G4HCofThisEvent *HCE = evt->GetHCofThisEvent();
	if(HCE)
	{
        THC = (SToGS::TrackerHitsCollection *)(HCE->GetHC(colltrackerID));
        if ( colltrackerID > -1 ) {
            nb_hits += THC->entries();
        }
		CHC = (SToGS::CopClusterHitsCollection*)(HCE->GetHC(collcaloID));
        if ( collcaloID > -1 ) {
            nb_hits += CHC->entries();
        }
  	}
    /*
    if ( fRecordOption == 1 ) {
        if ( nb_hits == 0 )
            return;
    }
    */
    // write primary part
    G4int istart = 0, iend = evt->GetNumberOfPrimaryVertex();
	
	for( G4int i = istart; i < iend ; i++ ){
        G4int jstart = 0, jend = evt->GetPrimaryVertex(i)->GetNumberOfParticle(); // to get the next vertex
		for( G4int j = jstart; j < jend; j++ ) {
            
			G4PrimaryParticle *prim = evt->GetPrimaryVertex(i)->GetPrimary(j); // get the next particle for vertex #i
			
			fOutputFile
                << std::setw(2)
                << "P "	<< " "
                << std::setw(4)
                << prim->GetTrackID() << " "
                << std::setw(3)
                << iend*jend << " "
                << std::setw(7)
                << (std::sqrt( prim->GetMomentum().mag2() + prim->GetMass() * prim->GetMass()) - prim->GetMass())/CLHEP::keV << " "
                << prim->GetPx() << " " << prim->GetPy() << " " << prim->GetPz()
                << std::endl;
		}
	}
	if(THC && colltrackerID != -1)
	{
		int n_hit = THC->entries();
		for (int i = 0 ; i < n_hit; i++) {
            fOutputFile
                << std::setw(2)
                << "H "	<< " "
                << std::setw(4)
                // << PrimaryID[(*THC)[i]->GetTrackID()]  << " " TODO through tracking info
                << 1 << " "
                << std::setw(3)
                << (*THC)[i]->GetDetID() << " "
                << std::setw(7)
                << (*THC)[i]->GetEdep()/CLHEP::keV << " "
                << (*THC)[i]->GetPos().x()/CLHEP::cm << " " << (*THC)[i]->GetPos().y()/CLHEP::cm << " " << (*THC)[i]->GetPos().z()/CLHEP::cm << " "
                << std::setw(2)
                << (*THC)[i]->GetToF()/CLHEP::ns
                << std::endl;
		}
  	}
	if(CHC && collcaloID != -1)
	{
		int n_hit = CHC->entries();
		for (int i = 0 ; i < n_hit; i++) {
            fOutputFile
                << std::setw(2)
                << "H "	<< " "
                << std::setw(4)
                << 1 << " "
                << std::setw(3)
                << (*CHC)[i]->GetDetID() << " "
                << std::setw(7)
                << (*CHC)[i]->GetEdep()/CLHEP::keV << " "
                << (*CHC)[i]->GetPos().x()/CLHEP::cm << " " << (*CHC)[i]->GetPos().y()/CLHEP::cm << " " << (*CHC)[i]->GetPos().z()/CLHEP::cm << " "
                << std::setw(2)
                << (*CHC)[i]->GetToF()/CLHEP::ns
                << std::endl;
		}
	}
}

SToGS::AsciiAction::AsciiAction(G4String conf):
    AllActions(conf),
    fOutputFile(),
    fPathToData("./"),
    fBaseName("SToGS_Ascii"),
    fMaxEvents(500000),
    fRecordOption(0),
    fPrintModulo(1000)
{
    // open the ascii file
    std::ifstream file(conf.data());
	if ( file.is_open() == false ) {
		G4cout << " ** SToGS::Ascii WARNING ** Cannot open file " << conf << ", default parameters to be used "<< G4endl;
    }
	else {
        std::string aline; getline(file,aline);
        while ( file.good() && !file.eof() ) {
            
            if ( aline[0] == '#' ){
                getline(file,aline);
                continue;
            } // this line is a comment
            
            std::string key, what, option; std::istringstream decode(aline);
            decode >> key;
            
            if ( key == "path:" ) {
                decode >> fPathToData;
            }
            if ( key == "basename:" ) {
                decode >> fBaseName;
            }
            if ( key == "modulo:" ) {
                decode >> fPrintModulo;
            }
            if ( key == "max_event:" ) {
                decode >> fMaxEvents;
            }
            if ( key == "record_option:" ) {
                decode >> fRecordOption;
            }
            getline(file,aline);
        }
        file.close();
    }
}

void SToGS::AsciiAction::OpenStream(G4int run_id)
{
    // for MT, thread id used in the name of the file to avoid pb
    G4int thread_id = 0;
#if G4MULTITHREADED
    thread_id = G4Threading::G4GetThreadId();
#endif
    std::ostringstream filename;
    filename.clear();
    filename << fPathToData << "/" << fBaseName << "_" << std::setfill('0') << std::setw(2) << thread_id << "_" << std::setw(4) << run_id << ".g4event";
    
    // close the file if already open
    if ( fOutputFile.is_open() )
    {
        fOutputFile.close();
        fOutputFile.clear();
    }
    fOutputFile.open(filename.str().data());
    //
    if ( fOutputFile ) {
		G4cout << " The File " << filename.str() << " is open to record data " << G4endl;
		
		fOutputFile << "#" << std::endl;
		fOutputFile << "# FORMAT: 'C' ID1  ID2  Energy(keV)    x(cm)     y(cm)     z(cm) " << std::endl;
		fOutputFile << "#    with X=P for primary gammas and in this case:" << std::endl;
		fOutputFile << "#        ID1 is the vertexID, ID2 the total number of primaries" << std::endl;
		fOutputFile << "#        E the energy of the emitted gamma " << std::endl;
		fOutputFile << "#        x,y,z represents the momentum" << std::endl;
		fOutputFile << "#    with X=H for a single hit and in this case:" << std::endl;
		fOutputFile << "#        ID1 is the vertexID it comes from, ID2 the detector number" << std::endl;
		fOutputFile << "#        E the energy of the impact " << std::endl;
		fOutputFile << "#        x,y,z represents its position" << std::endl;
		fOutputFile << "#" << std::endl;
	} else {
		G4cout << " *** ERROR *** cannot open " << filename.str() << " to record data " << G4endl;
	}
}

/*
std::pair<G4int, G4int> SToGS::AsciiAction::HitsinCollection(const G4Event)
{
    
}
*/

G4Run* SToGS::AsciiAction::GenerateRun()
{
    G4Run* therun = 0x0; G4cout << " In AsciiAction, Generate a new Run " << G4endl;

    // creates a new run. As the file is open bu AsciiRun no need to open something for master ... maybe one day to keep some globals ?
#if G4MULTITHREADED
    if ( G4Threading::IsWorkerThread() ) {
//    if ( IsMaster() ) {
        SToGS::AsciiRun *loc_therun = new SToGS::AsciiRun(fOutputFile);
        therun = loc_therun;
    }
    else therun = G4UserRunAction::GenerateRun();
#else
    therun = new SToGS::AsciiRun(fOutputFile);
#endif
    return therun;
}
    
void SToGS::AsciiAction::BeginOfRunAction(const G4Run *aRun)
{
	G4cout << " In AsciiAction, Begin of Run " << aRun->GetRunID()
           << ", Number of events to be simulated " << aRun->GetNumberOfEventToBeProcessed() << G4endl;
    //
    OpenStream(aRun->GetRunID());
}
void SToGS::AsciiAction::EndOfRunAction(const G4Run* aRun)
{
 	G4cout << " In AsciiAction, End of Run " << aRun->GetRunID() << " " << aRun->GetNumberOfEvent() << G4endl;
    //
    if ( fOutputFile.is_open() )
    {
        fOutputFile.close();
        fOutputFile.clear();
    }
}
void SToGS::AsciiAction::BeginOfEventAction(const G4Event *evt)
{
    G4int evtNb = evt->GetEventID();
    if ( evtNb % fPrintModulo == 0 ) {
    	G4cout << " In AsciiAction, Begin of event: " << evtNb << G4endl;
    }
}
void SToGS::AsciiAction::EndOfEventAction(const G4Event *evt)
{
    G4int evtNb = evt->GetEventID();
    if ( evtNb % fPrintModulo == 0 ) {
    	G4cout << " In AsciiAction, End of event: " << evtNb << G4endl;
    }
}

