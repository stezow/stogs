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

SToGS::AsciiRun::AsciiRun(G4String filename) :
    G4Run(),
    colltrackerID(-1),
    collcaloID(-1),
    fOutputFile(),
    fFilename(filename),
    fEventMarkOff('$','\0'),
    fMaxEvents(500000),
    fRecordOption(0)
{
   	G4SDManager* SDman = G4SDManager::GetSDMpointer();
    if ( SDman ) {
        // keep the collection ID associated to this collection
        if ( G4SDManager::GetSDMpointer()->FindSensitiveDetector("Tracker",false) )
            colltrackerID = SDman->GetCollectionID("TrackerHits");
        
        // keep the collection ID associated to this collection
        if ( G4SDManager::GetSDMpointer()->FindSensitiveDetector("CopCluster",false) )
            collcaloID = SDman->GetCollectionID("CopClusterHits");
    }
}
void SToGS::AsciiRun::OpenStream()
{
    // for MT, thread id used in the name of the file to avoid pb
    G4int thread_id = 0;
#if G4MULTITHREADED
    thread_id = G4Threading::G4GetThreadId();
#endif
    std::ostringstream filename;
    filename.clear();
    filename << fFilename << "_" << thread_id << "_" << GetRunID() << ".ascii";
    
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

SToGS::AsciiRun::~AsciiRun()
{
    if ( fOutputFile.is_open() )
    {
        fOutputFile.close();
        fOutputFile.clear();
    }
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
        if ( THC ) {
            nb_hits += THC->entries();
        }
		CHC = (SToGS::CopClusterHitsCollection*)(HCE->GetHC(collcaloID));
        if ( CHC ) {
            nb_hits += CHC->entries();
        }
  	}
    if ( fRecordOption == 1 ) {
        if ( nb_hits == 0 )
            return;
    }
    
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

SToGS::AsciiRunAction::AsciiRunAction(G4String path, G4String basename, G4int maxevent, G4int record_opt) :
    G4UserRunAction(),
    fPathToData(path),
    fBaseName(basename),
    fMaxEvents(maxevent),
    fRecordOption(record_opt)
{
    ;
}
SToGS::AsciiRunAction::~AsciiRunAction()
{
    ;
}
G4Run* SToGS::AsciiRunAction::GenerateRun()
{
    G4String filename = fPathToData;
    filename += fBaseName;
    
    // creates a new run. As the file is open bu AsciiRun no need to open something for master ... maybe one day to keep some globals ?
    SToGS::AsciiRun *localRun = new SToGS::AsciiRun(filename);
    if ( IsMaster() == false ) {
        localRun->OpenStream();
    }
    return localRun;
}
void SToGS::AsciiRunAction::BeginOfRunAction(const G4Run *aRun)
{
     if ( !IsMaster() ) {
         return;
     }
    
    const SToGS::AsciiRun *localRun = static_cast<const SToGS::AsciiRun *> (aRun);
	G4cout << " IN AsciiRunAction, Begin of Run " << aRun->GetRunID()
           << ", Nomber of events to be simulated " << localRun->GetNumberOfEventToBeProcessed() << G4endl;
}
void SToGS::AsciiRunAction::EndOfRunAction(const G4Run* aRun)
{
 	G4cout << " IN AsciiRunAction, End of Run " << aRun->GetRunID() <<  G4endl;
}
void SToGS::AsciiEventAction::BeginOfEventAction(const G4Event *evt)
{
    G4int evtNb = evt->GetEventID();
    
    if ( evtNb % fPrintModulo == 0) {
    	G4cout << " Begin of event: " << evtNb << G4endl;
    }
}
void SToGS::AsciiEventAction::EndOfEventAction(const G4Event *evt)
{
}
void SToGS::AsciiTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
}
void SToGS::AsciiTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
}
void SToGS::AsciiSteppingAction::UserSteppingAction(const G4Step* aStep)
{
}

SToGS::Ascii::Ascii(G4String conf):
    UserActionInitialization(),
    fConffile(conf),
    fPathToData("./"),
    fBaseName("SToGS_Ascii"),
    fMaxEvents(500000),
    fRecordOption(0),
    fPrintModulo(1000)
{
    // open the ascii file
    std::ifstream file(fConffile.data());
	if ( file.is_open() == false ) {
		G4cout << " ** SToGS::Ascii WARNING ** Cannot open file " << fConffile << ", default parameters to be used "<< G4endl;
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
SToGS::Ascii::~Ascii()
{
    ;
}
G4UserRunAction *SToGS::Ascii::GetRunAction() const
{
    return new SToGS::AsciiRunAction(fPathToData,fBaseName,fMaxEvents,fRecordOption);
}
G4UserEventAction *SToGS::Ascii::GetEventAction() const
{
    return new SToGS::AsciiEventAction(fPrintModulo);
}
G4UserTrackingAction *SToGS::Ascii::GetTrackingAction() const
{
    return new AsciiTrackingAction();
}
G4UserSteppingAction *SToGS::Ascii::GetSteppingAction() const
{
    return 0x0;
}

void SToGS::Ascii::BuildForMaster () const
{
#if G4VERSION_NUMBER < 1000
    G4cout << " *** ERROR *** SToGS::Ascii::BuildForMaster() should never be called by this version of Geant4 " << G4endl;
#else
    SetUserAction( new SToGS::AsciiRunAction(fPathToData,fBaseName,fMaxEvents,fRecordOption) );
#endif
}
void SToGS::Ascii::Build () const
{
#if G4VERSION_NUMBER < 1000
    G4cout << " *** ERROR *** SToGS::Ascii::BuildForMaster() should never be called by this version of Geant4 " << G4endl;
#else
    SetUserAction( GetGun(fWhichGenerator.first,fWhichGenerator.second) );
    SetUserAction( new SToGS::AsciiRunAction(fPathToData,fBaseName,fMaxEvents,fRecordOption) ) ;
    SetUserAction( new SToGS::AsciiEventAction(fPrintModulo) ) ;
    SetUserAction( new SToGS::AsciiTrackingAction() ) ;
    SetUserAction( new SToGS::AsciiSteppingAction() ) ;
#endif
}







