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

#include "SToGS_BaseROOT.hh"

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#if G4MULTITHREADED
#include "G4Threading.hh"
#endif
#ifdef G4MULTITHREADED
#include "G4AutoLock.hh"
namespace { G4Mutex buildMutex = G4MUTEX_INITIALIZER; }
#endif

#include "TFile.h"
#include "TTree.h"

#include <fstream>
#include <iomanip>

SToGS::BaseROOTTreeRun::BaseROOTTreeRun(TTree *tree) :
    G4Run(),
    colltrackerID(-1),
    collcaloID(-1)
//    fOutputFile(out)
{
   	G4SDManager* SDman = G4SDManager::GetSDMpointer();
    if ( SDman ) {
        // keep the collection ID associated to this collection
        colltrackerID = SDman->GetCollectionID("TrackerHits");
        
        // keep the collection ID associated to this collection
        collcaloID = SDman->GetCollectionID("CopClusterHits");
    }
    //SetBranc
}

void SToGS::BaseROOTTreeRun::RecordEvent(const G4Event* evt)
{
    
}

SToGS::BaseROOTTreeRun::~BaseROOTTreeRun()
{
}

SToGS::BaseROOTAction::BaseROOTAction(G4String conf):
    AllActions(conf),
    fRootFile(0x0),
    fPathToData("./"),
    fBaseName("SToGS_Ascii"),
    fMaxEvents(1900000000),
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
/*
G4Run* SToGS::BaseROOTAction::GenerateRun()
{
    G4Run* therun = 0x0; G4cout << " In BaseROOTAction, Generate a new Run " << G4endl;
    
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
 */

void SToGS::BaseROOTAction::OpenFile(G4int run_id)
{
#ifdef G4MULTITHREADED
    G4AutoLock lock(&buildMutex); // better to be protected since root is doing some stuff in global
#endif
    
    // for MT, thread id used in the name of the file to avoid pb
    G4int thread_id = 0;
#if G4MULTITHREADED
    thread_id = G4Threading::G4GetThreadId();
#endif
    std::ostringstream filename;
    filename.clear();
    filename << fPathToData << "/" << fBaseName << "_Thread" << std::setfill('0') << std::setw(2) << thread_id
             << "_Run" << std::setw(3) << run_id << ".root";
    
    fRootFile = new TFile(filename.str().data(),"UPDATE");
}
void SToGS::BaseROOTAction::CloseFile()
{
#ifdef G4MULTITHREADED
    G4AutoLock lock(&buildMutex); // better to be protected since root is doing some stuff in global
#endif
    if ( fRootFile )
    {
        fRootFile->Close();
        delete fRootFile; fRootFile = 0x0;
    }
}

void SToGS::BaseROOTAction::BeginOfRunAction(const G4Run *aRun)
{
	G4cout << " In BaseROOTAction, Begin of Run " << aRun->GetRunID()
            << ", Number of events to be simulated " << aRun->GetNumberOfEventToBeProcessed() << G4endl;
    //
    OpenFile(aRun->GetRunID());
}
void SToGS::BaseROOTAction::EndOfRunAction(const G4Run* aRun)
{
 	G4cout << " In BaseROOTAction, End of Run " << aRun->GetRunID()
            << ", Number of simulated events " << aRun->GetNumberOfEvent() << G4endl;
    //
    CloseFile();
}
void SToGS::BaseROOTAction::BeginOfEventAction(const G4Event *evt)
{
    G4int evtNb = evt->GetEventID();
    if ( evtNb % fPrintModulo == 0 ) {
    	G4cout << " In BaseROOTAction, Begin of event: " << evtNb << G4endl;
    }
}
void SToGS::BaseROOTAction::EndOfEventAction(const G4Event *evt)
{
    G4int evtNb = evt->GetEventID();
    if ( evtNb % fPrintModulo == 0 ) {
    	G4cout << " In BaseROOTAction, End of event: " << evtNb << G4endl;
    }
}

SToGS::BaseROOTTreeAction::BaseROOTTreeAction(G4String conffile) :
    SToGS::BaseROOTAction(conffile),
    fTree(0x0),
    fTreeName("SToGSTree"),
    fTreeTitle("Tree produced by SToGS")
{
    // open the ascii file
    std::ifstream file(conffile.data());
	if ( file.is_open() == false ) {
		G4cout << " ** SToGS::BaseROOTTreeAction WARNING ** Cannot open file " << conffile << ", default parameters to be used "<< G4endl;
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
            
            if ( key == "tree:" ) { // get name and title. Title could have spaces ...
                decode >> fTreeName ;
                
                size_t pos_title = aline.find(fTreeName);
                if ( pos_title != std::string::npos ) {
                    pos_title +=
                        fTreeName.length();
                    fTreeTitle = aline.substr(pos_title);
                }
            }
            getline(file,aline);
        }
        file.close();
    }
}

SToGS::BaseROOTTreeAction::~BaseROOTTreeAction()
{
    ;
}

void SToGS::BaseROOTTreeAction::OpenFile(G4int run_id)
{
#ifdef G4MULTITHREADED
    G4AutoLock lock(&buildMutex); 
    // better to be protected since root is doing some stuff in global
    // could not call BaseROOTAction::CloseFile() otherwise double lock !!!
#endif
    // for MT, thread id used in the name of the file to avoid pb
    G4int thread_id = 0;
#if G4MULTITHREADED
    thread_id = G4Threading::G4GetThreadId();
#endif
    if ( fTree == 0x0 ) {
        fTree = new TTree(fTreeName.data(),fTreeTitle.data());
        fTree->SetDirectory(0x0);
    }
    
    std::ostringstream filename;
    filename.clear();
    filename << fPathToData << "/" << fBaseName << "_Thread" << std::setfill('0') << std::setw(2) << thread_id
                << "_Run" << std::setw(3) << run_id << ".root";
    
    fRootFile = new TFile(filename.str().data(),"UPDATE");
    if ( fRootFile->IsOpen() ) {
        G4cout << " The File " << filename.str() << " is open to record data " << G4endl;
        fTree->SetDirectory(fRootFile);
    }
}
void SToGS::BaseROOTTreeAction::CloseFile()
{
#ifdef G4MULTITHREADED
    G4AutoLock lock(&buildMutex);
    // better to be protected since root is doing some stuff in global
    // could not call BaseROOTAction::CloseFile() otherwise double lock !!!
#endif
	TFile *current_root = fTree->GetCurrentFile();
	if ( current_root ) {
        G4cout << " Closing File " << current_root->GetName() << " " << G4endl;
		current_root->Write();
        fTree->SetDirectory(0x0);
        delete current_root; fRootFile = 0x0; // this is not an error !
	}
    delete fTree; fTree = 0x0;
}





