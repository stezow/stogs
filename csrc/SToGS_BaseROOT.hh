
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
//----------------------------------------------------------------------------

#ifndef SToGS_BaseROOT_h
#define SToGS_BaseROOT_h 1

#include "SToGS_UserActionInitialization.hh"
#include "G4Run.hh"

class TFile;
class TTree;

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //!
    /*!
     */
    class BaseROOTTreeRun : public G4Run
    {
    protected:
        G4int fRecordOption;
    protected:
        G4int colltrackerID;
        G4int collcaloID;
  
    public:
        BaseROOTTreeRun(G4int record_option = 0);
        virtual ~BaseROOTTreeRun();
        
        //virtual void RecordEvent(const G4Event* evt);
        //virtual void Merge(const G4Run*);
    };
    //! This class just print out once a new run begins/ends with the run number and the number of events to be simulated
    /*!
     \TODO
     */
    class BaseROOTAction : public AllActions
    {
    protected:
        //! the current root file
        TFile *fRootFile;
    protected:
        //! directory where to output data
        G4String fPathToData;
        //! base for all the files
        G4String fBaseName;
        //! max numer of event per files ... better to limit because of BaseROOT file could be uged !
        G4int fMaxEvents;
        //! 0 [default] means keep all, 1 only events which gives at least one hit in the full detector
        G4int fRecordOption;
        //! to print out status any fPrintModulo events
        G4int fPrintModulo;
    protected:
        //! Just check if there are collected hits in the collection
        // std::pair<G4int, G4int> HitsinCollection(const G4Event);
        
        //! Open the stream depending of the configuration
        virtual void OpenFile(G4int run_id);
        //! Make sure ths file is closed properly
        virtual void CloseFile();
        
    public:
        BaseROOTAction(G4String conffile = "setup/SToGS_root_actions.conf");
        virtual ~BaseROOTAction()
            {;}
                
        virtual void BeginOfRunAction(const G4Run * /*therun*/);
        virtual void EndOfRunAction(const G4Run * /*therun*/);
        virtual void BeginOfEventAction(const G4Event * /*event*/);
        virtual void EndOfEventAction(const G4Event * /*event*/);
        virtual void PreUserTrackingAction(const G4Track * /*atrack*/)
            {;}
        virtual void PostUserTrackingAction(const G4Track * /*atrack*/)
            {;}
        // virtual void UserSteppingAction(const G4Step * /*step*/);
    };
    class BaseROOTTreeAction : public BaseROOTAction
    {
    protected:
        //! the Tree to be filled
        TTree *fTree;
        
    protected:
        G4String fTreeName;
        G4String fTreeTitle;
        
    protected:
        //! Open the stream depending of the configuration and attach the Tree to the file
        virtual void OpenFile(G4int run_id);
        //! Make sure ths file is closed properly
        virtual void CloseFile();
        
    public:
        BaseROOTTreeAction(G4String conffile = "setup/SToGS_root_actions.conf");
        virtual ~BaseROOTTreeAction();
    };
    
    //! Extract informations from Geant4 using SToGS sensitives and write hits in root tree
    /*!
     */
    class BaseROOTTree : public  AllInOneUserActionInitialization<BaseROOTTreeAction>
    {
    public:
        BaseROOTTree(G4String conf = "setup/SToGS_root_actions.conf", G4String which_gene = "GPS", G4String which_gene_opt = "G4Macros/GPSPointLike.mac"):
        AllInOneUserActionInitialization<BaseROOTTreeAction>(conf,which_gene,which_gene_opt)
            {;}
        virtual ~BaseROOTTree()
            {;}
    };
    
} // SToGS Namespace


#endif
