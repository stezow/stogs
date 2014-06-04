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

#ifndef SToGS_BaseROOTEventsActions_h
#define SToGS_BaseROOTEventsActions_h 1

#include "SToGS_BaseROOT.hh"
#include "SToGS_BaseROOTEvents.h"

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! A Run
    /*!
     */
    class BaseROOTEventsRun : public BaseROOTTreeRun
    {
    protected:
        TTree *fTree;
    protected:
        // Events to be filled from G4 to ROOT
        SBRPEvent fPrimaryEvent;
        SBREvent  fEvent;
        
    public:
        BaseROOTEventsRun(TTree *tree);
        virtual ~BaseROOTEventsRun()
        {;}
        
        virtual void RecordEvent(const G4Event* evt);
    };

    //! The BaseROOTEventsUserAction Stores in a ROOT Tree BaseROOTEvents [SBREvent, SBRPEvent etc ... see analysis/SToGS/BaseROOTEvents.h file]
    /*!
     */
    class BaseROOTEventsUserAction : public SToGS::BaseROOTTreeAction
    {
    protected:
        //! true if one has also photons from scintillations
        G4int fisOptical;
        //! emitted optical photons
        SBROpticalEvent *fOpticalEventBeg;
        //! end of optical photons
        SBROpticalEvent *fOpticalEventEnd;
        
    public:
        BaseROOTEventsUserAction(G4String conffile = "setup/SToGS_Tree_actions.conf");
        virtual ~BaseROOTEventsUserAction()
        {
            delete fOpticalEventBeg;
            delete fOpticalEventEnd;
        }
        
        virtual G4Run* GenerateRun();
        
        //    virtual void BeginOfRunAction(const G4Run * /*therun*/);
        virtual void EndOfRunAction(const G4Run * /*therun*/);
        virtual void BeginOfEventAction(const G4Event * /*event*/);
        //    virtual void EndOfEventAction(const G4Event * /*event*/);
        virtual void PreUserTrackingAction(const G4Track * /*atrack*/);
        virtual void PostUserTrackingAction(const G4Track * /*atrack*/);
    };

    //! The BaseROOTEvents is use to init G4 kernel with the actions defined in BaseROOTEventsUserAction
    /*!
     */
    class BaseROOTEvents : public SToGS::AllInOneUserActionInitialization<BaseROOTEventsUserAction>
    {
    public:
        BaseROOTEvents(G4String conf = "setup/SToGS_tree_actions.conf", G4String which_gene = "GPS",
                                      G4String which_gene_opt = "G4Macros/GPSPointLike.mac"):
            AllInOneUserActionInitialization<BaseROOTEventsUserAction>(conf,which_gene,which_gene_opt)
        {;}
        virtual ~BaseROOTEvents()
        {;}
    };
    
} // SToGS Namespace

#endif

