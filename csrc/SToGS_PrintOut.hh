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

#ifndef SToGS_PrintOut_h
#define SToGS_PrintOut_h 1

#include "SToGS_UserActionInitialization.hh"
#include "G4Run.hh"

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //!
    /*!
     */
    class PrintOutRun : public G4Run
    {
    protected:
        G4int colltrackerID;
        G4int collcaloID;
        
    public:
        PrintOutRun();
        virtual ~PrintOutRun();
        
        virtual void RecordEvent(const G4Event* evt);
        //virtual void Merge(const G4Run*);
    };
    //! At any steps, it prints out some informations regarding run, event, step ...
    /*!
     */
    class PrintOutAction : public AllActions
    {
    public:
        PrintOutAction(G4String option = "") : AllActions(option)
            {;}
        virtual ~PrintOutAction()
            {;}
        
        virtual G4Run* GenerateRun();
        
        virtual void BeginOfRunAction(const G4Run * /*therun*/);
        virtual void EndOfRunAction(const G4Run * /*therun*/);
        virtual void BeginOfEventAction(const G4Event * /*event*/);
        virtual void EndOfEventAction(const G4Event * /*event*/);
        virtual void PreUserTrackingAction(const G4Track * /*atrack*/);
        virtual void PostUserTrackingAction(const G4Track * /*atrack*/);
        virtual void UserSteppingAction(const G4Step * /*step*/);
    };
    //! This is the ActionInitialization class to handle PrintOut at different levels (run, events ...) from a unique AllAction event Class
    /*!
     */
    class PrintOut : public AllInOneUserActionInitialization<PrintOutAction>
    {
    public:
        PrintOut(G4String option = "run;event;track;step", G4String which_gene = "GPS", G4String which_gene_opt = "G4Macros/GPSPointLike.mac") :
            AllInOneUserActionInitialization<PrintOutAction>(option,which_gene,which_gene_opt)
            {;}
        virtual ~PrintOut()
            {;}
        
        virtual G4UserRunAction *GetRunAction() const;
        virtual G4UserEventAction *GetEventAction() const;
        virtual G4UserTrackingAction *GetTrackingAction() const;
        virtual G4UserSteppingAction *GetSteppingAction() const;
        
        virtual void 	Build () const;
    };

    
} // SToGS Namespace
        
#endif   /* SToGS_PrintOut.hh */
        
        
        
        
