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
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"

#include <vector>

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
    //! This class just print out once a new run begins/ends with the run number and the number of events to be simulated
    /*!
     */
    class PrintOutRunAction : public G4UserRunAction
    {
    public:
        PrintOutRunAction();
        virtual ~PrintOutRunAction();
        
        virtual G4Run* GenerateRun()
        {
            return new PrintOutRun();
        }
        virtual void BeginOfRunAction(const G4Run *aRun);
        virtual void EndOfRunAction(const G4Run* aRun);
        
    };
    //! This class just print out once a new event begins/ends with the event
    /*!
     */
    class PrintOutEventAction : public G4UserEventAction
    {
    public:
        PrintOutEventAction()
            {;}
        virtual ~PrintOutEventAction()
            {;}
        
    public:
        virtual void BeginOfEventAction(const G4Event *event);
        virtual void EndOfEventAction(const G4Event *event);
    };
    //! This class just print out once a new track begins/ends with some informations PARTICLE, TRACKID, PARENTID, Total and kinetic energy
    /*!
     */
    class PrintOutTrackingAction : public G4UserTrackingAction
    {
    public:
        PrintOutTrackingAction()
            {;}
        virtual ~PrintOutTrackingAction()
            {;}
        
    public:
        virtual void PreUserTrackingAction(const G4Track* atrack);
        virtual void PostUserTrackingAction(const G4Track *atrack);
    };
    //! This class just print out information of step once going from one volume to another one
    /*!
     */
    class PrintOutSteppingAction : public G4UserSteppingAction
    {
    public:
        PrintOutSteppingAction() : G4UserSteppingAction()
            {;}
        virtual ~PrintOutSteppingAction()
            {;}
        
    public:
        virtual void UserSteppingAction(const G4Step *step);
    };
    
    //! Extract informations from Geant4 and print on the standard output hits informations
    /*!
     This class illustrates how to extract some informations from Geant4 and display them
     on the standard output.
     
     It shows first how to define the sensitive detectors (part of the detector for
     which collects hits in the geometry) and how to implements an EndOfEvent action.
     
     The sensitive part are defined in the constructor.
     */
    class PrintOut : public UserActionInitialization
    {
    protected:
        G4String fOption;
        
    public:
        PrintOut(G4String option = "run"): UserActionInitialization(), fOption(option)
            {;}
        virtual ~PrintOut();
        
        //! return a PrintOutRunAction action
        virtual G4UserRunAction *GetRunAction() const;
        virtual G4UserEventAction *GetEventAction() const;
        virtual G4UserTrackingAction *GetTrackingAction() const;
        virtual G4UserSteppingAction *GetSteppingAction() const;
        
        virtual void 	BuildForMaster () const;
        virtual void 	Build () const;
        };
        
        } // SToGS Namespace
        
#endif   /* SToGS_PrintOut.hh */
        
        
        
        
