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

#ifndef SToGS_Ascii_h
#define SToGS_Ascii_h 1

#include "SToGS_UserActionInitialization.hh"
#include "G4Run.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"

#include <vector>
#include <fstream>

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //!
    /*!
     */
    class AsciiRun : public G4Run
    {
    protected:
        G4int colltrackerID;
        G4int collcaloID;
    protected:
        //! current stream to output data
        std::ofstream &fOutputFile;
    protected:
        //! A new event in the file starts/ends with the following characters
        //std::pair < char, char >  fEventMarkOff;
        
    public:
        AsciiRun(std::ofstream &out);
        virtual ~AsciiRun();

        virtual void RecordEvent(const G4Event* evt);
        //virtual void Merge(const G4Run*);
    };
    //! This class just print out once a new run begins/ends with the run number and the number of events to be simulated
    /*!
     \TODO
     */
    class AsciiAction : public AllActions
    {
    protected:
        //! current stream to output data
        std::ofstream fOutputFile;
    protected:
        //! directory where to output data
        G4String fPathToData;
        //! base for all the files
        G4String fBaseName;
        //! max numer of event per files ... better to limit because of ascii file could be uged !
        G4int fMaxEvents;
        //! 0 [default] means keep all, 1 only events which gives at least one hit in the full detector
        G4int fRecordOption;
        //! to print out status any fPrintModulo events
        G4int fPrintModulo;
        
    protected:
        //! Just check if there are collected hits in the collection
        // std::pair<G4int, G4int> HitsinCollection(const G4Event);
        
        //! Open the stream depending of the configuration
        void OpenStream(G4int run_id);
        
    public:
        AsciiAction(G4String conffile = "setup/SToGS_ascii_actions.conf");
        virtual ~AsciiAction()
            {;}
        
        virtual G4Run* GenerateRun();
        
        virtual void BeginOfRunAction(const G4Run * /*therun*/);
        virtual void EndOfRunAction(const G4Run * /*therun*/);
        virtual void BeginOfEventAction(const G4Event * /*event*/);
        virtual void EndOfEventAction(const G4Event * /*event*/);
        // virtual void PreUserTrackingAction(const G4Track * /*atrack*/);
        // virtual void PostUserTrackingAction(const G4Track * /*atrack*/);
        // virtual void UserSteppingAction(const G4Step * /*step*/);
    };
    
    //! Extract informations from Geant4 using SToGS sensitives and write hits in ascii files
    /*!
     */
    class Ascii : public  AllInOneUserActionInitialization<AsciiAction>
    {
    public:
        Ascii(G4String conf = "setup/SToGS_ascii_actions.conf", G4String which_gene = "GPS", G4String which_gene_opt = "G4Macros/GPSPointLike.mac"):
            AllInOneUserActionInitialization<AsciiAction>(conf,which_gene,which_gene_opt)
            {;}
        virtual ~Ascii()
            {;}
        
        // virtual void Build () const;
    };
        
} // SToGS Namespace
        
#endif   /* SToGS_Ascii.hh */
        
        
        
        
