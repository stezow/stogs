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

#include "SToGS_UserActionManager.hh"
#include "G4Run.hh"

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    
    class PrintRun : public G4Run
    {
    public:
        PrintRun();
        virtual ~PrintRun();
        
        virtual void RecordEvent(const G4Event* evt);
        virtual void Merge(const G4Run*);
    };
    
    //! Extract informations from Geant4 and print on the standard output hits informations
    /*!
     This class illustrates how to extract some informations from Geant4 and display them
     on the standard output.
     
     It shows first how to define the sensitive detectors (part of the detector for
     which collects hits in the geometry) and how to implements an EndOfEvent action.
     
     The sensitive part are defined in the constructor.
     */
    class PrintOut : public UserActionManager
    {
    public:
        PrintOut(G4String conf = "setup/printout");
        virtual ~PrintOut();
        
        virtual void 	BuildForMaster () const;
        virtual void 	Build () const;
    };
    
} // SToGS Namespace

#endif   /* ParisPrintOut.hh */




