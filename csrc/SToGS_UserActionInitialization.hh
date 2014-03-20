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

#ifndef SToGS_UserActionInitialization_h
#define SToGS_UserActionInitialization_h 1

#include "G4Version.hh"
#include "G4String.hh"

#if G4VERSION_NUMBER < 1000
#else
#include "G4VUserActionInitialization.hh"
#endif

#include <utility>

class G4VSensitiveDetector;
class G4VUserPrimaryGeneratorAction;
class G4UserRunAction;
class G4UserEventAction;
class G4UserTrackingAction;
class G4UserSteppingAction;

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! Base class to manage SToGS user's actions.
    /*!
     
        This class is mandatory in case of Multi-Threading. It should implement Build and BuidForMaster methods for G4 Kernel.

        All classes that inherits from it should focus on the User's Actions i.e Run, Event, Stepping etc ...
        As Geant4.10 required to put toghether all these actions with the generator, 
        there are facilities in the base class to select a general pre-existing one (like GPS for instance)
        or the event generator could be proper to the class. In the case, just overwrite GetGun()
     
     */
#if G4VERSION_NUMBER < 1000
    class UserActionInitialization
#else
    class UserActionInitialization : public G4VUserActionInitialization
#endif
    {
    protected:
        //! to keep what generator to be used. Defaut is "-","-", means the concrete implementation provides its own generator
        std::pair < G4String, G4String > fWhichGenerator;

    protected:
        //! Find and get a gun in many possibilities
        virtual G4VUserPrimaryGeneratorAction *GetGun(G4String which, G4String opt) const;
        
    public:
        UserActionInitialization();
        virtual ~UserActionInitialization();
        
        //! once allocated the generator associated with a UserActionInititialization could be changed here before asking for one 
        /*!
         */
        std::pair < G4String, G4String > SetWhichGenerator(G4String which_gene, G4String option);
        
        //! to get a general SToGS tracker. In Multi-threading mode, return a new instance otherwise a global one
        static G4VSensitiveDetector *GetTrackerSD( G4String name = "/SToGS/Track" );
        //! to get a general SToGS Calo. In Multi-threading mode, return a new instance otherwise a global one
        static G4VSensitiveDetector *GetCaloSD( G4String name = "/SToGS/Calo" );
        
// INTERFACE provided for Geant4 < 10.0
    public:
        //! depending on one string, select a given gun
        /*!
            GPS: the gun is realized through the general GPS interface
         */
        virtual G4VUserPrimaryGeneratorAction *GetGun() const;
        //! depending on one string, select a given gun
        /*!
         */
        virtual G4UserRunAction *GetRunAction() const = 0;
        virtual G4UserEventAction *GetEventAction() const = 0;
        virtual G4UserTrackingAction *GetTrackingAction() const = 0;
        virtual G4UserSteppingAction *GetSteppingAction() const = 0;
        
// INTERFACE in principle required for Geant4 > 10.0
        virtual void 	BuildForMaster () const = 0;
        virtual void 	Build () const = 0;
    };
} // SToGS Namespace

#endif


