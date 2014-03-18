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

#ifndef SToGS_UserActionManager_h
#define SToGS_UserActionManager_h 1

// TOBE DONE : call it UserActionManager

//G4Version.hh
#include "G4String.hh"

#if G4VERSION_NUMBER < 1000
#else
#include "G4VUserUserActionManager.hh"
#endif

#include <utility>

class G4VSensitiveDetector;
class G4VUserPrimaryGeneratorAction;

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! Base class for SToGS user's actions.
    /*!
     */
#if G4VERSION_NUMBER < 1000
    class UserActionManager
#else
    class UserActionManager : public G4VUserUserActionManager
#endif
    {
    private:
        std::pair < G4String, G4String > fWhichGeometry;
        std::pair < G4String, G4String > fWhichPhysics;
        
        std::pair < G4String, G4String > fWhichGenerator;
        std::pair < G4String, G4String > fWhichActionManager;
        
    protected:
        //! depending on one string, select a given gun
        /*!
            GPS: the gun is realized through the general GPS interface
         */
        G4VUserPrimaryGeneratorAction *GetGun(G4String which = "GPS", G4String opt = "G4Macros/GPS_Cs137.mac");
        
    public:
        UserActionManager(G4String file = "setup/SToGS.global");
        virtual ~UserActionManager();
        
        virtual void 	BuildForMaster () const;
        virtual void 	Build () const;
        
    public:
        //! to get a general SToGS tracker. In Multi-threading mode, return a new instance otherwise a global one
        static G4VSensitiveDetector *GetTrackerSD( G4String name = "/SToGS/Track" );
        
        //! to get a general SToGS Calo. In Multi-threading mode, return a new instance otherwise a global one
        static G4VSensitiveDetector *GetCaloSD( G4String name = "/SToGS/Calo" );
    };
    //! Concrete Base class for SToGS user's actions.
    /*!
     */
    /*
    class CUserActionManager : public UserActionManager
    {
    private:
        UserActionManager *fImplementation;
        
    public:
        CUserActionManager(G4String file = "");
        virtual ~CUserActionManager();
        
        virtual void 	BuildForMaster () const;
        virtual void 	Build () const;
    };
*/
    
} // SToGS Namespace

#endif


