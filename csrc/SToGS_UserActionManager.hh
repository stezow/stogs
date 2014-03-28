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

#include "SToGS_UserActionInitialization.hh"

class G4VUserDetectorConstruction;
class G4VUserPhysicsList;

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //!  Manager which, based on a configuration file, manages user's actions + selection of the geometry and the physics list
    /*!
        It extends it also by adding the physics list and the geometry.
        It can be used in geant4.9.6 and geant > 10.0
     */
    class UserActionManager : public UserActionInitialization
    {
    private:
        SToGS::UserActionInitialization *fImplementation;
        
    protected:
        //! In principle in Geant4.10, geometry and physics are golbal so they are not provided by G4VUserActionInitialization
        /*!
            We are keeping however here some references in case the user would like to ask the manager to do the job for her/him
         */
        std::pair < G4String, G4String > fWhichGeometry;
        std::pair < G4String, G4String > fWhichPhysics;
        
        //! Real stuff the G4VUserActionInitialization should deal with
        std::pair < G4String, G4String > fWhichActionManager;
        
    protected:
        G4int fNbThreads;
        
    protected:
        SToGS::UserActionInitialization *GetUserActionInitialization();
        
    public:
        UserActionManager(G4String configurationfile = "");
        virtual ~UserActionManager()
            {;}
    public:
        //! allows the user to specify the number of threads
        G4int GetNbThread() const
        {
            return fNbThreads;
        }
        //! Return the detector depending of the configuration file
        virtual G4VUserDetectorConstruction *GetDetectorConstruction() const;
        //! Return the physics list depending of the configuration file
        virtual G4VUserPhysicsList *GetPhysicsList() const;
        
        //! Individual calls in case it is used for Geant4 < 10.0
        virtual G4UserRunAction *GetRunAction() const;
        virtual G4UserEventAction *GetEventAction() const;
        virtual G4UserTrackingAction *GetTrackingAction() const;
        virtual G4UserSteppingAction *GetSteppingAction() const;
        
        //! Geant4 > 10.0 interface to set user's local actions i.e Gun, Run, Event etc ... Action
        virtual void 	BuildForMaster () const;
        virtual void 	Build () const;
    };
} // SToGS Namespace

#endif


