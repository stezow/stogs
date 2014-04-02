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
#include "G4ios.hh"
#include "G4String.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4UserSteppingAction.hh"

#if G4VERSION_NUMBER < 1000
#else
#include "G4VUserActionInitialization.hh"
#endif

#include <vector>
#include <utility>

class G4VSensitiveDetector;
class G4VUserPrimaryGeneratorAction;


//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! Base class that regroups in the same space all user's hooks. Convenient for sharing similar data
    /*!
        Having such requires carefull manipulation since it trickes Geant4 which expect different objects.
        In partilcular, User's actions are deleting by the Geant4 kernel.
     
        This is why SToGS defines also SToGS::RunAction (etc ...) which are bridges like design patterns
     
        The option could be used to pass a filename that contains configuration fields
     */
    class AllActions : public G4UserRunAction, public G4UserEventAction, public G4UserTrackingAction, public G4UserSteppingAction
    {
    protected:
        G4String fOption;
        
    public:
        AllActions(G4String opt): G4UserRunAction(), G4UserEventAction(), G4UserTrackingAction(), G4UserSteppingAction(), fOption(opt)
            {;}
        virtual ~AllActions()
            {;}
        
        virtual G4Run* GenerateRun()
        {
            return 0x0;
        }
        virtual void BeginOfRunAction(const G4Run * /*therun*/)
        {
            ;
        }
        virtual void EndOfRunAction(const G4Run * /*therun*/)
        {
            ;
        }
        virtual void BeginOfEventAction(const G4Event * /*event*/)
        {
            ;
        }
        virtual void EndOfEventAction(const G4Event * /*event*/)
        {
            ;
        }
        virtual void PreUserTrackingAction(const G4Track * /*atrack*/)
        {
            ;
        }
        virtual void PostUserTrackingAction(const G4Track * /*atrack*/)
        {
            ;
        }
        virtual void UserSteppingAction(const G4Step * /*step*/)
        {
            ;
        }
    };
    //! Base class for a Run action that calls a concrete one.
    /*!
     */
    class RunAction : public G4UserRunAction
    {
    private:
        G4UserRunAction  *theRealAction;
    public:
        RunAction(G4UserRunAction *theaction = 0x0) : G4UserRunAction(), theRealAction(theaction)
            {;}
        virtual ~RunAction()
            {;}
        
    public:
        virtual G4Run* GenerateRun()
        {
            if (theRealAction)
                return theRealAction->GenerateRun();
            return G4UserRunAction::GenerateRun();
        }
        virtual void BeginOfRunAction(const G4Run *aRun)
        {
            if (theRealAction)
                theRealAction->BeginOfRunAction(aRun);
        }
        virtual void EndOfRunAction(const G4Run* aRun)
        {
            if (theRealAction)
                theRealAction->EndOfRunAction(aRun);
        }
        
    };
    //! Base class for a Event action that calls a concrete one.
    /*!
     */
    class EventAction : public G4UserEventAction
    {
    private:
        G4UserEventAction *theRealAction;
    public:
        EventAction(G4UserEventAction *theaction = 0x0) : G4UserEventAction(), theRealAction(theaction)
            {;}
        virtual ~EventAction()
            {;}
    public:
        virtual void BeginOfEventAction(const G4Event *event)
        {
            if (theRealAction)
                theRealAction->BeginOfEventAction(event);
        }
        virtual void EndOfEventAction(const G4Event *event)
        {
            if (theRealAction)
                theRealAction->EndOfEventAction(event);
        }
    };
    //! Base class for a Tracking action that calls a concrete one.
    /*!
     */
    class TrackingAction : public G4UserTrackingAction
    {
    private:
        G4UserTrackingAction *theRealAction;
    public:
        TrackingAction(G4UserTrackingAction *theaction = 0x0) : G4UserTrackingAction(), theRealAction(theaction)
            {;}
        virtual ~TrackingAction()
            {;}
    public:
        virtual void PreUserTrackingAction(const G4Track* atrack)
        {
            if (theRealAction)
                theRealAction->PreUserTrackingAction(atrack);
        }
        virtual void PostUserTrackingAction(const G4Track *atrack)
        {
            if (theRealAction)
                theRealAction->PostUserTrackingAction(atrack);
        }
    };
    //! Base class for a Steeping action that calls a concrete one.
    /*!
     */
    class SteppingAction : public G4UserSteppingAction
    {
    private:
        G4UserSteppingAction *theRealAction;
    public:
        SteppingAction(G4UserSteppingAction *theaction = 0x0) : G4UserSteppingAction(), theRealAction(theaction)
            {;}
        virtual ~SteppingAction()
            {;}
    public:
        virtual void UserSteppingAction(const G4Step *step)
        {
            if (theRealAction)
                theRealAction->UserSteppingAction(step);
        }
    };

    //! Base class to manage SToGS user's hooks + the generator
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
        //! list of allocated AllActions to delete them
        mutable std::vector < AllActions *> fAllUserAction;
        
    protected:
        //! action passed to the created user's actions if required
        G4String fUserActionOption;
        
    protected:
        //! to keep what generator to be used. Defaut is "-","-", means the concrete implementation provides its own generator
        std::pair < G4String, G4String > fWhichGenerator;

    protected:
        //! Find and get a gun in many possibilities
        virtual G4VUserPrimaryGeneratorAction *GetGun(G4String which, G4String opt) const;
        //! add the passed class to the list of allocated ones. It is required for a proper destruction
        /*!
            This call is thread locked
         */
        void BuildAndRegister ( AllActions* ) const;
        void Register ( AllActions* ) const;
     
    public:
        UserActionInitialization(G4String which_user_action_opt, G4String which_gene, G4String which_gene_opt);
        virtual ~UserActionInitialization();
        
        //! once allocated the generator associated with a UserActionInititialization could be changed here before asking for one 
        /*!
         */
        std::pair < G4String, G4String > SetWhichGenerator(G4String which_gene, G4String option);
        
        //! to get a general SToGS tracker. In Multi-threading mode, return a new instance otherwise a global one
        static G4VSensitiveDetector *GetTrackerSD( G4String name = "/SToGS/SD/Tracker" );
        //! to get a general SToGS Calorimeter. In Multi-threading mode, return a new instance otherwise a global one
        static G4VSensitiveDetector *GetCopClusterSD( G4String name = "/SToGS/SD/CopCluster" );
        
// INTERFACE provided for Geant4 < 10.0
    public:
        //! depending on one string, select a given gun
        /*!
            Ex GPS: the gun is realized through the general GPS interface
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
        //! provide default ... likely to be changed since RunManager may be different from Slave and Master
        virtual void 	BuildForMaster () const
        {
            G4cout << " ------ INF ------ from SToGS::AllInOneUserActionInitialization::BuildForMaster() " << G4endl;
#if G4VERSION_NUMBER < 1000
            G4cout << " *** ERROR *** SToGS::AllInOneUserActionInitialization::BuildForMaster() should never be called by this version of Geant4 " << G4endl;
#else
            // As real work is done by slaves, default is not RunManager for Master. Could be changed by overloading this method
            SetUserAction( new SToGS::RunAction(0x0) ) ;
#endif
            G4cout << " ------ END ------  from SToGS::AllInOneUserActionInitialization::BuildForMaster() " << G4endl;
        }
        virtual void 	Build () const = 0;
    };
    
    //! a G4 user's action manage by a single AllAction class
    template<typename _T> class AllInOneUserActionInitialization : public UserActionInitialization
    {
    public:
        AllInOneUserActionInitialization(G4String which_user_action_opt, G4String which_gene, G4String which_gene_opt) :
            UserActionInitialization(which_user_action_opt,which_gene,which_gene_opt)
        {
            ;
        }
        virtual ~AllInOneUserActionInitialization()
        {
            ;
        }
        virtual G4UserRunAction *GetRunAction() const
        {
            if ( fAllUserAction.size() == 0 ) {
                fAllUserAction.push_back( new _T(fUserActionOption) );
            }
            return new SToGS::RunAction(fAllUserAction[0]);
        }
        virtual G4UserEventAction *GetEventAction() const
        {
            if ( fAllUserAction.size() == 0 ) {
                fAllUserAction.push_back( new _T(fUserActionOption) );
            }
            return new SToGS::EventAction(fAllUserAction[0]);
        }
        virtual G4UserTrackingAction *GetTrackingAction() const
        {
            if ( fAllUserAction.size() == 0 ) {
                fAllUserAction.push_back( new _T(fUserActionOption) );
            }
            return new SToGS::TrackingAction(fAllUserAction[0]);
        }
        virtual G4UserSteppingAction *GetSteppingAction() const
        {
            if ( fAllUserAction.size() == 0 ) {
                fAllUserAction.push_back( new _T(fUserActionOption) );
            }
            return new SToGS::SteppingAction(fAllUserAction[0]);
        }
        virtual void Build () const
        {
            BuildAndRegister( new _T(fUserActionOption) );
        }
    };
    
    
} // SToGS Namespace

#endif


