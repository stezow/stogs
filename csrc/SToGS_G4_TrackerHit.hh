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

#ifndef SToGS_TrackerHit_h
#define SToGS_TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Types.hh"

#ifndef G4MULTITHREADED
#define G4ThreadLocal
#endif

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! Informations to be kept at each step if a positive energy is deposited in a sensitive detector
    /*!
     In G4, one have to define a class (a hit) that inherits from G4VHit to keep informations during the tracking.
     This one is dedicated for tracking since it kept infomations of any single impact in sensitive detectors. Here is
     the list of the available informations for each TrackerHit. Use Get/Set methods to obtain/change these values.
     
     - G4int trackID:   Track to which belongs this hit
     - G4int parentID:  Parent track of the track to which this hit belongs
     - G4int primaryID: Primary vertex from which this hit is coming
     
     - G4int detID: Detector ID in which the hit occured
     - G4int motherID: Mother ID of the detector in which the hit occured
     
     - G4double edep: Energy deposited
     - G4double ToF:  global time at which the interaction occured
     
     - G4ThreeVector pos: Position of the impact
     
     - G4String processName: Process that gives this hit
     - G4String particleName: Name of the particle that lost energie at that position
     - G4String detName: Name of the detector in which the hit occured
     - G4String motherDetName: CURRENTLY NOT AVAILABLE
     - G4int PDGcode: code of the particle
     (Nuclear codes are given as 10-digit numbers +-100ZZZAAAI.
     I gives the isomer level, with I = 0 corresponding
     to the ground state and I >0 to excitations
     for the elementary particle see http://www.slac.stanford.edu/BFROOT/www/Computing/Environment/NewUser/htmlbug/node51.html
     These codes are define in G4ParticleTable
     )
     */
    class TrackerHit : public G4VHit
    {
    public:
        TrackerHit();
        virtual ~TrackerHit();
        TrackerHit(const TrackerHit &right);
        
        const TrackerHit& operator=(const TrackerHit &right);
        G4int operator==(const TrackerHit &right) const;
        
        inline void *operator new(size_t);
        inline void operator delete(void *aHit);
        
    public:
// setters
        void SetTrackID(G4int id)
        {
            trackID = id;
        }
        void SetParentID(G4int id)
        {
            parentID = id;
        }
        void SetPrimaryID(G4int id)
        {
            primaryID = id;
        }
        
        void SetDetID  (G4int id)
        {
            detID = id;
        }
        void SetMotherID(G4int id)
        {
            motherID = id;
        }
        void SetEdep(G4double de)
        {
            edep = de;
        }
        void SetToF (G4double tf)
        {
            ToF = tf;
        }
        void SetPos(G4ThreeVector xyz)
        {
            pos = xyz;
        }
        void SetProcessName(G4String name)
        {
            processName = name;
        }
        void SetParticleName(G4String name)
        {
            particleName = name;
        }
        void SetDetName(G4String name)
        {
            detName = name;
        }
        void SetMotherDetName(G4String name)
        {
            motherDetName = name;
        }
        void SetPDGcode(G4int code)
        {
            PDGcode = code;
        }
 // getters
        G4int GetTrackID()
        {
            return trackID;
        }
        G4int GetParentID()
        {
            return parentID;
        }
        G4int GetPrimaryID()
        {
            return primaryID;
        }
        G4int GetDetID()
        {
            return detID;
        }
        G4int GetMotherID()
        {
            return motherID;
        }
        G4double GetEdep()
        {
            return edep;
        }
        G4double GetToF()
        {
            return ToF;
        }
        G4int  GetPDGcode()
        {
            return PDGcode;
        }
        G4ThreeVector &GetPos()
        {
            return pos;
        }
        G4String &GetProcessName()
        {
            return processName;
        }
        G4String &GetParticleName()
        {
            return particleName;
        }
        G4String &GetDetName()
        {
            return detName;
        }
        G4String &GetMotherDetName()
        {
            return motherDetName;
        }
        
// others
        void Draw();
        void Print();
        
    private:
        G4int trackID;	// current track ID to which this hit belongs
        G4int parentID;	//  parent track ID to which this hit belongs
        G4int primaryID;	// primary track ID to which this hit belongs
        
        G4int detID;
        G4int motherID;
        
        G4double edep;
        G4double ToF;
        
        G4ThreeVector pos;
        
        G4String processName;
        G4String particleName;
        G4String detName;
        G4String motherDetName;
        
        
        G4int PDGcode;	// PDGcode of particle which hit the detector
    };
    
    typedef G4THitsCollection<TrackerHit> TrackerHitsCollection;
    extern G4ThreadLocal G4Allocator<TrackerHit> *TrackerHitAllocator;
    
#ifdef G4MULTITHREADED
//    extern G4ThreadLocal G4Allocator<TrackerHit> *TrackerHitAllocator;
#else
//    extern G4Allocator<TrackerHit> *TrackerHitAllocator;
#endif
    
    inline void* TrackerHit::operator new(size_t)
    {
        if( !TrackerHitAllocator )
            TrackerHitAllocator = new G4Allocator<TrackerHit>;
        return TrackerHitAllocator->MallocSingle();
    }
    
    inline void TrackerHit::operator delete(void *aHit)
    {
        TrackerHitAllocator->FreeSingle((TrackerHit*) aHit);
    }
} // SToGS Namespace

#endif
