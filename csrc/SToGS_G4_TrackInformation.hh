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

#ifndef SToGS_G4_TrackInformation_h
#define SToGS_G4_TrackInformation_h 1

#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

#ifndef G4MULTITHREADED
#define G4ThreadLocal
#endif

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! User Track Information that kept the origine of any track
    /*!
     */
    class PrimaryTrackInformation : public G4VUserTrackInformation
    {
    private:
        //! track ID of the particle which is the origine of this track
        G4int fPrimaryID;
    public:
        PrimaryTrackInformation(G4int pid = 0) : G4VUserTrackInformation(), fPrimaryID(pid)
            {;}
        virtual ~PrimaryTrackInformation()
            {;}
        
    public:
        //! Set the primary ID
        void SetPrimaryID(G4int pid = 0)
        {
            fPrimaryID = pid;
        }
        //! Get the primary ID
        G4int GetPrimaryID() const
        {
            return fPrimaryID;
        }
        
    public:
        //! mandatory to have the allocator mechanism 
        inline void *operator new(size_t);
        inline void operator delete(void *aHit);
    };
    
    extern G4ThreadLocal G4Allocator<PrimaryTrackInformation> *aPrimaryTrackInformationAllocator;
    inline void* PrimaryTrackInformation::operator new(size_t)
    {
        if( !aPrimaryTrackInformationAllocator )
            aPrimaryTrackInformationAllocator = new G4Allocator<PrimaryTrackInformation>;
        return aPrimaryTrackInformationAllocator->MallocSingle();
    }
    inline void PrimaryTrackInformation::operator delete(void *aTrackinfo)
    {
        aPrimaryTrackInformationAllocator->FreeSingle((PrimaryTrackInformation*) aTrackinfo);
    }

} // SToGS Namespace

#endif
