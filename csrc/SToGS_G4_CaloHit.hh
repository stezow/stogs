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

#ifndef SToGS_CaloHit_h
#define SToGS_CaloHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Types.hh"

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! Informations to keep the energy is in a sensitive detector
    /*!
     In G4, one have to define a class (a hit) that inherits from G4VHit to keep informations during the tracking.
     This one is dedicated for tracking since it kept infomations of any single impact in sensitive detectors. Here is
     the list of the available informations for each ParisCaloHit. Use Get/Set methods to obtain/change these values.
     
     - G4int primaryID: Primary vertex from which this hit is coming
     
     - G4int detID: Detector ID in which the hit occured
     - G4int motherID: Mother ID of the detector in which the hit occured
     
     - G4double edep: Energy deposited
     - G4double ToF:  global time at which the interaction occured
     
     - G4ThreeVector pos: Position of the impact
     */
    class CaloHit : public G4VHit
    {
    public:
        CaloHit();
        virtual ~CaloHit();
        CaloHit(const CaloHit &right);
        
        const CaloHit& operator=(const CaloHit &right);
        G4int operator==(const CaloHit &right) const;
        
        inline void *operator new(size_t);
        inline void operator delete(void *aHit);
        
    public:
// setters
        void SetDetID  (G4int id)
        {
            detID = id;
        }
        void SetMotherID(G4int id)
        {
            motherID = id;
        }
        void SetEdep(G4double de = 0.0)
        {
            edep = de;
        }
        void SetToF (G4double tf = 0.0)
        {
            ToF = tf;
        }
        void SetPos(G4ThreeVector xyz = G4ThreeVector(0.,0.,0.))
        {
            pos = xyz;
        }
        void SetDetName(G4String name)
        {
            detName = name;
        }
        void SetMotherDetName(G4String name)
        {
            motherDetName = name;
        }
        void SetNbHits(G4int nb = 0)
		{
            NbHits = nb;
        }
// getters
        G4int GetDetID() const
		{
            return detID;
        }
        G4int GetMotherID() const
		{
            return motherID;
        }
        G4double GetEdep() const
		{
            return edep;
        }
        G4double GetToF()  const
		{
            return ToF;
        }
        const G4ThreeVector &GetPos()
        {
            return pos;
        }
        G4int GetNbHits() const
		{
            return NbHits;
        }
// others
        void AddOneHit()
		{
            NbHits++;
        }
        void AddEdep(G4double de)
		{
            edep += de;
        }
        void AddToF (G4double de, G4double t)
		{
            ToF += de*t;
        }
        void AddPos (G4double de, G4ThreeVector xyz)
		{
            pos.setX(pos.x()+de*xyz.x()); pos.setY(pos.y()+de*xyz.y()); pos.setZ(pos.z()+de*xyz.z());
        }
        //! to be applied at the end of an event to calculate properly the mean values
        void EndOfEvent();
        
        void Draw();
        void Print();
        
    private:
        G4int detID;
        G4int motherID;
        
        G4double edep;
        G4double ToF;
        
        G4ThreeVector pos;
        
        G4String detName;
        G4String motherDetName;
        
        G4int NbHits;
    };

    typedef G4THitsCollection<CaloHit> CaloHitsCollection;
    extern G4ThreadLocal G4Allocator<CaloHit> *CaloHitAllocator;
    
    inline void* CaloHit::operator new(size_t)
    {
        if( !CaloHitAllocator )
            CaloHitAllocator = new G4Allocator<CaloHit>;
        return CaloHitAllocator->MallocSingle();
    }
    inline void CaloHit::operator delete(void *aHit)
    {
        CaloHitAllocator->FreeSingle((CaloHit*) aHit);
    }
} // SToGS Namespace

#endif
