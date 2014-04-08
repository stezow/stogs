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

#ifndef SToGS_BaseROOTEvents_h
#define SToGS_BaseROOTEvents_h 1

#include "TObject.h"
#include "TClonesArray.h"

#include <list>

/*! SToGS Base Root Hit
 
 In case of primary particles (simulations)
 - fE : energy of the emitted particle
 - fPDG: PDG code of the emitted particle
 - fX : direction of the emitted particle
 - fY : direction of the emitted particle
 - fZ : direction of the emitted particle
 - fT : emission's time
 - fFlag : unique number in the list of emitted particles for an event
 - fUID  : user's id - unique detector number
 
 In case of a real detector
 - fE : energy of one cell [detector unit]
 - fPDG: PDG code of the reconstructed particle. Likely to be unknown or determined at the analysis level
 - fX : direction of the cell fired
 - fY : direction of the cell fired
 - fZ : direction of the cell fired
 - fT : detection's time
 - fFlag : flag concerning the cell, in principle default value
 - fUID  : user's id - unique detector number
 
 In case of simulated detector tacking mode [calo mode]
 - fE : energy of the impact [cell]
 - fPDG: PDG code of the particle giving the hit [unkown, several particle can deposited energy in a given cell]
 - fX : position of the impact [mean position in the cell]
 - fY : position of the impact [mean position in the cell]
 - fZ : position of the impact [mean position in the cell]
 - fT : time of the impact [mean time in the cell]
 - fFlag : unique number in the list of emitted particles for an event [number of hits in the cell]
 - fUID  : user's id - unique detector number
 
 In case of output of algorithms
 - fE : energy of the cluster
 - fPDG: unique id that depends of the type of the particle
 - fX : direction of the cluster
 - fY : direction of the cluster
 - fZ : direction of the cluster
 - fT : detection's time
 - fFlag : flag concerning the constructed cluster (depends on the algorithm)
 - fUID  : user's id - unique detector number
 */
class SBRHit : public TObject
{
public:
	Double32_t fE;  // energy
	Int_t fPDG;		// PDG of the particle giving this Hit
	Double32_t fX;	// X position
	Double32_t fY;	// Y position
	Double32_t fZ;	// Z position
	Double32_t fT;	// time of flight
	Int_t fFlag;	// a flag
	Int_t fUID;		// a universal ID
    
private:
	void Clear(Option_t *opt);
    
public:
	SBRHit();
	SBRHit(const SBRHit &);
	virtual ~SBRHit();
    
	ClassDef( SBRHit , 1 ) // SToGS Basic ROOR  Hit
};

/*! A SToGS Basic ROOT event is a list of SBRHits. It adds also the sum energy and event fold
 */
class SBREvent : public TObject
{
private:
	TClonesArray *fHits; //-> Collection of hits
	
	Double32_t fH;	// sum energy
	Int_t fK;		// fold
    
public:
	SBREvent();
	virtual ~SBREvent();
	
	//! helper function
	void CopyTo(std::vector <SBRHit *> &ordlist, Option_t *opt = "");
    
	// number of hits in that event
	Int_t GetNbHits() const
    { return fHits->GetEntries(); }
	
	SBRHit *IsUID(Int_t uid, Int_t *which = 0x0)
	{
		TClonesArray &ar = *fHits;
        
		SBRHit *ishit = 0x0; Int_t entries = fHits->GetEntries();
		for (Int_t j = 0; j < entries; j++) {
			SBRHit *hit = (SBRHit *)ar[j];
			if ( hit->fUID == uid ) {
				ishit = hit;
				if ( which )
					(*which) = j;
				break;
			}
		}
		return ishit;
	}
	
	//! to get a Hit
	SBRHit *GetHit(Int_t);
	//! add a hit to the current event
	SBRHit *AddHit();
	
	void AddHK(Double_t h, Int_t k = 1)
    { fH += h; fK += k; }
	
	void SetHK(Double_t h, Int_t k)
    { fH = h; fK = k; }
	
	Double_t GetH() const
    { return fH;}
	Int_t GetK() const
    { return fK; }
	
	//! clear the collection of hits, set H, K to 0
	void Clear(Option_t *opt);
	
	ClassDef( SBREvent , 1 ) // SToGS Basic ROOR Event
};

/*! A SToGS Basic ROOT event is a list of SBRHits. It adds also the sum energy and event fold
 */
class SBROpticalHit : public TObject
{
public:
	Double32_t fX;	// X position of the last point of the track
	Double32_t fY;	// Y position of the last point of the track
	Double32_t fZ;	// Z position of the last point of the track
	Double32_t fTA;	// arrival time for the last point of the track since event
	Double32_t fTL;	// arrival time for the last point of the track	since creation of the photon
	
	Double32_t fLength; // length of the trace
	
	Int_t fPrimaryID;	// Primary ID, in case more that one primary particle is giving scintillation
	Int_t fSecondaryID;	// detector ID where the optical photon is
	Int_t fNbSteps;		// number of steps between start and stop
	
private:
	void Clear(Option_t *opt);
	
public:
	SBROpticalHit();
	SBROpticalHit(const SBROpticalHit &);
	virtual ~SBROpticalHit();
	
	ClassDef( SBROpticalHit , 1 ) // Paris Optical Hit
};

/*! A SToGS Basic ROOT event is a list of SBRHits. It adds also the sum energy and event fold
 */
class SBROpticalEvent : public TObject
{
private:
	TClonesArray *fHits; //-> Collection of optical hits
	
public:
	SBROpticalEvent();
	virtual ~SBROpticalEvent();
    
	// number of hits in that event
	Int_t GetNbHits() const
    { return fHits->GetEntries(); }
	
	//! to get a Hit
	SBROpticalHit *GetHit(Int_t);
	//! add a hit to the current event
	SBROpticalHit *AddHit();
	
	void Clear(Option_t *opt);
    
	ClassDef( SBROpticalEvent , 1 ) // Paris Optical Event
};

#endif
