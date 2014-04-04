

#ifndef _ParisEvents
#define _ParisEvents 1

#include "TObject.h"
#include "TClonesArray.h"

#include <list>

/*! General Hit structure for Paris
 
 In case of primary particles (simulations)
	- fE : energy of the emitted particle
	- fID: unique id that depends of the type of the particle
	- fX : direction of the emitted particle
	- fY : direction of the emitted particle
	- fZ : direction of the emitted particle
	- fT : emission's time
	- fFlag : unique number in the list of emitted particles for an event
	- fUID  : user's id - unique detector number
 
 In case of a real detector
	- fE : energy of one cell
	- fID: type is in principle not known 
	- fX : direction of the cell fired
	- fY : direction of the cell fired
	- fZ : direction of the cell fired
	- fT : detection's time
	- fFlag : flag concerning the cell, in principle default value
	- fUID  : user's id - unique detector number

 In case of simulated detector tacking mode [calo mode]
	- fE : energy of the impact [cell]
	- fID: type of the primary particle that gives this impact [unkown, two primarires can deposited energy in a given cell]
	- fX : position of the impact [mean position in the cell]
	- fY : position of the impact [mean position in the cell]
	- fZ : position of the impact [mean position in the cell]
	- fT : time of the impact [mean time in the cell]
	- fFlag : unique number in the list of emitted particles for an event [number of hits in the cell]
	- fUID  : user's id - unique detector number
 
 In case of output of algorithms
	- fE : energy of the cluster
	- fID: unique id that depends of the type of the particle
	- fX : direction of the cluster
	- fY : direction of the cluster
	- fZ : direction of the cluster
	- fT : detection's time
	- fFlag : flag concerning the constructed cluster (depends on the algorithm)
	- fUID  : user's id - unique detector number
*/
class PHit : public TObject
{	
public:
	enum EParticle { kUnknown = -1, kGamma, kGDR, kNeutron };

public:
	Double32_t fE; // energy
	Int_t fID;		// an ID
	Double32_t fX;	// X position
	Double32_t fY;	// Y position
	Double32_t fZ;	// Z position
	Double32_t fT;	// time of flight
	Int_t fFlag;	// a flag
	Int_t fUID;		// a universal ID

private:
	void Clear(Option_t *opt);

public:
	PHit();
	PHit(const PHit &);
	virtual ~PHit();	
		
	ClassDef( PHit , 1 ) // Paris Hit
};

/*! A Paris Event
	A Paris event is composed of a collection of PHits with H, K associated
*/
class PEvent : public TObject
{
private:
	TClonesArray *fHits; //-> Collection of hits
	
	Double32_t fH;	// sum energy
	Int_t fK;		// fold

public:
	PEvent();
	virtual ~PEvent();	
	
	//! helper function 
	void CopyTo(std::vector <PHit *> &ordlist, Option_t *opt = "");

	// number of hits in that event
	Int_t GetNbHits() const 
		{ return fHits->GetEntries(); }
	
	PHit *IsUID(Int_t uid, Int_t *which = 0x0)
	{
		TClonesArray &ar = *fHits;

		PHit *ishit = 0x0; Int_t entries = fHits->GetEntries();
		for (Int_t j = 0; j < entries; j++) {
			PHit *hit = (PHit *)ar[j];
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
	PHit *GetHit(Int_t);
	//! add a hit to the current event
	PHit *AddHit();
	
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
	
	ClassDef( PEvent , 1 ) // Paris Event
};

class POpticalHit : public TObject
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
	POpticalHit();
	POpticalHit(const POpticalHit &);
	virtual ~POpticalHit();	
	
	ClassDef( POpticalHit , 1 ) // Paris Optical Hit
};

class POpticalEvent : public TObject
{
private:
	TClonesArray *fHits; //-> Collection of optical hits
	
public:
	POpticalEvent();
	virtual ~POpticalEvent();	

	// number of hits in that event
	Int_t GetNbHits() const 
		{ return fHits->GetEntries(); }
	
	//! to get a Hit
	POpticalHit *GetHit(Int_t);
	//! add a hit to the current event
	POpticalHit *AddHit();
	
	void Clear(Option_t *opt);

	ClassDef( POpticalEvent , 1 ) // Paris Optical Event
	
};

///// ======== Michal ========== ////////////////////////////////////////////

/*! 
*/
class cBasicREv {
public:
	
	Int_t fMult; // primary event multiplicity
	Int_t *fDetId; // [fMult] detector id's, one for each energy event
	Double_t *fEnergy; // [fMult] inetraction enery in each tuched detector
	
	cBasicREv();
	cBasicREv(Int_t Mult);
	~cBasicREv();
	
	Int_t SetMultiplicity(Int_t Mult);
	Int_t Fill( Int_t Mult, Int_t *DetId, Double_t *Energy);
	Int_t Fill2( Int_t NrInteraction, Int_t DetId, Double_t Energy);  
	Int_t Copy(cBasicREv BasicREv);
	Double_t SummEnergy();
	Double_t GetEnergy( Int_t DetId);
	ClassDef(cBasicREv,1)
};


class cRootEv {
public:
	Int_t fNrPGamma; // Nr of primary gammas in the generated event (cascade gamma)
	Double_t *fEnPGamma; //[fNrPGamma] Energies for the generated gammas
	cBasicREv *fDetEv; // [fNrPGamma] Events in the detectors gerated by the primary gammas
	cRootEv();
	cRootEv( Int_t NrPGamma, Double_t *EnPGamma=0x0, cBasicREv *DetEv=0x0);
	~cRootEv();
	
	//  Int_t  InitRootTree( char *treeName="Paris_Tree");
	Int_t SetMultiplicity( Int_t NrGamma, Int_t Mult);
	Int_t Fill( Int_t NrPGamma, Double_t *EnPGamma=0x0, cBasicREv *DetEv=0x0);
	Int_t DefineEventDep(  Int_t NrGamma, Int_t Mult, Int_t *DetId, Double_t *Energy);	  
	Int_t DefineEventDep2(  Int_t NrGamma, Int_t NrInteraction, Int_t DetId, Double_t Energy);  
	Int_t SetPrEnergy( Int_t NrGamma, Double_t PrEnergy);
	Double_t GetPrEnergy( Int_t NrGamma); // return energy of primary gamma
	Double_t SummEnergy(Int_t NrGamma); // return the total energy deposted in the detectors for the gamma nr NrPGamma
	Double_t SummEnergyInDet(Int_t NrGamma, Int_t DetId); //return  the total enrgy deposited in the detector DetId for the gamma Nr NrGamma
	Double_t SummEnergyDet( Int_t DetId); // return the total energy deposted in the detector nr DetId
  ClassDef(cRootEv,1)
};


///// ======== Michal ========== ////////////////////////////////////////////

#endif

