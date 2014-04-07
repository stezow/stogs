

#include "ParisEvents.h"
#include "TClass.h"

ClassImp(PHit)
ClassImp(PEvent)
ClassImp(POpticalEvent)

PHit::PHit() : TObject(),
	fE(0.0),
	fID(kUnknown),
	fX(0.0),
	fY(0.0),
	fZ(0.0),
	fT(0.0),
	fFlag(0),
	fUID(-1)
{
	PHit::Class()->IgnoreTObjectStreamer(); 
	Clear("");
	
//	printf("Const \n");
}

PHit::PHit(const PHit &phit) : TObject((TObject &)phit)
{
	fE		= phit.fE;
	fID	= phit.fID;
	fX		= phit.fX;
	fY		= phit.fY;	
	fZ		= phit.fZ;
	fT		= phit.fT;	
	fFlag	= phit.fFlag;
	fUID	= phit.fUID;
}


void PHit::Clear(Option_t * /*opt*/)
{
	fE  = 0;
	fID = kUnknown;
	fX  = 0.0;
	fY  = 0.0;
	fZ  = 0.0;
	fT  = 0.0;
	fFlag = 0;
	fUID  = -1;		
}

PHit::~PHit()
{
//	printf("Des \n");
}

POpticalHit::POpticalHit() : TObject(),
	fX(0.0),
	fY(0.0),
	fZ(0.0),
	fTA(0.0),
	fTL(0.0),
	fLength(0.0),
	fPrimaryID(0.0),
	fSecondaryID(0.0),
	fNbSteps(0.0)
{
	POpticalHit::Class()->IgnoreTObjectStreamer(); 
	Clear("");
	
	//	printf("Const \n");
}

POpticalHit::POpticalHit(const POpticalHit &phit) : TObject((TObject &)phit)
{	
	fX		= phit.fX;
	fY		= phit.fY;	
	fZ		= phit.fZ;
	fTA		= phit.fTA;	
	fTL		= phit.fTL;	
	fLength		= phit.fLength;

	fPrimaryID	= phit.fPrimaryID;
	fSecondaryID	= phit.fSecondaryID;
	fNbSteps	= phit.fNbSteps;
}


void POpticalHit::Clear(Option_t * /*opt*/)
{
	fX  = 0.0;
	fY  = 0.0;
	fZ  = 0.0;
	fTA  = 0.0;
	fTL  = 0.0;
	fLength	= 0.0;

	fPrimaryID  = -1;
	fSecondaryID = -1;
	fNbSteps = 0;
}

POpticalHit::~POpticalHit()
{
	//	printf("Des \n");
}

PEvent::PEvent() : TObject(),
	fHits(0x0),
	fH(0.0),
	fK(0)
{
	fHits = new TClonesArray("PHit",10000); 

	PEvent::Class()->IgnoreTObjectStreamer(); 
}

PEvent::~PEvent()
{
	delete fHits;
}

void PEvent::Clear(Option_t *opt)
{
	fHits->Clear(opt);
	
	fK = 0; 
	fH = 0.0;
}

PHit *PEvent::AddHit()
{
//	printf("entries IN %d \n",fHits->GetEntries());

	TClonesArray &ar = *fHits;
	PHit *p = new( ar[fHits->GetEntries()] ) PHit();
	
//	printf("entries OU %d \n",fHits->GetEntries());
	
	return p;
}

PHit *PEvent::GetHit(Int_t which)
{
	TClonesArray &ar = *fHits;

	if ( which < fHits->GetEntries() )
		return (PHit *)ar[which];
	return 0x0;
}

POpticalEvent::POpticalEvent() : TObject(),
	fHits(0x0)
{
	fHits = new TClonesArray("POpticalHit",100000); 
	
	POpticalEvent::Class()->IgnoreTObjectStreamer(); 
}

POpticalEvent::~POpticalEvent()
{
	delete fHits;
}

void POpticalEvent::Clear(Option_t *opt)
{
	fHits->Clear(opt);
}

POpticalHit *POpticalEvent::AddHit()
{
	//	printf("entries IN %d \n",fHits->GetEntries());
	
	TClonesArray &ar = *fHits;
	POpticalHit *p = new( ar[fHits->GetEntries()] ) POpticalHit();
	
	//	printf("entries OU %d \n",fHits->GetEntries());
	
	return p;
}

POpticalHit *POpticalEvent::GetHit(Int_t which)
{
	TClonesArray &ar = *fHits;
	
	if ( which < fHits->GetEntries() )
		return (POpticalHit *)ar[which];
	return 0x0;
}

#include <algorithm>

// function to have hit sorted with increasing energy
bool CompE_ASCENDING(PHit *hit1, PHit *hit2)
{
	if ( hit2->fE > hit1->fE )
		return true;
	
	return false;
}
// function to have hit sorted with decreasing energy
bool CompE_DESCENDING(PHit *hit1, PHit *hit2)
{
	if ( hit2->fE < hit1->fE )
		return true;
	
	return false;
}

void PEvent::CopyTo(std::vector <PHit *> &ordlist, Option_t *opt)
{
	TString o = opt; TClonesArray &ar = *fHits;

	for( Int_t i = 0; i < ar.GetEntries(); i++ ) {
		ordlist.push_back( (PHit *)ar[i] );
	}		
	if ( o == "e>" ) {
		std::sort( ordlist.begin(), ordlist.end(), CompE_ASCENDING );
	}
	if ( o == "e<" ) {
		std::sort( ordlist.begin(), ordlist.end(), CompE_DESCENDING );
	}
}






