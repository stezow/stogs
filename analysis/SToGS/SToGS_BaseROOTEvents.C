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
//

#include "SToGS_BaseROOTEvents.h"
#include "TClass.h"

ClassImp(SBRHit)
ClassImp(SBREvent)
ClassImp(SBRPHit)
ClassImp(SBRPEvent)
ClassImp(SBROpticalEvent)

SBRHit::SBRHit() : TObject(),
    fE(0.0),
    fPDG(-1),
    fX(0.0),
    fY(0.0),
    fZ(0.0),
    fT(0.0),
    fFlag(0),
    fUID(-1)
{
	SBRHit::Class()->IgnoreTObjectStreamer();
	Clear("");
	
    //	printf("Const \n");
}

SBRHit::SBRHit(const SBRHit &sbrhit) : TObject((TObject &)sbrhit)
{
	fE		= sbrhit.fE;
	fPDG    = sbrhit.fPDG;
	fX		= sbrhit.fX;
	fY		= sbrhit.fY;
	fZ		= sbrhit.fZ;
	fT		= sbrhit.fT;
	fFlag	= sbrhit.fFlag;
	fUID	= sbrhit.fUID;
}

SBRHit::~SBRHit()
{
    //	printf("Des \n");
}

SBRPHit::SBRPHit(const SBRPHit &sbrphit) : SBRHit((SBRPHit &)sbrphit)
{
	fDX		= sbrphit.fDX;
	fDY		= sbrphit.fDY;
	fDZ		= sbrphit.fDZ;
}

SBROpticalHit::SBROpticalHit() : TObject(),
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
	SBROpticalHit::Class()->IgnoreTObjectStreamer();
	Clear("");
	
	//	printf("Const \n");
}

SBROpticalHit::SBROpticalHit(const SBROpticalHit &SBRHit) : TObject((TObject &)SBRHit)
{
	fX		= SBRHit.fX;
	fY		= SBRHit.fY;
	fZ		= SBRHit.fZ;
	fTA		= SBRHit.fTA;
	fTL		= SBRHit.fTL;
	fLength		= SBRHit.fLength;
    
	fPrimaryID	= SBRHit.fPrimaryID;
	fSecondaryID	= SBRHit.fSecondaryID;
	fNbSteps	= SBRHit.fNbSteps;
}

void SBROpticalHit::Clear(Option_t * /*opt*/)
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

SBROpticalHit::~SBROpticalHit()
{
	//	printf("Des \n");
}

SBREvent::SBREvent() : TObject(),
    fHits(0x0),
    fETot(0.0),
    fMultTot(0)
{
	fHits = new TClonesArray("SBRHit",10000);
    
	SBREvent::Class()->IgnoreTObjectStreamer();
}

SBREvent::~SBREvent()
{
	delete fHits;
}

void SBREvent::Clear(Option_t *opt)
{
	fHits->Clear(opt);
	
	fMultTot = 0;
	fETot = 0.0;
}

SBRHit *SBREvent::AddHit()
{
    //	printf("entries IN %d \n",fHits->GetEntries());
    
	TClonesArray &ar = *fHits;
	SBRHit *p = new( ar[fHits->GetEntries()] ) SBRHit();
	
    //	printf("entries OU %d \n",fHits->GetEntries());
	
	return p;
}

SBRHit *SBREvent::GetHit(Int_t which)
{
	TClonesArray &ar = *fHits;
    
	if ( which < fHits->GetEntries() )
		return (SBRHit *)ar[which];
	return 0x0;
}

SBRPEvent::SBRPEvent() : TObject(),
    fHits(0x0),
    fETot(0.0),
    fMultTot(0)
{
	fHits = new TClonesArray("SBRPHit",10000);
    
	SBRPEvent::Class()->IgnoreTObjectStreamer();
}

SBRPEvent::~SBRPEvent()
{
	delete fHits;
}

void SBRPEvent::Clear(Option_t *opt)
{
	fHits->Clear(opt);
	
	fMultTot = 0;
	fETot = 0.0;
}

SBRPHit *SBRPEvent::AddHit()
{
    //	printf("entries IN %d \n",fHits->GetEntries());
    
	TClonesArray &ar = *fHits;
	SBRPHit *p = new( ar[fHits->GetEntries()] ) SBRPHit();
	
    //	printf("entries OU %d \n",fHits->GetEntries());
	
	return p;
}

SBRPHit *SBRPEvent::GetHit(Int_t which)
{
	TClonesArray &ar = *fHits;
    
	if ( which < fHits->GetEntries() )
		return (SBRPHit *)ar[which];
	return 0x0;
}


SBROpticalEvent::SBROpticalEvent() : TObject(),
    fHits(0x0)
{
	fHits = new TClonesArray("SBROpticalHit",100000);
	//
	SBROpticalEvent::Class()->IgnoreTObjectStreamer();
}

SBROpticalEvent::~SBROpticalEvent()
{
	delete fHits;
}

void SBROpticalEvent::Clear(Option_t *opt)
{
	fHits->Clear(opt);
}

SBROpticalHit *SBROpticalEvent::AddHit()
{
	//	printf("entries IN %d \n",fHits->GetEntries());
	
	TClonesArray &ar = *fHits;
	SBROpticalHit *p = new( ar[fHits->GetEntries()] ) SBROpticalHit();
	
	//	printf("entries OU %d \n",fHits->GetEntries());
	
	return p;
}

SBROpticalHit *SBROpticalEvent::GetHit(Int_t which)
{
	TClonesArray &ar = *fHits;
	
	if ( which < fHits->GetEntries() )
		return (SBROpticalHit *)ar[which];
	return 0x0;
}

#include <algorithm>

// function to have hit sorted with increasing energy
bool CompE_ASCENDING(SBRHit *hit1, SBRHit *hit2)
{
	if ( hit2->fE > hit1->fE )
		return true;
	
	return false;
}
// function to have hit sorted with decreasing energy
bool CompE_DESCENDING(SBRHit *hit1, SBRHit *hit2)
{
	if ( hit2->fE < hit1->fE )
		return true;
	
	return false;
}

void SBREvent::CopyTo(std::vector <SBRHit *> &ordlist, Option_t *opt)
{
	TString o = opt; TClonesArray &ar = *fHits;
    
	for( Int_t i = 0; i < ar.GetEntries(); i++ ) {
		ordlist.push_back( (SBRHit *)ar[i] );
	}
	if ( o == "e>" ) {
		std::sort( ordlist.begin(), ordlist.end(), CompE_ASCENDING );
	}
	if ( o == "e<" ) {
		std::sort( ordlist.begin(), ordlist.end(), CompE_DESCENDING );
	}
}





