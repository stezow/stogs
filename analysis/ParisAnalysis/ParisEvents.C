

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


///// ======== Michal ========== ////////////////////////////////////////////


ClassImp(cBasicREv)

cBasicREv::cBasicREv()
{
	fMult=0;
	fDetId=0x0;
	fEnergy=0x0;
}


cBasicREv:: cBasicREv(Int_t Mult)
{
	fMult=Mult;
	fDetId= new Int_t[fMult];
	fEnergy= new Double_t[fMult];
	
}

cBasicREv::  ~cBasicREv()
{
	
	if (fDetId!=0x0) delete [] fDetId;
	if (fEnergy!=0x0) delete [] fEnergy;
}


Int_t cBasicREv::SetMultiplicity(Int_t Mult)
{
	fMult=Mult;
	if (fDetId!=0x0) delete [] fDetId;
	if (fEnergy!=0x0) delete [] fEnergy;
	fDetId= new Int_t[fMult];
	fEnergy= new Double_t[fMult];
#ifdef DEBUG
	std::cout<<"\nIn cBasicREv::SetMultiplicity() \n fMult= "<<fMult<<std::flush;
#endif
	
	return 0;
}

Int_t cBasicREv::Fill( Int_t Mult, Int_t *DetId, Double_t *Energy)
{
	fMult=Mult;
	if (fDetId!=0x0) delete [] fDetId;
	if (fEnergy!=0x0) delete [] fEnergy;
	
	fDetId= new Int_t[fMult];
	fEnergy= new Double_t[fMult];
	
	memcpy( (void *)fDetId, (void *)DetId, fMult*sizeof(Int_t));
	memcpy( (void *)fEnergy, (void *)Energy, fMult*sizeof(Double_t));
	return 0;
}

Int_t cBasicREv::Fill2( Int_t NrInteraction, Int_t DetId, Double_t Energy) // NrInteraction starts from 1
{
	if( NrInteraction >=0 && NrInteraction <= fMult ){
		fDetId[NrInteraction] = DetId;
		fEnergy[NrInteraction] = Energy;
		return 0;
	}	 
	return -1;
}

Double_t cBasicREv::SummEnergy()
{
	Double_t Esum=0.;
	
	for( Int_t i=0;i<fMult;i++)
		Esum += fEnergy[i];
	return Esum;
}


Int_t cBasicREv::Copy(cBasicREv BasicREv)
{
	Fill(BasicREv.fMult,BasicREv.fDetId,BasicREv.fEnergy);
	return 0;
}

Double_t cBasicREv::GetEnergy( Int_t DetId)
{
	Double_t Esum=0.;
	for( Int_t i=0;i<fMult;i++)
		if( DetId== fDetId[i])
			Esum+= fEnergy[i];
	
	return Esum;
}


ClassImp(cRootEv)

cRootEv::cRootEv():fNrPGamma(0),fEnPGamma(0x0),fDetEv(0x0)
{
	
}

cRootEv::cRootEv( Int_t NrPGamma, Double_t *EnPGamma, cBasicREv *DetEv)
{
	fNrPGamma=NrPGamma;
	fEnPGamma= new Double_t[fNrPGamma];
	fDetEv = new cBasicREv[fNrPGamma];
	
	memcpy( (void *)fEnPGamma, (void *)EnPGamma, fNrPGamma*sizeof(Double_t));
	for( Int_t i=0;i<fNrPGamma;i++){
		fDetEv[i].Copy(DetEv[i]);
	}  
}

cRootEv:: ~cRootEv()
{
	if( fEnPGamma!=0x0) delete [] fEnPGamma;
	if( fDetEv!=0x0) delete [] fDetEv;
	
}

//  Int_t  InitRootTree( char *treeName="Paris_Tree");

Int_t cRootEv::Fill( Int_t NrPGamma, Double_t *EnPGamma, cBasicREv *DetEv)
{
	if( NrPGamma<1) return -1;
	
#ifdef DEBUG
	std::cout<<"\nIn cROOTEv::Fill() \n NRPGamma= "<<NrPGamma<<std::flush;
#endif
	
	fNrPGamma=NrPGamma;
	
	
	if( fEnPGamma!=0x0) delete [] fEnPGamma;
	if( fDetEv!=0x0) delete [] fDetEv;
	
	fEnPGamma= new Double_t[fNrPGamma];
	fDetEv = new cBasicREv[fNrPGamma];
	
#ifdef DEBUG
	std::cout<<"\n Definition of fEnPGamma... "<<std::flush;
#endif
	if( EnPGamma!=0x0)
		memcpy( (void *)fEnPGamma, (void *)EnPGamma, fNrPGamma*sizeof(Double_t));
	
#ifdef DEBUG
	std::cout<<"OK!\nDefinition of fDetEv..."<<std::flush;
#endif
	
	if( DetEv!=0x0)
		for( Int_t i=0;i<fNrPGamma;i++){
			fDetEv[i].Copy(DetEv[i]);
		}  
#ifdef DEBUG
	std::cout<<"\nExiting cROOTEv::Fill()... "<<std::flush;
#endif
	
	return 0;
}

Int_t cRootEv::DefineEventDep(  Int_t NrGamma, Int_t Mult, Int_t *DetId, Double_t *Energy)
{
	return   fDetEv[NrGamma].Fill(Mult,DetId,Energy);
	
}

Int_t cRootEv::DefineEventDep2(  Int_t NrGamma, Int_t NrInteraction, Int_t DetId, Double_t Energy)
{
	return fDetEv[NrGamma].Fill2(NrInteraction, DetId, Energy); 
}

Double_t cRootEv::SummEnergy(Int_t NrGamma)
{
	return fDetEv[NrGamma].SummEnergy();
	
}

Double_t cRootEv::SummEnergyDet( Int_t DetId)
{
	Double_t Esum=0.;
	
	for(Int_t i=0;i< fNrPGamma; i++){
		Esum+= fDetEv[i].GetEnergy( DetId);
		
	}
	return Esum;
}

Double_t cRootEv::SummEnergyInDet(Int_t NrGamma, Int_t DetId)
{
	Double_t Esum =0;
	Esum = fDetEv[NrGamma].GetEnergy(DetId);
	
	return Esum;
}
Int_t cRootEv::SetMultiplicity( Int_t NrGamma, Int_t Mult)
{
	fDetEv[NrGamma].SetMultiplicity(Mult);
	
	return 0;
}

Int_t cRootEv::SetPrEnergy( Int_t NrGamma, Double_t PrEnergy)
{
	if( NrGamma>=0 && NrGamma<fNrPGamma) {
		fEnPGamma[NrGamma]= PrEnergy ;
    return 0;
  }
  return -1;
}

Double_t cRootEv::GetPrEnergy( Int_t NrGamma)
{
  if(NrGamma>=0 && NrGamma<fNrPGamma){
    return fEnPGamma[NrGamma];	
}
  return -1;
}

///// ======== Michal ========== ////////////////////////////////////////////





