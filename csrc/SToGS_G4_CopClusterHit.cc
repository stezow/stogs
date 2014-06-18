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

#include "SToGS_G4_CopClusterHit.hh"

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"

#ifdef G4MULTITHREADED
G4ThreadLocal G4Allocator<SToGS::CopClusterHit> *SToGS::CopClusterHitAllocator = 0;
#else
G4Allocator<SToGS::CopClusterHit> *SToGS::CopClusterHitAllocator = 0;
#endif

SToGS::CopClusterHit::CopClusterHit():
    detID(-1),
    motherID(-1),
    edep(0.0),
    ToF(0.0),
    pos(0.0,0.0,0.0),
    detName("-"),
    motherDetName("-")
{
}

SToGS::CopClusterHit::~CopClusterHit()
{
    ;
}

SToGS::CopClusterHit::CopClusterHit(const SToGS::CopClusterHit &from) : G4VHit()
{
	if ( this == &from ) return;

	edep = from.edep;	
	pos  = from.pos;

	detID = from.detID;
	motherID = from.motherID;

	ToF = from.ToF;
		
	detName = from.detName;
	motherDetName = from.motherDetName;
    
    NbHits = from.NbHits;
}

const SToGS::CopClusterHit& SToGS::CopClusterHit::operator=(const SToGS::CopClusterHit &from)
{
	if ( this == &from )
        return (*this);
	
	edep = from.edep;
	pos  = from.pos;

	detID = from.detID;
	motherID = from.motherID;

	ToF = from.ToF;
	
	detName = from.detName;
	motherDetName = from.motherDetName;
		
    NbHits = from.NbHits;

	return (*this);
}

G4int SToGS::CopClusterHit::operator==(const SToGS::CopClusterHit &from) const
{
	return (this==&from) ? 1 : 0;
}

void SToGS::CopClusterHit::Draw()
{
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
	if(pVVisManager)
	{
		G4Circle circle(pos);
		circle.SetScreenSize(0.04); circle.SetFillStyle(G4Circle::filled);
		
		G4Colour colour(1.,0.,0.);
		G4VisAttributes attribs(G4Colour(1.,0.,0.)); circle.SetVisAttributes(attribs);
		pVVisManager->Draw(circle);
	}
}

void SToGS::CopClusterHit::Print()
{
	G4cout << "--------------- SToGS::CopClusterHit::Print() --------------------- " << G4endl;
	G4cout  << "detID: " << detID 
	 	<< ", detName: " << detName 
		<< ", motherID: " << motherID 
	 	<< ", motherDetName: " << motherDetName << G4endl; 
	G4cout << "edep: " << edep / CLHEP::keV << " keV"
	 	<< ", pos: " << pos 
	 	<< ", ToF: " << ToF  / CLHEP::ns << " ns" << G4endl; 
	G4cout << "--------------------------------------------------------------" << G4endl;
}

void SToGS::CopClusterHit::EndOfEvent()
{
	if ( edep == 0.0 )
        return;
    
    ToF /= edep;
    pos /= edep;
}








