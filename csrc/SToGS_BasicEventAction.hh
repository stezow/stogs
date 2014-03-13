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

#ifndef ParisBasicEventAction_h
#define ParisBasicEventAction_h 1

#include "G4UserEventAction.hh"
#include "ParisOutputManager.hh"

//! SToGS namespace to protect SToGS classes
namespace SToGS {
//! Basic event action consists just in calling ParisOutputManager event actions
/*!
	In principle no need to modify. The real work is done through the ouputmanager which is configured @ the beginning
*/
class ParisBasicEventAction : public G4UserEventAction
{
private:
	ParisOutputManager *fOutManager;
	
public:
	ParisBasicEventAction(ParisOutputManager *outmanager = 0x0);
	virtual ~ParisBasicEventAction() 
		{;}
	
public:
	virtual void BeginOfEventAction(const G4Event *event) 
	{ 
		fOutManager->BeginOfEventAction(event) ; 
	}
	virtual void EndOfEventAction(const G4Event *event)   
	{ 
		fOutManager->EndOfEventAction(event) ; 
	}
};
} // SToGS Namespace

#endif    
