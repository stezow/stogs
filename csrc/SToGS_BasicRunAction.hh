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


#ifndef ParisBasicRunAction_h
#define ParisBasicRunAction_h 1

#include "G4UserRunAction.hh"
#include "ParisOutputManager.hh"

class G4Run;

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! Basic event action consists just in calling ParisOutputManager run actions
/*!
 In principle no need to modify. The real work is done through the ouputmanager which is configured @ the beginning
*/
class ParisBasicRunAction : public G4UserRunAction
{
private:
	ParisOutputManager *fOutManager;
	
public:
	ParisBasicRunAction(ParisOutputManager *outmanager = 0x0);
	virtual ~ParisBasicRunAction() 
		{;}
	
public:
	virtual void BeginOfRunAction(const G4Run *therun) 
	{ 
		fOutManager->BeginOfRunAction(therun); 
	}
	virtual void EndOfRunAction(const G4Run *therun)   
	{ 
		fOutManager->EndOfRunAction(therun); 
	}
};

} // SToGS Namespace

#endif
