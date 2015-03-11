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
 
#ifndef toROOTGPSPrimaryGeneratorAction_h
#define toROOTGPSPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"
//
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"
#include "G4SPSRandomGenerator.hh"

#include "G4VUserPrimaryGeneratorAction.hh"

// to know SToGS events
#include "SToGS_BaseROOTEvents.h"


class G4Event;
class toROOTGPSPrimaryGeneratorActionMessanger;
class TChain;


//! first step to a general file (ROOT) GPS generator in which events are read from files and/or generatated through distributions
/*!
\code

\endcode	
*/
class toROOTGPSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	toROOTGPSPrimaryGeneratorAction();  
	toROOTGPSPrimaryGeneratorAction(G4String);      
	virtual ~toROOTGPSPrimaryGeneratorAction();

public:
	virtual void GeneratePrimaries(G4Event* anEvent);
  
public:     
	//! from a given file, it reads the characteristics of the gamma cascade 
	/*!
	*/	     
	void GetInputFiles(G4String filename = "setup/toROOTGPS");
    
private:
	/*
    G4SPSPosDistribution* posGenerator;
	G4SPSAngDistribution* angGenerator;
	G4SPSEneDistribution* eneGenerator;
	G4SPSRandomGenerator* biasRndm;
	//
	// Other particle properties
	G4int NumberOfParticlesToBeGenerated;
	G4ParticleDefinition * particle_definition;
	G4ParticleMomentum particle_momentum_direction;
	G4double particle_energy;
	G4double particle_charge;
	G4ThreeVector particle_position;
	G4double particle_time;
	G4ThreeVector particle_polarization;
	G4double particle_weight;
     */
    //! TChain containing ROOT events
    TChain *theChainOfPrimaryEvents;
    //! Current primary from TTree
    SBRPEvent *fPr;
    //! Current entry to be read
    Long64_t fCurrentEntry;
    Long64_t fEntries;

    //!
	toROOTGPSPrimaryGeneratorActionMessanger *theMessanger;
};

#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;

//! Messanger class for ParisBasicPrimaryGenerator
/*!
	The messanger class allows to reset the generator inside a Geant4 session with for instance:
	
	 \c /Paris/generator/BasicPrimaryGenerator \c setup/mynewbasic.gene
	
*/
class toROOTGPSPrimaryGeneratorActionMessanger: public G4UImessenger
{
public:
	toROOTGPSPrimaryGeneratorActionMessanger(toROOTGPSPrimaryGeneratorAction *);
	~toROOTGPSPrimaryGeneratorActionMessanger();
		
	void SetNewValue(G4UIcommand*, G4String);
	
private:
	toROOTGPSPrimaryGeneratorAction *theGenerator;
	
	G4UIdirectory *theDirectory;
	G4UIcmdWithAString *resetParametersCmd;
};
#endif


