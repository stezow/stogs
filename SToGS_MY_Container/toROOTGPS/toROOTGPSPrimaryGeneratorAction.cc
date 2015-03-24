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

#include "toROOTGPSPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "TChain.h"

// std includes
#include <iostream>
#include <sstream>
#include <fstream>

#if G4MULTITHREADED
#include "G4Threading.hh"
#endif
#ifdef G4MULTITHREADED
#include "G4AutoLock.hh"
namespace { G4Mutex buildMutex = G4MUTEX_INITIALIZER; G4Mutex runMutex = G4MUTEX_INITIALIZER; }
#endif

using namespace std;


//!
toROOTGPSPrimaryGeneratorActionMessanger *theMessanger;

toROOTGPSPrimaryGeneratorAction::toROOTGPSPrimaryGeneratorAction() :
    G4VUserPrimaryGeneratorAction(),
    theChainOfPrimaryEvents(0x0),
    fPr(0x0),
    fCurrentEntry(0L),
    fEntries(0L)
{
    fPr = new SBRPEvent(); fEntries = fCurrentEntry = 0L;
    /*
    
    NumberOfParticlesToBeGenerated = 1;
	particle_definition = G4Geantino::GeantinoDefinition();
	G4ThreeVector zero;
	particle_momentum_direction = G4ParticleMomentum(1, 0, 0);
	particle_energy = 1.0 * MeV;
	particle_position = zero;
	particle_time = 0.0;
	particle_polarization = zero;
	particle_charge = 0.0;
	particle_weight = 1.0;
    
	biasRndm = new G4SPSRandomGenerator();
	posGenerator = new G4SPSPosDistribution();
	posGenerator->SetBiasRndm(biasRndm);
	angGenerator = new G4SPSAngDistribution();
	angGenerator->SetPosDistribution(posGenerator);
	angGenerator->SetBiasRndm(biasRndm);
	eneGenerator = new G4SPSEneDistribution();
	eneGenerator->SetBiasRndm(biasRndm);
     
     */
    GetInputFiles();
    theMessanger = new toROOTGPSPrimaryGeneratorActionMessanger(this);
}

toROOTGPSPrimaryGeneratorAction::toROOTGPSPrimaryGeneratorAction(G4String filename) :
    G4VUserPrimaryGeneratorAction(),
    theChainOfPrimaryEvents(0x0),
    fPr(0x0),
    fCurrentEntry(0L),
    fEntries(0L)
{
    fPr = new SBRPEvent(); fEntries = fCurrentEntry = 0L;

    GetInputFiles(filename);
    theMessanger = new toROOTGPSPrimaryGeneratorActionMessanger(this);
}

toROOTGPSPrimaryGeneratorAction::~toROOTGPSPrimaryGeneratorAction()
{
    delete theMessanger;
}

void toROOTGPSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    G4ThreeVector particle_position(0,0,0), particle_momentum_direction;
    
    G4PrimaryParticle* particle = 0x0; G4PrimaryVertex* vertex = new G4PrimaryVertex(particle_position, 0);
    
    {
#if G4MULTITHREADED
        G4AutoLock lock(&runMutex); // better to be protected since root is doing some stuff in global
#endif
        if ( !(fCurrentEntry < fEntries) )
            return;
        theChainOfPrimaryEvents->GetEntry(fCurrentEntry++);
    }
    
    
    // create a new vertex
    for (G4int i = 0; i < fPr->GetNbHits(); i++) {
        
        // get next Pr Particle
        SBRPHit *aprimary = fPr->GetHit(i);
        particle =
            new G4PrimaryParticle(aprimary->fPDG);
        
        particle->SetKineticEnergy( aprimary->fE*CLHEP::keV );
        particle_momentum_direction.set(aprimary->fPX,aprimary->fPY,aprimary->fPZ);
        particle->SetMomentumDirection( particle_momentum_direction.unit() );

        //particle->SetPolarization(0,0,0);
        /*if (verbosityLevel > 1) {

            G4cout << "Particle name: " << particle_definition->GetParticleName() << G4endl;
            G4cout << "       Energy: " << particle_energy << G4endl;
            G4cout << "     Position: " << particle_position << G4endl;
            G4cout << "    Direction: " << particle_momentum_direction
            << G4endl;
        } */
        // Set bweight equal to the multiple of all non-zero weights
        // particle_weight = eneGenerator->GetWeight()*biasRndm->GetBiasWeight();
        // pass it to primary particle
        particle->SetWeight(1);
        
        vertex->SetPrimary(particle);
    }
    anEvent->AddPrimaryVertex(vertex);
    
    /*
    G4AutoLock l(&mutex);
    if (particle_definition == NULL)
        return;
    
    if (verbosityLevel > 1)
        G4cout << " NumberOfParticlesToBeGenerated: "
        << NumberOfParticlesToBeGenerated << G4endl;
    
    // Position stuff
    particle_position = posGenerator->GenerateOne();
    
    // create a new vertex
    G4PrimaryVertex* vertex = new G4PrimaryVertex(particle_position,
                                                  131       particle_time);
    
    for (G4int i = 0; i < NumberOfParticlesToBeGenerated; i++) {
        // Angular stuff
        particle_momentum_direction = angGenerator->GenerateOne();
        // Energy stuff
        particle_energy = eneGenerator->GenerateOne(particle_definition);
        
        if (verbosityLevel >= 2)
            G4cout << "Creating primaries and assigning to vertex" << G4endl;
        // create new primaries and set them to the vertex
        G4double mass = particle_definition->GetPDGMass();
        G4PrimaryParticle* particle =
        new G4PrimaryParticle(particle_definition);
        particle->SetKineticEnergy( particle_energy );
        particle->SetMass( mass );
        particle->SetMomentumDirection( particle_momentum_direction );
        particle->SetCharge( particle_charge );
        particle->SetPolarization(particle_polarization.x(), particle_polarization.y(),particle_polarization.z());
        if (verbosityLevel > 1) {
            G4cout << "Particle name: "
            << particle_definition->GetParticleName() << G4endl;
            G4cout << "       Energy: " << particle_energy << G4endl;
            G4cout << "     Position: " << particle_position << G4endl;
            G4cout << "    Direction: " << particle_momentum_direction
            << G4endl;
        }
        // Set bweight equal to the multiple of all non-zero weights
        particle_weight = eneGenerator->GetWeight()*biasRndm->GetBiasWeight();
        // pass it to primary particle
        particle->SetWeight(particle_weight);
        
        vertex->SetPrimary(particle);
        
    }
    // now pass the weight to the primary vertex. CANNOT be used here!
    //  vertex->SetWeight(particle_weight);
    evt->AddPrimaryVertex(vertex);
    if (verbosityLevel > 1)
        G4cout << " Primary Vetex generated !" << G4endl;
     */
    
}


void toROOTGPSPrimaryGeneratorAction::GetInputFiles(G4String filename)
{
	G4cout << G4endl << " ------ INFO ------ from toROOTGPSPrimaryGeneratorAction::GetInputFiles " << (void*)(this) << G4endl;
    
    // for MT, thread id used in the name of the file to avoid pb
    std::ostringstream fullfilename;
    fullfilename.clear();
    //
#if G4MULTITHREADED
    G4AutoLock lock(&buildMutex); // better to be protected since root is doing some stuff in global
    G4int thread_id = 0;
    thread_id = G4Threading::G4GetThreadId();
    fullfilename << filename << "_Thread" << std::setfill('0') << std::setw(2) << thread_id << ".gene";
#else
    fullfilename << filename << ".gene";
#endif
    //
	// open the ascii file
	ifstream file; file.open(fullfilename.str().data());
	if ( file.is_open() == false ) {
		G4cout << " ** WARNING ** cannot open file " << fullfilename << " (Default parameters are used) "<< G4endl;
		return;
	}
    
    // read the file and get the cascade
    std::string aline, key, rfile; std::getline (file,aline);
	while ( file.good() ) {
        
		if ( aline[0] == '#' ) {
            std::getline (file,aline);
            continue;
        } // this line is a comment
        
        // file to be added to the chain
        std::istringstream decode(aline);
        decode >> key >> rfile;
        
        if ( key == "tree_name" ) {
            theChainOfPrimaryEvents = new TChain(rfile.data());
            std::getline (file,aline);
            
            G4cout << " ==> TTree expected name is " << rfile.data() << G4endl;
            continue;
        }
        // default is ROOTGPS
        if ( theChainOfPrimaryEvents == 0x0 ) {
            theChainOfPrimaryEvents = new TChain("ROOTGPS");
            G4cout << " ==> TTree expected name is ROOTGPS " << G4endl;
        }
		Int_t added = theChainOfPrimaryEvents->AddFile(rfile.data());
        if ( added ) {
            G4cout << " ==> add the file to TChain " << rfile.data() << G4endl;
        }
        else
            G4cout << " !!! Cannot add the file to TChain " << rfile.data() << G4endl;

        std::getline (file,aline);

    }
    // theChainOfPrimaryEvents->Print();
    theChainOfPrimaryEvents->SetBranchAddress("Pr.",&fPr);
    fEntries = theChainOfPrimaryEvents->GetEntries();

	G4cout << " ------ END ------ from toROOTGPSPrimaryGeneratorAction::GetInputFiles " << G4endl << G4endl;
}

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

toROOTGPSPrimaryGeneratorActionMessanger::toROOTGPSPrimaryGeneratorActionMessanger(toROOTGPSPrimaryGeneratorAction *thegene): theGenerator(thegene)
{
	theDirectory = new G4UIdirectory("/toROOTGPS/generator/");
	theDirectory->SetGuidance("To modify generator's parameters");
	
	resetParametersCmd = new G4UIcmdWithAString("/toROOTGPS/generator/file", this);
	resetParametersCmd->SetGuidance("To load new parameters of toROOTGPSPrimaryGeneratorAction from a file");
	resetParametersCmd->SetGuidance("Required parameters: one valid filename (string)");
	resetParametersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

toROOTGPSPrimaryGeneratorActionMessanger::~toROOTGPSPrimaryGeneratorActionMessanger()
{
	delete theDirectory;
    delete resetParametersCmd;
}

void toROOTGPSPrimaryGeneratorActionMessanger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if( command == resetParametersCmd )
        theGenerator->GetInputFiles( newValue );
}


