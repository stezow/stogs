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

#include "SToGS_G4_GeneralPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4DecayPhysics.hh"
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4OpticalPhoton.hh"

SToGS::GeneralPhysics::GeneralPhysics(const G4String& name) : G4VPhysicsConstructor(name)
{
}

SToGS::GeneralPhysics::~GeneralPhysics()
{
}

void SToGS::GeneralPhysics::ConstructParticle()
{
    // pseudo-particles
    G4Geantino::GeantinoDefinition();
    G4ChargedGeantino::ChargedGeantinoDefinition();
    // to avoid tracking optical photon in SD detector, add this here event if optical physics not loaded
    // otherwise, crash in MT mode due to particle created out of pre_init state
    G4OpticalPhoton::OpticalPhotonDefinition();
}

void SToGS::GeneralPhysics::ConstructProcess()
{
    /*
    // Add Decay Process
    fDecayProcess = new G4Decay();
    //
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        if (fDecayProcess->IsApplicable(*particle)) {
            pmanager ->AddProcess(fDecayProcess);
            // set ordering for PostStepDoIt and AtRestDoIt
            pmanager ->SetProcessOrdering(fDecayProcess, idxPostStep);
            pmanager ->SetProcessOrdering(fDecayProcess, idxAtRest);
        }
    }
     */ 
}


