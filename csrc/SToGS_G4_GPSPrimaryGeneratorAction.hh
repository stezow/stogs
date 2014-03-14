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

#ifndef SToGS_G4_GPSPrimaryGeneratorAction_h
#define SToGS_G4_GPSPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4String.hh"

class G4GeneralParticleSource;
class G4Event;

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! This class is just calling the GPS interface - Could be exported to any other G4 Applications
    /*!
     */
    class GPSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
    {
    public:
        GPSPrimaryGeneratorAction(G4String mac = "");
        virtual ~GPSPrimaryGeneratorAction();
    public:
        virtual void GeneratePrimaries(G4Event* anEvent);
        
    private:
        G4GeneralParticleSource* particleGun;
    };
} // SToGS Namespace

#endif


