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

#include "SToGS_G4_GPSPrimaryGeneratorAction.hh"
#ifdef TO_ROOTGPS
#include "SToGS_G4_ROOTGeneralParticleSource.hh"
#endif

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4UImanager.hh"

// to protect the classic GPS during the time required to set
//#define TO_ROOTGPS 1

SToGS::GPSPrimaryGeneratorAction::GPSPrimaryGeneratorAction(G4String mac)
{
#ifdef TO_ROOTGPS
    // ROOTGPS + init à la main
    particleGun = new SToGS::ROOTGeneralParticleSource();
#else
    particleGun = new G4GeneralParticleSource();
    if ( mac != "" ) {
        G4bool is_mac = false; std::ifstream f(mac.c_str());
        if (f.good()) {
            f.close();
            is_mac = true;
        } else {
            f.close();
        }
        if ( is_mac ) {
            G4String command = "/control/execute ";
            command += mac;
            G4UImanager::GetUIpointer()->ApplyCommand(command);
        }
    }
#endif
}

SToGS::GPSPrimaryGeneratorAction::~GPSPrimaryGeneratorAction()
{
    delete particleGun;
}

void SToGS::GPSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    particleGun->GeneratePrimaryVertex(anEvent) ;
}






