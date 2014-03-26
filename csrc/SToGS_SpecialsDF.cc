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

// G4 includes
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"

// Project includes
#include "SToGS_SpecialsDF.hh"
//
#include "SToGS_TwoShellsDetectorConstruction.hh"


// list of all specific factories
namespace  {
    // all scintillators
    SToGS::SpecialsDF theGenericFactory("Generics/");
}

G4VPhysicalVolume *SToGS::SpecialsDF::Get(G4String basename, G4bool is_full)
{
    G4VUserDetectorConstruction *theConstructor = 0x0; G4VPhysicalVolume *theDetector = 0x0;
    
    // check is already loaded
    for (size_t i = 0; i < fLoadedPhysical.size(); i++) {
        if ( fLoadedPhysical[i].first == basename ) {
            theDetector = fLoadedPhysical[i].second;
            break;
        }
    }
    if ( theDetector ) {
        return theDetector;
    }
    // load using the g4 facility
    if ( basename.contains("TwoShells") ) {
        theConstructor = new SToGS::TwoShellsDetectorConstruction();
        theDetector = theConstructor->Construct();
    }
    if ( theDetector ) {
        std::pair < G4String, G4VPhysicalVolume *> p1(basename,theDetector); // add the new loaded detector to the list
        fLoadedPhysical.push_back(p1);
        std::pair < G4String, G4VUserDetectorConstruction *> p2(basename,theConstructor);
        fLoadedUserDetectorConstruction.push_back(p2);
        
        if ( is_full ) {
#if G4VERSION_NUMBER < 1000
#else
            theConstructor->ConstructSDandField();
#endif
        }
    }

    return theDetector;
}

void SToGS::SpecialsDF::GetAttributes(G4String basename)
{
    // check is already loaded
#if G4VERSION_NUMBER < 1000
#else
    for (size_t i = 0; i < fLoadedUserDetectorConstruction.size(); i++) {
        if ( fLoadedUserDetectorConstruction[i].first == basename ) {
            fLoadedUserDetectorConstruction[i].second->ConstructSDandField(basename);
            break;
        }
    }
#endif
}






