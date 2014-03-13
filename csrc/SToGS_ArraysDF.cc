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
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4Sphere.hh"

#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"

#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4SubtractionSolid.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"

// Project includes
#include "SToGS_ArraysDF.hh"

// list of all specific factories
namespace  {
    // all arrays
    SToGS::ArraysDF theArrays("Arrays/");
}

void SToGS::ArraysDF::MakeStore()
{
    //    ParisOutputManager::SetGCopyNb(0);
    //    MakeInStore("Chateau2Crystal","bare");
    ParisOutputManager::SetGCopyNb(0);
    MakeInStore("Chateau2Crystal","");
}

G4VPhysicalVolume *SToGS::ArraysDF::Make(G4String name, G4String version_string)
{
    G4VPhysicalVolume *theDetector = 0x0; G4String detname;

    if ( name == "Chateau2Crystal" ) {
        detname = GetDetName("Chateau2Crystal",version_string);
        theDetector = MakeAnArrayFromFactory("DetectorFactory/Arrays/Builders/Chateau2Crystal.dfb");
    }
    return theDetector;
}




