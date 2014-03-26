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


// Paris includes
#include "SemiConductorGeDF.hh"
#include "ParisMaterialConsultant.hh"
#include "ParisOutputManager.hh"

#include "G4GDMLParser.hh"

using namespace std;


// list of all specific factories
namespace  {
    // all SemiConductors
    SemiConductorGeDF theGeFactory("SemiConductors/Ge/");
}

/*
// Germanium isotopes
G4Isotope* Ge70 = new G4Isotope(name="Ge70", nprot=32, nnucl=70, a=69.9242*g/mole);
G4Isotope* Ge72 = new G4Isotope(name="Ge72", nprot=32, nnucl=72, a=71.9221*g/mole);
G4Isotope* Ge73 = new G4Isotope(name="Ge73", nprot=32, nnucl=73, a=72.9235*g/mole);
G4Isotope* Ge74 = new G4Isotope(name="Ge74", nprot=32, nnucl=74, a=73.9212*g/mole);
G4Isotope* Ge76 = new G4Isotope(name="Ge76", nprot=32, nnucl=76, a=75.9214*g/mole);
// germanium defined via its isotopes
G4Element* elGe = new G4Element(name="Germanium",symbol="Ge", nelements=5);
elGe->AddIsotope(Ge70, 0.2123);
elGe->AddIsotope(Ge72, 0.2766);
elGe->AddIsotope(Ge73, 0.0773);
elGe->AddIsotope(Ge74, 0.3594);
elGe->AddIsotope(Ge76, 0.0744);
myElements.push_back(elGe);
 */

G4VPhysicalVolume * SemiConductorGeDF::MakeEURO_PI(G4String detname, G4String opt)
{
    G4VPhysicalVolume *theDetector = 0x0; G4LogicalVolume *detlogicWorld; G4Box *detWorld; // G4bool do_caps = false, do_housing = false;

    // Option
    if ( opt != "bare" ) {
//            do_caps = true;
    }
    /*
    if ( housing_width != 0.0 ) {
        if ( opt != "bare" ) {
            do_housing = true;
        }
    }
     */
    
    // use a physical as a container to describe the detector
	detWorld= new G4Box(detname,10.*cm,10.*cm,50.*cm);
	detlogicWorld= new G4LogicalVolume(detWorld, ParisMaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
	
	detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	//  Must place the World Physical volume unrotated at (0,0,0).
	theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    detname,      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number

    // Fron drawings
    const G4double TronCrystalLength	= 70.0*mm; // cystal length

	const G4double TronCrystalEdgeDepth	= 40.0*mm; // depth at which starts TronCrystalRadiusMax
	const G4double TronCrystalRadiusMin	= 32.0*mm; // radius of the crystal at the entrance face
	const G4double TronCrystalRadiusMax	= 35.0*mm; // radius of the crystal at the back face
	
	const G4double TronCrystalHoleDepth	 = 15.0*mm; // depth at which starts the hole
	const G4double TronCrystalHoleRadius =  5.0*mm; // radius of the hole


	G4int nbSlice = 5;
	G4double zSlice[5] = {
        0.0*mm,
        TronCrystalHoleDepth, TronCrystalHoleDepth + 1.0*mm,
        TronCrystalEdgeDepth, TronCrystalLength };
    
	G4double InnRad[5] = {
        0.0*mm,
        0.0*mm, TronCrystalHoleRadius,
        TronCrystalHoleRadius, TronCrystalHoleRadius };
    
	G4double OutRad[5] = {
        TronCrystalRadiusMin,
        TronCrystalRadiusMin + (TronCrystalRadiusMax-TronCrystalRadiusMin)*TronCrystalHoleDepth/TronCrystalEdgeDepth,
        TronCrystalRadiusMin + (TronCrystalRadiusMax-TronCrystalRadiusMin)*(TronCrystalHoleDepth+1.0*mm)/TronCrystalEdgeDepth,
        TronCrystalRadiusMax, TronCrystalRadiusMax };
    
    // the Germanim crystal
    G4Polycone *crys =
        new G4Polycone("CrystalEUROGAM_PI",0.0*deg,360.0*deg,nbSlice,zSlice,InnRad,OutRad);
    // G4LogicalVolume *crys_logical  = new G4LogicalVolume( crys, matCrystal, G4String(sName), 0, 0, 0 );
    
    
    return theDetector;
}

G4VPhysicalVolume * SemiConductorGeDF::Make(G4String name, G4String version_string)
{
    G4VPhysicalVolume *theDetector = 0x0; G4String detname;
    
    // PARIS RELATED
    if ( name == "EuroPI" ) {
        
        detname = GetDetName("EuroPI",version_string);
        theDetector =
            MakeEURO_PI(detname,version_string);
    }
    
    return theDetector;
}

void SemiConductorGeDF::MakeStore()
{
// EUROGAM/BALL related
    ParisOutputManager::SetGCopyNb(0);
    //MakeInStore("ParisPW","2");
}



