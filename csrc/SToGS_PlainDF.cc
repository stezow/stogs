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
#include "G4GenericTrap.hh"
#include "G4Polycone.hh"

#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"

#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"

#include <stdlib.h>

// Project includes
#include "SToGS_PlainDF.hh"

#include "SToGS_DetectorFactory.hh"
#include "SToGS_MaterialConsultant.hh"
#include "SToGS_UserActionManager.hh"

using namespace std;

// list of all specific factories
namespace  {
    // all Plain Detectors
    SToGS::PlainDF thePlainFactory("PlainDetectors/");
}

G4VPhysicalVolume * SToGS::PlainDF::Make(G4String detname, G4String opt)
{
    
    // **************************************************************************
    // *                              the WORLD                                 *
    // **************************************************************************
    
    G4VPhysicalVolume *theDetector = 0x0;
    
    const G4double world_x = 10.*CLHEP::cm;
    const G4double world_y = 10.*CLHEP::cm;
    const G4double world_z = 10.*CLHEP::cm;
    
    // use a physical as a container to describe the detector
    
    G4Box *detWorld= new G4Box(detname,world_x,world_y,world_z);
    G4LogicalVolume *detlogicWorld= new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
    
    detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
    
    //  Must place the World Physical volume unrotated at (0,0,0).
    theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    detname,      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number
    
    // **************************************************************************
    // *                              the geometries                            *
    // **************************************************************************
    
    
    //cube geometry
    G4double cube_x = 1.*CLHEP::mm;
    G4double cube_y = 1.*CLHEP::mm;
    G4double cube_z = 1.*CLHEP::mm;
    
    //tube geometry
    G4double innerRadius = 10.*CLHEP::mm;
    G4double outerRadius = 20.*CLHEP::mm;
    G4double hz = 20.*CLHEP::mm;
    G4double startAngle = 0.*CLHEP::deg;
    G4double spanningAngle = 360.*CLHEP::deg;
    
    //sphere geometry
    G4double sphereRmin= 0.* CLHEP::mm;
    G4double sphereRmax= 100*CLHEP::mm;
    G4double sphereSPhi= 0.*CLHEP::mm;
    G4double sphereDPhi= 360.*CLHEP::deg;//2*Pi
    G4double sphereSTheta = 0 *CLHEP::deg;
    G4double sphereDTheta = 180 *CLHEP::deg;
    
    
    G4VSolid *plaindetector = 0x0; G4VPhysicalVolume *plain_physical = 0x0; //it means is a pointer
    G4String tmp;
    tmp = detname;
    tmp += "_Shape";
    if (detname=="PLAIN_cube")
    {
        plaindetector= new G4Box(tmp,cube_x,cube_y,cube_z);        
    }
    if (detname=="PLAIN_tube")
    {
        plaindetector = new G4Tubs(tmp,
                                   innerRadius,
                                   outerRadius,
                                   hz,
                                   startAngle,
                                   spanningAngle);
    }
    
    if (detname=="PLAIN_sphere")
    {
        plaindetector= new G4Sphere (tmp,
                                     sphereRmin,
                                     sphereRmax,
                                     sphereSPhi,
                                     sphereDPhi,
                                     sphereSTheta,
                                     sphereDTheta );
    }
    tmp = detname;
    tmp += "_LV";
    G4LogicalVolume *plain_logic= new G4LogicalVolume(plaindetector, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), tmp, 0, 0, 0);
    plain_logic->SetSensitiveDetector( SToGS::UserActionInitialization::GetTrackerSD() );
    
    G4VisAttributes *plain_logicVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
    plain_logic->SetVisAttributes(plain_logicVisAtt);
    
    //  Must place the Detector Physical volume unrotated at (0,0,0).
    plain_physical = new G4PVPlacement(0,         // no rotation
                                       G4ThreeVector(), // at (0,0,0)
                                       plain_logic,      // its logical volume
                                       "plain_P",      // its name
                                       detlogicWorld,               // its mother  volume
                                       false,           // no boolean operations
                                       0);              // copy number
    
    
    
    
    // **************************************************************************
    // *                              the CUBE                                 *
    // **************************************************************************
    
    
    /*
     //box geometry for test
     const G4double cube_x = 1.*CLHEP::mm;
     const G4double cube_y = 1.*CLHEP::mm;
     const G4double cube_z = 1.*CLHEP::mm;
     
     G4VPhysicalVolume *thecube = 0x0; //it means is a pointer
     
     // G4bool do_caps = false, do_housing = false;
     
     // use a physical as a container to describe the detector
     
     G4Box *detCube= new G4Box("detcube",cube_x,cube_y,cube_z);
     G4LogicalVolume *detlogicCube= new G4LogicalVolume(detCube, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), "det_log", 0, 0, 0);
     
     detlogicCube->SetSensitiveDetector( SToGS::UserActionInitialization::GetTrackerSD() );
     
     G4VisAttributes *detlogicCubeVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
     detlogicCube->SetVisAttributes(detlogicCubeVisAtt);
     
     //  Must place the Cube Physical volume unrotated at (0,0,0).
     thecube = new G4PVPlacement(0,         // no rotation
     G4ThreeVector(), // at (0,0,0)
     detlogicCube,      // its logical volume
     "det_P",      // its name
     detlogicWorld,               // its mother  volume
     false,           // no boolean operations
     0);              // copy number
     
     
     // **************************************************************************
     // *                              the TUBE                                 *
     // **************************************************************************
     
     
     //make a cylinder
     
     G4double innerRadius = 10.*CLHEP::mm;
     G4double outerRadius = 20.*CLHEP::mm;
     G4double hz = 20.*CLHEP::mm;
     G4double startAngle = 0.*CLHEP::deg;
     G4double spanningAngle = 360.*CLHEP::deg;
     
     G4VPhysicalVolume *theTube = 0x0;
     
     
     G4Tubs *detTube = new G4Tubs("dettube",
     innerRadius,
     outerRadius,
     hz,
     startAngle,
     spanningAngle);
     
     G4LogicalVolume *detlogicTube= new G4LogicalVolume(detTube, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), "tube_log", 0, 0, 0);
     
     
     G4VisAttributes *detlogicTubeVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0)); //cyan
     detlogicTube->SetVisAttributes(detlogicTubeVisAtt);
     detlogicTube->SetSensitiveDetector( SToGS::UserActionInitialization::GetTrackerSD() );
     
     theTube = new G4PVPlacement(0,         // no rotation
     G4ThreeVector(0,0,0), // at (0,0,0)
     detlogicTube,      // its logical volume
     "tube_P",      // its name
     detlogicWorld,               // its mother  volume
     false,           // no boolean operations
     0);              // copy number
     
     
     // **************************************************************************
     // *                              the SPHERE                                *
     // **************************************************************************
     
     //make a sphere
     G4double sphereRmin= 0.* CLHEP::mm;
     G4double sphereRmax= 100*CLHEP::mm;
     G4double sphereSPhi= 0.*CLHEP::mm;
     G4double sphereDPhi= 360.*CLHEP::deg;//2*Pi
     G4double sphereSTheta = 0 *CLHEP::deg;
     G4double sphereDTheta = 180 *CLHEP::deg;
     
     
     G4VPhysicalVolume *theSphere = 0x0;
     
     G4Sphere *detSphere= new G4Sphere ("detsphere",
     sphereRmin,
     sphereRmax,
     sphereSPhi,
     sphereDPhi,
     sphereSTheta,
     sphereDTheta );
     
     
     G4LogicalVolume *detlogicSphere= new G4LogicalVolume(detSphere, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), "sphere_log", 0, 0, 0);
     
     
     G4VisAttributes *detlogicSphereVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0)); //green
     detlogicSphere->SetVisAttributes(detlogicSphereVisAtt);
     detlogicSphere->SetSensitiveDetector( SToGS::UserActionInitialization::GetTrackerSD() );
     
     theSphere = new G4PVPlacement(0,         // no rotation
     G4ThreeVector(0,0,0), // at (0,0,0)
     detlogicSphere,      // its logical volume
     "sphere_P",      // its name
     detlogicWorld,               // its mother  volume
     false,           // no boolean operations
     0);              // copy number
     
     
     */
    
    return theDetector;
}
/*
 G4VPhysicalVolume * SToGS::PlainDF::Make(G4String name, G4String version_string)
 {
 G4VPhysicalVolume *theDetector = 0x0; G4String detname;
 
 if ( name == "PLAIN" ) {
 detname = GetDetName("PLAIN",version_string);
 theDetector =
 MakePLAIN(detname,version_string);
 }
 
 
 
 return theDetector;
 }
 */
void SToGS::PlainDF::MakeStore()
{
    
    
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("PLAIN_cube",""); //
    
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("PLAIN_tube","");//
    
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("PLAIN_sphere","");//
    
}




