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
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"


// Project includes
#include "SToGS_ScintillatorDF.hh"

#include "SToGS_DetectorFactory.hh"
#include "SToGS_MaterialConsultant.hh"
#include "SToGS_UserActionManager.hh"

// list of all specific factories
namespace  {
    // all scintillators
    SToGS::ScintillatorDF theScintillatorFactory("Scintillators/");
}

void SToGS::ScintillatorDF::MakeStore()
{
    // PARIS RELATED
    // SToGS::DetectorFactory::theMainFactory()->Clean();
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("ParisPW","2");
    // SToGS::DetectorFactory::theMainFactory()->Clean();
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("ParisPW","2-bare");
    //
    // SToGS::DetectorFactory::theMainFactory()->Clean();
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("ParisPW","3");
    // SToGS::DetectorFactory::theMainFactory()->Clean();
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("ParisPW","3-bare");
    //
    // SToGS::DetectorFactory::theMainFactory()->Clean();
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("CParisPW","2");
    // SToGS::DetectorFactory::theMainFactory()->Clean();
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("CParisPW","2-bare");
    
    // Chateau de Crystal
    // SToGS::DetectorFactory::theMainFactory()->Clean();
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("aChateau2Crystal","bare");
    // SToGS::DetectorFactory::theMainFactory()->Clean();
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("aChateau2Crystal","");
    
    // FATIMA
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("FATIMAM","bare");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("FATIMAQ","bare");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("FATIMAM","");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("FATIMAQ","");

    // EDEN
    // SToGS::DetectorFactory::theMainFactory()->Clean();
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("EDEN","v0");
    // SToGS::DetectorFactory::theMainFactory()->Clean();
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("EDEN","bare");
}

G4VPhysicalVolume * SToGS::ScintillatorDF::Make(G4String name, G4String version_string)
{
    G4VPhysicalVolume *theDetector = 0x0; G4String detname;
    
    // PARIS RELATED
    if ( name == "ParisPW" ) {
        
        detname = GetDetName("ParisPW",version_string);
        
        if ( version_string.contains("2") ) {
            if ( version_string == "2" ) {
                theDetector = MakePPW(detname,0.5*CLHEP::mm,0.0*CLHEP::mm,"");
            }
            if ( version_string.contains("2-bare") ) {
                theDetector = MakePPW(detname,0.5*CLHEP::mm,0.0*CLHEP::mm,"bare");
            }
        }
        if ( version_string.contains("3") ) {
            if ( version_string == "3" ) {
                theDetector = MakePPW(detname,1.0*CLHEP::mm,0.0*CLHEP::mm,"");
            }
            if ( version_string.contains("3-bare") ) {
                theDetector = MakePPW(detname,1.0*CLHEP::mm,0.0*CLHEP::mm,"bare");
            }
        }
    }
    if ( name == "CParisPW" ) {
        detname = GetDetName("CParisPW",version_string); // version string should be the same that ParisPW
        theDetector = MakeCPPW(detname,version_string);
    }
    // Chateau de Crystal
    if ( name == "aChateau2Crystal" ) {
        detname = GetDetName("aChateau2Crystal",version_string);
        theDetector = MakeCdC(detname,version_string);
    }
    
    // FATIMA
    if ( name.index("FATIMA") == 0 ) {
        if ( name == "FATIMAM" ) {
            detname = GetDetName("FATIMAM",version_string);
            theDetector = MakeFATIMAM(detname,version_string);
        }
        if ( name == "FATIMAQ" ) {
            detname = GetDetName("FATIMAQ",version_string);
            theDetector = MakeFATIMAQ(detname,version_string);
        }
    }

    // EDEN
    if ( name == "EDEN" ) {
        detname = GetDetName("EDEN",version_string);
        theDetector = MakeEDEN(detname,1.*CLHEP::mm, 10.*CLHEP::cm ,version_string);
    }
    
    return theDetector;
}


G4VPhysicalVolume *SToGS::ScintillatorDF::MakeCdC(G4String detname, G4String opt)
{
    G4VPhysicalVolume *theDetector = 0x0; G4LogicalVolume *detlogicWorld; G4Box *detWorld; G4bool do_caps = true;
    
    // from manifacture
	const G4double caps_width = 0.9*CLHEP::mm, teflon_width = 0.1*CLHEP::mm, crystal_radius = 5.0*CLHEP::cm, crystal_depth = 14*CLHEP::cm;
    const G4int numSide = 6, numZPlanes = 2;
    
	// for G4 radius is inner radius ==> ri = sqrt
	G4double inner_crystal_radius = crystal_radius * ::sqrt(3) / 2.;
	G4double
        z[]	     = { caps_width + teflon_width, crystal_depth + caps_width + teflon_width },
        rInner[] = { 0., 0. },
        rOuter[] = { inner_crystal_radius, inner_crystal_radius },
    
        z_caps_front[]	    = { 0.0, caps_width + teflon_width },
        rInner_caps_front[] = { 0., 0.},
        rOuter_caps_front[] = { inner_crystal_radius, inner_crystal_radius },
    
        z_caps[]      = { 0., crystal_depth + 2.0*(caps_width + teflon_width) },
        rInner_caps[] = { inner_crystal_radius, inner_crystal_radius },
        rOuter_caps[] = { inner_crystal_radius + teflon_width + caps_width, inner_crystal_radius + teflon_width + caps_width };
    
    // if option bare does not include the encapsulation
    if ( opt == "bare" ) {
        do_caps = false;
    }
    
    // use a physical as a container to describe the detector
	detWorld= new G4Box(detname,15.*CLHEP::cm,15.*CLHEP::cm,25.*CLHEP::cm);
	detlogicWorld= new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
	
	detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	//  Must place the World Physical volume unrotated at (0,0,0).
	theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    detname,      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number
    
	// crystal
	G4Polyhedra* a_solid
        = new G4Polyhedra("ShapeC2CCrys",0.,360.1*CLHEP::deg,numSide,numZPlanes,z,rInner,rOuter);
	//Setting a logical volume
	G4LogicalVolume *a_log
        = new G4LogicalVolume(a_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("BaF2"),"C2CCrysLV",0,0,0);
	
	G4VisAttributes *visatt = new G4VisAttributes( G4Colour(0.0, 0.6, 0.0) );
	visatt->SetVisibility(true);
	a_log->SetVisAttributes( visatt );
    a_log->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );
    
	// Rotation and translation to place each peace in the assembly
	G4ThreeVector Ta;
	
	// Add crystals Inner
	Ta.setX( 0.0 );
    Ta.setY( 0.0 );
    Ta.setZ( 0.0 );
    new G4PVPlacement(0,Ta,a_log,"C2CCrys",detlogicWorld,false,0);

    // encapsulation side
    G4VisAttributes *Capsule_visatt;
    a_solid = new G4Polyhedra("ShapeC2CCapsSide",0.,360.*CLHEP::deg,numSide,numZPlanes,z_caps,rInner_caps,rOuter_caps);
    if ( do_caps ) {
        a_log
            = new G4LogicalVolume(a_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("Al"),"C2CCapsSideLV",0,0,0);
        // grey for all passive part
        Capsule_visatt = new G4VisAttributes( G4Colour(0.8, 0.8, 0.8, 0.75) );
        Capsule_visatt->SetVisibility(true);
        a_log->SetVisAttributes( Capsule_visatt );
    }
    else {
        a_log
            = new G4LogicalVolume(a_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"),"C2CCapsSideLV",0,0,0);
        // white for all passive part if air
        Capsule_visatt = new G4VisAttributes( G4Colour(1, 1, 1, 0.5) );
        Capsule_visatt->SetVisibility(true);
        a_log->SetVisAttributes( Capsule_visatt );
    }
    //Setting a logical volume
    Ta.setX( 0.0 );
    Ta.setY( 0.0 );
    Ta.setZ( 0.0 );
    new G4PVPlacement(0,Ta,a_log,"C2CCapsSide",detlogicWorld,false,-1);
    
    // encapsulation front and back
    a_solid = new G4Polyhedra("C2CCapsFront",0.,360.*CLHEP::deg,numSide,numZPlanes,z_caps_front,rInner_caps_front,rOuter_caps_front);
    //Setting a logical volume
    if ( do_caps ) {
        a_log
            = new G4LogicalVolume(a_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("Al"),"C2CCapsFrontLV",0,0,0);
        // grey for all passive part
        a_log->SetVisAttributes( Capsule_visatt );
    }
    else {
        a_log
            = new G4LogicalVolume(a_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"),"C2CCapsFrontLV",0,0,0);
        // grey for all passive part
        a_log->SetVisAttributes( Capsule_visatt );
    }
    Ta.setX( 0.0 );
    Ta.setY( 0.0 );
    Ta.setZ( 0.0 );
    new G4PVPlacement(0,Ta,a_log,"C2CCapsFront",detlogicWorld,false,-1);
    Ta.setX( 0.0 );
    Ta.setY( 0.0 );
    Ta.setZ( crystal_depth + caps_width + teflon_width );
    new G4PVPlacement(0,Ta,a_log,"C2CCapsBack",detlogicWorld,false,-1);
	
    return theDetector;
}


// PARIS phoswitch - Version 0
// capsule = polyhedra closed at the front 
// cristals are attached to the world
// this shape may be the cause of singularities when the detector is used in simulation 
// (pathologic treatment of projectiles impacting at (x=0, y=0), 
// namely in the case of a centered beam or a point source emitting in z direction)

/*
G4VPhysicalVolume *SToGS::ScintillatorDF::MakePPW(G4String detname, G4double caps_width, G4double housing_width, G4String opt)
{
    G4VPhysicalVolume *theDetector = 0x0; G4LogicalVolume *detlogicWorld; G4Box *detWorld; G4bool do_caps = false, do_housing = false;
    G4RotationMatrix R;
	G4ThreeVector T;
    G4Transform3D Tr;
    
    // Option
    if ( caps_width != 0.0 ) {
        if ( opt != "bare" ) {
            do_caps = true;
        }
    }
    if ( housing_width != 0.0 ) {
        if ( opt != "bare" ) {
            do_housing = true;
        }
    }
    
    // use a physical as a container to describe the detector
	detWorld= new G4Box(detname,10.*CLHEP::cm,10.*CLHEP::cm,50.*CLHEP::cm);
	detlogicWorld= new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
	
	detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	//  Must place the World Physical volume unrotated at (0,0,0).
	theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    detname,      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number
    
	// material for first and second layer
	const char *matInner = "LaBr3", *matOuter = "NaI", *matCapsule = "Al", *matHousing = "Al";
	
	// full lenght in X,Y,Z
	G4double
        InnerX = 2*2.54*CLHEP::cm, InnerY = 2*2.54*CLHEP::cm, InnerZ = 2*2.54*CLHEP::cm,
        OuterX = 2*2.54*CLHEP::cm, OuterY = 2*2.54*CLHEP::cm, OuterZ = 6*2.54*CLHEP::cm;
    // in case of housing longer at the back of the PW
	G4double ExtraHousingBack = 20*CLHEP::mm, eps_caps =  caps_width / 100. , eps_housing = 0.05*CLHEP::mm; // eps remove in the capsule part to avoid overlapping
    // G4double r_caps_width = caps_width - eps_caps, r_housing_width = housing_width - eps_housing;
	
    
	// Stage 0
	G4Box *Inner_solid = new G4Box("ShapePW0", InnerX/2., InnerY/2., InnerZ/2.);
	G4LogicalVolume *Inner_logic =
        new G4LogicalVolume(Inner_solid, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matInner),"PWLV:0",0,0,0);
	G4VisAttributes *Inner_visatt = new G4VisAttributes( G4Colour(0.0, 0.0, 1.0) );
	Inner_visatt->SetVisibility(true);
	Inner_logic->SetVisAttributes( Inner_visatt );
    Inner_logic->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );
    
	T.setX( 0.0 );
    T.setY( 0.0 );
    T.setZ( housing_width + caps_width + InnerZ / 2. );
    
    new G4PVPlacement(0,T,Inner_logic,"PW:0:",detlogicWorld,false,0);
	
	// Stage 1
	G4Box *Outer_solid = new G4Box("ShapePW1", OuterX/2., OuterY/2., OuterZ/2.);
	G4LogicalVolume *Outer_logic =
        new G4LogicalVolume(Outer_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matOuter), "PWLV:1",0,0,0);
	G4VisAttributes *Outer_visatt = new G4VisAttributes( G4Colour(1.0, 0.0, 0.0) );
	Outer_visatt->SetVisibility(true);
	Outer_logic->SetVisAttributes( Outer_visatt );
    Outer_logic->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );
    
	T.setX( 0.0 );
    T.setY( 0.0 );
    T.setZ( housing_width + caps_width + InnerZ + OuterZ / 2. );
    new G4PVPlacement(0,T,Outer_logic,"PW:1:",detlogicWorld,false,1);
    
    // Encapsulation Al if on otherwise set as AIR ... useful to check overlapping
    G4double zplane[3], rinner[3], router[3];
    
    zplane[0] = 0.01*CLHEP::mm; zplane[1] = caps_width - eps_caps; zplane[2] = InnerZ + OuterZ + ExtraHousingBack;
    rinner[0] = 0.00*CLHEP::mm; rinner[1] = InnerX / 2 + eps_caps; rinner[2] = InnerX / 2 + eps_caps;
    router[0] = InnerX / 2 + caps_width ; router[1] = InnerX / 2 + caps_width ; router[2] = InnerX / 2 + caps_width;
    
    //rinner[0] = rinner[1] = rinner[2] = InnerX / 2 + eps_caps;
    //router[0] = router[1] = router[2] = InnerX / 2 + caps_width;
    
    G4Polyhedra *side_part =
        new G4Polyhedra("ShapePWCaps", 0.0, 361.0*CLHEP::deg, 4, 3, zplane, rinner, router);
    
    G4LogicalVolume *Capsule_logic_side; G4VisAttributes *Capsule_visatt;
    if ( do_caps ) {
        Capsule_logic_side =
            new G4LogicalVolume(side_part,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matCapsule),"PWCaps",0,0,0);
        
        Capsule_visatt = new G4VisAttributes( G4Colour(0.8, 0.8, 0.8, 0.75) );
        Capsule_visatt->SetVisibility(true);
        Capsule_logic_side->SetVisAttributes( Capsule_visatt );
    }
    else {
        Capsule_logic_side =
            new G4LogicalVolume(side_part,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"),"PWCaps",0,0,0);
        
        Capsule_visatt = new G4VisAttributes( G4Colour(1, 1, 1, 0.5) );
        Capsule_visatt->SetVisibility(true);
        Capsule_logic_side->SetVisAttributes( Capsule_visatt );
    }
    
    T.setX( 0.0 );
    T.setY( 0.0 );
    T.setZ( 0.0 );
    G4RotationMatrix *Ro = new G4RotationMatrix();
    Ro->rotateZ(45.0*CLHEP::deg);
    new G4PVPlacement(Ro,T,Capsule_logic_side,"PwCaps",detlogicWorld,false,-1);
    
    return theDetector;
}
*/

// PARIS phoswitch - Version 1
// capsule = box, with a box of air inside
// cristals are attached to the box of air inside the capsule
/*
G4VPhysicalVolume *SToGS::ScintillatorDF::MakePPW(G4String detname, G4double caps_width, G4double housing_width, G4String opt)
{
// Note : housing_width is not used in this version

    G4VPhysicalVolume *theDetector = 0x0; G4LogicalVolume *detlogicWorld; G4Box *detWorld; G4bool do_caps = false, do_housing = false;
    G4RotationMatrix R;
	G4ThreeVector T;
    G4Transform3D Tr;
    
    // Option
    if ( caps_width != 0.0 ) {
        if ( opt != "bare" ) {
            do_caps = true;
        }
    }
    if ( housing_width != 0.0 ) {
        if ( opt != "bare" ) {
            do_housing = true;
        }
    }
    
    // use a physical as a container to describe the detector
	detWorld= new G4Box(detname,10.*CLHEP::cm,10.*CLHEP::cm,50.*CLHEP::cm);
	detlogicWorld= new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
	
	detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	//  Must place the World Physical volume unrotated at (0,0,0).
	theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    detname,      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number
    
	// material for first and second layer
	const char *matInner = "LaBr3", *matOuter = "NaI", *matCapsule = "Al", *matHousing = "Al";
	
	// full lenght in X,Y,Z
	G4double
        InnerX = 2*2.54*CLHEP::cm, InnerY = 2*2.54*CLHEP::cm, InnerZ = 2*2.54*CLHEP::cm,
        OuterX = 2*2.54*CLHEP::cm, OuterY = 2*2.54*CLHEP::cm, OuterZ = 6*2.54*CLHEP::cm;
    // in case of housing longer at the back of the PW
	G4double 
		ExtraHousingBack = 20*CLHEP::mm, 
		eps_caps = 0.1*CLHEP::mm, //caps_width / 100., eps remove in the capsule part to avoid overlapping	
		Shift = 0.01*CLHEP::mm; // the detector entrance is not just at 0

    // Encapsulation Al if on, otherwise set as AIR ... useful to check overlapping

	G4Box *Capsule_solid = new G4Box("ShapePWCaps", 
					InnerX/2.+caps_width, 
					InnerY/2.+caps_width, 
					(InnerZ+OuterZ+caps_width+ExtraHousingBack)/2.);
    
    G4LogicalVolume *Capsule_logic; G4VisAttributes *Capsule_visatt;
    if ( do_caps ) {
        Capsule_logic =
            new G4LogicalVolume(Capsule_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matCapsule),"PWCaps",0,0,0);
        
        Capsule_visatt = new G4VisAttributes( G4Colour(0.8, 0.8, 0.8, 0.75) );
        Capsule_visatt->SetVisibility(true);
        Capsule_logic->SetVisAttributes( Capsule_visatt );
    }
    else {
        Capsule_logic =
            new G4LogicalVolume(Capsule_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"),"PWCaps",0,0,0);
        
        Capsule_visatt = new G4VisAttributes( G4Colour(1, 1, 1, 0.5) );
        Capsule_visatt->SetVisibility(true);
        Capsule_logic->SetVisAttributes( Capsule_visatt );
    }

    T.setX( 0.0 );
    T.setY( 0.0 );
    T.setZ( (InnerZ+OuterZ+caps_width+ExtraHousingBack)/2. + Shift);
	new G4PVPlacement(0,T,Capsule_logic,"PwCaps",detlogicWorld,false,-1);

    // Air box in the capsule box to make it hollow

	G4Box * Air_solid = new G4Box("ShapePWAirBox",
					InnerX/2.+eps_caps,
					InnerY/2.+eps_caps,
					(InnerZ+OuterZ+eps_caps+ExtraHousingBack)/2.);

	G4LogicalVolume *Air_logic = 
		new G4LogicalVolume(Air_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"),"PWAir",0,0,0);

	G4VisAttributes *Air_visatt = new G4VisAttributes( G4Colour(1, 1, 1, 0.5) );
	Air_visatt->SetVisibility(true);
    Air_logic->SetVisAttributes( Air_visatt );
	
    T.setX( 0.0 );
    T.setY( 0.0 );
    T.setZ( (caps_width-eps_caps)/2. );
	new G4PVPlacement(0,T,Air_logic,"PwAir",Capsule_logic,false,-1);

	// Stage 0
	G4Box *Inner_solid = new G4Box("ShapePW0", InnerX/2., InnerY/2., InnerZ/2.);
	G4LogicalVolume *Inner_logic =
        new G4LogicalVolume(Inner_solid, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matInner),"PWLV:0",0,0,0);
	G4VisAttributes *Inner_visatt = new G4VisAttributes( G4Colour(0.0, 0.0, 1.0) );
	Inner_visatt->SetVisibility(true);
	Inner_logic->SetVisAttributes( Inner_visatt );
    	Inner_logic->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );
    
	T.setX( 0.0 );
    T.setY( 0.0 );
    //T.setZ( housing_width + caps_width + InnerZ / 2. );
    //new G4PVPlacement(0,T,Inner_logic,"PW:0:",detlogicWorld,false,0);
	T.setZ( (eps_caps-OuterZ-ExtraHousingBack)/2. );	
	new G4PVPlacement(0,T,Inner_logic,"PW:0:",Air_logic,false,0);

	// Stage 1
	G4Box *Outer_solid = new G4Box("ShapePW1", OuterX/2., OuterY/2., OuterZ/2.);
	G4LogicalVolume *Outer_logic =
        new G4LogicalVolume(Outer_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matOuter), "PWLV:1",0,0,0);
	G4VisAttributes *Outer_visatt = new G4VisAttributes( G4Colour(1.0, 0.0, 0.0) );
	Outer_visatt->SetVisibility(true);
	Outer_logic->SetVisAttributes( Outer_visatt );
    	Outer_logic->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );
    
	T.setX( 0.0 );
    T.setY( 0.0 );
    //T.setZ( housing_width + caps_width + InnerZ + OuterZ / 2. );
    //new G4PVPlacement(0,T,Outer_logic,"PW:1:",detlogicWorld,false,1);
    T.setZ( (eps_caps+InnerZ-ExtraHousingBack)/2. );
    new G4PVPlacement(0,T,Outer_logic,"PW:1:",Air_logic,false,1);
    
    return theDetector;
}
*/

// PARIS phoswitch - Version 2
// capsule = polyhedra with a tap (box) in front
// cristals are attached to the world

G4VPhysicalVolume *SToGS::ScintillatorDF::MakePPW(G4String detname, G4double caps_width, G4double housing_width, G4String opt)
{
// Note : housing_width is not used in this version

    G4VPhysicalVolume *theDetector = 0x0; G4LogicalVolume *detlogicWorld; G4Box *detWorld; G4bool do_caps = false, do_housing = false;
    G4RotationMatrix R;
	G4ThreeVector T;
    G4Transform3D Tr;
    
    // Option
    if ( caps_width != 0.0 ) {
        if ( opt != "bare" ) {
            do_caps = true;
        }
    }
    if ( housing_width != 0.0 ) {
        if ( opt != "bare" ) {
            do_housing = true;
        }
    }
    
    // use a physical as a container to describe the detector
	detWorld= new G4Box(detname,10.*CLHEP::cm,10.*CLHEP::cm,50.*CLHEP::cm);
	detlogicWorld= new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
	
	detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	//  Must place the World Physical volume unrotated at (0,0,0).
	theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    detname,      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number
    
	// material for first and second layer
	const char *matInner = "LaBr3", *matOuter = "NaI", *matCapsule = "Al"; // *matHousing = "Al";
	
	// full lenght in X,Y,Z
	G4double
        InnerX = 2*2.54*CLHEP::cm, InnerY = 2*2.54*CLHEP::cm, InnerZ = 2*2.54*CLHEP::cm,
        OuterX = 2*2.54*CLHEP::cm, OuterY = 2*2.54*CLHEP::cm, OuterZ = 6*2.54*CLHEP::cm;
    // in case of housing longer at the back of the PW
	G4double 
		ExtraHousingBack = 20*CLHEP::mm, 
		eps_caps = 0.1*CLHEP::mm, //caps_width / 100., eps remove in the capsule part to avoid overlapping	
		Shift = 0.01*CLHEP::mm; // the detector entrance is not just at 0

    // Encapsulation Al if on, otherwise set as AIR ... useful to check overlapping

/*	G4Box *Capsule_solid = new G4Box("ShapePWCaps", 
					InnerX/2.+caps_width, 
					InnerY/2.+caps_width, 
					(InnerZ+OuterZ+caps_width+ExtraHousingBack)/2.);
*/

    G4double zplane[2], rinner[2], router[2];
    
	zplane[0] = Shift+caps_width-eps_caps; 
	zplane[1] = Shift+caps_width+InnerZ+OuterZ+ExtraHousingBack;
    	rinner[0] = InnerX / 2. + eps_caps; rinner[1] = InnerX / 2. + eps_caps;
    	router[0] = InnerX / 2. + caps_width ; router[1] = InnerX / 2. + caps_width;

	G4Polyhedra *Capsule_solid = new G4Polyhedra("ShapePWCaps", 
					45.0*CLHEP::deg,360.0*CLHEP::deg,
					4,2,zplane,rinner,router
					);
	G4Box *CapsTap_solid = new G4Box("ShapeCapsTap",
					InnerX/2.+caps_width,
					InnerY/2.+caps_width,
					caps_width-eps_caps);

	G4LogicalVolume *Capsule_logic; 
	G4LogicalVolume *CapsTap_logic;
	G4VisAttributes *Capsule_visatt;

    if ( do_caps ) {
        Capsule_logic =
            new G4LogicalVolume(Capsule_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matCapsule),"PWCaps",0,0,0);

	CapsTap_logic =
	    new G4LogicalVolume(CapsTap_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matCapsule),"PWCapsTaps",0,0,0);      

        Capsule_visatt = new G4VisAttributes( G4Colour(0.8, 0.8, 0.8, 0.75) );
        Capsule_visatt->SetVisibility(true);
        Capsule_logic->SetVisAttributes( Capsule_visatt );
        CapsTap_logic->SetVisAttributes( Capsule_visatt );
    }
    else {
        Capsule_logic =
            new G4LogicalVolume(Capsule_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"),"PWCaps",0,0,0);
        
	CapsTap_logic =
	    new G4LogicalVolume(CapsTap_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"),"PWCapsTaps",0,0,0);
  
        Capsule_visatt = new G4VisAttributes( G4Colour(1, 1, 1, 0.5) );
        Capsule_visatt->SetVisibility(true);
        Capsule_logic->SetVisAttributes( Capsule_visatt );
        CapsTap_logic->SetVisAttributes( Capsule_visatt );
    }

    T.setX( 0.0 );
    T.setY( 0.0 );
    T.setZ( 0.0 );
	new G4PVPlacement(0,T,Capsule_logic,"PwCaps",detlogicWorld,false,-1);

    T.setZ( Shift+(caps_width-eps_caps)/2. );
	new G4PVPlacement(0,T,CapsTap_logic,"PwCapsTaps",detlogicWorld,false,-1);

	// Stage 0
	G4Box *Inner_solid = new G4Box("ShapePW0", InnerX/2., InnerY/2., InnerZ/2.);
	G4LogicalVolume *Inner_logic =
        new G4LogicalVolume(Inner_solid, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matInner),"PWLV:0",0,0,0);
	G4VisAttributes *Inner_visatt = new G4VisAttributes( G4Colour(0.0, 0.0, 1.0) );
	Inner_visatt->SetVisibility(true);
	Inner_logic->SetVisAttributes( Inner_visatt );
    	Inner_logic->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );
    
	T.setX( 0.0 );
    	T.setY( 0.0 );
	T.setZ( Shift + caps_width + InnerZ/2. );	
	new G4PVPlacement(0,T,Inner_logic,"PW:0:",detlogicWorld,false,0);

	// Stage 1
	G4Box *Outer_solid = new G4Box("ShapePW1", OuterX/2., OuterY/2., OuterZ/2.);
	G4LogicalVolume *Outer_logic =
        new G4LogicalVolume(Outer_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matOuter), "PWLV:1",0,0,0);
	G4VisAttributes *Outer_visatt = new G4VisAttributes( G4Colour(1.0, 0.0, 0.0) );
	Outer_visatt->SetVisibility(true);
	Outer_logic->SetVisAttributes( Outer_visatt );
    	Outer_logic->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );
    
	T.setX( 0.0 );
    	T.setY( 0.0 );
	T.setZ( Shift + caps_width + InnerZ + OuterZ/2. );
    new G4PVPlacement(0,T,Outer_logic,"PW:1:",detlogicWorld,false,1);
    
    return theDetector;
}


G4VPhysicalVolume *SToGS::ScintillatorDF::MakeCPPW(G4String detname, G4String opt)
{
    G4VPhysicalVolume *theDetector = 0x0; G4LogicalVolume *detlogicWorld; G4Box *detWorld; G4double caps_width = 0.5*CLHEP::mm;
    G4RotationMatrix R;
	G4ThreeVector T;
    G4Transform3D Tr;
    
    // Should be at the beginning before def of new theDetector
    G4String base_element = "DetectorFactory/Scintillators/ParisPW_";
    base_element += opt;
    
    G4VPhysicalVolume    *a_pw = DetectorFactory::theMainFactory()->Get(base_element.data());
    G4AssemblyVolume *assembly = DetectorFactory::theMainFactory()->GetAssembly(base_element.data());
    
    // use a physical as a container to describe the detector
	detWorld= new G4Box(detname,20.*CLHEP::cm,20.*CLHEP::cm,25.*CLHEP::cm);
	detlogicWorld= new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
	
	detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	
	//  Must place the World Physical volume unrotated at (0,0,0).
	theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    detname,      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number
    
    if ( opt.contains("3") ) {
        caps_width = 1.0*CLHEP::mm; // version 3 have PW with caps of 1 mm
    }
    
    const G4double HalfInnerX = 2.54*CLHEP::cm + 1.001*caps_width, HalfInnerY = 2.54*CLHEP::cm + 1.001*caps_width; //
    
	G4double ShiftX[9];
	G4double ShiftY[9];
	
    // different order ?
    /*
     ShiftX[0] = -2*(HalfInnerX);  ShiftX[1] =  -2*(HalfInnerX); ShiftX[2] = -2*(HalfInnerX);
     ShiftX[3] =  0*(HalfInnerX);  ShiftX[4] =   0*(HalfInnerX); ShiftX[5] =  0*(HalfInnerX);
     ShiftX[6] =  2*(HalfInnerX);  ShiftX[7] =   2*(HalfInnerX); ShiftX[8] =  2*(HalfInnerX);
     
     ShiftY[0] = -2*(HalfInnerY);  ShiftY[1] =   0*(HalfInnerY); ShiftY[2] =  2*(HalfInnerY);
     ShiftY[3] = -2*(HalfInnerY);  ShiftY[4] =   0*(HalfInnerY); ShiftY[5] =  2*(HalfInnerY);
     ShiftY[6] = -2*(HalfInnerY);  ShiftY[7] =   0*(HalfInnerY); ShiftY[8] =  2*(HalfInnerY); */
    
    // 0 is center. then start top, left and turn clockwise .. to be done
    //       y  z
    //       | /
    //       |/
    //   x----
    //
    ShiftX[0] =  0*(HalfInnerX);  ShiftX[1] =  +2*(HalfInnerX); ShiftX[2] =  0*(HalfInnerX);
	ShiftX[3] = -2*(HalfInnerX);  ShiftX[4] =  -2*(HalfInnerX); ShiftX[5] = -2*(HalfInnerX);
	ShiftX[6] =  0*(HalfInnerX);  ShiftX[7] =   2*(HalfInnerX); ShiftX[8] =  2*(HalfInnerX);
	
	ShiftY[0] =  0*(HalfInnerY);  ShiftY[1] =  +2*(HalfInnerY); ShiftY[2] =  2*(HalfInnerY);
	ShiftY[3] =  2*(HalfInnerY);  ShiftY[4] =   0*(HalfInnerY); ShiftY[5] = -2*(HalfInnerY);
	ShiftY[6] = -2*(HalfInnerY);  ShiftY[7] =  -2*(HalfInnerY); ShiftY[8] =  0*(HalfInnerY);
    
    // Add crystals Inner
    for (G4int i = 0; i < 9; i++) {
        T.setX( ShiftX[i] );
        T.setY( ShiftY[i] );
        T.setZ( 0 );
        assembly->MakeImprint( detlogicWorld, T, &R );
    }
    // remap assembly
    DoMap(assembly,a_pw,SToGS::DetectorFactory::GetGCopyNb());
    
    return theDetector;
}

G4VPhysicalVolume *SToGS::ScintillatorDF::MakeFATIMAM(G4String detname, G4String opt)
{
    // Note : housing_width is not used in this version
    
    G4VPhysicalVolume *theDetector = 0x0; G4LogicalVolume *detlogicWorld; G4Box *detWorld; G4bool do_caps = true; // do_housing = false;
    G4RotationMatrix R;
	G4ThreeVector T;
    G4Transform3D Tr;
    
    // use a physical as a container to describe the detector
	detWorld= new G4Box(detname,10.*CLHEP::cm,10.*CLHEP::cm,50.*CLHEP::cm);
	detlogicWorld= new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
	
	detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	//  Must place the World Physical volume unrotated at (0,0,0).
	theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    detname,      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number
    
	// material for first and second layer
	G4String matPMT = "Al", matShielding = "Pb", matWindow = "SToGS_Quartz", matSD = "LaBr3";
    if ( opt.contains("bare") ) {
        matPMT = "AIR";
        matShielding = "AIR";
        matWindow = "AIR";
        
        do_caps = false;
    }
    
    // those numbers have been extracted from the AGATA code / Geometry of FATIMA implemented by M. Labiche, re-
    const double    Labr_Length = 45.75*CLHEP::mm;         //length of the crystal (+5mm = 2inch)
    const double    Labr_Rad    = 19.0*CLHEP::mm;          //radius of the crystal (*2 = 1.5 inch)
    const double    Window_Length = 5.0*CLHEP::mm;         //quartz window length behind LaBr crystal
    
    // to build the Al housing and lead shielding
    G4double zPlane1[16] =
    {0.0*CLHEP::mm,1.0*CLHEP::mm,1.0*CLHEP::mm,1.5*CLHEP::mm,1.5*CLHEP::mm,46.0*CLHEP::mm,46.0*CLHEP::mm,56.0*CLHEP::mm,
        56.0*CLHEP::mm,71.0*CLHEP::mm,71.0*CLHEP::mm,244.3*CLHEP::mm,244.3*CLHEP::mm,260.8*CLHEP::mm,260.8*CLHEP::mm,263.3*CLHEP::mm};
    G4double rInner1[16] =
    {19.0*CLHEP::mm,19.0*CLHEP::mm,0.0*CLHEP::mm,0.0*CLHEP::mm,Labr_Rad,Labr_Rad,Labr_Rad,Labr_Rad,
        32.0*CLHEP::mm,32.0*CLHEP::mm,32.0*CLHEP::mm,32.0*CLHEP::mm,32.0*CLHEP::mm,32.0*CLHEP::mm, 0.0*CLHEP::mm,0.0*CLHEP::mm};
    G4double rOuter1[16] =
    {22.75*CLHEP::mm,22.75*CLHEP::mm,22.75*CLHEP::mm,22.75*CLHEP::mm,22.75*CLHEP::mm,22.75*CLHEP::mm,37.5*CLHEP::mm,37.5*CLHEP::mm,
        37.5*CLHEP::mm,37.5*CLHEP::mm,35.0*CLHEP::mm,35.0*CLHEP::mm,37.5*CLHEP::mm,37.5*CLHEP::mm,37.5*CLHEP::mm,35.5*CLHEP::mm};
    
    G4Polycone *solid_PMT = new G4Polycone("ShapeFATIMAPMT", 0.0, 360.*CLHEP::deg, 16, zPlane1 ,rInner1, rOuter1);
    G4LogicalVolume *logicPMT = new G4LogicalVolume( solid_PMT, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matPMT), "FATIMAPMTLV", 0, 0, 0 );
    
    // to build the Lead shielding
    G4double zPlane2[5] = {0.0*CLHEP::mm  , 2.0*CLHEP::mm,   40.0*CLHEP::mm,  40.0*CLHEP::mm, 46.0*CLHEP::mm};
    G4double rInner2[5] = {22.8*CLHEP::mm , 22.8*CLHEP::mm,  22.8*CLHEP::mm,  22.8*CLHEP::mm, 22.8*CLHEP::mm};
    G4double rOuter2[5] = {25.75*CLHEP::mm, 27.75*CLHEP::mm, 27.75*CLHEP::mm, 37.5*CLHEP::mm, 37.5*CLHEP::mm};
    G4Polycone *solid_Shielding = new G4Polycone("ShapeFATIMAShielding", 0.0, 360.*CLHEP::deg, 5, zPlane2 ,rInner2, rOuter2);
    G4LogicalVolume *logicShielding = new G4LogicalVolume( solid_Shielding, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matShielding), "FATIMAShieldingLV", 0, 0, 0 );
    
    //Building the LaBr3 crystal + window:
    G4Tubs *solid_Scintillator = new G4Tubs("ShapeFATIMAScintillator", 0.0*CLHEP::mm, Labr_Rad, Labr_Length/2 ,0, 360.*CLHEP::deg);
    G4LogicalVolume *logicScintillator = new G4LogicalVolume( solid_Scintillator, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matSD), "FATIMAScintillatorLV", 0, 0, 0 );
    
    G4Tubs *solid_Window = new G4Tubs("ShapeFATIMAWindow", 0.0*CLHEP::mm, Labr_Rad, Window_Length/2 ,0, 360.*CLHEP::deg);
    G4LogicalVolume *logicWindow = new G4LogicalVolume( solid_Window, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matWindow), "FATIMAWindowLV", 0, 0, 0 );
    
    G4VisAttributes *visatt;
    
    T.setX( 0.0 );
    T.setY( 0.0 );
    T.setZ( 0.0 );
    if ( do_caps ) {
        visatt = new G4VisAttributes( G4Colour(0.0, 1.0, 1.0, 0.75) );
        visatt->SetVisibility(true);
        logicPMT->SetVisAttributes( visatt );
    }
    else {
        visatt = new G4VisAttributes( G4Colour(1.0, 1.0, 1.0, 0.5) );
        visatt->SetVisibility(true);
        logicPMT->SetVisAttributes( visatt );
    }
    new G4PVPlacement(0,T,logicPMT,"PMT",detlogicWorld,false,-1);
    
    T.setX( 0.0 );
    T.setY( 0.0 );
    T.setZ( 0.0 );
    if ( do_caps ) {
        visatt = new G4VisAttributes( G4Colour(1.0, 1.0, 0.0, 0.75) );
        visatt->SetVisibility(true);
        logicShielding->SetVisAttributes( visatt );
    }
    else {
        visatt = new G4VisAttributes( G4Colour(1.0, 1.0, 1.0, 0.5) );
        visatt->SetVisibility(true);
        logicShielding->SetVisAttributes( visatt );
    }
    new G4PVPlacement(0,T,logicShielding,"Shield",detlogicWorld,false,-1);
    
    T.setX( 0.0 );
    T.setY( 0.0 );
    T.setZ( 2.*CLHEP::mm + Labr_Length/2. );

    visatt = new G4VisAttributes( G4Colour(0.0, 1.0, 0.0, 1) );
    visatt->SetVisibility(true);
    logicScintillator->SetVisAttributes( visatt );
    logicScintillator->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );
    new G4PVPlacement(0,T,logicScintillator,"Sensor",detlogicWorld,false,0);
    
    T.setX( 0.0 );
    T.setY( 0.0 );
    T.setZ( 2.*CLHEP::mm + Labr_Length + Window_Length/2. );
    if ( do_caps ) {
        visatt = new G4VisAttributes( G4Colour(0.0, 0.0, 1.0, 0.75) );
        visatt->SetVisibility(true);
        logicWindow->SetVisAttributes( visatt );
    }
    else {
        visatt = new G4VisAttributes( G4Colour(1.0, 1.0, 1.0, 0.5) );
        visatt->SetVisibility(true);
        logicWindow->SetVisAttributes( visatt );
    }
    new G4PVPlacement(0,T,logicWindow,"Window",detlogicWorld,false,-1);

    
    return theDetector;
}


G4VPhysicalVolume *SToGS::ScintillatorDF::MakeFATIMAQ(G4String detname, G4String opt)
{

    G4VPhysicalVolume *theDetector = 0x0; G4LogicalVolume *detlogicWorld; G4Box *detWorld;
    G4RotationMatrix R;
	G4ThreeVector T;

    // Should be at the beginning before def of new theDetector
    G4String base_element = "DetectorFactory/Scintillators/FATIMAM";
    if ( opt.size() ) {
        base_element += "_";
        base_element += opt;
    }
    
    G4VPhysicalVolume    *a_it = DetectorFactory::theMainFactory()->Get(base_element.data());
    G4AssemblyVolume *assembly = DetectorFactory::theMainFactory()->GetAssembly(base_element.data());
    
    // use a physical as a container to describe the detector
	detWorld= new G4Box(detname,10.*CLHEP::cm,10.*CLHEP::cm,25.*CLHEP::cm);
	detlogicWorld= new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
	
	detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	
	//  Must place the World Physical volume unrotated at (0,0,0).
	theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    detname,      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number
    
    T.setX( +37.5*CLHEP::mm );
    T.setY( +37.5*CLHEP::mm  );
    T.setZ( 0.0 );
    assembly->MakeImprint( detlogicWorld, T, &R);
    T.setX( +37.5*CLHEP::mm );
    T.setY( -37.5*CLHEP::mm  );
    T.setZ( 0.0 );
    assembly->MakeImprint( detlogicWorld, T, &R);
    T.setX( -37.5*CLHEP::mm );
    T.setY( +37.5*CLHEP::mm  );
    T.setZ( 0.0 );
    assembly->MakeImprint( detlogicWorld, T, &R);
    T.setX( -37.5*CLHEP::mm );
    T.setY( -37.5*CLHEP::mm  );
    T.setZ( 0.0 );
    assembly->MakeImprint( detlogicWorld, T, &R);
    
    // remap assembly
    DoMap(assembly,a_it,SToGS::DetectorFactory::GetGCopyNb());
    
    return theDetector;
}



G4VPhysicalVolume *SToGS::ScintillatorDF::MakeEDEN(G4String detname, G4double caps_width, G4double extra_back, G4String opt)
{
    G4VPhysicalVolume *theDetector = 0x0;
    G4LogicalVolume *detlogicWorld;
    G4Box *detWorld;
    G4ThreeVector T;
    
    // WORLD //
    
    // use a physical as a container to describe the detector
	detWorld= new G4Box(detname,10.*CLHEP::cm,10.*CLHEP::cm,50.*CLHEP::cm);
	detlogicWorld= new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
	
	detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	//  Must place the World Physical volume unrotated at (0,0,0).
	theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    detname,      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number
    
    // EDEN //
	
	// material
	// version "v0"
	const char
	*matScint = "XYLENE",
	*matCapsule = "STAINLESS-STEEL",
	*matBack1 = "GLASS_PLATE", // Other possible options : air, vacuum, wet-air
	*matBack2 = "AIR";	// Other possible options : vacuum, wet-air
	
	// Sizes given in literature and fixed arbitrary quantities
	G4double
	r_cell = 10.*CLHEP::cm, //[Lau93]
	h_cell = 5.*CLHEP::cm, //[Lau93]
	front_width = 0.2*CLHEP::mm, //[Lau93]
	glass_width = 6.*CLHEP::mm, //[Lau93]
	eps = 0.1*CLHEP::mm, // Fixed distance at interfaces between materials
	Shift = 0.01*CLHEP::mm;	 // the detector entrance is not just at 0
	
    // version "v0"
    // Scintillator = cylindre, xylene
    // Capsule = tube, steel
    // Front = cylindre, steel
    // Back1 = cylindre, glass
    // Back2 = cylindre, air
    
	//Sizes used to define the shapes
	G4double
	// scintillator
	r_scint = r_cell,
	h_scint = h_cell,
	// capsule tube
	rint_caps = r_scint + eps,
	rext_caps = r_scint + caps_width,
	h_caps = h_scint + glass_width + 2.*eps + extra_back,
	// capsule front
	r_front = rext_caps,
	h_front = front_width,
	// back
	r_back = r_scint,
	h_back1 = glass_width,
	h_back2 = h_caps - h_scint - h_back1 - 2.*eps;
    
    G4Tubs * Scint_solid = new G4Tubs("ShapeEDENscint",0.*CLHEP::mm,r_scint,h_scint/2.,0.*CLHEP::deg,360.*CLHEP::deg);
    
    G4Tubs * Caps_solid = new G4Tubs("ShapeEDENcaps",rint_caps,rext_caps,h_caps/2.,0.*CLHEP::deg,360.*CLHEP::deg);
    
    G4Tubs * CapsFront_solid = new G4Tubs("ShapeEDENcapsfront",0.*CLHEP::mm,r_front,h_front/2.,0.*CLHEP::deg,360.*CLHEP::deg);
    
    G4Tubs * Back1_solid = new G4Tubs("ShapeEDENback1",0.*CLHEP::mm,r_back,h_back1/2.,0.*CLHEP::deg,360.*CLHEP::deg);
    
    G4Tubs * Back2_solid = new G4Tubs("ShapeEDENback2",0.*CLHEP::mm,r_back,h_back2/2.,0.*CLHEP::deg,360.*CLHEP::deg);
    
    // Logical volumes
    
    G4LogicalVolume *Scint_logic = new G4LogicalVolume(Scint_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matScint),"EDENscint",0,0,0);
    
    G4LogicalVolume *Caps_logic = new G4LogicalVolume(Caps_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matCapsule),"EDENcaps",0,0,0);
    
    G4LogicalVolume *CapsFront_logic = new G4LogicalVolume(CapsFront_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matCapsule),"EDENcapsfront",0,0,0);
    
    G4LogicalVolume *Back1_logic = new G4LogicalVolume(Back1_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matBack1),"EDENback1",0,0,0);
    
    G4LogicalVolume *Back2_logic = new G4LogicalVolume(Back2_solid,SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(matBack2),"EDENback2",0,0,0);
    
    
    // Visual attributes
    
    G4VisAttributes *Scint_visatt = new G4VisAttributes( G4Colour(0.0, 1.0, 0.0) );
    Scint_visatt->SetVisibility(true);
    Scint_logic->SetVisAttributes( Scint_visatt );
    
    G4VisAttributes *Capsule_visatt = new G4VisAttributes( G4Colour(0.8, 0.8, 0.8, 0.75) );
    Capsule_visatt->SetVisibility(true);
    Caps_logic->SetVisAttributes( Capsule_visatt );
    CapsFront_logic->SetVisAttributes( Capsule_visatt );
    
    G4VisAttributes *Back1_visatt = new G4VisAttributes( G4Colour(0.0, 0.0, 1.0, 0.5) );
    Back1_visatt->SetVisibility(true);
    Back1_logic->SetVisAttributes( Back1_visatt );
    
    G4VisAttributes *Back2_visatt = new G4VisAttributes( G4Colour(0.0, 0.0, 1.0, 0.25) );
    Back2_visatt->SetVisibility(true);
    Back2_logic->SetVisAttributes( Back2_visatt );
    
    
    // Placements
    
    T.setX( 0.0 );
    T.setY( 0.0 );
    
    if(opt=="v0")
	{
        T.setZ( Shift + h_front/2. );
        new G4PVPlacement(0,T,CapsFront_logic,"EdenCapsFront",detlogicWorld,false,-1);
        
        T.setZ( Shift + h_front + eps + h_scint/2. );
        new G4PVPlacement(0,T,Scint_logic,"EdenScint",detlogicWorld,false,0);
        
        T.setZ( Shift + h_front + eps + h_scint + eps + h_back1/2.);
        new G4PVPlacement(0,T,Back1_logic,"EdenBack1",detlogicWorld,false,-1);
        
        T.setZ( Shift + h_front + eps + h_scint + eps + h_back1 + eps + h_back2/2.);
        new G4PVPlacement(0,T,Back2_logic,"EdenBack2",detlogicWorld,false,-1);
        
        T.setZ( Shift + h_front + eps + h_caps/2.);
        new G4PVPlacement(0,T,Caps_logic,"EdenCaps",detlogicWorld,false,-1);
	}
    else if(opt=="bare")
	{
        T.setZ( Shift + h_front + eps + h_scint/2. );
        new G4PVPlacement(0,T,Scint_logic,"EdenScint",detlogicWorld,false,0);
	}
    
    // Sensitivity	
	
    //	Scint_logic->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );
	Scint_logic->SetSensitiveDetector( SToGS::UserActionInitialization::GetTrackerSD() );
	
	return theDetector;
}

