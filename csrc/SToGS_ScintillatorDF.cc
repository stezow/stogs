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
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("ParisPW","2");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("ParisPW","2-bare");
    //
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("ParisPW","3");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("ParisPW","3-bare");
    //
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("CParisPW","2");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("CParisPW","2-bare");
    
    // Chateau de Crystal
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("aChateau2Crystal","bare");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("aChateau2Crystal","");
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
	detlogicWorld= new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->GetMaterial("Air"), detname, 0, 0, 0);
	
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
        = new G4Polyhedra("C2CCrys",0.,360.1*CLHEP::deg,numSide,numZPlanes,z,rInner,rOuter);
	//Setting a logical volume
	G4LogicalVolume *a_log
        = new G4LogicalVolume(a_solid,SToGS::MaterialConsultant::theConsultant()->GetMaterial("BaF2"),"C2CCrys",0,0,0);
	
	G4VisAttributes *visatt = new G4VisAttributes( G4Colour(0.0, 0.6, 0.0) );
	visatt->SetVisibility(true);
	a_log->SetVisAttributes( visatt );
    a_log->SetSensitiveDetector( SToGS::UserActionManager::GetCaloSD() );
    
	// Rotation and translation to place each peace in the assembly
	G4ThreeVector Ta;
	
	// Add crystals Inner
	Ta.setX( 0.0 );
    Ta.setY( 0.0 );
    Ta.setZ( 0.0 );
    new G4PVPlacement(0,Ta,a_log,"C2CCrys",detlogicWorld,false,0);
    
    if ( do_caps ) {
        
        // encapsulation side
        a_solid = new G4Polyhedra("C2CCapsSide",0.,360.*CLHEP::deg,numSide,numZPlanes,z_caps,rInner_caps,rOuter_caps);
        //Setting a logical volume
        a_log
            = new G4LogicalVolume(a_solid,SToGS::MaterialConsultant::theConsultant()->GetMaterial("Aluminium"),"C2CCapsSide",0,0,0);
        
        // grey for all passive part
		G4VisAttributes *Capsule_visatt = new G4VisAttributes( G4Colour(0.8, 0.8, 0.8, 0.75) );
        Capsule_visatt->SetVisibility(true);
        a_log->SetVisAttributes( Capsule_visatt );
        
        Ta.setX( 0.0 );
        Ta.setY( 0.0 );
        Ta.setZ( 0.0 );
        new G4PVPlacement(0,Ta,a_log,"C2CCapsSide",detlogicWorld,false,-1);
        
        // encapsulation front and back
        a_solid = new G4Polyhedra("C2CCapsFront",0.,360.*CLHEP::deg,numSide,numZPlanes,z_caps_front,rInner_caps_front,rOuter_caps_front);
        //Setting a logical volume
        a_log
            = new G4LogicalVolume(a_solid,SToGS::MaterialConsultant::theConsultant()->GetMaterial("Aluminium"),"C2CCapsFront",0,0,0);
        
        // grey for all passive part
        a_log->SetVisAttributes( Capsule_visatt );
        
        Ta.setX( 0.0 );
        Ta.setY( 0.0 );
        Ta.setZ( 0.0 );
        new G4PVPlacement(0,Ta,a_log,"C2CCapsFront",detlogicWorld,false,-1);
        Ta.setX( 0.0 );
        Ta.setY( 0.0 );
        Ta.setZ( crystal_depth + caps_width + teflon_width );
        new G4PVPlacement(0,Ta,a_log,"C2CCapsBack",detlogicWorld,false,-1);
    }
	
    return theDetector;
}

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
	detlogicWorld= new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->GetMaterial("Air"), detname, 0, 0, 0);
	
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
	const char *matInner = "LaBr3", *matOuter = "NaI", *matCapsule = "Aluminium", *matHousing = "Aluminium";
	
	// full lenght in X,Y,Z
	G4double
    InnerX = 2*2.54*CLHEP::cm, InnerY = 2*2.54*CLHEP::cm, InnerZ = 2*2.54*CLHEP::cm,
    OuterX = 2*2.54*CLHEP::cm, OuterY = 2*2.54*CLHEP::cm, OuterZ = 6*2.54*CLHEP::cm;
    // in case of housing longer at the back of the PW
	G4double ExtraHousingBack = 20*CLHEP::mm, eps_caps =  caps_width / 100. , eps_housing = 0.05*CLHEP::mm; // eps remove in the capsule part to avoid overlapping
    // G4double r_caps_width = caps_width - eps_caps, r_housing_width = housing_width - eps_housing;
	
    
	// Stage 0
	G4Box *Inner_solid = new G4Box(matInner, InnerX/2., InnerY/2., InnerZ/2.);
	G4LogicalVolume *Inner_logic =
        new G4LogicalVolume(Inner_solid, SToGS::MaterialConsultant::theConsultant()->GetMaterial(matInner),"PW:0",0,0,0);
	G4VisAttributes *Inner_visatt = new G4VisAttributes( G4Colour(0.0, 0.0, 1.0) );
	Inner_visatt->SetVisibility(true);
	Inner_logic->SetVisAttributes( Inner_visatt );
    Inner_logic->SetSensitiveDetector( SToGS::UserActionManager::GetCaloSD() );
    
	T.setX( 0.0 );
    T.setY( 0.0 );
    T.setZ( housing_width + caps_width + InnerZ / 2. );
    
    new G4PVPlacement(0,T,Inner_logic,"PW:0",detlogicWorld,false,0);
	
	// Stage 1
	G4Box *Outer_solid = new G4Box(matOuter, OuterX/2., OuterY/2., OuterZ/2.);
	G4LogicalVolume *Outer_logic =
        new G4LogicalVolume(Outer_solid,SToGS::MaterialConsultant::theConsultant()->GetMaterial(matOuter), "PW:1",0,0,0);
	G4VisAttributes *Outer_visatt = new G4VisAttributes( G4Colour(1.0, 0.0, 0.0) );
	Outer_visatt->SetVisibility(true);
	Outer_logic->SetVisAttributes( Outer_visatt );
    Outer_logic->SetSensitiveDetector( SToGS::UserActionManager::GetCaloSD() );
    
	T.setX( 0.0 );
    T.setY( 0.0 );
    T.setZ( housing_width + caps_width + InnerZ + OuterZ / 2. );
    new G4PVPlacement(0,T,Outer_logic,"PW:1",detlogicWorld,false,1);
    
	if ( do_caps ) {
        
		// grey for all passive part
		G4VisAttributes *Capsule_visatt = new G4VisAttributes( G4Colour(0.8, 0.8, 0.8, 0.75) );
        Capsule_visatt->SetVisibility(true);
        
        G4double zplane[3], rinner[3], router[3];
        
        zplane[0] = 0.01*CLHEP::mm; zplane[1] = caps_width - eps_caps; zplane[2] = InnerZ + OuterZ + ExtraHousingBack;
        rinner[0] = 0.00*CLHEP::mm; rinner[1] = InnerX / 2 + eps_caps; rinner[2] = InnerX / 2 + eps_caps;
        router[0] = InnerX / 2 + caps_width ; router[1] = InnerX / 2 + caps_width ; router[2] = InnerX / 2 + caps_width;
        
        //rinner[0] = rinner[1] = rinner[2] = InnerX / 2 + eps_caps;
        //router[0] = router[1] = router[2] = InnerX / 2 + caps_width;
        
 		G4Polyhedra *side_part =
        new G4Polyhedra("PwCaps", 0.0, 361.0*CLHEP::deg, 4, 3, zplane, rinner, router);
        T.setX( 0.0 );
        T.setY( 0.0 );
        T.setZ( 0.0 );
        
        G4RotationMatrix *Ro = new G4RotationMatrix();
        Ro->rotateZ(45.0*CLHEP::deg);
        
 		G4LogicalVolume *Capsule_logic_side =
        new G4LogicalVolume(side_part,SToGS::MaterialConsultant::theConsultant()->GetMaterial(matCapsule),"PWCaps",0,0,0);
        
		Capsule_logic_side->SetVisAttributes( Capsule_visatt );
        new G4PVPlacement(Ro,T,Capsule_logic_side,"PwCaps",detlogicWorld,false,-1);
        
        /* OLD ... better to avoid boolean operation if possible
         // capsule
         G4Box *capsule_full =
         new G4Box("CAPSFULL", InnerX/2. + caps_width, InnerY/2. + caps_width, InnerZ/2. + OuterZ/2. + caps_width); //2*caps_width/2.
         // to remove the central part
         G4Box *tot_crys =
         new G4Box("TOTCRYS", InnerX/2. + Gap1 , InnerY/2. + Gap1, InnerZ/2. + OuterZ/2. + 2*caps_width ); // 2 so that the endcaps is open behind to let the possibility to add glass window. careful, ta has to be 3
         // needs to move a little bit in z direction to have the face in front of the source filled while back is empty
         // could be check by setting SD on Capsule_logic : DONE
         T.setX( 0.0 );
         T.setY( 0.0 );
         T.setZ( 3*caps_width );
         //
         Tr = G4Transform3D(R,T);
         //
         G4SubtractionSolid *Capsule_solid =
         new G4SubtractionSolid("PWCapsule", capsule_full, tot_crys, Tr);
         G4LogicalVolume *Capsule_logic =
         new G4LogicalVolume(Capsule_solid,SToGS::MaterialConsultant::theConsultant()->GetMaterial(matCapsule),"PWCapsule",0,0,0);
         
         // grey for all passive part
         G4VisAttributes *Capsule_visatt = new G4VisAttributes( G4Colour(0.8, 0.8, 0.8) );
         Capsule_visatt->SetVisibility(true);
         Capsule_logic->SetVisAttributes( Capsule_visatt );
         
         // now add it to the assembly
         T.setX( 0.0 );
         T.setY( 0.0 );
         T.setZ( housing_width + (InnerZ/2. + OuterZ/2. + caps_width)  );
         
         new G4PVPlacement(0,T,Capsule_logic,"PwCapsule",detlogicWorld,false,-1);
         */
	}
	
    // TO BE DONE !
	if ( do_housing ) {
		// housing
		G4Box *Housing_full =
            new G4Box("HOUSINGFULL", InnerX/2. + caps_width + housing_width, InnerY/2. + caps_width + housing_width, InnerZ/2. + OuterZ/2. + ExtraHousingBack/2. + caps_width + housing_width);
		// to remove the centrale part
		G4Box *tot_crys =
        new G4Box("HOUSINGCRYS", InnerX/2. + caps_width, InnerY/2. + caps_width, InnerZ/2. + OuterZ/2. + ExtraHousingBack/2. + caps_width + 2*housing_width ); // 2 so that the endcaps is open behind to let the possibility to add glass window. careful, ta has to be 3
		// needs to move a little bit in z direction to have the face in front of the source filled while back is empty
		// could be check by setting SD on Housing_logic : DONE
		T.setX( 0.0 );
        T.setY( 0.0 );
        T.setZ( 3*housing_width );
		//
		Tr = G4Transform3D(R,T);
		//
		G4SubtractionSolid *Housing_solid =
            new G4SubtractionSolid("PWHousing", Housing_full, tot_crys, Tr);
		G4LogicalVolume *Housing_logic =
            new G4LogicalVolume(Housing_solid,SToGS::MaterialConsultant::theConsultant()->GetMaterial(matHousing),"PWHousing",0,0,0);
        
		// grey for all passive part
		G4VisAttributes *Housing_visatt = new G4VisAttributes( G4Colour(0.8, 0.8, 0.8, 0.75) );
		Housing_visatt->SetVisibility(true);
		Housing_logic->SetVisAttributes( Housing_visatt );
		
		// now add it to the assembly
		T.setX( 0.0 );
        T.setY( 0.0 );
        T.setZ( InnerZ/2. + OuterZ/2. + ExtraHousingBack/2. + caps_width + housing_width);
        
        new G4PVPlacement(0,T,Housing_logic,"PWHousing",detlogicWorld,false,-1);
	}
    
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
	detlogicWorld= new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->GetMaterial("Air"), detname, 0, 0, 0);
	
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




