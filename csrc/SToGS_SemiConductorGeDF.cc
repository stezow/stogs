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
#include "SToGS_SemiConductorGeDF.hh"

#include "SToGS_DetectorFactory.hh"
#include "SToGS_MaterialConsultant.hh"
#include "SToGS_UserActionManager.hh"

using namespace std;

// list of all specific factories
namespace  {
    // all SemiConductors
    SToGS::SemiConductorGeDF theGeFactory("SemiConductors/Ge/");
}

G4VPhysicalVolume * SToGS::SemiConductorGeDF::MakeEXOCLOVER(G4String detname, G4String opt)
{
 
// **************************************************************************
// *                              the WORLD                                 *
// **************************************************************************  

  G4bool do_caps = true;
  // Option
  if ( opt.contains("bare" ) )
    {
      do_caps = false;
    }

 

  G4VPhysicalVolume *theDetector = 0x0; //it means is a pointer

  const G4double world_x = 10.*CLHEP::cm;
  const G4double world_y = 10.*CLHEP::cm;
  const G4double world_z = 10.*CLHEP::cm;

 
 
  // G4bool do_caps = false, do_housing = false;

  // use a physical as a container to describe the detector

  G4Box *detWorld= new G4Box(detname,world_x,world_y,world_z);
  G4LogicalVolume *detlogicWorld= new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
	
  //  detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
  
 G4VisAttributes *detlogicWorldVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0)); //turquoise
  detlogicWorld->SetVisAttributes(detlogicWorldVisAtt);

  //  Must place the World Physical volume unrotated at (0,0,0).
  theDetector = new G4PVPlacement(0,         // no rotation
				  G4ThreeVector(), // at (0,0,0)
				  detlogicWorld,      // its logical volume
				  detname,      // its name
				  0,               // its mother  volume
				  false,           // no boolean operations
				  -1);              // copy number



  //here is where you construct your clover EXOGAM

// **************************************************************************
// *                             CLOVER EXOGAM                              *
// **************************************************************************  
 
  const G4double CrystalLength       = 90.0*CLHEP::mm; // Ge crystal length
  const G4double CrystalHoleDepth    = 15.0*CLHEP::mm; // depth at which starts the hole

  const G4double CrystalOuterRadius  = 30.0*CLHEP::mm; // outer radius for crystal
  const G4double CrystalInnerRadius  =  5.0*CLHEP::mm; // inner radius for hole in crystal
  
  const G4double CrystalEdgeOffset1  = 26.0*CLHEP::mm; // distance of the edge from the center of the crystal
  const G4double CrystalEdgeOffset2  = 28.5*CLHEP::mm; // distance of the edge from the center of the crystal
	
  const G4double CrystalEdgeDepth    = 30.0*CLHEP::mm;  // depth to which the crystal is shaped
  const G4double CrystalEdgeAngle    = 22.5*CLHEP::deg; // bevel angle

  const G4double CapsuleWidth        = 1.5*CLHEP::mm;   // capsule width
  const G4double CapsuleLength       = 110.*CLHEP::mm;   // capsule length
  const G4double CapsuleEdgeDepth    = 3.3*CLHEP::cm;   // same as crystal !
  const G4double CrystalToCapsule    = 3.5*CLHEP::mm;   // to be adjusted ..
  
  const G4double BGOLength           = 120.0*CLHEP::mm;
  const G4double BGOWidth            = 25.0*CLHEP::mm;
  
  const G4double CsILength	     = 20.0*CLHEP::mm;

  const G4double Tolerance           = 0.1*CLHEP::mm; // distance between crystals
  const G4double Space               = 1.0*CLHEP::mm; // distance between Al capsule and BGO 
  	
  // define a coaxial shape that will be modify with SubstractSolid 

  G4VPhysicalVolume *CrystalA_phys = 0x0;
  G4VPhysicalVolume *CrystalB_phys = 0x0;
  G4VPhysicalVolume *CrystalC_phys = 0x0;
  G4VPhysicalVolume *CrystalD_phys = 0x0;
  G4VPhysicalVolume *Capsule_phys = 0x0;
  G4VPhysicalVolume *BGO_phys = 0x0;
  G4VPhysicalVolume *CsIBack_phys = 0x0;
  // the Ge crystal dimensions
  G4int nbZplanes = 4;
  G4double zPlaneGe[4] = { 0.0*CLHEP::mm,
			   CrystalHoleDepth,
			   CrystalHoleDepth + 3.0*CLHEP::mm,
			   CrystalLength};  // depth where is the hole
  G4double rInnerGe[4] = { 0.0*CLHEP::mm,
			   0.0*CLHEP::mm, 
			   CrystalInnerRadius, 
			   CrystalInnerRadius};       // to define the hole in the crystal
  G4double rOuterGe[4] = { CrystalOuterRadius, 
			   CrystalOuterRadius, 
			   CrystalOuterRadius,
			   CrystalOuterRadius};  // to define the external surface


  char sName[40]; // generic for named objects
  sprintf(sName, "Crystal");

  G4Polycone *detCrystal= new G4Polycone(G4String(sName),  //name
					 0.*CLHEP::deg,     //phi Start
					 360.*CLHEP::deg,   //phiTotal
					 nbZplanes, // number of sides
					 zPlaneGe,     // number of Z planes
					 rInnerGe,     // inner radius
					 rOuterGe);    // outer radius
 
  
  // box definition to remove some matter to the crystal
  
  G4double Edge[3];
 
  sprintf(sName, "LongEdge1");
  Edge[0] = (CrystalOuterRadius-CrystalEdgeOffset1);	// x half-width
  Edge[1] = 1.001*CrystalOuterRadius;			// y half-width
  Edge[2] = 1.001*CrystalLength/2.0;			// z half-width
  G4Box *cutEdge1  = new G4Box(G4String(sName),Edge[0],Edge[1],Edge[2]);

  sprintf(sName, "LongEdge2");
  Edge[0] = (CrystalOuterRadius-CrystalEdgeOffset2);	// x half-width
  Edge[1] = 1.001*CrystalOuterRadius;			// y half-width
  Edge[2] = 1.001*CrystalLength/2.0;			// z half-width
  G4Box *cutEdge2  = new G4Box(G4String(sName),Edge[0],Edge[1],Edge[2]);		

  sprintf(sName, "Bevel");
  Edge[0] = 1.001*CrystalOuterRadius; 
  Edge[1] = sin(CrystalEdgeAngle)*(CrystalEdgeDepth);		
  Edge[2] = 1.001*CrystalLength/2.0;
  G4Box *cutBevel = new G4Box(G4String(sName),Edge[0],Edge[1],Edge[2]); 
  
  
  // **************************************************************************
  // *                             SUBSTRACTIONS                              *
  // **************************************************************************  	

  // now remove previously defined box from coax. The box must be placed correctly before
  // since the box definition goes from negative to positive values. 	
		
   G4RotationMatrix rm; //  rm.SetName(G4String("RotationEdge"));	
  sprintf(sName, "coax_cut1_edge");
  G4SubtractionSolid *coax_cut1 
    = new G4SubtractionSolid (G4String(sName),  detCrystal, cutEdge1, &rm, G4ThreeVector(-CrystalOuterRadius,0.0,CrystalLength/2.0));

  sprintf(sName, "coax_cut2_edge"); 
  G4SubtractionSolid *coax_cut2 
    = new G4SubtractionSolid (G4String(sName), coax_cut1, cutEdge2, &rm, G4ThreeVector(CrystalOuterRadius,0.0,CrystalLength/2.0));

  sprintf(sName, "coax_cut3_edge");
  rm.rotateZ(90.0*CLHEP::deg);
  G4SubtractionSolid *coax_cut3 
    = new G4SubtractionSolid (G4String(sName), coax_cut2, cutEdge2, &rm, G4ThreeVector(0.0,CrystalOuterRadius,CrystalLength/2.0));

  sprintf(sName, "coax_cut4_edge"); 
  G4SubtractionSolid *coax_cut4 
    = new G4SubtractionSolid (G4String(sName), coax_cut3, cutEdge1, &rm, G4ThreeVector(0.0,-CrystalOuterRadius,CrystalLength/2.0));
  rm.rotateZ(-90.0*CLHEP::deg);

  sprintf(sName, "coax_cut5_edge");
   rm.rotateX(CrystalEdgeAngle);
  G4SubtractionSolid *coax_cut5 
    = new G4SubtractionSolid (G4String(sName), coax_cut4, cutBevel, &rm, G4ThreeVector(0.,CrystalEdgeOffset2,0.)); 
  // rm.rotateX(-CrystalEdgeAngle);
  //Bevel is already rotated X with CrystalEdgeAngle soyou have to put it in place to continu

  rm.rotateX(-CrystalEdgeAngle);
  sprintf(sName, "coax_cut6_edge");
   rm.rotateZ(90.0*CLHEP::deg);
   rm.rotateX(CrystalEdgeAngle);
   G4SubtractionSolid *coax_cut6
     = new G4SubtractionSolid (G4String(sName), coax_cut5, cutBevel, &rm, G4ThreeVector(CrystalEdgeOffset2,0.,0.)); 
   //rotation back  
   rm.rotateX(-CrystalEdgeAngle);
   rm.rotateZ(-90.0*CLHEP::deg);

   //end substraction		


   //Crystal A
   G4RotationMatrix* Crystal_90deg = new G4RotationMatrix();
   Crystal_90deg -> rotateZ(90*CLHEP::deg);

   sprintf(sName, "ShapedCrystalA");
   G4LogicalVolume *pCrystalA  = new G4LogicalVolume( coax_cut6, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("Ge"), G4String(sName), 0, 0, 0 );
   pCrystalA->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );	
  
   CrystalA_phys = new G4PVPlacement(Crystal_90deg,         // no rotation
				     G4ThreeVector(CrystalEdgeOffset1+Tolerance,-CrystalEdgeOffset1-Tolerance,CrystalToCapsule), // at (0,0,0)
				     pCrystalA,      // its logical volume
				     "CrystalA_P",      // its name
				     detlogicWorld,               // its mother  volume
				     false,           // no boolean operations
				     0);              // copy n
   G4VisAttributes *CrystalA_VisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0)); //blue
   pCrystalA ->SetVisAttributes(CrystalA_VisAtt);
 

  //Crystal B
   G4RotationMatrix* Crystal_180deg = new G4RotationMatrix();
   Crystal_180deg -> rotateZ(180*CLHEP::deg);

   sprintf(sName, "ShapedCrystalB");
   G4LogicalVolume *pCrystalB  = new G4LogicalVolume( coax_cut6, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("Ge"), G4String(sName), 0, 0, 0 );
   pCrystalB->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );	
  
   CrystalB_phys = new G4PVPlacement(Crystal_180deg,         // no rotation
				     G4ThreeVector(-CrystalEdgeOffset1-Tolerance,-CrystalEdgeOffset1-Tolerance,CrystalToCapsule), // at (-26*mm,-26*mm,3.5*mm)
				     pCrystalB,      // its logical volume
				     "CrystalB_P",      // its name
				     detlogicWorld,               // its mother  volume
				     false,           // no boolean operations
				     1);              // copy n
   G4VisAttributes *CrystalB_VisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0)); //green 
   pCrystalB ->SetVisAttributes(CrystalB_VisAtt);
 
   
   //Crystal C
   G4RotationMatrix* Crystal_270deg = new G4RotationMatrix();
   Crystal_270deg -> rotateZ(270*CLHEP::deg);

   sprintf(sName, "ShapedCrystalC");
   G4LogicalVolume *pCrystalC  = new G4LogicalVolume( coax_cut6, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("Ge"), G4String(sName), 0, 0, 0 );
   pCrystalC->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );	
   
   CrystalC_phys = new G4PVPlacement(Crystal_270deg,         // no rotation
				     G4ThreeVector(-CrystalEdgeOffset1-Tolerance,CrystalEdgeOffset1+Tolerance,CrystalToCapsule), // at (0,0,0)
				     pCrystalC,      // its logical volume
				     "CrystalC_P",      // its name
				     detlogicWorld,               // its mother  volume
				     false,           // no boolean operations
				     2);              // copy n
   G4VisAttributes *CrystalC_VisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red  
   pCrystalC ->SetVisAttributes(CrystalC_VisAtt);


   //Crystal D
   sprintf(sName, "ShapedCrystalD");
   G4LogicalVolume *pCrystalD  = new G4LogicalVolume( coax_cut6, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("Ge"), G4String(sName), 0, 0, 0 );
  
   pCrystalD->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );
	
   CrystalD_phys = new G4PVPlacement(0,         // no rotation
				     G4ThreeVector(CrystalEdgeOffset1+Tolerance,CrystalEdgeOffset1+Tolerance,CrystalToCapsule), // at (26*mm,26*mm,3.5*mm)
				     pCrystalD,      // its logical volume
				     "CrystalD_P",      // its name
				     detlogicWorld,               // its mother  volume
				     false,           // no boolean operations
				     3);              // copy n
   G4VisAttributes *CrystalD_VisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0)); //yellow
   pCrystalD ->SetVisAttributes(CrystalD_VisAtt);
  


  
  // **************************************************************************
  // *                               Al CAPSULE                               *
  // **************************************************************************  
       
    
    G4int nbslice = 7;
    const G4double widthface = 45.5*CLHEP::mm;
    G4double zSlice[7] = {  0.0*CLHEP::mm, 
			    CapsuleWidth-0.1*CLHEP::mm, 
			    CapsuleWidth, 
			    CapsuleEdgeDepth,
			    CapsuleLength-CapsuleWidth,
			    CapsuleLength-CapsuleWidth+0.1*CLHEP::mm,
			    CapsuleLength
    };  
    G4double InnRad[7] = {  0.00*CLHEP::mm, 
			    0.00*CLHEP::mm, 
			    widthface-CapsuleWidth,
			    CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule - CapsuleWidth, 
			    CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule - CapsuleWidth,
			    0.0*CLHEP::mm,
			    0.0*CLHEP::mm
    };      
    G4double OutRad[7] = {  widthface-1.5*CLHEP::mm, 
			    widthface, 
			    widthface,
			    CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule,
			    CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule,
			    CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule,
			    CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule
    };
  
  G4RotationMatrix* Cap_45deg = new G4RotationMatrix();
  Cap_45deg -> rotateZ(45*CLHEP::deg);
	  
   if ( do_caps ) 
     {	
  
      G4Polyhedra *caps = new G4Polyhedra(G4String("Capsule"),
					  0.*CLHEP::deg,
					  360.*CLHEP::deg,
					  4, 
					  nbslice,
					  zSlice,
					  InnRad,
					  OutRad);
      sprintf(sName, "Capsule");
      G4LogicalVolume * pCapsule  = new G4LogicalVolume( caps, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("Al"), G4String(sName), 0, 0, 0 );	

   
      Capsule_phys = new G4PVPlacement(Cap_45deg,         // no rotation
				       G4ThreeVector(), // at (0,0,0)
				       pCapsule,      // its logical volume
				       "Capsule_P",      // its name
				       detlogicWorld,               // its mother  volume
				       false,           // no boolean operations
				       -1);              // copy number


      G4VisAttributes *Capsule_VisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.75)); //grey
      pCapsule  ->SetVisAttributes(Capsule_VisAtt);
  }

    // **************************************************************************
    // *                            BGO AntiCompton1                            *
    // **************************************************************************  

    // define a coaxial shape that will be modify with SubstractSolid
    
    G4int numZplane = 3;
    G4double zSides[3] = { 0.0*CLHEP::mm,
			   0.0*CLHEP::mm,
			   BGOLength};  
    G4double rInnerBGO[3] = { CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule + Space,
			   CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule + Space ,
			   CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule + Space }; 
    G4double rOuterBGO[3] = { CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule + Space, 
			   CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule + BGOWidth,
			   CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule + BGOWidth}; 
	

    zSides[1] = BGOWidth / tan(CrystalEdgeAngle);
    
    G4Polyhedra *bgo = new G4Polyhedra(G4String("BGO"),  //pName
				       0.*CLHEP::deg,           //phiStart
				       360.*CLHEP::deg,         //phiTotal
				       4,                //numSide
				       numZplane,        //numZPlanes
				       zSides,           //zPlane[]
				       rInnerBGO,        //rInner[]
				       rOuterBGO);       //rOuter[]

// G4Polyhedra *bgo = new G4Polyhedra(G4String("BGO"), 0.*deg, 360.*deg, 4, numZplane, zSides, rInnerBGO, rOuterBGO);
    sprintf(sName, "BGORear");
    G4LogicalVolume *pBGO  = new G4LogicalVolume( bgo, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("SToGS_BGO"), G4String(sName), 0, 0, 0 );	
    pBGO->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );

    BGO_phys = new G4PVPlacement(Cap_45deg,         // no rotation
				 G4ThreeVector(0.0*CLHEP::mm,0.0*CLHEP::mm,CrystalEdgeDepth+CrystalToCapsule), // at (0,0,0)
				 pBGO,      // its logical volume
				 "BGO_P",      // its name
				 detlogicWorld,               // its mother  volume
				 false,           // no boolean operations
				 4);              // copy number


    G4VisAttributes *BGO_VisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0)); //cyan
    pBGO  ->SetVisAttributes(BGO_VisAtt);
   
    // **************************************************************************
    // *                        CsIBack Anticompton2                            *
    // **************************************************************************  	
  
    G4Tubs *hole= new G4Tubs(G4String("ColdFinger"),
			     0.0*CLHEP::mm,
			     15.0*CLHEP::mm, 
			     CsILength+0.3*CLHEP::mm,
			     0.0*CLHEP::deg,
			     360.*CLHEP::deg);
    
    G4Box  *fullcsi= new G4Box(G4String("FullCsIBack"),
			       CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule,
			       CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule,
			       CsILength);	
		
 	
    G4SubtractionSolid *hole_cut_csi= new G4SubtractionSolid (G4String("cut_csi"), fullcsi, hole, &rm,G4ThreeVector(0.,0.,0.)); 
    
    sprintf(sName, "CsIBack");		
    G4LogicalVolume *pCsIBack= new G4LogicalVolume( hole_cut_csi, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("SToGS_CsI"), G4String(sName), 0, 0, 0 );
    pCsIBack->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );	
    CsIBack_phys = new G4PVPlacement(0,         // no rotation
				     G4ThreeVector(0.,0.,CrystalLength+ 40 *CLHEP::mm+CrystalToCapsule), // at (0,0,0)
				     pCsIBack,      // its logical volume
				     "CsIBack_P",      // its name
				     detlogicWorld,               // its mother  volume
				     false,           // no boolean operations
				     5);              // copy number

    G4VisAttributes *CsIBack_VisAtt= new G4VisAttributes(G4Colour(1.0,0.0,1.0)); //magenta
    pCsIBack  ->SetVisAttributes(CsIBack_VisAtt);
    // CsIBack_VisAtt->SetForceWireframe(true);




    return theDetector;
}




G4VSolid *SToGS::SemiConductorGeDF::AGATAShaper(G4Polycone *polycone,
                                                 G4double *xfront, G4double *yfront, G4double *xback, G4double *yback,
                                                 G4double zback,
                                                 G4double added_dilatation)
{
    G4String name = polycone->GetName(), tmp;
    
    const G4int nb_edges = 6; G4double lxfront[nb_edges], lyfront[nb_edges], lxback[nb_edges], lyback[nb_edges];
    // apply scaling to the different shapes if required
    if ( added_dilatation == 0.0 ) {
        for (G4int i = 0; i < nb_edges; i++) {
            lxfront[i] = xfront[i];
            lyfront[i] = yfront[i];
            lxback[i]  = xback[i];
            lyback[i]  = yback[i];
        }
    }
    else {
        G4double D;
        for (G4int i = 0; i < nb_edges; i++) {
            D = std::sqrt( xfront[i]*xfront[i] + yfront[i]*yfront[i] );
            lxfront[i] = xfront[i] + added_dilatation*xfront[i]/D ;
            lyfront[i] = yfront[i] + added_dilatation*yfront[i]/D ;
            D = std::sqrt( xback[i]*xback[i] + yback[i]*yback[i] );
            lxback[i]  = xback[i] + added_dilatation*xback[i]/D ;
            lyback[i]  = yback[i] + added_dilatation*yback[i]/D ;
        }
    }
    // for each side of the hexagone, computes from the given point to remove an 'infinite' box [edge]
    G4double edge_length_x = 40*CLHEP::mm, edge_width_y = 10*CLHEP::mm, edge_depth_z = 1.2*zback/2.;
    G4int inext;
    
    G4VSolid *result = polycone;
    for (G4int i = 0; i < nb_edges; i++) {
        
        // a new box [edge] to cut the polygone
        tmp  = name; tmp += "_edge_";
        std::stringstream s1;
        s1 << i;
        tmp += s1.str();
        //
        G4Box *edge = new G4Box(tmp,edge_length_x,edge_width_y,edge_depth_z);
        
        // now defines the edges in 3D
        if ( i == nb_edges-1 ) {
            inext = 0;
        }
        else inext = i + 1;
        
        G4ThreeVector T;
        // Compute the 3D vector that goes from middle of front segment to middle of back segment
        // --> this gives the rotatation angle the edge should be rotated alongst X
        G4ThreeVector v_RotX_1((lxback[i]+lxback[inext])/2.-(lxfront[i]+lxfront[inext])/2.,
                               (lyback[i]+lyback[inext])/2.-(lyfront[i]+lyfront[inext])/2.,zback);
        
        // compute 2D vector at the center of the segment
        // --> it is used to determine the additional offset in X,Y
        //     in order to have the face of the box tangent to the surface extracted from the face of the hexagone
        G4ThreeVector v_RotZ_0(lxfront[inext]+lxfront[i],lyfront[inext]+lyfront[i] , 0 );
        // compute 2D vector going from i -> inext.
        // --> it is used to detemine the rotation angle along Z for the box
        G4ThreeVector v_RotZ_1(lxfront[inext]-lxfront[i],lyfront[inext]-lyfront[i] , 0 );
        
        // the total rotation matrix on the edge
        G4RotationMatrix R;
        R.set(0,0,0);
        R.rotateX(v_RotX_1.theta());
        R.rotateZ(v_RotZ_1.phi());
        // the translation to bring the edge at the center of the cylinder
        T.setX( (lxback[inext]+lxfront[i])/2. );
        T.setY( (lyback[inext]+lyfront[i])/2. );
        T.setZ( zback / 2. );
        // the additionnal offset to take into account the wdth of the edge
        G4ThreeVector offset( std::cos(v_RotZ_0.phi())*edge_width_y, std::sin(v_RotZ_0.phi())*edge_width_y , 0 );
        
#ifdef L_DEBUG
        G4cout << "Building Edges " << tmp << G4endl;
        G4cout << " hexogone def @ 0 and " << zback
        << " " << lxfront[i] << " " << lyfront[i] << " " << lxback[i] << " " << lyback[i] << G4endl;
        G4cout << " RotX[theta] " << v_RotX_1.theta()/CLHEP::deg << " RotZ[phi] "  << v_RotZ_1.phi()/CLHEP::deg << G4endl;
        G4cout << "Translation to bring to center  " << T << G4endl;
        G4cout << "Additionnal offset due to edge width " << offset << G4endl;
#endif
        
        T = T + offset;
        tmp  = name; tmp += "_step_";
        std::stringstream s2;
        s2 << i;
        tmp += s2.str();
        result = new G4SubtractionSolid(tmp, result, edge, G4Transform3D(R,T) );
        // used for visualization
        //result = new G4UnionSolid(tmp, result, edge, G4Transform3D((R),T) );
    }
    
    return result;
}
/*
G4LogicalVolume *SToGS::SemiConductorGeDF::MakeAGATACapsule(G4String detname, G4String opt)
{
    G4bool do_caps = true, do_passive = false; G4String tmp;
    
    // cotations
    G4double px[6], py[6], pX[6], pY[6], pz[6], pZ[6]; // definition of the front/back hexagone
    
    G4double HoleR;		// Coaxial hole with radius HoleR
	G4double HoleL; 	// Hole starts at the depth HoleL
	
	G4double CylR;		// Radius of the cylinder that is the base of the shaped crystal
	G4double CylL;		// Length of the cylinder that is the base of the shaped crystal
	G4double CylX;		// Additionnal offset
	G4double CylY;		// Additionnal offset
	G4double CylZ;		// Additionnal offset
	
	G4double ThickB;	// Thickness of the passive area (back)
	G4double ThickC; 	// Thickness of the passive area (coaxial)
	G4double CapS;		// The crystal-encapsulation spacing is capS
	G4double CapT;		// The capsule is capT thick
	G4double Tolerance;	// CapS+CapT+tol
	
	G4double ColX;		// RGB color components
	G4double ColY;		// RGB color components
	G4double ColZ;		// RGB color components
    
    G4double eps = 0.*CLHEP::mm;
    // to avoid overlapping in G4, add/remove atrificially this quantities at borders between caps and Ge
    
    // file to read cotations
    ifstream infil; infil.open("DetectorFactory/SemiConductors/Ge/Builders/agata_capsule.geo");
    if ( !infil.is_open() ) {
        G4cout << "[SToGS] *** Cannot open file " << "DetectorFactory/SemiConductors/Ge/Builders/agata_capsule.geo"
        << endl; // end read the file.
        return 0x0 ;
    }
    
    // Options: bare -> no capsules, passive -> add passive
    G4int which_id = 0;
    if ( detname.contains("Green") ) {
        which_id = 1;
//        innername = "BGreen";
    }
    if ( detname.contains("Blue") ) {
        which_id = 2;
//        innername = "CBlue";
    }
    //    if ( detname.contains("bare") )
    do_caps = false;
    //    if ( detname.contains("passive") )
    //        which_id = 2;
    
    int i1,i2,i3, nb_line = 0, nb_point = 0; double x,y,z,X,Y,Z; std::string line;
    while( infil.good() ) {
        
        getline(infil,line);
        //
        if ( line.size() < 2u )
            continue;
        if ( line[0] == '#' )
            continue;
        
        // decode the line
        if(sscanf(line.data(),"%d %d %d %lf %lf %lf %lf %lf %lf", &i1, &i2, &i3, &x, &y, &z, &X, &Y, &Z) != 9) {
            break;
        }
        
        if(which_id != i1) { // a new crystal is being defined, so init a new ASolid structure
            continue;
        }
        else nb_line++;
        
        if(i2==0 && i3==0) { // basic shape for the crystal
            HoleR = x * CLHEP::mm;
            CylR  = y * CLHEP::mm;
            CylL  = z * CLHEP::mm;
            CylX  = X * CLHEP::mm;
            CylY  = Y * CLHEP::mm;
            CylZ  = Z * CLHEP::mm;
        }
        else if(i2==0 && i3==1) { // passive and capsule
            HoleL  = x * CLHEP::mm;
            ThickB = y * CLHEP::mm;
            ThickC = z * CLHEP::mm;
            CapS   = X * CLHEP::mm;
            CapT   = Y * CLHEP::mm;
            Tolerance   = Z * CLHEP::mm;
        }
        else if(i2==0 && i3==2) { // colors
            ColX     = x;
            ColY     = y;
            ColZ     = z;
        }
        else { // a new point to define the ploygon
            
            px[nb_point] = x * CLHEP::mm;
            py[nb_point] = y * CLHEP::mm;
            pX[nb_point] = X * CLHEP::mm;
            pY[nb_point] = Y * CLHEP::mm;
            
            pz[nb_point] = z * CLHEP::mm;
            pZ[nb_point] = Z * CLHEP::mm;
            
            nb_point++;
        }
    } // while good
    infil.close();
    G4cout << " the file " << "DetectorFactory/SemiConductors/Ge/Builders/agata_capsule.geo"
            << " has been read  " << endl; // end read the file.
    
    // use a physical as a container to describe the detector
    G4RotationMatrix R; G4ThreeVector T; G4Transform3D Tr;
    
    // the coax part :
    // the raw crystal
    G4double *zSliceGe = new G4double[4];
    zSliceGe[0] = (-CylL/2.+eps)*CLHEP::mm;
    zSliceGe[1] = (-CylL/2.+HoleL)*CLHEP::mm;
    zSliceGe[2] = (-CylL/2.+HoleL+0.1)*CLHEP::mm;
    zSliceGe[3] = (+CylL/2.-eps)*CLHEP::mm;
    //
    G4double *InnRadGe = new G4double[4];
    InnRadGe[0] = 0.;
    InnRadGe[1] = 0.;
    InnRadGe[2] = HoleR*CLHEP::mm;
    InnRadGe[3] = HoleR*CLHEP::mm;
    //
    G4double *OutRadGe = new G4double[4];
    OutRadGe[0] = CylR*CLHEP::mm;
    OutRadGe[1] = CylR*CLHEP::mm;
    OutRadGe[2] = CylR*CLHEP::mm;
    OutRadGe[3] = CylR*CLHEP::mm;
    //
    tmp  = detname;
    tmp += "_coax";
    G4Polycone *coax = new G4Polycone(tmp, 0.*deg, 360.*deg, 4, zSliceGe, InnRadGe, OutRadGe ) ;
    // out of the shaper, the crystal with its complex shape
    tmp  = detname;
    tmp += "_crystal";
    G4VSolid *crystal = AGATAShaper(coax,px,py,pX,pY,CylL);
    crystal->SetName(tmp);
    
    // set attributes of the capsule and place it
    G4LogicalVolume *crystal_logic =
    new G4LogicalVolume(crystal,
                        SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("Ge"),tmp,0,0,0);
	G4VisAttributes *crystal_visatt = new G4VisAttributes( G4Colour(ColX, ColY, ColZ,1) );
	crystal_logic->SetVisAttributes( crystal_visatt );
    crystal_logic->SetSensitiveDetector( SToGS::UserActionInitialization::GetCopClusterSD() );
    // the coax part
    
    // encapsulation of the crystal
    G4LogicalVolume *capsule_logic = 0x0, *spacing_logic = 0x0;
    if (do_caps) {
        // spacing filled with air
        zSliceGe[0] = (-CylL/2. - CapT)*CLHEP::mm;
        zSliceGe[1] = (+CylL/2. + CapT)*CLHEP::mm;
        //
        InnRadGe[0] = 0.;
        InnRadGe[1] = 0.;
        //
        OutRadGe[0] = (CylR + CapT)*CLHEP::mm;
        OutRadGe[1] = (CylR + CapT)*CLHEP::mm;
        //
        tmp  = detname;
        tmp += "_coax_spacing";
        G4Polycone *coax_spacing_shape = new G4Polycone(tmp, 0.*deg, 360.*deg, 2, zSliceGe, InnRadGe, OutRadGe);
        //
        tmp  = detname;
        tmp += "_capsule_spacing";
        G4VSolid *spacing = AGATAShaper(coax_spacing_shape,px,py,pX,pY,CylL,CapT);
        spacing->SetName(tmp);
        //
        // set attributes
        spacing_logic = new G4LogicalVolume(spacing,
                                            SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"),tmp,0,0,0);
        //      G4VisAttributes *spacing_visatt = new G4VisAttributes( G4Colour(ColX, ColY, ColZ,0.4) );
        //      spacing_logic->SetVisAttributes( spacing_visatt );
        spacing_logic->SetVisAttributes( G4VisAttributes::Invisible );
        
        // spacing filled with air
        zSliceGe[0] = (-CylL/2. - CapT - CapS)*CLHEP::mm;
        zSliceGe[1] = (+CylL/2. + CapT + CapS)*CLHEP::mm;
        //
        InnRadGe[0] = 0.;
        InnRadGe[1] = 0.;
        //
        OutRadGe[0] = (CylR + CapT + CapS)*CLHEP::mm;
        OutRadGe[1] = (CylR + CapT + CapS)*CLHEP::mm;
        //
        tmp  = detname;
        tmp += "_coax_encapsulation";
        G4Polycone *coax_caps_shape = new G4Polycone(tmp, 0.*deg, 360.*deg, 2, zSliceGe, InnRadGe, OutRadGe);
        //
        tmp  = detname;
        tmp += "_capsule";
        G4VSolid *capsule = AGATAShaper(coax_caps_shape,px,py,pX,pY,CylL,CapT+CapS);
        capsule->SetName(tmp);
        //
        // set attributes
        capsule_logic = new G4LogicalVolume(capsule,
                                            SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("Al"),tmp,0,0,0);
        G4VisAttributes *capsule_visatt = new G4VisAttributes( G4Colour(0.3, 0.3, 0.3,0.8) );
        capsule_logic->SetVisAttributes( capsule_visatt );
    }
    
    if ( capsule_logic && spacing_logic ) {
        return capsule_logic;
    }
    
    return crystal_logic;
}
 */


G4VPhysicalVolume *SToGS::SemiConductorGeDF::MakeAGATACapsule(G4String detname, G4String opt, G4String geo_file)
{
    G4VPhysicalVolume *theDetector = 0x0; G4bool do_caps = true, do_passive = false;
    G4String tmp, crystal_name = "ARed", filename = geo_file; // file with definition
    
    // cotations
    G4double px[6], py[6], pX[6], pY[6], pz[6], pZ[6]; // definition of the front/back hexagone
    
    G4double HoleR;		// Coaxial hole with radius HoleR
	G4double HoleL; 	// Hole starts at the depth HoleL
	
	G4double CylR;		// Radius of the cylinder that is the base of the shaped crystal
	G4double CylL;		// Length of the cylinder that is the base of the shaped crystal
	G4double CylX;		// Additionnal offset
	G4double CylY;		// Additionnal offset
	G4double CylZ;		// Additionnal offset
	
	G4double ThickB;	// Thickness of the passive area (back)
	G4double ThickC; 	// Thickness of the passive area (coaxial)
	G4double CapS;		// The crystal-encapsulation spacing is capS
	G4double CapT;		// The capsule is capT thick
	G4double Tolerance;	// CapS+CapT+tol
	
	G4double ColX;		// RGB color components
	G4double ColY;		// RGB color components
	G4double ColZ;		// RGB color components
    
    G4double eps = 0.*CLHEP::mm;
    // to avoid overlapping in G4, add/remove atrificially this quantities at borders between caps and Ge
    
    // file to read cotations
    ifstream infil; infil.open(filename.data());
    if ( !infil.is_open() ) {
        G4cout << "[SToGS] *** Cannot open file " << filename.data() << endl; // end read the file.
        return 0x0 ;
    }
    
    // Options: bare -> no capsules, passive -> add passive
    G4int which_id = 0;
    if ( detname.contains("Green") ) {
        crystal_name = "BGreen";
        which_id = 1;
    }
    if ( detname.contains("Blue") ) {
        crystal_name = "CBlue";
        which_id = 2;
    }
    if ( opt.contains("bare") )
        do_caps = false;
    //    if ( detname.contains("passive") )
    //        which_id = 2;
    
// The detector
    G4LogicalVolume *detlogicWorld;
    G4Box *detWorld;
	detWorld = new G4Box(detname,10.*CLHEP::cm,10.*CLHEP::cm,25.*CLHEP::cm);
	detlogicWorld=
        new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
	detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	//  Must place the World Physical volume unrotated at (0,0,0).
	theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(0,0,0), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    detname,      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number

    
    int i1,i2,i3, nb_line = 0, nb_point = 0; double x,y,z,X,Y,Z; std::string line;
    while( infil.good() ) {
        
        getline(infil,line);
        //
        if ( line.size() < 2u )
            continue;
        if ( line[0] == '#' )
            continue;
        
        // decode the line
        if(sscanf(line.data(),"%d %d %d %lf %lf %lf %lf %lf %lf", &i1, &i2, &i3, &x, &y, &z, &X, &Y, &Z) != 9) {
            break;
        }
        
        if(which_id != i1) { // a new crystal is being defined, so init a new ASolid structure
            continue;
        }
        else nb_line++;
        
        if(i2==0 && i3==0) { // basic shape for the crystal
            HoleR = x * CLHEP::mm;
            CylR  = y * CLHEP::mm;
            CylL  = z * CLHEP::mm;
            CylX  = X * CLHEP::mm;
            CylY  = Y * CLHEP::mm;
            CylZ  = Z * CLHEP::mm;
        }
        else if(i2==0 && i3==1) { // passive and capsule
            HoleL  = x * CLHEP::mm;
            ThickB = y * CLHEP::mm;
            ThickC = z * CLHEP::mm;
            CapS   = X * CLHEP::mm;
            CapT   = Y * CLHEP::mm;
            Tolerance   = Z * CLHEP::mm;
        }
        else if(i2==0 && i3==2) { // colors
            ColX     = x;
            ColY     = y;
            ColZ     = z;
        }
        else { // a new point to define the ploygon
            
            px[i3] = x * CLHEP::mm;
            py[i3] = y * CLHEP::mm;
            pX[i3] = X * CLHEP::mm;
            pY[i3] = Y * CLHEP::mm;
            
            pz[i3] = z * CLHEP::mm;
            pZ[i3] = Z * CLHEP::mm;
            
            nb_point++;
        }
    } // while good
    infil.close();
    G4cout << " the file " << filename.data() << " has been read  " << endl; // end read the file.
    
    // use a physical as a container to describe the detector
    G4RotationMatrix R; G4ThreeVector T; G4Transform3D Tr;
    
    // the coax part :
    // the raw crystal
    G4double *zSliceGe = new G4double[4];
    zSliceGe[0] = (0)*CLHEP::mm;
    zSliceGe[1] = (HoleL)*CLHEP::mm;
    zSliceGe[2] = (HoleL+0.1)*CLHEP::mm;
    zSliceGe[3] = (+CylL)*CLHEP::mm;
    //
    G4double *InnRadGe = new G4double[4];
    InnRadGe[0] = 0.;
    InnRadGe[1] = 0.;
    InnRadGe[2] = HoleR*CLHEP::mm;
    InnRadGe[3] = HoleR*CLHEP::mm;
    //
    G4double *OutRadGe = new G4double[4];
    OutRadGe[0] = CylR*CLHEP::mm;
    OutRadGe[1] = CylR*CLHEP::mm;
    OutRadGe[2] = CylR*CLHEP::mm;
    OutRadGe[3] = CylR*CLHEP::mm;
    //
    tmp  = crystal_name;
    tmp += "ShapeCoax";
    G4Polycone *coax = new G4Polycone(tmp, 0.*CLHEP::deg, 360.*CLHEP::deg, 4, zSliceGe, InnRadGe, OutRadGe ) ;
    // out of the shaper, the crystal with its complex shape
    tmp  = crystal_name;
    tmp += "ShapeCrystal";
    G4VSolid *crystal = AGATAShaper(coax,px,py,pX,pY,CylL);
    crystal->SetName(tmp);
    // set attributes of the capsule and place it
    tmp  = crystal_name;
    tmp += "LV";
    G4LogicalVolume *crystal_logic =
    new G4LogicalVolume(crystal,
                        SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("Ge"),tmp,0,0,0);
	G4VisAttributes *crystal_visatt = new G4VisAttributes( G4Colour(ColX, ColY, ColZ,1) );
	crystal_logic->SetVisAttributes( crystal_visatt );
    crystal_logic->SetSensitiveDetector( SToGS::UserActionInitialization::GetTrackerSD() );
    
    // encapsulation of the crystal
    G4LogicalVolume *capsule_logic = 0x0, *spacing_logic = 0x0;
    // spacing filled with air
    zSliceGe[0] = (0.0)*CLHEP::mm;
    zSliceGe[1] = (CylL + 2*CapS)*CLHEP::mm;
    //
    InnRadGe[0] = 0.;
    InnRadGe[1] = 0.;
    //
    OutRadGe[0] = (CylR + CapS)*CLHEP::mm;
    OutRadGe[1] = (CylR + CapS)*CLHEP::mm;
    //
    tmp  = crystal_name;
    tmp += "ShapeCoaxSpacing";
    G4Polycone *coax_spacing_shape = new G4Polycone(tmp, 0.*CLHEP::deg, 360.*CLHEP::deg, 2, zSliceGe, InnRadGe, OutRadGe);
    //
    tmp  = crystal_name;
    tmp += "ShapeSpacing";
    G4VSolid *spacing = AGATAShaper(coax_spacing_shape,px,py,pX,pY,CylL,CapS);
    spacing->SetName(tmp);
    //
    // set attributes
    tmp  = crystal_name;
    tmp += "SpacingLV";
    spacing_logic = new G4LogicalVolume(spacing,
                                        SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"),tmp,0,0,0);
    G4VisAttributes *spacing_visatt = new G4VisAttributes( G4Colour(1, 1, 1,0.5) );
    spacing_logic->SetVisAttributes( spacing_visatt );
    
    // capsule filled with Al ... or AIR in case one would like to have only the Ge crystal ... to be used to know lost due to encapsulation
    zSliceGe[0] = (0.0)*CLHEP::mm;
    zSliceGe[1] = (CylL + 2*(CapT + CapS))*CLHEP::mm;
    //
    InnRadGe[0] = 0.;
    InnRadGe[1] = 0.;
    //
    OutRadGe[0] = (CylR + CapT + CapS)*CLHEP::mm;
    OutRadGe[1] = (CylR + CapT + CapS)*CLHEP::mm;
    //
    tmp  = crystal_name;
    tmp += "ShapeCoaxEncapsulation";
    G4Polycone *coax_caps_shape = new G4Polycone(tmp, 0.*CLHEP::deg, 360.*CLHEP::deg, 2, zSliceGe, InnRadGe, OutRadGe);
    //
    tmp  = crystal_name;
    tmp += "ShapeCapsule";
    G4VSolid *capsule = AGATAShaper(coax_caps_shape,px,py,pX,pY,CylL,(CapT+CapS));
    capsule->SetName(tmp);
    //
    // set attributes
    tmp  = crystal_name;
    tmp += "CapsuleLV";
    // in case the capsule is not there, just set it as air
    if (do_caps) {
        capsule_logic = new G4LogicalVolume(capsule,
                                            SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("Al"),tmp,0,0,0);
        G4VisAttributes *capsule_visatt = new G4VisAttributes( G4Colour(0.6, 0.6, 0.6, 0.75) );
        capsule_logic->SetVisAttributes( capsule_visatt );
    }
    else {
        capsule_logic = new G4LogicalVolume(capsule,
                                            SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"),tmp,0,0,0);
        G4VisAttributes *capsule_visatt = new G4VisAttributes( G4Colour(1, 1, 1, 0.5) );
        capsule_logic->SetVisAttributes( capsule_visatt );
    }
    
    // place volumes
    if ( capsule_logic && spacing_logic ) {
        tmp  = crystal_name;
        tmp += "Caps";
        T.setX( 0.0 );
        T.setY( 0.0 );
        T.setZ( 0.0 );
        new G4PVPlacement(0,T,capsule_logic,tmp,detlogicWorld,false,0);
        //
        tmp  = crystal_name;
        tmp += "Spacing";
        T.setZ( CapT );
        new G4PVPlacement(0,T,spacing_logic,tmp,capsule_logic,false,-1);
        //
        T.setZ( CapS + CapT );
        new G4PVPlacement(0,T,crystal_logic,crystal_name,spacing_logic,false,0);
    }

    return theDetector;    
}

G4VPhysicalVolume *SToGS::SemiConductorGeDF::MakeAGATACluster(G4String detname, G4String opt,G4String geo_file)
{
    G4VPhysicalVolume *theDetector = 0x0, *asubdetector = 0x0, *caps; std::vector < G4String > capsules(3);
    
    ifstream infil; infil.open(geo_file.data());
    if ( !infil.is_open() ) {
        G4cout << "[SToGS] *** Cannot open file " << geo_file.data()
        << endl; // end read the file.
        return 0x0 ;
    }
    
    G4LogicalVolume *detlogicWorld;
    G4Box *detWorld;
	detWorld = new G4Box(detname,30.*CLHEP::cm,30.*CLHEP::cm,50.*CLHEP::cm);
	detlogicWorld=
        new G4LogicalVolume(detWorld,
                SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
	detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	//  Must place the World Physical volume unrotated at (0,0,0).
	theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    detname,      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number

    // treatment option ... careful, the real detection volume could be inside a capsule ...
    // one should then look for it before setting copy number ...
    //
    if ( opt.contains("bare") && !opt.contains("bare1") ) {
        capsules[0] = "DetectorFactory/SemiConductors/Ge/AGATA-ARed_bare";
        capsules[1] = "DetectorFactory/SemiConductors/Ge/AGATA-BGreen_bare";
        capsules[2] = "DetectorFactory/SemiConductors/Ge/AGATA-CBlue_bare";

    }
    else { // Ge encapsulation on
        capsules[0] = "DetectorFactory/SemiConductors/Ge/AGATA-ARed";
        capsules[1] = "DetectorFactory/SemiConductors/Ge/AGATA-BGreen";
        capsules[2] = "DetectorFactory/SemiConductors/Ge/AGATA-CBlue";
    }
    // Remap capsuel to give copy number 0 1 2 to Red, Green Blue
    ReMap(Get(capsules[0]),0);
    ReMap(Get(capsules[1]),1);
    ReMap(Get(capsules[2]),2);

    G4int i1,i2,i3, nb_line = 0, which = 0; G4double x,y,z,ps, th, ph; std::string line;
    while( infil.good() ) {
        
        // get euler angles from files
        getline(infil,line);
        //
        if ( line.size() < 2u )
            continue;
        if ( line[0] == '#' )
            continue;
        if(sscanf(line.data(),"%d %d %d %lf %lf %lf %lf %lf %lf", &i1, &i2, &i3, &ps, &th, &ph, &x, &y, &z) != 9) {
            break;
        }
        
        G4ThreeVector T;G4RotationMatrix R; R.set(0,0,0);
        which = i2;
        //
        T.setX(x*CLHEP::mm);
        T.setY(y*CLHEP::mm);
        T.setZ(z*CLHEP::mm);
        //
        R.rotateZ(G4double(ps)*CLHEP::deg);
        R.rotateY(G4double(th)*CLHEP::deg);
        R.rotateZ(G4double(ph)*CLHEP::deg);
        //
        G4Transform3D Tr(R,T);
        Set(capsules[which],theDetector,0,&Tr);
    }
    infil.close();
    
    return theDetector;
}

/*
 
G4VPhysicalVolume * SemiConductorGeDF::MakeEURO_PI(G4String detname, G4String opt)
{
    G4VPhysicalVolume *theDetector = 0x0; G4LogicalVolume *detlogicWorld; G4Box *detWorld; // G4bool do_caps = false, do_housing = false;

    // Option
    if ( opt != "bare" ) {
//            do_caps = true;
    }
    if ( housing_width != 0.0 ) {
        if ( opt != "bare" ) {
            do_housing = true;
        }
    }

    
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
*/
 
G4VPhysicalVolume * SToGS::SemiConductorGeDF::Make(G4String name, G4String version_string)
{
    G4VPhysicalVolume *theDetector = 0x0; G4String detname;
    
    if ( name == "EXOCLOVER" ) {
      detname = GetDetName("EXOCLOVER",version_string);
      theDetector =
	MakeEXOCLOVER(detname,version_string);
    }

    //
    /*
    if ( name == "EXOGAM" ) {
        detname = GetDetName("EXOGAM",version_string);
        theDetector =
            MakeEXO_CLOVER(detname,version_string);
    }
    */
    if ( name.contains("AGATA") ) {
        if ( name.contains("ARed") || name.contains("BGreen") || name.contains("CBlue") ) {
            // AGATA asymmetric Capsules
            detname = GetDetName(name,version_string);
            theDetector =
                MakeAGATACapsule(detname,version_string);
        }
    }
    if ( name == "ATC"  ) { // ATC
        detname = GetDetName(name,version_string);
        theDetector =
        MakeAGATACluster(detname,version_string);
    }
    return theDetector;
}

void SToGS::SemiConductorGeDF::MakeStore()
{
// EUROGAM/BALL related
    
// EXOGAM
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("EXOCLOVER","A-bare"); //only Ge crystal

    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("EXOCLOVER","A-AC");//with Aluminium Capsule
    
// AGATA Capsult and Cluster. Full array in Array factory
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("AGATA-ARed","bare");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("AGATA-BGreen","bare");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("AGATA-CBlue","bare");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("ATC","bare");
    SToGS::DetectorFactory::theMainFactory()->Clean();

    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("AGATA-ARed","");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("AGATA-BGreen","");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("AGATA-CBlue","");
    SToGS::DetectorFactory::theMainFactory()->Clean();
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("ATC","");
}




