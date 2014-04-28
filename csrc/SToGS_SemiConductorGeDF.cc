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

G4VSolid *SToGS::SemiConductorGeDF::AGATAShaper(G4Polycone *polycone,
                                        G4double *xfront, G4double *yfront, G4double *xback, G4double *yback,
                                        G4double base_length,
                                        G4double added_dilatation)
{
    G4String name = polycone->GetName(), tmp;
    const G4int nb_edges = 6;
    G4double length, lxfront[nb_edges], lyfront[nb_edges], lxback[nb_edges], lyback[nb_edges];
    
    G4VSolid *result = polycone;
    
// apply scaling to the different shapes if required
    if ( added_dilatation == 0.0 ) {
        length = base_length;
        for (G4int i = 0; i < nb_edges; i++) {
            lxfront[i] = xfront[i];
            lyfront[i] = yfront[i];
            lxback[i]  = xback[i];
            lyback[i]  = yback[i];
        }
    }
    else {
        length = base_length + 2.*added_dilatation;
        for (G4int i = 0; i < nb_edges; i++) {
            G4double k = (base_length + added_dilatation) / base_length;
            lxfront[i] = k * xfront[i];
            lyfront[i] = k * yfront[i];
            lxback[i]  = k * xback[i];
            lyback[i]  = k * yback[i];
        }
    }
    // for each side of the hexagone, computes from the given point to remove an 'infinite' box [edge]
    G4double edge_width = 8*CLHEP::mm; G4int inext; G4ThreeVector T;
    for (G4int i = 0; i < nb_edges; i++) {
    
        if ( i == nb_edges-1 ) {
            inext = 0;
        }
        else inext = i + 1;
        
        tmp  = name; tmp += "_edge_";
        std::stringstream s;
        s << i;
        tmp += s.str();
        
        G4Box *edge = new G4Box(tmp,edge_width,40*CLHEP::mm,80*CLHEP::mm);
        
        // vectors used to get from the hexagonal points the different rotation angles
        G4ThreeVector v1(lxback[i]-lxfront[i],lyback[i]-lyfront[i],length);
        G4ThreeVector v2( (lxfront[inext]+lxfront[i]) , (lyfront[inext]+lyfront[i]) , 0);
        //
        G4ThreeVector offset(edge_width,0,0);
        //
        G4RotationMatrix rot;
        //
        rot.rotateY(v1.theta());
        rot.rotateZ(v2.phi());
        //
        offset = rot * offset; //

        T.setX( ((lxfront[inext]+lxfront[i])/2.) + offset.getX() );
        T.setY( ((lyfront[inext]+lyfront[i])/2.) + offset.getY() );
        T.setZ( -length / 2. + offset.getZ() );
        // G4cout << v1.theta() / CLHEP::deg << " " << v1.phi() / CLHEP::deg  << G4endl;
        
        result = new G4SubtractionSolid(tmp, result, edge, G4Transform3D(rot,T) );
    }

    return result;
}

G4VPhysicalVolume *SToGS::SemiConductorGeDF::MakeAGATACapsule(G4String detname, G4String opt)
{
    G4VPhysicalVolume *theDetector = 0x0; G4bool do_caps = true, do_passive = false; G4String tmp, innername("ARed");
    
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
    
    G4double eps = 0.0*CLHEP::mm;
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
        innername = "BGreen";
    }
    if ( detname.contains("Blue") ) {
        which_id = 2;
        innername = "CBlue";
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
    //
    G4LogicalVolume *detlogicWorld;
    G4Box *detWorld;
	detWorld = new G4Box(detname,10.*CLHEP::cm,10.*CLHEP::cm,25.*CLHEP::cm);
	detlogicWorld=
        new G4LogicalVolume(detWorld, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), detname, 0, 0, 0);
	detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	//  Must place the World Physical volume unrotated at (0,0,0).
	theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    detname,      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number
    
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
    
    // place volumes 
    if ( capsule_logic && spacing_logic ) {
        T.setX( 0.0 );
        T.setY( 0.0 );
        
        T.setZ( CylL / 2. + CapT + CapS );
        tmp  = detname;
        tmp += "_capsule";
        new G4PVPlacement(0,T,capsule_logic,tmp,detlogicWorld,false,-1);
        
        T.setZ( 0.0 );
        tmp  = detname;
        tmp += "_capsule_spacing";
        new G4PVPlacement(0,T,spacing_logic,tmp,capsule_logic,false,-1);
        tmp  = detname;
        tmp += "_crystal";
        new G4PVPlacement(0,T,crystal_logic,tmp,spacing_logic,false,AddGCopyNb());
    }
    else {
        T.setX( 0.0 );
        T.setY( 0.0 );
        T.setZ( CylL / 2. + CapT + CapS );
       // T.setZ( 0 );
      
        tmp  = detname;
        tmp += "_crystal";
        new G4PVPlacement(0,T,crystal_logic,tmp,detlogicWorld,false,AddGCopyNb());
    }

    return theDetector;
}

G4VPhysicalVolume *SToGS::SemiConductorGeDF::MakeAGATACluster(G4String detname, G4String opt)
{
    G4VPhysicalVolume *theDetector = 0x0, *caps; std::vector < G4LogicalVolume * > capsules(3);
    
    ifstream infil; infil.open("DetectorFactory/SemiConductors/Ge/Builders/agata_cluster.geo");
    if ( !infil.is_open() ) {
        G4cout << "[SToGS] *** Cannot open file " << "DetectorFactory/SemiConductors/Ge/Builders/agata_cluster.geo"
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

    
    // Should be at the beginning before def of new theDetector
    G4String base_element = "DetectorFactory/Scintillators/ParisPW_";
    base_element += opt;
    
    caps = DetectorFactory::theMainFactory()->Get("DetectorFactory/SemiConductors/Ge/AGATA-ARed_bare");
    capsules[0] = caps->GetLogicalVolume()->GetDaughter(0)->GetLogicalVolume();
    caps = DetectorFactory::theMainFactory()->Get("DetectorFactory/SemiConductors/Ge/AGATA-BGreen_bare");
    capsules[1] = caps->GetLogicalVolume()->GetDaughter(0)->GetLogicalVolume();
    caps = DetectorFactory::theMainFactory()->Get("DetectorFactory/SemiConductors/Ge/AGATA-CBlue_bare");
    capsules[2] = caps->GetLogicalVolume()->GetDaughter(0)->GetLogicalVolume();

    G4int i1,i2,i3, nb_line = 0, which = 0; G4double x,y,z,ps, th, ph; std::string line;
    while( infil.good() ) {
        
        getline(infil,line);
        //
        if ( line.size() < 2u )
            continue;
        if ( line[0] == '#' )
            continue;
        if(sscanf(line.data(),"%d %d %d %lf %lf %lf %lf %lf %lf", &i1, &i2, &i3, &ps, &th, &ph, &x, &y, &z) != 9) {
            break;
        }
        
        G4ThreeVector T(x*CLHEP::mm,y*CLHEP::mm,z*CLHEP::mm+100*CLHEP::mm); G4RotationMatrix R /* = new G4RotationMatrix() */;
        R.rotateZ(G4double(ps)*CLHEP::deg);
        R.rotateY(G4double(th)*CLHEP::deg);
        R.rotateZ(G4double(ph)*CLHEP::deg);

        G4Transform3D Tr(R,T);
        G4String tmp  = detname;
        tmp += "_Cluster";
        new G4PVPlacement(Tr,capsules[which],tmp,detlogicWorld,false,-1);
        
        which++;
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
        if ( name == "AGATA-TC"  ) { // AGATA Capsules
            detname = GetDetName(name,version_string);
            theDetector =
                MakeAGATACluster(detname,version_string);
        }
    }
    return theDetector;
}

void SToGS::SemiConductorGeDF::MakeStore()
{
// EUROGAM/BALL related
    
// EXOGAM
    SToGS::DetectorFactory::SetGCopyNb(0);
    //MakeInStore("EXOGAM","A-bare");
    
// AGATA Capsult and Cluster. Full array in Array factory
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("AGATA-ARed","bare");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("AGATA-BGreen","bare");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("AGATA-CBlue","bare");
    SToGS::DetectorFactory::SetGCopyNb(0);
    MakeInStore("AGATA-TC","bare");
}




