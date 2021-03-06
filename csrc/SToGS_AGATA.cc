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
#include "SToGS_AGATA.hh"

#include "SToGS_DetectorFactory.hh"
#include "SToGS_MaterialConsultant.hh"
#include "SToGS_UserActionManager.hh"

#define L_DEBUG 1

using namespace std;

G4VSolid *SToGS::AGATA::AGATAShaper(G4Polycone *polycone,
                                                G4double *xfront, G4double *yfront, G4double *xback, G4double *yback,
                                                G4double zback,
                                                G4double added_dilatation)
{
    G4String name = polycone->GetName(), tmp;
    
    const G4int nb_edges = 6;
    G4double lxfront[nb_edges], lyfront[nb_edges], lxback[nb_edges], lyback[nb_edges] /*offsetx, offsety, offsetz */;
    
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
        G4cout << "Additionnal offset dur to edge width " << offset << G4endl;
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

G4LogicalVolume *SToGS::AGATA::MakeAGATACapsule(G4String detname, G4String opt)
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
    ifstream infil; infil.open(faSolid.data());
    if ( !infil.is_open() ) {
        G4cout << "[SToGS] *** Cannot open file " << faSolid.data() << endl; // end read the file.
        return 0x0 ;
    }
    
    // Options: bare -> no capsules, passive -> add passive
    G4int which_id = 0;
    if ( detname.contains("Green") ) {
        which_id = 1;
    }
    if ( detname.contains("Blue") ) {
        which_id = 2;
    }
    if ( opt.contains("bare") )
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
    G4cout << " the file " << faSolid.data() << " has been read  " << endl; // end read the file.
    
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
    tmp  = detname;
    tmp += "_coax";
    G4Polycone *coax = new G4Polycone(tmp, 0.*CLHEP::deg, 360.*CLHEP::deg, 4, zSliceGe, InnRadGe, OutRadGe ) ;
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
        zSliceGe[0] = (0.0)*CLHEP::mm;
        zSliceGe[1] = (CylL + 2*CapT)*CLHEP::mm;
        //
        InnRadGe[0] = 0.;
        InnRadGe[1] = 0.;
        //
        OutRadGe[0] = (CylR + CapT)*CLHEP::mm;
        OutRadGe[1] = (CylR + CapT)*CLHEP::mm;
        //
        tmp  = detname;
        tmp += "_coax_spacing";
        G4Polycone *coax_spacing_shape = new G4Polycone(tmp, 0.*CLHEP::deg, 360.*CLHEP::deg, 2, zSliceGe, InnRadGe, OutRadGe);
        //
        tmp  = detname;
        tmp += "_capsule_spacing";
        G4VSolid *spacing = AGATAShaper(coax_spacing_shape,px,py,pX,pY,CylL,CapT);
        spacing->SetName(tmp);
        //
        // set attributes
        spacing_logic = new G4LogicalVolume(spacing,
                                            SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"),tmp,0,0,0);
        G4VisAttributes *spacing_visatt = new G4VisAttributes( G4Colour(1, 1, 1,0.1) );
        spacing_logic->SetVisAttributes( spacing_visatt );
        
        // capsule filled with Al
        zSliceGe[0] = (0.0)*CLHEP::mm;
        zSliceGe[1] = (2*(CapT + CapS))*CLHEP::mm;
        //
        InnRadGe[0] = 0.;
        InnRadGe[1] = 0.;
        //
        OutRadGe[0] = (CylR + CapT + CapS)*CLHEP::mm;
        OutRadGe[1] = (CylR + CapT + CapS)*CLHEP::mm;
        //
        tmp  = detname;
        tmp += "_coax_encapsulation";
        G4Polycone *coax_caps_shape = new G4Polycone(tmp, 0.*CLHEP::deg, 360.*CLHEP::deg, 2, zSliceGe, InnRadGe, OutRadGe);
        //
        tmp  = detname;
        tmp += "_capsule";
        G4VSolid *capsule = AGATAShaper(coax_caps_shape,px,py,pX,pY,CylL,(CapT+CapS));
        capsule->SetName(tmp);
        //
        // set attributes
        capsule_logic = new G4LogicalVolume(capsule,
                                            SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("Al"),tmp,0,0,0);
        G4VisAttributes *capsule_visatt = new G4VisAttributes( G4Colour(0.3, 0.3, 0.3, 0.6) );
        capsule_logic->SetVisAttributes( capsule_visatt );
    }
    
    if ( capsule_logic && spacing_logic ) {
        return capsule_logic;
    }
    
    return crystal_logic;
}

G4AssemblyVolume *SToGS::AGATA::MakeAGATACluster(G4String opt)
{
    G4AssemblyVolume *theCluster = 0x0; G4VPhysicalVolume *caps; std::vector < G4LogicalVolume * > capsules(3);
    
    ifstream infil; infil.open(faCluster.data());
    if ( !infil.is_open() ) {
        G4cout << "[SToGS] *** Cannot open file " << faCluster.data() << G4endl; // end read the file.
        return 0x0 ;
    }
    //
    capsules[0] = MakeAGATACapsule("AGATA-ARed","bare");
    capsules[1] = MakeAGATACapsule("AGATA-BGreen","bare");
    capsules[2] = MakeAGATACapsule("AGATA-CBlue","bare");
    //
    theCluster = new G4AssemblyVolume();
    
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
        G4String name; // name for each crystal
        if ( which == 0 )
            name = "ARed:0";
        if ( which == 1 )
            name = "BGreen:1";
        if ( which == 2 )
            name = "CBlue:2";
        
        //
        T.setX(x*CLHEP::mm);
        T.setY(y*CLHEP::mm);
        T.setZ(z*CLHEP::mm);
        //
        R.rotateZ(G4double(ps)*CLHEP::deg);
        R.rotateY(G4double(th)*CLHEP::deg);
        R.rotateZ(G4double(ph)*CLHEP::deg);
        //
        G4Transform3D Tr(R,T) ;
        theCluster->AddPlacedVolume( capsules[which] , Tr);
    }
    infil.close();
    G4cout << " the file " << faCluster.data() << " has been read  " << endl; // end read the file.

    return theCluster;
}

#ifdef L_DEBUG
#include "G4GDMLParser.hh"
#endif
G4VPhysicalVolume* SToGS::AGATA::Construct()
{
    G4VPhysicalVolume *theDetector = 0x0;
    
    // try to open the file
    ifstream infil; infil.open(faEuler.data());
    if ( !infil.is_open() )
        return 0x0 ;
    
    G4LogicalVolume *detlogicWorld;
    G4Box *detWorld;
	detWorld = new G4Box("AGATA",2.*CLHEP::m,2.*CLHEP::m,2.*CLHEP::m);
	detlogicWorld=
    new G4LogicalVolume(detWorld,
                        SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), "AGATA", 0, 0, 0);
	detlogicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	//  Must place the World Physical volume unrotated at (0,0,0).
	theDetector = new G4PVPlacement(0,         // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    detlogicWorld,      // its logical volume
                                    "AGATA",      // its name
                                    0,               // its mother  volume
                                    false,           // no boolean operations
                                    -1);              // copy number
    

    G4cout << "\n Reading description of AGATA from file " << faEuler.data() << " ..." << endl;
    
    G4AssemblyVolume *cluster = MakeAGATACluster("");
    
    G4int i1, i2, nclust = 0; G4double ps, th, ph, x, y, z; std::string line;
    while(infil.good())  {
        
        // get euler angles from files
        getline(infil,line);
        //
        if ( line.size() < 2u )
            continue;
        if ( line[0] == '#' )
            continue;
        if(sscanf(line.data(),"%d %d %lf %lf %lf %lf %lf %lf", &i1, &i2, &ps, &th, &ph, &x, &y, &z) != 8)
            break;
        
        G4ThreeVector T(x*CLHEP::mm,y*CLHEP::mm,z*CLHEP::mm); G4RotationMatrix R;
        R.rotateZ(ps*CLHEP::deg);
        R.rotateY(th*CLHEP::deg);
        R.rotateZ(ph*CLHEP::deg);

        G4Transform3D Tr(R,T);
        cluster->MakeImprint(detlogicWorld, Tr );
        nclust++;;
    }
    infil.close();
    
    cout << nclust << " cluster(s) composed this geometry " << endl;
    
#ifdef L_DEBUG
     G4GDMLParser parser;
     parser.Write("toto.gdml",detlogicWorld,false);
#endif
    
    return theDetector;
}




