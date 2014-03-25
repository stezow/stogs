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
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

// SToGS:: includes
#include "SToGS_ShellDetectorConstruction.hh"
#include "SToGS_MaterialConsultant.hh"
#include "SToGS_DetectorFactory.hh"

using namespace std;

SToGS::TwoShellsDetectorConstruction::aShell::aShell()
{
	Name		= "";
	MatName		= "";
	RMin		= 0.0;
	RMax		= 0.0;
	PhiStart	= 0.0;
	PhiDelta	= 0.0;
	ThetaStart	= 0.0;
	ThetaDelta	= 0.0;
	IsActive	= 0;
}

SToGS::TwoShellsDetectorConstruction::aShell::aShell(const aShell &from)
{
	Name		= from.Name;
	MatName		= from.MatName;
	RMin		= from.RMin;
	RMax		= from.RMax;
	PhiStart	= from.PhiStart;
	PhiDelta	= from.PhiDelta;
	ThetaStart	= from.ThetaStart;
	ThetaDelta	= from.ThetaDelta;
	IsActive	= from.IsActive;
}

void SToGS::TwoShellsDetectorConstruction::aShell::Print(std::ostream &out)
{
	out << Name << " " << MatName << " "
        << RMin/CLHEP::cm << " " << RMax/CLHEP::cm << " "
        << PhiStart<< " " << PhiDelta/CLHEP::deg << " " << ThetaStart/CLHEP::deg << " " << ThetaDelta/CLHEP::deg << " " << IsActive << endl;
}

SToGS::TwoShellsDetectorConstruction::TwoShellsDetectorConstruction():
    solidWorld(0),
    logicWorld(0),
    physiWorld(0)
{
	Inner.Name		= "Shell:0";
	Inner.MatName	= "NaI";
	
	Inner.RMin          = 10.*CLHEP::cm;
	Inner.RMax          = 15.*CLHEP::cm;
	Inner.PhiStart		= 0.*CLHEP::deg;
	Inner.PhiDelta		= 360.*CLHEP::deg;
	Inner.ThetaStart	= 0.*CLHEP::deg;
	Inner.ThetaDelta	= 180.*CLHEP::deg;
	Inner.IsActive		= 1;
	
	Outer.Name		= "Shell:1";
	Outer.MatName	= "BGO";
	
	Outer.RMin          = 25.*CLHEP::cm;
	Outer.RMax          = 40.*CLHEP::cm;
	Outer.PhiStart		= 0.*CLHEP::deg;
	Outer.PhiDelta		= 360.*CLHEP::deg;
	Outer.ThetaStart	= 0.*CLHEP::deg;
	Outer.ThetaDelta	= 180.*CLHEP::deg;
	Outer.IsActive		= 1;
	
	// reads parameters from an ascii file
	ComputeParameters();
}

SToGS::TwoShellsDetectorConstruction::TwoShellsDetectorConstruction(G4String filename):
    solidWorld(0),
    logicWorld(0),
    physiWorld(0)
{
	Inner.Name		= "Shell:0";
	Inner.MatName	= "NaI";
	
	Inner.RMin          = 10.*CLHEP::cm;
	Inner.RMax          = 15.*CLHEP::cm;
	Inner.PhiStart		= 0.0*CLHEP::deg;
	Inner.PhiDelta		= 360.*CLHEP::deg;
	Inner.ThetaStart	= 0.*CLHEP::deg;
	Inner.ThetaDelta	= 180.*CLHEP::deg;
	Inner.IsActive		= 1;
	
	Outer.Name		= "Shell:1";
	Outer.MatName	= "BGO";
	
	Outer.RMin          = 25.*CLHEP::cm;
	Outer.RMax          = 40.*CLHEP::cm;
	Outer.PhiStart		= 0.*CLHEP::deg;
	Outer.PhiDelta		= 360.*CLHEP::deg;
	Outer.ThetaStart	= 0.*CLHEP::deg;
	Outer.ThetaDelta	= 180.*CLHEP::deg;
	Outer.IsActive		= 1;
	
	// reads parameters from an ascii file
	ComputeParameters(filename);
}

SToGS::TwoShellsDetectorConstruction::~TwoShellsDetectorConstruction()
{
	for(unsigned int i = 0; i < otherShells.size() ; i++ ) {
        delete  otherShells[i];
    }
}

void SToGS::TwoShellsDetectorConstruction::ComputeParameters(G4String filename)
{
	G4cout << G4endl << " ------ INFO ------ from SToGS::TwoShellsDetectorConstruction::ComputeParameters with " << filename << G4endl;
	
	// init the MaterialConsultant to check the materials defined
	SToGS::MaterialConsultant *materialFactory = SToGS::MaterialConsultant::theConsultant();
	
	// open the ascii file
	ifstream file; file.open(filename.data());
	if ( file.is_open() == false ) {
		G4cout << " ** WARNING ** cannot open file " << filename << G4endl;
		G4cout << "  ==> Current parameters are used " << G4endl;
		Inner.Print(G4cout);
        Outer.Print(G4cout);
		return;
	}
	
	// read the file and init the Shell structures
	const G4int MAXWIDTH = 300; char aline[MAXWIDTH]; file.getline(aline,MAXWIDTH);
    aShell tmp;
	G4int nb_active = 0, nb_mandatory = 0; G4double rlastmin = 0.0, rlastmax = 0.0;
	while ( file.good() ) {
        
		if ( aline[0] == '#' ) {
            file.getline(aline,MAXWIDTH);
            continue;
        } // this line is a comment
		
		// from the line extract the sphere definition
		char name[30], mat[30];
		G4float tmp1, tmp2, tmp3, tmp4, tmp5, tmp6; G4int tmpi;
		
		G4int format = sscanf(aline,"%s %s %f %f %f %f %f %f %d",name,mat,&tmp1,&tmp2,&tmp3,&tmp4,&tmp5,&tmp6,&tmpi);
		
		tmp.Name = name; tmp.MatName = mat;
        
		tmp.RMin        = tmp1*CLHEP::cm;
		tmp.RMax        = tmp2*CLHEP::cm;
		tmp.PhiStart	= tmp3*CLHEP::deg;
		tmp.PhiDelta	= tmp4*CLHEP::deg;
		tmp.ThetaStart	= tmp5*CLHEP::deg;
		tmp.ThetaDelta	= tmp6*CLHEP::deg;
		tmp.IsActive	= tmpi;
        
		// check if these values are consistents ... if not the shell is not registered
		// geometrical considerations
		G4int check;
		
		check  = 0;
		if ( tmp.RMin > tmp.RMax  )
            check++;
		if (   tmp.PhiStart < 0.0 ||   tmp.PhiStart > 360.0/CLHEP::deg )
            check++;
		if (   tmp.PhiDelta < 0.0 ||   tmp.PhiDelta > 360.0/CLHEP::deg )
            check++;
		if ( tmp.ThetaStart < 0.0 || tmp.ThetaStart > 180.0/CLHEP::deg )
            check++;
		if ( tmp.ThetaDelta < 0.0 || tmp.ThetaDelta > 180.0/CLHEP::deg )
            check++;
		
		if ( check > 0 ) {
			G4cout << " ** WARNING ** the following shell has geometrical problems (removed from the final geometry) ";
            tmp.Print(G4cout);
		}
		// check if the material is known
		if ( materialFactory->GetMaterial(tmp.MatName) == NULL ) {
			G4cout << " ** WARNING ** the material of the following shell is unknown: "; tmp.Print(G4cout);
			G4cout << *(G4Material::GetMaterialTable()) << G4endl;
			check++;
		}
		if ( ( tmp.RMin > rlastmin && tmp.RMin < rlastmax ) || ( tmp.RMax > rlastmin && tmp.RMax < rlastmax ) ) {
			G4cout << " ** WARNING ** the following shell overlaps with a previous one: "; tmp.Print(G4cout);
			check++;
		}
		if ( check > 0 ) {
            file.getline(aline,MAXWIDTH);
            continue;
        }
        
		if ( format == 9 ) { // it has been read with the expected format = no errors
            
			rlastmin = tmp.RMin; rlastmax = tmp.RMax; // keep the radius of the last shell
			
			if ( tmp.Name == "Shell:0" ) {
				if ( tmp.MatName == "Air" )
					tmp.IsActive = 0; // Air cannot be a detector
				
				Inner = tmp; nb_mandatory++; if ( tmp.IsActive ) nb_active++;
			}
			else {
				if ( tmp.Name == "Shell:1" ) {
					if ( tmp.MatName == "Air" )
						tmp.IsActive = 0;
                    
					Outer = tmp; nb_mandatory++; if ( tmp.IsActive ) nb_active++;
				}
				else { // inactive shell
					tmp.IsActive = 0; // just to be sure
                    
					// add it to the collection of inactive shells
					aShell *ptshell = new aShell(); (*ptshell) = tmp; otherShells.push_back(ptshell);
				}
			}
		}
		else { G4cout << " ** WARNING ** the following line does not have the expected format:" << G4endl << aline << G4endl; }
		
		file.getline(aline,MAXWIDTH); // READ NEXT LINE
	}
	file.close();
	
	// now check out if something looks wrong (overlappings, etc ..)
	if ( nb_active == 0 ) { // there are no active shells ..
		G4cout << " ** WARNING **, there is no active shell defined " << G4endl;
	}
	if ( nb_mandatory != 2 ) { //
		G4cout << " ** WARNING **, " << 2 - nb_mandatory << " shell(s) missing " << G4endl;
	}
	// print oout the final geometry loaded
	G4cout << " Inner and outer shell definition: " << G4endl;
	G4cout << "\t"; Inner.Print(G4cout); G4cout << "\t"; Outer.Print(G4cout);
	if ( otherShells.size() ) {
		G4cout << " List of passive materials: " << G4endl;
		for(unsigned int i = 0; i < otherShells.size() ; i++ ) { G4cout << "\t"; otherShells[i]->Print(G4cout); }
	}
	G4cout << " ------ END ------ from SToGS::TwoShellsDetectorConstruction::ComputeParameters " << G4endl << G4endl;
}

G4VPhysicalVolume* SToGS::TwoShellsDetectorConstruction::Construct()
{
    // Clean old geometry, if any
    //
    // G4GeometryManager::GetInstance()->OpenGeometry();
    //  G4PhysicalVolumeStore::GetInstance()->Clean();
    //  G4LogicalVolumeStore::GetInstance()->Clean();
    //  G4SolidStore::GetInstance()->Clean();
    
    // build the world based on the outer shell
	G4double HalfWorldLength = 1.5*Outer.RMax;
	
	// just check out that a passive material has not been added after the outer shell
	for (unsigned int i = 0; i < otherShells.size() ; i++ ) {
		if ( otherShells[i]->RMax > Outer.RMax )
            HalfWorldLength = 1.5*otherShells[i]->RMax;
	}
	
	solidWorld= new G4Box("TwoShells",HalfWorldLength,HalfWorldLength,HalfWorldLength);
	logicWorld= new G4LogicalVolume(solidWorld, SToGS::MaterialConsultant::theConsultant()->GetMaterial("Air"), "TheWorld", 0, 0, 0);
	
	logicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	
	//  Must place the World Physical volume unrotated at (0,0,0).
	physiWorld = new G4PVPlacement(0,         // no rotation
                                   G4ThreeVector(), // at (0,0,0)
                                   logicWorld,      // its logical volume
                                   "TwoShells",      // its name
                                   0,               // its mother  volume
                                   false,           // no boolean operations
                                   -1);              // copy number
	
    // Inner and outer are built
	G4Sphere *asolidShell; G4LogicalVolume *alogicShell; G4VPhysicalVolume *aphysiShell; G4VisAttributes *visatt;
	
	asolidShell = new G4Sphere(Inner.Name, Inner.RMin, Inner.RMax, Inner.PhiStart, Inner.PhiDelta, Inner.ThetaStart, Inner.ThetaDelta);
	logicInner = new G4LogicalVolume(asolidShell,
                                      SToGS::MaterialConsultant::theConsultant()->GetMaterial(Inner.MatName),
                                      Inner.Name,0,0,0);
	visatt = new G4VisAttributes( G4Colour(0.0, 0.0, 1.0) );  visatt->SetVisibility(true);
  	logicInner->SetVisAttributes( visatt );
	aphysiShell = new G4PVPlacement(0,        // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    logicInner,     // its logical volume
                                    Inner.Name,      // its name
                                    logicWorld,      // its mother  volume
                                    false,           // no boolean operations
                                    0);              // copy number
    
    
	asolidShell = new G4Sphere(Outer.Name, Outer.RMin, Outer.RMax, Outer.PhiStart, Outer.PhiDelta, Outer.ThetaStart, Outer.ThetaDelta);
	logicOuter = new G4LogicalVolume(asolidShell,
                                      SToGS::MaterialConsultant::theConsultant()->GetMaterial(Outer.MatName),
                                      Outer.Name,0,0,0);
	visatt = new G4VisAttributes( G4Colour(1.0, 0.0, 0.) ); visatt->SetVisibility(true);
  	logicOuter->SetVisAttributes( visatt );
	aphysiShell = new G4PVPlacement(0,        // no rotation
                                    G4ThreeVector(), // at (0,0,0)
                                    logicOuter,     // its logical volume
                                    Outer.Name,      // its name
                                    logicWorld,      // its mother  volume
                                    false,           // no boolean operations
                                    1);              // copy number
    
    // do the same for other shells
	for (unsigned int i = 0; i < otherShells.size() ; i++ ) {
		aShell *tmp = otherShells[i];
		
		asolidShell = new G4Sphere(tmp->Name, tmp->RMin, tmp->RMax, tmp->PhiStart, tmp->PhiDelta, tmp->ThetaStart, tmp->ThetaDelta);
		alogicShell = new G4LogicalVolume(asolidShell,
                                          SToGS::MaterialConsultant::theConsultant()->GetMaterial(tmp->MatName),
                                          tmp->Name,0,0,0);
		visatt = new G4VisAttributes( G4Colour(0.3, 0.3, 0.3) ); visatt->SetVisibility(true);
  		alogicShell->SetVisAttributes( visatt );
		aphysiShell = new G4PVPlacement(0,        // no rotation
                                        G4ThreeVector(), // at (0,0,0)
                                        alogicShell,     // its logical volume
                                        tmp->Name,      // its name
                                        logicWorld,      // its mother  volume
                                        false,           // no boolean operations
                                        -1);              // copy number
	}
	
    // the world is returned
	return physiWorld;
}

void SToGS::TwoShellsDetectorConstruction::ConstructSDandField()
{
    // sensitive part. Because the sensitive part depends on what we would to extrac
	// for analysis, it asks th OutputManager to give it the sensitive part
    G4VSensitiveDetector *sd = SToGS::DetectorFactory::GetSD("/SToGS/Tracker");
	if ( sd && Inner.IsActive ) {
        //SetSensitiveDetector("Shell:0", sd);
  		logicInner->SetSensitiveDetector( sd );
	}
	// sensitive part. Because the sensitive part depends on what we would to extract
	// for analysis, it asks th OutputManager to give it the sensitive part
	if ( sd && Outer.IsActive ) {
        // SetSensitiveDetector("Shell:1", sd);
  		logicOuter->SetSensitiveDetector( sd );
	}
}





