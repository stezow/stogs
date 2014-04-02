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
// AUTHOR X. Fabian

// G4 includes
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"

#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"

#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "PW_OpticalDetectorConstruction.hh"
#include "SToGS_MaterialConsultant.hh"
#include "SToGS_DetectorFactory.hh"

using namespace std;

PW_OpticalDetectorConstruction::PW_OpticalDetectorConstruction(G4String filename):
    G4VUserDetectorConstruction(),
    logicWorld(0x0),
    logicDetector(0x0,0x0)
{
    if ( filename == "" ) {
        
    }
	else ReadParameters(filename);
}

PW_OpticalDetectorConstruction::~PW_OpticalDetectorConstruction()
{
}

#include "G4SubtractionSolid.hh"

G4VPhysicalVolume* PW_OpticalDetectorConstruction::Construct()
{
	G4Box*             solidWorld;  // pointer to the solid envelope 
	G4VPhysicalVolume* physiWorld;  // pointer to the physical envelope
	
	// build the world based on the outer shell
	G4double HalfWorldLength = 1. * CLHEP::m;
	
	solidWorld= new G4Box("TheWorld",HalfWorldLength,HalfWorldLength,HalfWorldLength);
	logicWorld= new G4LogicalVolume(solidWorld, SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), "TheWorld", 0, 0, 0);
	
	logicWorld->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	
	//  Must place the World Physical volume unrotated at (0,0,0).
	physiWorld = new G4PVPlacement(0,        // no rotation
				G4ThreeVector(), // at (0,0,0)
				logicWorld,      // its logical volume
				"TheWorld",      // its name
				0,               // its mother  volume - xfabian structure en arbre
				false,           // no boolean operations
				-1);              // copy number "detector number" : ID for each Geant4 object
	
	// inch to cm
	cube_side *= 2.54 * CLHEP::cm; cuboid_side *= 2.54 * CLHEP::cm;
	switch(whichGeometry)
	{
		case 0:
			Construct_cuboid(physiWorld);
			break;
		
		case 1:
			Construct_assembledPW(physiWorld);
			break;
		
		default:
			break;
	}
    
#ifdef G4MULTITHREADED
#else
    ConstructSDandField();
#endif
	
	// the world is returned
	return physiWorld;
}

#include "G4PhysicalVolumeStore.hh"
#include "G4SDManager.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

void PW_OpticalDetectorConstruction::Construct_assembledPW(G4VPhysicalVolume* physiWorld)
{
	
	// CONSTRUCT OPTICAL PROPERTIES
	SToGS::MaterialConsultant::theConsultant()->SetOpticalProperties(
                            SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(cube_material), G4String(cube_material));
	SToGS::MaterialConsultant::theConsultant()->SetOpticalProperties(
                            SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(cuboid_material), G4String(cuboid_material));
	
	G4OpticalSurface* cube_world_border = new G4OpticalSurface("cube_world_surface");
	G4OpticalSurface* cuboid_world_border = new G4OpticalSurface("cuboid_world_surface");
	G4OpticalSurface* cube_cuboid_border = new G4OpticalSurface("cube_cuboid_surface");
	
    // give the properties to the surfaces
	PW_SetScintOptSurfaces(cube_world_border, cuboid_world_border, cube_cuboid_border);
	
	// CONSTRUCT THE BLOCKS
	// logicPWPair.first = logical_cube ; logicPWPair.second = logical_cuboid
	std::pair<G4LogicalVolume*, G4LogicalVolume*> logicPWPair = PW_getScintSolids();
	std::pair<G4OpticalSurface*, G4LogicalVolume*> photocathodePair = getCubicPhotocathode();
	
	// Set physical volumes
	G4ThreeVector scint_cube_T(0 * CLHEP::cm, 0 * CLHEP::cm, cubePosition * CLHEP::cm);
	G4ThreeVector scint_cuboid_T(0 * CLHEP::cm, 0 * CLHEP::cm, cuboidPosition * CLHEP::cm);
	G4ThreeVector photocath_T(0 * CLHEP::cm, 0 * CLHEP::cm, cuboid_side / 2. + (cuboidPosition + 0.005) * CLHEP::cm);
	
	G4VPhysicalVolume *scint_cube_phys = new G4PVPlacement(0, scint_cube_T, logicPWPair.first, "scint_cube", logicWorld, false, 1);
	G4VPhysicalVolume *scint_cuboid_phys = new G4PVPlacement(0, scint_cuboid_T, logicPWPair.second, "scint_cuboid", logicWorld, false, 2); // (...,,,, DetectorID)
	G4VPhysicalVolume *photocath_phys = new G4PVPlacement(0, photocath_T, photocathodePair.second, "photocath", logicWorld, false, 3);
	
	// set optical borders (!= skins)
	new G4LogicalBorderSurface("cube_world_border", scint_cube_phys, physiWorld, cube_world_border);
	new G4LogicalBorderSurface("cuboid_world_border", scint_cuboid_phys, physiWorld, cuboid_world_border);
	new G4LogicalBorderSurface("cuboid_photocathode_border", scint_cuboid_phys, photocath_phys, photocathodePair.first);
	
	// set the optical border between the cube and the cuboid
	//new G4LogicalBorderSurface("cube_cuboid_border", scint_cube_phys, scint_cuboid_phys, cube_cuboid_border); // if defined, nothing gets through !
    
    logicDetector = logicPWPair;
}

// set the optical properties of the scintillating material borders, between them and with the world
void PW_OpticalDetectorConstruction::PW_SetScintOptSurfaces(G4OpticalSurface *cube_world_border, G4OpticalSurface *cuboid_world_border, G4OpticalSurface *cube_cuboid_border)
{
	G4MaterialPropertiesTable* cube_world_optPropTable = new G4MaterialPropertiesTable();
	G4MaterialPropertiesTable* cuboid_world_optPropTable = new G4MaterialPropertiesTable();
	G4MaterialPropertiesTable* cube_cuboid_optPropTable = new G4MaterialPropertiesTable();
	
	const G4int NUM = 2; G4double pp[NUM] = {3.0 * CLHEP::eV, 5.0 * CLHEP::eV};
	
	// parameters not explicited here are private attributes of this object
	cube_world_optPropTable->AddProperty("REFLECTIVITY", pp, scint_cube_reflectivity, NUM);
	cuboid_world_optPropTable->AddProperty("REFLECTIVITY", pp, scint_cuboid_reflectivity, NUM);
	cube_cuboid_optPropTable->AddProperty("REFLECTIVITY", pp, cube_cuboid_reflectivity, NUM);
	
	/*cube_world_optPropTable->AddProperty("EFFICIENCY", pp, scint_cube_efficiency, NUM);
	cuboid_world_optPropTable->AddProperty("EFFICIENCY", pp, scint_cuboid_efficiency, NUM);
	cube_cuboid_optPropTable->AddProperty("EFFICIENCY", pp, cube_cuboid_efficiency, NUM);*/
	
	//cube_cuboid_optPropTable->AddProperty("RINDEX", pp, cube_cuboid_rindex, NUM);
	
	cube_world_border->SetMaterialPropertiesTable(cube_world_optPropTable);
	cuboid_world_border->SetMaterialPropertiesTable(cuboid_world_optPropTable);
	cube_cuboid_border->SetMaterialPropertiesTable(cube_cuboid_optPropTable);
	
	WhichSurface(cube_world_border, whichOpticalSurface_cube);
	WhichSurface(cuboid_world_border, whichOpticalSurface_cuboid);
	//WhichSurface(cube_cuboid_border, whichOpticalSurface_cube_cuboid);
	
	return;
}

// Operates on the two pointers given to properly define the two scintillating materials logical volumes
// Returns : pair.first = cube logical volume and pair.second = cuboid logical volume
std::pair<G4LogicalVolume*, G4LogicalVolume*> PW_OpticalDetectorConstruction::PW_getScintSolids()
{
	std::pair<G4LogicalVolume*, G4LogicalVolume*> logicPWpair;
	
	// cf http://www.lcsim.org/software/geant4/doxygen/html/G4Box_8cc_source.html ; half-widths arguments
	G4Box* scint_cube = new G4Box("TestBox", cube_side / 2., cube_side / 2., cube_side / 2.);
	G4Box* scint_cuboid = new G4Box("TestBox", cube_side / 2., cube_side / 2., cuboid_side / 2.);
	
	logicPWpair.first = new G4LogicalVolume(scint_cube,
                                        SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(cube_material), "scint_cube", 0, 0, 0);
	logicPWpair.second = new G4LogicalVolume(scint_cuboid,
                                        SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(cuboid_material), "scint_cuboid", 0, 0, 0);
	
	// openGL visibility attributes
	G4VisAttributes *visatt_cube = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red cube
	G4VisAttributes *visatt_cuboid = new G4VisAttributes(G4Colour(1.0, 0.4, 0.0)); // orange cuboid
	visatt_cube->SetVisibility(true);
	visatt_cuboid->SetVisibility(true);
	logicPWpair.first->SetVisAttributes(visatt_cube);
	logicPWpair.second->SetVisAttributes(visatt_cuboid);
	
	return logicPWpair;
}

// Set the photocathode (the detection material at the rear of the PW)
// Returns : pair.first = the optical surface destined to be a border between the PC and th cuboid
// pair.second = the photocathode logical volume
std::pair<G4OpticalSurface*, G4LogicalVolume*> PW_OpticalDetectorConstruction::getCubicPhotocathode()
{
	std::pair<G4OpticalSurface*, G4LogicalVolume*> photocathodePair;
	
	// Surface properties
	const G4int NUM = 2;
	G4double photocath_EFF[NUM] = {1., 1.}; // Enables 'detection' of photons
	G4double photocath_ReR[NUM] = {1.92, 1.92};
	G4double photocath_ImR[NUM] = {1.69, 1.69};
	G4double Ephoton[NUM] = {3.0 * CLHEP::eV, 5.0 * CLHEP::eV}; // TODO prendre des vraies valeurs
	
	G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
	photocath_mt->AddProperty("EFFICIENCY", Ephoton, photocath_EFF, NUM);
	photocath_mt->AddProperty("REALRINDEX", Ephoton, photocath_ReR, NUM);
	photocath_mt->AddProperty("IMAGINARYRINDEX", Ephoton, photocath_ImR, NUM);
	
	photocathodePair.first = new G4OpticalSurface("photocath_opsurf", glisur, polished, dielectric_metal);
	photocathodePair.first->SetMaterialPropertiesTable(photocath_mt);
	
	// Physical properties
	G4Box* photocath = new G4Box("cubic_photocathode", cube_side / 2., cube_side / 2., 0.01 * CLHEP::cm);
	photocathodePair.second = new G4LogicalVolume(photocath, getPhotocathodeMaterial(), "photocath", 0, 0, 0);
	
	return photocathodePair;
}

void PW_OpticalDetectorConstruction::ReadParameters(G4String filename)
{
	G4cout << "------ INF ------ from PW_OpticalDetectorConstruction::ReadParamters() " << filename << G4endl;
	
	ifstream configFile(filename); string line, tempStr; G4double temp = 0; G4int tempInt = 0;
	if(configFile.is_open())
	{
		while(configFile.good())
		{
			std::getline(configFile, line);
			
			if(line.compare("#") == 0)
				continue;
			
			else if(line.compare("whichGeometry") == 0)
			{
				configFile >> tempInt;
				whichGeometry = tempInt;
				if(whichGeometry == 0)
					G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() whichGeometry set to a unique cube" << G4endl;
				if(whichGeometry == 1)
					G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() whichGeometry set to a cube/cuboid system" << G4endl;
				continue;
			}
			
			else if(line.compare("whichPhotocathMat") == 0)
			{
				configFile >> tempInt;
				whichPhotocathMat = tempInt;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() whichPhotocathMat set to " << whichPhotocathMat << G4endl;
				continue;
			}
			
			else if(line.compare("cube_side") == 0)
			{
				configFile >> temp;
				cube_side = temp;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() cube_side set to " << cube_side << G4endl;
				continue;
			}
			
			else if(line.compare("cuboid_side") == 0)
			{
				configFile >> temp;
				cuboid_side = temp;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() cuboid_side set to " << cuboid_side << G4endl;
				continue;
			}
			
			else if(line.compare("cube_position") == 0 && whichGeometry != 0)
			{
				configFile >> temp;
				cubePosition = temp;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() cube_position set to " << cubePosition << G4endl;
				continue;
			}
			
			else if(line.compare("cuboid_position") == 0)
			{
				configFile >> temp;
				cuboidPosition = temp;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() cuboid_position set to " << cuboidPosition << G4endl;
				continue;
			}
			
			else if(line.compare("scint_cube_reflectivity") == 0 && whichGeometry != 0)
			{
				configFile >> temp;
				scint_cube_reflectivity[0] = temp;
				configFile >> temp;
				scint_cube_reflectivity[1] = temp;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() scint_cube_reflectivity set to " << scint_cube_reflectivity[0] << " and " << scint_cube_reflectivity[1] << G4endl;
				continue;
			}
			
			else if(line.compare("scint_cuboid_reflectivity") == 0)
			{
				configFile >> temp;
				scint_cuboid_reflectivity[0] = temp;
				configFile >> temp;
				scint_cuboid_reflectivity[1] = temp;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() scint_cuboid_reflectivity set to " << scint_cuboid_reflectivity[0] << " and " << scint_cuboid_reflectivity[1] << G4endl;
				continue;
			}
			
			/*else if(line.compare("cube_cuboid_reflectivity") == 0 && whichGeometry != 0)
			{
				configFile >> temp;
				cube_cuboid_reflectivity[0] = temp;
				configFile >> temp;
				cube_cuboid_reflectivity[1] = temp;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() cube_cuboid_reflectivity set to " << cube_cuboid_reflectivity[0] << " and " << cube_cuboid_reflectivity[1] << G4endl;
				continue;
			}*/
			
			else if(line.compare("scint_cube_efficiency") == 0 && whichGeometry != 0)
			{
				configFile >> temp;
				scint_cube_efficiency[0] = temp;
				configFile >> temp;
				scint_cube_efficiency[1] = temp;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() scint_cube_efficiency set to " << scint_cube_efficiency[0] << " and " << scint_cube_efficiency[1] << G4endl;
				continue;
			}
			
			else if(line.compare("scint_cuboid_efficiency") == 0)
			{
				configFile >> temp;
				scint_cuboid_efficiency[0] = temp;
				configFile >> temp;
				scint_cuboid_efficiency[1] = temp;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() scint_cuboid_efficiency set to " << scint_cuboid_efficiency[0] << " and " << scint_cuboid_efficiency[1] << G4endl;
				continue;
			}
			
			/*else if(line.compare("cube_cuboid_efficiency") == 0 && whichGeometry != 0)
			{
				configFile >> temp;
				cube_cuboid_efficiency[0] = temp;
				configFile >> temp;
				cube_cuboid_efficiency[1] = temp;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() cube_cuboid_efficiency set to " << cube_cuboid_efficiency[0] << " and " << cube_cuboid_efficiency[1] << G4endl;
				continue;
			}*/
			
			else if(line.compare("cube_material") == 0 && whichGeometry != 0)
			{
				configFile >> tempStr;
				cube_material = tempStr;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() cube_material set to " << cube_material << G4endl;
				if(!(cube_material.compare("LaBr3") == 0 || cube_material.compare("NaI") == 0) || cuboid_material.compare("CsI") == 0)
					G4cout << " ** WARNING ** in PW_OpticalDetectorConstruction::ReadParameters() ; cube_material should be 'LaBr3', 'NaI' or 'CsI' !" << G4endl;
				continue;
			}
			
			else if(line.compare("cuboid_material") == 0)
			{
				configFile >> tempStr;
				cuboid_material = tempStr;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() cuboid_material set to " << cuboid_material << G4endl;
				if(!(cuboid_material.compare("LaBr3") == 0 || cuboid_material.compare("NaI") == 0 || cuboid_material.compare("CsI") == 0))
					G4cout << " ** WARNING ** in PW_OpticalDetectorConstruction::ReadParameters() ; cuboid_material should be 'LaBr3', 'NaI' or 'CsI' !" << G4endl;
				continue;
			}
			
			else if(line.compare("cube_cuboid_rindex") == 0 && whichGeometry != 0)
			{
				configFile >> temp;
				cube_cuboid_rindex[0] = temp;
				configFile >> temp;
				cube_cuboid_rindex[1] = temp;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() cube_cuboid_rindex set to " << cube_cuboid_rindex[0] << " and " << cube_cuboid_rindex[1] << G4endl;
				continue;
			}
			
			else if(line.compare("whichOpticalSurface_cube") == 0 && whichGeometry != 0)
			{
				configFile >> tempInt;
				whichOpticalSurface_cube = tempInt;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() whichOpticalSurface_cube set to " << whichOpticalSurface_cube << G4endl;
				continue;
			}
			
			else if(line.compare("whichOpticalSurface_cuboid") == 0)
			{
				configFile >> tempInt;
				whichOpticalSurface_cuboid = tempInt;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() whichOpticalSurface_cuboid set to " << whichOpticalSurface_cuboid << G4endl;
				continue;
			}
			
			else if(line.compare("whichOpticalSurface_cube_cuboid") == 0 && whichGeometry != 0)
			{
				configFile >> tempInt;
				whichOpticalSurface_cube_cuboid = tempInt;
				G4cout << " ** INFO ** in PW_OpticalDetectorConstruction::ReadParameters() whichOpticalSurface_cube_cuboid set to " << whichOpticalSurface_cube_cuboid << G4endl;
				continue;
			}
		}
		
		configFile.close();
	}
	
	else
		G4cout << " ** ERROR ** in PW_OpticalDetectorConstruction::set_constructedPW_properties() unable to open " << filename << G4endl;
	
	G4cout << "------ END ------ from PW_OpticalDetectorConstruction::ReadParamters() " << G4endl;
	
	return;
}

void PW_OpticalDetectorConstruction::WhichSurface(G4OpticalSurface* Surface, G4int which)
{
	// modify whichOpticalSurface in the .geo file to choose another surface
	switch(which)
	{
		case 0:
			Surface->SetModel(glisur); // model coming from GEANT3
			
			Surface->SetType(dielectric_metal);
			Surface->SetFinish(polished);
			break;
		case 1:
			Surface->SetModel(unified);
			
			Surface->SetType(dielectric_metal);
			Surface->SetFinish(polished);
			break;
		case 2:
			Surface->SetModel(unified);
			
			Surface->SetType(dielectric_dielectric);
			Surface->SetFinish(polishedteflonair);
			break;
		
		case 3:
			Surface->SetModel(unified);
			
			Surface->SetType(dielectric_dielectric);
			Surface->SetFinish(polished);
			break;
		
		// polished : too mirror ?
		case 4:
			Surface->SetModel(unified);
			Surface->SetType(dielectric_dielectric);
			Surface->SetFinish(etchedteflonair);
			break;
		
		case 5:
			// supposed to be pure lambertian, see http://hypernews.slac.stanford.edu/HyperNews/geant4/get/geometry/574/1.html
			Surface->SetModel(unified);
			Surface->SetType(dielectric_dielectric);
			Surface->SetFinish(groundfrontpainted);
			break;
		
		
		
		default:
			break;
	}
	/*
	 63    polished,                    // smooth perfectly polished surface
	 64    polishedfrontpainted,        // smooth top-layer (front) paint
	 65    polishedbackpainted,         // same is 'polished' but with a back-paint
	 66
	 67    ground,                      // rough surface
	 68    groundfrontpainted,          // rough top-layer (front) paint
	 69    groundbackpainted,           // same as 'ground' but with a back-paint
	 70
	 71    polishedlumirrorair,         // mechanically polished surface, with lumirror
	 72    polishedlumirrorglue,        // mechanically polished surface, with lumirror & meltmount
	 73    polishedair,                 // mechanically polished surface
	 74    polishedteflonair,           // mechanically polished surface, with teflon
	 75    polishedtioair,              // mechanically polished surface, with tio paint
	 76    polishedtyvekair,            // mechanically polished surface, with tyvek
	 77    polishedvm2000air,           // mechanically polished surface, with esr film
	 78    polishedvm2000glue,          // mechanically polished surface, with esr film & meltmount
	 79
	 80    etchedlumirrorair,           // chemically etched surface, with lumirror
	 81    etchedlumirrorglue,          // chemically etched surface, with lumirror & meltmount
	 82    etchedair,                   // chemically etched surface
	 83    etchedteflonair,             // chemically etched surface, with teflon
	 84    etchedtioair,                // chemically etched surface, with tio paint
	 85    etchedtyvekair,              // chemically etched surface, with tyvek
	 86    etchedvm2000air,             // chemically etched surface, with esr film
	 87    etchedvm2000glue,            // chemically etched surface, with esr film & meltmount
	 88
	 89    groundlumirrorair,           // rough-cut surface, with lumirror
	 90    groundlumirrorglue,          // rough-cut surface, with lumirror & meltmount
	 91    groundair,                   // rough-cut surface
	 92    groundteflonair,             // rough-cut surface, with teflon
	 93    groundtioair,                // rough-cut surface, with tio paint
	 94    groundtyvekair,              // rough-cut surface, with tyvek
	 95    groundvm2000air,             // rough-cut surface, with esr film
	 96    groundvm2000glue             // rough-cut surface, with esr film & meltmount
	 */
}

G4Material* PW_OpticalDetectorConstruction::getPhotocathodeMaterial()
{
	G4double a, z, density;	G4int natoms; G4String name; G4Material *pc_mat = 0x0;

	switch(whichPhotocathMat)
	{
		// Aluminium
		case 0:
		{
			density = 2.7 * CLHEP::g / CLHEP::cm3;
			a = 26.98 * CLHEP::g / CLHEP::mole;
			pc_mat = new G4Material(name="pc_aluminium", z=13., a, density);
			break;
		}
		
		//SbKCs (PM tubes book)
		case 1:
		{
			density = 2.7 * CLHEP::g / CLHEP::cm3; // TODO : ?
			
			pc_mat = new G4Material(name="pc_SbKCs", density, 3);
			pc_mat->AddElement(G4Element::GetElement("Sb"), natoms = 1);
			pc_mat->AddElement(G4Element::GetElement("K"), natoms = 1);
			pc_mat->AddElement(G4Element::GetElement("Cs"), natoms = 1);
			break;
		}
		
		default:
			break;
	}
	
	return pc_mat;
}

// Construct a cuboid with a photocathode on the highest z side
void PW_OpticalDetectorConstruction::Construct_cuboid(G4VPhysicalVolume* physiWorld)
{
	// CONSTRUCT OPTICAL PROPERTIES
	SToGS::MaterialConsultant::theConsultant()->SetOpticalProperties(
                                   SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(cuboid_material), G4String(cuboid_material));

	G4OpticalSurface* cuboid_world_border = new G4OpticalSurface("cuboid_world_surface");
	G4MaterialPropertiesTable* cuboid_world_optPropTable = new G4MaterialPropertiesTable();
	
	const G4int NUM = 2;
	G4double pp[NUM] = {3.0 * CLHEP::eV, 5.0 * CLHEP::eV};
	
	// parameters not explicited here are private attributes of this object
	cuboid_world_optPropTable->AddProperty("REFLECTIVITY", pp, scint_cuboid_reflectivity, NUM);
	cuboid_world_border->SetMaterialPropertiesTable(cuboid_world_optPropTable);
	WhichSurface(cuboid_world_border, whichOpticalSurface_cuboid);
	
	// CONSTRUCT THE BLOCKS
	// cf http://www.lcsim.org/software/geant4/doxygen/html/G4Box_8cc_source.html ; half-widths arguments
	G4Box *scint_cuboid = new G4Box("TestBox", cube_side / 2., cube_side / 2., cuboid_side / 2.);
	
	G4LogicalVolume *logic_cuboid = new G4LogicalVolume(scint_cuboid,
                                    SToGS::MaterialConsultant::theConsultant()->FindOrBuildMaterial(cuboid_material), "scint_cuboid", 0, 0, 0);
	
	// openGL visibility attributes
	G4VisAttributes *visatt_cuboid = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red cube
	visatt_cuboid->SetVisibility(true);
	logic_cuboid->SetVisAttributes(visatt_cuboid);
	
	std::pair<G4OpticalSurface*, G4LogicalVolume*> photocathodePair = getCubicPhotocathode();
	
	// Set physical volumes
	G4ThreeVector scint_cuboid_T(0 * CLHEP::cm, 0 * CLHEP::cm, cuboidPosition * CLHEP::cm);
	G4ThreeVector photocath_T(0 * CLHEP::cm, 0 * CLHEP::cm, cuboid_side / 2. + (cuboidPosition + 0.005) * CLHEP::cm);
	
	G4VPhysicalVolume *scint_cuboid_phys = new G4PVPlacement(0, scint_cuboid_T, logic_cuboid, "scint_cuboid", logicWorld, false, 1);
	G4VPhysicalVolume *photocath_phys = new G4PVPlacement(0, photocath_T, photocathodePair.second, "photocath", logicWorld, false, 3);
	
	// set optical borders (!= skins)
	new G4LogicalBorderSurface("cuboid_world_border", scint_cuboid_phys, physiWorld, cuboid_world_border);
	new G4LogicalBorderSurface("cuboid_photocathode_border", scint_cuboid_phys, photocath_phys, photocathodePair.first);
    
    // keep reference at the global level to add sensitivity
    logicDetector.first = logic_cuboid;
}

void PW_OpticalDetectorConstruction::ConstructSDandField()
{
    G4cout <<  " ------ INF ------ from PW_OpticalDetectorConstruction::ConstructSDandField() " << G4endl;

    G4VSensitiveDetector *sd = SToGS::DetectorFactory::GetSD("/SToGS/SD/Tracker");
    if ( logicDetector.first ) {
        logicDetector.first->SetSensitiveDetector( sd );
    }
    if ( logicDetector.second ) {
        logicDetector.second->SetSensitiveDetector( sd );
    }
    
    G4cout << " ------ END ------ from PW_OpticalDetectorConstruction::ConstructSDandField() " << G4endl;

}




