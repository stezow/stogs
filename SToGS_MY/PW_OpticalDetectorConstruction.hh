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
#ifndef PW_OpticalDetectorConstruction_h
#define PW_OpticalDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

// std includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

class G4Box;
class G4Sphere;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4AssemblyVolume;
class G4OpticalSurface;

class PW_OpticalDetectorConstruction : public G4VUserDetectorConstruction
{
private:
    //! pointer to the world logical envelope
	G4LogicalVolume* logicWorld;
    //! pointers to logical. In case of a unique cube, second is null otherwise two stage PW [init in Construct, used in ConstructSDandField]
    std::pair<G4LogicalVolume*, G4LogicalVolume*> logicDetector;
	
private:
    //! parameters for construction
	std::string cube_material, cuboid_material;
	G4int whichOpticalSurface_cube, whichOpticalSurface_cuboid, whichOpticalSurface_cube_cuboid;
	G4int whichGeometry, whichPhotocathMat;
	
	G4double cube_side, cuboid_side;
	G4double cubePosition, cuboidPosition;
	G4double scint_cube_reflectivity[2], scint_cuboid_reflectivity[2], cube_cuboid_reflectivity[2];
	G4double scint_cube_efficiency[2], scint_cuboid_efficiency[2];
//    G4double cube_cuboid_efficiency[2];
	G4double cube_cuboid_rindex[2];
	
private:
	void Construct_assembledPW(G4VPhysicalVolume* physiWorld);
	void PW_SetScintOptSurfaces(G4OpticalSurface *cube_world_border, G4OpticalSurface *cuboid_world_border, G4OpticalSurface *cube_cuboid_border);
	std::pair<G4LogicalVolume*, G4LogicalVolume*> PW_getScintSolids();
	
	void Construct_cuboid(G4VPhysicalVolume* physiWorld);
	
	std::pair<G4OpticalSurface*, G4LogicalVolume*> getCubicPhotocathode();
	G4Material* getPhotocathodeMaterial();
	
	void ReadParameters(G4String filename);
	void WhichSurface(G4OpticalSurface* Surface, G4int which);
	
public:
	PW_OpticalDetectorConstruction(G4String filename = "");
	virtual ~PW_OpticalDetectorConstruction();
	
	//! One of the mandatory class to be implemented in order to have G4 working properly
	virtual G4VPhysicalVolume* Construct();
    
    //! NEW G4.10
    virtual void ConstructSDandField();
};
#endif
