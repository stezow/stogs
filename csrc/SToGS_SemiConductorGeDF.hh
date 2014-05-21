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

#ifndef SToGS_SemiConductorGeDF_h
#define SToGS_SemiConductorGeDF_h 1

#include "SToGS_DetectorFactory.hh"

class G4Polycone;

//! SToGS namespace to protect SToGS classes
namespace SToGS {
//! Base classe to build all Ge like detectors
/*!
 */
class SemiConductorGeDF : public DetectorFactory
{
public:
    //! Should be implemented in any sub factory. It built (C++) a detector and return it
    virtual G4VPhysicalVolume * Make(G4String /* name */, G4String /* version_string */);
    
public:
    // From the two hexagones, the method shapes a polycone by substracting edges.
    // x/y_front/back are the points defining the hexagone at 0 and z.
    // The dilatation is used to build the spacing/capsule.
    // 1 mm means all points of the hexagones are moved by 1mm along their radius
    // IMPORTANT: the dilatation is not applied to the polycone, the user should add it to the radius and length
    G4VSolid *AGATAShaper(G4Polycone *polycone,
                          G4double *xfront, G4double *yfront, G4double *xback, G4double *yback,
                          G4double zback,
                          G4double added_dilatation = 0.0);
    //! From asolid file provided by AGATA, get a capsule
    /*!
        Available option: \n
     - bare --> no aluminum capsule. Default is capsule
     - passive --> add pasive part at back and around central contact
     */
    G4VPhysicalVolume *MakeAGATACapsule(G4String which_capsule, G4String opt,
                                        G4String geo_file = "DetectorFactory/SemiConductors/Ge/Builders/agata_capsule.geo");
    
    //! from acluster file, provided by AGATA, get a triple cluter
    /*!
     Available option: \n
     - bare1 --> no cryostat but Al around Ge crystals
     - bare2 --> no cryostat no Al
     - passive --> add pasive part at back and around central contact
     */
    G4VPhysicalVolume *MakeAGATACluster(G4String detname, G4String opt,
                                        G4String geo_file = "DetectorFactory/SemiConductors/Ge/Builders/agata_cluster.geo");

    //! make a eurogam p1 Ge detector ... to be finished 
    //G4VPhysicalVolume *MakeEURO_PI(G4String detname, G4String opt = "");
    //
  G4VPhysicalVolume *MakeEXOCLOVER(G4String detname, G4String opt = "");
public:
    SemiConductorGeDF(G4String path) : DetectorFactory(path)
        {;}
    virtual ~SemiConductorGeDF()
        {;}
    
    //! build the default store i.e. all the Ge detectors.
    virtual void MakeStore();
};
} // SToGS Namespace

#endif


