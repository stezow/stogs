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
// $Id: SToGS_AGATA.hh,v 1.8 2006/06/29 17:47:30 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
#ifndef SToGS_AGATA_h
#define SToGS_AGATA_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

// std includes
#include <iostream>
#include <sstream>
#include <fstream>

class G4Box;
class G4Sphere;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4Polycone;
class G4AssemblyVolume;


//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! Class to build AGATA from the same configuration files [asolid, acluster, aeuler] that the ones in the AGATA-G4 package
    /*!
     Probably never fully functionnal as it is mainly used for understanding/debugguing.
     
    In the SemiConductor/Ge Factory are defined the basic shapes i.e. the Red, Green and blue capsule as well as Triple Cluster
     
     */
    class AGATA : public G4VUserDetectorConstruction
    {
    private:
        //! path to the data files
        G4String faSolid;
        G4String faCluster;
        G4String faEuler;
        
    private:
        // From the two hexagones, the method shapes a polycone by substracting edges.
        // x/y_front/back are the points defining the hexagone at 0 and z.
        // The dilatation is used to build the spacing/capsule.
        // 1 mm means all points of the hexagones are moved by 1mm along their radius
        // IMPORTANT: the dilatation
        G4VSolid *AGATAShaper(G4Polycone *polycone,
                              G4double *xfront, G4double *yfront, G4double *xback, G4double *yback,
                              G4double zback,
                              G4double added_dilatation = 0.0);
        //! From asolid file provided by AGATA, get a capsule
        G4VPhysicalVolume *PlaceAGATACapsule(G4String which_capsule, G4String opt);
        G4LogicalVolume *MakeAGATACapsule(G4String which_capsule, G4String opt);
        
        //! from acluster file, provided by AGATA, get a triple cluter
        G4AssemblyVolume *MakeAGATACluster(G4String opt);
        
    public:
        //! constructor
        AGATA() :
            faSolid("DetectorFactory/SemiConductors/Ge/Builders/agata_capsule.geo"),
            faCluster("DetectorFactory/SemiConductors/Ge/Builders/agata_cluster.geo"),
            faEuler("DetectorFactory/SemiConductors/Ge/Builders/agata_euler.geo")
        {
            ;
        }
        //! TODO read where are builders from a single configuration file
        AGATA(G4String filename)
        {
            ;
        }
        virtual ~AGATA()
        {
            ;
        }
        //! One of the mandatory class to be implemented in order to have G4 working properly
        virtual G4VPhysicalVolume* Construct();
        
    };
} // SToGS Namespace

#endif
