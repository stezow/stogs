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
// $Id: SToGS_TwoShellsDetectorConstruction.hh,v 1.8 2006/06/29 17:47:30 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
#ifndef SToGS_TwoShellsDetectorConstruction_h
#define SToGS_TwoShellsDetectorConstruction_h 1

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

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! Detector composed of a set of concentric shells with two (possibly one) active shells
    /*!
     This is the simpliest detector composed of shells. The detector is initialised from an unique ascii file (\c setup/shells.geo)
     which is expected in the setup directory. An example of such a file is available in the package (shells.geo.demo). Here is a snapshot:
     \code
     #
     # Ascii file that described a collection of perfect shells for HermeTwoShellsDetectorConstruction
     # You must write from the smallest radius to the largest !!
     # Two shells are mandatory one named Inner and a second one named Outer.
     # In case you want to work with only one active shell the second is composed of AIR and set inative
     # (see Example 1)
     #
     # Full description of one shell per line with the following conventions:
     #   (lengths are in centimeter and angles in degres)
     #   (is_active is equal to 1 if it is an active shell 0 otherwise. active means a detector-like shell)
     #
     # name_of-the-shell  material rMin rMax starting_phi delta_phi starting_theta delta_theta is_active
     #
     # Example 1 - one unique shell (the second one is then composed of air and inactive)
     #	Inner	NaI 	10.	15.	0.	360.	0.	180.	1
     #	Outer	AIR 	25.	40.	0.	360.	0.	180.	0
     #
     # Example 2 - two active shells inner of NaI and outer BGO
     #	Inner	NaI 	10.	15.	0.	360.	0.	180.	1
     #	Outer	BGO 	25.	40.	0.	360.	0.	180.	1
     #
     # Example 3 - two active shells, the target chamber and a dead layer between the inner and outer
     #	Target	Al	8.	8.5	0.	360.	0.	180.	0
     #	Inner	NaI 	10.	15.	0.	360.	0.	180.	1
     #	Abs1	Pb	16.0	16.05	0.	360.	0.	180.	0
     #	Outer	BGO 	25.	40.	0.	360.	0.	180.	1
     #
     Inner	NaI 	10.	20.	0.	360.	0.	180.	1
     Outer	BGO 	25.	45.	0.	360.	0.	90.	1
     #
     #
     #
     \endcode
     As it can be seen, absorbers could be added easily. In case you would like to study only one shell, just set the
     material of the other one to Air (it will be automatically set to inactive).
     
     Here is picture with 2 actives shells and 2 inactives (Example 3).
     
     
     */
    class TwoShellsDetectorConstruction : public G4VUserDetectorConstruction
    {
    public:
        //! constructor
        /*!
         It reads the file named geometries/shells.geo. If not found, default are:
         \code
         #	Inner	NaI 	10.	15.	0.	360.	0.	180.	1
         #	Outer	BGO 	25.	40.	0.	360.	0.	180.	1
         \endcode
         */
        TwoShellsDetectorConstruction();
        TwoShellsDetectorConstruction(G4String filename);
        virtual ~TwoShellsDetectorConstruction();
        
        //! One of the mandatory class to be implemented in order to have G4 working properly
        virtual G4VPhysicalVolume* Construct();
        //! NEW G4.10 ... but also define for G4.9 except is is called explicitely at the end of Construct
        virtual void ConstructSDandField();
        
    private:
        //! Structure that defines a shell
        /*!
         A shell completely defines the parameters to build a G4Sphere
         */
        struct aShell
        {
            aShell();
            aShell(const aShell &);
            
            G4String Name;
            G4String MatName;
            
            G4double RMin;
            G4double RMax;
            G4double PhiStart;
            G4double PhiDelta;
            G4double ThetaStart;
            G4double ThetaDelta;
            
            G4int    IsActive;
            
            void Print(std::ostream &);
        };
        aShell Inner;
        aShell Outer;
        std::vector<aShell *> otherShells;
        
    public:
        //! from a given file, it computes the parameters of the two shells
        /*!
         The parameters are read from a file which is by default steup/shells.geo. In this file, the geometry
         of the two shells is completely defined. Some checkings are done at this level. For instance:
         - if the two shells overlap. In this case, the inner one is kept are the outer one truncated
         - if there are errors while reading the file. In this case, the properties are not modified	
         */	     
        void ComputeParameters(G4String filename = "DetectorFactory/Generics/TwoShells.geo");
        
    private:
        G4Box*             solidWorld;  // pointer to the solid envelope 
        G4LogicalVolume*   logicWorld;  // pointer to the logical envelope
        G4VPhysicalVolume* physiWorld;  // pointer to the physical envelope
        
        G4LogicalVolume*   logicInner;  // pointer to the logical envelope
        G4LogicalVolume*   logicOuter;  // pointer to the logical envelope
    };
} // SToGS Namespace

#endif
