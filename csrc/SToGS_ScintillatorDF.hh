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

#ifndef SToGS_ScintillatorDF_h
#define SToGS_ScintillatorDF_h 1

#include "SToGS_DetectorFactory.hh"

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! Base classe to build scintillators
    /*!
     */
    class ScintillatorDF : public DetectorFactory
    {
    public:
        // All individual methods to build some detector units
        // make PARIS PW with/without encapuslation, housing
        /*!
         Options \n
         bare : take into account caps and housing but they are not included in the PW. It allows studies of the impact of encapsulation
         */
        G4VPhysicalVolume *MakePPW(G4String detname, G4double caps_width = 0.0, G4double housing_width = 0.0, G4String opt = "bare");
        // make paris cluster from
        G4VPhysicalVolume *MakeCPPW(G4String detname, G4String opt = "bare");
        
        // make CHATEAU DE CRYTSAL module
        G4VPhysicalVolume *MakeCdC(G4String detname, G4String opt = "bare");
        
        // make a Fatima module
        G4VPhysicalVolume *MakeFATIMAM(G4String detname, G4String opt = "bare");
        G4VPhysicalVolume *MakeFATIMAQ(G4String detname, G4String opt = "bare");

        // make EDEN
        G4VPhysicalVolume *MakeEDEN(G4String detname, G4double caps_width = 0.0, G4double extra_back = 0.0, G4String opt = "v0");
        
    public:
        //! Should be implemented in any sub factory. It built (C++) a detector and return it
        virtual G4VPhysicalVolume *Make(G4String /* name */, G4String /* version_string */);
        
    public:
        ScintillatorDF(G4String path) : DetectorFactory(path)
            {;}
        virtual ~ScintillatorDF()
            {;}
        
        //! build the default store i.e. all the scintillators detectors.
        virtual void MakeStore();
    };
} // SToGS Namespace

#endif


