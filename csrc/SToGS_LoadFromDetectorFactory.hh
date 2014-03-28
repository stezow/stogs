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

#ifndef SToGS_LoadFromDetectorFactory_h
#define SToGS_LoadFromDetectorFactory_h 1

#include "G4VUserDetectorConstruction.hh"

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! to load a setup from the factory
    /*!
     */
    class LoadFromDetectorFactory : public G4VUserDetectorConstruction
    {
    public:
        LoadFromDetectorFactory(G4String name_in_factory = "") : G4VUserDetectorConstruction(), fNameInFactory(name_in_factory)
            {;}
        virtual ~LoadFromDetectorFactory();
        
        //! it built only the detector part i.e. load the xml file and the dmap
        virtual G4VPhysicalVolume *Construct();
        
        //! it set the sensitivity, colors and fields here
        virtual void ConstructSDandField();
        
    private:
        G4String fNameInFactory;
    };
    //! to load a setup from the factory
    /*!
     */
    class BuildFromDetectorFactory : public G4VUserDetectorConstruction
    {
    public:
        BuildFromDetectorFactory(G4String file_to_build = "default.setup") : G4VUserDetectorConstruction(), fInputFile(file_to_build)
            {;}
        virtual ~BuildFromDetectorFactory();
        
        //! it built only the detector part i.e. load the xml file and the amap, dmap
        virtual G4VPhysicalVolume *Construct();
        
    private:
        G4String fInputFile;
    };
} // SToGS Namespace


#endif
