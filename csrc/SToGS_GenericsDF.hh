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

#ifndef SToGS_GenericsDF_h
#define SToGS_GenericsDF_h 1

#include "SToGS_DetectorFactory.hh"

class G4VUserDetectorConstruction;

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! Factory which is an interface to G4VUserDetectorConstruction
    /*!
     */
    class GenericsDF : public DetectorFactory
    {
    protected:
        std::vector < std::pair < G4String, G4VUserDetectorConstruction *> > fLoadedUserDetectorConstruction;
        
    protected:
        //! return from the full path in the factory the configuration file which is following $
        /*!
            ex: DetectorFactory/Generics/TwoShell$toto.geo
         */
        G4String GetConf(G4String path_in_factory) const;

    public:
        GenericsDF(G4String path) : DetectorFactory(path), fLoadedUserDetectorConstruction()
            {;}
        virtual ~GenericsDF()
            {;}
        
        //! overwrite the Get method to retrieve the detector from the standard G4 way i.e. by calling Construct
        virtual G4VPhysicalVolume *Get(G4String basename);
        //! Read attrbiutes by calling COnstructSDandFileds
        virtual void GetAttributes(G4String basename, G4bool do_amap = true, G4bool do_dmap = false);

        //! build the default store i.e. nothing here
        virtual void MakeStore()
            {;}
    };
} // SToGS Namespace

#endif


