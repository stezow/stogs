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

#ifndef SToGS_SpecialsDF_h
#define SToGS_SpecialsDF_h 1

#include "SToGS_DetectorFactory.hh"

class G4VUserDetectorConstruction;

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! Factory which is an interface to G4VUserDetectorConstruction
    /*!
     */
    class SpecialsDF : public DetectorFactory
    {
    protected:
        std::vector < std::pair < G4String, G4VUserDetectorConstruction *> > fLoadedUserDetectorConstruction;

    public:
        SpecialsDF(G4String path) : DetectorFactory(path), fLoadedUserDetectorConstruction()
            {;}
        virtual ~SpecialsDF()
            {;}
        
        //! overwritte. 
        virtual G4VPhysicalVolume *Get(G4String basename, G4bool is_full = true);
        //! Read the amap file and apply atributes to the detector. if not found, it creates a deefault one from the sensitive detector founds
        virtual void GetAttributes(G4String basename);

        //! build the default store i.e. all the scintillators detectors.
        virtual void MakeStore()
            {;}
    };
} // SToGS Namespace

#endif


