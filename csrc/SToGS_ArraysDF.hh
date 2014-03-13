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

#ifndef SToGS_ArraysDF_h
#define SToGS_ArraysDF_h 1

// includes
#include "SToGS_DetectorFactory.hh"

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! base class for building standard arrays from pieces in other factories
    /*!
     */
    class ArraysDF : public DetectorFactory
    {
    protected:
        //! Should be implemented in any sub factory. It built (C++) a detector and return it
        virtual G4VPhysicalVolume * Make(G4String /* name */, G4String /* version_string */);
        
    public:
        ArraysDF(G4String path = "Arrays/") : DetectorFactory(path)
            {;}
        virtual ~ArraysDF()
            {;}
        
        //! to be filled with default arrays
        virtual void MakeStore();
    };
} // SToGS Namespace

#endif


