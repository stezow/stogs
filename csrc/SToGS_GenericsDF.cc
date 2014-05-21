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

// G4 includes
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4Version.hh"

// Project includes
#include "SToGS_GenericsDF.hh"
#include "SToGSConfig.hh"

// Specific to the Generic Factory
#include "SToGS_TwoShellsDetectorConstruction.hh"
#include "SToGS_AGATA.hh"

// list of all specific factories
namespace  {
    // all scintillators
    SToGS::GenericsDF theGenericFactory("Generics/");
}

G4String SToGS::GenericsDF::GetConf(G4String path_in_factory) const
{
    G4String result(""), tmp = path_in_factory;
    if ( !tmp.contains("$") ) {
        return result;
    }
    else { result = GetFactoryName(); }
            
    G4int start = tmp.last('$'); tmp.remove(0,start+1);
    result += tmp;

    return result;
}

G4VPhysicalVolume *SToGS::GenericsDF::Get(G4String basename, G4bool is_full)
{
    G4VUserDetectorConstruction *theConstructor = 0x0; G4VPhysicalVolume *theDetector = 0x0; G4String conffile;
    
    // check is already loaded
    for (size_t i = 0; i < fLoadedPhysical.size(); i++) {
        if ( fLoadedPhysical[i].first == basename ) {
            theDetector = fLoadedPhysical[i].second;
            break;
        }
    }
    if ( theDetector ) {
        return theDetector;
    }
    else { // get conf name. If null call default constructor otherwise call constructor with string parameter
        conffile = GetConf(basename);
    }
    // load using the g4 facility
    if ( basename.contains("/TwoShells") ) {
        if ( conffile == "" ) {
            theConstructor = new SToGS::TwoShellsDetectorConstruction(conffile);
        }
        else theConstructor = new SToGS::TwoShellsDetectorConstruction(conffile);
        theDetector = theConstructor->Construct();
    }
    if ( basename.contains("/AGATA") ) {
        theConstructor = new SToGS::AGATA();
        theDetector = theConstructor->Construct();
    }
#ifdef HAS_MYDET
    // based on My plugins defined in SToGSConfig, it build the detector name in factory
    if ( basename.contains(MYDET_) ) {
        if ( conffile == "" ) {
            theConstructor = new MYDET_CLASSTYPE();
        }
        else theConstructor = new MYDET_CLASSTYPE(conffile);
        theDetector = theConstructor->Construct();
    }
#endif
    if ( theDetector ) {
        std::pair < G4String, G4VPhysicalVolume *> p1(basename,theDetector); // add the new loaded detector to the list
        fLoadedPhysical.push_back(p1);
        std::pair < G4String, G4VUserDetectorConstruction *> p2(basename,theConstructor);
        fLoadedUserDetectorConstruction.push_back(p2);
        
        if ( is_full ) {
#if G4VERSION_NUMBER < 1000
#else
            theConstructor->ConstructSDandField();
#endif
        }
    }

    return theDetector;
}

void SToGS::GenericsDF::GetAttributes(G4String basename)
{
    // check is already loaded
#if G4VERSION_NUMBER < 1000
#else
    for (size_t i = 0; i < fLoadedUserDetectorConstruction.size(); i++) {
        if ( fLoadedUserDetectorConstruction[i].first == basename ) {
            fLoadedUserDetectorConstruction[i].second->ConstructSDandField();
            break;
        }
    }
#endif
}






