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

#include "SToGS_LoadFromDetectorFactory.hh"
#include "SToGS_DetectorFactory.hh"

SToGS::LoadFromDetectorFactory::~LoadFromDetectorFactory()
{
}

G4VPhysicalVolume *SToGS::LoadFromDetectorFactory::Construct()
{
    G4cout << G4endl << " ------ INF ------ from SToGS::LoadFromDetectorFactory::Construct " << fNameInFactory << G4endl;

    G4VPhysicalVolume *physiWorld = 0x0; SToGS::DetectorFactory *where_to_load = SToGS::DetectorFactory::GetFactory(fNameInFactory);
    if ( where_to_load ) {
        physiWorld = where_to_load->Get(fNameInFactory);
        if ( physiWorld == 0x0 ) {
            G4cout << "**** Cannot find setup " << fNameInFactory << " in Detector Factory " << where_to_load->GetFactoryName() << G4endl;
        }
    }
    else G4cout << "**** Cannot load setup " << fNameInFactory << " from Detector Factory " << G4endl;
    
    G4cout << " ------ END ------ from SToGS::LoadFromDetectorFactory::Construct " << G4endl;

    return physiWorld;
}

void SToGS::LoadFromDetectorFactory::ConstructSDandField()
{
    G4cout << " ------ INF ------ from SToGS::LoadFromDetectorFactory::ConstructSDandField " << fNameInFactory << G4endl;

    SToGS::DetectorFactory *where_to_load = SToGS::DetectorFactory::GetFactory(fNameInFactory);
    if ( where_to_load ) {
        where_to_load->GetAttributes(fNameInFactory,true,false); // assign amap, dmap
    }
    else
        G4cout << "**** Cannot load setup " << fNameInFactory << " from Detector Factory " << G4endl;
    
    G4cout << " ------ END ------ from SToGS::LoadFromDetectorFactory::ConstructSDandField " << G4endl;
}

SToGS::BuildFromDetectorFactory::~BuildFromDetectorFactory()
{
}

G4VPhysicalVolume *SToGS::BuildFromDetectorFactory::Construct()
{
    G4VPhysicalVolume *physiWorld = 0x0;
    
    SToGS::DetectorFactory *where_to_load = SToGS::DetectorFactory::GetFactory(fInputFile);
    if ( where_to_load == 0x0 ) {
        where_to_load = SToGS::DetectorFactory::GetFactory("DetectorFactory/MyStore/");
    }
    if ( where_to_load ) {
        physiWorld = where_to_load->MakeAnArrayFromFactory(fInputFile);
    }
    else
        G4cout << "**** Cannot open file " << fInputFile << " Detector Factory " << G4endl;
    
    return physiWorld;
}


