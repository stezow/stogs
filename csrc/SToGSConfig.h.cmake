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
//----------------------------------------------------------------------------
 
/** \file SToGSConfig.hh to set informations about the SToGS configuration
*
*
*/ 
#ifndef SToGS_Config_h
#define SToGS_Config_h

/* Define to the version of this package. */
#define PACKAGE_VERSION @PROJECT_VERSION@

/* Define to the full name of this package. */
#define PACKAGE_NAME @PROJECT_NAME@

/* Define MY */
// MyDetectorConstruction
#ifdef HAS_MYDET
#include "@MY_DET@DetectorConstruction.hh"
#define MYDET_          "@MY_DET@"
#define MYDET_CLASSTYPE  @MY_DET@DetectorConstruction
#endif

// MyPrimaryGeneratorAction
#ifdef HAS_MYPRI
#include "@MY_PRI@PrimaryGeneratorAction.hh"
#define MYPRI_          "@MY_PRI@"
#define MYPRI_CLASSTYPE  @MY_PRI@PrimaryGeneratorAction
#endif

// MyAction
#ifdef HAS_MYACT
#include "@MY_ACT@UserAction.hh"
#define MYACT_          "@MY_ACT@"
#define MYACT_CLASSTYPE  @MY_ACT@UserActionInitialization
#endif

#endif
