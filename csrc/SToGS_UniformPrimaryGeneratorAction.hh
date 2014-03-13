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
 
#ifndef ParisUniformPrimaryGeneratorAction_h
#define ParisUniformPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

// std includes 
#include <iostream>
#include <sstream>
#include <fstream>

class G4ParticleGun;
class G4Event;
class ParisUniformPrimaryGeneratorMessanger;


//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! It generates gamma-rays using an uniform distribution
/*!
	This generator read a setup file ( \c setup/uniform.gene) to generate a cascade of gamma-rays. A template
	of such a file is available is the setup directory (uniform.gene.demo). Here is a snapshot:	
\code
#
# Ascii file that described the particle cascade shoot in a given event for ParisUniformPrimaryGeneratorAction
# This is a basic UNIFORM photon generator !! I.e. only gammas are allowed and each of them is emitted from the same source point
# The number of gammas per cascade is free but can not exceed 100 (see below)
# The first line gives information on the source position (x, y, z) in the lab frame (experimental hall frame) as well as 
# its velocity components (vx,vy,vz) required for including Doppler effects (note that z is the direction of the beam).
# For sake of simplicity, the velocity components are given in units of the speed of light i.e. (beta_x, beta_y, beta_z).
# The second line gives the number of gamma to be shoot (multiplicity max of 100)
# The third line indicates the energy range over which the gammas will be randomly sampled as well as their emission cone  
#
# source position and velocity/c
0. 0. 0. 0. 0. 0.0
# number of particles
1
# cascade: energy (in KeV) min and max theta_min theta_max phi_min phi_max
0. 25000. 0. 180. 0. 360.
#
\endcode
	\warning The multiplicity of the cascade is restricted to 100.
*/ 
class ParisUniformPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
   ParisUniformPrimaryGeneratorAction();    
   ParisUniformPrimaryGeneratorAction(G4String);       
   virtual ~ParisUniformPrimaryGeneratorAction();

public:
    virtual void GeneratePrimaries(G4Event* anEvent);
  
public:     
     //! from a given file, it reads the characteristics of the gamma cascade 
     /*!
     	The parameters are read from a file which is by default setup/uniform.gene. In this file, the position
.	of the emitting gamma source, the angular coverage of the emission cone and the gamma energy can be modified
     */	     
	void ComputeParameters(G4String filename = "setup/uniform.gene");
	
	void ChangeMultiplicity(G4int smult);    

private:
  // Maximum gamma multiplicity in a cascade
    static const G4int MultMax = 100; 
    G4ParticleGun* particleGun;
    G4int mult;
    G4double Posx;
    G4double Posy;
    G4double Posz;
    G4double betax;
    G4double betay;
    G4double betaz;
    G4double Emin, Emax, Thetamin, Thetamax, Phimin, Phimax;
    G4double Egamma[MultMax];
    G4double EgammaDoppler[MultMax];
    
	ParisUniformPrimaryGeneratorMessanger *theMessanger;
};

#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//! Messanger class for ParisBasicPrimaryGenerator
/*!
 The messanger class allows to reset the generator inside a Geant4 session with for instance:
 
 \c /Paris/generator/BasicPrimaryGenerator \c setup/mynewbasic.gene
 
 */
class ParisUniformPrimaryGeneratorMessanger: public G4UImessenger
{
public:
	ParisUniformPrimaryGeneratorMessanger(ParisUniformPrimaryGeneratorAction *);
	~ParisUniformPrimaryGeneratorMessanger();
	
	void SetNewValue(G4UIcommand*, G4String);
private:
	ParisUniformPrimaryGeneratorAction *theGenerator;
	
	G4UIdirectory *theDirectory;
	G4UIcmdWithAString *resetParametersCmd;
	G4UIcmdWithAnInteger *multCmd;
};

} // SToGS Namespace

#endif
