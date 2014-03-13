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
 
#ifndef ParisBasicPrimaryGeneratorAction_h
#define ParisBasicPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

// std includes 
#include <iostream>
#include <sstream>
#include <fstream>

class G4ParticleGun;
class G4Event;
class ParisBasicPrimaryGeneratorActionMessanger;
 
//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! Basic generator to generate a cascade of discrete gamma-rays
/*!
	This generator read a setup file to generate a cascade of discrete gamma-rays. 
	Here is a snapshot of such a file:
\code
#
# Ascii file that described the particle cascade shoot in a given event for HermeBasicPrimaryGeneratoraction
# This is a BASIC photon generator !! I.e. only gammas are allowed and each of them is emitted from the same source point
# The number of gammas per cascade is free but can not exceed 100
# The first line gives information on the source position (x, y, z) in the lab frame (experimental hall frame) as well as 
# its velocity components (vx,vy,vz) required for including Doppler effects (note that z is the direction of the beam).
# For sake of simplicity, the velocity components are given in units of the speed of light i.e. (beta_x, beta_y, beta_z).
# The second line is a flag number depending on the information given below per particle
# 	- 0 : the emission cone is given for the first gamma only (similar to all gammas)
#       - 1 : gamma energy and angular emission cone are indicated for each particle
# NB: the above second option might be useful in case the shells of the detector do not cover 4pi but also for investigating the influence
#     of overlapping - or not - shells depending on the gamma energy (e.g. low vs. high energy gammas)  
#
# source position and velocity/c
0. 0. 0.  0. 0. 0.000
# file format
0
# cascade: energy (in KeV) theta_min theta_max phi_min phi_max (depending on the format)
25000. 0. 180. 0. 360. 
300 
400. 
500.
600. 
700 
800. 
900  
#
\endcode	
	\warning The multiplicity of the cascade is restricted to 100.
*/ 
class ParisBasicPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	ParisBasicPrimaryGeneratorAction();  
	ParisBasicPrimaryGeneratorAction(G4String);      
	virtual ~ParisBasicPrimaryGeneratorAction();

public:
	virtual void GeneratePrimaries(G4Event* anEvent);
  
public:     
	//! from a given file, it reads the characteristics of the gamma cascade 
	/*!
     	The parameters are read from a file which is by default setup/basic.gene. In this file, the position
.	of the emitting gamma source, the angular coverage of the emission cone and the gamma energy can be modified
	*/	     
	void ComputeParameters(G4String filename = "setup/basic.gene");
    
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
	G4double Egamma[MultMax];
	G4double EgammaDoppler[MultMax];
	G4double Theta_min[MultMax];
	G4double Theta_max[MultMax];
	G4double Phi_min[MultMax];
	G4double Phi_max[MultMax];    
	
	ParisBasicPrimaryGeneratorActionMessanger *theMessanger;
};

#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;

//! Messanger class for ParisBasicPrimaryGenerator
/*!
	The messanger class allows to reset the generator inside a Geant4 session with for instance:
	
	 \c /Paris/generator/BasicPrimaryGenerator \c setup/mynewbasic.gene
	
*/
class ParisBasicPrimaryGeneratorActionMessanger: public G4UImessenger
{
public:
	ParisBasicPrimaryGeneratorActionMessanger(ParisBasicPrimaryGeneratorAction *);
	~ParisBasicPrimaryGeneratorActionMessanger();
		
	void SetNewValue(G4UIcommand*, G4String);
	
private:
	ParisBasicPrimaryGeneratorAction *theGenerator;
	
	G4UIdirectory *theDirectory;
	G4UIcmdWithAString *resetParametersCmd;
};
    
} // SToGS Namespace

#endif


