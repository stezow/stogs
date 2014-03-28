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

#include "ParisBasicPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

using namespace std; 

ParisBasicPrimaryGeneratorAction::ParisBasicPrimaryGeneratorAction()
{

// restriction to gamma particles and 1 gamma per gun 
  particleGun = new G4ParticleGun(1);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));

// initialization
  mult = MultMax;
  Posx = 0.;
  Posy = 0.;
  Posz = 0.;
  betax = 0.;
  betay = 0.;
  betaz = 0.;
  for (G4int imult = 0; imult < mult; imult++ ) {
  	Egamma[imult] = 0.;
	EgammaDoppler[imult] = 0.;
  	Theta_min[imult] = 0.;
  	Theta_max[imult] = 0.;
  	Phi_min[imult] = 0.;
  	Phi_max[imult] = 0.;
  }  
    
// default values (isotropic emission in 4pi of one gamma with energy 200keV)
  mult = 1;
  Posx = 0. *cm;
  Posy = 0. *cm;
  Posz = 0. *cm;
  betax = 0.;
  betay = 0.;
  betaz = 0.;
  Egamma[0] = 200. *keV;
  Theta_min[0] = 0. *deg;
  Theta_max[0] = 180.*deg;
  Phi_min[0] = 0. *deg;
  Phi_max[0] = 360. *deg;
  
// particle cascade read from file  
  ComputeParameters();   
    
  theMessanger = new ParisBasicPrimaryGeneratorActionMessanger(this);
}

ParisBasicPrimaryGeneratorAction::ParisBasicPrimaryGeneratorAction(G4String filename)
{

// restriction to gamma particles and 1 gamma per gun 
  particleGun = new G4ParticleGun(1);
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="gamma"));

// initialization
  mult = MultMax;
  Posx = 0.;
  Posy = 0.;
  Posz = 0.;
  betax = 0.;
  betay = 0.;
  betaz = 0.;
  for (G4int imult = 0; imult < mult; imult++ ) {
  	Egamma[imult] = 0.;
	EgammaDoppler[imult] = 0.;
  	Theta_min[imult] = 0.;
  	Theta_max[imult] = 0.;
  	Phi_min[imult] = 0.;
  	Phi_max[imult] = 0.;
  }  
    
// default values (isotropic emission in 4pi of one gamma with energy 200keV)
  mult = 1;
  Posx = 0. *cm;
  Posy = 0. *cm;
  Posz = 0. *cm;
  betax = 0.;
  betay = 0.;
  betaz = 0.;
  Egamma[0] = 200. *keV;
  Theta_min[0] = 0. *deg;
  Theta_max[0] = 180.*deg;
  Phi_min[0] = 0. *deg;
  Phi_max[0] = 360. *deg;
  
// particle cascade read from file  
  ComputeParameters(filename);   
    
  theMessanger = new ParisBasicPrimaryGeneratorActionMessanger(this);
}

ParisBasicPrimaryGeneratorAction::~ParisBasicPrimaryGeneratorAction()
{
  delete particleGun; delete theMessanger;
}

void ParisBasicPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 

// temporary file for sake of check and test
//   ofstream sortie_test("test.dat",ios::out);
   
// generation of the primary particles
  particleGun->SetParticlePosition(G4ThreeVector(Posx, Posy, Posz));  
  for (G4int imult = 0; imult < mult; imult++ ) {  

	// randomization of the particle direction : uniforn distribution within the solid angle sustained by the particle
	G4double phi = Phi_min[imult] + (Phi_max[imult]-Phi_min[imult]) * G4UniformRand();
	G4double cosThetamin = std::cos(Theta_min[imult]), cosThetamax = std::cos(Theta_max[imult]);
	G4double cosTheta = cosThetamin + (cosThetamax-cosThetamin) * G4UniformRand();  
// recording in test output file
//		sortie_test << phi/deg << "    " << cosTheta << endl;
        G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
        G4double ux = sinTheta*std::cos(phi),
                 uy = sinTheta*std::sin(phi),
                 uz = cosTheta;
	// particle energy accouting for Doppler effect
	G4double beta=0.;
	G4double bgamma = 1.;
	beta = std::sqrt(betax*betax + betay*betay + betaz*betaz);
	G4double ThetaCNgamma =  std::acos(cosTheta);
	if (beta > 0.) {
		ThetaCNgamma = std::acos(cosTheta) - std::acos(betaz/beta);
	} 
	bgamma = 1. / (std::sqrt( (1.-beta*beta) ) );
	EgammaDoppler[imult] = Egamma[imult]/bgamma/(1.-beta*(std::cos(ThetaCNgamma)));
	
	particleGun->SetParticleEnergy(EgammaDoppler[imult]);
// example for printing out 
// 	G4cout <<  imult+1 << " " <<  Egamma[imult]/keV << " " << EgammaDoppler[imult]/keV << G4endl;
        particleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));
        particleGun->GeneratePrimaryVertex(anEvent);

  } 
 
 }


void ParisBasicPrimaryGeneratorAction::ComputeParameters(G4String filename)
{
	G4cout << G4endl << " ------ INF ------ from ParisBasicPrimaryGeneratorAction::ComputeParameters " << G4endl;
	
	// open the ascii file 
	ifstream file; 
	G4int fileformat;
	
	file.open(filename.data());
	if ( file.is_open() == false ) { 
		G4cout << " ** WARNING ** cannot open file " << filename << " (Default parameters are used) "<< G4endl;
		G4cout << "  --> Default parameters are used for the generator" << G4endl;
		G4cout << mult << " " << particleGun->GetParticleDefinition()->GetParticleName() << " with energy " << Egamma[0] << " keV emitted in (theta_min,theta_max,phi_min,phi_max) " << Theta_min[0] << " " << Theta_max[0] << " "  << Phi_min[0] << " " << Phi_max[0] << G4endl;
		return;
	}
	
	// read the file and get the cascade
	const G4int MAXWIDTH = 300; char aline[MAXWIDTH];
	G4float tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
	mult = 0;
	file.getline(aline,MAXWIDTH);
	G4int ilect = 0;
	while ( file.good() ) {

		if ( aline[0] == '#' ) { file.getline(aline,MAXWIDTH); continue; } // this line is a comment 
		
		if (ilect == 0) {
			// from the line extract the source position
			sscanf(aline, "%f %f %f %f %f %f", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6);
			Posx = tmp1 *cm;
			Posy = tmp2 *cm;
			Posz = tmp3 *cm;
			betax = tmp4;
			betay = tmp5;
			betaz = tmp6;
			ilect += 1;
			file.getline(aline,MAXWIDTH); 
			continue;
		}
		if (ilect == 1) {
			// from the line extract the format of following information
			sscanf(aline, "%d", &fileformat);
			ilect += 1;
			file.getline(aline,MAXWIDTH); 
			continue;
		}
			// from the line extract the characteristics of the gamma of the cascade
  		switch( fileformat ){
  		case 0:
			if (mult == 0) {
				sscanf(aline, "%f %f %f %f %f", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5);
				Egamma[mult] = tmp1 *keV;
				Theta_min[mult] = tmp2 *deg;
				Theta_max[mult] = tmp3 *deg;
				Phi_min[mult] = tmp4 *deg;
				Phi_max[mult] = tmp5 *deg;
			}
			else {
				sscanf(aline, "%f", &tmp1);
				Egamma[mult] = tmp1 *keV;
				Theta_min[mult] = Theta_min[0];
				Theta_max[mult] = Theta_max[0];
				Phi_min[mult] = Phi_min[0];
				Phi_max[mult] = Phi_max[0];
			}
			mult += 1;
			break;	
  		case 1:
			sscanf(aline, "%f %f %f %f %f", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5);
			Egamma[mult] = tmp1 *keV;
			Theta_min[mult] = tmp2 *deg;
			Theta_max[mult] = tmp3 *deg;
			Phi_min[mult] = tmp4 *deg;
			Phi_max[mult] = tmp5 *deg;
			mult += 1;
			break;
		default:
			G4cout << " ** WARNING ** empty generator file " << G4endl;
			break;	
  		}	
		
		file.getline(aline,MAXWIDTH);
		
	}
	file.close();	
	
	// check reading of the input file 
  	G4cout << "Multiplicity = " <<  mult << G4endl;
  	G4cout << "Source position = " <<  Posx/cm << " " << Posy/cm << " " << Posz/cm  << G4endl;
  	G4cout << "Source velocity/c = " <<  betax << " " << betay << " " << betaz  << G4endl;
  	G4cout << "Particle characteristics (energy, angles) = " << G4endl;
  	for (G4int imult = 0; imult < mult; imult++ ) {
  		G4cout <<  imult+1 << " " <<  Egamma[imult]/keV << G4endl;
  		G4cout <<  Theta_min[imult]/deg << " " << Theta_max[imult]/deg << " " << Phi_min[imult]/deg << " " << Phi_max[imult]/deg << " " <<G4endl;
  	}
	
	G4cout << " ------ END ------ from ParisBasicPrimaryGeneratorAction::ComputeParameters " << G4endl << G4endl;
}

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

ParisBasicPrimaryGeneratorActionMessanger::ParisBasicPrimaryGeneratorActionMessanger(ParisBasicPrimaryGeneratorAction *thegene): theGenerator(thegene)
{	
	theDirectory = new G4UIdirectory("/Paris/generator/");
	theDirectory->SetGuidance("To modify generator's parameters");
	
	resetParametersCmd = new G4UIcmdWithAString("/Paris/generator/file", this);
	resetParametersCmd->SetGuidance("To load new parameters of ParisBasicPrimaryGeneratorAction from a file");
	resetParametersCmd->SetGuidance("Required parameters: one valid filename (string)");
	resetParametersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);	
}

ParisBasicPrimaryGeneratorActionMessanger::~ParisBasicPrimaryGeneratorActionMessanger()
{
	delete theDirectory; delete resetParametersCmd;
}

void ParisBasicPrimaryGeneratorActionMessanger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if( command == resetParametersCmd ) theGenerator->ComputeParameters( newValue );
}


