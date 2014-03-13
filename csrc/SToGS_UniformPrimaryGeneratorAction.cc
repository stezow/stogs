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

#include "ParisUniformPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

using namespace std; 

ParisUniformPrimaryGeneratorAction::ParisUniformPrimaryGeneratorAction()
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
  Emin = 0.;
  Emax = 0.;
  Thetamin = 0.;
  Thetamax = 0.;
  Phimin = 0.;
  Phimax = 0.;
  for (G4int imult = 0; imult < mult; imult++ ) {
  	Egamma[imult] = 0.;
	EgammaDoppler[imult] = 0.;
  }  
    
// default values (isotropic emission in 4pi of 1 gamma whose energy is sampled in a flat distribution between 100keV and 10MeV)
  mult = 1;
  Posx = 0. *cm;
  Posy = 0. *cm;
  Posz = 0. *cm;
  betax = 0.;
  betay = 0.;
  betaz = 0.;
  Emin = 100. *keV;
  Emax = 10000. *keV;
  Thetamin = 0. *deg;
  Thetamax = 180.*deg;
  Phimin = 0. *deg;
  Phimax = 360. *deg;
  
// particle cascade read from file  
  ComputeParameters();   
    
	theMessanger = new ParisUniformPrimaryGeneratorMessanger(this);
	
}  
 
ParisUniformPrimaryGeneratorAction::ParisUniformPrimaryGeneratorAction(G4String filename)
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
  Emin = 0.;
  Emax = 0.;
  Thetamin = 0.;
  Thetamax = 0.;
  Phimin = 0.;
  Phimax = 0.;
  for (G4int imult = 0; imult < mult; imult++ ) {
  	Egamma[imult] = 0.;
	EgammaDoppler[imult] = 0.;
  }  
    
// default values (isotropic emission in 4pi of 1 gamma whose energy is sampled in a flat distribution between 100keV and 10MeV)
  mult = 1;
  Posx = 0. *cm;
  Posy = 0. *cm;
  Posz = 0. *cm;
  betax = 0.;
  betay = 0.;
  betaz = 0.;
  Emin = 100. *keV;
  Emax = 10000. *keV;
  Thetamin = 0. *deg;
  Thetamax = 180.*deg;
  Phimin = 0. *deg;
  Phimax = 360. *deg;
  
// particle cascade read from file  
  ComputeParameters(filename);   
	
	theMessanger = new ParisUniformPrimaryGeneratorMessanger(this);
}    

ParisUniformPrimaryGeneratorAction::~ParisUniformPrimaryGeneratorAction()
{
  delete particleGun;
    delete theMessanger;
}

void ParisUniformPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{ 

// temporary file for sake of check and test
//   ofstream sortie_test("test.dat",ios::out);
   
// generation of the primary particles
  particleGun->SetParticlePosition(G4ThreeVector(Posx, Posy, Posz));  
  for (G4int imult = 0; imult < mult; imult++ ) {  

	// randomization of the particle direction : uniforn distribution within the solid angle sustained by the particle
	G4double phi = Phimin + (Phimax-Phimin) * G4UniformRand();
	G4double cosThetamin = std::cos(Thetamin), cosThetamax = std::cos(Thetamax);
	G4double cosTheta = cosThetamin + (cosThetamax-cosThetamin) * G4UniformRand();  
// recording in test output file
//		sortie_test << phi/deg << "    " << cosTheta << endl;
        G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
        G4double ux = sinTheta*std::cos(phi),
                 uy = sinTheta*std::sin(phi),
                 uz = cosTheta;
	
	// sampling of the particle energy : uniforn distribution between Emin and Emax
	Egamma[imult] = Emin + (Emax-Emin) * G4UniformRand(); 
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
//		sortie_test << phi/deg << "  " << cosTheta << "  " << Egamma[imult] << endl;

  } 
 
 }


void ParisUniformPrimaryGeneratorAction::ComputeParameters(G4String filename)
{
	
	G4cout << G4endl << " ------ INFO ------ from ParisUniformPrimaryGeneratorAction::ComputeParameters " << G4endl;
	
	// open the ascii file 
	ifstream file;
	
	file.open(filename.data());
	if ( file.is_open() == false ) { 
		G4cout << " ** WARNING ** cannot open file " << filename << " (Default parameters are used) "<< G4endl;
		G4cout << "  --> Default parameters are used for the generator" << G4endl;
		G4cout << mult << " " << particleGun->GetParticleDefinition() << " with energy between" << Emin << " keV and " << Emax << "keV in (thetamin,thetamax,phimin,phimax) " << Thetamin << " " << Thetamax << " "  << Phimin << " " << Phimax << G4endl;
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
			// from the line extract the number of particles to be shoot
			sscanf(aline, "%d", &mult);
			ilect += 1;
			file.getline(aline,MAXWIDTH); 
			continue;
		}
			// from the line extract the characteristics of the energy and angular window
  		if (ilect == 2) {
			// energy and angular range to be sampled
			sscanf(aline, "%f %f %f %f %f %f", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6);
			Emin = tmp1 *keV;
			Emax = tmp2 *keV;
			Thetamin = tmp3 *deg;
			Thetamax = tmp4 *deg;
			Phimin = tmp5 *deg;
			Phimax = tmp6 *deg;
		}
				
		file.getline(aline,MAXWIDTH);
		
	}
	file.close();
	
	// check reading of the input file 
  	G4cout << "Multiplicity = " <<  mult << G4endl;
  	G4cout << "Source position = " <<  Posx/cm << " " << Posy/cm << " " << Posz/cm  << G4endl;
  	G4cout << "Source velocity/c = " <<  betax << " " << betay << " " << betaz  << G4endl;
  	G4cout << "Particle characteristics (energy, angles) = " << G4endl;
	G4cout << "Energy range in keV = " <<  Emin/keV << "   " << Emax/keV << G4endl;
	G4cout << "Angular range in deg = " <<  Thetamin/deg << "  " << Thetamax/deg << "  " << Phimin/deg << "  " << Phimax/deg << G4endl;
	
	G4cout << " ------ END ------ from ParisUniformPrimaryGeneratorAction::ComputeParameters " << G4endl << G4endl;
}

void ParisUniformPrimaryGeneratorAction::ChangeMultiplicity(G4int smult)
{
	G4cout << " ------ INFO ------ Changing multiplicity in ParisUniformPrimaryGeneratorAction " << G4endl;
	mult = smult;	
}   

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

ParisUniformPrimaryGeneratorMessanger::ParisUniformPrimaryGeneratorMessanger(ParisUniformPrimaryGeneratorAction *thegene): theGenerator(thegene)
{	
	theDirectory = new G4UIdirectory("/Paris/generator/");
	theDirectory->SetGuidance("To modify the parameters");
	
	resetParametersCmd = new G4UIcmdWithAString("/Paris/generator/file", this);
	resetParametersCmd->SetGuidance("To load new parameters of ParisRandomCascade from a file");
	resetParametersCmd->SetGuidance("Required parameters: one valid filename (string)");
	resetParametersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
	
	multCmd = new G4UIcmdWithAnInteger("/Paris/generator/mult", this);
	multCmd->SetGuidance("To change the multiplicity of the generated cascade");
	multCmd->SetGuidance("Required parameters: an integer");
	multCmd->AvailableForStates(G4State_PreInit,G4State_Idle);	
	
}

ParisUniformPrimaryGeneratorMessanger::~ParisUniformPrimaryGeneratorMessanger()
{
	delete theDirectory; delete resetParametersCmd; delete multCmd;
}
void ParisUniformPrimaryGeneratorMessanger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if( command == resetParametersCmd ) 
		theGenerator->ComputeParameters( newValue );
	if( command == multCmd ) 
		theGenerator->ChangeMultiplicity( multCmd->GetNewIntValue(newValue) );	
	
}




