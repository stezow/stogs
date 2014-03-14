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
//----------------------------------------------------------------------------------

#include "SToGS_MaterialConsultant.hh"
#include "G4MaterialTable.hh"
#include "globals.hh"

SToGS::MaterialConsultant *SToGS::MaterialConsultant::theMaterialConsultant = 0x0;

SToGS::MaterialConsultant::MaterialConsultant()
{
    G4double z,density; G4String name,symbol; G4int nel,natoms;
    
    //------------
    // Registering of all the needed element. Don't forget to do it if you use one
    // element to build a composite material
    //------------
    G4double a;
    
    a=1.01*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Hydrogen",symbol="H",z=1.,a));
    a=2.01*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Deuterium",symbol="D",z=1.,a));
    a=4.*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Helium",symbol="He",z=2.,a));
    a=6.94*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Lithium",symbol="Li",z=3.,a));
    a=9.01*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Berillium",symbol="Be",z=4.,a));
    a=12.01*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Carbon",symbol="C",z=6.,a));
    a=14.00674*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Nitrogen",symbol="N",z=7.,a));
    a=15.9994*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Oxygen",symbol="O",z=8.,a));
    a=20.18*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Neon",symbol="Ne",z=10.,a));
    a=22.99*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Sodium",symbol="Na",z=11.,a));
    a=24.305*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Magnesium",symbol="Mg",z=12.,a));
    a=26.98154*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Aluminium",symbol="Al",z=13.,a));
    a=28.085*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Silicon",symbol="Si",z=14.,a));
    a=35.453*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Chlorine",symbol="Cl",z=17.,a));
    // a=40.08*CLHEP::g/CLHEP::mole;
    // theElements.push_back(new G4Element(name="Calcium",symbol="Ca",z=20.,a));
    // a=54.938*CLHEP::g/CLHEP::mole;
    // theElements.push_back(new G4Element(name="Manganese",symbol="Mn",z=25.,a));
    a=55.850*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Iron",symbol="Fe",z=26.,a));
    a=63.546*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Copper",symbol="Cu",z=29.,a));
    a=137.327*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Barium",symbol="Ba",z=56.,a));
    a=18.9984032*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Fluorine",symbol="F",z=9.,a));
    
    // Germanium isotopes
    G4Isotope* Ge70 = new G4Isotope(name="Ge70", 32, 70, 69.9242*CLHEP::g/CLHEP::mole);
    G4Isotope* Ge72 = new G4Isotope(name="Ge72", 32, 72, 71.9221*CLHEP::g/CLHEP::mole);
    G4Isotope* Ge73 = new G4Isotope(name="Ge73", 32, 73, 72.9235*CLHEP::g/CLHEP::mole);
    G4Isotope* Ge74 = new G4Isotope(name="Ge74", 32, 74, 73.9212*CLHEP::g/CLHEP::mole);
    G4Isotope* Ge76 = new G4Isotope(name="Ge76", 32, 76, 75.9214*CLHEP::g/CLHEP::mole);
    // germanium defined via its isotopes
    G4Element* elGe = new G4Element(name="Germanium",symbol="Ge", 5);
    elGe->AddIsotope(Ge70, 0.2123);
    elGe->AddIsotope(Ge72, 0.2766);
    elGe->AddIsotope(Ge73, 0.0773);
    elGe->AddIsotope(Ge74, 0.3594);
    elGe->AddIsotope(Ge76, 0.0744);
    theElements.push_back(elGe);
    a=79.904*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Bromine",symbol="Br",z=35.,a));
    a=126.90477*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Iodine",symbol="I",z=53.,a));
    a=132.90545*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Cesium",symbol="Cs",z=55.,a));
    a=138.9055*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Lanthanum",symbol="La",z=57.,a));
    a=183.85*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Tungsten",symbol="W",z=74.,a));
    a=207.19*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Lead",symbol="Pb",z=82.,a));
    a=208.98038*CLHEP::g/CLHEP::mole;
    theElements.push_back(new G4Element(name="Bismuth",symbol="Bi",z=83.,a));
    // a=238.03*CLHEP::g/CLHEP::mole;
    // theElements.push_back(new G4Element(name="Uranium",symbol="U",z=92.,a));
    
    //-------------------
    // simple materials, mainly for absorbers
    //-------------------
    
    density = 2.7*CLHEP::g/CLHEP::cm3;
    a = 26.98*CLHEP::g/CLHEP::mole;
    theMaterials.push_back(new G4Material(name="Aluminium",z=13.,a,density));
    density = 7.87*CLHEP::g/CLHEP::cm3;
    a = 55.85*CLHEP::g/CLHEP::mole;
    theMaterials.push_back(new G4Material(name="Iron",z=26.,a,density));
    density = 8.96*CLHEP::g/CLHEP::cm3;
    a = 63.54*CLHEP::g/CLHEP::mole;
    theMaterials.push_back(new G4Material(name="Copper",z=29.,a,density));
    density = 19.3*CLHEP::g/CLHEP::cm3;
    a = 183.85*CLHEP::g/CLHEP::mole;
    theMaterials.push_back(new G4Material(name="Tungsten",z=74.,a,density));
    density = 11.35*CLHEP::g/CLHEP::cm3;
    a = 207.19*CLHEP::g/CLHEP::mole;
    theMaterials.push_back(new G4Material(name="Lead",z=82.,a,density));
    density = 8.96*CLHEP::g/CLHEP::cm3;
    a = 58.69*CLHEP::g/CLHEP::mole;
    theMaterials.push_back(new G4Material(name="Nickel",z=28.,a,density));
    
    //  density = 1.4*CLHEP::g/CLHEP::cm3;
    //  a = 39.95*CLHEP::g/CLHEP::mole;
    //  theMaterials.push_back(new G4Material(name="LiquidArgon",z=18.,a,density));
    //  density = 0.002*CLHEP::g/CLHEP::cm3;
    //  a = 39.95*CLHEP::g/CLHEP::mole;
    //  theMaterials.push_back(new G4Material(name="ArgonGas",z=18.,a,density));
    
    //------------------
    // mixtures
    //------------------
    
    G4Material *mix;
    
    density = 1.290*CLHEP::g/CLHEP::cm3; nel = 2; mix = new G4Material(name="Air",density,nel);
    mix->AddElement(GetElement("Nitrogen"), 70);
    mix->AddElement(GetElement("Oxygen"), 30);
    theMaterials.push_back(mix);
    
    density = 3.67*CLHEP::g/CLHEP::cm3, nel = 2; mix = new G4Material(name="NaI",density,nel);
    mix->AddElement(GetElement("Sodium"), natoms = 1);
    mix->AddElement(GetElement("Iodine"), natoms = 1);
    theMaterials.push_back(mix);
    
    density = 3.79*CLHEP::g/CLHEP::cm3, nel = 2; mix = new G4Material(name="LaCl3",density,nel);
    mix->AddElement(GetElement("Lanthanum"), natoms = 1);
    mix->AddElement(GetElement("Chlorine"), natoms = 3);
    theMaterials.push_back(mix);
    
    density  = 4.51*CLHEP::g/CLHEP::cm3, nel = 2; mix = new G4Material(name="CsI", density, nel);
    mix->AddElement(GetElement("Cesium"), natoms = 1);
    mix->AddElement(GetElement("Iodine"), natoms = 1);
    theMaterials.push_back(mix);
    
    density = 5.29*CLHEP::g/CLHEP::cm3, nel = 2; mix = new G4Material(name="LaBr3",density,nel);
    mix->AddElement(GetElement("Lanthanum"), natoms = 1);
    mix->AddElement(GetElement("Bromine"), natoms = 3);
    theMaterials.push_back(mix);
    
    density = 7.13*CLHEP::g/CLHEP::cm3, nel = 3; mix = new G4Material(name="BGO", density, nel);
    mix->AddElement(GetElement("Bismuth"), natoms = 4);
    mix->AddElement(GetElement("Germanium"), natoms = 3);
    mix->AddElement(GetElement("Oxygen"), natoms = 12);
    theMaterials.push_back(mix);
    
    density = 4.89*CLHEP::g/CLHEP::cm3, nel = 2; mix = new G4Material(name="BaF2", density, nel);
    mix->AddElement(GetElement("Barium"), natoms = 1);
    mix->AddElement(GetElement("Fluorine"), natoms = 2);
    
    // must have the right composition for stainless steel
    
    //  density = 8.96*CLHEP::g/CLHEP::cm3;
    //  StainlessSteel = new G4Material(name="StainlessSteel",density,nel=1);
    //  StainlessSteel->AddElement(elO, fractionmass = 1.);
    
    //  density              = 1.e-5*CLHEP::g/CLHEP::cm3;
    //  G4double pressure    = 2.e-2*bar;
    //  G4double temperature = STP_Temperature;         //from PhysicalConstants.h
    //  G4Material* Vacuum = new G4Material(name="Vacuum", density, nel=1,
    //				    kStateGas,temperature,pressure);
    //  Vacuum->AddMaterial(Air, fractionmass=1.);
    
}

SToGS::MaterialConsultant *SToGS::MaterialConsultant::theConsultant()
{
    if (theMaterialConsultant == 0x0) {
        theMaterialConsultant = new SToGS::MaterialConsultant();
    }
    return theMaterialConsultant;
}

G4Material *SToGS::MaterialConsultant::GetMaterial(G4String what) const
{
    G4Material* material = G4Material::GetMaterial(what);
    /*
     for(unsigned int i = 0; i < theMaterials.size() ; i++ ){
     if ( what == theMaterials[i]->GetName() ) { material = theMaterials[i]; break; }
     } */
    return material;
}

G4Element *SToGS::MaterialConsultant::GetElement(G4String what) const
{
    G4Element *element = G4Element::GetElement(what);
    /*
     for(unsigned int i = 0; i < theElements.size() ; i++ ){
     if ( what == theElements[i]->GetName() ) { element = theElements[i]; break; }
     } */
    return element;
}

void SToGS::MaterialConsultant::SetOpticalProperties(G4Material *mat, G4String which_properties)
{
	G4MaterialPropertiesTable *mat_mt = 0x0;
	//
	if ( which_properties == G4String("LaBr3") ) {
		
		const G4int NUMENTRIES = 3; G4double Energy[NUMENTRIES]    = { 2.0*CLHEP::eV , 5*CLHEP::eV, 9*CLHEP::eV };
		
		G4double SCINT[NUMENTRIES] = { 1.0, 1.0, 1.0 };
		G4double RIND[NUMENTRIES]  = { 1.9 , 1.9, 1.9 };
		
		G4double ABSL[NUMENTRIES]  = { 400.*CLHEP::cm, 400.*CLHEP::cm, 400.*CLHEP::cm};
		
        //		G4double ABSL[NUMENTRIES]  = { 0.35*CLHEP::cm, 0.35*CLHEP::cm, 0.35*CLHEP::cm};
        
		mat_mt = new G4MaterialPropertiesTable();
		mat_mt->AddProperty("FASTCOMPONENT", Energy, SCINT, NUMENTRIES);
		mat_mt->AddProperty("SLOWCOMPONENT", Energy, SCINT, NUMENTRIES);
		mat_mt->AddProperty("RINDEX",        Energy, RIND,  NUMENTRIES);
		mat_mt->AddProperty("ABSLENGTH",     Energy, ABSL,  NUMENTRIES); // mean free path before being absorbed
		mat_mt->AddConstProperty("SCINTILLATIONYIELD",20000./CLHEP::MeV);
        //		mat_mt->AddConstProperty("SCINTILLATIONYIELD",63000./CLHEP::MeV);
		mat_mt->AddConstProperty("FASTSCINTILLATIONRISETIME",1.*CLHEP::ns);
		mat_mt->AddConstProperty("SLOWSCINTILLATIONRISETIME",1.*CLHEP::ns);
		mat_mt->AddConstProperty("FASTTIMECONSTANT",20.*CLHEP::ns);
		mat_mt->AddConstProperty("SLOWTIMECONSTANT",20.*CLHEP::ns);
		
		mat_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
		mat_mt->AddConstProperty("YIELDRATIO",1.0); //relative strenght of the fast against total
		
		// Set the Birks Constant for the LXe scintillator
		// mat->GetIonisation()->SetBirksConstant(0.126*CLHEP::mm/CLHEP::MeV);
	}
	if ( which_properties == G4String("NaI") ) {
		
		const G4int NUMENTRIES = 3; G4double Energy[NUMENTRIES]    = { 2.0*CLHEP::eV , 5*CLHEP::eV, 9*CLHEP::eV };
		
		G4double SCINT[NUMENTRIES] = { 1.0, 1.0, 0.1 };
		G4double RIND[NUMENTRIES]  = { 1.85 , 1.85, 1.85 };
		G4double ABSL[NUMENTRIES]  = { 500.*CLHEP::cm, 500.*CLHEP::cm, 500.*CLHEP::cm };
		
		mat_mt = new G4MaterialPropertiesTable();
		mat_mt->AddProperty("FASTCOMPONENT", Energy, SCINT, NUMENTRIES);
		mat_mt->AddProperty("SLOWCOMPONENT", Energy, SCINT, NUMENTRIES);
		mat_mt->AddProperty("RINDEX",        Energy, RIND,  NUMENTRIES);
		mat_mt->AddProperty("ABSLENGTH",     Energy, ABSL,  NUMENTRIES); // mean free path before being absorbed
		mat_mt->AddConstProperty("SCINTILLATIONYIELD",3800./CLHEP::MeV);
        //		mat_mt->AddConstProperty("SCINTILLATIONYIELD",38000./CLHEP::MeV);
		mat_mt->AddConstProperty("FASTTIMECONSTANT",230.*CLHEP::ns);
		mat_mt->AddConstProperty("SLOWTIMECONSTANT",2300.*CLHEP::ns);
		mat_mt->AddConstProperty("FASTSCINTILLATIONRISETIME",1.*CLHEP::ns);
		mat_mt->AddConstProperty("SLOWSCINTILLATIONRISETIME",1.*CLHEP::ns);
		
		mat_mt->AddConstProperty("RESOLUTIONSCALE",4.0);
		mat_mt->AddConstProperty("YIELDRATIO",0.5);
		
		// Set the Birks Constant for the LXe scintillator
		//mat->GetIonisation()->SetBirksConstant(0.126*CLHEP::mm/CLHEP::MeV);
	}
	//
	if ( which_properties == G4String("Air") || which_properties == G4String("Vacuum") ) {
		
		const G4int NUMENTRIES = 3; G4double Energy[NUMENTRIES]={2.0*CLHEP::eV,7.0*CLHEP::eV,9*CLHEP::eV};
		
		G4double RIND[NUMENTRIES]={1.,1.,1.};
		G4double ABSL[NUMENTRIES]={0.1*CLHEP::mm,0.1*CLHEP::mm,0.1*CLHEP::mm};
		mat_mt = new G4MaterialPropertiesTable();
		mat_mt->AddProperty("ABSLENGTH",Energy,ABSL,NUMENTRIES);
		mat_mt->AddProperty("RINDEX", Energy, RIND,NUMENTRIES);
	}
	//
	if ( which_properties == G4String("Glass") ) {
		
		const G4int NUMENTRIES = 3; G4double Energy[NUMENTRIES]    = { 7.0*CLHEP::eV , 7.07*CLHEP::eV, 7.14*CLHEP::eV };
		G4double RIND[NUMENTRIES]={1.49,1.49,1.49};
		G4double ABSL[NUMENTRIES]={420.*CLHEP::cm,420.*CLHEP::cm,420.*CLHEP::cm};
		mat_mt = new G4MaterialPropertiesTable();
		mat_mt->AddProperty("ABSLENGTH",Energy,ABSL,NUMENTRIES);
		mat_mt->AddProperty("RINDEX",Energy,RIND,NUMENTRIES);
	}
	
	if ( mat_mt ) {
		mat->SetMaterialPropertiesTable(mat_mt);		
	}
}



