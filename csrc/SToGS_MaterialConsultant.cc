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
#include "G4NistManager.hh"

#include "globals.hh"

SToGS::MaterialConsultant *SToGS::MaterialConsultant::theMaterialConsultant = 0x0;

SToGS::MaterialConsultant::MaterialConsultant() :
    theElements(),
    theMaterials(),
    fIsSToGDBSSearchedFirst(true)
{
}

G4Element *SToGS::MaterialConsultant::BuildElement(G4String name_element)
{
    G4double z, a; G4String stogsname = name_element, name, symbol;
    
    if ( name_element.index("SToGS_") != 0 ) {
        stogsname = "SToGS_" + name_element;
    }
    if ( !stogsname.contains("_El") ) {
        stogsname += "_El";
    }
    
    // first search the already created elements
    G4Element *element = 0x0;
    for (size_t i = 0; i < theElements.size(); i++) {
        if ( theElements[i]->GetName() == stogsname ) {
            return theElements[i];
        }
    }
    size_t size_before = theElements.size();
    
    // not found so add it
    if ( stogsname == "SToGS_H_El" ) {
        a=1.01*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_H_El",symbol="H",z=1.,a));
    }
    if ( stogsname == "SToGS_Deuterium_El" ) {
        a=2.01*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Deuterium_El",symbol="D",z=1.,a));
    }
    if ( stogsname == "SToGS_D_El" ) {
        a=2.01*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_D_El",symbol="D",z=1.,a));
    }
    if ( stogsname == "SToGS_He_El" ) {
        a=4.*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_He_El",symbol="He",z=2.,a));
    }
    if ( stogsname == "SToGS_Li_El" ) {
        a=6.94*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Li_El",symbol="Li",z=3.,a));
    }
    if ( stogsname == "SToGS_Be_El" ) {
        a=9.01*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Be_El",symbol="Be",z=4.,a));
    }
    if ( stogsname == "SToGS_C_El" ) {
        a=12.01*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_C_El",symbol="C",z=6.,a));
    }
    if ( stogsname == "SToGS_N_El" ) {
        a=14.00674*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_N_El",symbol="N",z=7.,a));
    }
    if ( stogsname == "SToGS_O_El" ) {
        a=15.9994*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_O_El",symbol="O",z=8.,a));
    }
    if ( stogsname == "SToGS_Ne_El" ) {
        a=20.18*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Ne_El",symbol="Ne",z=10.,a));
    }
    if ( stogsname == "SToGS_Na_El" ) {
        a=22.99*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Na_El",symbol="Na",z=11.,a));
    }
    if ( stogsname == "SToGS_Mg_El" ) {
        a=24.305*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Mg_El",symbol="Mg",z=12.,a));
    }
    if ( stogsname == "SToGS_Al_El" ) {
        a=26.98154*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Al_El",symbol="Al",z=13.,a));
    }
    if ( stogsname == "SToGS_Si_El" ) {
        a=28.085*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Si_El",symbol="Si",z=14.,a));
    }
    if ( stogsname == "SToGS_Cl_El" ) {
        a=35.453*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Cl_El",symbol="Cl",z=17.,a));
    }
    if ( stogsname == "SToGS_Fe_El" ) {
        a=55.850*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Fe_El",symbol="Fe",z=26.,a));
    }
    if ( stogsname == "SToGS_Cu_El" ) {
        a=63.546*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Cu_El",symbol="Cu",z=29.,a));
    }
    if ( stogsname == "SToGS_Ba_El" ) {
        a=137.327*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Ba_El",symbol="Ba",z=56.,a));
    }
    if ( stogsname == "SToGS_F_El" ) {
        a=18.9984032*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_F_El",symbol="F",z=9.,a));
    }
    if ( stogsname == "SToGS_Ge_El" ) {
        // Germanium isotopes
        G4Isotope* Ge70 = new G4Isotope(name="SToGS_Ge70", 32, 70, 69.9242*CLHEP::g/CLHEP::mole);
        G4Isotope* Ge72 = new G4Isotope(name="SToGS_Ge72", 32, 72, 71.9221*CLHEP::g/CLHEP::mole);
        G4Isotope* Ge73 = new G4Isotope(name="SToGS_Ge73", 32, 73, 72.9235*CLHEP::g/CLHEP::mole);
        G4Isotope* Ge74 = new G4Isotope(name="SToGS_Ge74", 32, 74, 73.9212*CLHEP::g/CLHEP::mole);
        G4Isotope* Ge76 = new G4Isotope(name="SToGS_Ge76", 32, 76, 75.9214*CLHEP::g/CLHEP::mole);
        // germanium defined via its isotopes
        G4Element* elGe = new G4Element(name="SToGS_Ge_El",symbol="Ge", 5);
        elGe->AddIsotope(Ge70, 0.2123);
        elGe->AddIsotope(Ge72, 0.2766);
        elGe->AddIsotope(Ge73, 0.0773);
        elGe->AddIsotope(Ge74, 0.3594);
        elGe->AddIsotope(Ge76, 0.0744);
        theElements.push_back(elGe);
    }
    if ( stogsname == "SToGS_Br_El" ) {
        a=79.904*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Br_El",symbol="Br",z=35.,a));
    }
    if ( stogsname == "SToGS_I_El" ) {
        a=126.90477*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_I_El",symbol="I",z=53.,a));
    }
    if ( stogsname == "SToGS_Cs_El" ) {
        a=132.90545*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Cs_El",symbol="Cs",z=55.,a));
    }
    if ( stogsname == "SToGS_La_El" ) {
        a=138.9055*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_La_El",symbol="La",z=57.,a));
    }
    if ( stogsname == "SToGS_W_El" ) {
        a=183.85*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_W_El",symbol="W",z=74.,a));
    }
    if ( stogsname == "SToGS_Pb_El" ) {
        a=207.19*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Pb_El",symbol="Pb",z=82.,a));
    }
    if ( stogsname == "SToGS_Bi_El" ) {
        a=208.98038*CLHEP::g/CLHEP::mole;
        theElements.push_back(new G4Element(name="SToGS_Bi_El",symbol="Bi",z=83.,a));
    }
    // to be checked ..
    if ( stogsname == "SToGS_Calcium_El" ) {
        // a=40.08*CLHEP::g/CLHEP::mole;
        // theElements.push_back(new G4Element(name="SToGS_Calcium",symbol="Ca",z=20.,a));
    }
    if ( stogsname == "SToGS_Manganese_El" ) {
        // a=54.938*CLHEP::g/CLHEP::mole;
        // theElements.push_back(new G4Element(name="SToGS_Manganese",symbol="Mn",z=25.,a));
    }
    if ( stogsname == "SToGS_Uranium_El" ) {
        // a=238.03*CLHEP::g/CLHEP::mole;
        // theElements.push_back(new G4Element(name="SToGS_Uranium",symbol="U",z=92.,a));
    }
    if ( theElements.size() != size_before ) {
        element = theElements.back();
    }
    return element;
}

G4Material *SToGS::MaterialConsultant::BuildSimpleMaterial(G4String name_material, G4double density)
{
    G4String name_el = name_material + "_El"; G4Element *el = BuildElement(name_el); G4Material *mat = 0x0;
    if ( el == 0x0 )
        return 0x0;
    
    G4int nel = 1; mat = new G4Material(name_material,density,nel);
    mat->AddElement(el, 1);
    
    return mat;
}

G4Material *SToGS::MaterialConsultant::BuildMaterial(G4String name_material)
{
    G4double density; G4int nel,natoms; G4String name; G4String stogsname = name_material, symbol;
    
    if ( name_material.index("SToGS_") != 0 ) {
        stogsname = "SToGS_" + name_material;
    }
    
    // first search the already created elements
    G4Material *mat = 0x0;
    for (size_t i = 0; i < theMaterials.size(); i++) {
        if ( theMaterials[i]->GetName() == stogsname ) {
            return theMaterials[i];
        }
    }
    size_t size_before = theMaterials.size();
    
    //-------------------
    // simple materials
    //-------------------
    if ( stogsname == "SToGS_Al" ) {
        //density = 2.7*CLHEP::g/CLHEP::cm3;
        //a = 26.98*CLHEP::g/CLHEP::mole;
        //theMaterials.push_back(new G4Material(name="SToGS_Al",z=13.,a,density));
        
        mat = BuildSimpleMaterial(name_material,density = 2.6989*CLHEP::g/CLHEP::cm3);
        
    }
    if ( stogsname == "SToGS_Fe" ) {
        //density = 7.87*CLHEP::g/CLHEP::cm3;
        //a = 55.85*CLHEP::g/CLHEP::mole;
        //theMaterials.push_back(new G4Material(name="SToGS_Fe",z=26.,a,density));
        
        mat = BuildSimpleMaterial(name_material,density = 7.874*CLHEP::g/CLHEP::cm3);

    }
    if ( stogsname == "SToGS_Ni" ) {
        //density = 8.96*CLHEP::g/CLHEP::cm3;
        //a = 58.69*CLHEP::g/CLHEP::mole;
        //theMaterials.push_back(new G4Material(name="SToGS_Ni",z=28.,a,density));
        
        mat = BuildSimpleMaterial(name_material,density = 8.902*CLHEP::g/CLHEP::cm3);

    }
    if ( stogsname == "SToGS_Cu" ) {
        //density = 8.96*CLHEP::g/CLHEP::cm3;
        //a = 63.54*CLHEP::g/CLHEP::mole;
        //theMaterials.push_back(new G4Material(name="SToGS_Cu",z=29.,a,density));
        
        mat = BuildSimpleMaterial(name_material,density = 8.96*CLHEP::g/CLHEP::cm3);

    }
    if ( stogsname == "SToGS_W" ) {
        //density = 19.3*CLHEP::g/CLHEP::cm3;
        //a = 183.85*CLHEP::g/CLHEP::mole;
        //theMaterials.push_back(new G4Material(name="SToGS_W",z=74.,a,density));
        
        mat = BuildSimpleMaterial(name_material,density = 19.3*CLHEP::g/CLHEP::cm3);

    }
    if ( stogsname == "SToGS_Pb" ) {
        //density = 11.35*CLHEP::g/CLHEP::cm3;
        //a = 207.19*CLHEP::g/CLHEP::mole;
        //theMaterials.push_back(new G4Material(name="SToGS_Pb",z=82.,a,density));
        
        density = 11.35*CLHEP::mg/CLHEP::cm3; nel = 1; mat = new G4Material(name="SToGS_Pb",density,nel);
        mat->AddElement(BuildElement("SToGS_Pb"), 1);
        theMaterials.push_back(mat);
        
        mat = BuildSimpleMaterial(name_material,density = 11.35*CLHEP::mg/CLHEP::cm3);

    }

    //------------
    // Registering of all the needed element. Don't forget to do it if you use one
    // element to build a composite material
    //------------
    if ( stogsname == "SToGS_AIR" ) {
        density = 1.290*CLHEP::mg/CLHEP::cm3; nel = 2; mat = new G4Material(name="SToGS_AIR",density,nel);
        mat->AddElement(BuildElement("SToGS_N"), 70);
        mat->AddElement(BuildElement("SToGS_O"), 30);
        theMaterials.push_back(mat);
    }
    if ( stogsname == "SToGS_Ge" ) {
      density = 5.323*CLHEP::g/CLHEP::cm3; nel = 1; mat = new G4Material(name="SToGS_Ge",density,nel);
       mat->AddElement(BuildElement("SToGS_Ge"), 100*CLHEP::perCent);
      theMaterials.push_back(mat);
    }
    if ( stogsname == "SToGS_NaI" ) {
      density = 3.67*CLHEP::g/CLHEP::cm3, nel = 2; mat = new G4Material(name="SToGS_NaI",density,nel);
      mat->AddElement(BuildElement("SToGS_Na"), natoms = 1);
      mat->AddElement(BuildElement("SToGS_I"), natoms = 1);
      theMaterials.push_back(mat);
    }
    if ( stogsname == "SToGS_LaCl3" ) {
      density = 3.79*CLHEP::g/CLHEP::cm3, nel = 2; mat = new G4Material(name="SToGS_LaCl3",density,nel);
        mat->AddElement(BuildElement("SToGS_La"), natoms = 1);
        mat->AddElement(BuildElement("SToGS_Cl"), natoms = 3);
        theMaterials.push_back(mat);
    }
    if ( stogsname == "SToGS_CsI" ) {
        density  = 4.51*CLHEP::g/CLHEP::cm3, nel = 2; mat = new G4Material(name="SToGS_CsI", density, nel);
        mat->AddElement(BuildElement("SToGS_Cs"), natoms = 1);
        mat->AddElement(BuildElement("SToGS_I"), natoms = 1);
        theMaterials.push_back(mat);
    }
    if ( stogsname == "SToGS_LaBr3" ) {
        density = 5.29*CLHEP::g/CLHEP::cm3, nel = 2; mat = new G4Material(name="SToGS_LaBr3",density,nel);
        mat->AddElement(BuildElement("SToGS_La"), natoms = 1);
        mat->AddElement(BuildElement("SToGS_Br"), natoms = 3);
        theMaterials.push_back(mat);
    }
    if ( stogsname == "SToGS_BGO" ) {
        density = 7.13*CLHEP::g/CLHEP::cm3, nel = 3; mat = new G4Material(name="SToGS_BGO", density, nel);
        mat->AddElement(BuildElement("SToGS_Bi"), natoms = 4);
        mat->AddElement(BuildElement("SToGS_Ge"), natoms = 3);
        mat->AddElement(BuildElement("SToGS_O"), natoms = 12);
        theMaterials.push_back(mat);
    }
    if ( stogsname == "SToGS_BaF2" ) {
        density = 4.89*CLHEP::g/CLHEP::cm3, nel = 2; mat = new G4Material(name="SToGS_BaF2", density, nel);
        mat->AddElement(BuildElement("SToGS_Ba"), natoms = 1);
        mat->AddElement(BuildElement("SToGS_F"), natoms = 2);
        theMaterials.push_back(mat);
    }
    if ( stogsname == "SToGS_Quartz" ) {
        density = 2.648*CLHEP::g/CLHEP::cm3, nel = 2; mat = new G4Material(name="SToGS_Quartz", density, nel);
        mat->AddElement(BuildElement("SToGS_Si"), natoms = 1);
        mat->AddElement(BuildElement("SToGS_O") , natoms = 2);
        theMaterials.push_back(mat);
    }
    if ( name_material == "SToGS_Vacuum_Ge" ) {
         density = 2.376e-15*CLHEP::g/CLHEP::cm3, nel = 2; mat = new G4Material(name="SToGS_Vacuum_Ge", density, nel);
         mat->AddElement(BuildElement("SToGS_N"), 70);
         mat->AddElement(BuildElement("SToGS_O"), 30);
         theMaterials.push_back(mat);
     }
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
    
    if ( theMaterials.size() != size_before ) {
        mat = theMaterials.back();
    }
    return mat;
}


SToGS::MaterialConsultant *SToGS::MaterialConsultant::theConsultant()
{
    if (theMaterialConsultant == 0x0) {
        theMaterialConsultant = new SToGS::MaterialConsultant();
    }
    return theMaterialConsultant;
}

G4Material *SToGS::MaterialConsultant::FindOrBuildMaterial(G4String what)
{
    // to be sure the inner table is built
    G4Material* material = G4Material::GetMaterial(what,false); G4String lwhat = what;
    
    if ( material == 0x0 ) { // try to look for in Nist or by adding G4_ or StoGS_ to material
        if ( fIsSToGDBSSearchedFirst ) {
            material = SToGS::MaterialConsultant::theConsultant()->BuildMaterial(lwhat);
            if ( material == 0x0 ) {
                lwhat  = "SToGS_";
                lwhat += what;
                material = SToGS::MaterialConsultant::theConsultant()->BuildMaterial(lwhat);
                if ( material == 0x0 ) {
                    lwhat = what;
                    material = G4NistManager::Instance()->FindOrBuildMaterial(lwhat);
                    if ( material == 0x0 ) {
                        lwhat  = "G4_";
                        lwhat += what;
                        material = G4NistManager::Instance()->FindOrBuildMaterial(lwhat);
                    }
                }
            }
        }
        else {
            material = G4NistManager::Instance()->FindOrBuildMaterial(lwhat);
            if ( material == 0x0 ) {
                lwhat  = "G4_";
                lwhat += what;
                material = G4NistManager::Instance()->FindOrBuildMaterial(lwhat);
                if ( material == 0x0 ) {
                    lwhat = what;
                    SToGS::MaterialConsultant::theConsultant()->BuildMaterial(lwhat);
                    if ( material == 0x0 ) {
                        lwhat  = "SToGS_";
                        lwhat += what;
                        material = SToGS::MaterialConsultant::theConsultant()->BuildMaterial(lwhat);
                    }
                }
            }
        }
    }
    if ( material && what != lwhat )
        G4cout << "[W] in SToGS::MaterialConsultant::FindOrBuildMaterial " << what << " replaced by " << lwhat << G4endl;
    if ( material == 0x0 )
        G4cout << "[E] in SToGS::MaterialConsultant::FindOrBuildMaterial, material " << what << " not found " << G4endl;
    
    return material;
}


/*
G4Material *SToGS::MaterialConsultant::GetMaterial(G4String what) const
{
    // to be sure the inner table is built
    SToGS::MaterialConsultant::theConsultant(); G4Material* material = G4Material::GetMaterial(what,false); G4String lwhat = what;
    
    if ( material == 0x0 ) { // try NIST first else add
        material = G4NistManager::Instance()->FindOrBuildMaterial(lwhat);
        
        if ( material == 0x0 ) {
            G4String tmp = "SToGS_";
            tmp+=what;
            material = G4Material::GetMaterial(tmp);
        }
    }
    
    return material;
}
*/

G4Element *SToGS::MaterialConsultant::GetElement(G4String what) const
{
    G4Element *element = G4Element::GetElement(what);
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
	if ( which_properties == G4String("AIR") || which_properties == G4String("Vacuum") ) {
		
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



