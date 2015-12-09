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

#ifndef SToGS_G4_MaterialConsultant_h
#define SToGS_G4_MaterialConsultant_h 1

#include "G4Material.hh"
#include "G4Element.hh"
#include "globals.hh"

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! Ask this class to get a material (kind of material factory)
    /*!
     This class is the central point to define all the materials that will be needed to build
     your detector. This is a kind of material factory. To initialised the factory:
     
     \code MaterialConsultant *materialconsultant = MaterialConsultant::theConsultant();
     
     Once the factory is initialised, you can ask for a specific material with:
     
     \code G4Material *mymat = materialconsultant->GetMaterial("AIR")
     
     if this material does not exist, a WARNING is issued and mymat is equal to NULL.
     
     The full list of materials is constructed in the constructor of this class. So if you need to add a
     new material, you should add it there.
     */
    class MaterialConsultant
    {
    private:
        //! Private to make it a singleton
        /*!
         The list of all needed materials is defined there. This is the place to add a new one, just follow
         the way it is done for already defined ones.
         */
        MaterialConsultant();
        
    public:
        //! No need to delete the created materials/elements, Geant makes it !
        virtual ~MaterialConsultant()
        {
            theMaterialConsultant = 0x0;
        }
    public:
        //! It returns the unique consultant (singleton)
        static MaterialConsultant *theConsultant();
        
    public:
        //! attached to the given material some optical properties
        void SetOpticalProperties(G4Material *, G4String which_properties);
        
        //! It looks into the SToGS DATABASE and returns a G4Element. The name to be given is XX and the element is SToGS_XX_El
        G4Element  *GetElement(G4String)  const;
        //! It looks into the different DATABASE and returns a G4Material
        G4Material *FindOrBuildMaterial(G4String what);
        
    private:
        static SToGS::MaterialConsultant *theMaterialConsultant;
        
    protected:
        //! methods to be completed in case one would like to define differently some elements/materials. Cenventions 
        G4Element  *BuildElement(G4String name);
        //! methods to be completed in case one would like to define differently some elements/materials
        G4Material *BuildMaterial(G4String name);
       
    private:
        std::vector<G4Element*>  theElements;    // save pointers here to avoid
        std::vector<G4Material*> theMaterials;   // warnings of unused variables
        
        //! if true (default) try and look in SToGS DB for new objects in FindOrBuidMaterial. Otherwise is is NIST used first
        G4bool fIsSToGDBSSearchedFirst;
        
    protected:
        //! this is used to build a simple material using a single element
        // ! this is a work around a possible big in gdml ... 
        G4Material *BuildSimpleMaterial(G4String name, G4double d);
    };
} // SToGS Namespace

#endif   /* MaterialConsultant.hh */






