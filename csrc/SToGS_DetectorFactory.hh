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

#ifndef SToGS_DetectorFactory_h
#define SToGS_DetectorFactory_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"

// std includes 
#include <iostream>
#include <sstream>
#include <fstream>

class G4Box;
class G4Sphere;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4AssemblyVolume;
class G4VSensitiveDetector;

//! SToGS namespace to protect SToGS classes
namespace SToGS {
    
    class DFMessenger;

    //! Base Factory. This is a container of sub-factories
    /*!
     All operation should be done through this interface using theMainFactory (singleton).
     The main factory is in charge to find the correct sub factory and sent the work to it.
     
     New factories are autmoatically added with a global static
     */
    class DetectorFactory
    {
    /// DATA MEMBERS ///
    private:
        //! the main factory [NOT YET deleted ... should be done one day through a singleton destructor]
        static DetectorFactory *pMainFactory;
    private:
        //! use for messanger action
        static DFMessenger *pMainMessanger;
    private:
        //! list of all registered factories
        /*!
            A list of all registered factory is used.
         */
        static std::vector < DetectorFactory * > fSubFactory;
    private:
        //! global used to set a unique number to each sensitive volume
        static G4int gCopyNb;
        
    protected:
        //! Name of the factory, also used for the name of the directory containing the detectors handled by the factory
        G4String fPath;
    protected:
        //! in case a gdml file has already by loaded, it allows to get the assembly quickly
        std::vector < std::pair < G4String, G4AssemblyVolume *> > fLoadedAssembly;
        //! in case it has been loaded, keep a trace 
        std::vector < std::pair < G4String, G4VPhysicalVolume *> > fLoadedPhysical;
        
    /// METHODS ///
    private:
        //! Register a factory to the main factory
        /*!
         Warning ! nothing is checked so several factories could be registered with the same name ....
         In case of search in the list of subfactories only the first one could be returned ...
         
         So better to have one unique name per factory
         */
        void Register(DetectorFactory *df)
        {
            fSubFactory.push_back(df);
        }
        
    protected:
        //! used only for the main, does not register this in the factory
        DetectorFactory():
            fPath("DetectorFactory/"),
            fLoadedAssembly(0),
            fLoadedPhysical(0)
            {;}
        //! used for inheriting classes, it registers this in the factory
        DetectorFactory(G4String path):
            fPath(path),
            fLoadedAssembly(1000),
            fLoadedPhysical(1000)
        {
            fPath.prepend(theMainFactory()->GetFactoryName().data());
            theMainFactory()->Register(this) ;
        }
        
    public:
// Main factory facilities
        //! global copy number, used to set copy number once building detector. Not protected for MT since in principle init of geom done without MT
        static G4int AddGCopyNb()
        {
            G4int val = gCopyNb++;
            return val;
        }
        static G4int SetGCopyNb(G4int val)
        {
            G4int current = gCopyNb;
            gCopyNb = val;
            return current;
        }
        static G4int GetGCopyNb()
        {
            return gCopyNb;
        }
        //! to get the main factory
        static DetectorFactory *theMainFactory();
        
        //! get the specific factory for a full name to a file
        static DetectorFactory *GetFactory(G4String fullname);
        
        //! make a cubic world
        static G4VPhysicalVolume *MakeVCR(
                            G4String name,
                            G4double HalfX = 5.0*CLHEP::m, G4double HalfY = 5.0*CLHEP::m, G4double HalfZ = 5.0*CLHEP::m,
                            G4int copy = -1);
        
        //! Set to the given logical volume some attributs (sensitive, color) depending of the given option. used @ the import level
        /*!
         option has the following format : Volumename|material|sensitivedetectorname|color \n
         with color red;green;blue;alpha
         */
        static void SetActive(G4LogicalVolume *pv, const G4String &option = "*|*|*|0.5;0.5;0.5;0.5");
        
        //! Change some SD detector into an other type of SD. Opt is used eventually to filter on the name of the volume and its material name
        /*!
         if top is 0x0, tries with world
         */
        static void ChangeSD(G4String opt, G4VPhysicalVolume *top = 0x0);
        //! Get a particular SD. S means a SD while s is for Scorers
        static G4VSensitiveDetector * GetSD(G4String opt, const char sd_type = 'S');
        
        //! From a top volume, it collects into collection all logical and physical volumes
        static void CollectVolumes(G4VPhysicalVolume *theDetector,
                                   std::vector<G4LogicalVolume *> &logical_stored, std::vector<G4VPhysicalVolume *> &phycical_stored);
        
        protected:
// Methods to be overwritten in particular factory
        //! Should be implemented in any sub factory. It built (C++) a detector and return it
        virtual G4VPhysicalVolume * Make(G4String /* name */, G4String /* version_string */)
        {
            return 0x0;
        }
        //! recursively called in the tree struture
        void StoreMap(std::ofstream &amap, std::ofstream &dmap,
                      std::vector<G4LogicalVolume *> &logical_stored,
                      std::vector<G4VPhysicalVolume *> &phycical_stored,
                      G4VPhysicalVolume *theDetector /*,
                                                      const G4String &mother_name */);
        
        //! extrac from the name of the detector the touchable part
        void StreamTouchable(std::ofstream &dmap, G4String);
        
        //! import a gdml file (extention to ascii ?) into factory. See Active for opt_amap option
        /*!
         if opt_dmap == "T" touchable is done i.e. the physical name is changed to include the copy numnber PNAME:COPY#
         */
        G4VPhysicalVolume *Import(G4String gdmlfile, G4String detname, const G4String &opt_amap, const G4String &opt_dmap);
        
        
    public:
        virtual ~DetectorFactory()
        { /* what to delete ? ... the main factory .. */ ;}
        
        G4String GetFactoryName()
        {
            return fPath;
        }
        
    public:
        //! Get detector name from full path to the gdml file
        /*!
         It removes the path and the gdml extention. This name is supposed to be equal to the one of the detector inside the gdml file
         */
        G4String GetDetName(const G4String &path) const;
        G4String GetDetName(G4String name, G4String version) const;
        
        //! compute name of the file to store/read information
        G4String GetFileName(G4String detname, G4String extention);
        
        //! search for a detector in DetectorFactory
        /*!
         Get calls GetGeometry + GetAttributes. Separated since geant4.10 requires separation of geometrical part and attributes for Multi-Threading
         
         it has to start with the subfactory name and then the name of the file. It adds autmoatically .gdml
         Ex:
         
            Scintillators/ParisPW_2 --> DetectorFactory/Scintillators/ParisPW_2.gdml
         */
        virtual G4VPhysicalVolume *Get(G4String basename, G4bool is_full = true);
        G4AssemblyVolume  *GetAssembly(G4String basename);
        //! Read the amap file and apply atributes to the detector. if not found, it creates a deefault one from the sensitive detector founds
        virtual void GetAttributes(G4String basename);
        
        //! to be used once a detector is fully contructed to simply place it the the world
        G4bool Set(G4String basename, G4VPhysicalVolume *world, G4int copy_number_offset, G4ThreeVector *T = 0x0, G4RotationMatrix *R = 0x0);

        //! clear factory i.e. in memory collections of physicals and assemblies
        void Clean();

    protected:
        //! build in store a particular detector from its names and version. i.e. call th Make method of the sub factory and then Store
        /*!
            The most important method is Make which, depending of the name and version calls the correct C++ methods with the right parameters.
            This methods just check is has not yet been made in store and if not call Make then Store
         */
        void MakeInStore(G4String name, G4String version_string);
        
    public:
        //! build the default store i.e. all the detectors. Here It calls MakeStore for all registered sub factories
        virtual void MakeStore();
        
        //! Store in the sub-factory the given physical volume
        virtual void Store(G4VPhysicalVolume *);
        
        //! rename physical volume produced by the assembly
        void AssemblyRename(G4AssemblyVolume *assembly, const G4String &volume_to_find) const;
        
        //! remap all the physical volumes produced by the Imprint method.
        /*
         It loops over the produced physical volumes,
         changes their names by adding the copy number to them taking into account i.e. mame -> name:#:
         
         In principle then at the end an active volume should have a name XXX:A:B:C: ...
         where XXX is a generic string and A B C are numbers corresponding to touchable history as introduced by Geant4
         
         As well, if a volume is an active detector, its copy number is modified starting from copy_number_offset
         
         It returns the nmuber of active volume counted in the full assembly
         
         opt = down means in the copy number is set going into the detector tree struture from top
         if imprint is the snapshot number, volumeid is the position of an active volume in one snapshot and tot_imprint and tot_volumeid the total number of them,
         then \n
         - if opt = down, a volume as for copy id  copy_number_offset + volumeid + tot_volumeid
         - otherwise copy_number_offset + volumeid + tot_imprint
         
         if a detector is composed of several elements in case of down they are going to have consecutive copy numbers
         while with other option,
         
         */
        G4int DoMap(G4AssemblyVolume *assembly, G4VPhysicalVolume *volume_used_to_built_assembly, G4int copy_number_offset = 0, G4String opt = "->") const;
        
        //! built an array from the factory using the given input file
        virtual G4VPhysicalVolume * MakeAnArrayFromFactory(G4String input_file);
    };
} // SToGS Namespace

    
#include "G4UImessenger.hh"
    
class G4UIdirectory;
class G4UIcmdWithAString;
    
//! SToGS namespace to protect SToGS classes
namespace SToGS {
    //! the Messenger for the Detector Factory
    /*!
     */
    class DFMessenger: public G4UImessenger
    {
    public:
        DFMessenger();
        virtual ~DFMessenger();
        
        void SetNewValue(G4UIcommand*, G4String);
        
    private:
        //! current director in /Paris/DetectorFactory
        G4UIdirectory *theDirectory;
        
        //! the commands
        G4UIcmdWithAString *cmdToChangeSD;
        G4UIcmdWithAString *cmdToSaveInFactory;
    };
    
} // SToGS Namespace

#endif


