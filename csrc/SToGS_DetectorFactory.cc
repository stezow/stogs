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
 
// G4 includes
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4Sphere.hh"

#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"

#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4SubtractionSolid.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"


// Paris includes
#include "SToGS_DetectorFactory.hh"
#include "SToGS_MaterialConsultant.hh"
#include "SToGS_UserActionInitialization.hh"

#include "G4GDMLParser.hh"

using namespace std;
using namespace SToGS;

// mother of all factories
SToGS::DetectorFactory *SToGS::DetectorFactory::pMainFactory = 0x0; DFMessenger *SToGS::DetectorFactory::pMainMessanger = 0x0;
// 
G4int SToGS::DetectorFactory::gCopyNb = 0;

// containing the sub factories ...
std::vector < SToGS::DetectorFactory * > SToGS::DetectorFactory::fSubFactory;


SToGS::DetectorFactory *SToGS::DetectorFactory::theMainFactory()
{
    if ( pMainFactory == 0x0 ) {
        pMainFactory = new SToGS::DetectorFactory();
        pMainMessanger = new DFMessenger();
        
        SToGS::MaterialConsultant::theConsultant(); // just to be sure SToGS materials are built 
    }
    return pMainFactory;
}

SToGS::DetectorFactory *SToGS::DetectorFactory::GetFactory(G4String path)
{
    SToGS::DetectorFactory *sub_factory = 0x0;
    
    for (std::vector<SToGS::DetectorFactory *>::iterator it = fSubFactory.begin() ; it != fSubFactory.end(); ++it){
        G4String sub_name = (*it)->GetFactoryName();
        if ( path.contains(sub_name) ) {
            sub_factory = (*it);
            break;
        }
    }
    
    return sub_factory;
}

void SToGS::DetectorFactory::ChangeSD(G4String opt, G4VPhysicalVolume *top)
{
    G4VPhysicalVolume *l_top = top;
    if ( top == 0x0
        && G4TransportationManager::GetTransportationManager()
        && G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()) {
        l_top = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
    }
    if ( l_top == 0x0 )
        return;

    // collect from top volume all the logical volumes and modify SD accordingly
    std::vector<G4LogicalVolume *> logical_stored; std::vector<G4VPhysicalVolume *> phycical_stored;
    //
    CollectVolumes(l_top, logical_stored, phycical_stored);
    //
    for (size_t i = 0; i < logical_stored.size(); i++) {
//        G4cout << "Changing volume of " << logical_stored[i]->GetName() << G4endl;
        if ( logical_stored[i]->GetSensitiveDetector() )
            SetActive(logical_stored[i], opt);
    }
}


void SToGS::DetectorFactory::MakeStore()
{
    for (std::vector<SToGS::DetectorFactory *>::iterator it = fSubFactory.begin() ; it != fSubFactory.end(); ++it) {
        G4cout << "[+..] BUILDING Store " << (*it)->GetFactoryName() << G4endl;
        (*it)->MakeStore();
        G4cout << "[..+] BUILDING Store " << (*it)->GetFactoryName() << G4endl;        
    }
}

void SToGS::DetectorFactory::MakeInStore(G4String name, G4String version_string)
{
    // check if already in store    
    std::ifstream gdml_in(GetFileName(GetDetName(name,version_string), "gdml").data());
    if ( gdml_in.is_open() ) {
        gdml_in.close();
        G4cout << "[===] Detector " << GetDetName(name,version_string) << " already in store " << G4endl;
        return;
    }
    G4cout << "[+..] Building Detector " << GetDetName(name,version_string) << G4endl;
    G4VPhysicalVolume *theDetector = Make(name, version_string);
    if (theDetector) {
        Store(theDetector);
    }
    G4cout << "[..+] Building Detector " << GetDetName(name,version_string) << G4endl;
}

G4String SToGS::DetectorFactory::GetFileName(G4String detname, G4String extention)
{
    G4String result;

    result = GetFactoryName(); result += detname; result += "."; result += extention;    
    return result;
}

#include <sstream>

void SToGS::DetectorFactory::StreamTouchable(std::ofstream &dmap, G4String dname)
{
    size_t where = dname.find(":");
    if ( where != std::string::npos ) {
        std::string tmp = dname.substr( where );
        dmap << " " << tmp;
    }
}

void SToGS::DetectorFactory::StoreMap(std::ofstream &amap, std::ofstream &dmap,
                               std::vector<G4LogicalVolume *> &logical_stored, std::vector<G4VPhysicalVolume *> &phycical_stored,
                               G4VPhysicalVolume *theDetector /*,
                               const G4String &mother_name */)
{
    G4LogicalVolume *alogical = theDetector->GetLogicalVolume(); G4VPhysicalVolume *aphysical = 0x0;
    if ( alogical == 0x0 )
        return;
    
    // check first if already in the list --> means already written in the amap file
    G4bool in_the_list = false;
    for (size_t i = 0; i < logical_stored.size() ; i++) {
        if (alogical == logical_stored[i] )
            in_the_list = true;
    }
    if ( !in_the_list ) {
        
        G4String sd("S ");
        
        if ( alogical->GetSensitiveDetector() ) {
            sd += (alogical->GetSensitiveDetector())->GetFullPathName();
        }
        else
            sd += "-";
        
        const G4VisAttributes *att = alogical->GetVisAttributes(); G4Color color(0.8,0.8,0.8,0);
        if ( att )
            color = att->GetColor();
        
        amap << alogical->GetName()
            << "\t"
            << "C "
            << setprecision(2)
            << color.GetRed()
            << " "
            << color.GetGreen()
            << " "
            << color.GetBlue()
            << " "
            << color.GetAlpha()
            << "\t"
            << sd
            << " F -"
            << endl;
        
        logical_stored.push_back(alogical);
    }
    
    if ( theDetector->GetLogicalVolume()->GetSensitiveDetector() ) {
        
        dmap << setw(5) << setfill('0')
            << theDetector->GetCopyNo() << "\t"
            << theDetector->GetName() << "\t";
        StreamTouchable(dmap, theDetector->GetName()); dmap << "\t";
 //       dmap << theDetector->GetName() << "\t";
        
        if ( theDetector->GetTranslation().getX() < 0.0 )
            dmap << theDetector->GetTranslation().getX()/CLHEP::cm << " ";
        else
            dmap << "+" << theDetector->GetTranslation().getX()/CLHEP::cm << " ";
        
        if ( theDetector->GetTranslation().getY() < 0.0 )
            dmap << theDetector->GetTranslation().getY()/CLHEP::cm << " ";
        else
            dmap << "+" << theDetector->GetTranslation().getY()/CLHEP::cm << " ";
        
        if ( theDetector->GetTranslation().getZ() < 0.0 )
                dmap << theDetector->GetTranslation().getZ()/CLHEP::cm << " ";
            else
                dmap << "+" << theDetector->GetTranslation().getZ()/CLHEP::cm << " ";
    
        dmap << " cm\t";
        dmap << " "<< std::endl;
        // 
    }

    for (G4int i = 0; i < alogical->GetNoDaughters(); i++) {
        aphysical = alogical->GetDaughter(i);
        StoreMap(amap, dmap, logical_stored, phycical_stored, aphysical /*, theDetector->GetName() */);
    }
}

void SToGS::DetectorFactory::CollectVolumes(G4VPhysicalVolume *theDetector,
                            std::vector<G4LogicalVolume *> &logical_stored, std::vector<G4VPhysicalVolume *> &phycical_stored)
{
    G4LogicalVolume *alogical = theDetector->GetLogicalVolume(); G4VPhysicalVolume *aphysical = 0x0; G4bool in_the_list;
    
    in_the_list = false;
    for (size_t i = 0; i < phycical_stored.size() ; i++) {
        if ( aphysical == phycical_stored[i] ) {
            in_the_list = true;
            break;
        }
    }
    if ( !in_the_list ) {
        phycical_stored.push_back(theDetector);
    }
    in_the_list = false;
    for (size_t i = 0; i < logical_stored.size() ; i++) {
        if ( alogical == logical_stored[i] ) {
            in_the_list = true;
            break;
        }
    }
    if ( !in_the_list ) {
        logical_stored.push_back(alogical);
    }
    
    // recursive call
    for (G4int i = 0; i < alogical->GetNoDaughters(); i++) {
        aphysical = alogical->GetDaughter(i);
        CollectVolumes(aphysical,logical_stored, phycical_stored);
    }
}

#include <sstream>

void SToGS::DetectorFactory::Store(G4VPhysicalVolume *theDetector)
{
    G4String filename = GetFileName(theDetector->GetName(), "gdml");

    std::ifstream gdml_in(filename.data());
    if ( gdml_in.is_open() ) {
        gdml_in.close();
        G4cout << " file already exiting " << filename.data() << G4endl;
        return;
    }
    else {
        // store geometry in gdml
        G4GDMLParser parser;
        parser.Write(filename.data(),theDetector,false);
    }

    // store attribute
    filename = GetFileName(theDetector->GetName(), "amap"); std::ofstream amap(filename.data());
    filename = GetFileName(theDetector->GetName(), "dmap"); std::ofstream dmap(filename.data());
    
    if ( amap.is_open() && dmap.is_open() ) {
        
        std::vector<G4LogicalVolume *> logical_stored; std::vector<G4VPhysicalVolume *> phycical_stored;
        
        for (G4int i = 0; i < theDetector->GetLogicalVolume()->GetNoDaughters(); i++) {
            StoreMap(amap, dmap, logical_stored, phycical_stored, theDetector->GetLogicalVolume()->GetDaughter(i)
                     /*,theDetector->GetName() */);
        }
    }
}

G4String SToGS::DetectorFactory::GetDetName(const G4String &path) const
{
    G4String tmp = path; 
    
    // replace slashes by space
    G4int start = path.last('/');
    tmp.remove(0,start+1);
    start = tmp.last('.');
    tmp.remove(start);

    return tmp;
}

G4String SToGS::DetectorFactory::GetDetName(G4String name, G4String version_string) const
{
    G4String detname = name;
    
    if ( version_string != "" ) {
        detname += "_";
        detname += version_string;
    }
    
    return detname;
}

void SToGS::DetectorFactory::AssemblyRename(G4AssemblyVolume *assembly, const G4String &volume_to_find) const
{
	std::string tmp; 
	
	std::vector < G4VPhysicalVolume* >::iterator vol = assembly->GetVolumesIterator();
	for (size_t i = 0; i < assembly->TotalImprintedVolumes(); i++) {
		G4VPhysicalVolume *an_element = *vol++;
		
		tmp = (an_element)->GetName();
		if ( tmp.find(volume_to_find) != std::string::npos ) { // ok, rename volume
			//
			G4int av, impr, pv, nbread; char *tmpname = new char[tmp.size()];
            
			// necessary otherwise sscanf does not work and set string up to the end
			tmp.replace(tmp.find("_pv"), 3, std::string(" _pv"));
			nbread = sscanf(tmp.data(),"av_%d_impr_%d_%s _pv_%d",&av,&impr,tmpname,&pv);
			if ( nbread == 4 ) {
				G4cout << " Volume " <<  (an_element)->GetName() << ", with copy number " << an_element->GetCopyNo() << " " << (void*)(an_element->GetLogicalVolume()) << G4endl;;
				sprintf(tmpname,"%s_%d",tmpname,--impr);
				(an_element)->SetName(tmpname);
				G4cout << "  |--> becomes " <<  (an_element)->GetName() << G4endl;
			}
			delete [] tmpname;
		}
	}    
}

G4int SToGS::DetectorFactory::DoMap(G4AssemblyVolume *assembly, G4VPhysicalVolume *volume_used_to_built_assembly, G4int copy_number_offset, G4String opt) const
{

    G4cout << "[+] Results of snapshots of " <<  volume_used_to_built_assembly->GetName() << endl;

	std::string tmp; G4int item = 0, nb_sd; std::vector <G4int> items; // for each imprinted volumes, current number
    // from the detector used to make the assembly, collect logical and physical to remap the imprinted volumes
    std::vector<G4LogicalVolume *> logical_stored; std::vector<G4VPhysicalVolume *> phycical_stored; std::vector<G4int > volume_counter;
    //
    CollectVolumes(volume_used_to_built_assembly, logical_stored, phycical_stored);
    //
    
    // count number of SD in primary detector
    volume_counter.resize(phycical_stored.size(),0); items.resize( assembly->GetImprintsCount() , 0 );
    nb_sd = 0;
    for (size_t i = 0; i < logical_stored.size(); i++) {
        if ( logical_stored[i]->GetSensitiveDetector() ) {
            nb_sd++;
        }
    }
    
    std::vector < G4VPhysicalVolume* >::iterator vol = assembly->GetVolumesIterator();
    for (size_t i = 0; i < assembly->TotalImprintedVolumes(); i++) {
        
        G4VPhysicalVolume *an_element = *vol++;
        tmp = (an_element)->GetName();
        
        G4int av, impr, pv, nbread; char *tmpname = new char[tmp.size()];
        
        // necessary otherwise sscanf does not work and set string up to the end
        tmp.replace(tmp.find("_pv"), 3, std::string(" _pv"));
        nbread = sscanf(tmp.data(),"av_%d_impr_%d_%s _pv_%d",&av,&impr,tmpname,&pv);

        impr--;

        // look for the physical volume this imprinted volumes refered to, extract the name and add :impr
        G4int keep_j = 0;
        for (size_t j = 1; j < phycical_stored.size(); j++) {
            if ( an_element->GetLogicalVolume() == phycical_stored[j]->GetLogicalVolume() ) {
                
                if ( impr == volume_counter[j] ) {
                    keep_j = j;
                    break;
                }
            }
        }
        ostringstream hname, tmpsize; tmpsize << assembly->GetImprintsCount();
        
        hname << phycical_stored[keep_j]->GetName() << ":" << setw(tmpsize.str().size()) << setfill('0') << impr;
       
        (an_element)->SetName(hname.str());
        volume_counter[keep_j]++;
        //
        // set new copy number depending of option and if SD
        if ( an_element->GetLogicalVolume()->GetSensitiveDetector() ) {
            if ( opt == "->" )
               (an_element)->SetCopyNo(copy_number_offset + item );
            else
                (an_element)->SetCopyNo(copy_number_offset + items[impr]*assembly->GetImprintsCount() + impr );
            
            items[impr]++;
            item++;
        }
        else
            (an_element)->SetCopyNo(-1);
        //
        G4cout << " + " << (an_element)->GetName() << " with copy number " << (an_element)->GetCopyNo() << G4endl;

        
        delete [] tmpname;
    }
    
    G4cout << "[-] Results of snapshots of " <<  volume_used_to_built_assembly->GetName()  << endl;
    
    return item;
    
}

void SToGS::DetectorFactory::SetActive(G4LogicalVolume *pv, const G4String &opt)
{
    //
	std::string lopt = opt, tmp; G4bool ok = false;
//    if ( lopt.find('^') == std::string::npos ) { // only one string, force to finish it with separtor ^
 //       lopt += "^";
  //  }
    lopt += "^";
    
	std::vector < std::string > all_opt; std::vector < size_t > pos; size_t start = 0;
	//
	for (size_t i = 0; i < lopt.size(); i++) { // search for one sequence i.e a string stopped by ^
        
		if ( lopt[i] == '^' ) {
			
			tmp = lopt.substr(start,i-start); start = i+1;
            
			pos.clear();
			for (size_t j = 0; j < tmp.size(); j++) { // check the sequence contains 3 times the character |
				if ( tmp[j] == '|' ) {
					pos.push_back ( j );
				}
			}
            if ( pos.size() == 2 ) {
                tmp += "|";
                pos.push_back( tmp.size() - 1 ); // in case color is not defined just set pos[2] to end of string
            }
			if ( pos.size() == 3 ) { // valid expression
				
				G4String volumename = tmp.substr(0,pos[0]);
				G4String matname    = tmp.substr(pos[0]+1,pos[1]-pos[0]-1);
				G4String sdname     = tmp.substr(pos[1]+1,pos[2]-pos[1]-1);
				G4String color      = tmp.substr(pos[2]+1);
				
                //				G4cout << " dec (" << tmp << ") " << volumename << " " << matname << " " << sdname << " " << color << G4endl;
                //				G4cout << " from " << pv->GetName() << " " <<  pv->GetLogicalVolume()->GetMaterial()->GetName() << G4endl;
				
				G4int nb_ok = 0;
				// compare volume name, matname,
				if ( volumename.find("*") == std::string::npos ) { // a given type of volume
					if ( pv->GetName().find(volumename) != std::string::npos )
						nb_ok++;
				}
				else nb_ok++;
				
				if ( matname.find("*") == std::string::npos ) {// a given type of material
					if ( pv->GetMaterial()->GetName().find(matname) != std::string::npos )
						nb_ok++;
				}
				else nb_ok++;
                // decode color. if not given set to default only if color not yet defined
                G4VisAttributes *visatt = 0x0; G4int valid_col = 0; G4double r = 0.5, g = 0.5, b = 0.5, a = 0.5;
                for (size_t k = 0; k < color.size(); k++) {
                    if ( color[k] == ';' ) {
                        valid_col++;
                        color[k] = ' ';
                    }
                }
                if ( valid_col == 3 ) {
                    istringstream col(color);
                    col >> r >> g >> b >> a;
                    
                    visatt = new G4VisAttributes( G4Colour(r, g, b, a) );
                    pv->SetVisAttributes( visatt );
                }
                if ( pv->GetSensitiveDetector() == 0x0 ) {
                    visatt = new G4VisAttributes( G4Colour(r, g, b, a) );
                    pv->SetVisAttributes( visatt );
                }
                
				// check sensitivity. extention list of name separted by ; for scorer ?
				G4VSensitiveDetector *sensdet = GetSD(sdname);
				
				if ( nb_ok == 2 ) { // should be added to the list of sensitive detectors.
                    if ( sensdet ) {
     //                   if ( G4SDManager::GetSDMpointer()->FindSensitiveDetector(sdname) ) {
                        ok = true;
                    }
				}
				if ( ok ) {
					pv->SetSensitiveDetector(sensdet);
                    G4cout << pv->GetName() << " is now sensitive of type " << sensdet->GetFullPathName() << G4endl;
				}
				if ( ok )
                    break;
			}
		}
	}
}

G4VSensitiveDetector *SToGS::DetectorFactory::GetSD(G4String sdname, const char sd_type)
{
    G4VSensitiveDetector *sensdet = 0x0;
    
    // test first if found in SDManager
    sensdet = G4SDManager::GetSDMpointer()->FindSensitiveDetector(sdname,false);
    if ( sensdet == 0x0 ) {
        if ( sd_type == 'S' ) {
            if ( sdname == "/SToGS/Calo")
               sensdet = SToGS::UserActionInitialization::GetCaloSD(sdname);
            if ( sdname == "/SToGS/Track")
                sensdet = SToGS::UserActionInitialization::GetTrackerSD(sdname);
        }
        if ( sd_type == 's' ) {
        }
        
        if (sensdet) {
            G4cout << "[i] " << sensdet->GetFullPathName() << " has been added to the list of SD detector" << G4endl;
        }
    }
    return sensdet;
}

G4VPhysicalVolume *SToGS::DetectorFactory::Import(G4String gdmlfile, G4String detname, const G4String &opt_amap, const G4String &opt_dmap)
{
    G4VPhysicalVolume *theDetector = 0x0; // G4String detname, fullname;

    // check if file exist
    std::ifstream isgdml_in(gdmlfile.data());
    if ( isgdml_in.is_open() ) {
        isgdml_in.close();
    }
    else {
        return theDetector;
    }
    
    cout << "[+..] Importing " << gdmlfile << " in " << GetFactoryName() << detname << endl;
    
    G4GDMLParser parser;
    parser.Read(gdmlfile,false);
    
    theDetector = parser.GetWorldVolume();
    if ( theDetector ) {
        
        // change top names
        theDetector->SetName(detname);
        theDetector->GetLogicalVolume()->SetName(detname);
        
        // collect logical to eventual set sensitivity using opt
        std::vector<G4LogicalVolume *> logical_stored; std::vector<G4VPhysicalVolume *> phycical_stored;
        //
        CollectVolumes(theDetector, logical_stored, phycical_stored);
 
        for (size_t i = 0; i < logical_stored.size(); i++) {
            SetActive(logical_stored[i], opt_amap);
        }
        for (size_t i = 0; i < phycical_stored.size(); i++) {

            if ( phycical_stored[i]->GetLogicalVolume()->GetSensitiveDetector() ) {
                phycical_stored[i]->SetCopyNo( SToGS::DetectorFactory::AddGCopyNb() );
            }
            else
                phycical_stored[i]->SetCopyNo(-1);
            //
            if ( opt_dmap == "T" && phycical_stored[i] != theDetector &&  phycical_stored[i]->GetCopyNo() >=0 ) {
                ostringstream tmp;
                tmp << phycical_stored[i]->GetName() << ":" << phycical_stored[i]->GetCopyNo();
                
                phycical_stored[i]->SetName(tmp.str());
            }
        }
        Store(theDetector);
    }
    cout << "[..+] Importing " << gdmlfile << " in store " << GetFactoryName() << " with name " <<  detname << endl;
    
    return theDetector;
}

/*
G4VPhysicalVolume *SToGS::DetectorFactory::Get(G4String basename)
{
    G4VPhysicalVolume *theDetector = GetGeometry(basename);
    if ( theDetector ) {
        GetAttributes(basename);
    }
    return theDetector;
}
 */

G4VPhysicalVolume *SToGS::DetectorFactory::Get(G4String basename, G4bool is_full)
{
    G4VPhysicalVolume *theDetector = 0x0; G4String detname, fullname;
    
    fullname  = basename;
    fullname += ".gdml";
    
    // check is already loaded
    for (size_t i = 0; i < fLoadedPhysical.size(); i++) {
        if ( fLoadedPhysical[i].first == fullname ) {
            theDetector = fLoadedPhysical[i].second;
            break;
        }
    }
    if ( theDetector ) {
        return theDetector;
    }
    
    std::ifstream isgdml_in(fullname.data());
    if ( isgdml_in.is_open() ) {
        isgdml_in.close();
    }
    else {
        return theDetector;
    }

    cout << "[+..] Loading from store " << basename << endl;
    
    G4GDMLParser parser;
    parser.Read(fullname,false);
    //
    detname = GetDetName(fullname); theDetector = parser.GetWorldVolume(detname); // load geometry from gdml
    if ( theDetector == 0x0 ) {
        // should add a warning
        return theDetector;
    }

    std::pair < G4String, G4VPhysicalVolume *> p(fullname,theDetector); // add the new loaded detector to the list 
    fLoadedPhysical.push_back(p);
    
    if ( is_full ) {
        GetAttributes(basename);
    }
    
    cout << "[..+] Loading from store " << basename << endl;

    return theDetector;
}

void SToGS::DetectorFactory::GetAttributes(G4String basename)
{
    G4VPhysicalVolume *theDetector; G4String detname, fullname;
    
    theDetector = Get(basename,false);
    if ( theDetector == 0x0 )
        return;
    
    cout << "[+..] Loading Attributes from store " << basename << endl;
    
    // now apply map
    fullname  = basename; fullname += ".dmap"; std::ifstream dmap(fullname.data());
    fullname  = basename; fullname += ".amap"; std::ifstream amap(fullname.data());
    
    std::vector<G4LogicalVolume *> logical_stored; std::vector<G4VPhysicalVolume *> phycical_stored; G4String aline; G4bool has_done_amap = false;
    //
    CollectVolumes(theDetector, logical_stored, phycical_stored);
    //
    if ( amap.is_open() ) {
        
        // logicals
        getline(amap,aline);
        while ( amap.good() ) {
            
            istringstream decode(aline);
            G4String vname, key_sd, sd, key_color, touchable; G4double r,g,b,a;
            
            decode >> vname
            >> key_color
            >> r
            >> g
            >> b
            >> a
            >> key_sd
            >> sd
            ;
            
            for (size_t i = 0; i < logical_stored.size(); i++) {
                if ( vname == logical_stored[i]->GetName() ) {
                    cout << " Change attributes [colors] of " << vname << " to " << G4Color(r,g,b,a) << endl;
                    //
                    G4VisAttributes visatt(G4Color(r,g,b,a));
                    if ( key_color == "C" )
                        visatt.SetVisibility(true);
                    else
                        visatt.SetVisibility(false);
                    
                    logical_stored[i]->SetVisAttributes(visatt);
                    
                    // SD
                    if ( sd != "-" ) { // means a sensitive detector so look for it and assign to the volume
                        G4VSensitiveDetector * sensdet = GetSD(sd);
                        if (sensdet) {
                            cout << " Change attribute [sensitive] of " << vname << " to " << sd << endl;
                            logical_stored[i]->SetSensitiveDetector(sensdet);
                        }
                    }
                    // others ... field etc ...
                    
                    break;
                }
            }
            
            getline(amap,aline);
        }
    }
    else { // do an amap file with all the sensitive detectors
        fullname  = basename; fullname += ".amap"; std::ofstream amap_(fullname.data());
        
        if ( amap_.is_open() ) {
            has_done_amap = true;
            for (size_t i = 0; i < logical_stored.size(); i++) {
                
                //               if ( logical_stored[i] == theDetector->GetLogicalVolume() ) { // mother is not a
                //                   continue;
                //               }
                
                amap_ << logical_stored[i]->GetName()
                << "\t"
                << "C "
                << setprecision(2)
                << "0.5"
                << " "
                << "0.5"
                << " "
                << "0.5"
                << " "
                << "0.5"
                << "\t"
                << " S -"
                << " F -"
                << endl;
            }
        }
    }
    
    //
    if ( dmap.is_open() ) {
        
        for (size_t i = 0; i < phycical_stored.size(); i++) {
            phycical_stored[i]->SetCopyNo(-1);
        }
        
        // physicals
        getline(dmap,aline);
        while ( dmap.good() ) {
            
            G4String pname, unit, touchable; G4double x,y,z; G4int uid;
            istringstream decode(aline);
            
            decode >> uid
            >> pname
            >> touchable
            >> x
            >> y
            >> z
            >> unit
            ;
            //               cout << " Change Copy number of " << pname << phycical_stored.size() << endl;
            
            for (size_t i = 0; i < phycical_stored.size(); i++) {
                if ( pname == phycical_stored[i]->GetName() ) {
                    cout << " Load Copy number of " << pname << " : " << uid << endl;
                    //
                    phycical_stored[i]->SetCopyNo(uid);
                    break;
                }
            }
            //
            getline(dmap,aline);
        }
    }
    else {
        for (size_t i = 0; i < phycical_stored.size(); i++) {
            phycical_stored[i]->SetCopyNo(-1);
        }
        
        if ( !has_done_amap ) { // already an amap, change copy # of active volumes and store results in dmap
            
            fullname  = basename; fullname += ".dmap"; std::ofstream dmap_(fullname.data());
            
            for (size_t i = 0; i < phycical_stored.size(); i++) {
                if ( phycical_stored[i]->GetLogicalVolume()->GetSensitiveDetector() ) {
                    cout << " Change Name and Copy number of "
                    << phycical_stored[i]->GetName() << " -> "
                    << phycical_stored[i]->GetName() << ": " << SToGS::DetectorFactory::GetGCopyNb() << " "
                    << SToGS::DetectorFactory::GetGCopyNb() << endl;
                    //
                    ostringstream tmp; tmp << phycical_stored[i]->GetName(); tmp << ":"; tmp << SToGS::DetectorFactory::GetGCopyNb();
                    
                    phycical_stored[i]->SetName(tmp.str());
                    phycical_stored[i]->SetCopyNo( SToGS::DetectorFactory::AddGCopyNb() );
                    
                    if ( dmap_.is_open() ) {
                        dmap_ << setw(5) << setfill('0')
                        << phycical_stored[i]->GetCopyNo() << "\t"
                        << phycical_stored[i]->GetName() << "\t";
                        StreamTouchable(dmap_, phycical_stored[i]->GetName()); dmap_ << "\t";
                        
                        if ( phycical_stored[i]->GetTranslation().getX() < 0.0 )
                            dmap_ << phycical_stored[i]->GetTranslation().getX()/CLHEP::cm << " ";
                        else
                            dmap_ << "+" << phycical_stored[i]->GetTranslation().getX()/CLHEP::cm << " ";
                        
                        if ( phycical_stored[i]->GetTranslation().getY() < 0.0 )
                            dmap_ << phycical_stored[i]->GetTranslation().getY()/CLHEP::cm << " ";
                        else
                            dmap_ << "+" << phycical_stored[i]->GetTranslation().getY()/CLHEP::cm << " ";
                        
                        if ( phycical_stored[i]->GetTranslation().getZ() < 0.0 )
                            dmap_ << phycical_stored[i]->GetTranslation().getZ()/CLHEP::cm << " ";
                        else
                            dmap_ << "+" << phycical_stored[i]->GetTranslation().getZ()/CLHEP::cm << " ";
                        
                        dmap_ << " cm\t";
                        dmap_ << " "<< std::endl;
                    }
                }
            }
        }
    }
    
    cout << "[..+] Loading Attributes from store " << basename << endl;
}


G4bool SToGS::DetectorFactory::Set(G4String basename, G4VPhysicalVolume *mother, G4int copy_number_offset, G4ThreeVector *T, G4RotationMatrix *R)
{
    G4VPhysicalVolume *thefullDetector = 0x0, *subdetector ; G4LogicalVolume *volume_to_move;
    
    thefullDetector = Get(basename); // load from factory the detector, remove the envelop and move the content into mother with all its attribute
    if ( thefullDetector == 0x0 ) {
        return false;
    }
    else { volume_to_move = thefullDetector->GetLogicalVolume(); /* thefullDetector->SetLogicalVolume(0x0); */ }
    
    // set copy number as known by full detector
    for (G4int i = 0; i < volume_to_move->GetNoDaughters(); i++) {
        subdetector = volume_to_move->GetDaughter(i);
//        cout << " + " << subdetector->GetName() << " " << subdetector->GetCopyNo() << endl;
        
        G4ThreeVector *T_ = new G4ThreeVector(); G4RotationMatrix *R_ = new G4RotationMatrix();
        
        (*T_) = subdetector->GetObjectTranslation();
        if ( subdetector->GetRotation() ) {
            (*R_) = (*subdetector->GetRotation());
        }
        if ( T )
            (*T_) += (*T);
        if ( R )
            (*R_) = (*R) * (*R_);
        
        new G4PVPlacement(R_,(*T_),
                            subdetector->GetLogicalVolume(),
                            subdetector->GetName(),
                            mother->GetLogicalVolume(),
                            false,
                            copy_number_offset + subdetector->GetCopyNo());
    }
    
    return true;
}

G4AssemblyVolume *SToGS::DetectorFactory::GetAssembly(G4String basename)
{
    G4AssemblyVolume *theAssembly = 0x0;  G4VPhysicalVolume *theDetector = 0x0, *subdetector;  G4String detname, fullname;
    
    fullname  = basename;
    fullname += ".gdml";
    
    // check is already loaded
    for (size_t i = 0; i < fLoadedAssembly.size(); i++) {
        if ( fLoadedAssembly[i].first == fullname ) {
            theAssembly = fLoadedAssembly[i].second;
            break;
        }
    }
    
    if ( theAssembly == 0x0 ) {
        theDetector = Get(basename);
        if (theDetector == 0x0 ) {
            return 0x0;
        }
        theAssembly = new G4AssemblyVolume();
        //
        std::pair < G4String, G4AssemblyVolume *> p(fullname,theAssembly); // add the new assembly to the list of assembly
        fLoadedAssembly.push_back(p);
        
        // rebuilt the detector as an assembly
        for (G4int i = 0; i < theDetector->GetLogicalVolume()->GetNoDaughters(); i++) {
            subdetector = theDetector->GetLogicalVolume()->GetDaughter(i);
  //          cout << " + " << subdetector->GetName() << " " << subdetector->GetCopyNo() << endl;
            
            G4Transform3D Tr = G4Transform3D(subdetector->GetObjectRotationValue(),subdetector->GetObjectTranslation());
            theAssembly->AddPlacedVolume( subdetector->GetLogicalVolume(), Tr );
        }

    }
    return theAssembly;
}

void SToGS::DetectorFactory::Clean()
{    
    if ( this == DetectorFactory::theMainFactory() ) {
        for (size_t i = 0; i < fSubFactory.size(); i++) {
            fSubFactory[i]->Clean();
        }
    }
    else {
        // clean the two inner collection of this factory and call G4 store manager to clean physical and volumes
        for (size_t i = 0; i < fLoadedPhysical.size(); i++) {
            fLoadedPhysical[i].second = 0x0;
        }
        fLoadedPhysical.resize(0);
        for (size_t i = 0; i < fLoadedAssembly.size(); i++) {
            delete fLoadedAssembly[i].second; fLoadedAssembly[i].second = 0x0;
        }
        fLoadedAssembly.resize(0);
        //
        G4LogicalVolumeStore::Clean(); G4PhysicalVolumeStore::Clean();
    }
}

G4VPhysicalVolume *SToGS::DetectorFactory::MakeVCR(G4String name, G4double HalfX, G4double HalfY, G4double HalfZ, G4int copy_number)
{
    G4VPhysicalVolume *theRoom = 0x0; G4LogicalVolume *logicRoom; G4Box *room;
    G4RotationMatrix R;
	G4ThreeVector T;
    G4Transform3D Tr;
    
    // use a physical as a container to describe the detector
	room= new G4Box(name,HalfX,HalfY,HalfZ);
	logicRoom= new G4LogicalVolume(room, MaterialConsultant::theConsultant()->FindOrBuildMaterial("AIR"), name, 0, 0, 0);
	
	logicRoom->SetVisAttributes(G4VisAttributes::Invisible); // hide the world
	
	//  Must place the World Physical volume unrotated at (0,0,0).
	theRoom = new G4PVPlacement(0,         // no rotation
                                G4ThreeVector(), // at (0,0,0)
                                logicRoom,      // its logical volume
                                name,      // its name
                                0,               // its mother  volume
                                false,           // no boolean operations
                                copy_number);              // copy number
    
    return theRoom;
}


G4VPhysicalVolume * SToGS::DetectorFactory::MakeAnArrayFromFactory(G4String input_file)
{
    G4VPhysicalVolume *theDetector = 0x0; std::vector< std::pair< G4AssemblyVolume *,G4VPhysicalVolume * > > all_assembly;

    // open the file containing the array to be built
    std::ifstream g4map(input_file.data());
    //
    if ( ! g4map.is_open() ) {
        return theDetector;
    }
    // first read the
    
    G4String key, detector_name, subdetector_name, aline, unit1, what, option_copy_number; G4double X, Y, Z, r_value;
    
    // logicals
    getline(g4map,aline);
    while ( g4map.good() ) {
        
        istringstream decode(aline); decode.clear();
        
        decode >> key ;
        
        if ( key == "#" ) {
            getline(g4map,aline);
            continue;
        }
        // change global copy number
        if ( key == "@" ) {
            G4int offset; decode >> offset >> option_copy_number;

            SToGS::DetectorFactory::SetGCopyNb(offset);
            
            getline(g4map,aline);
            continue;
        }
        // imports an xml file, adds attributes and save
        if ( key == "i" ) {
            G4String gdmlfile, fullfactoryname, option_amap("*|*|*|0.5;0.5;0.5;0.5"), option_dmap("");
            decode >> gdmlfile >> fullfactoryname >> option_amap >> option_dmap;
            
            SToGS::DetectorFactory *factory = SToGS::DetectorFactory::theMainFactory()->GetFactory(fullfactoryname);
            if ( factory ) {
                detector_name = SToGS::DetectorFactory::theMainFactory()->GetDetName(fullfactoryname);
                theDetector = factory->Import(gdmlfile, detector_name, option_amap, option_dmap);
            }
            getline(g4map,aline);
            continue;
        }
        
        decode >> subdetector_name ;
        
        if ( key == "w" && decode.good() ) { // this is the name of the new setup/detector
                        
            // use a physical as a container to describe the detector ... should be enough !
            theDetector = SToGS::DetectorFactory::theMainFactory()->Get(subdetector_name.data());
            if ( theDetector == 0x0 ) {
                
                decode >> X >> Y >> Z >> unit1 ;
                G4double do_unit =
                    G4UnitDefinition::GetValueOf(unit1);

                theDetector = SToGS::DetectorFactory::MakeVCR(subdetector_name,do_unit*X,do_unit*Y,do_unit*Z,-1);
            }
        }
        if ( key == "+" && decode.good() && theDetector ) { // simply load the detector and add it
            
            G4RotationMatrix *R = new G4RotationMatrix();

            decode >> X >> Y >> Z >> unit1 ;
            G4double do_unit =
                G4UnitDefinition::GetValueOf(unit1);
            G4ThreeVector T(do_unit*X,do_unit*Y,do_unit*Z);

            while ( !decode.eof() ) {
                decode >> what;
                if ( what == "Rx" ) {
                    decode >> r_value;
                    R->rotateX(r_value*CLHEP::deg);
                }
                if ( what == "Ry" ) {
                    decode >> r_value;
                    R->rotateY(r_value*CLHEP::deg);
                }
                if ( what == "Rz" ) {
                    decode >> r_value;
                    R->rotateZ(r_value*CLHEP::deg);
                }
                if ( what == "Rt" ) {
                    T = (*R) * T;
                }
                what = "";
            }
            
            Set(subdetector_name, theDetector, SToGS::DetectorFactory::GetGCopyNb(), &T, R);
        }
        if ( key == "*" && decode.good() && theDetector ) { // this is a detector going to be replicated in space using the assembly mechanism
            
            G4VPhysicalVolume *detector = SToGS::DetectorFactory::theMainFactory()->Get(subdetector_name.data());
            G4AssemblyVolume  *assembly = SToGS::DetectorFactory::theMainFactory()->GetAssembly(subdetector_name.data());
            
            if ( assembly ) {
                
                decode >> X >> Y >> Z >> unit1 ;
                G4double do_unit =
                    G4UnitDefinition::GetValueOf(unit1);
                G4ThreeVector T(do_unit*X,do_unit*Y,do_unit*Z);
                
                // keep list of assemble required to build this new detector/setup
                G4bool in = false;
                for (size_t i = 0; i < all_assembly.size(); i++ ) {
                    if ( assembly == all_assembly[i].first ) {
                        in = true;
                        break;
                    }
                }
                if ( !in ) {
                    all_assembly.push_back( std::pair< G4AssemblyVolume *, G4VPhysicalVolume *>(assembly,detector) );
                }
                
                G4RotationMatrix *R =new G4RotationMatrix();
                
                while ( !decode.eof() ) {
                    decode >> what;
                    if ( what == "Rx" ) {
                        decode >> r_value;
                        R->rotateX(r_value*CLHEP::deg);
                    }
                    if ( what == "Ry" ) {
                        decode >> r_value;
                        R->rotateY(r_value*CLHEP::deg);
                    }
                    if ( what == "Rz" ) {
                        decode >> r_value;
                        R->rotateZ(r_value*CLHEP::deg);
                    }
                    if ( what == "Rt" ) {
                        T = (*R) * T;
                    }
                    what = "";
                }
                
                assembly->MakeImprint( theDetector->GetLogicalVolume(), T, R );
            }
        }

        getline(g4map,aline);
    }
    
    // treated @ the end i.e. after already built detector even if not in order in the input file
    for (size_t i = 0; i < all_assembly.size(); i++ ) {
        G4int nb_active = DoMap(all_assembly[i].first,all_assembly[i].second,SToGS::DetectorFactory::GetGCopyNb());
        SToGS::DetectorFactory::SetGCopyNb(SToGS::DetectorFactory::GetGCopyNb() + nb_active);
    }
    
    return theDetector;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

DFMessenger::DFMessenger():
    theDirectory(0x0),
    cmdToChangeSD(0x0),
    cmdToSaveInFactory(0x0)
{
	theDirectory = new G4UIdirectory("/DetectoryFactory/");
	theDirectory->SetGuidance("To manage the DetectoryFactory");
	
	cmdToChangeSD = new G4UIcmdWithAString("/DetectoryFactory/changeSD", this);
	cmdToChangeSD->SetGuidance("Change Sensitivity of Sensitive detectors ");
	cmdToChangeSD->SetGuidance("Requires one string having the form XXX:YYY:ZZZ  ");
	cmdToChangeSD->SetGuidance(" \t XXX is the volume name, could be * for all");
	cmdToChangeSD->SetGuidance(" \t YYY is the mat name, could be *  ");
	cmdToChangeSD->SetGuidance(" \t ZZZ is the name of the new SD  ");
   
	cmdToChangeSD->AvailableForStates(G4State_PreInit,G4State_Idle);
    
	cmdToSaveInFactory = new G4UIcmdWithAString("/DetectoryFactory/Save", this);
	cmdToSaveInFactory->SetGuidance("To save the current setup in MyStore ");
	cmdToSaveInFactory->SetGuidance("Requires one string which is the name of the setup to be stored");
	cmdToSaveInFactory->AvailableForStates(G4State_PreInit,G4State_Idle);
}

DFMessenger::~DFMessenger()
{
	delete theDirectory;
    // delete exportGEOMCmd;  delete exportMAPCmd;
}

void DFMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if( command == cmdToChangeSD ) {
        std::string modifier(newValue);
        if ( modifier[modifier.length()-1] != '^' )
            modifier += '^';
        SToGS::DetectorFactory::ChangeSD(modifier);
	}
    
	if( command == cmdToSaveInFactory ) {
        SToGS::DetectorFactory *where_to_store = SToGS::DetectorFactory::GetFactory("DetectorFactory/MyStore/");
        if ( where_to_store ) {
            G4VPhysicalVolume *world = 0x0;
            if ( G4TransportationManager::GetTransportationManager()
                 && G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()) {
                world = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
            }
            if ( world ) {
                world->SetName(newValue); world->GetLogicalVolume()->SetName(newValue);
                where_to_store->Store(world);
            }
            else {
                G4cout << " Cannot save world in MyStore ! " << G4endl;
            }
        }
	}
    
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








