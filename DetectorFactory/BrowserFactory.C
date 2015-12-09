

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"

//! From the Full geometry, extract the logical and physical [nodes] volumes
void CollectVolumes(TGeoNode *theDetector, TObjArray &logical_stored, TObjArray &phycical_stored)
{
    
    TGeoVolume *alogical = theDetector->GetVolume(); TGeoNode *aphysical = 0x0; Bool_t in_the_list;
    
    in_the_list = false;
    for (size_t i = 0; i < phycical_stored.GetEntries() ; i++) {
        if ( aphysical == phycical_stored[i] ) {
            in_the_list = true;
            break;
        }
    }
    if ( !in_the_list ) {
        phycical_stored.Add(theDetector);
    }
    in_the_list = false;
    for (size_t i = 0; i < logical_stored.GetEntries() ; i++) {
        if ( alogical == logical_stored[i] ) {
            in_the_list = true;
            break;
        }
    }
    if ( !in_the_list ) {
        logical_stored.Add(alogical);
    }
    
    // recursive call
    for (size_t i = 0; i < alogical->GetNdaughters(); i++) {
        aphysical = alogical->GetNode(i);
        CollectVolumes(aphysical,logical_stored, phycical_stored);
    }
}

//! Read the detector in the factory and show it in open GL window
TGeoVolume *ShowDetector(const char *basename = "DetectorFactory/Scintillators/CParisPW_2", Option_t *opt_draw = "ogl")
{
    TString detname, fullname; TObjArray logical_stored; TObjArray phycical_stored;
    
    fullname  = basename;
    fullname += ".gdml";
    cout << fullname << endl;
    
    TGeoManager::Import(fullname.Data()); world = gGeoManager->GetTopVolume();
    //
    if ( world == 0x0 ) {
        std::cout << "\n ***** The detector " << basename << " has not been found in the factory ! ***** " << std::endl;
        return 0x0;
    }
    
    // collect logical and physical (nodes) volumes
    CollectVolumes(gGeoManager->GetTopNode(),logical_stored,phycical_stored);
    
    fullname  = basename; fullname += ".amap"; std::ifstream amap(fullname.Data());
    if ( amap.is_open() ) {
        
        std::string touchable, aline;
        
        // logicals
        getline(amap,aline);
        while ( amap.good() ) {
            
            istringstream decode(aline); TString vname, key_sd, sd, key_color; Float_t r,g,b,a;
            
            decode >> vname
                >> key_color
                >> r
                >> g
                >> b
                >> a
                >> key_sd
                >> sd
                ;
            
            for (size_t i = 0; i < logical_stored.GetEntries(); i++) {
                TString tmp = logical_stored[i]->GetName();
                if ( tmp == vname ) {
                    cout << " **** Set Attributes to: " << tmp << " using " << aline << endl;

                    TGeoVolume *alogical = (TGeoVolume *)logical_stored[i];
                    
                    alogical->SetFillColor( TColor::GetColor(r,g,b) );
                    alogical->SetLineColor( TColor::GetColor(r,g,b) );
                    alogical->SetTransparency(100*(1-a));
                }
            }
            
            getline(amap,aline);
        }
     }

    if ( opt_draw != "-" ) {
        world->Draw(opt_draw);
    }
    gGeoManager->SetVisOption(0);

    return world;
}

//! Read the detector in the factory and show it in Eve
void ShowInEve(const char *basename = "DetectorFactory/Scintillators/CParisPW_2")
{
    TEveManager::Create();
    
    TGeoVolume *world = ShowDetector(basename,"-");
    TEveGeoTopNode* inn = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
    gEve->AddGlobalElement(inn);

    gEve->FullRedraw3D(kTRUE);
    
    // EClipType not exported to CINT (see TGLUtil.h):
    // 0 - no clip, 1 - clip plane, 2 - clip box
    TGLViewer *v = gEve->GetDefaultGLViewer();
    v->GetClipSet()->SetClipType(1);
    v->RefreshPadEditor(v);
    
    v->CurrentCamera().RotateRad(-.7, 0.5);
    v->DoDraw();
}

//! rename the three files
void RenameSToGSDetector(const char *from, const char *to)
{
    TString t1,t2;
    
    t1  = from;
    t1 += ".gdml";
    t2  = to;
    t2 += ".gdml";
    gSystem->Rename(t1.Data(),t2.Data());
    t1  = from;
    t1 += ".amap";
    t2  = to;
    t2 += ".amap";
    gSystem->Rename(t1.Data(),t2.Data());
    t1  = from;
    t1 += ".dmap";
    t2  = to;
    t2 += ".dmap";
    gSystem->Rename(t1.Data(),t2.Data());
}



/*
//! Read the detector in the factory and show it in open GL window
TGeoVolume *GetDetector(const char *basename = "DetectorFactory/Scintillators/CParisPW_2")
{
    TString detname, fullname; TObjArray logical_stored; TObjArray phycical_stored;
    
    fullname  = basename;
    fullname += ".gdml";
    cout << fullname << endl;
    
    if ( gGeoManager == 0x0 ) {
        new TGeoManager("GDMLImport", "Geometry imported from GDML");
        gGeoManager->SetVisOption(0);
    }
    
    TString cmd = TString::Format("TGDMLParse::StartGDML(\"%s\")", fullname.Data());
    TGeoVolume* world = (TGeoVolume*)gROOT->ProcessLineFast(cmd);
    if ( world == 0x0 ) {
        return 0x0;
    }
    
    // collect logical and physical (nodes) volumes
    CollectVolumes(world->GetNode(0),logical_stored,phycical_stored);
    
    fullname  = basename; fullname += ".amap"; std::ifstream amap(fullname.Data());
    if ( amap.is_open() ) {
        
        std::string touchable, aline;
        
        // logicals
        getline(amap,aline);
        while ( amap.good() ) {
            
            istringstream decode(aline); TString vname, key_sd, sd, key_color; Float_t r,g,b,a;
            
            decode >> vname
            >> key_color
            >> r
            >> g
            >> b
            >> a
            >> key_sd
            >> sd
            ;
            
            for (size_t i = 0; i < logical_stored.GetEntries(); i++) {
                TString tmp = logical_stored[i]->GetName();
                if ( tmp == vname ) {
                    cout << " **** Set Attributes to: " << tmp << " using " << aline << endl;
                    
                    TGeoVolume *alogical = (TGeoVolume *)logical_stored[i];
                    
                    alogical->SetFillColor( TColor::GetColor(r,g,b) );
                    alogical->SetLineColor( TColor::GetColor(r,g,b) );
                    alogical->SetTransparency(100*(1-a));
                }
            }
            
            getline(amap,aline);
        }
    }
    
    return world;
}

void ShowInEveT()
{
    TEveManager::Create();
    
    TGeoVolume *det1 = GetDetector("DetectorFactory/SemiConductors/Ge/ATC");
    TEveGeoTopNode* inn1 = new TEveGeoTopNode(gGeoManager, det1->GetNode(0));
    gEve->AddGlobalElement(inn1);
  
    TGeoVolume *det2 = GetDetector("DetectorFactory/SemiConductors/Ge/EXOCLOVER_A-AC");
    TEveGeoTopNode* inn2 = new TEveGeoTopNode(gGeoManager, det2->GetNode(0));
    gEve->AddGlobalElement(inn2);
    
    gEve->FullRedraw3D(kTRUE);
    
    // EClipType not exported to CINT (see TGLUtil.h):
    // 0 - no clip, 1 - clip plane, 2 - clip box
    TGLViewer *v = gEve->GetDefaultGLViewer();
    v->GetClipSet()->SetClipType(1);
    v->RefreshPadEditor(v);
    
    v->CurrentCamera().RotateRad(-.7, 0.5);
    v->DoDraw();
}
*/

/*
void ExportToROOT(const Char_t *rfilename = "DetectorFactory.root")
{
    // first canvas to draw
    TFile *f = new TFile(rfilename,"RECREATE");
    f->cd();
    
    // loop and draw
    TString d_name = allDetectors[0], tmp; Int_t i = 1;
    do {
        
        tmp  = d_name;
        tmp += ".gdml";
        
        if ( gSystem->AccessPathName(tmp.Data()) ) { // does not exist
            d_name = allDetectors[i++];
            cout << " **** This file does not exist " << d_name << endl;
            continue;
        }

        TGeoVolume *det = LoadOneDetector(d_name.Data());
        if ( det ) {
            det->Write();
        }
        
        // next in factory
        d_name = allDetectors[i++];
        
        delete gGeoManager;
        
    } while ( d_name != "end" );
    
}

 */



