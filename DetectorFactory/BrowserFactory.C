

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"

//void CollectVolumes(TGeoNode *theDetector, std::vector<TGeoVolume*> &logical_stored, std::vector<TGeoNode *> &phycical_stored)
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

Bool_t ShowDetector(const char *basename = "DetectorFactory/Scintillators/CParisPW_2", Option_t *opt_draw = "ogl")
{
    TString detname, fullname; TObjArray logical_stored; TObjArray phycical_stored;
    
    fullname  = basename;
    fullname += ".gdml";
    cout << fullname << endl;
    
    TGeoManager::Import(fullname.Data());
    world = gGeoManager->GetTopVolume();
    
    // collect logical and physical (nodes) volumes
    CollectVolumes(gGeoManager->GetTopNode(),logical_stored,phycical_stored);
    
    fullname  = basename; fullname += ".amap"; std::ifstream amap(fullname.Data());
    if ( amap.is_open() ) {
        
        std::string touchable, aline;
        
        // logicals
        getline(amap,aline);
        while ( amap.good() ) {
            
            istringstream decode(aline);
            TString vname, key_sd, sd, key_color; Float_t r,g,b,a;
            
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
                    cout << " **** Set Attributes to " << tmp << " " << aline << endl;

                    TGeoVolume *alogical = (TGeoVolume *)logical_stored[i];
                    
                  //  alogical->SetFillColor( TColor::GetColor(r,g,b) );
                    alogical->SetLineColor( TColor::GetColor(r,g,b) );
                    alogical->SetTransparency(100*(1-a));
                }
            }
            
            getline(amap,aline);
        }
     }
    
    gGeoManager->CloseGeometry();
    world->Draw(opt_draw);
    
    return true;
}

/*
void ShowInEve(const char *basename = "DetectorFactory/Scintillators/CParisPW_2")
{
    TEveManager::Create();
    

    gGeoManager = gEve->GetGeometry(basename);

    
    gGeoManager->DefaultColors();
    
    //TEveGeoTopNode* inn = new TEveGeoTopNode(gGeoManager, world);
    //gEve->AddGlobalElement(inn);
    
    
    
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



