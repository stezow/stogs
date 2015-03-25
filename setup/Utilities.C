

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

#include "TSystemDirectory.h"
#include "TList.h"
#include "TRegexp.h"

// for a given root file (name of it), it check whether or not it contains the stequence Thread_whichthread
bool guess_thread(const std::string &rfile, size_t whichthread)
{
    std::ostringstream o;
    o.clear();
    o << "_Thread" << std::setfill('0') << std::setw(2) << whichthread ;
    
    if ( rfile.find(o.str()) != std::string::npos) {
        return true;
    }
    return false;
}


//! function to split a list of ROOT file into configuration files used by toROOTGPS
void DOtoROOTGPS_gene( std::vector<std::string> &list_of_files, const char *treename, const char *confbasename, size_t nb_thread )
{
    size_t size_of_all_gene = 0;
    if ( nb_thread == 0 )
        size_of_all_gene = 1;
    else
        size_of_all_gene = nb_thread;
    //
    std::vector< std::ofstream * > all_gene(size_of_all_gene); // stack of configuration files

    if ( nb_thread == 0 ) { // means no Mthread ... i.e. prepare a single file for non-MT Geant4
        
        std::string fullfilename = confbasename; fullfilename += ".gene";
        all_gene[0] = new std::ofstream( fullfilename.data() );
        if ( !all_gene[0]->is_open() ) {
            std::cout << "Cannot open file " << fullfilename << std::endl;
            return;
        }
        else {
            (*all_gene[0]) << "# it starts with the name of the tree inside the root file " << std::endl;
            (*all_gene[0]) << "tree_name " << treename << std::endl;
            (*all_gene[0]) << "# then + to add a new file followed by the file name " << std::endl;

        }
    }
    else {
        for (size_t i = 0; i < nb_thread; i++) { // open all conf files and add tree definition
            std::ostringstream fullfilename;
            fullfilename.clear();
            fullfilename << confbasename << "_Thread" << std::setfill('0') << std::setw(2) << i << ".gene";
            
            all_gene[i] = new std::ofstream( fullfilename.str().data() );
            if ( !all_gene[i]->is_open() ) {
                std::cout << "Cannot open file " << fullfilename.str().data() << std::endl;
                return;
            }
            else {
                (*all_gene[i]) << "# it starts with the name of the tree inside the root file " << std::endl;
                (*all_gene[i]) << "tree_name " << treename << std::endl;
                (*all_gene[i]) << "# then + to add a new file followed by the file name " << std::endl;

            }
        }
    }
    
    // now adds the root files to the .gene. Try to see if the root file has already been labelled by a thread. Otherwise dispach on all one by one
    size_t modulo_on_thread = 0;
    for (size_t i = 0; i < list_of_files.size(); i++) {
        bool is_in = false;
        for (size_t which_thread = 0; which_thread < all_gene.size(); which_thread++) {
            if ( guess_thread(list_of_files[i],which_thread)  ) {
                (*all_gene[which_thread]) << "+ " << list_of_files[i] << std::endl;
                is_in = true;
                break;
            }
        }
        if ( !is_in ) {
            (*all_gene[modulo_on_thread++]) << "+ " << list_of_files[i] << std::endl;
            if ( modulo_on_thread == size_of_all_gene ) {
                modulo_on_thread = 0;
            }
        }
    }
}

//! function to split a list of ROOT file ()from the current directory into configuration files used by toROOTGPS
void DOtoROOTGPS_gene( const char *path, const char *rootpattern, const char *treename, const char *confbasename, size_t nb_thread )
{
    // vector of file extracted fron the path and following the rootpattern
    std::vector<std::string> list_of_files;

    TSystemDirectory dir_run; dir_run.SetDirectory(path); TList *list; TObject *entry; TString opt(rootpattern);
    //
    list = dir_run.GetListOfFiles();
    if ( list == 0x0 ) {
        std::cout << path <<  " Is Empty " << std::endl;
        return ; // no files
    }
    else {
        std::cout << list->GetEntries() << " files found in " << path << std::endl;
    }

    
    TIter iter_file(list); // start collecting .adf file
    //
    while ( (entry = iter_file()) )
    {
        TString fname = entry->GetName();
        
        Bool_t wildcard = false;
        if ( opt.Contains("*") )
            wildcard = true;
        TRegexp all("*",kTRUE), pattern(opt.Data(),wildcard);
        
        if ( pattern.Status() != TRegexp::kOK )
        {
            if ( gDebug > 0 )
                std::cerr << "pattern [*] intead of " << opt.Data()<< std::endl;
            pattern = all;
        }
        else
        {
            if ( gDebug > 0 ) std::cerr << "pattern " << opt.Data() << std::endl ;
        }
        if ( entry->InheritsFrom("TSystemFile") )
        {
            if ( fname.Contains( pattern ) ) {
                std::string tmp = path;
                if ( tmp[tmp.size()-1] != '/' )
                    tmp += "/";
                tmp += fname.Data();
                list_of_files.push_back( tmp );
                std::cout << " + ADD " << fname.Data() << " to the list of files to be treated " <<  std::endl;
            }
        }
    }
    if ( list_of_files.size() ) {
        DOtoROOTGPS_gene(list_of_files,treename,confbasename,nb_thread);
    }
}

void TEST()
{
    std::vector<std::string> list_of_files;
    /*
    list_of_files.push_back( "toROOTGPS_Thread00_Run000.root"  );
    list_of_files.push_back( "toROOTGPS_Thread01_Run000.root"  );
    list_of_files.push_back( "toROOTGPS_Thread02_Run000.root"  );
    list_of_files.push_back( "toROOTGPS_Thread03_Run000.root"  );
     */
    
    list_of_files.push_back( "../Run000.root"  );
    list_of_files.push_back( "../Run001.root"  );
    list_of_files.push_back( "../Run002.root"  );
    list_of_files.push_back( "../Run003.root"  );
    
    list_of_files.push_back( "../Run004.root"  );
    list_of_files.push_back( "../Run005.root"  );
    list_of_files.push_back( "../Run006.root"  );
    list_of_files.push_back( "../Run007.root"  );
    
    DOtoROOTGPS_gene(list_of_files,"TEST","SToGS",0);
}


