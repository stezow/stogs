

//! read aeuler file as given in the AGATA G4 package and creates a .dfb that can be used in SToGS
void ConvertAeulerTodfb(const Char_t *aeuler, const Char_t *dfb)
{
    ifstream infil; infil.open(aeuler);
    if ( !infil.is_open() ) {
        G4cout << "[ConvertAeulerTodfb] *** Cannot open input file " << aeuler << endl; // end read the file.
        return 0x0 ;
    }
    ofstream oufil; oufil.open(dfb);
    if ( !oufil.is_open() ) {
        G4cout << "[ConvertAeulerTodfb] *** Cannot open input file " << dfb << endl; // end read the file.
        return 0x0 ;
    }
    std::string line;  Int_t lline, i1, i2, nclust = 0; Double_t ps, th, ph, x, y, z; std::string sign_p[3], sign_r[3];
    
    //
    oufil << "# " << endl;
    oufil << "# First the world / detector envelop. if detector does not exist, box with dim given otherwise dim ignored " << endl;
    oufil << "# " << endl;
    oufil << "w AGATA +2. +2. +2. m " << endl;
    oufil << "# " << endl;
    oufil << "@ 0 -> " << endl;
    oufil << "# " << endl;
    oufil << "# " << endl;

    // read the input file
    while( infil.good() ) {
        
        getline(infil,line);
        //
        if ( line.size() < 2u )
            continue;
        if ( line[0] == '#' )
            continue;
        
        // decode the line
        if(sscanf(line.data(),"%d %d %lf %lf %lf %lf %lf %lf", &i1, &i2, &ps, &th, &ph, &x, &y, &z) != 8) {
            break;
        }
        if ( ps < 0.0 ) {
            sign_r[0] = "";
        }
        else sign_r[0] = "+";
        if ( th < 0.0 ) {
            sign_r[1] = "";
        }
        else sign_r[1] = "+";
        if ( ph < 0.0 ) {
            sign_r[2] = "";
        }
        else sign_r[2] = "+";
        //
        if ( x < 0.0 ) {
            sign_p[0] = "";
        }
        else sign_p[0] = "+";
        if ( y < 0.0 ) {
            sign_p[1] = "";
        }
        else sign_p[1] = "+";
        if ( z < 0.0 ) {
            sign_p[2] = "";
        }
        else sign_p[2] = "+";
        
        oufil << "* DetectorFactory/SemiConductors/Ge/AGATA-TC\t"
                << sign_p[0] << x << "\t" << sign_p[1] << y << "\t" << sign_p[2] << z << "\tmm"
                << "\tRz " << sign_r[0] << ps << "\tRy " << sign_r[1] << th << "\tRz " << sign_r[2] << ph << "" << endl;
    }
    oufil << "# " << endl;
    oufil << "# end " << endl;
    
}