// xGEMS is a C++ and Python library for thermodynamic modeling by Gibbs energy minimization
//
// Copyright (C) 2018 Allan Leal, Dmitrii Kulik
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License version 2.1 
// as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

// C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

// xGEMS includes
#include <xGEMS/ChemicalEngine.hpp>
#include "xGEMS/GEMSEngine.hpp"
using namespace xGEMS;

// Reading the whole file (given by path) into a string str
// and skipping the header e.g. {"dch": and the very last } at the end of file
int getJsonFile( const string path, string& str )
{
    ifstream f(path); //taking file as inputstream
    // string str;
    if(f) {
      ostringstream ss;
      ss << f.rdbuf(); // reading data
    //   // Looking for the first '{' to strip top-level [] 
    //   size_t from = ss.str().find_first_of('{'); 
    //   // Looking for the last '}'
    //   size_t until = ss.str().find_last_of('}');
    //   // Getting json string from until-1
    //   str = ss.str().substr(from, until-from+1);
      str = ss.str();
      // cout << str.substr(0,80) << endl;
      return 0;
    }
    cout << "Error reading file " << path << " into text string";
    return 1;
}

int main(int argc, char **argv)
{
    auto gemsd  = GEMSEngine("resources/Thermo-time-in/series1-dat.lst");
    Vector bj = gemsd.b_amounts;

    std::cout << "Vector bj: " << bj.transpose() << std::endl;
    return 0;
    // Testing with input of small keyvalue GEMS3K files
    // ChemicalEngine chemicalengine("resources/CalciteIC/CalciteIC-dat.lst");
    // Vector b = chemicalengine.elementAmounts();
    // std::cout << "Vector b: " << b.transpose() << std::endl;
    // std::cout << chemicalengine << std::endl;

//  Testing with input of large keyvalue GEMS3K files
    // ChemicalEngine chemicalengine1("resources/CemGEMS-keyvalue/CemHyds-dat.lst");
    // Vector b1 = chemicalengine1.elementAmounts();
    // std::cout << "\nVector b1: " << b1.transpose() << std::endl;
    // std::cout << chemicalengine1 << std::endl;
    
    // Test of input from JSON GEMS3K files
//    ChemicalEngine chemicalengine2("resources/CemGEMS-formatted/CemHyds-dat.lst");
    // ChemicalEngine chemicalengine2("resources/CemGEMS-condensed/CemHyds-dat.lst");
    // Vector b2 = chemicalengine2.elementAmounts();
    // std::cout << "\nVector b2: " << b2.transpose() << std::endl;
    // std::cout << chemicalengine2 << std::endl;

    // Testing= with input from JSON strings
    ChemicalEngine engine;
    int f1, f2, f3; 
    std::string dch_json = "";
    std::string ipm_json = "";
    std::string dbr_json = "";

//    f1 = getJsonFile( "resources/Test-json-strings/dch_test.json", dch_json );
//    f2 = getJsonFile( "resources/Test-json-strings/ipm_test.json", ipm_json );
//    f3 = getJsonFile( "resources/Test-json-strings/dbr_test.json", dbr_json );

//    f1 = getJsonFile( "resources/Kaolinite-formatted/pHtitrKaS-dch.json", dch_json );
//    f2 = getJsonFile( "resources/Kaolinite-formatted/pHtitrKaS-ipm.json", ipm_json );
//    f3 = getJsonFile( "resources/Kaolinite-formatted/pHtitrKaS-dbr-0-0000.json", dbr_json );

//    f1 = getJsonFile( "resources/Kaolinite-condensed/pHtitrKaS-dch.json", dch_json );
//    f2 = getJsonFile( "resources/Kaolinite-condensed/pHtitrKaS-ipm.json", ipm_json );
//    f3 = getJsonFile( "resources/Kaolinite-condensed/pHtitrKaS-dbr-0-0000.json", dbr_json );

//    f1 = getJsonFile( "resources/Export-formatted/Export-dch.json", dch_json );
//    f2 = getJsonFile( "resources/Export-formatted/Export-ipm.json", ipm_json );  
//    f3 = getJsonFile( "resources/Export-formatted/Export-dbr-0-0000.json", dbr_json ); 

    // f1 = getJsonFile( "resources/Export-condensed/Export-dch.json", dch_json );
    // f2 = getJsonFile( "resources/Export-condensed/Export-ipm.json", ipm_json );
    // f3 = getJsonFile( "resources/Export-condensed/Export-dbr-0-0000.json", dbr_json );

   f1 = getJsonFile( "resources/CemGEMS-formatted/CemHyds-dch.json", dch_json );
   f2 = getJsonFile( "resources/CemGEMS-formatted/CemHyds-ipm.json", ipm_json );
   f3 = getJsonFile( "resources/CemGEMS-formatted/CemHyds-dbr-0-0000.json", dbr_json );
    
    // f1 = getJsonFile( "resources/CemGEMS-condensed/CemHyds-dch.json", dch_json );
    // f2 = getJsonFile( "resources/CemGEMS-condensed/CemHyds-ipm.json", ipm_json );
    // f3 = getJsonFile( "resources/CemGEMS-condensed/CemHyds-dbr-0-0000.json", dbr_json );

//    f1 = getJsonFile( "resources/GEMSW-from-DB/GEMSW-dch.json", dch_json );
//    f2 = getJsonFile( "resources/GEMSW-from-DB/GEMSW-ipm.json", ipm_json );
//    f3 = getJsonFile( "resources/GEMSW-from-DB/GEMSW-dbr-0-0000.json", dbr_json );
    
    if( f1 == 0 && f2 == 0 && f3 == 0 )
    {
        engine.initializeFromJsonStrings(dch_json, ipm_json, dbr_json);

        Vector bj = engine.elementAmounts();
    
        std::cout << "Vector bj: " << bj.transpose() << std::endl;

        //std::cout << engine << std::endl;

        auto T = engine.temperature();
        auto P = engine.pressure();
        auto b = engine.elementAmounts();


        engine.equilibrate(T, P, b);
        return 0;   
    }    
    return 1;
}
