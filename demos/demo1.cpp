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
      // Looking for the first ':'
      size_t from = ss.str().find_first_of(':')+1; 
      // Looking for the last '}'
      size_t until = ss.str().find_last_of('}');
      // Getting json string from until-1
      str = ss.str().substr(from, until-from);
//      cout << str.substr(0,80) << endl;
      return 0;
    }
    cout << "Error reading file " << path << " into text string";
    return 1;
}

int main(int argc, char **argv)
{
    // ChemicalEngine chemicalengine("resources/CalciteIC-dat.lst");
    // Vector b = chemicalengine.elementAmounts();
    // std::cout << "Vector b: " << b.transpose() << std::endl;
    // std::cout << chemicalengine << std::endl;
    // chemicalengine.~ChemicalEngine();
    // ChemicalEngine chemicalengine2("resources/CemHyds-dat.lst");
    // Vector b2 = chemicalengine2.elementAmounts()
    // std::cout << "\nVector b2: " << b2.transpose() << std::endl;
    // std::cout << chemicalengine2 << std::endl;
    // Test input from JSON documents

    ChemicalEngine engine;
    
    int f1, f2, f3; 
    std::string dch_json = "";
    f1 = getJsonFile( "resources/dch_test.json", dch_json );

    std::string ipm_json = "";
    f2 = getJsonFile( "resources/ipm_test.json", ipm_json );

    std::string dbr_json = "";
    f3 = getJsonFile( "resources/dbr_test.json", dbr_json );
    
    if( f1 == 0 && f2 == 0 && f3 == 0 )
    {
        engine.initializeJson(dch_json, ipm_json, dbr_json);

        Vector bj = engine.elementAmounts();
    
        std::cout << "Vector bj: " << bj.transpose() << std::endl;

        std::cout << engine << std::endl;

        // engine.equilibrate(  );

        return 0;   
    }    
    return 1;
}
