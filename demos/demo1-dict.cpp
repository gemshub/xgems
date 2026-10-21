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

// xGEMS includes
#include "xGEMS/ChemicalEngineMaps.hpp"
#include <xGEMS/ChemicalEngine.hpp>

template <class T>
void print_map(const std::string title, const std::map<std::string, T>& data)
{
    std::cout << "\n" << title << std::endl;
    for(const auto& sp: data) {
        std::cout<< sp.first << ":" << sp.second << std::endl;
    }
}
static void demo_dict();
static void demo_warnings_map();
static void demo_warnings();
static void check_err();


int main(int argc, char **argv)
{
    //demo_dict();
    //demo_warnings_map();
    //demo_warnings();
    check_err();
    return 0;
}

void demo_warnings_map()
{

    auto engine = xGEMS::ChemicalEngineMaps("resources/solvus2/series1-dat.lst");
    std::cout << engine.equilibrate() << std::endl;

    print_map("engine.species_molar_mass()", engine.species_molar_mass());

    print_map("engine.element_molar_masses()", engine.element_molar_masses());

    print_map("engine.phases_moles()", engine.phases_moles());

}

void demo_warnings()
{

    xGEMS::ChemicalEngine engine("resources/solvus2/series1-dat.lst");

    auto T = engine.temperature();
    auto P = engine.pressure();
    auto b = engine.elementAmounts();
    engine.equilibrate(T, P, b);

    // Test setPT
    bool no_error = engine.setPT(101325, 298.15);
    std::cout << (no_error ? "PT set correctly" : "PT error") << std::endl;
    no_error = engine.setPT(50000000, 673.15);
    std::cout << (no_error ? "PT set correctly" : "PT error") << std::endl;

    //Test long names
    auto ispecies1 = engine.indexSpecies("Anorthite_long_substance_name");
    auto ispecies2 = engine.indexSpecies("Anorthite_long_substance_name", "Plagioclase");
    auto ispecies3 = engine.indexSpecies("Anorthite");
    std::cout << ispecies1 << " " << ispecies2 << " " <<ispecies3 << std::endl;
    std::cout << engine.speciesName(ispecies1) << " " << engine.speciesName(ispecies2) << " " <<engine.speciesName(ispecies3) << std::endl;

    auto iphase1= engine.indexPhase("Kaolinite_long_phase_name");
    std::cout << iphase1 << " " << engine.phaseName(iphase1) << std::endl;
    auto iphase2= engine.indexPhase("Kaolinite");
    std::cout << iphase2 << " " << engine.phaseName(iphase2) << std::endl;

    //indexSpecies dependen on phase optional
    std::cout << engine.speciesAmount(ispecies1) << " " << engine.speciesAmount(ispecies2) << std::endl;
    engine.setSpeciesAmount("Anorthite_long_substance_name", 10);
    engine.setSpeciesAmount("Anorthite_long_substance_name", 20, "Plagioclase");
    std::cout << engine.speciesAmount(ispecies1) << " " << engine.speciesAmount(ispecies2) << std::endl;

    print_map("Albite indexes", engine.indexSpeciesMap("Albite"));
    print_map("Quartz indexes", engine.indexSpeciesMap("Quartz"));

    // warning messages required
    std::cout << engine.phaseSpecificVolume(2) << std::endl;
    std::cout << engine.phaseSpecificVolume(30) << std::endl;
    std::cout << engine.phaseSpecificVolume(-1) << std::endl;

}


void demo_dict()
{
    auto engine = xGEMS::ChemicalEngineMaps("resources/CalciteBC/CalciteBC-dat.lst");

    //------------ calculate

    engine.T = 298.15;
    engine.P = 100000.0;
    xGEMS::ValuesMap  bulk_composition = { {"C", 1e-08}, {"Ca", 1e-08}, {"Cl", 0.002},
                                         {"H", 111.016746657646}, {"Mg", 0.001}, {"O", 55.5083933588231},
                                         {"Sn", 130.841288437146}, {"Zz", 0.0} };
    engine.set_bulk_composition(bulk_composition);
    std::cout << engine.equilibrate() << std::endl;

    //------------ results

    print_map("bulk_composition", engine.bulk_composition());

    std::cout << "pH: " << engine.pH()
              << ", pe: " << engine.pE()
              << ", Ionic strength: " << engine.ionic_strength()
              << ", System volume: " << engine.system_volume()
              << ", System mass: " << engine.system_mass()
              << std::endl;

    auto values = engine.element_molar_masses();
    std::cout << "Molar mass of 'O': " << values["O"] << std::endl;

    auto charge = engine.species_charges();
    std::cout << "Charge of 'OH-': " << charge["OH-"] << std::endl;

    auto values1 = engine.species_molar_mass();
    std::cout << "Molar mass of 'H2O@': " << values1["H2O@"] << std::endl;

    auto values2 = engine.species_molar_volumes();
    std::cout << "Molar volume of 'H2O@': " << values2["H2O@"] << std::endl;

    auto volumes = engine.phases_molar_volume();
    std::cout << "Molar volume of 'aq_gen': " << volumes["aq_gen"] << std::endl;

    auto sat_indices = engine.phase_sat_indices();
    std::cout << "Saturation indice of 'aq_gen': " << sat_indices["aq_gen"] << std::endl;

    auto molarity = engine.aq_elements_molarity();
    std::cout << "Molarity of 'O': " << molarity["O"] << std::endl;

    auto molality = engine.aq_elements_molality();
    std::cout << "Molality of 'O': " << molality["O"] << std::endl;

    auto molarity2 = engine.aq_species_molarity();
    std::cout << "Molarity of 'H2O@': " << molarity2["H2O@"] << std::endl;

    auto molality2 = engine.aq_species_molality();
    std::cout << "Molality of 'Mg(CO3)@': " << molality2["Mg(CO3)@"] << std::endl;

    auto amount = engine.aq_elements_moles();
    std::cout << "Amount of 'O' in the aqueous phase: " << amount["O"] << std::endl;

    auto amount2 = engine.solids_elements_moles();
    std::cout << "Amount of 'Sn' in the solids phases: " << amount2["Sn"] << std::endl;

    auto phase_el_moles = engine.phases_elements_moles();
    std::cout << "Amount of 'Ca' in 'aq_gen' phase: " << phase_el_moles["aq_gen"]["Ca"] << std::endl;

    auto amounts = engine.phases_moles();
    std::cout << "Amount of 'aq_gen' phase: " << amounts["aq_gen"] << std::endl;

    auto amounts2 = engine.species_moles();
    std::cout << "Amount of 'H2O@': " << amounts2["H2O@"] << std::endl;

    auto lnActivities = engine.species_ln_activities();
    std::cout << "ln Activity of 'H2O@': " << lnActivities["H2O@"] << std::endl;

    auto lnActCoeff = engine.species_ln_activity_coefficients();
    std::cout << "ln Activity coeff of 'H2O@': " << lnActCoeff["H2O@"] << std::endl;

    auto amounts3 = engine.phase_species_moles("aq_gen");
    std::cout << "Amount of 'H2O@' in 'aq_gen' phase: " << amounts3["H2O@"] << std::endl;

    auto mas_frac = engine.solids_mass_frac();
    std::cout << "Mass fraction of 'Tin' phase: " << mas_frac["Tin"] << std::endl;

    auto vol_frac = engine.solids_volume_frac();
    std::cout << "Volume fraction of 'Tin' phase from total system volume: " << vol_frac["Tin"] << std::endl;

    std::cout << "Aq volume fraction from total system volume: " << engine.aq_volume_frac() << std::endl;

    auto volumes4 = engine.phases_volume();
    std::cout << "Volume of 'Tin' phase: " << volumes4["Tin"] << std::endl;

    auto masses = engine.phases_mass();
    std::cout << "Mass of 'Tin' phase: " << masses["Tin"] << std::endl;

    auto volumes5 = engine.phases_volume_frac();
    std::cout << "Volume fraction of 'aq_gen' phase from total system volume: " << volumes5["aq_gen"] << std::endl;

    // ------------------- changes

    engine.add_species_amt("H2O@", 0.01, "kg");
    engine.add_multiple_species_amt({ {"Cl-",0.01}, {"H2@",2} }, "moles");
    engine.add_element_amt("Ca", 0.3, "moles");
    engine.add_multiple_elements_amt({ {"Cl",1.013077}, {"Sn",1.013077} }, "moles");
    engine.add_amt_from_formula( { {"Cl",2}, {"O",1} }, 4.108*1e-3, "kg");
    auto vect = engine.get_b_from_formula( { {"H",2},{"O",1}}, 0.1, "kg");

    engine.set_species_lower_bound( 8, 400, "moles");
    engine.set_species_upper_bound( 8, 900, "kg");
    engine.set_species_lower_bound( "Ca(HCO3)+", 200, "moles");
    engine.set_species_upper_bound( "CaOH+", 500, "kg");
    engine.set_multiple_species_lower_bound( {{"Mg(CO3)@",30}, {"Mg(HCO3)+",40}, {"Mg+2",50}});
    engine.set_multiple_species_upper_bound( {{"Mg(CO3)@",300}, {"Mg(HCO3)+",400}, {"Mg+2",500}});

    auto upperLimits = engine.species_upper_bounds();
    std::cout << "Upper bounds of 'Mg(HCO3)+': " << upperLimits["Mg(HCO3)+"] << std::endl;

    auto lowerLimits = engine.species_lower_bounds();
    std::cout << "Lower bounds of 'Mg(HCO3)+': " << lowerLimits["Mg(HCO3)+"] << std::endl;

    engine.suppress_phase("gas_gen");
    engine.suppress_multiple_phases({"Dolomite-dis", "Tin"});
    engine.suppress_species("Ca(CO3)@");
    engine.suppress_multiple_species({"ClO4-", "Cl-"});

    engine.activate_phase("gas_gen");
    engine.activate_multiple_phases({"Dolomite-dis", "Tin"});
    engine.activate_species("Ca(CO3)@");
    engine.activate_multiple_species({"Ca(HCO3)+", "CaOH+", "Mg(CO3)@", "Mg(HCO3)+", "Mg+2", "ClO4-", "Cl-"});


    //--------------------- phase data

    auto amounts11 = engine.phase_species_moles();
    std::cout << "Amount of 'CaOH+' in 'aq_gen' phase: " << amounts11["aq_gen"]["CaOH+"] << std::endl;

    auto lnActivities11 = engine.phase_species_ln_activities();
    std::cout << "ln activities of 'CaOH+' in 'aq_gen' phase: " << lnActivities11["aq_gen"]["CaOH+"] << std::endl;

    auto lnActCoeff11 = engine.phase_species_ln_activity_coefficients();
    std::cout << "ln activities coeff of 'CaOH+' in 'aq_gen' phase: " << lnActCoeff11["aq_gen"]["CaOH+"] << std::endl;

    auto upperLimits11 = engine.phase_species_upper_bounds();
    std::cout << "Upper limits of 'CaOH+' in 'aq_gen' phase: " << upperLimits11["aq_gen"]["CaOH+"] << std::endl;

    auto lowerLimits11 = engine.phase_species_lower_bounds();
    std::cout << "Lower limits of 'CaOH+' in 'aq_gen' phase: " << lowerLimits11["aq_gen"]["CaOH+"] << std::endl;

    print_map("engine.species_molar_mass()", engine.species_molar_mass());

    print_map("engine.element_molar_masses()", engine.element_molar_masses());

    print_map("engine.phases_moles()", engine.phases_moles());

}


void check_err()
{
    std::string path{"resources/gems_files/M1-6-dat.lst"};

    xGEMS::ChemicalEngine xgEngine(path);

    auto T = xgEngine.temperature();
    auto P = xgEngine.pressure();
    xGEMS::Vector b0 = xgEngine.elementAmounts();

    std::cout << "Task " << path << " T " << T<< " P " << P << std::endl;
    std::cout << b0 << std::endl;
    std::cout<< "reequilibrate " << xgEngine.reequilibrate() << std::endl;


    std::map<std::string, int>   element_indexes = {
        {"Pb", xgEngine.indexElement("Pb")},
        {"I", xgEngine.indexElement("I")},
        {"O", xgEngine.indexElement("O")},
        {"He", xgEngine.indexElement("He")}
    };

    std::map<std::string, int>   species_indexes = {
        {"Pb", xgEngine.indexSpecies("Pb")},
        {"I2", xgEngine.indexSpecies("I2")},
        {"I(g)", xgEngine.indexSpecies("I(g)")},
        {"I2(g)", xgEngine.indexSpecies("I2(g)")},
        {"Pb(g)", xgEngine.indexSpecies("Pb(g)")},
        {"Pb2(g)", xgEngine.indexSpecies("Pb2(g)")},
        {"PbI(g)", xgEngine.indexSpecies("PbI(g)")},
        {"PbI2(g)", xgEngine.indexSpecies("PbI2(g)")},
        {"PbI2(s)", xgEngine.indexSpecies("PbI2(s)")}
    };


    b0[element_indexes["Pb"]] = 1.61e-5;
    b0[element_indexes["I"]]  = 2 * b0[element_indexes["Pb"]];
    b0[element_indexes["He"]] = 0.75;
    b0[element_indexes["O"]]  = 2e-10;

    std::cout << "b0\n" << b0 << std::endl;
    P = 1e5;

    for(double t=25.; t<=1000.; t=t+1.) {
     auto code = xgEngine.equilibrate(t + 273.15, P, b0);
     std::cout << " T " << t << " code " << code << std::endl;


     if(code!=2 && code !=6) {
         std::cout <<"Warning: equilibrium not fully converged at " << t << " °C (code " << code << ")\n";
     }
    }
    std::cout <<"Finish" << std::endl;

}
