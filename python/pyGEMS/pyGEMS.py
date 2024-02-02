"""
A simple object oriented interface to control gems in python with pythonic 
naming convention and dictionaries as a more convenient wrapper on pyxGEMS.
Created by R.A.Patel (c) 2021, extended by G.D.Miron and D.A.Kulik (c) November 2021.
Updated/extended by D.A.Kulik (c) December 2022
"""
#%%import libraries
import xgems           # import pyxgems ?
import numpy as np
#%%

#class GEMSParameter:
#    def __init__(self, symbol, type_, value, error, distribution='normal'):
#        self.symbol = symbol
#        self.value = value
#        self.type_ = type_
#        if distribution=='normal':
#            self.distribution = cp.Normal(value, error)

class GEMS(object):
    """
    Gems interface in calculator format for easy Python using dictionaries
    """
    def __init__(self,input_file, reset_calc=True, cold_start=True):
        """
        Initialization of the pyxGEMS calculator
        reset_calc: True will reset b vector to zero
        """
        self.input_file = input_file
        self.gem = xgems.ChemicalEngine(input_file)
        self.T = self.gem.temperature() # K
        self.P = self.gem.pressure() # Pa
        self.b = np.array(self.gem.elementAmounts()) # moles 
        if cold_start: self.cold_start()
        self.equilibrate()
        self.nelements = self.gem.numElements() # number of elements in the system
        self.nphases = self.gem.numPhases() # number of phases in the system
        self.nspecies = self.gem.numSpecies() # number of species in the system
        self.element_names = [] # 
        self.element_molar_masses = {} # dictionary containing elements and their molar masses
        elemolarmass = self.gem.elementMolarMasses()
        self.aq_phase_symbol = self.gem.phaseName(0)
        #self.aq_phase_symbol = 'aq_gen'               Potential problem if no aq phase in system
        #self.gas_phase_symbol = self.gem.phaseName(1) # use index 0 if no aq phase in the system, or
        self.gas_phase_symbol = 'gas_gen'             # !!!!!!!!!!! Potential problem if no gas phase in system
                                                        #             or it has a different name
        for i in range(self.nelements):
            self.element_names.append(self.gem.elementName(i))
            self.element_molar_masses[self.gem.elementName(i)] = elemolarmass[i]
        self.species_names =[]
        for i in range(self.nspecies):
            self.species_names.append(self.gem.speciesName(i))
        self.phase_names = []
        for i in range(self.nphases):
            self.phase_names.append(self.gem.phaseName(i))
        self.species_in_phase = {} # dictionary containing species in phase
        for i in range(self.nphases):
            if self.gem.numSpeciesInPhase(i) > 0: 
                  self.species_in_phase[self.phase_names[i]] = self.species_names[self.gem.indexFirstSpeciesInPhase(i):
                      self.gem.indexFirstSpeciesInPhase(i)+self.gem.numSpeciesInPhase(i)]
        self.species_charges={}
        for i in range(self.nspecies):
            self.species_charges[self.species_names[i]] = self.gem.speciesCharge(i)
        molar_mass = self.gem.speciesMolarMasses()
        self.species_molar_mass={} # dictionary containing species and their molar masses
        for i in range(self.nspecies):
            self.species_molar_mass[self.species_names[i]] = molar_mass[i]
        self.species_molar_volumes = {} # dictionary containing species and their molar volumes
        for i in range(self.nspecies):
            self.species_molar_volumes[self.species_names[i]] = self.gem.standardMolarVolume(i) 
        self.formulaMatrix=self.gem.formulaMatrix().T

        if reset_calc:self.clear()

    def clear(self,cvalue=1e-15):
        if cvalue > 0:
            self.b[:] = cvalue
        if 'Zz' in self.element_names:
            self.b[self.nelements-1] = 0.0

    def set_gems_parameters(self, parameters):
        for p in parameters:
            if p.type_ == 'G0':
                self.set_species_G0(p.symbol, p.value)

    def set_species_G0(self, symbol, value):
        self.gem.setStandardMolarGibbsEnergy(symbol, value)

    @property
    def bulk_composition(self):
        """
        return input bulk elemental composition (vector b) in moles
        """
        out = {}
        for i in range(self.nelements):
            out[self.element_names[i]] =self.b[i]
        return out
    vector_b = bulk_composition   # alias     

    @property
    def aq_elements_molarity(self):
        """
        aq solution composition in mol/L aq solution
        """
        out = {}
        aq_index = self.gem.indexPhase(self.aq_phase_symbol)
        moles_elements = self.gem.elementAmountsInPhase(aq_index)
        for i, v in enumerate(moles_elements):
            out[self.gem.elementName(i)] = v / (self.gem.phaseVolume(self.phase_names.index(self.aq_phase_symbol))*1000) 
            # volume from m3 to L
        return out
    aq_el_M = aq_elements_molarity        
     
    @property
    def aq_elements_molality(self):      # def aq_species_composition(self):
        """
        aq solution elemental composition in mol/kgH2O
        """
        aq_index = self.gem.indexPhase(self.aq_phase_symbol)
        H2Oindex = self.gem.numSpeciesInPhase(aq_index)-1
        H2Oamount = self.gem.speciesAmounts()[H2Oindex]
        H2Ommass = self.gem.speciesMolarMasses()[H2Oindex]
        H2Omass = H2Oamount*H2Ommass/1000 # in kg
        out = {}
        moles_elements = self.gem.elementAmountsInPhase(aq_index)
        for i, v in enumerate(moles_elements):
            out[self.gem.elementName(i)] = v / H2Omass # volume from m3 to L
        return out
    aq_el_my = aq_elements_molality     
    aq_species_composition = aq_elements_molality  # this alias is misleading!   
    
    @property
    def aq_species_molarity(self):
         """
         aq solution composition in mol/L of aqueous species
         """
         out = {}
         aq_index = self.gem.indexPhase(self.aq_phase_symbol)
         moles_species = self.gem.speciesAmounts()
         for i, v in enumerate(moles_species):
             out[self.gem.speciesName(i)] = v / (self.gem.phaseVolume(self.phase_names.index(self.aq_phase_symbol))*1000) 
             #volume from m3 to L
         return out 
    aq_sp_M = aq_species_molarity     
    
    @property
    def aq_species_molality(self):
        """
        aq solution composition in mol/kg H2O of aqueous species (speciation)
        """
        out = {}
        molalities =  self.gem.speciesMolalities()
        aq_species_names=self.species_in_phase[self.aq_phase_symbol][:-1]
        for name in aq_species_names:
            out[name] =molalities[self.species_names.index(name)]
        return out
    aq_sp_my = aq_species_molality
    aq_composition = aq_species_molality   # this alias is misleading!       
    
    @property
    def aq_elements_moles(self):
        """
        aq solution elements amount in moles
        """
        out = {}
        aq_index = self.gem.indexPhase(self.aq_phase_symbol)
        moles_elements = self.gem.elementAmountsInPhase(aq_index)
        for i, v in enumerate(moles_elements):
            out[self.gem.elementName(i)] = v
        return out
    aq_elements_amounts = aq_elements_moles    

    def set_bulk_composition(self, b_input): 
        """
        set input bulk elemental composition (vector b) in moles
        """
        for i in range(self.nelements):
            if self.element_names[i] in b_input.keys():
                self.b[i] = b_input[self.element_names[i]]
                if self.b[i] < 1e-15:
                    self.b[i] = 1e-15
                if self.element_names[i] == 'Zz':
                    self.b[i] = 0.0
    set_vector_b = set_bulk_composition                

    def reset_aq_composition(self):  # Not quite clear what this access method really does (DK)
        """
        Removes bulk elemental aqueous solution composition from vector b
        be careful as this will also remove water i.e H+ and OH-
        """
        peamt = self.gem.phaseAmounts()
        aqupx = self.gem.indexPhase(self.aq_phase_symbol)
        b_aqup = self.b
        # b_aqup[:]=0.0
        if aqupx < self.nphases:
            if peamt[aqupx] > 1e-12:
                b_aqup = self.gem.elementAmountsInPhase(aqupx)
                # b_aqup = np.array(self.gem.elementAmountsInPhase(aqupx))
                self.b -= b_aqup
        for i in range(self.nelements):
            if self.b[i] < 1e-12:
                self.b[i] = 1e-12
    clear_b_aq_part = reset_aq_composition            

    @property
    def solids_elements_moles(self):
        """
        return a dictionary containing mole amounts of elements in all solids together
        """
        b_solid = np.array(self.gem.elementAmounts())
        b_aqup = self.b
        b_gasp = self.b
        peamt = np.array(self.gem.phaseAmounts()) # self.gem.phaseAmounts()
        aqupx = self.gem.indexPhase(self.aq_phase_symbol)
        if aqupx < self.nphases:
            if peamt[aqupx] > 1e-12:
                b_aqup = self.gem.elementAmountsInPhase(aqupx)
                b_solid -= b_aqup
        gaspx = self.gem.indexPhase(self.gas_phase_symbol)
        if gaspx < self.nphases:
            if peamt[gaspx] > 1e-12:
                b_gasp = self.gem.elementAmountsInPhase(gaspx)
                b_solid -= b_gasp
        out = {}
        for i in range(self.nelements):
            if b_solid[i] < 1e-14:
                b_solid[i] = 0.0
            out[self.element_names[i]] = b_solid[i]
        return out
    solid_elements_amounts = solids_elements_moles

    @property
    def phases_elements_moles(self):
        """
        return a dictionary (table) containing amounts of elements in phases in moles
        """
        out = {}
        for k in range(self.nphases):
            peamt = self.gem.elementAmountsInPhase(k)
            dictelems = {}
            for i in range(self.nelements):
                dictelems[self.element_names[i]] = peamt[i]
            out[self.phase_names[k]] = dictelems
#     for pname in self.phase_names:
#         peamt = self.gem.elementAmountsInPhase(self.phase_names.index(pname))
#         dictelems = {}
#         for ename in self.element_names:
#             dictelems[ename] = peamt[self.element_names.index(ename)]
#         out[pname] = dictelems
        return out
    phase_elements_amounts = phases_elements_moles

    @property
    def phases_moles(self):
        """
        return phases amounts in moles
        """
        peamt = self.gem.phaseAmounts()
        out = {}
        for name in self.phase_names:
            out[name]=peamt[self.phase_names.index(name)]
        return out
    phase_amounts = phases_moles

    @property
    def species_moles(self):
        """
        returns all species amounts in moles
        """
        speamt = self.gem.speciesAmounts()
        out = {}
        for name in self.species_names:
            out[name]=speamt[self.species_names.index(name)]
        return out
    species_amounts = species_moles
     
    @property
    def species_ln_activities(self):
        """
        returns species ln(activities)
        """
        speamt = self.gem.lnActivities()
        out = {}
        for name in self.species_names:
            out[name]=speamt[self.species_names.index(name)]
        return out
    
    @property
    def species_ln_activity_coefficients(self):
        """
        returns species ln(activity_coefficient)
        """
        speamt = self.gem.lnActivityCoefficients()
        out = {}
        for name in self.species_names:
            out[name]=speamt[self.species_names.index(name)]
        return out
    
    def phase_species_moles(self, phase_symbol):
        """
        returns species in phase in moles
        """
        out = {}
        amounts =  self.gem.speciesAmounts()
        species_names=self.species_in_phase[phase_symbol]
        for name in species_names:
            out[name] =amounts[self.species_names.index(name)]
        return out
    phase_species_amounts = phase_species_moles

    @property
    def solids_mass_frac(self):
        """
        mass(phase)/mass(system) ratios for [solid] phases
        """
        mfrac = self.gem.phaseMasses()/np.sum(self.gem.phaseMasses())
        out = {}
        for name in self.phase_names:
            out[name]=mfrac[self.phase_names.index(name)]
        return out
    solid_mass_frac = solids_mass_frac

    @property
    def solids_volume_frac(self):
        """
        volume(phase)/volume(total) ratio for solid phases
        """
        out = self.phases_volume_frac
        del out[self.aq_phase_symbol]
        gaspx = self.gem.indexPhase(self.gas_phase_symbol)
        if gaspx < self.nphases:
            del out[self.gas_phase_symbol]
        return out
    solid_volume_frac = solids_volume_frac

    @property
    def aq_volume_frac(self):
        """
        Volume fraction of aqueous phase in the system
        """
        return self.phase_volume_frac[self.aq_phase_symbol]

    @property
    def phases_volume(self):
        """
        returns a dict. with phases and their absolute volume in m3
        """
        v = self.gem.phaseVolumes()
        out = {}
        for name in self.phase_names:
            out[name]=v[self.phase_names.index(name)]
        return out
    phase_volumes = phases_volume    

    @property
    def phases_mass(self):
        """
        returns a dict. with phases and their mass in kg
        """
        mass = self.gem.phaseMasses()
        out = {}
        for name in self.phase_names:
            out[name]=mass[self.phase_names.index(name)]
        return out
    phase_masses = phases_mass    

    @property
    def phases_volume_frac(self):
        """
        returns a dict. with phases and their volume fractions in the system
        """
        vfrac = self.gem.phaseVolumes()/self.system_volume
        out = {}
        for name in self.phase_names:
            out[name]=vfrac[self.phase_names.index(name)]
        return out
    phase_volume_frac = phases_volume_frac    

    def equilibrate(self):
        """
        runs GEM equilibriation of the current (internally set) system
        """
        outcode= self.gem.equilibrate(self.T,self.P,self.b)
        return self._status_encoder[outcode]

    def cold_start(self):
        self.gem.setColdStart()

    def warm_start(self):
        self.gem.setWarmStart()

    def add_multiple_species_amt(self,input_dict,units="moles"):
        """
        add species amount in the system useful for adding aqueous solution composition
        units= moles, kg, m3
        """
        for name, val in input_dict.items():
            self.add_species_amt(name,val,units)

    def add_species_amt(self,species,val,units="moles"):
        """
        add species amount in the system useful for adding aqueous solution composition
        units= moles, kg, m3
        """
        species_idx = self.species_names.index(species)
        if units == "kg":
            val/=self.species_molar_mass[species]
        if units == "m3":
            val/=self.species_molar_volumes[species]
        self.b += self.formulaMatrix[species_idx]*val

    def add_element_amt(self,element_name,val, units = "moles"):
        """
        add element amount in the system
        units = moles, kg
        """
        if units  == "kg":
            val /= self.element_molar_masses[element_name]
        self.b[self.element_names.index(element_name)]+=val

    def add_multiple_elements_amt(self,input_dict,units="moles"):
        """
        add elements amount in the system useful for adding aqueous solution composition
        units= moles,kg
        """
        for name, val in input_dict.items():
            self.add_element_amt(name,val,units)

    def add_amt_from_formula(self,formula, val, units="moles"):
        """
        add element amount using user defined formula
        units = moles,kg
        """
        if units == "kg":
            molarmass =0
            for element in formula.keys():
                molarmass += formula[element] * self.element_molar_masses[element]
            val/=molarmass
        for element in formula.keys():
            self.add_element_amt(element, val * formula[element])

    def get_b_from_formula(self,formula, val=1, units="moles"):
        """
        returns a bulk vector b from user-defined formula (as dict. {"H":2,"O":1} ) 
        and amount of the formula [object] in units of 'moles' or 'kg'
        """
        bx = [v for v in self.b]
        #bx[:] = 0.0
        if units == "kg":
            molarmass =0
            for element in formula.keys():
                molarmass += formula[element] * self.element_molar_masses[element]
            val/=molarmass
        for element in formula.keys():
            bx[self.element_names.index(element)]+=val * formula[element]
        return np.array(bx)
    vector_b_from_formula = get_b_from_formula    

    def set_multiple_species_lower_bound(self, input_dict,units="moles"):
        """
        constrain species amount to a specified lower bound
        units= moles,kg,m3
        """
        for species,val in input_dict.items():
            self.set_species_lower_bound(species,val,units)
    multiple_species_lower_bound = set_multiple_species_lower_bound        

    def set_multiple_species_upper_bound(self, input_dict,units="moles"):
        """
        constrain species amount to a specified lower bound
        units= moles,kg,m3
        """
        for species,val in input_dict.items():
            self.set_species_upper_bound(species,val,units)
    multiple_species_upper_bound = set_multiple_species_upper_bound        

    def set_species_lower_bound(self, species,val,units="moles"):
        """
        constrain species amount to a specified lower bound
        units= moles,kg,m3
        """
        if units == "kg":
            val/=self.species_molar_mass[species]
        if units == "m3":
            val/=self.species_molar_volumes[species]
        species_idx = self.species_names.index(species)
        self.gem.setSpeciesLowerLimit(species_idx,val)
    species_lower_bound = set_species_lower_bound    

    def set_species_upper_bound(self,species,val,units="moles"):
        """
        constrain species amount to a specified upper bound
        units= moles,kg,m3
        """
        if units == "kg":
            val/=self.species_molar_mass[species]
        if units == "m3":
            val/=self.species_molar_volumes[species]
        self.gem.setSpeciesUpperLimit(species,val)
    species_upper_bound = set_species_upper_bound    

    def supress_phase(self,phase_name):
        """
        supresses a phase in GEM calculation
        """
        for species in self.species_in_phase[phase_name]:
            self.supress_species(species)

    def supress_multiple_phases(self,phase_name_list):
        """
        supresses multiple phase in calculation as given in phase names list
        """
        for phase in phase_name_list:
            self.supress_phase(phase)

    def supress_species(self,species_name):
        """
        supresses species in calculation
        """
        self.set_species_lower_bound(species_name,0)
        self.set_species_upper_bound(species_name,1e-15)
        
    def supress_multiple_species(self,species_list):
        """
        supresses multiple species in GEM calculation as given in species name list
        """
        for species in species_list:
            self.supress_species(species)

    def activate_phase(self,phase_name):
        """
        activate supressed phase
        """
        for species in self.species_in_phase[phase_name]:
            self.activate_species(species)

    def activate_multiple_phases(self,phase_name_list):
        """
        activate multiple supressed phases given in list
        """
        for phase in phase_name_list:
            self.activate_phase(phase)

    def activate_multiple_species(self,species_list):
        """
        activate multiple supressed species given in the list
        """
        for species in species_list:
            self.activate_species(species)

    def activate_species(self,species_name):
        """
        activate a supressed species in phase
        """
        self.set_species_lower_bound(species_name,0)
        self.set_species_upper_bound(species_name,1e6)

    @property
    def pH(self):
        """
        returns pH of the solution
        """
        return self.gem.pH()

    @property
    def pE(self):
        """
        returns pE of the solution
        """
        return self.gem.pe()

    @property
    def ionic_strength(self):
        """
        returns ionic strength of the solution
        """
        return self.gem.ionicStrength()
    IS = ionic_strength    

    @property
    def system_volume(self):
        """
        returns volume of the system in m3
        """
        return self.gem.systemVolume()

    @property
    def system_mass(self):
        """
        returns mass of the system in kg
        """
        return self.gem.systemMass()

    @property
    def phases_molar_volume(self):
        """
        returns molar volume of phases in m3/mol
        """
        phase_mvol = {}
        for i in range(self.nphases):
            phase_mvol[self.phase_names[i]]=self.gem.phaseMolarVolume(i)
        return phase_mvol
    phase_molar_volume = phases_molar_volume   

    @property
    def phase_sat_indices(self):
        """
        returns saturation indices of phases
        """
        phaseSatIndices= self.gem.phaseSatIndices()
        satIndices = {}
        for i in range(self.nphases):
            satIndices[self.phase_names[i]]=phaseSatIndices[i]
        return satIndices
    phases_sat_index = phase_sat_indices
    phases_SI = phase_sat_indices 

    @property
    def _status_encoder(self):
        out = {}
        out[0]="No GEM re-calculation needed"
        out[1]= "Need GEM calculation with LPP (automatic) initial approximation (AIA)"
        out[2]="OK after GEM calculation with LPP AIA"
        out[3]="Bad (not fully trustful) result after GEM calculation with LPP AIA"
        out[4]="Failure (no result) in GEM calculation with LPP AIA"
        out[5]="Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation (full DATABR lists only)"
        out[6]="OK after GEM calculation with SIA"
        out[7]="Bad (not fully trustful) result after GEM calculation with SIA"
        out[8]="Failure (no result) in GEM calculation with SIA"
        out[9]="Terminal error has occurred in GEMS3K (e.g. memory corruption). Restart is required."
        return out        

# End of class GEMS

