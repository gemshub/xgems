# -*- coding: utf-8 -*-
"""
Simple demo to illustrate use of xGems api
"""
from __future__ import division,print_function
import os,sys
mygemspath = '../build/lib'#change this to path where xgems is compiled. If directory structure is same as in repo dont change this
sys.path.append(mygemspath)
import xgems
input_file =os.path.join('resources','CalciteIC/CalciteIC-dat.lst')
gem = xgems.ChemicalEngine(input_file)
#%% equilibriate system
T = gem.temperature() #get temprature
P = gem.pressure() #get pressure
b = gem.elementAmounts() #get element amount
gem.equilibrate(T, P, b) #equilibriate
print ("number of iterations to reach equilibrium: %s"%gem.numIterations())
#%% print some system information
nelements = gem.numElements()
nphases = gem.numPhases()
nspecies = gem.numSpecies()
print ("number of elements in the defined system: %s" %nelements)
print ("number of species in the defined system: %s"%nspecies)
print("number of phases in the defined system: %s"%nphases)
print ("temprature: %s pressure: %s volume: %s pH: %s pE: %s" 
%(gem.temperature(),gem.pressure(),gem.systemVolume(),gem.pH(),gem.pe()))

#%%
print ("lets get name of the elements and corresponding indexes...")
element_names = []
for i in range(nelements):
    print("index:%s,name:%s"%(i,gem.elementName(i)))
    element_names.append(gem.elementName(i))
#%%
print ("lets get name of the phases and corresponding indexes...")
phase_names =[]
for i in range(nphases):
    print("index:%s,name:%s"%(i,gem.phaseName(i)))
    phase_names.append(gem.phaseName(i))
#%%
print ("lets get name of the species and corresponding indexes...")
species_names =[]
for i in range(nspecies):
    print("index:%s,name:%s"%(i,gem.speciesName(i)))
    species_names.append(gem.speciesName(i))
#%%
print ("printing element amounts...")
element_amounts = gem.elementAmounts()
print("Name \t Amounts (in moles)")
for i in range(nelements):    
    print("%s \t %s"%(element_names[i],element_amounts[i]))
#%%
print ("printing phases amounts...")
phase_amounts = gem.phaseAmounts()
print("Name \t Amounts (in moles)")
for i in range(nphases):
    print("%s \t %s"%(phase_names[i],phase_amounts[i]))
#%%
print ("printing amount of elements in each phases")
for i in range(nphases):
    temp = gem.elementAmountsInPhase(i)
    print("------------------")
    print("%s"%phase_names[i])
    print("------------------")
    for i in range(nelements):    
        print("%s \t %s"%(element_names[i],temp[i]))
#%%
print ("printing species amounts...")
species_amounts = gem.speciesAmounts()
print("Name \t \t Amounts (in moles)")
for i in range(nspecies):    
    print("%s \t \t %s"%(species_names[i],species_amounts[i]))
#%%
print("obtaining index of element from element name...")
print("index of %s is %s"%("Ca",gem.indexElement("Ca")))

#%%
print ("printing phase amounts...")
temp = gem.phaseAmounts()
for i,name in enumerate(phase_names):
    print("%s: %s"%(name,temp[i]))
    
print ("printing phase volumes...")
temp = gem.phaseVolumes()
for i,name in enumerate(phase_names):
    print("%s: %s"%(name,temp[i]))

print ("printing phase volumes fraction...")
temp = gem.phaseVolumes()
for i,name in enumerate(phase_names):
    print("%s: %s"%(name,temp[i]/gem.systemVolume()))
