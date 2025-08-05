TEMPLATE = app
LANGUAGE = C++
TARGET = demo

CONFIG -= qt
CONFIG -= warn_on
CONFIG += debug
CONFIG += console
CONFIG += c++17

XGEMS_CPP = ..
XGEMS_H   = $$XGEMS_CPP
DEPENDPATH += $$XGEMS_H
INCLUDEPATH += $$XGEMS_H

DEFINES += IPMGEMPLUGIN
DEFINES += NODEARRAYLEVEL
DEFINES += USE_THERMOFUN
DEFINES += USE_THERMO_LOG
#DEFINES += OVERFLOW_EXCEPT  #compile with nan inf exceptions

#GEMS3K_CPP = ../../standalone/GEMS3K
#GEMS3K_H   = $$GEMS3K_CPP
#DEPENDPATH += $$GEMS3K_H
#INCLUDEPATH += $$GEMS3K_H
#GEM2MT_CPP = ../../standalone/nodearray-gem
#GEM2MT_H   = $$GEM2MT_CPP
#DEPENDPATH += $$GEM2MT_H
#INCLUDEPATH += $$GEM2MT_H

#include($$GEMS3K_CPP/gems3k.pri)
#include($$GEM2MT_CPP/gem2mt.pri)

OBJECTS_DIR = obj

LIBS += -lGEMS3K
contains(DEFINES, USE_THERMOFUN) {
    LIBS += -lThermoFun -lChemicalFun
} ## end USE_THERMOFUN

HEADERS	 += \
    $$XGEMS_H/xGEMS/ChemicalEngine.hpp \
    $$XGEMS_H/xGEMS/ChemicalEngineMaps.hpp \
    $$XGEMS_H/xGEMS/Eigen.hpp \
    $$XGEMS_H/xGEMS/Index.hpp \
    $$XGEMS_H/xGEMS/Interface.hpp
SOURCES	+=  demo1-dict.cpp \
    $$XGEMS_CPP/xGEMS/ChemicalEngine.cpp \
    $$XGEMS_CPP/xGEMS/ChemicalEngineMaps.cpp



