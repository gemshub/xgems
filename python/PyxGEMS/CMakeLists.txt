# Collect the source files of PyxGEMS
file(GLOB_RECURSE HEADER_FILES *.hpp)
file(GLOB_RECURSE SOURCE_FILES *.cpp)

# Include the ThermoFun/python directory
include_directories(${PROJECT_SOURCE_DIR}/python)
include_directories(${CONDA_AWARE_PREFIX}/include/eigen3)
include_directories(/usr/local/include/eigen3)

# Create the PyThermoFun library
pybind11_add_module(PyxGEMS ${SOURCE_FILES} NO_EXTRAS) # NO_EXTRAS used to decrease linking time, but causes larger module file size

# Link PyThermoFun against ThermoFun library
target_link_libraries(PyxGEMS LINK_PRIVATE XGEMS_STATIC GEMS3K-static ThermoFun ChemicalFun )
