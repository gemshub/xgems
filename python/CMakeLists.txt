# Find pybind11
#find_package(pybind11 REQUIRED)

# Build the Python wrappers of xGEMS library
add_subdirectory(PyxGEMS)

# Build the Python package for ThermoFun
add_subdirectory(xGEMS)
