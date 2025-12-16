message(STATUS "Building python library ")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

# Collect the source files of PyChemicalData
file(GLOB_RECURSE HEADER_FILES *.hpp)
file(GLOB_RECURSE SOURCE_FILES *.cpp)

# Include the xGEMS/python directory
include_directories(${PROJECT_SOURCE_DIR}/python)
include_directories(${CONDA_AWARE_PREFIX}/include/eigen3)
include_directories(/usr/local/include/eigen3)

# Create the PyxGEMS library
pybind11_add_module(PyxGEMS ${SOURCE_FILES} NO_EXTRAS) # NO_EXTRAS used to decrease linking time, but causes larger module file size

# Link PyxGEMS against xGEMS and GEMS3K library
target_link_libraries(PyxGEMS PUBLIC xGEMS-static)
