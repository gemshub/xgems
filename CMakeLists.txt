# Set cmake version requirement
cmake_minimum_required(VERSION 3.0)

# Set the name of the project
project(xGEMS)

# Set the cmake module path of the project
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Include this module to enable testing with CTest and dashboard submissions to CDash
include(CTest)

# Include the cmake variables with values for installation directories
include(GNUInstallDirs)

# Ensure proper configuration if in a conda environment
include(CondaAware)

# Set the version of the project
set(XGEMS_VERSION_MAJOR "1")
set(XGEMS_VERSION_MINOR "0")
set(XGEMS_VERSION_MICRO "0")
set(XGEMS_VERSION "${XGEMS_VERSION_MAJOR}.${XGEMS_VERSION_MINOR}.${XGEMS_VERSION_MICRO}")

# Set the output directories of the built libraries and binaries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Define which components of Reaktoro to build
option(XGEMS_BUILD_ALL    "Build everything." OFF)
option(XGEMS_BUILD_PYTHON "Build the python bindings." ON)
option(XGEMS_BUILD_DEMOS  "Build demos." OFF)
option(XGEMS_BUILD_DOCS   "Build documentation." OFF)
option(XGEMS_BUILD_TESTS  "Build tests." OFF)

# Modify the XGEMS_BUILD_X variables according to XGEMS_BUILD_ALL
if(XGEMS_BUILD_ALL MATCHES ON)
    set(XGEMS_BUILD_DEMOS  ON)
    set(XGEMS_BUILD_DOCS   ON)
    set(XGEMS_BUILD_PYTHON ON)
    set(XGEMS_BUILD_TESTS  ON)
endif()

# Set the CTest option BUILD_TESTING to the value of XGEMS_BUILD_TESTS
set(BUILD_TESTING XGEMS_BUILD_TESTS)

# Define which types of libraries to build
option(XGEMS_BUILD_SHARED_LIBS "Build shared libraries." ON)
option(XGEMS_BUILD_STATIC_LIBS "Build static libraries." ON)

# Set libraries to be compiled with position independent code mode (i.e., fPIC option in GNU compilers)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Cannot currently produce a dynamic library using MSVC, only a static library
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    set(XGEMS_BUILD_SHARED_LIBS OFF)
endif()

# Set the default build type to Release
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# add_compile_options(-std=c++11)
add_compile_options(-std=c++17)

# Set the list of compiler flags for GNU compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
#   add_compile_options(-std=c++11 -Wall -Wno-misleading-indentation -Wno-ignored-attributes -Wno-pedantic -Wno-variadic-macros -Wno-deprecated)
    add_compile_options(-Wall -Wno-misleading-indentation -Wno-ignored-attributes -Wno-pedantic -Wno-variadic-macros -Wno-deprecated)
endif()

# Set the list of compiler flags for Clang compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
   add_compile_options(-Wall -Wno-ignored-attributes -Wno-pedantic -Wno-variadic-macros -Wno-deprecated)
#  add_compile_options(-std=c++11 -Wall -Wno-ignored-attributes -Wno-pedantic -Wno-variadic-macros -Wno-deprecated)
endif()

# Set the list of compiler flags for Intel compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
    add_compile_options(-Wall -Wno-variadic-macros -Wno-deprecated)
endif()

# Set the list of compiler flags for MSVC compiler
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    add_compile_options("/W0 -D_SCL_SECURE_NO_WARNINGS /MP4")
endif()

# Find Eigen library and include its directory
#find_package(Eigen3 3.3.7 REQUIRED)
#message(STATUS "Found Eigen3 v${Eigen3_VERSION}: ${Eigen3_INCLUDE_DIRS}")

# Find the Python interpreter
if(XGEMS_BUILD_PYTHON)
    if(NOT PYBIND11_PYTHON_VERSION)
        set(PYBIND11_PYTHON_VERSION 2.7)
    endif()
    
    #find_package(pybind11 2.0.0 REQUIRED)
    find_package(pybind11 2.2.0 REQUIRED)
    message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
    
    if(NOT pybind11_FOUND)
        set(XGEMS_BUILD_PYTHON OFF)
        message(WARNING "pybind11 has not been found. "
            "Cannot build python bindings for XGEMS.")
    endif()
endif()

# Set the project directory names
set(XGEMS_SOURCE_DIR ${CMAKE_SOURCE_DIR}/xGEMS)
set(XGEMS_DEMOS_DIR  ${CMAKE_SOURCE_DIR}/demos)
set(XGEMS_DOCS_DIR   ${CMAKE_SOURCE_DIR}/docs)
set(XGEMS_TESTS_DIR  ${CMAKE_SOURCE_DIR}/tests)

# Set the include directories
include_directories(${CMAKE_SOURCE_DIR})
#include_directories(${EIGEN3_INCLUDE_DIR})


# Set the local directory where the third-party libraries are installed
set(THIRDPARTY_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/thirdparty/install)

# Set the include directories for third-party libraries
include_directories(${THIRDPARTY_INSTALL_PREFIX}/include)

# Set the link directories for third-party libraries
link_directories(${THIRDPARTY_INSTALL_PREFIX}/lib)

# Include the ExternalProject module
include(ExternalProject)

# Set the common cmake arguments to all external projects
set(GEMS3K_CMAKE_ARGS
    -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
    -DCMAKE_INSTALL_PREFIX=${THIRDPARTY_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
#    -DCMAKE_C_FLAGS=${CFLAGS}
#    -DCMAKE_CXX_FLAGS=${CXXFLAGS}
#    -DCMAKE_INSTALL_INCLUDEDIR=include
#    -DCMAKE_INSTALL_LIBDIR=lib
#    -DCMAKE_INSTALL_BINDIR=bin
#    -DBUILD_SHARED_LIBS=OFF
     # Disables the cmake warning message `Manually-specified
     # variables were not used by the project` because of
     # unused arguments above
    --no-warn-unused-cli
    )

# Clone GEMS3K from git repository and build it 
ExternalProject_Add(EXTERNAL_PROJECT_GEMS3K
    PREFIX build
#    GIT_REPOSITORY git@bitbucket.org:gems4/gems3k.git
    GIT_REPOSITORY https://bitbucket.org/gems4/gems3k.git
    GIT_TAG trunk
#    GIT_TAG arango    
    GIT_SHALLOW TRUE
    SOURCE_DIR build/src/gems3k
    BINARY_DIR build/src/gems3k-build
    STAMP_DIR build/src/gems3k-stamp
    CMAKE_ARGS ${GEMS3K_CMAKE_ARGS}
    )

 set(EIGEN_CMAKE_ARGS
    -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
    -DCMAKE_INSTALL_PREFIX=${THIRDPARTY_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DPYBIND11_TEST=OFF
    )
    # Clone pybind11 from git repository and build it 
 ExternalProject_Add(EXTERNAL_PROJECT_EIGEN
    PREFIX build
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG before-git-migration
    GIT_SHALLOW TRUE
    SOURCE_DIR build/src/eigen
   BINARY_DIR build/src/eigen-build
    STAMP_DIR build/src/eigen-stamp
    CMAKE_ARGS ${EIGEN_CMAKE_ARGS})

# Build the C++ library and the Python bindings
add_subdirectory(xGEMS)

# Build the demos
if(XGEMS_BUILD_DEMOS)
    add_subdirectory(demos)
else()
    add_subdirectory(demos EXCLUDE_FROM_ALL)
endif()

# Build the docs
if(XGEMS_BUILD_DOCS)
    add_subdirectory(docs)
else()
    add_subdirectory(docs EXCLUDE_FROM_ALL)
endif()

# Build the tests
if(XGEMS_BUILD_TESTS)
    add_subdirectory(tests)
else()
    add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()

# Add target "demos" for manual building of demos, as `make demos`, if XGEMS_BUILD_DEMOS is OFF
add_custom_target(demos
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/demos")

# Add target "tests" for manual building of tests, as `make tests`, if XGEMS_BUILD_TESTS is OFF
add_custom_target(tests
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests")
