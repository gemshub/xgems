# Require a certain version of cmake
cmake_minimum_required(VERSION 3.0)

# Set the name of the project
project(GEMS3K)

# Recursively collect all header files from the current directory
file(GLOB HEADER_FILES *.h)

# Recursively collect all source files from the current directory
file(GLOB SOURCE_FILES *.cpp)

# Set some necessary definitions
add_definitions(-DIPMGEMPLUGIN)
add_definitions(-DNOPARTICLEARRAY)

# Create the static library
add_library(GEMS3K STATIC ${HEADER_FILES} ${SOURCE_FILES})