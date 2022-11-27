# Collect all hpp files from the current directory
file(GLOB HPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)

# Collect all cpp files from the current directory
file(GLOB CPP_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# Compile the python binding files into a Python module called xgems
if(XGEMS_BUILD_PYTHON)
# Collect all cpp files from the current directory
    file(GLOB CXX_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cxx)
endif()

# Compile the C++ source code into object files
add_library(XGEMS_OBJECT OBJECT ${HPP_FILES} ${CPP_FILES})

target_compile_definitions(XGEMS_OBJECT PUBLIC -DSPDLOG_FMT_EXTERNAL=yes) 

# Specify the GEMS3K dependency to ensure proper build sequence
#add_dependencies(XGEMS_OBJECT EXTERNAL_PROJECT_GEMS3K EXTERNAL_PROJECT_EIGEN)

#message(${CMAKE_CURRENT_BINARY_DIR}/../thirdparty/install/include/eigen3)
#include_directories(${CMAKE_CURRENT_BINARY_DIR}/../thirdparty/install/include/eigen3)

#target_link_libraries(XGEMS_OBJECT
#    PUBLIC ChemicalFun::ChemicalFun
#    PUBLIC ThermoFun::ThermoFun
#    PRIVATE spdlog::spdlog
#    PRIVATE nlohmann_json::nlohmann_json
#)

# Check if a shared library must be built
if(XGEMS_BUILD_SHARED_LIBS)
    add_library(XGEMS_SHARED SHARED $<TARGET_OBJECTS:XGEMS_OBJECT>)
    target_link_libraries(XGEMS_SHARED GEMS3K)
    set_target_properties(XGEMS_SHARED PROPERTIES OUTPUT_NAME xGEMS)
    install(TARGETS XGEMS_SHARED DESTINATION "lib" COMPONENT libraries)
endif()

# Check if a static library must be built
if(XGEMS_BUILD_STATIC_LIBS)
    add_library(XGEMS_STATIC STATIC $<TARGET_OBJECTS:XGEMS_OBJECT>)
    target_link_libraries(XGEMS_STATIC GEMS3K-static)
    set_target_properties(XGEMS_STATIC PROPERTIES OUTPUT_NAME xGEMS)
    install(TARGETS XGEMS_STATIC DESTINATION "lib" COMPONENT libraries)
endif()

# Create an install target for the header files
install(DIRECTORY ${XGEMS_SOURCE_DIR}
    DESTINATION "include" COMPONENT headers
    FILES_MATCHING PATTERN "*.hpp")
