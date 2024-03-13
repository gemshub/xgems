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

# Check if a shared library must be built
if(XGEMS_BUILD_SHARED_LIBS)
    add_library(XGEMS_SHARED SHARED $<TARGET_OBJECTS:XGEMS_OBJECT>)
    target_link_libraries(XGEMS_SHARED PUBLIC GEMS3K::GEMS3K)
    set_target_properties(XGEMS_SHARED PROPERTIES OUTPUT_NAME xGEMS)
    install(TARGETS XGEMS_SHARED DESTINATION "lib" COMPONENT libraries)
endif()

# Check if a static library must be built
if(XGEMS_BUILD_STATIC_LIBS)
    add_library(XGEMS_STATIC STATIC $<TARGET_OBJECTS:XGEMS_OBJECT>)
    target_link_libraries(XGEMS_STATIC PUBLIC GEMS3K::GEMS3K-static)
    set_target_properties(XGEMS_STATIC PROPERTIES OUTPUT_NAME xGEMS)
    install(TARGETS XGEMS_STATIC DESTINATION "lib" COMPONENT libraries)
endif()


# Create an install target for the header files
install(DIRECTORY ${XGEMS_SOURCE_DIR}
    DESTINATION "include" COMPONENT headers
    FILES_MATCHING PATTERN "*.hpp")
