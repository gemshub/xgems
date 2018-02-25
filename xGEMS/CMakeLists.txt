# Compile the GEMS3K kernel code
add_subdirectory(GEMS3K)

# Recursively collect all header files from the current directory
file(GLOB HEADER_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp)

# Recursively collect all source files from the current directory
file(GLOB SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# Compile the C++ source code into object files
add_library(XGEMS_OBJECT OBJECT ${HEADER_FILES} ${SOURCE_FILES})

# Specify the GEMS3K dependency to ensure proper build sequence
add_dependencies(XGEMS_OBJECT GEMS3K)

# Check if a shared library must be built
if(XGEMS_BUILD_SHARED_LIBS)
    add_library(XGEMS_SHARED SHARED $<TARGET_OBJECTS:XGEMS_OBJECT>)
    target_link_libraries(XGEMS_SHARED ${GEMS3K})
    set_target_properties(XGEMS_SHARED PROPERTIES OUTPUT_NAME xGEMS)
    install(TARGETS XGEMS_SHARED DESTINATION "lib" COMPONENT libraries)
endif()

# Check if a static library must be built
if(XGEMS_BUILD_STATIC_LIBS)
    add_library(XGEMS_STATIC STATIC $<TARGET_OBJECTS:XGEMS_OBJECT>)
    target_link_libraries(XGEMS_STATIC ${GEMS3K})
    set_target_properties(XGEMS_STATIC PROPERTIES OUTPUT_NAME xGEMS)
    install(TARGETS XGEMS_STATIC DESTINATION "lib" COMPONENT libraries)
endif()

# Create an install target for the header files
install(DIRECTORY ${XGEMS_SOURCE_DIR}
    DESTINATION "include" COMPONENT headers
    FILES_MATCHING PATTERN "*.hpp")
