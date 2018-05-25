# Copy resources dir for the demos to the binary build dir along with the demos
file(COPY resources DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

file(GLOB DEMO_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(DEMO_FILE ${DEMO_FILES})
    get_filename_component(DEMO_NAME ${DEMO_FILE} NAME_WE)
    add_executable(${DEMO_NAME} ${DEMO_FILE})
    target_link_libraries(${DEMO_NAME} XGEMS_STATIC GEMS3K-static)
endforeach()
