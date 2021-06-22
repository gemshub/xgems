# Copy the xGEMS directory to the binary directory
file(COPY . DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Define a variable with the name of the compiled python module file as defined by pybind11 (which is used in the file setup.py.in)
set(XGEMS_PYTHON_MODULE_FILENAME PyxGEMS${PYTHON_MODULE_EXTENSION})

# Configure the setup.py file
set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
set(SETUP_PY ${CMAKE_CURRENT_BINARY_DIR}/setup.py)

# Configure the setup.py script with the provided CMake variables
configure_file(${SETUP_PY_IN} ${SETUP_PY})

# Create a custom target to build xGEMS python package during build stage
add_custom_target(xgems ALL
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PyxGEMS> 
        ${CMAKE_CURRENT_BINARY_DIR}/$<TARGET_FILE_NAME:PyxGEMS>
    COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} --quiet install --prefix=${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Set dependencies of xGEMS target
add_dependencies(xgems PyxGEMS)

# Set the path where the python package is installed to CMAKE_INSTALL_PREFIX if not given
if(NOT DEFINED XGEMS_PYTHON_INSTALL_PREFIX)
    set(XGEMS_PYTHON_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()
# If the path is already in Windows format (with backslashes), it can't be added directly
# to the string below, otherwise CMake will later complain about "Invalid escape sequence".
file(TO_CMAKE_PATH "${XGEMS_PYTHON_INSTALL_PREFIX}" XGEMS_PYTHON_INSTALL_PREFIX)

# Install the xGEMS python package using setuptools
install(CODE
"
    file(TO_NATIVE_PATH \"${XGEMS_PYTHON_INSTALL_PREFIX}\" XGEMS_PYTHON_INSTALL_PREFIX_NATIVE)

    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/../../lib/PyxGEMS.pdb)
        string(REPLACE .pyd .pdb XGEMS_PDB_FILENAME \"${XGEMS_PYTHON_MODULE_FILENAME}\")

        execute_process(
            COMMAND \${CMAKE_COMMAND} -E copy ../../lib/PyxGEMS.pdb \${XGEMS_PDB_FILENAME}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif()

    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=\${XGEMS_PYTHON_INSTALL_PREFIX_NATIVE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
"
)
