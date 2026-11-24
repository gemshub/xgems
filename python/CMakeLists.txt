message(STATUS "Building python")

# Configure the pyproject.toml script with the provided CMake variables
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml.in ${CMAKE_CURRENT_BINARY_DIR}/pyproject.toml)

add_subdirectory(src/xgems)

# Building a wheel with scikit-build-core (pyproject.toml at the repo root):
# install the module and python files straight into the wheel and skip pip below.
if(SKBUILD)
    install(TARGETS PyxGEMS LIBRARY DESTINATION xgems COMPONENT python
                            RUNTIME DESTINATION xgems COMPONENT python)
    install(FILES src/xgems/__init__.py src/xgems/PyxGEMS.pyi src/xgems/py.typed
            DESTINATION xgems COMPONENT python)
    return()
endif()

# Install where the python package is installed to CMAKE_INSTALL_PREFIX if not given
if(NOT DEFINED XGEMS_PYTHON_INSTALL_PREFIX)

# Install the xgems python package using setuptools
install(CODE
"
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PyxGEMS>
                ${CMAKE_CURRENT_BINARY_DIR}/src/xgems/$<TARGET_FILE_NAME:PyxGEMS>
      )

   execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pip install ${CMAKE_CURRENT_BINARY_DIR} -vv
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
"
)

else()

# If the path is already in Windows format (with backslashes), it can't be added directly
# to the string below, otherwise CMake will later complain about "Invalid escape sequence".
file(TO_CMAKE_PATH "${XGEMS_PYTHON_INSTALL_PREFIX}" XGEMS_PYTHON_INSTALL_PREFIX)
message(STATUS "XGEMS_PYTHON_INSTALL_PREFIX: ${XGEMS_PYTHON_INSTALL_PREFIX}")

# Install the xgems python package using setuptools
install(CODE
"
   file(TO_NATIVE_PATH \"${XGEMS_PYTHON_INSTALL_PREFIX}\" XGEMS_PYTHON_INSTALL_PREFIX_NATIVE)

    execute_process(
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:PyxGEMS>
                ${CMAKE_CURRENT_BINARY_DIR}/src/xgems/$<TARGET_FILE_NAME:PyxGEMS>
        )

    execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pip install --prefix \${XGEMS_PYTHON_INSTALL_PREFIX_NATIVE} ${CMAKE_CURRENT_BINARY_DIR} -vv
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
"
)

endif()

