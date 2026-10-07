file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/gems3k" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/test1.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/test_systems.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/test_solver_modes.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/test_logging.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

enable_testing()

add_test(NAME test_systems
COMMAND ${Python3_EXECUTABLE} "-m" "unittest" "-v" test_systems.py)

add_test(NAME test_solver_modes
COMMAND ${Python3_EXECUTABLE} "-m" "unittest" "-v" test_solver_modes.py)

add_test(NAME test_logging
COMMAND ${Python3_EXECUTABLE} "-m" "unittest" "-v" test_logging.py)
