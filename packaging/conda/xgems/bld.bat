cmake -S . -B build -G Ninja ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
      -DCMAKE_PREFIX_PATH="%LIBRARY_PREFIX%" ^
      -DPYTHON_EXECUTABLE:FILEPATH="%PYTHON%" ^
      -DXGEMS_PYTHON_INSTALL_PREFIX="%PREFIX%" ^
      -DXGEMS_BUILD_SHARED_LIBS=OFF ^
      -DXGEMS_BUILD_DEMOS=OFF ^
      -DCMAKE_CXX_FLAGS="/utf-8"
if errorlevel 1 exit 1
cmake --build build
if errorlevel 1 exit 1
cmake --install build
if errorlevel 1 exit 1
