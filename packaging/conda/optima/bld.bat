cmake -S . -B build -G Ninja ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
      -DBUILD_SHARED_LIBS=OFF ^
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON ^
      -DOPTIMA_BUILD_PYTHON=OFF ^
      -DOPTIMA_BUILD_DEMOS=OFF ^
      -DOPTIMA_BUILD_DOCS=OFF ^
      -DOPTIMA_BUILD_BENCH=OFF
if errorlevel 1 exit 1
cmake --build build
if errorlevel 1 exit 1
cmake --install build
if errorlevel 1 exit 1
