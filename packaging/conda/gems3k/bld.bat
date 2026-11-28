:: 1. Optima from the gemshub fork (source folder optima-src), static, into a folder of this package.
set OPT_SUB=lib/gems3k-optima
cmake -S optima-src -B build-optima -G Ninja ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%/%OPT_SUB%" ^
      -DCMAKE_INSTALL_LIBDIR=lib ^
      -DCMAKE_PREFIX_PATH="%LIBRARY_PREFIX%" ^
      -DBUILD_SHARED_LIBS=OFF ^
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON ^
      -DOPTIMA_BUILD_PYTHON=OFF ^
      -DOPTIMA_BUILD_DEMOS=OFF ^
      -DOPTIMA_BUILD_DOCS=OFF ^
      -DOPTIMA_BUILD_BENCH=OFF
if errorlevel 1 exit 1
cmake --build build-optima
if errorlevel 1 exit 1
cmake --install build-optima
if errorlevel 1 exit 1
:: 2. GEMS3K against that Optima; GEMS3KConfig.cmake then finds it for programs that use GEMS3K.
cmake -S . -B build -G Ninja ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
      -DCMAKE_INSTALL_LIBDIR=lib ^
      -DCMAKE_PREFIX_PATH="%LIBRARY_PREFIX%/%OPT_SUB%;%LIBRARY_PREFIX%" ^
      -DOptima_DIR="%LIBRARY_PREFIX%/%OPT_SUB%/lib/cmake/Optima" ^
      -DGEMS3K_OPTIMA_SUBDIR="%OPT_SUB%" ^
      -DBUILD_SOLMOD_PYTHON=OFF ^
      -DCMAKE_CXX_FLAGS="/utf-8" ^
      -DUSE_OPTIMA_SOLVER=ON
if errorlevel 1 exit 1
cmake --build build
if errorlevel 1 exit 1
cmake --install build
if errorlevel 1 exit 1
