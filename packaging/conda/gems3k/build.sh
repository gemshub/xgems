#!/bin/bash
set -eu
# 1. Optima from the gemshub fork (source folder optima-src), static, into a folder of this package.
OPT_SUB=lib/gems3k-optima
cmake -S optima-src -B build-optima -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="$PREFIX/$OPT_SUB" \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_PREFIX_PATH="$PREFIX" \
      -DBUILD_SHARED_LIBS=OFF \
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
      -DOPTIMA_BUILD_PYTHON=OFF \
      -DOPTIMA_BUILD_DEMOS=OFF \
      -DOPTIMA_BUILD_DOCS=OFF \
      -DOPTIMA_BUILD_BENCH=OFF
cmake --build build-optima -j"${CPU_COUNT}"
cmake --install build-optima
# 2. GEMS3K against that Optima; GEMS3KConfig.cmake then finds it for programs that use GEMS3K.
cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_PREFIX_PATH="$PREFIX/$OPT_SUB;$PREFIX" \
      -DOptima_DIR="$PREFIX/$OPT_SUB/lib/cmake/Optima" \
      -DGEMS3K_OPTIMA_SUBDIR="$OPT_SUB" \
      -DBUILD_SOLMOD_PYTHON=OFF \
      -DUSE_OPTIMA_SOLVER=ON
cmake --build build -j"${CPU_COUNT}"
cmake --install build
