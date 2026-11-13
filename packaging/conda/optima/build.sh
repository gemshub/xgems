#!/bin/bash
set -eu
cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DBUILD_SHARED_LIBS=OFF \
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
      -DOPTIMA_BUILD_PYTHON=OFF \
      -DOPTIMA_BUILD_DEMOS=OFF \
      -DOPTIMA_BUILD_DOCS=OFF \
      -DOPTIMA_BUILD_BENCH=OFF
cmake --build build -j"${CPU_COUNT}"
cmake --install build
