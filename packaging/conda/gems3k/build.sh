#!/bin/bash
set -eu
cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_PREFIX_PATH="$PREFIX" \
      -DBUILD_SOLMOD_PYTHON=OFF \
      -DUSE_OPTIMA_SOLVER=ON
cmake --build build -j"${CPU_COUNT}"
cmake --install build
