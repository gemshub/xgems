#!/bin/bash
set -eu
cmake -S . -B build -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCMAKE_PREFIX_PATH="$PREFIX" \
      -DPYTHON_EXECUTABLE:FILEPATH="$PYTHON" \
      -DPython3_EXECUTABLE:FILEPATH="$PYTHON" \
      -DPython3_ROOT_DIR="$PREFIX" \
      -DPython3_FIND_FRAMEWORK=NEVER \
      -DPython3_FIND_STRATEGY=LOCATION \
      -DXGEMS_PYTHON_INSTALL_PREFIX="$PREFIX" \
      -DXGEMS_BUILD_SHARED_LIBS=OFF \
      -DXGEMS_BUILD_DEMOS=OFF
cmake --build build -j"${CPU_COUNT}"
cmake --install build
