#!/bin/sh
# Put the xGEMS version from CMakeLists.txt into docs/rtd-environment.yml.
# Read the Docs runs this before it creates the conda environment.
set -e
XV=$(sed -n 's/^project(xGEMS VERSION \([0-9.]*\).*/\1/p' CMakeLists.txt)
[ -n "$XV" ] || { echo "xGEMS version not found in CMakeLists.txt"; exit 1; }
sed -i "s/@XGEMS_VERSION@/$XV/" docs/rtd-environment.yml
