#!/usr/bin/env bash
# Build the C++ dependencies of xgems as static libraries into $XGEMS_DEPS_PREFIX,
# so the Python wheel is self-contained (nothing to install besides the wheel).
# Called from pypi.yml / cibuildwheel (before-all) on Linux, macOS and Windows (Git Bash).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="${XGEMS_DEPS_PREFIX:?set XGEMS_DEPS_PREFIX}"
WORK="${XGEMS_DEPS_WORK:-$(mktemp -d)}"
JOBS="${JOBS:-4}"

# Pinned versions/branches (edit here)
EIGEN_TAG=3.4.0
EIGEN_URL="${EIGEN_URL:-https://gitlab.com/libeigen/eigen.git}"
JSON_TAG=v3.11.3
SPDLOG_TAG=v1.15.3   # >=1.15 needed: bundled fmt of older versions fails on recent Apple clang
CHEMICALFUN_REF=master
THERMOFUN_REF=master
GEMS3K_REF=master

mkdir -p "$PREFIX" "$WORK"
cd "$WORK"

# On Windows use the Visual Studio generator for x64 (same as scikit-build-core).
GEN_ARGS=()
if [[ "${OS:-}" == "Windows_NT" ]]; then
  GEN_ARGS=(-A x64)
fi

common=(
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_INSTALL_PREFIX="$PREFIX"
  -DCMAKE_PREFIX_PATH="$PREFIX"
  -DCMAKE_INSTALL_LIBDIR=lib
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON
)
# Header-only deps are found through include paths (CMake adds CXXFLAGS to its defaults).
export CXXFLAGS="-I$PREFIX/include -I$PREFIX/include/eigen3 ${CXXFLAGS:-}"
[[ -n "${MACOSX_DEPLOYMENT_TARGET:-}" ]] && common+=(-DCMAKE_OSX_DEPLOYMENT_TARGET="$MACOSX_DEPLOYMENT_TARGET")
[[ -n "${CMAKE_OSX_ARCHITECTURES:-}" ]] && common+=(-DCMAKE_OSX_ARCHITECTURES="$CMAKE_OSX_ARCHITECTURES")

build() {  # build <src-dir> [extra cmake args...]
  local src="$1"; shift
  cmake -S "$src" -B "$src/_build" ${GEN_ARGS[@]+"${GEN_ARGS[@]}"} "${common[@]}" "$@"
  cmake --build "$src/_build" --config Release --parallel "$JOBS"
  cmake --install "$src/_build" --config Release
}

fetch() {  # fetch <url> <ref> <dir>
  [[ -d "$3" ]] || git clone --depth 1 --branch "$2" "$1" "$3"
}

# ChemicalFun/ThermoFun/GEMS3K install $<TARGET_PDB_FILE:...> under MSVC, which CMake
# rejects for static libraries ("allowed only for targets with linker created artifacts").
# Drop those install() blocks; we never need the .pdb files.
strip_pdb_install() {  # strip_pdb_install <src-dir>
  find "$1" -name CMakeLists.txt -print0 | xargs -0 perl -0pi -e \
    's/install\s*\(\s*FILES\s+\$<TARGET_PDB_FILE:[^>]+>[^)]*\)//g'
}

# --- header-only libraries ---------------------------------------------------
fetch "$EIGEN_URL" "$EIGEN_TAG" eigen
build eigen -DBUILD_TESTING=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_PKGCONFIG=OFF

fetch https://github.com/nlohmann/json.git "$JSON_TAG" json
build json -DJSON_BuildTests=OFF -DJSON_MultipleHeaders=ON

fetch https://github.com/gabime/spdlog.git "$SPDLOG_TAG" spdlog
build spdlog -DSPDLOG_BUILD_EXAMPLE=OFF -DSPDLOG_BUILD_TESTS=OFF

# --- compiled libraries (static, no python bindings) -------------------------
fetch https://github.com/thermohub/chemicalfun.git "$CHEMICALFUN_REF" chemicalfun
strip_pdb_install chemicalfun
build chemicalfun -DBUILD_SHARED_LIBS=OFF -DCHEMICALFUN_BUILD_PYTHON=OFF \
  -DCHEMICALFUN_BUILD_TESTS=OFF -DCHEMICALFUN_BUILD_EXAMPLES=OFF

fetch https://github.com/thermohub/thermofun.git "$THERMOFUN_REF" thermofun
strip_pdb_install thermofun
build thermofun -DBUILD_SHARED_LIBS=OFF -DTFUN_BUILD_PYTHON=OFF -DTFUN_BUILD_TESTS=OFF

fetch https://github.com/gemshub/GEMS3K.git "$GEMS3K_REF" GEMS3K
strip_pdb_install GEMS3K
build GEMS3K -DCMAKE_PROJECT_INCLUDE="$HERE/preload-deps.cmake" -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DUSE_THERMOFUN=ON \
  -DBUILD_TOOLS=OFF -DBUILD_SOLMOD=OFF -DBUILD_SOLMOD_PYTHON=OFF

echo "Dependencies installed into $PREFIX"
