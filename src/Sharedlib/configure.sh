#!/bin/bash

# Adding the path to Trilinos 
if [ -z "$TRILINOS_ROOT" ]; then
    echo "Error: TRILINOS_ROOT environment variable is not set."
    echo "Please set it before running this script. For example:"
    echo "export TRILINOS_ROOT=/path/to/TrilinosInstall"
    exit 1
fi

echo "Using Trilinos installation at: $TRILINOS_ROOT"

SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SOURCE_DIR}/../.." && pwd)"
BUILD_DIR="${SOURCE_DIR}/build"

julia --project="${PROJECT_DIR}" -e '
using MPIPreferences
if MPIPreferences.binary != "system"
    error(
        "MPI.jl is configured to use $(MPIPreferences.binary), but GridapTrilinos " *
        "must use the same MPI as Trilinos. Run `MPIPreferences.use_system_binary(...)`, " *
        "restart Julia, and rebuild.",
    )
end
'

trilinos_lib_paths=()
[ -d "${TRILINOS_ROOT}/lib" ] && trilinos_lib_paths+=("${TRILINOS_ROOT}/lib")
[ -d "${TRILINOS_ROOT}/lib64" ] && trilinos_lib_paths+=("${TRILINOS_ROOT}/lib64")

if [ ${#trilinos_lib_paths[@]} -gt 0 ]; then
    trilinos_ld_path=$(IFS=:; echo "${trilinos_lib_paths[*]}")
    export LD_LIBRARY_PATH="${trilinos_ld_path}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
fi

export CMAKE_PREFIX_PATH="${TRILINOS_ROOT}${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}"

if [ -d "${TRILINOS_ROOT}/lib/cmake/Trilinos" ]; then
    export Trilinos_DIR="${TRILINOS_ROOT}/lib/cmake/Trilinos"
elif [ -d "${TRILINOS_ROOT}/lib64/cmake/Trilinos" ]; then
    export Trilinos_DIR="${TRILINOS_ROOT}/lib64/cmake/Trilinos"
fi

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}" || exit 1

cmake \
    "${SOURCE_DIR}"

cmake --build "${BUILD_DIR}"
