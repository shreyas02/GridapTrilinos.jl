#!/usr/bin/env bash
set -euo pipefail

BUILD_TYPE="${BUILD_TYPE:-Release}"
TRILINOS_SOURCE="${TRILINOS_SOURCE:-.ci/TrilinosSrc}"
TRILINOS_BUILD="${TRILINOS_BUILD:-.ci/TrilinosBuild}"
TRILINOS_INSTALL="${TRILINOS_INSTALL:-.ci/trilinos-install}"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_dir="${repo_root}/${TRILINOS_SOURCE}"
build_dir="${repo_root}/${TRILINOS_BUILD}"
install_dir="${repo_root}/${TRILINOS_INSTALL}"

if [[ ! -d "${source_dir}" ]]; then
  echo "Trilinos source directory not found: ${source_dir}" >&2
  exit 1
fi

mkdir -p "${build_dir}" "${install_dir}"

cmake -S "${source_dir}" -B "${build_dir}" -G Ninja \
  -D CMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -D CMAKE_INSTALL_PREFIX:PATH="${install_dir}" \
  -D BUILD_SHARED_LIBS:BOOL=ON \
  -D CMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_Fortran_COMPILER=mpif90 \
  -D Trilinos_ENABLE_Fortran:BOOL=ON \
  -D Trilinos_ENABLE_TESTS:BOOL=OFF \
  -D Trilinos_TRACE_ADD_TEST:BOOL=OFF \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
  -D Trilinos_ASSERT_MISSING_PACKAGES:BOOL=ON \
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
  -D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
  -D Trilinos_ENABLE_Amesos2:BOOL=ON \
  -D Trilinos_ENABLE_Belos:BOOL=ON \
  -D Trilinos_ENABLE_Galeri:BOOL=ON \
  -D Trilinos_ENABLE_MueLu:BOOL=ON \
  -D Trilinos_ENABLE_ShyLU_DD:BOOL=ON \
  -D Trilinos_ENABLE_ShyLU_DDFROSch:BOOL=ON \
  -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
  -D Trilinos_ENABLE_Thyra:BOOL=ON \
  -D Trilinos_ENABLE_Tpetra:BOOL=ON \
  -D Trilinos_ENABLE_Xpetra:BOOL=ON \
  -D Kokkos_ENABLE_SERIAL:BOOL=ON \
  -D Tpetra_ENABLE_DEPRECATED_CODE:BOOL=ON \
  -D Xpetra_ENABLE_DEPRECATED_CODE:BOOL=ON \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D TPL_ENABLE_BLAS:BOOL=ON \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
  -D TPL_ENABLE_UMFPACK:BOOL=ON \
  -D UMFPACK_INCLUDE_DIRS:PATH=/usr/include/suitesparse \
  -D UMFPACK_LIBRARY_DIRS:PATH=/usr/lib/x86_64-linux-gnu \
  -D TPL_ENABLE_Boost:BOOL=OFF \
  -D TPL_ENABLE_Netcdf:BOOL=OFF \
  -D TPL_ENABLE_X11:BOOL=OFF \
  -D Galeri_ENABLE_TESTS:BOOL=OFF \
  -D Kokkos_ENABLE_TESTS:BOOL=OFF \
  -D ShyLU_DDBDDC_ENABLE_TESTS:BOOL=OFF \
  -D ShyLU_DD_ENABLE_TESTS:BOOL=OFF

cmake --build "${build_dir}" --target install --parallel "$(nproc)"
