#!/bin/bash
set -euo pipefail

if [ -z "${TRILINOS_ROOT:-}" ]; then
  echo "Error: TRILINOS_ROOT environment variable is not set."
  echo "Please set it before running this script. For example:"
  echo "export TRILINOS_ROOT=/path/to/TrilinosInstall"
  exit 1
fi

echo "Using Trilinos installation at: ${TRILINOS_ROOT}"

SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SOURCE_DIR}/../.." && pwd)"

julia --project="${PROJECT_DIR}" "${PROJECT_DIR}/deps/build.jl"
