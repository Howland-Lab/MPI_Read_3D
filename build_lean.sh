#!/usr/bin/env bash

# Unified builder for mpi_diag_lean.F90.
#
# Usage:
#   ./build_lean.sh anvil gcc
#   ./build_lean.sh anvil impi
#   ./build_lean.sh archer2
#   ./build_lean.sh stampede3 gcc
#   ./build_lean.sh stampede3 impi
#
# The setup files load the platform/toolchain modules and define FC,
# COMPILER_ID, and DECOMP_PATH.  NetCDF is part of the setup because
# netcdf.mod must be compatible with the active Fortran compiler.

set -euo pipefail

platform="${1:-}"
toolchain="${2:-}"

if [[ -z "${platform}" ]]; then
  echo "Usage: $0 <anvil|archer2|stampede3> [gcc|impi]" >&2
  exit 2
fi

case "${platform}:${toolchain}" in
  anvil:gcc)
    setup_file="setup/SetupEnv_Anvil_gcc.sh"
    ;;
  anvil:gcc_debug)
    setup_file="setup/SetupEnv_Anvil_gcc_debug.sh"
    ;;
  anvil:impi)
    setup_file="setup/SetupEnv_Anvil_impi.sh"
    ;;
  archer2:|archer2:gcc)
    setup_file="setup/SetupEnv_Archer2.sh"
    ;;
  stampede3:gcc)
    setup_file="setup/SetupEnv_Stampede3_gcc.sh"
    ;;
  stampede3:impi)
    setup_file="setup/SetupEnv_Stampede3_impi.sh"
    ;;
  *)
    echo "Unsupported platform/toolchain: ${platform} ${toolchain}" >&2
    echo "Valid examples: anvil gcc, anvil impi, archer2, stampede3 gcc, stampede3 impi" >&2
    exit 2
    ;;
esac

source "${setup_file}"

if [[ -z "${FC:-}" ]]; then
  echo "ERROR: setup file did not define FC." >&2
  exit 3
fi

if [[ -z "${DECOMP_PATH:-}" ]]; then
  echo "ERROR: setup file did not define DECOMP_PATH." >&2
  exit 3
fi

if [[ ! -d "${DECOMP_PATH}/include" || ! -d "${DECOMP_PATH}/lib" ]]; then
  echo "ERROR: DECOMP_PATH does not contain include/ and lib/: ${DECOMP_PATH}" >&2
  exit 3
fi

if ! command -v nf-config >/dev/null 2>&1; then
  echo "ERROR: nf-config not found after loading ${setup_file}." >&2
  echo "Load a NetCDF Fortran module compatible with this compiler/MPI stack or update the setup file." >&2
  exit 3
fi

source_file="${SOURCE_FILE:-mpi_diag_lean.F90}"
case "${platform}" in
  archer2)
    default_build_dir="build"
    ;;
  *)
    default_build_dir="build/${toolchain}"
    ;;
esac
build_dir="${BUILD_DIR:-${default_build_dir}}"
program="${PROGRAM:-MPIR3D_lean}"
program_path="${build_dir}/${program}"

mkdir -p "${build_dir}"

common_flags=("-O3" "-g")

case "${COMPILER_ID:-}" in
  GNU)
    fcflags=("${common_flags[@]}" "-Wall" "-fopenmp" "-ffree-line-length-none" "-fallow-argument-mismatch")
    mod_flags=("-J${build_dir}")
    ;;
  Intel)
    fcflags=("${common_flags[@]}" "-traceback" "-warn" "all" "-qopenmp")
    mod_flags=("-module" "${build_dir}")
    ;;
  *)
    fcflags=("${common_flags[@]}")
    mod_flags=("-J${build_dir}")
    ;;
esac

"${FC}" "${fcflags[@]}" \
  "${mod_flags[@]}" \
  -I"${DECOMP_PATH}/include" $(nf-config --fflags) \
  "${source_file}" \
  -L"${DECOMP_PATH}/lib" -l2decomp_fft \
  $(nf-config --flibs) \
  -o "${program_path}"

echo "Built ${program_path} using ${setup_file}"
