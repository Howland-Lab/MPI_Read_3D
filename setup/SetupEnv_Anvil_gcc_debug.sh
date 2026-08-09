#!/usr/bin/env bash

# Source the GNU/OpenMPI PadeOps stack on Anvil, then load the extra
# NetCDF Fortran module required by mpi_diag_lean.

export PADEOPS_ROOT="${PADEOPS_ROOT:-/anvil/projects/x-atm170028/karim/PadeOps}"

pushd "${PADEOPS_ROOT}" >/dev/null
source "${PADEOPS_ROOT}/setup/SetupEnv_Anvil_gcc.sh"

# Override DECOMP_PATH to use the debug build of 2decomp_fft.
export DECOMP_PATH=/anvil/projects/x-atm170028/karim/PadeOps/dependencies/gcc-debug/2decomp_fft
popd >/dev/null

module load netcdf-fortran
module list
