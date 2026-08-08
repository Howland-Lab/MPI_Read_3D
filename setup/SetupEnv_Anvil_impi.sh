#!/usr/bin/env bash

# Source the Intel/IMPI PadeOps stack on Anvil, then load the extra
# NetCDF Fortran module required by mpi_diag_lean.

export PADEOPS_ROOT="${PADEOPS_ROOT:-/anvil/projects/x-atm170028/karim/PadeOps}"

pushd "${PADEOPS_ROOT}" >/dev/null
source "${PADEOPS_ROOT}/setup/SetupEnv_Anvil_impi.sh"
popd >/dev/null

module load netcdf-fortran
module list
