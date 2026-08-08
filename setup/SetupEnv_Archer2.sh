#!/usr/bin/env bash

# Source the Archer2 PadeOps stack, then load the extra NetCDF module
# required by mpi_diag_lean.

export PADEOPS_ROOT="${PADEOPS_ROOT:-/mnt/lustre/a2fs-work3/work/e773/e773/pounds/PadeOps}"

pushd "${PADEOPS_ROOT}" >/dev/null
source "${PADEOPS_ROOT}/setup/SetupEnv_Archer2.sh"
popd >/dev/null

module load cray-netcdf-hdf5parallel/4.9.0.7
module list
