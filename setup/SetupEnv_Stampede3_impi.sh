#!/usr/bin/env bash

# Source the Intel/IMPI PadeOps stack on Stampede3, then load the extra
# NetCDF modules used by the MPI_Read_3D diagnostics.

export PADEOPS_ROOT="${PADEOPS_ROOT:-/work2/10829/kali/stampede3/PadeOps}"

pushd "${PADEOPS_ROOT}" >/dev/null
source "${PADEOPS_ROOT}/setup/SetupEnv_Stampede3_impi.sh"
popd >/dev/null

module load netcdf
module load zlib
module load hdf5
module list
