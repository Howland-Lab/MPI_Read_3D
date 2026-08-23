#!/usr/bin/env bash

# Source the Intel/IMPI PadeOps stack on Anvil, then load the extra
# NetCDF Fortran module required by mpi_diag_lean.

export PADEOPS_ROOT="${PADEOPS_ROOT:-/anvil/projects/x-atm170028/karim/PadeOps}"

pushd "${PADEOPS_ROOT}" >/dev/null
source "${PADEOPS_ROOT}/setup/SetupEnv_Anvil_impi.sh"
popd >/dev/null

if [[ -n "${NETCDF_FORTRAN_ROOT:-}" ]]; then
  export PATH="${NETCDF_FORTRAN_ROOT}/bin:${PATH}"
  export LD_LIBRARY_PATH="${NETCDF_FORTRAN_ROOT}/lib:${LD_LIBRARY_PATH:-}"
else
  module load netcdf-fortran/4.5.3
fi
module list

if ! command -v nf-config >/dev/null 2>&1; then
  echo "ERROR: nf-config not found. Set NETCDF_FORTRAN_ROOT to a NetCDF Fortran install built with ${FC}." >&2
  return 1 2>/dev/null || exit 1
fi

netcdf_fc="$(nf-config --fc 2>/dev/null || true)"
case "${netcdf_fc}" in
  *ifx*|*ifort*|*mpiifx*|*mpiifort*)
    ;;
  *)
    echo "ERROR: incompatible NetCDF Fortran compiler: ${netcdf_fc:-unknown}" >&2
    echo "The Anvil netcdf-fortran/4.5.3 module is GCC-built and cannot be used with the PadeOps Intel/IMPI stack." >&2
    echo "Build or load NetCDF Fortran with ${FC}, then rerun with NETCDF_FORTRAN_ROOT=/path/to/netcdf-fortran." >&2
    return 1 2>/dev/null || exit 1
    ;;
esac
