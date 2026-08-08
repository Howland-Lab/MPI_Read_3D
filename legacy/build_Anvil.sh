#!/usr/bin/env bash

# --- Modules ---
module purge 
module load intel/19.0.5.281
module load mvapich2/2.3.6
module load netcdf-c
module load netcdf-fortran
module load zlib
module load hdf5
module list  

# --- Compilers (use Cray wrappers) ---
export COMPILER_ID=Intel
export FC=mpiifort 
export CC=mpiicc
export CXX=mpiicpc

export DECOMP_DIR="/anvil/projects/x-atm170028/padeops_setup/dependencies/2decomp_fft"
export DECOMP2D_INC="${DECOMP_DIR}/include"
export DECOMP2D_LIB="${DECOMP_DIR}/lib"
export source="mpi_diag.F90"
export prog="MPIR3D"

# INTEL
# ${FC} -O3 -traceback -g -warn all -qopenmp -I${DECOMP2D_INC} ${source} -L${DECOMP2D_LIB} -l2decomp_fft -o ${prog}

${FC} -O3 -traceback -g -warn all -qopenmp -I${DECOMP2D_INC} $(nf-config --fflags) \
  ${source} \
  -L${DECOMP2D_LIB} -l2decomp_fft \
  $(nf-config --flibs) \
  -o ${prog}