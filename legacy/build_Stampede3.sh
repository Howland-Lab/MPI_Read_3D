#!/usr/bin/env bash

# --- Modules ---
module purge
module load intel impi netcdf zlib hdf5
module list

# --- Compilers (use Cray wrappers) ---
export COMPILER_ID=Intel
export FC=mpiifort 
export CC=mpiicc
export CXX=mpiicpc

export DECOMP_DIR="/work2/10829/kali/stampede3/PadeOps/dependencies/2decomp_fft"
export DECOMP2D_INC="${DECOMP_DIR}/include"
export DECOMP2D_LIB="${DECOMP_DIR}/lib"
export source="mpi_diag.F90"
export prog="MPIR3D"

${FC} -O3 -traceback -g -warn all -qopenmp -I${DECOMP2D_INC} $(nf-config --fflags) \
  ${source} \
  -L${DECOMP2D_LIB} -l2decomp_fft \
  $(nf-config --flibs) \
  -o ${prog}