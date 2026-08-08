#!/usr/bin/env bash

module purge
module load StdEnv
module load gcc/12.2.0
module load community-modules
module load cmake/3.27.9
module load openmpi/4.1.4
module load netcdf-c
module load netcdf-fortran
module load hdf5

module list

export COMPILER_ID=GNU
export FC=mpifort
export CC=mpicc
export CXX=mpicxx

export DECOMP_DIR="/home/karimali/work/PadeOps/dependencies/2decomp_fft_gcc"
export DECOMP2D_INC="${DECOMP_DIR}/include"
export DECOMP2D_LIB="${DECOMP_DIR}/lib"
export source="mpi_diag.F90"
export prog="MPIR3D"

NF_PREFIX=$(nf-config --prefix)
NC_PREFIX=$(nc-config --prefix)

${FC} -O3 -g -fopenmp -ffree-line-length-none \
  -I${DECOMP2D_INC} $(nf-config --fflags) \
  ${source} \
  -L${DECOMP2D_LIB} -l2decomp_fft \
  -L${NF_PREFIX}/lib -L${NF_PREFIX}/lib64 \
  -L${NC_PREFIX}/lib -L${NC_PREFIX}/lib64 \
  $(nf-config --flibs) \
  -o ${prog}