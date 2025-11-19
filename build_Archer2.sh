#!/usr/bin/env bash

# --- Modules ---
module purge
module load PrgEnv-gnu
module load craype-x86-rome  

# --- Compilers (use Cray wrappers) ---
export COMPILER_ID=GNU
export CC=cc
export CXX=CC
export FC=ftn

export DECOMP_DIR="/mnt/lustre/a2fs-work3/work/e773/e773/pounds/PadeOps/dependencies/2decomp_fft"
export DECOMP2D_INC="${DECOMP_DIR}/include"
export DECOMP2D_LIB="${DECOMP_DIR}/lib"
export source="mpi_diag.F90"
export prog="MPIR3D"

# GNU + OpenMPI
mpif90 -O3 -fbacktrace -g -Wall -fopenmp -I${DECOMP2D_INC} ${source} -L${DECOMP2D_LIB} -l2decomp_fft -o ${prog}