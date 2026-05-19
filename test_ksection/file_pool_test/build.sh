#!/bin/bash
# Build script for Phase D file_pool unit test.
# Uses ifx + Intel MPI + parallel HDF5 (built with mpiifx).
set -e

HDF5_ROOT=/home/kjhan/local/hdf5
FC=mpiifx
FCFLAGS="-O2 -fpp -DHDF5 -I${HDF5_ROOT}/include"
LDFLAGS="-L${HDF5_ROOT}/lib -Wl,-rpath,${HDF5_ROOT}/lib -lhdf5_fortran -lhdf5"

SRC_DIR=$(cd "$(dirname "$0")" && pwd)
RAMSES_IO=/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/patch/cuda/ramses_hdf5_io.f90

cd "${SRC_DIR}"
rm -f *.o *.mod file_pool_test

echo ">>> Compiling amr_parameters_stub.f90"
${FC} ${FCFLAGS} -c amr_parameters_stub.f90

echo ">>> Compiling ramses_hdf5_io.f90 (with -DHDF5)"
${FC} ${FCFLAGS} -c ${RAMSES_IO} -o ramses_hdf5_io.o

echo ">>> Compiling file_pool_test.f90"
${FC} ${FCFLAGS} -c file_pool_test.f90

echo ">>> Linking file_pool_test"
${FC} ${FCFLAGS} -o file_pool_test \
    file_pool_test.o ramses_hdf5_io.o amr_parameters_stub.o \
    ${LDFLAGS}

echo ">>> Build OK: $(ls -l file_pool_test | awk '{print $5, $9}')"
