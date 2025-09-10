#!/bin/bash

# CUDA NP 8 Resonant Box
mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H0_P2/3D_Resonant_Box_H0_P2.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H0_P3/3D_Resonant_Box_H0_P3.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H0_P4/3D_Resonant_Box_H0_P4.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P2/3D_Resonant_Box_BLEG_H0_P2.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P3/3D_Resonant_Box_BLEG_H0_P3.json \
  --device cuda
# Doesn't fit into VRAM memory
mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P4/3D_Resonant_Box_BLEG_H0_P4.json \
  --device cuda
# MPI NP 8 Resonant Box
mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H0_P2/3D_Resonant_Box_H0_P2.json \

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H0_P3/3D_Resonant_Box_H0_P3.json \

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H0_P4/3D_Resonant_Box_H0_P4.json \

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P2/3D_Resonant_Box_BLEG_H0_P2.json \

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P3/3D_Resonant_Box_BLEG_H0_P3.json \

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P4/3D_Resonant_Box_BLEG_H0_P4.json \

## Single Core Resonant Box
./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H0_P2/3D_Resonant_Box_H0_P2.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H0_P3/3D_Resonant_Box_H0_P3.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_H0_P4/3D_Resonant_Box_H0_P4.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P2/3D_Resonant_Box_BLEG_H0_P2.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P3/3D_Resonant_Box_BLEG_H0_P3.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box/3D_Resonant_Box_BLEG_H0_P4/3D_Resonant_Box_BLEG_H0_P4.json \