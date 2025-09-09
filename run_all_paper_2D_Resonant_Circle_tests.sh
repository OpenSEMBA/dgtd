#!/bin/bash

# CUDA NP 8 Resonant Circle
#mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_H1_P2/2D_Resonant_Circle_H1_P2.json \
#  --device cuda
#
#mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_H1_P3/2D_Resonant_Circle_H1_P3.json \
#  --device cuda
#
#mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_H1_P4/2D_Resonant_Circle_H1_P4.json \
#  --device cuda
#
#mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_BLEG_H1_P2/2D_Resonant_Circle_BLEG_H1_P2.json \
#  --device cuda
#
#mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_BLEG_H1_P3/2D_Resonant_Circle_BLEG_H1_P3.json \
#  --device cuda
#
#mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_BLEG_H1_P4/2D_Resonant_Circle_BLEG_H1_P4.json \
#  --device cuda

# MPI NP 8 Resonant Circle
#mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_H1_P2/2D_Resonant_Circle_H1_P2.json \
#
#mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_H1_P3/2D_Resonant_Circle_H1_P3.json \
#
#mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_H1_P4/2D_Resonant_Circle_H1_P4.json \
#
#mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_BLEG_H1_P2/2D_Resonant_Circle_BLEG_H1_P2.json \
#
#mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_BLEG_H1_P3/2D_Resonant_Circle_BLEG_H1_P3.json \
#
#mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_BLEG_H1_P4/2D_Resonant_Circle_BLEG_H1_P4.json \

# Single Core Resonant Circle
./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_H1_P2/2D_Resonant_Circle_H1_P2.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_H1_P3/2D_Resonant_Circle_H1_P3.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_H1_P4/2D_Resonant_Circle_H1_P4.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_BLEG_H1_P2/2D_Resonant_Circle_BLEG_H1_P2.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_BLEG_H1_P3/2D_Resonant_Circle_BLEG_H1_P3.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Circle/2D_Resonant_Circle_BLEG_H1_P4/2D_Resonant_Circle_BLEG_H1_P4.json \