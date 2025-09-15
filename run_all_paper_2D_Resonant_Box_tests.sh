#!/bin/bash

# CUDA NP 8 Resonant Box

# mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H1_P1/2D_Resonant_Box_H1_P1.json \
#   --device cuda

# mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H1_P2/2D_Resonant_Box_H1_P2.json \
#   --device cuda

# mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H1_P3/2D_Resonant_Box_H1_P3.json \
#   --device cuda

# mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P1/2D_Resonant_Box_H2_P1.json \
#   --device cuda

# mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P2/2D_Resonant_Box_H2_P2.json \
#   --device cuda

# mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P3/2D_Resonant_Box_H2_P3.json \
#   --device cuda

# mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H3_P1/2D_Resonant_Box_H3_P1.json \
#   --device cuda

# mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H3_P2/2D_Resonant_Box_H3_P2.json \
#   --device cuda

# mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H3_P3/2D_Resonant_Box_H3_P3.json \
#   --device cuda

# OPEN BASIS

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H1_P1/2D_Resonant_Box_BLEG_H1_P1.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H1_P2/2D_Resonant_Box_BLEG_H1_P2.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H1_P3/2D_Resonant_Box_BLEG_H1_P3.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H2_P1/2D_Resonant_Box_BLEG_H2_P1.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H2_P2/2D_Resonant_Box_BLEG_H2_P2.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H2_P3/2D_Resonant_Box_BLEG_H2_P3.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H3_P1/2D_Resonant_Box_BLEG_H3_P1.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H3_P2/2D_Resonant_Box_BLEG_H3_P2.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H3_P3/2D_Resonant_Box_BLEG_H3_P3.json \
  --device cuda

# MPI NP 8 Resonant Box
# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H1_P1/2D_Resonant_Box_H1_P1.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H1_P2/2D_Resonant_Box_H1_P2.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H1_P3/2D_Resonant_Box_H1_P3.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P1/2D_Resonant_Box_H2_P1.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P2/2D_Resonant_Box_H2_P2.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P3/2D_Resonant_Box_H2_P3.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H3_P1/2D_Resonant_Box_H3_P1.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H3_P2/2D_Resonant_Box_H3_P2.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H3_P3/2D_Resonant_Box_H3_P3.json \

# OPEN BASIS

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H1_P1/2D_Resonant_Box_BLEG_H1_P1.json \

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H1_P2/2D_Resonant_Box_BLEG_H1_P2.json \

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H1_P3/2D_Resonant_Box_BLEG_H1_P3.json \

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H2_P1/2D_Resonant_Box_BLEG_H2_P1.json \

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H2_P2/2D_Resonant_Box_BLEG_H2_P2.json \

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H2_P3/2D_Resonant_Box_BLEG_H2_P3.json \

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H3_P1/2D_Resonant_Box_BLEG_H3_P1.json \

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H3_P2/2D_Resonant_Box_BLEG_H3_P2.json \

mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H3_P3/2D_Resonant_Box_BLEG_H3_P3.json \


# Single Core Resonant Box
# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H1_P1/2D_Resonant_Box_H1_P1.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H1_P2/2D_Resonant_Box_H1_P2.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H1_P3/2D_Resonant_Box_H1_P3.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P1/2D_Resonant_Box_H2_P1.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P2/2D_Resonant_Box_H2_P2.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H2_P3/2D_Resonant_Box_H2_P3.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H3_P1/2D_Resonant_Box_H3_P1.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H3_P2/2D_Resonant_Box_H3_P2.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/2D_Resonant_Box/2D_Resonant_Box_H3_P3/2D_Resonant_Box_H3_P3.json \

# OPEN BASIS

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H1_P1/2D_Resonant_Box_BLEG_H1_P1.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H1_P2/2D_Resonant_Box_BLEG_H1_P2.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H1_P3/2D_Resonant_Box_BLEG_H1_P3.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H2_P1/2D_Resonant_Box_BLEG_H2_P1.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H2_P2/2D_Resonant_Box_BLEG_H2_P2.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H2_P3/2D_Resonant_Box_BLEG_H2_P3.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H3_P1/2D_Resonant_Box_BLEG_H3_P1.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H3_P2/2D_Resonant_Box_BLEG_H3_P2.json \

./build/gnu-release-mpi/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/2D_Resonant_Box_BLEG/2D_Resonant_Box_BLEG_H3_P3/2D_Resonant_Box_BLEG_H3_P3.json \

