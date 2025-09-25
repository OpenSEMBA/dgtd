maxwellInputs#!/bin/bash

# CUDA NP 8 Resonant Box

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H0_P1/3D_Resonant_Box_TM55_H0_P1.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H0_P2/3D_Resonant_Box_TM55_H0_P2.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H0_P3/3D_Resonant_Box_TM55_H0_P3.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H1_P1/3D_Resonant_Box_TM55_H1_P1.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H1_P2/3D_Resonant_Box_TM55_H1_P2.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H1_P3/3D_Resonant_Box_TM55_H1_P3.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H2_P1/3D_Resonant_Box_TM55_H2_P1.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H2_P2/3D_Resonant_Box_TM55_H2_P2.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H2_P3/3D_Resonant_Box_TM55_H2_P3.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H3_P1/3D_Resonant_Box_TM55_H3_P1.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H4_P1/3D_Resonant_Box_TM55_H4_P1.json \
  --device cuda


mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H0_P1/3D_Resonant_Box_BLEG_TM55_H0_P1.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H0_P2/3D_Resonant_Box_BLEG_TM55_H0_P2.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H0_P3/3D_Resonant_Box_BLEG_TM55_H0_P3.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H1_P1/3D_Resonant_Box_BLEG_TM55_H1_P1.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H1_P2/3D_Resonant_Box_BLEG_TM55_H1_P2.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H1_P3/3D_Resonant_Box_BLEG_TM55_H1_P3.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H2_P1/3D_Resonant_Box_BLEG_TM55_H2_P1.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H2_P2/3D_Resonant_Box_BLEG_TM55_H2_P2.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H2_P3/3D_Resonant_Box_BLEG_TM55_H2_P3.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H3_P1/3D_Resonant_Box_BLEG_TM55_H3_P1.json \
  --device cuda

mpirun -np 8 ./build/gnu-release-cuda/bin/opensemba_dgtd \
  -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H4_P1/3D_Resonant_Box_BLEG_TM55_H4_P1.json \
  --device cuda

# MPI NP 8 Resonant Box
# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H0_P1/3D_Resonant_Box_TM55_H0_P1.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H0_P2/3D_Resonant_Box_TM55_H0_P2.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H0_P3/3D_Resonant_Box_TM55_H0_P3.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H1_P1/3D_Resonant_Box_TM55_H1_P1.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H1_P2/3D_Resonant_Box_TM55_H1_P2.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H1_P3/3D_Resonant_Box_TM55_H1_P3.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H2_P1/3D_Resonant_Box_TM55_H2_P1.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H2_P2/3D_Resonant_Box_TM55_H2_P2.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H2_P3/3D_Resonant_Box_TM55_H2_P3.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H0_P1/3D_Resonant_Box_BLEG_TM55_H0_P1.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H0_P2/3D_Resonant_Box_BLEG_TM55_H0_P2.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H0_P3/3D_Resonant_Box_BLEG_TM55_H0_P3.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H1_P1/3D_Resonant_Box_BLEG_TM55_H1_P1.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H1_P2/3D_Resonant_Box_BLEG_TM55_H1_P2.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H1_P3/3D_Resonant_Box_BLEG_TM55_H1_P3.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H2_P1/3D_Resonant_Box_BLEG_TM55_H2_P1.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H2_P2/3D_Resonant_Box_BLEG_TM55_H2_P2.json \

# mpirun -np 8 ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H2_P3/3D_Resonant_Box_BLEG_TM55_H2_P3.json \

# ## Single Core Resonant Box
# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H0_P1/3D_Resonant_Box_TM55_H0_P1.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H0_P2/3D_Resonant_Box_TM55_H0_P2.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H0_P3/3D_Resonant_Box_TM55_H0_P3.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H1_P1/3D_Resonant_Box_TM55_H1_P1.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H1_P2/3D_Resonant_Box_TM55_H1_P2.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H1_P3/3D_Resonant_Box_TM55_H1_P3.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H2_P1/3D_Resonant_Box_TM55_H2_P1.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H2_P2/3D_Resonant_Box_TM55_H2_P2.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_TM55/3D_Resonant_Box_TM55_H2_P3/3D_Resonant_Box_TM55_H2_P3.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H0_P1/3D_Resonant_Box_BLEG_TM55_H0_P1.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H0_P2/3D_Resonant_Box_BLEG_TM55_H0_P2.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H0_P3/3D_Resonant_Box_BLEG_TM55_H0_P3.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H1_P1/3D_Resonant_Box_BLEG_TM55_H1_P1.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H1_P2/3D_Resonant_Box_BLEG_TM55_H1_P2.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H1_P3/3D_Resonant_Box_BLEG_TM55_H1_P3.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H2_P1/3D_Resonant_Box_BLEG_TM55_H2_P1.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H2_P2/3D_Resonant_Box_BLEG_TM55_H2_P2.json \

# ./build/gnu-release-mpi/bin/opensemba_dgtd \
#   -i ./testData/maxwellInputs/3D_Resonant_Box_BLEG_TM55/3D_Resonant_Box_BLEG_TM55_H2_P3/3D_Resonant_Box_BLEG_TM55_H2_P3.json \
