#!/bin/bash

# This script is run by mpirun on each rank, launching gdbserver on a unique port

# MPI rank from OpenMPI
RANK=${OMPI_COMM_WORLD_RANK:-0}
PORT=$((20000 + RANK))  # unique port per rank

echo "Rank $RANK launching gdbserver on port $PORT"

# Start gdbserver, waiting for VS Code debugger to attach before running your test with filter
exec gdbserver :$PORT ./build/gnu-debug-parallel/bin/mfem_tests --gtest_filter=ParallelSpaces
