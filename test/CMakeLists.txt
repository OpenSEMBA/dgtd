# Find dependencies
find_package(GTest CONFIG REQUIRED)
find_package(MPI REQUIRED)

include_directories(
    ./
    ${maxwell_INCLUDE_DIRS}
    ${GTEST_INCLUDE_DIRS}
)

add_library(TestLauncher OBJECT TestLauncher.cpp)
target_include_directories(TestLauncher PRIVATE ${MFEM_INCLUDE_DIRS})
target_link_libraries(TestLauncher PRIVATE mfem)

add_subdirectory(mfem)
add_subdirectory(hesthavenComparison)
add_subdirectory(maxwell)
add_subdirectory(cases)
add_subdirectory(rcs)

add_test(NAME mfem
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUM_MPI_PROCS} ${CMAKE_BINARY_DIR}/bin/mfem_tests)

add_test(NAME hesthavenComparison
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUM_MPI_PROCS} ${CMAKE_BINARY_DIR}/bin/hesthavenComparison_tests)

add_test(NAME maxwell
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUM_MPI_PROCS} ${CMAKE_BINARY_DIR}/bin/maxwell_tests)

add_test(NAME rcs
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUM_MPI_PROCS} ${CMAKE_BINARY_DIR}/bin/rcs_tests)

add_test(NAME maxwell_solver
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUM_MPI_PROCS} ${CMAKE_BINARY_DIR}/bin/maxwell_solver_tests)

add_test(NAME cases
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUM_MPI_PROCS} ${CMAKE_BINARY_DIR}/bin/cases_tests)
