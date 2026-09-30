cmake_minimum_required(VERSION 3.25.2)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(dgtd)

set(CMAKE_CXX_STANDARD 17)

option(SEMBA_DGTD_ENABLE_MPI "Enable MPI support" OFF )
option(SEMBA_DGTD_ENABLE_OPENMP "Enable OpenMP support" OFF )
option(SEMBA_DGTD_ENABLE_CUDA "Enable CUDA support" OFF )
option(SEMBA_DGTD_ENABLE_PYTHON_BINDINGS "Enable Python bindings" OFF )
option(SEMBA_DGTD_ENABLE_MFEM_AS_SUBDIRECTORY "Use MFEM as a subdirectory" ON )

option(SEMBA_DGTD_ENABLE_EXTENSIVE_CASE_TESTS "Enable extensive case tests" ON )
option(SEMBA_DGTD_ENABLE_EXTENSIVE_SOLVER_TESTS "Enable extensive solver tests" ON )
option(SEMBA_DGTD_ENABLE_EXTENSIVE_RCS_TESTS "Enable extensive RCS tests" ON )
option(SEMBA_DGTD_ENABLE_L2_ANALYSIS_TESTS "Enable L2 Analysis Tests" OFF)

option(SEMBA_DGTD_ENABLE_TIMER_INFORMATION "Enable timer information" ON)

include_directories(
	src/ 
	src/maxwell/
)

# MFEM enables CUDA in add_subdirectory(); arch must be fixed in cache first
# (stale CMAKE_CUDA_ARCHITECTURES from a prior machine/toolkit can break nvcc).
if (SEMBA_DGTD_ENABLE_CUDA)
	if (CUDA_ARCH)
		string(REGEX REPLACE "^sm_" "" _dgtd_cuda_arch_num "${CUDA_ARCH}")
		if (NOT CUDA_ARCH STREQUAL "sm_${_dgtd_cuda_arch_num}")
			message(FATAL_ERROR "Unknown CUDA_ARCH: ${CUDA_ARCH} (expected sm_XX)")
		endif()
		set(CMAKE_CUDA_ARCHITECTURES "${_dgtd_cuda_arch_num}" CACHE STRING
			"CUDA architectures to build for." FORCE)
	elseif (NOT CMAKE_CUDA_ARCHITECTURES)
		message(FATAL_ERROR
			"SEMBA_DGTD_ENABLE_CUDA=ON requires CUDA_ARCH (e.g. sm_89) or CMAKE_CUDA_ARCHITECTURES")
	endif()
	message(STATUS "CUDA arch for MFEM + dgtd: sm_${CMAKE_CUDA_ARCHITECTURES}")
endif()

if (SEMBA_DGTD_ENABLE_MFEM_AS_SUBDIRECTORY)
	add_subdirectory(external/mfem-geg)
	# Blackwell (sm_120): cicc can take hours/OOM on heavily templated MFEM
	# EA/TMOP/LOR TUs at -O3 (mfem/mfem#5363). Keep arch=120 but drop opt on
	# those sources only; other CUDA TUs stay at Release -O3.
	if (SEMBA_DGTD_ENABLE_CUDA AND "${CMAKE_CUDA_ARCHITECTURES}" STREQUAL "120")
		set(_dgtd_mfem_dir "${CMAKE_SOURCE_DIR}/external/mfem-geg")
		# MFEM 4.9 moved TMOP PA kernels under fem/tmop/{pa,assemble,mult,tools,metrics}/
		file(GLOB _dgtd_mfem_cuda_o0
			"${_dgtd_mfem_dir}/fem/integ/bilininteg_*_ea.cpp"
			"${_dgtd_mfem_dir}/fem/integ/bilininteg_*_kernels.cpp"
			"${_dgtd_mfem_dir}/fem/lor/lor_batched.cpp"
			"${_dgtd_mfem_dir}/fem/dgmassinv.cpp"
			"${_dgtd_mfem_dir}/fem/tmop.cpp"
			"${_dgtd_mfem_dir}/fem/tmop_tools.cpp"
			"${_dgtd_mfem_dir}/fem/tmop_amr.cpp"
			"${_dgtd_mfem_dir}/fem/qinterp/*.cpp"
		)
		file(GLOB_RECURSE _dgtd_mfem_tmop_cuda_o0
			"${_dgtd_mfem_dir}/fem/tmop/*.cpp"
		)
		list(APPEND _dgtd_mfem_cuda_o0 ${_dgtd_mfem_tmop_cuda_o0})
		# Properties must be set in mfem's directory scope (CMake 3.18+).
		set_source_files_properties(${_dgtd_mfem_cuda_o0}
			DIRECTORY "${_dgtd_mfem_dir}"
			TARGET_DIRECTORY mfem
			PROPERTIES COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:CUDA>:-O0>")
		list(LENGTH _dgtd_mfem_cuda_o0 _dgtd_mfem_cuda_o0_n)
		message(STATUS
			"CUDA sm_120: compiling ${_dgtd_mfem_cuda_o0_n} heavy MFEM TUs at -O0 (cicc workaround)")
	endif()
else ()
	find_package(MFEM CONFIG REQUIRED)
	include_directories(${MFEM_INCLUDE_DIRS})
endif()

if(SEMBA_DGTD_ENABLE_PYTHON_BINDINGS)
	add_subdirectory(pythonBindings)
endif()

macro(files_to_cuda EXE_SRCS)
    foreach(SRC_FILE IN LISTS ${EXE_SRCS})
        if (SEMBA_DGTD_ENABLE_CUDA)
            set_source_files_properties(${SRC_FILE} PROPERTIES LANGUAGE CUDA)
        endif()
    endforeach()
endmacro()

if(SEMBA_DGTD_ENABLE_CUDA)
    message(STATUS "Project files will be compiled with CUDA support")
    file(GLOB_RECURSE SRC_FILES "src/*.cpp" "test/*.cpp")
    files_to_cuda(SRC_FILES)
    
    if (NOT CMAKE_CUDA_HOST_COMPILER)
        set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    endif()
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD} CACHE STRING
        "CUDA standard to use.")
    set(CMAKE_CUDA_STANDARD_REQUIRED ON CACHE BOOL
        "Force the use of the chosen CUDA standard.")
    set(CMAKE_CUDA_EXTENSIONS OFF CACHE BOOL "Enable CUDA standard extensions.")
    set(CUDA_FLAGS "--expt-extended-lambda")
    if (CMAKE_VERSION VERSION_LESS 3.18.0)
        set(CUDA_FLAGS "-arch=${CUDA_ARCH} ${CUDA_FLAGS}")
    elseif (NOT CMAKE_CUDA_ARCHITECTURES)
        string(REGEX REPLACE "^sm_" "" ARCH_NUMBER "${CUDA_ARCH}")
        if ("${CUDA_ARCH}" STREQUAL "sm_${ARCH_NUMBER}")
        set(CMAKE_CUDA_ARCHITECTURES "${ARCH_NUMBER}")
        else()
        message(FATAL_ERROR "Unknown CUDA_ARCH: ${CUDA_ARCH}")
        endif()
    else()
        set(CUDA_ARCH "CMAKE_CUDA_ARCHITECTURES: ${CMAKE_CUDA_ARCHITECTURES}")
    endif()
    message(STATUS "Using CUDA architecture: ${CUDA_ARCH}")
    if (CMAKE_VERSION VERSION_LESS 3.12.0)
        # CMake versions 3.8 and 3.9 require this to work; 3.10 and 3.11 are not
        # tested and may not actually need this (but should be ok to keep).
        set(CUDA_FLAGS "-ccbin=${CMAKE_CXX_COMPILER} ${CUDA_FLAGS}")
        set(CMAKE_CUDA_HOST_LINK_LAUNCHER ${CMAKE_CXX_COMPILER})
    endif()
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${CUDA_FLAGS}")
    find_package(CUDAToolkit REQUIRED)
    set(CUSPARSE_FOUND TRUE)
    get_target_property(CUSPARSE_LIBRARIES CUDA::cusparse LOCATION)

endif()

add_subdirectory(src)

enable_testing()

add_subdirectory(test/)

