message(STATUS "Creating build system for maxwell-evolution")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

set(SOURCES
    "MaxwellEvolution.cpp"
    "MaxwellEvolutionMethods.cpp"
    "HesthavenEvolutionMethods.cpp"
    "HesthavenEvolution.cpp"
    "GlobalEvolution.cpp"
)

if (SEMBA_DGTD_ENABLE_CUDA)
    enable_language(CUDA)

    list(APPEND SOURCES "GlobalEvolution.cu")

    set_source_files_properties(GlobalEvolution.cu PROPERTIES LANGUAGE CUDA)

    message(STATUS "CUDA is enabled: Compiling with GlobalEvolution.cu")
else()
    message(STATUS "CUDA is disabled: Skipping GlobalEvolution.cu")
endif()


add_library(maxwell-evolution STATIC ${SOURCES})

if (SEMBA_DGTD_ENABLE_CUDA)
    set_target_properties(maxwell-evolution PROPERTIES
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
    )
endif()

get_filename_component(PARENT_DIR ../ ABSOLUTE)
target_include_directories(maxwell-evolution PRIVATE ${PARENT_DIR})

target_link_libraries(maxwell-evolution 
    maxwell-components 
    mfem 
    maxwell-mfemExtension 
    Eigen3::Eigen
)

if (SEMBA_DGTD_ENABLE_OPENMP)
    message(STATUS "Linking maxwell-evolution with OpenMP libraries")
    target_link_libraries(maxwell-evolution OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()