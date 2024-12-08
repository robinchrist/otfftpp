cmake_minimum_required(VERSION 3.27)
project(otfftpp CXX)

option(OTFFTPP_BUILD_TESTS "Build tests" OFF)
option(OTFFTPP_WITH_OPENMP "Use OpenMP" OFF)

add_library(
    otfftpp
)
add_library(otfftpp::otfftpp ALIAS otfftpp)

target_sources(
    otfftpp
    PRIVATE
        src/otfft.cpp
        src/otfft_algo.cpp
        src/otfft_algo_dif16.cpp
        src/otfft_algo_dif4.cpp
        src/otfft_algo_dif8.cpp
        src/otfft_algo_dit16.cpp
        src/otfft_algo_dit4.cpp
        src/otfft_algo_dit8.cpp
        src/otfft_algo_mixed.cpp
        src/otfft_algo_sixstep.cpp
    PUBLIC
        FILE_SET 
            otfftpp_public_headers
        TYPE
            HEADERS
        BASE_DIRS
            include
        FILES
            include/otfftpp/otfft_complex.h
            include/otfftpp/otfft_fwd.h
            include/otfftpp/otfft.h
    PRIVATE
        FILE_SET
            otfftpp_private_headers
        TYPE
            HEADERS
        BASE_DIRS
            include
)

if(OTFFTPP_WITH_OPENMP)
    find_package (OpenMP REQUIRED)
    target_link_libraries (otfftpp PUBLIC OpenMP::OpenMP_CXX)
else()
    target_compile_definitions(otfftpp PRIVATE "DO_SINGLE_THREAD=1")
endif()

# Required for legacy reasons
target_compile_definitions(otfftpp PRIVATE "USE_UNALIGNED_MEMORY=1")
# This should handle some edge cases that would otherwise prevent optimization
target_compile_options(otfftpp PRIVATE "-fno-math-errno")

if(OTFFT_BUILD_TESTS)
    message(STATUS "Testing is enabled!")
    enable_testing()
    add_subdirectory(tests)
endif()


install(
    TARGETS 
        otfftpp
    FILE_SET
        otfftpp_public_headers
)