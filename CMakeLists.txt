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
            inc
)

if(OTFFTPP_WITH_OPENMP)
    find_package (OpenMP REQUIRED)
    target_link_libraries (otfftpp PUBLIC OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(otfftpp PRIVATE "USE_AVX2" "USE_AVX" "USE_SSE2" "OTFFT_NAMESPACE=OTFFT_AVX2")
