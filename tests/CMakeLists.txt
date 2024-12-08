#
# otfftpp unit tests

find_package(Boost
    REQUIRED 
    COMPONENTS
        unit_test_framework
)

add_executable(otfftpp-tests)

target_sources(otfftpp-tests
    PRIVATE
        test_module.cpp
        otfft_basic_test.cpp
        otfft_bluestein_test.cpp
        otfft_dct_test.cpp
        otfft_find_peak_test.cpp
        otfft_optimization_test.cpp
        otfft_transform_test.cpp
)

target_link_libraries(otfftpp-tests
    PRIVATE
        otfftpp::otfftpp
        Boost::unit_test_framework
)


add_test(NAME otfftpp-tests COMMAND otfftpp-tests)
