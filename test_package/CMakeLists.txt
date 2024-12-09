cmake_minimum_required(VERSION 3.15)
project(PackageTest CXX)

find_package(otfftpp CONFIG REQUIRED)



add_executable(otfftpp_test_package src/otfftpp_test_package.cpp)
target_link_libraries(otfftpp_test_package otfftpp::otfftpp)
