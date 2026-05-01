cmake_minimum_required(VERSION 3.15)

project(AMBER_ff LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/d72f9c8aea6817cdf1ca0ac10887f328de7f3da2.zip
)
FetchContent_MakeAvailable(googletest)

find_package(Armadillo REQUIRED)

add_library(amber_ff_lib 
    src/molecule.cpp
    src/read_xyz.cpp
    src/read_params.cpp
    src/ff.cpp
    src/util.cpp
    src/atom_typing.cpp
)
target_include_directories(amber_ff_lib PUBLIC src ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(amber_ff_lib PUBLIC ${ARMADILLO_LIBRARIES})

add_executable(amber_ff src/main.cpp)
target_link_libraries(amber_ff PRIVATE amber_ff_lib)

# Test Cases
enable_testing()
include(GoogleTest)

add_executable(test_hello tests/test_hello.cpp)
target_link_libraries(test_hello PRIVATE amber_ff_lib GTest::gtest_main)
gtest_discover_tests(test_hello)

add_executable(test_xyz tests/test_xyz.cpp)
target_link_libraries(test_xyz PRIVATE amber_ff_lib GTest::gtest_main)
gtest_discover_tests(test_xyz)

add_executable(test_params tests/test_params.cpp)
target_link_libraries(test_params PRIVATE amber_ff_lib GTest::gtest_main)
gtest_discover_tests(test_params)

add_executable(test_calc_ff tests/test_calc_ff.cpp)
target_link_libraries(test_calc_ff PRIVATE amber_ff_lib GTest::gtest_main)
gtest_discover_tests(test_calc_ff)

