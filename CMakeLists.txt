cmake_minimum_required(VERSION 3.15)

project(AMBER_ff LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Armadillo REQUIRED)

add_executable(amber_ff
    src/main.cpp
    src/molecule.cpp
    src/read_xyz.cpp
    src/read_params.cpp
    src/ff.cpp
    src/util.cpp
)

target_include_directories(amber_ff PRIVATE src)
target_link_libraries(amber_ff PRIVATE ${ARMADILLO_LIBRARIES})
target_include_directories(amber_ff PRIVATE ${ARMADILLO_INCLUDE_DIRS})
