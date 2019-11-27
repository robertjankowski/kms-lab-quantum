cmake_minimum_required(VERSION 3.15)
project(kms_lab_quantum)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic -O3")
set(source_dir "${PROJECT_SOURCE_DIR}/src/")

add_executable(quantum "src/main.cpp" "src/simulation.cpp")