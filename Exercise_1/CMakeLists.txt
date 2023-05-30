# complie mpic++ -std=c++17 -O3 -Wall -pedantic -march=native -ffast-math main.cpp -o jacobiMPI
cmake_minimum_required(VERSION 3.12)

project(mpiJacobi
        VERSION 0.1
        LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS true)

set(CMAKE_CXX_STANDARD 17)

add_executable(mpiJacobi src/jacobiMPI.cpp)

target_include_directories(mpiJacobi PRIVATE include)