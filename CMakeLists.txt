cmake_minimum_required(VERSION 3.14)
project(aks_project)

set(CMAKE_CXX_STANDARD 17)

find_package(MPI REQUIRED)

add_executable(aks_project main.cpp inc/random.hpp lib/random.cpp inc/ga_int.hpp lib/ga_int.cpp inc/ga_double.hpp lib/ga_double.cpp
        inc/ga_double_bit.hpp lib/ga_double_bit.cpp inc/functions.hpp lib/functions.cpp)


find_package(Threads REQUIRED)
target_link_libraries(aks_project Threads::Threads)
target_include_directories(aks_project PRIVATE ${MPI_CXX_INCLUDE_PATH})
target_compile_options(aks_project PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(aks_project ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})