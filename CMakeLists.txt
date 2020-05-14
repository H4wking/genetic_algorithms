cmake_minimum_required(VERSION 3.14)
project(aks_project)

set(CMAKE_CXX_STANDARD 17)

add_executable(aks_project main.cpp inc/random.hpp lib/random.cpp inc/ga_int.hpp lib/ga_int.cpp inc/ga_double.hpp lib/ga_double.cpp
        inc/ga_double_bit.hpp lib/ga_double_bit.cpp)


find_package(Threads REQUIRED)
target_link_libraries(aks_project Threads::Threads)