cmake_minimum_required(VERSION 3.15)
project(aks_project)

set(CMAKE_CXX_STANDARD 17)

add_executable(aks_project main.cpp)

find_package(Threads REQUIRED)
target_link_libraries(aks_project Threads::Threads)