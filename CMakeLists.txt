cmake_minimum_required(VERSION 3.14)
project(CppNavierStokesSolver)

set(CMAKE_CXX_STANDARD 17)

# ADIOS2
set(CMAKE_PREFIX_PATH "${ADIOS2_DIR}" ${CMAKE_PREFIX_PATH})
find_package(adios2 REQUIRED)

# MGARD
set(MGARD_DIR "/home/adios/local/MGARD_install")
include_directories(${MGARD_DIR}/include)
link_directories(${MGARD_DIR}/lib)

# Set compiler flags based on build type
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Configuring Debug build")
    add_compile_options(-O0 -g -Wall -Wextra)
else()
    message(STATUS "Configuring Release build")
    add_compile_options(-O3 -march=native -ffast-math -fopenmp -fopenmp-simd -Wall -Wextra)
endif()

add_executable(main main.cpp)

target_link_libraries(main PRIVATE
    adios2::adios2
    mgard
)
