cmake_minimum_required(VERSION 3.14)
project(CppNavierStokesSolver)

set(CMAKE_CXX_STANDARD 17)

# ADIOS2
set(CMAKE_PREFIX_PATH "${ADIOS2_DIR}" ${CMAKE_PREFIX_PATH})
find_package(adios2 REQUIRED)


# MGARD (manually set up include & lib paths)
set(MGARD_DIR "/home/adios/local/MGARD_install")
include_directories(${MGARD_DIR}/include)
link_directories(${MGARD_DIR}/lib)

add_executable(main main.cpp)

target_link_libraries(main PRIVATE
  adios2::adios2
  mgard
)

