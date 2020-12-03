cmake_minimum_required(VERSION 3.13)

project(tonberry VERSION 1.1)

# Uses C++ 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmakefiles")

add_executable(tonberry "")

# Source code
add_subdirectory(src)
add_subdirectory(src/mesh)
add_subdirectory(src/hydro)

find_package("NetCDF" REQUIRED)
find_package("SILO" REQUIRED)
find_package("Iris" REQUIRED)
find_package("log4cxx" REQUIRED)

target_link_libraries(tonberry netcdf)
target_link_libraries(tonberry silo)
target_link_libraries(tonberry iris)
target_link_libraries(tonberry log4cxx)
