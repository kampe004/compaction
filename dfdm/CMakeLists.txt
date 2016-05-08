cmake_minimum_required(VERSION 2.8)
project(dfdm)
add_definitions(-Wall)

SET( GCC_NETCDF_INCLUDE "-I/opt/local/include")
SET( GCC_NETCDF_LIBRARY "-L/opt/local/lib -lnetcdf_c++4")
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_NETCDF_LIBRARY}" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_NETCDF_INCLUDE} -std=c++14")

#set(CMAKE_BUILD_TYPE debug)
# or add to cmake command: -DCMAKE_BUILD_TYPE=Debug
set(CMAKE_C_FLAGS_DEBUG "-g -O0 ")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 ")

add_subdirectory(src)
add_executable(dfdm.exe ${SOURCE})
