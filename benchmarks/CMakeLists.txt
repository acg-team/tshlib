cmake_minimum_required(VERSION 3.5)
project(tshlib_benchmarks)

set(CMAKE_CXX_STANDARD 11)

SET(CMAKE_CXX_FLAGS_DEBUG  "-Wall -Wshadow -Weffc++ -O0 -g")
SET(CMAKE_C_FLAGS_DEBUG "-Wall -Wshadow -O0 -g")

SET(CMAKE_CXX_FLAGS_RELEASE  "-O3")
SET(CMAKE_C_FLAGS_RELEASE "-O3")

SET(CMAKE_CXX_FLAGS_INTEL  "-O3")
SET(CMAKE_C_FLAGS_INTEL "-O3")

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(WINDOWS TRUE)
    message(STATUS "Compilation will be performed under Windows")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LINUX TRUE)
    message(STATUS "Compilation will be performed under Linux")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    message(STATUS "Compilation will be performed under Apple MacOS")
    find_package (glog 0.3.5 REQUIRED)
    find_package(gflags REQUIRED)

endif ()

if (UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lstdc++ -ldl -pthread ")
endif (UNIX)


find_package("Eigen3")
if (EIGEN3_FOUND)
    message("Found Eigen3 in: ${EIGEN3_INCLUDE_DIR}")
endif (EIGEN3_FOUND)
include_directories(${EIGEN3_INCLUDE_DIR})


include_directories("../src")

set(SOURCES
        main.cpp
        main.hpp)

add_executable(tshlib_benchmarks ${SOURCES})
target_link_libraries(tshlib_benchmarks tshlib)
target_link_libraries(tshlib_benchmarks glog)
target_link_libraries(tshlib_benchmarks gflags)


