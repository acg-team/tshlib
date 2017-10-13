cmake_minimum_required(VERSION 3.5)
project(tshexe)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-Wall -Wshadow -Weffc++")


if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(WINDOWS TRUE)
    message(STATUS "Compilation will be performed under Windows")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LINUX TRUE)
    message(STATUS "Compilation will be performed under Linux")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    message(STATUS "Compilation will be performed under Apple MacOS")
endif ()

if (UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")
endif (UNIX)

add_definitions(-DVERBOSE=0)


find_package("Eigen3")
if (EIGEN3_FOUND)
    message("Found Eigen3 in: ${EIGEN3_INCLUDE_DIR}")
endif (EIGEN3_FOUND)
include_directories(${EIGEN3_INCLUDE_DIR})

#add_subdirectory("lib/googletest-master")

add_subdirectory("src")
add_subdirectory("test")
include_directories("src")


set(SOURCES
        main.cpp
        main.hpp)

add_executable(tshexe ${SOURCES})
target_link_libraries(tshexe tshlib)
