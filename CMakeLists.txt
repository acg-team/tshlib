# Licensed Materials - Property of Lorenzo Gatti & Massimo Maiolo
# Copyright (C) 2015-2019 by Lorenzo Gatti & Massimo Maiolo
# -------------------------------------------------------------------------
# This file is part of Tree Search Heuristics Library (TSH-LIB)
#
# TSH-LIB is a C++ Library whose purpose is generate alternative tree
# toplogies applying NNI, SPR, and F-VFNI rearrangements on an initial tree
# #
# TSH-LIB is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.
#
# TSH-LIB is a free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# You should have received a copy of the GNU Lesser General Public
# License along with TSH-LIB. If not, see <http://www.gnu.org/licenses/>.
# -------------------------------------------------------------------------


cmake_minimum_required(VERSION 3.5)

set(TSHLIB_SOFTWARENAME "tshlib")
set(TSHLIB_DESCRIPTION "Tree Search Heuristic Library")

project(${TSHLIB_SOFTWARENAME}
        DESCRIPTION ${TSHLIB_DESCRIPTION}
        LANGUAGES "CXX")

## Store the git hash of the current head
if(EXISTS "${PROJECT_SOURCE_DIR}/.git/HEAD")
    file(READ "${PROJECT_SOURCE_DIR}/.git/HEAD"
            PROJECT_SOURCE_VERSION)
    if("${PROJECT_SOURCE_VERSION}" MATCHES "^ref:")
        string(REGEX REPLACE "^ref: *([^ \n\r]*).*" "\\1"
                PROJECT_GIT_REF "${PROJECT_SOURCE_VERSION}")
        file(READ "${PROJECT_SOURCE_DIR}/.git/${PROJECT_GIT_REF}" PROJECT_SOURCE_VERSION)
    endif()
    string(STRIP "${PROJECT_SOURCE_VERSION}" PROJECT_SOURCE_VERSION)
endif()

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

# Store the build date
if(WIN32)
    execute_process(COMMAND "cmd" " /c date /t"
            OUTPUT_VARIABLE DATE)
    string(REGEX REPLACE "[^0-9]*(..).*" "\\1" MONTH "${DATE}")
    set(MONTHS ""
            "Jan" "Feb" "Mar" "Apr" "May" "Jun"
            "Jul" "Aug" "Sep" "Oct" "Nov" "Dec")
    list(GET MONTHS "${MONTH}" MONTH)
    string(REGEX REPLACE "[^/]*/(..)/(....).*" "\\1 ${MONTH} \\2" PROJECT_BUILD_DATE "${DATE}")
    execute_process(COMMAND "cmd" " /c echo %TIME%" OUTPUT_VARIABLE TIME)
    string(REGEX REPLACE "[^0-9]*(..:..:..).*" "\\1" PROJECT_BUILD_TIME "${TIME}")
    execute_process(COMMAND git describe --abbrev=0 --tags
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_TAG_VERSION)
    string(REGEX REPLACE "\n" "" GIT_TAG_VERSION_STRIPPED "${GIT_TAG_VERSION}")
else()
    execute_process(COMMAND "date" "+%d %b %Y/%H:%M:%S" OUTPUT_VARIABLE DATE_TIME)
    string(REGEX REPLACE "([^/]*)/.*" "\\1" PROJECT_BUILD_DATE "${DATE_TIME}")
    string(REGEX REPLACE "[^/]*/([0-9:]*).*" "\\1" PROJECT_BUILD_TIME "${DATE_TIME}")
    execute_process(COMMAND git describe --abbrev=0 --tags
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_TAG_VERSION)
    string(REGEX REPLACE "\n" "" GIT_TAG_VERSION_STRIPPED "${GIT_TAG_VERSION}")
endif()

set(TSHLIB_RELTIME ${PROJECT_BUILD_TIME})
set(TSHLIB_RELDATE ${PROJECT_BUILD_DATE})
set(TSHLIB_VERSION ${GIT_TAG_VERSION_STRIPPED})

# Set library package version from git repository release tag version
set(TSHLIB_PACKAGE_VERSION ${TSHLIB_VERSION})
string(REPLACE "." ";" TSHLIB_PACKAGE_VERSION_LIST ${TSHLIB_PACKAGE_VERSION})
list (GET TSHLIB_PACKAGE_VERSION_LIST 0 MAJOR_PACKAGE_VERSION)

project(${TSHLIB_SOFTWARENAME} VERSION ${TSHLIB_VERSION})

add_definitions(-DPRJ_GITBRANCH=\"${PROJECT_GIT_REF}\"
        -DPRJ_VERSION=\"${MINIJATI_VERSION}\"
        -DPRJ_GITREF=\"${PROJECT_SOURCE_VERSION}\"
        -DPRJ_DESC=\"${MINIJATI_DESCRIPTION}\"
        -DPRJ_NAME=\"${MINIJATI_SOFTWARENAME}\"
        -DPRJ_DATE=\"${MINIJATI_RELDATE}\"
        -DPRJ_TIME=\"${MINIJATI_RELTIME}\" )

message("-- Compilation will be performed with the following release of the software:
\tbranch   ${PROJECT_GIT_REF}
\tref      ${PROJECT_SOURCE_VERSION}
\ttime     ${TSHLIB_RELDATE} ${TSHLIB_RELTIME}
\tversion  ${TSHLIB_VERSION}")


set(CMAKE_CXX_STANDARD 11)


SET(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wshadow -O0 -g -DNDEBUG=1 ")
SET(CMAKE_C_FLAGS_DEBUG "-Wall -Wshadow -O0 -g -DNDEBUG=1 ")

SET(CMAKE_CXX_FLAGS_RELEASE "-O3  -g0 -DNDEBUG=0 ")
SET(CMAKE_C_FLAGS_RELEASE "-O3  -g0 -DNDEBUG=0 ")

SET(CMAKE_CXX_FLAGS_INTEL "-O3  -g0 -DNDEBUG=0")
SET(CMAKE_C_FLAGS_INTEL "-O3  -g0 -DNDEBUG=0")

if(WIN32 AND NOT CYGWIN)
    set(DEF_INSTALL_CMAKE_DIR CMake)
else()
    set(DEF_INSTALL_CMAKE_DIR lib/tshlib)
endif()

set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH  "Installation directory for CMake files")


# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
    set(var INSTALL_${p}_DIR)
    if(NOT IS_ABSOLUTE "${${var}}")
        set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
    endif()
endforeach()


message(STATUS "Project source is in ${PROJECT_SOURCE_DIR}")
message(STATUS "Project binary is in ${PROJECT_BINARY_DIR}")

# set up include-directories
include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")

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

# Dependencies not covered by find package should be found in the following directories
#if (${CMAKE_PREFIX_PATH})
#    include_directories("${CMAKE_PREFIX_PATH}/include")
#    LINK_DIRECTORIES("${CMAKE_PREFIX_PATH}}/lib")
#    LINK_DIRECTORIES("${CMAKE_PREFIX_PATH}/lib64")
#    message(STATUS "Looking for libraries in the following directory: ${CMAKE_PREFIX_PATH}/lib")
#    message(STATUS "Looking for libraries in the following directory: ${CMAKE_PREFIX_PATH}/lib64")
#    message(STATUS "Looking for headers in the following directory: ${CMAKE_PREFIX_PATH}/include")
#endif ()

find_package(glog 0.3.5 REQUIRED)
if(GLOG_FOUND)
    message(STATUS "Found GLOG in ${GLOG_INCLUDE_DIRS}")
    include_directories(${GLOG_INCLUDE_DIRS})
endif()

set(Boost_USE_STATIC_LIBS        OFF) # only find static libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost REQUIRED COMPONENTS system unit_test_framework)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    message("\t static lib:  ${Boost_LIBRARY_DIRS}" )
    message("\t include dir: ${Boost_INCLUDE_DIRS}" )
endif()


# Add local include folder
include_directories("${PROJECT_SOURCE_DIR}/include")

# Files to compile
file(GLOB SOURCES src/*.cpp )
file(GLOB HEADERS include/*.hpp )

# Create the shared library
add_library(tshlib SHARED ${SOURCES} ${HEADERS})
target_link_libraries(tshlib glog::glog)
set_target_properties(tshlib PROPERTIES PUBLIC_HEADER "${HEADERS}")
set_target_properties(tshlib PROPERTIES VERSION ${GIT_TAG_VERSION_STRIPPED} SOVERSION ${MAJOR_PACKAGE_VERSION})


# Create the static library
add_library(tshlib-static STATIC ${SOURCES} ${HEADERS})
target_link_libraries(tshlib-static glog::glog)
set_target_properties(tshlib-static PROPERTIES PUBLIC_HEADER "${HEADERS}")
set_target_properties(tshlib-static PROPERTIES OUTPUT_NAME tshlib)
set_target_properties(tshlib-static PROPERTIES VERSION ${GIT_TAG_VERSION_STRIPPED} SOVERSION ${MAJOR_PACKAGE_VERSION})

# Test target
enable_testing()
add_subdirectory (tests)

# Install shared library
install(TARGETS tshlib
        # IMPORTANT: Add the tshlib library to the "export-set"
        EXPORT TshLibTargets
        RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
        LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
        PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/tshlib"
        COMPONENT dev)

# Install static library
install(TARGETS tshlib-static
        EXPORT  TshLibTargets
        RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
        LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
        ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT dev
        PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/tshlib"
        COMPONENT dev)

# Add all targets to the build-tree export set
export(TARGETS tshlib FILE "${PROJECT_BINARY_DIR}/TshLibTargets.cmake")
export(TARGETS tshlib-static FILE "${PROJECT_BINARY_DIR}/TshLibTargets.cmake")
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE tshlib)
export(PACKAGE tshlib-static)

# Create the FooBarConfig.cmake and FooBarConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
#set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
set(CONF_INCLUDE_DIRS "${INSTALL_INCLUDE_DIR}/tshlib")
configure_file(TshLibConfig.cmake.in "${PROJECT_BINARY_DIR}/TshLibConfig.cmake" @ONLY)
# ... for the install tree
#set(CONF_INCLUDE_DIRS "\${TSHLIB_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(TshLibConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TshLibConfig.cmake" @ONLY)
# ... for both
configure_file(TshLibConfigVersion.cmake.in "${PROJECT_BINARY_DIR}/TshLibConfigVersion.cmake" @ONLY)

# Install the FooBarConfig.cmake and FooBarConfigVersion.cmake
install(FILES
        "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TshLibConfig.cmake"
        "${PROJECT_BINARY_DIR}/TshLibConfigVersion.cmake"
        DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# Install the export set for use with the install-tree
install(EXPORT TshLibTargets DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# Add target for CLION installation
add_custom_target(install_${PROJECT_NAME}
        "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target install
        DEPENDS ${PROJECT_NAME}
        COMMENT "Installing ${PROJECT_NAME}")
