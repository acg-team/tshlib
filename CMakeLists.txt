# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.


cmake_minimum_required(VERSION 3.8.2)

project(tshlib)
set(THSLIB_MAJOR_VERSION 0)
set(THSLIB_MINOR_VERSION 1)
set(THSLIB_PATCH_VERSION 0)
set(THSLIB_VERSION ${FOOBAR_MAJOR_VERSION}.${FOOBAR_MINOR_VERSION}.${FOOBAR_PATCH_VERSION})

#project(tshlib THSLIB_VERSION DESCRIPTION "Tree Search Heuristic Library")

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

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


message("Project source is in ${PROJECT_SOURCE_DIR}")
message("Project binary is in ${PROJECT_BINARY_DIR}")


# set up include-directories
include_directories(
        "${PROJECT_SOURCE_DIR}"   # to find
        "${PROJECT_BINARY_DIR}")  # to find



# The interesting stuff goes here
# ===============================

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
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    message(STATUS "Compilation will be performed under Apple MacOS")
    find_package(Boost 1.65.1 REQUIRED COMPONENTS date_time filesystem iostreams)
    find_package (glog 0.3.5 REQUIRED)
    find_package(gflags REQUIRED)
endif ()

find_package("Eigen3")
if (EIGEN3_FOUND)
    message("Found Eigen3 in: ${EIGEN3_INCLUDE_DIR}")
    include_directories(${EIGEN3_INCLUDE_DIR})
endif (EIGEN3_FOUND)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

file(GLOB SOURCES src/*.cpp )
file(GLOB HEADERS include/*.hpp )

add_library(tshlib SHARED ${SOURCES} ${HEADERS})
target_link_libraries(tshlib ${EIGEN3_INCLUDE_DIR})
target_link_libraries(tshlib glog)
target_link_libraries(tshlib gflags)

set_target_properties(tshlib PROPERTIES PUBLIC_HEADER "${HEADERS}")

install(TARGETS tshlib
        # IMPORTANT: Add the foo library to the "export-set"
        EXPORT TshLibTargets
        RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
        LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
        PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/tshlib"
        COMPONENT dev)

# Add all targets to the build-tree export set
export(TARGETS tshlib FILE "${PROJECT_BINARY_DIR}/TshLibTargets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE tshlib)

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

#set_target_properties(tshlib PROPERTIES
#        VERSION ${PROJECT_VERSION}
#        SOVERSION 1
#        PUBLIC_HEADER ${HEADERS})

