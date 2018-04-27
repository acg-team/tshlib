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

project(${TSHLIB_SOFTWARENAME} VERSION ${GIT_TAG_VERSION_STRIPPED})

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
\ttime     ${MINIJATI_RELDATE} ${MINIJATI_RELTIME}
\tcurrent  ${MINIJATI_VERSION} (latest version)")


set(CMAKE_CXX_STANDARD 11)

SET(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wshadow -Weffc++ -O0 -g -DNDEBUG=1")
SET(CMAKE_C_FLAGS_DEBUG "-Wall -Wshadow -O0 -g -DNDEBUG=1")

SET(CMAKE_CXX_FLAGS_RELEASE "-O3  -g0 -DNDEBUG=0")
SET(CMAKE_C_FLAGS_RELEASE "-O3  -g0 -DNDEBUG=0")

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
if (${CMAKE_PREFIX_PATH})
    include_directories("${CMAKE_PREFIX_PATH}/include")
    LINK_DIRECTORIES("${CMAKE_PREFIX_PATH}}/lib")
    LINK_DIRECTORIES("${CMAKE_PREFIX_PATH}/lib64")
    message(STATUS "Looking for libraries in the following directory: ${CMAKE_PREFIX_PATH}/lib")
    message(STATUS "Looking for libraries in the following directory: ${CMAKE_PREFIX_PATH}/lib64")
    message(STATUS "Looking for headers in the following directory: ${CMAKE_PREFIX_PATH}/include")
endif ()

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

find_package(Eigen3 REQUIRED)
if (EIGEN3_FOUND)
    message(STATUS "Found Eigen3")
        message("\t include dir: ${EIGEN3_INCLUDE_DIR}")
    include_directories(${EIGEN3_INCLUDE_DIR})
endif ()


# Add local include folder
include_directories("${PROJECT_SOURCE_DIR}/include")

# Files to compile
file(GLOB SOURCES src/*.cpp )
file(GLOB HEADERS include/*.hpp )

# Create the shared library
add_library(tshlib SHARED ${SOURCES} ${HEADERS})
target_link_libraries(tshlib glog::glog)
set_target_properties(tshlib PROPERTIES PUBLIC_HEADER "${HEADERS}")

# Create the static library
add_library(tshlib-static STATIC ${SOURCES} ${HEADERS})
target_link_libraries(tshlib-static glog::glog)
set_target_properties(tshlib-static PROPERTIES PUBLIC_HEADER "${HEADERS}")
set_target_properties(tshlib-static PROPERTIES OUTPUT_NAME tshlib)

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
