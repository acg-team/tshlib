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
project(tshlib)

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


#SET(EIGEN_INCLUDE_DIR "/usr/include/eigen3" CACHE STRING "path to libEigen3 headers")
find_package("Eigen3")
if (EIGEN3_FOUND)
    message("Found Eigen3 in: ${EIGEN3_INCLUDE_DIR}")
endif (EIGEN3_FOUND)
################################################################################
#                 _                        _   _                 
#   ___ __ _  ___| |__   ___    ___  _ __ | |_(_) ___  _ __  ___ 
#  / __/ _` |/ __| '_ \ / _ \  / _ \| '_ \| __| |/ _ \| '_ \/ __|
# | (_| (_| | (__| | | |  __/ | (_) | |_) | |_| | (_) | | | \__ \
#  \___\__,_|\___|_| |_|\___|  \___/| .__/ \__|_|\___/|_| |_|___/
#                                   |_|                          
################################################################################

#SET(WITH_DEBUG OFF CACHE BOOL "enable debug output")
#SET(WITH_SSE ON CACHE BOOL "enable code using SSE2 instructions")
#SET(WITH_SSE OFF CACHE BOOL "enable code using SSE2 instructions")
#SET(WITH_LS_LOG ON CACHE BOOL "use fast approximation of logarithm")
#SET(WITH_CODON ON CACHE BOOL "enable codon alignment")
#SET(WITH_DNA ON CACHE BOOL "enable DNA alignment")
#SET(WITH_AA ON CACHE BOOL "enable AA alignment")

#temporary flags
#SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")
#set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
#SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#SET(BUILD_SHARED_LIBRARIES OFF)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")

if (UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")
endif ()

################################################################################
#                                       ___     _                       _       
#  ___  ___  _   _ _ __ ___ ___  ___   ( _ )   | |_ __ _ _ __ __ _  ___| |_ ___ 
# / __|/ _ \| | | | '__/ __/ _ \/ __|  / _ \/\ | __/ _` | '__/ _` |/ _ \ __/ __|
# \__ \ (_) | |_| | | | (_|  __/\__ \ | (_>  < | || (_| | | | (_| |  __/ |_\__ \
# |___/\___/ \__,_|_|  \___\___||___/  \___/\/  \__\__,_|_|  \__, |\___|\__|___/
#                                                            |___/              
################################################################################

INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(SOURCES
        main.cpp
        newick.cpp
        PhyTree.cpp)

ADD_EXECUTABLE(tshlib ${SOURCES})
target_link_libraries(tshlib ${EIGEN3_INCLUDE_DIR})

################################################################################
#                                                              
#  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___  ___  _ __ 
# | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|/ _ \| '__|
# | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \ (_) | |   
# | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/\___/|_|   
# |_|            |_|                                           
#      _       __ _       _ _   _                 
#   __| | ___ / _(_)_ __ (_) |_(_) ___  _ __  ___ 
#  / _` |/ _ \ |_| | '_ \| | __| |/ _ \| '_ \/ __|
# | (_| |  __/  _| | | | | | |_| | (_) | | | \__ \
#  \__,_|\___|_| |_|_| |_|_|\__|_|\___/|_| |_|___/
#                                                 
################################################################################

#IF(WITH_DEBUG)
#   ADD_DEFINITIONS(-DDEBUG)
#ENDIF(WITH_DEBUG)

#IF(WITH_SSE)
#   #TODO icc does only recognize -msse3
#   ADD_DEFINITIONS(-msse2)
#ENDIF(WITH_SSE)

#IF(WITH_LS_LOG)
#   ADD_DEFINITIONS(-DUSE_LS_LOG)
#ENDIF(WITH_LS_LOG)


################################################################################
#  _           _        _ _       _   _             
# (_)_ __  ___| |_ __ _| | | __ _| |_(_) ___  _ __  
# | | '_ \/ __| __/ _` | | |/ _` | __| |/ _ \| '_ \ 
# | | | | \__ \ || (_| | | | (_| | |_| | (_) | | | |
# |_|_| |_|___/\__\__,_|_|_|\__,_|\__|_|\___/|_| |_|
#                                                   
################################################################################

#INSTALL(TARGETS NNI_SPR DESTINATION bin)

#SET(CPACK_GENERATOR "STGZ;TGZ;TBZ2;DEB")
#SET(CPACK_PACKAGE_CONTACT "mmaiolo@gmail.com")
#INCLUDE(CPack)
