cmake_minimum_required(VERSION 3.5)
project(tshlib_tests)

include(CTest)
enable_testing()

#set(gtest_disable_pthreads on) #needed in MinGW


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    include(../lib/DP/DownloadProject.cmake)
    download_project(
            PROJ googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG master
            UPDATE_DISCONNECTED 1)

    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()


# Include TSHLIB
include_directories(../src)

#set(test_SRCS ADD ALL TEST SOURCE FILES HERE)
add_executable(runUnitTests basic_tests.cpp ${test_SRCS})
target_link_libraries(runUnitTests gtest gmock)
target_link_libraries(runUnitTests tshlib)
#add_test(runUnitTests runUnitTests) #included in all tutorials but I don't know what it actually does.
