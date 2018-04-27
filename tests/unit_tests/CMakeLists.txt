include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

file(GLOB_RECURSE UnitTests_sources *.cpp)

add_executable(Unit_Tests_run ${UnitTests_sources})

target_link_libraries(Unit_Tests_run gtest gtest_main)
target_link_libraries(Unit_Tests_run tshlib)