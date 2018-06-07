include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

file(GLOB_RECURSE UnitTests_sources *.cpp)

add_executable(Unit_Tests_run ${UnitTests_sources} testData.hpp)
target_link_libraries(Unit_Tests_run gtest gtest_main)
target_link_libraries(Unit_Tests_run tshlib)

#add_test(NAME Unit_Tests_run
#        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests/unit_tests/
#        COMMAND ${PROJECT_BINARY_DIR}/tests/unit_tests/Unit_Tests_run})

add_test(NAME UNIT_ComputeDepthNode
        COMMAND Unit_Tests_run "--gtest_output=xml:${PROJECT_BINARY_DIR}/unit_tests.xml")