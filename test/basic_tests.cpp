//
// Created by Lorenzo Gatti on 08.10.17.
//

#include "gtest/gtest.h"

namespace {
// The fixture for testing class TSHLIB_LowLevel.
    class TSHLIB_LowLevel : public ::testing::Test {
    protected:
        // You can remove any or all of the following functions if its body
        // is empty.

        TSHLIB_LowLevel() {
            // You can do set-up work for each test here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:

        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        ~TSHLIB_LowLevel() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // Objects declared here can be used by all tests in the test case for Foo.
    };
    /*
    // Tests that the Foo::Bar() method does Abc.
    TEST_F(FooTest, MethodBarDoesAbc) {
        const string input_filepath = "this/package/testdata/myinputfile.dat";
        const string output_filepath = "this/package/testdata/myoutputfile.dat";
        Foo f;
        EXPECT_EQ(0, f.Bar(input_filepath, output_filepath));
    }
    */

    TEST_F(TSHLIB_LowLevel, test_eq) {

        EXPECT_EQ(1, 1);

    }

    TEST_F(TSHLIB_LowLevel, test_neq) {
        EXPECT_NE(1, 0);

    }

}  // namespace


namespace {
// The fixture for testing class TSHLIB_LowLevel.
    class TSHLIB_MiddleLevel : public ::testing::Test {
    protected:
        // You can remove any or all of the following functions if its body
        // is empty.

        TSHLIB_MiddleLevel() {
            // You can do set-up work for each test here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:

        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        ~TSHLIB_MiddleLevel() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // Objects declared here can be used by all tests in the test case for Foo.
    };
    /*
    // Tests that the Foo::Bar() method does Abc.
    TEST_F(FooTest, MethodBarDoesAbc) {
        const string input_filepath = "this/package/testdata/myinputfile.dat";
        const string output_filepath = "this/package/testdata/myoutputfile.dat";
        Foo f;
        EXPECT_EQ(0, f.Bar(input_filepath, output_filepath));
    }
    */

    TEST_F(TSHLIB_MiddleLevel, test_eq) {

        EXPECT_EQ(1, 1);

    }

    TEST_F(TSHLIB_MiddleLevel, test_neq) {
        EXPECT_NE(1, 1);

    }

}  // namespace


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}