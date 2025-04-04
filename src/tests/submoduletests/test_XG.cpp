#include "XG.hpp"
#include "gtest/gtest.h"

/*
 * XGTest: Wrapper class for tests
 * XGTest is effectively a mock of Infil_All module but done indirectly. Due to the complexity of the module classes, it was easier to write this.  
 * The member variables with the _ prefix are inputs that are supplied to the constructor of Crack.
 * Default values are given and used for most tests.
 * Other member variables are parameters that are also supplied to Crack unless it has the const specifier, then it is just useful for these tests.
 * Member functions are just tools to enable the tests.
 * Initialization of XGTest assumes that the frozen period has just begun. 
 * 
 */
class XGTest : public testing::Test
{
protected:

    XGTest()
    {
    };

};


TEST_F(XGTest, PutNameHere) {
    XG xg();

    xg.run();

    xg.get_thaw_front_depth();

    xg.get_freeze_front_depth();

    xg.get_first_front_depth();
};

