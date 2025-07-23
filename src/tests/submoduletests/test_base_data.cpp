#include <gtest/gtest.h>
#include "data_base.hpp"
#include <limits>
#include <memory>

// Mock CacheType for testing
struct MockCache : cache_base {
    double test_value = std::numeric_limits<double>::quiet_NaN();
};

class data : public data_base<MockCache>
{
public:
    double& test_value()
    {
        return cache_->test_value;
    };

    void update_test1()
    {
        update_field(test_value(), []{ return 42.0; });
    };

    void update_test2()
    {
        update_field(test_value(), []{ return 200.0; });
    };

    size_t get_timestep()
    {
        return cache_->last_timestep;
    };

};

// Fixture
class DataBaseTest : public ::testing::Test {
protected:
    void SetUp() override {
        mock_global = boost::make_shared<global>();
        mock_global->timestep_counter = 0;
        d.set_global(mock_global);
    }

    data d;
    boost::shared_ptr<global> mock_global;
};

// Tests
TEST_F(DataBaseTest, InitCacheOnFirstUse) {
    d.reset_cache();
    d.update_test1();
    EXPECT_EQ(d.get_timestep(), 0);
    EXPECT_DOUBLE_EQ(d.test_value(), 42.0);
}

TEST_F(DataBaseTest, CachePersistsSameTimestep) {
    d.update_test1();
    EXPECT_DOUBLE_EQ(d.test_value(),42.0);
    d.update_test2();
    EXPECT_DOUBLE_EQ(d.test_value(), 42.0);
}

TEST_F(DataBaseTest, CacheResetsOnNewTimestep) {
    d.update_test1();

    EXPECT_DOUBLE_EQ(d.test_value(), 42.0);
    d.reset_cache();
    d.update_test2();

    EXPECT_DOUBLE_EQ(d.test_value(), 200.0);
}

