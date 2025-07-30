#include <gtest/gtest.h>
#include "data_base.hpp"
#include <limits>
#include <memory>

// Mock CacheType for testing
struct MockCache : public cache_base {
    double value = std::numeric_limits<double>::quiet_NaN();
};

class DataBaseTest : public ::testing::Test {
protected:
    void SetUp() override {
        mock_global = boost::make_shared<global>();
    }

    mesh_elem mock_face;
    boost::shared_ptr<global> mock_global;
    pt::ptree mock_cfg;

    static inline auto A = [](){ return 42.0;};
    static inline auto B = [](){ return 123.0;};
};

class data : public data_base<MockCache>
{
public:
    data(mesh_elem face, boost::shared_ptr<global> param, pt::ptree& cfg)
        : data_base<MockCache>(face,param,cfg,true) {};
    ~data() {};

    template<typename V,typename L>
    void update(V&& v, L&& l)
    {
        update_field(v, l);
    };
    
    template<typename T>
    void output(T& output, const T& t)
    {
       set_output(output,t);
    }; 

};

TEST_F(DataBaseTest, CacheInitializesOnFirstUpdate) {
    data db(mock_face, mock_global, mock_cfg);
    double test_value = std::numeric_limits<double>::quiet_NaN();

    db.update(test_value, [](){ return 42.0; });

    EXPECT_FALSE(std::isnan(test_value));
    EXPECT_EQ(test_value, 42.0);
}

TEST_F(DataBaseTest, CacheRespectsTimestepChanges) {
    data db(mock_face, mock_global, mock_cfg);
    double test_value = std::numeric_limits<double>::quiet_NaN();
    
    // First call at timestep 0
    mock_global->timestep_counter = 0;
    db.update(test_value, A);

    // Second call at same timestep - should use cached value
    db.update(test_value,B);
    EXPECT_EQ(test_value, A());  // No update to B

    // Force cache reset by changing timestep
    mock_global->timestep_counter = 1;
    db.update(test_value, B);
    EXPECT_EQ(test_value, B()); // Updates
}

TEST_F(DataBaseTest, ManualCacheResetWorks) {
    data db(mock_face, mock_global, mock_cfg);
    double test_value = std::numeric_limits<double>::quiet_NaN();

    db.update(test_value, A);
    db.reset_cache();

    db.update(test_value, B);
    EXPECT_EQ(test_value, B());
}

TEST_F(DataBaseTest, SetOutputUpdatesValue) {
    data db(mock_face, mock_global, mock_cfg);
    double output = 0.0;
    
    db.output(output, 3.14);
    EXPECT_EQ(output, 3.14);
}

TEST_F(DataBaseTest, OnlyUpdatesNaNValues) {
    data db(mock_face, mock_global, mock_cfg);
    double test_value = 10.0;  // Not NaN

    db.update(test_value, A);
    EXPECT_EQ(test_value, 10.0);  // Should remain unchanged
}

TEST_F(DataBaseTest, CacheInitializesViaBothMethods) {
    data db(mock_face, mock_global, mock_cfg);
    
    // Test via update_field
    double test_value = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(db.get_cache().has_value());  // Cache should start empty
    
    mock_global->timestep_counter = 99;

    db.update(test_value, A);
    EXPECT_TRUE(db.get_cache().has_value());  // Cache should now exist
    EXPECT_EQ(db.get_cache()->last_timestep, mock_global->timestep_counter);

    // Reset and test via set_output
    db.reset_cache();
    double output_value = 0.0;
    
    db.output(output_value, 3.14);
    EXPECT_TRUE(db.get_cache().has_value());  // Cache should be recreated
    EXPECT_EQ(db.get_cache()->last_timestep, -1);            
}

// Mock CacheType for testing
//struct MockCache : cache_base {
//    double test_value = std::numeric_limits<double>::quiet_NaN();
//};
//
//class data : public data_base<MockCache>
//{
//public:
//    double& test_value()
//    {
//        return cache_->test_value;
//    };
//
//    void update_test1()
//    {
//        update_field(test_value(), []{ return 42.0; });
//    };
//
//    void update_test2()
//    {
//        update_field(test_value(), []{ return 200.0; });
//    };
//
//    size_t get_timestep()
//    {
//        return cache_->last_timestep;
//    };
//
//};
//
//// Fixture
//class DataBaseTest : public ::testing::Test {
//protected:
//    void SetUp() override {
//        mock_global = boost::make_shared<global>();
//        mock_global->timestep_counter = 0;
//        d.set_global(mock_global);
//    }
//
//    data d;
//    boost::shared_ptr<global> mock_global;
//};
//
//// Tests
//TEST_F(DataBaseTest, InitCacheOnFirstUse) {
//    d.reset_cache();
//    d.update_test1();
//    EXPECT_EQ(d.get_timestep(), 0);
//    EXPECT_DOUBLE_EQ(d.test_value(), 42.0);
//}
//
//TEST_F(DataBaseTest, CachePersistsSameTimestep) {
//    d.update_test1();
//    EXPECT_DOUBLE_EQ(d.test_value(),42.0);
//    d.update_test2();
//    EXPECT_DOUBLE_EQ(d.test_value(), 42.0);
//}
//
//TEST_F(DataBaseTest, CacheResetsOnNewTimestep) {
//    d.update_test1();
//
//    EXPECT_DOUBLE_EQ(d.test_value(), 42.0);
//    d.reset_cache();
//    d.update_test2();
//
//    EXPECT_DOUBLE_EQ(d.test_value(), 200.0);
//}
//
