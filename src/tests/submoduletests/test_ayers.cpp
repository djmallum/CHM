#include <gtest/gtest.h>
#include "Ayers_old.hpp"  // new version
#include "submodules/Ayers.hpp" // old version
#include "Soil.h" // your soil data class

class AyersComparisonTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common test data
        soil_data_obj = std::make_shared<Soil::soils_na>();
        
        // Test cases covering different scenarios
        test_cases = {
            {"coarse_over_coarse", "bare_soil", 10.0, 0.0},
            {"medium_over_medium", "row_crop", 5.0, 0.0},
            {"fine_over_fine", "good_pasture", 3.0, 0.0},
            {"soil_over_bedrock", "forested", 1.0, 0.0},
            {"coarse_over_coarse", "forested", 100.0, 0.0}, // high rainfall
            {"medium_over_medium", "bare_soil", 0.0, 5.0}, // snowmelt only
            {"fine_over_fine", "poor_pasture", 2.0, 3.0}, // mixed
            {"coarse_over_coarse", "row_crop", 0.0, 0.0} // zero input
        };
    }

    std::shared_ptr<Soil::soils_na> soil_data_obj;
    std::vector<std::tuple<std::string, std::string, double, double>> test_cases;
};

// Test that texture maps are identical
TEST_F(AyersComparisonTest, TextureMapsAreIdentical) {
    AyersTextureLookup new_lookup;
    
    // Test all combinations in the original map
    std::vector<std::string> textures = {
        "coarse_over_coarse", "medium_over_medium", "medium_over_fine", 
        "fine_over_fine", "soil_over_bedrock"
    };
    
    std::vector<std::string> covers = {
        "bare_soil", "row_crop", "poor_pasture", "small_grains", 
        "good_pasture", "forested"
    };
    
    for (const auto& texture : textures) {
        for (const auto& cover : covers) {
            double new_value = new_lookup.get_max_infiltration(texture, cover);
            double old_value = (soil_data_obj.get()->*&Soil::soils_na::ayers_texture)(texture, cover);
            
            EXPECT_DOUBLE_EQ(new_value, old_value) 
                << "Mismatch for texture: " << texture << ", cover: " << cover;
        }
    }
}

// Test individual cases with rainfall only
TEST_F(AyersComparisonTest, RainfallOnlyCases) {
    for (const auto& [texture, cover, rainfall, snowmelt] : test_cases) {
        // Old version
        Ayers_old<Soil::soils_na, &Soil::soils_na::ayers_texture> old_ayers(
            rainfall, snowmelt, texture, cover, *soil_data_obj);
        old_ayers.run();
        
        // New version - create mock data object
        struct TestData {
            double _rainfall, _snowmelt;
            std::string _texture, _ground_cover;
            double _runoff, _inf, _snow_inf;
            
            double rainfall() const { return _rainfall; }
            double snowmelt() const { return _snowmelt; }
            std::string texture() const { return _texture; }
            std::string ground_cover() const { return _ground_cover; }
            void runoff(double value) { _runoff = value; }
            void inf(double value) { _inf = value; }
            void snow_inf(double value) { _snow_inf = value; }
        };
        
        TestData test_data{rainfall, snowmelt, texture, cover, 0.0, 0.0, 0.0};
        Ayers<TestData> new_ayers;
        new_ayers.execute(test_data);
        
        // Compare results
        EXPECT_NEAR(test_data._runoff, old_ayers.get_runoff(), 1e-12)
            << "Runoff mismatch for texture: " << texture << ", cover: " << cover;
        EXPECT_NEAR(test_data._inf, old_ayers.get_inf(), 1e-12)
            << "Infiltration mismatch for texture: " << texture << ", cover: " << cover;
        EXPECT_NEAR(test_data._snow_inf, old_ayers.get_snow_inf(), 1e-12)
            << "Snow infiltration mismatch for texture: " << texture << ", cover: " << cover;
    }
}

// Test edge cases
TEST_F(AyersComparisonTest, EdgeCases) {
    // Test with very small values that should round to zero
    struct TestData {
        double _rainfall = 1e-13; // Very small rainfall
        double _snowmelt = 0.0;
        std::string _texture = "coarse_over_coarse";
        std::string _ground_cover = "bare_soil";
        double _runoff, _inf, _snow_inf;
        
        double rainfall() const { return _rainfall; }
        double snowmelt() const { return _snowmelt; }
        std::string texture() const { return _texture; }
        std::string ground_cover() const { return _ground_cover; }
        void runoff(double value) { _runoff = value; }
        void inf(double value) { _inf = value; }
        void snow_inf(double value) { _snow_inf = value; }
    };
    
    TestData test_data;
    
    // Old version
    Ayers_old<Soil::soils_na, &Soil::soils_na::ayers_texture> old_ayers(
        test_data._rainfall, test_data._snowmelt, test_data._texture, 
        test_data._ground_cover, *soil_data_obj);
    old_ayers.run();
    
    // New version
    Ayers<TestData> new_ayers;
    new_ayers.execute(test_data);
    
    EXPECT_DOUBLE_EQ(test_data._runoff, old_ayers.get_runoff());
    EXPECT_DOUBLE_EQ(test_data._inf, old_ayers.get_inf());
    EXPECT_DOUBLE_EQ(0.0, old_ayers.get_runoff()); // Should be exactly zero
}

// Test multiple executions with same instance (new version advantage)
TEST_F(AyersComparisonTest, MultipleExecutions) {
    struct TestData {
        double _rainfall, _snowmelt;
        std::string _texture, _ground_cover;
        double _runoff, _inf, _snow_inf;
        
        double rainfall() const { return _rainfall; }
        double snowmelt() const { return _snowmelt; }
        std::string texture() const { return _texture; }
        std::string ground_cover() const { return _ground_cover; }
        void runoff(double value) { _runoff = value; }
        void inf(double value) { _inf = value; }
        void snow_inf(double value) { _snow_inf = value; }
    };
    
    Ayers<TestData> new_ayers; // Single instance
    
    // Execute multiple times with different data
    std::vector<TestData> multiple_cases = {
        {10.0, 0.0, "coarse_over_coarse", "bare_soil", 0.0, 0.0, 0.0},
        {5.0, 2.0, "medium_over_medium", "row_crop", 0.0, 0.0, 0.0},
        {0.0, 3.0, "fine_over_fine", "good_pasture", 0.0, 0.0, 0.0}
    };
    
    for (auto& test_data : multiple_cases) {
        // Old version needs new instance each time
        Ayers_old<Soil::soils_na, &Soil::soils_na::ayers_texture> old_ayers(
            test_data._rainfall, test_data._snowmelt, test_data._texture,
            test_data._ground_cover, *soil_data_obj);
        old_ayers.run();
        
        // New version reuses same instance
        new_ayers.execute(test_data);
        
        EXPECT_NEAR(test_data._runoff, old_ayers.get_runoff(), 1e-12);
        EXPECT_NEAR(test_data._inf, old_ayers.get_inf(), 1e-12);
        EXPECT_NEAR(test_data._snow_inf, old_ayers.get_snow_inf(), 1e-12);
    }
}

// Test concept enforcement
TEST_F(AyersComparisonTest, ConceptEnforcement) {
    // This should compile - valid data object
    struct ValidData {
        double _rainfall = 5.0;
        std::string _texture = "coarse_over_coarse";
        // ... other required members
        
        double rainfall() const { return _rainfall; }
        double snowmelt() const { return 0.0; }
        std::string texture() const { return _texture; }
        std::string ground_cover() const { return "bare_soil"; }
        void runoff(double) {}
        void inf(double) {}
        void snow_inf(double) {}
    };
    
    // This should work
    Ayers<ValidData> ayers;
    ValidData data;
    EXPECT_NO_THROW(ayers.execute(data));
    
    // Note: Invalid data objects would fail to compile due to concept
    // This is tested by attempting to compile with missing methods
}

// Performance test - compare creation and execution time
TEST_F(AyersComparisonTest, PerformanceComparison) {
    const int iterations = 1000;
    
    auto start_old = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        Ayers_old<Soil::soils_na, &Soil::soils_na::ayers_texture> old_ayers(
            10.0, 0.0, "coarse_over_coarse", "bare_soil", *soil_data_obj);
        old_ayers.run();
        // Destructor called each iteration
    }
    auto end_old = std::chrono::high_resolution_clock::now();
    
    struct TestData {
        double _rainfall = 10.0;
        double _snowmelt = 0.0;
        std::string _texture = "coarse_over_coarse";
        std::string _ground_cover = "bare_soil";
        double _runoff, _inf, _snow_inf;
        
        double rainfall() const { return _rainfall; }
        double snowmelt() const { return _snowmelt; }
        std::string texture() const { return _texture; }
        std::string ground_cover() const { return _ground_cover; }
        void runoff(double value) { _runoff = value; }
        void inf(double value) { _inf = value; }
        void snow_inf(double value) { _snow_inf = value; }
    };
    
    auto start_new = std::chrono::high_resolution_clock::now();
    Ayers<TestData> new_ayers; // Single instance
    for (int i = 0; i < iterations; ++i) {
        TestData data;
        new_ayers.execute(data);
    }
    auto end_new = std::chrono::high_resolution_clock::now();
    
    auto duration_old = std::chrono::duration_cast<std::chrono::microseconds>(end_old - start_old);
    auto duration_new = std::chrono::duration_cast<std::chrono::microseconds>(end_new - start_new);
    
    std::cout << "Old version: " << duration_old.count() << " microseconds" << std::endl;
    std::cout << "New version: " << duration_new.count() << " microseconds" << std::endl;
    
    // New version should be faster due to single instance creation
    EXPECT_LT(duration_new.count(), duration_old.count());
}
/*
 * CrackTest: Wrapper class for tests
 * CrackTest is effectively a mock of Infil_All module but done indirectly. Due to the complexity of the module classes, it was easier to write this.  
 * The member variables with the _ prefix are inputs that are supplied to the constructor of Crack.
 * Default values are given and used for most tests.
 * Other member variables are parameters that are also supplied to Crack unless it has the const specifier, then it is just useful for these tests.
 * Member functions are just tools to enable the tests.
 * Initialization of CrackTest assumes that the frozen period has just begun. 
 * 
 */
//class AyersTest : public testing::Test
//{
//protected:
//
//    AyersTest()
//    {
//    };
//
//	typedef const Soil::soils_na S;
//
//	typedef Ayers<S,&S::ayers_texture> MyAyers;
//
//    double _snowmelt = 1.0;
//    double _rainfall = 0.0;
//	std::string _ground_cover = "bare_soil";
//	std::string _texture = "coarse_over_coarse";		
//
//	S& soils = Soil::get_soil_obj<S>();
//
//	MyAyers DoAyers(double& snowmelt, double& rainfall)
//	{
//		MyAyers ayers(rainfall,snowmelt,_texture,_ground_cover, soils);
//
//        ayers.run();
//		return ayers;
//	};
//
//	void DoAssert(MyAyers& ayers,const double inf,const double runoff,const double snowinf)
//	{
//		ASSERT_EQ(ayers.get_inf(),inf);
//		ASSERT_EQ(ayers.get_runoff(),runoff);
//		ASSERT_EQ(ayers.get_snow_inf(),snowinf);
//	};
//};
//
//TEST_F(AyersTest, ZeroInputs)
//{
//	_rainfall = 0.0;
//	_snowmelt = 0.0;
//
//	MyAyers ayers = DoAyers(_snowmelt, _rainfall);
//
//	DoAssert(ayers,0.0,0.0,0.0);
//};
//
//TEST_F(AyersTest, NonZeroRainfallSmall)
//{
//	_rainfall = 1e-6;
//	_snowmelt = 0.0;
//
//	MyAyers ayers = DoAyers(_snowmelt,_rainfall);
//
//	DoAssert(ayers,_rainfall,0.0,0.0);
//
//};
//
//TEST_F(AyersTest, NonZeroBigRainfall)
//{
//	_rainfall = 1e4;
//	_snowmelt = 0.0;
//
//	MyAyers ayers = DoAyers(_snowmelt,_rainfall);
//
//	DoAssert(ayers,7.6,_rainfall - 7.6,0.0);
//
//};
//
//TEST_F(AyersTest, NonZeroSnowMelt)
//{
//	_rainfall = 0.0;
//	_snowmelt = 1.3e2;
//
//	MyAyers ayers1 = DoAyers(_snowmelt,_rainfall);
//
//	DoAssert(ayers1,_snowmelt,0.0,_snowmelt);
//
//	_snowmelt = 3.14159;
//
//	MyAyers ayers2 = DoAyers(_snowmelt,_rainfall);
//
//	DoAssert(ayers2,_snowmelt,0.0,_snowmelt);
//
//};
//
//TEST_F(AyersTest, MeltAndRain)
//{
//	_rainfall = 99;
//	_snowmelt = 13;
//	
//	_ground_cover = "small_grains";
//	_texture = "medium_over_medium";
//
//	MyAyers ayers = DoAyers(_snowmelt,_rainfall);
//	double maxinfil = 10.2;
//	DoAssert(ayers,maxinfil+_snowmelt,_rainfall - maxinfil,_snowmelt);
//
//};
//
//
