#include <concepts>
#include <gtest/gtest.h>
#include "Ayers_infiltration.hpp"

class AyersTest : public ::testing::Test
{
protected:

};

TEST_F(AyersTest,InfiltrationMap)
{
    // Confirms that the class doesn't actually hold the map, just a reference
    EXPECT_LE(sizeof(Ayers::infiltration_map), sizeof(void*) *3 );

    Ayers::infiltration_map map;

    auto output = map.max_infil_lookup("coarse_over_coarse","bare_soil");
    EXPECT_DOUBLE_EQ(output,7.6);

    output = map.max_infil_lookup("medium_over_fine","good_pasture");
    EXPECT_DOUBLE_EQ(output,5.1);

    output = map.max_infil_lookup("medium_over_medium","row_crop");
    EXPECT_DOUBLE_EQ(output,5.1);

};

TEST_F(AyersTest,RainInfiltrarion)
{
    Ayers::rain_infiltration rain_inf;

    // small value
    auto rainfall = 0.0001;
    auto result = rain_inf.compute("coarse_over_coarse","forested",rainfall);
    EXPECT_DOUBLE_EQ(result.runoff,0.0);
    EXPECT_DOUBLE_EQ(result.inf,rainfall);

    // big value
    rainfall = 1000.0;
    result = rain_inf.compute("medium_over_medium","bare_soil",rainfall); 
    Ayers::infiltration_map map;

    auto output = map.max_infil_lookup("medium_over_medium","bare_soil");
    EXPECT_DOUBLE_EQ(result.runoff,rainfall - output);
    EXPECT_DOUBLE_EQ(result.inf,output);

    rainfall = output + 1e-13;
    result = rain_inf.compute("medium_over_fine","poor_pasture",rainfall);
    output = map.max_infil_lookup("medium_over_fine","poor_pasture");

    EXPECT_DOUBLE_EQ(result.inf,output);
    EXPECT_EQ(result.runoff,0.0);
    
};

TEST_F(AyersTest,AyersModel)
{

    class data
    {
    public:
        double rainfall() { return 5.5; };
        double snowmelt() { return 1.2; };
        const std::string& texture() { 
            static constexpr std::string s = "coarse_over_coarse";
            return s;
        };
        const std::string& ground_cover() {
            static constexpr std::string s = "forested";
            return s;
        };

        void runoff(const double in) { _runoff = in;};
        void snow_infiltrated(const double in) {_snow_infiltrated = in;};
        void infiltrated(const double in) { _infiltrated = in;};

        double _runoff = 0.0;
        double _snow_infiltrated = 0.0;
        double _infiltrated = 0.0;
    } d;

    Ayers::Model<data> ayers;

    ayers.execute(d);

    EXPECT_DOUBLE_EQ(d._runoff + d._infiltrated, d.rainfall() + d.snowmelt());

};
