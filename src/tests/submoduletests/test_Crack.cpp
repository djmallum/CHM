#include <gtest/gtest.h>
#include <span>
#include <stdexcept>
#include "Crack.hpp"

class CrackTest_TestOffDay;
//class CrackTest;

class Data
{
public:
    
    FRIEND_TEST(CrackTest,TestOffDay);
    FRIEND_TEST(CrackTest,HasLens);
    FRIEND_TEST(CrackTest,Accumulate);
    FRIEND_TEST(CrackTest,SoilSaturatedAtFreezeResultsInZeroInfiltration);
    FRIEND_TEST(CrackTest,PriorMajorMelt);
    FRIEND_TEST(CrackTest,Restricted);
    FRIEND_TEST(CrackTest,Unlimited);
    FRIEND_TEST(CrackTest,EquationCalculation);
    FRIEND_TEST(CrackTest,Limited);
    friend class CrackTest;
private:    
    Crack::State _state;

    bool _new_day;
    double _snowmelt;
    double _rainfall;
    double _swe;
    double _t;

    double _inf;
    double _runoff;
    double _snow_inf;
    double _melt_runoff;
    double _rain_on_snow;

public:
    Crack::State& get_state() { return _state;};
    const bool is_newday() { return _new_day;};
    double snowmelt() { return _snowmelt;};
    double rainfall() { return _rainfall;};
    double swe() { return _swe;};
    double air_temperature() { return _t;};

    void infiltrated(const double out) { _inf = out;};
    void runoff(const double out) { _runoff = out;};
    void snow_infiltrated(const double out) { _snow_inf = out;};
    void melt_runoff(const double out) { _melt_runoff = out;};
    void rain_on_snow(const double out) { _rain_on_snow = out;};
};

using namespace Crack;

class CrackTest : public ::testing::Test
{
protected:
    Model<Data> crack;
    Data d;
    State* s;

    static constexpr size_t HOURS_PER_DAY = static_cast<size_t>(24);
    static constexpr size_t DAYS = static_cast<size_t>(3);
    void SetUp() override
    {
        s = &d.get_state(); 
        Params::major_melt_threshold = 5.0;
        Params::infDays = 6;
        Params::lenstemp = -10.0;
        Params::allow_early_inf = true;
        
    }

    void assign_current(std::span<double> v)
    {
        s->current_inf = v[0];
        s->current_runoff = v[1];
        s->current_snow_inf = v[2];
        s->current_melt_runoff = v[3];    
    };
    
    bool compare_current(std::span<double,4> v)
    {

        return false; 
    };
    
    void wind_up_limited(const double saturation)
    {
        d._swe = 100.0;
        d._snowmelt = 0.0;
        d._t = -5.0;
        s->soil_saturation_at_freeze = saturation;
        constexpr size_t NUM_DAYS = 15;

        for (size_t day = 0; day < NUM_DAYS; ++day)
        {
            for (size_t hour = 0; hour < HOURS_PER_DAY; ++hour)
            {
                d._new_day = hour == static_cast<size_t>(0);
                crack.execute(d);
            };
        };
    }

    void bind()
    {
        s->daily_melt_total.bind_target(d._snowmelt);
        s->daily_rain_total.bind_target(d._rainfall);
    };
}; 

TEST_F(CrackTest,TestOffDay)
{
    d._new_day = false;
    EXPECT_THROW(
            {
            crack.execute(d);
            },std::logic_error);
    bind();
    
    std::vector<double> v{1.0,2.0,8.5,1.3}; 
    assign_current(v);
    crack.execute(d); 

    EXPECT_EQ(d._inf,v[0]);
    EXPECT_EQ(d._runoff,v[1]);
    EXPECT_EQ(d._snow_inf,v[2]);
    EXPECT_EQ(d._melt_runoff,v[3]);

};

TEST_F(CrackTest,Accumulate)
{
    d._snowmelt = 1.0;
    d._rainfall = 3.0;
    bind();
    d._new_day = false;

    crack.execute(d);
    crack.execute(d);
    crack.execute(d);
    
    EXPECT_DOUBLE_EQ(s->daily_melt_total.get_yesterday(),0.0);
    EXPECT_DOUBLE_EQ(s->daily_rain_total.get_yesterday(),0.0);
    d._new_day = true;

    crack.execute(d);
    
    constexpr size_t days = static_cast<size_t>(3);
    EXPECT_DOUBLE_EQ(s->daily_melt_total.get_yesterday(),d._snowmelt * days);
    EXPECT_DOUBLE_EQ(s->daily_rain_total.get_yesterday(),d._rainfall * days);
};

TEST_F(CrackTest,HasLens)
{
    d._new_day = false;
    s->soil_saturation_at_freeze = 50.0;
    s->major_melt_count = 3;
    bind();
    d._t = -30.0;
    // daily_max_temp is set to _t on the first new day, so
    // Should rewrite this test so that _t sets it properly
    s->daily_max_temp = -30.0;
    auto melt = 5.5;
    auto rain = 2.0;
    d._snowmelt = melt;
    d._rainfall = rain;
    crack.execute(d); //Once to accumulate
    d._new_day = true;
    crack.execute(d); //Once to activate 

    EXPECT_DOUBLE_EQ(d._runoff,melt + rain);
    EXPECT_TRUE(s->major_melt_count > Params::infDays);
};



TEST_F(CrackTest, SoilSaturatedAtFreezeResultsInZeroInfiltration)
{
    // Setup: Soil is saturated at freeze (100% saturation)
    s->soil_saturation_at_freeze = 100.0;

    // Initialize other state variables
    s->current_inf = 0.0;
    s->current_runoff = 0.0;
    s->current_snow_inf = 0.0;
    s->current_melt_runoff = 0.0;
    s->major_melt_count = 0;
    s->index = 0.0;
    s->max_major_per_melt = 0.0;
    s->init_SWE = 0.0;
    s->daily_max_temp = -5.0;  // Below freezing

    // Bind accumulators to data members
    bind();

    // Test parameters
    const double snowmelt_rate = 3.0;  // mm/hour
    const int calls_before_new_day = 5;  // N calls before new day
    d._new_day = false;
    d._snowmelt = snowmelt_rate;
    d._rainfall = 0.0;
    d._swe = 10.0;
    d._t = -2.0;  // Below freezing

    double total_snowmelt_before_new_day = 0.0;

    // Call execute N times before a new day
    for (int i = 0; i < calls_before_new_day; i++) {
        // Execute the crack model
        crack.execute(d);

        // Accumulate snowmelt for verification
        total_snowmelt_before_new_day += snowmelt_rate;

        // Verify infiltration is zero on each call
        EXPECT_EQ(d._inf, 0.0) << "Infiltration should be 0.0 on call " << i;
        EXPECT_EQ(d._snow_inf, 0.0) << "Snow infiltration should be 0.0 on call " << i;
        
        EXPECT_EQ(d._runoff, 0.0)
            << "Runoff should equal 0.0 on call before first new day " << i;
        EXPECT_EQ(d._melt_runoff, 0.0)
            << "Melt runoff should equal 0.0 on call before first new day " << i;

        // Verify runoff equals accumulated snowmelt
        //EXPECT_EQ(d._runoff, total_snowmelt_before_new_day)
        //    << "Runoff should equal total snowmelt on call " << i;
        //EXPECT_EQ(d._melt_runoff, total_snowmelt_before_new_day)
        //    << "Melt runoff should equal total snowmelt on call " << i;
    }

    // Now simulate a new day
    d._new_day = true;

    // Execute once more with new day flag
    crack.execute(d);

    // After new day processing, verify final state
    EXPECT_EQ(d._inf, 0.0) << "Infiltration should still be 0.0 after new day";
    EXPECT_EQ(d._snow_inf, 0.0) << "Snow infiltration should still be 0.0 after new day";

    // Runoff should include the additional snowmelt from the new day call
    EXPECT_EQ(d._runoff, total_snowmelt_before_new_day)
        << "Runoff should equal total snowmelt across all calls";
    EXPECT_EQ(d._melt_runoff, total_snowmelt_before_new_day)
        << "Melt runoff should equal total snowmelt across all calls";

    // Verify that soil_saturation_at_freeze didn't change
    EXPECT_EQ(s->soil_saturation_at_freeze, 100.0)
        << "soil_saturation_at_freeze should remain 100.0";

    // Verify state was properly reset on new day
    EXPECT_EQ(s->current_inf, 0.0) << "current_inf should be zero because no inf possible";
    EXPECT_EQ(s->current_runoff, total_snowmelt_before_new_day) << "current_runoff should equal to the total_snowmelt_before_new_day";
    EXPECT_EQ(s->current_snow_inf, 0.0) << "current_snow_inf should be reset to 0.0";
    EXPECT_EQ(s->current_melt_runoff, total_snowmelt_before_new_day) << "current_melt_runoff should be equal to the total_snowmelt_before_new_day";
}

TEST_F(CrackTest,PriorMajorMelt)
{
    d._snowmelt = 0.01;
    d._swe = 100.0;
    bind();

    d._new_day = true;

    ASSERT_TRUE(d._snowmelt * HOURS_PER_DAY < Params::major_melt_threshold);

    s->soil_saturation_at_freeze = 75.0;
    for (size_t ii = 0; ii < HOURS_PER_DAY; ++ii)
    {
        crack.execute(d);
        d._new_day = false;

    };

    d._new_day = true;
    crack.execute(d);

    EXPECT_EQ(s->major_melt_count,0);
    EXPECT_TRUE(s->daily_melt_total.get_yesterday() < Params::allow_early_inf);
    EXPECT_DOUBLE_EQ(s->daily_melt_total.get_yesterday(),d._snowmelt * HOURS_PER_DAY);
};

TEST_F(CrackTest,Restricted)
{
    constexpr double INITIAL_SATURATION = 100.0;

    bind();
    wind_up_limited(INITIAL_SATURATION);
    
    d._snowmelt = 0.6;
    
    ASSERT_TRUE(d._snowmelt * HOURS_PER_DAY > Params::major_melt_threshold);
    
    
    for (size_t day = 0; day < DAYS; ++day)
    {
        for (size_t hour = 0; hour < HOURS_PER_DAY; ++hour)
        {
            d._new_day = hour == static_cast<size_t>(0);
            crack.execute(d);
        };
        switch (day)
        {
            case 0:
                EXPECT_DOUBLE_EQ(d._inf,0.0);
                EXPECT_DOUBLE_EQ(d._runoff,0.0);
                EXPECT_DOUBLE_EQ(d._snow_inf,0.0);
                EXPECT_DOUBLE_EQ(d._melt_runoff,0.0);
                break;
            case 1:
            case 2:
                EXPECT_DOUBLE_EQ(d._inf,0.0);
                EXPECT_DOUBLE_EQ(d._runoff,d._snowmelt * HOURS_PER_DAY);
                EXPECT_DOUBLE_EQ(d._snow_inf,0.0);
                EXPECT_DOUBLE_EQ(d._melt_runoff,d._snowmelt * HOURS_PER_DAY);
                break;
            default:
                ASSERT_TRUE(false) << "Shouldn't be reachable...";
                break;
        }
    };
};

TEST_F(CrackTest,Unlimited)
{
    constexpr double INITIAL_SATURATION = 0.0;

    bind();
    wind_up_limited(INITIAL_SATURATION);
    
    d._snowmelt = 0.6;
    
    
    ASSERT_TRUE(d._snowmelt * HOURS_PER_DAY > Params::major_melt_threshold);
    
    for (size_t day = 0; day < DAYS; ++day)
    {
        for (size_t hour = 0; hour < HOURS_PER_DAY; ++hour)
        {
            d._new_day = hour == static_cast<size_t>(0);
            crack.execute(d);
        };
        switch (day)
        {
            case 0:
                EXPECT_DOUBLE_EQ(d._inf,0.0);
                EXPECT_DOUBLE_EQ(d._runoff,0.0);
                EXPECT_DOUBLE_EQ(d._snow_inf,0.0);
                EXPECT_DOUBLE_EQ(d._melt_runoff,0.0);
                break;
            case 1:
            case 2:
                EXPECT_DOUBLE_EQ(d._inf,d._snowmelt * HOURS_PER_DAY);
                EXPECT_DOUBLE_EQ(d._runoff,0.0);
                EXPECT_DOUBLE_EQ(d._snow_inf,d._snowmelt * HOURS_PER_DAY);
                EXPECT_DOUBLE_EQ(d._melt_runoff,0.0);
                break;
            default:
                ASSERT_TRUE(false) << "Shouldn't be reachable...";
                break;
        }
    };
};

TEST_F(CrackTest,Limited)
{
    constexpr double INITIAL_SATURATION = 50.0;

    bind();
    wind_up_limited(INITIAL_SATURATION);

    d._snowmelt = 0.25;

    
    ASSERT_TRUE(d._snowmelt * HOURS_PER_DAY > Params::major_melt_threshold);
    auto INF = [](double swe) { 
        return 5 * ( 1 - INITIAL_SATURATION/100.0) * std::pow(swe,0.584); };
    ASSERT_TRUE(d._snowmelt * HOURS_PER_DAY < INF(d._swe) / Params::infDays) << d._snowmelt * HOURS_PER_DAY;
    
    struct Last
    {
        double index;
        double init_SWE;
        double max_major_per_melt;
    } last;
    for (size_t day = 0; day < DAYS + 1; ++day)
    {
        for (size_t hour = 0; hour < HOURS_PER_DAY; ++hour)
        {
            d._new_day = hour == static_cast<size_t>(0);
            crack.execute(d);
        };
        
        switch (day)
        {
            case 0:
                EXPECT_DOUBLE_EQ(d._inf,0.0);
                EXPECT_DOUBLE_EQ(d._runoff,0.0);
                EXPECT_DOUBLE_EQ(d._snow_inf,0.0);
                EXPECT_DOUBLE_EQ(d._melt_runoff,0.0);
                break;
            case 1:
            case 2:
                EXPECT_DOUBLE_EQ(s->init_SWE,d._swe);
                EXPECT_DOUBLE_EQ(s->index,INF(s->init_SWE)/s->init_SWE);
                EXPECT_DOUBLE_EQ(s->max_major_per_melt,s->index * d._swe / Params::infDays);

                EXPECT_DOUBLE_EQ(d._inf,d._snowmelt * HOURS_PER_DAY * s->index);
                EXPECT_DOUBLE_EQ(d._runoff,d._snowmelt * HOURS_PER_DAY - d._inf);
                EXPECT_DOUBLE_EQ(d._snow_inf,d._inf);
                EXPECT_DOUBLE_EQ(d._melt_runoff,d._runoff);

                if (day == 1)
                {
                    last = Last{s->index,s->init_SWE,s->max_major_per_melt};
                }
                else
                {
                    EXPECT_DOUBLE_EQ(s->index,last.index);
                    EXPECT_DOUBLE_EQ(s->init_SWE,last.init_SWE);
                    EXPECT_DOUBLE_EQ(s->max_major_per_melt,last.max_major_per_melt);
                    d._swe = s->init_SWE * 1.01;
                };
                break;
            case 3:
                EXPECT_NE(s->index,last.index);
                EXPECT_DOUBLE_EQ(s->index,INF(last.init_SWE*1.01)/(last.init_SWE*1.01));
                EXPECT_NE(s->init_SWE,last.init_SWE);
                EXPECT_DOUBLE_EQ(s->init_SWE,last.init_SWE * 1.01);
                EXPECT_NE(s->max_major_per_melt,last.max_major_per_melt);
                EXPECT_DOUBLE_EQ(s->max_major_per_melt,s->index * d._swe / Params::infDays);
                break;
            default:
                ASSERT_TRUE(false) << "Shouldn't be reachable...";
                break;
        }
    };
};
