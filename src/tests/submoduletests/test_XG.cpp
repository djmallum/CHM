#include "I_freeze_thaw_depths.hpp"
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
        set_default_S();
        set_default_P();
    };

    double surface_temperature = 0.0;
    XG_algorithm::param set_default_P();
    XG_algorithm::state set_default_S();
    XG_algorithm::state S = set_default_S();
    XG_algorithm::param P = set_default_P();
    void set_default_S();

    void distribute_moisture();
};


TEST_F(XGTest, NoRunDefaultOutput) {

    XG_algorithm xg(surface_temperature,S,P);
    
    double thaw = xg.get_thaw_front_depth();

    double freeze = xg.get_freeze_front_depth();

    double first = xg.get_first_front_depth();

    ASSERT_EQ(thaw,0.0);
    ASSERT_EQ(freeze,0.0);
    ASSERT_EQ(first,0.0);
};


XG_algorithm::state XGTest::set_default_S(void)
{
    XG_algorithm::state S();

    S.nfront = 0;
    S.Fz_low = 1;
    S.Th_low = 1;
    S.t_trend = 0.0;
    S.Zd_front.assign(P.n,0.0);
    S.XG_moist_d = 0.0;
    S.XG_rechr_d = 0.0;
    S.tc_composite.assign(P.n,0.0);
    S.tc_composite2.assign(P.n,0.0);
    S.XG_max.assign(P.n);
        
    std::transform(P.por.begin(),P.por.end(),
            P.depths.begin(),
            S.XG_max.begin(), 
            [](double x,double y) {return x* y;}
            );
    
    S.theta = P.theta_default;



};

XG_algorithm::param XGTest::set_default_P(void)
{
    const int n = 10;
    std::vector<double> depths(n,0.5);
    std::vector<double> por(n,0.5);
    std::vector<double> theta_default(n,0.5);
    std::vector<double> soil_km(n,2.5);
    std::vector<double> soil_km_ki(n,1.98);
    soil_km_ki[0] = 1.55;
    std::vector<double> soil_km_kw(n,1.67);
    soil_km_kw[0] = 0.8;
    XG_algorithm::param P(100,depths,por,n,theta_default,0.001,soil_km,soil_km_ki,soil_km_kw,0.35,1,1,1,250.0,750.0);

    return P;
};

void XGTest::distribute_moisture()
{
    
};

