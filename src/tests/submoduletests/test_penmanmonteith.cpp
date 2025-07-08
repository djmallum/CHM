#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "CSVreader.hpp"
#include "PenmanMonteith.hpp"
#define diff 1e-4
class PenmanMonteithTest : public ::testing::Test
{
protected:
    
    double leaf_area_index = 1.0;
    //double LAImin;
    //double seasonal_growth;
    double leaf_area_index_max =  0.5;
    double Veg_height = 1.5;
    double wind_measurement_height = 10.0; // This one might be uniform...
    double stomatal_resistance_min = 25.0; // Also might be domain wide
    double soil_depth = 1.5;
    double Frac_to_ground = 0.1;
    const double heat_capacity_air = 1005;
    const double kappa = 0.4; // also might be domain wide
    const double air_entry_tension = 0.786;
    const double pore_size_dist = 5.3; 
    const double wilt_point = 0.13;
    const double porosity = 0.485;
    double wind_speed, short_wave_in, all_wave_net, t, soil_storage,
              vapour_pressure,saturated_vapour_pressure,P_atm;

    PM_output out;
    const int i = 420;
    PM_vars getPMvar(const int ii) {
        wind_speed = reader.getValue<double>("hru_u",ii);
        short_wave_in = reader.getValue<double>("Qsi",ii);
        all_wave_net = reader.getValue<double>("Rn",ii);
        t = reader.getValue<double>("hru_t",ii);
        soil_storage = reader.getValue<double>("soil_moist",ii-1);
        vapour_pressure = reader.getValue<double>("hru_ea",ii);
        saturated_vapour_pressure = estar(t);
        P_atm = getP(); 
        
        PM_vars var(wind_speed,short_wave_in,all_wave_net,t,soil_storage,vapour_pressure,saturated_vapour_pressure,P_atm);
       return var; 
    }

    double estar(double t)
    {
       if (t > 0.0)
            return 0.611 * exp(17.27 * t / (t + 237.3));
       else
            return 0.611 * exp(21.88 * t / (t + 265.5));
    }; 
    double const elevation = 1845;
    double getP()
    {
        return 101.3 * pow( ( 293.0 - 0.0065 * elevation)/293.0,5.26);
    }; 

};

TEST_F(PenmanMonteithTest, CompareWithCRHMTest)
{
    PenmanMonteith PM(leaf_area_index,leaf_area_index_max,Veg_height,
            wind_measurement_height,stomatal_resistance_min,soil_depth,
            Frac_to_ground,heat_capacity_air,kappa,air_entry_tension,
            pore_size_dist,wilt_point,porosity);
    PM_vars var = getPMvar(i);
    

    PM.CalcEvapT(var,out);

    double rc = reader.getValue<double>("rc",i);
    double ET = reader.getValue<double>("hru_evap",i);

    EXPECT_NEAR(out.ET,ET,diff);
    EXPECT_NEAR(out.stomatal_resistance,rc,diff);
};
