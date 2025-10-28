#include <gtest/gtest.h>
#include <iomanip> 
#include "Atmosphere.h"
#include "Penman_Montieth.hpp"

namespace legacy_data
{
    constexpr size_t N_COMPARES = 6, N_BOUNDARIES = 3;
    constexpr std::array<double,N_COMPARES> R_C_Comparison{133.33333333333334,5000.0,666.66666666666674,133.33333333333334,5000.0,5000.0};
    constexpr std::array<double,N_COMPARES> E_T_Comparison{0.095635101141937504,0.0030968023432009018,0.09712545157162103,0.095635101141937504,0.0037307822786114184,0.0029488454061942023};
    constexpr std::array<double,N_BOUNDARIES> R_C_Boundaries{5000.0,5000.0,5000.0};
    constexpr std::array<double,N_BOUNDARIES> E_T_Boundaries{0.0049307977258960124,0.00070336782232421284,0.0};
};
class MockPenmanData 
{
public:
    double wind_measurement_height() const;
    double d() const;
    double Z0() const;
    double kappa() const;
    double wind_speed() const;
    double stomatal_resistance_min() const;
    bool has_vegetation() const;
    double Veg_height() const;
    double leaf_area_index_max() const;
    double short_wave_in() const;
    double saturated_vapour_pressure() const;
    double vapour_pressure() const;
    double air_entry_tension() const;
    double porosity() const;
    double volumetric_moisture_content() const;
    double pore_size_dist() const;
    double air_temperature() const;
    double delta() const;
    double Q_net() const;
    double Q_g() const;
    double air_density() const;
    double heat_capacity_air() const;
    double aero_resistance() const;
    double gamma() const;
    double stomatal_resistance() const;
    double lambda() const;
    int s_per_time_step() const;
    void stomatal_resistance(const double T) const;
    void ET(const double T) const;

    // Test parameters
    double LAI, LAImax, veg_Ht, wind_height, stomatal_res_min;
    double soil_d, F_to_g, _s_per_step, Cp, K, tension;
    double pore_sz, theta_pwp, phi;
    double Qsw, soil_storage, ea, ea_star, P_atm;
    // Added member variables to support the functions
    double _Q_net, _air_temperature, _wind_speed;
    mutable double _stomatal_resistance, _ET;
};
// Test fixture for PenmanMonteith tests
class PenmanMonteithTest : public ::testing::Test {
protected:
        MockPenmanData data;
    void SetUp() override {
        // Common test parameters
        data.LAI = 2.5;
        data.LAImax = 3.0;
        data.veg_Ht = 1.5;
        data.wind_height = 2.0;
        data.stomatal_res_min = 100.0;
        data.soil_d = 0.75;
        data.F_to_g = 0.2;
        data._s_per_step = 3600.0;
        data.Cp = 1013.0;//Atmosphere::Cp;
        data.K = 0.41;
        data.tension = 0.2;
        data.pore_sz = 0.5;
        data.theta_pwp = 0.15;
        data.phi = 0.45;
    
        data._wind_speed = 2.0;
        data.Qsw = 500.0;
        data._Q_net = 100.0;
        data._air_temperature = 20.0;
        data.soil_storage = 300;
        data.ea = 1.5;
        data.ea_star = 2.0;
        data.P_atm = 101.3;
        data._stomatal_resistance = 0.0;
        data._ET = 0.0;

        std:: cout << std::setprecision(17); 
    };
        
    Penman_monteith<MockPenmanData> refactored_pm;
};
    
// Comparison test for stomatal resistance under various conditions
TEST_F(PenmanMonteithTest, StomatalResistanceComparison) {
    struct TestCase {
        double short_wave;
        double vapour_pressure_diff;
        double soil_storage;
        double temperature;
        std::string description;
    };
    
    std::array<TestCase,legacy_data::N_COMPARES> test_cases = {
        TestCase{600.0, 0.5, 0.25, 25.0, "normal_conditions"},
        TestCase{0.0, 0.5, 0.25, 25.0, "night_time"},
        TestCase{600.0, 2.5, 0.25, 25.0, "high_vpd"},
        TestCase{600.0, 0.5, 0.1, 25.0, "dry_soil"},
        TestCase{600.0, 0.5, 0.25, 45.0, "high_temperature"},
        TestCase{600.0, 0.5, 0.25, 2.0, "low_temperature"}
    };
    
    Penman_monteith<MockPenmanData>::stomatal_resistance_jarvis r_s;
    int i = 0;
    for (const auto& test_case : test_cases) {
        // Set up test conditions using mock data
        data.Qsw = test_case.short_wave;
        data._air_temperature = test_case.temperature;
        data.soil_storage = test_case.soil_storage;
        data.ea_star = data.ea + test_case.vapour_pressure_diff;
        
        r_s.calculate(data); 

        // Compare results
        EXPECT_NEAR(legacy_data::R_C_Comparison[i], data.stomatal_resistance(), 1e-12)
            << "Stomatal resistance mismatch for: " << test_case.description;
        EXPECT_NEAR(legacy_data::R_C_Comparison[i], data._stomatal_resistance, 1e-12)
            << "Stomatal resistance mismatch for: " << test_case.description;
        
        refactored_pm.execute(data); 
        // Compare results
        EXPECT_NEAR(legacy_data::R_C_Comparison[i], data.stomatal_resistance(), 1e-12)
            << "Stomatal resistance mismatch for (full system): " << test_case.description;
        EXPECT_NEAR(legacy_data::R_C_Comparison[i], data._stomatal_resistance, 1e-12)
            << "Stomatal resistance mismatch for (full system): " << test_case.description;
        EXPECT_NEAR(legacy_data::E_T_Comparison[i], data._ET, 1e-12)
            << "ET mismatch for (full system): " << test_case.description;
        ++i;
    }
    
}

// Comparison test for boundary conditions
TEST_F(PenmanMonteithTest, BoundaryConditionsComparison) {
    struct ExtremeCase {
        double wind_speed, Qsw, Qnet, temp, soil, ea, ea_star, P;
        std::string description;
    };
    
    std::array<ExtremeCase,legacy_data::N_BOUNDARIES> extreme_cases = {
        ExtremeCase{100.0, 1000.0, 500.0, 50.0, 0.5, 3.0, 4.0, 110.0, "high_values"},
        {0.1, 10.0, 10.0, -10.0, 0.01, 0.1, 0.2, 90.0, "low_values"},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 101.3, "zeros"},
    };
    size_t i = 0; 
    for (const auto& test_case : extreme_cases) {
        // Set extreme values in mock data
        data._wind_speed = test_case.wind_speed;
        data.Qsw = test_case.Qsw;
        data._Q_net = test_case.Qnet;
        data._air_temperature = test_case.temp;
        data.soil_storage = test_case.soil;
        data.ea = test_case.ea;
        data.ea_star = test_case.ea_star;
        data.P_atm = test_case.P;
        
        // Test refactored version
        refactored_pm.execute(data); 
        // Compare behavior under extreme conditions
        EXPECT_NEAR(legacy_data::E_T_Boundaries[i], data._ET, 1e-6)
            << "ET mismatch under boundary conditions: " << test_case.description;
        EXPECT_NEAR(legacy_data::R_C_Boundaries[i], data._stomatal_resistance, 1e-6)
            << "Stomatal resistance mismatch under boundary conditions: " << test_case.description;
        ++i;
    }
    
    // Restore original values
    data._wind_speed = 2.0;
    data.Qsw = 500.0;
    data._Q_net = 100.0;
    data._air_temperature = 20.0;
    data.soil_storage = 0.3;
    data.ea = 1.5;
    data.ea_star = 2.0;
    data.P_atm = 101.3;
}

// Performance comparison test
TEST_F(PenmanMonteithTest, PerformanceComparison) {
    constexpr auto LEGACY_TIME = std::chrono::nanoseconds(200000);
    const int iterations = 1000;
    
    // Time refactored implementation
    auto start_refactored = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        refactored_pm.execute(data);
    }
    auto end_refactored = std::chrono::high_resolution_clock::now();
    
    auto original_duration = LEGACY_TIME;
    auto refactored_duration = end_refactored - start_refactored;
    EXPECT_LE(refactored_duration, original_duration) 
        << "Refactored version should not be significantly slower\n"
        << "It can sometimes fail because of an unlucky run, always run again";
}
