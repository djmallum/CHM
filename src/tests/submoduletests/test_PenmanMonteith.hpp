#include <gtest/gtest.h>
#include "PenmanMonteith.hpp"
#include "Penman_Montieth.hpp"
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
        data.Cp = 1013.0;
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

        // Create original implementation
        original_pm = std::make_unique<PenmanMonteith>(
            data.LAI, data.LAImax, data.veg_Ht, data.wind_height, data.stomatal_res_min,
            data.soil_d, data.F_to_g, data._s_per_step, data.Cp, data.K, data.tension, data.pore_sz,
            data.theta_pwp, data.phi
        );
        

        // Create refactored implementation
        // refactored_pm = std::make_unique<RefactoredPenmanMonteith>(
        //     LAI, LAImax, veg_Ht, wind_height, stomatal_res_min,
        //     soil_d, F_to_g, s_per_step, Cp, K, tension, pore_sz,
        //     theta_pwp, phi
        // );
    }
    
    void TearDown() override {
        original_pm.reset();
        //refactored_pm.reset();
    }
    
    std::unique_ptr<PenmanMonteith> original_pm;
    Penman_monteith<MockPenmanData> refactored_pm;
    // std::unique_ptr<RefactoredPenmanMonteith> refactored_pm;
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
    
    std::vector<TestCase> test_cases = {
        {600.0, 0.5, 0.25, 25.0, "normal_conditions"},
        {0.0, 0.5, 0.25, 25.0, "night_time"},
        {600.0, 2.5, 0.25, 25.0, "high_vpd"},
        {600.0, 0.5, 0.1, 25.0, "dry_soil"},
        {600.0, 0.5, 0.25, 45.0, "high_temperature"},
        {600.0, 0.5, 0.25, 2.0, "low_temperature"}
    };
    
    Penman_monteith<MockPenmanData>::stomatal_resistance_jarvis r_s;

    for (const auto& test_case : test_cases) {
        // Set up test conditions using mock data
        data.Qsw = test_case.short_wave;
        data._air_temperature = test_case.temperature;
        data.soil_storage = test_case.soil_storage;
        data.ea_star = data.ea + test_case.vapour_pressure_diff;
        
        PM_vars vars(data.wind_speed(), data.Qsw, data.Q_net(), data.air_temperature(),
                    data.soil_storage, data.ea, data.ea_star, data.P_atm);
        PM_output original_output;
        original_pm->CalcEvapT(vars, original_output);
        
        r_s.calculate(data); 

        // Compare results
        EXPECT_NEAR(original_output.stomatal_resistance, data.stomatal_resistance(), 1e-12)
            << "Stomatal resistance mismatch for: " << test_case.description;
        EXPECT_NEAR(original_output.stomatal_resistance, data._stomatal_resistance, 1e-12)
            << "Stomatal resistance mismatch for: " << test_case.description;
        
        refactored_pm.execute(data); 
        // Compare results
        EXPECT_NEAR(original_output.stomatal_resistance, data.stomatal_resistance(), 1e-12)
            << "Stomatal resistance mismatch for (full system): " << test_case.description;
        EXPECT_NEAR(original_output.stomatal_resistance, data._stomatal_resistance, 1e-12)
            << "Stomatal resistance mismatch for (full system): " << test_case.description;
        EXPECT_NEAR(original_output.ET, data._ET, 1e-12)
            << "ET mismatch for (full system): " << test_case.description;
    }
    
}

// Comparison test for boundary conditions
TEST_F(PenmanMonteithTest, BoundaryConditionsComparison) {
    struct ExtremeCase {
        double wind_speed, Qsw, Qnet, temp, soil, ea, ea_star, P;
        std::string description;
    };
    
    std::vector<ExtremeCase> extreme_cases = {
        {100.0, 1000.0, 500.0, 50.0, 0.5, 3.0, 4.0, 110.0, "high_values"},
        {0.1, 10.0, 10.0, -10.0, 0.01, 0.1, 0.2, 90.0, "low_values"},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 101.3, "zeros"},
    };
    
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
        
        PM_vars vars(data.wind_speed(), data.Qsw, data.Q_net(), data.air_temperature(),
                    data.soil_storage, data.ea, data.ea_star, data.P_atm);
        PM_output original_output;
        original_pm->CalcEvapT(vars, original_output);
        
        // Test refactored version
        refactored_pm.execute(data); 
        // Compare behavior under extreme conditions
        EXPECT_NEAR(original_output.ET, data._ET, 1e-6)
            << "ET mismatch under boundary conditions: " << test_case.description;
        EXPECT_NEAR(original_output.stomatal_resistance, data._stomatal_resistance, 1e-6)
            << "Stomatal resistance mismatch under boundary conditions: " << test_case.description;
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
    const int iterations = 1000;
    PM_vars test_vars(data.wind_speed(), data.Qsw, data.Q_net(), data.air_temperature(),
                     data.soil_storage, data.ea, data.ea_star, data.P_atm);
    
    // Time original implementation
    auto start_original = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        PM_output output;
        original_pm->CalcEvapT(test_vars, output);
    }
    auto end_original = std::chrono::high_resolution_clock::now();
    
    // Time refactored implementation
    auto start_refactored = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        PM_output output;
        refactored_pm.execute(data);
    }
    auto end_refactored = std::chrono::high_resolution_clock::now();
    
    auto original_duration = end_original - start_original;
    auto refactored_duration = end_refactored - start_refactored;
    
    EXPECT_LE(refactored_duration, original_duration * 1.01) 
        << "Refactored version should not be significantly slower";
}
