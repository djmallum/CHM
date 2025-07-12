#include <gtest/gtest.h>
#include "net_radiation.hpp"

// Mock class that satisfies the NetRadiationData concept
class MockNetRadiationData {
public:
    double max_sun_hours_val = 0.0;
    double air_temp_val = 0.0;
    double vapour_pressure_val = 0.0;
    double bright_sun_ratio_val = 0.0;
    double actual_sun_hours_val = 0.0;
    double incident_short_wave_clear_val = 0.0;
    double diffuse_short_wave_clear_val = 0.0;
    double albedo_val = 0.0;
    double net_all_wave_output = 0.0;

    double max_sun_hours() const { return max_sun_hours_val; }
    double air_temperature() const { return air_temp_val; }
    double vapour_pressure() const { return vapour_pressure_val; }
    double bright_sun_ratio() const { return bright_sun_ratio_val; }
    double actual_sun_hours() const { return actual_sun_hours_val; }
    double incident_short_wave_clear() const { return incident_short_wave_clear_val; }
    double diffuse_short_wave_clear() const { return diffuse_short_wave_clear_val; }
    double albedo() const { return albedo_val; }
    void net_all_wave(double value) { net_all_wave_output = value; }
};

class NetRadiationTest : public ::testing::Test {
protected:
    MockNetRadiationData mock_data;
    net_radiation<MockNetRadiationData> net_rad{mock_data};

    void SetUp() override {
        // Default values that can be overridden in individual tests
        mock_data.max_sun_hours_val = 12.0;
        mock_data.air_temp_val = 20.0;
        mock_data.vapour_pressure_val = 10.0;
        mock_data.bright_sun_ratio_val = 0.5;
        mock_data.actual_sun_hours_val = 6.0;
        mock_data.incident_short_wave_clear_val = 100.0;
        mock_data.diffuse_short_wave_clear_val = 50.0;
        mock_data.albedo_val = 0.3;
    }
};

TEST_F(NetRadiationTest, ExecuteSetsNetAllWave) {
    net_rad.execute();
    EXPECT_NE(mock_data.net_all_wave_output, 0.0);
}

TEST_F(NetRadiationTest, StefanBoltzmannLaw) {
    // Test with 0°C (273.15K)
    double result = net_rad.stefan_boltzmann_law(0.0);
    const double sigma = 5.670374418e-8;
    double expected = sigma * std::pow(273.15, 4);
    EXPECT_DOUBLE_EQ(result, expected);

    // Test with 20°C (293.15K)
    result = net_rad.stefan_boltzmann_law(20.0);
    expected = sigma * std::pow(293.15, 4);
    EXPECT_DOUBLE_EQ(result, expected);
}

TEST_F(NetRadiationTest, NetClearLongWave) {
    double result = net_rad.net_clear_long_wave();
    // Expected value calculated manually based on the formula
    double sigma = 5.670374418e-8;
    double T_kelvin = 20.0 + 273.15;
    double expected = sigma * std::pow(T_kelvin, 4) * (-0.39 + 0.093 * std::sqrt(10.0));
    EXPECT_DOUBLE_EQ(result, expected);
}

TEST_F(NetRadiationTest, CloudCoverLong) {
    mock_data.bright_sun_ratio_val = 0.5;
    double result = net_rad.cloud_cover_long();
    EXPECT_DOUBLE_EQ(result, 0.26 + 0.81 * 0.5);

    mock_data.bright_sun_ratio_val = 0.0;
    result = net_rad.cloud_cover_long();
    EXPECT_DOUBLE_EQ(result, 0.26);

    mock_data.bright_sun_ratio_val = 1.0;
    result = net_rad.cloud_cover_long();
    EXPECT_DOUBLE_EQ(result, 0.26 + 0.81);
}

TEST_F(NetRadiationTest, DirectRadiation) {
    mock_data.bright_sun_ratio_val = 0.5;
    double result = net_rad.direct_radiation();
    EXPECT_DOUBLE_EQ(result, 0.974 * std::pow(0.5, 1.35));
}

TEST_F(NetRadiationTest, DiffuseRadiation) {
    mock_data.bright_sun_ratio_val = 0.5;
    double result = net_rad.diffuse_radiation();
    EXPECT_DOUBLE_EQ(result, 2.2 * 0.5 - 3.85 * std::pow(0.5, 2));
}

TEST_F(NetRadiationTest, GetLongWaveWithSun) {
    mock_data.max_sun_hours_val = 12.0;
    double result = net_rad.get_long_wave();
    double net_clear = net_rad.net_clear_long_wave();
    double cloud_cover = net_rad.cloud_cover_long();
    double expected = -0.85 + 0.97 * net_clear * cloud_cover;
    EXPECT_DOUBLE_EQ(result, expected);
}

TEST_F(NetRadiationTest, GetLongWaveNoSun) {
    mock_data.max_sun_hours_val = 0.0;
    double result = net_rad.get_long_wave();
    EXPECT_DOUBLE_EQ(result, -0.85); // Just the C1 constant
}

TEST_F(NetRadiationTest, GetShortWaveWithSun) {
    mock_data.actual_sun_hours_val = 6.0;
    mock_data.max_sun_hours_val = 12.0;
    double result = net_rad.get_short_wave();
    
    double direct = net_rad.direct_radiation();
    double diffuse = net_rad.diffuse_radiation();
    double expected = (0.024 + direct) * mock_data.incident_short_wave_clear_val + 
                     (2.68 + diffuse) * mock_data.diffuse_short_wave_clear_val;
    
    EXPECT_DOUBLE_EQ(result, expected);
}

TEST_F(NetRadiationTest, GetShortWaveNoSun) {
    mock_data.actual_sun_hours_val = 0.0;
    double result = net_rad.get_short_wave();
    double expected = 0.024 * mock_data.incident_short_wave_clear_val + 
                     2.68 * mock_data.diffuse_short_wave_clear_val;
    EXPECT_DOUBLE_EQ(result, expected);
}

TEST_F(NetRadiationTest, SetNetAllWave) {
    net_radiation<MockNetRadiationData>::Net net;
    net.long_wave = 50.0;
    net.short_wave = 100.0;
    mock_data.albedo_val = 0.2;
    
    net_rad.set_net_all_wave(net);
    
    double expected = 100.0 * (1.0 - 0.2) + 50.0;
    EXPECT_DOUBLE_EQ(mock_data.net_all_wave_output, expected);
}

TEST_F(NetRadiationTest, NetStructGet) {
    net_radiation<MockNetRadiationData>::Net net;
    net.long_wave = 30.0;
    net.short_wave = 70.0;
    
    double result = net.get(0.5);
    EXPECT_DOUBLE_EQ(result, 70.0 * 0.5 + 30.0);
}
