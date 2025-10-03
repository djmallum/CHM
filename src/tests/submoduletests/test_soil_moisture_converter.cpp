#include <gtest/gtest.h>
#include "soil_moisture_converter.hpp"

// Mock data class that satisfies the SoilMoistureConverterData concept
class MockSoilData {
public:
    bool storage_is_total_moisture_return = false;
    double fractional_cutoff_return = 0.0;
    double soil_storage_return = 0.0;
    double soil_storage_max_return = 0.0;
    double porosity_return = 0.0;
    double last_volumetric_content = 0.0;
    double last_saturation = 0.0;

    bool storage_is_total_moisture() { return storage_is_total_moisture_return; }
    double fractional_cutoff() { return fractional_cutoff_return; }
    double soil_storage() { return soil_storage_return; }
    double soil_storage_max() { return soil_storage_max_return; }
    double porosity() { return porosity_return; }
    double volumetric_moisture_content() { return last_volumetric_content; }
    void volumetric_moisture_content(const double value) { last_volumetric_content = value; }
    void saturation(const double value) { last_saturation = value; }
};

class SoilMoistureConverterTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Default test values
        mock_data.porosity_return = 0.5;
        mock_data.soil_storage_max_return = 100.0;
    }

    MockSoilData mock_data;
    soil_moisture_converter<MockSoilData> converter;
};

TEST_F(SoilMoistureConverterTest, Execute_TotalMoisturePath) {
    mock_data.storage_is_total_moisture_return = true;
    mock_data.soil_storage_return = 50.0;  // 50% of max storage

    converter.execute(mock_data);

    // Verify volumetric moisture content calculation
    EXPECT_DOUBLE_EQ(mock_data.last_volumetric_content, 0.25);  // 0.5 * (50/100)

    // Verify saturation calculation
    EXPECT_DOUBLE_EQ(mock_data.last_saturation, 0.5);  // 0.25 / 0.5
}

TEST_F(SoilMoistureConverterTest, Execute_FractionalMoisturePath) {
    mock_data.storage_is_total_moisture_return = false;
    mock_data.fractional_cutoff_return = 0.1;  // Wilt point at 10%
    mock_data.soil_storage_return = 45.0;      // 50% of available moisture (90 units available)

    converter.execute(mock_data);

    // Verify calculation: 0.1 + (45/100)*(0.5-0.1)
    EXPECT_DOUBLE_EQ(mock_data.last_volumetric_content, 0.1 + 0.45*0.4);

    // Verify saturation calculation
    EXPECT_DOUBLE_EQ(mock_data.last_saturation, mock_data.last_volumetric_content / mock_data.porosity_return);
}

TEST_F(SoilMoistureConverterTest, Execute_EmptySoil) {
    mock_data.soil_storage_return = 0.0;
    mock_data.fractional_cutoff_return = 0.1; 
    converter.execute(mock_data);
    
    EXPECT_DOUBLE_EQ(mock_data.last_volumetric_content, 0.1);  // Just the fractional cutoff
    EXPECT_DOUBLE_EQ(mock_data.last_saturation, 0.1 / 0.5);
}

TEST_F(SoilMoistureConverterTest, Execute_FullSoil) {
    mock_data.soil_storage_return = 100.0;
    mock_data.fractional_cutoff_return = 0.1;
    
    converter.execute(mock_data);
    
    EXPECT_DOUBLE_EQ(mock_data.last_volumetric_content, 0.5);  // Full porosity
    EXPECT_DOUBLE_EQ(mock_data.last_saturation, 1.0);          // Fully saturated
}

TEST_F(SoilMoistureConverterTest, CheckCutoffValidity_ValidRange) {
    // Should not throw for values in [0, 1]
    EXPECT_NO_THROW(converter.execute(mock_data));
    
    mock_data.fractional_cutoff_return = 1.0;
    EXPECT_NO_THROW(converter.execute(mock_data));
}

TEST_F(SoilMoistureConverterTest, CheckCutoffValidity_InvalidLow) {
    mock_data.fractional_cutoff_return = -0.1;
    EXPECT_THROW(converter.execute(mock_data), std::logic_error);
}

TEST_F(SoilMoistureConverterTest, CheckCutoffValidity_InvalidHigh) {
    mock_data.fractional_cutoff_return = 1.1;
    EXPECT_THROW(converter.execute(mock_data), std::logic_error);
}

