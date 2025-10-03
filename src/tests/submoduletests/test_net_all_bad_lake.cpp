#include "net_all_bad_lake.hpp"
#include "gtest/gtest.h"

// Mock class satisfying NetAllData concept
struct MockNetAllData {
    double albedo_value;
    double short_wave_value;
    double net_all_wave_result;

    double albedo() const { return albedo_value; }
    double incoming_short_wave() const { return short_wave_value; }
    void net_all_wave(const double value) { net_all_wave_result = value; }
};

// Test fixture
class NetAllBadLakeTest : public ::testing::Test {
protected:
    net_all_bad_lake<MockNetAllData> step;
    MockNetAllData mock_data;
};

// Test case: Verify execute() computes net_all_wave correctly
TEST_F(NetAllBadLakeTest, ComputesNetAllWaveCorrectly) {
    // Setup
    mock_data.albedo_value = 0.1;          // 10% albedo
    mock_data.short_wave_value = 100.0;    // 100 W/m² incoming shortwave

    // Action
    step.execute(mock_data);

    // Verification
    const double expected = -2.24 + 0.651 * 100.0 * (1 - 0.1);
    EXPECT_DOUBLE_EQ(mock_data.net_all_wave_result, expected);
}

// Test case: Edge case (zero incoming shortwave)
TEST_F(NetAllBadLakeTest, HandlesZeroShortWave) {
    mock_data.albedo_value = 0.3;
    mock_data.short_wave_value = 0.0;

    step.execute(mock_data);

    EXPECT_DOUBLE_EQ(mock_data.net_all_wave_result, -2.24);  // Only 'a' term remains
}

