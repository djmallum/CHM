#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "surface_temperature.hpp"
#define diff 0.0001

// This is a fake mock
// Actual mocks require virtual functions but we want to avoid this overhead
// instead just overwrite the functions
class MockDataSurfaceTemperature : 
	public surface_temperature::data
{
public:
	
};

class SurfaceTemperatureTest : public ::testing::Test
{
protected:
	typedef surface_temperature_calculator<MockDataSurfaceTemperature>
		calculator;
	void SetUp override 
	{
		calculation = new calculator(d);	
	};
	
	MockDataSurfaceTemperature d;
	calculator* calculation;
	
};

TEST_F(SurafaceTemperatureTest,BareGroundTest)
{
	d
};




