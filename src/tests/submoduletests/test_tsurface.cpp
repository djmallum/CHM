#include <gmock/gmock.h>
#include <gtest/gtest.h>
#define diff 0.0001

class SurfaceTemperatureTest : public ::testing::Test
{
protected:
	void SetUp override 
	{

	};

};

class MockDataSurfaceTemperature : 
	public surface_temperature::data
{
public:




