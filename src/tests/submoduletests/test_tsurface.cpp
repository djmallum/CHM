#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "surface_temperature.hpp"
#include <cmath>
#include "CSVreader.hpp"
#define diff 0.0001

class data
{
public:
	//state variables
	double thaw_front_depth_last = 0.0;
	double snow_temperature_daily_mean = 0.0;

	//TODO Add this struct to the actual submodule class
	//struct W
	//{
	//	constexpr static double a = 0.77;
	//	constexpr static double b = 0.02;
	//	constexpr static double c = 7.0;
	//	constexpr static double d = 0.03;
	//};
	//[[no_unique_address]] W w;
	

	double swe() const
	{
		return _swe;
	};
	double thaw_front_depth() const
	{
		return _thaw_front_depth;
	};
	double air_temperature() const
	{
		return _air_temperature;
	};
	double net_radiation() const
	{
		return _net_radiation;
	};
	double snow_density() const
	{
		return _snow_density;
	};
	double ground_heat_flux() const
	{
		return _ground_heat_flux;
	};
	double snow_depth() const
	{
		return _snow_depth;
	};
	double snow_temperature() const
	{
		return _snow_temperature;
	};
	double surface_temperature() const
	{
		return _surface_temperature;	
	};
	double snow_thermal_conductivity() const
	{
		return _snow_thermal_conductivity;
	};
	
	//setters just to init these tests
	void swe(double& in)
	{
		_swe = in;
	};
	void thaw_front_depth(double& in)
	{
		_thaw_front_depth = in;
	};
	void air_temperature(double& in)
	{
		_air_temperature = in;
	};
	void net_radiation(double& in)
	{
		_net_radiation = in;
	};
	void snow_density(double& in)
	{
		_snow_density = in;
	};
	void ground_heat_flux(double& in)
	{
		_ground_heat_flux = in;
	};
	void snow_depth(double& in)
	{
		_snow_depth = in;
	};
	void snow_temperature(double& in)
	{
		_snow_temperature = in;
	};

	//actual outputs
	void surface_temperature(double& T)
	{
		_surface_temperature = T;
	};
	void snow_thermal_conductivity(double& Tc)
	{
		_snow_thermal_conductivity = Tc;
	};
private:
	//inputs
	double _swe;
	double _thaw_front_depth;
	double _air_temperature;
	double _net_radiation;
	double _snow_density;
	double _ground_heat_flux;
	double _snow_depth;
	double _snow_temperature;
	
	//outputs
	double _surface_temperature;	
	double _snow_thermal_conductivity;

};

class SurfaceTemperatureTest : public ::testing::Test
{
protected:
	void SetUp() override 
	{
		luce_tarboton = new luce_tarboton_surface_temperature(d);
		rcc = new RCC_surface_temperature(d);
	};
	
	void TearDown() override
	{
		delete luce_tarboton;
		delete rcc;
	};

	data d;
	RCC_surface_temperature<data>* rcc;
	luce_tarboton_surface_temperature<data>* luce_tarboton;
	
	void do_snow_test()
	{
		d.swe(#value);
		d.thaw_front_depth(#value);
		d.air_temperature(#value);
		d.net_radiation(#value);
		d.snow_density(#value);
		d.ground_heat_flux(#value);
		d.snow_depth(#value);
		d.snow_temperature(#value);

		calculation.execute();
	}

	void do_bare_ground_test()
	{
		d.thaw_front_depth(#value);
		d.air_temperature(#value);
		d.net_radiation(#value);

		rcc.execute();
	};
};

TEST_F(SurfaceTemperatureTest,BareGroundZeroThawDepthTest)
{
	d.thaw_front_depth(0.0);
	d.air_temperature(3.0);
	using RCC_surface_temperature;
	double CRHM_tsurf = 
		(a*d.air_temperature() + b*d.net_radiation())
		* std::atan(c*(d.thaw_front_depth_last + d)) 
		* 2.0/3.14159265;
	ASSERT_DOUBLE_EQ(d.surface_temperature(),CRHM_tsurf);
};

TEST_F(SurfaceTemperatureTest,BareGroundSmallThawDepthTest)
{
	d.thaw_front_depth_last = 1.0;
	d.thaw_front_depth(0.8);
	d.air_temperature(3.0);
	using RCC_surface_temperature;
	double CRHM_tsurf = 
		(a*d.air_temperature() + b*d.net_radiation())
		* std::atan(c*(d.thaw_front_depth_last + d)) 
		* 2.0/3.14159265;
	ASSERT_DOUBLE_EQ(d.surface_temperature(),CRHM_tsurf);
};

TEST_F(SurfaceTemperatureTest,BareGroundLargeThawDepthTest)
{
	d.thaw_front_depth_last = 1.0;
	d.thaw_front_depth(1.8);
	d.air_temperature(3.0);
	using RCC_surface_temperature;
	double CRHM_tsurf = 
		(a*d.air_temperature() + b*d.net_radiation())
		* std::atan(c*(d.thaw_front_depth() + d)) 
		* 2.0/3.14159265;
	ASSERT_DOUBLE_EQ(d.surface_temperature(),CRHM_tsurf);
};

TEST_F(SurfaceTemperatureTest,SnowCoveredLargeDensityLargeDailyTempTest)
{
	struct dense_const_test
	{
		constexpr static double a = 0.138;
		constexpr static double b = 1.01;
		constexpr static double c = 3.233;
		constexpr static double density_threshold_test = 156.0;
	};
	constexpr static double kg_per_m3_to_g_per_m3_test = 1.0 / 1000.0;
	constexpr static double temperaure_threshold_test = -70.0;

	using luce_tarboton_surface_temperature;	
	
	static_assert(kg_per_m3_to_g_per_m3_test == kg_per_m3_to_g_per_m3);
	static_assert(dense_const_test::density_threshold_test ==
			snow_density_threshold);
	static_assert(temperature_threshold_test == temperature_threshold);
	static_assert(dense_const_test::a == dense_const::a);
	static_assert(dense_const_test::b == dense_const::b);
	static_assert(dense_const_test::c == dense_const::c);

	d.air_temperature(-20.0);
	// temperature threshold is negative, so halfing makes it "large"
	d.snow_temperature_daily_mean(temperature_threshold/2.0);
	d.snow_density(snow_density_threshold*2.0);
	d.snow_depth(3.2);
	d.ground_heat_flux(32.333);

	luce_tarboton.execute();

	double CRHM_tsurf = 
		d.snow_temperature_daily_mean + 
		0.5*d.ground_heat_flux()*d.snow_depth() / 
		d.snow_thermal_conductivity();

	ASSERT_DOUBLE_EQ(d.surface_temperature(),CRHM_tsurf);

	double CRHM_tcond = 
		dense_const::a - 
		dense_const::b * d.snow_density() * kg_per_m3_to_g_per_cm3 + 
		dense_const::c * std::pow(d.snow_density() * kg_per_m3_to_g_per_cm3);

	ASSERT_DOUBLE_EQ(d.snow_thermal_conductivity(),CRHM_tcond);
};

TEST_F(SurfaceTemperatureTest,SnowCoveredSmallDensityLargeDailyTempTest)
{
	struct sparse_const_test
	{
		constexpr static double a = 0.023;
		constexpr static double b = 0.234;
		constexpr static double density_threshold_test = 156.0;
	};
	constexpr static double kg_per_m3_to_g_per_m3_test = 1.0 / 1000.0;
	constexpr static double temperaure_threshold_test = -70.0;

	using luce_tarboton_surface_temperature;	
	
	static_assert(kg_per_m3_to_g_per_m3_test == kg_per_m3_to_g_per_m3);
	static_assert(dense_const_test::density_threshold_test ==
			snow_density_threshold);
	static_assert(temperature_threshold_test == temperature_threshold);
	static_assert(sparse_const_test::a == sparse_const::a);
	static_assert(sparse_const_test::b == sparse_const::b);

	d.air_temperature(-20.0);
	// temperature threshold is negative, so halfing makes it "large"
	d.snow_temperature_daily_mean(temperature_threshold/2.0);
	d.snow_density(snow_density_threshold/2.0);
	d.snow_depth(3.2);
	d.ground_heat_flux(32.333);

	luce_tarboton.execute();

	double CRHM_tsurf = 
		d.snow_temperature_daily_mean + 
		0.5*d.ground_heat_flux()*d.snow_depth() / 
		d.snow_thermal_conductivity();

	ASSERT_DOUBLE_EQ(d.surface_temperature(),CRHM_tsurf);

	double CRHM_tcond = 
		sparse_const::a - 
		sparse_const::b * d.snow_density() * kg_per_m3_to_g_per_cm3;

	ASSERT_DOUBLE_EQ(d.snow_thermal_conductivity(),CRHM_tcond);
};

TEST_F(SurfaceTemperatureTest,SnowCoveredSmallDensitySmallDailyTempTest)
{
	struct sparse_const_test
	{
		constexpr static double a = 0.023;
		constexpr static double b = 0.234;
		constexpr static double density_threshold_test = 156.0;
	};
	constexpr static double kg_per_m3_to_g_per_m3_test = 1.0 / 1000.0;
	constexpr static double temperaure_threshold_test = -70.0;

	using luce_tarboton_surface_temperature;	
	
	static_assert(kg_per_m3_to_g_per_m3_test == kg_per_m3_to_g_per_m3);
	static_assert(dense_const_test::density_threshold_test ==
			snow_density_threshold);
	static_assert(temperature_threshold_test == temperature_threshold);
	static_assert(sparse_const_test::a == sparse_const::a);
	static_assert(sparse_const_test::b == sparse_const::b);

	d.air_temperature(-20.0);
	// temperature threshold is negative, so halfing makes it "large"
	d.snow_temperature_daily_mean(temperature_threshold*2.0);
	d.snow_density(snow_density_threshold/2.0);
	d.snow_depth(3.2);
	d.ground_heat_flux(32.333);

	luce_tarboton.execute();

	double CRHM_tsurf = d.air_temperature(); 

	ASSERT_DOUBLE_EQ(d.surface_temperature(),CRHM_tsurf);

	double CRHM_tcond = 
		sparse_const::a - 
		sparse_const::b * d.snow_density() * kg_per_m3_to_g_per_cm3;

	ASSERT_DOUBLE_EQ(d.snow_thermal_conductivity(),CRHM_tcond);
};

TEST_F(SurfaceTemperatureTest,TemperatureAccumulationTest)
{
	size_t N = 48;
	double period = 24.0;
	std::array<double,N> temperatures;
	size_t start = 0;

	bool is_new_day(size_t& ind) { return i % 24 == 0; };

	for (size_t i = 0; i < N; ++i)
	{
		if (is_new_day(i))
		{
			d.temperature_accumulated = 0.0;
			start = i;
		};
		double phase = 2.0 * 3.14159265 * 
			i / period;

		temperatures[i] = std::pow(std::sin(phase),2);
		d.air_temperature(temperatures[i]);
		accumulate_temperature.execute();
		double total = 0.0;
		for (size_t j = 0; j < i + 1; ++j)
		{
			total += temperatures[j];
		}
		ASSERT_DOUBLE_EQ(d.temperature_accumulated,total);
	}
};

TEST_F(SurfaceTemperatureTest,CRHMCompareTest)
{
	// TODO look at CRHM data, then find a point that works.
	// Hard code so not reliant on reader
	
};
