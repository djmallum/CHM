#pragma once
#include "base_step.hpp"

template <typename D>
concept NetRadiationData = requires(D d) {
    // Check required member functions
	{ d.max_sun_hours() } -> std::convertible_to<double>;
	{ d.air_temperature() } -> std::convertible_to<double>;
    { d.vapour_pressure() } -> std::convertible_to<double>;
    { d.bright_sun_ratio() } -> std::convertible_to<double>;
    { d.actual_sun_hours() } -> std::convertible_to<double>;
    { d.incident_short_wave_clear() } -> std::convertible_to<double>;
    { d.diffuse_short_wave_clear() } -> std::convertible_to<double>;
    { d.albedo() } -> std::convertible_to<double>;

    // Check net_all_wave is callable with a double
    { d.net_all_wave(std::declval<double>()) } -> std::same_as<void>;
};

template<NetRadiationData data>
class net_radiation : public base_step
{
public:
	explicit net_radiation(data& _d, double& sb) : 
		base_step(_d),
		stefan_boltzmann(sb)
		{};

	~net_radiation() {};

	void override execute();
private:

	//empirical constants
	struct longwave
	{
		static constexpr double a = -0.85;
		static constexpr double b = 0.97;
	};

	struct brunt
	{
		static constexpr double a = -0.39;
		static constexpr double b = 0.093;
	};

	struct cloud_cover_long
	{
		static constexpr double a = 0.26;
		static constexpr double b = 0.81;
	};

	struct cloud_cover_short
	{
		static constexpr double a = 0.024;
		static constexpr double b = 0.974;
		static constexpr double c = 1.35;
	};

	struct Net
	{
		double long_wave = 0.0;
		double short_wave = 0.0;
	
		double get(double& albedo)
		{
			return short_wave*(1.0 - albedo) + long_wave;
		};
	};

	double& stefan_boltzmann;
};

/*
 * TODO
 *
 * - auto tabbing on new line not working sometime in shortwave
 * - create helper functions
 * - decide on consistent approach for static constants while also considering the helper functions
 * - Consider unit conversions found in CRHM for shortwave.
 *		- Consider that modules that use it (PenmanMonteith) likely expects the mm version not the W or MJ versions. 
 *		- How easy would a refactor be?
 */

void net_radiation::execute()
{
	Net net;

	if (this->d.max_sun_hours() > 0.0)
	{
		net.long_wave = longwave::a + longwave::b * stefan_boltzmann * std::pow(this->d.air_temperature() + 273.0,4) * (brunt::a + brunt::b * std::sqrt(this->d.vapour_pressure()) * (cloud_cover_long::a + cloud_cover_long::b * (this->d.bright_sun_ratio()));
	}
	else
	{
		net.long_wave = longwave::a;
	}

	if (this->d.actual_sun_hours() > 0.0 && this->d.max_sun_hours() > 0.0)
	{
		net.short_wave = (cloud_cover_short::a + cloud_cover_short::b * std::pow(this->d.bright_sun_ratio(),cloud_cover_short::c)) * this->d.incident_short_wave_clear() + (2.68 + 2.2*(this->d.bright_sun_ratio()) - 3.85*std::pow(this->d.bright_sun_ratio(),2)) * this->d.diffuse_short_wave_clear();
	}
	else
	{
		net.short_wave = (cloud_cover_short::a * this->d.incident_short_wave_clear() + 2.68 * this->d.diffuse_short_wave_clear());
	}

	this->d.net_all_wave(net.get(this->d.albedo()));

};


				


};

