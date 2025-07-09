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
	explicit net_radiation(data& _d) : base_step(_d) {};

	~net_radiation() {};

	void execute() final;
private:

	struct Net
	{
		double long_wave = 0.0;
		double short_wave = 0.0;
	
		double get(const double& albedo) const
		{
			return short_wave*(1.0 - albedo) + long_wave;
		};
	};

	double get_long_wave() const;
	double get_short_wave() const;
	void set_net_all_wave(Net& net);

	//long-wave
	double stefan_boltzmann_law(const double& T) const;
	double incoming_long_wave() const;
	double cloud_cover_long() const;

	//short-wave
	double direct_radiation() const;
	double diffuse_radiation() const;

};

/*
 * TODO
 *
 * - Consider unit conversions found in CRHM for shortwave.
 *		- Consider that modules that use it (PenmanMonteith) likely expects the mm version not the W or MJ versions. 
 *		- How easy would a refactor be?
 */

void net_radiation::execute()
{
	Net net;
	
	net.long_wave = get_long_wave();

	net.short_wave = get_short_wave();

	set_net_all_wave();

};

double net_radiation::get_long_wave() const
{
	static const double C = -0.85;
	static const double brunt_a = -0.39;
	static const double brunt_b = 0.093;
	static const double cloud_cover_a = 0.26;
	static const double cloud_cover_b = 0.81;
	double long_wave = 0.0;

	if (this->d.max_sun_hours() > 0.0)
	{
		long_wave = C + incoming_long_wave() * cloud_cover_long();
	}
	else
	{
		long_wave = C;
	}

	return long_wave;
};

double net_radiation::get_short_wave() const
{
	static const double a_direct 0.024;
	static const double a_diffuse = 2.68;
	double short_wave = 0.0;

	if (this->d.actual_sun_hours() > 0.0 && this->d.max_sun_hours() > 0.0)
	{
		short_wave = (a_direct + direct_radiation()) * this->d.incident_short_wave_clear() 
			+ (a_diffuse + diffuse_radiation()) * this->d.diffuse_short_wave_clear();
	}
	else
	{
		short_wave = a_direct * this->d.incident_short_wave_clear()
			+ a_diffuse * this->d.diffuse_short_wave_clear();
	}

	return short_wave;
};

void net_radiation::set_net_all_wave(Net& net)
{	
	this->d.net_all_wave(net.get(this->d.albedo()));
};

double net_radiation::stefan_boltzmann_law(const double& T) const
{
	// Stefan-Boltzmann constant
	static const double sigma = 5.670374418e-8;
	static const double emissivity = 0.97;
	static const double convert_to_kelvin = 273.15;

	return emissivity * sigma * 
		std::pow(T + convert_to_kelvin,4);
};
	
double net_radiation::incoming_long_wave() const
{
	// brunt equation to compute long-wave
	// from the atmosphere
	static const double a = -0.39;
	static const double b = 0.093;

	return stefan_boltzmann_law(this->d.air_temperature()) * 
		(a + b * std::sqrt(this->d.vapour_pressure()));
};

double net_radiation::cloud_cover_long() const
{
	//computes effect of cloud cover on long wave radiation
	
	static const double a = 0.26;
	static const double b = 0.81;

	return a + b * this->d.bright_sun_ratio();	
};

double net_radiation::direct_radiation() const
{
	/* 
	 * Q_{d}/Q_{d0} - a
	 *
	 * Q_d: direct short-wave radiation
	 * Q_d0: direct short-wave radiation with clear skies
	 *
	 * This function returns an empirical formula for this ratio
	 * 
	 * For convenience, the constant a is defined in get_short_wave()
	 * as a_direct.
	 */

	static const double b = 0.974;
	static const double c = 1.35;

	return b * 
		std::pow(this->d.bright_sun_ratio(),c);
};

double net_radiation::diffuse_radiation() const
{
	// See comment above for direct version
	// mention a_diffuse
	static const double b = 2.2;
	static const double c = 3.85;

	return 2.2 * this->d.bright_sun_ratio()
		- 3.85 * std::pow(this->d.bright_sun_ratio(),2);
};
