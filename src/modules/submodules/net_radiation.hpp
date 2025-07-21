#pragma once
#include "base_step.hpp"

template <typename D>
concept NetRadiationData = requires(D d) {
    // Check required member functions
	// Maximum sunlight hours for the day
	{ d.max_sun_hours() } -> std::convertible_to<double>;
	
	// Air temperature (Celsius)
	{ d.air_temperature() } -> std::convertible_to<double>;
    
	// Vapour pressure in mb
	{ d.vapour_pressure() } -> std::convertible_to<double>;
    
	// actual_sun_hours/max_sun_hours
	{ d.bright_sun_ratio() } -> std::convertible_to<double>;
    
	// Real sun hours for that day
	{ d.actual_sun_hours() } -> std::convertible_to<double>;
    
	// Incident short wave radiation if there were no clouds
	{ d.direct_short_wave_clear() } -> std::convertible_to<double>;
    
	// Diffuse radiation without clouds
	{ d.diffuse_short_wave_clear() } -> std::convertible_to<double>;
    
	// Surface albedo
	{ d.albedo() } -> std::convertible_to<double>;
	
	// Check net_all_wave is callable with a double
	// This sets the output
    { d.net_all_wave(std::declval<double>()) } -> std::same_as<void>;
};

template<NetRadiationData data>
class net_radiation : public base_step<data>
{
public:
	explicit net_radiation() {};
	friend class NetRadiationTest;
	~net_radiation() {};

	void execute(data& d) final;
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

	double get_long_wave(data& d) const;
	double get_short_wave(data& d) const;
	void set_net_all_wave(data& d,Net& net);

	//long-wave
	double stefan_boltzmann_law(const double& T) const;
	double net_clear_long_wave(data& d) const;
	double cloud_cover_long(data& d) const;

	//short-wave
	double direct_radiation(data& d) const;
	double diffuse_radiation(data& d) const;

};

/*
 * TODO
 *
 * - Consider unit conversions found in CRHM for shortwave.
 *		- Consider that modules that use it (PenmanMonteith) likely expects the mm version not the W or MJ versions. 
 *		- How easy would a refactor be?
 */

template<NetRadiationData data>
void net_radiation<data>::execute(data& d)
{
	/*
	 * Computes long-wave and short-wave radiation and outputs.
	 *
	 * get_long_wave() is in units of MJ/m^2/day.
	 *
	 * MJ_per_day_to_W converts these units to W/m^2.
	 *
	 * 1e6 is J/MJ
	 *
	 * 86400 is seconds/day
	 */ 
	static const double MJ_per_day_to_W = 1e6/86400;
	Net net;
	
	net.long_wave = get_long_wave(data& d) * 
		MJ_per_day_to_W;

	net.short_wave = get_short_wave(data& d);

	set_net_all_wave(net);
};

template<NetRadiationData data>
double net_radiation<data>::get_long_wave(data& d) const
{
	/*
	 * Computes long-wave radiation based on a semi-empirical equation 
	 * for net long-wave radiation on clear skies following
	 *
	 * Q_ln/Q_ln0 = cloud_cover_long()
	 * 
	 * and Q_ln0 is 
	 *
	 * net_clear_long_wave()
	 *
	 * so long_wave is cloud_cover_long()*Q_ln0 with best-fit constants 
	 */

	// Best-fit coefficients
	static const double C1 = -0.85;
	static const double C2 = 0.97; //Not an emissivity
								   //Looks like one in CRHM, but isn't
	double long_wave = 0.0;

	if (d.max_sun_hours() > 0.0)
	{
		long_wave = C1 + C2*net_clear_long_wave(data& d) * cloud_cover_long(data& d);
	}
	else
	{
		// C is not obviously defined in the literature
		// I believe it is an intercept from linear fitting
		// Therefire, this is the long_wave when there is no sun.
		long_wave = C1;
	}

	return long_wave;
};

template<NetRadiationData data>
double net_radiation<data>::get_short_wave(data& d) const
{
	/*
	 * Computes short-wave radiation sub-divided into direct and diffuse radiation
	 * 
	 * Diffuse and direct obey the following form:
	 *
	 * Q_df/Q_df0 = diffuse_radiation() + a_diffuse
	 *
	 * Q_dr/Q_dr0 = direct_radiation() + a_direct
	 *
	 * where Q_dx is the actual diffuse (f) or direct (r) radiation, and Q_dx0 is
	 * the equivalent assuming clear skies.
	 *
	 * We divide a_diffuse and a_direct constants from the functions diffuse/direct_radiation()
	 * for convenience in this if-statement and to avoid repeating a_diffuse definitions.
	 *
	 * It should be known that these constants should be thought of as part of the functions.
	 *
	 *
	 */ 
	static const double a_direct 0.024;
	static const double a_diffuse = 2.68;
	double short_wave = 0.0;

	if (d.actual_sun_hours() > 0.0 && d.max_sun_hours() > 0.0)
	{
		short_wave = (a_direct + direct_radiation(data& d)) * d.direct_short_wave_clear() 
			+ (a_diffuse + diffuse_radiation(data& d)) * d.diffuse_short_wave_clear();
	}
	else
	{
		short_wave = a_direct * d.direct_short_wave_clear()
			+ a_diffuse * d.diffuse_short_wave_clear();
	}

	return short_wave;
};

template<NetRadiationData data>
void net_radiation<data>::set_net_all_wave(data& d,Net& net)
{
	double net_all = net.get(d.albedo()); 
	
	d.net_all_wave(net_all);
};

template<NetRadiationData data>
double net_radiation<data>::stefan_boltzmann_law(const double& T) const
{
	/*
	 * Computes \sigma T^4
	 *
	 * Input: double T (Celius)
	 *
	 * T is converted to Kelvin in the equation.
	 */ 
	// Stefan-Boltzmann constant
	static const double sigma = 5.670374418e-8;
	static const double convert_to_kelvin = 273.15;

	return sigma * std::pow(T + convert_to_kelvin,4);
};

template<NetRadiationData data>
double net_radiation<data>::net_clear_long_wave(data& d) const
{
	/*
	 * Net long-wave radiation
	 *
	 * Derived assuming the Brunt equation for downwards radiation
	 * and treating the snow as a black-body.
	 *
	 * Then, assuming the snow temperature at the surface is the 
	 * air temperature, those terms are combined.
	 *
	 * The constant b is empirical.
	 *
	 * The constant a is semi-empirical. It is the difference of 
	 * empirical constant 0.58 and the approximate emissivity of snow.
	 *
	 * Details in 
	 *
	 * Gray, Donald Maurice, and P. G. Landine. Development and 
	 * performance evaluation of an energy budget snowmelt model. 
	 * Research Management Division, Alberta Environment, 1988.
	 *
	 * NOTE: Some of the values are not perfect there, so also 
	 * consider 
	 *
	 * Gray, D. M., and P. G. Landine. "An energy-budget snowmelt 
	 * model for the Canadian Prairies." Canadian Journal of Earth 
	 * Sciences 25.8 (1988): 1292-1303.
	 *
	 * Both are important for a complete understanding.
	 */ 
	// brunt equation to compute long-wave
	// from the atmosphere
	static const double a = -0.39;
	static const double b = 0.093;

	return stefan_boltzmann_law(d.air_temperature()) * 
		(a + b * std::sqrt(d.vapour_pressure()));
};

template<NetRadiationData data>
double net_radiation<data>::cloud_cover_long(data& d) const
{
	/*
	 * Equal to 
	 *
	 * Q_ln/Q_ln0
	 */ 
	
	static const double a = 0.26;
	static const double b = 0.81;

	return a + b * d.bright_sun_ratio();	
};

template<NetRadiationData data>
double net_radiation<data>::direct_radiation(data& d) const
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
		std::pow(d.bright_sun_ratio(),c);
};

template<NetRadiationData data>
double net_radiation<data>::diffuse_radiation(data& d) const
{
	/* 
	 * Q_{d}/Q_{d0} - a
	 *
	 * Q_d: diffuse short-wave radiation
	 * Q_d0: diffuse short-wave radiation with clear skies
	 *
	 * This function returns an empirical formula for this ratio
	 * 
	 * For convenience, the constant a is defined in get_short_wave()
	 * as a_direct.
	 */
	static const double b = 2.2;
	static const double c = 3.85;

	return b * d.bright_sun_ratio()
		- c * std::pow(d.bright_sun_ratio(),2);
};
