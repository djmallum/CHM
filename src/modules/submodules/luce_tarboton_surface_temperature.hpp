#pragma once
#include "base_step.hpp"
#include <cmath>
#include <concepts>

/*
 * Equation (11) of [1], solved for T_s
 *
 * [1] C. H. Luce and D. G. Tarboton, “Evaluation of alternative formulae for calculation of surface temperature in snowmelt models using frequency analysis of temperature observations,” Hydrol. Earth Syst. Sci., 2010.
 *
 * Small note: the equations slightly differ. Instead of d, a damping depth, we have snow_depth/2. Not sure why, but this is what CRHM does.
 */  
template<class T>
concept luce_tarboton_data = requires(T& t)
{
    // Inputs	
    { t.air_temperature() } -> std::floating_point;

	{ t.snow_depth() } -> std::floating_point;

	{ t.snow_density() } -> std::floating_point;

	{ t.ground_heat_flux() } -> std::floating_point;

    { t.daily_mean_temperature() } -> std::floating_point;

    // Outputs
	{ t.surface_temperature(std::declval<const double>())} -> std::same_as<void>;

	{ t.snow_thermal_conductivity(std::declval<const double>()) } -> std::same_as<void>;
};

template<luce_tarboton_data data>
class luce_tarboton_surface_temperature : public base_step<data>
{
public:
	explicit luce_tarboton_surface_temperature() {};
	~luce_tarboton_surface_temperature() {};

	void execute(data& d) override final;

private:

    static constexpr inline auto kg_per_m3_to_g_per_cm3 = 1.0 / 1000.0;
    double low_density_thermal_conductivity(data& d);  
    double high_density_thermal_conductivity(data& d);
};

template<luce_tarboton_data data>
void luce_tarboton_surface_temperature<data>::execute(data& d)
{
    static constexpr auto density_threshold = 156.0;
    static constexpr auto temperature_threshold = -70.0;

	double tc,T;
	if (d.snow_density() < density_threshold)
	{
        tc = low_density_thermal_conductivity(d);
	}
	else
	{
		tc = high_density_thermal_conductivity(d); 
	}

	d.snow_thermal_conductivity(tc);

	if (d.daily_mean_temperature() < temperature_threshold)
		T = d.air_temperature();
	else
	{
		T = d.daily_mean_temperature() + 
			d.ground_heat_flux() * 0.5 * d.snow_depth() / tc;
	}

	d.surface_temperature(T);
};

template<luce_tarboton_data data>
double luce_tarboton_surface_temperature<data>::low_density_thermal_conductivity(data& d)
{
    static constexpr auto a = 0.023;
    static constexpr auto b = 0.234;

    return a + b * d.snow_density() * kg_per_m3_to_g_per_cm3;
};

template<luce_tarboton_data data>
double luce_tarboton_surface_temperature<data>::high_density_thermal_conductivity(data& d)
{
    static constexpr auto a = 0.138;
    static constexpr auto b = 1.01;
    static constexpr auto c = 3.233;

    double snow_density = d.snow_density() * kg_per_m3_to_g_per_cm3;

    return a - b * snow_density + c * std::pow(snow_density,2);
};

