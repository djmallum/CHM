#pragma once
#include "base_step.hpp"
#include <cmath>

template<class T>
concept luce_tarboton_data = requires(T& t,const double& out)
{
	// gets the air temperature
	{ t.air_temperature() } -> std::convertible_to<double&>;

	{ t.snow_depth() } -> std::convertible_to<double&>;

	{ t.snow_density() } -> std::convertible_to<double&>;

	{ t.ground_heat_flux() } -> std::convertible_to<double&>;

    { t.daily_mean_temperature() } -> std::convertible_to<const double>;

	{ t.surface_temperature(out) } -> std::same_as<void>;

	{ t.snow_thermal_conductivity(out) } -> std::same_as<void>;
};

template<luce_tarboton_data data>
class luce_tarboton_surface_temperature : public base_step<data>
{
public:
	explicit luce_tarboton_surface_temperature() {};
	~luce_tarboton_surface_temperature() {};

	void execute(data& d) override final;

private:

    static constexpr inline double kg_per_m3_to_g_per_cm3 = 1.0 / 1000.0;
    double low_density_thermal_conductivity(data& d);  
    double high_density_thermal_conductivity(data& d);
};

template<luce_tarboton_data data>
void luce_tarboton_surface_temperature<data>::execute(data& d)
{
    static const double density_threshold = 156.0;
    static const double temperature_threshold = -70.0;

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
    static const double a = 0.023;
    static const double b = 0.234;

    return a + b * d.snow_density() * kg_per_m3_to_g_per_cm3;
};

template<luce_tarboton_data data>
double luce_tarboton_surface_temperature<data>::high_density_thermal_conductivity(data& d)
{
    static const double a = 0.138;
    static const double b = 1.01;
    static const double c = 3.233;

    double snow_density = d.snow_density() * kg_per_m3_to_g_per_cm3;

    return a - b * snow_density + c * std::pow(snow_density,2);
};

