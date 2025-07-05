#pragma once
#include "base_step.hpp"
#include "daily_accumulator.hpp"
#include <cmath>

template<class T>
concept luce_tarboton_data = requires(T& t,const double& out)
{
	// gets the air temperature
	{ t.air_temperature() } -> std::convertible_to<double&>;

	{ t.snow_depth() } -> std::convertible_to<double&>;

	{ t.snow_density() } -> std::convertible_to<double&>;

	{ t.ground_heat_flux() } -> std::convertible_to<double&>;

	{ t.surface_temperature(out) } -> std::same_as<void>;

	{ t.snow_thermal_conductivity(out) } -> std::same_as<void>;
};



template<luce_tarboton_data data>
class luce_tarboton_surface_temperature : public base_step<data>
{
public:
	explicit luce_tarboton_surface_temperature(data& _d);
	~luce_tarboton_surface_temperature() {};

	void execute() override final;

private:
	struct Constants
	{
		constexpr static double density_threshold = 156.0;
		constexpr static double temperature_threshold = -70.0;
		constexpr static double kg_per_m3_to_g_per_m3 = 1.0 / 1000.0;
	};

	struct dense_const
	{
		constexpr static double a = 0.138;
		constexpr static double b = 1.01;
		constexpr static double c = 3.233;
	};

	struct sparse_const
	{
		constexpr static double a = 0.023;
		constexpr static double b = 0.234;
	};

	daily_accumulator<data> temperature_accumulator;

};

template<luce_tarboton_data data>
luce_tarboton_surface_temperature<data>::luce_tarboton_surface_temperature(data& _d) 
    : base_step<data>(_d), temperature_accumulator(_d)
{
	temperature_accumulator.bind_to_var(this->d.air_temperature());
};



template<luce_tarboton_data data>
void luce_tarboton_surface_temperature<data>::execute()
{
	double tc,T;
	if (this->d.snow_density() < Constants::densty_threshold)
	{
		tc = sparse_const::a + this->d.sparse_const::b * this->d.snow_density() * Constants::kg_per_m3_to_g_per_cm3;
	}
	else
	{
		double snow_density = this->d.snow_density() * 
			Constants::kg_per_m3_to_g_per_cm3;
		tc = dense_const::a - dense_const::b * snow_density + dense_const::c * std::pow(snow_density,2);
	}

	this->d.snow_thermal_conductivity(tc);

	if (temperature_accumulator.get_last_mean() < Constants::temperature_threshold)
		T = this->d.air_temperature();
	else
	{
		T = temperature_accumulator.get_last_mean() + 
			this->d.ground_heat_flux() * 0.5 * this->d.snow_depth() / tc;
	}

	this->d.surface_temperature(T);

	temperature_accumulator.execute();
};
