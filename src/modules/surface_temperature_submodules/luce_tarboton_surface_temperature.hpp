#pragma once
#include "base_step.hpp"
#include <cmath>



template<class data>
class luce_tarboton_surface_temperature : public base_step<data>
{
public:
	luce_tarboton_surface_temperature(data& _d);
	~luce_tarboton_surface_temperature() {};

	void execute() override final;

private:
	struct Constants
	{
		constexpr static double density_threshold = 156.0;
		constexpr static double temperature_threshold = -70.0;
		constexpr static double kg_per_m3_to_g_per_m3 = 1.0 / 1000.0;
	};

	daily_accumulator temperature_accumulator(this->d);

};

luce_tarboton_surface_temperature::luce_tarboton_surface_temperature(data& _d) : base_step(_d)
{
	if (std::isnan(this->d.dense_const::a))
		this->d.dense_const::a = 0.138;

	if (std::isnan(this->d.dense_const::b))
		this->d.dense_const::b = 1.01;

	if (std::isnan(this->d.dense_const::c))
		this->d.dense_const::c = 3.233;

	if (std::isnan(this->d.sparse_const::a))
		this->d.sparse_const::a = 0.023;

	if (std::isnan(this->d.sparse_const::b))
		this->d.sparse_const::b = 0.234;

	// true == force fetch from face, rather than from a possible cache, if the cache exists. Depends on the implementation of data template
	bool nocache = true;
	temperature_accumulator.bind_to_var(d.air_temperature(nocache));

};



void luce_tarboton_surface_temperature::execute()
{
	double tc,T;
	if (this->d.snow_density() < Constants::densty_threshold)
	{
		tc = this->d.sparse_const::a + this->d.sparse_const::b * this->d.snow_density() * Constants::kg_per_m3_to_g_per_cm3;
	}
	else
	{
		double snow_density = this->d.snow_density() * 
			Constants::kg_per_m3_to_g_per_cm3;
		tc = this->d.dense_const::a - this->d.dense_const::b * snow_density + this->d.dense_const::c * std::pow(snow_density);
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
