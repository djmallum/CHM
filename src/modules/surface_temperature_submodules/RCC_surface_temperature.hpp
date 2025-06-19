#pragma once
#include "base_step.hpp"
#include <cmath>
struct Constants
{
	constexpr static double a = 0.77;
	constexpr static double b = 0.02;
	constexpr static double c = 7.0;
	constexpr static double d = 0.03;
};

template<class data>
class RCC_surface_temperature : public base_step<data>
{
public:
	RCC_surface_temperature(data& _d) : base_step(_d) {};
	~RCC_surface_temperature() {};

	void execute() override final;
};

void RCC_surface_temperature::execute()
{
	if (this->d.thaw_front_depth() > this->d.thaw_front_depth_last)
		this->d.thaw_front_depth_last = this->d.thaw_front_depth();

	double T = 
		(Constants::a * this->d.air_temperature() + Constants::b * this->d.net_radiation()) * 
		std::atan(Constants::c * (this->d.thaw_front_depth_last + Constants::d)) * 2.0 / 3.14156265;
	this->d.surface_temperature(T);
};
