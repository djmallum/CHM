#pragma once
#include "base_step.hpp"
#include <cmath>

template<typename T>
concept RCC_data = requires(T& t,const double& out)
{
	{ t.thaw_front_depth() } -> std::convertible_to<double&>;

	{ t.air_temperature() } -> std::convertible_to<double&>;

	{ t.net_radiation() } -> std::convertible_to<double&>;
	
	{ t.surface_temperature(out) } -> std::same_as<void>;
};


template<RCC_data data>
class RCC_surface_temperature : public base_step<data>
{
public:
	explicit RCC_surface_temperature(data& _d) : base_step<data>(_d) {};
	~RCC_surface_temperature() {};

	void execute() override final;

private:
    double thaw_front_depth_last = 0.0;

	struct Constants
	{
		constexpr static double a = 0.77;
		constexpr static double b = 0.02;
		constexpr static double c = 7.0;
		constexpr static double d = 0.03;
	};
};


template<RCC_data data>
void RCC_surface_temperature<data>::execute()
{
    thaw_front_depth_last = std::max(thaw_front_depth_last,
            this->d.thaw_front_depth());

	double T = 
		(Constants::a * this->d.air_temperature() + Constants::b * this->d.net_radiation()) * 
		std::atan(Constants::c * (thaw_front_depth_last + Constants::d)) * 2.0 / 3.14156265;
	this->d.surface_temperature(T);
};
