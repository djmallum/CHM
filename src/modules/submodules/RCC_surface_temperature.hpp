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

    { t.thaw_front_depth_last } -> std::convertible_to<double>;
};


template<RCC_data data>
class RCC_surface_temperature : public base_step<data>
{
public:
	explicit RCC_surface_temperature() {};
	~RCC_surface_temperature() {};

	void execute(data& d) override final;

private:
	struct Constants
	{
			};
};


template<RCC_data data>
void RCC_surface_temperature<data>::execute()
{
    constexpr static double a = 0.77;
    constexpr static double b = 0.02;
    constexpr static double c = 7.0;
    constexpr static double d = 0.03;

    d.thaw_front_depth_last = std::max(d.thaw_front_depth_last,
            this->d.thaw_front_depth());

	double T = 
		(a * this->d.air_temperature() + b * this->d.net_radiation()) * 
		std::atan(c * (thaw_front_depth_last + d)) * 2.0 / 3.14156265;
	this->d.surface_temperature(T);
};
