#pragma once
#include "base_step.hpp"
#include <algorithm>
#include <concepts>

/*
 * Equation (7) from [1].
 * 
 * Computes the temperature of the ground surface using the Radiative-Conductive-Convection Approach, used to estimate an empiracal formula using linear regression.
 *
 * [1] T. J. Williams, J. W. Pomeroy, J. R. Janowicz, S. K. Carey, K. Rasouli, and W. L. Quinton, “A radiative–conductive–convective approach to calculate thaw season ground surface temperatures for modelling frost table dynamics,” Hydrological Processes, vol. 29, no. 18, pp. 3954–3965, Aug. 2015, doi: 10.1002/hyp.10573.
 */ 
template<typename T>
concept RCC_data = requires(T& t)
{
    // Inputs
	{ t.thaw_front_depth() } -> std::floating_point;

	{ t.air_temperature() } -> std::floating_point; 

	{ t.net_radiation() } -> std::floating_point; 
    
    { t.thaw_front_depth_last() } -> std::floating_point;
    
    // Outputs    
	{ t.surface_temperature(std::declval<const double>()) } -> std::same_as<void>;

    { t.thaw_front_depth_last(std::declval<const double>()) } -> std::same_as<void>;
};

template<RCC_data data>
class RCC_surface_temperature : public base_step<data>
{
public:
	explicit RCC_surface_temperature() {};
	~RCC_surface_temperature() {};

	void execute(data& d) override final;

private:
};


template<RCC_data data>
void RCC_surface_temperature<data>::execute(data& d)
{
    static constexpr auto a = 0.77;
    static constexpr auto b = 0.02;
    static constexpr auto c = 7.0;
    static constexpr auto e = 0.03; // d is taken by input argument
                                      //

    // arctan is often used because it varies smoothly from -pi/2 to pi/2
    // this value makes it from -1 to 1.
    static constexpr auto ARCTAN_NORMALIZER = 2.0 / 3.14156265; 
    d.thaw_front_depth_last(std::max(d.thaw_front_depth_last(),
            d.thaw_front_depth()));
    

	auto T = 
		(a * d.air_temperature() + b * d.net_radiation()) * 
		std::atan(c * (d.thaw_front_depth_last() + e)) * ARCTAN_NORMALIZER;
	d.surface_temperature(T);
};

