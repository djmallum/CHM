#pragma once

#include "base_step.hpp"
#include "double_range.hpp"
#include <concepts>
#include <utility>

/*
 * Source:
 * 
 * Equation 23 of:
 *
 * Granger, R.J. and Gray, D.M. 1990. A net radiation model for calculating daily snowmelt in open environments. Nordic Hydrology 21: 217-234.
 *
 * Computes the best-fit, linear regression of the net radiation from just the incoming short-wave
 */ 

template<typename T>
concept NetAllData = requires(T& t)
{
    // Inputs
    { t.incoming_short_wave() } -> std::floating_point;

    { t.albedo() } -> std::floating_point;


    // Outputs
    { t.net_all_wave(std::declval<const double>()) } -> std::same_as<void>;


};

template<NetAllData data>
class net_all_bad_lake : public base_step<data>
{
public:
    explicit net_all_bad_lake() {};
    ~net_all_bad_lake() {};

    void execute(data& d) override final;

private:
    static inline constexpr double a = -2.24;
    static inline constexpr double b = 0.651;
};

template<NetAllData data>
void net_all_bad_lake<data>::execute(data& d)
{
    double net_all_wave;
    double_range::Unit albedo = d.albedo();

    net_all_wave = a + b*d.incoming_short_wave() * ( 1 - albedo );

    d.net_all_wave(net_all_wave);
};

