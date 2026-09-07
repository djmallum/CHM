//
// Created by Allum, Donovan on 2026-09-07.
//
import module physics_functions;
export module priestly_taylor;
struct Output {
    double evapotranspiration;
};
struct Input {
    double net_radiation;
    double air_pressure;
    double air_temperature;

};
struct Params {
    double percent_ground_flux;
    double correction_factor;
    size_t seconds_per_step;
    double heat_capacity_air;
};

Output calc_evapotranspiration(const Input& input, const Params& parameters)
{
    static constexpr double WATER_DENSITY = 1000.0;
    static constexpr double MM_PER_M = 1000.0;
    double Q = input.net_radiation * (1 - parameters.percent_ground_flux);
    Output output;
    output.evapotranspiration = parameters.correction_factor * delta(input.air_temperature) * Q / (delta(input.air_temperature) + gamma(input.air_pressure,input.air_temperature,parameters.heat_capacity_air) );
    // TODO PM and PT methods should both stop at W/m^2 and conversion done in module
    // Justification for unit conversions in PenmanMonteith
    output.evapotranspiration *= 1.0 / (WATER_DENSITY * lambda(input.air_temperature));

    output.evapotranspiration *= MM_PER_M * parameters.seconds_per_step;
    return output;
}
