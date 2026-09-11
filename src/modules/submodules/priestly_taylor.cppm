//
// Created by Allum, Donovan on 2026-09-07.
//
export module priestly_taylor;
import physics_functions;
struct Output {
    double evapotranspiration;
};
struct Input {
    WattsPerM2 net_radiation;
    KiloPascals air_pressure;
    Celsius air_temperature;

};

struct Params {
    Percent ground_flux;
    double correction_factor;
    size_t seconds_per_step;
    JoulesPerKgPerKelvin heat_capacity_air;
};

Output calc_evapotranspiration(const Input& input, const Params& parameters)
{
    static constexpr double WATER_DENSITY = 1000.0;
    static constexpr double MM_PER_M = 1000.0;
    const double Q = input.net_radiation * (Percent{1.0} - parameters.ground_flux);
    Output output;
    const auto delta = saturation_vapour_pressure_slope(input.air_temperature);
    output.evapotranspiration = parameters.correction_factor * delta * Q / (delta + psychrometric_constant(input.air_temperature,input.air_pressure,parameters.heat_capacity_air) );
    // TODO PM and PT methods should both stop at W/m^2 and conversion done in module
    // Justification for unit conversions in PenmanMonteith
    output.evapotranspiration *= 1.0 / (WATER_DENSITY * (input.air_temperature));

    output.evapotranspiration *= MM_PER_M * parameters.seconds_per_step;
    return output;
}
