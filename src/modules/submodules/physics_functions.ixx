//
// Created by Allum, Donovan on 2026-09-07.
//

export module physics_functions;
#include <cmath>

struct Units {
    double value{};
};
struct Celsius : public Units{};
struct KiloPascals : public Units{};
struct JoulesPerKgPerKelvin : public Units{};
struct KpAPerDegree : public Units {};
struct KgPerM3 : public Units {};
struct JoulePerKg : public Units {};

double delta(const Celsius air_temperature);
double air_density();
double gamma(const Celsius air_temperature, const KiloPascals air_pressure, const JoulesPerKgKelvin heat_capacity_air);
double lambda(const Celsius air_temperature);

// Implementations
KpAPerDegree delta(const Celsius air_temperature)
{
    auto result = KpAPerDegree{};
    if (const double t = air_temperature.value; t > 0.0)
        result.value = (2504.0*exp(17.27 * t/(t+237.3)) / pow(t+237.3,2));
    else
        result.value = (3549.0*exp( 21.88 * t/(t+265.5)) / pow(t+265.5,2));

    return result;
}

KgPerM3 air_density()
{
    static_assert(false,"Function not implemented");
    return KgPerM3{};
    // TODO
}

KpAPerDegree gamma(const Celsius air_temperature, const KiloPascals air_pressure, const JoulesPerKgKelvin heat_capacity_air )
{
    auto result = KpAPerDegree{};
    result.value = heat_capacity_air.value * air_pressure.value / (0.622 * lambda(air_temperature)); // lambda (J/kg)
    return result;
}

JoulePerKg lambda(const AirTemperature air_temperature)
{
    const double t = air_temperature.value;
    static constexpr double J_PER_MJ = 1e6;
    auto result = JoulePerKg{};
    result.value = (2.501 - 0.002361 * t) * J_PER_MJ;
    return result; // original is MegaJoules/kg, 1e6 returns it to joules/kg
}
