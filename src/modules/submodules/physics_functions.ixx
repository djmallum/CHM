//
// Created by Allum, Donovan on 2026-09-07.
//

export module physics_functions;
#include <cmath>

struct Units {
    double value{};

    constexpr double operator+(const Units & other) const {
        return value + other.value;
    }
    constexpr double operator-(const Units & other) const
    {
        return value - other.value;
    }
    constexpr bool operator==(const Units& other) const
    {
        return value == other.value;
    };
};

static_assert(Units{1.0} - Units{0.75} == 0.25);
static_assert(Units{0.75} - Units{1.0} == -0.25);

export struct Celsius : Units{};
export struct KiloPascals : Units{};
export struct JoulesPerKgPerKelvin : Units{};
export struct KpAPerKelvin : Units {};
export struct KgPerM3 : Units {};
export struct JoulePerKg : Units {};
export struct WattsPerM2 : Units {};
export struct MetresPerSecond : Units {};
export struct SecondsPerMetre : Units {};
export struct Metre : Units {};
export struct MiliMetre : Units {};

export KpAPerKelvin delta(Celsius air_temperature);
export KgPerM3 air_density(Celsius air_temperature,KiloPascals air_pressure,KiloPascals vapour_pressure);
export KpAPerKelvin gamma(Celsius air_temperature,KiloPascals air_pressure,JoulesPerKgPerKelvin heat_capacity_air);
export JoulePerKg lambda(Celsius air_temperature);

// Implementations
KpAPerKelvin delta(const Celsius air_temperature)
{
    auto result = KpAPerKelvin{};
    if (const double t = air_temperature.value; t > 0.0)
        result.value = (2504.0*exp(17.27 * t/(t+237.3)) / pow(t+237.3,2));
    else
        result.value = (3549.0*exp( 21.88 * t/(t+265.5)) / pow(t+265.5,2));

    return result;
}

KgPerM3 air_density(const Celsius air_temperature,const KiloPascals air_pressure, const KiloPascals vapour_pressure)
{
    constexpr double R0 = 2870;
    const double result =  (1E4*air_pressure.value /(R0*( 273.15 + air_temperature.value))*(1.0 - 0.379*(vapour_pressure.value/air_pressure.value))); //
    return KgPerM3{result};
}

KpAPerKelvin gamma(const Celsius air_temperature, const KiloPascals air_pressure, const JoulesPerKgPerKelvin heat_capacity_air )
{
    auto result = KpAPerKelvin{};
    result.value = heat_capacity_air.value * air_pressure.value / (0.622 * lambda(air_temperature)); // lambda (J/kg)
    return result;
}

JoulePerKg lambda(const Celsius air_temperature)
{
    const double t = air_temperature.value;
    static constexpr double J_PER_MJ = 1e6;
    auto result = JoulePerKg{};
    result.value = (2.501 - 0.002361 * t) * J_PER_MJ;
    return result; // original is MegaJoules/kg, 1e6 returns it to joules/kg
}
