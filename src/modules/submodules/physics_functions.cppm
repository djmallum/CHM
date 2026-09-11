//
// Created by Allum, Donovan on 2026-09-07.
//

export module physics_functions;
#include <cmath>
#include <compare>

/*
 **NewType** A design pattern where instead of passing around raw doubles, you pass around a struct with a single double
 member.

 The `Units` struct is the base struct for the NewType pattern implemention. It includes a + and - operator, returning
 double.

 It is not intended to be used for unit conversion, just for strong, explicit units enforced by type safety. The
 compiler does a little bit of work by ensuring that the right type is passed. However, this still requires the
 programmer to do the work themselves. Nothing is stopping them from accidentally passing something that has units of
 Metres and initializing a MiliMetre type. But what this does is allow readers to know what the units are supposed to
 be.
 */
export struct Units {
    double value{};

    constexpr double operator+(const Units & other) const {
        return value + other.value;
    }
    constexpr double operator-(const Units & other) const
    {
        return value - other.value;
    }
    constexpr double operator*(const Units & other) const
    {
        return value * other.value;
    }
    constexpr double operator/(const Units & other) const
    {
        return value / other.value;
    }
    constexpr double operator*(const double & other) const
    {
        return value * other;
    }
    constexpr double operator/(const double & other) const
    {
        return value / other;
    }
    constexpr bool operator==(const Units& other) const
    {
        return value == other.value;
    };
};
export double operator*(const double lhs, const Units& rhs)
{
    return lhs * rhs.value;
};
export double operator/(const double lhs, const Units& rhs)
{
    return lhs / rhs.value;
};
export constexpr std::partial_ordering operator<=>(const Units& lhs, const Units& rhs)
{
    const auto result = lhs.value <=> rhs.value;
    return result;
}

static_assert(Units{1.0} - Units{0.75} == 0.25);
static_assert(Units{0.75} - Units{1.0} == -0.25);

export struct Celsius : Units{};
export struct KiloPascals : Units{};
export struct JoulesPerKgPerKelvin : Units{};
/*
 Kilopascals per Kelvin
 */
export struct KpAPerKelvin : Units {};
/*
 Kg/m^3
 */
export struct KgPerM3 : Units {};
export struct JoulePerKg : Units {};
/*
 W/m^2
 */
export struct WattsPerM2 : Units {};
export struct MetresPerSecond : Units {};
export struct SecondsPerMetre : Units {};
export struct Metre : Units {};
export struct MiliMetre : Units {};
export struct Percent : Units
{
    explicit Percent(const double v)
    {
        value = v;
        // TODO error if v<0 or v>1
    }
};

export KpAPerKelvin saturation_vapour_pressure_slope(Celsius air_temperature);
export KgPerM3 air_density(Celsius air_temperature,KiloPascals air_pressure,KiloPascals vapour_pressure);
export KpAPerKelvin psychrometric_constant(Celsius air_temperature,KiloPascals air_pressure,JoulesPerKgPerKelvin heat_capacity_air);
export JoulePerKg latent_heat_vapour(Celsius air_temperature);

// Implementations
KpAPerKelvin saturation_vapour_pressure_slope(const Celsius air_temperature)
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

KpAPerKelvin psychrometric_constant(const Celsius air_temperature, const KiloPascals air_pressure, const JoulesPerKgPerKelvin heat_capacity_air )
{
    auto result = KpAPerKelvin{};
    result.value = heat_capacity_air.value * air_pressure.value / (0.622 * latent_heat_vapour(air_temperature)); // lambda (J/kg)
    return result;
}

JoulePerKg latent_heat_vapour(const Celsius air_temperature)
{
    const double t = air_temperature.value;
    static constexpr double J_PER_MJ = 1e6;
    auto result = JoulePerKg{};
    result.value = (2.501 - 0.002361 * t) * J_PER_MJ;
    return result; // original is MegaJoules/kg, 1e6 returns it to joules/kg
}
