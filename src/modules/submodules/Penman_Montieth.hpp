#include "base_step.hpp"
#include <algorithm>
#include <concepts>

template<typename T>
concept Penman_data = requires(const T& t)
{
    // Inputs
    { t.wind_measurement_height() } -> std::floating_point;

    { t.d() } -> std::floating_point;

    { t.Z0() } -> std::floating_point;

    { t.kappa() } -> std::floating_point;

    { t.wind_speed() } -> std::floating_point;

    { t.stomatal_resistance_min() } -> std::floating_point;

    { t.has_vegetation() } -> std::same_as<bool>;

    { t.Veg_height() } -> std::floating_point;

    { t.leaf_area_index_max() } -> std::floating_point;

    { t.short_wave_in() } -> std::floating_point;

    { t.saturated_vapour_pressure() } -> std::floating_point;

    { t.vapour_pressure() } -> std::floating_point;

    { t.air_entry_tension() } -> std::floating_point;

    { t.porosity() } -> std::floating_point;

    { t.volumetric_moisture_content() } -> std::floating_point;

    { t.pore_size_dist() } -> std::floating_point;

    { t.air_temperature() } -> std::floating_point;

    { t.delta() } -> std::floating_point;

    { t.Q_net() } -> std::floating_point;
    
    { t.Q_g() } -> std::floating_point;

    { t.air_density() } -> std::floating_point;

    { t.heat_capacity_air() } -> std::floating_point;

    { t.gamma() } -> std::floating_point;

    { t.stomatal_resistance() } -> std::floating_point;

    { t.lambda() } -> std::floating_point;

    { t.s_per_time_step() } -> std::integral;


    // Outputs
    { t.stomatal_resistance(std::declval<const double>()) } -> std::same_as<void>;

    { t.ET(std::declval<const double>()) } -> std::same_as<void>;


};



template<Penman_data data>
class Penman_monteith : public base_step<data>
{
public:
    explicit Penman_monteith() {};
    ~Penman_monteith() {};

    void execute(data& d) final;

    class aerodynamic_resistance_calculator;
    class stomatal_resistance_jarvis;
private:
    aerodynamic_resistance_calculator aero_resistance;
    stomatal_resistance_jarvis stomatal_resistance;
};

template<Penman_data data>
class Penman_monteith<data>::aerodynamic_resistance_calculator {
public:
    double calculate(data& d) const;
};

template<Penman_data data>
double Penman_monteith<data>::aerodynamic_resistance_calculator::calculate(data& d) const 
{
    if (d.wind_measurement_height() - d.d() > 0) {
        return pow(log((d.wind_measurement_height() - d.d()) / d.Z0()), 2) / 
               (pow(d.kappa(), 2) * d.wind_speed());
    } else {
        return 0; // I don't know if this is right, but it is at least... safe.
    }
};

template<Penman_data data>
class Penman_monteith<data>::stomatal_resistance_jarvis
{
private:
    static constexpr double UPPER_LIMIT = 5000.0;

public:
    // Interface for the data required by the calculation

    void calculate(const data& d) const; 

private:
    double CalculateF1(const data& d) const;

    double CalculateF2(const data& d) const;

    double CalculateF3(const data& d) const;

    double CalculateF4(const data& d) const;
};

template<Penman_data data>
void Penman_monteith<data>::stomatal_resistance_jarvis::calculate(const data& d) const 
{
    double rcstar = d.stomatal_resistance_min();

    // In CRHM, the below calculation is an option, for now just use the minimum option.
    if (d.has_vegetation()) {
        double LAI = d.Veg_height() / 2.0 * d.leaf_area_index_max(); // TODO ad hoc for test
        rcstar *= d.leaf_area_index_max() / LAI;
        // rcstar = stomatal_resistance_min * leaf_area_index_max / leaf_area_index; TODO commented for test
    }
    
    double f1 = CalculateF1(d);
    double f2 = CalculateF2(d);
    double f3 = CalculateF3(d);
    double f4 = CalculateF4(d);
    
    if (d.short_wave_in() <= 0) {
        d.stomatal_resistance(UPPER_LIMIT);
    } else {
        d.stomatal_resistance(std::min(rcstar * f1 * f2 * f3 * f4, UPPER_LIMIT));
    }

};

template<Penman_data data>
double Penman_monteith<data>::stomatal_resistance_jarvis::CalculateF1(const data& d) const 
{
    if (d.short_wave_in() > 0.0) {
        return std::max(1.0, 500.0 / d.short_wave_in() - 1.5);
    }
    return 1.0;
};

template<Penman_data data>
double Penman_monteith<data>::stomatal_resistance_jarvis::CalculateF2(const data& d) const 
{
    return std::max(1.0, 2.0 * (d.saturated_vapour_pressure() - d.vapour_pressure()));
};

template<Penman_data data>
double Penman_monteith<data>::stomatal_resistance_jarvis::CalculateF3(const data& d) const 
{
    double p = d.air_entry_tension() * pow(d.porosity() / d.volumetric_moisture_content(), d.pore_size_dist());
    return std::max(1.0, p / 40.0);
};

template<Penman_data data>
double Penman_monteith<data>::stomatal_resistance_jarvis::CalculateF4(const data& d) const {
    double t = d.air_temperature();
    if (t < 5.0 || t > 40.0) {
        return UPPER_LIMIT / d.stomatal_resistance_min();
    }
    return 1.0;
}

template<Penman_data data>
void Penman_monteith<data>::execute(data& d)
{
	constexpr double MM_PER_M = 1000.0; // mm/m
    constexpr double WATER_DENSITY = 1000.0; // kg/m^3
    
	double r_a = aero_resistance.calculate(d);
	stomatal_resistance.calculate(d);

	double radiation = d.delta() * (d.Q_net() - d.Q_g());  //Units: kPa/K * W/m^2 (in order, left to right)
	
    double mass = d.air_density() * d.heat_capacity_air() * 
        (d.saturated_vapour_pressure() - d.vapour_pressure())/ r_a;

	double ET = (radiation + mass) / 
		( d.delta() + d.gamma() * ( 1 + d.stomatal_resistance() / r_a));
		// Units are W/m^2

	ET *= 1.0 / (WATER_DENSITY * d.lambda()); // Converts units to m/s 
	
	ET *= MM_PER_M * d.s_per_time_step(); // Converts from m/s to mm/(time step)

    d.ET(ET);

};
