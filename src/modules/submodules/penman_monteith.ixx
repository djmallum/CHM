export module penman_monteith;
import physics_functions;
#include <cmath>

export struct GridParams
{
    Metre vegetation_height;
    Metre wind_measurement_height;
};
export struct Input
{
    Celsius air_temperature;
    KiloPascals vapour_pressure;
    KiloPascals saturated_vapour_pressure;
    KiloPascals air_pressure;
    WattsPerM2 net_radiation;
    MetresPerSecond wind_speed;
    const GridParams* const params;
    Metre z_0() const {return Metre{params->vegetation_height.value/7.6};}
    Metre d() const {return Metre{params->vegetation_height.value*0.67};}
    Metre wind_measurement_height() const { return params->wind_measurement_height;}
    Input() = delete;
};
export class Params
{
public:
    KpAPerKelvin heat_capacity_air;
    WattsPerM2 ground_flux;
    size_t seconds_per_step;
    double kappa;
    Params() = delete;
};

/*
    Numerator struct because units of `Components::radiation` and `Components::mass` are W*kPa / (m^2*K) and would
    be ugly to name a struct for this, instead Numerator is used.
 */
struct Numerator : Units {};

class Components
{
public:
    Components(const Input*,const Params*);
    SecondsPerMetre aero_resistance;
    SecondsPerMetre stomatal_resistance;
    Numerator radiation;
    Numerator mass;
    KpAPerKelvin gamma_;
    KpAPerKelvin delta_;
    JoulePerKg lambda_;
    size_t s_per_time_step() const;
private:
    const Input* const _input;
    const Params* const _parameters;
};
export struct Output
{
    void _compute_ET(const Components& components);
    explicit Output(const Components&);
    Output() = delete;
    Output(double,double) = delete;

    MiliMetre evapotranspiration; // TODO make unitful
    SecondsPerMetre stomatal_resistance;
};
SecondsPerMetre get_aero_resistance      (const Input&);
SecondsPerMetre get_stomatal_resistance  (const Input&);
export Output calc_evapotranspiration(const Params& parameters, const Input& input)
{

    const auto components = Components(&input,&parameters);

    const auto output = Output(components);
    return output;

};

size_t Components::s_per_time_step() const
{
    return _parameters->seconds_per_step;
}
void Output::_compute_ET(const Components& components)
{
    static constexpr auto WATER_DENSITY = 1000.0;
    static constexpr auto MM_PER_M = 1000.0;
    double ET = (components.radiation + components.mass) /
                (components.delta_ + components.gamma_.value *
                                         (1 + components.stomatal_resistance.value / components.aero_resistance.value));

    // TODO PM and PT methods should both stop at W/m^2 and conversion done in module

    ET *= 1.0 / (WATER_DENSITY * components.lambda_.value); // Converts units to m/s

    ET *= MM_PER_M * components.s_per_time_step(); // Converts from m/s to mm/(time step)

    evapotranspiration = MiliMetre{ET};
}
Output::Output(const Components& components)
{
    _compute_ET(components);

    stomatal_resistance = components.stomatal_resistance;
}
SecondsPerMetre get_aero_resistance(const Input& input, const Params& parameters)
{
    if (input.wind_measurement_height() - input.d() > 0)
    {
        return SecondsPerMetre{std::pow( std::log((input.wind_measurement_height() - input.d())/input.z_0().value),2) /
            (std::pow(parameters.kappa,2) * input.wind_speed.value)};
    }
    else
    {
        return SecondsPerMetre{0.0}; // I don't know if this is right, but it is at least... safe.
    }
}

Components::Components(const Input* input,const Params* parameters) : _input(input), _parameters(parameters)
{
    aero_resistance = get_aero_resistance(*_input,*_parameters);
    stomatal_resistance = get_stomatal_resistance(*_input);
    this->delta_ = delta(_input->air_temperature);
    this->lambda_ = lambda(_input->air_temperature);
    radiation = Numerator{this->delta_.value * _input->net_radiation.value * (1 - _parameters->ground_flux.value)};  //Units: W/m^2 * kPa/K

    // Units: W/m^2 * kPa/K
    mass = Numerator{air_density(_input->air_temperature,_input->air_pressure,_input->vapour_pressure).value * _parameters->heat_capacity_air.value *
    (_input->saturated_vapour_pressure - _input->vapour_pressure)/ aero_resistance};
}
