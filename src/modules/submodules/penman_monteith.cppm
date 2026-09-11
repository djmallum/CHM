module;
#include <algorithm>
export module penman_monteith;
import physics_functions;

export struct GridParams
{
};
export struct Input
{
    Celsius air_temperature;
    KiloPascals vapour_pressure;
    KiloPascals saturated_vapour_pressure;
    KiloPascals air_pressure;
    WattsPerM2 net_radiation;
    MetresPerSecond wind_speed;
    WattsPerM2 short_wave_in;
    MiliMetre soil_storage;
    Input() = delete;
};
export struct DomainParameters
{
    KpAPerKelvin heat_capacity_air;
    Percent ground_flux;
    size_t seconds_per_step;
    SecondsPerMetre stomatal_resistance_min;
    Metre wind_measurement_height;
    double kappa;
    /*
     LAI can be greater than one so leaving unitless
     */
    double leaf_area_index_max;
};
export class Params
{
public:
    Metre vegetation_height;
    MiliMetre soil_depth;
    Percent porosity;
    Percent wilt_point;
    // TODO Unitful!
    double air_entry_tension{};
    double pore_size_dist{};
    [[nodiscard]] Metre wind_measurement_height() const { return domain_parameters->wind_measurement_height;};
    [[nodiscard]] KpAPerKelvin heat_capacity_air() const { return domain_parameters->heat_capacity_air;}
    [[nodiscard]] Percent ground_flux() const { return domain_parameters->ground_flux;}
    [[nodiscard]] size_t seconds_per_step() const { return domain_parameters->seconds_per_step;}
    [[nodiscard]] SecondsPerMetre stomatal_resistance_min() const { return domain_parameters->stomatal_resistance_min;}
    [[nodiscard]] double kappa() const { return domain_parameters->kappa;}
    [[nodiscard]] double leaf_area_index_max() const { return domain_parameters->leaf_area_index_max;}
    [[nodiscard]] explicit Params(const DomainParameters* domain_parameters)
        : porosity(), wilt_point(), domain_parameters(domain_parameters)
    {
    }
    [[nodiscard]] Metre z_0() const {return Metre{vegetation_height/7.6};}
    [[nodiscard]] Metre d() const {return Metre{vegetation_height*0.67};}
    [[nodiscard]] bool has_vegetation() const {return vegetation_height>Metre{0.0};}
private:
    const DomainParameters* const domain_parameters;
};

/*
    Numerator struct because units of `Components::radiation` and `Components::mass` are W*kPa / (m^2*K) and would be ugly to name a struct for this, instead Numerator is used. Did not leave it as double so that users can see the type and redirect themselves to this definition.
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
    [[nodiscard]] size_t s_per_time_step() const;
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

    MiliMetre evapotranspiration;
    SecondsPerMetre stomatal_resistance;
};
SecondsPerMetre get_aero_resistance      (const Input&);
SecondsPerMetre get_stomatal_resistance  (const Input&,const Params& parameters);
export Output calc_evapotranspiration(const Params& parameters, const Input& input)
{
    const auto components = Components(&input,&parameters);

    const auto output = Output(components);
    return output;
};

size_t Components::s_per_time_step() const
{
    return _parameters->seconds_per_step();
}
void Output::_compute_ET(const Components& components)
{
    static constexpr auto WATER_DENSITY = 1000.0;
    static constexpr auto MM_PER_M = 1000.0;
    double ET = (components.radiation + components.mass) /
                (components.delta_.value + components.gamma_.value *
                                         (1 + components.stomatal_resistance.value / components.aero_resistance.value));

    // TODO PM and PT methods should both stop at W/m^2 and conversion done in module

    ET *= 1.0 / (WATER_DENSITY * components.lambda_.value); // Converts units to m/s

    ET *= MM_PER_M * static_cast<double>(components.s_per_time_step()); // Converts from m/s to mm/(time step)

    evapotranspiration = MiliMetre{ET};
}
Output::Output(const Components& components)
{
    _compute_ET(components);

    stomatal_resistance = components.stomatal_resistance;
}
SecondsPerMetre get_stomatal_resistance(const Input& input,const Params& parameters)
{
    SecondsPerMetre stomatal_res_min = parameters.stomatal_resistance_min();
    // In CRHM, the below calculation is an option, for now just use the minimum option.
    if (parameters.has_vegetation())
    {

        const double LAI = parameters.vegetation_height/2.0*parameters.leaf_area_index_max(); //TODO ad hoc for test
        stomatal_res_min.value = parameters.stomatal_resistance_min() * parameters.leaf_area_index_max() / LAI;
        // rcstar = stomatal_resistance_min * leaf_area_index_max / leaf_area_index; TODO commented for test
    }
    // TODO check units. for example, short_wave_in - 1.5 is suspect
    double f1 = 1.0;
    if (input.short_wave_in > WattsPerM2{0.0})
        f1 = std::max(1.0, 500.0/input.short_wave_in - 1.5);

    const double f2 = std::max(1.0, 2.0 * (input.saturated_vapour_pressure - input.vapour_pressure) );

    const double p = parameters.air_entry_tension * pow(parameters.porosity / (input.soil_storage/parameters.soil_depth + parameters.wilt_point.value), parameters.pore_size_dist);
    const double f3 = std::max(1.0, p/40.0);

    double f4 = 1.0;
    if (input.air_temperature < Celsius{5.0} || input.air_temperature > Celsius{40.0})
        f4 = 5000.0/parameters.stomatal_resistance_min();

    if (input.net_radiation <= WattsPerM2{0.0})
        return SecondsPerMetre{5000.0};
    else
    {
        return SecondsPerMetre{std::min(parameters.stomatal_resistance_min() * f1 * f2 * f3 * f4, 5000.0)};
    }
}
SecondsPerMetre get_aero_resistance(const Input& input, const Params& parameters)
{
    if (parameters.wind_measurement_height() - parameters.d() > 0)
    {
        return SecondsPerMetre{std::pow( std::log((parameters.wind_measurement_height() - parameters.d())/parameters.z_0().value),2) /
            (std::pow(parameters.kappa(),2) * input.wind_speed.value)};
    }
    else
    {
        return SecondsPerMetre{0.0}; // I don't know if this is right, but it is at least... safe.
    }
}

Components::Components(const Input* input,const Params* parameters) : _input(input), _parameters(parameters)
{
    aero_resistance = get_aero_resistance(*_input,*_parameters);
    stomatal_resistance = get_stomatal_resistance(*_input,*_parameters);
    this->delta_ = saturation_vapour_pressure_slope(_input->air_temperature);
    this->lambda_ = latent_heat_vapour(_input->air_temperature);
    radiation = Numerator{this->delta_.value * _input->net_radiation.value * (1 - _parameters->ground_flux().value)};  //Units: W/m^2 * kPa/K

    // Units: W/m^2 * kPa/K
    mass = Numerator{air_density(_input->air_temperature,_input->air_pressure,_input->vapour_pressure).value * _parameters->heat_capacity_air().value *
    (_input->saturated_vapour_pressure - _input->vapour_pressure)/ aero_resistance};
}
