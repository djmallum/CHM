//
// Created by Allum, Donovan on 2026-09-03.
//
export module penman_monteith;

export struct Input
{
    double temperature;
    double vapour_pressure;
    double saturated_vapour_pressure;
    double air_pressure;
    double net_radiation;
};
export class Params
{
public:
    double heat_capacity_air;
    double ground_flux;
    size_t seconds_per_step;
    bool first_run() const;
    void set_first_run(bool cond) const;
    void heights() const;
private:
    // TODO smells bad, change so that it is just set on construction
    mutable bool _first_run = true;
    mutable double _veg_height = 0.0;
};

class Components
{
public:
    Components(const Input*,const Params*);
    double aero_resistance;
    double stomatal_resistance;
    double radiation;
    double mass;
    double gamma;
    double delta;
    double lambda;
    double s_per_time_step() const;
private:
    const Input* const _input;
    const Params* const _parameters;
};
export struct Output
{
    explicit Output(const Components&);
    Output() = delete;
    Output(double,double) = delete;

    double evapotranspiration;
    double stomatal_resistance;
};
double delta                (const Input&);
double gamma                (const Input&,const Params&);
double air_density          (const Input&);
double get_aero_resistance      (const Input&);
double get_stomatal_resistance  (const Input&);

export Output calc_evapotranspiration(const Params& parameters, const Input& input)
{
    if (parameters.first_run())
    {
        // TODO put heights inside set_first_run
        parameters.heights();
        parameters.set_first_run(false);
    }

    const auto components = Components(&input,&parameters);

    const auto output = Output(components);
    return output;

};

double Components::s_per_time_step() const
{
    return _parameters->seconds_per_step;
}
Output::Output(const Components& components)
{
    static constexpr auto WATER_DENSITY = 1000.0;
    static constexpr auto MM_PER_M = 1000.0;
    evapotranspiration = (components.radiation + components.mass) /
            ( components.delta + components.gamma
                * ( 1 + components.stomatal_resistance / components.aero_resistance ));
    // Units are W/m^2

    // TODO PM and PT methods should both stop at W/m^2 and conversion done in module

    evapotranspiration *= 1.0 / (WATER_DENSITY * components.lambda); // Converts units to m/s

    evapotranspiration *= MM_PER_M * components.s_per_time_step(); // Converts from m/s to mm/(time step)

    stomatal_resistance = components.stomatal_resistance;
}

Components::Components(const Input* input,const Params* parameters) : _input(input), _parameters(parameters)
{
    aero_resistance = get_aero_resistance(*_input);
    stomatal_resistance = get_stomatal_resistance(*_input);
    this->delta = delta(*_input);
    radiation = delta(*_input) * _input->net_radiation * (1 - _parameters->ground_flux);  //Units: W/m^2 * kPa/K

    // Units: W/m^2 * kPa/K
    mass = air_density(*_input) * _parameters->heat_capacity_air *
    (_input->saturated_vapour_pressure - _input->vapour_pressure)/ aero_resistance;
}
