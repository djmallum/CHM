#include "PriestleyTaylor.hpp"

PriestleyTaylor::PriestleyTaylor(const double& alpha_const, const double& Cp, const int dt) : alpha(alpha_const), heat_capacity_air(Cp), s_per_time_step(dt)
{

}

PriestleyTaylor::~PriestleyTaylor()
{
    // Do nothing
}

void PriestleyTaylor::CalcEvapT(var_base& basevar, model_output& output)
{
    const PT_vars& var = static_cast<const PT_vars&>(basevar);
    constexpr double water_density = 1000.0;
    constexpr double mm_per_m = 1000.0;
    double Q = var.all_wave_net * (1 - Frac_to_ground);
   
    output.ET = alpha * delta(var.air_temperature) * Q / (delta(var.air_temperature) + gamma(var.P_atm,var.air_temperature,heat_capacity_air) );
    // TODO PM and PT methods should both stop at W/m^2 and conversion done in module
    // Justification for unit conversions in PenmanMonteith
    output.ET *= 1.0 / (water_density * lambda(var.air_temperature));

    output.ET *= mm_per_m * s_per_time_step;

}
