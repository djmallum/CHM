#include "PriestleyTaylor.hpp"

PriestleyTaylor::PriestleyTaylor(const double& alpha_const) : alpha(alpha_const)
{

}

PriestleyTaylor::~PriestleyTaylor()
{
    // Do nothing
}

void PriestleyTaylor::CalcEvapT(var_base& basevar, model_output& output)
{
    const PT_vars& var = static_cast<const PT_vars&>(basevar);

    double Q = var.all_wave_net * (1 - Frac_to_ground);

    output.ET = alpha * delta(var.air_temperature) * Q / (delta(var.air_temperature) + gamma(var.P_atm,var.air_temperature) );
}
