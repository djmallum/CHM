#include "PenmanMonteith.hpp"


PenmanMonteith::PenmanMonteith(double& LAI, double& LAImax, const double& Cp, const double& K, const double& tension, const double& pore_sz, const double& theta_pwp, const double& phi) 
    : leaf_area_index(LAI), leaf_area_index_max(LAImax), heat_capacity(Cp), kappa(K), air_entry_tension(tension), pore_size_dist(pore_sz), wilt_point(theta_pwp), porosity(phi)
{
    if (leaf_area_index == 0)
        has_vegetation = false;
    else
        has_vegetation = true;
}

PenmanMonteith::~PenmanMonteith()
{
    //Do nothing
}

void PenmanMonteith::CalcHeights()
{
    Z0 = Veg_height/7.6;
    d = Veg_height*0.67;
}

double PenmanMonteith::CalcAeroResistance(const PM_vars& var)
{
    if (wind_measurement_height - d > 0)
    {
        return pow( log((wind_measurement_height - d)/Z0),2) / 
            (pow(kappa,2) * var.wind_speed);
    }
    else
    {
        return 0; // I don't know if this is right, but it is at least... safe.
    }
}

double PenmanMonteith::CalcStomatalResistance(const PM_vars& var)
{
    double rcstar = stomatal_resistance_min;

    // In CRHM, the below calculation is an option, for now just use the minimum option.
    if (has_vegetation)
    {
        rcstar = stomatal_resistance_min * LAImax / leaf_area_index;
    }
    // TODO check units. for example, short_wave_in - 1.5 is suspect
    double f1 = 1.0;
    if (var.short_wave_in > 0.0)
       f1 = std::max(1.0, 500.0/(var.short_wave_in - 1.5));  
//max <double> (1.0, 500.0/(var.short_wave_in - 1.5));  
    
    double f2 = std::max(1.0, 2.0 * (var.saturated_vapour_pressure - var.vapour_pressure) );
//<double> (1.0, 2.0 * (var.saturated_vapour_pressure - var.vapour_pressure) );

    double p = air_entry_tension * pow(porosity / (var.soil_storage/soil_depth + wilt_point), pore_size_dist);  
    
    double f3 = std::max(1.0, p/40.0);

    double f4 = 1.0;
    if (var.t < 5.0 || var.t > 40.0)
        f4 = 5000/50;

    if (var.short_wave_in <= 0)
        return 5000;
    else
    {
        return std::min(rcstar * f1 * f2 * f3 * f4, 5000.0);
    }

}

void PenmanMonteith::CalcEvapT(var_base& basevar, model_output& output)
{
    const PM_vars & var = static_cast<const PM_vars&>(basevar);
    
    double Q =  var.all_wave_net * (1 - Frac_to_ground);

    if (IsFirstRun)
    {
        CalcHeights();
        IsFirstRun = false;
    }

    double aero_resistance = CalcAeroResistance(var);
    double stomatal_resistance = CalcStomatalResistance(var);

    output.ET = ( delta(var.t) * Q + AirDensity(var.t,var.vapour_pressure,var.P_atm) * heat_capacity_air / (lambda(var.t)*1e3) * ( var.saturation_vapour_pressure - var.vapour_pressure )/ aero_resistance )
       / ( delta(var.t) + gamma(var.P_atm, var.t) * (1 + stomatal_resistance / aero_resistance ) ); 
}

// TODO make delta, gamma, density fucntions

double PenmanMonteith::AirDensity(double& t, double& ea, double& Pa) // atmospheric density (kg/m^3)
{
  const double R0 = 2870;
   return (1E4*Pa /(R0*( 273.15 + t))*(1.0 - 0.379*(ea/Pa))); //
}

