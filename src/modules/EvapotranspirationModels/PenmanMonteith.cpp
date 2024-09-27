#include "PenmanMonteith.hpp"


PenmanMonteith::PenmanMonteith(const double& Cp, const double& K, const double& tension, const double& pore_sz) 
    : heat_capacity_air(Cp), kappa(K), air_entry_tension(tension), pore_size(pore_sz)
{
    
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
     return pow( log((wind_measurement_height - d)/Z0),2) / 
        (pow(kappa,2) * var.wind_speed);
}

double PenmanMonteith::CalcStomatalResistance(const PM_vars& var)
{
    double rcstar = stomatal_resistance_min;

    double leaf_area_index = Veg_height/Veg_height_max * 
        (LAImin + seasonal_growth * (LAImax - LAImin));
    rcstar = stomatal_resistance_min * LAImax / leaf_area_index;

    double f1 = 1.0;
    if (var.short_wave_in > 0.0)
       f1 = std::max(1.0, 500.0/(var.short_wave_in - 1.5));  
//max <double> (1.0, 500.0/(var.short_wave_in - 1.5));  
    
    double f2 = std::max(1.0, 2.0 * (var.saturated_vapour_pressure - var.vapour_pressure) );
//<double> (1.0, 2.0 * (var.saturated_vapour_pressure - var.vapour_pressure) );

    double volumetric_soil_storage = (var.soil_storage/soil_depth + 
           SetSoilProperties[1]) / SetSoilProperties[3]; // TODO Name is wrong, see equation 6-12 in Dingman, I should derive this myself so I fully understand it. And then make list of things to include 
  
    double p = air_entry_tension * pow(1.0 / volumetric_soil_storage, pore_size);  // TODO soilproperties is disconnected, will add it as a member function to triangulation, and likely handle it in the main evap class
    
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

