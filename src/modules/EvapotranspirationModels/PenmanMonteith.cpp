#include "PenmanMonteith.hpp"


PenmanMonteith::PenmanMonteith(const double& Cp, const double& K, const double& tension, const double& pore_sz, const double& theta_pwp, const double& phi) 
    : heat_capacity(Cp), kappa(K), air_entry_tension(tension), pore_size(pore_sz), wilt_point(theta_pwp), porosity(phi)
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

double PenmanMonteith::CalcAeroResistance(PM_var& var)
{
     return pow( log((wind_measurement_height - d)/Z0),2) / 
        (pow(kappa,2) * var.wind_speed);
}

double PenmanMonteith::CalcStomatalResistance(PM_var& var)
{
    double rcstar = stomatal_resistance_min;

    double leaf_area_index = Veg_height/Veg_height_max * 
        (LAImin + seasonal_growth * (LAImax - LAImin));
    rcstar = stomatal_resistance_min * LAImax / leaf_area_index;

    double f1 = 1.0;
    if (var.ShortWave_in > 0.0)
       f1 = max <double> (1.0, 500.0/(var.ShortWave_in - 1.5));  
    
    double f2 = max <double> (1.0, 2.0 * (var.saturated_vapour_pressure - var.vapour_pressure) );

    double volumetric_soil_storage = (var.soil_storage/soil_depth + 
           wilt_point) / porosity; // TODO SetSoilProperties is currently disconnected. Chris says to add it to a single file and share it to relevant modules.
                                                       
  
    double p = air_entry_tension * pow(1.0 / volumetric_soil_storage, pore_size_dist);  // TODO soilproperties is disconnected, will add it as a member function to triangulation, and likely handle it in the main evap class
    
    double f3 = max <double> (1.0, p/40.0);

    double f4 = 1.0;
    if (var.t < 5.0 || var.t > 40.0)
        f4 = 5000/50;

    if (var.ShortWave_in <= 0)
        return 5000;
    else
    {
        return min <double> (rcstar * f1 * f2 * f3 * f4, 5000.0);
    }

}

void PenmanMonteith::CalcEvapT(var_base& basevar, model_output& output)
{
    const PM_var & var = static_cast<const PM_var&>(basevar);
    
    double Q =  var.Rnet * (1 - Frac_to_ground);

    if (IsFirstRun)
    {
        CalcHeights();
        IsFirstRun = false;
    }

    double aero_resistance = CalcAeroResistance(var);
    double stomatal_resistance = CalcStomatalResistance(var);

    output.ET = ( delta(var.t) * Q + AirDensity(var.t,var.vapour_pressure,var.P_atm) * heat_capacity / (lambda(var.t)*1e3) * ( var.saturation_vapour_pressure - var.vapour_pressure )/ aero_resistance )
       / ( delta(var.t) + gamma(var.P_atm, var.t) * (1 + stomatal_resistance / aero_resistance ) ); 
}

// TODO make delta, gamma, density fucntions

double PenmanMonteith::AirDensity(double& t, double& ea, double& Pa) // atmospheric density (kg/m^3)
{
  const double R0 = 2870;
   return (1E4*Pa /(R0*( 273.15 + t))*(1.0 - 0.379*(ea/Pa))); //
}
