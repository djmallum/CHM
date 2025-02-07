#include "PenmanMonteith.hpp"


PenmanMonteith::PenmanMonteith(double& LAI, double& LAImax, double& veg_Ht, double& wind_height, double& stomatal_res_min, double& soil_d, double& F_to_g, const double& Cp, const double& K, const double& tension, const double& pore_sz, const double& theta_pwp, const double& phi) 
    : leaf_area_index(LAI), leaf_area_index_max(LAImax), Veg_height(veg_Ht), wind_measurement_height(wind_height), stomatal_resistance_min(stomatal_res_min), soil_depth(soil_d), Frac_to_ground(F_to_g), heat_capacity_air(Cp), kappa(K), air_entry_tension(tension), pore_size_dist(pore_sz), wilt_point(theta_pwp), porosity(phi)
{
    if (leaf_area_index == 0)
        has_vegetation = false;
    else
        has_vegetation = true;
     
    SPDLOG_DEBUG(air_entry_tension);
    SPDLOG_DEBUG(pore_size_dist);
    SPDLOG_DEBUG(wilt_point);
    SPDLOG_DEBUG(porosity); 

    // TODO other checks on values might be good
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
    //SPDLOG_DEBUG("CalcAeroResistance");
    //SPDLOG_DEBUG(wind_measurement_height);
    //SPDLOG_DEBUG(d);
    //SPDLOG_DEBUG(Z0);
    //SPDLOG_DEBUG(kappa);
    //SPDLOG_DEBUG(var.wind_speed);
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
    SPDLOG_DEBUG("StomatalREsistance");
    //SPDLOG_DEBUG(stomatal_resistance_min);
    //SPDLOG_DEBUG(has_vegetation);
    //SPDLOG_DEBUG(leaf_area_index_max);
    //SPDLOG_DEBUG(var.short_wave_in);
    //SPDLOG_DEBUG(air_entry_tension);
    //SPDLOG_DEBUG(porosity);
    double rcstar = stomatal_resistance_min;

    // In CRHM, the below calculation is an option, for now just use the minimum option.
    if (has_vegetation)
    {
        SPDLOG_DEBUG("has_veg");

        double LAI = Veg_height/2.0*leaf_area_index_max; //TODO ad hoc for test
        SPDLOG_DEBUG(LAI);
        rcstar = stomatal_resistance_min * leaf_area_index_max / LAI;
        SPDLOG_DEBUG(rcstar);
        // rcstar = stomatal_resistance_min * leaf_area_index_max / leaf_area_index; TODO commented for test
    }
    // TODO check units. for example, short_wave_in - 1.5 is suspect
    double f1 = 1.0;
    if (var.short_wave_in > 0.0)
       f1 = std::max(1.0, 500.0/(var.short_wave_in - 1.5));  
//max <double> (1.0, 500.0/(var.short_wave_in - 1.5));  
    
    double f2 = std::max(1.0, 2.0 * (var.saturated_vapour_pressure - var.vapour_pressure) );
//<double> (1.0, 2.0 * (var.saturated_vapour_pressure - var.vapour_pressure) );

    double p = air_entry_tension * pow(porosity / (var.soil_storage/soil_depth/1000.0 + wilt_point), pore_size_dist);  
    SPDLOG_DEBUG("p");
    SPDLOG_DEBUG(p);
    SPDLOG_DEBUG(air_entry_tension);
    SPDLOG_DEBUG(var.soil_storage);
    SPDLOG_DEBUG(soil_depth);
    SPDLOG_DEBUG(wilt_point);
    SPDLOG_DEBUG(pore_size_dist);
    SPDLOG_DEBUG(porosity / (var.soil_storage/soil_depth/1000.0 + wilt_point));
    double f3 = std::max(1.0, p/40.0);

    double f4 = 1.0;
    if (var.t < 5.0 || var.t > 40.0)
        f4 = 5000/50;
    
    SPDLOG_DEBUG(f1);
    SPDLOG_DEBUG(f2);
    SPDLOG_DEBUG(f3);
    SPDLOG_DEBUG(f4);

    if (var.short_wave_in <= 0)
        return 5000;
    else
    {
        return std::min(rcstar * f1 * f2 * f3 * f4, 5000.0);
    }

}

void PenmanMonteith::CalcEvapT(var_base& basevar, model_output& output)
{
    SPDLOG_DEBUG("Calc");
    SPDLOG_DEBUG(air_entry_tension);
    SPDLOG_DEBUG(pore_size_dist);
    SPDLOG_DEBUG(wilt_point);
    SPDLOG_DEBUG(porosity); 

    const PM_vars & var = static_cast<const PM_vars&>(basevar);
    
    double Q =  var.all_wave_net * (1 - Frac_to_ground);

    if (IsFirstRun)
    {
        CalcHeights();
        IsFirstRun = false;
    }

    double aero_resistance = CalcAeroResistance(var);
    double stomatal_resistance = CalcStomatalResistance(var);
    
    SPDLOG_DEBUG("final calc"); 
    SPDLOG_DEBUG(var.t);
    SPDLOG_DEBUG(delta(var.t));
    SPDLOG_DEBUG(Q);
    SPDLOG_DEBUG(var.vapour_pressure);
    SPDLOG_DEBUG(var.P_atm);
    SPDLOG_DEBUG(lambda(var.t));
    SPDLOG_DEBUG(var.saturated_vapour_pressure);
    SPDLOG_DEBUG(gamma(var.P_atm, var.t));
    SPDLOG_DEBUG(stomatal_resistance);
    SPDLOG_DEBUG(aero_resistance);
    output.ET = ( delta(var.t) * Q + AirDensity(var.t,var.vapour_pressure,var.P_atm) * heat_capacity_air / (lambda(var.t)*1e3) * ( var.saturated_vapour_pressure - var.vapour_pressure )/ aero_resistance )
       / ( delta(var.t) + gamma(var.P_atm, var.t) * (1 + stomatal_resistance / aero_resistance ) );
    SPDLOG_DEBUG("Output");
    SPDLOG_DEBUG(output.ET); 
}

// TODO make delta, gamma, density fucntions

double PenmanMonteith::AirDensity(double& t, double& ea, double& Pa) // atmospheric density (kg/m^3)
{
  const double R0 = 2870;
   return (1E4*Pa /(R0*( 273.15 + t))*(1.0 - 0.379*(ea/Pa))); //
}

