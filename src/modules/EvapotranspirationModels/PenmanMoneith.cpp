#include "PenmanMonteith.hpp"


PenmanMonteith::PenmanMonteith(PM_param& param)
{
    
    CalcHeights(param);
}

PenmanMonteith::~PenmanMonteith()
{
    //Do nothing
}

void PenmanMonteith::CalcHeights(PM_param& param)
{


    param.Z0 = param.Veg_height/7.6;
    param.d = param.Veg_height*0.67;
}

void PenmanMonteith::CalcAeroResistance(PM_param& param, PM_var& var)
{
    var.aero_resistance = pow( log((param.wind_measurement_height - param.d)/param.Z0),2) / 
        (pow(param.kappa,2) * var.wind_speed);
}

// TODO no such function CalcSaturationVapourPressure, put it here for now
//

double PenmanMonteith::CalcSaturationVapourPressure(PM_var& var)
{
    if (var.t > 0.0)
	{
		return 0.611 * exp(17.27 * var.t / (var.t + 237.3));
	}
	else
	{
		return 0.611 * exp(21.88 * var.t / (var.t + 265.5));
	}
}

void PenmanMonteith::CalcStomatalResistance(PM_param& param, PM_var& var)
{
    double rcstar = param.rc_min;

    double leaf_area_index = param.Veg_height/param.Veg_height_max * 
        (param.LAImin + param.seasonal_growth * (param.LAImax - param.LAImin));
    rcstar = param.rc_min * param.LAImax / leaf_area_index;

    double f1 = 1.0;
    if (var.ShortWave_in > 0.0)
       f1 = max <double> (1.0, 500.0/(var.ShortWave_in - 1.5));  
    
    double f2 = max <double> (1.0, 2.0 * ( CalcSaturationVapourPressure(t) - var.vapour_pressure) );

    double volumetric_soil_storage = (var.soil_storage/param.soil_depth + 
           SetSoilProperties[1]) / SetSoilProperties[3]; // TODO SetSoilProperties is currently disconnected. Chris says to add it to a single file and share it to relevant modules.
                                                       
  
    double p = soilproperties[AIRENT] * pow(1.0 / volumetric_soil_storage, soilproperties[PORESZ]);  // TODO soilproperties is disconnected, will add it as a member function to triangulation, and likely handle it in the main evap class
    
    double f3 = max <double> (1.0, p/40.0);

    double f4 = 1.0;
    if (var.t < 5.0 || var.t > 40.0)
        f4 = 5000/50;

    if (var.ShortWave_in <= 0)
        var.stomatal_resistance = 5000;
    else
    {
        var.stomatal_resistance = min <double> (rcstar * f1 * f2 * f3 * f4, 
                5000.0);
    }

}

void PenmanMonteith::CalcEvapT(param_base& baseparam, var_base& basevar, model_output& output) const
{
    const PM_param & param = static_cast<const PM_param&>(baseparam);
    const PM_var & var = static_cast<const PM_var&>(basevar);
    double Q =  var.Rnet * (1 - param.Frac_to_ground);

    CalcAeroResistance(param,var);
    CalcStomatalResistance(param,var);

    output.ET = ( delta(var.t) * Q + AirDensity(var.t,var.vapour_pressure,var.P_atm) * param.Cp / (lambda(var.t)*1e3) * ( CalcSaturationVapourPressure(var.t) - var.vapour_pressure )/ var.aero_resistance )
       / ( delta(var.t) + gamma(var.P_atm, var.t) * (1 + stomatal_resistance / aero_resistance ) ); 
    return;
}

// TODO make delta, gamma, density fucntions

double PenmanMonteith::AirDensity(double& t, double& ea, double& Pa) // atmospheric density (kg/m^3)
{
  const double R0 = 2870;
   return (1E4*Pa /(R0*( 273.15 + t))*(1.0 - 0.379*(ea/Pa))); //
}

