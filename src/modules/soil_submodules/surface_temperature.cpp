#include "surface_temperature.hpp"



void surface_temperature::run()
{
    if (DTO.swe <= 0)
        set_tsurface_bare(); // 
    else if (DTO.t > 0.0) // snowcovered but warm out
        DTO.t_surface = 0.0;
    else
        set_tsurface_snow();

    if (DTO.is_day_start())
        set_daily_tsurface();

    increment_daily_counters();

};

void surface_temperature::set_tsurface_base(void)
{
    if (DTO.Zdt > Zdt_last)
        Zdt_last = DTO.Zdt;
		
	double Qn = DTO.netD; // W/m^2, CRHM converts, maybe I should too
	DTO.t_surface  = (DTO.W_a * DTO.t + DTO.W_b*Qn) * atan(DTO.W_c * (Zdt_last + DTO.W_d)) * 2.0 / 3.14159265; //TODO check if M_PI is CRHM is just pi.
};

void surface_temperature::set_tsurface_snow(void)
{

    double SWE_tc = 0.0;
	if (DTO.snow_density < 156) // Sturm et al. 1997. The thermal conductivity of seasonal snow
	    SWE_tc = 0.023 + 0.234*DTO.snow_density/1000.0; // TODO check if snow density /1000 makes sense, acutally might be dividied by reference density.
	else 
		SWE_tc = 0.138 - 1.01 * DTO.snow_density / 1000.0 + 3.233*pow(DTO.snow_density/1000.0,2.0);

    if (yesterday_snow_temp < -70.0)
        DTO.t_surface = DTO.t;
    else
        DTO.t_surface = yesterday_snow_temp + DTO.Ground_heat_flux*0.5*DTO.snow_depth/SWE_tc; 
 
};

void surface_temperature::set_daily_tsurface()
{
    int steps_per_day_fraction = 1 / (86400 / get_dt(DTO));
    if (steps_per_day_fraction > 0)
        DTO.call_model_error();
            
    daily_t_surface = t_surface_total * steps_per_day_fraction;

    yesterday_snow_temp = snow_temp_total * steps_per_daiy_fraction;
        
    t_surface_total = 0.0;
    snow_temp_total = 0.0;

};
void surface_temperature::increment_daily_counters()
{
    t_surface_total += DTO.t_surface;
    snow_temp_total += DTO.snow_temp;
};
