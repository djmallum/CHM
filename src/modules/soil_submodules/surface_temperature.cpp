void surface_temperature::run()
{
	if (DTO.swe <= 0) // bare ground
	{
		if (DTO.Zdt > Zdt_last)
			Zdt_last = DTO.Zdt;
		
		double Qn = DTO.netD; // W/m^2, CRHM converts, maybe I should too
		DTO.t_surface  = (DTO.W_a * DTO.t + DTO.W_b*Qn) * atan(DTO.W_c * (Zdt_last + DTO.W_d)) * 2.0 / 3.14159265; //TODO check if M_PI is CRHM is just pi.
	}
	else //snowcover
	{
		double SWE_tc = 0.0;
		else if (DTO.snow_density < 156) // Sturm et al. 1997. The thermal conductivity of seasonal snow
			SWE_tc = 0.023 - 1.01 * DTO.snow_density / 1000.0 + 0.234*pow(snow_density/1000.0,2.0); // TODO check if snow density /1000 makes sense, acutally might be dividied by reference density.
		else 
			SWE_tc = 0.138 - 1.01 * DTO.snow_density / 1000.0 + 3.233*pow(snow_density/1000.0,2.0);
		DTO.t_surface = yesterday_snow_temp - DTO.G*0.5*DTO.snow_depth/SWE_tc;
	
	
		if (DTO.t > 0.0) // If snowcovered, ignore positive temperatures.
			DTO.t_surface = 0.0;
		else
			DTO.t_surface = DTO.t; // TODO looks like a bug, looks like t_surface is never set.
								   // compared with CRHM and email logan
								   // I should also look at CRHM in vim just to make sure I got the {} 
								   // bounds correct
	}
	
	// TODO add the rest here, stuff about daily, is the daily value used in XG?
}