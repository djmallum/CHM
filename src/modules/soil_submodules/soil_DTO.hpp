struct two_layer_DTO
{
    double thaw_front_depth = 0.0; //in - XG
    double freeze_front_depth = 0.0; //in - XG
    double potential_ET = 0.0; //in - evapT
    double swe = 0.0; //in - snobal TODO should this include input from PBSM3D?
    double infil = 0.0; //in - infil_all
    double runoff = 0.0; //in - infil_all 
    double routing_residual = 0.0; //in - NetRoute

    double condensation = 0.0; //out (and used)
    double actual_ET = 0.0; //out
    double soil_excess_to_runoff = 0.0; //out (and used)
    double soil_excess_to_gw = 0.0; //out (and used)
    double ground_water_outflow = 0.0; //out (and used)
    double soil_to_ssr = 0.0; //out (and used)
    double excess = 0.0;

			// All K_ stuff are in
    double K_soil_to_gw = 0.0; // All these K's could be swapped to held from previous and could be global/system wide
    double K_rechr_to_ssr = 0.0;
    double K_lower_to_ssr = 0.0;
    double K_detention_snow_to_runoff = 0.0;
    double K_detention_organic_to_runoff = 0.0;
    double K_depression_to_ssr = 0.0;

    // held from previous 
    double soil_storage_max = 0.0;
    double soil_storage = 0.0;
    double soil_rechr_storage = 0.0;
    double thaw_fraction_rechr = 0.0;
    double thaw_fraction_lower = 0.0;
    double soil_rechr_max = 0.0;
    double soil_rechr_storage = 0.0;
    bool excess_to_ssr = true; 
    double detention_max = 0.0;
    double detention_snow_max = 0.0;
    double detention_organic_max = 0.0;
    double detention_snow_init = 0.0;
    double detention_organic_init = 0.0;
    double depression_max = 0.0;
    double depression_storage = 0.0;
    double ground_water_storage = 0.0;
    double ground_water_max = 0.0;
    double ground_cover_type = 0.0;
    int soil_type = 0;
    
};

