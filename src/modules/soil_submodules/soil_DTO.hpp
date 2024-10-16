class two_layer_DTO
{
public:
    /*
********************************************************************
*                             GETTERS                              *
********************************************************************
    */
    double get_thaw_front_depth() {return thaw_front_depth}; //in - XG
    double get_freeze_front_depth() {return freeze_front_depth}; //in - XG
    double get_potential_ET() {return potential_ET}; //in - evapT
    double get_swe() {return swe}; //in - snobal TODO should this include input from PBSM3D?
    double get_infil() {return infil}; //in - infil_all
    double get_runoff() {return runoff}; //in - infil_all 
    double get_routing_residual() {return routing_residual}; //in - NetRoute

    double get_condensation() { return condensation}; //out (and used)
    double get_actual_ET() { return actual_ET}; //out
    double get_soil_excess_to_runoff() { return soil_excess_to_runoff}; //out (and used)
    double get_soil_excess_to_gw() { return soil_excess_to_gw}; //out (and used)
    double get_ground_water_outflow() { return ground_water_outflow}; //out (and used)
    double get_soil_to_ssr() { return soil_to_ssr}; //out (and used)
			// All K_ stuff are in
    double get_K_soil_to_gw() { return K_soil_to_gw}; // All these K's could be swapped to held from previous and could be global/system wide
    double get_K_rechr_to_ssr() { return K_rechr_to_ssr};
    double get_K_lower_to_ssr() { return K_lower_to_ssr};
    double get_K_detention_snow_to_runoff() { return K_detention_snow_to_runoff};
    double get_K_detention_soil_to_runoff() { return K_detention_soil_to_runoff};
    double get_K_depression_to_ssr() { return K_depression_to_ssr};

    // held from previous 
    double get_soil_storage_max() { return soil_storage_max};
    double get_soil_storage() { return soil_storage};
    double get_soil_rechr_storage() { return soil_rechr_storage};
    std::vector<double> layer_thaw_fraction; //2
    double get_soil_rechr_max() { return soil_rechr_max};
    double get_soil_rechr_storage() { return soil_rechr_storage};
    bool get_excess_to_ssr() { return excess_to_ssr}; 
    double get_detention_max() { return detention_max};
    double get_detention_snow_max() { return detention_snow_max};
    double get_detention_organic_max() { return detention_organic_max};
    double get_detention_snow_init() { return detention_snow_init};
    double get_detention_organic_init() { return detention_organic_init};
    double get_depression_max() { return depression_max};
    double get_depression_storage() { return depression_storage};
    double get_ground_water_storage() { return ground_water_storage};
    double get_ground_water_max() { return ground_water_max};
    double get_ground_cover_type() { return ground_cover_type};
    int get_soil_type() { return soil_type};

    /*
********************************************************************
*                             SETTERS                              *
********************************************************************
    */

    void set_thaw_front_depth(double value) {thaw_front_depth = value}; //in - XG
    void set_freeze_front_depth(double value) {freeze_front_depth = value}; //in - XG
    void set_potential_ET(double value) {potential_ET = value}; //in - evapT
    void set_swe(double value) {swe = value}; //in - snobal TODO should this include input from PBSM3D?
    void set_infil(double value) {infil = value}; //in - infil_all
    void set_runoff(double value) {runoff = value}; //in - infil_all 
    void set_routing_residual(double value) {routing_residual = value}; //in - NetRoute

    void set_condensation(double value) {condensation = value}; //out (and used)
    void set_actual_ET(double value) {actual_ET = value}; //out
    void set_soil_excess_to_runoff(double value) {soil_excess_to_runoff = value}; //out (and used)
    void set_soil_excess_to_gw(double value) {soil_excess_to_gw = value}; //out (and used)
    void set_ground_water_outflow(double value) {ground_water_outflow = value}; //out (and used)
    void set_soil_to_ssr(double value) {soil_to_ssr = value}; //out (and used)
			// All K_ stuff are in
    void set_K_soil_to_gw(double value) {K_soil_to_gw = value}; // All these K's could be swapped to held from previous and could be global/system wide
    void set_K_rechr_to_ssr(double value) {K_rechr_to_ssr = value};
    void set_K_lower_to_ssr(double value) {K_lower_to_ssr = value};
    void set_K_detention_snow_to_runoff(double value) {K_detention_snow_to_runoff = value};
    void set_K_detention_soil_to_runoff(double value) {K_detention_soil_to_runoff = value};
    void set_K_depression_to_ssr(double value) {K_depression_to_ssr = value};

    // held from previous 
    void set_soil_storage_max(double value) {soil_storage_max = value};
    void set_soil_storage(double value) {soil_storage = value};
    void set_soil_rechr_storage(double value) {soil_rechr_storage = value};
    void set_thaw_fraction_rechr(double value) {thaw_fraction_rechr = value};
    void set_thaw_fraction_lower(double value) {
    std::vector<double> layer_thaw_fraction; //2
    void set_soil_rechr_max(double value) {soil_rechr_max = value};
    void set_soil_rechr_storage(double value) {soil_rechr_storage = value};
    void set_excess_to_ssr(bool value) {excess_to_ssr = value}; 
    void set_detention_max(double value) {detention_max = value};
    void set_detention_snow_max(double value) {detention_snow_max = value};
    void set_detention_organic_max(double value) {detention_organic_max = value};
    void set_detention_snow_init(double value) {detention_snow_init = value};
    void set_detention_organic_init(double value) {detention_organic_init = value};
    void set_depression_max(double value) {depression_max = value};
    void set_depression_storage(double value) {depression_storage = value};
    void set_ground_water_storage(double value) {ground_water_storage = value};
    void set_ground_water_max(double value) {ground_water_max = value};
    void set_ground_cover_type(double value) {ground_cover_type = value};
    void set_soil_type(int value) {soil_type = value};

private:
    double get_thaw_front_depth; //in - XG
    double freeze_front_depth; //in - XG
    double potential_ET; //in - evapT
    double swe; //in - snobal TODO should this include input from PBSM3D?
    double infil; //in - infil_all
    double runoff; //in - infil_all 
    double routing_residual; //in - NetRoute

    double condensation; //out (and used)
    double actual_ET; //out
    double soil_excess_to_runoff; //out (and used)
    double soil_excess_to_gw; //out (and used)
    double ground_water_outflow; //out (and used)
    double soil_to_ssr; //out (and used)
			// All K_ stuff are in
    double K_soil_to_gw; // All these K's could be swapped to held from previous and could be global/system wide
    double K_rechr_to_ssr;
    double K_lower_to_ssr;
    double K_detention_snow_to_runoff;
    double K_detention_soil_to_runoff;
    double K_depression_to_ssr;

    // held from previous 
    double soil_storage_max;
    double soil_storage;
    double soil_rechr_storage;
    double thaw_fraction_rechr;
    double thaw_fraction_lower;
    double soil_rechr_max;
    double soil_rechr_storage;
    bool excess_to_ssr; 
    double detention_max;
    double detention_snow_max;
    double detention_organic_max;
    double detention_snow_init;
    double detention_organic_init;
    double depression_max;
    double depression_storage;
    double ground_water_storage;
    double ground_water_max;
    double ground_cover_type;
    int soil_type;
};

