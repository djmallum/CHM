#pragma once

#include "soil_base.hpp"
#include "soil.hpp"


class two_layer_DTO
{
public:
    // note does not include actual ET except for main variables
    // new each time step
    // TODO this is memory inefficient, could store in/out variables as unique_ptr that is started when ready.
    // TODO split soil modules into submodules (one for gw, one for settting up ssr, one for layers, etc.)
    // TODO after splitting, make a separate DTO which is then combined into a single DTO (this one).
    double thaw_front_depth; //in - XG
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
    std::vector<double> layer_thaw_fraction; //2
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

    // global constants
    double num_layers;
};
// Temporary forward declaration
// TODO add include after writing this
class k_estimate;

struct input_soil_two_layer : input_base
{
    double& thaw_front_depth;
    double& freeze_front_depth;
    double& potential_ET;
    double& ET;    
    double& swe;
    double& infil;
    double& runoff;
};

struct output_soil_two_layer : output_base
{
    double condensation = 0.0;
    double actual_ET 0.0;
    double soil_excess_to_runoff = 0.0;
    double soil_excess_to_gw = 0.0;
    double gw_flow = 0.0;
    double soil_to_ssr = 0.0;
}

class soil_two_layer : public soil_base
{
public: // methods
    soil_two_layer(int& _num_layers, soil_module::data& _d);
    ~soil_two_layer();

    void run(input_base& input, output_base& output) override;

private: //methods
    
    
    int _find_thaw_front_layer(const std::vector<double>& depths, double& value);
    void set_K_values();
    void set_thaw_fraction(double value, int layer); // set the new value in a layer, ensure values are between 0 and 1.

public: //members
    // change if modified by other classes
    double detention_storage;
    double depression_storage;
    double soil_storage;
    double soil_rechr_storage;
    double ground_water_storage;
    const double routing_residual = 0.0;


private: //members
    int& num_layers;
    double& soil_storage_max;
    double& soil_rechr_max;
    double& detention_storage_max;
    double& detention_max;
    double& detention_snow_max;
    double& detention_soil_max;
    double& depression_max;
    double& detention_snow_init; // Might make these passed through the constructor jsut for initialization of detention storage
    double& detention_soil_init; //same as above
    double& ground_water_max;
    double& ground_cover_type;
    int soil_type;
    double K_depression_to_ssr;
    double K_detention_snow_to_runoff;
    double K_detention_soil_to_runoff;
    double K_rechr_to_ssr;
    double K_lower_to_ssr;
    double K_soil_to_gw;
    bool& excess_to_ssr; //global
    
    k_estimate k_estimator; 
    
    std::vector<double> layer_thaw_fraction; // thaw_layers_lay
    std::vector<double> 

};
