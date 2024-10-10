#pragma once

#include "soil_base.hpp"
#include "soil.hpp"

// Temporary forward declaration
// TODO add include after writing this
class k_estimate;

struct input_soil_layers : input_base
{
    double& Zdt;
    double& Zdf;
    double& potential_ET;
    double& ET;    
    double& swe;
    double& infil;
    double& runoff;
};

struct output_soil_layers : output_base
{
    double condensation;
    double actual_ET;
    double soil_excess_to_runoff;
    double soil_excess_to_gw;
    double gw_flow;
    double soil_to_ssr;
    double ET;
}

class soil_layers : public soil_base
{
public: // methods
    soil_layers(int& _num_layers, soil_module::data& _d);
    ~soil_layers();

    void run (input_base& input, output_base& output) override;

private: //methods
    
    
    int _find_thaw_front_layer(const std::vector<double>& depths, double& value);
    void set_K_values();
    void set_thaw_fraction(double value, int layer); // set the new value in a layer, ensure values are between 0 and 1.

public: //members
    // change if modified by other classes
    double detention_storage;
    double depression_storage;
    double soil_lower;
    double soil_storage;
    double soil_rechr;
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
