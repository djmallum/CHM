#pragma once

#include "soil_two_layer.hpp"

soil_two_layer::soil_two_layer(double& _num_layers) : num_layers(_num_layers)
{
    thaw_layer_fraction(num_layers);        
};

soil_two_layer::~soil_two_layer(void) 
{

};

void soil_two_layer::run(input_base& input, output_base& output) 
{
    const input_soil_two_layer& _input = static_cast<const input_soil_two_layer&>(input);
    const output_soil& _output = static_cast<const output_soil&>(output);

    set_layer_thaw_fraction(_input);

    set_condensation(_input.ET,_output.condensation);

    organize_soil_layers(_input,_output);

    set_detention_storage(_input,_output);

    set_depression_storage(_input,_output);

    set_ground_water(_input,_output);

    set_subsurface_runoff(_input,_output);

    set_soil_ET(_input,_output);

};

void soil_two_layer::organize_soil_layers(input_soil_two_layer& _input,output_soil_two_layer& _output);
{
    double excess = 0.0; // mm
    if (soil_storage_max > 0.0)
    {
        double soil_lower_storage = soil_storage - soil_rechr_storage;

        double potential = _input.infil + _output.condensation;

        double possible = layer_thaw_fraction[0] * (soil_rechr_max - soil_rechr_storage);

        if (possible > potential)
            possible = potential;
        else
            soil_excess_to_runoff = potential - possible;

        soil_rechr_storage += possible;

        if (soil_rechr_storage > soil_rechr_max)
            _push_excess_down(soil_rechr_storage,soil_rechr_max,soil_lower_storage);

        soil_storage = soil_lower_storage + soil_rechr_storage;

        if (soil_storage > soil_storage_max)
            _push_excess_down(soil_storage,soil_storage_max,soil_excess_to_gw);

        // TODO add code for rechr_ssr, see line 555 of soilX, in the crhmcomments branch. I asked logan why this is there. 

        if (soil_excess_to_gw > K_soil_to_gw * layer_thaw_fraction[1])
        {
            double excess_to_gw_max = K_soil_to_gw * layer_thaw_fraction[1];
            _push_excess_down(soil_excess_to_gw,excess_to_gw_max,excess);
        }

        // Line 607 of SoilX crhmcommetns branch, comment says upper layer but code says lower-layer, ask logan about this.
        if (excess_to_ssr && excess > 0.0)
        {
            double excess_to_ssr_max = excess * (1.0 -  layer_thaw_fraction[1]);
            _push_excess_down(excess,excess_to_ssr_max,soil_to_ssr);
        }

        
    }
    else
        excess = _input.infil + _output.condensation;

    
};

void soil_two_layer::_push_excess_down(double& layer_storage, double& layer_max, double& layer_down)
{
    double excess = layer_storage - layer_max;

    layer_storage = layer_max;

    layer_down += excess;
}
    


void soil_two_layer::set_layer_thaw_fraction()
{
    std::fill(layer_thaw_fraction.begin(),layer_thaw_fraction.end(),0.0);


    for (int layer = 0; layer < num_layers; ++layer)
    {
        if (thaw_front_depth == 0.0 && freeze_front_depth == 0.0)
            layer_thaw_fraction[layer] = 1.0; 
        else 
        {
            
            layer_thaw_fraction[layer] = _find_thaw_fraction(layer);
            if (layer_thaw_fraction[layer] < 1.0)
                break;
        }
    }

    
};

double soil_two_layers::_find_thaw_fraction(int& layer)
{
    double layer_start = 0.0 // mm
    if (depth[layer] >= thaw_front_depth)
    {
        layer_start = std::accumulate(depth.begin(), depth.begin + layer + 1, 0.0);
        return (thaw_front_depth - layer_start) / depth[layer];
    else
        return 1.0;
    }
};

void soil_two_layers::set_condensation(double& ET, double& condense)
{
    if (ET < 0.0)
    {
        condense = -1.0 * ET;
        ET = 0.0;
    }
};
void soil_two_layer::set_data(soil_module::data& _d)
{
    d = &_d;
};




