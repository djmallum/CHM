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
        
        if (_input.swe == 0.0) // if there is no snowcover
        {
            rechr_to_ssr = soil_rechr_storage / soil_rechr_max * K_rechr_to_ssr * layer_thaw_fraction[0];

            if (rechr_to_ssr > soil_rechr_storage * layer_thaw_fraction[0])
                rechr_to_ssr = soil_rechr_storage * layer_thaw_fraction[0];

            soil_rechr_storage = std::max(0.0, soil_rechr_storage - rechr_to_ssr);

            soil_storage -= rechr_to_ssr;
            soil_to_ssr = rechr_to_ssr;

        }

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

        
        
    else
        excess = _input.infil + _output.condensation;

    
};

void soil_two_layer::manage_detention(input_soil_two_layers& _input, output_soil_two_layer& _output,double& excess)
{
    _output.soil_excess_to_runoff += _input.runoff + excess + routing_residual / face_area; // routing_residual comes from the crhm varaible redirected_residual which has units of mm*km^2/int (not sure why), so face_area is there for now for consistency.
    // TODO is the above in _output?                                   
        
    routing_residual = 0.0;

    if (_output.soil_excess_to_runoff > 0.0)
    {
        if (_input.swe <= snow_covered_threshold)
            detention_max = detention_snow_max;
        else
            detention_max = detention_organic_max;

        double detention_space = detention_max - detention_storage;

        if (detention_space > 0.0)
        {
            if (_output.soil_excess_to_runoff > detention_space)
            {
                _output.soil_excess_to_runoff = std::max(0.0,_output.soil_excess_to_runoff - detention_space); 
                detention_storage += detention_space;
            }
            else
            {
                detention_storage += _output.soil_excess_to_runoff;
                soil_excess_to_runoff = 0.0;
            }
        }
    }

    if (_input.swe <= snow_covered_threshold) // default will be zero
        K_detention_to_runoff = K_detention_organic_to_runoff;
    else
        K_detention_to_runoff = K_detention_snow_to_runoff;

    if (detention_storage > 0.0 && K_detention_to_runoff > 0.0)
    {
        double transfer = std::min(detention_storage,K_detention_storage);
        _output.soil_excess_to_runoff += transfer;
        detention_storage -= transfer;
    }

    if (detention_storage < 0.0001) // from CRHM, for safety and to drop any Floating-point errors
        detention_storage = 0.0;

};

void soil_two_layer::manage_depression(input_soil_two_layers& _input, output_soil_two_layers& _output)
{
    if (_output.soil_excess_to_runoff > 0.0 && depression_storage > 0.0)
    {
        double exponent = -1.0 * std::min(12.0,_output.soil_excess_to_runoff / depression_max);

        double depression_space = (depression_max - depression_storage) * (1 - exp(exponent));        

        if (soil_storage_max == 0.0)
            depression_space = depression_max - depression_storage;

        if (depression_space > 0.0)
        {
            if (_output.soil_excess_to_runoff > depression_space)
            {
                _output.soil_excess_to_runoff = std::max(0.0,_output_soil_excess_to_runoff - depression_space);
                depression_storage += depression_space;
                // TODO add total tracker
            }
            else
            {
                depression_storage += _output.soil_excess_to_runoff;
                // TODO add total tracker
                _output.soil_excess_to_runoff = 0.0;
            }

        }
    }

    if (depression_storage > 0.0 && K_depression_to_gw > 0.0)
    {
        double value = transfer_min(depression_storage,K_depression_to_gw);
        soil_depresssion_to_gw += value;
        depression_storage -= value;
        if (depression_storage < 0.0) // floatig point error safety?
            depression_storage = 0.0;
    }


};

double soil_two_layer::transfer_min(double& val1, double& val2)
{
    return std::min(val1,val2);
};

void soil_two_layer::manage_groundwater()
{
    ground_water_storage += soil_depression_to_gw;
    ground_water_flow = 0.0;

    if (ground_water_storage > ground_water_max)
    {
        _push_excess_down(ground_water_storage,ground_water_max,ground_water_flow);
    }

    if (ground_water_max > 0.0) // divide by zero safety
    {
        double spilled = ground_water_storage / ground_water_max * K_ground_water_out;
        ground_water_storage -= spill;
        ground_water_flow += spill;
    }

    // TODO daily and simulation totals, like in CRHM?

};

void soil_two_layer::manage_subsurface_runoff(input_soil_two_layers& _input, output_soil_two_layers& _output)
{
    if (depression_storage > 0.0 && K_depression_to_ssr > 0.0)
    {
        double value = transfer_min(depression_storage,K_depression_to_ssr);
        soil_to_ssr += value;
        depression_storage -= value;
        if (depression_storage < 0.0)
            depression_storage = 0.0;
    }

    if (K_lower_to_ssr > 0.0)
    {
        double available = soil_storage - soil_rechr;
        double value = transfer_min(K_lower_to_ssr * layer_thaw_fraction[1],available);
        soil_storage -= value;
        soil_to_ssr += value;
    }

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




