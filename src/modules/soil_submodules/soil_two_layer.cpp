#pragma once

#include "soil_two_layer.hpp"

soil_layers::soil_layers(double& _num_layers) : num_layers(_num_layers)
{
    
};

soil_layers::~soil_layers(void) 
{

};

void soil_layers::run(input_base& input, output_base& output) 
{
    const input_soil_layers& _input = static_cast<const input_soil_layers&>(input);
    const output_soil& _output = static_cast<const output_soil&>(output);

    // daily step
    //
    // condensation
    //
    // infiltrate soil and organize
    //
    // detention
    //
    // depression
    //
    // groundwater
    //
    // subsurface runoff
    //
    // ET
    //

};

void soil_layers::set_data(soil_module::data& _d)
{
    d = &_d;
};




