#include "soil_two_layer.hpp"

soil_two_layer::soil_two_layer(two_layer_DTO& DTO,I_K_estimate& estimate) : 
    init(DTO),
    condensation(DTO), 
    infiltrate(DTO),
    detention(DTO),
    depression(DTO),
    groundwater(DTO),
    ssr(DTO),
    k_estimator(estimate)
{
};
void soil_two_layer::run() 
{
    
    init.zero_single_step_vars();
    init.layer_thaw_fraction();

    condensation.set();

    infiltrate.distribute();

    detention.manage();

    depression.manage();

    groundwater.manage();

    ssr.manage();
    
    k_estimator.run();
};
