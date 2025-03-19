#pragma once

#include "I_K_estimate.hpp"
#include "soil_DTO.hpp"
#include <memory>

class surface_temperature
{
public:
    surface_temperature(two_layer_DTO& _DTO) : DTO(_DTO) 
    {

    };
    ~surface_temperature() {};

    void run(void);
private:

    two_layer_DTO& DTO;
	double Zdt_last = 0.0;
    void set_tsurface_bare(void);
    void set_tsurface_snow(void);
    void set_daily_tsurface(void);
    void increment_daily_counters();


     
};
