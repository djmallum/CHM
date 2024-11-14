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
    ~K_estimate() {};

    void run(void);

private:

    two_layer_DTO& DTO;
	double Zdt_last = 0.0;


     
};
