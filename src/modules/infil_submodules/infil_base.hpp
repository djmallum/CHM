#pragma once

#include <submodule_base.hpp>

struct output : _output_base
{
    double inf = 0.0;
    double snowinf = 0.0;
    double runoff = 0.0;
    double melt_runoff = 0.0;
    double rain_on_snow = 0.0;
};

struct _soil_input : _input_base
{
    // something
}

class _soil_base : public _submodule_base
{
public:
    _soil_base();
    ~_soil_base();

    virtual void run(_input_base& input, _output_base& output) override = 0;

};
        
