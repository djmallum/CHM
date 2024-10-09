#pragma once

#include <cmath>

struct _input_base
{

};

struct _output_base
{
    
};

class _submodule_base
{
public:
    virtual ~_submodule_base = default;

    virtual void run(_input_base& input, _output_base& output) = 0;
    
};

