#pragma once

#include <cmath>

struct input_base
{

};

struct output_base
{
    
};

class submodule_base
{
public:
    virtual ~submodule_base = default;

    virtual void run(input_base& input, output_base& output) = 0;
    
};

