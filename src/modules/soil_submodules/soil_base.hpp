#pragma once

#include <submodule_base.hpp>

struct input_soil_base : input_base
{
};

struct output_soil : output_base
{
    double soil_storage;
};

class soil_base : public submodule_base
{
public:
    virtual soil_base = default;
    virtual ~soil_base = default;

    virtual void run(input_base& input, output_base& output) = 0;

};

