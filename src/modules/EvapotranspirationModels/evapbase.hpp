#pragma once

#include <cmath>
#include "logger.hpp"

struct var_base
{

};

struct model_output
{
    double ET;
    double rc;
};

class evapT_base
{
public:
    virtual ~evapT_base() = default;

    virtual void CalcEvapT(var_base& basevar, model_output& output) = 0;
    
    static double delta(const double t);
    static double gamma(const double P_atm, const double t, const double c_a);
    static double lambda(const double t);
};

