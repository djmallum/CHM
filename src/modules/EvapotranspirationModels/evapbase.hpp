#pragma once

#include <cmath>

struct var_base
{

};

struct model_output
{
    double ET;
};

class evapT_base
{
public:
    virtual ~evapT_base() = default;

    virtual void CalcEvapT(var_base& basevar, model_output& output) = 0;
    
    double delta(double& t);
    double gamma(double& P_atm, double& t);
    double lambda(double& t);
};

