#pragma once

#include <cmath>


class submodule_base
{
public:
    virtual ~submodule_base() = default;

    virtual void run() = 0;
    
};

