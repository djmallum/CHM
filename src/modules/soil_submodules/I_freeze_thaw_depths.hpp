#pragma once
class I_freeze_thaw_depths
{
public:
    virtual ~I_freeze_thaw_depths() = default;
    virtual void run(void) = 0;

    double thaw_front_depth = 0.0;
    double freeze_front_depth = 0.0;

    double get_freeze_front() { return freeze_front_depth;};

    double get_thaw_depth() { return thaw_front_depth; };
};
