#pragma once
class I_freeze_thaw_depths
{
public:
    virtual ~I_freeze_thaw_depths() = default;
    virtual void run(void) = 0;
};
