#pragma once
#include <stdexcept>

class daily_accumulator
{
private:
    double* _target_var = nullptr;
    double _accumulator = 0.0;
    double _mean_value = 0.0;
    
public:
    daily_accumulator() {};

    void bind_target(double& target)
    {
        _target_var = &target;
    };

    const double& get_last_mean() const { return _mean_value; };

    void accumulate(bool is_new_day, double steps_per_day)
    {
        if (!_target_var) 
            throw std::logic_error("_target_var should be bound before calling execute()");

        if (is_new_day)
        {
            _mean_value = _accumulator / steps_per_day;
            _accumulator = 0.0;
        };

        _accumulator += *_target_var;
    };

};
