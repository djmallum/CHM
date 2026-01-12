#pragma once
#include <stdexcept>
class daily_accumulator
{
private:
    const double* _target_var = nullptr;
    double _accumulator = 0.0;
    double _yesterday_value = 0.0;
    
public:
    daily_accumulator() {};

    void bind_target(const double& target)
    {
        _target_var = &target;
    };

    const double get_yesterday() const { return _yesterday_value; };

    void accumulate(bool is_new_day, double steps_per_day)
    {
        if (!_target_var) 
            throw std::logic_error("_target_var should be bound before calling accumulate(bool,double)");

        if (is_new_day)
        {
            [[unlikely]] _yesterday_value = _accumulator / steps_per_day;
            _accumulator = *_target_var;
        }
        else
            _accumulator += *_target_var;
    };

    void accumulate(bool is_new_day)
    {
        if (!_target_var)
            throw std::logic_error("_target_var should be bound before calling accumulate(bool)");

        if (is_new_day) 
        {
            [[unlikely]] _yesterday_value = _accumulator;
            _accumulator = *_target_var;
        }
        else
            _accumulator += *_target_var;
    };

    bool operator==(const daily_accumulator& other) const
    {
        return _yesterday_value == other._yesterday_value;
    };
};
