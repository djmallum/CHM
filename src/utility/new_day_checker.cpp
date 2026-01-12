#include "new_day_checker.hpp"

bool new_day_checker::check() const
{
    // TODO This has hard coded elements, Chris suggested something different here: https://godbolt.org/z/3c51T1avT
	auto td = _global->posix_time().time_of_day().total_seconds();
    auto time_to_midnight = SECONDS_PER_DAY - td;
    if (td >= 0 && td < _global->dt()) //(time_to_midnight >= _global->dt())
    {
        return true;
    }
    else
        return false;
};

void new_day_checker::set_global(std::shared_ptr<global> _g) { 
    _global = _g;
};
