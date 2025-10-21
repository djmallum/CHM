#include "double_range.hpp"
#include <stdexcept>

namespace double_range
{
    Unit::Unit(const double c) : _c(c)
    {
        if (c < 0.0 || c > 1.0)
            throw std::logic_error("cutoff must be between (inclusive) 0 and 1");
    };
    Unit::operator double() const { return _c; };
};
