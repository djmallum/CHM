#include "global.hpp"
#include <memory>

class new_day_checker
{
public:
    bool check() const;

    void set_global(std::shared_ptr<global> _g);

private:
    std::shared_ptr<global> _global;
    static constexpr size_t SECONDS_PER_DAY = 86400;
};
