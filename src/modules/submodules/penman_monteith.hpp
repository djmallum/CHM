#include "base_step.hpp"
#include <concepts>

template<class T>
concept penman_data = requires(T& t)
{
    {t.leaf_area_index()} -> std::convertible_to<double>;
    
    {t.leaf_area_index_max()} -> std::convertible_to<double>;

    {t.vegetation_height()} -> std::convertible_to<double>;

    {t.wind_measurement_height()} -> std::convertible_to<double>;

    {t.stomatal_resistance_min()} -> std::convertible_to<double>;

    {t.soil_depth()} -> std::convertible_to<double>;    

    {t.radiation_to_ground()} -> std::convertible_to<double>;

    {t.s_per_step()} -> std::convertible_to<double>;

    {t.heat_capacity_air()} -> std::convertible_to<double>;

    {t.kappa()} -> std::convertible_to<double>;

    {t.air_entry_tension()} -> std::convertible_to<double>;

    {t.pore_size_dist()} -> std::convertible_to<double>;

    {t.wilt_point()} -> std::convertible_to<double>;

    {t.porosity()} -> std::convertible_to<double>;

    // might want these to be const doubles, so I could add the requires above...
    {t.stomatal_resistance(0.0)} -> std::convertible_to<void>;

    {t.ET(0.0)} -> std::convertible_to<void>;
};

template<penman_data data>
class penman_monteith : public base_step<data>
{
public:
    explicit penman_monteith() {};
    ~penman_monteith() {};

    void execute(data& d) final;
};

template<penman_data data>
void penman_monteith<data>::execute(data& d)
{

};
