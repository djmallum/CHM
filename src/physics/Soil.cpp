#include "Soil.h"

namespace Soil
{

    soils_na::soils_na()
    {
        _make_hash();
    }

    soils_na::~soils_na()
    {

    }

    double soils_na::porosity(std::string soil_type)
    {
        return _porosity[soil_type];
    }

    double soils_na::pore_size_dist(std::string soil_type)
    {
        return _pore_size_dist[soil_type];
    }

    double soils_na::pore_space(std::string soil_type)
    {
        return _pore_space[soil_type];
    }

    double soils_na::wilt_point(std::string soil_type)
    {
        return _wilt_point[soil_type];
    }

    double soils_na::air_entry_tension(std::string soil_type)
    {
        return _air_entry_tension[soil_type];
    }

    double soils_na::saturated_conductivity(std::string soil_type)
    {
        return _saturated_conductivity[soil_type];
    }

    double soils_na::_make_hash()
    {
        // illegal to not include this step for the dense_hash_map
        _porosity.set_empty_key("");
        _pore_size_dist.set_empty_key("");
        _pore_space.set_empty_key("");
        _wilt_point.set_empty_key("");
        _air_entry_tension.set_empty_key("");
        _saturated_conductivity.set_empty_key("");



        
    }

















    double soils_na::get_textureproperties(const std::size_t texture_type, const std::size_t groundcover_type)
    {
        return textureproperties[texture_type][groundcover_type];
    }
}
