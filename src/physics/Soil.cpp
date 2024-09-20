#include "physics/Soil.h"

namespace Soil
{
    double SoilData::get_soilproperties(const std::size_t property_type, const std::size_t soil_type)
    {
        return soilproperties[property_type][soil_type];
    }


    double SoilData::get_textureproperties(const std::size_t texture_type, const std::size_t groundcover_type)
    {
        return textureproperties[texture_type][groundcover_type];
    }
}
