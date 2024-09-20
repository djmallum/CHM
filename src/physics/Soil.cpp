#include "physics/Soil.h"

namespace Soil
{

    SoilData::SoilData()
    {

    }

    SoilData::~SoilData()
    {

    }

    double SoilData::get_soilproperties(const std::size_t soil_type, const std::size_t property_type)
    {
        return soilproperties[soil_type][property_type];
    }


    double SoilData::get_textureproperties(const std::size_t texture_type, const std::size_t groundcover_type)
    {
        return textureproperties[texture_type][groundcover_type];
    }
}
