#include "Ayers.hpp"

Ayers::Ayers(const double& _rainfall, const double& _snowmelt, const std::string& _texture, const std::string& _ground_cover, const soil_data& _soil_data) : rainfall(_rainfall), snowmelt(_snowmelt), texture(_texture), ground_cover(_ground_cover), soils(_soil_data)
{
};

void Ayers::run()
{
    double maxinfil = (soils.*max_infil)(texture,ground_cover); 
    if (maxinfil > rainfall)
    {
        inf = rainfall;
    }
    else
    {
        inf = maxinfil;
        runoff = rainfall - maxinfil;
    }
}
