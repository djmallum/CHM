#include "freeze_thaw_soils.hpp"

REGISTER_MODULE_CPP(soil_module);

freeze_thaw_soils::soil_module(config_file cfg) : module_base("soil_module", parallel::data, cfg)
{
    depends("surface_temp");

    provides("thaw_front_depth");
    provides("freeze_front_depth");
    provides("first_front_depth");
};

freeze_thaw_soils::~soil_module()
{
    
};

void freeze_thaw_soils::init(mesh& domain)
{

    SoilDataObj = Soil::get_soil_obj<Soil::soils_na>();//std::make_unique<Soil::soils_na>();

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        
    }
};

void freeze_thaw_soils::run(mesh_elem& face)
{

    auto& d = face->make_module_data<freeze_thaw_soils::data>(ID);
    
    

};


