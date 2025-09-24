#include "soil_water_unit_conversion.hpp"
#include <stdexcept>

soil_water_unit_converter::soil_water_unit_converter(config_file cfg)
    : module_base("soil_water_unit_converter", parallel::data, cfg)
{
    depends("soil_storage");

    provides("soil_volumetric_content");
    provides("soil_saturation");
    provides("soil_saturation_at_freeze");
};

double soil_water_unit_converter::data::soil_storage() 
{
    update_field([this]() -> auto& { return cache_->soil_storage;},
           [this]() { return (*face)["soil_storage"_s]; } );

    return cache_->soil_storage;
};

void soil_water_unit_converter::init(mesh& domain)
{
    for (size_t i = 0; i < domain->size_local_faces(); i++)
    {
        // Doing nothing except creating the data objects
        auto face = domain->face(i);
        auto& d = face->make_module_data<soil_water_unit_converter::data>(ID, face, global_param, cfg);
    };
};

void soil_water_unit_converter::run(mesh_elem& face)
{
    // Constructor is empty because only fetching...
    auto& d = face->get_module_data<soil_water_unit_converter::data>(ID);

    converter.execute(d);

    d.reset_local_cache();
};

bool soil_water_unit_converter::data::storage_is_total_moisture() { return false; };

double soil_water_unit_converter::data::fractional_cutoff()
{
    return SoilDataObj.wilt_point(get_soil_type());
};

const double soil_water_unit_converter::data::soil_storage_max()
{
    if (!soil_storage_max_)
    {
        soil_storage_max_.emplace(face->soil_attribute<double>("soil_storage_max"));
    }
    return *soil_storage_max_;
};

double soil_water_unit_converter::data::porosity()
{
    return SoilDataObj.porosity(get_soil_type());
};

void soil_water_unit_converter::data::volumetric_moisture_content(const double out)
{
    set_output([this]() -> auto& { return cache_->volumetric_moisture_content; },out); 
};

double soil_water_unit_converter::data::volumetric_moisture_content()
{
    auto static err = 
        []() {  throw std::runtime_error("volumetric_moisture_content accessed before being set");
                return 1.0; };
    
    update_field([this]() -> auto& { return cache_->volumetric_moisture_content; },err);

    return cache_->volumetric_moisture_content;
};

void soil_water_unit_converter::data::saturation(const double out) 
{
    set_output([this]() -> auto& { return cache_->saturation; },out); 
};

double soil_water_unit_converter::data::saturation()
{
    auto static err = 
       []() { throw std::runtime_error("saturation accessed before being set");
              return 1.0;   };

    update_field([this]() -> auto& { return cache_->saturation; },err);

    return cache_->saturation;
};

std::string& soil_water_unit_converter::data::get_soil_type()
{
    if (!soil_type)
    {
        soil_type = cfg_.get<std::string>("soil_type");
    }

    return *soil_type;
};

void soil_water_unit_converter::data::reset_local_cache()
{
    (*face)["soil_volumetric_content"_s] = volumetric_moisture_content();
    (*face)["soil_saturation"_s] = saturation();
    reset_cache();
};

