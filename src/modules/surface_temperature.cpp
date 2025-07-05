#include "surface_temperature.hpp"

REGISTER_MODULE_CPP(surface_temperature);

surface_temperature::surface_temperature(config_file cfg) : module_base("surface_temperature", parallel::data,cfg)
{
    depends("swe");
    depends("air_temperature");
    depends("snow_depth");
    depends("ground_heat_flux");
    depends("snow_density");

    provides("surface_temperature");
    provides("snow_thermal_conductivity");
    
};
void surface_temperature::init(mesh& domain)
{
	// enforces all access to data via cache is through face rather than the cache.
	// automatically ends this option when it goes out of scope
	auto guard = init_phase();

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<surface_temperature::data>(ID);

        d.RCC = std::make_unique<RCC_surface_temperature<surface_temperature::API>>(d.api);

        d.luce_tarboton = std::make_unique<luce_tarboton_surface_temperature<surface_temperature::API>>(d.api);
        
        d.api.set_face(face);
        d.api.set_global(global_param_);
    };
};

void surface_temperature::run(mesh_elem& face)
{
    auto& d = face->make_module_data<surface_temperature::data>(ID);

    double swe = (*face)["swe"_s];

    if (swe == 0.0)
        d.RCC->execute();
    else
        d.luce_tarboton->execute();

    d.api.set_outputs_to_face();
};

double& surface_temperature::API::air_temperature()
{
	return Cache_.get_field(
			[this]() -> double& { return (*face)["air_temperature"_s]; },
			[this]() -> double& { return Cache_.air_temperature; }
			);
};

double& surface_temperature::API::snow_depth()
{
	return Cache_.get_field(
			[this]() -> double& { return (*face)["snow_depth"_s]; },
			[this]() -> double& { return Cache_.snow_depth; }
			);
};

double& surface_temperature::API::ground_heat_flux()
{
	return Cache_.get_field(
			[this]() -> double& { return (*face)["ground_heat_flux"_s]; },
			[this]() -> double& { return Cache_.ground_heat_flux; }
			);
};

void surface_temperature::API::surface_temperature(const double& value)
{
	Cache_.set_field(
			[this]() -> double& { return Cache_.surface_temperature; }, value);
};

void surface_temperature::API::snow_thermal_conductivity(const double& value)
{
	Cache_.set_field(
			[this]() -> double& { return Cache_.snow_thermal_conductivity; }, value);
};

void surface_temperature::API::set_outputs_to_face()
{
    TempCache& c = Cache_.get_cache();
	(*face)["surface_temperature"_s] = c.surface_temperature;
	(*face)["snow_thermal_conductivity"_s] = c.snow_thermal_conductivity;
	c.reset();
};

bool surface_temperature::API::is_new_day()
{
    int td = global_param->posix_time().time_of_day().total_seconds();
    int time_to_midnight = 86400 - td;
    if (td >= 0 && td < global_param->dt()) //(time_to_midnight >= global_param->dt())
    {
        return true;
    }
    else
        return false;
};
