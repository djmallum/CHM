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

double& surface_temperature::API::air_temperature() const
{
	if (std::isnan(air_temperature)) 
		air_temperature = (*face)["air_temperature"_s];

	return air_temperature;
};

double& surface_temperature::API::thaw_front_depth() const
{
	if (std::isnan(thaw_front_depth)) 
		thaw_front_depth = (*face)["thaw_front_depth"_s];

	return thaw_front_depth;
};

double& surface_temperature::API::snow_depth() const
{
	if (std::isnan(snow_depth)) 
		snow_depth = (*face)["snow_depth"_s];

	return snow_depth;
};

double& surface_temperature::API::snow_density() const
{
	if (std::isnan(snow_density)) 
		snow_density = (*face)["snow_density"_s];

	return snow_density;
};

double& surface_temperature::API::ground_heat_flux() const
{
	if (std::isnan(ground_heat_flux)) 
		ground_heat_flux = (*face)["ground_heat_flux"_s];

	return ground_heat_flux;
};

double& surface_temperature::API::net_radiation() const
{
	if (std::isnan(net_radiation)) 
		net_radiation = (*face)["net_radiation"_s];

	return net_radiation;
}


void surface_temperature::API::surface_temperature(const double& value)
{
	surface_temperature = value;
};

void surface_temperature::API::snow_thermal_conductivity(const double& value)
{
	snow_thermal_conductivity = value;
};

void surface_temperature::API::set_outputs_to_face()
{
	(*face)["surface_temperature"_s] = surface_temperature;
	
	(*face)["snow_thermal_conductivity"_s] = snow_thermal_conductivity;
	
	reset_cache();
};

void surface_temperature::API::reset_cache();
{
	air_temperature = std::numeric_limits<double>::quiet_nan();
	thaw_front_depth = std::numeric_limits<double>::quiet_nan();
	snow_depth = std::numeric_limits<double>::quiet_nan();
	snow_density = std::numeric_limits<double>::quiet_nan();
	ground_heat_flux = std::numeric_limits<double>::quiet_nan();
	net_radiation = std::numeric_limits<double>::quiet_nan();

	surface_temperature = 0.0;
	snow_thermal_conductivity = 0.0;
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

int surface_temperature::API::steps_per_day()
{
	static const auto [S,D] = [this]() ->std::pair<int, double>
	{
		const int s = 86400;
		const int d = global_param->dt();

		if (d > s)
		{
			CHM_THROW_EXCEPTION(module_error, "Surface_temerature: time step must be less than a day");
		}

		return {s,d};
	}();

	return s / d;
};
