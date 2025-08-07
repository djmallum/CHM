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
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<surface_temperature::data>(ID,face,global_param,cfg);
        
        d.mean_temperature.bind_target((*face)["air_temperature"_s]);
        d.set_face(face);
        d.set_global(global_param_);
    };
};

void surface_temperature::run(mesh_elem& face)
{
    auto& d = face->get_module_data<surface_temperature::data>(ID);

    double swe = (*face)["swe"_s];

    if (swe == 0.0)
        RCC.execute(d);
    else
        luce_tarboton.execute(d);

    d.mean_temperature.accumulate(is_new_day(),steps_per_day());

    d.set_outputs_to_face();
};

double& surface_temperature::data::air_temperature() const
{
    update_field(cache_->air_temperature,
            [this]() { return (*face)["air_temperature"_s]; }
            );

    return cache_->air_temperature;
};

double& surface_temperature::data::thaw_front_depth() const
{
    update_field(cache_->thaw_front_depth,
            [this]() { return (*face)["thaw_front_depth"_s]; }
            );

    return cache_->thaw_front_depth;
};

double& surface_temperature::data::snow_depth() const
{
    update_field(cache_->snow_depth,
            [this]() { return (*face)["snow_depth"_s]; }
            );

    return cache_->snow_depth;
};

double& surface_temperature::data::snow_density() const
{
    update_field(cache_->snow_density,
            [this]() { return (*face)["snow_density"_s]; }
            );

    return cache_->snow_density;
};

double& surface_temperature::data::ground_heat_flux() const
{
    update_field(cache_->ground_heat_flux,
            [this]() { return (*face)["ground_heat_flux"_s]; }
            );

    return cache_->ground_heat_flux;
};

double& surface_temperature::data::daily_mean_temperature() const
{
    return mean_temperature.get_last_mean();
};

double& surface_temperature::data::net_radiation() const
{
    update_field(cache_->net_radiation,
            [this]() { return (*face)["net_radiation"_s]; }
            );

    return cache_->net_radiation;
}

void surface_temperature::data::surface_temperature(const double& value)
{
    set_output(cache_->surface_temperature,value);
};

void surface_temperature::data::snow_thermal_conductivity(const double& value)
{
    set_output(cache_->snow_thermal_conductivity, value);
};

void surface_temperature::data::set_outputs_to_face()
{
	(*face)["surface_temperature"_s] = surface_temperature;
	
	(*face)["snow_thermal_conductivity"_s] = snow_thermal_conductivity;
	
	reset_cache();
};

bool surface_temperature::is_new_day()
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

int surface_temperature::steps_per_day()
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
    
	return S / D;
};
