#include "net_all.hpp"
net_all::net_all(config_file cfg) : module_base("net_all", parallel::data, cfg)
{
	depends("t");
	depends("rh");
	depends("iswr_direct");
	depends("iswr_diffuse");

    provides("net_all_wave");

};

void net_all::init(mesh& domain)
{
	for (size_t i = 0; i < domain->size_local_faces(); i++)
	{
		auto face = domain->face(i);
		auto& d = face->make_module_data<net_all::data>(ID);

		d.set_face(face);
		d.set_global(global_param);
	}
};

void net_all::run(mesh_elem& face)
{
	auto& d = face->make_module_data<net_all::data>(ID);

	net_rad.execute(d);

	d.set_outputs_to_face();
};

double& net_all::data::air_temperature() const
{
    update_field(
            cache_->air_temperature,
            [this]() { return (*face)["t"_s]; });

	if (cache_->air_temperature > 150.0)
		CHM_THROW_EXCEPTION(module_error, "net_all: Airtemperature too large, likely unphysical or Kelvin");

	return cache_->air_temperature;
};

double& net_all::data::vapour_pressure() const
{
    update_field(
            cache_->vapour_pressure,
            [this]() {
            const static double CtoKelvin = 273.15;
            double relative_humidity = (*face)["rh"_s];
            return relative_humidity * Atmosphere::saturatedVapourPressure(air_temperature() + CtoKelvin);
            }
            );

	return cache_->vapour_pressure;
};

double& net_all::data::max_sun_hours() const
{
    update_field(
            cache_->max_sun_hours,
            [this]() { return (*face)["max_sun_hours"_s]; }
            );

	return cache_->max_sun_hours;;
};

double& net_all::data::actual_sun_hours() const
{
    update_field(
            cache_->actual_sun_hours,
            [this]() { return (*face)["actual_sun_hours"_s]; }
            );

	return cache_->actual_sun_hours;
};

double& net_all::data::bright_sun_ratio() const
{
    update_field(
            cache_->bright_sun_ratio,
            [this]() { return actual_sun_hours() / max_sun_hours(); }
            );
             
	return cache_->bright_sun_ratio;
};

double& net_all::data::direct_short_wave_clear() const
{
    update_field(
            cache_->direct_short_wave_clear,
            [this]() { return (*face)["direct_short_wave_clear"_s]; }
            );

	return cache_->direct_short_wave_clear;
};

double& net_all::data::diffuse_short_wave_clear() const
{
    update_field(
            cache_->diffuse_short_wave_clear,
            [this]() { return (*face)["diffuse_short_wave_clear"_s]; }
            );

	return cache_->diffuse_short_wave_clear;
};

double& net_all::data::albedo() const
{
    update_field(
            cache_->albedo,
            [this]() { return (*face)["albedo"_s]; }
            );

	return cache_->albedo;
};

void net_all::data::net_all_wave(const double& out)
{
    set_output(cache_->net_all_wave, out);
};

void net_all::data::set_outputs_to_face()
{
	(*face)["net_all_wave"_s] = cache_->net_all_wave;

    reset_cache();    
};
