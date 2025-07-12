net_all::net_all()
{
	// depends and provides
};

void net_all::init(mesh& domain)
{

};

void net_all::run(mesh_elem& face)
{
	auto& d = face->make_module_data<net_all::data>(ID);

	net_rad.execute(d);

	d.set_outputs_to_face();
};

double& net_all::max_sun_hours() const
{
	if (std::isnan(max_sun_hours))
		max_sun_hours = (*face)["max_sun_hours"_s];

	return max_sun_hours;
};

double& net_all::air_temperature() const
{
	if (std::isnan(air_temperature))
		air_temperature = (*face)["air_temperature"_s];

	return air_temperature
};

double& net_all::vapour_pressure() const
{
	if (std::isnan(vapour_pressure))
		vapour_pressure = (*face)["vapour_pressure"_s];

	return vapour_pressure
};

double& net_all::actual_sun_hours() const
{
	if (std::isnan(actual_sun_hours))
		actual_sun_hours = (*face)["actual_sun_hours"_s];

	return actual_sun_hours
};

double& net_all::direct_short_wave_clear() const
{
	if (std::isnan(direct_short_wave_clear))
		direct_short_wave_clear = (*face)["direct_short_wave_clear"_s];

	return direct_short_wave_clear
};

double& net_all::diffuse_short_wave_clear() const
{
	if (std::isnan(diffuse_short_wave_clear))
		diffuse_short_wave_clear = (*face)["diffuse_short_wave_clear"_s];

	return diffuse_short_wave_clear
};

double& net_all::albedo() const
{
	if (std::isnan(albedo))
		albedo = (*face)["albedo"_s];

	return albedo
};

void net_all::net_all_wave(const double& out)
{
	net_all_wave = out;
};

void net_all::set_outputs_to_face()
{
	(*face)["net_all_wave"_s] = net_all_wave;
};

void net_all::data::reset_cache()
{
	max_sun_hours = std::numeric_limits<double>::quiet_nan();
	air_temperature = std::numeric_limits<double>::quiet_nan();
	vapour_pressure = std::numeric_limits<double>::quiet_nan();
	actual_sun_hours = std::numeric_limits<double>::quiet_nan();
	direct_short_wave_clear = std::numeric_limits<double>::quiet_nan();
	diffuse_short_wave_clear = std::numeric_limits<double>::quiet_nan();
	albedo = std::numeric_limits<double>::quiet_nan();

	net_all_wave = 0.0;
}
