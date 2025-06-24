void surface_temperature::init(mesh& domain)
{
	// enforces all access to data via cache is through face rather than the cache.
	// automatically ends this option when it goes out of scope
	auto guard = cache_handler<TempCache>::scoped_init();
};

double& surface_temperature::data::air_temperature()
{
	return Cache.get_field(
			[this]() -> double& { return (*face)["air_temperature"_s]; },
			[](TempCache& c) -> double& { return c.air_temperature; }
			)
};

double& surface_temperature::data::snow_depth()
{
	return Cache.get_field(
			[this]() -> double& { return (*face)["snow_depth"_s]; },
			[](TempCache& c) -> double& { return c.snow_depth; }
			)
};

double& surface_temperature::data::ground_heat_flux()
{
	return Cache.get_field(
			[this]() -> double& { return (*face)["ground_heat_flux"_s]; },
			[](TempCache& c) -> double& { return c.ground_heat_flux; }
			)
};

void surface_temperature::data::surface_temperature(double& value)
{
	Cache.set_field(
			[](TempCache& c) -> double& { return c.surface_temperature; }, value);
};

void surface_temperature::data::snow_thermal_conductivity(double& value)
{
	Cache.set_field(
			[](TempCache& c) -> double& { return c.snow_thermal_conductivity; }, value);
};

void surface_temperature::data::set_outputs_to_face(TempCache&& c)
{
	(*face)["surface_temperature"_s] = c.surface_temperature;
	(*face)["snow_thermal_conductivity"_s] = c.snow_thermal_conductivity;
	c.reset();
};
