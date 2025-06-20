double surface_temperature::data::air_temperature()
{
	if (!is_initialized)
	{
		return (*face)["air_temperature"_s];
	};

	if (!Cache && face)
	{
		Cache = make_unique<TempCache>();
		Cache.air_temperature = 
			(*face)["air_temperature"_s];
	}
	else if (!face)
	{
		CHM_THROW_EXCEPTION(module_error,
				"pointer to face uninstantiated in surface_temperature::data");
	}

	return Cache.air_temperature;
};

double surface_temperature::data::air_temperature()
{
	if (!Cache && face)
	{
		Cache = make_unique<TempCache>();
		Cache.air_temperature = 
			(*face)["air_temperature"_s];
	}
	else if (!face)
	{
		CHM_THROW_EXCEPTION(module_error,
				"pointer to face uninstantiated in surface_temperature::data");
	}

	return Cache.air_temperature;
};

double surface_temperature::data::air_temperature()
{
	if (!Cache && face)
	{
		Cache = make_unique<TempCache>();
		Cache.air_temperature = 
			(*face)["air_temperature"_s];
	}
	else if (!face)
	{
		CHM_THROW_EXCEPTION(module_error,
				"pointer to face uninstantiated in surface_temperature::data");
	}

	return Cache.air_temperature;
};

double surface_temperature::data::air_temperature()
{
	if (!Cache && face)
	{
		Cache = make_unique<TempCache>();
		Cache.air_temperature = 
			(*face)["air_temperature"_s];
	}
	else if (!face)
	{
		CHM_THROW_EXCEPTION(module_error,
				"pointer to face uninstantiated in surface_temperature::data");
	}

	return Cache.air_temperature;
};
