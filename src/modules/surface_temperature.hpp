#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "surface_temperature_submodules/RCC_surface_temperature.hpp"
#include "surface_temperature_submodules/luce_tarboton_surface_temperature.hpp"

class surface_temperature : public module_base
{
REGISTER_MODULE_HPP(surface_temperature);
public:
	surface_temperature(config_file cfg);
	~surface_temperature();

	void init(mesh& domain);
	void run(mesh_elem& face);

	class data : public face_info
	{
	public:
		double air_temperature();
		
		double snow_depth();
		double ground_heat_flux();
		
		void surface_temperature(const double& in);
		void snow_thermal_conductivity(const double& in);
		
		void done_init()
		{ is_inialized = true; };
	private:
		mesh_elem* face = nullptr;
		std::unqiue_ptr<TempCache> Cache = nullptr;

		struct TempCache
		{
			double air_temperature = 0.0;
			double snow_depth = 0.0;
			double snow_density = 0.0;
			double ground_heat_flux = 0.0;
			double surface_temperature = 0.0;
			double snow_thermal_conductivity = 0.0;
		};


		static bool is_initialized = false;
	};

};

