#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "RCC_surface_temperature.hpp"
#include "luce_tarboton_surface_temperature.hpp"
#include "cache_handler.hpp"
#include <memory>

class surface_temperature : public module_base
{
REGISTER_MODULE_HPP(surface_temperature);
public:
	surface_temperature(config_file cfg);
	~surface_temperature();

	void init(mesh& domain);
	void run(mesh_elem& face);

    class API
    {
        mesh_elem face = nullptr;
        global* global_param = nullptr;

    public:
        struct TempCache
        {
            double air_temperature = 0.0;
            double thaw_front_depth = 0.0;
            double snow_depth = 0.0;
            double snow_density = 0.0;
            double ground_heat_flux = 0.0;
            double net_radiation = 0.0;

            double surface_temperature = 0.0;
            double snow_thermal_conductivity = 0.0;
        };
    private:
		cache_handler<TempCache> Cache_;
    
    public:
        //inputs
        double& air_temperature();
	    double& thaw_front_depth();
		double& snow_depth();
	    double& snow_density();
        double& ground_heat_flux();
	    double& net_radiation();

        //outputs
		void surface_temperature(const double& in);
		void snow_thermal_conductivity(const double& in); // just if snow-covered
        
        // used by module    
        void set_outputs_to_face();   
        void set_face(mesh_elem& face_in)
        { face = *face_in; };
        void set_global(global* global_param_)
        { global_param = global_param_; };

        // daily accumulator
        bool is_new_day();
        double steps_per_day();
    };

	class data : public face_info
	{
	public:
        std::unique_ptr<RCC_surface_temperature<API>> RCC;
        std::unique_ptr<luce_tarboton_surface_temperature<API>> luce_tarboton;

        API api;
	};

    auto init_phase()
    {
        return cache_handler<API::TempCache>::scoped_init(); 
    };

};

