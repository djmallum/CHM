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
	
		// Inputs
		mutable double air_temperature = std::numeric_limits<double>::quiet_nan();
		mutable double thaw_front_depth = std::numeric_limits<double>::quiet_nan();
		mutable double snow_depth = std::numeric_limits<double>::quiet_nan();
		mutable double snow_density = std::numeric_limits<double>::quiet_nan();
		mutable double ground_heat_flux = std::numeric_limits<double>::quiet_nan();
		mutable double net_radiation = std::numeric_limits<double>::quiet_nan();
		mutable double daily_mean_temperature = std::numeric_limits<double>::quiet_nan();

		// Outputs
		double surface_temperature = 0.0; 
		double snow_thermal_conductivity = 0.0;
    
    public:
        //inputs
        double& air_temperature() const;
	    double& thaw_front_depth() const;
		double& snow_depth() const;
	    double& snow_density() const;
        double& ground_heat_flux() const;
        double& daily_mean_temperature() const;
	    double& net_radiation() const;

        //outputs
		void surface_temperature(const double& in);
		void snow_thermal_conductivity(const double& in); // just if snow-covered
       
        // used by module    
        void set_outputs_to_face();   
        void set_face(mesh_elem& face_in)
        { face = *face_in; };
        void set_global(global* global_param_)
        { global_param = global_param_; };
		void reset_cache();

        // daily accumulator
        bool is_new_day();
        double steps_per_day();
    };

	class data : public face_info
	{
	public:
        std::unique_ptr<RCC_surface_temperature<API>> RCC;
        
        API api;
	};

    auto init_phase()
    {
        return cache_handler<API::TempCache>::scoped_init(); 
    };

private:
    luce_tarboton_surface_temperature<API> luce_tarboton;

};

