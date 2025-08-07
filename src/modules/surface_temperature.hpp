#pragma once

#include "daily_accumulator.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "data_base.hpp"
#include "RCC_surface_temperature.hpp"
#include "luce_tarboton_surface_temperature.hpp"

class surface_temperature : public module_base
{
REGISTER_MODULE_HPP(surface_temperature);
public:
	surface_temperature(config_file cfg);
	~surface_temperature() = default;

	void init(mesh& domain) override;
	void run(mesh_elem& face) override;

private:
    struct Cache :  public cache_base
    {
        // Inputs
		mutable double air_temperature = std::numeric_limits<double>::quiet_NaN();
		mutable double thaw_front_depth = std::numeric_limits<double>::quiet_NaN();
		mutable double snow_depth = std::numeric_limits<double>::quiet_NaN();
		mutable double snow_density = std::numeric_limits<double>::quiet_NaN();
		mutable double ground_heat_flux = std::numeric_limits<double>::quiet_NaN();
		mutable double net_radiation = std::numeric_limits<double>::quiet_NaN();
		mutable double daily_mean_temperature = std::numeric_limits<double>::quiet_NaN();

		// Outputs
		double surface_temperature = 0.0; 
		double snow_thermal_conductivity = 0.0;
    };

	class data : public face_info, public data_base<Cache>
	{
	public:
          explicit data(mesh_elem& face_in, boost::shared_ptr<global> param, config_file cfg);
          double& air_temperature() const;
          double& thaw_front_depth() const;
          double& snow_depth() const;
          double& snow_density() const;
          double& ground_heat_flux() const;
          const double daily_mean_temperature() const;
          double& net_radiation() const;

          double thaw_front_depth_last = 0.0;
          // outputs
          void surface_temperature(const double& in);
          void snow_thermal_conductivity(const double& in); // just if snow-covered

          // used by module
          void set_outputs_to_face();
          
          daily_accumulator mean_temperature;
	};

private:
    luce_tarboton_surface_temperature<data> luce_tarboton;
    RCC_surface_temperature<data> RCC;

    // daily accumulator
    bool is_new_day();
    int steps_per_day();

};

