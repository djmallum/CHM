//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#pragma once

#include "triangulation.hpp"
#include "module_base.hpp"
#include "Soil.h"
#include <cstdlib>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>
#include "EvapotranspirationModels/evapbase.hpp"
#include "EvapotranspirationModels/PenmanMonteith.hpp"
#include "EvapotranspirationModels/PriestleyTaylor.hpp"
#include "data_base.hpp"
#include "net_all_bad_lake.hpp"

/**
 * \ingroup modules exp evap
 * @{
 * \class Evapotranspiration_All

 * Calculates evapo-transpiration via Penman-Monteith.
 *
 * Not currently maintained.
 *
 * Depends:
 * - Incoming shortwave radaition "iswr" [\f$ W \cdot m^{-2} \f$ ]
 * - Incoming longwave radiation "ilwr" [\f$ W \cdot m^{-2} \f$ ]
 * - Relative humidity "rh" [%]
 * - Windspeed at 2m "U_2m_above_srf" [\f$ m \cdot s^{-1} \f$ ]
 *
 * Provides:
 * - Evapotranspiration "ET" [\f$ mm \cdot dt^{-1} \f$ ]
 *
 * @}
 */
class Evapotranspiration_All : public module_base
{
REGISTER_MODULE_HPP(Evapotranspiration_All);
public:
    Evapotranspiration_All(config_file cfg);
    ~Evapotranspiration_All();
    void init(mesh& domain);
    void run(mesh_elem& face);

    class Cache : public cache_base
    {
    public:
        //Inputs
        Input<double> incoming_short_wave;
        Input<double> wind_measurement_height;
        Input<double> d;
        Input<double> Z0;
        Input<double> kappa;
        Input<double> wind_speed;
        Input<double> stomatal_resistance_min;
        Input<bool> has_vegetation;
        Input<double> Veg_height;
        Input<double> leaf_area_index_max;
        Input<double> short_wave_in;
        Input<double> saturated_vapour_pressure;
        Input<double> vapour_pressure;
        Input<double> air_entry_tension;
        Input<double> porosity;
        Input<double> volumetric_moisture_content;
        Input<double> pore_size_dist;
        Input<double> air_temperature;
        Input<double> delta;
        Input<double> Q_net;
        Input<double> Q_g;
        Input<double> air_density;
        Input<double> heat_capacity_air;
        Input<double> gamma;
        Input<double> stomatal_resistance;
        Input<double> lambda;
        Input<int> s_per_time_step;
        Input<double> P_atm;

        // Output variables
        Output<double> ET;
        Output<double> net_all_wave;
    };

    class data : public face_info, public data_base<Cache>
    {
    public:
        std::unique_ptr<evapT_base> MyPenmanMonteith; 
        std::unique_ptr<evapT_base> MyPriestleyTaylor;
        
        double LAI;
        double LAImax;
        double vegetation_height;
        double soil_depth;

        double albedo();
        double incoming_short_wave();
        double net_all_wave();
        double wind_measurement_height() const;
        double d();
        double Z0() const;
        double kappa() const;
        double wind_speed() const;
        double stomatal_resistance_min() const;
        bool has_vegetation() const;
        double Veg_height() const;
        double leaf_area_index_max() const;
        double short_wave_in() const;
        double saturated_vapour_pressure() const;
        double vapour_pressure() const;
        double air_entry_tension() const;
        double porosity() const;
        double volumetric_moisture_content() const;
        double pore_size_dist() const;
        double air_temperature() const;
        double delta() const;
        double Q_net() const;
        double Q_g() const;
        double air_density() const;
        double heat_capacity_air() const;
        double gamma() const;
        double stomatal_resistance() const;
        double lambda() const;
        int s_per_time_step() const;

        // Output setters
        void stomatal_resistance(const double value);
        void ET(const double value);
        void net_all_wave(const double val);

        data(const mesh_elem& face_in, const boost::shared_ptr<global> param,
                const config_file cfg) : data_base<Cache>(face_in,param,cfg) {};
        ~data() {}; 
    };

private:

    const Soil::soils_na& SoilDataObj = Soil::get_soil_obj<Soil::soils_na>();
 
    // PriestleyTaylor
    double alpha;
    // TODO add PT methods here 

    void init_PriestleyTaylor(Evapotranspiration_All::data& d, const double alpha,const int dt);
    PT_vars set_PriestleyTaylor_vars(mesh_elem& face,data& d);
   
    // PenmanMonteith
    double stomatal_resistance_min;
    double Frac_to_ground;
    
    void init_PenmanMonteith(Evapotranspiration_All::data& d, mesh_elem& face, const double wind_height, 
            const double stomatal_resistance_min, const double Frac_to_ground);
    
    PM_vars set_PenmanMonteith_vars(mesh_elem& face, const double t, 
            const double saturated_vapour_pressure, const double vapour_pressure, data& d);
    
    const double get_dt();

    net_all_bad_lake<data> net_radiation; 
};

