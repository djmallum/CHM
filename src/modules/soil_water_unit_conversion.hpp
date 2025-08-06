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
#include "data_base.hpp"
#include "Soil.h"
#include "soil_moisture_converter.hpp"
#include <limits>

/*
 * Documentation here
 */ 
class soil_water_unit_converter : public module_base
{
REGISTER_MODULE_HPP(soil_water_unit_converter)
public:
    soil_water_unit_converter(config_file cfg);

    ~soil_water_unit_converter() {};

    void init(mesh& domain);
    void run(mesh_elem& face);

    struct Cache : public cache_base
    {
        double soil_storage = std::numeric_limits<double>::quiet_NaN();
        double volumetric_moisture_content = std::numeric_limits<double>::quiet_NaN();
        double saturation = std::numeric_limits<double>::quiet_NaN();
    };

    class data : public face_info, public data_base<Cache>
    {
        static const inline Soil::soils_na& SoilDataObj = 
            Soil::get_soil_obj<const Soil::soils_na>();
        std::optional<std::string> soil_type;
        std::optional<const double> soil_storage_max_;
        std::string& get_soil_type();

      public:
        explicit data(mesh_elem& face_in, const boost::shared_ptr<global> param, const config_file& cfg) 
            : data_base(face_in,param,cfg) {};

        bool storage_is_total_moisture();
        double fractional_cutoff();
        double soil_storage() const;
        const double soil_storage_max();
        double porosity();
        void volumetric_moisture_content(const double out) const;
        double volumetric_moisture_content() const;
        void saturation(const double out) const;
        double saturation() const;
        void reset_local_cache();
    };
    
private:

    soil_moisture_converter<data> converter;
};

