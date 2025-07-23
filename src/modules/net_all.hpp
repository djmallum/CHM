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

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "data_base.hpp"
#include "net_radiation.hpp"
#include "TPSpline.hpp"
#include "Atmosphere.h"
#include <cmath>


/**
 * \ingroup modules infil soils exp
 * @{
 * \class net_all
 *
 *
 * @}
 */
class net_all : public module_base
{
REGISTER_MODULE_HPP(net_all)
public:
	explicit net_all(config_file cfg);

    ~net_all() {};

    void run(mesh_elem &face) override;
    void init(mesh& domain) override;
    
private:

    struct Cache : public cache_base
    {
        //inputs
        double max_sun_hours = std::numeric_limits<double>::quiet_NaN();
        double air_temperature = std::numeric_limits<double>::quiet_NaN();
        double vapour_pressure = std::numeric_limits<double>::quiet_NaN();
        double actual_sun_hours = std::numeric_limits<double>::quiet_NaN();
        double bright_sun_ratio = std::numeric_limits<double>::quiet_NaN();
        double direct_short_wave_clear = std::numeric_limits<double>::quiet_NaN();
        double diffuse_short_wave_clear = std::numeric_limits<double>::quiet_NaN();
        double albedo = std::numeric_limits<double>::quiet_NaN();

        //outputs
        double net_all_wave = 0.0;	
    };

    class data : public face_info, public data_base<Cache>
    {
    public:
        double& max_sun_hours() const;
        double& air_temperature() const;
        double& vapour_pressure() const;
        double& actual_sun_hours() const;
        double& bright_sun_ratio() const;
        double& direct_short_wave_clear() const;
        double& diffuse_short_wave_clear() const;
        double& albedo() const;

        void net_all_wave(const double& out);

        void set_outputs_to_face();
    };

private:

	net_radiation<data> net_rad;
};
