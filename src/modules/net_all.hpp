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
#include "net_radiation.hpp"
#include "TPSpline.hpp"
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

    class data : public face_info
    {
        
        struct Cache
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

            size_t last_timestep = -1;

            bool is_stale(int tn)
            { return last_timestep != tn; };
        }; 
        mutable std::optional<Cache> cache_;

        void init_cache() const
        {
            size_t current = global_param->timestep_counter;
            if (!cache_ || cache_->is_stale(current))
            {
                cache_.emplace();
                cache_->last_timestep = current;
            };
        };

        template<typename Fetch>
        void update_field(double& value, Fetch fetch) const;

        mesh_elem face{nullptr};
        boost::shared_ptr<global> global_param;
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
        void set_face(mesh_elem& face_in)
        { face = face_in; };
        void set_global(boost::shared_ptr<global> param)
        { global_param = param; };
        void reset_cache();
    };

private:

	net_radiation net_rad;
};
