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
#include "TPSpline.hpp"
#include <cmath>
#include "Soil.h"
#include "soil_two_layer.hpp"
#include "soil_ET.hpp"
#include "soil_DTO.hpp"

/**
 * \ingroup modules infil soil_module exp
 * @{
 * \class soil_module *
 * 
 * Organizes soil process models (separate classes)
 *      - estimating freeze/thaw fronts using the XG algorithm: XG-freeze_thaw
 *      - layer per-time-step outflow: k_estimate
 *      - two layer soil model: soil_two_layer
 *      - soil Evapotranspiration: soil_ET
 *
 *      Each are found in soil_submodules/
 *
 * **Depends:**
 * - Snow water equivalent "swe" [mm]
 * - thaw front in soil "thaw_front_depth" [mm]
 * - freeze front in soil "freeze_front_depth" [mm]
 * - potential evapotranspiration "potential_ET" [mm]
 * - infiltrated moisture "inf" [mm]
 * - runoff/excess from infiltration scheme "runoff" [mm]
 * - leftover from routing between triangles "routing_residual" [mm]             
 *
 * **Provides:** TODO match provides with actual .cpp file
 * - Infiltration "inf" [\f$mm \cdot dt^{-1}\f$]
 * - Total infiltration "total_inf" [mm]
 * - Total infiltration excess "total_excess" [mm]
 * - Total runoff "runoff" [mm]
 * - Total soil storage "soil_storage"
 * - Potential infiltration "potential_inf"
 * - Opportunity time for infiltration to occur "opportunity_time"
 * - Available storage for water of the soil "available_storage"
 *
 * \rst
 * .. note:: TODO add any notes
 *    Has hardcoded soil parameters that need to be read from the mesh parameters.
 *
 * \endrst
 *
 * **References:** TODO add any references
 * - Gray, D., Toth, B., Zhao, L., Pomeroy, J., Granger, R. (2001). Estimating areal snowmelt infiltration into frozen soils
 * Hydrological Processes  15(16), 3095-3111. https://dx.doi.org/10.1002/hyp.320
 * @}
 */

class soil_module : public module_base
{
REGISTER_MODULE_HPP(soil_module)
public:
    soil_module(config_file cfg);

    ~soil_module();

    void run(mesh_elem& face);
    void init(mesh& domain);

    class data : public face_info, public two_layer_DTO, public soil_ET_DTO
    {
    public:
        std::unique_ptr<soil_base> soil_layers;
        std::unique_ptr<soil_base> ET;
        std::shared_ptr<mesh_elem> my_face;

        // overridden
        bool is_lake(soil_ET_DTO& DTO) override;
        
        struct my_module
        {
            soil_module& my_soil;

            my_module(soil_module& _my_soil) : my_soil(_my_soil) {};
        };
        std::shared_ptr<soil_module> local_module;
        //my_module* local_module;
    private:


    };

    void set_local_module(soil_module::data& d);


private:

    void get_soil_inputs(mesh_elem& face, data& d);
    void set_soil_outputs(mesh_elem& face, data& d);
    void set_soil_params(soil_module::data& d);
    void set_ET_params(soil_module::data& d);
    void initial_soil_conditions(soil_module::data& d);

    

};
