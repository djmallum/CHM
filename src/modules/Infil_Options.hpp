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
#include "Crack.hpp"
#include "Ayers_infiltration.hpp"
#include "data_base.hpp"
#include <string_view>
/**
 * \ingroup modules infil soils exp
 * @{
 * \class Infil_All
 *
 * TODO some of this is correct, but likely will be outdated eventually.
 * Estimates areal snowmelt infiltration into frozen soils for:
 *    a) Restricted -  Water entry impeded by surface conditions
 *    b) Limited - Capiliary flow dominates and water flow influenced by soil physical properties
 *    c) Unlimited - Gravity flow dominates
 *
 * **Depends:**
 * - Snow water equivalent "swe" [mm]
 * - Snow melt for interval "snowmelt_int" [\f$mm \cdot dt^{-1}\f$]
 *
 * **Provides:**
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
 * .. note::
 *    Has hardcoded soil parameters that need to be read from the mesh parameters.
 *
 * \endrst
 *
 * **References:**
 * - Gray, D., Toth, B., Zhao, L., Pomeroy, J., Granger, R. (2001). Estimating areal snowmelt infiltration into frozen soils
 * Hydrological Processes  15(16), 3095-3111. https://dx.doi.org/10.1002/hyp.320
 * @}
 */
class Infil_Options : public module_base
{
REGISTER_MODULE_HPP(Infil_Options)

    enum class Status
    {
        FROZEN,
        THAWED,
        TO_FROZEN,
        TO_THAWED
    };
public:
    Infil_Options(config_file cfg);

    ~Infil_Options();

    void run(mesh_elem &face);
    void init(mesh& domain);
    
    struct Cache : public cache_base
    {
        // Shared
        double snowmelt = default_value<double>();
        double rainfall = default_value<double>();
        // Crack
        double swe = default_value<double>();
        double air_temperature = default_value<double>();

        double infiltrated = 0.0;
        double runoff = 0.0;
        double snow_infiltrated = 0.0;
        double melt_runoff = 0.0; 
        double rain_on_snow = 0.0;
    };

    class data : public data_base<Cache>,face_info
    {
        friend class Infil_Options;
        Crack::State state;
        std::string _texture;
        std::string _ground_cover;
        bool saturation_set;
        void set_outputs() const;
    public:
        Crack::State& get_state();
        bool is_newday() const;
        double snowmelt();
        double rainfall();
        double swe();
        double air_temperature();
        std::string_view texture();
        std::string_view ground_cover();

        void infiltrated(const double);
        void runoff(const double);
        void snow_infiltrated(const double);
        void melt_runoff(const double);
        void rain_on_snow(const double);
        
        struct CrackCHMDetails
        {
            Status status = Status::THAWED;
            bool end_freeze_tomorrow = false;
        } crack_details;
    };

private:
     
	bool is_new_day(void);

    Crack::Model<data> crack;
    Ayers::Model<data> ayers;
    
    Status get_status(data& d,const Status Old);

    struct DomainConstantsCrack
    {
        double min_swe_to_freeze;
        size_t day_of_year_to_freeze;
    };

    std::unique_ptr<const DomainConstantsCrack> constants_crack;
};

