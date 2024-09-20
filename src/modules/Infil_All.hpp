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
#include "physics/soil.hpp"

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
class Infil_All : public module_base
{
REGISTER_MODULE_HPP(Infil_All)
public:
    Infil_All(config_file cfg);

    ~Infil_All();

    void run(mesh_elem &face);
    void init(mesh& domain);

    class data : public face_info
    {
    public:
        struct tempvars
        {
            double intensity;
            double soil_storage_deficit;
            double capillary_suction;
            double initial_rate;
            double initial_storage;
            double final_rate;
            double final_storage;
            double pond;
            double storage_at_ponding;
            double time_to_ponding;
        };
        double total_inf;
        double total_excess;
        double total_snowinf;
        double total_meltexcess;
        double total_rain_on_snow;
        
        // Crack
        bool frozen;
        double index;
        double max_major_per_melt;
        double init_SWE;
        unsigned int major_melt_count;
        
        // GreenAmpt
        double soil_storage;
        std::unique_ptr<tempvars> GA_temp{nullptr};
        
    };

private:
    
    // Soil Data
    Soil::SoilData SoilDataObj;

    // Crack
    double major;
    double min_swe_to_freeze;
    unsigned int infDays;
    bool AllowPriorInf;   
    double lenstemp;

    // General, thawed soil
    enum ThawOptions { AYERS, GREENAMPT};
    unsigned int ThawType;

    // Ayers
    unsigned int texture;
    unsigned int groundcover;
    

    // GreenAmpt
    unsigned int soil_type;
    double max_soil_storage;
    double soil_depth;
    double porosity;
    enum GATable {PSI, KSAT, WILT, FCAP, PORG, PORE, AIENT, PORESZ, AVAIL}; // Used for mapping the soil table, PSI and KSAT are used, the others are unused but may but used in the future or other modules.    
    double ksaturated;
    enum GAVars {TOTINF, RATEINF, SUCTION, THETA};
    //double soilproperties[][9];
    //double textureproperties[][6];
     

    // General Functions
    void Increment_Totals(data &d, double &runoff, double &melt_runoff, double &inf, double &snowinf, double &rain_on_snow);

    // Crack Functions
    void Calc_Index(data &d, double &swe, double &theta);
    double Calc_Actual_Inf(data &d, double &melt);
    void Check_for_ice_lens(data &d, double &t); 

    // Green-Ampt Functions
    double convert_to_rate_hourly(double &rainfall); 
    bool is_space_in_dry_soil(double &moist, double &max, double &rainfall); 
    void Initialize_GA_Variables(data &d); 
    void initialize_ponding_vars(std::unique_ptr<data::tempvars> &GA); 
    void find_final_storage(std::unique_ptr<data::tempvars> &GA, \
        double &initial_storage, double &dt); 
    double calc_GA_infiltration_rate(std::unique_ptr<data::tempvars> &GA, double &F);





};
