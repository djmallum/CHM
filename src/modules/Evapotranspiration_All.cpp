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


#include "Evapotranspiration_All.hpp"
#include "EvapotranspirationModels/evapbase.hpp"
#include "EvapotranspirationModels/PenmanMoneith.hpp"
// PT here #include "EvapotranspirationModels/PriestelyTaylor.hpp"

REGISTER_MODULE_CPP(Evapotranspiration_All);

Evapotranspiration_All::Evapotranspiration_All(config_file cfg)
        :module_base("Evapotranspiration_All", parallel::data, cfg)
{
    // TODO Constructor is not properly editted with all new inputs (see set vars function at the end)
    depends("iswr");
    depends("ilwr");
    depends("rh");
    depends("t");
    depends("U_2m_above_srf"); // 
    depends("snow_albedo"); // named inaccurate to this module, 
    depends("soil_storage");                      // but is how albedo is used in CHM as of Sept, 2024

    provides("ET");


}

void Evapotranspiration_All::init(mesh& domain)
{
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<Infil_All::data>(ID);
        
        // Consider if an if statement is necessary.

        init_PenmanMonteith(d);

        // put PristelyTaylor parameters here.
    }

}

void Evapotranspiration_All::run(mesh_elem& face)
{

    auto& d = face->get_module_data<Evapotranspiration_All::data>(ID);

    if (is_water(face))
    {
        // Do PriestlyTaylor
    }
    // else if (If wetlands)
    // {
    //     Also Do PriestlyTaylor
    // }
    else
    {
        // All members of PM_vars are references and must be set at initialization
        // Therefore we need a copy of SVP to reference
        // t is made its own copy to avoid dereferencing face for "t" twice

        double t = (*face)["t"_s];
        double SVP = Atmosphere::saturatedVapourPressure(t);  
        PM_vars my_PM_vars = set_PenmanMoneith_vars(face,t,AVP);
    
        model_output output;
    
        d.MyPenmanMonteith->CalcEvapT(my_PM_vars,output);
    }
    
    
    (*face)["ET"_s] = output.ET;
    
}

Evapotranspiration_All::~Evapotranspiration_All()
{

}

void Evapotranspiration_All::init_PenmanMonteith(Evapotransporation_All::data& d)
{
    
    const std::string soil_type = cfg.get("soil_type","sand");
   
    SoilDataObj = std::make_unique<Soil::soils_na>(); 
    d.MyPenmanMonteith = std::make_unique<PenmanMonteith>(Atmosphere::Cp, Atmosphere::kappa,
            SoilDataObj.air_entry_tension(soil_type),
            SoilDataObj.pore_size_dist(soil_type),SoilDataObj.wilt_point(soil_type),SoilDataObj.porosity(soil_type));


    d.MyPenmanMonteith->Veg_height = cfg.get("Veg_height",1);
    d.MyPenmanMonteith->Veg_height_max = cfg.get("Veg_height_max",1);
    d.MyPenmanMonteith->wind_measurement_height = cfg.get("wind_measurement_height",2); // 2 is default and should be used if U_2m_above_srf is used
    d.MyPenmanMonteith->stomatal_resistance_min = cfg.get("stomatal_resistance_min",1);
    d.MyPenmanMonteith->soil_depth = cfg.get("soil_depth",1);
    d.MyPenmanMonteith->Frac_to_ground = cfg.get("Frac_to_ground",1);
    d.MyPenmanMonteith->LAImin = cfg.get("Leaf_area_index_min",1);
    d.MyPenmanMonteith->LAImax = cfg.get("Lead_area_index_max",1);
    d.MyPenmanMonteith->seasonal_growth = cfg.get("seasonal_growth",1);

}

PM_vars Evapotranspiration_All::set_PenmanMonteith_vars(mesh_elem& face,double& t, double& saturated_vapour_pressure)
{
    PM_vars vars((*face)["U_2m_above_srf"_s],(*face)["iswr"_s],(*face)["netall"_s],t,(*face)["soil_storage"_s],saturated_vapour_pressure,(*face)["P_atm"]);


    return vars;
}

