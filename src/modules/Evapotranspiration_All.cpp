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
    alpha = cfg.get("alpha_PriestelyTaylor",1.26);
    wind_height = cfg.get("wind_measurement_height,2);
    stomatal_resistance_min = cfg.get("stomatal_resistance_min",62);
    soil_depth = cfg.get("soil_depth",1);
    Frac_to_ground = cfg.get("Frac_to_ground",1); // iswr_subcanopy exists.
                                                  
    SoilDataObj = std::make_unique<Soil::soils_na>();

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<Infil_All::data>(ID);
        
        // Consider if an if statement is necessary.


        init_PenmanMonteith(d, face, wind_height, stomatal_resistance_min, soil_depth, Frac_to_ground);
        
        init_PriestelyTaylor(d,alpha);

        // put PristelyTaylor parameters here.
    }

}

void Evapotranspiration_All::run(mesh_elem& face)
{

    auto& d = face->get_module_data<Evapotranspiration_All::data>(ID);

    if (is_water(face))
    {
        PT_vars my_PT_vars = set_PriestelyTaylor(face);

        model_output output;

        d.MyPriestelyTaylor->CalcEvapT(my_PT_vars,output);
        // Do PriestlyTaylor
    }
    // TODO wetlands or saturated soils, need input from soils module.
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
        double VP = (*face)["rh"_s] * SVP;
        PM_vars my_PM_vars = set_PenmanMoneith_vars(face,t,SVP,VP);
    
        model_output output;
    
        d.MyPenmanMonteith->CalcEvapT(my_PM_vars,output);
    }
    
    // TODO total_ET, as well as PT ET and PM ET as separate. 
    (*face)["ET"_s] = output.ET;
    
}

Evapotranspiration_All::~Evapotranspiration_All()
{

}

void Evapotranspiration_All::init_PriestelyTaylor(Evapotranspiration_All::data& d,double& alpha)
{
    d.MyPriestelyTaylor = std::make_unique<PriestelyTaylor>(alpha);
}

void Evapotranspiration_All::init_PenmanMonteith(Evapotransporation_All::data& d,mesh_elem& face, double& wind_height, double& stomatal_resistance_min, double& soil_depth, double& Frac_to_ground)
{
    
    const std::string soil_type = cfg.get("soil_type","sand");

   
    
    const double& Cp = Atmoshpere::Cp;
    const double& kappa = Atmosphere::kappa;
    const double& air_entry_tension = SoilDataObj.air_entry_tension(soil_type);
    const double& pore_size_dist = SoilDataObj.pore_size_dist(soil_type);
    const double& wilt_point = SoilDataObj.wilt_point(soil_type);
    const double& porosity = SoilDataObj.porosity(soil_type);

    // Leaf area index is not used if no vegetation, but LAI and LAImax are references in the PenmanMonteith model, therefore values are needed for initialization. It is ok if these values go out of scope as long as there is no vegetation. 
    if (face->has_vegetation())
    {
        double& LAI = face->veg_attribute("LAI");
        double& LAImax = face->veg_attribute("LAImax");

        d.MyPenmanMonteith = std::make_unique<PenmanMonteith>(LAI, LAImax, Cp, kappa, 
                air_entry_tension, pore_size_dist, wilt_point, porosity);

    }
    else
    {
        d.MyPenmanMonteith = std::make_unique<PenmanMonteith>(dummyvar, dummyvar, Cp, kappa, 
                air_entry_tension, pore_size_dist, wilt_point, porosity);
    }
    
    // Currently, one copy per triangle because that may be true later.
    // For now it is uniform, so set a single copy per module instance (one per mpi process). 
    d.MyPenmanMonteith->wind_measurement_height = wind_height // 2 is default and should be used if U_2m_above_srf is used
    d.MyPenmanMonteith->stomatal_resistance_min = stomatal_resistance_min; // default based on Armstrong (2011) value, page 56.
    d.MyPenmanMonteith->soil_depth = soil_depth;
    d.MyPenmanMonteith->Frac_to_ground = Frac_to_ground; // iswr_subcanopy exists.

}

PM_vars Evapotranspiration_All::set_PenmanMonteith_vars(mesh_elem& face,double& t, double& saturated_vapour_pressure,double& vapour_pressure)
{
    PM_vars vars((*face)["U_2m_above_srf"_s],(*face)["iswr"_s],(*face)["netall"_s],t,(*face)["soil_storage"_s],vapour_pressure,saturated_vapour_pressure,(*face)["P_atm"]); 


    return vars;
}

PT_vars Evapotransporaqtion_All:set_PriestelyTaylor_vars(mesh_elem& face)
{
    // TODO P_atm, is a state variable and should have a copy local to the run function 
    // because it is calculated from the Atmosphere namespace, not done yet
    PT_vars vars((*face)["netall"_s],(*face)["P_atm"]);

    return vars;
}
