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

#include "Infil_All.hpp"
REGISTER_MODULE_CPP(Infil_All);

Infil_All::Infil_All(config_file cfg)
        : module_base("Infil_All", parallel::data, cfg)
{

    depends("swe");
    depends("snowmelt_int");

    provides("inf");
    provides("total_inf");
    provides("total_excess");
    provides("runoff");
    provides("soil_storage");
    provides("potential_inf");
    provides("opportunity_time");
    provides("available_storage");

}

Infil_All::~Infil_All()
{

}

void Infil_All::init(mesh& domain)
{

    //store all of snobals global variables from this timestep to be used as ICs for the next timestep
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<Infil_All::data>(ID);

        d.soil_depth = 400;
        d.porosity = .4;
        d.max_storage =  d.soil_depth * d.porosity;
        //        d.storage =  d.max_storage  - (1 - d.max_storage * face->parameter("sm"_s)/100.);
        d.storage =  d.max_storage * face->parameter("sm"_s)/100.;

        d.last_ts_potential_inf = 0;
        d.opportunity_time=0.;
        d.total_inf = 0.;
        d.total_excess = 0.;
        d.total_rain_on_snow = 0.;


        d.frozen = false; // TODO not real currently, crhm equivalent is crackon in PrairieInfiltration
        d.frozen_phase = 0; // TODO  not real currently, crhm equivalent is crackstat in PrairieInfiltration

    }
}
void Infil_All::run(mesh_elem &face)
{
    if(is_water(face))
    {
        set_all_nan_on_skip(face);
        return;
    }


    auto& d = face->get_module_data<Infil_All::data>(ID);

    auto id = face->cell_local_id;
    // TODO These are old and for the parametric equation
    double C = 2.;
    double S0 = 1;
    double SI = face->get_initial_condition("sm")/100.;

    double TI = 272.;
    // CRHM does total infil for snow and total infil separately, wonder if I should do this
    double runoff = 0.;
    double inf = 0.;

    double snowmelt = (*face)["snowmelt_int"_s];
    double rainfall = (*face)["snowmelt_int"_s]; // TODO Get correct input, likely from observations
    double swe = (*face)["swe"_s];

    double potential_inf = d.last_ts_potential_inf;
    double avail_storage = (d.max_storage - d.storage);

    if (swe > 25.0 && !d.frozen)
    {
        d.frozen = true; // Initiate frozen soil at 25 mm depth (as in CRHM)
        d.frozen_phase = 0;
    }
    
    if (d.frozen)
    {
        if (rainfall > 0.0)
        {
            d.total_rain_on_snow += rainfall;
        }

        if (snowmelt > 0.0)
        {
            if (d.frozen_phase == 0)
            {
                inf += snowmelt;
                d.crackstatus = 1; 
            }
            else if (d.frozen_phase == 1)
                if (snowmelt > = Major || d.crackstatus >= 1)
                    if (swe > 
        }
    }
    else if (d.ThawType == 0) // Ayers
    {

    }
    else if (d.ThawType == 1) // GreenAmpt
    {

    }
    


    if(snowmelt > 0 || rainfall > 0)
    {
        
        // Note: The potential INF uses the parametric equation which John advises against. So I will need to change this
        // Also doesn't appear to consider RESTRICTED, LIMITED and UNLIMITED specifiers
        // d.opportunity_time += global_param->dt() / 3600.;
        double t0 = d.opportunity_time;
        potential_inf = C * pow(S0,2.92) * pow((1. - SI),1.64) * pow((273.15 - TI) / 273.15, -0.45) * pow(t0,0.44);


        //cap the total infiltration to be no more than our available storage
        if (potential_inf > avail_storage )
        {
            potential_inf = avail_storage;
        }

        d.last_ts_potential_inf = potential_inf;

        if( d.total_inf + snowmelt > potential_inf)
        {
            runoff = (d.total_inf + snowmelt) - potential_inf;
        }


        //infiltrate everything else
        inf    = snowmelt - runoff;


        d.storage += inf;

        if (d.storage > d.max_storage)
        {
            d.storage = d.max_storage;
        }


        d.total_inf += inf;
        d.total_excess += runoff;

    }


    if( !is_nan(face->get_initial_condition("sm")))
    {
        (*face)["total_excess"_s]=d.total_excess;
        (*face)["total_inf"_s]=d.total_inf;

        (*face)["runoff"_s]=runoff;
        (*face)["inf"_s]=inf;
        (*face)["potential_inf"_s]=potential_inf;
        (*face)["soil_storage"_s]= d.storage;
        (*face)["opportunity_time"_s]=d.opportunity_time;
        (*face)["available_storage"_s]=avail_storage;
    }
    else
    {
        set_all_nan_on_skip(face);
    }

}
