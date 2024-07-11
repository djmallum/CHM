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
    depends("rainfall_int"); // NEW

    provides("inf");
    provides("total_inf");
    provides("snowinf"); // NEW
    provides("total_snowinf"); // NEW
    provides("total_excess");
    provides("total_meltexcess"); // NEW
    provides("runoff");
    provides("melt_runoff"); // NEW
    provides("total_rain_on_snow") // NEW

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

        d.soil_depth = 400; // TODO uniform soil depth is probably not ideal but ok for now
        d.porosity = .4; // TODO likewise with  uniform porosity. 
                         // Future version: user supplies data to CHM depending on soil type to get porosity
                         // soil depth data could also be supplied similarly
        d.max_storage =  d.soil_depth * d.porosity;
        //        d.storage =  d.max_storage  - (1 - d.max_storage * face->parameter("sm"_s)/100.);
        d.storage =  d.max_storage * face->parameter("sm"_s)/100.; // TODO what is sm? Soil moisture? why is this a parameter? Also _s vs without

        d.total_inf = 0.;
        d.total_snowinf = 0.; // NEW
        d.total_excess = 0.;
        d.total_meltexcess = 0.; // NEW
        d.total_rain_on_snow = 0.; // NEW
        
        d.frozen = false; // NEW, Maybe initial condition, not always necessary because of SWE check to freeze the ground
        d.crackstatus = 0; // NEW, For Gray frozen soil routine, counts number of major melts
        d.frozen_phase = 0; // NEW, Set by the Volumetric module in crhm, once per season. Crhm uses a % system, I am using a bool system (with three options, 0, 1 and 2)
        d.Major = 5; // NEW, default is 5 mm/day in crhm 
        d.Xinfil = new double*[3]; // TODO This has some utility for GA, I might find a better version later. 
                                   // Xinfil[0] is INF/SWE, Xinfil[1] is INF
    


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
    // double C = 2.;
    // double S0 = 1;
    // double SI = face->get_initial_condition("sm")/100.;
    // double TI = 272.;

    // CRHM does total infil for snow and total infil separately, wonder if I should do this
    double runoff = 0.;
    double inf = 0.;
    double snowinf = 0.;
    double melt_runoff = 0.;
        

    double snowmelt = (*face)["snowmelt_int"_s];
    double rainfall = (*face)["rainfall_int"_s]; // NEW
    double swe = (*face)["swe"_s]; 
    double Major = (*face)["Major"_s]; // NEW

    double potential_inf = d.last_ts_potential_inf;
    double avail_storage = (d.max_storage - d.storage);

    if (swe > 25.0 && !d.frozen)
    {
        d.frozen = true; // Initiate frozen soil at 25 mm depth (as in CRHM)
        d.frozen_phase = 0;
    }
    
    if (d.frozen) // Gray's infiltration, 1985
    {
        if (rainfall > 0.0)
        {
            d.rain_on_snow += rainfall;
        }
        
        if (snowmelt > 0.0)
        {
            if (d.frozen_phase == 0) // Unlimited
            {
                inf += snowmelt;
                d.crackstatus = 1; 
            }
            // TODO  Xinfil is not defined
            else if (d.frozen_phase == 1) // Limited
            {
                if (snowmelt >= Major || d.crackstatus >= 1)
                {
                    if (swe > Xinfil[2] && snowmelt >= Major)
                    {
                        // TODO Add Gray equation here
                    }
                    if (snowmelt >= Major)
                    {
                        if (d.crackstatus <= 0)
                        {
                            d.crackstatus = 1;
                        }
                        else
                        {  
                            d.crackstatus += 1;
                        }    
                        
                        snowinf += snowmelt * Xinfil[0];

                        if (snowinf > Xinfil[1])
                        {
                            snowinf = Xinfil[1];
                        }
                    }
                    else
                    {
                        snowinf += snowmelt * Xinfil[0];
                    }

                    if (d.crackstatus > infDays)
                    {
                        snowinf = 0;
                    }
                }
                else
                {
                    if (d.Prior_Inf)
                    {
                        snowinf = snowmelt;
                    }
                }
            }


            else if (d.frozen_phase == 2) // Restricted
            {
                snowinf = 0.;
                d.crackstatus = 1;
            }

            meltrunoff = snowmelt - snowinf;

            if (snowinfil > 0.0)
            {
                snowinf += d.rain_on_snow;
            }
            else
            {
                meltrunoff = d.rain_on_snow;
            }

            // TODO total values incremented here
            // snowinfil, meltrunoff, rainonsnow
            // Zero somethings

        }
    }
    else if (d.ThawType == 0) // if not frozen, do Ayers
    {
        if (rainfall > 0.0)
            double maxinfil = d.texture[texture][groundcover]; // TODO, completey pseduocode, might need units altered
            if (maxinfil > rainfall)
            {
                inf = rainfall;
            }
            else
            {
                inf = maxinfil;
                runoff = rainfall - maxinfil;
            }

            // Increment totals
    }
    else if (d.ThawType == 1) // if not frozen, do GreenAmpt
    {

            inf = 0.0;
            runoff = 0.0;
            snowinf = 0.0;
            meltrunoff = 0.0;
            F1 = 0.0;

            double melt = snowmelt/Global::Freq; // TODO FREQ is from CRHM, need to make this work for CHM. Freq is timesteps/day, so snowmelt/freq converts from snowmelt/day to snowmelt/timestep (this has units of mm, but often written as mm/interval to reduce confusion (which creates some in me...)
            double All = rainfall + melt;

            // TODO remove [hh], define variables like pond, look for lines labeled NEW.

            if(All > 0.0) {

              // TODO Old line, CRHM converted All to mm/h for consistency, not necessarily needed except to ensure units are correct, FIX for CHM
                intensity = All*Global::Freq/24.0; // TODO convert to mm/h, FIX for CHM
                bool is_space_in_dry_soil = soil_moist <= 0.0 && soil_moist_max <= All;
                if(soil_type == 12){ // handle pavement separately
                  runoff = All;
                }
                else if(soil_type[hh] == 0 || is_space_in_dry_soil){ // handle water separately, CRHM version of this simply to check if location is water OR if the soil is dry. I know the soil model handles whether infiltrated actually gets into the soil, but here I changed it to add an upper limit imposed by the space in the soil.
                  infil[hh] =  All;
                }
                else {
                  F1[hh] = soil_moist; // TODO, need to add soil model which tracks this, total moisture in soil, also not declared
                  dthbot    = (1.0 - soil_moist/soil_moist_max); // TODO Not yet declared
                  psidthbot = soilproperties[soil_type[hh]][PSI]*dthbot[hh]; // Not yet declared, also Soil_type not catalogued anywhere
                  f1[hh] = calcf1(F1[hh], psidthbot[hh])*Global::Interval*24.0; // infiltrate first interval rainfall, TODO, not declared

                  infiltrate(); // TODO add function

                  infil = F1 - F0; // TODO F0 not yet declared
        //          cuminfil[hh] += infil[hh];

                  if(pond > 0.0){

                    runoff = pond; //TODO pond not yet declared
                  }
                }
              

              if(melt >= infil[hh]){
                snowinfil[hh] = melt;
                infil[hh] = 0.0;
              }
              else if(melt > 0.0){
                snowinfil[hh] = melt;
                infil[hh] -= snowinfil[hh];
              }
              else
                snowinfil[hh] = 0.0;

              if(melt - snowinfil[hh] >= pond){
                meltrunoff[hh] = melt - snowinfil[hh];
                runoff[hh] = 0.0;
              }
              else if(melt - snowinfil[hh] > 0.0){
                meltrunoff[hh] = melt - snowinfil[hh];
                runoff[hh] = pond - meltrunoff[hh];
              }
              else{
                meltrunoff[hh] = 0.0;
                runoff[hh] = pond;
              }

              cuminfil[hh] += infil[hh];
              cumrunoff[hh] += runoff[hh];

              cumsnowinfil[hh] += snowinfil[hh];
              cummeltrunoff[hh] += meltrunoff[hh];

            } // if(net_rain[hh] + net_snow[hh] > 0.0) greenampt routine
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



