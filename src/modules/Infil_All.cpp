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
    depends("soil_storage_at_freeze"); // NEW, depends on Volumetric model, equivalent to fallstat in crhm

    provides("inf");
    provides("total_inf");
    provides("snowinf"); // NEW
    provides("total_snowinf"); // NEW
    provides("total_excess");
    provides("total_meltexcess"); // NEW
    provides("runoff");
    provides("melt_runoff"); // NEW
    provides("total_rain_on_snow") // NEW
    provides("rain_on_snow") // NEW

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

        // Triangle Totals
        d.total_inf = 0.;
        d.total_snowinf = 0.; // NEW
        d.total_excess = 0.;
        d.total_meltexcess = 0.; // NEW
        d.total_rain_on_snow = 0.; // NEW

        // Triangle variables
        d.frozen = false; // NEW, Maybe initial condition, not always necessary because of SWE check to freeze the ground
	    d.major_melt_count = 0; // NEW, For Gray frozen soil routine, counts number of major melts
        d.index = 0;
        d.max_major_per_day = 0.;
        d.init_SWE = 0.;
        
        // Model Parameters
        infDays = cfg.get("max_inf_days",6);
        min_swe_to_freeze = cfg.get("min_swe_to_freeze",25);
        major = cfg.get("major",5); 

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

    // CRHM does total infil for snow and total infil separately, wonder if I should do this
    double runoff = 0.;
    double inf = 0.;
    double snowinf = 0.;
    double melt_runoff = 0.;
    double rain_on_snow = 0.;


    double snowmelt = (*face)["snowmelt_int"_s];
    double rainfall = (*face)["rainfall_int"_s]; // NEW
    double swe = (*face)["swe"_s]; 
    double soil_storage_at_freeze = (*face)["soil_storage_at_freeze"_s];
    double avail_storage = (d.max_storage - d.storage);

    if (swe > d.min_swe_to_freeze && !d.frozen)
    {
        d.frozen = true; // Initiate frozen soil at 25 mm depth (as in CRHM)
        d.frozen_phase = 0;

        d.index = 0.;
        d.max_major_per_day = 0.;
        d.init_SWE = 0.;
    }

    if (d.frozen) // Gray's infiltration, 1985
    {
        // TODO this version is not optimal for CHM because it is hourly and this model is daily only.
        // Idea for new version:
        // At the end of every day, check if it is Major and increment the count
        // problem: How to deal with the issue where if the melt is not major, it all infiltrates?
        // Suggestion, every hour we examine the rate and if its on pace for major: it obeys mI/S 
        // and if it is below pace for major, all infiltrates
        // new problem: in the old system
        if (rainfall > 0.0)
        {
            rain_on_snow = rainfall;
        }

        if (snowmelt > 0.0)
        {
            if (soil_storage_at_freeze == 0) // Unlimited
            {
                inf += snowmelt;
                d.major_melt_count = 1; 
            }
            else if (soil_storage_at_freeze > 0 && soil_storage_at_freeze < 100) // Limited
            {
		        if ((d.major_melt_count == 0 & snowmelt >= major) || swe >= d.init_SWE) {
                    // TODO Add Gray eqn call here
                    //
                    // calculate snowinf here, make it a function
                    d.major_melt_count += 1;
                }
                else if (d.major_melt_count > 0 && major_melt_count < infDays) {
                    // TODO calculate snowinf here, make it a function
                    d.major_melt_count += 1;
                }
                else if (d.major_melt_count == 0 and AllowPriorInf) {
                    snowinf = snowmelt;
                }

            }
            // Stopped here aug 14
            else if (soil_storage_at_freeze == 100) // Restricted
            {
                snowinf = 0.;
                d.major_melt_count = 1;
            }

            melt_runoff = snowmelt - snowinf;

            if (snowinf > 0.0)
            {
                snowinf += d.rain_on_snow;
            }
            else
            {
                melt_runoff = d.rain_on_snow;
            }

            // TODO total values incremented here
            // snowinfil, melt_runoff, rainonsnow
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
        melt_runoff = 0.0;
        F1 = 0.0;

        double All = rainfall + snowmelt;

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


            if(snowmelt >= infil[hh]){
                snowinf[hh] = snowmelt;
                infil[hh] = 0.0;
            }
            else if(snowmelt > 0.0){
                snowinf[hh] = snowmelt;
                infil[hh] -= snowinf[hh];
            }
            else
                snowinf[hh] = 0.0;

            if(snowmelt - snowinf[hh] >= pond){
                melt_runoff[hh] = snowmelt - snowinf[hh];
                runoff[hh] = 0.0;
            }
            else if(snowmelt - snowinf[hh] > 0.0){
                melt_runoff[hh] = snowmelt - snowinf[hh];
                runoff[hh] = pond - melt_runoff[hh];
            }
            else{
                melt_runoff[hh] = 0.0;
                runoff[hh] = pond;
            }

            // Increment totals

        } // if(net_rain[hh] + net_snow[hh] > 0.0) greenampt routine
    }  



    // TODO increment totals, everywhere, maybe do once
    (*face)["total_excess"_s]=d.total_excess;
    (*face)["total_meltexcess"_s]=d.total_meltexcess;
    (*face)["total_inf"_s]=d.total_inf;
    (*face)["total_snownf"_s]=d.total_snowinf;

    (*face)["runoff"_s]=runoff;
    (*face)["inf"_s]=inf;
    (*face)["rain_on_snow"_s]=rain_on_snow;
    (*face)["potential_inf"_s]=potential_inf;
    (*face)["soil_storage"_s]= d.storage;
    (*face)["opportunity_time"_s]=d.opportunity_time;
    (*face)["available_storage"_s]=avail_storage;

}

// Crack Functions
bool Infil_All::Do_Initial_INF_Calculation(const double& melt, const double& swe, auto& d) {
    return (melt >= d.major & d.major_melt_count = 0) || swe >= d.init_SWE 
}

double Infil_All::Compute_inf_for_single_melt(const double& melt, const double& swe, auto& d ) {
    double inf = 0;

    if ( Do_Initial_INF_Calculation(melt, swe, d) ) {
        // calc inf and swe

        inf = 

// Green-Ampt Functions
void Infil_All::infiltrate(void){

    F0[hh] = F1[hh];
    f0[hh] = f1[hh];

    if(soil_type[hh] == 0) { // water!
        pond += garain;
        return;
    }
    pond = 0.0;

    f0[hh] = calcf1(F0[hh], psidthbot[hh])*Global::Interval*24.0;

    if(intensity > f0[hh]) {
        ponding();
        return;
    }

    F1[hh] = F0[hh] + garain;

    f1[hh] = calcf1(F1[hh], psidthbot[hh])*Global::Interval*24.0;

    if(intensity > f1[hh])
        startponding();
}

void Infil_All::ponding(void){

    F1[hh] = F0[hh] + garain;

    howmuch(F0[hh], Global::Interval*24.0);

    pond    += garain - (F1[hh] - F0[hh]);

    howmuch(F0[hh], Global::Interval*24.0);

}

void Infil_All::startponding(void){

    double Fp = k[hh]*psidthbot[hh]/(intensity - k[hh]);
    double dt = (Fp - F0[hh])/intensity;

    howmuch(Fp, Global::Interval*24.0 - dt);

    pond    += garain - (F1[hh] - F0[hh]);

}

void Infil_All::howmuch(double F0, double dt) {

    double LastF1;

    do {
        LastF1 = F1[hh];
        // This equation is just equation (6-44) in Dingman. The basic GA solution.
        F1[hh] = F0 + k[hh]*dt + psidthbot[hh]*log((F1[hh] + psidthbot[hh])/(F0 + psidthbot[hh]));
    } while(fabs(LastF1 - F1[hh]) > 0.001);
}

double Infil_All::calcf1(double F, double psidth){

    return k[hh]*(psidth/F + 1.0);

}

//End Green-Ampt functions


// TODO add Gray's 1985 equation here


