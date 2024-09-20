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
    depends("t");

    provides("inf");
    provides("total_inf");
    provides("snowinf"); // NEW
    provides("total_snowinf"); // NEW
    provides("total_excess");
    provides("total_meltexcess"); // NEW
    provides("runoff");
    provides("melt_runoff"); // NEW
    provides("total_rain_on_snow"); // NEW
    provides("rain_on_snow"); // NEW

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
        d.max_major_per_melt = 0.;
        d.init_SWE = 0.;
        d.soil_storage = 0.;
            
        // Model Parameters
        infDays = cfg.get("max_inf_days",6);
        min_swe_to_freeze = cfg.get("min_swe_to_freeze",25);
        major = cfg.get("major",5); 
        AllowPriorInf = cfg.get("AllowPriorInf",true);
        ThawType = cfg.get("ThawType",0); // Default is Ayers
        texture = cfg.get("soil_texture",0);
        groundcover = cfg.get("soil_groundcover",0);
        lenstemp = cfg.get("temperature_ice_lens",-10.0);
        soil_type = cfg.get("soil_type",0); // default is sand
        porosity = cfg.get("soil_porosity",0.5);
        soil_depth = cfg.get("soil_depth",1); // metres, default 1 m
        max_soil_storage = porosity * soil_depth;
        ksaturated = SoilDataObj.get_soilproperties(soil_type,KSAT);


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
    double melt_runoff = 0.;
    double inf = 0.;
    double snowinf = 0.;
    double rain_on_snow = 0.;

    double snowmelt = (*face)["snowmelt_int"_s];
    double rainfall = (*face)["rainfall_int"_s]; // NEW
    double swe = (*face)["swe"_s]; 
    double soil_storage_at_freeze = (*face)["soil_storage_at_freeze"_s];
    double airtemp = (*face)["t"_s];

    if (swe > min_swe_to_freeze && !d.frozen)
    {
        d.frozen = true; // Initiate frozen soil at 25 mm depth (as in CRHM)

        d.index = 0.;
        d.max_major_per_melt = 0.;
        d.init_SWE = 0.;
    }

    if (d.frozen) // Gray's infiltration, 1985
    {
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
		        
                Check_for_ice_lens(d,airtemp);
                
                if ((d.major_melt_count == 0 & snowmelt >= major) || swe >= d.init_SWE) {
                    Calc_Index(d,swe,soil_storage_at_freeze);
                    
                    snowinf = Calc_Actual_Inf(d,snowmelt);
                    
                    d.major_melt_count += 1;
                }
                else if (d.major_melt_count > 0 && d.major_melt_count < infDays) {
                    snowinf = Calc_Actual_Inf(d,snowmelt);

                    d.major_melt_count += 1;
                }
                else if (d.major_melt_count == 0 and AllowPriorInf) {
                    snowinf = snowmelt;
                }

            }
            else if (soil_storage_at_freeze == 100) // Restricted
            {
                snowinf = 0.;
                d.major_melt_count = 1;
            }

            melt_runoff = snowmelt - snowinf;

            // melt_runoff and snowinf only track melt related quantities
            // total runoff and infiltrated amounts from ANY source are stored in
            // inf and runoff
            
            // This is a weird function, if there is any snowinf, then the rain on the snow also infiltrates
            // this is ported directly from CRHM module crack
            if (snowinf > 0.0)
            {
                inf += rain_on_snow;
            }
            else
            {
                runoff += rain_on_snow;
            }
            runoff += melt_runoff;
            inf += snowinf;

            Increment_Totals(d,runoff,melt_runoff,inf,snowinf,rain_on_snow);
          

        }
    }
    else if (ThawType == AYERS) // if not frozen, do Ayers
    {
        if (rainfall > 0.0)
        {
            double maxinfil = SoilDataObj.get_textureproperties(texture,groundcover); // TODO Currently texture properties is assumed uniform, later make this triangle specific.
            if (maxinfil > rainfall)
            {
                inf = rainfall;
            }
            else
            {
                inf = maxinfil;
                runoff = rainfall - maxinfil;
            }
        }

        // Increment totals
        Increment_Totals(d,runoff,melt_runoff,inf,snowinf,rain_on_snow);
    }
    else if (ThawType == GREENAMPT) // if not frozen, do GreenAmpt
    {
        d.GA_temp = std::make_unique<data::tempvars>();

        if(rainfall > 0.0) {
            d.GA_temp->intensity = convert_to_rate_hourly(rainfall);

            if(soil_type == 12){ // handle pavement separately
                runoff = rainfall;
            }
            else if(is_space_in_dry_soil(d.soil_storage,max_soil_storage,rainfall)){
                inf =  rainfall;
            }
            else {
                
                //double F0 = d.soil_storage;
                Initialize_GA_Variables(d);

                // TODO what about ponding
                if (d.GA_temp->intensity > d.GA_temp->initial_rate) { // ponding is ongoing

                    d.GA_temp->final_storage = d.GA_temp->initial_storage + rainfall;
                    
                    double ponding_time = global_param->dt() / 3600.0;
                    
                    find_final_storage(d.GA_temp,d.GA_temp->initial_storage,ponding_time);    
                    
                    d.GA_temp->pond = rainfall - (d.GA_temp->final_storage - d.GA_temp->initial_storage); 

                    // TODO CRHM version calls find_final_storage again here, but appears unchanged. 
                    // Note that during tests that this could cause a difference.
                }
                else {

                    d.GA_temp->final_storage = d.GA_temp->initial_storage + rainfall;
                    d.GA_temp->final_rate = calc_GA_infiltration_rate(d.GA_temp,d.GA_temp->final_storage); //TODO calcf1 not a function anymore

                    if (d.GA_temp->intensity > d.GA_temp->final_rate) { // ponding starts midway through the time step
                        initialize_ponding_vars(d.GA_temp);

                        double ponding_time = global_param->dt() / 3600.0 - d.GA_temp->time_to_ponding;
                        
                        find_final_storage(d.GA_temp,d.GA_temp->storage_at_ponding,ponding_time);
                        
                        d.GA_temp->pond = rainfall - (d.GA_temp->final_storage - d.GA_temp->initial_storage); 
                    }

                }


                inf = d.GA_temp->final_storage - d.soil_storage; 
                if(d.GA_temp->pond > 0.0){
                    runoff = d.GA_temp->pond; 
                }
            }


            // Increment totals
            Increment_Totals(d,runoff,melt_runoff,inf,snowinf,rain_on_snow);
            d.soil_storage += d.GA_temp->final_storage;  
            d.GA_temp.reset();
        } // if(net_rain[hh] + net_snow[hh] > 0.0) greenampt routine
    }  



    // TODO increment totals, everywhere, maybe do once
    (*face)["total_excess"_s]=d.total_excess;
    (*face)["total_meltexcess"_s]=d.total_meltexcess;
    (*face)["total_inf"_s]=d.total_inf;
    (*face)["total_snownf"_s]=d.total_snowinf;
    (*face)["total_rain_on_snow"_s]=d.total_rain_on_snow;

    (*face)["runoff"_s]=runoff;
    (*face)["inf"_s]=inf;
    (*face)["rain_on_snow"_s]=rain_on_snow;
    (*face)["snowinf"_s]=snowinf;
    (*face)["melt_runoff"_s]=melt_runoff;
}

//General Functions
void Infil_All::Increment_Totals(Infil_All::data &d, double &runoff, double &melt_runoff, double &inf, double &snowinf, double &rain_on_snow) {
    d.total_inf += inf;
    d.total_excess += runoff;
    d.total_snowinf += snowinf;
    d.total_meltexcess += melt_runoff;
    d.total_rain_on_snow += rain_on_snow;
}      


// Crack Functions
void Infil_All::Calc_Index(Infil_All::data &d, double &swe, double &theta) {
    d.index = 5 * (1 - theta/100.0) * std::pow(swe,0.584);
    // d.major_major_per_melt is obtained by dividing d.index by the 
    // total number of time steps to get to d.index
    // This only works if 86400 / dt is a fraction which turns infDays into an integer
    // Example: if dt is 345600 (4 days in seconds) and infDays is 6 days. 
    // the denominator is 1.5, which then requires 2 major melts for it to stop, not 1.5
    // this is actually OK behaviour, but difficult to understand.
     
    d.max_major_per_melt = d.index / (infDays * 86400.0 / global_param->dt() );
    d.index = d.index / swe;
    d.init_SWE = swe;
}

double Infil_All::Calc_Actual_Inf(Infil_All::data &d, double &melt) {
    double inf = melt * d.index;
    if (inf > d.max_major_per_melt) {
        inf = d.max_major_per_melt;
    }
    return inf;
}


void Infil_All::Check_for_ice_lens(Infil_All::data &d, double &t) {
    if (d.major_melt_count > 0 && t < lenstemp) {
        d.major_melt_count = infDays + 1;
    }
}

// Ayers

// Green-Ampt Functions

double Infil_All::convert_to_rate_hourly(double &rainfall) {
    return rainfall / (global_param->dt() / 3600.0);
}

bool Infil_All::is_space_in_dry_soil(double &moist, double &max, double &rainfall) {
    return moist == 0.0 && max >= rainfall;
}

void Infil_All::Initialize_GA_Variables(Infil_All::data &d) {
    // This function requires d.soil_storage so the full object d must be passed.
    // For simple reading, defined GA pointer to be consistent with other functions that use GA 
    // rather than GA_temp
    std::unique_ptr<Infil_All::data::tempvars> &GA = d.GA_temp;
    
    GA->soil_storage_deficit = (1.0 - d.soil_storage/max_soil_storage); 
    GA->initial_rate = calc_GA_infiltration_rate(GA,d.soil_storage);
    GA->initial_storage = d.soil_storage;
    GA->final_storage = GA->initial_storage;
    GA->final_rate = GA->initial_rate;
}

void Infil_All::initialize_ponding_vars(std::unique_ptr<Infil_All::data::tempvars> &GA) {
    GA->storage_at_ponding = ksaturated * GA->capillary_suction / (GA->intensity - ksaturated); 
    GA->time_to_ponding = (GA->storage_at_ponding - GA->initial_storage)/GA->intensity;
}

void Infil_All::find_final_storage(std::unique_ptr<Infil_All::data::tempvars> &GA, \
        double &initial_storage, double &dt) {
    
    double LastF1;

    do {
    
        LastF1 = GA->final_storage;
    
        GA->final_storage = initial_storage + ksaturated*dt + GA->capillary_suction * \
                            log((GA->final_storage + GA->capillary_suction) \
                            / (initial_storage + GA->capillary_suction));
    
    } while(fabs(LastF1 - GA->final_storage) > 0.001);

}

double Infil_All::calc_GA_infiltration_rate(std::unique_ptr<Infil_All::data::tempvars> &GA, double &F){

    return ksaturated*(GA->capillary_suction/F + 1.0);

}

//End Green-Ampt functions
