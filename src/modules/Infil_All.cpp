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
    depends("t")

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
        d.soil_storage = 0.;
            
        // Model Parameters
        infDays = cfg.get("max_inf_days",6);
        min_swe_to_freeze = cfg.get("min_swe_to_freeze",25);
        major = cfg.get("major",5); 
        AllowPriorInf = cfg.get("AllowPriorInf",true);
        ThawType = cfg.get("ThawType",0); // Default is Ayers
        texture = cfg.get("soil_texture",0);
        groundcover = cfg.get("soil_groundcover",0);
        soil_type = cfg.get("soil_type",0); // default is sand
        porosity = cfg.get("soil_porosity",0.5);
        soil_depth = cfg.get("soil_depth",1); // metres, default 1 m
        max_soil_storage = porosity * soil_depth;
        ksaturated = soilproperties[soil_type][KSAT];
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
        d.frozen_phase = 0;

        d.index = 0.;
        d.max_major_per_day = 0.;
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
		        
                Check_for_ice_lens(d,soil_storage_at_freeze,airtemp);
                
                if ((d.major_melt_count == 0 & snowmelt >= major) || swe >= d.init_SWE) {
                    Calc_Index(d,swe);
                    
                    snowinf = Calc_Actual_Inf(d,snowmelt);
                    
                    d.major_melt_count += 1;
                }
                else if (d.major_melt_count > 0 && major_melt_count < infDays) {
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
            double maxinfil = textureproperties[texture][groundcover]; // Currently texture properties is assumed uniform, later make this triangle specific.
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
        Increment_Totals(d,runoff,melt_runoff,inf,snowinf,rain_on_snow);
    }
    else if (ThawType == GREENAMPT) // if not frozen, do GreenAmpt
    {
        d.GA_temp = std::make_unique<tempvars>();

        if(rainfall > 0.0) {
            d.GA_temp->intensity = convert_to_rate_hourly(rainfall);

            if(soil_type == 12){ // handle pavement separately
                runoff = rainfall;
            }
            else if(is_space_in_dry_soil(soil_storage,max_soil_storage,rainfall)){
                inf =  rainfall;
            }
            else {
                
                //double F0 = soil_storage;
                Initialize_GA_Variables(d.GA_temp);

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

                    if (intensity > d.GA_temp->final_rate) { // ponding starts midway through the time step
                        initialize_ponding_vars(d.GA_temp);

                        double ponding_time = global_param->dt() / 3600.0 - d.GA_temp->time_to_ponding;
                        
                        find_final_storage(d.GA_temp,d.GA_temp->storage_at_ponding,ponding_time);
                        
                        d.GA_temp->pond = rainfall - (d.GA_temp->final_storage - d.GA_temp->initial_storage); 
                    }

                }


                inf = d.GA_temp->final_storage - soil_storage; 
                if(d.GA_temp->pond > 0.0){
                    runoff = d.GA_temp->pond; 
                }
            }


            // Increment totals
            Increment_Totals(d,runoff,melt_runoff,inf,snowinf,rain_on_snow);
            
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
void Infil_All::Calc_Index(Infil_All::data &d, const double &swe) {
    d.index = 5 * (1 - theta) * swe ** (0.584);
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

double Infil_All::Calc_Actual_Inf(Infil_All::data &d, const double &melt) {
    double inf = melt * d.index;
    if (inf > d.max_major_per_melt) {
        inf = d.max_major_per_melt;
    }
    return inf;
}


void Infil_All::Check_for_ice_lens(Infil_All::data &d,double &soil_storage_at_freeze, double &t) {
    if (d.major_melt_count > 0 && t < lenstemp) {
        d.major_melt_count = infDays + 1;
    }
}

// Ayers
static double textureproperties[][6] = { // mm/hour
  {7.6, 12.7, 15.2, 17.8, 25.4, 76.2},  // coarse over coarse
  {2.5,  5.1,  7.6, 10.2, 12.7,  15.2}, // medium over medium
  {1.3,  1.8,  2.5,  3.8,  5.1,  6.4},  // medium/fine over fine
  {0.5,  0.5,  0.5,  0.5,  0.5,  0.5}   // soil over shallow bedrock
};


// Green-Ampt Functions

static double soilproperties[][9] = {
  { 0.0,  999.9, 0.000, 0.00, 1.100,  1.000,	0.000,	0.0,  4},  //      0  water
  { 49.5, 117.8, 0.020, 0.10, 0.437,  0.395,	0.121,	4.05, 1},  //      1  sand
  { 61.3,  29.9, 0.036, 0.16, 0.437,  0.41 ,	0.09,	4.38, 4},  //      2  loamsand
  {110.1,  10.9, 0.041, 0.23, 0.453,  0.435,	0.218,	4.9,  2},  //      3  sandloam
  { 88.9,   3.4, 0.029, 0.26, 0.463,  0.451,	0.478,	5.39, 2},  //      4  loam
  {166.8,   6.5, 0.045, 0.38, 0.501,  0.485,	0.786,	5.3,  2},  //      5  siltloam
  {218.5,   1.5, 0.068, 0.38, 0.398,  0.420,	0.299,	7.12, 3},  //      6  saclloam
  {208.8,   1.0, 0.155, 0.39, 0.464,  0.476,	0.63,	8.52, 2},  //      7  clayloam
  {273.3,   1.0, 0.039, 0.40, 0.471,  0.477,	0.356,	7.75, 2},  //      8  siclloam
  {239.0,   0.6, 0.110, 0.41, 0.430,  0.426,	0.153,	10.4, 3},  //      9  sandclay
  {292.2,   0.5, 0.056, 0.43, 0.479,  0.492,	0.49,	10.4, 3},  //      10 siltclay
  {316.3,   0.3, 0.090, 0.46, 0.475,  0.482,	0.405,	11.4, 3},  //      11 clay
  {  0.0,   0.0, 0.000, 0.00, 0.000,  0.000,	0.0,	 0.0, 4}   //      12 pavement. Values not used
};


double Infil_All::convert_to_rate_hourly(double &rainfall) {
    return rainfall / (global_param->dt() / 3600.0);
}

bool Infil_All::is_space_in_dry_soil(double &moist, double &max, double &rainfall) {
    return moist == 0.0 && max >= rainfall;
}

void Infil_All::Initialize_GA_Variables(std::unique_ptr<Infil_All::data::tempvars &GA) {
    GA->soil_storage_deficit = (1.0 - soil_storage/max_soil_storage); 
    GA->initial_rate = calc_GA_infiltration_rate(GA,soil_storage);
    GA->final_storage = GA->initial_storage;
    GA->final_rate = GA->initial_rate;
}

void Infil_All::initialize_ponding_vars(std::unique_ptr<Infil_All::data::tempvars> &GA) {
    GA->storage_at_ponding = ksaturated * GA->capillary_suction / (GA->intensity - ksaturated); 
    GA->time_at_ponding = (GA->storage_at_ponding - GA->initial_storage)/GA->intensity;
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

double Infil_All::calc_GA_infiltration_rate(std::unique_ptr<Infil_All::data::tempvars &GA, double &F){

    return ksaturated*(GA->capillary_suction/F + 1.0);

}

//End Green-Ampt functions
