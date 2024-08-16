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
        double Vars[4] = {}; // {TOTINF, RATEINF, SUCTION, THETA};

        if(rainfall > 0.0) {
            double intensity = Convert_To_Rate_Hourly(rainfall);

            if(soil_type == 12){ // handle pavement separately
                runoff = rainfall;
            }
            else if(soil_type == 0 || is_space_in_dry_soil(soil_storage,max_soil_storage,rainfall)){ // handle water separately, CRHM version of this simply to check if location is water OR if the soil is dry. I know the soil model handles whether infiltrated actually gets into the soil, but here I changed it to add an upper limit imposed by the space in the soil.
                inf =  rainfall;
            }
            else {
                
                //double F0 = soil_storage;
                double soil_storage_deficit = (1.0 - soil_storage/max_soil_storage); 
                double capillary_suction = soilproperties[soil_type][PSI]*theta; 
                double initial_rate = calc_GA_infiltration_rate(soil_storage, soil_storage_deficit, capillary_suction)*Global::Interval*24.0; // TODO FActor to the right is leftover from CRHM 
                
                if (intensity > initial_rate) {
                    Find_ponding( // Stopped here, August 19, I'm copying the infiltrate() function here. My thinking is that too much is going on in the infiltrate function, which makes it hard to know how things work, so put "everything" here instead and make a few small functions that are more descriptive.


                infiltrate(); // TODO add function

                inf = F1 - soil_storage; 
                soil_storage = F1;
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

double Calc_Actual_Inf(Infil_All::data &d, const double &melt) {
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


double Infil_All::Convert_To_Rate_Hourly(double &rainfall) {
    return rainfall / (global_param->dt() / 3600.0);
}

bool is_space_in_dry_soil(double &moist, double &max, double &rainfall) {
    return moist == 0.0 && max >= rainfall;
}

void Initialize_GA_Variables(double &F1, double &f1, double &theta, double &suction) {
    F1 = soil_storage; 
    theta = (1.0 - soil_storage/max_soil_storage); 
    suction = soilproperties[soil_type][PSI]*Vars[THETA]; 
    f1 = calcf1(F1, soil_storage_deficit, capillary_suction)*Global::Interval*24.0;
}

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


