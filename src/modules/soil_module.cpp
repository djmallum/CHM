#include "soil_module.hpp"

soil_module::soil_module(config_file cfg) : module_base("soil_module", parallel::data, cfg)
{
    depends("swe");
    depends("thaw_front_depth"); 
    depends("freeze_front_depth");
    depends("potential_ET");
    depends("inf");
    depends("runoff");
    depends("routing_residual");

    provides("condensation");
    provides("actual_ET");
    provides("soil_excess_to_runoff");
    provides("soil_excess_to_gw");
    provides("ground_water_out");
    provides("soil_to_ssr");
    provides("soil_storage");
    provides("soil_rechr_storage");
    provides("depression_storage");
    provides("ground_water_storage");
};

soil_module::~soil_module()
{

};

void soil_module::init(mesh& domain)
{
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<soil_module::data>(ID);
        d.my_face = std::make_shared<mesh_elem>(face);

        // it might be wise to actually have soil have an instance of ET inside of it, rather than separate here.
        //
        // while the processes are not tightly coupled, they operate on the same construct 
        d.soil_layers = std::make_unique<soil_two_layer>(d);
        d.ET = std::make_unique<soil_ET>(d);  
        if (d.soil_layers)
            set_soil_params(d);
        if (d.ET)
            set_ET_params(d);
    
        initial_soil_conditions(d);

        set_local_module(d);

    }
};

void soil_module::run(mesh_elem& face)
{

    auto& d = face->make_module_data<soil_module::data>(ID);
    
    get_soil_inputs(face,d);


    // order of operations here is hard coded, but it wouldn't be physically wrong to impose ET before soil
    // This is fine because this soil module is run in this order.

    d.soil_layers->run();
    
    d.ET->run();

    set_soil_outputs(face,d);
};

void soil_module::get_soil_inputs(mesh_elem& face,soil_module::data& d)
{
    d.swe = (*face)["swe"_s];
    d.thaw_front_depth = 0.0; //(*face)["thaw_front_depth"_s];
    d.freeze_front_depth = 0.0; //(*face)["freeze_front_depth"_s];
    d.potential_ET = (*face)["ET"_s];
    d.infil = (*face)["inf"_s];
    d.runoff = (*face)["runoff"_s];
    d.routing_residual = 0.0; //(*face)["routine_residual"_s];
};

void soil_module::set_soil_outputs(mesh_elem& face,soil_module::data& d)
{
    (*face)["condensation"_s] = d.condensation;
    (*face)["actual_ET"_s] = d.actual_ET; 
    (*face)["soil_excess_to_runoff"_s] = d.soil_excess_to_runoff; 
    (*face)["soil_excess_to_gw"_s] = d.soil_excess_to_gw; 
    (*face)["ground_water_out"_s] = d.ground_water_out; 
    (*face)["soil_to_ssr"_s] = d.soil_to_ssr;
    (*face)["soil_storage"_s] = d.soil_storage;
    (*face)["soil_rechr_storage"_s] = d.soil_rechr_storage;
    (*face)["depression_storage"_s] = d.depression_storage;
    (*face)["ground_water_storage"_s] = d.ground_water_storage;
};

void soil_module::set_soil_params(soil_module::data& d)
{
    // TODO actually connect to stuff
    d.soil_storage_max = 0.0;
    d.soil_rechr_max = 0.0;
    d.excess_to_ssr = true; 
    d.detention_max = 0.0;
    d.detention_snow_max = 0.0;
    d.detention_organic_max = 0.0;
    d.depression_max = 0.0;
    d.ground_water_max = 0.0;
    d.ground_cover_type = 0.0;
    d.soil_type_rechr = 0;
    d.soil_type_lower = 0;
};


void soil_module::set_ET_params(soil_module::data& d)
{
    d.ground_cover_type = 1; //TODO access ground cover type from face and then figure out how to get it right
    d.soil_type_rechr = 1; //TODO same as above. 
    d.soil_type_lower = 1;
};


void soil_module::initial_soil_conditions(soil_module::data& d)
{
    // TODO actually connect to stuff
    // requires MESHER or at least data for one station
    d.soil_storage = 0.0;
    d.soil_rechr_storage = 0.0;
    d.thaw_fraction_rechr = 0.0;
    d.thaw_fraction_lower = 0.0;
    d.soil_rechr_storage = 0.0;
    d.detention_snow_init = 0.0;
    d.detention_organic_init = 0.0;
    d.depression_storage = 0.0;
    d.ground_water_storage = 0.0;
};

bool soil_module::data::is_lake(soil_ET_DTO& DTO)
{
    soil_module::data& d = static_cast<soil_module::data&>(DTO);

    return d.local_module->is_water(*d.my_face);
};

void soil_module::set_local_module(soil_module::data& d)
{
    d.local_module = std::make_shared<soil_module>(*this);
    //d.local_module = new soil_module::data::my_module(*this);
};
    // TODO this is just a copy of is_water and this is a bad practice but currently the is_water function is not accessible by the data class. Fix: create a separate object taht module_base inherits that contains these functions. face_info will also inherit these functions.
 //   soil_module::data& d = static_cast<soil_module::data&>(DTO);
 //   bool is = false;

 //   if(d.my_face->has_parameter("landcover"_s))
 //   {
 //       int LC = face->parameter("landcover"_s);
 //       is = global_param->parameters.get<bool>("landcover." + std::to_string(LC) + ".is_water",false);
 //   }
 //   return is
      
//}; 
