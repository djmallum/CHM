soil_module::soil_module(config_file cfg)
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
    provides("ground_water_outflow");
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
    for (size_t = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<soil_module::data>(ID);

        d.soil = make_unique<soil_layers>(d);
    
        set_soil_params(d);
        initial_soil_conditions(d);
    }
};

void soil_module::run(mesh_elem& face)
{
    get_soil_inputs(face);

    d.soil.run();

    set_soil_outputs(face);
};

void soil_module::get_soil_inputs(mesh_elem& face)
{
    auto& d->get_module_data<soil_module::data>(ID);

    d.swe = (*face)["swe"_s];
    d.thaw_front_depth = 0.0; //(*face)["thaw_front_depth"_s];
    d.freeze_front_depth = 0.0; //(*face)["freeze_front_depth"_s];
    d.potential_ET = (*face)["ET"_s];
    d.infil = (*face)["inf"_s];
    d.runoff = (*face)["runoff"_s];
    d.routing_residual = 0.0; //(*face)["routine_residual"_s];
};

void soil_module::set_soil_outputs(mesh_elem& face)
{
    (*face)["condensation"_s] = d.condensation;
    (*face)["actual_ET"_s] = d.actual_ET; 
    (*face)["soil_excess_to_runoff"_s] = d.soil_excess_to_runoff; 
    (*face)["soil_excess_to_gw"_s] = d.soil_excess_to_gw; 
    (*face)["ground_water_outflow"_s] = d.ground_water_outflow; 
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
    d.soil_type = 0;

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
