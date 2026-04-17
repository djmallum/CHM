#include "katabatic_routing_glacier.hpp"
#include "Atmosphere.h"
#include "Glacier.hpp"
#include "PhysConst.h"
#include "melt_routing_glacier.hpp"
#include <deque>

// data and view constructors
REGISTER_MODULE_CPP(katabatic_routing_glacier);

katabatic_routing_glacier::katabatic_routing_glacier(config_file cfg)
	: module_base("katabatic_routing_glacier", parallel::data, cfg)
{
	depends("Pa");
	depends("t");
	depends("rh");
	depends("t_lapse_rate");	
	depends("snowmelt_int");
    depends("p_subcanopy");
	depends("swe");
    depends("iswr_net");
	depends("iswr_subcanopy");
	depends("ilwr_subcanopy");
    depends("T_rain");

	// All submodule outputs
	provides("glacier_water_equivalent");
	provides("total_depth");
	provides("firnmelt");
	provides("icemelt");
	provides("latent_heat");
	provides("sensible_heat");
	provides("snowmelt_delayed");
	provides("firnmelt_delayed");
	provides("icemelt_delayed");
	provides("total_delayed");
	provides("firn");
	provides("ice");

};

katabatic_routing_glacier::~katabatic_routing_glacier() {};

void katabatic_routing_glacier::init(mesh& domain)
{
    construct_data_objects(domain);
    
    init_glacier();

    init_katabatic();
    
    init_routing(); 
	
};

void katabatic_routing_glacier::run(mesh_elem& face)
{
	
    do_katabatic(face);

	if (is_new_day())
    {
        do_glacier(face);
    }

    do_routing(face);
    
    auto& d = face->get_module_data<data>(ID);
	d.reset_cache();

};

void katabatic_routing_glacier::construct_data_objects(mesh& domain) {
    for (size_t i = 0; i < domain->size_local_faces(); ++i)
    {
        auto face = domain->face(i);
        auto& p_glacier = glacier.get_params();
        auto& p_routing = routing.get_params();
        auto& d = face->make_module_data<data>(ID,face,global_param,&cfg,&p_glacier,&p_routing);
    }
};

void katabatic_routing_glacier::init_glacier()
{
    using namespace Glacier;
    auto& p = glacier.get_params();

    auto densify_version = cfg.get("densify_version","Linear");
    if (densify_version == "Linear")
        p.densify_version = DensifyVersion::Linear;
    else if (densify_version == "HerronLangway")
        p.densify_version = DensifyVersion::HerronLangway;
    else 
    {
        p.densify_version = DensifyVersion::HerronLangway;
        SPDLOG_DEBUG("Unknown densify version specified. Using HerronLangway as the default");
    }

    switch (p.densify_version)
    {
        case DensifyVersion::Linear:
            p.small_increment = Units::DensitySI{cfg.get("small_increment",25.0)};
            p.big_increment = Units::DensitySI{cfg.get("big_increment",50.0)};
            break;
        case DensifyVersion::HerronLangway:
            break;
    };

    p.critical_density = Units::DensitySI{cfg.get("critical_density",550.0)};
    p.firn_to_ice_density = Units::DensitySI{cfg.get("firn_to_ice_density",830.0)};
    
    p.thermal_factor = cfg.get<double>("thermal_factor",0.95);
    p.seconds_per_step = global_param->dt();
};

void katabatic_routing_glacier::init_katabatic() {
    auto& p = katabatic.get_params();
    p.prandtl = cfg.get<double>("Pr",5.0);
    p.k = cfg.get<double>("k",4e-4);
    p.k2 = cfg.get<double>("k2",1.0);
    p.seconds_per_step = global_param->dt();
}

double katabatic_routing_glacier::rain_sun_energy(const mesh_elem& face,data& d)
{
    using namespace PhysConst;
	auto iswr = (*face)["iswr_subcanopy"_s];
	auto ilwr = (*face)["ilwr_subcanopy"_s];
	const auto& firn = d.glacier_state.firn;
	const auto& ice = d.glacier_state.ice;
	auto emissivity = 0.0;
	if (d.swe().value > 0.0)
		return 0.0;
	else if (firn.water_equivalent().value)
		emissivity = d.firn_emissivity;
	else if (ice.water_equivalent().value)
		emissivity = d.ice_emissivity;
	else
		return 0.0;
	
	auto olwr = PhysConst::sbc() * emissivity
		* std::pow(d.air_temperature().value,4.0);
	auto oswr = (*face)["glacier_albedo"_s] * iswr;
	auto Qsun = (ilwr - olwr) + (iswr - oswr);

    constexpr auto M_PER_MM = 1 / 1000.0;
	auto Qrain = Cw() * water_reference_density() * (*face)["p_subcanopy"_s] * M_PER_MM
        * (d.air_temperature().value - d.glacier_temperature().value) / global_param->dt();
	return Qsun + Qrain;
};



void katabatic_routing_glacier::do_katabatic(mesh_elem& face) {
    auto& d = face->get_module_data<data>(ID);
    
    katabatic_view data_k(d,face);
    
    katabatic.execute(data_k);
    
    const auto& cache = d.get_cache();
    
    d.total_energy += rain_sun_energy(face,d) + cache->latent_heat +
        cache->sensible_heat;
};

void katabatic_routing_glacier::do_glacier(mesh_elem& face) {
    auto& d = face->get_module_data<data>(ID);
    
    glacier_view data_g(d,face);
    
    glacier.execute(data_g);

    d.total_energy = 0.0;
};

void katabatic_routing_glacier::do_routing(mesh_elem& face) {

    auto& d = face->get_module_data<data>(ID);
    
    routing_view data_r(d,face);
    
    routing.execute(data_r);
};

bool katabatic_routing_glacier::is_new_day()
{
	// TODO This has hard coded elements, Chris suggested something different here: https://godbolt.org/z/3c51T1avT
    int td = global_param->posix_time().time_of_day().total_seconds();
    int time_to_midnight = 86400 - td;
    if (td >= 0 && td < global_param->dt()) //(time_to_midnight >= global_param->dt())
    {
        return true;
    }
    else
        return false;

};

katabatic_routing_glacier::katabatic_view::~katabatic_view()
{	
	auto& cache = d.get_cache();

	(*face)["latent_heat"_s] = cache->latent_heat;	
	(*face)["sensible_heat"_s] = cache->sensible_heat;
};

katabatic_routing_glacier::routing_view::~routing_view()
{
	auto& cache = d.get_cache();

	(*face)["snowmelt_delayed"_s] = cache->snowmelt_delayed;
	(*face)["firnmelt_delayed"_s] = cache->firnmelt_delayed;
	(*face)["icemelt_delayed"_s] = cache->icemelt_delayed;
	(*face)["total_delayed"_s] = cache->snowmelt_delayed + cache->firnmelt_delayed + cache->icemelt_delayed;
};

katabatic_routing_glacier::glacier_view::~glacier_view()
{
	auto& cache = d.get_cache();
	(*face)["glacier_water_equivalent"] = d.glacier_state.total_water_equiv().value;
	(*face)["total_depth"] = d.glacier_state.total_depth().value;
	(*face)["firnmelt"] = cache->firnmelt;
	(*face)["icemelt"] = cache->icemelt;

	(*face)["firn"] = d.glacier_state.firn.water_equivalent().value;
	(*face)["ice"] = d.glacier_state.ice.water_equivalent().value;
};

static Glacier::Ice construct_ice(const Glacier::Params* p,const mesh_elem& face)
{

	using namespace Glacier;
	auto init = Units::Milimetres{face->land_attribute("InitialGlacierIce")};
	Ice ice(p,init);

	return ice;
};

static Glacier::LayeredFirn construct_firn(const Glacier::Params* p,const mesh_elem& face)
{
    // Only supporting single initial firn layer start up
    
    //  TODO Expand this for multiple layers either via linear assumption
    // for initial layering, or using a "wind-up" approach where intiial firn 
    // is divided up and densified.
    // Possible source of error: if density is too large after wind-up, then firn converts to
    // ice by accident.
	using namespace Glacier;
    auto init = Units::Milimetres{face->land_attribute("InitialGlacierFirn")};
    std::deque<Layer> layers;
    layers.emplace_back(init);

    LayeredFirn firn(p,layers);

	return firn;
};
katabatic_routing_glacier::data::data(const mesh_elem& face_in, std::shared_ptr<global> global_param, const config_file* cfg,
		const Glacier::Params* p_g, const GlacierRouting::Params* p_r) 
	: data_base(face_in,global_param,cfg), routing_state(p_r), glacier_state(construct_firn(p_g,face),construct_ice(p_g,face)) {};

Units::Kelvin katabatic_routing_glacier::data::glacier_temperature()
{
    static const auto T = get_domain_param("glacier_temperature",273.15);

	return Units::Kelvin{T};
};
Units::Kelvin katabatic_routing_glacier::data::air_temperature()
{
	update_value( [this]() -> auto& { return cache_->air_temperature; },
			[this]() { return (*face)["t"_s]; } );

	return Units::Kelvin{cache_->air_temperature + 273.15};
};
Units::Pa katabatic_routing_glacier::data::air_pressure()
{
	update_value( [this]() -> auto& { return cache_->air_pressure; },
			[this]() { return (*face)["Pa"_s]; } );

	return Units::Pa{cache_->air_pressure};
};
Units::Pa katabatic_routing_glacier::data::vapour_pressure()
{
	update_value( [this]() -> auto& { return cache_->vapour_pressure; },
			[this]() { 
                auto rh = (*face)["rh"_s];
                auto output = rh * Atmosphere::saturatedVapourPressure(air_temperature().value);
                return output; }
                );

	return Units::Pa{cache_->vapour_pressure};
};
Units::Pa katabatic_routing_glacier::data::vapour_pressure_surface()
{
    static const auto value = get_domain_param("vapour_pressure_surface",611.3);

	return Units::Pa{cache_->vapour_pressure_surface};
};
Units::LapseRateSI katabatic_routing_glacier::data::lapse_rate()
{
	update_value( [this]() -> auto& { return cache_->lapse_rate; },
			[this]() { return (*face)["t_lapse_rate"_s]; } );

	return Units::LapseRateSI{cache_->lapse_rate};
};
void katabatic_routing_glacier::data::latent_heat(const double v)
{
	set_output( [this]() -> auto& { return cache_->latent_heat; },
			v);
};
void katabatic_routing_glacier::data::sensible_heat(const double v)
{
	set_output( [this]() -> auto& { return cache_->sensible_heat; },
			v);
};

const Units::Milimetres katabatic_routing_glacier::data::snowmelt()
{
	update_value( [this]() -> auto& { return cache_->snowmelt; },
			[this]() { return (*face)["snowmelt_int"_s]; } );

	return Units::Milimetres{cache_->snowmelt};
};
const Units::Milimetres katabatic_routing_glacier::data::firnmelt()
{
	if (!cache_ || std::isnan(cache_->firnmelt))
		return Units::Milimetres{0.0};
	
	return Units::Milimetres{cache_->firnmelt};
};
const Units::Milimetres katabatic_routing_glacier::data::icemelt()
{
	if (!cache_ || std::isnan(cache_->icemelt))
		return Units::Milimetres{0.0};
	
	return Units::Milimetres{cache_->icemelt};
};
void katabatic_routing_glacier::data::snowmelt_delayed(double v)
{
	set_output( [this]() -> auto& { return cache_->snowmelt_delayed; },
			v);
};
void katabatic_routing_glacier::data::firnmelt_delayed(double v)
{
	set_output( [this]() -> auto& { return cache_->firnmelt_delayed; },
			v);

};
void katabatic_routing_glacier::data::icemelt_delayed(double v)
{
	set_output( [this]() -> auto& { return cache_->icemelt_delayed; },
			v);
};
void katabatic_routing_glacier::data::total_delayed(double v)
{
	set_output( [this]() -> auto& { return cache_->total_delayed; },
			v);
};

const Units::Milimetres katabatic_routing_glacier::data::swe()
{
	update_value( [this]() -> auto& { return cache_->swe; },
			[this]() { return (*face)["swe"_s]; } );

	return Units::Milimetres{cache_->swe};
};
const Units::Watts_per_m2 katabatic_routing_glacier::data::melt_energy()
{
    static const auto firn_emissivity = cfg->get<double>("firn_emissivity");
    static const auto ice_emissivity = cfg->get<double>("ice_emissivity");
	auto iswr = (*face)["iswr_subcanopy"_s];
	auto ilwr = (*face)["ilwr_subcanopy"_s];
	const auto& firn = glacier_state.firn;
	const auto& ice = glacier_state.ice;
	auto emissivity = 0.0;
	if (swe().value > 0.0)
		return Units::Watts_per_m2{0.0};
	else if (firn.water_equivalent().value)
		emissivity = firn_emissivity;
	else if (ice.water_equivalent().value)
		emissivity = ice_emissivity;
	else
		return Units::Watts_per_m2{0.0};
	

	auto olwr = PhysConst::sbc() * emissivity
		* std::pow(air_temperature().value,4.0);
	auto oswr = (*face)["glacier_albedo"_s] * iswr;
	auto Qsun = (ilwr - olwr) + (iswr - oswr);
	auto Qrain = (*face)["Qrain"_s];
	// Means that katabatic_melt_energy didn't run
	if (!cache_)
		return Units::Watts_per_m2{Qsun + Qrain};

	auto Q_sensible = 0.0;	
	auto Q_latent = 0.0;

	if ( !std::isnan(cache_->latent_heat) )
		Q_latent = cache_->latent_heat;

	if ( !std::isnan(cache_->sensible_heat) )
		Q_sensible = cache_->sensible_heat;

	return Units::Watts_per_m2{Qsun + Qrain + Q_sensible + Q_latent};
};
bool katabatic_routing_glacier::data::update_now()
{
	static const auto day = get_domain_param("change_day"_s,300u);
	
	auto p = global_param->posix_time();
	auto d = p.date();

	auto jd = d.day_of_year();
	std::cout << jd << std::endl;
	
	return day == jd; 
};

void katabatic_routing_glacier::data::glacier_water_equivalent(const double out)
{
	set_output( [this]() -> auto& { return cache_->glacier_water_equivalent; },
			out);
};
void katabatic_routing_glacier::data::total_depth(const double out)
{
	set_output( [this]() -> auto& { return cache_->total_depth; },
			out);
};
void katabatic_routing_glacier::data::firn_melt(const double out)
{
	set_output( [this]() -> auto& { return cache_->firnmelt; },
			out);
};
void katabatic_routing_glacier::data::ice_melt(const double out)
{
	set_output( [this]() -> auto& { return cache_->icemelt; },
			out);
};


