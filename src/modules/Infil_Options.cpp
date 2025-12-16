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

#include "Infil_Options.hpp"
#include <string_view>
#include <cassert>

REGISTER_MODULE_CPP(Infil_Options);


Infil_Options::Infil_Options(config_file cfg) : module_base("Infil_Options", parallel::data, cfg)
{

    depends("swe");
    depends("snowmelt_int");
    depends("rainfall_int");
    depends("t");

    provides("infiltrated");
    provides("runoff");
    provides("snow_infiltrated");
    provides("melt_runoff");

};

Infil_Options::~Infil_Options() {};

void Infil_Options::init(mesh& domain)
{
    auto major_melt_threshold = cfg.get("major",5);
    auto infDays = cfg.get("max_inf_days",6);
    auto lenstemp = cfg.get("temperature_ice_lens",-10.0);
    auto allow_early_inf = cfg.get("AllowPriorInf",true);

    Crack::Params::set(major_melt_threshold,infDays,lenstemp,allow_early_inf);

    auto min_swe_to_freeze = cfg.get("min_swe_to_freeze",25.0);
    auto day_of_year_to_freeze = static_cast<size_t>(cfg.get("day_of_year_to_freeze",300));
    constants_crack = std::make_unique<DomainConstantsCrack>(min_swe_to_freeze,day_of_year_to_freeze);

    for (size_t i = 0; i < domain->size_local_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<Infil_Options::data>(ID,face,global_param,cfg);

        d._texture = face->soil_attribute<std::string>("soil_texture","soils");
        d._ground_cover = face->soil_attribute<std::string>("soil_ground_cover","soils");
        const auto& melt_per_step = (*face)["snowmelt_int"_s];
        const auto& rain_per_step = (*face)["rainfall_int"_s];
        d.state.daily_melt_total.bind_target(melt_per_step);
        d.state.daily_rain_total.bind_target(rain_per_step);
    };
};

void Infil_Options::run(mesh_elem& face)
{

    auto& d = face->get_module_data<Infil_Options::data>(ID);
    auto& s = d.get_state();
    constexpr double DECIMAL_TO_PERCENT = 100.0;

    if (d.is_newday()) [[unlikely]] 
    {
        d.crack_details.status = get_status(d,d.crack_details.status);
        if (global_param->day() == constants_crack->day_of_year_to_freeze) [[unlikely]]
        {
            s.soil_saturation_at_freeze = 
                (*face)["soil_saturation"_s] * DECIMAL_TO_PERCENT;
            d.saturation_set = true;
        }    
    }

    switch(d.crack_details.status) 
    {
        case Status::THAWED:
            ayers.execute(d);
            break;
        case Status::FROZEN:
            crack.execute(d);
            break;
        case Status::TO_THAWED:
            if (d.is_newday())
            {
                s.soil_saturation_at_freeze = 0.0; 
                d.saturation_set = false;
            }
            crack.execute(d);
            break;
        case Status::TO_FROZEN:
            if (!d.saturation_set)
            {
                s.soil_saturation_at_freeze = 
                    (*face)["soil_saturation"_s] * DECIMAL_TO_PERCENT;
                d.saturation_set = true;
            }
            crack.execute(d);
            break;
    }

    d.set_outputs();

};

void Infil_Options::data::set_outputs() const
{
    auto& c = get_cache();
    if (!c)
    {
        (*face)["infiltrated"_s] = 0.0;
        (*face)["snow_infiltrated"_s] = 0.0;
        (*face)["runoff"_s] = 0.0;
        (*face)["melt_runoff"_s] = 0.0;
        (*face)["rain_on_snow"_s] = 0.0;
    }
    else
    {
        (*face)["infiltrated"_s] = c->infiltrated;
        (*face)["snow_infiltrated"_s] = c->snow_infiltrated;
        (*face)["runoff"_s] = c->runoff;
        (*face)["melt_runoff"_s] = c->melt_runoff;
        (*face)["rain_on_snow"_s] = c->rain_on_snow;
    }
};

bool Infil_Options::data::is_newday() const
{
    // TODO This has hard coded elements, Chris suggested something different here: https://godbolt.org/z/3c51T1avT
	auto td = global_param->posix_time().time_of_day().total_seconds();
    auto time_to_midnight = 86400 - td;
    if (td >= 0 && td < global_param->dt()) //(time_to_midnight >= global_param->dt())
    {
        return true;
    }
    else
        return false;
};

double Infil_Options::data::snowmelt()
{
    update_value(
            [this]() -> auto& { return cache_->snowmelt;},
            [this]() { return (*face)["snowmelt_int"_s];}
            );

    return cache_->snowmelt;
};

double Infil_Options::data::rainfall()
{
    update_value(
            [this]() -> auto& { return cache_->rainfall;},
            [this]() { return (*face)["rainfall_int"_s];}
            );

    return cache_->rainfall;
};

double Infil_Options::data::swe()
{
    update_value(
            [this]() -> auto& { return cache_->swe;},
            [this]() { return (*face)["swe"_s];}
            );

    return cache_->swe;
};

double Infil_Options::data::air_temperature()
{
    update_value(
            [this]() -> auto& { return cache_->air_temperature;},
            [this]() { return (*face)["t"_s];}
            );

    return cache_->air_temperature;
};

std::string_view Infil_Options::data::texture()
{
    assert(!_texture.empty() && "Infil_Options: texture much be allocated");

    return _texture;
};

std::string_view Infil_Options::data::ground_cover()
{
    assert(!_ground_cover.empty() && "Infil_Options: ground_cover much be allocated");

    return _ground_cover;
};

void Infil_Options::data::infiltrated(const double out)
{
    set_output(
            [this]() -> auto& { return cache_->infiltrated;},
            out
            );
};

void Infil_Options::data::runoff(const double out)
{
    set_output(
            [this]() -> auto& { return cache_->runoff;},
            out
            );
};

void Infil_Options::data::snow_infiltrated(const double out)
{
    set_output(
            [this]() -> auto& { return cache_->snow_infiltrated;},
            out
            );
};

void Infil_Options::data::melt_runoff(const double out)
{
    set_output(
            [this]() -> auto& { return cache_->melt_runoff;},
            out
            );
};

void Infil_Options::data::rain_on_snow(const double out)
{
    set_output(
            [this]() -> auto& { return cache_->rain_on_snow;},
            out
            );
};

Crack::State& Infil_Options::data::get_state()
{
    return state;
};

Infil_Options::Status Infil_Options::get_status(data& d,const Status Old)
{
    switch(Old)
    {
        case Status::TO_FROZEN:
            return Status::FROZEN;

        case Status::TO_THAWED:
            return Status::THAWED;

        case Status::FROZEN:
            if (d.swe() > 0.0)
                return Status::FROZEN;
            else
                return Status::TO_THAWED;

        case Status::THAWED:
            if (d.swe() <= constants_crack->min_swe_to_freeze)
                return Status::THAWED;
            else
                return Status::TO_FROZEN;
    }

};
