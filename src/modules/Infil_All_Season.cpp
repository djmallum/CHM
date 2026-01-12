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

#include "Infil_All_Season.hpp"
#include <cassert>

REGISTER_MODULE_CPP(Infil_All_Season);


Infil_All_Season::Infil_All_Season(config_file cfg) : module_base("Infil_All_Season", parallel::data, cfg)
{

    depends("swe");
    depends("snowmelt_int");
    depends("rainfall_int");
    depends("t");

    provides("rain_on_snow");
    provides("infiltrated");
    provides("runoff");
    provides("snow_infiltrated");
    provides("melt_runoff");

};

Infil_All_Season::~Infil_All_Season() {};

void Infil_All_Season::init(mesh& domain)
{
    Crack::State::major_melt_threshold = cfg.get("major",5);
    Crack::State::infDays = cfg.get("max_inf_days",6);
    Crack::State::lenstemp = cfg.get("temperature_ice_lens",-10.0);
    Crack::State::allow_early_inf = cfg.get("AllowPriorInf",true);

    new_day.set_global(global_param);
    
    auto min_swe_to_freeze = cfg.get("min_swe_to_freeze",25.0);
    
    // TODO this is now broken I think, because algorithm_selector is nowa private member of data
    // Think about this.
    algorithm_selector.set_min_swe(min_swe_to_freeze);
    day_of_year_to_freeze = static_cast<size_t>(cfg.get("day_of_year_to_freeze",300));

    for (size_t i = 0; i < domain->size_local_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<Infil_All_Season::data>(ID,face,global_param,cfg);

        d.ayers_params.texture = face->soil_attribute<std::string>("soil_texture","soils");
        d.ayers_params.ground_cover = face->soil_attribute<std::string>("soil_ground_cover","soils");
        const auto& melt_per_step = (*face)["snowmelt_int"_s];
        const auto& rain_per_step = (*face)["rainfall_int"_s];
        auto& s = d.get_state();
        s.daily_melt_total.bind_target(melt_per_step);
        s.daily_rain_total.bind_target(rain_per_step);
    };

};

void Infil_All_Season::run(mesh_elem& face)
{

    auto& d = face->get_module_data<Infil_All_Season::data>(ID);

    if (new_day.check()) [[unlikely]] 
    {
        d.crack_details.status = algorithm_selector.get(
                [&d]() {return d.swe(); });

        auto& s = d.get_state();
        if (global_param->day() == day_of_year_to_freeze) [[unlikely]]
        {
            s.soil_saturation_at_freeze = 
                (*face)["soil_saturation"_s] * Infil_All_Season::DECIMAL_TO_PERCENT;
            d.crack_details.saturation_set = true;
        }    
    }

    algorithm_runner.run(new_day,d,face);

    set_outputs(face,d);

};

template<typename T>
void runner<T>::run(const new_day_checker& new_day,T& d,mesh_elem& face)
{
    auto& s = d.get_state();
    auto& cr = d.crack_details;
    switch(d.crack_details.status) 
    {
        case Status::TO_THAWED:
            if (d.is_newday())
            {
                s.soil_saturation_at_freeze = 0.0; 
                cr.saturation_set = false;
            }
        case Status::THAWED:
            ayers.execute(d);
            break;
        case Status::TO_FROZEN:
            if (!cr.saturation_set)
            {
                s.soil_saturation_at_freeze = 
                    (*face)["soil_saturation"_s] * Infil_All_Season::DECIMAL_TO_PERCENT;
                cr.saturation_set = true;   
            }
        case Status::FROZEN:
            crack.execute(d);
            break;
    }
};

void Infil_All_Season::set_outputs(mesh_elem& face,const data& d) const
{
    auto& c = d.get_cache();
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

Infil_All_Season::data::data(mesh_elem& face_in,std::shared_ptr<global> param, config_file cfg)
    : data_base<Cache>(face_in, param, cfg) {};

bool Infil_All_Season::data::is_newday() const
{
    return new_day.check();
};

double Infil_All_Season::data::snowmelt()
{
    update_value(
            [this]() -> auto& { return cache_->snowmelt;},
            [this]() { return (*face)["snowmelt_int"_s];}
            );

    return cache_->snowmelt;
};

double Infil_All_Season::data::rainfall()
{
    update_value(
            [this]() -> auto& { return cache_->rainfall;},
            [this]() { return (*face)["rainfall_int"_s];}
            );

    return cache_->rainfall;
};

double Infil_All_Season::data::swe()
{
    update_value(
            [this]() -> auto& { return cache_->swe;},
            [this]() { return (*face)["swe"_s];}
            );

    return cache_->swe;
};

double Infil_All_Season::data::air_temperature()
{
    update_value(
            [this]() -> auto& { return cache_->air_temperature;},
            [this]() { return (*face)["t"_s];}
            );

    return cache_->air_temperature;
};

const std::string& Infil_All_Season::data::texture()
{
    return ayers_params.texture;
};

const std::string& Infil_All_Season::data::ground_cover()
{
    return ayers_params.ground_cover;
};

void Infil_All_Season::data::infiltrated(const double out)
{
    set_output(
            [this]() -> auto& { return cache_->infiltrated;},
            out
            );
};

void Infil_All_Season::data::runoff(const double out)
{
    set_output(
            [this]() -> auto& { return cache_->runoff;},
            out
            );
};

void Infil_All_Season::data::snow_infiltrated(const double out)
{
    set_output(
            [this]() -> auto& { return cache_->snow_infiltrated;},
            out
            );
};

void Infil_All_Season::data::melt_runoff(const double out)
{
    set_output(
            [this]() -> auto& { return cache_->melt_runoff;},
            out
            );
};

void Infil_All_Season::data::rain_on_snow(const double out)
{
    set_output(
            [this]() -> auto& { return cache_->rain_on_snow;},
            out
            );
};

Crack::State& Infil_All_Season::data::get_state()
{
    return state;
};

template<typename swe_getter>
Status infil_chooser::get(swe_getter&& swe)
{
    switch(status)
    {
        case Status::TO_FROZEN:
            return Status::FROZEN;

        case Status::TO_THAWED:
            return Status::THAWED;

        case Status::FROZEN:
            if (swe() > 0.0)
                return Status::FROZEN;
            else
                return Status::TO_THAWED;

        case Status::THAWED:
            if (swe() <= min_swe_to_freeze)
                return Status::THAWED;
            else
                return Status::TO_FROZEN;
    }
};

void infil_chooser::set_min_swe(const double val) {
    min_swe_to_freeze = val;
};
