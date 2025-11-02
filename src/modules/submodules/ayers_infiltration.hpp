#pragma once

#include <concepts>
#include <string>
#include "base_step.hpp"

// Concept for data getters/setters
template<typename T, typename Data>
concept HasDataAccessors = requires(Data d) {
    { d.rainfall() } -> (std::floating_point<> || std::integral<>);
    { d.snowmelt() } -> (std::floating_point<> || std::integral<>);
    { d.texture() } -> std::same_as<std::string>;
    { d.ground_cover() } -> std::same_as<std::string>;
    { d.runoff(std::declval<double>()) } -> std::same_as<void>;
    { d.inf(std::declval<double>()) } -> std::same_as<void>;
    { d.snow_inf(std::declval<double>()) } -> std::same_as<void>;
};

// Container for ayers soil infiltration measurements
// Based on texture and ground cover.
class Ayers_texture_lookup
{
public:
    Ayers_texture_lookup()
    {
        _texture_map.set_empty_key("");
        
        std::string c1 = "bare_soil";
        std::string c2 = "row_crop";
        std::string c3 = "poor_pasture";
        std::string c4 = "small_grains";
        std::string c5 = "good_pasture";
        std::string c6 = "forested";
        
        google::dense_hash_map<std::string, double> inner_map;
        inner_map.set_empty_key("");

        // coarse_over_coarse
        inner_map[c1] = 7.6;
        inner_map[c2] = 12.7;
        inner_map[c3] = 15.2;
        inner_map[c4] = 17.8;
        inner_map[c5] = 25.4;
        inner_map[c6] = 76.2;
        _texture_map["coarse_over_coarse"] = inner_map;

        // medium over medium
        inner_map.clear();
        inner_map.set_empty_key("");
        inner_map[c1] = 2.5;
        inner_map[c2] = 5.1;
        inner_map[c3] = 7.6;
        inner_map[c4] = 10.2;
        inner_map[c5] = 12.7;
        inner_map[c6] = 15.2;
        _texture_map["medium_over_medium"] = inner_map;

        // medium/fine over fine        
        inner_map.clear();
        inner_map.set_empty_key("");
        inner_map[c1] = 1.3;
        inner_map[c2] = 1.8;
        inner_map[c3] = 2.5;
        inner_map[c4] = 3.8;
        inner_map[c5] = 5.1;
        inner_map[c6] = 6.4;
        _texture_map["medium_over_fine"] = inner_map;
        _texture_map["fine_over_fine"] = inner_map;

        // soil over shallow bedrock
        inner_map.clear();
        inner_map.set_empty_key("");
        inner_map[c1] = 0.5;
        inner_map[c2] = 0.5;
        inner_map[c3] = 0.5;
        inner_map[c4] = 0.5;
        inner_map[c5] = 0.5;
        inner_map[c6] = 0.5;
        _texture_map["soil_over_bedrock"] = inner_map;
    }

    double get_max_infiltration(const std::string& texture, const std::string& ground_cover) const
    {
        auto texture_it = _texture_map.find(texture);
        if (texture_it != _texture_map.end()) {
            auto cover_it = texture_it->second.find(ground_cover);
            if (cover_it != texture_it->second.end()) {
                return cover_it->second;
            }
        }
        return 0.0; // or throw an exception, or return a default
    }

private:
    google::dense_hash_map<std::string, google::dense_hash_map<std::string, double>> _texture_map;
};



template<typename Data>
requires HasDataAccessors<double, Data>
class Ayers : public base_step<Data>
{
public:
    Ayers() = default;
    ~Ayers() = default;

    void execute(Data& d) override
    {
        auto rainfall = d.rainfall();
        auto snowmelt = d.snowmelt();
        auto texture = d.texture();
        auto ground_cover = d.ground_cover();

        if (rainfall == 0.0 && snowmelt == 0.0) 
        {
            d.runoff(0.0);
            d.inf(0.0);
            d.snow_inf(0.0);
            return;
        }
        
        double runoff = 0.0;
        double inf = 0.0;
        double snow_inf = 0.0;

        if (rainfall > 0.0)    
        {
            auto maxinfil = _texture_lookup.get_max_infiltration(texture, ground_cover);
            inf = std::min(maxinfil, rainfall);
            runoff = rainfall - inf;        
            if (runoff < 1e-12)
                runoff = 0.0;  
        }

        if (snowmelt > 0.0)    
        {
            inf += snowmelt;   
            snow_inf = snowmelt; 
        }
        
        d.runoff(runoff);
        d.inf(inf);
        d.snow_inf(snow_inf);
    }

private:
    Ayers_texture_lookup _texture_lookup;
};
