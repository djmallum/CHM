#include "Ayers_infiltration.hpp"

namespace Ayers
{
    static const MaxInfiltrationMap initialize_max_infil_map()
    {
        struct Values
        {
            double bare_soil;
            double row_crop;
            double poor_pasture;
            double small_grains;
            double good_pasture;
            double forested;
        };

        auto set_inner = [](Values&& v) -> InnerMap {
            InnerMap inner_map;
            inner_map.set_empty_key("");

            constexpr std::string bare_soil = "bare_soil";
            constexpr std::string row_crop = "row_crop";
            constexpr std::string poor_pasture = "poor_pasture";
            constexpr std::string small_grains = "small_grains";
            constexpr std::string good_pasture = "good_pasture";
            constexpr std::string forested = "forested";

            inner_map[bare_soil] = v.bare_soil;
            inner_map[row_crop] = v.row_crop;
            inner_map[poor_pasture] = v.poor_pasture;
            inner_map[small_grains] = v.small_grains;
            inner_map[good_pasture] = v.good_pasture;
            inner_map[forested] = v.forested;

            return inner_map;
        };

        constexpr std::string coarse_over_coarse = "coarse_over_coarse";
        constexpr std::string medium_over_medium = "medium_over_medium";
        constexpr std::string medium_over_fine = "medium_over_fine";
        constexpr std::string fine_over_fine = "fine_over_fine";
        constexpr std::string soil_over_bedrock = "soil_over_bedrock";

        MaxInfiltrationMap outer_map;
        outer_map.set_empty_key("");

        outer_map[coarse_over_coarse] = 
            set_inner(Values{7.6,12.7,15.2,17.8,25.4,76.2});
        
        outer_map[medium_over_medium] = 
            set_inner(Values{2.5,5.1,7.6,10.2,12.7,15.2});

        outer_map[medium_over_fine] = 
            set_inner(Values{1.3,1.8,2.5,3.8,5.1,6.4});

        outer_map[fine_over_fine] = outer_map[medium_over_fine]; 

        outer_map[soil_over_bedrock] = 
            set_inner(Values{0.5,0.5,0.5,0.5,0.5,0.5});
        
        return outer_map;
    };

    infiltration_map::infiltration_map() : _max_infil(get_max_infil_map()) {};

    const MaxInfiltrationMap& infiltration_map::get_max_infil_map()
    {
        static const MaxInfiltrationMap value = initialize_max_infil_map();
        return value;
    };

    double infiltration_map::max_infil_lookup(const std::string& outer_key,const std::string& inner_key) const
    {
        auto inner_map = Soil::lookup(_max_infil,outer_key);
        return Soil::lookup(inner_map,inner_key);
    };

    Result rain_infiltration::compute(const std::string& texture, const std::string& ground_cover, const double rainfall) const
    {
        Result result;
        auto max_infiltration = _max_infil.max_infil_lookup(texture,ground_cover);
        result.inf = std::min(max_infiltration,rainfall);

        result.runoff = rainfall - result.inf;
        if (result.runoff < 1e-12)
            result.runoff = 0.0;

        return result;
    };
};
