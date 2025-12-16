#include "base_step.hpp"
#include "Soil.h"
#include <string_view>
#include <concepts>
#include <sparsehash/dense_hash_map>

namespace Ayers
{
    template<class T>
    concept AyersData = requires(T& t)
    {
        {t.rainfall()} -> std::floating_point;
        {t.snowmelt()} -> std::floating_point;
        {t.texture()} -> std::same_as<std::string_view>;
        {t.ground_cover()} -> std::same_as<std::string_view>;
        
        {t.runoff(submodules::test<double>)} -> std::same_as<void>;
        {t.snow_infiltrated(submodules::test<double>)} -> std::same_as<void>;
        {t.infiltrated(submodules::test<double>)} -> std::same_as<void>;
    };

    template<AyersData data>
    class Model : public submodules::base_step<Model<data>,data>
    {
    public:
        Model();
        ~Model() {};

        void execute_impl(data& d) const;

    private:    
        struct Result
        {
            double inf;
            double runoff;
        };
        Result do_rainfall(data& d, double) const;
        
        using InnerMap = Soil::string_map<double>;
        using MaxInfiltrationMap = Soil::string_map<Soil::string_map<double>>;//google::dense_hash_map<std::string,InnerMap>;
        static const MaxInfiltrationMap _max_infil; 
        static MaxInfiltrationMap initialize_max_infil_map();
        double max_infil_lookup(const std::string& outer_key, const std::string& inner_key);
    };

    template<AyersData data>
    void Model<data>::execute_impl(data& d) const
    {
        auto rainfall = d.rainfall();
        auto snowmelt = d.snowmelt();

        Result result;
        if (rainfall > 0.0)
        {
            result = do_rainfall(d,rainfall);
        }

        if (snowmelt > 0.0)
        {
            result.inf += snowmelt;
            d.snow_inf(snowmelt);
        }

        d.inf(result.inf);
        d.runoff(result.runoff);
    };

    template<AyersData data>
    Model<data>::Result Model<data>::do_rainfall(data& d,double rainfall) const
    {
        auto texture = d.texture();
        auto ground_cover = d.ground_cover();
        Result result;
        auto max_infiltration = _max_infil(texture,ground_cover);
        result.inf = std::min(max_infiltration,rainfall);

        result.runoff = rainfall - result.inf;
        if (result.runoff < 1e-12)
            result.runoff = 0.0;

        return result;
    };

    template<AyersData data>
    double Model<data>::max_infil_lookup(const std::string& outer_key,const std::string& inner_key)
    {
        auto inner_map = Soil::lookup(_max_infil,outer_key);
        return Soil::lookup(inner_map,inner_key);
    };

    template<AyersData data>
    const Model<data>::MaxInfiltrationMap Model<data>::_max_infil = 
        Model::initialize_max_infil_map();

    template<AyersData data>
    Model<data>::MaxInfiltrationMap Model<data>::initialize_max_infil_map()
    {
        std::string c1 = "bare_soil";
        std::string c2 = "row_crop";
        std::string c3 = "poor_pasture";
        std::string c4 = "small_grains";
        std::string c5 = "good_pasture";
        std::string c6 = "forested";
        
        struct Values
        {
            double bare_soil;
            double row_crop;
            double poor_pasture;
            double small_grains;
            double good_pasture;
            double forested;
        };

        auto set_inner = [](Values&& v) -> Model::InnerMap {
            Model::InnerMap inner_map;
            inner_map.set_empty_key("");

            inner_map["bare_soil"] = v.bare_soil;
            inner_map["row_crop"] = v.row_crop;
            inner_map["poor_pasture"] = v.poor_pasture;
            inner_map["small_grains"] = v.small_grains;
            inner_map["good_pasture"] = v.good_pasture;
            inner_map["forested"] = v.forested;

            return inner_map;
        };

        Model::MaxInfiltrationMap outer_map;
        outer_map.set_empty_key("");

        outer_map["coarse_over_coarse"] = 
            set_inner(Values{7.6,12.7,15.2,17.8,25.4,76.2});
        
        outer_map["medium_over_medium"] = 
            set_inner(Values{2.5,5.1,7.6,10.2,12.7,15.2});

        outer_map["medium_over_fine"] = 
            set_inner(Values{1.3,1.8,2.5,3.8,5.1,6.4});

        outer_map["fine_over_fine"] = 
            Soil::lookup(outer_map,"medium_over_fine");

        outer_map["soil_over_bedrock"] = 
            set_inner(Values{0.5,0.5,0.5,0.5,0.5,0.5});
        
        return outer_map;
    };
};
