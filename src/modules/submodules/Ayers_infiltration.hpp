#include "base_step.hpp"
#include "Soil.h"
#include <concepts>
#include <sparsehash/dense_hash_map>

namespace Ayers
{

    /* Convenience types for infiltration map */
    using InnerMap = Soil::string_map<double>;
    using MaxInfiltrationMap = Soil::string_map<InnerMap>;

    /* Convenience struct for outputs */
    struct Result
    {
        double inf;
        double runoff;
    };

    /* Concept to enforce getters/setters required by Ayers::Model */
    template<class T>
    concept AyersData = requires(T& t)
    {
        {t.rainfall()} -> std::floating_point;
        {t.snowmelt()} -> std::floating_point;
        {t.texture()} -> std::same_as<const std::string&>;
        {t.ground_cover()} -> std::same_as<const std::string&>;
        
        {t.runoff(submodules::test<double>)} -> std::same_as<void>;
        {t.snow_infiltrated(submodules::test<double>)} -> std::same_as<void>;
        {t.infiltrated(submodules::test<double>)} -> std::same_as<void>;
    };

    /* Class containing the ayers infiltration table, stored as a MaxInfiltrationMap 
     * Uses lazy initialization so every instance gets a reference to the same copy */
    class infiltration_map
    {
    public:
        double max_infil_lookup(const std::string& outer_key,const std::string& inner_key) const;

        infiltration_map();
        ~infiltration_map() {};
    private:
        const MaxInfiltrationMap& _max_infil;
        const MaxInfiltrationMap& get_max_infil_map();

    };

    /* Engine of the Ayers::Model class. Converts rainfall into Result instance */
    class rain_infiltration
    {
    public:
        Result compute(const std::string& texture, const std::string& ground_cover, const double rainfall) const;
        rain_infiltration() {};
    private:
        const infiltration_map _max_infil;
    };

    /* Actual Model. Stores instance of rain_infiltration privately, passes data from data object to members and calculations */
    template<AyersData data>
    class Model : public submodules::base_step<Model<data>,data>
    {
    public:
        Model() {};
        ~Model() {};

        void execute_impl(data& d) const;

    private:    

        rain_infiltration rain_inf;
         
    };
    
    /* public interface of Model class */
    template<AyersData data>
    void Model<data>::execute_impl(data& d) const 
    {
        auto rainfall = d.rainfall();
        auto snowmelt = d.snowmelt();

        Result result;
        if (rainfall > 0.0)
        {
            const auto& texture = d.texture();
            const auto& ground_cover = d.ground_cover();
            result = rain_inf.compute(texture,ground_cover,rainfall);
        }

        if (snowmelt > 0.0)
        {
            result.inf += snowmelt;
            d.snow_infiltrated(snowmelt);
        }

        d.infiltrated(result.inf);
        d.runoff(result.runoff);
    };

};
