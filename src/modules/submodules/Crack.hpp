#include "base_step.hpp"
#include <concepts>
#include <cstddef>

namespace Crack
{
    struct State
    {
        // Commented out are part of the module state, not the submodule state, Separate!
        // bool frozen = false;
        size_t major_melt_count = 0;
        double index = 0.0;
        double max_major_per_melt = 0.0;
        double init_SWE = 0.0;
        double daily_melt_total = 0.0;
        double daily_rain_total = 0.0;
        double tmax = 0.0;
        double current_inf = 0.0;
        double current_snow_inf = 0.0;
        double current_runoff = 0.0;
        double current_melt_runoff = 0.0;
        double yesterday_melt = 0.0;
        double soil_saturation_at_freeze;
        //bool end_freeze_tomorrow = false; 
    };

    template<class T>
    concept CrackSetUpData = requires(T& t)
    {
        {t.get_state()} -> std::same_as<Crack::State&>;
        {t.is_new_day()} -> std::convertible_to<const bool>;
        {t.snowmelt()} -> std::floating_point;
        {t.rainfall()} -> std::floating_point;
        {t.swe()} -> std::floating_point;
        {t.air_temperature()} -> std::floating_point;
        
        {t.infiltrated(std::declval<double>())} -> std::same_as<void>;
        {t.runoff(std::declval<double>())} -> std::same_as<void>;
        {t.snow_infiltrated(std::declval<double>())} -> std::same_as<void>;
        {t.melt_runoff(std::declval<double>())} -> std::same_as<void>;
    };

    template<typename Data>
    class Model : public submodules::base_step<Model<Data>,Data>
    {
    public:

        void execute_impl(Data& d);
        
    };
};



