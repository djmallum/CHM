#include "base_step.hpp"
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <utility>
#include "daily_accumulator.hpp"

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
        daily_accumulator daily_melt_total;
        daily_accumulator daily_rain_total;
        double daily_max_temp = 0.0;
        double current_inf = 0.0;
        double current_snow_inf = 0.0;
        double current_runoff = 0.0;
        double current_melt_runoff = 0.0;
        double soil_saturation_at_freeze;
                //bool end_freeze_tomorrow = false; 

    };

    enum class InfilPhase;
    

    InfilPhase determine_infiltration_phase(const State&,const double swe);
    double get_limited_inf(State&,const double swe);
    double calc_index_inf(State&);    
    void process_new_day(State&,const double swe);

    template<class T>
    concept CrackData = requires(T& t)
    {
        {t.get_state()} -> std::same_as<Crack::State&>;
        {t.is_newday()} -> std::convertible_to<const bool>;
        {t.snowmelt()} -> std::floating_point;
        {t.rainfall()} -> std::floating_point;
        {t.swe()} -> std::floating_point;
        {t.air_temperature()} -> std::floating_point;
        
        {t.infiltrated(std::declval<double>())} -> std::same_as<void>;
        {t.runoff(std::declval<double>())} -> std::same_as<void>;
        {t.snow_infiltrated(std::declval<double>())} -> std::same_as<void>;
        {t.melt_runoff(std::declval<double>())} -> std::same_as<void>;
        {t.rain_on_snow(std::declval<double>())} -> std::same_as<void>;
    };

    template<CrackData Data>
    class Model : public submodules::base_step<Model<Data>,Data>
    {
    public:

        void execute_impl(Data& d) const;

    };
    
    struct Params
    { 
        static inline double major_melt_threshold;
        static inline size_t infDays;
        static inline double lenstemp;
        static inline bool allow_early_inf;

        static void set(const double major_melt_threshold_, const size_t infDays_,
                const double lenstemp_, const bool allow_early_inf_)
        {
            major_melt_threshold = major_melt_threshold_;
            infDays = infDays_;
            lenstemp = lenstemp_;
            allow_early_inf = allow_early_inf_;
        };
    };
};


template<Crack::CrackData Data>
void Crack::Model<Data>::execute_impl(Data& d) const
{
    auto& s = d.get_state();

    auto is_newday = d.is_newday();
    
    if (!is_newday) [[likely]]
    {
        auto t = d.air_temperature();
        s.daily_max_temp = std::max(s.daily_max_temp, t);

        d.runoff(s.current_runoff);
        d.melt_runoff(s.current_melt_runoff);
        d.infiltrated(s.current_inf);
        d.snow_infiltrated(s.current_snow_inf);
        constexpr bool new_day = false;
        s.daily_rain_total.accumulate(new_day);
        s.daily_melt_total.accumulate(new_day);

        return;
    }
    
    const auto swe = d.swe();
    process_new_day(s,swe);
    
    s.daily_max_temp = d.air_temperature();

    auto yesterday_rain = s.daily_rain_total.get_yesterday();

    d.rain_on_snow(yesterday_rain);
    d.runoff(s.current_runoff);
    d.melt_runoff(s.current_melt_runoff);
    d.infiltrated(s.current_inf);
    d.snow_infiltrated(s.current_snow_inf);


};





