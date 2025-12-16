#include "base_step.hpp"
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <stdexcept>
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

        
    private:
        void process_new_day(Data&,State&) const;
        double get_limited_inf(Data&,State&) const;
        double calc_index_inf(State&) const;

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
        s.daily_rain_total.accumulate(is_newday);
        s.daily_melt_total.accumulate(is_newday);

        return;
    }

    process_new_day(d,s);

};

template<Crack::CrackData Data>
void Crack::Model<Data>::process_new_day(Data& d, State& s) const
{
    auto inf = 0.0;
    auto runoff = 0.0;
    auto melt_runoff = 0.0;
    auto snow_inf = 0.0;
    
    constexpr bool new_day = true;
    s.daily_rain_total.accumulate(new_day);
    s.daily_melt_total.accumulate(new_day);
    
    double yesterday_melt = s.daily_melt_total.get_yesterday();
    double yesterday_rain = s.daily_rain_total.get_yesterday();
    // TODO Profile, consider removing unpredictable branch
    if (yesterday_melt > 0.0)
    {
        if (s.soil_saturation_at_freeze > 0.0 &&
                s.soil_saturation_at_freeze < 100.0 ) [[likely]]
        {
            inf = get_limited_inf(d,s);
        }
        else [[unlikely]]
        {
            if (s.soil_saturation_at_freeze == 0.0)
            {
                inf = yesterday_melt;
                s.major_melt_count = 1;
            }
            else if (s.soil_saturation_at_freeze == 100.0)
            {
                s.major_melt_count = 0;
            }
            else
                throw std::logic_error("soil_saturation_at_freeze not bounded between 0 and 100");
        }

        runoff = yesterday_melt - inf;
        // Moved the following lines before the if, because daily_rain_total shouldn't 
        // be added tomelt_runoff or snow_inf (which track only melt related quantities, not rain).
        melt_runoff = runoff;
        snow_inf = inf;
        if (inf > 0.0)
            inf += yesterday_rain;
        else
            runoff += yesterday_rain;

    } // if

    s.current_inf = inf;
    s.current_snow_inf = snow_inf;
    s.current_runoff = runoff;
    s.current_melt_runoff = melt_runoff;

    
    // has to be new day, therefore, maximum temperature is the current temperature 
    // for this new day
    s.daily_max_temp = d.air_temperature();

    d.rain_on_snow(yesterday_rain);
    d.runoff(s.current_runoff);
    d.melt_runoff(s.current_melt_runoff);
    d.infiltrated(s.current_inf);
    d.snow_infiltrated(s.current_snow_inf);
    
    

    
};

enum class Crack::InfilPhase
{
    FIRST_MAJOR,
    RESET,
    LIMITED_PHASE,
    PRIOR_INFILTRATION,
    NONE
};

template<Crack::CrackData Data>
double Crack::Model<Data>::get_limited_inf(Data& d,State& s) const
{
    const auto swe = d.swe();

    InfilPhase phase = determine_infiltration_phase(s,swe);

    switch(phase)
    {
        // returns behave as "breaks;" commands
        case InfilPhase::FIRST_MAJOR:            
            // FIRST_MAJOR and RESET have same outcome
        case InfilPhase::RESET:
            s.index = 5 * (1 - s.soil_saturation_at_freeze/100.0) * std::pow(swe,0.584);
            
            s.max_major_per_melt = s.index / Params::infDays;
            s.index = std::min(s.index/swe,1.0);
            s.init_SWE = swe;
            // RESET is same as LIMITED_PHASE but with a calculation of index
            // No return so go to next
        case InfilPhase::LIMITED_PHASE:
            return calc_index_inf(s);
        case InfilPhase::PRIOR_INFILTRATION:
            return s.daily_melt_total.get_yesterday(); 
        case InfilPhase::NONE:
            return 0.0;
        default:
            throw std::logic_error("Must be a phase from enum InfilPhase");
    }    
};

template<Crack::CrackData Data>
double Crack::Model<Data>::calc_index_inf(State& s) const
{
    return std::min(s.daily_melt_total.get_yesterday() * s.index,
            s.max_major_per_melt);
};
