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
        double soil_saturation_at_freeze;
            
        struct Current
        {
            double inf = 0.0;
            double snow_inf = 0.0;
            double runoff = 0.0;
            double melt_runoff = 0.0;
            double rain = 0.0;
        } current;

        static inline size_t infDays = 6;
        static inline double major_melt_threshold = 25.0;
        static inline double lenstemp = -10.0;
        static inline bool allow_early_inf = true;
                //bool end_freeze_tomorrow = false; 

        State() {};
    private:
        /* private constructor for tests to wind up daily_melt_total and 
         * daily_rain_total */
        State(const double& melt, const double& rain);
    public:
        /* public factory for test constructor */
        static State construct_test_state(const double& melt, const double& rain)
        {
            return State(melt,rain);
        };

    };

    enum class LimitedPhase
    {
        FIRST_MAJOR, // First Major melt of the season
        RESET, // Reset index method
        NORMAL, // Compute infiltration from index
        PRIOR_INFILTRATION, // Before the first major melt
        NONE // Limited phase blocked
    };

    class phase_checker
    {
    public:
        LimitedPhase get(const State&, const double swe) const;
    };

    class limited_inf
    {
    public:
        double get(State&,const double swe) const;

    private:
        void set_index(State&,const double swe) const;
        double inf_from_index(State&) const;
        const phase_checker phase;

    };

    class inf_calculator
    {
    public:
        double limited(State&,const double swe) const;
        double unlimited(State&) const;
        double restricted(State&) const;
    private:
         const limited_inf inf;
    };

    class Processor
    {
    public:
        void new_day(State&,const double swe) const;

    private:
        inf_calculator inf_calc;
    };

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
        Processor process;
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

    s.daily_melt_total.accumulate(is_newday);
    s.daily_rain_total.accumulate(is_newday);
    
    const auto yesterday_rain = s.daily_rain_total.get_yesterday();

    if (!is_newday) [[likely]]
    {
        auto t = d.air_temperature();
        s.daily_max_temp = std::max(s.daily_max_temp, t);

        // Note confident on this rain_on_snow... shouldn't it be based on rain right now? Not yesterday?
        // Check with CRHM version.
        d.rain_on_snow(s.current.rain);
        d.runoff(s.current.runoff);
        d.melt_runoff(s.current.melt_runoff);
        d.infiltrated(s.current.inf);
        d.snow_infiltrated(s.current.snow_inf);
        return;
    }
    
    const auto swe = d.swe();
    const auto yesterday_melt = s.daily_melt_total.get_yesterday();
    process.new_day(s,swe);    

    s.daily_max_temp = d.air_temperature();

    d.rain_on_snow(s.current.rain);
    d.runoff(s.current.runoff);
    d.melt_runoff(s.current.melt_runoff);
    d.infiltrated(s.current.inf);
    d.snow_infiltrated(s.current.snow_inf);


};





