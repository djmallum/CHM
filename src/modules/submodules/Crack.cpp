#include "Crack.hpp"
#include <stdexcept>

namespace Crack
{
    enum class InfilPhase
    {
        FIRST_MAJOR,
        RESET,
        LIMITED_PHASE,
        PRIOR_INFILTRATION,
        NONE
    };

    InfilPhase determine_infiltration_phase(const Crack::State& s,const double swe)
    {
        const bool count_nonzero = s.major_melt_count > 0;
        const bool has_ice_lens = count_nonzero && s.daily_max_temp < Params::lenstemp
            && s.soil_saturation_at_freeze > 0.0;

        if (has_ice_lens)
        {   
            State& s_temp = const_cast<State&>(s);
            s_temp.major_melt_count = Params::infDays + 4;
            return InfilPhase::NONE;
        }

        const bool before_first_major = s.major_melt_count == 0;
        const bool is_major_melt = s.daily_melt_total.get_yesterday() > Params::major_melt_threshold;

        if (before_first_major)
        {
            if ( is_major_melt && swe > 0.0) 
                return InfilPhase::FIRST_MAJOR;
            else   
                return InfilPhase::PRIOR_INFILTRATION;
        }

        const bool count_below_max = s.major_melt_count <= Params::infDays;

        if (count_nonzero && count_below_max)
        {
            const bool reset_index = is_major_melt
            && swe > s.init_SWE;
            
            if (reset_index)
                return InfilPhase::RESET;
                
            return InfilPhase::LIMITED_PHASE;
        }

        return InfilPhase::NONE; 
        
    };

    double get_limited_inf(State& s,const double swe)
    {

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

    double calc_index_inf(State& s)
    {
        return std::min(s.daily_melt_total.get_yesterday() * s.index,
                s.max_major_per_melt);
    };

    void process_new_day(State& s, const double swe)
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
                inf = get_limited_inf(s,swe);
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

        
        
    };
};
