#include "Crack.hpp"
#include <stdexcept>

namespace Crack
{

    State::State(const double& melt, const double& rain)
    {
        daily_melt_total.bind_target(melt);
        daily_rain_total.bind_target(rain);

        daily_melt_total.accumulate(false);
        daily_rain_total.accumulate(false);

        daily_melt_total.accumulate(true);
        daily_rain_total.accumulate(true);
    };
    
    LimitedPhase phase_checker::get(const Crack::State& s,const double swe) const
    {
        const bool count_nonzero = s.major_melt_count > 0;
        const bool has_ice_lens = count_nonzero && s.daily_max_temp < State::lenstemp
            && s.soil_saturation_at_freeze > 0.0;

        if (has_ice_lens)
        {   
            State& s_temp = const_cast<State&>(s);
            s_temp.major_melt_count = State::infDays + 4;
            return LimitedPhase::NONE;
        }

        const bool before_first_major = s.major_melt_count == 0;
        const bool is_major_melt = s.daily_melt_total.get_yesterday() > State::major_melt_threshold;

        if (before_first_major)
        {
            if (is_major_melt) 
                return LimitedPhase::FIRST_MAJOR;
            else   
                return LimitedPhase::PRIOR_INFILTRATION;
        }

        const bool count_below_max = s.major_melt_count <= State::infDays;

        if (count_nonzero && count_below_max)
        {
            const bool reset_index = is_major_melt
            && swe > s.init_SWE;
            
            if (reset_index)
                return LimitedPhase::RESET;
                
            return LimitedPhase::NORMAL;
        }

        return LimitedPhase::NONE; 
        
    };

    double limited_inf::get(State& s,const double swe) const
    {
        switch(phase.get(s,swe))
        {
            // returns behave as "breaks;" commands
            case LimitedPhase::FIRST_MAJOR:            
                // ensure major_melt_count starts at 0, will be incremented in RESET
                // FIRST_MAJOR and RESET have same outcome otherwise
                s.major_melt_count = 0;
            case LimitedPhase::RESET:
                set_index(s,swe); 
                // RESET is same as NORMAL but with a calculation of index
                // No return so go to next
                ++s.major_melt_count;
            case LimitedPhase::NORMAL:
                return inf_from_index(s);
            case LimitedPhase::PRIOR_INFILTRATION:
                return s.daily_melt_total.get_yesterday(); 
            case LimitedPhase::NONE:
                return 0.0;
        }    
    };

    void limited_inf::set_index(State& s, const double swe) const
    {
        s.index = 5 * (1 - s.soil_saturation_at_freeze/100.0) * std::pow(swe,0.584);
                
        s.max_major_per_melt = s.index / State::infDays;
        s.index = std::min(s.index/swe,1.0);
        s.init_SWE = swe;
    };

    double limited_inf::inf_from_index(State& s) const
    {
        return std::min(s.daily_melt_total.get_yesterday() * s.index,
                s.max_major_per_melt);
    };

    double inf_calculator::limited(State& s, const double swe) const
    {
        return inf.get(s,swe);
    };

    double inf_calculator::restricted(State& s) const
    {
        s.major_melt_count = 0;
        return 0.0;
    };
    
    double inf_calculator::unlimited(State& s) const
    {
        s.major_melt_count = 1;
        return s.daily_melt_total.get_yesterday();
    };
    
    void Processor::new_day(State& s, const double swe) const
    {
        auto inf = 0.0;
        auto runoff = 0.0;
        auto melt_runoff = 0.0;
        auto snow_inf = 0.0;
        const auto rain_yesterday = s.daily_rain_total.get_yesterday();
        
        // TODO Profile, consider removing unpredictable branch
        if (s.daily_melt_total.get_yesterday() > 0.0)
        {
            if (s.soil_saturation_at_freeze > 0.0 &&
                    s.soil_saturation_at_freeze < 100.0 ) [[likely]]
            {
                inf = inf_calc.limited(s, swe);
            }
            else [[unlikely]]
            {
                if (s.soil_saturation_at_freeze == 0.0)
                {
                    inf = inf_calc.unlimited(s);
                }
                else if (s.soil_saturation_at_freeze == 100.0)
                {
                    inf = inf_calc.restricted(s);
                }
                else
                    throw std::logic_error("soil_saturation_at_freeze not bounded between 0 and 100");
            }

            runoff = s.daily_melt_total.get_yesterday() - inf;
            // Moved the following lines before the if, because daily_rain_total shouldn't 
            // be added tomelt_runoff or snow_inf (which track only melt related quantities, not rain).
            melt_runoff = runoff;
            snow_inf = inf;
            if (inf > 0.0)
                inf += rain_yesterday;
            else
                runoff += rain_yesterday;

            s.current.rain = 0.0;

        } // if

        s.current.inf = inf;
        s.current.snow_inf = snow_inf;
        s.current.runoff = runoff;
        s.current.melt_runoff = melt_runoff;
        s.current.rain += rain_yesterday;

    };
};
