#include "Crack.hpp"

namespace Crack
{
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

};
