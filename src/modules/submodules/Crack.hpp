#pragma once
#include "base_step.hpp"
#include <concepts>
#include <spdlog/spdlog.h>

// Concept for Crack data accessors
template<typename Data>
concept HasCrackDataAccessors = requires(Data d) {
    // Input getters
    { d.snowmelt() } -> std::floating_point<>;
    { d.rainfall() } -> std::floating_point<>;
    { d.swe() } -> std::floating_point<>;
    { d.soil_saturation_at_freeze() } -> std::floating_point<>;
    { d.airtemp() } -> std::floating_point<>;
    { d.is_newday() } -> std::same_as<bool>;
    
    // Parameter getters
    { d.major() } -> std::floating_point<>;
    { d.min_swe_to_freeze() } -> std::floating_point<>;
    { d.infDays() } -> std::integral<>;
    { d.AllowPriorInf() } -> std::same_as<bool>;
    { d.lenstemp() } -> std::floating_point<>;
    { d.steps_per_day() } -> std::floating_point<>;
    
    // Output setters
    { d.runoff(std::declval<double>()) } -> std::same_as<void>;
    { d.melt_runoff(std::declval<double>()) } -> std::same_as<void>;
    { d.inf(std::declval<double>()) } -> std::same_as<void>;
    { d.snow_inf(std::declval<double>()) } -> std::same_as<void>;
    { d.rain_on_snow(std::declval<double>()) } -> std::same_as<void>;
    
    // State accessors
    { d.crack_info() } -> std::same_as<typename Crack<Data>::info&>;
};

template<typename Data>
requires HasCrackDataAccessors<Data>
class Crack : public base_step<Data> {
public:
    struct info {
        bool frozen = false;
        unsigned int major_melt_count = 0;
        double index = 0.0;
        double max_major_per_melt = 0.0;
        double init_SWE = 0.0;
        double daily_melt_total = 0.0;
        double daily_rain_total = 0.0;
        bool current_day_is_major = false;
        double tmax = 0.0;
        double current_inf = 0.0;
        double current_snow_inf = 0.0;
        double current_runoff = 0.0;
        double current_melt_runoff = 0.0;
        double yesterday_melt = 0.0;
        bool end_freeze_tomorrow = false;
        
        void init() {
            frozen = false;
            major_melt_count = 0;
            index = 0.0;
            max_major_per_melt = 0.0;
            init_SWE = 0.0;
            daily_melt_total = 0.0;
            daily_rain_total = 0.0;
            current_day_is_major = false;
            tmax = 0.0;
            current_inf = 0.0;
            current_snow_inf = 0.0;
            current_runoff = 0.0;
            current_melt_runoff = 0.0;
            yesterday_melt = 0.0;
            end_freeze_tomorrow = false;
        }

        void begin_freeze() {
            frozen = true;
            index = 0.0;
            max_major_per_melt = 0.0;
            init_SWE = 0.0;
            current_inf = 0.0;
            current_runoff = 0.0;
            yesterday_melt = 0.0;
            end_freeze_tomorrow = false;
        }

        void end_freeze() {
            frozen = false;
            major_melt_count = 0;
        }
    };

    Crack() = default;
    ~Crack() = default;

    void execute(Data& d) override {
        auto& crack_info = d.crack_info();
        
        crack_info.daily_rain_total += d.rainfall();
        
        if (d.is_newday()) {
            process_new_day(d, crack_info);
        } else {
            carry_over_previous_values(d, crack_info);
        }
        
        if (!d.is_CRHM_compare_test()) {
            crack_info.daily_melt_total += d.snowmelt();
        }
        update_t_max(d, crack_info);
    }

private:
    void process_new_day(Data& d, info& crack_info) {
        if (crack_info.daily_melt_total > 0.0) {
            double soil_sat = d.soil_saturation_at_freeze();
            
            if (soil_sat == 0) { // Unlimited
                process_unlimited_case(d, crack_info);
            } else if (soil_sat > 0 && soil_sat < 100) { // Limited
                process_limited_case(d, crack_info);
            } else if (soil_sat == 100) { // Restricted
                process_restricted_case(d, crack_info);
            }
            
            // Update rain handling
            update_rain_handling(d, crack_info);
            
            // Store current state
            update_state_variables(d, crack_info);
        }
        
        // Reset daily totals
        reset_daily_totals(crack_info);
    }

    void process_unlimited_case(Data& d, info& crack_info) {
        d.inf(crack_info.daily_melt_total);
        crack_info.major_melt_count = 1;
    }

    void process_limited_case(Data& d, info& crack_info) {
        increment_major_count(d, crack_info);
        check_for_ice_lens(d, crack_info);

        if (is_first_major(d, crack_info)) {
            SPDLOG_DEBUG("First Major");
            calc_index(d, crack_info);
            calc_actual_inf(d, crack_info);
        } else if (is_limited_phase(crack_info)) {
            SPDLOG_DEBUG("Limited Phase");
            calc_actual_inf(d, crack_info);
        } else if (is_prior_first_major(crack_info) && d.AllowPriorInf()) {
            SPDLOG_DEBUG("Prior");
            d.inf(crack_info.daily_melt_total);
        }
    }

    void process_restricted_case(Data& d, info& crack_info) {
        d.inf(0.0);
        crack_info.major_melt_count = 1;
    }

    void update_rain_handling(Data& d, info& crack_info) {
        double runoff = crack_info.daily_melt_total - d.inf();
        
        if (d.inf() > 0.0) {
            d.inf(d.inf() + crack_info.daily_rain_total);
        } else {
            runoff += crack_info.daily_rain_total;
        }
        
        d.runoff(runoff);
        d.melt_runoff(runoff);
        d.snow_inf(d.inf());
    }

    void update_state_variables(Data& d, info& crack_info) {
        crack_info.yesterday_melt = crack_info.daily_melt_total;
        crack_info.current_inf = d.inf();
        crack_info.current_snow_inf = d.snow_inf();
        crack_info.current_runoff = d.runoff();
        crack_info.current_melt_runoff = d.melt_runoff();
        d.rain_on_snow(crack_info.daily_rain_total);
    }

    void reset_daily_totals(info& crack_info) {
        crack_info.daily_melt_total = 0.0;
        crack_info.daily_rain_total = 0.0;
    }

    void carry_over_previous_values(Data& d, info& crack_info) {
        d.runoff(crack_info.current_runoff);
        d.melt_runoff(crack_info.current_melt_runoff);
        d.inf(crack_info.current_inf);
        d.snow_inf(crack_info.current_snow_inf);
    }

    void calc_index(Data& d, info& crack_info) {
        double soil_sat = d.soil_saturation_at_freeze();
        crack_info.index = 5 * (1 - soil_sat/100.0) * std::pow(d.swe(), 0.584);
        crack_info.max_major_per_melt = crack_info.index / d.infDays();
        crack_info.index = std::min(crack_info.index / d.swe(), 1.0);
        crack_info.init_SWE = d.swe();
    }

    void calc_actual_inf(Data& d, info& crack_info) {
        double calculated_inf = crack_info.daily_melt_total * crack_info.index;
        if (calculated_inf > crack_info.max_major_per_melt && is_major_melt(crack_info, d.major())) {
            calculated_inf = crack_info.max_major_per_melt;
        }
        d.inf(calculated_inf);
    }

    void update_t_max(Data& d, info& crack_info) {
        if (d.is_newday()) {
            crack_info.tmax = d.airtemp();
        } else {
            crack_info.tmax = std::max(crack_info.tmax, d.airtemp());
        }
    }

    void check_for_ice_lens(Data& d, info& crack_info) {
        if (crack_info.major_melt_count > 0 && crack_info.tmax < d.lenstemp()) {
            SPDLOG_DEBUG("Ice lens found");
            crack_info.major_melt_count = d.infDays() + 4;
        }
        crack_info.tmax = 0.0;
    }

    bool is_first_major(Data& d, info& crack_info) {
        return is_major_melt(crack_info, d.major()) && d.swe() >= crack_info.init_SWE && 
               (is_prior_first_major(crack_info) || is_limited_phase(crack_info));
    }

    bool is_major_melt(info& crack_info, double major_threshold) {
        return crack_info.daily_melt_total > major_threshold;
    }

    void increment_major_count(Data& d, info& crack_info) {
        if (is_major_melt(crack_info, d.major())) {
            crack_info.major_melt_count++;
        }
    }

    bool is_limited_phase(info& crack_info) {
        return crack_info.major_melt_count > 0 && crack_info.major_melt_count <= infDays;
    }

    bool is_prior_first_major(info& crack_info) {
        return crack_info.major_melt_count == 0;
    }
};
