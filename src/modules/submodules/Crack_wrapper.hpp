#pragma once
#include "Crack.hpp"

// Extended concept for wrapper
template<typename Data>
concept HasCrackWrapperDataAccessors = HasCrackDataAccessors<Data> && requires(Data d) {
    // Additional methods needed by wrapper
    { d.is_new_day() } -> std::same_as<bool>;
    { d.dt() } -> std::floating_point;
    { d.soil_saturation_at_freeze_available() } -> std::same_as<bool>;
    { d.get_soil_saturation() } -> std::floating_point;
    { d.set_soil_saturation_at_freeze(std::declval<double>()) } -> std::same_as<void>;
    { d.set_steps_per_day(std::declval<double>()) } -> std::same_as<void>;
    { d.set_final_runoff(std::declval<double>()) } -> std::same_as<void>;
    { d.set_final_melt_runoff(std::declval<double>()) } -> std::same_as<void>;
    { d.set_final_inf(std::declval<double>()) } -> std::same_as<void>;
    { d.set_final_snowinf(std::declval<double>()) } -> std::same_as<void>;
    { d.set_final_rain_on_snow(std::declval<double>()) } -> std::same_as<void>;
    { d.reset_soil_saturation_at_freeze() } -> std::same_as<void>;
};

template<typename Data>
class CrackWrapper {
public:
    CrackWrapper() = default;
    
    void execute_if_frozen(Data& d) {
        auto& crack_info = d.crack_info();
        
        // Check if we should start the crack model
        if (d.swe() > d.min_swe_to_freeze() && !crack_info.frozen && d.is_new_day()) {
            crack_info.begin_freeze();
            crack_info.end_freeze_tomorrow = false;
            
            // Set soil saturation if not already set
            if (!d.soil_saturation_at_freeze_available()) {
                d.set_soil_saturation_at_freeze(d.get_soil_saturation());
            }
        }
        
        // Run crack model if frozen
        if (crack_info.frozen) {
            double steps_per_day = 86400.0 / d.dt();
            d.set_steps_per_day(steps_per_day);
            
            // Update daily melt total (scaled by steps per day)
            crack_info.daily_melt_total = d.snowmelt() * steps_per_day;
            
            // Execute crack model
            crack_model.execute(d);
            
            // Scale outputs back
            double runoff = d.runoff() / steps_per_day;
            double melt_runoff = d.melt_runoff() / steps_per_day;
            double inf = d.inf() / steps_per_day;
            double snowinf = d.snow_inf() / steps_per_day;
            double rain_on_snow = d.rain_on_snow();
            
            // Update final outputs
            d.set_final_runoff(runoff);
            d.set_final_melt_runoff(melt_runoff);
            d.set_final_inf(inf);
            d.set_final_snowinf(snowinf);
            d.set_final_rain_on_snow(rain_on_snow);
            
            // Check if we should end freeze
            if (d.is_new_day() && d.swe() <= 0.0 && crack_info.major_melt_count > 0) {
                crack_info.end_freeze();
                crack_info.end_freeze_tomorrow = true;
                d.reset_soil_saturation_at_freeze();
            }
        }
    }

private:
    Crack<Data> crack_model;
};
