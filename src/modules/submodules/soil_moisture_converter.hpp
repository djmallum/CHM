#pragma once
#include "base_step.hpp"
#include <concepts>
#include <utility>

template<typename T>
concept SoilMoistureConverterData = requires(T& t) {
    { t.storage_is_total_moisture() } -> std::convertible_to<bool>;
    { t.fractional_cutoff() } -> std::convertible_to<double>;
    { t.soil_storage() } -> std::convertible_to<double>;
    { t.soil_storage_max() } -> std::convertible_to<double>;
    { t.porosity() } -> std::convertible_to<double>;
    { t.volumetric_moisture_content(std::declval<const double>()) } -> std::same_as<void>;
    { t.volumetric_moisture_content() } -> std::convertible_to<double>;
    { t.saturation(std::declval<const double>()) } -> std::same_as<void>;
};

template<SoilMoistureConverterData data>
class soil_moisture_converter : public base_step<data>
{
public:
    explicit soil_moisture_converter() {};
    ~soil_moisture_converter() {};

    void execute(data& d) override final;

    // Currently constants (e.g., wilt_point) and inputs (e.g., CurrentDay) are indistinguishable
    // TODO Consider in the future if having them be different is necessary.
private:
    double total_volumetric_moisture(data& d) const;
    double fractional_volumetric_moisture(const double lower_bound_fraction, data& d) const;
    void check_cutoff_validity(double c) const;
};

template<SoilMoistureConverterData data>
void soil_moisture_converter<data>::execute(data& d)
{
    double volumetric_moisture;

    const double cutoff = d.fractional_cutoff();

    check_cutoff_validity(cutoff);
   
    if (d.storage_is_total_moisture())
        volumetric_moisture = total_volumetric_moisture(d); 
    else
        volumetric_moisture = fractional_volumetric_moisture(cutoff,d);
    
    d.volumetric_moisture_content(volumetric_moisture);   

    d.saturation(volumetric_moisture / d.porosity()); 
};

template<SoilMoistureConverterData data>
void soil_moisture_converter<data>::check_cutoff_validity(double c) const
{
    if (c < 0.0 || c > 1.0)
        throw std::logic_error("cutoff must be between (inclusive) 0 and 1");
};

/*
 * It's common for the soil moisture storage (a depth) to not actually include all moisture in the soil in
 * the variable. Depending on the processes, not all moisture in all pores can be included in that process. 
 * For example, transpiration does not occur below the wilt point and so if the important process is not 
 * about transpiration or there is no transpiration, the soil storage might only track moisture above the 
 * wilt point.
 *
 * Therefore there are two functions to compute the volumetric moisture content:
 *
 * 1. total_volumetric_moisture(data& d);
 * 
 * This function computes assuming that ALL moisture in the soil is counted in the soil storage varaible 
 * (here we mean d.soil_storage()).
 *
 * 2. fractional_volumetric_moisture(double lower_bound_fraction, data& d);
 *
 * This function assumes that a percentage of the soil moisture content is always full and never used by the
 * calculation that computes d.soil_storage() and accounts for that to compute the actual volumetric soil 
 * moisture content. lower_bound_fraction could be the field capacity or wilt point. User who writes the data
 * class will use their knowledge of where the d.soil_storage() output comes from to determine which version
 * to use by setting d.storage_is_total_moisture() as true (total_volumetric_moisture) or false (this function).
 */

template<SoilMoistureConverterData data>
double soil_moisture_converter<data>::total_volumetric_moisture(data& d) const
{
    /*
     * Volumetric moisture content if d.soil_storage() is all of the moisture in the soil.
     */
    return d.soil_storage()/d.soil_storage_max() * d.porosity();
};

template<SoilMoistureConverterData data>
double soil_moisture_converter<data>::fractional_volumetric_moisture(const double lower_bound_fraction, data& d) const
{
    /*
     * Volumetric moisture content for a case where d.soil_storage() is not the actual total moisture
     * in the soil but rather it is the moisture above a specific threshold, set by lower_bound_fraction.
     *
     * NOTE: this function reduces to total_volumetric_moisture if lower_bound_fraction = 0
     */ 

    return lower_bound_fraction + d.soil_storage()/d.soil_storage_max() 
        * (d.porosity() - lower_bound_fraction);
};

