#pragma once
#include "base_step.hpp"
#include <concepts>
#include <stdexcept>
#include <utility>

/*
 * Submodule to compute volumetric moisture content and degree of saturation assuming 
 */ 

template<typename T>
concept SoilMoistureConverterData = requires(T& t) {
    // Inputs
    { t.fractional_cutoff() } -> std::floating_point;
    { t.soil_storage() } -> std::floating_point;
    { t.soil_storage_max() } -> std::floating_point;
    { t.porosity() } -> std::floating_point;
    { t.volumetric_moisture_content() } -> std::floating_point;
    
    // Outputs
    { t.saturation(std::declval<const double>()) } -> std::same_as<void>;
    { t.volumetric_moisture_content(std::declval<const double>()) } -> std::same_as<void>;
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
    double volumetric_moisture(data& d) const;
    class Cutoff; // Forward declaration
};

template<SoilMoistureConverterData data>
void soil_moisture_converter<data>::execute(data& d)
{
    double value;

    value = volumetric_moisture(d);
    
    d.volumetric_moisture_content(value);   

    d.saturation(value / d.porosity()); 
};

/*
 * As described below, the output of d.fraction_cutoff() must be between [0,1] for physical constraints. 
 * I created this simple type 'Cutoff' which performs automatic conversions to double while checking the
 * value upon construction.
 */
template<SoilMoistureConverterData data>
class soil_moisture_converter<data>::Cutoff
{
    double _c;
public:
    Cutoff(const double c) : _c(c)
    {
        if (c < 0.0 || c > 1.0)
            throw std::logic_error("cutoff must be between (inclusive) 0 and 1");
    };
    ~Cutoff() = default;
    operator double() const { return _c; }
};

/*
 * It's common for the soil moisture storage (a depth) to not actually include all moisture in the soil in
 * the variable. Depending on the processes, not all moisture in all pores can be included in that process. 
 * For example, transpiration does not occur below the wilt point and so if the important process is not 
 * about transpiration or there is no transpiration, the soil storage might only track moisture above the 
 * wilt point.
 *
 * Therefore, volumetric moisture accouts for it by using the following equation:
 *
 * $$\theta = \gamma + S/S_{\max} * (\phi - \gamma)$$
 *
 * where $\theta$ is the volumetric moisture content, $\gamma$ is the lower bound percentage of moisture 
 * that is included in $S$ (soil moisture storage), $\phi$ is the porosity (e.g., volumetric moisture content 
 * at saturation).
 *
 * This equation basically says that $S$ and $S_{\max}$ only capture soil soil moisture content above the 
 * percentage set by $\gamma$.
 */
template<SoilMoistureConverterData data>
double soil_moisture_converter<data>::volumetric_moisture(data& d) const
{

    /*
     * Volumetric moisture content for a case where d.soil_storage() is not the actual total moisture
     * in the soil but rather it is the moisture above a specific threshold, set by fractional_cutoff().
     */ 
    
    Cutoff lower_bound_fraction = d.fractional_cutoff();

    return lower_bound_fraction + d.soil_storage()/d.soil_storage_max() 
        * (d.porosity() - lower_bound_fraction);
};
