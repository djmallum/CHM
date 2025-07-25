#pragma once
#include "base_step.hpp"

template<class data>
class volumetric : public base_step
{
public:
    explicit volumetric(Soil::soils_na& S) : SoilDataObj(S) {};
    ~volumetric() {};

    void execute(data& d) override final;

    // Currently constants (e.g., wilt_point) and inputs (e.g., CurrentDay) are indistinguishable
    // TODO Consider in the future if having them be different is necessary.
private:
    std::optional<size_t> CalcDay_;

    void set_CalcDay(const size_t& D)
    {
        if (CalcDay_.has_value())
            throw std::runtime_error("CalcDay already set");

        CalcDay_ = D;
    };

    void get_CalcDay()
    {
        if (!CalcDay_.has_value())
            throw std::runtime_error("CalcDay not yet set");

        return *CalcDay_;
    };
};

template<class data>
void volumetric<data>::execute(data& d)
{
    if (get_CalcDay() != d.CurrentDay())
        return;
   
    // TODO Verify equation. CRHM has two formulas. One which uses recharge moisture and the other uses the full column.
    // In CRHM runs I was given, the full soil column was used, not the recharge layer. That said, recharge might be more accurate since infiltration is injected into the top layer, nowhere else. 
    double volumetric_moisture = recharge_volumetric_moisture(); 

    d.volumetric_moisture_content(volumetric_moisture);
        
        (d.soil_recharge_storage() / d.soil_depth() + d.wilt_point()); // TODO Verify units compared to CRHM use
    double saturation_fraction = get_saturation_fraction();
    if (d.porosity() > 0.0)
    {
        soil_storage = volumetric_moisture / d.porosity(); // TODO Verify units compared to CRHM use      
    };

    d.degree_of_saturation(saturation_fraction); // TODO verify units as used in CHM
};

template<class data>
double& volumetric<data>::recharge_volumetric_moisture(data& d)
{
    /*
     * Soil moisture s_r in the recharge layer is typically transported as the amount of moisture in mm
     * above the wilt point. Volumetric moisture content does not make a distinction between above
     * or below the wilt point or field capacity. Therefore, we must first convert the soil moisture 
     * content from a depth above the wilt point in mm to a percentage, and then convert this percentage
     * to the actual moisture that it represents in a volume fraction.
     *
     * TODO Finish note
     */ 
    double available_moisture_max = d.porosity() - d.wilt_point();
    double available_moisture = d.soil_recharge_storage() / d.soil_recharge_depth();

    return actual_moisture * available_moisture_max + d.wilt_point(); 


    return d.soil_recharge_storage() / d.soil_recharge_depth();   
};

template<class data>
double& volumetric<data>::get_saturation_fraction(data& d)
{
    
};



