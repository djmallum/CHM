#include "evapbase.hpp"
#include <algorithm>

// Implementation:
// include the header
// std::unique_ptr<evapT_base> mymodel = std::make_unqiue<PenmanMonteith>();
// std::unique_ptr<PM_param> param = std::make_unique<PM_param>(); TODO may not need to make param a pointer because it doesnt need to be released on each timestep (in fact, it shouldn't)
// std::unique_ptr<PM_var> var = std::make_unique<PM_var>();
// 
// var must be remade because values are gathered from other modules.
// set param and var members. param never changes, var changes every timestep and should be released.


class PenmanMonteith : public evapT_base
{
public:

    PenmanMonteith(const double& Cp, const double& K, const double& tension, const double& pore_sz, const double& theta_pwp, const double& phi); 
                     
    ~PenmanMonteith(void) override; // Deconstructor
                     

    void CalcEvapT(var_base& vars, model_output& output) override;
    
    // TODO These should be references
    double Veg_height;
    double Veg_height_max;
    double wind_measurement_height; // This one might be uniform...
    double stomatal_resistance_min; // Also might be domain wide
    double soil_depth;
    double Frac_to_ground;
    double LAImin;
    double seasonal_growth;
    double LAImax;
    const double& heat_capcity_air;
    const double& kappa; // also might be domain wide
    const double& air_entry_tension;
    const double& pore_size; 
    const double& wilt_point;
    const double& porosity; 
private:

    // dont delete
    void CalcHeights(void);
    double CalcAeroResistance(const PM_vars& var);
    double CalcStomatalResistance(const PM_vars& var);
    double AirDensity(double& t, double& ea, double& Pa);

    double Z0;
    double d;
    bool IsFirstRun = true;
}


struct PM_vars : public var_base
{
    double wind_speed;
    double short_wave_in;
    double all_wave_net;
    double t;
    double soil_storage;
    double vapour_pressure;
    double saturated_vapour_pressure; 

    // TODO Add a constructor, make these references.
    // maybe not needed, remember it is normal to make a local copy of values on faces, which these are.
}


    // Make this a variable passed to evap? double ShortWave_in;
    // Same as above double t;


    
