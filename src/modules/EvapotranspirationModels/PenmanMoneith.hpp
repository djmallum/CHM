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


class PenmanMoneith : public evapT_base
{
public:

    PenmanMoneith(PM_param& param); // TODO Add inputs here to constructor, just param)
                     
    ~PenmanMoneith(void) override; // Deconstructor
                     

    double CalcEvapT(param_base& param, var_base& var) const override;

private:

    // dont delete
    void CalcHeights(param_base& param, var_base& var);
    void CalcAeroResistance(param_base& param, var_base& var);
    void CalcSaturationVapourPressure(param_base& param, var_base& var);
    void CalcStomatalResistance(param_base& param, var_base& var);

}


struct PM_param : public param_base 
{
    double Veg_height;
    double Veg_height_max;
    double wind_measurement_height; // This one might be uniform...
    double kappa; // also might be domain wide
    double stomatal_resistance_min; // Also might be domain wide
    double soil_depth;
    double Frac_to_ground;
    double Z0;
    double d;
    double Cp;
} // TODO need a separate container for parameters which never change and are global.


struct PM_var : public var_base
{
    double wind_speed;
    double ShortWave_in;
    double Rnet;
    double t;
    double soil_storage;
    double aero_resistance;
    double stomatal_resistance;
}


    // Make this a variable passed to evap? double ShortWave_in;
    // Same as above double t;


    
