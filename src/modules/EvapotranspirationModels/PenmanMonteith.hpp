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

    PenmanMonteith(); // TODO Add inputs here to constructor, just param)
                     
    ~PenmanMonteith(void) override; // Deconstructor
                     

    void CalcEvapT(PM_vars& vars, model_output& output) override;
    
    // TODO some of these are global, might need to pass them as arguments and store more efficiently. 
    double Veg_height;
    double Veg_height_max;
    double wind_measurement_height; // This one might be uniform...
    double kappa; // also might be domain wide
    double stomatal_resistance_min; // Also might be domain wide
    double soil_depth;
    double Frac_to_ground;
    double Cp;
    double LAImin;
    double seasonal_growth;
    double LAImax;
    double air_entry_tension;
    double pore_size; 
private:

    // dont delete
    void CalcHeights(void);
    double CalcAeroResistance(var_base& var);
    double CalcStomatalResistance(var_base& var);
    double AirDensity(double& t, double& ea, double& Pa);

    double Z0;
    double d;
    bool IsFirstRun = true;
}


struct PM_var : public var_base
{
    double wind_speed;
    double ShortWave_in;
    double Rnet;
    double t;
    double soil_storage;
    double vapour_pressure;
    double saturated_vapour_pressure; 
}


    // Make this a variable passed to evap? double ShortWave_in;
    // Same as above double t;


    
