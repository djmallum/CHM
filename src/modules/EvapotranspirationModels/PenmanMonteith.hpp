#pragma once

#include "evapbase.hpp"
// Implementation:
// include the header
// std::unique_ptr<evapT_base> mymodel = std::make_unqiue<PenmanMonteith>();
// std::unique_ptr<PM_param> param = std::make_unique<PM_param>(); TODO may not need to make param a pointer because it doesnt need to be released on each timestep (in fact, it shouldn't)
// std::unique_ptr<PM_var> var = std::make_unique<PM_var>();
// 
// var must be remade because values are gathered from other modules.
// set param and var members. param never changes, var changes every timestep and should be released.

struct PM_vars;
struct PM_output;

class PenmanMonteith : public evapT_base
{
public:

    PenmanMonteith(const double& LAI, const double& LAImax, const double& veg_Ht, const double& wind_height, 
            const double& stomatal_res_min, const double& soil_d, const double& F_to_g, const double& s_per_step,
			const double Cp, const double K, const double tension, const double pore_sz, 
            const double theta_pwp, const double phi); 
   
                     
    ~PenmanMonteith(void) override; // Deconstructor
                     

    void CalcEvapT(var_base& vars, model_output& output) override;
    
private:

    // dont delete
    void CalcHeights(void);
    double CalcAeroResistance(const PM_vars& var);
    double CalcStomatalResistance(const PM_vars& var);
    double AirDensity(const double t, const double ea, const double Pa);

    double Z0;
    double d;
    bool has_vegetation;
    bool IsFirstRun = true;
    static constexpr double water_density = 1000; //kg/m^3

    const double& leaf_area_index;
    const double& leaf_area_index_max;
    const double& Veg_height;
    const double& wind_measurement_height; // This one might be uniform...
    const double& stomatal_resistance_min; // Also might be domain wide
    const double& soil_depth;
    const double& Frac_to_ground;
    const double& s_per_time_step;  
    const double heat_capacity_air;
    const double kappa; // also might be domain wide
    const double air_entry_tension;
    const double pore_size_dist; 
    const double wilt_point;
    const double porosity;
};


struct PM_vars : public var_base
{
    const double& wind_speed;
    const double& short_wave_in;
    const double& all_wave_net;
    const double& t;
    const double& soil_storage;
    const double& vapour_pressure;
    const double& saturated_vapour_pressure; 
    const double& P_atm;

    PM_vars(const double& U, const double& Qsw, const double& Qnet, const double& temp, const double& soil, const double& ea, const double& ea_star, const double& P) 
        : wind_speed(U), short_wave_in(Qsw), all_wave_net(Qnet), t(temp), soil_storage(soil), vapour_pressure(ea), saturated_vapour_pressure(ea_star), P_atm(P) {};
};


struct PM_output : public model_output
{
    double stomatal_resistance;
};
