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


class PriestleyTaylor : public evapT_base
{
public:

    PriestleyTaylor(double& alpha_const); 
                     
    ~PriestleyTaylor(void) override; // Deconstructor
                     

    void CalcEvapT(var_base& vars, model_output& output) override;

    double Frac_to_ground;
    const double& alpha; 
    
private:


};


struct PT_vars : public var_base
{
    double& all_wave_net;  
    double& P_atm;

    PT_vars(double& Q, double& P) : all_wave_net(Q), P_atm(P) {};
};

    // Make this a variable passed to evap? double ShortWave_in;
    // Same as above double t;


    
