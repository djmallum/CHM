#include "I_freeze_thaw_depths.hpp"
#include <cmath>

class XG_algorithm : public I_freeze_thaw_depths
{
public:
    ~XG_algorithm(double& t, double& _SWE, state& _S, param& _P) : surface_temp(t), SWE(_SWE), S(_S), P(_P) {};
    XG_algorithm() {};

    virtual void run() override;

    double& surface_temp;
    double& SWE;
    class state;
    class param;
    state& S;
    param& P;
    
    class state
    {
    public:
        // Depth-related
        double Zdf;                // (m) depth of freezing front  
        double Zdt;                // (m) depth of thawing front  
        std::vector<double> Zd_front; // (m) depths of all freezing/thawing fronts (thaw: +, freeze: -)  
        int Th_low;            // lowest thawed layer  
        int Fz_low;            // lowest frozen layer  
        int nfront;                // number of freezing/thawing fronts  

        // Degree-day counters  
        double Bfr;                // (ºC*d) freeze degree days  
        double Bth;                // (ºC*d) thaw degree days  

        // Layer ratios  
        std::vector<double> pf;    // soil layers freezing ratios  
        std::vector<double> pt;    // soil layers thawing ratios  

        // Thermal properties  
        std::vector<double> ttc;   // (W/(m*K)) thawing thermal conductivity  
        std::vector<double> ftc;   // (W/(m*K)) freezing thermal conductivity  
        std::vector<double> tc_composite;  // (W/(m*K)) composite ftc/ttc  
        std::vector<double> tc_composite2; // (W/(m*K)) composite2 ftc/ttc  

        // Moisture/water content  
        std::vector<double> theta;           // (m³/m³) layer theta  
        std::vector<double> layer_h2o;       // (kg/m³) layer water content  
        double XG_moist_d;      // (m) layer depth of soil moisture in XG  
        double XG_rechr_d;      // (m) layer depth of soil recharge in XG  
        std::vector<double> XG_max;          // (mm) layer max soil moisture  
        std::vector<double> XG_moist;        // (mm) layer moisture content  
        double check_XG_moist;     // (mm) sum of XG soil_moist  

        // Local/temporary  
        double B;                  // (ºC*d) interval degree-day sum  
        double TrigAcc;            // (ºC*d) freeze/thaw cycle detection  
        int TrigState;             // 1/0/-1 → thaw/idle/freeze  
        std::vector<int> ttc_contents;          // 0/1 → thaw/freeze  
        std::vector<int> ftc_contents;          // 0/1 → freeze/thaw  
        double t_trend;            // (°C) temperature long-term trend  

        // Fractions  
        std::vector<double> rechr_fract;   // fraction of layer (soil_rechr_max)  
        std::vector<double> moist_fract;   // fraction of layer (soil_moist_max)  
        std::vector<double> default_fract; // fraction of layer (theta_default)  

    // Variation #1 only  
    //int n_factor_T;            // days after start of thaw  
    //double n_factor;           // calculated n_factor value
    };

    class params
    {
    public:
    // Core parameters  
        const double Trigthrhld;           // (ºC*d) Trigger reference level  
        const std::vector<double> depths;  // (m) soil layer thicknesses  
        const std::vector<double> por;     // soil porosity  
        const int N_Soil_layers;           // number of soil layers (≤ nlay)  
        const double theta_default;        // (m³/m³) default theta  
        const double theta_min;            // (m³/m³) minimum theta  
        const std::vector<double> soil_solid_km;        // (W/(m*K)) dry soil conductivity  
        const std::vector<double> soil_solid_km_ki;     // (W/(m*K)) saturated frozen conductivity  
        const std::vector<double> soil_solid_km_kw;     // (W/(m*K)) saturated unfrozen conductivity  
        const double SWE_k;                // (W/(m*K)) snow thermal conductivity UNUSED IN CRHM 
        //Next 3 are initial conditions
        //const double Zdf_init;             // (m) initial freezing front depth  
        //const double Zdt_init;             // (m) initial thawing front depth  
        //const double Zpf_init;             // (m) initial permafrost depth  
        const bool freeze_kw_ki_update;     // update kw→ki behind freeze front  
        const bool thaw_ki_kw_update;       // update ki→kw behind thaw front  
        const int k_update;                // 0=never, 1=post-layer, 2=continuous  
        const double soil_rechr_max;       // (mm) max recharge zone capacity  
        const double soil_moist_max;       // (mm) max rooting zone capacity  
        const bool is_newday;
        const double time_step_per_day;
    };
    // Variation #1 parameters  
    //const double& n_factor_a;           // surface-to-air temp ratio  
    //const double& n_factor_b;           // surface-to-air temp ratio  
    //const double& n_factor_c;           // surface-to-air temp ratio  
};
