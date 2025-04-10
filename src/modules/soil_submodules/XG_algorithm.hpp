#include "I_freeze_thaw_depths.hpp"
#include <cmath>

class XG_algorithm : public I_freeze_thaw_depths
{
public:
    ~XG_algorithm(double& t, state& _S, param& _P) : surface_temp(t), S(_S), P(_P) {};
    XG_algorithm() {};

    virtual void run() override;

    double& surface_temp;
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
        // Constructor that initializes vector sizes
        explicit state(int N) :
            Zd_front(N),
            pf(N),
            pt(N),
            ttc(N),
            ftc(N),
            tc_composite(N),
            tc_composite2(N),
            theta(N),
            layer_h2o(N),
            XG_max(N),
            XG_moist(N),
            ttc_contents(N),
            ftc_contents(N),
            rechr_fract(N),
            moist_fract(N),
            default_fract(N)
        {
            // Initialize other members
            Zdf = 0.0;
            Zdt = 0.0;
            Th_low = 1;
            Fz_low = 1;
            nfront = 0;
            Bfr = 0.0;
            Bth = 0.0;
            XG_moist_d = 0.0;
            XG_rechr_d = 0.0;
            check_XG_moist = 0.0;
            B = 0.0;
            TrigAcc = 0.0;
            TrigState = 0;
            t_trend = 0.0;
        }

        void set_XG_max(params& P);
        void set_theta(params& P);
        void distriute_moisture(params& P);
        void set_thermal_conductivities(params& P);
        void set_freezethaw_ratios(params& P);

    };

    class params
    {
    public:
    // Core parameters  
        const double Trigthrhld;           // (ºC*d) Trigger reference level  
        const std::vector<double> depths;  // (m) soil layer thicknesses  
        const std::vector<double> por;     // soil porosity  
        const int N_Soil_layers;           // number of soil layers (≤ nlay)  
        const std::vector<double> theta_default;        // (m³/m³) default theta  
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

        Param(
            const double& t, const std::vector<double>& d, const std::vector<double>& p,
            const int& n, const double& td, const double& tm, const double& skm,
            const double& ski, const double& skw, const double& swk, const double& zdi,
            const double& zti, const double& zpi, const int& fku, const int& tku,
            const int& ku, const double& srm, const double& smm,
        ) : 
            Trigthrhld(t), depths(d), por(p), N_Soil_layers(n), theta_default(td),
            theta_min(tm), soil_solid_km(skm), soil_solid_km_ki(ski),
            soil_solid_km_kw(skw), SWE_k(swk), freeze_kw_ki_update(fku), 
            thaw_ki_kw_update(tku), k_update(ku), soil_rechr_max(srm), 
            soil_moist_max(smm), 
        {}
    };
    // Variation #1 parameters  
    //const double& n_factor_a;           // surface-to-air temp ratio  
    //const double& n_factor_b;           // surface-to-air temp ratio  
    //const double& n_factor_c;           // surface-to-air temp ratio  
};

class StateBuilder {
public:
    explicit StateBuilder(int num_layers) {
        state_ = std::make_unique<state>();
        initialize_vectors(num_layers);
    }

    StateBuilder& size_check() 
    {
        //check sizes
    };

    StateBuilder& set_initial_freezethaw_depths(const double Zdf, const double Zdt)
    {
        state_->Zdf = Zdf;
        state_->Zdt = Zdt;

        return *this;
    };

    StateBuilder& set_XG_max(std::vector<double> por,std::vector<double> depths) {

        for (auto&& [m,p,d] : std::views::zip(state_->XG_max,por,depths))
            m = p * d * 1000.0;

        return *this;
    };

    StateBuilder& set_theta(std::vector<double> theta_default) {
        state_->theta = P.theta_default;

        return *this;
    };

    StateBuilder& set_layer_moisture_maximums() {
        double sum = 0.0;
        for (double val : P.depths)
        {
            sum += val;
        }
        if (sum < state_->Zdf || sum < state_->Zdt)
        {
            //TODO Add A CHM exception to say that the total soil depth is less than initial Zdt,Zdf
        } 
        double rechrmax = P.soil_rechr_max;
        double soilmax = P.soil_moist_max;

        for (int layer; layer < P.N_Soil_layers; ++layer)
        {
            state_->XG_max[layer] = P.por[layer] * P.depths[layer] * 1000.0;
            state_->theta[layer] = P.theta_default[layer];

            if (rechrmax > 0.0)
            {
                if (rechrmax > state_->XG_max[layer])
                {
                    state_->XG_rechr_d += P.depths[layer];
                    state_->rechr_fract[layer] = 1.0;
                    rechrmax -= state_->XG_max[layer];
                }
                else
                {
                    const double amount = rechrmax / state_->XG_max[layer];
                    state_->rechr_fract[layer] = rechrmax / state_->XG_max[layer];

                    state_->XG_rechr_d += P.depths[layer] * amount;
                    const double amount_remaining = 1.0 - amount;
                    if (soilmax >= state_->XG_max[layer]*amount_remaining)
                    {
                        state_->moist_fract[layer] = amount_remaining;
                        soilmax -= state_->XG_max[layer] * amount_remaining;
                        state_->XG_moist_d[layer] += P.depths[layer];
                    }
                    else
                    {
                        state_->moist_fract[layer] = (soilmax -  rechrmax) / state_->XG_max[layer];
                        const double used = state_->rechr_fract[layer] + state_->moist_fract[layer];
                        state_->default_fract[layer] = 1.0 - used;
                        state_->XG_moist_d += state_->XG_rechr_d[layer] + P.depths[layer] * used;
                        soilmax = 0.0;
                    }
                    rechrmax = 0.0;
                }
            }
            else if (soilmax > 0.0)
            {
                if (soilmax >= state_->XG_max[layer]) {
                    state_->XG_moist_d += P.depths[layer];
                    state_->moist_fract[layer] = 1.0;
                    soilmax -= state_->XG_max[layer];
                }
                else
                {
                    const double amount = soilmax / state_->XG_max[layer];
                    state_->XG_moist_d[layer] += P.depths[layer] * amount;
                    state_->moist_fract[layer] = amount;
                    state_->default_fract[layer] = 1.0 - amount;
                    soilmax = 0.0;
                }
            }
            else
            {
                state_->default_fract[layer] = 1.0;
            }
        }

        if (rechrmax != 0.0 || soilmax != 0.0)
        {
            // put CHM exception here
        }
        return *this;
    };


    StateBuilder& set_thermal_conductivities(const XG_algorithm::params& P) {
        state_->set_thetmal_conductivities(P);
        return *this;
    };

    StateBuilder& set_freezethaw_ratios(const XG_algorithm::params& P) {
        state_->set_freezethaw_ratios(P); 
        return *this;
    };

    std::unique_ptr<state> build() {
        return std::move(state_);
    };

private:
    void initialize_vectors(int num_layers) {
        state_->theta.resize(num_layers);
        state_->XG_max.resize(num_layers);
        state_->tc_composite.assign(num_layers,0.0);
        state_->tc_composite2.assign(num_layers,0.0);
       // ... resize all other vectors ...
    }

    std::unique_ptr<XG_algorithm::state> state_;
};
