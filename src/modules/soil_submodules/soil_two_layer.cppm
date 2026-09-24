module;
#include <algorithm>
#include <cmath>
#include <cstddef>

export module soil_two_layer;

namespace soil_two_layer
{
export struct Params;
export struct Input;
export struct Output;
export struct ThawFractions
{
    double rechr = 0.0;
    double lower = 0.0;
    void recompute(const Input& input,const Params& params);
};
export struct State {
    double soil_storage = 0.0;
    double soil_rechr_storage = 0.0;
    double depression_storage = 0.0;
    double soil_excess_to_gw = 0.0;
    ThawFractions thaw_fraction{};
    double detention_storage = 0.0;
    double detention_max = 0.0; // This is the only storage parameter that can change
    double ground_water_storage = 0.0;

    Output step_soil(
        const Input& input, 
        const Params& params
    );
private:    
    void distribute_infiltration(input, K);
    void manage_detention(input, params, K);
    void manage_depression(input, params);
    void manage_groundwater(input, params, K);
    void manage_subsurface_runoff(input);
    void remove_tiny_moisutre();
};

export struct Output {
    double condensation = 0.0;
    double excess = 0.0;
    double depression_to_gw = 0.0;
    double soil_excess_to_runoff = 0.0;
    double runoff_to_depression = 0.0;
    double ground_water_out = 0.0;
    double soil_to_ssr = 0.0;
    double rechr_to_ssr = 0.0;
    double remaining_ET = 0.0;
}

export struct Input
{
    double swe = 0.0;
    bool new_day = false;
    double infil = 0.0;
    double runoff = 0.0;
    double actual_ET = 0.0;
    double thaw_front_depth = 0.0;
    double freeze_front_depth = 0.0;
    double freeze_thaw_first_front = 0.0;
    double routing_residual = 0.0;
};

export struct SaturatedConductivity
{
    double Ksaturated_rechr = 0.0;
    double Ksaturated_lower = 0.0;
    double Ksaturated_ground_water = 0.0;
    double Ksaturated_snow = 0.0;
    double Ksaturated_organic = 0.0;
};

export struct Params
{
    double soil_storage_max = 0.0;
    double soil_rechr_max = 0.0;
    double depression_max = 0.0;
    double detention_snow_max = 0.0;
    double detention_organic_max = 0.0;
    double ground_water_max = 0.0;
    double local_slope = 0.0;
    double pore_size_dist = 0.0;
    double pore_size_dist_organic = 0.0;
    double soil_index = 0.0;
    double snow_grain_diameter = 0.0;
    double detention_snow_init = 0.0; 
    double detention_organic_init = 0.0; 
    double snow_covered_threshold = 0.0;
    double porosity = 0.0;
    SaturatedConductivity Ksaturated{};
    bool allow_runoff_from_infiltration = false;
    int seconds_per_step = 0;
    bool excess_to_ssr = true;
};

namespace
{
void push_excess_down(double& storage, double maximum, double& downstream)
{
    downstream += storage - maximum;
    storage = maximum;
};

void ThawFractions::recompute(const Input& input, const Params& params)
{
    rechr = 0.0;
    lower = 0.0;

    if (params.soil_storage_max == 0.0 || params.porosity == 0.0)
        return;

    const double recharge_depth = params.soil_rechr_max / params.porosity / 1000.0;
    const double soil_depth = params.soil_storage_max / params.porosity / 1000.0;
    if (!params.allow_runoff_from_infiltration || input.freeze_thaw_first_front == 0.0)
        rechr = lower = 1.0;

    if (params.allow_runoff_from_infiltration && input.freeze_thaw_first_front > 0.0)
    {
        if (input.thaw_front_depth < recharge_depth)
            rechr = input.thaw_front_depth / recharge_depth;
        else if (input.thaw_front_depth < soil_depth)
        {
            rechr = 1.0;
            lower = (input.thaw_front_depth - recharge_depth) /
                                         (soil_depth - recharge_depth);
        }
        else
            rechr = lower = 1.0;
    }
};

void State::distribute_infiltration(
    Output& output,
    const Input& input, 
    const UnsaturatedConductivity& K
)
{
    if (params.soil_storage_max <= 0.0)
    {
        excess = input.infil + condensation;
        return;
    }

    double lower_storage = soil_storage - soil_rechr_storage;
    const double potential = input.infil + condensation;
    double possible = thaw_fraction.rechr *
                      (params.soil_rechr_max - params.soil_rechr_storage);
    if (possible > potential || !params.allow_runoff_from_infiltration)
        possible = potential;
    else
        soil_excess_to_runoff = potential - possible;

    soil_rechr_storage += possible;
    if (soil_rechr_storage > params.soil_rechr_max)
        push_excess_down(soil_rechr_storage, params.soil_rechr_max, lower_storage);
    soil_storage = lower_storage + soil_rechr_storage;
    if (soil_storage > params.soil_storage_max)
        push_excess_down(soil_storage, params.soil_storage_max,
                         soil_excess_to_gw);

    if (input.swe == 0.0 && params.soil_rechr_max > 0.0)
    {
        rechr_to_ssr = soil_rechr_storage / params.soil_rechr_max *
                             K.rechr_to_ssr * thaw_fration.rechr;
        output.rechr_to_ssr = std::min(output.rechr_to_ssr,
                                      soil_rechr_storage * thaw_fraction.rechr);
        soil_rechr_storage = std::max(0.0,
                                            soil_rechr_storage - output.rechr_to_ssr);
        state.soil_storage -= output.rechr_to_ssr;
        output.soil_to_ssr = output.rechr_to_ssr;
    }

    const double groundwater_limit = K.soil_to_gw * thaw_fraction.lower;
    if (soil_excess_to_gw > groundwater_limit)
        push_excess_down(soil_excess_to_gw, groundwater_limit, excess);
    if (params.excess_to_ssr && excess > 0.0)
    {
        const double excess_to_ssr_value = excess * (1.0 - thaw_fraction.lower);
        push_excess_down(excess, excess_to_ssr_value, output.soil_to_ssr);
    }
}

enum class DetentionRegime {
    Snow,
    Organic
};

constexpr DetentionRegime detention_regime_for(double swe)
{
    return swe > 0.0 ? DetentionRegime::Snow : DetentionRegime::Organic;
}

void manage_detention(
    const Output& output,
    const Input& input, 
    const Params& params,
    const UnsaturatedConductivity& K
)
{
    const auto regime = detention_regime_for(input.swe);
    state.detention_max = regime == DetentionRegime::Organic 
        ? params.detention_organic_max : params.detention_snow_max;

    output.soil_excess_to_runoff += input.runoff + output.excess + input.routing_residual; // TODO Deal with routing Residual

    if (state.soil_excess_to_runoff > 0.0)
    {
        const double space = state.detention_max - state.detention_storage;
        if (space > 0.0)
        {
            const double stored = std::min(output.soil_excess_to_runoff, space);
            state.detention_storage += stored;
            state.soil_excess_to_runoff -= stored;
        }
    }
    if (state.detention_storage > 0.0 && K.detention_to_runoff > 0.0)
    {
        const double transfer = std::min(state.detention_storage,
                                          K.detention_to_runoff);
        state.detention_storage -= transfer;
        state.soil_excess_to_runoff += transfer;
        if (state.detention_storage < 0.0001)
            state.detention_storage = 0.0;
    }
}

void manage_depression(State& state, const Input& input, const Params& params)
{
    if (state.soil_excess_to_runoff > 0.0 && input.depression_max > 0.0)
    {
        double space = (input.depression_max - state.depression_storage) *
                       (1.0 - std::exp(-std::min(12.0,
                                                  state.soil_excess_to_runoff /
                                                  input.depression_max)));
        if (input.soil_storage_max == 0.0)
            space = input.depression_max - state.depression_storage;
        if (space > 0.0)
        {
            const double stored = std::min(state.soil_excess_to_runoff, space);
            state.depression_storage += stored;
            state.soil_excess_to_runoff -= stored;
            state.runoff_to_depression += stored;
        }
    }
    if (state.depression_storage > 0.0 && K.depression_to_gw > 0.0)
    {
        const double transfer = std::min(state.depression_storage,
                                          K.depression_to_gw);
        state.depression_storage -= transfer;
        state.depression_to_gw += transfer;
    }
}

void manage_groundwater(
    State& state, 
    const Input& input, 
    const Params& params,
    const UnsaturatedConductivity& K
)
{
    state.soil_excess_to_gw += state.depression_to_gw;
    state.depression_to_gw = 0.0;
    state.ground_water_storage += state.soil_excess_to_gw;
    if (state.ground_water_storage > params.ground_water_max)
        push_excess_down(state.ground_water_storage, params.ground_water_max,
                         state.ground_water_out);
    if (input.ground_water_max > 0.0)
    {
        const double spilled = state.ground_water_storage / params.ground_water_max *
                               K.ground_water_out;
        state.ground_water_storage -= spilled;
        state.ground_water_out += spilled;
    }
}

void manage_subsurface_runoff(State& state, const Input& input)
{
    if (state.depression_storage > 0.0 && K.depression_to_ssr > 0.0)
    {
        const double transfer = std::min(state.depression_storage,
                                          K.depression_to_ssr);
        state.depression_storage -= transfer;
        state.soil_to_ssr += transfer;
    }
    if (K.lower_to_ssr > 0.0)
    {
        const double available = state.soil_storage - state.soil_rechr_storage;
        const double transfer = std::min(K.lower_to_ssr * state.thaw_fraction_lower,
                                          available);
        state.soil_storage -= transfer;
        state.soil_to_ssr += transfer;
    }
}

struct UnsaturatedConductivity
{
    double soil_to_gw = 0.0;
    double rechr_to_ssr = 0.0;
    double lower_to_ssr = 0.0;
    double detention_to_runoff = 0.0;
    double depression_to_ssr = 0.0;
    double depression_to_gw = 0.0;
    double ground_water_out = 0.0;
    UnsaturatedConductivity(const State& state);
}
static void remove_tiny_moisutre(State& state)
{
    if (state.soil_storage < 0.0001)
        state.soil_storage = 0.0;
    if (state.soil_rechr_storage < 0.0001)
        state.soil_rechr_storage = 0.0;
    if (state.detention_storage < 0.0001)
        state.detention_storage = 0.0;
    if (state.depression_storage < 0.0001)
        state.depression_storage = 0.0;
    if (state.ground_water_storage < 0.0001)
        state.ground_water_storage = 0.0;
}
UnsaturatedConductivity::UnsaturatedConductivity(const State& state,
                                                 const Params& params))
{
    if (state.soil_storage_max <= 0.0) return;

    const double exponent = 3.0 + 2.0 / params.pore_size_dist;
    const double organic_exponent = 3.0 + 2.0 / params.pore_size_dist_organic;
    const double factor = 1000.0 * 9.8 * 0.001787;
    double lateral_lower = params.Ksaturated_lower *
        std::pow((state.soil_storage - state.soil_rechr_storage) /
                 (params.soil_storage_max - params.soil_rechr_max), exponent) *
        std::tan(input.local_slope);
    double vertical = params.Ksaturated_lower *
        std::pow(state.soil_storage / params.soil_storage_max, exponent);
    double lateral_groundwater = params.Ksaturated_ground_water * std::tan(state.local_slope);
    double lateral_detention = 0.0;
    if (state.detention_max > 0.0)
    {
        if (state.swe > 0.0 && state.snow_density > 100.0)
        {
            const double saturated_snow = 0.077 * std::pow(state.snow_grain_diameter / 1000.0, 2.0) *
                                          std::exp(-7.8 * (state.snow_density / 1000.0)) * factor;
            lateral_detention = saturated_snow *
                std::pow(state.detention_storage / state.detention_max, state.soil_index) *
                std::sin(state.local_slope);
        }
        else
            lateral_detention = params.Ksaturated_organic *
                std::pow(state.detention_storage / state.detention_max, organic_exponent) *
                std::tan(state.local_slope);
    }
    const double unit_changer = state.seconds_per_step * 1000.0;
    const double lateral_units = unit_changer / 1000.0;
    
    rechr_to_ssr = input.swe > 0.0 ? 0.0 :
            params.Ksaturated_rechr * std::pow(state.soil_rechr_storage /
            state.soil_rechr_max, exponent) * std::tan(state.local_slope) *
            state.soil_rechr_max * lateral_units;
    lower_to_ssr = lateral_lower * (state.soil_storage_max - state.soil_rechr_max) * lateral_units;
    depression_to_ssr = lateral_lower * state.soil_storage_max * lateral_units;
    depression_to_gw = input.swe > 0.0 ? 0.0 : vertical * unit_changer;
    soil_to_gw = input.swe > 0.0 ? vertical * unit_changer : vertical * unit_changer;
    ground_water_out = lateral_groundwater * state.ground_water_storage * lateral_units;
    detention_to_runoff = lateral_detention * state.detention_max * lateral_units;
};

export Output step_soil(
    State& state, 
    const Input& input, 
    const Params& params
)
{
    Output output{};
    if (params.soil_storage_max == 0.0)
    {
        output.excess = input.infil + state.condensation;
        return output;
    }

    thaw_fraction.recompute(input, params);

    const auto K = UnsaturatedConductivity(*this, params);

    if (input.actual_ET < 0.0 && input.swe == 0.0)
    {
        output.condensation = -state.actual_ET;
        output.remaining_ET = 0.0;
    }

    distribute_infiltration(input, K);
    manage_detention(input, params, K);
    manage_depression(input, params);
    manage_groundwater(input, params, K);
    manage_subsurface_runoff(input);
    remove_tiny_moisutre();

    return output;
}
}