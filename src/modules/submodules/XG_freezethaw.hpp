#pragma once
#include "base_step.hpp"
#include "Params.hpp"
#include "State.hpp"
#include "Front.hpp"
#include <concepts>
#include <vector>
#include <cmath>

namespace freeze_thaw {

class State;
class Params;

template<typename T>
concept XGData = requires(T t) {
    { t.get_freezethaw_state() } -> std::same_as<State&>;
    { t.get_surface_temp() } -> std::convertible_to<double>;
    { t.get_soil_moist() } -> std::convertible_to<double>;
    { t.get_soil_rechr() } -> std::convertible_to<double>;
    { t.set_thaw_depth(0.0) } -> std::same_as<void>;
    { t.set_freeze_depth(0.0) } -> std::same_as<void>;
};

template<XGData Data>
class Model : public base_step<Model<Data>,Data> {
public:
    // Constructor: takes const reference to Params
    explicit Model(const Params& p) : params(p) {}
    
    virtual ~Model() = default;

    // base_step interface
    virtual void execute() override {}

    // Execute freeze-thaw algorithm for a single grid point
    void execute_impl(Data&);

    // ========== CORE ALGORITHM METHODS ==========
    void freeze(State&, const std::vector<double>& depths);
    void thaw(State&, const std::vector<double>& depths);
    
    double stefan_equation(double& surface_index, double& thermal_conductivity, 
                          size_t& layer, State& state, const std::vector<double>& depths);
    
    double Interpolated_ttc(double Za, size_t layer, State& state, 
                           const std::vector<double>& depths);
    double Interpolated_ftc(double Za, size_t layer, State& state, 
                           const std::vector<double>& depths);

    // ========== INITIALIZATION METHODS ==========
    void find_thaw_D(double dt, State& state, const std::vector<double>& depths);
    void find_freeze_D(double df, State& state, const std::vector<double>& depths);
    void init_freezethaw_degreedays(const double& Zdf_init, const double& Zdt_init, 
                                   State& state, const std::vector<double>& depths);

    // ========== THERMAL PROPERTY METHODS ==========
    double get_ftc(size_t layer, State& state,
                   const double& layer_h2o, const double& theta);
    double get_ttc(size_t layer, State& state,
                   const double& layer_h2o, const double& theta);

    // ========== STATE QUERY METHODS ==========
    bool freezing(const State& state) const;
    bool thawing(const State& state) const;
    bool net_negative_degree_days(State& state, double B) const;
    bool thaw_front_overtaken(const State& state) const;
    bool freeze_front_overtaken(const State& state) const;
    bool exist_excess_fronts(const State& state) const;

    // ========== STATE MANAGEMENT METHODS ==========
    void determine_freeze_thaw_idle(State& state, double B);
    void accumulate_degree_days(State& state, double surface_temp);

    // ========== FRONT MERGE METHODS ==========
    void merge_freeze_fronts(State& state, const std::vector<double>& depths);
    void merge_thaw_fronts(State& state, const std::vector<double>& depths);

private:
    const Params& params;  // Const reference: safe, efficient, non-null
    
    static constexpr double tolerance = 0.000001;
    static constexpr double L = 335000.0;  // latent heat of fusion (J/kg)
};

} // namespace freeze_thaw
