// ClarkUnitHydrograph.cpp
#include "ClarkUnitHydrograph.hpp"
#include <algorithm>
#include <vector>

ClarkUnitHydrograph::ClarkUnitHydrograph(StorageView view) 
    : view_(std::move(view)) 
{
    initialize_lag_array();
}

void ClarkUnitHydrograph::execute() {
    auto& state = view_.get_state();
    auto& lag_array = view_.get_lag_array();
    
    // Add current input and any pending release to lag array
    lag_array[state.current_lag_index] = view_.get_input() + state.pending_release;
    state.pending_release = 0.0;
    
    // Advance circular buffer
    state.current_lag_index = (state.current_lag_index + 1) % state.lag_steps;
    
    // Calculate Clark unit hydrograph output
    const double delayed_input = lag_array[state.current_lag_index];
    const auto [c1, c2] = calculate_coefficients();
    
    const double output = c1 * (delayed_input + state.last_input) + c2 * state.last_output;
    
    // Update outputs and state
    view_.get_output() = output;
    state.last_input = delayed_input;
    state.last_output = output;
}

double ClarkUnitHydrograph::change_storage_coefficient(double new_storage_coefficient) {
    auto& state = view_.get_state();
    const auto [current_c1, current_c2] = calculate_coefficients();
    
    if (current_c2 >= 1.0) return 0.0; // No delay case
    
    // Calculate current storage
    const double current_storage = (1.0 / (1.0 - current_c2)) * 
                                 (current_c1 * state.last_input + current_c2 * state.last_output);
    
    if (current_storage <= 0.0) return 0.0;
    
    // Recalculate with new coefficient
    const double new_c2 = calculate_c2(new_storage_coefficient);
    const double new_last_output = (current_storage * (1.0 - new_c2) - 
                                  calculate_c1(new_storage_coefficient) * state.last_input) / new_c2;
    
    state.last_output = new_last_output;
    return current_storage;
}

double ClarkUnitHydrograph::change_lag_time(double new_lag_hours) {
    auto& state = view_.get_state();
    auto& lag_array = view_.get_lag_array();
    
    const long new_lag_steps = calculate_lag_steps(new_lag_hours);
    double total_lag_storage = 0.0;
    
    // Calculate current water in lag system
    for (long i = 1; i < state.lag_steps; ++i) {
        total_lag_storage += lag_array[(state.current_lag_index + i) % state.lag_steps];
    }
    
    if (new_lag_steps == state.lag_steps) {
        return total_lag_storage;
    }
    
    // Redistribute storage when lag time changes
    redistribute_lag_storage(new_lag_steps, total_lag_storage);
    
    return total_lag_storage;
}

double ClarkUnitHydrograph::total_water_in_system() const {
    const auto& state = view_.get_state();
    const auto& lag_array = view_.get_lag_array();
    
    double lag_storage = 0.0;
    for (long i = 1; i < state.lag_steps; ++i) {
        lag_storage += lag_array[(state.current_lag_index + i) % state.lag_steps];
    }
    
    const auto [c1, c2] = calculate_coefficients();
    if (c2 >= 1.0) return 0.0; // No storage case
    
    const double linear_storage = (1.0 / (1.0 - c2)) * (c1 * state.last_input + c2 * state.last_output);
    
    return lag_storage + linear_storage;
}

void ClarkUnitHydrograph::initialize_lag_array() {
    auto& state = view_.get_state();
    state.lag_steps = calculate_lag_steps(view_.get_lag_time());
    state.current_lag_index = 0;
    state.pending_release = 0.0;
    
    auto& lag_array = view_.get_lag_array();
    lag_array.resize(state.lag_steps);
    std::fill(lag_array.begin(), lag_array.end(), 0.0);
}

long ClarkUnitHydrograph::calculate_lag_steps(double lag_hours) const {
    // Using a simplified calculation - in practice this would use the model frequency
    const double default_frequency = 24.0; // Assuming hourly time steps
    return static_cast<long>(std::max(lag_hours, 0.0) / 24.0 * default_frequency + 1.1);
}

std::pair<double, double> ClarkUnitHydrograph::calculate_coefficients() const {
    const double storage_coeff = view_.get_storage_coefficient();
    const double default_interval = 1.0 / 24.0; // Assuming daily time step in days
    const double denominator = storage_coeff + default_interval * 0.5;
    
    const double c1 = default_interval * 0.5 / denominator;
    const double c2 = (storage_coeff - default_interval * 0.5) / denominator;
    
    return {c1, c2};
}

double ClarkUnitHydrograph::calculate_c1(double storage_coefficient) const {
    const double default_interval = 1.0 / 24.0; // Assuming daily time step in days
    return default_interval * 0.5 / (storage_coefficient + default_interval * 0.5);
}

double ClarkUnitHydrograph::calculate_c2(double storage_coefficient) const {
    const double default_interval = 1.0 / 24.0; // Assuming daily time step in days
    return (storage_coefficient - default_interval * 0.5) / (storage_coefficient + default_interval * 0.5);
}

void ClarkUnitHydrograph::redistribute_lag_storage(long new_lag_steps, double total_storage) {
    auto& state = view_.get_state();
    auto& lag_array = view_.get_lag_array();
    
    if (new_lag_steps == 1) {
        // Release all storage immediately
        state.pending_release = total_storage;
        lag_array.resize(1);
        lag_array[0] = 0.0;
    } else if (state.lag_steps > 1 && total_storage > 0.0) {
        // Accumulate current distribution
        std::vector<double> cumulative_distribution(state.lag_steps);
        cumulative_distribution[0] = 0.0;
        for (long i = 1; i < state.lag_steps; ++i) {
            cumulative_distribution[i] = cumulative_distribution[i-1] + 
                                       lag_array[(state.current_lag_index + i) % state.lag_steps];
        }
        
        // Create new distribution
        lag_array.resize(new_lag_steps);
        std::fill(lag_array.begin(), lag_array.end(), 0.0);
        
        double last_cumulative = 0.0;
        for (long new_step = 1; new_step < new_lag_steps - 1; ++new_step) {
            const double scaled_position = double(new_step) / (new_lag_steps - 1) * (state.lag_steps - 1);
            const long original_step = static_cast<long>(scaled_position);
            const double fraction = scaled_position - original_step;
            
            const double new_cumulative = cumulative_distribution[original_step] + 
                                        fraction * (cumulative_distribution[original_step + 1] - 
                                                   cumulative_distribution[original_step]);
            
            lag_array[new_step] = new_cumulative - last_cumulative;
            last_cumulative = new_cumulative;
        }
        
        lag_array[new_lag_steps - 1] = cumulative_distribution[state.lag_steps - 1] - last_cumulative;
    }
    
    state.lag_steps = new_lag_steps;
    state.current_lag_index = 0;
}
