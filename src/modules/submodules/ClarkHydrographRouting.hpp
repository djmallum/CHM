// ClarkUnitHydrograph.hpp
#pragma once
#include <functional>
#include <vector>
#include <memory>

class ClarkUnitHydrograph {
public:
    struct State {
        double last_input = 0.0;
        double last_output = 0.0;
        double pending_release = 0.0;
        long current_lag_index = 0;
        long lag_steps = 1;
    };

    struct StorageView {
        std::function<double()> get_input;
        std::function<double()> get_storage_coefficient;
        std::function<double()> get_lag_time;
        std::function<double&()> get_output;
        std::function<State&()> get_state;
        std::function<std::vector<double>&()> get_lag_array;
    };

    explicit ClarkUnitHydrograph(StorageView view);
    
    void execute();
    double change_storage_coefficient(double new_storage_coefficient);
    double change_lag_time(double new_lag_hours);
    double total_water_in_system() const;

private:
    StorageView view_;

    void initialize_lag_array();
    long calculate_lag_steps(double lag_hours) const;
    std::pair<double, double> calculate_coefficients() const;
    double calculate_c1(double storage_coefficient) const;
    double calculate_c2(double storage_coefficient) const;
    void redistribute_lag_storage(long new_lag_steps, double total_storage);
};
