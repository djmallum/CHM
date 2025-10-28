#include "test_PenmanMonteith.hpp"
#include "Atmosphere.h"
#include <iostream>

//double MockPenmanData::leaf_area_index() { return LAI; }
double MockPenmanData::leaf_area_index_max() const { return LAImax; }
double MockPenmanData::short_wave_in() const { return Qsw; }
double MockPenmanData::Veg_height() const { return veg_Ht; }
double MockPenmanData::wind_measurement_height() const { return wind_height; }
double MockPenmanData::stomatal_resistance_min() const { return stomatal_res_min; }
bool MockPenmanData::has_vegetation() const { return Veg_height() > 0.0; }
double MockPenmanData::Q_net() const { return _Q_net; }
double MockPenmanData::Q_g() const { return _Q_net * F_to_g; }
int MockPenmanData::s_per_time_step() const { return _s_per_step; }
double MockPenmanData::heat_capacity_air() const { return Cp; }
double MockPenmanData::kappa() const { return K; }
double MockPenmanData::air_entry_tension() const { return tension; }
double MockPenmanData::pore_size_dist() const { return pore_sz; }
double MockPenmanData::porosity() const { return phi; }
double MockPenmanData::air_temperature() const { return _air_temperature; }
double MockPenmanData::d() const { return Veg_height()*0.67; }
double MockPenmanData::Z0() const { return Veg_height()/7.6; }
double MockPenmanData::wind_speed() const { return _wind_speed; }
double MockPenmanData::saturated_vapour_pressure() const { return ea_star;};
double MockPenmanData::vapour_pressure() const { return ea;};
double MockPenmanData::volumetric_moisture_content() const { return soil_storage/soil_d/1000.0 + theta_pwp;};
double MockPenmanData::delta() const { return evapT_base::delta(air_temperature()); };
double MockPenmanData::air_density() const { return Atmosphere::air_density(air_temperature(),vapour_pressure(),P_atm);};
double MockPenmanData::gamma() const { return evapT_base::gamma(P_atm,air_temperature(),heat_capacity_air());};
double MockPenmanData::lambda() const { return  evapT_base::lambda(air_temperature());};
double MockPenmanData::stomatal_resistance() const { return _stomatal_resistance; };
void MockPenmanData::stomatal_resistance(const double T) const { _stomatal_resistance = T; }
void MockPenmanData::ET(const double T) const { _ET = T; }
