#include "evapbase.hpp"
#include "../module_base.hpp"

class PenmanMoneith : public evapbase
{
public:

    PenmanMoneith(config_file cfg); // TODO Add inputs here to constructor
                     
    ~PenmanMoneith(); // Deconstructor
                     

    virtual double CalcEvap() const override;

private:

    double Veg_height;
    double Veg_heigh_max;
    double wind_speed;
    double wind_measurement_height;
    double kappa;
    double stomatal_resistance_min;
    double ShortWave_in;
    double t;
    double vapour_pressure;
    double soil_storage;
    double soil_depth;
    double Rnet;
    double Frac_to_ground;
    double aero_resistance;
    double saturation_vapour_pressure;
    double stomatal_resistance;

    void CalcHeights();
    void CalcAeroResistance();
    void CalcSaturationVapourPressure();
    void CalcStomatalResistance();

}


