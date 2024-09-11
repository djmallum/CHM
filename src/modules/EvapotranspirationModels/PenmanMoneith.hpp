#include "evapbase.hpp"

class PenmanMoneith : public evapbase
{
public:

    PenmanMoneith(); // TODO Add inputs here to constructor
                     
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
}


