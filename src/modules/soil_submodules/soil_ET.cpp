soil_ET::soil_ET(soil_ET_DTO& _DTO) : DTO(_DTO);
{

};

soil_ET::~soil_ET()
{

};

void soil_ET::run()
{
    double available_to_evap = DTO.potential_ET;
    if (DTO.depression_storage + DTO.soil_storage > 0.0)
    {
        available_to_evap *= DTO.depression_storage / 
                (DTO.depression_storage + DTO.soil_storage);
    }
    else
        available_to_evap = 0.0;

    if (DTO.depression_storage > 0.0 && available_to_evap > 0.0)
    {
        if (DTO.depression_storage > available_to_evap)
            DTO.depression_storage = std::max(0.0,DTO.depression_storage - available_to_evap);
        else
        {
            available_to_evap = DTO.depression_storage;
            depression_storage = 0.0;
        }
        DTO.actual_ET = available_to_evap;
    }
    else
        available_to_evap = 0.0;

    available_to_evap = DTO.potential_ET - available_to_evap;

    if (available_to_evap > 0.0 && DTO.soil_storage > 0.0 && DTO.ground_cover_type > 0)
    {
        double percent_available_lower;
        double percent_available_rechr;
        double ET_lower;
        double ET_rechr;
        double soil_lower_storage = DTO.soil_storage_max - DTO.soil_rechr_max;
        double soil_lower_max = DTO.soil_storage_max - DTO.soil_rechr_max;

        if ( soil_lower_max > 0.0 ) // soil_lower > 0.0
            pctl = soil_lower_storage / soil_lower_max;
        else
            pctl = 0.0;

        pctr = DTO.soil_rechr_storage / DTO.soil_rechr_max;

        etr = available_to_evap;

        etr = set_ET_layer(etr,percent,soil_type_rechr);

        if (etr > available_to_evap) 
        {
            etl = 0.0;
            etr = available_to_evap;
        }
        else
            etl = available_to_evap - etr;

        if (etl > 0.0)
        {
            etl = set_ET_layer(etl, pctl, soil_type_lower);
        }

        double ET = 0.0;
        
        switch (ground_cover_type)
        {
        case 0: // bare soil, no evap
            break; 
        case 1: // recharge layer only, crops 
            if (etr > DTO.soil_rechr_storage)
            {
                ET = DTO.soil_rechr_storage;
                DTO.soil_rechr_storage = 0.0;
            }
            else
            {
                DTO.soil_rechr_storage -= etr;
                ET = etr;
            }
            DTO.soil_storage -= etr;
            break;
        case 2: // all soil moisture, grasses & shrubs
            if (etr + etl >= DTO.soil_storage)
            {
                ET = DTO.soil_storage;
                DTO.soil_storage = 0.0;
                DTO.soil_rechr_storage = 0.0;
            }
            else
            {
                ET = etr + etl;
                DTO.soil_storage -= ET;

                DTO.soil_rechr_storage = std::max(DTO.soil_rechr_storage - etr, 0.0);
            }
            break;
        }

        DTO.actual_ET += ET;

        if (DTO.is_lake(DTO))
            DTO.actual_ET = DTO.potential_ET;
};

double soil_ET::set_ET_layer(double et, double& percent, int& type)
{
    switch type
    {
    case 1:
        if (percent < 0.25)
            et = 0.5 * percent * et;
        break;
    case 2:
        if (percent < 0.5)
            et = percent * et;
        break;
    case 3:
        if (percent <= 0.33)
            et = 0.5 * percent * et;
        else if (percent < 0.67)
            et = percent * et;
        break;
    default:
        // do nothing, use default above
        break;
    }

    return et;
};
