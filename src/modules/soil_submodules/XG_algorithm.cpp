#include "XG_algorithm.hpp"

void XG_algorithm::run()
{

    if (DTO.is_day_start(DTO))
        B = 0.0;

    B += DTO.t_surface;

    if (STO.is_day_start(DTO))
    {
        // idle to freeze/thaw or freeze/thaw to idle
        
        // Keep T_total within +/- T_threshold
        // TODO might not be needed, might be diagnositic like B
        if (T_total > T_threshold)
            T_total = T_threshold;
        else if (T_total < -T_threshold)
            T_total = -T_threshold;

        if (std::abs(T_total) >= T_threshold/2.0 && freeze_thaw_state == 0)
        {
            if (freeze_front_depth < 0.0 || num_fronts)
                freeze_thaw_state = 1;
            else
                freeze_thaw_state = -1;

            temp_trend = 0.0;
        }

        if (freeze_thaw_state && std::abs(T_total) >= T_threshold/2.0 && std::abs(temp_trend) > 0.0)
        {
            freeze_thaw_state = 0;

            if (freeze_front_depth > 0.0 && thaw_front_depth > 0.0)
            {
                if ( is_freezing() )
                    switch_to_idle(thaw_front_depth);
                else
                    switch_to_idle(-freeze_front_depth);
            }
        }

        // Calculate thermal conductivities

        // Issues 

        thaw_low_layer = 1;
        freeze_low_layer = 1; 
        
        thaw_k_T_rechr = 
        freeze_k_T_rechr 
        for (int ii = 0; ii < N_layers; ++ii)
        {
            if (DTO.soil_storage_max > 0.0)
            {
                 
            }

        }    
    }






            
};

bool XG_algorithm::is_freezing()
{
    if (freeze_thaw_state == -1 && T_total > 0.0 && temp_trend > 0.0)
        return true;
    else if (freeze_thaw_state == 1 && T_total < 0.0 && temp_trend < 0.0)
        return false;
}

void XG_algorithm::switch_to_idle(double front_depth)
{
    push_front(front_depth); // TODO write this function
    // TODO in crhm both parts modify Zd_front_array (something I've decided not to include in XG) come back if necessary
    if (front_depth > 0.0 && thaw_front_depth > freeze_front_depth)
    {
        thaw_front_depth = 0.0;
        thaw_degree_days = 0.0;
    }
    else (front_depth < 0.0 && freeze_front_depth > thaw_front_depth)
    {
        freeze_front_depth = 0.0;
        freeze_degree_days = 0.0;
    }

};




