#include "XG_algorithm.hpp"


void XG_algorithm::run()
{
    double L = 335000;
    
    if (P.is_newday)
       S.B = 0.0;

    S.B += surface_temp / P.time_step_per_day;

    S.TrigAcc += B;

    S.t_trend -= S.t_trend / 192;
    S.t_trend += B/192;

    // ensure that Zd_front is the right size here

    if (P.is_newday)
    {
        if (S.TrigAcc > P.Trigthrhld)
            S.TrigAcc = P.Trigthrhld;
        else if (S.TrigAcc < -P.Trigthrhld)
            S.TrigAcc = -P.Trigthrhold;

        if (S.TrigAcc >= P.Trigthrhld / 2.0 && S.TrigState == 0 && (S.Zdf > 0.0 || nfront > 0))
        {
            S.TrigState = 1;
            S.Zd_front[1] = -S.Zdf;
            S.t_trend = 0.0;
        }            

        if (S.TrigAcc <= -P.Trigthrhld / 2.0 && S.TrigState == 0)
        {
            S.TrigState = -1;
            S.Zd_front[1] = S.Zdt;
            S.t_trend = 0.0;
        }

        if (S.TrigState == -1 && S.TrigAcc >= P.Trigthrhld/2.0 && S.t_trend > 0.0)
        {
            S.TrigState = 0;

            if (S.Zdt > 0.0 && S.Zdf > 0.0)
            {
                if (S.Zdt > S.Zdf)
                {
                    push_front(S.Zdt);
                    S.Zdt = 0.0;
                    S.Bth = 0.0;
                    S.Zd_front[0] = 0.0;
                    S.Zd_front[1] = -S.Zdf;
                }
            }
        }

        if (S.TrigState == 1 && S.TrigAcc <= -P.Trigthrhld / 2.0 && S.t_trend < 0.0)
        {
            S.TrigState = 0;

            if (S.Zdf > 0.0 && S.Zdt > 0.0)
            {
                if (S.Zdf > S.Zdt)
                {
                    push_front(-S.Zdf);
                    S.Zdf = 0.0;
                    S.Bfr = 0.0;
                    S.Zd_front[0] = 0.0;
                    S.Zd_front[1] = S.Zdt;
                }
            }
        }

        // Calculate Thermal Conductivities

        S.Th_low = 1;
        S.Fz_low = 1;
        S.Check_XG_moist = 0.0;

        for (int layer = 0; layer < S.Zd_front.size(); ++layer)
        {
            if (P.soil_moist_max > 0.0)
                S.XG_moist[layer] = S.rechr_fract[layer] * S.XG_max[layer] * P.soil_rechr / P.soil_rechr_max
                    + S.moist_fract[layer] * S.XG_max[layer] * P.soil_lower / P.soil_lower_max;
            else
                S.XG_moist[layer] = 0.0;

            S.check_XG_moist += S.XG_moist[layer];

            S.XG_moist[layer] += S.default_fract[layer] * S.XG_max[layer] * P.theta_default[layer];

            S.theta[layer] = S.XG_moist[layer] / S.XG_max[layer];

            if (S.theta[layer] < P.theta_min)
                S.theta[layer] = P.theta_min;

            S.h2o[layer] = S.theta[layer] * S.por[layer] * 1000;
            
            if(P.k_update)
            { // change all layers to dynamic
                if(S.ftc_contents[layer] == 1) // always set somewhere else
                  S.ftc[layer] = get_ttc(layer);
                else
                  S.ftc[layer] = get_ftc(layer);

                if(ttc_contents[layer] == 1) // always set somewhere else
                  S.ttc[layer] = get_ftc(layer);
                else
                  S.ttc[layer] = get_ttc(layer);
            }
            else
            {
                S.ftc[layer] = get_ftc(layer);
                S.ttc[layer] = get_ttc(layer);

                S.ftc_contents[layer] = 0;
                S.ttc_contents[layer] = 0;
            }
        }
    
        for (long layer = 1; layer < S.Zd_front.size(); ++layer) 
        {
            S.pf[layer] = std::sqrt(P.ftc[layer-1]*S.h2o[layer]/(P.ftc[layer]*S.h2o[layer-1])); // water kg/m3
            S.pt[layer] = std::sqrt(P.ttc[layer-1]*S.h2o[layer]/(P.ttc[layer]*S.h2o[layer-1])); // water kg/m3
        }

        if(S.TrigState < 0.0) // handle freezing
        {
            S.Bfr -= S.B;  // Calculate the value of Degree-day for every day
            
            if(S.Bfr > 0.0)
            {
                freeze(); // XG-Algorithm - Freezing
                
                // check for thaw front
                if(S.Zdt > 0.0 && S.Zdf >= S.Zdt)
                {
                    if(S.nfront > 0)
                    {
                        double Last = last_front();
                        
                        if(Last < 0.0) // frozen front
                        {
                            S.Zdf = pop_front();
                            find_freeze_D(S.Zdf);
                            double Last = last_front();
                            
                            if(Last > 0.0) // thaw front
                            {
                                S.Zdt = pop_front();
                                find_thaw_D(S.Zdf);
                                S.Zd_front[1] = S.Zdt;
                            }
                            else if(Last < 0.0) // never two frozen fronts together
                            {
                               //CRHM has thos throw an error, never two frozen fronts 
                            }
                            else // no thaw front
                            {
                                S.Zdt = 0.0;
                                S.Bth = 0.0;
                                S.Zd_front[1] = 0.0;
                            }
                        }
                        else if(Last < 0.0) // never two freeze fronts together
                        {
                            // CRHM had another throw here, but I dont think its possible to reach this ever.
                        }
                        else // no thaw layer
                        {
                            S.Zdt = 0.0;
                            S.Bth = 0.0;
                            S.Zd_front[1] = 0.0;
                        }
                    }
                    else // no fronts
                    {
                        S.Zdt = 0.0;
                        S.Bth = 0.0;
                        S.Zd_front[1] = 0.0;
                    }
                }
            }
            S.Zd_front[0] = -S.Zdf;
        } // freezing handled
        else if(S.TrigState > 0) // Surface thawing lower ground frozen
        {
            S.Bth += S.B;  // Accumulate thawing degree-days

            if(S.Bth <= 0.0)
            {
                S.Bfr = S.Bth;
                S.Zdt = 0.0;
                S.Bth = 0.0;
            }
            else
            {
                thaw(); // XG-Algorithm - Thawing
                
                // check for freeze front
                if(S.Zdf > 0.0 && S.Zdt >= S.Zdf)
                {
                    if(S.nfront > 0)
                    {
                        double Last = last_front();
                        
                        if(Last > 0.0) // thaw front
                        {
                            S.Zdt = pop_front();
                            find_thaw_D(S.Zdt);
                            Last = last_front();
                            
                            if(Last < 0.0) // frozen front
                            {
                                S.Zdf = pop_front();
                                find_freeze_D(S.Zdf);
                                S.Zd_front[1] = -S.Zdf;
                            }
                            else if(Last > 0.0) // never two thaw fronts together
                            {
                                //CRHM throws exception here
                            }
                            else // no frozen front
                            {
                                S.Zdf = 0.0; // no frozen layers
                                S.Bfr = 0.0;
                                S.Zd_front[1] = 0.0;
                            }
                        }
                        else if(Last < 0.0) // never two freeze fronts together
                        {
                            //Again, another throw that cannot be accessed, can remove this i think                        
                        }
                        else // no thaw layer
                        {
                            S.Zdf = 0.0;
                            S.Bfr = 0.0;
                            S.Zd_front[1] = 0.0;
                        }
                    }
                    else // no fronts
                    {
                        S.Zdf = 0.0; // no frozen layer
                        S.Bfr = 0.0;
                        S.Zdt = 0.0;
                        S.Bth = 0.0;
                        S.TrigState = 0;
                        S.Zd_front[1] = 0.0;
                    }
                } // if check freeze front
            } // else thaw processing
            
            S.Zd_front[0] = S.Zdt;
        } // thawing handled
    }
};

//void XG_algorithm::run()
//{
//
//    if (DTO.is_day_start(DTO))
//        B = 0.0;
//
//    B += DTO.t_surface;
//
//    if (STO.is_day_start(DTO))
//    {
//        // idle to freeze/thaw or freeze/thaw to idle
//        
//        // Keep T_total within +/- T_threshold
//        // TODO might not be needed, might be diagnositic like B
//        if (T_total > T_threshold)
//            T_total = T_threshold;
//        else if (T_total < -T_threshold)
//            T_total = -T_threshold;
//
//        if (std::abs(T_total) >= T_threshold/2.0 && freeze_thaw_state == 0)
//        {
//            if (freeze_front_depth < 0.0 || num_fronts)
//                freeze_thaw_state = 1;
//            else
//                freeze_thaw_state = -1;
//
//            temp_trend = 0.0;
//        }
//
//        if (freeze_thaw_state && std::abs(T_total) >= T_threshold/2.0 && std::abs(temp_trend) > 0.0)
//        {
//            freeze_thaw_state = 0;
//
//            if (freeze_front_depth > 0.0 && thaw_front_depth > 0.0)
//            {
//                if ( is_freezing() )
//                    switch_to_idle(thaw_front_depth);
//                else
//                    switch_to_idle(-freeze_front_depth);
//            }
//        }
//
//        // Calculate thermal conductivities
//
//        // Issues 
//
//        thaw_low_layer = 1;
//        freeze_low_layer = 1; 
//         
//        double thaw_k_T_rechr = get_k_T(DTO.soil_rechr_storage);
//        double freeze_k_T_rechr = ;
//        double thaw_k_T_lower = ;
//        double freeze_k_T_lower = ;
//         
//        for (int ii = 0; ii < N_layers; ++ii)
//        {
//            if (DTO.soil_storage_max > 0.0)
//            {
//                 
//            }
//
//        }    
//    }
//
//
//
//
//
//
//            
//};
//
//bool XG_algorithm::is_freezing()
//{
//    if (freeze_thaw_state == -1 && T_total > 0.0 && temp_trend > 0.0)
//        return true;
//    else if (freeze_thaw_state == 1 && T_total < 0.0 && temp_trend < 0.0)
//        return false;
//}
//
//void XG_algorithm::switch_to_idle(double front_depth)
//{
//    push_front(front_depth);// TODO write this function
//    // TODO in crhm both parts modify Zd_front_array (something I've decided not to include in XG) come back if necessary
//    if (front_depth > 0.0 && thaw_front_depth > freeze_front_depth)
//    {
//        thaw_front_depth = 0.0;
//        thaw_degree_days = 0.0;
//    }
//    else (front_depth < 0.0 && freeze_front_depth > thaw_front_depth)
//    {
//        freeze_front_depth = 0.0;
//        freeze_degree_days = 0.0;
//    }
//
//};
//
//
//
//
