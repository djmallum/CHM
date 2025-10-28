//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#include "physics/Atmosphere.h"

namespace Atmosphere
{
    // Logrithmic, assuming no snow cover and no canopy bewteen Z_in and Z_out
    // Because filter is before runtime, we do not know what the snowdepth will be, thus this
    // introduces some error latter when wind is scaled down taking into account the snowdepth.
    double log_scale_wind(double u, double Z_in, double Z_out, double snowdepthavg, double z0)
    {
        double Z0_SNOW  = z0; //Snow::Z0_SNOW; // Snow roughness (m)

        u = u * log((Z_out - (snowdepthavg + Z0_SNOW)) / Z0_SNOW) / log((Z_in - (snowdepthavg + Z0_SNOW)) / Z0_SNOW);
        return u;
    }

    // Exponential following Inoue (1963)
    double exp_scale_wind(double u, double Z_in, double Z_out, const double alpha)
    {

        u = u * exp(alpha*(Z_out/Z_in-1));
        return u;
    }

    // Correct precipitation input using triangle slope when input preciptation are given for the horizontally projected area.
    // See Fig 1 and 2 of Kienzle (2010, Hydrological Processes)  
    // Slope in radian
    double corr_precip_slope(double p, double slope)
    {
        p = p * cos(slope);
        return p;
    }

    /**
    * @brief Saturated vapour pressure
    * @param T air temperature (K)
    * @return Saturated vapour pressure (Pa)
    */
    double saturatedVapourPressure(const double T)
    {
        double TA = T - 273.15;
        double Es, E, Rhi, Rhw, Rh;                         //saturation and current water vapro pressure
        const double Aw = 611.21, Bw = 17.502, Cw = 240.97; //parameters for water
        const double Ai = 611.15, Bi = 22.452, Ci = 272.55; //parameters for ice
        const double Tfreeze = 0.;                          //freezing temperature

        if (T >= Tfreeze )
        {
            //above freezing point, water
            Es = Aw * exp((Bw * TA) / (Cw + TA));
        }
        else
        {
            Es = Ai * exp( (Bi * TA) / (Ci + TA) );
        }
        return Es;
    }

    /**
     * @brief Slope of saturated vapour pressure vs temperature
     * @param T air temperature (DEGREE_CELSIUS)
     * @return Delta (kPa/DEGREE_CELSIUS)
     */  
    double saturatedVapourPressure_slope(const double T) // Slope of sat vap p vs t, kPa/DEGREE_CELSIUS
    {
        if (T > 0.0)
            return(2504.0*exp(17.27 * T/(T+237.3)) / pow(T+237.3,2));
        else
            return(3549.0*exp( 21.88 * T/(T+265.5)) / pow(T+265.5,2));
    }

    /**
     * @brief latent heat of vaporization of water
     * @param T air temperature (DEGREE_CELSIUS)
     * @return latent heat of vaporization of water (J/kg)
     */ 
    double latent_heat_vapour_air(const double T) // Latent heat of vaporization (J/kg)
    {
        // Equation 7-8 Dingman (2002) Second Edition
        return (2.501 - 0.002361 * T) * 1e6; // original is MegaJoules/kg, 1e6 returns it to joules/kg
    }

    /**
     * @brief Psychrometric constant
     * @param P_a air pressure (kPa)
     * @param T air temperature (DEGREE_CELSIUS)
     * @return gamma (kPa/DEGREE_CELSIUS)
     */ 
    double psychrometric_constant(const double P_a, const double T, const double c_air) // Psychrometric constant (kPa/DEGREE_CELSIUS)
    {
        // Equation 7-13 Dingman Second Edition 2002
        return c_air * P_a / (0.622 * latent_heat_vapour_air(T)); // lambda (J/kg)
    }

    /**
     * @brief Density of air
     * @param T air_temperature (DEGREE_CELSIUS)
     * @param e_a vapour pressure (kPa)
     * @param P_a air pressure (kPa)
     */  
    double air_density(const double T,const double e_a,const double P_a)
    {
        static constexpr double R0 = 2870.0;
        return 1E4 * P_a / (R0 * (273.15 + T)) * (1.0 - 0.379*(e_a/P_a));
    }
};
