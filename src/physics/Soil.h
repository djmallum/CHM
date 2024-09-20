/* * Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
 * modular unstructured mesh based approach for hydrological modelling
 * Copyright (C) 2018 Christopher Marsh
 *
 * This file is part of Canadian Hydrological Model.
 *
 * Canadian Hydrological Model is free software: you can redistribute it and/or
 * modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Canadian Hydrological Model is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Canadian Hydrological Model.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

//
// Created by nicway on 12/07/16.
//

#pragma once

#include <iostream>

namespace Soil {
    /********* Soil ************/

    enum GATable {PSI, KSAT, WILT, FCAP, PORG, PORE, AIRENT, PORESZ, AVAIL}; // Used for mapping the soil table, PSI and KSAT are used, the others are unused but may but used in the future or other modules.    
    // see GreenAmpt module in CRHM wiki for details.

    class SoilData
    {
    public:

        SoilData();
        ~SoilData();

        double get_soilproperties(const std::size_t soil_type, const std::size_t property_type);
        double get_textureproperties(const std::size_t texture_type, const std::size_t groundcover_type);
    private:
            
        static constexpr double soilproperties[13][9] = {
            { 0.0,  999.9, 0.000, 0.00, 1.100,  1.000,	0.000,	0.0,  4},  //      0  water
            { 49.5, 117.8, 0.020, 0.10, 0.437,  0.395,	0.121,	4.05, 1},  //      1  sand
            { 61.3,  29.9, 0.036, 0.16, 0.437,  0.41 ,	0.09,	4.38, 4},  //      2  loamsand
            {110.1,  10.9, 0.041, 0.23, 0.453,  0.435,	0.218,	4.9,  2},  //      3  sandloam
            { 88.9,   3.4, 0.029, 0.26, 0.463,  0.451,	0.478,	5.39, 2},  //      4  loam
            {166.8,   6.5, 0.045, 0.38, 0.501,  0.485,	0.786,	5.3,  2},  //      5  siltloam
            {218.5,   1.5, 0.068, 0.38, 0.398,  0.420,	0.299,	7.12, 3},  //      6  saclloam
            {208.8,   1.0, 0.155, 0.39, 0.464,  0.476,	0.63,	8.52, 2},  //      7  clayloam
            {273.3,   1.0, 0.039, 0.40, 0.471,  0.477,	0.356,	7.75, 2},  //      8  siclloam
            {239.0,   0.6, 0.110, 0.41, 0.430,  0.426,	0.153,	10.4, 3},  //      9  sandclay
            {292.2,   0.5, 0.056, 0.43, 0.479,  0.492,	0.49,	10.4, 3},  //      10 siltclay
            {316.3,   0.3, 0.090, 0.46, 0.475,  0.482,	0.405,	11.4, 3},  //      11 clay
            {  0.0,   0.0, 0.000, 0.00, 0.000,  0.000,	0.0,	 0.0, 4}   //      12 pavement. Values not used
            };

        static constexpr double textureproperties[4][6] = { // mm/hour
            {7.6, 12.7, 15.2, 17.8, 25.4, 76.2},  // coarse over coarse
            {2.5,  5.1,  7.6, 10.2, 12.7,  15.2}, // medium over medium
            {1.3,  1.8,  2.5,  3.8,  5.1,  6.4},  // medium/fine over fine
            {0.5,  0.5,  0.5,  0.5,  0.5,  0.5}   // soil over shallow bedrock
            };

    };
}

