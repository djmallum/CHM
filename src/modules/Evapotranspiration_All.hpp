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

#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "Atmosphere.h"
#include "Soil.h"
#include <cstdlib>
#include <string>
#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>
#include "EvapotranspirationModels/evapbase.hpp"
#include "EvapotranspirationModels/PenmanMonteith.hpp"

/**
 * \ingroup modules exp evap
 * @{
 * \class Evapotranspiration_All

 * Calculates evapo-transpiration via Penman-Monteith.
 *
 * Not currently maintained.
 *
 * Depends:
 * - Incoming shortwave radaition "iswr" [\f$ W \cdot m^{-2} \f$ ]
 * - Incoming longwave radiation "ilwr" [\f$ W \cdot m^{-2} \f$ ]
 * - Relative humidity "rh" [%]
 * - Windspeed at 2m "U_2m_above_srf" [\f$ m \cdot s^{-1} \f$ ]
 *
 * Provides:
 * - Evapotranspiration "ET" [\f$ mm \cdot dt^{-1} \f$ ]
 *
 * @}
 */
class Evapotranspiration_All : public module_base
{
REGISTER_MODULE_HPP(Evapotranspiration_All);
public:
    Evapotranspiration_All(config_file cfg);
    ~Evapotranspiration_All();
    void init(mesh& domain)
    void run(mesh_elem& face);

    class data : public face_info
    {
        public:
            std::unique_ptr<evapT_base> MyPenmanMonteith; //TODO should be a pointer to the base class
            // PriestleyTaylor here

    }

private:

    void init_PenmanMonteith(Evapotranspiration_All::data& d);
    PM_vars set_PenmanMonteith_vars(mesh_elem& face);
    std::unique_ptr<Soils::_soils_base> SoilDataObj;

};