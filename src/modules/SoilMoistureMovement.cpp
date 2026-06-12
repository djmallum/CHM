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

#include "SoilMoistureMovement.hpp"

#include "SoilMoistureSolverCore/Details.hpp"
#include <boost/xpressive/detail/utility/traits_utils.hpp>
REGISTER_MODULE_CPP(SoilMoistureMovement)

SoilMoistureMovement::SoilMoistureMovement(const config_file& cfg)
    : module_base("SoilMoistureMovement", parallel::domain, cfg),
    _solver(ID,SoilMoistureSolver::param_builder(cfg,*global_param))
{
    // TODO Add depends/provides

    auto _provides = [this](HashName s)
    {
        this->provides(s);
    };

    auto _depends = [this](HashName s)
    {
        this->depends(s);
    };

    Solver::depends(_depends);
    Solver::provides(_provides);
}

void SoilMoistureMovement::init(mesh& domain)
{
    _solver.init(domain);
}

void SoilMoistureMovement::run(mesh& domain)
{
    // TODO maybe include the following if necessary
    // if(is_water(face))
    // {
    //     set_all_nan_on_skip(face);
    //     return;
    // }
    _solver.run(domain);
}
