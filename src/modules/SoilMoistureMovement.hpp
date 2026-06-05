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

#include "LinearAlgebra.hpp"
#include "module_base.hpp"
#include "submodules/SoilMoistureSolverCore/SoilMoistureSolverCore.hpp"
#include "triangulation.hpp"

/**
 * \ingroup TODO KEYWORDS HERE
 * @{
 * \class SoilMoistureMovement
 *
 * TODO DESCRIPTION HERE
 *
 * **Depends:**
 * TODO DEPENDS HERE
 *
 * **Provides:**
 * TODO PROVIDES HERE (units in [])
 *
 * \rst
 * .. note::
 *  TODO ANY NOTES HERE
 *
 * \endrst
 *
 * **References:**
 * TODO REFERENCES AS NEEDED
 *
 * @}
 */
class SoilMoistureMovement : public module_base
{
REGISTER_MODULE_HPP(SoilMoistureMovement)
public:
    explicit SoilMoistureMovement(const config_file& cfg);

    ~SoilMoistureMovement() override = default;

    void run(mesh &domain) override;
    void init(mesh& domain) override;

private:
    SoilMoistureSolver::SoilMoistureSolverCore<mesh> _solver;
};
