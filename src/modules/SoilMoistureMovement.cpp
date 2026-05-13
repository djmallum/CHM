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
REGISTER_MODULE_CPP(SoilMoistureMovement);

SoilMoistureMovement::SoilMoistureMovement(config_file cfg) : module_base("SoilMoistureMovement", parallel::domain, cfg)
{
    // TODO Add depends/provides
}

SoilMoistureMovement::~SoilMoistureMovement()
{
    // Do nothing
}

void SoilMoistureMovement::build_matrix(mesh& domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_local_faces(); i++)
    {
        for (size_t z = 0; z < sizes.vert_layers; z++)
        {
            auto face = domain->face(i);
            auto& d = face->get_module_data<data>(ID);

            solverData solver_data(d, face, z, sizes);

            math::LinearAlgebra::lateral_neighbours(*moisture_content_solver, solver_data);

            math::LinearAlgebra::vertical_neighbours(*moisture_content_solver, solver_data);
        }
    }
}


void SoilMoistureMovement::run(mesh& domain)
{
    // TODO maybe include the following if necessary
    // if(is_water(face))
    // {
    //     set_all_nan_on_skip(face);
    //     return;
    // }
    build_matrix(domain);

    auto results = try_solution(domain);
    SPDLOG_DEBUG("Moisture Content. Iterations: {}, residual: {}", results.numIters, results.residual);

    auto runoff_sol = moisture_content_solver->getSolutionView();

    write_output(domain, runoff_sol);
}
void SoilMoistureMovement::init(mesh& domain)
{
    /*
     * Since this is made to work directly with a soil module with specific parameters,
     * vert_layers is a constant
     */
    sizes.local = domain->size_local_faces();
    sizes.global = domain->size_global_faces();
    sizes.vert_layers = 3;//cfg.get<int>("num_soil_layers");
    double recharge_depth = cfg.get<double>("recharge_depth");
    double lower_depth = cfg.get<double>("lower_depth");
    double detention_depth = cfg.get<double>("detention_depth");
    moisture_content_solver.emplace(domain, sizes.vert_layers);

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_local_faces(); i++)
    {
        auto face = domain->face(i);
        for (size_t z = 0; z < sizes.vert_layers; z++)
        {
            // TODO set centre to face and centre to centre distances for each neighbour
            // TODO loop over neighbours
            std::array<std::optional<faceInterpolator::Pair>,NUM_NEIGHBOURS> neighbours;

            for (size_t nn = 0; nn < NUM_NEIGHBOURS; nn++)
            {
                using Pair = faceInterpolator::Pair;
                auto layer_half_depth = layer_depth(z) / 2.0;
                auto cell_centre = face->center();
                cell_centre.z() -= layer_half_depth;
                // TODO handle faces at the edge of their process,
                // currently just leaving the std::optional in neighbours uninitialized
                switch (static_cast<Neighbour>(nn))
                {
                case Neighbour::Lateral_0:
                case Neighbour::Lateral_1:
                case Neighbour::Lateral_2:
                {
                    auto edge = face->edge_midpoint<Point_3>(nn);
                    edge.z() -= layer_half_depth;
                    double owner = CGAL::sqrt(CGAL::squared_distance(cell_centre,edge));

                    auto centre_neighbour = face->neighbor(nn)->center();
                    centre_neighbour.z() -= layer_half_depth;

                    double neighbour = CGAL::sqrt(CGAL::squared_distance(centre_neighbour,edge));
                    Pair pair;
                    pair.owner = owner;
                    pair.neighbour = neighbour;
                    if (neighbours.size() > nn)
                    {
                        std::string err = std::format("Error: Allocation of pairs proceeding in the wrong order"
                            "at nearest neighbour {} in layer {}",nn,z);
                        CHM_THROW_EXCEPTION(module_error,err);
                    }
                    neighbours[nn].emplace(pair);
                    break;
                }
                case Neighbour::Top:
                    if (z != sizes.vert_layers - 1)
                    {
                        neighbours[nn].emplace(Pair(layer_half_depth,layer_half_depth));
                    }
                    break;
                case Neighbour::Bottom:
                    if (z != 0)
                    {
                        neighbours[nn].emplace(Pair(layer_half_depth,layer_half_depth));
                    }

                    break;
                //default:
                //    std::string err = std::format("Index {} exceeds the number of neighbours"
                //        "in a mesh of triangular prisms, stacked vertically.",nn);
                //    CHM_THROW_EXCEPTION(module_error,err);
                }

            }
        }

        // TODO convert std::array<Pairs,NUM_NEIGHBOURS> to faceInterpolator for data constructor
        auto& d = face->make_module_data<data>(ID);



        // TODO set per-triangle parameters
    }

    // TODO set domain parameters
}
SoilMoistureMovement::faceInterpolator::faceInterpolator(const Pair p, double d)
    : distance_to_face(p), cell_centre_distance(d) {}

double SoilMoistureMovement::faceInterpolator::interp(const Pair p) const
{
    return cell_centre_distance * p.owner * p.neighbour /
        (distance_to_face.owner * p.neighbour + distance_to_face.neighbour * p.owner);
}

SoilMoistureMovement::solverData::solverData(data& d, mesh_elem& face, const int z, const Sizes size)
    : d(d), face(face), z_idx(z), sizes(size)
{
}

size_t SoilMoistureMovement::solverData::idx() const { return sizes.global * z_idx + face->cell_global_id; }
bool SoilMoistureMovement::solverData::has_neighbour(int f) const { return d.has_neighbour(f); }
size_t SoilMoistureMovement::solverData::neighbour_idx(int f) const { return d.neighbour_idx(f); }
double SoilMoistureMovement::solverData::diagonal(int f) const
{
    return 1.0 / NUM_NEIGHBOURS + d.alpha[f] / soil_moisture_capacity(f) * K_unsat(f); /* TODO everything about geometry or constant
                                                                      *  in time goes in alpha, could make it a type
                                                                      */
}
double SoilMoistureMovement::solverData::off_diagonal(int f) const
{
    return -\alpha[f] / soil_moisture_capacity(f) * K_unsat(f); // TODO see diagonal
    // TODO off diagonal contributions
}
double SoilMoistureMovement::solverData::boundary_diagonal(int f) const
{
    return 1.0 / NUM_NEIGHBOURS; /*
                                  * no flux, so contribution from psi[n+1] and psi_j[n+1] terms are zero
                                  * Remaining contribution is from the psi[n+1] from the time derivative
                                  */
}
math::LinearAlgebra::VertLayer SoilMoistureMovement::solverData::layer() const { return z_idx; }
double SoilMoistureMovement::solverData::bottom_boundary_diagonal() const
{
    // TODO another boundary condition at the bottom
}
double SoilMoistureMovement::solverData::bottom_boundary_rhs() const
{
    // TODO Dirichlet boundary condition
}
double SoilMoistureMovement::solverData::top_boundary_rhs_in() const
{
    // TODO again, Dirichlet boundary condition
}
