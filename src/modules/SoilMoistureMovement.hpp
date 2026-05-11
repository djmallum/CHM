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
#include "StencilAssembly.hpp"
#include "module_base.hpp"
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
template<typename T>
concept Indexable = requires(T obj, size_t idx) {
    { obj[idx] } -> std::convertible_to<size_t>;  // Return type convertible to size_t
};
class SoilMoistureMovement : public module_base
{
REGISTER_MODULE_HPP(SoilMoistureMovement)
public:
    SoilMoistureMovement(config_file cfg);

    ~SoilMoistureMovement() override;

    void build_matrix(mesh& domain);
    void run(mesh &domain) override;
    void init(mesh& domain) override;

private:

    class faceInterpolator
    {
    public:
        struct Pair
        {
            double owner{};
            double neighbour{};
        };
        faceInterpolator(Pair,double d);
        ~faceInterpolator() = default;

        [[nodiscard]] double interp(Pair) const;
    private:
        const Pair distance_to_face;
        const double cell_centre_distance;
    };

    static constexpr auto NUM_NEIGHBOURS = 5;
    class data : public face_info
    {
        const std::array<faceInterpolator,NUM_NEIGHBOURS> interpolator;
    public:
        data(const std::array<faceInterpolator,NUM_NEIGHBOURS>& interp);
        // TODO Persistent data here
        bool has_neighbour(int f);
        int neighbour_idx(int f);
    };

    struct Sizes
    {
        size_t local;
        size_t global;
        size_t vert_layers = 1;
    } sizes{};
    enum class Neighbour
    {
        Lateral_0 = 0,
        Lateral_1,
        Lateral_2,
        Top,
        Bottom
    };
    // Do this to double check Neighbour will cast correctly
    static_assert(Neighbour::Bottom == static_cast<Neighbour>(4) &&
        static_cast<int>(Neighbour::Bottom) == 4,
        "Enum Class Neighbour, defined for clarity, not casting to the correct values");
    class solverData
    {
        mesh_elem& face;
        data& d;
        const int z_idx;
        const Sizes sizes;
    public:
        solverData(data& d, mesh_elem& face, int z, Sizes sizes);

        [[nodiscard]] size_t idx() const ;
        [[nodiscard]] bool has_neighbour(int f) const;
        [[nodiscard]] size_t neighbour_idx(int f) const;
        [[nodiscard]] double diagonal(int f) const;
        [[nodiscard]] double off_diagonal(int f) const;
        [[nodiscard]] double boundary_diagonal(int f) const;

        [[nodiscard]] math::LinearAlgebra::VertLayer layer() const;
        [[nodiscard]] static constexpr int top_face() { return static_cast<int>(Neighbour::Top); }
        [[nodiscard]] static constexpr int bottom_face() { return static_cast<int>(Neighbour::Bottom); }
        [[nodiscard]] double bottom_boundary_diagonal() const;
        [[nodiscard]] double bottom_boundary_rhs() const;
        [[nodiscard]] double top_boundary_rhs_in() const;
    };
    std::optional<math::LinearAlgebra::NearestNeighborProblem> moisture_content_solver;

    auto try_solution(const mesh& domain)
    {
        try
        {
            return moisture_content_solver->Solve();
        }
        catch (const Belos::StatusTestError& e)
        {
            auto rank = 0;
#ifdef USE_MPI
            rank = domain->_comm_world.rank();
#endif
            SPDLOG_ERROR("Rank {}");
            std::string prefix = "Runoff.rank" + std::to_string(rank);
            moisture_content_solver->writeSolutionMatrixMarket(prefix);
            SPDLOG_ERROR(e.what());
            CHM_THROW_EXCEPTION(module_error, e.what());
        }
    }

    template<Indexable T>
    void write_output(mesh& domain, T& runoff_sol)
    {
#pragma omp parallel for
        for (size_t i = 0; i < domain->size_local_faces(); i++)
        {
            auto face = domain->face(i);
            auto& d = face->get_module_data<data>(ID);

            for (int z = 0; z < sizes.vert_layers; ++z)
            {
                auto theta = runoff_sol[sizes.global * z + face->cell_local_id];

                std::string output_name = "theta" + std::to_string(z);
                (*face)[output_name] = theta;
                domain->ghost_neighbors_communicate_variable(output_name);
            }
        }
    }


};
