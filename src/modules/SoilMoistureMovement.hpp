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
    explicit SoilMoistureMovement(const config_file& cfg);

    ~SoilMoistureMovement() override = default;

    void build_matrix(const mesh& domain);
    void run(mesh &domain) override;
    void init(mesh& domain) override;

    static constexpr auto NUM_NEIGHBOURS = 5;
    static constexpr auto NUM_LAYERS = 3;

    static std::string get_theta_name(const int layer) noexcept { return "theta" + std::to_string(layer); }
    static std::string get_psi_name(const int layer) noexcept { return "psi" + std::to_string(layer); }

private:
    struct Pair
    {
        double owner{};
        double neighbour{};
    };

    struct Geometry
    {
        Pair to_face;
        Pair elevation;
        double cell_centre_distance;
    };
    class faceInterpolator
    {
    public:
        explicit faceInterpolator(const Geometry&);
        ~faceInterpolator() = default;

        void set_geometry(const Geometry&);
        [[nodiscard]] double interp(Pair) const;
        [[nodiscard]] Pair elevation() const;
    private:
        const Geometry geometry_;
    };

    struct GeoHelper;
    class Boundary
    {
    public:
        const double DeltaZ;
        const double distance_to_face;
        explicit Boundary(double dz,double df);
    };
    class Interior
    {
    public:
        explicit Interior(const mesh_elem&,GeoHelper,const config_file&);
        const Geometry geometry;
        faceInterpolator interpolator{geometry};


    private:
        // Static Helpers
        static Geometry set_interior_geometry(const mesh_elem& face, GeoHelper geo_help, const config_file& cfg);
        static Geometry set_bottom_geo(const Point_3& cell_centre, Pair depth);
        static Geometry lower_layer_boundary(const mesh_elem& face, GeoHelper geo_help, const Point_3& face_centre,
                                         double lower_depth, double recharge_depth);
        static Geometry recharge_layer_boundary(const mesh_elem& face,
                                          GeoHelper geo_help,
                                          const Point_3& face_centre, double lower_depth,
                                          double recharge_depth, double detention_depth);
        static Geometry detention_layer_boundary(const mesh_elem& face,
                                                                       GeoHelper geo_help,
                                                                       const Point_3& face_centre, double recharge_depth,
                                                                       double detention_depth);
        static Geometry set_top_geo(const Point_3& cell_centre, Pair depth);
        static Geometry set_lateral_geo(const mesh_elem& face, GeoHelper geo_help,
                                                       const Point_3& cell_centre,
                                                       double centre_depth);
    };

    class data : public face_info
    {
    public:
        template<typename T>
        using opt = std::optional<T>;
        template<typename T>
        using LayerNeighbourArray = std::array<std::array<T,NUM_NEIGHBOURS>,NUM_LAYERS>;
        using faceType = std::variant<Boundary,Interior>;

        const LayerNeighbourArray<opt<faceInterpolator>> interp_to_face;
        const LayerNeighbourArray<faceType> cell_geometry;
        const LayerNeighbourArray<double> alpha;
        const LayerNeighbourArray<opt<int>> neighbour_idx_;
        explicit data(const LayerNeighbourArray<opt<faceInterpolator>>& interp,
            const LayerNeighbourArray<faceType>& cell,
            const LayerNeighbourArray<double>& alpha,
            const LayerNeighbourArray<opt<int>>& neighbour_idx,
            const config_file& cfg);

        double K_saturated;
        double air_entry_tension;
        double pore_size_dist_index;
        double theta(int) const;
        std::array<double,NUM_LAYERS> psi_n;
    };
    template <size_t NumLayers, size_t... Is>
    std::array<data::faceType, NumLayers> build_layers(mesh_elem&, size_t i, std::index_sequence<Is...>);

    template <size_t NumNeighbours = NUM_NEIGHBOURS, size_t NumLayers = NUM_LAYERS, size_t... Ns>
    std::array<std::array<data::faceType, NumNeighbours>, NumLayers>
    build_geometry(mesh_elem& face,std::index_sequence<Ns...>);

    struct Sizes
    {
        size_t local;
        size_t global;
        size_t vert_layers = 1;
    } sizes{};
    enum class Neighbour
    {
        Lateral_0 = 0, // switch cases will receive a number from 0-4, even if lateral cases are the same
        Lateral_1,
        Lateral_2,
        Top,
        Bottom
    };
    // Do this to double-check Neighbour will cast correctly
    static_assert(Neighbour::Bottom == static_cast<Neighbour>(4) &&
        static_cast<int>(Neighbour::Bottom) == 4,
        "Enum Class Neighbour, defined for clarity, not casting to the correct values");
    class solverData
    {
        mesh_elem& face;
        data& d;
        const int z_idx;
        const Sizes sizes;
        const std::string_view ID;

        double soil_water_capacity(int) const;
        double K_unsaturated(int) const;
    public:
        solverData(data& d, mesh_elem& face, int z, const Sizes& sizes, std::string_view);

        [[nodiscard]] size_t idx() const ;
        [[nodiscard]] bool has_neighbour(int f) const;
        [[nodiscard]] size_t neighbour_idx(int f) const;

        [[nodiscard]] static constexpr int top_face() { return static_cast<int>(Neighbour::Top); }
        [[nodiscard]] static constexpr int bottom_face() { return static_cast<int>(Neighbour::Bottom); }

        // without tags
        using donor_choice = math::without_donor_tag;
        using boundary_donor_choice = math::without_boundary_donor_tag;
        using boundary_choice = math::without_boundary_tag;

        // Normal faces
        using rhs_choice = math::with_rhs_tag;
        [[nodiscard]] double rhs(int) const;
        [[nodiscard]] double diagonal(int f) const;
        [[nodiscard]] double off_diagonal(int f) const;

        // Boundaries
        using side_boundary_choice = math::with_side_boundary_tag;
        [[nodiscard]] static double side_boundary_diagonal(int);
        [[nodiscard]] static double side_boundary_off_diagonal(int);
        [[nodiscard]] double side_boundary_rhs(int) const;

        using bottom_boundary_choice = math::with_bottom_boundary_tag;
        [[nodiscard]] static double bottom_boundary_diagonal(int);
        [[nodiscard]] static double bottom_boundary_off_diagonal(int);
        [[nodiscard]] double bottom_boundary_rhs(int) const;

        using top_boundary_choice = math::with_top_boundary_tag;
        [[nodiscard]] static double top_boundary_diagonal(int);
        [[nodiscard]] static double top_boundary_off_diagonal(int);
        [[nodiscard]] double top_boundary_rhs(int) const;

    };
    std::optional<math::LinearAlgebra::NearestNeighborProblem> moisture_content_solver;

    data::faceType get_lateral_boundary(const mesh_elem& face, size_t nn) const;
    data::faceType get_geometry(const mesh_elem& face,GeoHelper geo_helper) const;

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


    template <Indexable T>
    void write_output(mesh& domain, T& runoff_sol)
    {
#pragma omp parallel for
        for (size_t i = 0; i < domain->size_local_faces(); i++)
        {
            auto face = domain->face(i);
            auto& d = face->get_module_data<data>(ID);

            for (int layer = 0; layer < NUM_LAYERS; ++layer)
            {
                auto psi = runoff_sol[sizes.global * layer + face->cell_local_id];
                d.psi_n[layer] = psi;
                const std::string theta_name = get_theta_name(layer);
                (*face)[theta_name] = d.theta(layer);
                const std::string psi_name = get_psi_name(layer);
                (*face)[psi_name] = psi;
            }
        }

        for (int layer = 0; layer < NUM_LAYERS; ++layer)
        {
            // Don't need to communicate theta
            const std::string psi_name = get_psi_name(layer);
            domain->ghost_neighbors_communicate_variable(psi_name);
        }
    }


};
