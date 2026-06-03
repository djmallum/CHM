//
// Created by Allum, Donovan on 2026-06-03.
//

#pragma once
#include "LinearAlgebra.hpp"
#include "StencilAssembly.hpp"
#include "Concepts.hpp"
#include "Details.hpp"
#include <spdlog/spdlog.h>
#include "Data.hpp"

namespace SoilMoistureSolver
{
namespace detail
{
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
}

using namespace detail;
template<Mesh mesh,Element mesh_elem>
class SoilMoistureSolverCore
{
public:
    explicit SoilMoistureSolverCore(const std::string& id,const Params& p) : _ID(id),_params(p) {}
    void run(mesh& domain);
    void init(mesh& domain);
private:
    const std::string _ID;
    const Params _params;
    std::optional<math::LinearAlgebra::NearestNeighborProblem> moisture_content_solver;
    Sizes sizes{};

    void build_matrix(mesh&);
    auto try_solution(const mesh& domain);
    template<Indexable T>
    void write_output(mesh& domain, T& runoff_sol);
};

template<Mesh mesh,Element mesh_elem>
void SoilMoistureSolverCore<mesh,mesh_elem>::build_matrix(mesh& domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_local_faces(); i++)
    {
        for (size_t z = 0; z < sizes.vert_layers; z++)
        {
            auto face = domain->face(i);
            auto& d = face->template get_module_data<data>(_ID);

            solverData solver_data(d, face, z, sizes,_ID);


            math::LinearAlgebra::assemble_all_neighbours<NUM_NEIGHBOURS>(*moisture_content_solver, solver_data);
        }
    }
}

template<Mesh mesh,Element mesh_elem>
auto SoilMoistureSolverCore<mesh,mesh_elem>::try_solution(const mesh& domain)
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
        const std::string prefix = "Runoff.rank" + std::to_string(rank);
        moisture_content_solver->writeSolutionMatrixMarket(prefix);
        SPDLOG_ERROR(e.what());
        CHM_THROW_EXCEPTION(module_error, e.what());
    }
}

static std::string get_theta_name(const int layer) noexcept { return "theta" + std::to_string(layer); }
static std::string get_psi_name(const int layer) noexcept { return "psi" + std::to_string(layer); }

template<Mesh mesh, Element mesh_elem>
template <Indexable T>
void SoilMoistureSolverCore<mesh, mesh_elem>::write_output(mesh& domain, T& runoff_sol)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_local_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->template get_module_data<data>(_ID);

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
template <Mesh mesh, Element mesh_elem>
void SoilMoistureSolverCore<mesh, mesh_elem>::run(mesh& domain)
{
    // TODO maybe include the following if necessary
    // if(is_water(face))
    // {
    //     set_all_nan_on_skip(face);
    //     return;
    // }
    build_matrix(domain);

    auto [numIters, residual] = try_solution(domain);
    const std::string warn = std::format("Richards Eqn Solve. Iterations: {}, residual: {}", numIters, residual);
    SPDLOG_DEBUG(warn);

    auto runoff_sol = moisture_content_solver->getSolutionView();

    write_output(domain, runoff_sol);
}

template <Mesh mesh, Element mesh_elem>
void SoilMoistureSolverCore<mesh, mesh_elem>::init(mesh& domain)
{
    /*
     * Since this is made to work directly with a soil module with specific parameters,
     * vert_layers is a constant
     */
    sizes.local = domain->size_local_faces();
    sizes.global = domain->size_global_faces();
    sizes.vert_layers = NUM_LAYERS;//cfg.get<int>("num_soil_layers");
    moisture_content_solver.emplace(domain, sizes.vert_layers);

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_local_faces(); i++)
    {
        auto face = domain->face(i);
        auto geometry = build_geometry<NUM_NEIGHBOURS,NUM_LAYERS>(face,
                std::make_index_sequence<NUM_LAYERS>{});

        LayerNeighbourArray<opt<faceInterpolator>> face_interp;
        LayerNeighbourArray<double> alpha{};
        LayerNeighbourArray<double> face_area{};
        LayerNeighbourArray<opt<int>> neighbour_idx{};

        static_assert(geometry.size() == NUM_LAYERS && geometry[0].size() == NUM_NEIGHBOURS,
            "geometry matrix built improperly");

        std::array<double,NUM_LAYERS> volume;
        const double tri_area = face->get_area();
        const std::array depths{
            _params.lower_soil_depth,
            _params.recharge_soil_depth,
            _params.detention_depth
        };

        for (size_t layer = 0; layer < NUM_LAYERS; layer++)
        {
            volume.at(layer) = tri_area * depths.at(layer);
        }

        for (const auto neighbour : all_neighbours) // size_t nn = 0; nn < NUM_NEIGHBOURS; ++nn)
        {
            for (size_t layer = 0; layer < NUM_LAYERS; ++layer)
            {
                const auto nn = static_cast<size_t>(neighbour);
                std::visit([&face_interp,nn,layer]<typename T0>(T0&& arg)
                {
                    using T = std::decay_t<T0>;
                    if constexpr (std::is_same_v<T, Interior>)
                    {
                        face_interp[layer][nn].emplace(arg.geometry);
                    }

                }, geometry[layer][nn]);


                switch (neighbour)
                {
                case Neighbour::Lateral_0:
                case Neighbour::Lateral_1:
                case Neighbour::Lateral_2:
                {
                    const auto side_length = face->edge_length(nn);
                    face_area[layer][nn] = side_length * depths[layer];
                    //neighbour
                    neighbour_idx[layer][nn] =
                        std::visit([layer,nn,this,&face]<typename T0>(T0&& arg) -> std::optional<int>
                        {
                            using T = std::decay_t<T0>;
                            if constexpr(std::is_same_v<T, Interior>)
                            {
                                return std::optional<int>(sizes.global * layer + face->neighbor(nn)->cell_global_id);
                            }
                            else
                                return std::nullopt;
                        },geometry[layer][nn]);
                    break;
                }
                case Neighbour::Top:
                    face_area[layer][nn] = tri_area;
                    //neighbour
                    neighbour_idx[layer][nn] =
                        std::visit([layer,this,&face]<typename T0>(T0&& arg) -> std::optional<int>
                        {
                            using T = std::decay_t<T0>;
                            if constexpr(std::is_same_v<T, Interior>)
                            {
                                if (layer == NUM_LAYERS - 1)
                                {
                                    const std::string err = std::format("geometry object is the Interior faceType variant but"
                                        "layer was allowed to be at the top: {}."
                                        "Should be lower than {}",layer,NUM_LAYERS - 1);
                                    CHM_THROW_EXCEPTION(module_error,err);
                                }

                                return std::optional<int>(this->sizes.global * (layer + 1) + face->cell_global_id);
                            }
                            else
                                return std::nullopt;
                        },geometry[layer][nn]);
                    break;
                case Neighbour::Bottom:
                    face_area[layer][nn] = tri_area;
                    //neighbour
                    neighbour_idx[layer][nn] =
                        std::visit([layer,this,&face]<typename T0>(T0&& arg) -> std::optional<int>
                        {
                            using T = std::decay_t<T0>;
                            if constexpr(std::is_same_v<T, Interior>)
                            {
                                if (layer == 0)
                                {
                                    const std::string err = std::format("geometry object is the Interior faceType variant but"
                                        "layer was allowed to be at the bottom: {}."
                                        "Should be higher than 0.",layer);
                                    CHM_THROW_EXCEPTION(module_error,err);
                                }
                                return std::optional<int>(this->sizes.global * (layer - 1) + face->cell_global_id);
                            }
                            else
                                return std::nullopt;
                        },geometry[layer][nn]);
                    break;
                }

                alpha[layer][nn] = _params.time_step_seconds * face_area[layer][nn] / volume[layer];
                alpha[layer][nn] *= std::visit([]<typename T0>(T0&& arg)
                {
                    using T = std::decay_t<T0>;
                    if constexpr (std::is_same_v<T, Interior>)
                    {
                        return arg.geometry.cell_centre_distance;
                    }
                    else if constexpr (std::is_same_v<T, Boundary>)
                    {
                        // There exists a mythical neighbour past the boundary
                        return 2.0 * arg.distance_to_face;
                    }
                },geometry[layer][nn]);
            }
        }

        face->template make_module_data<data>(_ID,face_interp, geometry, alpha, neighbour_idx, _params);


    }

    // TODO set domain parameters
}
} // namespace SoilMoistureSolver