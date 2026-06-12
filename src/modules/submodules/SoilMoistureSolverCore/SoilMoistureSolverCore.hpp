#pragma once
#include "LinearAlgebra.hpp"
#include "StencilAssembly.hpp"
#include "Concepts.hpp"
#include "Details.hpp"
#include <spdlog/spdlog.h>

#include <utility>
#include "Data.hpp"
#include "Adaptor.hpp"

#include <array>

namespace SoilMoistureSolver
{


using namespace detail;

static Params param_builder(const config_file& cfg, const global& g)
{
    Params p{};
    p.recharge_soil_depth = cfg.get<double>("recharge_soil_depth");
    p.lower_soil_depth = cfg.get<double>("lower_soil_depth");
    p.detention_depth = cfg.get<double>("detention_depth");
    p.time_step_seconds = g.dt();

    return p;
}

template<MeshInterface M>
class SoilMoistureSolverCore
{
public:
    explicit SoilMoistureSolverCore(std::string id,const Params& p) : ID(std::move(id)),_params(p) {}
    void run(M& domain);
    //template <ElementInterface E> cellInfo<cellFacesAndVerticalLayers> build_cell_geometry(E& face);
    void init(M& domain);

    template<typename F> requires std::invocable<F,HashName>
    static void depends(F&& f);

    template<typename F> requires std::invocable<F,HashName>
    static void provides(F&& f);
private:
    const std::string ID;
    const Params _params;
    std::optional<math::LinearAlgebra::NearestNeighborProblem> moisture_content_solver;
    Sizes sizes{};

    void build_matrix(M&);
    auto try_solution(const M& domain);
    template<Indexable T>
    void write_output(M& domain, T& runoff_sol);
};

template<MeshInterface M>
void SoilMoistureSolverCore<M>::build_matrix(M& domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_local_faces(); i++)
    {
        for (size_t z = 0; z < sizes.vert_layers; z++)
        {
            auto face = domain->face(i);
            auto& d = face->template get_module_data<data>(ID);

            solverData solver_data(d, face, z, sizes,ID);

            math::LinearAlgebra::assemble_all_neighbours<cellFacesAndVerticalLayers.face>(*moisture_content_solver, solver_data);
        }
    }
}

template<MeshInterface M>
auto SoilMoistureSolverCore<M>::try_solution(const M& domain)
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

struct Name
{
    HashName psi;
    HashName theta;
};
static const std::array<Name, cellFacesAndVerticalLayers.layer>& build_output_names()
{
    static constexpr std::array<Name,cellFacesAndVerticalLayers.layer> names = []
    {
        std::array n{
            Name{"theta_lower"_s,"psi_lower"_s},
            Name{"theta_recharge"_s,"psi_recharge"_s},
            Name{"theta_detention"_s,"psi_detention"_s}
        };
        return n;
    }();

    return names;
}

template<MeshInterface M>
template <Indexable T>
void SoilMoistureSolverCore<M>::write_output(M& domain, T& runoff_sol)
{
    const auto& output_names = build_output_names();
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_local_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->template get_module_data<data>(ID);

        for (int layer = 0; layer < output_names.size(); ++layer)
        {
            auto psi = runoff_sol[sizes.global * layer + face->cell_local_id];
            d.psi_n[layer] = psi;
            (*face)[output_names.at(layer).theta] = d.theta(layer);
            (*face)[output_names.at(layer).psi] = psi;
        }
    }

    for (int layer = 0; layer < cellFacesAndVerticalLayers.layer; ++layer)
    {
        // Don't need to communicate theta
        domain->ghost_neighbors_communicate_variable(output_names.at(layer).psi);
    }
}
template <MeshInterface M>
void SoilMoistureSolverCore<M>::run(M& domain)
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

//template <MeshInterface M>
//template <ElementInterface E> cellInfo<dim_size> SoilMoistureSolverCore<M>::build_cell_geometry(E& face)
//{
//    auto geometry =
//        build_geometry<E, dim_size.face, dim_size.layer>(face, _params, std::make_index_sequence<dim_size.layer>{});
//
//    LayerNeighbourArray<opt<faceInterpolator>> face_interp;
//    LayerNeighbourArray<double> alpha{};
//    LayerNeighbourArray<double> face_area{};
//    LayerNeighbourArray<opt<int>> neighbour_idx{};
//
//    static_assert(geometry.size() == dim_size.layer && geometry[0].size() == dim_size.face,
//                  "geometry matrix built improperly");
//
//    std::array<double, dim_size.layer> volume;
//    const double tri_area = face->get_area();
//    const std::array depths{_params.lower_soil_depth, _params.recharge_soil_depth, _params.detention_depth};
//
//    for (size_t layer = 0; layer < dim_size.layer; layer++)
//    {
//        volume.at(layer) = tri_area * depths.at(layer);
//    }
//
//    for (const auto neighbour : all_neighbours) // size_t nn = 0; nn < dim_size.face; ++nn)
//    {
//        for (size_t layer = 0; layer < dim_size.layer; ++layer)
//        {
//            const auto nn = static_cast<size_t>(neighbour);
//            std::visit(
//                [&face_interp, nn, layer]<typename T>(T&& arg)
//                {
//                    if constexpr (std::is_same_v<T, Interior>)
//                    {
//                        face_interp[layer][nn].emplace(arg.geometry);
//                    }
//                },
//                geometry[layer][nn]);
//
//            switch (neighbour)
//            {
//            case Neighbour::Lateral_0:
//            case Neighbour::Lateral_1:
//            case Neighbour::Lateral_2:
//            {
//                const auto side_length = face->edge_length(nn);
//                face_area[layer][nn] = side_length * depths[layer];
//                // neighbour
//                neighbour_idx[layer][nn] = std::visit(
//                    [layer, nn, this, &face]<typename T>(T&&) -> std::optional<int>
//                    {
//                        if constexpr (std::is_same_v<T, Interior>)
//                        {
//                            return std::optional<int>(sizes.global * layer + face->neighbor(nn)->cell_global_id);
//                        }
//                        else
//                            return std::nullopt;
//                    },
//                    geometry[layer][nn]);
//                break;
//            }
//            case Neighbour::Top:
//                face_area[layer][nn] = tri_area;
//                // neighbour
//                neighbour_idx[layer][nn] = std::visit(
//                    [layer, this, &face]<typename T>(T&&) -> std::optional<int>
//                    {
//                        if constexpr (std::is_same_v<T, Interior>)
//                        {
//                            if (layer == dim_size.layer - 1)
//                            {
//                                const std::string err =
//                                    std::format("geometry object is the Interior faceType variant but"
//                                                "layer was allowed to be at the top: {}."
//                                                "Should be lower than {}",
//                                                layer, dim_size.layer - 1);
//                                CHM_THROW_EXCEPTION(module_error, err);
//                            }
//
//                            return std::optional<int>(this->sizes.global * (layer + 1) + face->cell_global_id);
//                        }
//                        else
//                            return std::nullopt;
//                    },
//                    geometry[layer][nn]);
//                break;
//            case Neighbour::Bottom:
//                face_area[layer][nn] = tri_area;
//                // neighbour
//                neighbour_idx[layer][nn] = std::visit(
//                    [layer, this, &face]<typename T>(T&&) -> std::optional<int>
//                    {
//                        if constexpr (std::is_same_v<T, Interior>)
//                        {
//                            if (layer == 0)
//                            {
//                                const std::string err =
//                                    std::format("geometry object is the Interior faceType variant but"
//                                                "layer was allowed to be at the bottom: {}."
//                                                "Should be higher than 0.",
//                                                layer);
//                                CHM_THROW_EXCEPTION(module_error, err);
//                            }
//                            return std::optional<int>(this->sizes.global * (layer - 1) + face->cell_global_id);
//                        }
//                        else
//                            return std::nullopt;
//                    },
//                    geometry[layer][nn]);
//                break;
//            }
//
//            alpha[layer][nn] = _params.time_step_seconds * face_area[layer][nn] / volume[layer];
//            alpha[layer][nn] *= std::visit(
//                []<typename T>(T&& arg) -> double
//                {
//                    if constexpr (std::is_same_v<T, Interior>)
//                    {
//                        return arg.geometry.cell_centre_distance;
//                    }
//                    else if constexpr (std::is_same_v<T, Boundary>)
//                    {
//                        // There exists a mythical neighbour past the boundary
//                        return 2.0 * arg.distance_to_face;
//                    }
//
//                    CHM_THROW_EXCEPTION(module_error, "Visiting to invalid type, should be unreachable");
//                },
//                geometry[layer][nn]);
//        }
//    }
//    return cellInfo{face_interp, geometry, alpha, neighbour_idx};
//}
template <MeshInterface M>
void SoilMoistureSolverCore<M>::init(M& domain)
{
    /*
     * Since this is made to work directly with a soil module with specific parameters,
     * vert_layers is a constant
     *
     * TODO Extract this data into a separate class and then this goes into the constructor
     * Maybe not allll of it.
     */
    sizes.local = domain->size_local_faces();
    sizes.global = domain->size_global_faces();
    sizes.vert_layers = cellFacesAndVerticalLayers.layer;//cfg.get<int>("num_soil_layers");
    moisture_content_solver.emplace(domain, sizes.vert_layers);

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_local_faces(); i++)
    {
        constexpr auto OP = orderedPair{.layer=cellFacesAndVerticalLayers.layer,.face=cellFacesAndVerticalLayers.face};
        auto face = domain->face(i);

        cellInfo<OP> cell_info = cellInfo<OP>::build(face,_params,sizes);

        face->template make_module_data<data>(ID,cell_info,face);
    }

    // TODO set domain parameters
}
template <MeshInterface M>
template <typename F>
    requires std::invocable<F,HashName>
void SoilMoistureSolverCore<M>::depends(F&&)
{
    //depends on nothing
}

template <MeshInterface M>
template <typename F>
    requires std::invocable<F,HashName>
void SoilMoistureSolverCore<M>::provides(F&& f)
{
    for (const auto& names = build_output_names(); const auto& [psi, theta] : names)
    {
        f(psi);
        f(theta);
    }
}
} // namespace SoilMoistureSolver
