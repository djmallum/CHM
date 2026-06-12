#pragma once

#include <utility>

#include "triangulation.hpp"
#include "Details.hpp"
#include "Soil.h"


namespace SoilMoistureSolver::detail
{
template<typename T>
using opt = std::optional<T>;
template<typename T,orderedPair P = cellFacesAndVerticalLayers>
using LayerNeighbourArray = std::array<std::array<T,P.face>,P.layer>;
using faceType = std::variant<Boundary,Interior>;

template<orderedPair P>
struct cellInfo
{
    template<class E> requires ElementInterface<E>
    static cellInfo build(E& face,const Params&,Sizes);

    [[nodiscard]] double get_coefficient(const orderedPair&) const;
    template <class F> [[nodiscard]] bool is_faceType(const orderedPair& p) const;
    template <class F> [[nodiscard]] const F& get_faceType(const orderedPair& p) const;
    [[nodiscard]] double value_at_face(const orderedPair&, const Pair& p) const;
    [[nodiscard]] std::optional<int> get_neighbour_idx(const orderedPair&) const;

private:

    const LayerNeighbourArray<opt<faceInterpolator>,P> interp_to_face;
    const LayerNeighbourArray<faceType,P> cell_geometry;
    const LayerNeighbourArray<double,P> coefficient;
    const LayerNeighbourArray<opt<int>,P> neighbour_idx_;

    cellInfo(
    const LayerNeighbourArray<opt<faceInterpolator>>& interp,
    const LayerNeighbourArray<faceType>& geom,
    const LayerNeighbourArray<double>& a,
    const LayerNeighbourArray<opt<int>>& neighbours
) : interp_to_face(interp),
    cell_geometry(geom),
    coefficient(a),
    neighbour_idx_(neighbours)
    {}
};

class data : public face_info
{
public:
    template<ElementInterface E>
    explicit data(cellInfo<cellFacesAndVerticalLayers>  cell,
        const E& face);

    cellInfo<cellFacesAndVerticalLayers> cell_info;
    double K_saturated;
    double air_entry_tension;
    double pore_size_dist_index;
    [[nodiscard]] double theta(int) const;
    std::array<double,cellFacesAndVerticalLayers.layer> psi_n{0.0};
};

inline double data::theta(const int i) const
{
    const auto result =std::pow(psi_n.at(i)/ air_entry_tension,-1/pore_size_dist_index);

    return result;
}

template <ElementInterface E>
data::data(cellInfo<cellFacesAndVerticalLayers> cell,
    const E& face)
: cell_info(std::move(cell))
{
    const auto& soil_param = Soil::get_soil_obj<Soil::soils_na>();
    const auto soil_type = face->template soil_attribute<std::string>("soil_type");
    pore_size_dist_index = soil_param.pore_size_dist(soil_type);
    air_entry_tension = soil_param.air_entry_tension(soil_type);
    K_saturated = soil_param.saturated_conductivity(soil_type);

    const auto initial_psi = face->template soil_attribute<double>("initial_pressure_head");
    // TODO this should be in theta, converted to psi, with a block for 0 saturation
    std::ranges::fill(psi_n, initial_psi);
}

template<ElementInterface E>
static faceType get_lateral_boundary(const E& face,const size_t nn)
{
    const auto dir_vec = face->downslope_dir();
    const auto edge_normal = face->template edge_unit_normal<Vector_3>(nn);
    const auto face_centre = face->center();

    if (CGAL::scalar_product(dir_vec, edge_normal) >= 0.0)
    {
        // Downslope out the boundary
        const auto dist_to_face =
            CGAL::sqrt(CGAL::squared_distance(face_centre, face->template edge_midpoint<Point_3>(static_cast<int>(nn))));
        const auto DeltaZ = -2 * dist_to_face / std::tan(face->slope());
        return Boundary(DeltaZ, dist_to_face, face_index_to_Neighbour(nn));
    }
    // Upslope out the boundary
    const auto dist_to_face =
        CGAL::sqrt(CGAL::squared_distance(face_centre, face->template edge_midpoint<Point_3>(static_cast<int>(nn))));
    const auto DeltaZ = 2 * dist_to_face / std::tan(face->slope());

    return Boundary(DeltaZ, dist_to_face, face_index_to_Neighbour(nn));
}

template<ElementInterface E>
static faceType get_geometry(const E& face, const Params& p ,const orderedPair& op)
{

    // TODO handle faces at the edge of their process,
    // currently just leaving the std::optional in neighbours uninitialized

    if (op.top_or_bottom_boundary()
        || face->neighbor(op.face) == nullptr)
    {
        // boundary conditions
        switch (face_index_to_Neighbour(op))
        {
        case Neighbour::Lateral_0:
        case Neighbour::Lateral_1:
        case Neighbour::Lateral_2:
        {
            return get_lateral_boundary(face, op.face);
        }
        case Neighbour::Top:
        {
            const auto detention_depth = p.detention_depth;
            return Boundary(0.0,detention_depth / 2.0,face_index_to_Neighbour(op));
        }
        case Neighbour::Bottom:
        {
            const auto lower_depth = p.lower_soil_depth;
            const auto dist_to_face = lower_depth / 2.0;
            const auto DeltaZ = 2.0 * dist_to_face;
            return Boundary(DeltaZ,dist_to_face,face_index_to_Neighbour(op));
        }
        }
    }

    return Interior(face, op, p);
}

template <ElementInterface E, size_t NumNeighbours, size_t... Is>
std::array<faceType, NumNeighbours> build_layers(E& face, const Params& p, const size_t i, std::index_sequence<Is...>) {
    std::array<std::optional<faceType>, NumNeighbours> scratch{};
    for (const auto neighbour : all_neighbours)
    {
        orderedPair op{.layer =i,.face = static_cast<size_t>(Neighbour_to_face_index(neighbour))};
        scratch[i].emplace(get_geometry(face, p, op));
    }
    return { *std::move(scratch[Is])... };
}


template <ElementInterface E, orderedPair P = cellFacesAndVerticalLayers, size_t... Ns>
std::array<std::array<faceType, P.face>, P.layer>
build_geometry(E& face, const Params& p, std::index_sequence<Ns...>) {
    return { build_layers<E,P.face>(face, p, Ns, std::make_index_sequence<P.layer>{}) ... };
}

template<orderedPair P>
template<class E>
    requires ElementInterface<E>
cellInfo<P> cellInfo<P>::build(E& face,const Params& _params,const Sizes sizes)
{
    auto geometry =
        build_geometry<E, P>(face, _params, std::make_index_sequence<P.layer>{});

    LayerNeighbourArray<opt<faceInterpolator>,P> face_interp;
    LayerNeighbourArray<double,P> alpha{};
    LayerNeighbourArray<double,P> face_area{};
    LayerNeighbourArray<opt<int>,P> neighbour_idx{};

    static_assert(geometry.size() == P.layer && geometry[0].size() == P.face,
                  "geometry matrix built improperly");

    std::array<double, cellFacesAndVerticalLayers.layer> volume;
    const double tri_area = face->get_area();
    const std::array depths{_params.lower_soil_depth, _params.recharge_soil_depth, _params.detention_depth};

    for (size_t layer = 0; layer < cellFacesAndVerticalLayers.layer; layer++)
    {
        volume.at(layer) = tri_area * depths.at(layer);
    }

    for (const auto neighbour : all_neighbours) // size_t nn = 0; nn < cellFacesAndVerticalLayers.face; ++nn
    {
        for (size_t layer = 0; layer < cellFacesAndVerticalLayers.layer; ++layer)
        {
            const auto nn = static_cast<size_t>(neighbour);
            std::visit(
                [&face_interp, nn, layer]<typename T>(T&& arg)
                {
                    if constexpr (std::is_same_v<T, Interior>)
                    {
                        face_interp[layer][nn].emplace(arg.geometry);
                    }
                },
                geometry[layer][nn]);

            switch (neighbour)
            {
            case Neighbour::Lateral_0:
            case Neighbour::Lateral_1:
            case Neighbour::Lateral_2:
            {
                const auto side_length = face->edge_length(nn);
                face_area[layer][nn] = side_length * depths[layer];
                // neighbour
                neighbour_idx[layer][nn] = std::visit(
                    [layer, nn,sizes, &face]<typename T>(T&&) -> std::optional<int>
                    {
                        if constexpr (std::is_same_v<T, Interior>)
                        {
                            return std::optional<int>(sizes.global * layer + face->neighbor(nn)->cell_global_id);
                        }
                        else
                            return std::nullopt;
                    },
                    geometry[layer][nn]);
                break;
            }
            case Neighbour::Top:
                face_area[layer][nn] = tri_area;
                // neighbour
                neighbour_idx[layer][nn] = std::visit(
                    [layer,sizes, &face]<typename T>(T&&) -> std::optional<int>
                    {
                        if constexpr (std::is_same_v<T, Interior>)
                        {
                            if (layer == cellFacesAndVerticalLayers.layer - 1)
                            {
                                const std::string err =
                                    std::format("geometry object is the Interior faceType variant but"
                                                "layer was allowed to be at the top: {}."
                                                "Should be lower than {}",
                                                layer, cellFacesAndVerticalLayers.layer - 1);
                                CHM_THROW_EXCEPTION(module_error, err);
                            }

                            return std::optional<int>(sizes.global * (layer + 1) + face->cell_global_id);
                        }
                        else
                            return std::nullopt;
                    },
                    geometry[layer][nn]);
                break;
            case Neighbour::Bottom:
                face_area[layer][nn] = tri_area;
                // neighbour
                neighbour_idx[layer][nn] = std::visit(
                    [layer,sizes, &face]<typename T>(T&&) -> std::optional<int>
                    {
                        if constexpr (std::is_same_v<T, Interior>)
                        {
                            if (layer == 0)
                            {
                                const std::string err =
                                    std::format("geometry object is the Interior faceType variant but"
                                                "layer was allowed to be at the bottom: {}."
                                                "Should be higher than 0.",
                                                layer);
                                CHM_THROW_EXCEPTION(module_error, err);
                            }
                            return std::optional<int>(sizes.global * (layer - 1) + face->cell_global_id);
                        }
                        else
                            return std::nullopt;
                    },
                    geometry[layer][nn]);
                break;
            }

            alpha[layer][nn] = _params.time_step_seconds * face_area[layer][nn] / volume[layer];
            alpha[layer][nn] *= std::visit(
                []<typename T>(T&& arg) -> double
                {
                    if constexpr (std::is_same_v<T, Interior>)
                    {
                        return arg.geometry.cell_centre_distance;
                    }
                    else if constexpr (std::is_same_v<T, Boundary>)
                    {
                        // There exists a mythical neighbour past the boundary
                        return 2.0 * arg.distance_to_face;
                    }

                    CHM_THROW_EXCEPTION(module_error, "Visiting to invalid type, should be unreachable");
                },
                geometry[layer][nn]);
        }
    }
    return cellInfo{face_interp, geometry, alpha, neighbour_idx};
}
template <orderedPair P> double cellInfo<P>::get_coefficient(const orderedPair& p) const
{
    return coefficient[p.layer][p.face];
}
template <orderedPair P>
template <typename F>
bool cellInfo<P>::is_faceType(const orderedPair& p) const
{
    return std::holds_alternative<F>(cell_geometry[p.layer][p.face]);
}
template <orderedPair P> template<class F> const F& cellInfo<P>::get_faceType(const orderedPair& p) const
{
    return std::get<F>(cell_geometry[p.layer][p.face]);
}
template <orderedPair P> double cellInfo<P>::value_at_face(const orderedPair& op, const Pair& p) const
{
    return interp_to_face.at(op.layer).at(op.face)->interp(p);
}
template <orderedPair P> std::optional<int> cellInfo<P>::get_neighbour_idx(const orderedPair& p) const
{
    return neighbour_idx_.at(p.layer).at(p.face);
}
} // namespace SoilMoistureSolver::detail
