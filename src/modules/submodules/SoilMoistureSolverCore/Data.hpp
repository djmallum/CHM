//
// Created by Allum, Donovan on 2026-06-03.
//

#pragma once
#include "triangulation.hpp"
#include "Details.hpp"
#include "Soil.h"

namespace SoilMoistureSolver::detail
{

template<typename T>
using opt = std::optional<T>;
template<typename T>
using LayerNeighbourArray = std::array<std::array<T,NUM_NEIGHBOURS>,NUM_LAYERS>;
using faceType = std::variant<Boundary,Interior>;

struct cellInfo
{
    const LayerNeighbourArray<opt<faceInterpolator>> interp_to_face;
    const LayerNeighbourArray<faceType> cell_geometry;
    const LayerNeighbourArray<double> alpha;
    const LayerNeighbourArray<opt<int>> neighbour_idx_;
};

class data : public face_info
{
public:
    template<ElementInterface E>
    explicit data(const cellInfo& cell,
        const E& face);

    cellInfo cell_info;
    double K_saturated;
    double air_entry_tension;
    double pore_size_dist_index;
    double theta(int) const;
    std::array<double,NUM_LAYERS> psi_n{0.0};
};

inline double data::theta(const int i) const
{
    const auto result =std::pow(psi_n.at(i)/ air_entry_tension,-1/pore_size_dist_index);

    return result;
}

template <ElementInterface E>
data::data(const cellInfo& cell,
    const E& face)
: cell_info(cell)
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
        return Boundary(DeltaZ, dist_to_face);
    }
    // Upslope out the boundary
    const auto dist_to_face =
        CGAL::sqrt(CGAL::squared_distance(face_centre, face->template edge_midpoint<Point_3>(static_cast<int>(nn))));
    const auto DeltaZ = 2 * dist_to_face / std::tan(face->slope());

    return Boundary(DeltaZ, dist_to_face);
}

template<ElementInterface E>
static faceType get_geometry(const E& face, const Params& p ,const GeoHelper geo_helper)
{

    // TODO handle faces at the edge of their process,
    // currently just leaving the std::optional in neighbours uninitialized

    if (const size_t nn = static_cast<const size_t>(geo_helper.nn); geo_helper.top_or_bottom_boundary()
        || face->neighbor(nn) == nullptr)
    {
        // boundary conditions
        switch (geo_helper.nn)
        {
        case Neighbour::Lateral_0:
        case Neighbour::Lateral_1:
        case Neighbour::Lateral_2:
        {
            return get_lateral_boundary(face, nn);
        }
        case Neighbour::Top:
        {
            const auto detention_depth = p.detention_depth;
            return Boundary(0.0,detention_depth / 2.0);
        }
        case Neighbour::Bottom:
        {
            const auto lower_depth = p.lower_soil_depth;
            const auto dist_to_face = lower_depth / 2.0;
            const auto DeltaZ = 2.0 * dist_to_face;
            return Boundary(DeltaZ,dist_to_face);
        }
        }
    }

    return Interior(face, geo_helper, p);
}

template <ElementInterface E, size_t NumNeighbours, size_t... Is>
std::array<faceType, NumNeighbours> build_layers(E& face, const Params& p, const size_t i, std::index_sequence<Is...>) {
    std::array<std::optional<faceType>, NumNeighbours> scratch{};
    for (const auto neighbour : all_neighbours)
    {
        GeoHelper geo_helper;
        geo_helper.nn = neighbour;
        geo_helper.layer = i;
        scratch[i].emplace(get_geometry(face, p, geo_helper));
    }
    return { *std::move(scratch[Is])... };
}


template <ElementInterface E, size_t NumNeighbours, size_t NumLayers, size_t... Ns>
std::array<std::array<faceType, NumNeighbours>, NumLayers>
build_geometry(E& face, const Params& p, std::index_sequence<Ns...>) {
    return { build_layers<E,NumNeighbours>(face, p, Ns, std::make_index_sequence<NumNeighbours>{}) ... };
}
}
