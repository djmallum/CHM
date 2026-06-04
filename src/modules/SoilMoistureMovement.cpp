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

#include "Soil.h"
#include "SoilMoistureSolverCore/Details.hpp"
#include <boost/xpressive/detail/utility/traits_utils.hpp>
REGISTER_MODULE_CPP(SoilMoistureMovement)

SoilMoistureMovement::SoilMoistureMovement(const config_file& cfg)
    : module_base("SoilMoistureMovement", parallel::domain, cfg),
    _solver(ID,SoilMoistureSolver::param_builder(cfg,*global_param))
{
    // TODO Add depends/provides
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

auto translate_down(const double centre_depth, const Point_3& original_point)
{
    return Point_3(original_point.x(), original_point.y(), original_point.z() - centre_depth);
}
auto translate_up(const double centre_depth, const Point_3& original_point)
{
    return translate_down(-centre_depth, original_point);
}
SoilMoistureMovement::Geometry SoilMoistureMovement::Interior::set_lateral_geo(const mesh_elem& face, const GeoHelper geo_help,
                                               const Point_3& cell_centre,
                                               const double centre_depth)
{
    Geometry g;
    const auto edge = translate_down(centre_depth,face->edge_midpoint<Point_3>(static_cast<int>(geo_help.nn)));
    g.to_face.owner = CGAL::sqrt(CGAL::squared_distance(cell_centre, edge));

    const auto centre_neighbour = translate_down(centre_depth,face ->neighbor(static_cast<int>(geo_help.nn))->center());

    g.to_face.neighbour = CGAL::sqrt(CGAL::squared_distance(centre_neighbour, edge));

    g.elevation.owner = cell_centre.z();
    g.elevation.neighbour = centre_neighbour.z();

    g.cell_centre_distance = CGAL::sqrt(CGAL::squared_distance(cell_centre, centre_neighbour));
    return g;
}
SoilMoistureMovement::Geometry SoilMoistureMovement::Interior::set_top_geo(const Point_3& cell_centre, const Pair depth)
{
    const auto lower_depth = depth.owner;
    const auto upper_depth = depth.neighbour;

    Geometry g;
    g.to_face.owner = lower_depth / 2.0;
    g.to_face.neighbour = upper_depth / 2.0;

    g.elevation.owner = cell_centre.z();
    g.elevation.neighbour = cell_centre.z() + g.to_face.owner + g.to_face.neighbour;

    g.cell_centre_distance = upper_depth - lower_depth;
    return g;
}
SoilMoistureMovement::Geometry SoilMoistureMovement::Interior::set_bottom_geo(const Point_3& cell_centre,
                                         const Pair depth)
{
    const auto upper_depth = depth.owner;
    const auto lower_depth = depth.neighbour;
    Geometry g;
    g.to_face.owner = upper_depth / 2.0;
    g.to_face.neighbour = lower_depth / 2.0;

    g.cell_centre_distance = g.to_face.owner + g.to_face.neighbour;

    g.elevation.owner = cell_centre.z();
    g.elevation.neighbour = g.elevation.owner - g.cell_centre_distance;
    return g;
}
SoilMoistureMovement::Geometry SoilMoistureMovement::Interior::lower_layer_boundary(const mesh_elem& face,
                                                          const GeoHelper geo_help,
                                                          const Point_3& face_centre, const double lower_depth,
                                                          const double recharge_depth)
{
    const auto centre_depth = lower_depth / 2.0 + recharge_depth;
    const auto cell_centre = translate_down(centre_depth, face_centre);

    switch (static_cast<Neighbour>(geo_help.nn))
    {
    case Neighbour::Lateral_0:
    case Neighbour::Lateral_1:
    case Neighbour::Lateral_2:
    {
        return set_lateral_geo(face, geo_help, cell_centre, centre_depth);
    }
    case Neighbour::Top:
    {
        Pair p;
        p.owner = lower_depth;
        p.neighbour = recharge_depth;
        return set_top_geo(cell_centre, p);
    }
    case Neighbour::Bottom:
    {
        const std::string err = std::format("Not Possible, interior boundaries only for this type");
        CHM_THROW_EXCEPTION(module_error, err);
    }
    }
    return Geometry{};
}
SoilMoistureMovement::Geometry
SoilMoistureMovement::Interior::recharge_layer_boundary(const mesh_elem& face, const GeoHelper geo_help,
                                                        const Point_3& face_centre, const double lower_depth,
                                                        const double recharge_depth, const double detention_depth)
{
    Geometry g;
    const auto centre_depth = recharge_depth / 2.0;
    const auto cell_centre = translate_down(centre_depth, face_centre);

    switch (static_cast<Neighbour>(geo_help.nn))
    {
    case Neighbour::Lateral_0:
    case Neighbour::Lateral_1:
    case Neighbour::Lateral_2:
        g = set_lateral_geo(face, geo_help, cell_centre, centre_depth);
        return g;
    case Neighbour::Top:
    {
        Pair p;
        p.owner = recharge_depth;
        p.neighbour = detention_depth;
        g = set_top_geo(cell_centre, p);
        return g;
    }
    case Neighbour::Bottom:
    {
        Pair p;
        p.owner = recharge_depth;
        p.neighbour = lower_depth;
        g = set_bottom_geo(cell_centre, p);
        return g;
    }
    }
    return g;
}
SoilMoistureMovement::Geometry SoilMoistureMovement::Interior::detention_layer_boundary(
    const mesh_elem& face, const GeoHelper geo_help, const Point_3& face_centre,
    const double recharge_depth, const double detention_depth)
{
    Geometry g;
    const auto centre_depth = detention_depth / 2.0;
    const auto cell_centre = translate_up(centre_depth, face_centre);

    switch (static_cast<Neighbour>(geo_help.nn))
    {
    case Neighbour::Lateral_0:
    case Neighbour::Lateral_1:
    case Neighbour::Lateral_2:
        g = set_lateral_geo(face, geo_help, cell_centre, centre_depth);
        return g;
    case Neighbour::Top:
    {
        const std::string err = std::format("Not Possible, interior boundaries only for this type");
        CHM_THROW_EXCEPTION(module_error, err);
    }
    case Neighbour::Bottom:
    {
        Pair p;
        p.owner = detention_depth;
        p.neighbour = recharge_depth;
        g = set_bottom_geo(cell_centre, p);
        return g;
    }
    }
    return g;
}
SoilMoistureMovement::Geometry SoilMoistureMovement::Interior::set_interior_geometry(const mesh_elem& face, const GeoHelper geo_help, const config_file& cfg)
{
    const auto face_centre = face->center();
    const auto lower_depth = cfg.get<double>("lower_soil_depth");
    const auto recharge_depth = cfg.get<double>("recharge_soil_depth");
    const auto detention_depth = cfg.get<double>("detention_layer_depth");
    switch (geo_help.layer)
    {
    case 0: // lower soil layer
    {
        return lower_layer_boundary(face, geo_help, face_centre, lower_depth, recharge_depth);
    }
    case 1: // recharge depth
    {
        return recharge_layer_boundary(face, geo_help, face_centre, lower_depth, recharge_depth, detention_depth);
    }
    case 2:
    {
        return detention_layer_boundary(face, geo_help, face_centre, recharge_depth, detention_depth);
    }
    default:
    {
        const std::string err = "There are only 3 vertical layers, exceeded in geometry setup";
        CHM_THROW_EXCEPTION(module_error,err);
    }
    }
}
SoilMoistureMovement::Interior::Interior(const mesh_elem& face, const GeoHelper gh, const config_file& cfg) : geometry(set_interior_geometry(face,gh,cfg))
{

}
SoilMoistureMovement::Boundary::Boundary(const double dz, const double df) : DeltaZ(dz), distance_to_face(df) {};
SoilMoistureMovement::data::faceType SoilMoistureMovement::get_lateral_boundary(const mesh_elem& face,
                                                                                const size_t nn) const
{
    const auto dir_vec = face->downslope_dir();
    const auto edge_normal = face->edge_unit_normal<Vector_3>(nn);
    const auto face_centre = face->center();

    if (CGAL::scalar_product(dir_vec, edge_normal) >= 0.0)
    {
        // Downslope out the boundary
        const auto dist_to_face =
            CGAL::sqrt(CGAL::squared_distance(face_centre, face->edge_midpoint<Point_3>(static_cast<int>(nn))));
        const auto DeltaZ = -2 * dist_to_face / std::tan(face->slope());
        return Boundary(DeltaZ, dist_to_face);
    }
    // Upslope out the boundary
    const auto dist_to_face =
        CGAL::sqrt(CGAL::squared_distance(face_centre, face->edge_midpoint<Point_3>(static_cast<int>(nn))));
    const auto DeltaZ = 2 * dist_to_face / std::tan(face->slope());

    return Boundary(DeltaZ, dist_to_face);
}
SoilMoistureMovement::data::faceType
SoilMoistureMovement::get_geometry(const mesh_elem& face, const GeoHelper geo_helper)
{

    // TODO handle faces at the edge of their process,
    // currently just leaving the std::optional in neighbours uninitialized
    if (const auto nn = geo_helper.nn; face->neighbor(nn) == nullptr)
    {
        // boundary conditions
        switch (static_cast<Neighbour>(nn))
        {
        case Neighbour::Lateral_0:
        case Neighbour::Lateral_1:
        case Neighbour::Lateral_2:
        {
            return get_lateral_boundary(face, nn);
        }
        break;
        case Neighbour::Top:
        {
            const auto detention_depth = cfg.get<double>("detention_soil_depth");
            return Boundary(0.0,detention_depth / 2.0);
        }
        case Neighbour::Bottom:
        {
            const auto lower_depth = cfg.get<double>("lower_soil_depth");
            const auto dist_to_face = lower_depth / 2.0;
            const auto DeltaZ = 2.0 * dist_to_face;
            return Boundary(DeltaZ,dist_to_face);
        }
        }
    }
    else
    {
        return Interior(face, geo_helper,cfg);
    }
}
double SoilMoistureMovement::data::theta(const int i) const
{

    const auto result =std::pow(psi_n.at(i)/ air_entry_tension,-1/pore_size_dist_index);

    return result;
}
// --- helpers to fill array without default-constructing faceType ---
template <size_t NumNeighbours, size_t... Is>
std::array<SoilMoistureMovement::data::faceType, NumNeighbours> SoilMoistureMovement::build_layers(mesh_elem& face, const size_t i, std::index_sequence<Is...>) {
    std::array<std::optional<data::faceType>, NumNeighbours> scratch{};
    for (size_t nn = 0; nn < NumNeighbours; ++nn)
    {
        GeoHelper geo_helper;
        geo_helper.nn = nn;
        geo_helper.layer = i;
        scratch[i].emplace(get_geometry(face, geo_helper));
    }
    return { *std::move(scratch[Is])... };
}

template <size_t NumNeighbours, size_t NumLayers, size_t... Ns>
std::array<std::array<SoilMoistureMovement::data::faceType, NumNeighbours>, NumLayers>
SoilMoistureMovement::build_geometry(mesh_elem& face,std::index_sequence<Ns...>) {
    return { build_layers<NumNeighbours>(face, Ns, std::make_index_sequence<NumNeighbours>{}) ... };
}
void SoilMoistureMovement::init(mesh& domain)
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
        mesh_elem face = domain->face(i);
        auto geometry = build_geometry<NUM_NEIGHBOURS,NUM_LAYERS>(face,
                std::make_index_sequence<NUM_LAYERS>{});

        data::LayerNeighbourArray<data::opt<faceInterpolator>> face_interp;
        data::LayerNeighbourArray<double> alpha{};
        data::LayerNeighbourArray<double> face_area{};
        data::LayerNeighbourArray<data::opt<int>> neighbour_idx{};

        static_assert(geometry.size() == NUM_LAYERS && geometry[0].size() == NUM_NEIGHBOURS,
            "geometry matrix built improperly");

        std::array<double,NUM_LAYERS> volume;
        const double tri_area = face->get_area();
        const std::array depths{
            cfg.get<double>("lower_soil_depth"),
            cfg.get<double>("recharge_soil_depth"),
            cfg.get<double>("detention_depth")
        };

        for (size_t layer = 0; layer < NUM_LAYERS; layer++)
        {
            volume.at(layer) = tri_area * depths.at(layer);
        }

        for (size_t nn = 0; nn < NUM_NEIGHBOURS; ++nn)
        {
            for (size_t layer = 0; layer < NUM_LAYERS; ++layer)
            {
                std::visit([&face_interp,nn,layer]<typename T0>(T0&& arg)
                {
                    using T = std::decay_t<T0>;
                    if constexpr (std::is_same_v<T, Interior>)
                    {
                        face_interp[layer][nn].emplace(arg.geometry);
                    }

                }, geometry[layer][nn]);


                switch (static_cast<Neighbour>(nn))
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

                alpha[layer][nn] = global_param->dt() * face_area[layer][nn] / volume[layer];
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

        face->make_module_data<data>(ID,face_interp, geometry, alpha, neighbour_idx, cfg);


    }

    // TODO set domain parameters
}
SoilMoistureMovement::faceInterpolator::faceInterpolator(const Geometry& g)
    : geometry_(g) {}

double SoilMoistureMovement::faceInterpolator::interp(const Pair p) const
{
    return geometry_.cell_centre_distance * p.owner * p.neighbour /
        (geometry_.to_face.owner * p.neighbour + geometry_.to_face.neighbour * p.owner);
}
SoilMoistureMovement::Pair SoilMoistureMovement::faceInterpolator::elevation() const
{
    return geometry_.elevation;
}
SoilMoistureMovement::data::data(const LayerNeighbourArray<opt<faceInterpolator>>& interp,
                                 const LayerNeighbourArray<faceType>& cell,
                                 const LayerNeighbourArray<double>& alpha,
                                 const LayerNeighbourArray<opt<int>>& neighbour_idx,
                                 const config_file& cfg)
: interp_to_face(interp), cell_geometry(cell), alpha(alpha), neighbour_idx_(neighbour_idx)
{
    const auto& soil_param = Soil::get_soil_obj<Soil::soils_na>();
    const auto soil_type = cfg.get<std::string>("soil_type"); // Should be triangle specific in future versions
    pore_size_dist_index = soil_param.pore_size_dist(soil_type);
    air_entry_tension = soil_param.air_entry_tension(soil_type);
    K_saturated = soil_param.saturated_conductivity(soil_type);

    const auto initial_psi = cfg.get<double>("initial_pressure_head");
    // TODO this should be in theta, converted to psi, with a block for 0 saturation
    std::ranges::fill(psi_n, initial_psi);

}
// SoilMoistureMovement::Pair SoilMoistureMovement::data::get_elevation(const int f) const
// {
//     return interpolator[f].elevation();
// }
