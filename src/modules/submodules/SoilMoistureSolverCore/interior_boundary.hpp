#pragma once

#include "Details.hpp"
#include "Concepts.hpp"

namespace SoilMoistureSolver::detail
{
class Interior
{
public:
    template<ElementInterface E>
    explicit Interior(const E&, const orderedPair&, const Params&);
    const Geometry geometry;
    faceInterpolator interpolator{geometry};
};


template<ElementInterface E>
Interior::Interior(const E& face, const orderedPair& op, const Params& p) : geometry(set_interior_geometry(face,op,p))
{

}

template<ElementInterface E>
static Geometry set_interior_geometry(const E& face, const orderedPair& op, const Params& p )
{
    const auto face_centre = face->center();
    //const auto lower_depth = cfg.get<double>("lower_soil_depth");
    //const auto recharge_depth = cfg.get<double>("recharge_soil_depth");
    //const auto detention_depth = cfg.get<double>("detention_layer_depth");
    switch (op.layer)
    {
    case 0: // lower soil layer
    {
        return lower_layer_boundary(face, op, face_centre, p.lower_soil_depth, p.recharge_soil_depth);
    }
    case 1: // recharge depth
    {
        return recharge_layer_boundary(face, op, face_centre, p.lower_soil_depth, p.recharge_soil_depth, p.detention_depth);
    }
    case 2:
    {
        return detention_layer_boundary(face, op, face_centre, p.recharge_soil_depth, p.detention_depth);
    }
    default:
    {
        const std::string err = "There are only 3 vertical layers, exceeded in geometry setup";
        CHM_THROW_EXCEPTION(module_error,err);
    }
    }
}

Geometry set_top_geo(const Point_3& cell_centre, Pair depth);
Geometry set_bottom_geo(const Point_3& cell_centre, Pair depth);

template<ElementInterface E>
static Geometry detention_layer_boundary(
    const E& face, const orderedPair& op, const Point_3& face_centre,
    const double recharge_depth, const double detention_depth)
{
    const auto centre_depth = detention_depth / 2.0;
    const auto cell_centre = translate_up(centre_depth, face_centre);

    switch (static_cast<Neighbour>(op.face))
    {
    case Neighbour::Lateral_0:
    case Neighbour::Lateral_1:
    case Neighbour::Lateral_2:
        return set_lateral_geo(face, op, cell_centre, centre_depth);
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
        return set_bottom_geo(cell_centre, p);
    }
    default:
    {
        const std::string err = "There are only 3 vertical layers, exceeded in geometry setup";
        CHM_THROW_EXCEPTION(module_error,err);
    }
    }
}

template<ElementInterface E> static Geometry
recharge_layer_boundary(const E& face, const orderedPair& op,
                                                        const Point_3& face_centre, const double lower_depth,
                                                        const double recharge_depth, const double detention_depth)
{
    const auto centre_depth = recharge_depth / 2.0;
    const auto cell_centre = translate_down(centre_depth, face_centre);

    switch (static_cast<Neighbour>(op.face))
    {
    case Neighbour::Lateral_0:
    case Neighbour::Lateral_1:
    case Neighbour::Lateral_2:
        return set_lateral_geo(face, op, cell_centre, centre_depth);
    case Neighbour::Top:
    {
        Pair p;
        p.owner = recharge_depth;
        p.neighbour = detention_depth;
        return set_top_geo(cell_centre, p);
    }
    case Neighbour::Bottom:
    {
        Pair p;
        p.owner = recharge_depth;
        p.neighbour = lower_depth;
        return set_bottom_geo(cell_centre, p);
    }
    default:
    {
        const std::string err = "Unreachable path, reached."
            "This path should not be possible but added to silence compiler warning";
        CHM_THROW_EXCEPTION(module_error,err);
    }
    }
}

template<ElementInterface E>
static Geometry lower_layer_boundary(const E& face,
                                                          const orderedPair& op,
                                                          const Point_3& face_centre, const double lower_depth,
                                                          const double recharge_depth)
{
    const auto centre_depth = lower_depth / 2.0 + recharge_depth;
    const auto cell_centre = translate_down(centre_depth, face_centre);

    switch (static_cast<Neighbour>(op.face))
    {
    case Neighbour::Lateral_0:
    case Neighbour::Lateral_1:
    case Neighbour::Lateral_2:
    {
        return set_lateral_geo(face, op, cell_centre, centre_depth);
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
    default:
    {
        const std::string err = "Unreachable path, reached."
            "This path should not be possible but added to silence compiler warning";
        CHM_THROW_EXCEPTION(module_error,err);
    }
    }
}

template<ElementInterface E>
static Geometry set_lateral_geo(const E& face, const orderedPair& op,
                                               const Point_3& cell_centre,
                                               const double centre_depth)
{
    Geometry g;
    const auto edge = translate_down(centre_depth,face->template edge_midpoint<Point_3>(static_cast<int>(op.face)));
    g.to_face.owner = CGAL::sqrt(CGAL::squared_distance(cell_centre, edge));

    const auto centre_neighbour = translate_down(centre_depth,face->neighbor(static_cast<int>(op.face))->center());

    g.to_face.neighbour = CGAL::sqrt(CGAL::squared_distance(centre_neighbour, edge));

    g.elevation.owner = cell_centre.z();
    g.elevation.neighbour = centre_neighbour.z();

    g.cell_centre_distance = CGAL::sqrt(CGAL::squared_distance(cell_centre, centre_neighbour));
    return g;
}
}
