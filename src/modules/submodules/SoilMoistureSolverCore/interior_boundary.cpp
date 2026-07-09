//
// Created by Allum, Donovan on 2026-07-09.
//

#include "interior_boundary.hpp"

namespace SoilMoistureSolver::detail
{
Geometry set_top_geo(const Point_3& cell_centre, const Pair depth)
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

Geometry set_bottom_geo(const Point_3& cell_centre,
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
}