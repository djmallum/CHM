#include "Details.hpp"
#include "triangulation.hpp"

namespace SoilMoistureSolver::detail
{

faceInterpolator::faceInterpolator(const Geometry& g)
    : geometry_(g) {}

double faceInterpolator::interp(const Pair p) const
{
    return geometry_.cell_centre_distance * p.owner * p.neighbour /
        (geometry_.to_face.owner * p.neighbour + geometry_.to_face.neighbour * p.owner);
}
Pair faceInterpolator::elevation() const
{
    return geometry_.elevation;
}
Boundary::Boundary(const double dz, const double df, const Neighbour type) : DeltaZ(dz), distance_to_face(df), neighbour_type(type) {}
} // namespace SoilMoistureSolver::detail
