//
// Created by Allum, Donovan on 2026-06-03.
//

#pragma once

#include "Concepts.hpp"

namespace SoilMoistureSolver::detail
{
static inline constexpr size_t NUM_LAYERS = 3;
static inline constexpr size_t NUM_NEIGHBOURS = 5;

struct Params
{
    double lower_soil_depth;
    double recharge_soil_depth;
    double detention_depth;
    size_t time_step_seconds;
};

struct Sizes
{
    size_t local;
    size_t global;
    size_t vert_layers = 1;
};

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

enum class Neighbour
{
    // WARNING
    // DO NOT CHANGE THESE VALUES
    // Array indexing relies on this
    Lateral_0 = 0, // switch cases will receive a number from 0-4, even if lateral cases are the same
    Lateral_1 = 1,
    Lateral_2 = 2,
    Top = 3,
    Bottom = 4
};

// Do this to double-check Neighbour will cast correctly
static_assert(Neighbour::Bottom == static_cast<Neighbour>(4) &&
    static_cast<int>(Neighbour::Bottom) == 4,
    "Enum Class Neighbour, defined for clarity, not casting to the correct values");

constexpr std::array all_neighbours = {
    Neighbour::Lateral_0,
    Neighbour::Lateral_1,
    Neighbour::Lateral_2,
    Neighbour::Top,
    Neighbour::Bottom
};

static_assert(all_neighbours.size() == NUM_NEIGHBOURS,"Length of all_neighbours array does not match NUM_NEIGHBOURS");
struct GeoHelper
{
    Neighbour nn;
    size_t layer;

    bool top_or_bottom_boundary() const
    {
        switch (nn)
        {
        case Neighbour::Top:
            if (layer == NUM_LAYERS) return true;
            return false;
        case Neighbour::Bottom:
            if (layer == 0) return true;
            return false;
        default:
            return false;
        }
    }
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

class Boundary
{
  public:
    const double DeltaZ;
    const double distance_to_face;
    explicit Boundary(double dz, double df);
};
class Interior
{
  public:
    explicit Interior(const mesh_elem&, GeoHelper, const Params&);
    const Geometry geometry;
    faceInterpolator interpolator{geometry};

  private:
    // Static Helpers
    static Geometry set_interior_geometry(const mesh_elem& face, GeoHelper geo_help, const Params& cfg);
    static Geometry set_bottom_geo(const Point_3& cell_centre, Pair depth);
    static Geometry lower_layer_boundary(const mesh_elem& face, GeoHelper geo_help, const Point_3& face_centre,
                                         double lower_depth, double recharge_depth);
    static Geometry recharge_layer_boundary(const mesh_elem& face, GeoHelper geo_help, const Point_3& face_centre,
                                            double lower_depth, double recharge_depth, double detention_depth);
    static Geometry detention_layer_boundary(const mesh_elem& face, GeoHelper geo_help, const Point_3& face_centre,
                                             double recharge_depth, double detention_depth);
    static Geometry set_top_geo(const Point_3& cell_centre, Pair depth);
    static Geometry set_lateral_geo(const mesh_elem& face, GeoHelper geo_help, const Point_3& cell_centre,
                                    double centre_depth);
};

} // namespace SoilMoistureSolverCore::detail