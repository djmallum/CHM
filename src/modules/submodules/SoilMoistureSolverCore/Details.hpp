
#pragma once

#include "Concepts.hpp"

namespace SoilMoistureSolver::detail
{
struct orderedPair {
    const size_t layer{};
    const size_t face{};

    auto operator<=>(const orderedPair&) const = default;
    [[nodiscard]] constexpr bool top_or_bottom_boundary() const;
};

static inline constexpr size_t NUM_SIDES = 3; // Because triangle
constexpr auto cellFacesAndVerticalLayers = orderedPair{.layer = 3,.face = 5};


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
    Lateral_0, // switch cases will receive a number from 0-4, even if lateral cases are the same
    Lateral_1,
    Lateral_2,
    Top,
    Bottom
};

constexpr Neighbour face_index_to_Neighbour(const orderedPair& op)
{
    switch (op.face)
    {
    case 0:
        return Neighbour::Lateral_0;
    case 1:
        return Neighbour::Lateral_1;
    case 2:
        return Neighbour::Lateral_2;
    case 3:
        return Neighbour::Bottom;
    case 4:
        return Neighbour::Top;
    default:
        CHM_THROW_EXCEPTION(module_error,"Index out of bounds");
    }
}

constexpr Neighbour face_index_to_Neighbour(const size_t nn)
{
    const orderedPair op{.layer = 0,.face = nn};
    return face_index_to_Neighbour(op);
}

constexpr size_t Neighbour_to_face_index(const Neighbour N)
{
    switch (N) {
    case Neighbour::Lateral_0:
        return 0u;
    case Neighbour::Lateral_1:
        return 1u;
    case Neighbour::Lateral_2:
        return 2u;
    case Neighbour::Bottom:
        return 3u;
    case Neighbour::Top:
        return 4u;
    default:
        CHM_THROW_EXCEPTION(module_error,"Neighbour out of bounds");
    }
}

constexpr bool orderedPair::top_or_bottom_boundary() const
{
    return (layer == 0 && face_index_to_Neighbour(*this) == Neighbour::Bottom) ||
        (layer == cellFacesAndVerticalLayers.layer - 1 && face_index_to_Neighbour(*this) == Neighbour::Top);
}

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

static_assert(all_neighbours.size() == cellFacesAndVerticalLayers.face,"Length of all_neighbours array does not match dim_size.face");

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
    const Neighbour neighbour_type;
    explicit Boundary(double dz, double df,Neighbour type);
};

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

inline auto translate_down(const double centre_depth, const Point_3& original_point)
{
    return Point_3(original_point.x(), original_point.y(), original_point.z() - centre_depth);
}
inline auto translate_up(const double centre_depth, const Point_3& original_point)
{
    return translate_down(-centre_depth, original_point);
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
} // namespace SoilMoistureSolverCore::detail
