
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

constexpr bool operator==(Neighbour lhs, const size_t rhs) {
    return static_cast<int>(lhs) == rhs;
}
constexpr bool operator==(const size_t lhs,const Neighbour rhs) {
    return rhs == lhs;
}

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

    [[nodiscard]] double interp(Pair) const;
    [[nodiscard]] Pair elevation() const;

  private:
    const Geometry geometry_;
};

struct Boundary
{
    const double DeltaZ;
    const double distance_to_face;
    const Neighbour neighbour_type;
    explicit Boundary(double dz, double df,Neighbour type);
};


inline auto translate_down(const double centre_depth, const Point_3& original_point)
{
    return Point_3(original_point.x(), original_point.y(), original_point.z() - centre_depth);
}
inline auto translate_up(const double centre_depth, const Point_3& original_point)
{
    return translate_down(-centre_depth, original_point);
}

} // namespace SoilMoistureSolverCore::detail
