#include "Details.hpp"
#include "gtest/gtest.h"

using namespace SoilMoistureSolver::detail;
class TestFaceInterpolator : public ::testing::Test {
protected:
    Geometry geometry  {.to_face = Pair{.owner = 1.0,.neighbour = 2.5},
                        .elevation = Pair{.owner = 1500.0,.neighbour = 1550.0},
                        .cell_centre_distance = 3.0}; //to_face and cell_centre_distance are NOT invariants, they don't have to be equal

    faceInterpolator face_interpolator{geometry};
};

TEST_F(TestFaceInterpolator, ElevationGetter)
{
    const auto [owner, neighbour] = face_interpolator.elevation();

    EXPECT_DOUBLE_EQ(owner,geometry.elevation.owner);
    EXPECT_DOUBLE_EQ(neighbour,geometry.elevation.neighbour);
}

TEST_F(TestFaceInterpolator,ResultDoesNotNeedToBeBoundedByInputs)
{
    Geometry geometry  {.to_face = Pair{.owner = 1.0,.neighbour = 2.5},
                        .elevation = Pair{.owner = 1500.0,.neighbour = 1550.0},
                        .cell_centre_distance = 3.0}; //to_face and cell_centre_distance are NOT invariants, they don't have to be equal

    faceInterpolator face_interpolator{geometry};

    constexpr auto p_greater = Pair{.owner = 100.0,.neighbour = 60.0};
    ASSERT_GT(geometry.to_face.owner * p_greater.neighbour / p_greater.owner + geometry.to_face.neighbour,geometry.cell_centre_distance)
        << "This assertion is necessary for this test to be valid"
        << "\nBecause the grid is non-orthogonal, it is not the case that the sum of the distances to face"
        << "\nmust be equal to the cell centred distance, in face it is often smaller"
        << "\nWhen this assertion holds, we expected that the interpolated value will be less than p.neighbour"
        << "\nThis can be shown by setting the equation to be less than p.neighbour and rearranging for geometry.cell_centre_distance";
    constexpr auto p_equal = Pair{.owner = 100.0,.neighbour = 50.0};
    ASSERT_DOUBLE_EQ(geometry.to_face.owner * p_equal.neighbour / p_equal.owner + geometry.to_face.neighbour,geometry.cell_centre_distance)
        << "This assertion is necessary for this test to be valid"
        << "\nBecause the grid is non-orthogonal, it is not the case that the sum of the distances to face"
        << "\nmust be equal to the cell centred distance, in face it is often smaller"
        << "\nWhen this assertion holds, we expected that the interpolated value will be equal to p.neighbour"
        << "\nThis can be shown by setting the equation to be less than p.neighbour and rearranging for geometry.cell_centre_distance";

    const auto interpolated_value_greater_than = face_interpolator.interp(p_greater);
    const auto interpolated_value_equal = face_interpolator.interp(p_equal);

    EXPECT_GT(p_greater.neighbour,interpolated_value_greater_than);
    EXPECT_DOUBLE_EQ(p_equal.neighbour,interpolated_value_equal);
}

TEST_F(TestFaceInterpolator, ResultIsHarmonicMeanForSpecialCase)
{
    constexpr Geometry geometry {.to_face = Pair{.owner = 1.0,.neighbour = 1.0},
                        .elevation = Pair{.owner = 1500.0,.neighbour = 1600.0},
                        .cell_centre_distance = 2.0};

    const faceInterpolator face_interpolator {geometry};
    const auto p = Pair{.owner = 3.0, .neighbour = 1.0};
    const auto harmonic_mean = 2 * p.owner * p.neighbour / (p.owner + p.neighbour);

    const auto result = face_interpolator.interp(p);


    EXPECT_DOUBLE_EQ(result,harmonic_mean);
}

TEST_F(TestFaceInterpolator, CompareAgainstManualCompute)
{
    const auto p = Pair{.owner = 1.0,.neighbour = 10.0};
    const auto expected = geometry.cell_centre_distance * (p.owner * p.neighbour) /
        (geometry.to_face.owner * p.neighbour + geometry.to_face.neighbour * p.owner);

    const auto result = face_interpolator.interp(p);

    EXPECT_DOUBLE_EQ(expected,result);
}