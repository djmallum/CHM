#include "SoilMoistureSolverCore/SoilMoistureSolverCore.hpp"
#include "gtest/gtest.h"

using namespace SoilMoistureSolver::detail;

class TestSoilMoistureHelperFunctions : public ::testing::Test
{
    protected:
    static constexpr auto op_top    = orderedPair{.layer = cellFacesAndVerticalLayers.layer - 1, .face = cellFacesAndVerticalLayers.face - 1};
    static constexpr auto op_bottom = orderedPair{.layer = 0u, .face = 3u};
    static constexpr auto op_interior = orderedPair{.layer = 1u, .face = 0u};

    static_assert(Neighbour::Lateral_0 == 0u);
    static_assert(Neighbour::Lateral_1 == 1u);
    static_assert(Neighbour::Lateral_2 == 2u);
    static_assert(Neighbour::Top == 3u);
    static_assert(Neighbour::Bottom == 4u);

};

TEST_F(TestSoilMoistureHelperFunctions, FaceIndexToNeighbourFromOrderedPair)
{
    using namespace SoilMoistureSolver::detail;

    constexpr auto target = Neighbour::Top;

    EXPECT_EQ(face_index_to_Neighbour(op_top), target);
    EXPECT_EQ(face_index_to_Neighbour(op_bottom), Neighbour::Bottom);
    EXPECT_EQ(face_index_to_Neighbour(op_interior), Neighbour::Lateral_0);
}

TEST_F(TestSoilMoistureHelperFunctions, NeighbourToFaceIndex)
{
    using namespace SoilMoistureSolver::detail;

    EXPECT_EQ(Neighbour_to_face_index(Neighbour::Top), op_top.face);
    EXPECT_EQ(Neighbour_to_face_index(Neighbour::Bottom), op_bottom.face);
    EXPECT_EQ(Neighbour_to_face_index(Neighbour::Lateral_0), op_interior.face);
}

TEST_F(TestSoilMoistureHelperFunctions, FaceIndexToNeighbourFromFaceIndex)
{
    using namespace SoilMoistureSolver::detail;

    EXPECT_EQ(face_index_to_Neighbour(op_top.face), Neighbour::Top);
    EXPECT_EQ(face_index_to_Neighbour(op_bottom.face), Neighbour::Bottom);
    EXPECT_EQ(face_index_to_Neighbour(op_interior.face), Neighbour::Lateral_0);
}

TEST_F(TestSoilMoistureHelperFunctions, FaceIndexToNeighbourAndBack)
{
    using namespace SoilMoistureSolver::detail;

    EXPECT_EQ(Neighbour_to_face_index(face_index_to_Neighbour(op_top)), op_top.face);
    EXPECT_EQ(Neighbour_to_face_index(face_index_to_Neighbour(op_bottom)), op_bottom.face);
    EXPECT_EQ(Neighbour_to_face_index(face_index_to_Neighbour(op_interior)), op_interior.face);
}
