#include "SoilMoistureSolverCore/Data.hpp"
#include "SoilMoistureSolverCore/Details.hpp"
#include "gtest/gtest.h"

class MockData {};
template<typename T>
struct MockFace
{
    MockData data;
    template<typename... Args>
    T& make_module_data(const std::string&, Args&&... args) { return data;}
    T& get_module_data(const std::string&) {return data;}

    Vector_3 edge_unit_normal(size_t) { return Vector_3{};}
    Point_3 edge_midpoint(size_t) { return Point_3{};}
    Vector_3 downslope_dir() {return Vector_3{};}

    T neighbor(size_t) { return MockFace{};}
    size_t cell_global_id() { return 1u;}

};

using namespace SoilMoistureSolver::detail;
class TestBuildCellGeometry : public ::testing::Test {
protected:
    std::unique_ptr<MockFace<MockData>> face;
    Params p{   .lower_soil_depth = 3.0,
                .recharge_soil_depth = 2.0,
                .detention_depth = 1.0,
                .time_step_seconds = 3600};
    static constexpr auto Dims = orderedPair{.layer = 4u, .face = 2u};
};
template <ElementInterface E, orderedPair Dims = cellFacesAndVerticalLayers, size_t... LayerIndices>
struct foo{};
TEST_F(TestBuildCellGeometry,DoStuff)
{

    const auto geo = build_cell_geometry(face,  p, std::make_index_sequence<Dims.layer>());


}
