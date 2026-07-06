#include "SoilMoistureSolverCore/SoilMoistureSolverCore.hpp"
#include "StencilAssembly.hpp"

#include "gtest/gtest.h"
#include <iomanip>

using namespace SoilMoistureSolver::detail;

class TestHelperFunctions : public ::testing::Test
{
    protected:
    static constexpr auto op_top    = orderedPair{.layer = cellFacesAndVerticalLayers.layer - 1, .face = cellFacesAndVerticalLayers.face - 1};
    static constexpr auto op_bottom = orderedPair{.layer = 0u, .face = 3u};
    static constexpr auto op_interior = orderedPair{.layer = 1u, .face = 0u};

};

TEST_F(TestHelperFunctions, FaceIndexToNeighbourFromOrderedPair)
{
    using namespace SoilMoistureSolver::detail;

    constexpr auto target = Neighbour::Top;

    EXPECT_EQ(face_index_to_Neighbour(op_top), target);
    EXPECT_EQ(face_index_to_Neighbour(op_bottom), Neighbour::Bottom);
    EXPECT_EQ(face_index_to_Neighbour(op_interior), Neighbour::Lateral_0);
}

TEST_F(TestHelperFunctions, NeighbourToFaceIndex)
{
    using namespace SoilMoistureSolver::detail;

    EXPECT_EQ(Neighbour_to_face_index(Neighbour::Top), op_top.face);
    EXPECT_EQ(Neighbour_to_face_index(Neighbour::Bottom), op_bottom.face);
    EXPECT_EQ(Neighbour_to_face_index(Neighbour::Lateral_0), op_interior.face);
}

TEST_F(TestHelperFunctions, FaceIndexToNeighbourFromFaceIndex)
{
    using namespace SoilMoistureSolver::detail;

    EXPECT_EQ(face_index_to_Neighbour(op_top.face), Neighbour::Top);
    EXPECT_EQ(face_index_to_Neighbour(op_bottom.face), Neighbour::Bottom);
    EXPECT_EQ(face_index_to_Neighbour(op_interior.face), Neighbour::Lateral_0);
}

TEST_F(TestHelperFunctions, RoundTripFaceIndexToNeighbourAndBack)
{
    using namespace SoilMoistureSolver::detail;

    EXPECT_EQ(Neighbour_to_face_index(face_index_to_Neighbour(op_top)), op_top.face);
    EXPECT_EQ(Neighbour_to_face_index(face_index_to_Neighbour(op_bottom)), op_bottom.face);
    EXPECT_EQ(Neighbour_to_face_index(face_index_to_Neighbour(op_interior)), op_interior.face);
}

struct MatrixSizes
{
    size_t rows;
    size_t columns;
};
class TestLinearSystem
{
    MatrixSizes _s;
    std::vector<double> _d;
    std::vector<double> _rhs;
public:
    explicit TestLinearSystem(const MatrixSizes&& s) : _s(s), _d(_s.rows * _s.columns,0.0), _rhs(_s.rows,0.0) {};
    void matrixSumIntoGlobalValues(const int i, const int j, const double v)
    {
        _d.at(i + j *_s.rows) = v;
    };
    void rhsSumIntoGlobalValue(const int i, const double v)
    {
        _rhs.at(i) = v;
    };

    friend std::ostream& operator<<(std::ostream& os, const TestLinearSystem& sys);

    const std::vector<double>& get_matrix() { return _d; }
    const std::vector<double>& get_rhs() { return _rhs; }
};

std::ostream& operator<<(std::ostream& os, const TestLinearSystem& sys)
{
    os << "RHS Vector:\n";
    for (size_t row = 0; row < sys._rhs.size(); ++row)
        os << "  [" << std::setw(3) << row << "] " << std::setw(12) << sys._rhs[row] << "\n";

    os << "\nMatrix (" << sys._s.rows << " x " << sys._s.columns << "):\n";
    os << std::scientific << std::setprecision(6);
    for (size_t row = 0; row < sys._s.rows; ++row) {
        os << "  ";
        for (size_t col = 0; col < sys._s.columns; ++col) {
            os << std::setw(14) << sys._d[row + col * sys._s.rows];
        }
        os << "\n";
    }

    return os;
}

static_assert(math::optin::LinearSystem<TestLinearSystem>);

class TestStencilAssembly : public ::testing::Test {
    protected:
    static constexpr auto ROWS = 10u;
    static constexpr auto COLUMNS = 8u;
    std::unique_ptr<TestLinearSystem> _t = std::make_unique<TestLinearSystem>(MatrixSizes{.rows=ROWS,.columns= COLUMNS});
    TestStencilAssembly() = default;
};

TEST_F(TestStencilAssembly, LinearSystemZeroed)
{
    const auto& matrix = _t->get_matrix();
    const auto& rhs = _t->get_rhs();

    for (const auto m : matrix)
        EXPECT_EQ(m,0.0);

    for (const auto r : rhs)
        EXPECT_EQ(r,0.0);
}

TEST_F(TestStencilAssembly, LinearSystemPutToMatrix)
{
    constexpr auto value = 2.5;
    const auto row = 3u;
    const auto column = 6u;
    static_assert(row <= ROWS, "row constant violating definition");
    static_assert(column <= COLUMNS, "column constant violating definition");

    _t->matrixSumIntoGlobalValues(row, column, value);

    const auto& matrix = _t->get_matrix();
    const auto& rhs = _t->get_rhs();

    auto count = 0u;
    for (const auto m : matrix)
    {
        if (count == row + ROWS * column)
            EXPECT_EQ(m,value) << *_t;
        else
            EXPECT_EQ(m,0.0) << *_t;
        count++;
    }

    for (const auto r : rhs)
        EXPECT_EQ(r,0.0) << *_t;
}

TEST_F(TestStencilAssembly, LinearSystemPutToRHS)
{
    constexpr auto value = 2.5;
    constexpr auto row = 3u;
    _t->rhsSumIntoGlobalValue(row,value);
    const auto& matrix = _t->get_matrix();
    const auto& rhs = _t->get_rhs();

    auto count = 0u;
    for (const auto m : matrix)
    {
        ASSERT_EQ(m,0.0) << *_t << "\nCount: " << count << "\nm: " << m;
        count++;
    }

    count = 0u;
    for (const auto r : rhs)
    {
        if (count == row)
            EXPECT_EQ(r,value) << *_t;
        else
            EXPECT_EQ(r,0.0) << *_t;
        count++;
    }
};
