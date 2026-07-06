#include "StencilAssembly.hpp"

#include "gtest/gtest.h"
#include <iomanip>
#include <vector>
#include <numeric>

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
    static constexpr auto COLUMNS = 10u;
    static constexpr auto face_number = 4u;
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

class SolverData
{
public:
    SolverData() = default;

    size_t idx() const {}
    bool has_neighbour(const size_t f) const {}
    size_t neighbour_idx(const size_t f) const {}
    double diagonal(const size_t f) const {}
    double off_diagonal(const size_t f) const {}

    constexpr size_t top_face() const { return 3u;}
    constexpr size_t bottom_face() const { return 4u;}

    using donor_choice = math::without_donor_tag;
    using boundary_donor_choice = math::without_boundary_donor_tag;
    using boundary_choice = math::without_boundary_tag;

    using rhs_choice = math::with_rhs_tag;

};

template<size_t T>
static void assign_neighbours(TestLinearSystem& ls, const size_t idx, const std::array<unsigned, T> neigh)
{
    for (const auto n : neigh)
    {
        ls.matrixSumIntoGlobalValues(idx, n, n);
    }
}
static void build_expected_matrix(TestLinearSystem& ls)
{
    /*
     * For testing, construct the following system
     *
     * 4 triangle system with 3 layers
     *
     * Imagine a single, central, vertically stacked column of three triangular prisms
     *
     * With a similar stacking of triangular prisms at each of the three side faces of the central column.
     *
     * Each vertical layer has 4 triangles, and with three layers for a total of 12 cells.
     *
     * Each cell has 5 faces, only the middle triangle in the central stack has no impact of boundary conditions.
     *
     * Indexing is as follows: 0 for bottom centre, then 1, 2, 3 in a counter-clockwise ordering.
     *
     * Layer 2 has 4, then 5, 6, 7.
     *
     * Layer 3 has 8, then 9, 10, 11
     *
     * Final Matrix is 12x12, 12 equations per cell and 12 possibly contributing cells. Only neighbours will contribute.
     */

    constexpr auto num_cells = 12u;
    constexpr auto num_neighbours = 3u;
    const std::vector<std::vector<size_t>> real_neighbors = {{
        {1, 2, 3, 4},       // Cell 0
        {0, 5},             // Cell 1
        {0, 6},             // Cell 2
        {0, 7},             // Cell 3
        {0, 5, 6, 7, 8},    // Cell 4
        {4, 1, 9},          // Cell 5
        {4, 2, 10},         // Cell 6
        {4, 3, 11},         // Cell 7
        {4, 9, 10, 11},     // Cell 8
        {8, 5},             // Cell 9
        {8, 6},             // Cell 10
        {8, 7},             // Cell 11
    }};
    constexpr auto indices = []() {
        std::array<size_t,num_cells> arr;
        std::iota(arr.begin(), arr.end(), 0u);
        return arr;
    }();
    static_assert(indices[0] == 0);
    static_assert(indices[11] == 11);

    for (const auto idx : indices)
    {
        // diagonal
        ls.matrixSumIntoGlobalValues(idx,idx,idx);
        const auto neighbours = real_neighbors[idx];
        for (const auto neighbour : neighbours)
        {

        }
    }

    // Cell 0: bottom layer, central
    size_t idx = indices[0];
    ls.rhsSumIntoGlobalValue(idx,-static_cast<int>(idx));
    assign_neighbours(ls, idx, std::array{1u,2u,3u,4u});

    // Cell 1: Bottom layer, edge
    size_t idx = indices[1];
    ls.rhsSumIntoGlobalValue(idx,-3u*static_cast<int>(idx));
    assign_neighbours(ls, idx, std::array{0u,5u});


}

TEST_F(TestStencilAssembly, BuildMatrixRhsThroughPublicInterface)
{
    //
    auto solver_data = SolverData{};
    auto linear_system = TestLinearSystem{*_t};

    build_expected_matrix(linear_system);



    math::LinearAlgebra::assemble_all_neighbours<face_number>(*_t, solver_data);




}
