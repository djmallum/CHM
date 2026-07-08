#include "StencilAssembly.hpp"

#include "gtest/gtest.h"
#include <boost/math/policies/policy.hpp>
#include <boost/multi_array/base.hpp>
#include <iomanip>
#include <vector>
#include <optional>


struct TestLinearSystem {
    // Recorded calls
    struct MatrixCall { size_t i, j; double v; };
    struct RhsCall    { size_t i; double v; };

    std::vector<MatrixCall> matrix_calls;
    std::vector<RhsCall>    rhs_calls;

    // Concept interface (just records)
    void matrixSumIntoGlobalValues(const size_t i, const size_t j, const double v) {
        matrix_calls.push_back({i, j, v});
    }

    void rhsSumIntoGlobalValue(const size_t i, const double v) {
        rhs_calls.push_back({i, v});
    }
};

static_assert(math::optin::LinearSystem<TestLinearSystem>);

static constexpr size_t num_neighbours= 3;

template<typename T>
struct CountedValue {
    T value{};
    mutable size_t access_count = 0;

    CountedValue() = default;
    CountedValue(T v) : value(v) {}  // Add this
    operator T() const {
        ++access_count;
        return value;
    }

    CountedValue& operator=(const T& v) {
        value = v;
        return *this;
    }
};

struct SolverData
{
    CountedValue<size_t> _idx;
    std::array<std::optional<size_t>,num_neighbours> _neighbours;
    CountedValue<double> _diagonal;
    CountedValue<double> _off_diagonal;
    CountedValue<double> _boundary_diagonal;
    CountedValue<double> _boundary_rhs;
    CountedValue<double> _rhs;
    CountedValue<size_t> _top_face;
    CountedValue<size_t> _bottom_face;

    using donor_choice = math::without_donor_tag;
    using boundary_donor_choice = math::without_boundary_donor_tag;
    using boundary_choice = math::with_boundary_tag;
    using rhs_choice = math::with_rhs_tag;

    size_t idx() const { return _idx; }
    size_t neighbour_idx(const size_t f) const { return _neighbours[f].value(); }
    size_t top_face() const { return _top_face;}
    size_t bottom_face() const { return _bottom_face;}

    bool has_neighbour(const size_t f) const {
        return _neighbours[f].has_value();
    }

    double diagonal(const size_t) const { return _diagonal; }
    double off_diagonal(const size_t) const { return _off_diagonal; }
    double boundary_diagonal(size_t) const { return _boundary_diagonal;}
    double boundary_rhs(size_t) const { return _boundary_rhs;}
    double rhs(const size_t) const { return _rhs; }

};

class TestStencilAssembly : public ::testing::Test {
protected:
    static auto order_lambda()
    {
        return [](auto const& a, auto const& b)
        {
            if (a.i != b.i) return a.i < b.i;
            if constexpr (requires { a.j; }) {
                if (a.j != b.j) return a.j < b.j;
            }
            return a.v < b.v;
        };
    }
    using MatrixCall = TestLinearSystem::MatrixCall;
    using RhsCall = TestLinearSystem::RhsCall;
    TestLinearSystem linear_system{};
    SolverData solver_data{};
    struct Expected
    {
        std::vector<MatrixCall> matrix;
        std::vector<RhsCall>    rhs;
    };
    Expected expected;
    void enable_no_boundaries()
    {
        auto count = 1u;
        for (auto& neighbour : solver_data._neighbours)
        {
            neighbour.emplace(count);
            ++count;
        }

        build_expected_matrix_rhs();
    }

    void enable_with_boundaries()
    {
        solver_data._neighbours[1u].emplace(2u);

        build_expected_matrix_rhs();
    }

    /**
     * Build expected struct, matrix and rhs, representing A and b in Ax=b,
     *
     * Core Assumptions
     *
     * 1. Only looking at cell with index 0.
     * 2. Access members of solver_data with .value to avoid incremented the access counter and interfering with test isolation
     */
    void build_expected_matrix_rhs()
    {
        for (const auto neighbour : solver_data._neighbours)
        {
            if (neighbour)
            {
                expected.matrix.push_back(MatrixCall{.i = 0, .j = 0, .v = solver_data._diagonal.value});
                expected.matrix.push_back(MatrixCall{.i = 0, .j = *neighbour, .v = solver_data._off_diagonal.value});
                expected.rhs.push_back(RhsCall{.i = 0, .v = solver_data._rhs.value});
            }
            else
            {
                expected.matrix.push_back(MatrixCall{.i = 0, .j = 0, .v = solver_data._boundary_diagonal.value});
                expected.rhs.push_back(RhsCall{.i = 0, .v = solver_data._boundary_rhs.value});
            }
        }

        std::ranges::sort(expected.matrix, order_lambda());
        std::ranges::sort(expected.rhs, order_lambda());
    }
    void do_setup()
    {
        solver_data = SolverData{._idx = 0,
                                 ._diagonal = {2.5},
                                 ._off_diagonal = {1.0},
                                 ._boundary_diagonal = {-9999.0},
                                 ._boundary_rhs = {-3333.0},
                                 ._rhs = {33.0},
                                 ._top_face = {2},
                                 ._bottom_face = {3}};

    }

    std::pair<std::vector<MatrixCall>,std::vector<RhsCall>> get_finished_matrix()
    {
        math::LinearAlgebra::assemble_all_neighbours<num_neighbours>(linear_system, solver_data);

        std::ranges::sort(linear_system.matrix_calls,order_lambda());
        std::ranges::sort(linear_system.rhs_calls,order_lambda());

        const auto& matrix = linear_system.matrix_calls;
        const auto& rhs = linear_system.rhs_calls;

        return {std::move{matrix},std::move{rhs}};
    }
    void SetUp() override { do_setup(); };
};

TEST_F(TestStencilAssembly, NoBoundaryCheckAccessNumbers)
{
    enable_no_boundaries();

    math::LinearAlgebra::assemble_all_neighbours<num_neighbours>(linear_system, solver_data);

    const auto expected_access_count = solver_data._neighbours.size();
    EXPECT_EQ(solver_data._idx.access_count,expected_access_count)
        << "_idx is accessed directly via idx() once per neighbour";
    EXPECT_EQ(solver_data._diagonal.access_count,expected_access_count)
        << "diagonal() should be called once per neighbour";
    EXPECT_EQ(solver_data._off_diagonal.access_count,expected_access_count)
        << "off_diagonal() should be called once per neighbour";
    EXPECT_EQ(solver_data._rhs.access_count,expected_access_count)
        << "rhs() should be called once per neighbour";
    EXPECT_EQ(solver_data._boundary_diagonal.access_count,0u)
        << "Because _neighbours.size() == num_neighbours, none of the faces are boundaries and so this should never be accessed";
    EXPECT_EQ(solver_data._boundary_rhs.access_count,0u)
        << "Because _neighbours.size() == num_neighbours, none of the faces are boundaries and so this should never be accessed";
    EXPECT_EQ(solver_data._top_face.access_count,0u)
        << "top_face() is only accessed on boundary faces, and there are no boundary faces for this test";
    EXPECT_EQ(solver_data._bottom_face.access_count,0u)
        << "bottom_face() is only accessed on boundary faces, and there are no boundary faces for this test";
}

TEST_F(TestStencilAssembly, WithBoundaryCheckAccessNumbers)
{
    enable_with_boundaries();

    math::LinearAlgebra::assemble_all_neighbours<num_neighbours>(linear_system, solver_data);

    EXPECT_EQ(solver_data._boundary_diagonal.access_count,2u)
        << "2 faces of " << num_neighbours << " should be boundaries. Check set up to confirm";
    EXPECT_EQ(solver_data._boundary_rhs.access_count,2u)
        << "2 faces of " << num_neighbours << " should be boundaries. Check set up to confirm";
}

TEST_F(TestStencilAssembly, NoBoundaryCheckValuesSet)
{
    enable_no_boundaries();

    auto [matrix,rhs] = get_finished_matrix();

    ASSERT_EQ(matrix.size(),expected.matrix.size());
    ASSERT_EQ(rhs.size(),expected.rhs.size());

    for (size_t i = 0; i < expected.matrix.size(); ++i)
    {
        EXPECT_EQ(expected.matrix[i].i,matrix[i].i) << "Entry: " << i;
        EXPECT_EQ(expected.matrix[i].j,matrix[i].j) << "Entry: " << i;
        EXPECT_DOUBLE_EQ(expected.matrix[i].v,matrix[i].v) << "Entry: " << i;
    }
    for (size_t i = 0; i < expected.rhs.size(); ++i)
    {
        EXPECT_EQ(expected.rhs[i].i,rhs[i].i) << "Entry: " << i;
        EXPECT_DOUBLE_EQ(expected.rhs[i].v,rhs[i].v) << "Entry: " << i;
    }
}

TEST_F(TestStencilAssembly, WithBoundaryCheckValuesSet)
{
    enable_with_boundaries();

    auto [matrix,rhs] = get_finished_matrix();

    ASSERT_EQ(matrix.size(),expected.matrix.size());
    ASSERT_EQ(rhs.size(),expected.rhs.size());

    for (size_t i = 0; i < expected.matrix.size(); ++i)
    {
        EXPECT_EQ(expected.matrix[i].i,matrix[i].i) << "Entry: " << i;
        EXPECT_EQ(expected.matrix[i].j,matrix[i].j) << "Entry: " << i;
        EXPECT_DOUBLE_EQ(expected.matrix[i].v,matrix[i].v) << "Entry: " << i;
    }
    for (size_t i = 0; i < expected.rhs.size(); ++i)
    {
        EXPECT_EQ(expected.rhs[i].i,rhs[i].i) << "Entry: " << i;
        EXPECT_DOUBLE_EQ(expected.rhs[i].v,rhs[i].v) << "Entry: " << i;
    }
}
TEST_F(TestStencilAssembly, WithBoundaryNotFirstRowCheckValuesSet)
{
    enable_with_boundaries();

    auto [matrix,rhs] = get_finished_matrix();

    ASSERT_EQ(matrix.size(),expected.matrix.size());
    ASSERT_EQ(rhs.size(),expected.rhs.size());

}

