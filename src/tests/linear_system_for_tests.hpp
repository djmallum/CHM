#pragma once

#include <vector>

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
