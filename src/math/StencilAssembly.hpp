// math/StencilAssembly.hpp
//
// Generic free functions that assemble a nearest-neighbour stencil into 
// a matrix object that satisfies the LinearSystem concept, for a cell-data 
// object that satisfies the CellStencil concept.
//
// Designed for a two- or three-dimensional grid aligned with the CHM mesh, 
// consisting of triangular prisms, vertically stacked into columns. The code 
// here is modelled after the PBSM3D mesh allocation. Additional new models
// that get added to CHM over time could use this or refactor to accommodate
// more general schemes as necessary.
//
// Donor (upwind) advection is OPTIONAL: if the cell-data type does not
// model DonorScheme, the advection term is silently skipped via
// `if constexpr`. Implementers who don't need it just don't write it.

#pragma once

#include <concepts>
#include <cstddef>

namespace math::LinearAlgebra
{
	// container that holds a matrix, consistent with the NearestNeighborProblem
	// class syntax
    template <class S>
    concept LinearSystem = requires(S& s, std::size_t i, std::size_t j, double v) {
        { s.matrixSumIntoGlobalValues(i, j, v) };
        { s.rhsSumIntoGlobalValue(i, v) };
    };

    // f indexes a face of the cell (0..N-1). For a triangular prism:
    //   0..2 = lateral, 3 = top, 4 = bottom.
    template <class C>
    concept CellStencil = requires(const C& c, int f) {
        { c.idx() }                  -> std::convertible_to<std::size_t>;
        { c.has_neighbour(f) }       -> std::convertible_to<bool>;
        { c.neighbour_idx(f) }       -> std::convertible_to<std::size_t>;
        { c.diagonal(f) }            -> std::convertible_to<double>;
        { c.off_diagonal(f) }        -> std::convertible_to<double>;
        { c.boundary_diagonal(f) }   -> std::convertible_to<double>;
    };

    // OPTIONAL: Donor scheme for terms.
    //
    // Implementers who want contributions that are added to either the 
	// diagonal or off-diagonal based on a time-step-dependent check:
    //   - donor_on_diag(f): true when the term belongs on the
    //                       diagonal, false when it belongs on
    //                       the off-diagonal.
    //   - advect(f): the signed advection coefficient ( e.g. -A[f]*udotm[f] ).
    template <class C>
    concept DonorScheme = requires(const C& c, int f) {
        { c.donor_on_diag(f) }       -> std::convertible_to<bool>;
        { c.donor_term(f) }          -> std::convertible_to<double>;
    };

    // ---- 4. OPTIONAL: vertical extension for stacked / extruded meshes.
    //
    // Layer category lets the assembler know whether the cell sits at
    // the bottom, top, or middle of a stacked column. Top/bottom faces
    // are reported via top_face()/bottom_face() so the implementer
    // chooses their indexing convention.
    enum class VertLayer { Bottom, Middle, Top };  // _____ of the stack

    template <class C>
    concept VerticallyStacked = requires(const C& c) {
        { c.layer() }                       -> std::same_as<VertLayer>;
        { c.top_face() }                    -> std::convertible_to<int>;
        { c.bottom_face() }                 -> std::convertible_to<int>;
        { c.bottom_boundary_diagonal() }    -> std::convertible_to<double>;
        { c.bottom_boundary_rhs() }         -> std::convertible_to<double>;
        { c.top_boundary_rhs_in() }         -> std::convertible_to<double>; // outflow case
        requires requires { { c.top_boundary_rhs_out() }        -> std::convertible_to<double>;} == DonorScheme<C>; // inflow case
    };

    // ============================================================
    //  Helpers
    // ============================================================

    namespace detail
    {
        // Assemble one face of the upwind FVM stencil into (S, c).
        // Pure: same shape for lateral or vertical face -- the only
        // thing that varies is which face index `f` is.
        template <LinearSystem S, CellStencil C>
        void assemble_face(S& sys, const C& c, int f)
        {
            const auto i = c.idx();

            if (!c.has_neighbour(f)) {
                sys.matrixSumIntoGlobalValues(i, i, c.boundary_diagonal(f));
                return;
            }

            const auto j = c.neighbour_idx(f);

            if constexpr (DonorScheme<C>) {
                const double a = c.donor_term(f);
                if (c.donor_on_diag(f)) {
                    sys.matrixSumIntoGlobalValues(i, i, c.diagonal(f) + a);
                    sys.matrixSumIntoGlobalValues(i, j, c.off_diagonal(f));
                } else {
                    sys.matrixSumIntoGlobalValues(i, i, c.diagonal(f));
                    sys.matrixSumIntoGlobalValues(i, j, c.off_diagonal(f) + a);
                }
            } else {
                sys.matrixSumIntoGlobalValues(i, i, c.diagonal(f));
                sys.matrixSumIntoGlobalValues(i, j, c.off_diagonal(f));
            }
        }
    }

    // ============================================================
    //  Public free functions
    // ============================================================

    // Assemble all lateral faces (0 .. NLateral-1) of a cell.
    template <LinearSystem S, CellStencil C, int NLateral = 3>
    void lateral_neighbours(S& sys, const C& c)
    {
        for (int f = 0; f < NLateral; ++f)
            detail::assemble_face(sys, c, f);
    }

    // Assemble vertical (top + bottom) faces with layer-aware boundary
    // conditions. Only enabled when C models VerticallyStacked.
    template <LinearSystem S, CellStencil C>
        requires VerticallyStacked<C>
    void vertical_neighbours(S& sys, const C& c)
    {
        const auto i  = c.idx();
        const int top = c.top_face();
        const int bot = c.bottom_face();

        switch (c.layer()) {
        case VertLayer::Bottom:
            // Bottom face is a Dirichlet/flux boundary, top face is interior.
            sys.matrixSumIntoGlobalValues(i, i, c.bottom_boundary_diagonal());
            sys.rhsSumIntoGlobalValue   (i,    c.bottom_boundary_rhs());
            detail::assemble_face(sys, c, top);
            break;

        case VertLayer::Top:
            // Top face is a (precip) boundary, bottom face is interior.
            if constexpr (DonorScheme<C>) {
                if (c.donor_on_diag(top)) {
                    sys.matrixSumIntoGlobalValues(i, i, c.diagonal(top) + c.donor_term(top));
                    sys.rhsSumIntoGlobalValue   (i,    c.top_boundary_rhs_in());
                } else {
                    sys.matrixSumIntoGlobalValues(i, i, c.diagonal(top));
                    sys.rhsSumIntoGlobalValue   (i,    c.top_boundary_rhs_out());
                }
            } else {
                sys.matrixSumIntoGlobalValues(i, i, c.diagonal(top));
                sys.rhsSumIntoGlobalValue   (i,    c.top_boundary_rhs_in());
            }
            detail::assemble_face(sys, c, bot);
            break;

        case VertLayer::Middle:
            detail::assemble_face(sys, c, top);
            detail::assemble_face(sys, c, bot);
            break;
        }
    }

} // namespace math::LinearAlgebra
