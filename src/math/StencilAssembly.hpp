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
    };
}

// -----------------------------------------------------------------------------
// Opt-in / Opt-out generator
// -----------------------------------------------------------------------------
//
// Usage: DEFINE_OPTIONAL_FEATURE(Name, MEMBER_CHECK)
//   - Generates NameOptIn<NameOptOut, NameMember, NameChoiceMade>
//   - MEMBER_CHECK is a requires-expression fragment of the form:
//       (const C& c, int f) { { c.boundary_diagonal(f) } -> std::convertible_to<double>; }
//   - For features with different signatures, write the MEMBER_CHECK accordingly.
//
// The generator uses nested-type tokens by default (e.g. `using BoundaryDiagonalOptIn = void;`).
// It is straightforward to extend detection to boolean flags if you prefer.
//
// -----------------------------------------------------------------------------

#define DEFINE_OPTIONAL_FEATURE(NAME, MEMBER_CHECK)                                  \
namespace math::optin {                                                              \
                                                                                     \
    template <class C>                                                               \
    concept NAME##OptOut = requires { typename C::NAME##OptOut; };                   \
                                                                                     \
    template <class C>                                                               \
    concept NAME##Member = requires MEMBER_CHECK;                                    \
                                                                                     \
    template <class C>                                                               \
    concept NAME##ChoiceMade = NAME##Member<C> || NAME##OptOut<C>;                    \
                                                                                     \
} /* namespace math::optin */

// boundary_diagonal(int) -> double
DEFINE_OPTIONAL_FEATURE(BoundaryDiagonal,
    (const C& c, int f) { { c.boundary_diagonal(f) } -> std::convertible_to<double>; })

// boundary_off_diagonal(int) -> double
DEFINE_OPTIONAL_FEATURE(BoundaryOffDiagonal,
    (const C& c, int f) { { c.boundary_off_diagonal(f) } -> std::convertible_to<double>; })

// boundary_rhs(int) -> double
DEFINE_OPTIONAL_FEATURE(BoundaryRHS,
    (const C& c, int f) { { c.boundary_rhs(f) } -> std::convertible_to<double>; })

// rhs(int) -> double  (inhomogeneous RHS for interior faces)
DEFINE_OPTIONAL_FEATURE(RHS,
    (const C& c, int f) { { c.rhs(f) } -> std::convertible_to<double>; })

DEFINE_OPTIONAL_FEATURE(DonorScheme,
    (const C& c, int f) {
    { c.donor_on_diag(f) } -> std::convertible_to<bool>;
    { c.donor_term(f) }    -> std::convertible_to<double>;
    }
    )

// // Donor scheme: donor_on_diag(int) -> bool, donor_term(int) -> double
// namespace math::optin {
// template <class C>
// concept DonorSchemeMember = requires(const C& c, int f) {
//     { c.donor_on_diag(f) } -> std::convertible_to<bool>;
//     { c.donor_term(f) }    -> std::convertible_to<double>;
// };
//
// template <class C>
// concept DonorSchemeOptIn = requires { typename C::DonorSchemeOptIn; };
//
// template <class C>
// concept DonorSchemeOptOut = requires { typename C::DonorSchemeOptOut; };
//
// template <class C>
// concept DonorSchemeChoiceMade = DonorSchemeOptIn<C> || DonorSchemeOptOut<C>;
// }

    // template <class C>
    // concept BoundaryDiagonalSpecial = CellStencil<C> && requires (const C& c, int f)
    // {
    //     { c.boundary_diagonal(f) }   -> std::convertible_to<double>;
    // };
    //
    // template <class C>
    // concept BoundaryOffDiagonalSpecial = CellStencil<C> && requires (const C& c, int f)
    // {
    //     { c.boundary_off_diagonal(f) } -> std::convertible_to<double>;
    // };
    //
    // template <class C>
    // concept BoundaryRHSSpecial = CellStencil<C> && requires (const C& c, int f)
    // {
    //     { c.boundary_rhs(f) } -> std::convertible_to<double>;
    // };
    //
    // template <class C>
    // concept InHomogeneous = CellStencil<C> && requires (const C& c, int f)
    // {
    //     { c.rhs(f) } -> std::convertible_to<double>;
    // };

    // template <class C>
    // concept Homogeneous = !InHomogeneous<C>;

// Undefine helper macro to avoid leaking
#undef DEFINE_OPTIONAL_FEATURE

// -----------------------------------------------------------------------------
// Re-exported convenience concepts for library users
// -----------------------------------------------------------------------------
namespace math::LinearAlgebra
{
    // Opt-in/out/member/choice concepts are available under math::optin namespace.
    // For convenience, provide short aliases used in the assembler below.

    template <class C>
    concept BoundaryDiagonalOptOut = optin::BoundaryDiagonalOptOut<C>;

    template <class C>
    concept BoundaryDiagonalMember = optin::BoundaryDiagonalMember<C>;

    template <class C>
    concept BoundaryDiagonalChoiceMade = optin::BoundaryDiagonalChoiceMade<C>;

    template <class C>
    concept BoundaryOffDiagonalOptOut = optin::BoundaryOffDiagonalOptOut<C>;

    template <class C>
    concept BoundaryOffDiagonalMember = optin::BoundaryOffDiagonalMember<C>;

    template <class C>
    concept BoundaryOffDiagonalChoiceMade = optin::BoundaryOffDiagonalChoiceMade<C>;

    template <class C>
    concept BoundaryRHSOptOut = optin::BoundaryRHSOptOut<C>;

    template <class C>
    concept BoundaryRHSMember = optin::BoundaryRHSMember<C>;

    template <class C>
    concept BoundaryRHSChoiceMade = optin::BoundaryRHSChoiceMade<C>;

    template <class C>
    concept RHSOptOut = optin::RHSOptOut<C>;

    template <class C>
    concept RHSMember = optin::RHSMember<C>;

    template <class C>
    concept RHSChoiceMade = optin::RHSChoiceMade<C>;

    template <class C>
    concept DonorSchemeOptOut = optin::DonorSchemeOptOut<C>;

    template <class C>
    concept DonorSchemeMember = optin::DonorSchemeMember<C>;

    template <class C>
    concept DonorSchemeChoiceMade = optin::DonorSchemeChoiceMade<C>;

}

    // OPTIONAL: Donor scheme for terms.
    //
    // Implementers who want contributions that are added to either the 
	// diagonal or off-diagonal based on a time-step-dependent check:
    //   - donor_on_diag(f): true when the term belongs on the
    //                       diagonal, false when it belongs on
    //                       the off-diagonal.
    // //   - advect(f): the signed advection coefficient ( e.g. -A[f]*udotm[f] ).
    // template <class C>
    // concept DonorScheme = requires(const C& c, int f) {
    //     { c.donor_on_diag(f) }       -> std::convertible_to<bool>;
    //     { c.donor_term(f) }          -> std::convertible_to<double>;
    // };

    // ---- 4. OPTIONAL: vertical extension for stacked / extruded meshes.
    //
    // Layer category lets the assembler know whether the cell sits at
    // the bottom, top, or middle of a stacked column. Top/bottom faces
    // are reported via top_face()/bottom_face() so the implementer
    // chooses their indexing convention.
namespace math::LinearAlgebra {
    enum class VertLayer { Bottom, Middle, Top };  // _____ of the stack

    template <class C>
    concept VerticallyStacked = requires(const C& c) {
        { c.layer() }                       -> std::same_as<VertLayer>;
        { c.top_face() }                    -> std::convertible_to<int>;
        { c.bottom_face() }                 -> std::convertible_to<int>;
        { c.bottom_boundary_diagonal() }    -> std::convertible_to<double>;
        { c.bottom_boundary_off_diagonal() }    -> std::convertible_to<double>;
        { c.bottom_boundary_rhs() }         -> std::convertible_to<double>;
        { c.top_boundary_rhs_in() }         -> std::convertible_to<double>; // outflow case
        requires requires { { c.top_boundary_rhs_out() }        -> std::convertible_to<double>;} == !DonorSchemeOptOut<C>; // inflow case
    };

    // ============================================================
    //  Helpers
    // ============================================================

    namespace detail
    {
    template <LinearSystem S, CellStencil C> void normal_face(S& sys, const C& c, int f)
    {
        const auto i = c.idx();
        const auto j = c.neighbour_idx(f);
        sys.matrixSumIntoGlobalValues(i, i, c.diagonal(f)); // diagonal so row = col
        sys.matrixSumIntoGlobalValues(i, j, c.off_diagonal(f)); // off diagonal so row != col
        if constexpr (RHSMember<C>)
        {
            sys.rhsSumIntoGlobalValue(i, c.rhs(f));
        }
    }
    template <LinearSystem S, CellStencil C> void boundary_conditions(S& sys, const C& c, int f)
    {
        const auto i = c.idx();

        if constexpr (BoundaryDiagonalMember<C>)
        {
            sys.matrixSumIntoGlobalValues(i, i, c.boundary_diagonal(f));
        }

        // Boundary off-diagonal handling
        if constexpr (BoundaryOffDiagonalMember<C>)
        {
            const auto j = c.neighbour_idx(f);
            sys.matrixSumIntoGlobalValues(i, j, c.boundary_off_diagonal(f));
        }

        // Boundary RHS handling
        if constexpr (BoundaryRHSMember<C>)
        {
            sys.rhsSumIntoGlobalValue(i, c.boundary_rhs(f));
        }
    }
    template <LinearSystem S, CellStencil C> void do_donor(S& sys, const C& c, int f)
    {
        const double a = c.donor_term(f);
        const auto i = c.idx();
        const auto j = c.neighbour_idx(f);
        if (c.donor_on_diag(f))
        {
            sys.matrixSumIntoGlobalValues(i, i, c.diagonal(f) + a);
            sys.matrixSumIntoGlobalValues(i, j, c.off_diagonal(f));
        }
        else
        {
            sys.matrixSumIntoGlobalValues(i, i, c.diagonal(f));
            sys.matrixSumIntoGlobalValues(i, j, c.off_diagonal(f) + a);
        }
    }
#define FEATURE_HOWTO \
"To opt out, add `using FeatureOptOut = void;` to your cell type. " \
"Otherwise provide the required member(s)."
    template <CellStencil C> static void check_opt_ins()
    {
        // Ensure the implementer explicitly chose opt-in or opt-out for all features
        static_assert(BoundaryDiagonalChoiceMade<C>, FEATURE_HOWTO);
        static_assert(BoundaryOffDiagonalChoiceMade<C>, FEATURE_HOWTO);
        static_assert(BoundaryRHSChoiceMade<C>, FEATURE_HOWTO);
        static_assert(RHSChoiceMade<C>, FEATURE_HOWTO);
        static_assert(DonorSchemeChoiceMade<C>, FEATURE_HOWTO);
    }

    template <LinearSystem S, CellStencil C> void interior_face(S& sys, const C& c, int f)//int f, const auto i, const auto j)
    {
        if constexpr (DonorSchemeMember<C>)
        {
            do_donor(sys, c, f);
        }
        else
        {
            normal_face(sys, c, f);
        }
    }
    // Assemble one face of the upwind FVM stencil into (S, c).
        // Pure: same shape for lateral or vertical face -- the only
        // thing that varies is which face index `f` is.
        template <LinearSystem S, CellStencil C>
        void assemble_face(S& sys, const C& c, int f)
        {
            if (!c.has_neighbour(f)) {
                boundary_conditions(sys, c, f);
                return;
            }

            interior_face(sys, c, f);
        }
    }

    // ============================================================
    //  Public free functions
    // ============================================================

    // Assemble all lateral faces (0 .. NLateral-1) of a cell.
    template <size_t Neighbours = 3, LinearSystem S, CellStencil C>
    void lateral_neighbours(S& sys, const C& c)
    {
        check_opt_ins<C>();

        for (size_t f = 0; f < Neighbours; ++f)
            detail::assemble_face(sys, c,
                static_cast<int>(f));
    }

    // Assemble vertical (top + bottom) faces with layer-aware boundary
    // conditions. Only enabled when C models VerticallyStacked.
    template <LinearSystem S, CellStencil C>
        requires VerticallyStacked<C>
    void vertical_neighbours(S& sys, const C& c)
    {
        check_opt_ins<C>();

        // Now do runtime selection of which face is boundary vs interior.
        const int top = c.top_face();
        const int bot = c.bottom_face();

        switch (c.layer()) {
        case VertLayer::Bottom:
            // bottom is boundary, top is interior
            detail::boundary_conditions(sys, c, bot);
            detail::interior_face(sys, c, top);
            break;

        case VertLayer::Top:
            // top is boundary, bottom is interior
            detail::boundary_conditions(sys, c, top);
            detail::interior_face(sys, c, bot);
            break;

        case VertLayer::Middle:
            // both faces are interior
            detail::interior_face(sys, c, bot);
            detail::interior_face(sys, c, top);
            break;
        }
        // static_assert(DonorSchemeChoiceMade<C>,FEATURE_HOWTO);
        //
        // const int top = c.top_face();
        // const int bot = c.bottom_face();
        //
        // switch (c.layer()) {
        //     case VertLayer::Bottom:
        //         detail::boundary_conditions(sys, c, bot);
        //         detail::interior_face(sys, c, top);
        //         break;
        //     case VertLayer::Top:
        //         detail::boundary_conditions(sys, c, top);
        //         detail::interior_face(sys, c, bot);
        //         break;
        //     case VertLayer::Middle:
        //         detail::interior_face(sys, c, bot);
        //         detail::interior_face(sys, c, top);
        //         break;
        // }
    }

#undef FEATURE_HOWTO
} // namespace math::LinearAlgebra
// -----------------------------------------------------------------------------
// Just keeping this here for full version
// // Public free functions (lateral + vertical assembly)
// // -----------------------------------------------------------------------------
// namespace math::LinearAlgebra
// {
//     // Assemble all lateral faces (0 .. NLateral-1) of a cell.
//     template <LinearSystem S, CellStencil C, int NLateral = 3>
//     void lateral_neighbours(S& sys, const C& c)
//     {
//         for (int f = 0; f < NLateral; ++f)
//             detail::assemble_face(sys, c, f);
//     }
//
//     // Assemble vertical (top + bottom) faces with layer-aware boundary
//     // conditions. Only enabled when the cell type explicitly opts in to vertical stacking.
//     template <LinearSystem S, CellStencil C>
//         requires VerticallyStackedChoiceMade<C>
//     void vertical_neighbours(S& sys, const C& c)
//     {
//         // Ensure the implementer explicitly chose opt-in or opt-out for vertical stacking
//         static_assert(VerticallyStackedChoiceMade<C>,
//             "Cell type must declare either VerticallyStackedOptIn or VerticallyStackedOptOut");
//
//         if constexpr (!VerticallyStackedOptIn<C>) {
//             // If the cell explicitly opted out, this function should not be used.
//             // Provide a helpful static_assert to catch misuse at compile time.
//             static_assert(VerticallyStackedOptIn<C>,
//                 "vertical_neighbours called for a cell that opted out of vertical stacking");
//         }
//
//         // From here on we know VerticallyStackedOptIn<C> is true
//         static_assert(VerticallyStackedMember<C>,
//             "Cell declared VerticallyStackedOptIn but is missing required vertical members");
//
//         // If DonorScheme is also opted-in, ensure the vertical-specific outflow member exists
//         if constexpr (DonorSchemeOptIn<C>) {
//             static_assert(requires(const C& cc) { { cc.top_boundary_rhs_out() } -> std::convertible_to<double>; },
//                 "Cell declared DonorSchemeOptIn and VerticallyStackedOptIn but is missing top_boundary_rhs_out()");
//         }
//
//         const auto i  = c.idx();
//         const int top = c.top_face();
//         const int bot = c.bottom_face();
//
//         switch (c.layer()) {
//         case VertLayer::Bottom:
//             // Bottom face is a Dirichlet/flux boundary, top face is interior
//             sys.matrixSumIntoGlobalValues(i, i, c.bottom_boundary_diagonal());
//             sys.rhsSumIntoGlobalValue   (i,    c.bottom_boundary_rhs());
//             detail::assemble_face(sys, c, top);
//             break;
//
//         case VertLayer::Top:
//             // Top face is a (precip) boundary, bottom face is interior.
//             if constexpr (DonorSchemeOptIn<C>) {
//                 // donor scheme handles inflow/outflow split
//                 if (c.donor_on_diag(top)) {
//                     sys.matrixSumIntoGlobalValues(i, i, c.diagonal(top) + c.donor_term(top));
//                     sys.rhsSumIntoGlobalValue   (i,    c.top_boundary_rhs_in());
//                 } else {
//                     sys.matrixSumIntoGlobalValues(i, i, c.diagonal(top));
//                     sys.rhsSumIntoGlobalValue   (i,    c.top_boundary_rhs_out());
//                 }
//             } else {
//                 sys.matrixSumIntoGlobalValues(i, i, c.diagonal(top));
//                 sys.rhsSumIntoGlobalValue   (i,    c.top_boundary_rhs_in());
//             }
//             detail::assemble_face(sys, c, bot);
//             break;
//
//         case VertLayer::Middle:
//             detail::assemble_face(sys, c, top);
//             detail::assemble_face(sys, c, bot);
//             break;
//         }
//     }
//
// } // namespace math::LinearAlgebra