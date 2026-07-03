#pragma once

#include "Concepts.hpp"
#include "Data.hpp"
#include "Details.hpp"
#include "StencilAssembly.hpp"
#include "core.hpp"

namespace SoilMoistureSolver::detail
{
template<ElementInterface E>
class solverData
{
    E& face;
    data& d;
    const size_t z_idx;
    const Sizes sizes;
    const std::string_view ID;

    [[nodiscard]] double soil_water_capacity() const;
    double _K_unsat_lateral_face(size_t i, orderedPair op) const;
    double _K_unsat_vertical_Face(size_t i, orderedPair op) const;
    [[nodiscard]] double K_unsaturated(size_t) const;
public:
    solverData(data& d, E& face, size_t layer, const Sizes& sizes, std::string_view);

    [[nodiscard]] size_t idx() const ;
    [[nodiscard]] bool has_neighbour(size_t f) const;
    [[nodiscard]] size_t neighbour_idx(size_t f) const;

    static constexpr size_t top_face() { return static_cast<size_t>(Neighbour::Top); }
    static constexpr size_t bottom_face() { return static_cast<size_t>(Neighbour::Bottom); }

    // without tags
    using donor_choice = math::without_donor_tag;
    using boundary_donor_choice = math::without_boundary_donor_tag;
    using boundary_choice = math::without_boundary_tag;

    // Normal faces
    using rhs_choice = math::with_rhs_tag;
    [[nodiscard]] double rhs(size_t) const;
    [[nodiscard]] double diagonal(size_t f) const;
    [[nodiscard]] double off_diagonal(size_t f) const;

    // Boundaries
    using side_boundary_choice = math::with_side_boundary_tag;
    [[nodiscard]] static double side_boundary_diagonal(size_t);
    [[nodiscard]] static double side_boundary_off_diagonal(size_t);
    [[nodiscard]] double side_boundary_rhs(size_t) const;

    using bottom_boundary_choice = math::with_bottom_boundary_tag;
    [[nodiscard]] static double bottom_boundary_diagonal(size_t);
    [[nodiscard]] static double bottom_boundary_off_diagonal(size_t);
    [[nodiscard]] double bottom_boundary_rhs(size_t) const;

    using top_boundary_choice = math::with_top_boundary_tag;
    [[nodiscard]] static double top_boundary_diagonal(size_t);
    [[nodiscard]] static double top_boundary_off_diagonal(size_t);
    [[nodiscard]] double top_boundary_rhs(size_t) const;

};

template<ElementInterface E>
double solverData<E>::soil_water_capacity() const
{
    /*
     * Derivative of theta(psi) with respect to psi. Comes from the chain rule of the time derivative
     * Pointed out everywhere that doing this with no approximation or iteration will result
     * in poor physics, particularly near the wetting front. We are not resolving the wetting front.
     */
    //    const auto result =-std::pow(d.psi_n / d.air_entry_tension,-1/d.pore_size_dist_index) / (d.pore_size_dist_index * d.psi_n);

    return -std::pow(d.psi_n.at(z_idx) / d.air_entry_tension,-1/d.pore_size_dist_index) / (d.pore_size_dist_index * d.psi_n.at(z_idx)); //may need to correct this near saturation...
}

static double K_unsaturated_Campbell(const data& d, const size_t layer)
{
    return d.K_saturated * std::pow(d.air_entry_tension / d.psi_n.at(layer),2+3/d.pore_size_dist_index);
}
template <ElementInterface E> double solverData<E>::_K_unsat_lateral_face(size_t i, const orderedPair op) const
{
    const auto& neigh = face->neighbor(i);
    if (has_neighbour(i))
    {
        const auto& d_neigh = neigh->template get_module_data<data>(ID.data());

        // TODO testing this equation for accuracy AND behaviour near saturation and dry soil
        // Source is Campbell (1974)
        // Good source is also Deb and Shukla (2012)

        const Pair pair_to_interp{.owner = K_unsaturated_Campbell(d, op.layer),
                                  .neighbour = K_unsaturated_Campbell(d_neigh, op.layer)};

        return d.cell_info.value_at_face(op, pair_to_interp);
    }



    return K_unsaturated_Campbell(d,op.layer);

}
template <ElementInterface E> double solverData<E>::_K_unsat_vertical_Face(const size_t i, const orderedPair op) const
{
    if (has_neighbour(i))
    {
        size_t layer_to{};
        if (face_index_to_Neighbour(i) == Neighbour::Top)
        {
            layer_to = z_idx + 1;
        }
        else if (face_index_to_Neighbour(i) == Neighbour::Bottom)
        {
            layer_to = z_idx - 1;
        }
        else
        {
            CHM_THROW_EXCEPTION(module_error, "Should only be top or bottom face in this part of the code");
        };

        if (z_idx == 0 || z_idx == cellFacesAndVerticalLayers.layer - 1)
        {
            const std::string err = std::format("z_idx value indicates that this should be a boundary cell"
                                                " but the code reached here. z_idx = {}",
                                                z_idx);
            CHM_THROW_EXCEPTION(module_error, err);
        }

        const Pair pair_to_interp{.owner = K_unsaturated_Campbell(d, op.layer),
                                  .neighbour = K_unsaturated_Campbell(d, layer_to)};

        return d.cell_info.value_at_face(op, pair_to_interp);
    }

    if (const auto& boundary = d.cell_info.get_faceType<Boundary>(op); boundary.neighbour_type == Neighbour::Top)
    {
        constexpr std::string_view err = "Top boundary condition does not require K_saturated and yet it was requested";
        CHM_THROW_EXCEPTION(module_error, err.data());
    }

    return K_unsaturated_Campbell(d, op.layer);
}
template <ElementInterface E>
double solverData<E>::K_unsaturated(size_t i) const
{
    /* Campbell (1974)
     * Also see: Deb and Shukla (2012) for long list
     */
    // TODO K_unsaturated should also handle vertical neighbours too!!
    // TODO Must handle boundary as well

    const orderedPair op{.layer = z_idx, .face = static_cast<size_t>(i)};

    if (i < NUM_SIDES)
    {
        return _K_unsat_lateral_face(i, op);
    }

    return _K_unsat_vertical_Face(i, op);

}

template<ElementInterface E>
solverData<E>::solverData(data& d, E& face, const size_t layer, const Sizes& sizes, const std::string_view id)
    : face(face), d(d), z_idx(layer), sizes(sizes), ID(id)
{
}


template<ElementInterface E>
size_t solverData<E>::idx() const { return sizes.global * z_idx + face->cell_global_id; }

template<ElementInterface E>
bool solverData<E>::has_neighbour(const size_t f) const { return d.cell_info.get_neighbour_idx(orderedPair{.layer=z_idx,.face=f}).has_value(); }

template<ElementInterface E>
size_t solverData<E>::neighbour_idx(const size_t f) const
{
    const auto i = d.cell_info.get_neighbour_idx(orderedPair{.layer=z_idx,.face=f});
    if (!i)
    {
        const std::string err = std::format("neighbour_idx invoked for a face without a neighbour at index {} of face {}",
            f,face->cell_global_id);
        CHM_THROW_EXCEPTION(module_error,err);
    }
    return *i;
}

template<ElementInterface E>
double solverData<E>::rhs(const size_t f) const
{
    // const auto elevation = d.get_elevation(f);
    const auto& interior_face = d.cell_info.get_faceType<Interior>(orderedPair{.layer=z_idx,.face=f});

    // TODO neighbour.z - z = the distance between the neighbour centre and the centre of the current cell.
    // Therefore, when one defines alpha for this term, as written it must include the distance between
    return d.psi_n.at(z_idx) / cellFacesAndVerticalLayers.face + d.cell_info.get_coefficient(orderedPair{.layer=z_idx,.face=f}) / soil_water_capacity() *
        K_unsaturated(f) * (interior_face.geometry.elevation.neighbour - interior_face.geometry.elevation.owner);
}

template<ElementInterface E>
double solverData<E>::diagonal(const size_t f) const
{
    return 1.0 / cellFacesAndVerticalLayers.face + d.cell_info.get_coefficient(orderedPair{.layer=z_idx,.face=f}) / soil_water_capacity() * K_unsaturated(f); /* TODO everything about geometry or constant
                                                                      *  in time goes in alpha, could make it a type
                                                                      */
}

template<ElementInterface E>
double solverData<E>::off_diagonal(const size_t f) const
{
    return -d.cell_info.get_coefficient(orderedPair{.layer=z_idx,.face=f}) / soil_water_capacity() * K_unsaturated(f); // TODO see diagonal
    // TODO off diagonal contributions
}

template<ElementInterface E>
double solverData<E>::side_boundary_diagonal(size_t) {
    return 1.0 / cellFacesAndVerticalLayers.face;
}

template<ElementInterface E>
double solverData<E>::side_boundary_off_diagonal(const size_t)
{
    return 0.0;
}

template<ElementInterface E>
double solverData<E>::side_boundary_rhs(const size_t f) const
{
    if (f > 2)
        CHM_THROW_EXCEPTION(module_error,"side boundary must only be faces 0, 1, or 2");

    const auto& boundary_face = d.cell_info.get_faceType<Boundary>(orderedPair{.layer=z_idx,.face=f});
    const auto elevation_change = boundary_face.DeltaZ;

    return d.psi_n.at(z_idx) / cellFacesAndVerticalLayers.face + d.cell_info.get_coefficient(orderedPair{.layer=z_idx,.face=f}) / soil_water_capacity() * K_unsaturated(f) * elevation_change;
}

template<ElementInterface E>
double solverData<E>::bottom_boundary_diagonal([[maybe_unused]] size_t) {
    return 1.0 / cellFacesAndVerticalLayers.face;
}

template<ElementInterface E>
double solverData<E>::bottom_boundary_off_diagonal([[maybe_unused]] size_t)
{
    return 0.0;
}

template<ElementInterface E>
double solverData<E>::bottom_boundary_rhs(const size_t f) const
{
    // WARNING: alpha for the bottom boundary must not include the distance to the neighbour
    // cell centre. Impossible to obtain since it doesn't exist. But its worth being aware.
    return d.psi_n.at(z_idx) / cellFacesAndVerticalLayers.face + d.cell_info.get_coefficient(orderedPair{.layer=z_idx,.face=f}) / soil_water_capacity() * K_unsaturated(f);
}

template<ElementInterface E>
double solverData<E>::top_boundary_diagonal(size_t)
{
    return 1.0 / cellFacesAndVerticalLayers.face; // TODO add to comment here what kind of BC this represents
}

template<ElementInterface E>
double solverData<E>::top_boundary_off_diagonal(size_t)
{
    // top boundary is no flux
    // off-diagonal terms require that psi_j - psi != 0
    return 0.0;
}

template<ElementInterface E>
double solverData<E>::top_boundary_rhs(size_t) const
{
    return d.psi_n.at(z_idx) / cellFacesAndVerticalLayers.face;
}
}