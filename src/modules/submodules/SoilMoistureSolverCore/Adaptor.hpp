#pragma once

#include "Concepts.hpp"
#include "Data.hpp"
#include "Details.hpp"
#include "StencilAssembly.hpp"

namespace SoilMoistureSolver::detail
{
template<ElementInterface E>
class solverData
{
    E& face;
    data& d;
    const int z_idx;
    const Sizes sizes;
    const std::string_view ID;

    [[nodiscard]] double soil_water_capacity() const;
    double _K_unsat_lateral_face(int i, orderedPair op) const;
    double _K_unsat_vertical_Face(int i, orderedPair op) const;
    [[nodiscard]] double K_unsaturated(int) const;
public:
    solverData(data& d, E& face, int layer, const Sizes& sizes, std::string_view);

    [[nodiscard]] size_t idx() const ;
    [[nodiscard]] bool has_neighbour(int f) const;
    [[nodiscard]] size_t neighbour_idx(int f) const;

    static constexpr int top_face() { return static_cast<int>(Neighbour::Top); }
    static constexpr int bottom_face() { return static_cast<int>(Neighbour::Bottom); }

    // without tags
    using donor_choice = math::without_donor_tag;
    using boundary_donor_choice = math::without_boundary_donor_tag;
    using boundary_choice = math::without_boundary_tag;

    // Normal faces
    using rhs_choice = math::with_rhs_tag;
    [[nodiscard]] double rhs(int) const;
    [[nodiscard]] double diagonal(int f) const;
    [[nodiscard]] double off_diagonal(int f) const;

    // Boundaries
    using side_boundary_choice = math::with_side_boundary_tag;
    [[nodiscard]] static double side_boundary_diagonal(int);
    [[nodiscard]] static double side_boundary_off_diagonal(int);
    [[nodiscard]] double side_boundary_rhs(int) const;

    using bottom_boundary_choice = math::with_bottom_boundary_tag;
    [[nodiscard]] static double bottom_boundary_diagonal(int);
    [[nodiscard]] static double bottom_boundary_off_diagonal(int);
    [[nodiscard]] double bottom_boundary_rhs(int) const;

    using top_boundary_choice = math::with_top_boundary_tag;
    [[nodiscard]] static double top_boundary_diagonal(int);
    [[nodiscard]] static double top_boundary_off_diagonal(int);
    [[nodiscard]] double top_boundary_rhs(int) const;

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
template <ElementInterface E> double solverData<E>::_K_unsat_lateral_face(int i, const orderedPair op) const
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
template <ElementInterface E> double solverData<E>::_K_unsat_vertical_Face(const int i, const orderedPair op) const
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
double solverData<E>::K_unsaturated(int i) const
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
solverData<E>::solverData(data& d, E& face, const int layer, const Sizes& sizes, const std::string_view id)
    : face(face), d(d), z_idx(layer), sizes(sizes), ID(id)
{
}


template<ElementInterface E>
size_t solverData<E>::idx() const { return sizes.global * z_idx + face->cell_global_id; }

template<ElementInterface E>
bool solverData<E>::has_neighbour(const int f) const { return d.cell_info.neighbour_idx_[f][z_idx].has_value(); }

template<ElementInterface E>
size_t solverData<E>::neighbour_idx(const int f) const
{
    const auto i = d.cell_info.get_neighbour_idx(orderedPair{.layer=z_idx,.face=static_cast<size_t>(f)});
    if (!i)
    {
        const std::string err = std::format("neighbour_idx invoked for a face without a neighbour at index {} of face {}",
            f,face->cell_global_id);
        CHM_THROW_EXCEPTION(module_error,err);
    }
    return *i;
}

template<ElementInterface E>
double solverData<E>::rhs(int f) const
{
    // const auto elevation = d.get_elevation(f);
    const auto* face_geometry = std::get_if<Interior>(&d.cell_info.cell_geometry[f][z_idx]);
    if (!face_geometry)
    {
        const std::string err = std::format("Boundary face detected in non-boundary function. At triangle {}, face {}, and layer {}",
            this->face->cell_global_id,f,z_idx);
        CHM_THROW_EXCEPTION(module_error,err);
    }

    // TODO neighbour.z - z = the distance between the neighbour centre and the centre of the current cell.
    // Therefore, when one defines alpha for this term, as written it must include the distance between
    return d.psi_n.at(z_idx) / cellFacesAndVerticalLayers.face + d.cell_info.coefficient.at(f).at(z_idx) / soil_water_capacity() *
        K_unsaturated(f) * (face_geometry->geometry.elevation.neighbour - face_geometry->geometry.elevation.owner);
}

template<ElementInterface E>
double solverData<E>::diagonal(const int f) const
{
    return 1.0 / cellFacesAndVerticalLayers.face + d.cell_info.coefficient.at(f).at(z_idx) / soil_water_capacity() * K_unsaturated(f); /* TODO everything about geometry or constant
                                                                      *  in time goes in alpha, could make it a type
                                                                      */
}

template<ElementInterface E>
double solverData<E>::off_diagonal(const int f) const
{
    return -d.cell_info.coefficient.at(f).at(z_idx) / soil_water_capacity() * K_unsaturated(f); // TODO see diagonal
    // TODO off diagonal contributions
}

template<ElementInterface E>
double solverData<E>::side_boundary_diagonal(int) {
    return 1.0 / cellFacesAndVerticalLayers.face;
}

template<ElementInterface E>
double solverData<E>::side_boundary_off_diagonal(const int)
{
    return 0.0;
}

template<ElementInterface E>
double solverData<E>::side_boundary_rhs(const int f) const
{
    if (f > 2)
        CHM_THROW_EXCEPTION(module_error,"side boundary must only be faces 0, 1, or 2");

    const auto* face_geo = std::get_if<Boundary>(&d.cell_info.cell_geometry[z_idx][f]);
    if (!face_geo)
    {
        const std::string err = std::format("Interior face detected in a boundary function. At triangle {}, face {}, and layer {}",
            this->face->cell_global_id,f,z_idx);
        CHM_THROW_EXCEPTION(module_error,err);
    }
    const auto elevation_change = face_geo->DeltaZ;

    return d.psi_n.at(z_idx) / cellFacesAndVerticalLayers.face + d.cell_info.coefficient.at(f).at(z_idx) / soil_water_capacity() * K_unsaturated(f) * elevation_change;
}

template<ElementInterface E>
double solverData<E>::bottom_boundary_diagonal([[maybe_unused]] int) {
    return 1.0 / cellFacesAndVerticalLayers.face;
}

template<ElementInterface E>
double solverData<E>::bottom_boundary_off_diagonal([[maybe_unused]] int)
{
    return 0.0;
}

template<ElementInterface E>
double solverData<E>::bottom_boundary_rhs(const int f) const
{
    // WARNING: alpha for the bottom boundary must not include the distance to the neighbour
    // cell centre. Impossible to obtain since it doesn't exist. But its worth being aware.
    return d.psi_n.at(z_idx) / cellFacesAndVerticalLayers.face + d.cell_info.coefficient.at(f).at(z_idx) / soil_water_capacity() * K_unsaturated(f);
}

template<ElementInterface E>
double solverData<E>::top_boundary_diagonal(int)
{
    return 1.0 / cellFacesAndVerticalLayers.face; // TODO add to comment here what kind of BC this represents
}

template<ElementInterface E>
double solverData<E>::top_boundary_off_diagonal(int)
{
    // top boundary is no flux
    // off-diagonal terms require that psi_j - psi != 0
    return 0.0;
}

template<ElementInterface E>
double solverData<E>::top_boundary_rhs(int) const
{
    return d.psi_n.at(z_idx) / cellFacesAndVerticalLayers.face;
}
}