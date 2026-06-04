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

    double soil_water_capacity() const;
    double K_unsaturated(int) const;
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

template<ElementInterface E>
double solverData<E>::K_unsaturated(int i) const
{
    /* Campbell (1974)
     * Also see: Deb and Shukla (2012) for long list
     */
    // TODO K_unsaturated should also handle vertical neighbours too!!
    const auto neigh = face->neighbor(i);
    const auto d_neigh = neigh->template get_module_data<data>(ID.data());

    // TODO testing this equation for accuracy AND behaviour near saturation and dry soil
    // Source is Campbell (1974)
    // Good source is also Deb and Shukla (2012)
    // TODO i is face number not layer!!
    auto campbell_eqn = [layer = this->z_idx](const data& d)
    {
        return d.K_saturated * std::pow(d.air_entry_tension / d.psi_n.at(layer),2+3/d.pore_size_dist_index);
    };

    Pair pair;
    pair.owner = campbell_eqn(d);
    pair.neighbour = campbell_eqn(d_neigh);

    return d.cell_info.interp_to_face[z_idx][i]->interp(pair);
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
    const auto i = d.cell_info.neighbour_idx_[f][z_idx];
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
    return d.psi_n.at(z_idx) / NUM_NEIGHBOURS + d.cell_info.alpha.at(f).at(z_idx) / soil_water_capacity() *
        K_unsaturated(f) * (face_geometry->geometry.elevation.neighbour - face_geometry->geometry.elevation.owner);
}

template<ElementInterface E>
double solverData<E>::diagonal(const int f) const
{
    return 1.0 / NUM_NEIGHBOURS + d.cell_info.alpha.at(f).at(z_idx) / soil_water_capacity() * K_unsaturated(f); /* TODO everything about geometry or constant
                                                                      *  in time goes in alpha, could make it a type
                                                                      */
}

template<ElementInterface E>
double solverData<E>::off_diagonal(const int f) const
{
    return -d.cell_info.alpha.at(f).at(z_idx) / soil_water_capacity() * K_unsaturated(f); // TODO see diagonal
    // TODO off diagonal contributions
}

template<ElementInterface E>
double solverData<E>::side_boundary_diagonal(int) {
    return 1.0 / NUM_NEIGHBOURS;
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

    return d.psi_n.at(z_idx) / NUM_NEIGHBOURS + d.cell_info.alpha.at(f).at(z_idx) / soil_water_capacity() * K_unsaturated(f) * elevation_change;
}

template<ElementInterface E>
double solverData<E>::bottom_boundary_diagonal([[maybe_unused]] int) {
    return 1.0 / NUM_NEIGHBOURS;
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
    return d.psi_n.at(z_idx) / NUM_NEIGHBOURS + d.cell_info.alpha.at(f).at(z_idx) / soil_water_capacity() * K_unsaturated(f);
}

template<ElementInterface E>
double solverData<E>::top_boundary_diagonal(int)
{
    return 1.0 / NUM_NEIGHBOURS; // TODO add to comment here what kind of BC this represents
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
    return d.psi_n.at(z_idx) / NUM_NEIGHBOURS;
}
}