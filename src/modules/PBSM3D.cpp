//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a
// novel modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#include "PBSM3D.hpp"

#include "StencilAssembly.hpp"

REGISTER_MODULE_CPP(PBSM3D);

struct my_fill_topo_params
{
    double a1;
    double b1;
    double a2;
    double b2;
    double m_tpi;
    double s_tpi;
};

double my_fill_topo(double x, void* p)
{
    struct my_fill_topo_params* params = (struct my_fill_topo_params*)p;
    double a1 = (params->a1);
    double b1 = (params->b1);
    double a2 = (params->a2);
    double b2 = (params->b2);
    double m_tpi = (params->m_tpi);
    double s_tpi = (params->s_tpi);

    double x2 = x + 0.25;
    if (x2 <= 0)
        return (1. - a1 * tanh(b1 * x2)) * gsl_ran_gaussian_pdf(x - m_tpi, s_tpi);
    else
        return (1. - a2 * tanh(b2 * x2)) * gsl_ran_gaussian_pdf(x - m_tpi, s_tpi);
}

struct my_fill_topo2_params
{
    double a1;
    double b1;
    double a2;
    double b2;
    double m_tpi;
    double s_tpi;
    double hs;
    double ii;
};

double my_fill_topo2(double x, void* p)
{
    struct my_fill_topo2_params* params = (struct my_fill_topo2_params*)p;
    double a1 = (params->a1);
    double b1 = (params->b1);
    double a2 = (params->a2);
    double b2 = (params->b2);
    double m_tpi = (params->m_tpi);
    double s_tpi = (params->s_tpi);
    double hs = (params->hs);
    double ii = (params->ii);

    double x2 = x + 0.25;
    if (x2 <= 0)
        return (1. - a1 * tanh(b1 * x2)) * hs / ii * gsl_ran_gaussian_pdf(x - m_tpi, s_tpi);
    else
        return (1. - a2 * tanh(b2 * x2)) * hs / ii * gsl_ran_gaussian_pdf(x - m_tpi, s_tpi);
}

struct my_fill_topo3_params
{
    double m_tpi;
    double s_tpi;
    double fac_fill;
};

double my_fill_topo3(double x, void* p)
{
    struct my_fill_topo3_params* params = (struct my_fill_topo3_params*)p;
    double m_tpi = (params->m_tpi);
    double s_tpi = (params->s_tpi);
    double fac_fill = (params->fac_fill);

    return -fac_fill * x * gsl_ran_gaussian_pdf(x - m_tpi, s_tpi);
}

enum Subgrid
{
    DoNotUse,
    SubgridV1,
    SubgridV2
};

PBSM3D::PBSM3D(config_file cfg) : module_base("PBSM3D", parallel::domain, cfg)
{
    depends("U_2m_above_srf");
    depends("vw_dir");
    depends("swe");
    depends("t");
    depends("rh");
    depends("U_R");

    provides("pbsm_more_than_avail");

    provides("global_cell_id");
    //    depends("p_snow");
    //    depends("p");

    // Determine if we want fetch, and if so, which one. Default is to use Pomeroy
    // Tanh
    // if we use Liston fetch, we don't need hours since snowfall. If we use
    // Essery 1999, then we need hours since snowfall we have to do this here so
    // we can correctly setup the depends for module linking.
    use_exp_fetch = cfg.get("use_exp_fetch", false);
    use_tanh_fetch = cfg.get("use_tanh_fetch", true);
    use_PomLi_probability = cfg.get("use_PomLi_probability", false);
    z0_ustar_coupling = cfg.get("z0_ustar_coupling", false);

    // Determine if we account for sub-grid topography impact on snow redistribution
    subgrid_topo = static_cast<Subgrid>(cfg.get<int>("use_subgrid_topo",0));

    if (use_exp_fetch && use_tanh_fetch)
    {
        CHM_THROW_EXCEPTION(module_error, "PBSM3d: Cannot specify both exp_fetch and tanh_fetch");
    };

    if (use_exp_fetch || use_tanh_fetch)
        depends("fetch");
    else
        depends("p_snow_hours");
    use_R94_lambda = cfg.get("use_R94_lambda", true);


    provides("blowingsnow_probability");
    debug_output = cfg.get("debug_output", false);

    if (debug_output)
    {
        for (int i = 0; i < sizes.vert_layers; ++i)
        {
            provides("K" + std::to_string(i));
            provides("c" + std::to_string(i));
            provides("cz" + std::to_string(i));
            provides("rm" + std::to_string(i));

            provides("csubl" + std::to_string(i));
            provides("settling_velocity" + std::to_string(i));
            provides("u_z" + std::to_string(i));

            provides_vector("uvw"+std::to_string(i));
        }

        provides("is_drifting");
        provides("Km_coeff");
        provides("Qsusp_pbsm");
        provides("height_diff");
        provides("suspension_mass");
        provides("saltation_mass");
        //        provides("Ti");

        provides("w");
        provides("hs");
        provides("ustar");
        provides("l");
        provides("z0");
        provides("lambda");

        provides("U_10m");

        provides("csalt");
        provides("csalt_orig");
        provides("csalt_reset");
        provides("mass_qsalt");
        provides("c_salt_fetch_big");

        provides("u*_th");
        provides("u*_n");
        provides("tau_n_ratio");

        provides("dm/dt");
        provides("mm");
    }
    provides("Qsubl");
    provides("Qsubl_mass");
    provides("sum_subl");

    provides("drift_mass"); // kg/m^2
    provides("Qsusp");
    provides("Qsalt");

    provides("sum_drift");

    switch (subgrid_topo)
    {
    case Subgrid::DoNotUse:
        break;
    case Subgrid::V1:
        provides("frac_contrib");
        provides("frac_contrib_nosnw");
        provides("hold_topo");
        break;
    case Subgrid::V2:
        provides("frac_contrib");
        provides("hold_topo");
        provides("test_int");
        provides("test_err");
        provides("tpi_lim");
        break;
    }
}

void PBSM3D::init(mesh& domain)
{
    sizes.local = domain->size_local_faces();
    sizes.global = domain->size_global_faces();
    sizes.vert_layers = cfg.get("sizes.vert_layers", 10);

    susp_depth = 5;                      // 5m as per pomeroy
    v_edge_height = susp_depth / sizes.vert_layers; // height of each vertical prism
    l__max = 40;                         // mixing length for diffusivity calculations

    do_fixed_settling = cfg.get("do_fixed_settling", false);

    // settling_velocity is used if the user chooses a fixed settling velocity
    // (do_fixed_settling=true)
    settling_velocity = cfg.get("settling_velocity",
            0.5); // m/s, Lehning, M., H. Löwe, M. Ryser, and N. Raderschall
    // (2008), Inhomogeneous precipitation distribution and snow
    // transport in steep terrain, Water Resour. Res., 44(7),
    // 1–19, doi:10.1029/2007WR006545.

    if (settling_velocity < 0)
    {
        CHM_THROW_EXCEPTION(module_error, "PBSM3D settling velocity must be positive");
    };

    do_sublimation = cfg.get("do_sublimation", true);
    do_lateral_diff = cfg.get("do_lateral_diff", true);
    eps = cfg.get("smooth_coeff", 820);
    min_sd_trans = cfg.get("min_sd_trans", 0.1); // m Snow holding capacity for flat terrain

    cutoff = cfg.get("cutoff",
            0.3); // cutoff veg-snow diff (m) that we inhibit saltation entirely

    snow_diffusion_const = cfg.get("snow_diffusion_const",
            0.3); // Beta * K, this is beta and scales the eddy diffusivity
    rouault_diffusion_coeff = cfg.get("rouault_diffusion_coef", false);

    enable_veg = cfg.get("enable_veg", true);

    iterative_subl = cfg.get("iterative_subl", false);

    if (rouault_diffusion_coeff)
    {
        SPDLOG_WARN( "rouault_diffusion_coef overrides const snow_diffusion_const values.");
    }

    n_non_edge_tri = 0;

    SPDLOG_DEBUG("#face={}",sizes.local);

    // **************************************************************
    // **************************************************************
    // TODO can this loop be combined with the trilinos GrsGraph creations?
    // - it's easier to follow along with separated loops, but likely less efficient
    // **************************************************************
    // **************************************************************
#pragma omp parallel for
    for (size_t i = 0; i < sizes.local; i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<data>(ID);

        if (!face->has_vegetation() && enable_veg)
        {
            SPDLOG_ERROR("Vegetation is enabled, but no vegetation parameter was found.");
        }
        if (face->has_vegetation() && enable_veg)
        {
            d.CanopyHeight = face->veg_attribute("CanopyHeight");

            // only grab LAI if we are using the R90 lambda formulation
            if (use_R94_lambda)
            {
                d.LAI = face->veg_attribute("LAI");
                d.N = 0;
                d.dv = 0;
            }
            else // use stalk formulation
            {
                d.LAI = 0;

                try{
                    d.N =  face->veg_attribute("stalk_number");
                    d.dv = face->veg_attribute("stalk_diameter");
                }
                catch(module_error& e)
                {
                    // default values that work well for scrubby alpine terrain as per Macdonald's work in the Arctic
                    d.N = 1;
                    d.dv = 0.8;
                }

            }

        }
        else
        {
            d.CanopyHeight = 0;
            d.LAI = 0;
            d.N = 0;
            d.dv = 0;
            enable_veg = false;
        }

        d.F_fill.function = &my_fill_topo;
        d.F_fill2.function = &my_fill_topo2;
        d.F_fill3.function = &my_fill_topo3;

        // pre alloc for the windpseeds
        d.u_z_susp.resize(sizes.vert_layers);

        auto& m = d.m;
        // edge unit normals
        m[0].set_size(3);
        m[0](0) = face->edge_unit_normal(0).x();
        m[0](1) = face->edge_unit_normal(0).y();
        m[0](2) = 0;

        m[1].set_size(3);
        m[1](0) = face->edge_unit_normal(1).x();
        m[1](1) = face->edge_unit_normal(1).y();
        m[1](2) = 0;

        m[2].set_size(3);
        m[2](0) = face->edge_unit_normal(2).x();
        m[2](1) = face->edge_unit_normal(2).y();
        m[2](2) = 0;

        // top
        m[3].set_size(3);
        m[3].fill(0);
        m[3](2) = 1;

        // bottom
        m[4].set_size(3);
        m[4].fill(0);
        m[4](2) = -1;

        // face areas
        for (int j = 0; j < 3; ++j)
            d.A[j] = face->edge_length(j) * v_edge_height;

        // top, bottom
        d.A[3] = d.A[4] = face->get_area();

        d.is_edge = false;
        // which faces have neighbors? Ie, are we an edge?
        for (int a = 0; a < 3; ++a)
        {
            auto neigh = face->neighbor(a);
            if (neigh == nullptr)
            {
                d.face_neigh[a] = false;
                d.is_edge = true;
            }
            else
            {
                d.face_neigh[a] = true;
            }
        }
        if (!d.is_edge)
        {
            d.cell_local_id = n_non_edge_tri;
            ++n_non_edge_tri;
        }

        d.sum_drift = 0;
        d.sum_subl = 0;
        d.csubl.resize(sizes.vert_layers);
        (*face)["sum_drift"_s]=0;

    }

    suspension_NNP.reset(new math::LinearAlgebra::NearestNeighborProblem(domain,ID,sizes.vert_layers));
    deposition_NNP.reset(new math::LinearAlgebra::NearestNeighborProblem(domain,ID));

}

struct iterHelpers
{
    uintmax_t max_iter = 500;
    double rho_p = PhysConst::rho_ice;
    static bool tol(double a, double b)  { return fabs(a - b) < 1e-8; } ;
};
struct suspensionParams
{
    double fetch = 1000.0;
    double uref = 0.0;
    double snow_depth = 0.0;
    double u2 = 0.0;
    double u10 = 0.0;
};
suspensionParams PBSM3D::get_suspension_params(mesh_elem face)
{
    suspensionParams p;
    if (use_exp_fetch || use_tanh_fetch)
        p.fetch = (*face)["fetch"_s];

    p.uref = (*face)["U_R"_s];
    p.snow_depth = (*face)["snowdepthavg"_s];
    p.snow_depth = is_nan(p.snow_depth) ? 0 : p.snow_depth;

    p.u2 = (*face)["U_2m_above_srf"_s];
    double z10; // 10-m height above the snow surface
    z10 = 10. + p.snow_depth;

    if (z10 < Atmosphere::Z_U_R)
    {
        // u10 is used by the pom probability forumuation, so don't hide behide debug output
        p.u10 = Atmosphere::log_scale_wind(p.uref, Atmosphere::Z_U_R, z10, p.snow_depth);
    }
    else
    {
        p.u10 = p.uref; // Extreme case (avalanche gone crazy case)
    }
    return p;
}
double PBSM3D::do_topo_v1(mesh_elem face, suspensionParams p)
{
    double frac_contrib = 1.0;
    double frac_contrib_nosnw = 1.0;
    double min_sd_trans_avg = min_sd_trans;
    if (!is_nan(face->parameter("TPI_neg_frac"_s))) // Fraction of negative TPI is defined
    {
        double frac_neg = face->parameter("TPI_neg_frac"_s); // fraction of the grid covered by negative TPI
        if (frac_neg > 0.)
        {
            double moy_tpi_neg = std::max(-10.0, std::min(-0.05, face->parameter("TPI_neg_mean"_s)));
            double std_tpi_neg = std::min(5.0, std::max(0.1, face->parameter("TPI_neg_std"_s)));

            //              // Derive parameters of the gamma distribution representing the distrib of TPI
            double shape_gam = pow(moy_tpi_neg, 2.0) / pow(std_tpi_neg, 2.0);
            double scale_gam = -pow(std_tpi_neg, 2.0) / moy_tpi_neg;
            //
            frac_contrib = (1. - frac_neg) + // f_TPI>0.
                           frac_neg * gsl_cdf_gamma_P(p.snow_depth, shape_gam, scale_gam);

            frac_contrib_nosnw = (1. - frac_neg);

            if (p.snow_depth > min_sd_trans)
            {

                //   double fac = shape_gam/(gsl_sf_gamma(shape_gam)*pow(scale_gam,shape_gam));
                double hint = shape_gam * scale_gam -
                              scale_gam * (p.snow_depth * gsl_ran_gamma_pdf(p.snow_depth, shape_gam, scale_gam) -
                                           min_sd_trans * gsl_ran_gamma_pdf(min_sd_trans, shape_gam, scale_gam) +
                                           shape_gam * gsl_cdf_gamma_P(min_sd_trans, shape_gam, scale_gam) +
                                           shape_gam * gsl_cdf_gamma_Q(p.snow_depth, shape_gam, scale_gam));

                min_sd_trans_avg =
                    (1. - frac_neg) * min_sd_trans +
                    frac_neg * (min_sd_trans * gsl_cdf_gamma_P(min_sd_trans, shape_gam, scale_gam) + hint +
                                p.snow_depth * gsl_cdf_gamma_Q(p.snow_depth, shape_gam, scale_gam));
            }
        }
    }
    (*face)["frac_contrib"_s] = frac_contrib;
    (*face)["frac_contrib_nosnw"_s] = frac_contrib_nosnw;
    (*face)["hold_topo"_s] = min_sd_trans_avg;
    return frac_contrib;
}
double PBSM3D::do_topo_v2(iterHelpers helpers, mesh_elem face, suspensionParams p)
{
    auto& d = face->get_module_data<data>(ID);
    double min_sd_trans_avg = min_sd_trans;
    double frac_contrib = 1.0;
    // Default values for the TPI threshold above which gullies are filled.
    double tpi_lim = -min_sd_trans;

    if (!is_nan(face->parameter("TPI_std"_s)))
    { // Std value of TPI is defined

        double moy_tpi = std::max(-5.0, std::min(5.0, face->parameter("TPI_mean"_s)));
        double std_tpi = std::min(5.0, std::max(0.1, face->parameter("TPI_std"_s)));

        if (p.snow_depth > min_sd_trans)
        {
            double fac_fill = 0.8;

            // Coefficient for the filling function that give normalized snow depth as a function of TPI

            double a1 = 1.5;
            const double b1 = 0.3;
            const double a2 = 0.6;
            double b2 = 0.55;
            if (p.snow_depth > 0.75 and p.snow_depth < 1.25)
            {
                a1 = 1.15;
            }
            else if (p.snow_depth < 1.75)
            {
                a1 = 1.1;
                b2 = 0.4;
            }
            else if (p.snow_depth < 2.25)
            {
                a1 = 0.9;
                b2 = 0.4;
            }
            else if (p.snow_depth < 2.75)
            {
                a1 = 0.85;
                b2 = 0.4;
            }
            else if (p.snow_depth < 3.25)
            {
                a1 = 0.75;
                b2 = 0.4;
            }
            else
            {
                a1 = 0.6;
                b2 = 0.35;
            }

            // Compute normalization factor
            struct my_fill_topo_params params = {a1, b1, a2, b2, moy_tpi, std_tpi};
            d.F_fill.params = &params;

            gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
            double result, error;
            int code = gsl_integration_qags(&d.F_fill, -50, 50, 0, 1e-7, 1000, w, &result, &error);
            gsl_integration_workspace_free(w);

            (*face)["test_int"_s] = result;

            // Determine TPI threshold above which gullies are considered as filled.
            auto frootFn = [&](double xx) -> double
            { return (1 - a1 * tanh(b1 * (xx + 0.25))) * p.snow_depth / result + fac_fill * xx; };
            try
            {
                auto r = boost::math::tools::bracket_and_solve_root(frootFn, -1.0, 1.0, true, iterHelpers::tol,
                                                                    helpers.max_iter);
                tpi_lim = r.first + (r.second - r.first) / 2.0;
            }
            catch (...)
            {
                // Didn't converge
            }

            // Determine area-averaged snow depth which is stored in the non-filled gullies
            struct my_fill_topo2_params params2 = {a1, b1, a2, b2, moy_tpi, std_tpi, p.snow_depth, result};
            d.F_fill2.params = &params2;

            gsl_integration_workspace* w2 = gsl_integration_workspace_alloc(1000);
            double h1;
            int code2 = gsl_integration_qags(&d.F_fill2, -50, tpi_lim, 0, 1e-7, 1000, w2, &h1, &error);
            gsl_integration_workspace_free(w2);

            // Determine area-averaged snow depth which is stored in the filled gullies
            struct my_fill_topo3_params params3 = {moy_tpi, std_tpi, fac_fill};
            d.F_fill3.params = &params3;

            gsl_integration_workspace* w3 = gsl_integration_workspace_alloc(1000);
            double h2;
            int code3 =
                gsl_integration_qags(&d.F_fill3, tpi_lim, -min_sd_trans / fac_fill, 0, 1e-7, 1000, w3, &h2, &error);
            gsl_integration_workspace_free(w3);

            // Determine area-averaged snow depth hold in the area of positive TPI
            double h3 = min_sd_trans * gsl_cdf_gaussian_Q(-min_sd_trans / fac_fill - moy_tpi, std_tpi);

            // Compute total holding capacity
            min_sd_trans_avg = std::min(h1 + h2 + h3, p.snow_depth);
        }

        // Determine fraction of the triangle that contributes to snow transport
        frac_contrib = gsl_cdf_gaussian_Q(tpi_lim - moy_tpi, std_tpi);
    }
    (*face)["tpi_lim"_s] = tpi_lim;
    (*face)["frac_contrib"_s] = frac_contrib;
    (*face)["hold_topo"_s] = min_sd_trans_avg;
    return frac_contrib;
}


template<size_t N>
static std::array<double,N> get_outward_normal_component(arma::vec (&m)[N], arma::vec uvw)
{
    std::array<double,N> udotm = {};
    for (int j = 0; j < N; ++j)
    {
        udotm[j] = arma::dot(uvw, m[j]);
    }
    return udotm;
}
struct surfaceParams
{
    double u_star_saltation_threshold;
    double hs;
    double c_salt;
    double t;
    double rh;
    double es;
    double ea;
    double T;
    double height_diff;
    double z0 = 0.0;
    double ustar = 1.3;

    void sanity_check(mesh_elem face, bool debug_output);
};
void surfaceParams::sanity_check(mesh_elem face,const bool debug_output)
{
    z0 = std::max(Snow::Z0_SNOW, z0);
    ustar = std::max(0.01, ustar);
    if (debug_output)
        (*face)["ustar"_s] = ustar;
    if (debug_output)
        (*face)["z0"_s] = z0;
}

surfaceParams PBSM3D::set_surfaceParams(const iterHelpers helpers, mesh_elem face, const suspensionParams p, const double frac_contrib)
{
    // The strategy here is as follows:
    // 0) If the exposed vegetation is above $cutoff, inhibit saltation and
    // use the classical vegheight*0.12=z0 and use that for calculating u*
    // 1) Wait for vegetation to fill up until it is within $cutoff of the
    // top of the veg
    //     then use Pomeroy & Li 2000 eqn 4 to calculate an iterative
    //     solution to u* under blowing snow conditions This effectively
    //     allows wind to blow snow out of the vegetation
    // 2) Calculate a u* that uses the blowing snow z0. Then test this
    // against the blowing snow u* threshold 3) If blowing snow isn't
    // happening, recalculate u* using the normal z0.
    auto& d = face->get_module_data<PBSM3D::data>(ID);

    surfaceParams sp{};
    double min_sd_trans_avg = min_sd_trans;
    double swe = (*face)["swe"_s];
    swe = is_nan(swe) ? 0 : swe;   // handle the first timestep where swe won't have been
    // updated if we override the module order

    // height difference between snowcover and veg
    sp.height_diff = std::max(0.0, d.CanopyHeight - p.snow_depth);
    if (!enable_veg)
        sp.height_diff = 0;
    if (debug_output)
        (*face)["height_diff"_s] = sp.height_diff;
    // This is the lamdba from Li and Pomeroy eqn 4 that is used to include
    // exposed vegetation w/ the z0 estimate
    double lambda = 0;

    d.saltation = false; // default case

    // threshold friction velocity. Compute here as it's used below as well
    // Pomeroy and Li, 2000
    // Eqn 7
    sp.T = (*face)["t"_s];
    sp.u_star_saltation_threshold = 0.35 + (1.0 / 150.0) * sp.T + (1.0 / 8200.0) * sp.T * sp.T;
    if (debug_output)
        (*face)["u*_th"_s] = sp.u_star_saltation_threshold;

    // we don't have too high of veg. Check for blowing snow
    if (sp.height_diff <= cutoff && p.snow_depth >= min_sd_trans_avg && !is_water(face))
    {

        // lambda -> 0 when height_diff ->, such as full or no veg
        if (use_R94_lambda)
            // LAI/2.0 suggestion from Raupach 1994 (DOI:10.1007/BF00709229)
            // Section 3(a)
            lambda = 0.5 * d.LAI * sp.height_diff;
        else
            lambda = d.N * d.dv * sp.height_diff; // Pomeroy formulation

        if (debug_output)
            (*face)["lambda"_s] = lambda;

        if (z0_ustar_coupling)
        {
            // Calculate the new value of z0 to take into account partially filled
            // vegetation and the momentum sink
            auto ustarFn = [&](double ustar) -> double
            {
                // Li and Pomeroy 2000, eqn 5.
                // This formulation has the following coeffs built in
                // c_2 = 1.6;
                // c_3 = 0.07519;
                // c_4 = 0.5;
                // g   = 9.81;

                return p.u2 * PhysConst::kappa / log(2.0 / (0.6131702345e-2 * ustar * ustar + .5 * lambda)) - ustar;
            };
            try
            {
                auto max_iter = helpers.max_iter;
                auto r = boost::math::tools::bracket_and_solve_root(ustarFn, 1.0, 1.0, false, iterHelpers::tol,
                                                                    max_iter);
                sp.ustar = r.first + (r.second - r.first) / 2.0;
            }
            catch (...)
            {
                // Didn't converge
                d.saltation = false;
            }
        }
        else
        {
            // follow PBSM (Pom & Li 2000; Alpine3D) and don't calculate the feedback of z0 on u*
            sp.ustar = p.u2 * PhysConst::kappa / log(2.0 / 0.0002);
        }

        if (sp.ustar >= sp.u_star_saltation_threshold)
        {
            d.saltation = true;

            if (z0_ustar_coupling)
            {
                // Update z0 for blowing snow conditions
                // Li and Pomeroy 2000, eqn 5.
                // This formulation has the following coeffs built in
                // c_2 = 1.6;
                // c_3 = 0.07519;
                // c_4 = 0.5;
                // g   = 9.81;
                sp.z0 =
                    0.6131702345e-2 * sp.ustar * sp.ustar + .5 * lambda; // pom and li 2000, eqn 4
            }
            else
            {
                sp.z0 = Snow::Z0_SNOW;
            }
        }
    }

    if (!d.saltation)
    {
        // we still need a u* for spatial K estimation later
        sp.z0 = Snow::Z0_SNOW;
        sp.ustar = std::max(0.01, PhysConst::kappa * p.uref / log(Atmosphere::Z_U_R / sp.z0));
    }

    sp.sanity_check(face,debug_output);

    sp.hs = 0;
    if (d.saltation)
        sp.hs = 0.08436 * pow(sp.ustar, 1.27); // pomeroy

    d.hs = sp.hs;
    if (debug_output)
        (*face)["hs"_s] = sp.hs;
    if (debug_output)
        (*face)["is_drifting"_s] = 0;
    if (debug_output)
        (*face)["Qsusp_pbsm"_s] = 0; // for santiy checks against pbsm

    double Qsalt = 0;
    sp.c_salt = 0;
    double T_kelvin = (*face)["t"_s] + 273.15;

    // Check if we can blow snow in this triagnle
    // Are we above saltation threshold?
    // Do we have enough mass in this triangle?
    // Has saltation been disabled because there is too much veg?
    if (d.saltation)
    {

        double rho_f = mio::Atmosphere::stdDryAirDensity(face->get_z(),
                                                         T_kelvin); // air density kg/m^3, comment in mio is wrong.1.225;

        if (debug_output)
            (*face)["blowingsnow_probability"_s] = 0; // default to 0%

        if (debug_output)
        {
            double pbsm_qsusp = pow(p.u10, 4.13) / 674100.0;
            (*face)["Qsusp_pbsm"_s] = pbsm_qsusp;
        }

        if (debug_output)
            (*face)["is_drifting"_s] = 1;

        // Pomeroy and Li 2000, eqn 8
        double Beta = 202.0; // 170.0;
        double m = 0.16;

        // tau_n_ratio = ustar_n^2 / ustar^2 from MacDonald 2009 eq 3;
        // we need (ustar_n / ustar)^2 which the original derivation gives
        // so this is is correctly squared
        double tau_n_ratio = (m * Beta * lambda) / (1.0 + m * Beta * lambda);

        if (debug_output)
            (*face)["tau_n_ratio"_s] = tau_n_ratio;

        // Pomeroy 1992, eqn 12, see note above for ustar_n calc, but ustar_n
        // is correctly squared already
        sp.c_salt = rho_f / (3.29 * sp.ustar) *
                 (1.0 - tau_n_ratio -
                  (sp.u_star_saltation_threshold * sp.u_star_saltation_threshold) / (sp.ustar * sp.ustar));

        // occasionally happens to happen at low wind speeds where the
        // parameterization breaks.
        if (sp.c_salt < 0 || std::isnan(sp.c_salt))
        {
            sp.c_salt = 0;
            d.saltation = false;
        }

        if (debug_output)
            (*face)["c_salt_fetch_big"_s] = sp.c_salt;

        // exp decay of Liston, eq 10
        // 95% of max saltation occurs at fetch = 500m
        // Liston, G., & Sturm, M. (1998). A snow-transport model for complex
        // terrain. Journal of Glaciology.
        if (use_exp_fetch && p.fetch < 500)
        {
            double fetch_ref = 500;
            double mu = 3.0;
            sp.c_salt *= 1.0 - exp(-mu * p.fetch / fetch_ref);
        }
        else if (use_tanh_fetch && p.fetch <= 300.) // use Pomeroy & Male 1986 tanh fetch
        {
            double fetch_ref = 300;
            double Lc = 0.5 * tanh(0.1333333333e-1 * fetch_ref - 2.0) + 0.5;

            sp.c_salt *= Lc;
        }

        // consider the temporal non-steady effects
        if (use_PomLi_probability) // Pomeroy and Li 2000 upscaled
                                   // probability
        {
            //    1.Essery, R., Li, L. & Pomeroy, J. A distributed model of blowing snow over complex terrain.
            //    Hydrological Processes 13, 2423–2438 (1999).
            // Probability of blowing snow
            double A = (*face)["p_snow_hours"_s];                              // hours since last snowfall
            double u_mean = 11.2 + 0.365 * sp.T + 0.00706 * sp.T * sp.T + 0.9 * log(A); // eqn 10  T -> air temp, degC
            double delta = 0.145 * sp.T + 0.00196 * sp.T * sp.T + 4.3;                  // eqn 11

            double z0v = (d.N * d.dv * sp.height_diff) / 2.0; // eqn 14
            double us = p.u10 / sqrt((1 + 340.0 * z0v));   // eqn 13

            double Pu10 = 1.0 / (1.0 + exp((sqrt(M_PI) * (u_mean - us)) / delta)); // eqn 12
            (*face)["blowingsnow_probability"_s] = Pu10;

            // decrease the saltation by the probability amount
            sp.c_salt *= Pu10;
        }

        // consider subgrid topographic effect
        if (subgrid_topo != Subgrid::DoNotUse)
        {
            sp.c_salt *= frac_contrib;
        }

        // wind speed in the saltation layer Pomeroy and Gray 1990
        double uhs = 2.8 * sp.u_star_saltation_threshold; // eqn 7

        // kg/(m*s)
        Qsalt = sp.c_salt * uhs * sp.hs; // integrate over the depth of the saltation layer, kg/(m*s)

        double mass = 0;
        double phi = (*face)["vw_dir"_s];
        Vector_2 v = -math::gis::bearing_to_cartesian(phi);

        // setup wind vector
        arma::vec uvw(3);
        uvw(0) = v.x(); // U_x
        uvw(1) = v.y(); // U_y
        uvw(2) = 0;
        double V = face->get_area();
        double udotm[3] = {0, 0, 0};
        double E[3] = {0, 0, 0};

        // compute a divergence mass flux of saltation assuming our neighbors are 0 flux
        // however, we can't compute it w/ an upwind scheme like we do for the true solution later
        // as we don't know the neighbors values (might not have been computed yet). So this is just an
        for (int j = 0; j < 3; ++j)
        {
            udotm[j] = arma::dot(uvw, d.m[j]);
            E[j] = face->edge_length(j);
            mass += -E[j] * Qsalt * udotm[j];
        }

        mass = mass / V * global_param->dt();

        if (debug_output)
        {
            (*face)["csalt_orig"_s] = sp.c_salt;
            (*face)["mass_qsalt"_s] = mass;
        }

        if (mass < 0 && std::fabs(mass) > swe)
        {
            sp.c_salt = 0;
            //-swe*V/(hs*uhs*(E[0]*udotm[0]+E[1]*udotm[1]+E[2]*udotm[2])*global_param->dt());
            // kg/(m*s)
            Qsalt = sp.c_salt * uhs * sp.hs; // integrate over the depth of the saltation layer, kg/(m*s)

            if (debug_output)
                (*face)["csalt_reset"_s] = sp.c_salt;
        }
    }

    if (debug_output)
        (*face)["csalt"_s] = sp.c_salt;

    (*face)["Qsalt"_s] = Qsalt;

    sp.rh = (*face)["rh"_s] / 100.;
    sp.es = Atmosphere::saturatedVapourPressure(sp.T);
    sp.ea = sp.rh * sp.es / 1000.;
    return sp;
}

math::LinearAlgebra::VertLayer get_layer_category(const size_t z_idx,const size_t top_idx)
{
    using namespace math::LinearAlgebra;
    if (z_idx == 0)
    {
        return VertLayer::Bottom;
    }
    else if (z_idx == top_idx)
    {
        return VertLayer::Top;
    }
    else
        return VertLayer::Middle;
};
#define SUSP_WARN "Return value must be used"
struct suspensionData
{
    mesh_elem face;
    const PBSM3D::data& d;
    const math::LinearAlgebra::Sizes& sizes;
    const surfaceParams& surf_param;
    size_t z{};
    double csubl{};
    double _edge_height{};

    // holds A_ * K_ / h_
    // _0 -> _2 are the horizontal sides
    // _3 -> is the top of the prism
    // _4 -> is the bottom the prism
    std::array<double,5> K{};
    std::array<double,5> alpha{};
    std::array<double,5> udotm{};
    size_t _idx{};
    size_t idx() const { return _idx;}
    double V{};

    [[nodiscard(SUSP_WARN)]]
    double diagonal(int f) const
    {
        return V * csubl - alpha[f];
    }

    [[nodiscard(SUSP_WARN)]]
    double off_diagonal(int f) const
    {
        return alpha[f];
    }

    [[nodiscard(SUSP_WARN)]]
    double boundary_diagonal(int f) const
    {
        return -0.1e-1 * alpha[f] - .99 * d.A[f] * udotm[f] + csubl * V;
    }

    [[nodiscard(SUSP_WARN)]]
    size_t neighbour_idx(int f) const
    {
        return sizes.global * z + face->neighbor(f)->cell_global_id;
    }

    [[nodiscard(SUSP_WARN)]]
    bool has_neighbour(int f) const
    {
        return d.face_neigh[f];
    };

    [[nodiscard(SUSP_WARN)]]
    double edge_height() const
    {
        return _edge_height;
    };

    [[nodiscard(SUSP_WARN)]]
    bool donor_on_diag(int f) const
    {
        return udotm[f] > 0;
    };

    [[nodiscard(SUSP_WARN)]]
    double donor_term(int f) const
    {
        return - d.A[f] * udotm[f];
    };

    [[nodiscard(SUSP_WARN)]]
    math::LinearAlgebra::VertLayer layer() const
    {
        return get_layer_category(z,sizes.vert_layers - 1);
    };

    [[nodiscard(SUSP_WARN)]]
    static constexpr int top_face()
    {
        return 3;
    };

    [[nodiscard(SUSP_WARN)]]
    static constexpr int bottom_face()
    {
        return 4;
    };

    double alpha_bottom(const size_t face_idx) const
    {
        double alpha_special = -d.A[face_idx] * K[face_idx] / (surf_param.hs / 2.0 + _edge_height / 2.0);
        return alpha_special;
    };

    [[nodiscard(SUSP_WARN)]]
    double bottom_boundary_diagonal() const
    {
        constexpr size_t face_idx = bottom_face();

        double alpha_special = alpha_bottom(face_idx);

        return V * csubl - d.A[face_idx] * udotm[face_idx] - alpha_special;
    };

    [[nodiscard(SUSP_WARN)]]
    double bottom_boundary_rhs() const
    {
        return -alpha_bottom(bottom_face()) * surf_param.c_salt;
    };

    [[nodiscard(SUSP_WARN)]]
    double top_boundary_rhs_in() const
    {
        return 0.0;
    };

    [[nodiscard(SUSP_WARN)]]
    double top_boundary_rhs_out() const
    {
        return 0.0;
    };

    suspensionData(mesh_elem& face, const PBSM3D::data& d, const math::LinearAlgebra::Sizes& s, const surfaceParams& sp) : face(face), d(d), sizes(s), surf_param(sp) {};
};
static_assert(math::LinearAlgebra::DonorScheme<suspensionData>, "Should satisfy DonorScheme<T>");

// Step 1 of the refactoring is to put everything in member functions,
// using structs to clean up parameters.
//
// Step 2 is to extract functions into free functions in math namespace using concepts
// to enforce member functions so that this can work with generic code.

// void PBSM3D::lateral_neighbours(math::LinearAlgebra::NearestNeighborProblem& NNP,const suspensionData& sd)
// {
//     if constexpr(DonorSchemeImpl<suspensionData>)
//     {
//         for (int f = 0; f < 3; f++)
//         {
//             if (sd.donor_on_diag()) // (sd.udotm[f] > 0)
//             {
//
//                 if (sd.d.face_neigh[f])
//                 {
//                     // TODO check these cell_global_id calls or uses of neighbour_idx for id local_id or global_id is used!!
//                     // int nidx = sd.sizes.global * sd.z + face->neighbor(f)->cell_global_id;
//
//                     // Diagonal value
//                     if constexpr (DonorSchemeImpl<suspensionData>)
//                         NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(f) - sd.d.A[f] * sd.udotm[f]);
//
//                     // Off diagonal value
//                     NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(f), sd.off_diagonal(f));
//                 }
//                 else // missing neighbor case
//                 {
//                     // no mass in
//                     //                            elements[ idx_idx_off ] += V*csubl-d.A[f]*udotm[f]-alpha[f];
//
//                     // allow mass into the domain from ghost cell
//                     NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.boundary(f));
//                 }
//             }
//             else
//             {
//                 if (sd.d.face_neigh[f])
//                 {
//                     // Diagonal entry
//                     NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(f));
//                     // Off diagonal entry
//                     NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(f), -sd.d.A[f] * sd.udotm[f] + sd.off_diagonal(f));
//                 }
//                 else
//                 {
//                     // No mass in
//                     //                            elements[ idx_idx_off ] +=  V*csubl-alpha[f];
//
//                     // allow mass in
//                     NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.boundary(f));
//                 }
//             }
//         }
//     }
// }
// void PBSM3D::donor_scheme(math::LinearAlgebra::NearestNeighborProblem& NNP,const suspensionData& sd)
// {
//     if (sd.donor_condition() && sd.)
//     {
//
//
//     }
// };
// enum class VertLayer
// {
//     Bottom,
//     Middle,
//     Top
// };
//
// double suspensionData::vertical_boundary_diagonal() const
// {
//     constexpr size_t face_idx = 4;
//
//     double alpha = - d.A[face_idx] * K[face_idx] / (surf_param.hs / 2.0 + _edge_height / 2.0);
//
//     return V * csubl - d.A[face_idx] * udotm[face_idx] - alpha;
// }
// void PBSM3D::vertical_neighbours(math::LinearAlgebra::NearestNeighborProblem& NNP, const suspensionData& sd)
// {
//     const size_t top_layer = sd.sizes.vert_layers - 1;
//     const auto& surf_param = sd.surf_param;
//
//     switch (VertLayer layer = get_layer_category(sd.z, top_layer))
//     {
//         case VertLayer::Bottom:
//         {
//             double alpha4 = sd.d.A[4] * sd.K[4] / (surf_param.hs / 2.0 + sd.edge_height() / 2.0);
//
//             // bottom face, only turbulent diffusion
//             //              elements[idx_idx_off] += V * csubl - alpha4;
//
//             // includes advection term
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx,sd.vertical_boundary_diagonal());
//             // RHS
//             // TODO generalize this and extract to be done later
//             double val = -alpha4 * surf_param.c_salt;
//             NNP.rhsSumIntoGlobalValue(sd.idx, val);
//
//             if (sd.udotm[3] > 0)
//             {
//                 // Diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(3) - d.A[3] * sd.udotm[3]);
//                 // Off diagonal
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(3), sd.off_diagonal(3));
//             }
//             else
//             {
//                 // Diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(3));
//                 // Off diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(3),
//                                                           -sd.d.A[3] * sd.udotm[3] + sd.off_diagonal(3));
//             }
//             break;
//         }
//         case VertLayer::Top:
//         {
//             //(kg/m^2/s)/(m/s)  ---->  kg/m^3
//             double cprecip = 0; //(*face)["p_snow"_s]/global_param->dt()/w;
//
//             // (*face)["p_snow"_s]=0;
//             // (*face)["p"_s]=0;
//
//             if (sd.udotm[3] > 0)
//             {
//                 // Diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(3) - sd.d.A[3] * sd.udotm[3]);
//                 // RHS
//                 const double val = -sd.alpha[3] * cprecip;
//                 NNP.rhsSumIntoGlobalValue(sd.idx, val);
//             }
//             else
//             {
//                 // Diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(3));
//                 // RHS
//                 double val = sd.d.A[3] * cprecip * sd.udotm[3] - sd.alpha[3] * cprecip;
//                 NNP.rhsSumIntoGlobalValue(sd.idx, val);
//             }
//
//             if (sd.udotm[4] > 0)
//             {
//                 // Diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(4) - sd.d.A[4] * sd.udotm[4]);
//
//                 // Off diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(4), sd.off_diagonal(4));
//             }
//             else
//             {
//                 // Diagonal entry
//                 NNP.matrixSumIntoGlobalValues(s.idx, s.idx, sd.diagonal(4));
//                 // Off diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(4),
//                                                           -sd.d.A[4] * sd.udotm[4] + sd.off_diagonal(4));
//             }
//             break;
//         }
//         case VertLayer::Middle:
//         {
//             if (sd.udotm[3] > 0)
//             {
//                 // Diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(3) - sd.d.A[3] * sd.udotm[3]);
//                 // Off diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(3), sd.off_diagonal(3));
//             }
//             else
//             {
//                 // Diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(3));
//                 // Off diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(3), -sd.d.A[3] * sd.udotm[3] + sd.diagonal(3));
//             }
//
//             if (sd.udotm[4] > 0)
//             {
//                 // Diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(4) - sd.d.A[4] * sd.udotm[4]);
//                 // Off diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(4), sd.off_diagonal(4));
//             }
//             else
//             {
//                 // Diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(4));
//                 // Off diagonal entry
//                 NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(4),
//                                                           -sd.d.A[4] * sd.udotm[4] + sd.off_diagonal(4));
//             }
//             break;
//         }
//     }
//
//
//
//
//     }
//     if (sd.z == 0)
//     {
//         double alpha4 = sd.d.A[4] * sd.K[4] / (surf_param.hs / 2.0 + v_edge_height / 2.0);
//
//         // bottom face, only turbulent diffusion
//         //              elements[idx_idx_off] += V * csubl - alpha4;
//
//         // includes advection term
//         NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx,
//                                                   sd.V * sd.csubl - sd.d.A[4] * sd.udotm[4] - alpha4);
//         // RHS
//         double val = -alpha4 * surf_param.c_salt;
//         NNP.rhsSumIntoGlobalValue(sd.idx, val);
//
//         if (sd.udotm[3] > 0)
//         {
//             // Diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(3) - d.A[3] * sd.udotm[3]);
//             // Off diagonal
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(3), sd.off_diagonal(3));
//         }
//         else
//         {
//             // Diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(3));
//             // Off diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(3),
//                                                       -sd.d.A[3] * sd.udotm[3] + sd.off_diagonal(3));
//         }
//     }
//     else if (sd.z == sizes.vert_layers - 1) // top z layer
//     {
//         //(kg/m^2/s)/(m/s)  ---->  kg/m^3
//         double cprecip = 0; //(*face)["p_snow"_s]/global_param->dt()/w;
//
//         // (*face)["p_snow"_s]=0;
//         // (*face)["p"_s]=0;
//
//         if (sd.udotm[3] > 0)
//         {
//             // Diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(3) - sd.d.A[3] * sd.udotm[3]);
//             // RHS
//             const double val = -sd.alpha[3] * cprecip;
//             NNP.rhsSumIntoGlobalValue(sd.idx, val);
//         }
//         else
//         {
//             // Diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(3));
//             // RHS
//             double val = sd.d.A[3] * cprecip * sd.udotm[3] - sd.alpha[3] * cprecip;
//             NNP.rhsSumIntoGlobalValue(sd.idx, val);
//         }
//
//         if (sd.udotm[4] > 0)
//         {
//             // Diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(4) - sd.d.A[4] * sd.udotm[4]);
//
//             // Off diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(4), sd.off_diagonal(4));
//         }
//         else
//         {
//             // Diagonal entry
//             NNP.matrixSumIntoGlobalValues(s.idx, s.idx, sd.diagonal(4));
//             // Off diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(4),
//                                                       -sd.d.A[4] * sd.udotm[4] + sd.off_diagonal(4));
//         }
//     }
//     else // middle layers
//     {
//         if (sd.udotm[3] > 0)
//         {
//             // Diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(3) - sd.d.A[3] * sd.udotm[3]);
//             // Off diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(3), sd.off_diagonal(3));
//         }
//         else
//         {
//             // Diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(3));
//             // Off diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(3), -sd.d.A[3] * sd.udotm[3] + sd.diagonal(3));
//         }
//
//         if (sd.udotm[4] > 0)
//         {
//             // Diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(4) - sd.d.A[4] * sd.udotm[4]);
//             // Off diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(4), sd.off_diagonal(4));
//         }
//         else
//         {
//             // Diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.idx, sd.diagonal(4));
//             // Off diagonal entry
//             NNP.matrixSumIntoGlobalValues(sd.idx, sd.neighbour_idx(4),
//                                                       -sd.d.A[4] * sd.udotm[4] + sd.off_diagonal(4));
//         }
//     }
// }
void PBSM3D::setup_suspension_sys(mesh& domain)
{
    // Helpers for the u* iterative solver
    // - needs to be here in thread pool, otherwise there are thread consistency
    // issues with the solver
    iterHelpers helpers;
#pragma omp for
    for (size_t i = 0; i < sizes.local; i++)
    {

        auto face = domain->face(i);

        auto& d = face->get_module_data<data>(ID);
        auto& m = d.m;

        auto p = get_suspension_params(face);

        double min_sd_trans_avg = min_sd_trans; // Grid-averaged value for the topographic subgrid holding capacity

        if (debug_output)
            (*face)["U_10m"_s] = p.u10;

        //double swe = (*face)["swe"_s]; // mm   -->    kg/m^2
        //swe = is_nan(swe) ? 0 : swe;   // handle the first timestep where swe won't have been
        // updated if we override the module order


        double frac_contrib;
        switch (subgrid_topo)
        {
        case Subgrid::DoNotUse:
            frac_contrib = 1.0; // Default value for the fraction of the grid contributing to snow transport
            break;
        case Subgrid::V1:
            // Topographic holding capacity associated with subgrid topographic features
            // This method uses the distribution of negative TPI and the mean snow depth to determine the
            // fraction of the triangle which contributes to snow transport and the associated topographic
            // holding capacity
            frac_contrib = do_topo_v1(face, p);
            break;
        case Subgrid::V2:
            // Topographic holding capacity associated with subgrid topographic features
            // This method combines the distribution of TPI with a filling function to obtain
            // an estimation of the subgrid snow depth distribution per triangle
            // A filling criteria is then used to detertmine which fraction of the triangle
            // contributes to snow transport (area of positive TPI + filled gullies)
            // and which fraction does not (gullies which are not filled).

            // Areas with negative TPI are assumed to be filled when SD = fac_fill * TPI.

            frac_contrib = do_topo_v2(helpers, face, p);
            //func here
            break;
        }

        auto surf_param = set_surfaceParams(helpers, face, p, frac_contrib);

        // iterate over the vertical layers
        for (int z = 0; z < sizes.vert_layers; ++z)
        {
            // height in the suspension layer, floats above the snow surface
            double cz = z * v_edge_height + surf_param.hs + v_edge_height / 2.; // cell center height

            // compute new U_z at this height in the suspension layer
            double u_z = 0;

            // Height above the ground (snow+free) of the suspension layer
            double hz = cz + p.snow_depth;

            // the suspension layer discretization 'floats' on top of the snow
            // surface so height_diff = d.CanopyHeight - snowdepth which is
            // looking to see if cz is within this part of the canopy
            if (d.saltation && cz < surf_param.height_diff)
            {
                // saltating so used the z0 with veg, but we are in the canopy so
                // use the saltation vel

                // wind speed in the saltation layer Pomeroy and Gray 1990
                u_z = 2.8 * surf_param.u_star_saltation_threshold; // eqn 7
            }
            else if (cz < surf_param.height_diff)
            {
                // we/re in a canopy, but not saltating, just do nothing
                u_z = 0.01; // essentially do nothing when we are in sub canopy

                //                // LAI used as attenuation coefficient introduced
                //                by Inoue (1963) and increases with canopy density
                //                double LAI =
                //                std::max(0.01,face->veg_attribute("LAI"));
                //                //bring wind down to canopy top
                //                double u_cantop = std::max(0.01,
                //                Atmosphere::log_scale_wind(suspension_params.uref,
                //                Atmosphere::Z_U_R, d.CanopyHeight, 0 , veg_params.z0));
                //
                //                u_z = Atmosphere::exp_scale_wind(u_cantop,
                //                d.CanopyHeight, cz, LAI);
            }
            else
            {
                if (hz < Atmosphere::Z_U_R)
                {
                    u_z = std::max(0.01, Atmosphere::log_scale_wind(p.uref, Atmosphere::Z_U_R, hz, p.snow_depth, surf_param.z0));
                }
                else
                {
                    u_z = std::max(0.01, p.uref);
                }
            }

            d.u_z_susp.at(z) = u_z;

            // calculate dm/dt from
            // equation 13 from Pomeroy and Li 2000
            // To do so, use equations 12 - 16 in Pomeroy et al 2000
            // Pomeroy, J. W., and L. Li (2000), Prairie and arctic areal snow
            // cover mass balance using a blowing snow model, J. Geophys. Res.,
            // 105(D21), 26619–26634, doi:10.1029/2000JD900149. [online] Available
            // from: http://www.agu.org/pubs/crossref/2000/2000JD900149.shtml

            // these are from
            // Pomeroy, J. W., D. M. Gray, and P. G. Landine (1993), The prairie
            // blowing snow model: characteristics, validation, operation, J.
            // Hydrol., 144(1–4), 165–192.

            // eqn 18, mean particle radius
            // This is 'r_r' in Pomeroy and Gray 1995, eqn 53
            double rm = 4.6e-5 * pow(cz, -0.258);
            if (debug_output)
            {
                (*face)["rm" + std::to_string(z)] = rm;
                (*face)["cz" + std::to_string(z)] = cz;
            }

            // calculate mean mass, eqn 23, 24 in Pomeroy 1993 (PBSM)
            // 52, 53 P&G 1995
            double mm_alpha = 4.08 + 12.6 * cz; // 24
            double mm = 4. / 3. * M_PI * helpers.rho_p * rm * rm * rm *
                        (1.0 + 3.0 / mm_alpha + 2. / (mm_alpha * mm_alpha));
            // mean mass, eqn 23

            // mean radius of mean mass particle
            double r_z = pow((3.0 * mm) / (4 * M_PI * helpers.rho_p), 0.3333333); // 50 in p&g 1995
            if (debug_output)
                (*face)["mm"_s] = mm;

            double xrz = 0.005 * pow(u_z, 1.36); // eqn 16

            double omega = settling_velocity;
            ; // Settling velocity
            if (!do_fixed_settling)
            {
                omega = 1.1e7 * pow(r_z, 1.8); // eqn 15 settling velocity
            }

            if (debug_output)
                (*face)["settling_velocity" + std::to_string(z)] = omega;
            double Vr = omega + 3.0 * xrz * cos(M_PI / 4.0); // eqn 14

            double v = 1.88e-5; // kinematic viscosity of air, below eqn 13 in
            // Pomeroy 1993
            double Re = 2.0 * r_z * Vr / v; // eqn  55 in p&g 1995

            double Nu, Sh;
            Nu = Sh = 1.79 + 0.606 * pow(Re, 0.5); // eqn 12

            // define above, T is in C, t is in K

            // (A.6)
            double D = 2.06e-5 * pow(surf_param.t / 273.15,
                                     1.75); // diffusivity of water vapour in air, t in K,
            // eqn A-7 in Liston 1998 or Harder 2013 A.6

            // (A.9)
            double lambda_t = 0.000063 * surf_param.t + 0.00673; //  thermal conductivity, user Harder 2013 A.9, Pomeroy's
            //  is off by an order of magnitude, this matches this
            //  https://www.engineeringtoolbox.com/air-properties-d_156.html

            // Standard constant value, e.g.,
            // https://link.springer.com/referenceworkentry/10.1007%2F978-90-481-2642-2_329
            //          double L = 2.38e6; // Latent heat of sublimation, J/kg
            double L = 2.838e6; // Latent heat of sublimation, J/kg, Corrected value

            double dmdtz = 0;

            // use Pomeroy and Li 2000 iterative sol'n for Schmidt's equation
            if (iterative_subl)
            {
                /*
                 * The *1000 and /1000 are important unit conversions. Doesn't quite
                 * match the harder paper, but Phil assures me it is correct.
                 */
                double mw = 0.01801528 * 1000.0; //[kg/mol]  ---> g/mol
                double R = 8.31441 / 1000.0;     // [J mol-1 K-1]

                double rho = (mw * surf_param.ea) / (R * surf_param.t);

                // use Harder 2013 (A.5) Formulation, but Pa formulation for e
                auto fx = [=](double Ti)
                {
                    return boost::math::make_tuple(
                        surf_param.T +
                            D * L *
                                (rho / (1000.0) -
                                 .611 * mw * exp(17.3 * Ti / (237.3 + Ti)) / (R * (Ti + 273.15) * (1000.0))) /
                                lambda_t -
                            Ti,
                        D * L *
                                (-0.6110000000e-3 * mw * (17.3 / (237.3 + Ti) - 17.3 * Ti / pow(237.3 + Ti, 2)) *
                                     exp(17.3 * Ti / (237.3 + Ti)) / (R * (Ti + 273.15)) +
                                 0.6110000000e-3 * mw * exp(17.3 * Ti / (237.3 + Ti)) / (R * pow(Ti + 273.15, 2))) /
                                lambda_t -
                            1);
                };

                double guess = surf_param.T;
                double min = -100;
                double max = 50;
                int digits = 6;

                double Ti = boost::math::tools::newton_raphson_iterate(fx, guess, min, max, digits);
                double Ts = Ti + 273.15; // dmdtz expects in K

                // now use equation 13 with our solved Ts to compute dm/dt(z)
                dmdtz = 2.0 * M_PI * rm * lambda_t / L * Nu * (Ts - (surf_param.t + 273.15)); // eqn 13 in Pomeroy and Li 2000
            }
            else // Use PBSM. Eqn 11 Pomeroy, Gray, Ladine, 1993  "﻿The
                 // prairie blowing snow model: characteristics, validation,
                 // operation"
            {
                double M = 18.01; // molecular weight of water kg kmol-1
                double R = 8313;  // universal fas constant J mol-1 K-1

                double sigma = (surf_param.rh - 1.0) * (1.019 + 0.027 * log(cz)); // undersaturation, Pomeroy and Li 2000, eqn 14

                double rho = (M * surf_param.es) / (R * surf_param.t); // saturation vapour density at t
                // radiative energy absorebed by the particle -- take from CRHM's
                // PBSM implimentation
                double Qr = 0.9 * M_PI * rm * rm * 120.0; // 120.0 = PBSM_constants::Qstar (Solar Radiation
                // Input), 0.9 comes from Schmidt (1972) assuming a
                // snow particle albedo of 0.5 and a snow surface
                // albedo of 0.8
                // rm ise used here as in Liston and Sturm (1998)

                // eqn 11 in PGL 1993, r_z is used here as in PG95 and Liston ans Sturm (1998)
                dmdtz = Sh * rho * D *
                        (6.283185308 * Nu * R * r_z * sigma * surf_param.t * surf_param.t * lambda_t - L * M * Qr + Qr * R * surf_param.t) /
                        (D * L * Sh * (L * M - R * surf_param.t) * rho + lambda_t * surf_param.t * surf_param.t * Nu * R);
            }
            if (debug_output)
                (*face)["dm/dt"_s] = dmdtz;

            if (debug_output)
                (*face)["mm"_s] = mm;

            suspensionData susp_data(face,d,sizes,surf_param);
            susp_data._edge_height = v_edge_height;
            susp_data.csubl = dmdtz / mm; // EQN 21 POMEROY 1993 (PBSM)
            susp_data._idx = sizes.global * z + face->cell_global_id;
            // eddy diffusivity (m^2/s)
            // 0,1,2 will all be K = 0, as no horizontal diffusion process
            double K[5] = {0, 0, 0, 0, 0};

            // compute alpha and K for edges
            if (do_lateral_diff)
            {
                for (int a = 0; a < 3; ++a)
                {
                    // auto neigh = face->neighbor(a);
                    susp_data.alpha[a] = d.A[a];

                    // do just very low horz diffusion for numerics
                    K[a] = 0.00001;
                    susp_data.alpha[a] *= K[a];
                }
            }
            // Li and Pomeroy 2000
            double l = PhysConst::kappa * (cz + surf_param.z0) * l__max / (PhysConst::kappa * (cz + surf_param.z0) + l__max);
            if (debug_output)
                (*face)["l"_s] = l;

            double w = omega; // settling_velocity;
            if (debug_output)
                (*face)["w"_s] = w;

            double diffusion_coeff = snow_diffusion_const; // snow_diffusion_const is a shared param so
            // need a seperate copy here we can
            // overwrite
            if (rouault_diffusion_coeff)
            {
                double c2 = 1.0;
                double dc = 1.0 / (1.0 + (c2 * w * w) / (1.56 * surf_param.ustar * surf_param.ustar));
                diffusion_coeff = dc; // nope, snow_diffusion_const is shared, use a new
            }
            if (debug_output)
                (*face)["Km_coeff"_s] = diffusion_coeff;

            // snow_diffusion_const is pretty much a calibration constant. At 1 it
            // seems to over predict transports.
            // with pomeroy fall velocity, 0.3 gives good agreement w/ published
            // Qsusp values. Low value compensates for low fall velocity
            K[3] = K[4] = diffusion_coeff * surf_param.ustar * l;

            if (debug_output)
                (*face)["K" + std::to_string(z)] = K[3];
            // top
            susp_data.alpha[3] = d.A[3] * K[3] / v_edge_height;
            // bottom
            susp_data.alpha[4] = d.A[4] * K[4] / v_edge_height;

            double phi = (*face)["vw_dir"_s]; // wind direction
            Vector_2 vwind = -math::gis::bearing_to_cartesian(phi);

            // setup wind vector
            arma::vec uvw(3);
            uvw(0) = vwind.x(); // U_x
            uvw(1) = vwind.y(); // U_y
            uvw(2) = 0;

            // above we have just the direction so it's unit vector. scale it to
            // have the same magnitude as u_z
            uvw *= u_z / arma::norm(uvw, 2);

            // now we can add in the settling_velocity
            uvw(2) = -w;

            if (debug_output)
                (*face)["u_z" + std::to_string(z)] = u_z;

            // negate as direction it's blowing instead of where it is from!!
            if (debug_output)
            {
                Vector_3 v3(-uvw(0), -uvw(1), uvw(2));
                face->set_face_vector("uvw" + std::to_string(z), v3);
            }

            susp_data.udotm = get_outward_normal_component(m, uvw);
            // lateral

            // TODO move this to a separate step, like with donor scheme
            susp_data.V = face->get_area() * v_edge_height;
            // the sink term is added on for each edge check, which isn't right
            // and ends up 5x counting it so / by 5 for V so it's
            // not 5x counted.
            susp_data.V /= 5.0;

            if (!do_sublimation)
            {
                susp_data.csubl = 0.0;
            }

            d.csubl[z] = susp_data.csubl;


            //lateral_neighbours(*suspension_NNP, susp_data);
            math::LinearAlgebra::lateral_neighbours(*suspension_NNP,susp_data);

            // vertical_neighbours(*suspension_NNP, susp_data);
            math::LinearAlgebra::vertical_neighbours(*suspension_NNP, susp_data);

        } // end z iter

    } // end face iter
}
bool PBSM3D::do_suspension_solve(mesh& domain)
{
    bool suspension_present = false;
    auto suspension_rhs_max = suspension_NNP->getRhsMax();
    if (suspension_rhs_max > suspension_present_threshold)
    {
        suspension_present = true;
    }

    if (suspension_present)
    {

        try
        {
            auto suspension_results = suspension_NNP->Solve();
            SPDLOG_DEBUG("  suspension (isolated) iterations: {} residual: {} ", suspension_results.numIters,
                         suspension_results.residual);
        }
        catch (const Belos::StatusTestError& e)
        {
            int rank = 0;
#ifdef USE_MPI
            rank = domain->_comm_world.rank();
#endif
            SPDLOG_ERROR("Rank {} suspension_rhs_max={} and suspension_present_threshold={}", rank, suspension_rhs_max,
                         suspension_present_threshold);
            //            std::string prefix = "suspension.rank" + std::to_string(rank);
            //            suspension_NNP->writeSystemMatrixMarket(prefix);
            SPDLOG_ERROR(e.what());
            suspension_present = false;
            CHM_THROW_EXCEPTION(module_error, e.what());
        }

        ////////////////////////////////////////////////////////////////////////////
        // Write solution
        //        suspension_NNP->writeSolutionMatrixMarket(prefix);
        ////////////////////////////////////////////////////////////////////////////

    } // if suspension_present fails
    else
    {
        SPDLOG_DEBUG("  No suspended snow.");
    }

    // Note we still have to do the following if there is no suspended snow.
    // - In that case suspension_solution should be the 0 vector it was initialized to for this iteration

    // auto suspension_sol_array = suspension_solution->get1dView();
    auto suspension_sol_array = suspension_NNP->getSolutionView();

#pragma omp parallel for
    for (size_t i = 0; i < sizes.local; i++)
    {
        auto face = domain->face(i);
        auto& d = face->get_module_data<data>(ID);
        double Qsusp = 0;

        double Qsubl = 0;
        for (int z = 0; z < sizes.vert_layers; ++z)
        {
            double c = suspension_sol_array[sizes.local * z + face->cell_local_id];
            c = c < 0 || is_nan(c) ? 0 : c; // harden against some numerical issues that
            // occasionally come up for unknown reasons.

            double u_z = d.u_z_susp.at(z);

            Qsusp += c * u_z * v_edge_height; /// kg/m^3 ---->  kg/(m.s)

            if (debug_output)
            {
                (*face)["c" + std::to_string(z)] = c;
                (*face)["csubl" + std::to_string(z)] = d.csubl[z];
                // This is an approximation as it uses after transport concentrations.
                // However this will have already taken into account sublimation during the coupled transport phase
                // Eqn 20 Pomeroy 1993
            }
            Qsubl += d.csubl[z] * c * v_edge_height; //  kg/(m^2 *s)=> per unit area of snowcover
        }
        (*face)["Qsusp"_s] = Qsusp;

        (*face)["Qsubl"_s] = Qsubl;
        (*face)["Qsubl_mass"_s] = Qsubl * global_param->dt(); // kg/m^2 or mm
        d.sum_subl += (*face)["Qsubl_mass"_s];
        (*face)["sum_subl"_s] = d.sum_subl;
    }
    return suspension_present;
}
void PBSM3D::setup_deposition_sys(mesh& domain)
{

    /*
       Setup and solve the linear system for deposition
    */

#pragma omp parallel for
    for (size_t i = 0; i < sizes.local; i++)
    {
        auto face = domain->face(i);
        auto& d = face->get_module_data<data>(ID);
        auto& m = d.m;

        double phi = (*face)["vw_dir"_s];
        Vector_2 v = -math::gis::bearing_to_cartesian(phi);

        // setup wind vector
        arma::vec uvw(3);
        uvw(0) = v.x(); // U_x
        uvw(1) = v.y(); // U_y
        uvw(2) = 0;

        double udotm[3] = {0, 0, 0};    // Qt is in direction u_hat
        double E[3] = {0, 0, 0};        // edge lengths b/c 2d now
        double dx[3] = {2.0, 2.0, 2.0}; // cell centre distances

        double V = face->get_area(); // V for consistency but actually an area

        int global_row, local_col, global_col;
        global_row = static_cast<int>(face->cell_global_id);

        if(is_nan(V))
        {
            SPDLOG_DEBUG("Triangle {} area is nan", face->cell_global_id);
        }
        // Diagonal element
        deposition_NNP->matrixReplaceGlobalValues(global_row, global_row, V);

        // iterate over edges
        for (int j = 0; j < 3; j++)
        {
            // just unit vectors as qsusp/qsalt flux has magnitude
            udotm[j] = arma::dot(uvw, m[j]);
            E[j] = face->edge_length(j);

            double Qtj = 0;
            double Qsj = 0;

            // Pointing same way, we advect downwind
            if (udotm[j] > 0)
            {
                Qtj = (*face)["Qsusp"_s];
                Qsj = (*face)["Qsalt"_s];

                if(is_nan(Qtj))
                {
                    SPDLOG_DEBUG("udotm >0 Qtj is nan");
                }
                if(is_nan(Qsj))
                {
                    SPDLOG_DEBUG("udotm >0 Qsj is nan");
                }
            }
            else // pointing into the wind, use upwind as donor
            {
                if (d.face_neigh[j])
                {
                    auto neigh = face->neighbor(j);

                    // if (neigh->_is_ghost) {
                    //   LOG_DEBUG << "Ghost global_id " << neigh->cell_global_id << " ghost_Qsusp: " << (*neigh)["Qsusp"_s];
                    //   LOG_DEBUG << "Ghost global_id " << neigh->cell_global_id << " ghost_Qsalt: " << (*neigh)["Qsalt"_s];
                    // }

                    Qtj = (*neigh)["Qsusp"_s];
                    Qsj = (*neigh)["Qsalt"_s];
                    if(is_nan(Qsj))
                    
                    {
                        //todo: remove this
                        Qsj = 0;
                        SPDLOG_DEBUG("udotm <0 neigh Qsj is nan");
                    }
                }
                else
                {
                    // neighbor doesn't exist, i.e., we're on an actual boundary, treat as duplicate node outside of domain
                    Qtj = (*face)["Qsusp"_s];
                    Qsj = (*face)["Qsalt"_s];

                    if(is_nan(Qsj))
                    {
                        SPDLOG_DEBUG("udotm <0 non neigh Qsj is nan");
                    }
                }
            }

            // we now have an edge Qtj & Qsj estimate
            // build up our neighbors
            if (d.face_neigh[j])
            {
                auto neigh = face->neighbor(j);
                global_col = static_cast<int>(neigh->cell_global_id);
                dx[j] = math::gis::distance(face->center(), neigh->center());

                if(is_nan(eps))
                {
                    SPDLOG_DEBUG("eps is nan!");
                }
                if(is_nan(dx[j]))
                {
                    SPDLOG_DEBUG("dx is nan!");
                }
                // diagonal entry
                deposition_NNP->matrixSumIntoGlobalValues(global_row, global_row, eps * E[j] / dx[j]);

                // off diagonal entry
                deposition_NNP->matrixSumIntoGlobalValues(global_row, global_col, -eps * E[j] / dx[j]);
            }

            // RHS
            double val = -E[j] * (Qtj + Qsj) * udotm[j];
            if (is_nan(val))
            {
                SPDLOG_DEBUG("Detected val is nan:");
                domain->print_ghost_neighbor_info();

                SPDLOG_DEBUG("\tSusp: {} salt: {}", Qtj, Qsj);

                SPDLOG_DEBUG("\ttri global id: {}", face->cell_global_id);
                (*face)["global_cell_id"_s] = face->cell_global_id;
                SPDLOG_DEBUG("\tlocal_cell_id: {}", face->cell_local_id);
                SPDLOG_DEBUG("\tis_ghost: {}", face->is_ghost);
                SPDLOG_DEBUG("\towner: {}", face->owner);
                for (int i = 0; i < 3; i++)
                {

                    SPDLOG_DEBUG("\t Neigh {} information:", i);
                    SPDLOG_DEBUG("\t\tcell_global_id: {}", face->neighbor(i)->cell_global_id);

                    SPDLOG_DEBUG("\t\tqsalt: {}", face->neighbor(i)->operator[]("Qsalt"_s));
                    SPDLOG_DEBUG("\t\tis_ghost: {}", face->neighbor(i)->is_ghost);
                    SPDLOG_DEBUG("\t\towner: {}", face->neighbor(i)->owner);
                }
                SPDLOG_DEBUG("-------------------------------------------------");
            }
            deposition_NNP->rhsSumIntoGlobalValue(global_row, val);
        }
    } // end face iteration
}
struct Present
{
    bool suspension = false;
    bool deposition = false;
};

void PBSM3D::do_deposition_solve(mesh& domain, Present present)
{
    // Check if we exceed the threshold for blowing snow
    auto deposition_rhs_max = deposition_NNP->getRhsMax();
    if (present.suspension && deposition_rhs_max > deposition_present_threshold)
    {
        present.deposition = true;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Printing out the mat/rhs
    ////////////////////////////////////////////////////////////////////////////
    // std::string deposition_file_prefix="Deposition_";
    // deposition_file_prefix += std::to_string(count);
    // deposition_NNP->writeSystemMatrixMarket(deposition_file_prefix);
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    if (present.deposition)
    {

        ////////////////////////////////////////////////////////////////////////////
        // Printing out the solution
        ////////////////////////////////////////////////////////////////////////////
        // deposition_NNP->writeSolutionMatrixMarket(deposition_file_prefix);
        // count++;
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////

        try
        {
            auto deposition_results = deposition_NNP->Solve();
            SPDLOG_DEBUG("  deposition (isolated) iterations: {} residual: {}", deposition_results.numIters, deposition_results.residual);
        } catch(Belos::StatusTestError& e)
        {
            int rank = 0;
#ifdef USE_MPI
            rank = domain->_comm_world.rank();
#endif
            SPDLOG_ERROR("Rank {} deposition_rhs_max={} and deposition_present_threshold={}", rank, deposition_rhs_max, deposition_present_threshold);

//            std::string prefix = "deposition.rank" + std::to_string(rank);
//            deposition_NNP->writeSystemMatrixMarket(prefix);
            SPDLOG_ERROR(e.what());

//            return;
            CHM_THROW_EXCEPTION(module_error, e.what());
        }

        // auto deposition_sol_array = deposition_solution->get1dView();
        auto deposition_sol_array = deposition_NNP->getSolutionView();

#pragma omp parallel for
        for (size_t i = 0; i < sizes.local; i++)
        {
            auto face = domain->face(i);
            double qdep = is_nan(deposition_sol_array[i]) ? 0 : deposition_sol_array[i];

            double mass = 0;

            //     take one FE integration step to get the total mass (SWE) that is eroded or deposited
            mass = qdep * global_param->dt(); // kg/m^2*s *dt -> kg/m^2

            // could we have eroded more mass than what exists? cap it
             double swe = (*face)["swe"_s]; // mm   -->    kg/m^2
             swe = is_nan(swe) ? 0 : swe;   // handle the first timestep where swe won't have been
            // updated if we override the module order
            if( mass < 0 && std::fabs(mass) > swe )
            {
                (*face)["pbsm_more_than_avail"_s] = 1;
                mass = -swe;
            }

            // if this triangle did not saltate, it is not a candidate for removal
            auto& d = face->get_module_data<data>(ID);
            if (mass < 0 && !d.saltation)
            {
                mass = 0;
            }

            (*face)["drift_mass"_s] = mass;
            (*face)["sum_drift"_s] += mass;
        }

    } // if present.deposition fails
    else
    {
        SPDLOG_DEBUG("  No deposited snow.");
    }
}

void PBSM3D::comm_to_ghosts(mesh& domain)
{
    /*
           Communicate necessary (neighbor) vars for deposition linear system setup
           */
    // LOG_DEBUG << "Qsusp"_s << "     " << "Qsalt"_s;
    domain->ghost_neighbors_communicate_variable("Qsusp"_s);
    domain->ghost_neighbors_communicate_variable("Qsalt"_s);
}
void PBSM3D::run(mesh& domain)
{

    SPDLOG_DEBUG("PBSM: ");

    suspension_NNP->zeroSystem();
    deposition_NNP->zeroSystem();

    // Set these flags if the RHS of the suspension system is ever nonzero
    // Thread-safe because it is only ever switched in one direction
    Present present;

#pragma omp parallel
    {
        setup_suspension_sys(domain);

    }



    present.suspension = do_suspension_solve(domain);

    comm_to_ghosts(domain);

    setup_deposition_sys(domain);

    do_deposition_solve(domain, present);

}

PBSM3D::~PBSM3D() {
}

void PBSM3D::checkpoint(mesh& domain,  netcdf& chkpt)
{
    chkpt.create_variable1D("PBSM3D:sum_drift", sizes.local);

    for (size_t i = 0; i < sizes.local; i++)
    {
        auto face = domain->face(i);
        chkpt.put_var1D("PBSM3D:sum_drift", i,
                        (*face)["sum_drift"_s]);
    }

}

void PBSM3D::load_checkpoint(mesh& domain,  netcdf& chkpt)
{
    for (size_t i = 0; i < sizes.local; i++)
    {
        auto face = domain->face(i);
        (*face)["sum_drift"_s] = chkpt.get_var1D("PBSM3D:sum_drift", i);
    }
}
