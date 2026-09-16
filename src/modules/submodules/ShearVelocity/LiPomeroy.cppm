//
// Created by dallumu on 9/15/26.
//
module;
#include "physics/PhysConst.h"
#include <boost/math/tools/roots.hpp>
export module ShearVelocity.LiPomeroy;

export namespace ShearVelocity::LiPomeroy
{
struct Input
{
    double lambda;
    double u2;
    uintmax_t max_iter;
};

struct Output
{
    double ustar{};
    double z0{};
    bool saltation = true;
};

Output friction_velocity(Input& input)
{
    // Calculate the new value of z0 to take into account partially filled
    // vegetation and the momentum sink
    auto tol = [](const double a, const double b) -> bool { return fabs(a - b) < 1e-8; };
    Output output;
    auto ustarFn = [&](const double ustar) -> double
    { return input.u2 * PhysConst::kappa / log(2.0 / z0(ustar,input.lambda)) - ustar; };

    try
    {
        const auto r = boost::math::tools::bracket_and_solve_root(ustarFn, 1.0, 1.0, false, tol, input.max_iter);
        output.ustar = r.first + (r.second - r.first) / 2.0;
    }
    catch (...)
    {
        // Didn't converge
        output.saltation = false;
    }
    return output;
}
double z0(const double ustar, const double lambda)
{
    return 0.6131702345e-2 * ustar * ustar + .5 * lambda;
}
} // namespace ShearVelocity::LiPomeroy