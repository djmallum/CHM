//
// Created by dallumu on 9/15/26.
//

export module ShearVelocity.YuPomeroy;

export namespace ShearVelocity::YuPomeroy
{
struct Input {};
struct Output
{
    double ustar{};
};
Output friction_velocity(const Input& input);
double z0();
}


namespace ShearVelocity::YuPomeroy
{
Output friction_velocity(const Input& input)
{
    Output output;
    return output;
}
}