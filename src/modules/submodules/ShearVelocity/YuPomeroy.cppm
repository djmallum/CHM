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
Output z0(const Input& input);
}


namespace ShearVelocity::YuPomeroy
{
Output z0(const Input& input)
{
    Output output;
    return output;
}
}