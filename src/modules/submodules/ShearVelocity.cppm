//
// Created by dallumu on 9/15/26.
//

export module ShearVelocity;

export import ShearVelocity.LiPomeroy;
export import ShearVelocity.YuPomeroy;

export namespace ShearVelocity
{
enum class Type
{
    LiPomeroy,
    YuPomeroy,
    Off,
};
}