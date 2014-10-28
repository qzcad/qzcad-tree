#include "femforce.h"

FEMForce::FEMForce()
{
    forceType_ = NODAL_FORCE;
}

FEMForce::FEMForce(ForceType forceType)
{
    forceType_ = forceType;
}

FEMForce::FEMForce(const FEMForce &force)
{
    forceType_ = force.forceType_;
}
ForceType FEMForce::forceType() const
{
    return forceType_;
}

void FEMForce::setForceType(const ForceType &forceType)
{
    forceType_ = forceType;
}

