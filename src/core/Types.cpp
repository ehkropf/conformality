#include "Types.h"

#include "../domains/Domain.h"

MappingType determineMappingType(const Domain& source, const Domain& target)
{
    bool sourceHasInfinity = source.isUnbounded();
    bool targetHasInfinity = target.isUnbounded();

    if (!sourceHasInfinity && !targetHasInfinity)
    {
        return MappingType::INTERIOR_TO_INTERIOR;
    }
    else if (sourceHasInfinity && !targetHasInfinity)
    {
        return MappingType::EXTERIOR_TO_INTERIOR;
    }
    else if (!sourceHasInfinity && targetHasInfinity)
    {
        return MappingType::INTERIOR_TO_EXTERIOR;
    }
    else
    {
        return MappingType::EXTERIOR_TO_EXTERIOR;
    }
}