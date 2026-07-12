
#include "surfTorqueNone.hpp"

pFlow::coupling::surfTorqueNone::surfTorqueNone
(
    const unresolvedCouplingSystem &uCS, 
    const porosity &prsty
)
:
    surfaceRotationTorque(uCS, prsty)
{

}


void pFlow::coupling::surfTorqueNone::calculateSurfaceTorque
(
    const Foam::volVectorField&     U,
    const Plus::realx3ProcCMField&  parVel,
    const Plus::realx3ProcCMField&  parRotVel,
    const Plus::realProcCMField&    diameter,
    Plus::realx3ProcCMField&        particleTorque
)
{

}