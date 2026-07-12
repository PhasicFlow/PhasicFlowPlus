
#include "lowReynolds.hpp"
#include "lift.hpp"

pFlow::coupling::lowReynolds::lowReynolds
(
    const unresolvedCouplingSystem &uCS, 
    const porosity &prsty
)
:
    surfaceRotationTorque(uCS, prsty)
{

}


void pFlow::coupling::lowReynolds::calculateSurfaceTorque
(
    const Foam::volVectorField&     U,
    const Plus::realx3ProcCMField&  parVel,
    const Plus::realx3ProcCMField&  parRotVel,
    const Plus::realProcCMField&    diameter,
    Plus::realx3ProcCMField&        particleTorque
)
{

    const size_t nPar = diameter.size();
    const auto& parCellInd =  this->parCellIndex();
    const auto& nu = this->mesh().template lookupObject<Foam::volScalarField>("nu");
    const auto& rho = this->mesh().template lookupObject<Foam::volScalarField>("rho");
    
    auto curlUPtr = lift::fluidVorticity(U);
    const auto& curlU = curlUPtr.ref();

    #pragma omp parallel for schedule(dynamic)
    for(size_t i=0; i<nPar; ++i)
    {
        const Foam::label cellI = parCellInd[i];
        
        if(cellI == -1 )continue;

        const Foam::vector up{
            parVel[i].x(), 
            parVel[i].y(), 
            parVel[i].z()};

        const Foam::vector wp{
            parRotVel[i].x(),
            parRotVel[i].y(),
            parRotVel[i].z()
        };

        const Foam::scalar dp = diameter[i];
        
        const Foam::scalar f_spin = 1.0;
        const Foam::scalar f_shear = 1.0;
            
        Foam::vector torque = - Foam::constant::mathematical::pi * rho[cellI] * nu[cellI]*
            dp*dp*dp*(f_spin*wp - 0.5*f_shear* curlU[cellI]);

        particleTorque[i] += realx3(torque.x(), torque.y(), torque.z());

       
    }

}