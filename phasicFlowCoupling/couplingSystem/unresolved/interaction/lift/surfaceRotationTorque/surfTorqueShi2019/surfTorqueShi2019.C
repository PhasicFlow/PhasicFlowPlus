
#include "surfTorqueShi2019.hpp"
#include "lift.hpp"

pFlow::coupling::surfTorqueShi2019::surfTorqueShi2019
(
    const unresolvedCouplingSystem &uCS, 
    const porosity &prsty
)
:
    surfaceRotationTorque(uCS, prsty),
    residualRe_(this->dict().get<Foam::scalar>("residualRe"))
{

}


void pFlow::coupling::surfTorqueShi2019::calculateSurfaceTorque
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

        const Foam::vector uf = U[cellI];
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

        const Foam::vector uRelVec = up - uf;
        const Foam::scalar uRel = Foam::mag(uRelVec);

        const Foam::scalar Rep = Foam::max(uRel*dp/nu[cellI], residualRe_);
        const Foam::scalar Rew = Foam::max(Foam::mag(wp)*dp*dp/nu[cellI], Foam::SMALL);
        const Foam::scalar Res = Foam::max(Foam::mag(curlU[cellI])*dp*dp/nu[cellI], Foam::SMALL);
       
        
        const Foam::scalar f_spin = 1.0 + 5.0/(64.0*Foam::constant::mathematical::pi)*Foam::pow(Rew,0.6);
        const Foam::scalar f_shear = f_spin * 
            (
                (1 + 0.4 * (Foam::exp(-0.0135*Res) - 1)) * 
                (1-0.0702*Foam::pow(Rep, 0.455))
            );
            
        Foam::vector torque = - Foam::constant::mathematical::pi * rho[cellI] * nu[cellI]*
            dp*dp*dp*(f_spin*wp - 0.5*f_shear* curlU[cellI]);

        particleTorque[i] += realx3(torque.x(), torque.y(), torque.z());

       
    }

}