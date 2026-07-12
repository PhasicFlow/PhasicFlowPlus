#include "noneVirtualMass.hpp"
#include "unresolvedCouplingSystem.hpp"

pFlow::coupling::noneVirtualMass::noneVirtualMass
(
    const unresolvedCouplingSystem& uCS,
    const porosity& prsty
)
:
    virtualMass(uCS, prsty)
{
    tmpVirtualMassForce_ = Foam::tmp<Foam::volVectorField>::New
    (
        Foam::IOobject
        (
            "virtualMassForce",
            Foam::timeName(this->mesh().time()),
            this->mesh(),
            Foam::IOobject::NO_READ,
            Foam::IOobject::NO_WRITE
        ),
        this->mesh(),
        Foam::dimensionedVector
        (
            "virtualMassForce",
            Foam::dimensionSet(1,-2,-2,0,0),
            Foam::vector(0,0,0)
        )
    );
}

void pFlow::coupling::noneVirtualMass::calculateVirtualMassForce
(
    const Foam::volVectorField& U,
    const Plus::realx3ProcCMField& parAcc,
    const Plus::realProcCMField& diameter,
    Plus::realx3ProcCMField& particleForce
)
{}

