#ifndef __noneVirtualMass_hpp__
#define __noneVirtualMass_hpp__

#include "virtualMass.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class noneVirtualMass : public virtualMass
{
private:

    tmp<Foam::volVectorField> tmpVirtualMassForce_;

public:

    TypeInfo("none");

    noneVirtualMass(const unresolvedCouplingSystem& uCS, const porosity& prsty);

    virtual ~noneVirtualMass() = default;

    add_vCtor
    (
        virtualMass,
        noneVirtualMass,
        couplingSystem
    );

    void calculateVirtualMassForce
    (
        const Foam::volVectorField& U,
        const Plus::realx3ProcCMField& parAcc,
        const Plus::realProcCMField& diameter,
        Plus::realx3ProcCMField& particleForce
    ) override;

    Foam::tmp<Foam::volVectorField> virtualMassForce()const override
    {
        return tmpVirtualMassForce_;
    }
};

} // pFlow::coupling


#endif
