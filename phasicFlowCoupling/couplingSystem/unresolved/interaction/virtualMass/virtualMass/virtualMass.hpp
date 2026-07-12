#ifndef __virtualMass_hpp__
#define __virtualMass_hpp__

#include "OFCompatibleHeader.hpp"
#include "virtualConstructor.hpp"
#include "procCMFields.hpp"
#include "porosity.hpp"

namespace pFlow::coupling
{

class unresolvedCouplingSystem;

class virtualMass
{
private:
    
    /// Reference to porosity object for momentum coupling
    const porosity&         porosity_;
    
    /// Flag to enable/disable virtual mass force field to the output file
    Foam::Switch         printVirtualMass_;

protected:

    bool printVirtualMass()const
    {
        return printVirtualMass_;
    }

public:

    TypeInfo("virtualMass");

    virtualMass(const unresolvedCouplingSystem& uCS, const porosity& prsty);

    virtual ~virtualMass() = default;

    create_vCtor
    (
        virtualMass,
        couplingSystem,
        (
            const unresolvedCouplingSystem& uCS,
            const porosity& prsty
        ),
        (uCS, prsty)
    );

    virtual
    void calculateVirtualMassForce
    (
        const Foam::volVectorField& U,
        const Plus::realx3ProcCMField& parAcc,
        const Plus::realProcCMField& diameter,
        Plus::realx3ProcCMField& particleForce
    ) = 0;

    virtual
    Foam::tmp<Foam::volVectorField> virtualMassForce()const = 0;

    const porosity& Porosity()const
    {
        return porosity_;
    }

    const auto& parCellIndex()const
    {
        return porosity_.parCellIndex();
    }

    const Foam::fvMesh& mesh()const
    {
        return porosity_.mesh();
    }

    const Foam::volScalarField& alpha()const
    {
        return static_cast<const Foam::volScalarField&>(porosity_);
    }

    const Foam::dictionary& dict()const;

    static
    const Foam::dictionary& getDict(const unresolvedCouplingSystem& uCS);

    static
    uniquePtr<virtualMass> create
    (
        const unresolvedCouplingSystem& uCS,
        const porosity& prsty
    );
};

} // pFlow::coupling

#endif

