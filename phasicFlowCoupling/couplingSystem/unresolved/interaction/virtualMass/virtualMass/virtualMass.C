#include "virtualMass.hpp"
#include "unresolvedCouplingSystem.hpp"

namespace pFlow::coupling
{

virtualMass::virtualMass(const unresolvedCouplingSystem& uCS, const porosity& prsty)
:
    porosity_(prsty),
    printVirtualMass_(this->dict().getOrDefault("printVirtualMass", Foam::Switch(true)))
{}

const Foam::dictionary& virtualMass::dict()const
{
    return virtualMass::getDict(porosity_.uCS());
}

const Foam::dictionary& virtualMass::getDict(const unresolvedCouplingSystem& uCS)
{
    return uCS.unresolvedDict().subDict("momentumInteraction").subDict("virtualMass");
}

uniquePtr<virtualMass> virtualMass::create
(
    const unresolvedCouplingSystem& uCS,
    const porosity& prsty
)
{
    const auto& vmDict = virtualMass::getDict(uCS);

    auto vmType = vmDict.getOrDefault<Foam::word>("model", "none");

    if(couplingSystemvCtorSelector_.search(vmType))
    {
        Foam::Info<<"    Creating virtual mass force "<<Green_Text(vmType)<<" ...\n\n";
        return couplingSystemvCtorSelector_[vmType](uCS, prsty);
    }

    if(Plus::processor::isMaster())
    {
        printKeys
        (
            fatalErrorInFunction << "Ctor Selector "<< vmType << " does not exist"
            " for virtualMass method in "<< vmDict.name()
            <<"\nAvailable ones are: \n",
            couplingSystemvCtorSelector_
        )<<endl;
    }
    Plus::processor::abort(0);

    return nullptr;
}

} // pFlow::coupling
