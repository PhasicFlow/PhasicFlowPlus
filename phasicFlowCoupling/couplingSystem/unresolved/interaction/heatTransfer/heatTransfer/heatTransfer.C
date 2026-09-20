/*------------------------------- phasicFlow ---------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow
------------------------------------------------------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
------------------------------------------------------------------------------
Licence:
  This file is part of phasicFlow code. It is a free software for simulating
  granular and multiphase flows. You can redistribute it and/or modify it under
  the terms of GNU General Public License v3 or any other later versions.

  phasicFlow is distributed to help others in their research in the field of
  granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

-----------------------------------------------------------------------------*/

#include "heatTransfer.hpp"
#include "unresolvedCouplingSystem.hpp"

namespace pFlow
{
namespace coupling
{

//----------------------------- protected methods -----------------------------

void heatTransfer::setSuSpToZero()
{
    // Resetting source terms to zero prevents indefinite accumulation 
    // from previous time steps. Note: In OpenFOAM, the assignment 
    // operator '=' is strictly required here; the '==' operator creates 
    // an implicit fvMatrix equation and does not assign values.
    heatSu_ = Foam::dimensionedScalar("zero", heatSu_.dimensions(), 0.0);
    heatSp_ = Foam::dimensionedScalar("zero", heatSp_.dimensions(), 0.0);
}

//----------------------------- constructors ----------------------------------

heatTransfer::heatTransfer(
    const unresolvedCouplingSystem& uCS, 
    const porosity&                 prsty)
: 
    heatSu_(
        Foam::IOobject(
            "heatSu",
            prsty.mesh().time().timeName(),
            prsty.mesh(),
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::AUTO_WRITE),
        prsty.mesh(),
        // Dimension: Mass=1, Length=-1, Time=-3 -> [kg/(m.s^3)] = [W/m^3]
        Foam::dimensionedScalar(
            "heatSu",
            Foam::dimensionSet(1, -1, -3, 0, 0, 0, 0), 
            0.0)),
    heatSp_(
        Foam::IOobject(
            "heatSp",
            prsty.mesh().time().timeName(),
            prsty.mesh(),
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::AUTO_WRITE),
        prsty.mesh(),
        // Dimension: Mass=1, Length=-1, Time=-3, Temp=-1 -> [W/(m^3.K)]
        Foam::dimensionedScalar(
            "heatSp", 
            Foam::dimensionSet(1, -1, -3, -1, 0, 0, 0), 
            0.0)),
    porosity_(prsty) 
{
    // Initialization fully handled via initializer list
}

//---------------------------- public methods ---------------------------------

const Foam::dictionary& heatTransfer::getDict(
    const unresolvedCouplingSystem& uCS)
{
    return uCS.unresolvedDict()
              .subDict("heatInteraction")
              .subDict("heatTransferModel");
}

const Foam::dictionary& heatTransfer::dict() const
{
    return heatTransfer::getDict(porosity_.uCS());
}

uniquePtr<heatTransfer> heatTransfer::create(
    const unresolvedCouplingSystem& uCS, 
    const porosity&                 prsty)
{
    const auto& htDict = heatTransfer::getDict(uCS);
    
    // Extract base geometry type (e.g., "sphere") and the model (e.g., "Gunn")
    auto shapeName = uCS.shapeTypeName();
    auto htType    = lookupDict<Foam::word>(htDict, "model");
    
    // Construct the fully qualified typename expected in the registry 
    // Example: "<sphereHeatTransfer>Gunn"
    Foam::word fullType = angleBracketsNames(
        shapeName + "HeatTransfer", 
        htType);
    
    // Search the factory registry and instantiate if found
    if (couplingSystemvCtorSelector_.search(fullType))
    {
        Foam::Info << "    Creating heat transfer model " 
                   << Green_Text(fullType) << " ...\n\n";
                   
        return couplingSystemvCtorSelector_[fullType](uCS, prsty);
    }
    else
    {
        // Graceful failure: List available models to the user
        if (Plus::processor::isMaster())
        {
            printKeys( 
                fatalErrorInFunction 
                    << "Ctor Selector " << fullType << " does not exist.\n",
                couplingSystemvCtorSelector_) << '\n';
        }
        
        Plus::processor::abort(0);
    }
    
    return nullptr;
}

//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

} // coupling
} // pFlow



