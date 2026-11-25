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

#include "massTransfer.hpp"
#include "unresolvedCouplingSystem.hpp"

namespace pFlow
{

namespace coupling
{

//----------------------------- constructors ----------------------------------

massTransfer::massTransfer(
    const unresolvedCouplingSystem& uCS,
    const porosity&                 prsty)
:
    porosity_(prsty)
{}

//---------------------------- public methods ---------------------------------

const Foam::dictionary& massTransfer::dict() const
{
    return massTransfer::getDict(porosity_.uCS());
}

//---------------------------- static members ---------------------------------

uniquePtr<massTransfer> massTransfer::create(
    const unresolvedCouplingSystem& uCS,
    const porosity&                 prsty)
{
    const auto& mtDict = massTransfer::getDict(uCS);

    auto shapeName = uCS.shapeTypeName();
    auto mtType    = lookupDict<Foam::word>(mtDict, "model");

    // Example: "sphereMassTransfer<Gunn>"
    Foam::word fullType = angleBracketsNames(
        shapeName + "MassTransfer",
        mtType);

    if (couplingSystemvCtorSelector_.search(fullType))
    {
        Foam::Info << "    Creating mass transfer model "
            << Green_Text(fullType) << " ...\n\n";

        return couplingSystemvCtorSelector_[fullType](uCS, prsty);
    }
    else
    {
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

bool massTransfer::found(const unresolvedCouplingSystem& uCS)
{
    return uCS.unresolvedDict()
        .subDict("reactionInteraction")
        .found("massTransferModel");
}

const Foam::dictionary& massTransfer::getDict(
    const unresolvedCouplingSystem& uCS)
{
    return uCS.unresolvedDict()
        .subDict("reactionInteraction")
        .subDict("massTransferModel");
}

//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

} // coupling

} // pFlow
