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

#include "reactionTransfer.hpp"
#include "unresolvedCouplingSystem.hpp"

namespace pFlow::coupling
{

//----------------------------- protected methods ------------------------------

void reactionTransfer::setSuSpToZero()
{
    forAll(speciesSu_, speciesI)
    {
        auto& field = speciesSu_[speciesI];
        forAll(field, celli)
        {
            field[celli] = 0.0;
        }
    }

    forAll(speciesSp_, speciesI)
    {
        auto& field = speciesSp_[speciesI];
        forAll(field, celli)
        {
            field[celli] = 0.0;
        }
    }
}

//----------------------------- constructors ----------------------------------

reactionTransfer::reactionTransfer
(
    const unresolvedCouplingSystem& uCS,
    const porosity&                 prsty
)
:
    porosity_(prsty)
{
    const Foam::dictionary& rxnDict = getDict(uCS);

    gasSpeciesNames_ = rxnDict.get<Foam::wordList>("productGas");

    speciesSu_.setSize(gasSpeciesNames_.size());
    speciesSp_.setSize(gasSpeciesNames_.size());

    forAll(gasSpeciesNames_, i)
    {
        // Explicit mass source: [kg/(m^3.s)]
        speciesSu_.set
        (
            i,
            new Foam::volScalarField
            (
                Foam::IOobject
                (
                    "Su_" + gasSpeciesNames_[i],
                    prsty.mesh().time().timeName(),
                    prsty.mesh(),
                    Foam::IOobject::NO_READ,
                    Foam::IOobject::AUTO_WRITE
                ),
                prsty.mesh(),
                Foam::dimensionedScalar
                (
                    "SuZero",
                    Foam::dimensionSet(1, -3, -1, 0, 0, 0, 0),
                    0.0
                )
            )
        );

        // Implicit sink coefficient: [kg/(kmol.s)]
        speciesSp_.set
        (
            i,
            new Foam::volScalarField
            (
                Foam::IOobject
                (
                    "Sp_" + gasSpeciesNames_[i],
                    prsty.mesh().time().timeName(),
                    prsty.mesh(),
                    Foam::IOobject::NO_READ,
                    Foam::IOobject::AUTO_WRITE
                ),
                prsty.mesh(),
                Foam::dimensionedScalar
                (
                    "SpZero",
                    Foam::dimensionSet(1, 0, -1, 0, -1, 0, 0),
                    0.0
                )
            )
        );
    }

    Foam::Info << "    Registered Reaction Source Fields for species: "
               << gasSpeciesNames_ << Foam::endl;
}

//---------------------------- public methods ---------------------------------

const Foam::dictionary& reactionTransfer::dict() const
{
    return reactionTransfer::getDict(porosity_.uCS());
}

uniquePtr<reactionTransfer> reactionTransfer::create
(
    const unresolvedCouplingSystem& uCS,
    const porosity&                 prsty
)
{
    const auto& dict = getDict(uCS);

    auto shapeName = uCS.shapeTypeName();
    auto rtType    = lookupDict<Foam::word>(dict, "model");

    Foam::word fullType = angleBracketsNames
    (
        shapeName + "ReactionTransfer",
        rtType
    );

    if (couplingSystemvCtorSelector_.search(fullType))
    {
        Foam::Info << "    Creating reaction transfer model "
                   << Green_Text(fullType) << " ...\n";

        return couplingSystemvCtorSelector_[fullType](uCS, prsty);
    }
    else
    {
        if (Plus::processor::isMaster())
        {
            printKeys
            (
                fatalErrorInFunction
                    << "Ctor Selector " << fullType << " does not exist.\n",
                couplingSystemvCtorSelector_
            ) << '\n';
        }

        Plus::processor::abort(0);
    }

    return nullptr;
}

const Foam::dictionary& reactionTransfer::getDict
(
    const unresolvedCouplingSystem& uCS
)
{
    return uCS.unresolvedDict()
              .subDict("reactionInteraction")
              .subDict("reactionTransferModel");
}

//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

} // namespace pFlow::coupling
