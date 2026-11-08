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

#include "reactionInteraction.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "porosity.hpp"

namespace pFlow::coupling
{

reactionInteraction::reactionInteraction
(
    const unresolvedCouplingSystem& uCS,
    const porosity&                 prsty
)
:
    porosity_(prsty),
    reactionInteractionTimer_("reactionInteraction", &uCS.couplingTimers())
{
    auto rxnExch = dict().get<Foam::word>("reactionExchange");

    if (rxnExch == "distribution")
    {
        reactionExchangeDistribute_ = true;
    }
    else if (rxnExch == "cell")
    {
        reactionExchangeDistribute_ = false;
    }
    else
    {
        Foam::Info
            << "Unknown reaction exchange method: " << rxnExch
            << " in " << dict().name() << Foam::endl;
        Plus::processor::abort(0);
    }

    reactionTransfer_ = reactionTransfer::create(uCS, porosity_);
    requireCellDistribution_ = reactionExchangeDistribute_;

    if (!reactionExchangeDistribute_)
    {
        noDistribution_ = makeUnique<PCM>(uCS.cMesh(), uCS.centerMass());
    }
}

const unresolvedCouplingSystem& reactionInteraction::uCS() const
{
    return porosity_.uCS();
}

const Foam::dictionary& reactionInteraction::dict() const
{
    return reactionInteraction::getDict(uCS());
}

const Foam::dictionary& reactionInteraction::getDict
(
    const unresolvedCouplingSystem& uCS
)
{
    return uCS.unresolvedDict().subDict("reactionInteraction");
}

void reactionInteraction::calculateCoupling
(
    const std::vector<real>&     gasMassSource,
    const std::vector<real>&     gasMassSourceSp
)
{
    reactionInteractionTimer_.start();

    reactionTransfer_->calculateReactionTransfer
    (
        gasMassSource,
        gasMassSourceSp,
        reactionExchangeDistribute_
            ? uCS().distribution()
            : noDistribution_()
    );

    reactionInteractionTimer_.end();

    Foam::Info
        << Blue_Text("Reaction interaction time: ")
        << Yellow_Text(reactionInteractionTimer_.lastTime())
        << Yellow_Text(" s") << Foam::endl;
}

} // namespace pFlow::coupling
