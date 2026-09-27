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

#include "heatInteraction.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "porosity.hpp"
#include "fluidAveraging.hpp"
#include "solidAveraging.hpp"

namespace pFlow
{
namespace coupling
{

//----------------------------- constructors ----------------------------------

heatInteraction::heatInteraction(
    const unresolvedCouplingSystem& uCS, 
    const porosity&                 prsty)
:
    porosity_(prsty),
    heatInteractionTimer_(
        "heatInteraction", 
        &uCS.couplingTimers())
{
    // Extract the user's chosen mapping method from the simulation dictionary.
    // No "similarToMomentum" option here - always one of "distribution"/"cell".
    auto heatExch = dict().get<Foam::word>("heatSourceExchange");

    if (heatExch == "distribution")
    {
        // Advanced mapping: Heat source is spread across multiple fluid cells
        heatExchangeDistribute_ = true;
    }
    else if (heatExch == "cell")
    {
        // Point-Center mapping: Heat source is dumped into a single fluid cell
        heatExchangeDistribute_ = false;
    }
    else
    {
        Foam::Info << "Unknown heatSourceExchange method: " << heatExch 
                   << " in " << dict().name() << Foam::endl;
        Plus::processor::abort(0);
    }

    // fluidVelocity / solidVelocity: reuse momentum's averaging when asked
    // to ("similarToMomentum"), otherwise build a dedicated one for heat.
    auto fldAvr = dict().get<Foam::word>("fluidVelocity");
    if (fldAvr != "similarToMomentum")
    {
        fluidAveraging_ = 
            fluidAveraging::create(fldAvr, uCS, "fluidVelocity_heat");
    }
    else
    {
        Foam::Info << "    heatInteraction: fluidVelocity is "
                   << "similarToMomentum - reusing momentum's averaging."
                   << Foam::endl;
    }

    auto sldVel = dict().get<Foam::word>("solidVelocity");
    if (sldVel != "similarToMomentum")
    {
        solidAveraging_ = 
            solidAveraging::create(sldVel, uCS, porosity_, "solidVelocity_heat");
    }
    else
    {
        Foam::Info << "    heatInteraction: solidVelocity is "
                   << "similarToMomentum - reusing momentum's averaging."
                   << Foam::endl;
    }

    // Instantiate the physical calculation model (delegation)
    heatTransfer_ = heatTransfer::create(uCS, porosity_);

    requireCellDistribution_ = 
        heatExchangeDistribute_ ||
        (fluidAveraging_ && fluidAveraging_->requireCellDistribution()) ||
        (solidAveraging_ && solidAveraging_->requireCellDistribution());

    // heatExchangeDistribute_ == false (dictionary value "cell") is the
    // branch that constructs and uses a PCM mapper below. This mirrors
    // momentumInteraction's (corrected) convention: "distribution" means
    // the real distributionBase is used directly, "cell" means PCM.
    if (!heatExchangeDistribute_)
    {
        noDistribution_ = makeUnique<PCM>(uCS.cMesh(), uCS.centerMass());
    }
}

//---------------------------- public methods ---------------------------------

const unresolvedCouplingSystem& heatInteraction::uCS() const 
{ 
    return porosity_.uCS(); 
}

const Foam::dictionary& heatInteraction::dict() const 
{ 
    return heatInteraction::getDict(uCS()); 
}

const Foam::dictionary& heatInteraction::getDict(
    const unresolvedCouplingSystem& uCS)
{
    return uCS.unresolvedDict().subDict("heatInteraction");
}

void heatInteraction::calculateCoupling(
    const Foam::volVectorField&     U,
    const Plus::realx3ProcCMField&  vp,
    const fluidAveraging&           fluidVelocity,
    const solidAveraging&           parVelocity,
    const Plus::realProcCMField&    dp,
    const Plus::realProcCMField&    Tp,
    Plus::realProcCMField&          Qp,
    const Plus::realProcCMField&    emissivity,
    const Plus::realProcCMField&    radSumTemp,
    const Plus::uint32ProcCMField&  radNumPrt,
    Plus::realProcCMField&          QpRad)
{
    heatInteractionTimer_.start();

    // Populates heat's own averaging objects, when it has them.
    if (fluidAveraging_)
    {
        fluidAveraging_->calculate(U);
    }

    if (solidAveraging_)
    {
        solidAveraging_->calculate(vp);
    }

    // Use heat's own averaging object when the dictionary asked for one;
    // otherwise fall back to whatever was passed in (momentum's).
    const fluidAveraging& fldVel = 
        fluidAveraging_ ? fluidAveraging_() : fluidVelocity;

    const solidAveraging& sldVel = 
        solidAveraging_ ? solidAveraging_() : parVelocity;

    // Delegate the actual Nusselt number and Stefan-Boltzmann calculations
    // to the heatTransfer_ object, passing the appropriate distribution mapper.
    heatTransfer_->calculateHeatTransfer(
        fldVel,
        sldVel,
        dp, 
        Tp, 
        heatExchangeDistribute_ ? uCS().distribution() : noDistribution_(),
        Qp,
        emissivity,
        radSumTemp,
        radNumPrt,
        QpRad);

    heatInteractionTimer_.end();

    // Print profiling information to the standard output
    Foam::Info << Blue_Text("Heat interaction time: ") 
               << Yellow_Text(heatInteractionTimer_.lastTime())
               << Yellow_Text(" s") << Foam::endl;
}

//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

} // coupling
} // pFlow


