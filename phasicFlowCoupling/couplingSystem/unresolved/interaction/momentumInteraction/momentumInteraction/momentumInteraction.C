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

#include "momentumInteraction.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "porosity.hpp"

namespace pFlow::coupling
{


momentumInteraction::momentumInteraction
(
    const unresolvedCouplingSystem &uCS, 
    const porosity &prsty
)
:
    porosity_(prsty),
    momentumInteractionTimer_
    (
        "momentumInteraction", 
        &uCS.couplingTimers()
    )
{

    auto momExch = dict().get<Foam::word>("momentumExchange");

    if(momExch == "distribution")
    {
        momentumExchangeDistribute_ = true;
    }
    else if(momExch == "cell")
    {
        momentumExchangeDistribute_ = false;
    }
    else
    {
        Foam::Info
            << "Unknown momentum exchange method: " << momExch 
            << " in "<< dict().name() << Foam::endl;
        Plus::processor::abort(0);
    }

    drag_ = drag::create(uCS, porosity_);

    lift_ = lift::create(uCS, porosity_);

    auto fluidAveragingType = dict().template get<Foam::word>("fluidVelocity");

    fluidAveraging_ = fluidAveraging::create(fluidAveragingType, uCS, "fluidVelocity");

    auto solidAveragingType = dict().template get<Foam::word>("solidVelocity");
    solidAveraging_ = solidAveraging::create(solidAveragingType, uCS, porosity_, "Us");

    requireCellDistribution_ = momentumExchangeDistribute_ ||
        fluidAveraging_ -> requireCellDistribution() ||
        solidAveraging_ -> requireCellDistribution();

    if(requireCellDistribution_)
    {
        INFORMATION<<"Cell distribution is active."<<END_INFO;
    }
    
    if(momentumExchangeDistribute_)
    {
        noDistribution_ = makeUnique<PCM>(uCS.cMesh(), uCS.centerMass());
    }

}

const unresolvedCouplingSystem& momentumInteraction::uCS()const
{
    return porosity_.uCS();
}


const Foam::dictionary& momentumInteraction::dict()const
{
    return  momentumInteraction::getDict(uCS());

}

void momentumInteraction::calculateCoupling
(
    const Foam::volVectorField&     U,
    const Plus::realx3ProcCMField&  vp,
    const Plus::realx3ProcCMField&  wp,
    Plus::realx3ProcCMField&        fluidForce,
    Plus::realx3ProcCMField&        fluidTorque
)
{

    momentumInteractionTimer_.start();

    // calculate fluid averaging
    fluidAveraging_->calculate(U);

    // calculate solid averaging
    solidAveraging_->calculate(vp);

    // calculate drag force
    drag_->calculateDragForce
    (
        *fluidAveraging_,
        *solidAveraging_,
        porosity_.particleDiameter(),
        momentumExchangeDistribute_ ? noDistribution_() : uCS().distribution(),
        fluidForce);

    lift_->calculateLiftForceTorque(
        U,
        vp,
        wp,
        porosity_.particleDiameter(),
        fluidForce,
        fluidTorque);

    momentumInteractionTimer_.end();

    Foam::Info<<Blue_Text("Momentum interaction time: ")<< 
                Yellow_Text(momentumInteractionTimer_.lastTime())<<
                Yellow_Text(" s")<<Foam::endl;
}


const Foam::dictionary& momentumInteraction::getDict(const unresolvedCouplingSystem& uCS)
{
    return uCS.unresolvedDict().subDict("momentumInteraction");	
}


} // pFlow::coupling