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

#include "momentumSphereUnresolvedCouplingSystem.hpp"

pFlow::coupling::momentumSphereUnresolvedCouplingSystem::momentumSphereUnresolvedCouplingSystem
(
    word shapeTypeName,
    word couplingSystemType, 
    Foam::fvMesh& mesh,
    int argc, 
    char* argv[]
)
:
    unresolvedCouplingSystem
    (
        shapeTypeName, 
        couplingSystemType, 
        mesh, 
        argc, 
        argv
    ),
    porosity_
    (
        porosity::create
        (
            *this, 
            this->cMesh(),
            this->particleDiameter()
        )
    ),
    momentumInteraction_
    (
        *this, 
        porosity_()
    ),
    porosityTimer_
    (
        "porosity", 
        &this->couplingTimers()
    )
{
    requiresDistribution_ = 
        porosity_().requireCellDistribution()|| momentumInteraction_.requireCellDistribution();
}


void pFlow::coupling::momentumSphereUnresolvedCouplingSystem::calculatePorosity()
{
    // update coupling mesh and map particles 
    this->cMesh().update();

    porosityTimer_.start();

    // update weights for distribution (if coupling requires it)
    if(requiresDistribution_)
        this->updateDistributionWeights();

    // calculate porosity 
    porosity_->calculatePorosity();

    porosityTimer_.end();

    Foam::Info<<Blue_Text("Porosity time: ")<< 
                Yellow_Text(porosityTimer_.lastTime())<<
                Yellow_Text(" s")<<Foam::endl;
}


void pFlow::coupling::momentumSphereUnresolvedCouplingSystem::calculateMomentumCoupling()
{
    const auto& U = this->cMesh().mesh().template lookupObject<Foam::volVectorField>("U");

    const auto& vp = this->particleVelocity();

    const auto& wp = this->particleRVelocity();

    auto& fluidForce = this->fluidForce();

    auto& fluidTorque = this->fluidTorque();

    std::fill(fluidForce.begin(), fluidForce.end(), 0.0);
    std::fill(fluidTorque.begin(), fluidTorque.end(), 0.0);

    momentumInteraction_.calculateCoupling(U, vp, wp, fluidForce, fluidTorque);
}

void pFlow::coupling::momentumSphereUnresolvedCouplingSystem::calculateHeatCoupling()
{
    notImplementedFunction;
}

void pFlow::coupling::momentumSphereUnresolvedCouplingSystem::calculateMassCoupling()
{
    notImplementedFunction;
}

Foam::tmp<Foam::volScalarField> 
    pFlow::coupling::momentumSphereUnresolvedCouplingSystem::Sp()const
{
    return Foam::tmp<Foam::volScalarField>(momentumInteraction_.Sp());
}

Foam::tmp<Foam::volVectorField> 
    pFlow::coupling::momentumSphereUnresolvedCouplingSystem::Su()const
{

    const auto& SU = momentumInteraction_.Su();
    auto tmpLift = momentumInteraction_.liftForce();
    const auto& lift = tmpLift.ref();

    auto SUall = Foam::tmp<Foam::volVectorField>::New
    (
        Foam::IOobject
        (
            "SUall",
            Foam::timeName(this->cMesh().mesh().time()),
            this->cMesh().mesh(),
            Foam::IOobject::NO_READ,
            Foam::IOobject::NO_WRITE,
            false
        ),
        SU + lift
    );
    return SUall;
}

