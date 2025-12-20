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

#include "Saffmann.hpp"
#include "unresolvedCouplingSystem.hpp"

pFlow::coupling::Saffmann::Saffmann
(
    const unresolvedCouplingSystem &uCS, 
    const porosity &prsty
)
:
    lift(uCS, prsty),
    residualRe_
    (
        this->dict().template get<Foam::scalar>("residualRe")
    )
{
    this->setLiftActive(true);

    tmpLiftForce_ = Foam::tmp<Foam::volVectorField>::New
    (
        Foam::IOobject
        (
            "liftForce",
            Foam::timeName(this->mesh().time()),
            this->mesh(),
            Foam::IOobject::READ_IF_PRESENT,
            (this->printLift()?Foam::IOobject::AUTO_WRITE:Foam::IOobject::NO_WRITE)
        ),
        this->mesh(),
        Foam::dimensionedVector
        (
            "liftForce",
            Foam::dimensionSet(1,-2,-2,0,0),
            Foam::vector(0,0,0)
        )
    );
}


void pFlow::coupling::Saffmann::calculateLiftForce
(
    const Foam::volVectorField&     U,
    const Plus::realx3ProcCMField&  parVel,
    const Plus::realx3ProcCMField&  parRotVel,
    const Plus::realProcCMField&    diameter,
    Plus::realx3ProcCMField&        particleForce,
    Plus::realx3ProcCMField&        particleTorque
) 
{
    auto& liftForce = tmpLiftForce_.ref();

    // initialize all terms to zero
	forAll(liftForce, celli)
	{
		liftForce[celli] = Foam::vector(0,0,0);
	}

    const size_t nPar = diameter.size();
    const auto& parCellInd =  this->parCellIndex();
    const auto& nu = this->mesh().template lookupObject<Foam::volScalarField>("nu");
    const auto& rho = this->mesh().template lookupObject<Foam::volScalarField>("rho");
    
    auto curlUPtr = this->fluidVorticity(U);
    const auto& curlU = curlUPtr.ref();

    #pragma omp parallel for schedule(dynamic)
    for(size_t i=0; i<nPar; ++i)
    {
        const Foam::label cellI = parCellInd[i];
        
        if(cellI == -1 )continue;

        const Foam::vector uf = U[cellI];
        const Foam::vector up{
            parVel[i].x(), 
            parVel[i].y(), 
            parVel[i].z()};

        const Foam::vector wp{
            parRotVel[i].x(),
            parRotVel[i].y(),
            parRotVel[i].z()
        };

        const Foam::scalar dp = diameter[i];

        const Foam::vector uRelVec = up - uf;
        const Foam::scalar uRel = Foam::mag(uRelVec);

        const Foam::scalar Rep = Foam::max(uRel*dp/nu[cellI], residualRe_);
        const Foam::scalar Rew = Foam::max(Foam::mag(wp)*dp*dp/nu[cellI], Foam::SMALL);
        const Foam::scalar Res = Foam::max(Foam::mag(curlU[cellI])*dp*dp/nu[cellI], Foam::SMALL);
        const Foam::scalar Sr = Res/Rep;
        const Foam::scalar Rr = Rew/Rep;
        
        // spin lift coefficient
        const Foam::scalar Cl_spin = Rr ;
                    
        // Fluid shear lift coefficient
        const Foam::scalar J = 2.255;
        const Foam::scalar Cl_shear = (18.0/Foam::constant::mathematical::pi*Foam::constant::mathematical::pi)*
                Foam::sqrt(Sr/Rep)*J - 11.0/8.0*Sr*Foam::exp(-0.5*Rep);

        
        Foam::vector dirVectorShear = curlU[cellI] ^ uRelVec;
        Foam::vector dirVectorSpin = wp ^ uRelVec;

        dirVectorShear /= Foam::max(Foam::mag(dirVectorShear), Foam::SMALL);
        dirVectorSpin /= Foam::max(Foam::mag(dirVectorSpin), Foam::SMALL);

        Foam::vector lF = (Foam::constant::mathematical::pi/8.0*dp*dp*rho[cellI]) * uRel * uRel *
            ( 
                Cl_shear *  Foam::vector(dirVectorShear) + // fluid shear contribution
                Cl_spin  *  Foam::vector(dirVectorSpin)    // particle rotation contribution
            );
            
        
        particleForce[i] += realx3(lF.x(), lF.y(), lF.z());

        #pragma omp atomic
            liftForce[cellI].x() += lF.x();
        
        #pragma omp atomic
            liftForce[cellI].y() += lF.y();

        #pragma omp atomic
            liftForce[cellI].z() += lF.z();
    }

    const auto& Vcells = this->mesh().V();

    forAll(Vcells, celli)
    {
        liftForce[celli] /= Vcells[celli];
    }

    liftForce.correctBoundaryConditions();

}