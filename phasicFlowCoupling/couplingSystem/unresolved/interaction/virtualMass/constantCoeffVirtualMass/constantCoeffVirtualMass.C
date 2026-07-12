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

#include "fvc.H"

#include "constantCoeffVirtualMass.hpp"
#include "unresolvedCouplingSystem.hpp"

pFlow::coupling::constantCoeffVirtualMass::constantCoeffVirtualMass
(
    const unresolvedCouplingSystem& uCS,
    const porosity& prsty
)
:
    virtualMass(uCS, prsty),
    Cvm_(this->dict().getOrDefault<Foam::scalar>("Cvm", 0.5)),
    isCompressible_
    (
      this->mesh().template lookupObject<Foam::volScalarField>("p").dimensions() 
      == Foam::dimPressure
    )

    
{
    tmpVirtualMassForce_ = Foam::tmp<Foam::volVectorField>::New
    (
        Foam::IOobject
        (
            "virtualMassForce",
            Foam::timeName(this->mesh().time()),
            this->mesh(),
            Foam::IOobject::READ_IF_PRESENT,
            (this->printVirtualMass()?Foam::IOobject::AUTO_WRITE:Foam::IOobject::NO_WRITE)
        ),
        this->mesh(),
        Foam::dimensionedVector
        (
            "virtualMassForce",
            Foam::dimensionSet(1,-2,-2,0,0),
            Foam::vector(0,0,0)
        )
    );
}

void pFlow::coupling::constantCoeffVirtualMass::calculateVirtualMassForce
(
    const Foam::volVectorField& U,
    const Plus::realx3ProcCMField& parAcc,
    const Plus::realProcCMField& diameter,
    Plus::realx3ProcCMField& particleForce
)
{
    auto& vMF = tmpVirtualMassForce_.ref();

    forAll(vMF, celli)
    {
        vMF[celli] = Foam::vector(0,0,0);
    }

    const auto& rho = this->mesh().template lookupObject<Foam::volScalarField>("rho");
    const auto& parCellInd = this->parCellIndex();
    const auto& Vcells = this->mesh().V();
    const auto& phi = this->mesh().template lookupObject<Foam::surfaceScalarField>("phi");
    
    Foam::tmp<Foam::volVectorField> tDDtU;
    
    /// Calculate the material derivative of fluid velocity 
    if(isCompressible_)
    {
        tDDtU  = (fvc::ddt(rho, U) + fvc::div(phi, U))/rho;
    }
    else
    {
        tDDtU = fvc::ddt(U) + fvc::div(phi, U);
    }

    const auto& DDtU = tDDtU();
    
    const size_t nPar = diameter.size();

    #pragma omp parallel for schedule(dynamic)
    for(size_t i=0; i<nPar; ++i)
    {
        const Foam::label cellI = parCellInd[i];

        if(cellI < 0) continue;

        const Foam::scalar dp = diameter[i];
        const Foam::scalar Vp = Foam::constant::mathematical::pi/6.0 * Foam::pow(dp,3.0);
        
        /// calculation for particle acceleration
        const Foam::vector ap(parAcc[i].x(), parAcc[i].y(), parAcc[i].z());
               
   	/// Calculate virtual mass force
        const Foam::vector vmForce = Cvm_ * rho[cellI] * Vp * (DDtU[cellI] - ap);
        
    
        particleForce[i] += realx3(vmForce.x(), vmForce.y(), vmForce.z());

        #pragma omp atomic
        vMF[cellI].x() += vmForce.x();

        #pragma omp atomic
        vMF[cellI].y() += vmForce.y();

        #pragma omp atomic
        vMF[cellI].z() += vmForce.z();
    }

    forAll(Vcells, celli)
    {
        vMF[celli] /= Vcells[celli];
    }

    vMF.correctBoundaryConditions();
}
