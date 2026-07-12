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

#include "distributionSolidAveraging.hpp"
#include "porosity.hpp"
#include "unresolvedCouplingSystem.hpp"

pFlow::coupling::distributionSolidAveraging::distributionSolidAveraging
(
    const word&                     type,
    const unresolvedCouplingSystem& uCS,
    const porosity&                 prsty,
    const word&                     name
)
:
    solidAveraging(type, uCS, prsty, name),
    cellAvField_
    (
        Foam::IOobject
	    (
	        name,
	        Foam::timeName(prsty.mesh().time()),
	        prsty.mesh(),
	        Foam::IOobject::NO_READ,
	        Foam::IOobject::AUTO_WRITE
	    ),
    	prsty.mesh(),
        Foam::dimensionedVector(Foam::dimless, Foam::vector(0,0,0))
    ),
    calcParAvField_( 
        "calcParAvField."+name,
        prsty.centerMass()
    )
{}

void pFlow::coupling::distributionSolidAveraging::calculate
(
    const Plus::realx3ProcCMField& particleField
)
{
    const auto&  parCellInd = porosity_.parCellIndex();
    const size_t numPar = parCellInd.size();        
    const auto&  parDiam =  porosity_.particleDiameter();
    const auto&  cellVol = porosity_.mesh().V();
    const auto&  alpha = porosity_.alpha();
    const auto& distributor = uCS_.distribution();
    
    forAll(cellAvField_,celli)
    {
        cellAvField_[celli] = Foam::Zero;
    }

    #pragma omp parallel for schedule (dynamic)
    for(size_t i=0; i<numPar; i++)
    {
        Foam::scalar pVol = pFlow::Pi/6 *
                Foam::pow(parDiam[i], static_cast<real>(3.0));

        const Foam::label cellId = parCellInd[i];
        if( cellId >= 0 )
        {
            Foam::vector upv{ 
                pVol*particleField[i].x(), 
                pVol*particleField[i].y(), 
                pVol*particleField[i].z()};

            distributor.distributeValue_OMP(i, cellId, cellAvField_, upv);				
        }
    }

    forAll(cellAvField_,celli)
    {
        cellAvField_[celli] /= Foam::max( (1-alpha[celli])*cellVol[celli], Foam::SMALL);
    }

    distributor.smoothenField(cellAvField_);
    cellAvField_.correctBoundaryConditions();

    #pragma omp parallel for schedule (dynamic)
    for(size_t i=0; i<numPar; i++)
    {
        calcParAvField_[i] = realx3(
            cellAvField_[parCellInd[i]][0],
            cellAvField_[parCellInd[i]][1],
            cellAvField_[parCellInd[i]][2]);
    }

    parAvFieldRef_ = &calcParAvField_;
}

void pFlow::coupling::distributionSolidAveraging::calculateNumberBased
(
    const Plus::realx3ProcCMField& particleField
)
{
    const auto&  parCellInd = porosity_.parCellIndex();
    const size_t numPar = parCellInd.size();        
    const auto&  parDiam =  porosity_.particleDiameter();
    const auto& distributor = uCS_.distribution();
    
    forAll(cellAvField_,celli)
    {
        cellAvField_[celli] = Foam::Zero;
    }

    #pragma omp parallel for schedule (dynamic)
    for(size_t i=0; i<numPar; i++)
    {
        Foam::scalar pVol = pFlow::Pi/6 *
                Foam::pow(parDiam[i], static_cast<real>(3.0));

        const Foam::label cellId = parCellInd[i];
        if( cellId >= 0 )
        {
            Foam::vector upv{ 
                pVol*particleField[i].x(), 
                pVol*particleField[i].y(), 
                pVol*particleField[i].z()};

            distributor.distributeValue_OMP(i, cellId, cellAvField_, upv);				
        }
    }

    
    distributor.smoothenField(cellAvField_);
    cellAvField_.correctBoundaryConditions();

    #pragma omp parallel for schedule (dynamic)
    for(size_t i=0; i<numPar; i++)
    {
        calcParAvField_[i] = realx3(
            cellAvField_[parCellInd[i]][0],
            cellAvField_[parCellInd[i]][1],
            cellAvField_[parCellInd[i]][2]);
    }

    parAvFieldRef_ = &calcParAvField_;
}
