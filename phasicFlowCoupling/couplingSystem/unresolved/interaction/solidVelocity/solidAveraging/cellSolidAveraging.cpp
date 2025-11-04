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

template<typename DistributionType>
pFlow::coupling::cellSolidAveraging<DistributionType>::cellSolidAveraging
(
    const word& type,
    const couplingMesh& cMesh,
    Plus::realx3ProcCMField& particleField
)
:
    solidAveraging(type, cMesh, particleField),
    cellAvField_
    (
        Foam::IOobject
	    (
	        "cellAvField."+particleField.name(),
	        Foam::timeName(cMesh.mesh().time()),
	        cMesh.mesh(),
	        Foam::IOobject::NO_READ,
	        Foam::IOobject::AUTO_WRITE
	    ),
    	cMesh.mesh(),
        Foam::dimensionedVector(Foam::dimless, Foam::vector(0,0,0))
    ),
    calcParAvField_( 
        "calcParAvField_"+particleField.name(),
        cMesh.couplingMeshProcCMField().centerMass()
    )
{}

template<typename DistributionType>
void pFlow::coupling::cellSolidAveraging<DistributionType>::calculate
(
    const Plus::realx3ProcCMField& particleField
    const DistributorType& cellDistributor,
    const porosity& prsty
)
{
    const auto&  parCellInd = prsty.parCellIndex();
    const size_t numPar = parCellInd.size();        
    const auto&  parDiam =  prsty.particleDiameter();
    const auto&  cellVol = prsty.mesh().V();
    const auto&  alpha = prsty.alpha();
    
    
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

    forAll(Us_,celli)
    {
        cellAvField_[celli] /= Foam::max( (1-alpha[celli])*cellVol[celli], Foam::SMALL);
    }

    distributor.smoothenField(cellAvField_);
    cellAvField_.correctBoundaryConditions();

    #pragma omp parallel for schdule (dynamic)
    for(size_t i=0; i<numPar; i++)
    {
        calcParAvField_[i] = cellAvField_[parCellInd[i]];
    }

    parAvFieldRef_ = calcParAvField_;
}

template<typename DistributionType>
void pFlow::coupling::cellSolidAveraging<DistributionType>::calculateNumberBased
(
    const Plus::realx3ProcCMField& particleField
    const DistributorType& cellDistributor,
    const porosity& prsty
)
{}
