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

template<typename DistributorType, typename DragClosureType, bool useCellDistribution>
void pFlow::coupling::grainDrag<DistributorType, DragClosureType, useCellDistribution>::
	calculateDragForce
	(
		const DragClosureType& 			dragClosure,
		const Foam::volVectorField& 	U,
		const Plus::realx3ProcCMField& 	parVelocity,
		const Plus::realProcCMField& 	diameter,
		Plus::realx3ProcCMField& 		particleForce
	)
{

	if(!fluidVelocityPtr_)
	{
		fluidVelocityPtr_ = makeUnique<fluidVelocity>(
			fVelocityType_,
			U,
			this->cMesh()
		);
	}

	if(!solidVelocityPtr_)
	{
		solidVelocityPtr_ = makeUnique<solidVelocity>(
			sVelocityType_,
			parVelocity,
			this->cMesh()
		);
	}

	auto& fluidVelocity = fluidVelocityPtr_();
	fluidVelocity.interpolate(this->cMesh());

	auto& solidVelocity = solidVelocityPtr_();
	solidVelocity.average(
		this->cellDistribution(),
		this->Porosity());

	self selfDictribution;
	const auto& cellDist = this->cellDistribution();

	const auto& parCells =  this->parCellIndex();
	const auto& nu = this->mesh().template lookupObject<Foam::volScalarField>("nu");
	const auto& rho = this->mesh().template lookupObject<Foam::volScalarField>("rho");
	const auto& Vcells = this->mesh().V();
	const auto& alpha = this->alpha();
	auto& Su = this->Su();
	auto& Sp = this->Sp();

	// gets pressure gradient 
	auto pGradPtr = this->pressureGradient(rho);
	const auto& pGrad = pGradPtr();

	size_t numPar = parCells.size();

	#pragma ParallelRegion
	for(size_t parIndx=0; parIndx<numPar; parIndx++)
	{
		auto cellIndx = parCells[parIndx];

		if(cellIndx >= 0 )
		{
			auto rhoi = rho[cellIndx];
			auto mui = nu[cellIndx]* rhoi;
			auto ef = alpha[cellIndx];
			auto dp = diameter[parIndx];
			auto cgf = courseGrainFactor_[parIndx];
			auto dps = dp/cgf;
			auto vp = Foam::constant::mathematical::pi/6 * Foam::pow(dp,3.0);

			Foam::vector up = solidVelocity.vSolid(cellIndx, parIndx);

			Foam::vector ur = fluidVelocity.uFluid(cellIndx, parIndx)-up;
			Foam::scalar Res = ef * rhoi * Foam::mag(ur) * dps /mui ;

			Foam::scalar sp = 3 * Foam::pow(cgf,3) * Foam::constant::mathematical::pi * 
					mui * ef * dps * dragClosure.dimlessDrag(Res, ef);
			
			Foam::vector pf = static_cast<real>(sp)*ur - vp*pGrad[cellIndx];
			particleForce[parIndx] = realx3(pf.x(), pf.y(), pf.z());
			
			if constexpr (useCellDistribution)
			{
				cellDist.distributeValue_OMP(parIndx, cellIndx, Su, -(sp*up));
				cellDist.distributeValue_OMP(parIndx, cellIndx, Sp,   sp);
			}
			else
			{
				selfDictribution.distributeValue_OMP(parIndx, cellIndx, Su, -(sp*up));
				selfDictribution.distributeValue_OMP(parIndx, cellIndx, Sp,   sp);
			}
			
		}
	}

	forAll(Vcells, i)
	{
		Su[i] /= Vcells[i];
		Sp[i] /= Vcells[i];	
	}

	if constexpr (useCellDistribution)
	{
		cellDist.smoothenField(Su);
		cellDist.smoothenField(Sp);
	}
	else
	{
		selfDictribution.smoothenField(Su);
		selfDictribution.smoothenField(Sp);
	}

}


template<typename DistributorType, typename DragClosureType, bool useCellDistribution>
pFlow::coupling::grainDrag<DistributorType, DragClosureType, useCellDistribution>::grainDrag
(
	const UnresolvedCouplingSystem<DistributorType>& uCS, 
	const porosity& 								prsty
)
:
	DragType(uCS, prsty),
	dragClosure_(uCS.unresolvedDict().subDict("drag")),
	courseGrainFactor_(
		static_cast<
			const grainUnresolvedCouplingSystem<DistributorType>&>(uCS).
				courseGrainFactor())
{
	const auto& dict = uCS.unresolvedDict().subDict("drag");

	Foam::word fVelType(dict.lookup("fluidVelocity"));
	Foam::word sVelType(dict.lookup("solidVelocity"));

	if(fVelType == "particle" || fVelType == "cell" )
	{
		fVelocityType_ = fVelType;
	}
	else
	{
		fatalErrorInFunction
			<<"Valid options for fluidVelocity in dictionary "
			<<dict.name()
			<<" are: particles, cell"<<endl;
		Plus::processor::abort(0);
	}

	if(sVelType == "particle" || sVelType == "cell" )
	{
		sVelocityType_ = sVelType;
	}
	else
	{
		fatalErrorInFunction
			<<"Valid options for solidVelocity in dictionary "
			<<dict.name()
			<<" are: particles, cell"<<endl;
		Plus::processor::abort(0);
	}
}

template<typename DistributorType, typename DragClosureType, bool useCellDistribution>
void 
pFlow::coupling::grainDrag<DistributorType, DragClosureType, useCellDistribution>::
	calculateDragForce
	(
		const Foam::volVectorField& 	fluidVelocity,
		const Plus::realx3ProcCMField& 	parVelocity,
		const Plus::realProcCMField& 	diameter,
		Plus::realx3ProcCMField& 		particleForce
	)
{
	this->setSuSpToZero();
	particleForce = realx3(0,0,0);

	calculateDragForce(dragClosure_, fluidVelocity, parVelocity, diameter, particleForce);

	this->Sp().correctBoundaryConditions();
	this->Su().correctBoundaryConditions();
}
