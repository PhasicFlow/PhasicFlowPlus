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
void pFlow::coupling::sphereDrag<DistributorType, DragClosureType, useCellDistribution>::
	calculateDragForce
	(
		const DragClosureType& 			dragClosure,
		const Foam::volVectorField& 	fluidVelocity,
		const Plus::realx3ProcCMField& 	parVelocity,
		const Plus::realProcCMField& 	diameter,
		Plus::realx3ProcCMField& 		particleForce
	)
{

	self selfDictribution;
	const auto& cellDist = this->cellDistribution();

	const auto& parCells =  this->particleCellIndex();
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

#pragma omp parallel for
	for(size_t parIndx=0; parIndx<numPar; parIndx++)
	{
		auto cellIndx = parCells[parIndx];

		if(cellIndx >= 0 )
		{
			auto rhoi = rho[cellIndx];
			auto mui = nu[cellIndx]* rhoi;
			auto ef = alpha[cellIndx];
			auto dp = diameter[parIndx];
			auto vp = Foam::constant::mathematical::pi/6 * Foam::pow(dp,3.0);

			Foam::vector up = {
				parVelocity[parIndx].x(), 
				parVelocity[parIndx].y(),
				parVelocity[parIndx].z()};

			Foam::vector ur = fluidVelocity[cellIndx]-up;
			Foam::scalar Re = ef * rhoi * Foam::mag(ur) * dp /mui;

			Foam::scalar sp = 3 * Foam::constant::mathematical::pi * mui * ef * dp * dragClosure.dimlessDrag(Re, ef);
			
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
}


template<typename DistributorType, typename DragClosureType, bool useCellDistribution>
pFlow::coupling::sphereDrag<DistributorType, DragClosureType, useCellDistribution>::sphereDrag
(
	const UnresolvedCouplingSystem<DistributorType>& uCS, 
	const porosity& 								prsty
)
:
	DragType(uCS, prsty),
	dragClosure_(uCS.unresolvedDict().subDict("drag"))
{

}

template<typename DistributorType, typename DragClosureType, bool useCellDistribution>
void 
pFlow::coupling::sphereDrag<DistributorType, DragClosureType, useCellDistribution>::
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
}