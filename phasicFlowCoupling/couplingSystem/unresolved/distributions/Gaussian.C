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

#include "Gaussian.hpp"
#include "couplingMesh.hpp"
#include "schedule.hpp"

void pFlow::coupling::Gaussian::constructBoundaryLists
(
	const Foam::scalar searchLen
)
{
	if(!mesh_.hasCellCells())
	{
		mesh_.cellCells();
	}

	const Foam::vectorField& cellC = mesh_.cellCentres();
	const auto nCells = cellC.size();
	const auto& bndris = mesh_.boundary();
	
	// check for capacity and size 
	if(static_cast<size_t>(nCells) != boundaryCell_.size() )
	{	
		boundaryCell_.clear();
		boundaryCell_.reserve(nCells);
		boundaryCell_.resize(nCells);
	}

    // loop over all cells
	#pragma omp parallel for schedule (dynamic)
	for(Foam::label celli = 0; celli < nCells; celli++) 
	{
        const Foam::vector& ci = cellC[celli];
        
		Foam::label nghbrB = -1;
		Foam::label nghbrFaceIndex = -1;
		Foam::scalar minDistance = 1.0e15;
		for(Foam::label nB = 0; nB < bndris.size(); nB++)
		{
			//const labelUList& faceCells = bndris[nB].faceCells();
			const auto& bndry = bndris[nB];

			for(auto j = 0; j<bndry.size(); j++)
			{
				auto cellFaceDist = Foam::mag(ci - bndry.Cf()[j]);
				if( cellFaceDist < searchLen  && cellFaceDist<minDistance)
				{
					nghbrB = nB;
					nghbrFaceIndex = j;
					minDistance = cellFaceDist;		
				}
			}
		}
		
		boundaryCell_[celli] = {nghbrB, nghbrFaceIndex};
	}
}

pFlow::coupling::Gaussian::Gaussian
(
	Foam::dictionary 		dict, 
	const couplingMesh& 	cMesh,
	const Plus::centerMassField& centerMass
)
:
	distribution(dict, cMesh, centerMass),	
	standardDeviation_(lookupDict<Foam::scalar>(dict, "standardDeviation")),
	maxLayers_(lookupOrDefaultDict<Foam::label>(dict, "maxLayers", 2))
{
	constructLists(3*standardDeviation_, maxLayers_);
	constructBoundaryLists(1.5*standardDeviation_ );
}

void pFlow::coupling::Gaussian::updateWeights
(
	const Plus::procCMField<Foam::label> & parCellIndex,
	const Plus::procCMField<real> &
)
{

	const auto& centerMass = weights_.centerMass();
	const size_t numPar = centerMass.size();
	const Foam::scalar b2 = Foam::pow(standardDeviation_,2);
	const Foam::vectorField& allCellCntr = mesh_.cellCentres();
	
	// lambda for distribution function 
	auto distFunc  = [b2](const Foam::vector& x) -> Foam::scalar
	{
		return Foam::exp(-Foam::dot(x,x)/(2*b2));
	};

	auto distFunc2 = [b2](const Foam::vector& x, const Foam::scalar dist)-> Foam::scalar
	{
		auto ksi2 = Foam::dot(x,x)/(2*b2);
		auto shifted_ksi2 = (dot(x,x)+4*dist*dist)/(2*b2);
		return Foam::exp(-ksi2)+Foam::exp(-shifted_ksi2);
	};
	
	#pragma omp parallel for schedule (dynamic)
	for(size_t i=0; i<numPar; i++)
	{
		
		const Foam::label targetCellId = parCellIndex[i];
		auto& weightsPar = weights_[i];

		weightsPar.clear();
		if( targetCellId < 0 )continue;
		
		// get all the neighbors of cell 
		const auto& neighbors = neighborList_[targetCellId];
		
		// center of particle 
		const realx3& cp = centerMass[i];
		Foam::vector CP{cp.x(), cp.y(), cp.z()};	
		
		Foam::scalar pSubTotal = 0;
		
		auto [bndryIndex, faceIndex] = boundaryCell_[targetCellId];

		if( bndryIndex == -1 )
		{
			// this is not a boundary cell 
			
			for(auto cellId:neighbors)
			{	
				Foam::scalar f = distFunc(CP - allCellCntr[cellId]);
				weightsPar.push_back({cellId,f});
				pSubTotal += f;	
			}
		}
		else
		{
			// this is a boundary cell 
			const auto& bndry = mesh_.boundary()[bndryIndex];
			const Foam::vector parFace = CP - bndry.Cf()[faceIndex];
			const Foam::vector normal = Foam::normalised(bndry.Sf()[faceIndex]);
			Foam::scalar parFaceDist = std::abs(normal & parFace);

			for(auto cellId:neighbors)
			{	
				Foam::scalar f;
				if(bndryIndex== boundaryCell_[cellId].first )
				{
					f = distFunc2(CP - allCellCntr[cellId], parFaceDist);
				}
				else
				{
					f = distFunc(CP - allCellCntr[cellId]);
				}	
				weightsPar.push_back({cellId, f});
				pSubTotal += f;
			}

		}
		
		pSubTotal = Foam::max(pSubTotal, static_cast<Foam::scalar>(1.0e-10));
		for(auto& [i, w]:weightsPar) w /= pSubTotal;
	}

}

