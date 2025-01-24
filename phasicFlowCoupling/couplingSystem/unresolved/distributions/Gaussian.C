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

// from OpenFOAM
#include "Gaussian.hpp"
#include "couplingMesh.hpp"

void pFlow::coupling::Gaussian::constructLists(const Foam::scalar searchLen)
{
	const Foam::vectorField& cellC = mesh_.cellCentres();
	const Foam::scalarField& cellV = mesh_.cellVolumes();
	const auto& bndris = mesh_.boundary();
	const auto nCells = cellC.size();
	
	// check for capacity and size 
	if(nCells != neighborList_.size() )
	{
		neighborList_.clear();
		neighborList_.setSize(nCells);
		boundaryCell_.clear();
		boundaryCell_.setSize(nCells);
	}
	
	// loop over all cells
	// TODO: later this should be recursive calls 
	#pragma omp parallel for
	for(Foam::label i=0; i<neighborList_.size(); i++)
	{
		const Foam::vector& ci = cellC[i];
		const Foam::scalar  lCell = 0.5* Foam::pow(cellV[i], 0.33333);
		std::vector<Foam::label> vN;
		vN.reserve(50);
		vN.clear();
		vN.push_back(i);

		for(Foam::label j=0; j< nCells; j++)
		{
			if( i == j)continue;

			auto dist = Foam::mag(cellC[j]-ci);
			if(dist < searchLen + lCell )
			{
				vN.push_back(j);
			}
		}

		neighborList_[i].setSize(vN.size());
		auto& lst = neighborList_[i];
		for(size_t jj=0ul; jj<vN.size(); jj++)
		{
			lst[jj] = vN[jj];
		}

		Foam::label nghbrB = -1;
		Foam::label nghbrFaceIndex = -1;
		Foam::scalar minDistance = 1.0e15;
		for(Foam::label nB = 0; nB < bndris.size(); nB++)
		{
			//const labelUList& faceCells = bndris[nB].faceCells();
			const auto& bndry = bndris[nB];

			for(auto j =0; j<bndry.size(); j++)
			{
				auto cellFaceDist = Foam::mag(ci - bndry.Cf()[j]);
				if( cellFaceDist < searchLen + lCell)
				{
					if(cellFaceDist<minDistance)
					{
						nghbrB = nB;
						nghbrFaceIndex = j;
						minDistance = cellFaceDist;
					}		
				}
			}
		}
		
		boundaryCell_[i] = std::pair<Foam::label, Foam::label>{nghbrB, nghbrFaceIndex};
	}

}


pFlow::coupling::Gaussian::Gaussian
(
	Foam::dictionary 		dict, 
	const couplingMesh& 	cMesh,
	const Plus::centerMassField& centerMass
)
:
	neighborLength_(dict.lookup<Foam::scalar>("neighborLength")),
	weights_("weights",centerMass),
	mesh_(cMesh.mesh())
{
	constructLists(3.0*neighborLength_);
}


void pFlow::coupling::Gaussian::updateWeights
(
	const Plus::procCMField<Foam::label> & parCellIndex
)
{

	const auto& centerMass = weights_.centerMass();
	const size_t numPar = centerMass.size();
	const Foam::scalar b2 = Foam::pow(neighborLength_,2);
	const Foam::vectorField& allCellCntr = mesh_.cellCentres();
	
	// lambda for distribution function 
	auto distFunc  = [b2](const Foam::vector& x) -> Foam::scalar
	{
		return Foam::exp(-dot(x,x)/b2);
	};

	auto distFunc2 = [b2](const Foam::vector& x, const Foam::scalar dist)-> Foam::scalar
	{
		auto ksi2 = dot(x,x)/b2;
		auto shifted_ksi2 = (dot(x,x)+4*dist*dist)/b2;
		return Foam::exp(-ksi2)+Foam::exp(-shifted_ksi2);
	};
	
	
	#pragma omp parallel for 
	for(size_t i=0; i<numPar; i++)
	{
		
		const Foam::label targetCellId = parCellIndex[i];
		auto& weightsI = weights_[i];

		weightsI.clear();
		if( targetCellId < 0 )continue;
		
		// get all the neighbors of cell 
		const Foam::labelList& neighbors = neighborList_[targetCellId];
		
		// center of particle 
		const realx3& cp = centerMass[i];
		Foam::vector CP{cp.x(), cp.y(), cp.z()};	

		weightsI.reserve(neighbors.size());
		
		Foam::scalar pSubTotal = 0;
		
		auto [bndryIndex, faceIndex] = boundaryCell_[targetCellId];

		if( bndryIndex == -1 )
		{
			// this is not a boundary cell 
			for(auto j:neighbors)
			{
				Foam::scalar k = distFunc(CP - allCellCntr[j]);
				weightsI.push_back(k);
				pSubTotal += k;
			}
		}
		else
		{
			// this is a boundary cell 
			const auto& bndry = mesh_.boundary()[bndryIndex];
			const Foam::vector parFace = CP - bndry.Cf()[faceIndex];
			const Foam::vector normal = Foam::normalised(bndry.Sf()[faceIndex]);
			Foam::scalar parFaceDist = std::abs(normal & parFace);

			for(auto j:neighbors)
			{	
				Foam::scalar k;
				if(bndryIndex== boundaryCell_[j].first )
				{
					k = distFunc2(CP - allCellCntr[j], parFaceDist);
				}
				else
					k = distFunc(CP - allCellCntr[j]);
				weightsI.push_back(k);
				pSubTotal += k;
			}

		}
		pSubTotal = Foam::max(pSubTotal, Foam::vSmall);
		for(auto& w:weightsI) w /= pSubTotal;
	}

}

/*bool pFlow::coupling::GaussianDistributionKernel::internalFieldUpdate()
{
	
	auto solidVoldTmp = Foam::volScalarField::Internal::New(
		"solidVol",
		this->mesh(),
		 Foam::dimensioned("solidVol", Foam::dimVol, Foam::scalar(0))
		 	);
	
	auto& solidVol = solidVoldTmp.ref();
	size_t numPar = centerMass_.size();
	const vectorField& allCellCntr = this->mesh().cellCentres();
	
	const auto b2 = Foam::pow(neighborLength(),2);

	// lambda for distribution function 
	auto distFunc  = [b2](const Foam::vector& x) -> Foam::scalar
	{
		return Foam::exp(-dot(x,x)/b2);
	};

	auto distFunc2 = [b2](const Foam::vector& x, const Foam::scalar dist)-> Foam::scalar
	{
		auto ksi2 = dot(x,x)/b2;
		auto shifted_ksi2 = (dot(x,x)+4*dist*dist)/b2;
		return Foam::exp(-ksi2)+Foam::exp(-shifted_ksi2);
	};
	
	
	#pragma omp parallel for 
	for(size_t i=0; i<numPar; i++)
	{
		real parVolume = static_cast<real>(3.14159265358979/6)*
				pFlow::pow(particleDiameter_[i], static_cast<real>(3.0));
		
		const Foam::label targetCellId = parCellIndex_[i];

		if( targetCellId < 0 )continue;
		
		// get all the neighbors of cell 
		const auto& nList = neighborList_[targetCellId];
		
		// center of particle 
		const realx3& cp = centerMass_[i];
		Foam::vector CP{cp.x(), cp.y(), cp.z()};	


		std::vector<Foam::scalar> ks; 
		ks.clear();
		ks.reserve(nList.size());
		
		Foam::scalar pSubTotal = 0;
		Foam::label bndryIndex = boundaryCell_[targetCellId].first;
		Foam::label faceIndex = boundaryCell_[targetCellId].second;
		

		if( bndryIndex == -1 )
		{
			// this is not a boundary cell 
			for(auto j:nList)
			{
				Foam::scalar k = distFunc(CP - allCellCntr[j]);
				ks.push_back(k);
				pSubTotal += k;
			}
		}
		else
		{
			// this is a boundary cell 
			const auto& bndry = cMesh_.mesh().boundary()[bndryIndex];
			const Foam::vector parFace = CP - bndry.Cf()[faceIndex];
			const Foam::vector normal = normalised(bndry.Sf()[faceIndex]);
			Foam::scalar parFaceDist = std::abs(normal & parFace);

			for(auto j:nList)
			{	
				Foam::scalar k;
				if(bndryIndex== boundaryCell_[j].first )
				{
					k = distFunc2(CP - allCellCntr[j], parFaceDist);
				}
				else
					k = distFunc(CP - allCellCntr[j]);
				ks.push_back(k);
				pSubTotal += k;
			}

		}
		
		size_t n = 0;
		pSubTotal = Foam::max(pSubTotal, Foam::vSmall);
		for(auto j:nList)
		{
			#pragma omp atomic
			solidVol[j] += (ks[n++] /pSubTotal * parVolume);
		}

	}
	
	this->ref() = Foam::max(
		1 - solidVol/this->mesh().V(), 
		static_cast<Foam::scalar>(this->alphaMin()));

	return true;
}*/



