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

#include <set>

// from OpenFOAM
#include "Gaussian2.hpp"
#include "couplingMesh.hpp"

void pFlow::coupling::Gaussian2::parseNeighbors(
	const Foam::label 		targetCelli,
	const Foam::vector& 	targetCellCentre, 
	const Foam::scalar 		searchLen,
	const Foam::scalar 		celli,
	std::set<Foam::label>& 	finalList,
	const Foam::label 		layerNumber)
{
	const Foam::labelList& neigborCells = mesh_.cellCells(celli);
	const Foam::vectorField& cellC = mesh_.cellCentres(); 
	const Foam::scalar cellV = mesh_.cellVolumes()[targetCelli];
	const Foam::scalar lCell = 0.5* Foam::pow(cellV, 0.33333);
	// first find all the cells around 
	int found = 0; 
	for(auto nghbrCelli:neigborCells)
	{
		if(nghbrCelli==targetCelli || finalList.count(nghbrCelli) == 1) continue;
		
		Foam::scalar dist = Foam::mag(cellC[nghbrCelli] - targetCellCentre);	
		if(dist <= searchLen + lCell)
		{
			found++;
			finalList.insert(nghbrCelli);
		}
	}

    if(layerNumber <=4 && found != 0)
    {
    	for(auto nghbrCelli:neigborCells)
    	{
    		parseNeighbors
    		(
    			targetCelli,
    			targetCellCentre,
    			searchLen,
    			nghbrCelli,
    			finalList,
    			layerNumber+1
    		);
    	}
    }
}


void pFlow::coupling::Gaussian2::constructLists(const Foam::scalar searchLen)
{

	if(!mesh_.hasCellCells())
	{
		mesh_.cellCells();
	}

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
	for(Foam::label celli = 0; celli < nCells; celli++) 
	{
		const Foam::vector& ci = cellC[celli];
		const Foam::scalar  lCell = 0.5* Foam::pow(cellV[celli], 0.33333);
		
		std::set<Foam::label> finalList;

		finalList.insert(celli);

        parseNeighbors(celli, ci, searchLen, celli, finalList, 1);        
         
        neighborList_[celli].setSize(finalList.size());
        Foam::label n = 0;
        for(auto nbr:finalList)
        {
        	neighborList_[celli][n]=nbr;
        	n++;
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
		
		boundaryCell_[celli] = std::pair<Foam::label, Foam::label>{nghbrB, nghbrFaceIndex};
	}

}


pFlow::coupling::Gaussian2::Gaussian2
(
	Foam::dictionary 		dict, 
	const couplingMesh& 	cMesh,
	const Plus::centerMassField& centerMass
)
:
	neighborLength_(dict.lookup<Foam::scalar>("neighborLength")),
	lengthExtent_(dict.lookupOrDefault<Foam::scalar>("lengthExtent", 3.0)),
	weights_("weights",centerMass),
	mesh_(cMesh.mesh())
{
	constructLists(lengthExtent_*neighborLength_);

}


void pFlow::coupling::Gaussian2::updateWeights
(
	const Plus::procCMField<Foam::label> & parCellIndex,
		const Plus::procCMField<real> &
)
{

	const auto& centerMass = weights_.centerMass();
	const size_t numPar = centerMass.size();
	const Foam::scalar b2 = Foam::pow(neighborLength_,2);
	const Foam::vectorField& allCellCntr = mesh_.cellCentres();
	
	// lambda for distribution function 
	auto distFunc  = [b2](const Foam::vector& x) -> Foam::scalar
	{
		return Foam::exp(-dot(x,x)/(2*b2));
	};

	auto distFunc2 = [b2](const Foam::vector& x, const Foam::scalar dist)-> Foam::scalar
	{
		auto ksi2 = dot(x,x)/(2*b2);
		auto shifted_ksi2 = (dot(x,x)+4*dist*dist)/(2*b2);
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

