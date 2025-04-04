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

#include "distribution.hpp"
#include "couplingMesh.hpp"

void pFlow::coupling::distribution::constructLists(
    const Foam::scalar searchLen,
    const Foam::label maxLayers)
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
	if(static_cast<size_t>(nCells) != neighborList_.size() )
	{
		neighborList_.clear();
		neighborList_.reserve(nCells);
        neighborList_.resize(nCells);
		boundaryCell_.clear();
		boundaryCell_.reserve(nCells);
        neighborList_.resize(nCells);
	}
	
	// loop over all cells
	#pragma omp parallel for
	for(Foam::label celli = 0; celli < nCells; celli++) 
	{
        const Foam::vector& ci = cellC[celli];
        const Foam::scalar  lCell = 0.5* Foam::pow(cellV[celli], 0.33333);
        
        std::set<Foam::label> finalList;

        finalList.insert(celli);

        parseNeighbors(celli, ci, searchLen, celli, finalList, 1, maxLayers);        
        neighborList_[celli].clear();
        for(const auto nbr:finalList)
        {
            neighborList_[celli].push_back(nbr);
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
		
		boundaryCell_[celli] = {nghbrB, nghbrFaceIndex};
	}
}

void pFlow::coupling::distribution::parseNeighbors
(
    const Foam::label targetCelli, 
    const Foam::vector &targetCellCentre, 
    const Foam::scalar searchLen, 
    const Foam::scalar celli, 
    std::set<Foam::label> &finalList, 
    const Foam::label layerNumber, 
    const Foam::label maxLayers
)
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

    if(layerNumber <=maxLayers && found != 0)
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
                layerNumber+1,
                maxLayers
            );
        }
    }
}

pFlow::coupling::distribution::distribution(
    const Foam::dictionary& 	 dict, 
    const couplingMesh &cMesh, 
    const Plus::centerMassField &centerMass
)
:
    weights_("weights",centerMass),
    mesh_(cMesh.mesh())
{}

