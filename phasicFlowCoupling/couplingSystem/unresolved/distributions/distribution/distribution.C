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
#include "streams.hpp"


void pFlow::coupling::distribution::parseNeighbors
(
    const Foam::label       targetCelli,
    const Foam::label       layerNumber,
    const Foam::label       maxLayers, 
    const Foam::label       celli, 
    std::set<Foam::label>&  finalList       
)
{
    const labelListList& cellPoints = mesh().cellPoints();
    const labelListList& pointCells = mesh().pointCells();

    const auto& points = cellPoints[celli];
    
    std::vector<Foam::label> thisNewCells;

    int nFound = 0;
    for(const auto& p:points)
    {
        const auto& cells = pointCells[p];
        for(auto c:cells)
        {
            if( c == targetCelli || finalList.count(c) == 1) continue;
            nFound++;
            thisNewCells.push_back(c);
            finalList.insert(c);
        }
    }

    if( layerNumber < maxLayers && nFound != 0 )
    {
        for(const auto& newCell: thisNewCells)
        {
            parseNeighbors
            (
                targetCelli,
                layerNumber+1,
                maxLayers,
                newCell,
                finalList
            );
        }
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
    const Foam::labelList& neigborCells = mesh().cellCells(celli);
    const Foam::vectorField& cellC = mesh().cellCentres(); 
    const Foam::scalar cellV = mesh().cellVolumes()[targetCelli];
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

    if(layerNumber < maxLayers && found != 0)
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

void pFlow::coupling::distribution::constructLists(
    const Foam::scalar searchLen,
    const Foam::label maxLayers)
{
    if(!mesh().hasCellCells())
	{
		mesh().cellCells();
	}

	const Foam::vectorField& cellC = mesh().cellCentres();
   	const auto nCells = cellC.size();

	
	// check for capacity and size 
	if(static_cast<size_t>(nCells) != neighborList_.size() )
	{
		neighborList_.clear();
		neighborList_.reserve(nCells);
        neighborList_.resize(nCells);
	}

    // loop over all cells
	#pragma omp parallel for schedule (dynamic)
	for(Foam::label celli = 0; celli < nCells; celli++) 
	{
        const Foam::vector& ci = cellC[celli];
        
        std::set<Foam::label> finalList;

        finalList.insert(celli);

        parseNeighbors(celli, ci, searchLen, celli, finalList, 1, maxLayers);        
        neighborList_[celli].clear();
        for(const auto nbr:finalList)
        {
            neighborList_[celli].push_back(nbr);
        }
	}
}

void pFlow::coupling::distribution::constructLists(const Foam::label maxLayers)
{
    if(!mesh().hasPointCells())
    {
        mesh().pointCells();
        mesh().cellPoints();
    }

    const auto nCells = mesh().nCells();
    
    // check for capacity and size 
	if(static_cast<size_t>(nCells) != neighborList_.size() )
	{
		neighborList_.clear();
		neighborList_.reserve(nCells);
        neighborList_.resize(nCells);
	}

    #pragma omp parallel for schedule (dynamic)
    for(Foam::label celli = 0; celli < nCells; celli++)
    {
        
        std::set<Foam::label> finalList;

        finalList.insert(celli);

        parseNeighbors(celli, 1, maxLayers, celli, finalList);

        neighborList_[celli].clear();
        for(const auto nbr:finalList)
        {
            neighborList_[celli].push_back(nbr);
        }
    }
}

pFlow::coupling::distribution::distribution(
    const Foam::dictionary& 	 parrentDict, 
    const couplingMesh& 	     cMesh,
    const Plus::centerMassField& centerMass
)
:
    distributionBase(true, parrentDict, cMesh, centerMass)
{}


void pFlow::coupling::distribution::smoothenField(Foam::volScalarField &field) const
{
}

void pFlow::coupling::distribution::smoothenField(Foam::volVectorField &field) const
{
}
