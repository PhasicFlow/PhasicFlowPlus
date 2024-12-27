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
#include "fvCFD.H"

#include "statistical.hpp"
#include "porosity.hpp"

bool pFlow::coupling::statistical::cellNeighborsSearch()
{
	
	if(!performCellNeighborSearch()) return true;


	const vectorField& cellC = cMesh_.mesh().cellCentres();
	const scalarField& cellV = cMesh_.mesh().cellVolumes();
	const auto& bndris = cMesh_.mesh().boundary();

	const size_t nCells = cMesh_.mesh().nCells();

	const Foam::scalar b = neighborLength_ * boundRatio();
	
	// check for capacity and size 
	if(nCells > neighborList_.capacity() )
	{
		neighborList_.clear();
		neighborList_.reserve(nCells);
		boundaryCell_.clear();
		boundaryCell_.reserve(nCells);
	}

	// clear the elements of each list 
	for(auto& nL: neighborList_)
	{
		nL.clear();
	}
	
	// check for size 
	if(nCells > neighborList_.size() )
	{
		neighborList_.resize(nCells);

		boundaryCell_.resize(nCells);
		/*std::fill(
			boundaryCell_.begin(), 
			boundaryCell_.end(),
			std::pair<Foam::label, Foam::label>{-1,-1});*/
	}

	
	
	// loop over all cells
	// TODO: later this should be recursive calls 
	//#pragma omp parallel for
	for(size_t i=0; i<20; i++)
	{
		const Foam::vector& ci = cellC[i];
		const Foam::scalar  lCell = 0.5* Foam::pow(cellV[i], 0.33333);
		
		neighborList_[i].push_back(i);

		for(size_t j=0; j< nCells; j++)
		{
			
			if( i == j)continue;

			auto dist = mag(cellC[j]-ci);
			if(dist < b + lCell )
			{
				neighborList_[i].push_back(j);
			}
		}
		
		Foam::Info<<"Size of "<< neighborList_[i].size()<<Foam::endl;
		Foam::Info<<cMesh_.findSphere(i, b)<<Foam::endl<<Foam::endl;
		for(auto j = 0; j <neighborList_[i].size(); j++)
		{
			Foam::Info<<neighborList_[i][j]<<Foam::endl;
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
				if( cellFaceDist < b + lCell)
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

	for(auto i=0; i<nCells; i++)
	{
		boundaryPatchNum_[i] = boundaryCell_[i].first;
	}

	listConstructed_ = true;
	return true;
}



pFlow::coupling::statistical::statistical(
	Foam::dictionary 		dict, 
	couplingMesh& 			cMesh, 
	Plus::centerMassField& 	centerMass, 
	Plus::realProcCMField& 	parDiam)
:
	porosity(dict, cMesh, centerMass, parDiam),
	neighborLength_(dict.lookup<Foam::scalar>("neighborLength")),
	neighborList_(cMesh.mesh().nCells()),
	boundaryCell_(cMesh.mesh().nCells(), std::pair<Foam::label, Foam::label>{-1,-1}),
	listConstructed_(false),
	boundaryPatchNum_
	(
		Foam::IOobject
	    (
	        "boundaryPatchNum",
	        cMesh.mesh().time().timeName(),
	        cMesh.mesh(),
	        Foam::IOobject::NO_READ,
	        Foam::IOobject::AUTO_WRITE
	    ),
   		cMesh.mesh(),
   		Foam::dimless
	)
{
}




