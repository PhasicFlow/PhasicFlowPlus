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




pFlow::coupling::statistical::statistical(
	Foam::dictionary 		dict, 
	couplingMesh& 			cMesh, 
	MPI::centerMassField& 	centerMass, 
	MPI::realProcCMField& 	parDiam)
:
	porosity(dict, cMesh, centerMass, parDiam),
	b_(dict.lookup<Foam::scalar>("b")),
	neighborList_(cMesh.mesh().nCells())
{
	
	const vectorField& cellC = cMesh.mesh().cellCentres();

	// loop over all cells
	for(size_t i=0; i<neighborList_.size(); i++)
	{
		Foam::vector ci = cellC[i];
		for(size_t j=0; j< neighborList_.size(); j++)
		{
			auto dist = mag(cellC[j]-ci);

			if(dist < 3*b_)
			{
				neighborList_[i].push_back(j);
			}
		}
	}

	
}


bool pFlow::coupling::statistical::internalFieldUpdate()
{
	
	auto solidVoldTmp = Foam::volScalarField::Internal::New(
		"solidVol",
		this->mesh(),
		 Foam::dimensioned("solidVol", Foam::dimVol, Foam::scalar(0))
		 	);
	
	auto& solidVol = solidVoldTmp.ref();
	
	size_t numPar = centerMass_.size();
	const vectorField& allCellCntr = this->mesh().cellCentres();
	
	
	//#pragma omp parallel for 
	for(size_t i=0; i<numPar; i++)
	{
		real parVolume = static_cast<real>(3.14159265358979/6)*
				pFlow::pow(particleDiameter_[i], static_cast<real>(3.0));

		const Foam::label targetCellId = parCellIndex_[i];
		
		if( targetCellId < 0 )continue;
		
		const auto& nList=neighborList_[targetCellId];
		realx3 cp = centerMass_[i];	

		std::vector<Foam::scalar> ks; 
		ks.reserve(nList.size());
		
		Foam::scalar pSubTotal = 0;
		
		for(auto j:nList)
		{
			auto c2 = allCellCntr[j];
			Foam::vector dx{cp.x()-c2.x(), cp.y()-c2.y(), cp.z()-c2.z()};

			Foam::scalar ksi2 = dot(dx,dx)/(b_*b_);
			
			Foam::scalar k = 0;
			if (ksi2<1)
			{
				k = Foam::pow(1-ksi2,4);
			}
			else
			{
				k = 0 ;
			}
			ks.push_back(k);
			pSubTotal += k;

		}

		size_t n = 0;
		for(auto j:nList)
		{
			solidVol[j] += ks[n++] /pSubTotal * parVolume;
		}

	} 
	
	this->ref() = Foam::max(
		1 - solidVol/this->mesh().V(), 
		static_cast<Foam::scalar>(this->alphaMin()) );

	return true;
}



